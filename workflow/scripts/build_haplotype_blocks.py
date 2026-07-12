#!/usr/bin/env python3

import argparse
from collections import Counter, defaultdict, deque
import itertools
import os

import pandas as pd

from genotype_utils import context_key_from_snp_row, safe_int, tsv_header

# Columns actually consumed (grouping/keys + per-SNP metadata + the ZT/ZG/ZN/ZM tags copied into
# the molecules output). Loading only these -- with repeated string columns as categoricals --
# avoids materializing all 21 object columns of the ~1.7 GB / 7.5M-row molecule_snps table.
WANTED_COLS = [
    "sample", "qname", "snp_id", "chrom", "pos1", "ref", "alt", "allele_class",
    "observed_base", "gene_names", "gene_ids", "metagene_indices", "ZT", "ZG", "ZN", "ZM",
]
# Categoricals ONLY for columns that are never a sort key and never grouped/pivoted with a
# behavior that depends on category order. chrom is EXCLUDED: it is part of the sort_values key,
# and a categorical sorts by category codes (appearance order) rather than lexicographically,
# which would change deterministic block numbering. allele_class stays object (it is .isin'd and
# grouped). snp_id stays object (it is grouped and its .index becomes keep_snps).
CATEGORICAL = {
    "ref": "category", "alt": "category",
    "gene_names": "category", "gene_ids": "category", "metagene_indices": "category",
    "observed_base": "category", "ZT": "category", "ZG": "category", "ZN": "category", "ZM": "category",
}


def parse_args():
    ap = argparse.ArgumentParser(description="Build local read-backed haplotype blocks from candidate SNP molecules.")
    ap.add_argument("--molecule-snps", required=True, help="Molecule SNP TSV")
    ap.add_argument("--out-blocks-tsv", required=True, help="Output haplotype block TSV")
    ap.add_argument("--out-molecules-tsv", required=True, help="Output molecule haplotype TSV")
    ap.add_argument("--min-alt-reads", type=int, default=4)
    ap.add_argument("--min-cocover-reads", type=int, default=4)
    ap.add_argument("--max-block-snps", type=int, default=4)
    ap.add_argument("--min-haplotype-reads", type=int, default=4)
    return ap.parse_args()

def split_component(snps_sorted, max_block_snps):
    if len(snps_sorted) <= max_block_snps:
        return [snps_sorted]
    return [snps_sorted[i:i + max_block_snps] for i in range(0, len(snps_sorted), max_block_snps)]


def block_context(chunk, snp_meta):
    """Gene + coordinate context for a haplotype block, derived from its member SNPs.
    Returns gene_names (unique, ';'-joined), region (chrom:start-end, 1-based),
    start1/end1/span_bp, and a readable per-SNP coordinate string (chrom:pos ref>alt)."""
    chrom = str(snp_meta[chunk[0]].get("chrom", ""))
    positions = [safe_int(snp_meta[s].get("pos1")) for s in chunk]
    start1 = min(positions) if positions else 0
    end1 = max(positions) if positions else 0
    genes = []
    for s in chunk:
        for tok in str(snp_meta[s].get("gene_names", "") or "").split(";"):
            tok = tok.strip()
            if tok and tok.lower() not in {"nan", "none", "null"} and tok not in genes:
                genes.append(tok)
    snp_coords = "; ".join(
        f"{snp_meta[s].get('chrom', '')}:{snp_meta[s].get('pos1', '')} "
        f"{snp_meta[s].get('ref', '')}>{snp_meta[s].get('alt', '')}"
        for s in chunk
    )
    return {
        "gene_names": ";".join(genes),
        "region": f"{chrom}:{start1}-{end1}" if chrom else "",
        "start1": start1,
        "end1": end1,
        "span_bp": end1 - start1,
        "snp_coords": snp_coords,
    }


def main():
    args = parse_args()
    header = tsv_header(args.molecule_snps)
    usecols = [c for c in WANTED_COLS if c in header]
    dtype = {c: t for c, t in CATEGORICAL.items() if c in usecols}
    df = pd.read_csv(args.molecule_snps, sep="\t", usecols=usecols, dtype=dtype, low_memory=False)
    if df.empty:
        pd.DataFrame(columns=["block_id", "context_key", "gene_names", "chrom", "region", "start1", "end1", "span_bp", "n_snps", "snp_ids", "snp_coords", "support_reads", "complete_reads", "haplotypes"]).to_csv(
            args.out_blocks_tsv, sep="\t", index=False
        )
        pd.DataFrame(columns=["sample", "qname", "block_id", "context_key", "chrom", "haplotype", "support_rank", "ZT", "ZG", "ZN", "ZM"]).to_csv(
            args.out_molecules_tsv, sep="\t", index=False
        )
        return

    df = df[df["allele_class"].isin(["ref", "alt"])].copy()
    if df.empty:
        pd.DataFrame(columns=["block_id", "context_key", "gene_names", "chrom", "region", "start1", "end1", "span_bp", "n_snps", "snp_ids", "snp_coords", "support_reads", "complete_reads", "haplotypes"]).to_csv(
            args.out_blocks_tsv, sep="\t", index=False
        )
        pd.DataFrame(columns=["sample", "qname", "block_id", "context_key", "chrom", "haplotype", "support_rank", "ZT", "ZG", "ZN", "ZM"]).to_csv(
            args.out_molecules_tsv, sep="\t", index=False
        )
        return

    alt_support = df.loc[df["allele_class"] == "alt"].groupby("snp_id").size()
    keep_snps = set(alt_support[alt_support >= int(args.min_alt_reads)].index)
    df = df[df["snp_id"].isin(keep_snps)].copy()
    if df.empty:
        pd.DataFrame(columns=["block_id", "context_key", "gene_names", "chrom", "region", "start1", "end1", "span_bp", "n_snps", "snp_ids", "snp_coords", "support_reads", "complete_reads", "haplotypes"]).to_csv(
            args.out_blocks_tsv, sep="\t", index=False
        )
        pd.DataFrame(columns=["sample", "qname", "block_id", "context_key", "chrom", "haplotype", "support_rank", "ZT", "ZG", "ZN", "ZM"]).to_csv(
            args.out_molecules_tsv, sep="\t", index=False
        )
        return

    df["context_key"] = df.apply(context_key_from_snp_row, axis=1)
    # Order-invariance: block numbering (HAPBLOCK<i>) follows context_key first-appearance
    # and read iteration follows row order, so a deterministic sort makes the haplotype
    # blocks independent of upstream (BAM x chrom) shard completion order.
    df = df.sort_values(["context_key", "chrom", "pos1", "sample", "qname"], kind="stable").reset_index(drop=True)
    meta_cols = ["snp_id", "chrom", "pos1", "ref", "alt", "context_key"]
    for extra in ("gene_names", "gene_ids"):
        if extra in df.columns:
            meta_cols.append(extra)
    snp_meta = (
        df[meta_cols]
        .drop_duplicates("snp_id")
        .set_index("snp_id")
        .to_dict("index")
    )

    block_rows = []
    molecule_rows = []
    block_idx = 0

    for ctx, sub in df.groupby("context_key", sort=False):
        read_snps = defaultdict(dict)
        for row in sub.itertuples(index=False):
            read_snps[(row.sample, row.qname)][row.snp_id] = row.observed_base

        edge_counts = Counter()
        for snp_map in read_snps.values():
            snps = sorted(snp_map)
            for a, b in itertools.combinations(snps, 2):
                edge_counts[(a, b)] += 1

        adjacency = defaultdict(set)
        for (a, b), count in edge_counts.items():
            if count >= int(args.min_cocover_reads):
                adjacency[a].add(b)
                adjacency[b].add(a)

        seen = set()
        components = []
        component_ids = set(adjacency.keys())
        for vals in adjacency.values():
            component_ids.update(vals)
        for snp_id in sorted(component_ids):
            if snp_id in seen:
                continue
            queue = deque([snp_id])
            comp = []
            seen.add(snp_id)
            while queue:
                cur = queue.popleft()
                comp.append(cur)
                for nxt in adjacency.get(cur, []):
                    if nxt not in seen:
                        seen.add(nxt)
                        queue.append(nxt)
            components.append(comp)

        if not components:
            singleton_snps = sorted(sub["snp_id"].unique(), key=lambda x: (snp_meta[x]["chrom"], safe_int(snp_meta[x]["pos1"])))
            components = [[s] for s in singleton_snps if len(read_snps) >= int(args.min_cocover_reads)]

        for comp in components:
            comp = sorted(comp, key=lambda x: (snp_meta[x]["chrom"], safe_int(snp_meta[x]["pos1"])))
            for chunk in split_component(comp, int(args.max_block_snps)):
                if len(chunk) < 2:
                    continue
                block_idx += 1
                block_id = f"HAPBLOCK{block_idx}"
                chrom = snp_meta[chunk[0]]["chrom"]
                hap_counter = Counter()
                hap_members = []
                for (sample, qname), snp_map in read_snps.items():
                    if not all(s in snp_map for s in chunk):
                        continue
                    haplotype = "|".join(snp_map[s] for s in chunk)
                    hap_counter[haplotype] += 1
                    hap_members.append((sample, qname, haplotype))
                if not hap_counter:
                    continue
                keep_haps = {h for h, n in hap_counter.items() if n >= int(args.min_haplotype_reads)}
                complete_reads = sum(hap_counter.values())
                ctxinfo = block_context(chunk, snp_meta)
                block_rows.append({
                    "block_id": block_id,
                    "context_key": ctx,
                    "gene_names": ctxinfo["gene_names"],
                    "chrom": chrom,
                    "region": ctxinfo["region"],
                    "start1": ctxinfo["start1"],
                    "end1": ctxinfo["end1"],
                    "span_bp": ctxinfo["span_bp"],
                    "n_snps": len(chunk),
                    "snp_ids": ";".join(chunk),
                    "snp_coords": ctxinfo["snp_coords"],
                    "support_reads": sum(1 for v in read_snps.values() if any(s in v for s in chunk)),
                    "complete_reads": complete_reads,
                    "haplotypes": ";".join(f"{h}:{hap_counter[h]}" for h in sorted(hap_counter, key=lambda x: (-hap_counter[x], x))),
                })
                rank = {hap: i + 1 for i, (hap, _) in enumerate(sorted(hap_counter.items(), key=lambda x: (-x[1], x[0])))}
                for sample, qname, hap in hap_members:
                    if hap not in keep_haps:
                        hap = "OTHER"
                    first = sub[(sub["sample"] == sample) & (sub["qname"] == qname)].iloc[0]
                    molecule_rows.append({
                        "sample": sample,
                        "qname": qname,
                        "block_id": block_id,
                        "context_key": ctx,
                        "chrom": chrom,
                        "haplotype": hap,
                        "support_rank": rank.get(hap, 999),
                        "ZT": first.get("ZT", ""),
                        "ZG": first.get("ZG", ""),
                        "ZN": first.get("ZN", ""),
                        "ZM": first.get("ZM", ""),
                    })

    os.makedirs(os.path.dirname(args.out_blocks_tsv) or ".", exist_ok=True)
    block_df = pd.DataFrame(block_rows)
    if block_df.empty:
        block_df = pd.DataFrame(columns=["block_id", "context_key", "gene_names", "chrom", "region", "start1", "end1", "span_bp", "n_snps", "snp_ids", "snp_coords", "support_reads", "complete_reads", "haplotypes"])
    mol_df = pd.DataFrame(molecule_rows)
    if mol_df.empty:
        mol_df = pd.DataFrame(columns=["sample", "qname", "block_id", "context_key", "chrom", "haplotype", "support_rank", "ZT", "ZG", "ZN", "ZM"])
    block_df.to_csv(args.out_blocks_tsv, sep="\t", index=False)
    mol_df.to_csv(args.out_molecules_tsv, sep="\t", index=False)

if __name__ == "__main__":
    main()
