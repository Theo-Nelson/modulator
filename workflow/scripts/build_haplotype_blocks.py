#!/usr/bin/env python3

import argparse
from collections import Counter, defaultdict, deque
import itertools
import os

import pandas as pd


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


def context_key(row):
    mg = [x for x in str(row.get("metagene_indices", "")).split(";") if x]
    if len(set(mg)) == 1:
        return f"MG:{mg[0]}"
    genes = [x for x in str(row.get("gene_names", "")).split(";") if x]
    if len(set(genes)) == 1:
        return f"GENE:{genes[0]}"
    return f"CHR:{row['chrom']}"


def split_component(snps_sorted, max_block_snps):
    if len(snps_sorted) <= max_block_snps:
        return [snps_sorted]
    return [snps_sorted[i:i + max_block_snps] for i in range(0, len(snps_sorted), max_block_snps)]


def main():
    args = parse_args()
    df = pd.read_csv(args.molecule_snps, sep="\t", low_memory=False)
    if df.empty:
        pd.DataFrame(columns=["block_id", "context_key", "chrom", "n_snps", "snp_ids", "support_reads", "complete_reads", "haplotypes"]).to_csv(
            args.out_blocks_tsv, sep="\t", index=False
        )
        pd.DataFrame(columns=["sample", "qname", "block_id", "context_key", "chrom", "haplotype", "support_rank", "ZT", "ZG", "ZN", "ZM"]).to_csv(
            args.out_molecules_tsv, sep="\t", index=False
        )
        return

    df = df[df["allele_class"].isin(["ref", "alt"])].copy()
    if df.empty:
        pd.DataFrame(columns=["block_id", "context_key", "chrom", "n_snps", "snp_ids", "support_reads", "complete_reads", "haplotypes"]).to_csv(
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
        pd.DataFrame(columns=["block_id", "context_key", "chrom", "n_snps", "snp_ids", "support_reads", "complete_reads", "haplotypes"]).to_csv(
            args.out_blocks_tsv, sep="\t", index=False
        )
        pd.DataFrame(columns=["sample", "qname", "block_id", "context_key", "chrom", "haplotype", "support_rank", "ZT", "ZG", "ZN", "ZM"]).to_csv(
            args.out_molecules_tsv, sep="\t", index=False
        )
        return

    df["context_key"] = df.apply(context_key, axis=1)
    snp_meta = (
        df[["snp_id", "chrom", "pos1", "ref", "alt", "context_key"]]
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
                block_rows.append({
                    "block_id": block_id,
                    "context_key": ctx,
                    "chrom": chrom,
                    "n_snps": len(chunk),
                    "snp_ids": ";".join(chunk),
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
        block_df = pd.DataFrame(columns=["block_id", "context_key", "chrom", "n_snps", "snp_ids", "support_reads", "complete_reads", "haplotypes"])
    mol_df = pd.DataFrame(molecule_rows)
    if mol_df.empty:
        mol_df = pd.DataFrame(columns=["sample", "qname", "block_id", "context_key", "chrom", "haplotype", "support_rank", "ZT", "ZG", "ZN", "ZM"])
    block_df.to_csv(args.out_blocks_tsv, sep="\t", index=False)
    mol_df.to_csv(args.out_molecules_tsv, sep="\t", index=False)


def safe_int(x, default=0):
    try:
        return int(x)
    except Exception:
        return default


if __name__ == "__main__":
    main()
