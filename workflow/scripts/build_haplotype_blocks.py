#!/usr/bin/env python3

import argparse
from collections import Counter, defaultdict, deque
import itertools
import os

import numpy as np
import pandas as pd
from pyroaring import BitMap

from genotype_utils import open_chrom_table, safe_int, snp_context_keys_lists

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

def split_component(snps_sorted, max_block_snps, read_snps=None, min_reads=1, snp_reads=None):
    """Split a phased component into PHASEABLE blocks (>= min_reads reads co-cover every SNP in the
    block), each <= max_block_snps SNPs.

    The old blind fixed-index windowing had three failure modes: a chunk whose SNPs no single read
    fully covered was silently discarded (so a WIDER --max-block-snps LOST blocks); a size-1 remainder
    was dropped; and --max-block-snps 1 produced ZERO blocks (every chunk failed the >=2 rule). This
    read-aware greedy grows a block from the sorted SNPs while it stays within the cap AND enough reads
    still co-cover the whole block, closing it (and restarting) otherwise -- so every emitted block is
    actually phaseable and SNPs are grouped by real co-coverage rather than array index."""
    cap = max(2, int(max_block_snps))   # a haplotype block needs >= 2 SNPs; cap 1 would emit nothing
    if read_snps is None:
        # coverage-agnostic fallback (preserves old behaviour when no read map is supplied)
        return [snps_sorted[i:i + cap] for i in range(0, len(snps_sorted), cap)]

    if snp_reads is not None:
        # per-SNP read bitmaps: co-coverage of a chunk = size of the intersection (was a scan of EVERY
        # read of the context for every SNP added to a block -> O(SNPs x reads) per component)
        def n_complete(chunk):
            bm = snp_reads[chunk[0]]
            for s_ in chunk[1:]:
                bm = bm & snp_reads[s_]
            return len(bm)
    else:
        def n_complete(chunk):
            return sum(1 for v in read_snps.values() if all(s in v for s in chunk))

    blocks, cur = [], []
    for s in snps_sorted:
        trial = cur + [s]
        if len(trial) <= cap and (len(trial) < 2 or n_complete(trial) >= min_reads):
            cur = trial
        else:
            if len(cur) >= 2:
                blocks.append(cur)
            cur = [s]
    if len(cur) >= 2:
        blocks.append(cur)
    elif len(cur) == 1 and blocks:
        # DISJOINT rebalance: the greedy commits maximally leftward, so a SNP that only pairs rightward
        # (or a trailing n%cap==1 remainder) is stranded as a size-1 block and would be dropped -- and
        # a WIDER --max-block-snps then LOSES SNPs. Move the previous block's LAST SNP into a new
        # [prev[-1], orphan] block so no SNP is lost and blocks stay NON-overlapping (keeping the
        # haplotype-mod tests independent). Requires prev to remain >= 2 after the move and both blocks
        # to stay phaseable; otherwise the orphan cannot be recovered disjointly at this cap and is
        # dropped (unavoidable when the cap is too small to tile the component without overlap).
        prev = blocks[-1]
        newblk = [prev[-1], cur[0]]
        if len(prev) >= 3 and n_complete(newblk) >= min_reads:
            blocks[-1] = prev[:-1]
            blocks.append(newblk)
    return blocks


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


BLOCK_COLS = ["block_id", "context_key", "gene_names", "chrom", "region", "start1", "end1", "span_bp",
              "n_snps", "snp_ids", "snp_coords", "support_reads", "complete_reads", "haplotypes"]
MOL_COLS = ["sample", "qname", "block_id", "context_key", "chrom", "haplotype", "support_rank", "ZT", "ZG", "ZN", "ZM"]


def _write_empty(args):
    pd.DataFrame(columns=BLOCK_COLS).to_csv(args.out_blocks_tsv, sep="\t", index=False)
    pd.DataFrame(columns=MOL_COLS).to_csv(args.out_molecules_tsv, sep="\t", index=False)


def _blocks_for_chrom(df, args):
    """Haplotype blocks for ONE chromosome's molecule-SNP rows. Returns {context_key: [(block_row_without_id,
    [molecule_rows_without_id]), ...]} in the per-context block order. Block ids are assigned by the caller
    over the GLOBAL sorted context order, so the numbering is identical to the old whole-table pass."""
    out = {}
    df = df[df["allele_class"].isin(["ref", "alt"])].copy()
    if df.empty:
        return out
    alt_support = df.loc[df["allele_class"] == "alt"].groupby("snp_id").size()
    keep_snps = set(alt_support[alt_support >= int(args.min_alt_reads)].index)
    df = df[df["snp_id"].isin(keep_snps)].copy()
    if df.empty:
        return out

    # Fan a SNP out to ALL its metagene (MG:) contexts, then explode -- matching the snp_mod_assoc fix
    # (715916d). The old singular context_key_from_snp_row collapsed a multi-metagene SNP to a single
    # CHR:-keyed context, so those blocks never received a haplotype x mod test.
    df["context_key"] = snp_context_keys_lists(df)
    df = df.explode("context_key", ignore_index=True)
    df = df[df["context_key"].astype(str) != ""]
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

    for ctx, sub in df.groupby("context_key", sort=False):
        # first SNP id of this (context, chrom): orders chromosome groups of a context that spans
        # several contigs exactly as the old whole-table pass did (components sorted by snp_id string)
        _ctx_first_snp = str(sub["snp_id"].astype(str).min())
        read_snps = defaultdict(dict)
        # One-time (sample, qname) -> first row lookup for this context (column iteration, not
        # to_dict("records") -- 73 GiB on an 8-sample run). Output is identical.
        first_by_read = {}
        _zt = sub["ZT"].tolist() if "ZT" in sub.columns else [""] * len(sub)
        _zg = sub["ZG"].tolist() if "ZG" in sub.columns else [""] * len(sub)
        _zn = sub["ZN"].tolist() if "ZN" in sub.columns else [""] * len(sub)
        _zm = sub["ZM"].tolist() if "ZM" in sub.columns else [""] * len(sub)
        for sm, qn, zt, zg, zn, zm in zip(sub["sample"].tolist(), sub["qname"].tolist(), _zt, _zg, _zn, _zm):
            first_by_read.setdefault((sm, qn), {"ZT": zt, "ZG": zg, "ZN": zn, "ZM": zm})
        for sm, qn, sid, ob in zip(sub["sample"].tolist(), sub["qname"].tolist(),
                                   sub["snp_id"].tolist(), sub["observed_base"].tolist()):
            read_snps[(sm, qn)][sid] = ob
        # per-SNP bitmaps over the context's reads (read index = first-appearance order)
        snp_reads = defaultdict(BitMap)
        for ri, (_rk, snp_map) in enumerate(read_snps.items()):
            for sid in snp_map:
                snp_reads[sid].add(ri)

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
            components = [[s_] for s_ in singleton_snps if len(read_snps) >= int(args.min_cocover_reads)]

        ctx_blocks = []
        for comp in components:
            comp = sorted(comp, key=lambda x: (snp_meta[x]["chrom"], safe_int(snp_meta[x]["pos1"])))
            for chunk in split_component(comp, int(args.max_block_snps), read_snps=read_snps,
                                         min_reads=int(args.min_cocover_reads), snp_reads=snp_reads):
                if len(chunk) < 2:
                    continue
                chrom = snp_meta[chunk[0]]["chrom"]
                hap_counter = Counter()
                hap_members = []
                for (sample, qname), snp_map in read_snps.items():
                    if not all(s_ in snp_map for s_ in chunk):
                        continue
                    haplotype = "|".join(snp_map[s_] for s_ in chunk)
                    hap_counter[haplotype] += 1
                    hap_members.append((sample, qname, haplotype))
                if not hap_counter:
                    continue
                keep_haps = {h for h, n in hap_counter.items() if n >= int(args.min_haplotype_reads)}
                complete_reads = sum(hap_counter.values())
                ctxinfo = block_context(chunk, snp_meta)
                support = snp_reads[chunk[0]]
                for s_ in chunk[1:]:
                    support = support | snp_reads[s_]
                block_row = {
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
                    "support_reads": len(support),
                    "complete_reads": complete_reads,
                    "haplotypes": ";".join(f"{h}:{hap_counter[h]}" for h in sorted(hap_counter, key=lambda x: (-hap_counter[x], x))),
                }
                rank = {hap: i + 1 for i, (hap, _) in enumerate(sorted(hap_counter.items(), key=lambda x: (-x[1], x[0])))}
                mol_rows = []
                for sample, qname, hap in hap_members:
                    if hap not in keep_haps:
                        hap = "OTHER"
                    first = first_by_read.get((sample, qname), {})
                    mol_rows.append({
                        "sample": sample,
                        "qname": qname,
                        "context_key": ctx,
                        "chrom": chrom,
                        "haplotype": hap,
                        "support_rank": rank.get(hap, 999),
                        "ZT": first.get("ZT", ""),
                        "ZG": first.get("ZG", ""),
                        "ZN": first.get("ZN", ""),
                        "ZM": first.get("ZM", ""),
                    })
                ctx_blocks.append((block_row, mol_rows))
        if ctx_blocks:
            out[ctx] = (_ctx_first_snp, ctx_blocks)
    return out


def main():
    args = parse_args()
    tbl = open_chrom_table(args.molecule_snps)
    usecols = [c for c in WANTED_COLS if c in tbl.header_cols]
    dtype = {c: t for c, t in CATEGORICAL.items() if c in usecols}
    if not tbl.chroms:
        _write_empty(args)
        return

    # One chromosome at a time (a byte range of the chrom-sorted table; 73 GiB when loaded whole on an
    # 8-sample mosquito run). A context_key embeds its chromosome/metagene, so no context spans two
    # chromosomes and the per-chrom blocks are exactly the old per-context blocks; only the global
    # HAPBLOCK numbering is assigned afterwards, over the same sorted context order the old code used.
    # A metagene can span a primary chromosome AND its alt contig (the assembler keys genes by
    # name/id, not contig), so a context_key may appear in more than one chromosome: APPEND each
    # chromosome's block list and order the groups by their first snp_id string, which is the order
    # the old single pass (components sorted by snp_id) produced.
    by_ctx = {}
    n_rows = 0
    try:
        for chrom in tbl.chroms:
            df = tbl.read(chrom, usecols=usecols, dtype=dtype)
            n_rows += len(df)
            if df.empty:
                continue
            for ctx, (first_snp, blocks) in _blocks_for_chrom(df, args).items():
                by_ctx.setdefault(ctx, []).append((first_snp, blocks))
            del df
    finally:
        if hasattr(tbl, "close"):
            tbl.close()
    if n_rows == 0 or not by_ctx:
        _write_empty(args)
        return

    block_rows = []
    molecule_rows = []
    block_idx = 0
    for ctx in sorted(by_ctx):
        for _first_snp, blocks in sorted(by_ctx[ctx], key=lambda t: t[0]):
          for block_row, mol_rows in blocks:
            block_idx += 1
            block_id = f"HAPBLOCK{block_idx}"
            block_rows.append({"block_id": block_id, **block_row})
            for m in mol_rows:
                molecule_rows.append({"sample": m["sample"], "qname": m["qname"], "block_id": block_id,
                                      **{k: v for k, v in m.items() if k not in ("sample", "qname")}})

    os.makedirs(os.path.dirname(args.out_blocks_tsv) or ".", exist_ok=True)
    block_df = pd.DataFrame(block_rows)
    if block_df.empty:
        block_df = pd.DataFrame(columns=BLOCK_COLS)
    mol_df = pd.DataFrame(molecule_rows)
    if mol_df.empty:
        mol_df = pd.DataFrame(columns=MOL_COLS)
    block_df = block_df[BLOCK_COLS]
    mol_df = mol_df[MOL_COLS]
    block_df.to_csv(args.out_blocks_tsv, sep="\t", index=False)
    mol_df.to_csv(args.out_molecules_tsv, sep="\t", index=False)

if __name__ == "__main__":
    main()
