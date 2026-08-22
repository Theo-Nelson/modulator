#!/usr/bin/env python3
"""
Find SNPs that fall EXACTLY ON a modified base and report them separately, then
(optionally) FLAG those sites in the differential-modification result tables.

Why this exists
---------------
A segregating variant at the modified base is a confounder that the differential
tests (between-isoform `test_diffs`, between-condition `between_conditions`) do not
model: a base-changing SNP routes the alt-allele reads into modkit's `Ndiff`
bucket (not `Nvalid_cov`), so the modification fraction is computed on the
REFERENCE allele only, and a genotype difference can masquerade as a modification
difference. Rather than recalibrate the tests, this module (a) enumerates every
modified site that coincides with a candidate SNP, classifies the coincidence, and
(b) writes a boolean `snp_at_mod_base` flag into the differential tables so those
sites can be scrutinised or excluded downstream.

Classification at the modified base (distance 0)
------------------------------------------------
Because a SNP changes the base, at distance 0 the variant is never neutral:
  EDITING_SELF_REPORT   inosine (mod 17596), transcript-oriented A>G -- the edit
                        IS the "SNP"; the association is circular.
  PSEU_SELF_REPORT      pseudouridine (mod 17802), T>C -- U->C basecall error is
                        called as the variant at its own site; circular.
  MOD_BASE_CREATED      the ALT allele carries the modifiable base (ref does not);
                        the variant CREATES the substrate, and the modification is
                        read off the ALT reads. High-confidence cis class.
  MOD_BASE_ABLATED      the REF allele carries the base and the ALT removes it, so the
                        alt allele cannot carry the modification -- a genuine cis
                        genotype effect (fraction read off the REF reads).
  INCONSISTENT          neither reported allele is the modifiable base: a modification
                        called where no allele can carry it (a third allele under the
                        min-alt floor, or a strand/annotation mismatch) -- a data-
                        consistency flag, not a clean biological category.
  AT_BASE_UNKNOWN_MOD   the mod code has no known modifiable base in this table.

Inputs: candidate SNPs (`discover_candidate_snps.py`) + the ZN modification table
(`aggregate_zn <prefix>_FILTERED_sites_long.tsv`). Reference FASTA optional.
"""
from __future__ import annotations
import argparse
import csv
import os
import sys
from collections import defaultdict

# mod_code -> canonical (transcript-oriented) base it sits on
MOD_BASE = {"17596": "A", "a": "A", "m": "C", "17802": "T",
            "69426": "A", "19228": "C", "19229": "G", "19227": "T"}
# self-reporting signatures, transcript-oriented (mod_code -> (ref, alt, class))
SELF_REPORT = {"17596": ("A", "G", "EDITING_SELF_REPORT"),
               "17802": ("T", "C", "PSEU_SELF_REPORT")}
COMP = str.maketrans("ACGTN", "TGCAN")


def log(*a):
    print("[snp@mod]", *a, file=sys.stderr, flush=True)


def load_snps(path):
    """(chrom,pos0) -> dict(ref, alt, alt_frac, snp_id)."""
    out = {}
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            try:
                pos0 = int(row["pos1"]) - 1
            except (KeyError, ValueError):
                continue
            out[(row["chrom"], pos0)] = dict(
                ref=row.get("ref", ""), alt=row.get("alt", ""),
                alt_frac=row.get("alt_frac", ""), snp_id=row.get("snp_id", ""))
    return out


def load_mod_sites(path):
    """(chrom,start0,strand,mod_code) -> dict(gene_name, nmod, ncov, n_samples)."""
    agg = {}
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            try:
                key = (row["chrom"], int(row["start0"]), row["strand"], row["mod_code"])
                nmod = float(row.get("Nmod", 0) or 0)
                ncov = float(row.get("Nvalid_cov", 0) or 0)
            except (KeyError, ValueError):
                continue
            d = agg.setdefault(key, dict(gene_name=row.get("gene_name", ""),
                                         nmod=0.0, ncov=0.0, samples=set()))
            # the ZN long table is per (site x ZN-partition x sample); count DISTINCT samples, not rows
            # (else a single-sample run with 3 ZN partitions would report n_samples=3).
            d["nmod"] += nmod; d["ncov"] += ncov; d["samples"].add(str(row.get("sample", "")))
    return agg


def classify(mod_code, strand, ref, alt):
    tref = ref if strand == "+" else ref.translate(COMP)
    talt = alt if strand == "+" else alt.translate(COMP)
    sig = SELF_REPORT.get(str(mod_code))
    if sig and tref == sig[0] and talt == sig[1]:
        return sig[2], tref, talt
    base = MOD_BASE.get(str(mod_code))
    if base is None:
        return "AT_BASE_UNKNOWN_MOD", tref, talt     # no known modifiable base for this mod code
    if talt == base:
        # the ALT allele CARRIES the modifiable base (the variant creates it); modification is read
        # off the ALT reads, not the REF. Previously mislabelled "AT_BASE_OTHER" -- it is in fact the
        # highest-confidence cis class (the alt allele can carry the modification).
        return "MOD_BASE_CREATED", tref, talt
    if tref == base:
        # REF carries the base, ALT removes it -> the alt allele cannot carry the modification: a
        # genuine cis genotype effect (the fraction is read off the REF reads).
        return "MOD_BASE_ABLATED", tref, talt
    # neither reported allele is the modifiable base: a modification called where no allele can carry
    # it -> a data-consistency signal (e.g. a third allele under the min-alt floor, or a strand/annot
    # mismatch), NOT a clean ablation.
    return "INCONSISTENT", tref, talt


def detect(snps, mods):
    rows = []
    hits = {}  # (chrom,start0,mod_code) -> class  (for annotation)
    for (chrom, start0, strand, mod_code), d in mods.items():
        snp = snps.get((chrom, start0))
        if not snp:
            continue
        cls, tref, talt = classify(mod_code, strand, snp["ref"], snp["alt"])
        frac = (d["nmod"] / d["ncov"]) if d["ncov"] else 0.0
        rows.append(dict(
            chrom=chrom, pos0=start0, pos1=start0 + 1, strand=strand, mod_code=mod_code,
            canonical_base=MOD_BASE.get(str(mod_code), ""), ref=snp["ref"], alt=snp["alt"],
            tx_ref=tref, tx_alt=talt, alt_frac=snp["alt_frac"],
            snp_at_mod_base_class=cls, gene_name=d["gene_name"], n_samples=len(d["samples"]),
            total_Nmod=int(d["nmod"]), total_Nvalid_cov=int(d["ncov"]),
            pooled_frac_mod=round(frac, 4), snp_id=snp["snp_id"]))
        # key by STRAND too: two modifications at one genomic position on opposite strands (antisense
        # gene overlap) otherwise collide in this index and one classification is silently overwritten.
        hits[(chrom, start0, strand, str(mod_code))] = cls
    rows.sort(key=lambda r: (r["chrom"], r["pos0"], r["mod_code"]))
    return rows, hits


def annotate_tsv(path, hits):
    """Add `snp_at_mod_base` (0/1) + `snp_at_mod_base_class` columns in place
    (idempotent), joining on (chrom, mod-start0, mod-code).

    The mod-site position column is `start0` (mod-diff tables) or `mod_start0`
    (snp_mod_assoc); the mod-code column is `mod_code` or `target_mod_code`. This lets
    the same flag be surfaced onto snp_mod_assoc, whose rows would otherwise be silently
    skipped for lacking a literal `start0`/`mod_code` column."""
    try:
        with open(path) as fh:
            rr = csv.DictReader(fh, delimiter="\t")
            fields = list(rr.fieldnames or [])
            pos_col = "start0" if "start0" in fields else ("mod_start0" if "mod_start0" in fields else None)
            mc_col = "mod_code" if "mod_code" in fields else ("target_mod_code" if "target_mod_code" in fields else None)
            if "chrom" not in fields or pos_col is None:
                return False
            rows = list(rr)
    except (OSError, csv.Error):
        return False
    for c in ("snp_at_mod_base", "snp_at_mod_base_class"):
        if c not in fields:
            fields.append(c)
    for row in rows:
        ch = row.get("chrom"); s0 = _to_int(row.get(pos_col))
        st = row.get("strand", ""); mc = str(row.get(mc_col, "")) if mc_col else ""
        cls = hits.get((ch, s0, st, mc))
        if cls is None:
            # progressively looser match for differential tables missing strand and/or mod_code
            for (hch, hs0, hst, hmc), v in hits.items():
                if hch == ch and hs0 == s0 and (not mc or hmc == mc) and (not st or hst == st):
                    cls = v
                    break
        row["snp_at_mod_base"] = "1" if cls else "0"
        row["snp_at_mod_base_class"] = cls or ""
    # atomic write (.tmp + os.replace) so a crash/full-disk mid-write can't truncate a table that
    # took hours to produce -- matching the convention in discover_candidate_snps / build_molecule_snp_table.
    tmp = path + ".tmp"
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    os.replace(tmp, path)
    return True


def _to_int(x):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidate-snps", required=True)
    ap.add_argument("--mod-sites", required=True, help="aggregate_zn <prefix>_FILTERED_sites_long.tsv")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--annotate", nargs="*", default=[],
                    help="differential-result TSVs to add a snp_at_mod_base flag column to (in place)")
    ap.add_argument("--verbose", action="store_true")
    a = ap.parse_args()

    snps = load_snps(a.candidate_snps)
    mods = load_mod_sites(a.mod_sites)
    rows, hits = detect(snps, mods)

    cols = ["chrom", "pos0", "pos1", "strand", "mod_code", "canonical_base", "ref", "alt",
            "tx_ref", "tx_alt", "alt_frac", "snp_at_mod_base_class", "gene_name", "n_samples",
            "total_Nmod", "total_Nvalid_cov", "pooled_frac_mod", "snp_id"]
    with open(a.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    if a.verbose:
        by_cls = defaultdict(int)
        for r in rows:
            by_cls[r["snp_at_mod_base_class"]] += 1
        log(f"{len(rows)} modified site(s) with a SNP at the modified base "
            f"({len(snps)} candidate SNPs, {len(mods)} mod sites) -> {a.out_tsv}")
        for k, v in sorted(by_cls.items()):
            log(f"    {k}: {v}")

    for t in a.annotate:
        if annotate_tsv(t, hits) and a.verbose:
            log(f"flagged snp_at_mod_base in {t}")


if __name__ == "__main__":
    main()
