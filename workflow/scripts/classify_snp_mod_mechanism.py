#!/usr/bin/env python3
"""Why does this SNP change this modification? -- mechanistic classification of SNP x mod pairs.

test_snp_mod_assoc.py says a SNP's alleles carry different modification rates. This says WHY, on
three orthogonal axes, for every (SNP, mod-site) pair:

  A. POSITIONAL  -- where the SNP sits relative to the modified base:
       AT_MOD_BASE (d=0) > IN_MOTIF_CORE (d<=2, the 5-mer) > IN_MOTIF_EXTENDED (d<=4, the 9-mer)
       > PROXIMAL_CIS (d<=--proximal-bp) > DISTAL_CIS
  B. MOTIF EFFECT -- m6A only (DRACH is the only mod here with a real 5-mer consensus; pseU /
       inosine / 5mC are structure- or neighbor-driven, so they get k-mer context but no verdict):
       does the alt allele DISRUPT / CREATE / PRESERVE the DRACH consensus (or is it absent in both,
       which makes the m6A call itself suspect)?
  C. CONCORDANCE -- does the observed allelic direction match what the motif predicts? A disrupted
       DRACH should show LESS modification on alt. Agreement promotes the pair to a coherent causal
       cis variant; disagreement points at trans/structure effects or a mis-call.

Special case worth its own flag: if the SNP IS the modified base and the alt allele cannot carry the
modification (m6A needs an A), the "association" is DEFINITIONAL, not regulatory -- there is simply
no A to methylate on alt reads. Those are flagged ``mod_base_ablated`` so they can be excluded.

Primary key = ``{positional_class}__{motif_effect}``.
"""
import argparse
import json
import os
import re
import sys

import pandas as pd
import pysam

DRACH_RE = re.compile(r"^[AGT][AG]AC[ACT]$")   # D-R-A*-C-H, methylated A at index 2
M6A_CODE = "a"
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# SELF-REPORTING VARIANTS -- the "SNP" IS the modification, not a genomic variant.
# SNPs are discovered from RNA reads, so a modification that changes the basecall gets called as a
# variant at its own position, and the resulting SNP x mod association is perfectly circular
# (reads carrying "alt" ARE the modified reads; alt_mod_rate is 0 by construction).
#   * A-to-I editing: inosine is basecalled as G  -> transcript-orientation A>G at the inosine site
#   * pseudouridine : characteristic U-to-C error -> transcript-orientation T>C at the pseU site
# Keyed by mod code -> (transcript-orientation ref, alt, flag).
SELF_REPORT_SIGNATURE = {
    "17596": ("A", "G", "EDITING_SELF_REPORT"),
    "17802": ("T", "C", "PSEU_SELF_REPORT"),
}

OUT_COLS = ["snp_id", "mod_site_id", "chrom", "snp_pos1", "mod_pos1", "strand", "target_mod_code",
            "gene_names", "distance_bp", "positional_class", "snp_ref_tx", "snp_alt_tx",
            "ref_5mer", "alt_5mer", "ref_9mer", "alt_9mer",
            "motif_effect", "artifact_flag", "mod_base_ablated", "ref_mod_rate", "alt_mod_rate",
            "observed_direction", "predicted_direction", "direction_concordance", "class_key",
            "n_reads", "effect_abs_delta_mod_frac", "p_adj_bh"]


def revcomp(s):
    return s.translate(_COMP)[::-1]


def parse_args():
    ap = argparse.ArgumentParser(description="Mechanistic (positional + motif) classification of SNP x mod pairs.")
    ap.add_argument("--snp-mod-assoc", required=True, help="*_snp_mod_assoc.tsv")
    ap.add_argument("--reference-fa", required=True, help="Reference FASTA (indexed)")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--proximal-bp", type=int, default=50, help="max distance for PROXIMAL_CIS before DISTAL_CIS")
    ap.add_argument("--verbose", action="store_true")
    return ap.parse_args()


def _kmer(fa, chrom, mod_pos0, strand, half):
    """(2*half+1)-mer centred on the modified base, in TRANSCRIPT orientation.

    For an odd window centred on the base, revcomp keeps the modified base at the centre index.
    """
    try:
        s = fa.fetch(chrom, mod_pos0 - half, mod_pos0 + half + 1).upper()
    except Exception:
        return None
    if len(s) != 2 * half + 1:
        return None
    return s if strand == "+" else revcomp(s)


def _sub(kmer, half, genomic_offset, alt_base, strand):
    """Substitute the SNP's alt allele into the k-mer. Returns None if the SNP falls outside it."""
    idx = half + genomic_offset if strand == "+" else half - genomic_offset
    if idx < 0 or idx >= len(kmer):
        return None
    alt_t = alt_base if strand == "+" else alt_base.translate(_COMP)
    return kmer[:idx] + alt_t + kmer[idx + 1:]


def main():
    args = parse_args()
    try:
        df = pd.read_csv(args.snp_mod_assoc, sep="\t", low_memory=False)
    except Exception:
        df = pd.DataFrame()
    if df.empty or "snp_id" not in df.columns:
        pd.DataFrame(columns=OUT_COLS).to_csv(args.out_tsv, sep="\t", index=False)
        return

    fa = pysam.FastaFile(args.reference_fa)
    rows = []
    for r in df.itertuples(index=False):
        # snp_id = chrom:pos1:REF>ALT ; mod_site_id = chrom:start0-end0:strand:code
        try:
            # parse from the RIGHT: contig names can contain ':' (HLA / ALT / decoy contigs), which
            # a plain split(':') would over-split, silently dropping the pair via the except below.
            s_chrom, s_pos, s_alleles = str(r.snp_id).rsplit(":", 2)
            s_ref, s_alt = s_alleles.split(">")
            snp_pos1 = int(s_pos)
            m_parts = str(r.mod_site_id).rsplit(":", 3)   # -> [chrom, start0-end0, strand, code]
            mod_start0 = int(m_parts[1].split("-")[0])
            strand = m_parts[2]
            code = m_parts[3]
        except Exception:
            continue
        if strand not in ("+", "-"):
            continue
        mod_pos0 = mod_start0
        mod_pos1 = mod_pos0 + 1
        snp_pos0 = snp_pos1 - 1
        gen_off = snp_pos0 - mod_pos0            # genomic offset of the SNP from the modified base
        d = abs(gen_off)

        # ---- Axis A: positional ladder ----
        if d == 0:
            positional = "AT_MOD_BASE"
        elif d <= 2:
            positional = "IN_MOTIF_CORE"
        elif d <= 4:
            positional = "IN_MOTIF_EXTENDED"
        elif d <= int(args.proximal_bp):
            positional = "PROXIMAL_CIS"
        else:
            positional = "DISTAL_CIS"

        ref5 = _kmer(fa, s_chrom, mod_pos0, strand, 2)
        ref9 = _kmer(fa, s_chrom, mod_pos0, strand, 4)
        alt5 = _sub(ref5, 2, gen_off, s_alt, strand) if ref5 else None
        alt9 = _sub(ref9, 4, gen_off, s_alt, strand) if ref9 else None

        # SNP alleles are reported in genomic orientation; the motif/mod live in transcript orientation.
        ref_t = s_ref.translate(_COMP) if strand == "-" else s_ref
        alt_t = s_alt.translate(_COMP) if strand == "-" else s_alt

        # ---- Axis B: motif effect (m6A/DRACH only, and only when the SNP is inside the 5-mer) ----
        motif_effect, predicted = "NOT_APPLICABLE", None
        if code == M6A_CODE and ref5 and alt5 and d <= 2:
            ref_ok, alt_ok = bool(DRACH_RE.match(ref5)), bool(DRACH_RE.match(alt5))
            if ref_ok and not alt_ok:
                motif_effect, predicted = "MOTIF_DISRUPTED", "alt_lower"
            elif alt_ok and not ref_ok:
                motif_effect, predicted = "MOTIF_CREATED", "alt_higher"
            elif ref_ok and alt_ok:
                motif_effect = "MOTIF_PRESERVED"
            else:
                motif_effect = "MOTIF_ABSENT_BOTH"   # neither allele is DRACH -> m6A call suspect

        # ---- Self-reporting / definitional artifacts: the "SNP" and the modification are one event ----
        artifact, ablated = "NONE", False
        if d == 0:
            sig = SELF_REPORT_SIGNATURE.get(str(code))
            if sig and ref_t == sig[0] and alt_t == sig[1]:
                artifact = sig[2]                    # the variant IS the edit / the pseU basecall error
            elif code == M6A_CODE and alt5 and alt5[2] != "A":
                artifact, ablated = "MOD_BASE_ABLATED", True   # no A on alt -> nothing to methylate

        # ---- Axis C: observed direction + concordance ----
        ref_rate = alt_rate = None
        observed = "UNKNOWN"
        try:
            st = json.loads(r.per_state_json)
            rm, rn = float(st.get("ref_modified", 0)), float(st.get("ref_not_target", 0))
            am, an = float(st.get("alt_modified", 0)), float(st.get("alt_not_target", 0))
            if (rm + rn) > 0 and (am + an) > 0:
                ref_rate, alt_rate = rm / (rm + rn), am / (am + an)
                observed = "alt_lower" if alt_rate < ref_rate else ("alt_higher" if alt_rate > ref_rate else "equal")
        except Exception:
            pass
        if predicted is None or observed in ("UNKNOWN", "equal"):
            concordance = "NOT_TESTABLE"
        else:
            concordance = "CONCORDANT" if predicted == observed else "DISCORDANT"

        rows.append({
            "snp_id": r.snp_id, "mod_site_id": r.mod_site_id, "chrom": s_chrom,
            "snp_pos1": snp_pos1, "mod_pos1": mod_pos1, "strand": strand, "target_mod_code": code,
            "gene_names": getattr(r, "gene_names", ""), "distance_bp": d,
            "positional_class": positional, "snp_ref_tx": ref_t, "snp_alt_tx": alt_t,
            "ref_5mer": ref5 or "", "alt_5mer": alt5 or "",
            "ref_9mer": ref9 or "", "alt_9mer": alt9 or "",
            "motif_effect": motif_effect, "artifact_flag": artifact, "mod_base_ablated": ablated,
            "ref_mod_rate": round(ref_rate, 4) if ref_rate is not None else "",
            "alt_mod_rate": round(alt_rate, 4) if alt_rate is not None else "",
            "observed_direction": observed, "predicted_direction": predicted or "",
            "direction_concordance": concordance,
            "class_key": f"{positional}__{motif_effect}",
            "n_reads": getattr(r, "n_reads", ""),
            "effect_abs_delta_mod_frac": getattr(r, "effect_abs_delta_mod_frac", ""),
            "p_adj_bh": getattr(r, "p_adj_bh", ""),
        })
    fa.close()

    out = pd.DataFrame(rows)
    if out.empty:
        out = pd.DataFrame(columns=OUT_COLS)
    else:
        out = out[OUT_COLS].sort_values(["p_adj_bh", "distance_bp"]).reset_index(drop=True)
    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
    out.to_csv(args.out_tsv, sep="\t", index=False)
    if args.verbose and not out.empty:
        print(f"[snp_mod_mech] {len(out)} pairs; positional: " +
              ", ".join(f"{k}={v}" for k, v in out['positional_class'].value_counts().items()), flush=True)


if __name__ == "__main__":
    main()
