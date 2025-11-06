#!/usr/bin/env python3
import argparse, os, re, glob, sys
from collections import defaultdict, OrderedDict

def read_summary(summary_tsv):
    code_meta = {}
    with open(summary_tsv) as f:
        for line in f:
            if not line or line.startswith("#"): continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 18:  # guard
                continue
            code = parts[0]
            gname = parts[6]
            tid   = parts[7]
            code_meta[code] = (gname, tid)
    return code_meta

def parse_bedmethyl(path):
    """Yield tuples: (mod_code, Nvalid_cov, Nmod)."""
    with open(path) as f:
        for line in f:
            if not line or line[0] in "#tT":  # header lines may start with track/browser/#chrom...
                continue
            cols = line.strip().split()
            if len(cols) < 12:
                continue
            mod_code = cols[3]
            try:
                nvalid = int(cols[10])
                nmod   = int(cols[12])
            except Exception:
                continue
            yield (mod_code, nvalid, nmod)

def main():
    ap = argparse.ArgumentParser(description="Aggregate modkit bedMethyl partitioned by ZT tag (transcript code)")
    ap.add_argument("--modkit-dir", required=True, help="Directory containing SAMPLE_<ZT>.bed files from `modkit pileup --partition-tag ZT --prefix SAMPLE`")
    ap.add_argument("--summary-tsv", required=True, help="*_classification_summary.tsv from assembler v6")
    ap.add_argument("--mods", nargs="+", default=["17596","a","m","17802","69426","19228","19229","19227"],
                    help="Which mod codes to aggregate")
    ap.add_argument("--out-tsv", default="modkit_by_transcript.tsv")
    args = ap.parse_args()

    code_meta = read_summary(args.summary_tsv)  # code -> (gene_name, matched_tid)

    # collect all bed files
    beds = glob.glob(os.path.join(args.modkit_dir, "*.bed"))
    if not beds:
        sys.exit("No .bed files found in modkit-dir")

    # filename pattern: <SAMPLE>_<CODE>.bed  (CODE has no underscores in our assembler)
    # sample may contain underscores; take CODE as the suffix after last underscore
    recs = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: [0,0])))
    # recs[mod_code][code][sample] = [sum_nvalid, sum_nmod]

    samples = OrderedDict()
    codes_seen = set()

    for bed in beds:
        fn = os.path.basename(bed)
        if not fn.endswith(".bed"): continue
        root = fn[:-4]
        if "_" not in root: 
            # skip non-partitioned outputs
            continue
        code = root.split("_")[-1]
        sample = root[:-(len(code)+1)]
        samples[sample] = True
        if code not in code_meta:
            # still allow aggregation; gene/tid will be NA
            pass
        sums = defaultdict(lambda: [0,0])
        for mod_code, nvalid, nmod in parse_bedmethyl(bed):
            if mod_code not in args.mods: 
                continue
            sums[mod_code][0] += nvalid
            sums[mod_code][1] += nmod
        for mod_code, (nv, nm) in sums.items():
            recs[mod_code][code][sample][0] += nv
            recs[mod_code][code][sample][1] += nm
        if sums:
            codes_seen.add(code)

    samples = list(samples.keys())
    with open(args.out_tsv, "w") as out:
        header = ["code","gene_name","matched_tid","mod_code"]
        for s in samples:
            header += [f"{s}_frac", f"{s}_cov"]
        out.write("\t".join(header)+"\n")

        for mod_code in args.mods:
            if mod_code not in recs: 
                continue
            for code, by_sample in recs[mod_code].items():
                gname, tid = code_meta.get(code, ("NA","NA"))
                row = [code, gname, tid, mod_code]
                for s in samples:
                    nv, nm = by_sample.get(s, [0,0])
                    frac = (nm/nv) if nv>0 else "NA"
                    row += [f"{frac:.6f}" if isinstance(frac,float) else "NA", str(nv)]
                out.write("\t".join(row)+"\n")

    print(f"[OK] Wrote {args.out_tsv}", file=sys.stderr)

if __name__ == "__main__":
    main()

