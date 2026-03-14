#!/usr/bin/env python3

import argparse
import csv
import os
from collections import defaultdict


def parse_args():
    ap = argparse.ArgumentParser(
        description="Aggregate per-sample multigene scrap transcript counts into a transcript x sample matrix."
    )
    ap.add_argument("--counts", nargs="+", required=True, help="Per-sample scrap transcript count TSVs")
    ap.add_argument("--out", required=True, help="Output transcript x sample TSV")
    return ap.parse_args()


def main():
    args = parse_args()

    counts_by_code = defaultdict(dict)
    samples = set()

    for path in sorted(args.counts):
        if not os.path.exists(path):
            continue
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                sample = row.get("sample", "").strip()
                code = row.get("code", "").strip()
                if not sample or not code:
                    continue
                samples.add(sample)
                counts_by_code[code][sample] = int(row.get("scrapped_assigned_reads", 0) or 0)

    sample_cols = sorted(samples)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

    with open(args.out, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["code"] + sample_cols)
        for code in sorted(counts_by_code):
            writer.writerow([code] + [counts_by_code[code].get(sample, 0) for sample in sample_cols])


if __name__ == "__main__":
    main()
