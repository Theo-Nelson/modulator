#!/usr/bin/env python3
"""Convert a per-read table between TSV and parquet (one row group per chromosome block).

    convert_table_format.py --in <table.tsv|.parquet> --out <table.parquet|.tsv> [--chunk-rows N]

The parquet form is ~5-10x smaller (zstd + dictionary encoding) and every modulator consumer reads
either form through genotype_utils.open_chrom_table; the `.chromidx.tsv` sidecar is written for both.
"""
import argparse
import os
import sys

import pandas as pd

from genotype_utils import ParquetBlockWriter, is_parquet_path, open_chrom_table, write_chrom_index


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--chunk-rows", type=int, default=2_000_000)
    args = ap.parse_args()
    tbl = open_chrom_table(args.inp)
    cols = list(tbl.header_cols)
    if is_parquet_path(args.out):
        w = ParquetBlockWriter(args.out, cols)
        for chrom in tbl.chroms:
            for chunk in tbl.iter_chunks(chrom, chunksize=args.chunk_rows):
                w.write(chrom, chunk)
        w.close()
    else:
        tmp = args.out + ".tmp"
        blocks = []
        with open(tmp, "w") as fh:
            wrote = False
            if not tbl.chroms:
                pd.DataFrame(columns=cols).to_csv(fh, sep="\t", index=False)
            for chrom in tbl.chroms:
                off = fh.tell()
                for chunk in tbl.iter_chunks(chrom, chunksize=args.chunk_rows):
                    chunk.to_csv(fh, sep="\t", index=False, header=not wrote)
                    wrote = True
                fh.flush()
                blocks.append((chrom, off, fh.tell() - off))
        os.replace(tmp, args.out)
        write_chrom_index(args.out, blocks)
    print(f"[info] wrote {args.out} ({os.path.getsize(args.out) / 1e9:.2f} GB) from {args.inp} "
          f"({os.path.getsize(args.inp) / 1e9:.2f} GB)", file=sys.stderr)


if __name__ == "__main__":
    main()
