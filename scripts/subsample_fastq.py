#!/usr/bin/env python3
"""Write the first N FASTQ records to a new FASTQ file."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Subsample FASTQ records without randomization.")
    parser.add_argument("--fastq", required=True)
    parser.add_argument("--reads", type=int, required=True)
    parser.add_argument("--out", required=True)
    return parser.parse_args()


def open_read(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r", encoding="utf-8")


def open_write(path: str):
    return gzip.open(path, "wt") if path.endswith(".gz") else open(path, "w", encoding="utf-8")


def main() -> None:
    args = parse_args()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with open_read(args.fastq) as in_handle, open_write(str(out)) as out_handle:
        while written < args.reads:
            header = in_handle.readline()
            if not header:
                break
            seq = in_handle.readline()
            plus = in_handle.readline()
            qual = in_handle.readline()
            if not qual:
                break
            out_handle.write(header)
            out_handle.write(seq)
            out_handle.write(plus)
            out_handle.write(qual)
            written += 1
    print(f"subsampled_reads\t{written}")


if __name__ == "__main__":
    main()
