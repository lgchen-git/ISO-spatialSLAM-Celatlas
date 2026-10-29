#!/usr/bin/env python3
"""Convert LR per-read metadata into the compact SLAM barcode table."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build read_id/barcode/UMI/x/y table for SLAM calling.")
    parser.add_argument("--per-read-tsv", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_path = outdir / f"{args.sample_id}.slam_barcode_table.tsv.gz"
    qc = {"sample": args.sample_id, "total_reads": 0, "kept_reads": 0, "failed_poly_t": 0, "barcode_not_in_whitelist": 0, "reads_without_spatial_coordinates": 0}

    with gzip.open(args.per_read_tsv, "rt") as in_handle, gzip.open(out_path, "wt") as out_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        out_handle.write("read_id\tbarcode\tumi\tx\ty\n")
        for row in reader:
            qc["total_reads"] += 1
            qc["kept_reads"] += 1
            if int(row["x"]) < 0 or int(row["y"]) < 0:
                qc["reads_without_spatial_coordinates"] += 1
            out_handle.write(
                "\t".join(
                    [
                        row["read_id"],
                        row["barcode_corr"],
                        row["umi_seq"],
                        row["x"],
                        row["y"],
                    ]
                )
                + "\n"
            )

    (outdir / f"{args.sample_id}.barcode_qc.json").write_text(json.dumps(qc, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
