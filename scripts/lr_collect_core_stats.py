#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
from collections import Counter, defaultdict

import pysam


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Collect core per-read, per-barcode and alignment summaries.")
    parser.add_argument("--structure-tsv", required=True)
    parser.add_argument("--barcode-tsv", required=True)
    parser.add_argument("--per-read-tsv", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def read_gz_tsv(path: str):
    with gzip.open(path, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            yield row


def main() -> None:
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    structure_counts = Counter()
    barcode_counts = Counter()
    barcode_to_reads = Counter()
    barcode_to_source = {}
    confidence_counts = Counter()
    per_read = {}

    for row in read_gz_tsv(args.structure_tsv):
        structure_counts[row["scan_status"]] += 1

    for row in read_gz_tsv(args.barcode_tsv):
        barcode_counts[row["assign_status"]] += 1

    for row in read_gz_tsv(args.per_read_tsv):
        per_read[row["read_id"]] = row
        barcode = row["barcode_corr"]
        barcode_to_reads[barcode] += 1
        barcode_to_source.setdefault(barcode, row.get("structure_source", "."))
        confidence_counts[row.get("confidence_class", "reject")] += 1

    align_summary = Counter()
    barcode_mapped = Counter()
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for record in bam.fetch(until_eof=True):
            if record.is_unmapped:
                align_summary["unmapped"] += 1
                continue
            align_summary["mapped"] += 1
            if record.is_secondary:
                align_summary["secondary"] += 1
            if record.is_supplementary:
                align_summary["supplementary"] += 1
            read_id = record.query_name
            if read_id in per_read:
                barcode_mapped[per_read[read_id]["barcode_corr"]] += 1

    with open(os.path.join(args.outdir, "structure_summary.tsv"), "w", encoding="utf-8") as handle:
        handle.write("category\tcount\n")
        for key, value in sorted(structure_counts.items()):
            handle.write(f"{key}\t{value}\n")

    with open(os.path.join(args.outdir, "barcode_summary.tsv"), "w", encoding="utf-8") as handle:
        handle.write("category\tcount\n")
        for key, value in sorted(barcode_counts.items()):
            handle.write(f"{key}\t{value}\n")

    with open(os.path.join(args.outdir, "confidence_summary.tsv"), "w", encoding="utf-8") as handle:
        handle.write("confidence_class\tcount\n")
        for key, value in sorted(confidence_counts.items()):
            handle.write(f"{key}\t{value}\n")

    with open(os.path.join(args.outdir, "barcode_level_summary.tsv"), "w", encoding="utf-8") as handle:
        handle.write("barcode\treads_total\treads_mapped\tstructure_source\n")
        for barcode, total_reads in barcode_to_reads.most_common():
            source = barcode_to_source.get(barcode, ".")
            handle.write(f"{barcode}\t{total_reads}\t{barcode_mapped.get(barcode, 0)}\t{source}\n")

    with open(os.path.join(args.outdir, "alignment_summary.json"), "w", encoding="utf-8") as handle:
        json.dump(dict(align_summary), handle, indent=2, sort_keys=True)


if __name__ == "__main__":
    main()
