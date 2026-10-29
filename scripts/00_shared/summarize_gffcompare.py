#!/usr/bin/env python3
"""Convert gffcompare statistics into a compact report table."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


PATTERNS = [
    ("query_transcripts", r"Query mRNAs\s*:\s*(\d+)"),
    ("query_loci", r"Query mRNAs\s*:\s*\d+\s+in\s+(\d+) loci"),
    ("multi_exon_transcripts", r"\((\d+) multi-exon transcripts\)"),
    ("matching_intron_chains", r"Matching intron chains:\s*(\d+)"),
    ("matching_transcripts", r"Matching transcripts:\s*(\d+)"),
    ("matching_loci", r"Matching loci:\s*(\d+)"),
    ("novel_exons", r"Novel exons:\s*(\d+)/"),
    ("novel_introns", r"Novel introns:\s*(\d+)/"),
    ("novel_loci", r"Novel loci:\s*(\d+)/"),
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stats", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    text = Path(args.stats).read_text(encoding="utf-8")
    rows = []
    for metric, pattern in PATTERNS:
        match = re.search(pattern, text)
        if not match:
            raise SystemExit(f"Missing gffcompare metric: {metric}")
        rows.append((metric, match.group(1)))
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerows(rows)


if __name__ == "__main__":
    main()
