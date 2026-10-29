#!/usr/bin/env python3
"""Compute configurable quality gates from QC, extraction and primary-BAM outputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path


def metric(path: Path, key: str) -> int:
    with path.open(encoding="utf-8") as handle:
        values = dict(line.rstrip("\n").split("\t", 1) for line in handle.readlines()[1:] if "\t" in line)
    return int(float(values[key]))


def flagstat_total(path: Path) -> int:
    for line in path.read_text(encoding="utf-8").splitlines():
        if " in total " in line:
            return int(line.split("+", 1)[0].strip())
    raise ValueError(f"No total line in {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qc-metrics", required=True, type=Path)
    parser.add_argument("--structure-summary", required=True, type=Path)
    parser.add_argument("--extract-summary", required=True, type=Path)
    parser.add_argument("--primary-flagstat", required=True, type=Path)
    parser.add_argument("--per-read", required=True, type=Path)
    parser.add_argument("--chip-bounds", type=int, required=True)
    parser.add_argument("--min-structure-rate", type=float, required=True)
    parser.add_argument("--min-cdna-rate", type=float, required=True)
    parser.add_argument("--min-mapq20-rate", type=float, required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    raw = metric(args.qc_metrics, "n_reads")
    structure = json.loads(args.structure_summary.read_text(encoding="utf-8"))
    extract = json.loads(args.extract_summary.read_text(encoding="utf-8"))
    kept = int(extract.get("counts", {}).get("kept", 0))
    mapq20 = flagstat_total(args.primary_flagstat)
    coordinate_ok = True
    with gzip.open(args.per_read, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if not (0 <= int(row["x"]) <= args.chip_bounds and 0 <= int(row["y"]) <= args.chip_bounds):
                coordinate_ok = False
                break
    structure_rate = structure.get("candidate_found", 0) / raw if raw else 0.0
    cdna_rate = kept / raw if raw else 0.0
    mapq20_rate = mapq20 / kept if kept else 0.0
    checks = {
        "structure_rate": structure_rate >= args.min_structure_rate,
        "cdna_rate": cdna_rate >= args.min_cdna_rate,
        "mapq20_rate": mapq20_rate >= args.min_mapq20_rate,
        "coordinates_valid": coordinate_ok,
    }
    result = {
        "overall": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "raw_reads": raw,
        "kept_cdna_reads": kept,
        "primary_mapq20_reads": mapq20,
        "structure_rate": structure_rate,
        "cdna_rate": cdna_rate,
        "mapq20_rate": mapq20_rate,
        "coordinates_valid": coordinate_ok,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
