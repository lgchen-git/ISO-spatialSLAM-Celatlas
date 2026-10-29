#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter
from typing import Dict, Optional, Tuple

DNA_BASES = ("A", "C", "G", "T")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Assign parsed barcodes to whitelist positions.")
    parser.add_argument("--structure-tsv", required=True)
    parser.add_argument("--whitelist-csv", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--anchor-barcode-max-mismatch", type=int, default=1)
    parser.add_argument("--fallback-barcode-max-mismatch", type=int, default=1)
    return parser.parse_args()


def load_whitelist(path: str) -> Dict[str, Tuple[int, int]]:
    data: Dict[str, Tuple[int, int]] = {}
    with open(path, "r", encoding="utf-8") as handle:
      for row in csv.reader(handle):
          if len(row) < 3:
              continue
          bc = row[2].strip().upper()
          if not bc or bc == "barcode":
              continue
          data[bc] = (int(row[0]), int(row[1]))
    return data


def hamming_distance(a: str, b: str) -> int:
    return sum(1 for x, y in zip(a, b) if x != y)


def iter_hamming_neighbors(seq: str, max_mm: int):
    yield seq
    letters = list(seq)
    if max_mm >= 1:
        for i, base in enumerate(letters):
            for alt in DNA_BASES:
                if alt == base:
                    continue
                copy = letters.copy()
                copy[i] = alt
                yield "".join(copy)
    if max_mm >= 2:
        for i in range(len(letters)):
            for j in range(i + 1, len(letters)):
                for alt_i in DNA_BASES:
                    if alt_i == letters[i]:
                        continue
                    for alt_j in DNA_BASES:
                        if alt_j == letters[j]:
                            continue
                        copy = letters.copy()
                        copy[i] = alt_i
                        copy[j] = alt_j
                        yield "".join(copy)


def correct_barcode(barcode_raw: str, whitelist: Dict[str, Tuple[int, int]], max_mm: int) -> tuple[str, Optional[str], Optional[int], int, int]:
    if barcode_raw in whitelist:
        x, y = whitelist[barcode_raw]
        return "exact", barcode_raw, 0, x, y
    best_dist: Optional[int] = None
    best_candidates: set[str] = set()
    for candidate in iter_hamming_neighbors(barcode_raw, max_mm):
        if candidate == barcode_raw or candidate not in whitelist:
            continue
        dist = hamming_distance(barcode_raw, candidate)
        if best_dist is None or dist < best_dist:
            best_dist = dist
            best_candidates = {candidate}
        elif dist == best_dist:
            best_candidates.add(candidate)
    if best_dist is None:
        return "barcode_not_in_whitelist", None, None, -1, -1
    if len(best_candidates) > 1:
        return "ambiguous_barcode", None, None, -1, -1
    barcode = next(iter(best_candidates))
    x, y = whitelist[barcode]
    return "rescued", barcode, best_dist, x, y


def confidence_class(source: str, assign_status: str) -> str:
    if source == "anchor_5p":
        if assign_status == "exact":
            return "A"
        if assign_status == "rescued":
            return "B"
    if source == "linker1_fallback":
        if assign_status == "exact":
            return "C"
        if assign_status == "rescued":
            return "D"
    return "reject"


def main() -> None:
    args = parse_args()
    whitelist = load_whitelist(args.whitelist_csv)
    import os
    os.makedirs(args.outdir, exist_ok=True)
    out_path = f"{args.outdir}/barcode_assign.tsv.gz"
    summary_path = f"{args.outdir}/barcode_summary.json"
    counts = Counter()
    ed_counts = Counter()

    with gzip.open(args.structure_tsv, "rt") as in_handle, gzip.open(out_path, "wt") as out_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        out_handle.write(
            "read_id\tscan_status\tstructure_source\tbarcode_raw\tassign_status\tconfidence_class\tbarcode_corr\twhitelist_ed\tx\ty\n"
        )
        for row in reader:
            scan_status = row["scan_status"]
            source = row.get("structure_source", ".")
            if scan_status != "candidate_found":
                counts[scan_status] += 1
                out_handle.write(f"{row['read_id']}\t{scan_status}\t{source}\t{row['barcode_raw']}\t{scan_status}\treject\t.\t-1\t-1\t-1\n")
                continue
            max_mm = args.anchor_barcode_max_mismatch if source == "anchor_5p" else args.fallback_barcode_max_mismatch
            assign_status, barcode_corr, ed, x, y = correct_barcode(row["barcode_raw"], whitelist, max_mm)
            conf = confidence_class(source, assign_status)
            counts[assign_status] += 1
            if ed is not None:
                ed_counts[str(ed)] += 1
                if assign_status == "rescued":
                    counts[f"mismatch_{ed}"] += 1
            out_handle.write(
                f"{row['read_id']}\t{scan_status}\t{source}\t{row['barcode_raw']}\t{assign_status}\t{conf}\t{barcode_corr or '.'}\t{ed if ed is not None else -1}\t{x}\t{y}\n"
            )

    with open(summary_path, "w", encoding="utf-8") as handle:
        json.dump({"assign_counts": dict(counts), "edit_distance_counts": dict(ed_counts)}, handle, indent=2, sort_keys=True)


if __name__ == "__main__":
    main()
