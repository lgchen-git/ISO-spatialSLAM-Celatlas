#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter

RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract valid cDNA FASTQ after structure scan and barcode assignment.")
    parser.add_argument("--fastq", required=True)
    parser.add_argument("--structure-tsv", required=True)
    parser.add_argument("--barcode-tsv", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument(
        "--retain-all-barcodes",
        action="store_true",
        help="Retain candidate reads even when barcode is off-whitelist or ambiguous. Reads without valid spatial coordinates keep x/y=-1.",
    )
    return parser.parse_args()


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(RC_TABLE)[::-1]


def main() -> None:
    args = parse_args()
    import os
    os.makedirs(args.outdir, exist_ok=True)
    valid_out = f"{args.outdir}/{args.sample_id}.cdna.fastq.gz"
    per_read_out = f"{args.outdir}/per_read.tsv.gz"
    discard_out = f"{args.outdir}/discard_reason.tsv.gz"
    summary_out = f"{args.outdir}/extract_summary.json"
    counts = Counter()

    with open_text(args.fastq) as in_handle, gzip.open(args.structure_tsv, "rt") as structure_handle, gzip.open(args.barcode_tsv, "rt") as barcode_handle, gzip.open(valid_out, "wt") as fastq_handle, gzip.open(per_read_out, "wt") as per_handle, gzip.open(discard_out, "wt") as discard_handle:
        structure = csv.DictReader(structure_handle, delimiter="\t")
        barcode = csv.DictReader(barcode_handle, delimiter="\t")
        per_handle.write("read_id\tassign_status\tconfidence_class\tstructure_source\tbarcode_corr\tumi_seq\tstrand\tcdna_len\tx\ty\n")
        discard_handle.write("read_id\treason\n")
        while True:
            header = in_handle.readline()
            if not header:
                break
            seq = in_handle.readline().strip()
            in_handle.readline()
            qual = in_handle.readline().strip()
            read_id = header.strip()[1:].split()[0]
            try:
                srow = next(structure)
                brow = next(barcode)
            except StopIteration as exc:
                raise RuntimeError("Structure/barcode table ended before FASTQ") from exc
            if srow["read_id"] != read_id or brow["read_id"] != read_id:
                raise RuntimeError(
                    f"Read order mismatch: fastq={read_id}, structure={srow['read_id']}, barcode={brow['read_id']}"
                )
            if srow.get("scan_status") != "candidate_found":
                counts[srow.get("scan_status", "structure_not_found")] += 1
                discard_handle.write(f"{read_id}\t{srow.get('scan_status', 'structure_not_found')}\n")
                continue
            has_spatial_barcode = brow["assign_status"] in {"exact", "rescued"} and brow.get("confidence_class", "reject") != "reject"
            if not has_spatial_barcode and not args.retain_all_barcodes:
                counts[brow["assign_status"]] += 1
                discard_handle.write(f"{read_id}\t{brow['assign_status']}\n")
                continue
            strand = srow["strand"]
            source = srow.get("structure_source", ".")
            conf = brow.get("confidence_class", "reject")
            barcode_corr = brow["barcode_corr"] if has_spatial_barcode else brow.get("barcode_raw", ".")
            oriented_seq = seq if strand == "+" else reverse_complement(seq)
            oriented_qual = qual if strand == "+" else qual[::-1]
            start = int(srow["cdna_start"])
            end = int(srow["cdna_end"])
            if start < 0 or end <= start:
                counts["invalid_cdna_interval"] += 1
                discard_handle.write(f"{read_id}\tinvalid_cdna_interval\n")
                continue
            cdna_seq = oriented_seq[start:end]
            cdna_qual = oriented_qual[start:end]
            counts["kept"] += 1
            if not has_spatial_barcode:
                counts["kept_without_spatial_barcode"] += 1
            fastq_handle.write(f"@{read_id}\n{cdna_seq}\n+\n{cdna_qual}\n")
            per_handle.write(
                f"{read_id}\t{brow['assign_status']}\t{conf}\t{source}\t{barcode_corr}\t{srow['umi_seq']}\t{strand}\t{len(cdna_seq)}\t{brow['x']}\t{brow['y']}\n"
            )

    with open(summary_out, "w", encoding="utf-8") as handle:
        json.dump({"counts": dict(counts)}, handle, indent=2, sort_keys=True)


if __name__ == "__main__":
    main()
