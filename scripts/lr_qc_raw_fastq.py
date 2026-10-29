#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import statistics
import sys
from pathlib import Path

from lr_common import command_path, ensure_dir, run_command


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run raw FASTQ QC.")
    parser.add_argument("--fastq", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--skip-nanoplot", action="store_true")
    parser.add_argument("--skip-fastqc", action="store_true")
    return parser.parse_args()


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8")


def mean_q(qual: str) -> float:
    values = [ord(ch) - 33 for ch in qual.rstrip()]
    return round(sum(values) / len(values), 4) if values else 0.0


def n50(lengths: list[int]) -> int:
    ordered = sorted(lengths, reverse=True)
    half = sum(ordered) / 2 if ordered else 0
    running = 0
    for value in ordered:
        running += value
        if running >= half:
            return value
    return 0


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    ensure_dir(outdir)

    import subprocess

    if command_path("seqkit"):
        with (outdir / "seqkit.stats.tsv").open("w", encoding="utf-8") as handle:
            subprocess.run(["seqkit", "stats", "-Ta", args.fastq], check=True, stdout=handle)

    lengths: list[int] = []
    qualities: list[float] = []
    reads = 0
    bases = 0
    poly_t20 = 0
    with open_text(args.fastq) as in_handle, gzip.open(outdir / "per_read_qc.tsv.gz", "wt") as out_handle:
        out_handle.write("read_id\tread_len\tmean_q\tpoly_t20plus\n")
        while True:
            header = in_handle.readline()
            if not header:
                break
            seq = in_handle.readline().strip()
            in_handle.readline()
            qual = in_handle.readline().strip()
            read_id = header.strip()[1:].split()[0]
            read_len = len(seq)
            q = mean_q(qual)
            has_poly_t = int("T" * 20 in seq)
            out_handle.write(f"{read_id}\t{read_len}\t{q}\t{has_poly_t}\n")
            reads += 1
            bases += read_len
            poly_t20 += has_poly_t
            lengths.append(read_len)
            qualities.append(q)

    with (outdir / "qc_metrics.tsv").open("w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"n_reads\t{reads}\n")
        handle.write(f"bases\t{bases}\n")
        handle.write(f"mean_len\t{round(bases / reads, 4) if reads else 0}\n")
        handle.write(f"median_len\t{int(statistics.median(lengths)) if lengths else 0}\n")
        handle.write(f"N50\t{n50(lengths)}\n")
        handle.write(f"mean_q\t{round(statistics.mean(qualities), 4) if qualities else 0}\n")
        handle.write(f"polyT20plus_rate\t{round(poly_t20 / reads, 6) if reads else 0}\n")

    nanoplot = command_path("NanoPlot")
    if nanoplot and not args.skip_nanoplot:
        subprocess.run(
            [sys.executable, nanoplot, "--fastq", args.fastq, "--threads", str(args.threads), "--outdir", str(outdir / "NanoPlot"), "--plots", "hex", "dot", "--no_static", "--tsv_stats"],
            check=True,
        )
    if command_path("fastqc") and not args.skip_fastqc:
        ensure_dir(outdir / "fastqc")
        run_command(["fastqc", "-t", str(args.threads), "-o", str(outdir / "fastqc"), args.fastq])


if __name__ == "__main__":
    main()
