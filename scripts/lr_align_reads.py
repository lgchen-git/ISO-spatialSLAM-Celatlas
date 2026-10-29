#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path

from lr_common import ensure_dir, require_file


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Align valid cDNA reads with minimap2 and samtools.")
    parser.add_argument("--fastq", required=True)
    parser.add_argument("--mmi", required=True)
    parser.add_argument("--junc-bed", default="")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--workdir", required=True, help="Scratch directory for the unfiltered sorted BAM.")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--mm2-preset", default="splice")
    parser.add_argument("--mm2-extra", default="")
    parser.add_argument("--minimap2-bin", default="minimap2")
    parser.add_argument("--samtools-bin", default="samtools")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    require_file(args.fastq)
    require_file(args.mmi)
    if args.junc_bed:
        require_file(args.junc_bed)
    outdir = Path(args.outdir)
    workdir = Path(args.workdir)
    ensure_dir(outdir)
    ensure_dir(workdir)

    sorted_bam = workdir / f"{args.sample_id}.sorted.bam"
    primary_bam = outdir / f"{args.sample_id}.primary.q{args.min_mapq}.bam"

    mm2_cmd = [
        args.minimap2_bin,
        "-ax",
        args.mm2_preset,
        "--secondary=no",
        "--MD",
        "--cs=long",
        "-t",
        str(args.threads),
    ]
    if args.junc_bed:
        mm2_cmd.extend(["--junc-bed", args.junc_bed])
    if args.mm2_extra.strip():
        mm2_cmd.extend(args.mm2_extra.split())
    mm2_cmd.extend([args.mmi, args.fastq])

    sort_cmd = [args.samtools_bin, "sort", "-@", str(args.threads), "-o", str(sorted_bam), "-"]
    mm2 = subprocess.Popen(mm2_cmd, stdout=subprocess.PIPE)
    assert mm2.stdout is not None
    sort_proc = subprocess.Popen(sort_cmd, stdin=mm2.stdout)
    mm2.stdout.close()
    if sort_proc.wait() != 0:
        raise SystemExit("ERROR: samtools sort failed")
    if mm2.wait() != 0:
        raise SystemExit("ERROR: minimap2 failed")

    with (outdir / f"{args.sample_id}.flagstat.txt").open("w", encoding="utf-8") as handle:
        subprocess.run([args.samtools_bin, "flagstat", str(sorted_bam)], check=True, stdout=handle)
    with (outdir / f"{args.sample_id}.idxstats.tsv").open("w", encoding="utf-8") as handle:
        subprocess.run([args.samtools_bin, "idxstats", str(sorted_bam)], check=True, stdout=handle)
    with (outdir / f"{args.sample_id}.samtools.stats.txt").open("w", encoding="utf-8") as handle:
        subprocess.run([args.samtools_bin, "stats", str(sorted_bam)], check=True, stdout=handle)
    subprocess.run([args.samtools_bin, "index", "-@", str(args.threads), str(sorted_bam)], check=True)

    with primary_bam.open("wb") as handle:
        subprocess.run([args.samtools_bin, "view", "-b", "-q", str(args.min_mapq), "-F", "2304", str(sorted_bam)], check=True, stdout=handle)
    subprocess.run([args.samtools_bin, "index", "-@", str(args.threads), str(primary_bam)], check=True)
    with (outdir / f"{args.sample_id}.primary.flagstat.txt").open("w", encoding="utf-8") as handle:
        subprocess.run([args.samtools_bin, "flagstat", str(primary_bam)], check=True, stdout=handle)


if __name__ == "__main__":
    main()
