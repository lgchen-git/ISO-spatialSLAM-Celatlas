#!/usr/bin/env python3
"""Summarize strand-agnostic substitution channels in a primary-read BAM."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import pandas as pd
import pysam

BASES = "ACGT"


def load_vcf(path: str) -> set[tuple[str, int]]:
    sites: set[tuple[str, int]] = set()
    if not path:
        return sites
    opener = pysam.VariantFile(path)
    for record in opener:
        sites.add((record.contig, record.pos - 1))
    opener.close()
    return sites


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--exclude-vcf", default="")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--min-baseq", type=int, default=20)
    parser.add_argument("--max-reads", type=int, default=200000, help="Maximum passing primary reads; 0 uses all.")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    excluded = load_vcf(args.exclude_vcf)
    counts: Counter[str] = Counter()
    coverage: Counter[str] = Counter()
    qc: Counter[str] = Counter()
    fasta = pysam.FastaFile(args.fasta)
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            qc["bam_records"] += 1
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.mapping_quality < args.min_mapq:
                qc["skipped_records"] += 1
                continue
            if args.max_reads and qc["processed_primary_reads"] >= args.max_reads:
                break
            qc["processed_primary_reads"] += 1
            qualities = read.query_qualities
            try:
                aligned = read.get_aligned_pairs(matches_only=True, with_seq=True)
            except ValueError:  # MD tag absent: retain a reference-FASTA fallback.
                aligned = [(qpos, rpos, None) for qpos, rpos in read.get_aligned_pairs(matches_only=True)]
            for qpos, rpos, md_ref in aligned:
                if qpos is None or rpos is None or (qualities is not None and qualities[qpos] < args.min_baseq):
                    continue
                chrom = read.reference_name
                if (chrom, rpos) in excluded:
                    qc["excluded_snp_bases"] += 1
                    continue
                ref = md_ref.upper() if md_ref else fasta.fetch(chrom, rpos, rpos + 1).upper()
                query = read.query_sequence[qpos].upper()
                if ref not in BASES or query not in BASES:
                    continue
                coverage[ref] += 1
                if ref != query:
                    counts[f"{ref}>{query}"] += 1
    fasta.close()

    rows = []
    for ref in BASES:
        for alt in BASES:
            if ref == alt:
                continue
            channel = f"{ref}>{alt}"
            rows.append({"channel": channel, "count": counts[channel], "reference_coverage": coverage[ref], "rate": counts[channel] / coverage[ref] if coverage[ref] else 0.0})
    frame = pd.DataFrame(rows)
    frame.to_csv(outdir / "mismatch_spectrum.tsv", sep="\t", index=False)
    background_rows = frame[~frame["channel"].isin(["T>C", "A>G"])]
    background = max(1e-5, float(background_rows["rate"].max()) if not background_rows.empty else 1e-5)
    summary = {
        **dict(qc),
        "excluded_snp_sites": len(excluded),
        "variant_exclusion_mode": "provided_vcf" if args.exclude_vcf else "none",
        "variant_exclusion_vcf": str(Path(args.exclude_vcf).resolve()) if args.exclude_vcf else "",
        "background_rate": background,
        "background_definition": "maximum non-T>C/non-A>G substitution rate after provided VCF exclusion" if args.exclude_vcf else "maximum non-T>C/non-A>G substitution rate without external variant exclusion",
        "min_mapq": args.min_mapq,
        "min_baseq": args.min_baseq,
        "max_reads": args.max_reads,
    }
    (outdir / "mismatch_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
