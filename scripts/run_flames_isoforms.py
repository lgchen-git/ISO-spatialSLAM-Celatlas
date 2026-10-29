#!/usr/bin/env python3
"""Run FLAMES full-length isoform discovery on an existing long-read BAM."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--flames-python-dir", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--min-support", type=int, required=True)
    parser.add_argument("--min-count-fraction", type=float, required=True)
    parser.add_argument("--min-support-fraction", type=float, required=True)
    parser.add_argument("--strand-specific", type=int, choices=(0, 1, 2), required=True)
    parser.add_argument("--downsample-ratio", type=float, required=True)
    parser.add_argument("--threads", type=int, required=True)
    args = parser.parse_args()
    if (
        args.min_support <= 0
        or args.threads <= 0
        or not 0 <= args.min_count_fraction <= 1
        or not 0 <= args.min_support_fraction <= 1
        or not 0 < args.downsample_ratio <= 1
    ):
        raise SystemExit("Invalid FLAMES parameters")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, args.flames_python_dir)
    from find_isoform import find_isoform

    template = Path(args.flames_python_dir).parent / "extdata" / "template_config.json"
    config = json.loads(template.read_text(encoding="utf-8"))
    config["isoform_parameters"].update(
        {
            "generate_raw_isoform": True,
            "min_sup_cnt": args.min_support,
            "min_cnt_pct": args.min_count_fraction,
            "min_sup_pct": args.min_support_fraction,
            "strand_specific": args.strand_specific,
            "downsample_ratio": args.downsample_ratio,
        }
    )
    config["pipeline_parameters"]["bambu_isoform_identification"] = False
    (outdir / "flames_config.json").write_text(
        json.dumps(config, indent=2) + "\n", encoding="utf-8"
    )
    find_isoform(
        args.gtf,
        args.bam,
        str(outdir / "isoforms.gff3"),
        str(outdir / "tss_tes_stats.tsv"),
        args.genome,
        str(outdir / "isoforms.fa"),
        args.threads,
        config,
        str(outdir / "raw_isoforms.gff3"),
    )


if __name__ == "__main__":
    main()
