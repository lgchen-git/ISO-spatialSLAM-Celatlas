#!/usr/bin/env python3
"""Apply configured SLAM evidence classes with BH multiple-testing control."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def bh_adjust(values: pd.Series) -> np.ndarray:
    p = pd.to_numeric(values, errors="coerce").fillna(1.0).clip(0.0, 1.0).to_numpy()
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty(len(p), dtype=float)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--molecule-calls", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--min-t-covered", type=int, default=10)
    parser.add_argument("--min-tc-count", type=int, default=2)
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--max-other-mismatch-rate", type=float, default=0.10)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(args.molecule_calls, sep="\t", compression="infer")
    if frame.empty:
        frame["slam_qvalue"] = []
        frame["exploratory_class"] = []
    else:
        frame["slam_qvalue"] = bh_adjust(frame["binomial_p"])
        enough_t = frame["umi_t_covered"].astype(int) >= args.min_t_covered
        clean = pd.to_numeric(frame["umi_other_mismatch_rate"], errors="coerce").fillna(1.0) <= args.max_other_mismatch_rate
        supported = enough_t & clean & (frame["umi_tc_count"].astype(int) >= args.min_tc_count) & (frame["slam_qvalue"] < args.fdr)
        no_conversion = enough_t & (frame["umi_tc_count"].astype(int) == 0)
        frame["exploratory_class"] = np.select([supported, no_conversion], ["new_supported", "no_conversion_detected"], default="ambiguous")
    frame["assignment"] = frame["exploratory_class"]
    output = outdir / f"{args.sample}.slam_molecule_calls.tsv.gz"
    temporary = output.with_name(output.name + ".tmp")
    frame.to_csv(
        temporary,
        sep="\t",
        index=False,
        compression={"method": "gzip", "compresslevel": 1, "mtime": 0},
    )
    temporary.replace(output)

    counts = frame["exploratory_class"].value_counts().to_dict() if not frame.empty else {}
    total = len(frame)
    summary = {
        "sample": args.sample,
        "molecules_total": total,
        "class_counts": {key: int(value) for key, value in counts.items()},
        "new_supported_fraction_all": counts.get("new_supported", 0) / total if total else 0.0,
        "parameters": vars(args),
    }
    (outdir / f"{args.sample}.slam_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    legacy = outdir / f"{args.sample}.threshold_sensitivity.tsv"
    if legacy.exists():
        legacy.unlink()


if __name__ == "__main__":
    main()
