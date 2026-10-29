#!/usr/bin/env python3
"""Plot all 12 single-base substitution channels from mismatch QC output."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

_MATPLOTLIB_CACHE = Path(os.environ.get("TMPDIR", Path.cwd() / "tmp")) / "matplotlib-cache"
_MATPLOTLIB_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_MATPLOTLIB_CACHE.resolve()))

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Patch

from lr_figure_style import apply_style, clean_axis, save_figure

CHANNEL_ORDER = [
    "A>C", "A>G", "A>T",
    "C>A", "C>G", "C>T",
    "G>A", "G>C", "G>T",
    "T>A", "T>C", "T>G",
]
SLAM_CHANNELS = {"A>G", "T>C"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--font-family", default="DejaVu Sans")
    parser.add_argument("--base-size", type=float, default=10.5)
    parser.add_argument("--figure-only", action="store_true")
    args = parser.parse_args()
    apply_style(args.font_family, args.base_size)

    source = Path(args.spectrum).resolve()
    summary_source = Path(args.summary).resolve()
    outdir = Path(args.outdir)
    table_dir = outdir / "01_tables"
    figure_dir = outdir / "02_figures"
    metadata_dir = outdir / "05_metadata"
    for directory in (table_dir, figure_dir, metadata_dir):
        directory.mkdir(parents=True, exist_ok=True)

    frame = pd.read_csv(source, sep="\t")
    required = {"channel", "count", "reference_coverage", "rate"}
    if not required.issubset(frame.columns):
        raise SystemExit(f"Mismatch spectrum is missing columns: {sorted(required - set(frame.columns))}")
    if set(frame["channel"]) != set(CHANNEL_ORDER):
        raise SystemExit("Mismatch spectrum must contain exactly the 12 directed substitution channels")
    frame["channel"] = pd.Categorical(frame["channel"], categories=CHANNEL_ORDER, ordered=True)
    frame = frame.sort_values("channel").reset_index(drop=True)
    frame["rate_percent"] = frame["rate"] * 100.0
    frame["expected_slam_channel"] = frame["channel"].astype(str).isin(SLAM_CHANNELS)
    if not args.figure_only:
        frame.to_csv(table_dir / "base_substitution_comparison.tsv", sep="\t", index=False)

    colors = ["#E64B35" if value else "#4DBBD5" for value in frame["expected_slam_channel"]]
    x = list(range(len(frame)))
    fig, (axis_rate, axis_count) = plt.subplots(
        2, 1, figsize=(9.2, 6.2), sharex=True,
        gridspec_kw={"hspace": 0.10, "height_ratios": [1.15, 1.0]},
    )
    axis_rate.bar(x, frame["rate_percent"], color=colors, width=0.78)
    axis_rate.set_ylabel("Mismatch rate (%)")
    axis_rate.set_title("Single-base substitution comparison", loc="left")
    clean_axis(axis_rate, "y")
    for position, value in enumerate(frame["rate_percent"]):
        axis_rate.text(position, value, f"{value:.3f}", ha="center", va="bottom", fontsize=7, rotation=35)
    axis_rate.margins(y=0.20)

    axis_count.bar(x, frame["count"], color=colors, width=0.78)
    axis_count.set_yscale("log")
    axis_count.set_ylabel("Mismatch count (log scale)")
    axis_count.set_xticks(x, frame["channel"].astype(str), rotation=35, ha="right")
    clean_axis(axis_count, "y")
    fig.legend(
        handles=[
            Patch(facecolor="#E64B35", label="SLAM-associated A>G / T>C"),
            Patch(facecolor="#4DBBD5", label="Other substitutions"),
        ],
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.985),
        ncol=2,
    )
    fig.subplots_adjust(top=0.87, bottom=0.13, left=0.10, right=0.98)
    save_figure(fig, figure_dir / "base_substitution_comparison", args.dpi)

    if args.figure_only:
        return
    mismatch_summary = json.loads(summary_source.read_text(encoding="utf-8"))
    metadata = {
        "source_spectrum": str(source),
        "source_summary": str(summary_source),
        "channels": CHANNEL_ORDER,
        "slam_associated_channels": sorted(SLAM_CHANNELS),
        "rate_definition": "mismatch count divided by aligned reference-base coverage",
        "variant_exclusion_mode": mismatch_summary.get("variant_exclusion_mode", "unknown"),
        "display": "rates on a linear percent scale; counts on a logarithmic scale; values unchanged",
    }
    (metadata_dir / "base_substitution_comparison.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
