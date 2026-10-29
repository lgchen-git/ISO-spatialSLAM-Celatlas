#!/usr/bin/env python3
"""Summarize SLAM molecule classes by gene."""

from __future__ import annotations

import argparse
import json
import math
import os
from collections import Counter, defaultdict
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
_MATPLOTLIB_CACHE = Path(os.environ.get("TMPDIR", Path.cwd() / "tmp")) / "matplotlib-cache"
_MATPLOTLIB_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_MATPLOTLIB_CACHE.resolve()))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lr_figure_style import apply_style, clean_axis, save_figure

CLASSES = ("new_supported", "no_conversion_detected", "ambiguous")
CLASS_COLORS = {"new_supported": "#E64B35", "no_conversion_detected": "#4DBBD5", "ambiguous": "#B09C85"}


def wilson(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total <= 0:
        return math.nan, math.nan
    p = successes / total
    denominator = 1 + z * z / total
    center = (p + z * z / (2 * total)) / denominator
    half = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denominator
    return max(0.0, center - half), min(1.0, center + half)


def save(fig, stem: Path, dpi: int) -> None:
    save_figure(fig, stem, dpi)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--molecule-calls")
    parser.add_argument("--summary-table")
    parser.add_argument("--class-summary")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--min-informative", type=int, default=100)
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument("--chunk-rows", type=int, default=250000)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--font-family", default="DejaVu Sans")
    parser.add_argument("--base-size", type=float, default=10.5)
    args = parser.parse_args()
    if bool(args.molecule_calls) == bool(args.summary_table):
        raise SystemExit("Provide exactly one of --molecule-calls or --summary-table")
    if args.summary_table and not args.class_summary:
        raise SystemExit("--class-summary is required with --summary-table")
    apply_style(args.font_family, args.base_size)

    source = Path(args.molecule_calls).resolve() if args.molecule_calls else None
    outdir = Path(args.outdir)
    tables = outdir / "01_tables"
    figures = outdir / "02_figures"
    metadata = outdir / "05_metadata"
    for directory in (tables, figures, metadata):
        directory.mkdir(parents=True, exist_ok=True)

    if args.summary_table:
        frame = pd.read_csv(args.summary_table, sep="\t")
        totals = Counter(json.loads(Path(args.class_summary).read_text(encoding="utf-8")).get("class_counts", {}))
    else:
        gene_counts: dict[tuple[str, str], Counter[str]] = defaultdict(Counter)
        totals = Counter()
        for chunk in pd.read_csv(
            source,
            sep="\t",
            compression="infer",
            usecols=["gene_id", "gene_name", "exploratory_class"],
            chunksize=args.chunk_rows,
            dtype=str,
        ):
            chunk["gene_name"] = chunk["gene_name"].fillna("")
            grouped = chunk.groupby(["gene_id", "gene_name", "exploratory_class"], dropna=False).size()
            for (gene_id, gene_name, label), count in grouped.items():
                gene_counts[(str(gene_id), str(gene_name))][str(label)] += int(count)
                totals[str(label)] += int(count)

        rows = []
        for (gene_id, gene_name), counts in gene_counts.items():
            informative = counts["new_supported"] + counts["no_conversion_detected"]
            low, high = wilson(counts["new_supported"], informative)
            rows.append(
                {
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                    "new_supported": counts["new_supported"],
                    "no_conversion_detected": counts["no_conversion_detected"],
                    "ambiguous": counts["ambiguous"],
                    "informative_molecules": informative,
                    "new_supported_fraction": counts["new_supported"] / informative if informative else math.nan,
                    "wilson95_low": low,
                    "wilson95_high": high,
                    "passes_min_informative": informative >= args.min_informative,
                }
            )
        frame = pd.DataFrame(rows).sort_values(["informative_molecules", "gene_id"], ascending=[False, True])
        frame.to_csv(tables / "gene_new_supported_summary.tsv", sep="\t", index=False)

    plot = frame[frame["passes_min_informative"]].head(args.top).sort_values("new_supported_fraction")
    fig_height = max(5.2, 0.27 * len(plot) + 1.6)
    y = np.arange(len(plot))
    values = plot["new_supported_fraction"].to_numpy()
    errors = np.vstack([values - plot["wilson95_low"].to_numpy(), plot["wilson95_high"].to_numpy() - values])
    labels = [name if name else gene for name, gene in zip(plot["gene_name"], plot["gene_id"])]
    high = plot["wilson95_high"].to_numpy()
    sorted_values = np.sort(values)
    use_break = len(values) >= 2 and sorted_values[-1] > 0.45 and sorted_values[-1] - sorted_values[-2] > 0.20
    if use_break:
        split = (sorted_values[-1] + sorted_values[-2]) / 2
        fig, (axis_left, axis_right) = plt.subplots(
            1, 2, figsize=(9.0, fig_height), sharey=True,
            gridspec_kw={"width_ratios": [4.7, 1.25], "wspace": 0.06},
        )
        axes = (axis_left, axis_right)
        left_max = min(split, max(0.05, sorted_values[-2] * 1.16))
        right_min = max(0.0, sorted_values[-1] - max(0.035, sorted_values[-1] * 0.07))
        right_max = min(1.0, max(high[-1] if len(high) else sorted_values[-1], sorted_values[-1]) * 1.04)
        axis_left.set_xlim(0, left_max)
        axis_right.set_xlim(right_min, right_max)
        axis_left.spines["right"].set_visible(False)
        axis_right.spines["left"].set_visible(False)
        axis_right.tick_params(axis="y", left=False, labelleft=False)
        diagonal = 0.012
        kwargs = {"color": "#65717B", "clip_on": False, "linewidth": 1.0}
        axis_left.plot((1 - diagonal, 1 + diagonal), (-diagonal, +diagonal), transform=axis_left.transAxes, **kwargs)
        axis_left.plot((1 - diagonal, 1 + diagonal), (1 - diagonal, 1 + diagonal), transform=axis_left.transAxes, **kwargs)
        axis_right.plot((-diagonal, +diagonal), (-diagonal, +diagonal), transform=axis_right.transAxes, **kwargs)
        axis_right.plot((-diagonal, +diagonal), (1 - diagonal, 1 + diagonal), transform=axis_right.transAxes, **kwargs)
    else:
        fig, axis_left = plt.subplots(figsize=(8.4, fig_height))
        axes = (axis_left,)
        axis_left.set_xlim(0, max(0.05, min(1.0, float(plot["wilson95_high"].max()) * 1.12)) if len(plot) else 1.0)
    for axis in axes:
        axis.errorbar(
            values, y, xerr=errors, fmt="o", color="#E64B35", ecolor="#8491B4",
            capsize=2.5, markersize=4.8, elinewidth=1.2,
        )
        clean_axis(axis, "x")
    axis_left.set_yticks(y, labels)
    fig.supxlabel("New-supported fraction among informative molecules (95% Wilson CI)", fontsize=10)
    fig.suptitle(
        f"Gene-level SLAM signal · ≥{args.min_informative} informative molecules",
        fontsize=args.base_size + 1.0,
    )
    fig.subplots_adjust(left=0.18, right=0.98, bottom=0.10, top=0.91)
    save(fig, figures / "gene_new_supported_fraction", args.dpi)

    total = sum(totals.values())
    fig, ax = plt.subplots(figsize=(6.8, 3.7))
    labels = list(CLASSES)
    values = [totals[label] for label in labels]
    bars = ax.bar(labels, values, color=[CLASS_COLORS[label] for label in labels])
    ax.set_ylabel("Molecules")
    ax.set_title("SLAM molecule-class composition", loc="left")
    clean_axis(ax, "y")
    ax.tick_params(axis="x", rotation=12)
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, value, f"{value:,}\n({value / total:.1%})", ha="center", va="bottom", fontsize=8)
    ax.margins(y=0.18)
    save(fig, figures / "slam_class_composition", args.dpi)

    if args.summary_table:
        return
    (metadata / "gene_slam_summary.json").write_text(
        json.dumps(
            {
                "source_molecule_calls": str(source),
                "minimum_informative_molecules": args.min_informative,
                "genes": int(len(frame)),
                "genes_passing_threshold": int(frame["passes_min_informative"].sum()),
                "class_counts": dict(totals),
                "interpretation": "Gene-level conversion evidence summarized from configured molecule classes.",
            },
            indent=2,
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
