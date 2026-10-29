#!/usr/bin/env python3
"""Summarize and plot gffcompare transcript structure classes."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lr_figure_style import apply_style, clean_axis, save_figure


CLASS_LABELS = {
    "=": "Exact reference match",
    "c": "Contained in reference",
    "j": "Novel junction combination",
    "o": "Other same-strand overlap",
    "k": "Contains reference",
    "m": "Retained intron",
    "n": "Reference-intron overlap",
}
CLASS_ORDER = ["=", "c", "j", "o", "k", "m", "n"]
COLORS = ["#3C5488", "#4DBBD5", "#00A087", "#E64B35", "#8491B4", "#F39B7F", "#91D1C2"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tmap", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--font-family", default="DejaVu Sans")
    parser.add_argument("--base-size", type=float, default=10.5)
    args = parser.parse_args()
    apply_style(args.font_family, args.base_size)

    source = Path(args.tmap).resolve()
    outdir = Path(args.outdir)
    tables, figures, metadata = (
        outdir / "01_tables",
        outdir / "02_figures",
        outdir / "05_metadata",
    )
    for directory in (tables, figures, metadata):
        directory.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(source, sep="\t")
    required = {"class_code", "num_exons", "len"}
    if not required.issubset(frame.columns):
        raise SystemExit(f"TMAP is missing columns: {sorted(required - set(frame.columns))}")
    frame["class_code"] = frame["class_code"].astype(str)
    observed = [code for code in CLASS_ORDER if code in set(frame["class_code"])]
    observed += sorted(set(frame["class_code"]) - set(observed))
    counts = frame["class_code"].value_counts()
    summary = pd.DataFrame(
        {
            "class_code": observed,
            "description": [CLASS_LABELS.get(code, "Other gffcompare class") for code in observed],
            "transcripts": [int(counts.get(code, 0)) for code in observed],
        }
    )
    summary["fraction"] = summary["transcripts"] / max(summary["transcripts"].sum(), 1)
    summary.to_csv(tables / "isoform_class_summary.tsv", sep="\t", index=False)
    structure = (
        frame.assign(
            exon_group=np.where(
                pd.to_numeric(frame["num_exons"], errors="coerce").fillna(0) <= 1,
                "Single exon",
                "Multi exon",
            ),
            transcript_length=pd.to_numeric(frame["len"], errors="coerce"),
        )
        .groupby(["class_code", "exon_group"], observed=True)
        .agg(transcripts=("class_code", "size"), median_length=("transcript_length", "median"))
        .reset_index()
    )
    structure.to_csv(tables / "isoform_structure_summary.tsv", sep="\t", index=False)

    display = summary.sort_values("transcripts")
    color_map = dict(zip(CLASS_ORDER, COLORS))
    fig, ax = plt.subplots(figsize=(8.6, 5.0))
    labels = [f"{row.class_code}  {row.description}" for row in display.itertuples(index=False)]
    ax.barh(labels, display["transcripts"], color=[color_map.get(x, "#6E7F80") for x in display["class_code"]])
    ax.set_xlabel("Candidate transcripts")
    ax.set_title("Full-length transcript structural classes", loc="left")
    clean_axis(ax, "x")
    for index, value in enumerate(display["transcripts"]):
        ax.text(value, index, f" {value:,}", va="center", fontsize=8)
    fig.subplots_adjust(left=0.31, right=0.97, bottom=0.13, top=0.90)
    save_figure(fig, figures / "isoform_class_distribution", args.dpi)

    plot = structure.pivot(index="class_code", columns="exon_group", values="transcripts").fillna(0).reindex(observed)
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    left = np.zeros(len(plot))
    for label, color in (("Single exon", "#91D1C2"), ("Multi exon", "#3C5488")):
        values = plot.get(label, pd.Series(0, index=plot.index)).to_numpy()
        ax.bar(plot.index, values, bottom=left, label=label, color=color)
        left += values
    ax.set_xlabel("gffcompare class code")
    ax.set_ylabel("Candidate transcripts")
    ax.set_title("Transcript exon structure", loc="left")
    clean_axis(ax, "y")
    ax.legend(frameon=False, ncol=2, loc="upper right")
    fig.subplots_adjust(left=0.11, right=0.97, bottom=0.14, top=0.90)
    save_figure(fig, figures / "isoform_structure_distribution", args.dpi)

    (metadata / "isoform_class_summary.json").write_text(
        json.dumps(
            {
                "source_tmap": str(source),
                "transcripts": int(len(frame)),
                "class_counts": {str(row.class_code): int(row.transcripts) for row in summary.itertuples(index=False)},
                "interpretation": "Reference-comparison classes describe candidate transcript structures.",
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
