#!/usr/bin/env python3
"""Shared Nature/NPG-inspired plotting defaults for LR-Celatlas figures."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


NPG_CATEGORICAL_16 = [
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488",
    "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
    "#7E6148", "#B09C85", "#6F99AD", "#8A7AA9",
    "#C9A227", "#6E7F80", "#D4537E", "#378ADD",
]
NPG_GENERAL = ["#E5E7EB", "#4DBBD5", "#3C5488", "#E64B35"]
NPG_ENRICHMENT = ["#E5E7EB", "#91D1C2", "#00A087", "#DC0000"]
NPG_SIGNED = ["#3C5488", "#F7F7F7", "#E64B35"]
GENERAL_CMAP = LinearSegmentedColormap.from_list("npg_general", NPG_GENERAL)
ENRICHMENT_CMAP = LinearSegmentedColormap.from_list("npg_enrichment", NPG_ENRICHMENT)
SIGNED_CMAP = LinearSegmentedColormap.from_list("npg_signed", NPG_SIGNED)
INK = "#263238"
MUTED = "#65717B"
GRID = "#E5E7EB"


def apply_style(font_family: str = "DejaVu Sans", base_size: float = 10.5) -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [font_family, "DejaVu Sans", "Arial"],
            "font.size": base_size,
            "axes.titlesize": base_size + 1.5,
            "axes.titleweight": "semibold",
            "axes.labelsize": base_size,
            "axes.edgecolor": MUTED,
            "axes.linewidth": 0.8,
            "axes.facecolor": "white",
            "xtick.labelsize": base_size - 1.0,
            "ytick.labelsize": base_size - 1.0,
            "xtick.color": INK,
            "ytick.color": INK,
            "text.color": INK,
            "axes.labelcolor": INK,
            "legend.fontsize": base_size - 1.0,
            "legend.title_fontsize": base_size - 0.5,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.transparent": False,
        }
    )


def named_palette(levels: list[str]) -> dict[str, str]:
    if len(levels) > len(NPG_CATEGORICAL_16):
        raise ValueError("More than 16 categories require an explicit extended palette")
    return dict(zip(levels, NPG_CATEGORICAL_16))


def clean_axis(axis, grid_axis: str | None = None) -> None:
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(False)
    if grid_axis:
        axis.grid(axis=grid_axis, color=GRID, linewidth=0.65, zorder=0)
    axis.set_axisbelow(True)


def save_figure(fig, stem: Path | str, dpi: int, pdf: bool = True) -> None:
    stem = Path(stem)
    suffixes = ("png", "pdf") if pdf else ("png",)
    for suffix in suffixes:
        fig.savefig(
            stem.with_suffix(f".{suffix}"),
            dpi=dpi,
            bbox_inches="tight",
            facecolor="white",
        )
    plt.close(fig)
