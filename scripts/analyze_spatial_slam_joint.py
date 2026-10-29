#!/usr/bin/env python3
"""Join whole-transcriptome spatial clusters with SLAM molecule classes."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lr_figure_style import ENRICHMENT_CMAP, apply_style, clean_axis, named_palette, save_figure


CLASSES = ("new_supported", "no_conversion_detected", "ambiguous")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--molecule-calls")
    parser.add_argument("--clusters", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--bin-size", type=int, required=True)
    parser.add_argument("--min-informative", type=int, required=True)
    parser.add_argument("--gene-plot-top", type=int, required=True)
    parser.add_argument("--chunk-rows", type=int, default=250000)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--font-family", default="DejaVu Sans")
    parser.add_argument("--base-size", type=float, default=10.5)
    parser.add_argument("--render-only", action="store_true")
    return parser.parse_args()


def wilson(success: np.ndarray, total: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    z = 1.959963984540054
    total = total.astype(float)
    proportion = np.divide(success, total, out=np.zeros_like(total), where=total > 0)
    denominator = 1 + z * z / np.maximum(total, 1)
    center = (proportion + z * z / (2 * np.maximum(total, 1))) / denominator
    spread = z * np.sqrt(
        np.divide(proportion * (1 - proportion), np.maximum(total, 1))
        + z * z / (4 * np.maximum(total, 1) ** 2)
    ) / denominator
    return np.clip(center - spread, 0, 1), np.clip(center + spread, 0, 1)


def build_tables(args: argparse.Namespace, clusters: pd.DataFrame, tables: Path) -> None:
    if not args.molecule_calls:
        raise SystemExit("--molecule-calls is required unless --render-only is used")
    bin_to_cluster = clusters.set_index("bin_id")["cluster"].astype(str).to_dict()
    cluster_parts: list[pd.DataFrame] = []
    gene_parts: list[pd.DataFrame] = []
    input_molecules = assigned_molecules = 0
    usecols = ["x", "y", "gene_id", "gene_name", "exploratory_class"]
    for chunk in pd.read_csv(
        args.molecule_calls,
        sep="\t",
        compression="infer",
        usecols=usecols,
        chunksize=args.chunk_rows,
    ):
        input_molecules += len(chunk)
        x = pd.to_numeric(chunk["x"], errors="coerce")
        y = pd.to_numeric(chunk["y"], errors="coerce")
        valid = x.notna() & y.notna()
        chunk = chunk.loc[valid].copy()
        chunk["bin_id"] = (
            ((x.loc[valid].astype(int) // args.bin_size) * args.bin_size + args.bin_size // 2).astype(str)
            + "_"
            + ((y.loc[valid].astype(int) // args.bin_size) * args.bin_size + args.bin_size // 2).astype(str)
        )
        chunk["cluster"] = chunk["bin_id"].map(bin_to_cluster)
        chunk = chunk[chunk["cluster"].notna() & chunk["exploratory_class"].isin(CLASSES)]
        assigned_molecules += len(chunk)
        cluster_parts.append(
            chunk.groupby(["cluster", "exploratory_class"], observed=True).size().rename("molecules").reset_index()
        )
        gene_parts.append(
            chunk.groupby(["cluster", "gene_id", "gene_name", "exploratory_class"], observed=True, dropna=False)
            .size()
            .rename("molecules")
            .reset_index()
        )
    cluster_long = pd.concat(cluster_parts, ignore_index=True).groupby(
        ["cluster", "exploratory_class"], as_index=False, observed=True
    )["molecules"].sum()
    cluster_wide = cluster_long.pivot(
        index="cluster", columns="exploratory_class", values="molecules"
    ).fillna(0).reset_index()
    for label in CLASSES:
        if label not in cluster_wide:
            cluster_wide[label] = 0
        cluster_wide[label] = cluster_wide[label].astype(int)
    cluster_wide["informative_molecules"] = (
        cluster_wide["new_supported"] + cluster_wide["no_conversion_detected"]
    )
    cluster_wide["new_supported_fraction"] = np.divide(
        cluster_wide["new_supported"],
        cluster_wide["informative_molecules"],
        out=np.zeros(len(cluster_wide), dtype=float),
        where=cluster_wide["informative_molecules"] > 0,
    )
    low, high = wilson(
        cluster_wide["new_supported"].to_numpy(),
        cluster_wide["informative_molecules"].to_numpy(),
    )
    cluster_wide["wilson95_low"], cluster_wide["wilson95_high"] = low, high
    cluster_wide.sort_values("cluster", key=lambda x: x.astype(int)).to_csv(
        tables / "cluster_slam_summary.tsv", sep="\t", index=False
    )

    gene_long = pd.concat(gene_parts, ignore_index=True).groupby(
        ["cluster", "gene_id", "gene_name", "exploratory_class"],
        as_index=False,
        observed=True,
        dropna=False,
    )["molecules"].sum()
    gene_wide = gene_long.pivot_table(
        index=["cluster", "gene_id", "gene_name"],
        columns="exploratory_class",
        values="molecules",
        aggfunc="sum",
        fill_value=0,
    ).reset_index()
    for label in CLASSES:
        if label not in gene_wide:
            gene_wide[label] = 0
        gene_wide[label] = gene_wide[label].astype(int)
    gene_wide["total_molecules"] = gene_wide[list(CLASSES)].sum(axis=1)
    gene_wide["informative_molecules"] = (
        gene_wide["new_supported"] + gene_wide["no_conversion_detected"]
    )
    gene_wide["new_supported_fraction"] = np.divide(
        gene_wide["new_supported"],
        gene_wide["informative_molecules"],
        out=np.zeros(len(gene_wide), dtype=float),
        where=gene_wide["informative_molecules"] > 0,
    )
    gene_wide["passes_min_informative"] = gene_wide["informative_molecules"] >= args.min_informative
    low, high = wilson(
        gene_wide["new_supported"].to_numpy(),
        gene_wide["informative_molecules"].to_numpy(),
    )
    gene_wide["wilson95_low"], gene_wide["wilson95_high"] = low, high
    gene_wide.sort_values(
        ["cluster", "informative_molecules"], ascending=[True, False]
    ).to_csv(tables / "cluster_gene_slam.tsv", sep="\t", index=False)
    (tables.parent / "05_metadata" / "slam_joint_summary.json").write_text(
        json.dumps(
            {
                "sample": args.sample,
                "input_molecules": int(input_molecules),
                "molecules_joined_to_whole_transcriptome_bins": int(assigned_molecules),
                "clusters": int(cluster_wide["cluster"].nunique()),
                "genes": int(gene_wide["gene_id"].nunique()),
                "min_informative_molecules": args.min_informative,
                "joint_level": "whole-transcriptome spatial cluster and gene",
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def plot_results(args: argparse.Namespace, tables: Path, figures: Path) -> None:
    cluster = pd.read_csv(tables / "cluster_slam_summary.tsv", sep="\t")
    cluster["cluster"] = cluster["cluster"].astype(str)
    cluster = cluster.sort_values("cluster", key=lambda x: x.astype(int))
    palette = named_palette(list(cluster["cluster"]))
    positions = np.arange(len(cluster))

    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    totals = cluster[list(CLASSES)].sum(axis=1).replace(0, np.nan)
    bottom = np.zeros(len(cluster))
    labels = [
        ("new_supported", "New-supported", "#E64B35"),
        ("no_conversion_detected", "No conversion detected", "#3C5488"),
        ("ambiguous", "Ambiguous", "#B09C85"),
    ]
    for column, label, color in labels:
        values = (cluster[column] / totals).fillna(0).to_numpy()
        ax.bar(positions, values, bottom=bottom, color=color, label=label)
        bottom += values
    ax.set_xticks(positions, [f"C{x}" for x in cluster["cluster"]])
    ax.set_ylabel("Molecule fraction")
    ax.set_title("SLAM classes within whole-transcriptome clusters", loc="left")
    ax.set_ylim(0, 1)
    clean_axis(ax, "y")
    ax.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.17))
    fig.subplots_adjust(left=0.10, right=0.98, bottom=0.25, top=0.90)
    save_figure(fig, figures / "cluster_slam_composition", args.dpi)

    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    fraction = cluster["new_supported_fraction"].to_numpy()
    errors = np.vstack(
        [fraction - cluster["wilson95_low"].to_numpy(), cluster["wilson95_high"].to_numpy() - fraction]
    )
    ax.errorbar(
        positions,
        fraction,
        yerr=errors,
        fmt="none",
        ecolor="#6E7F80",
        elinewidth=1,
        capsize=2.5,
        zorder=1,
    )
    ax.scatter(positions, fraction, c=[palette[x] for x in cluster["cluster"]], s=42, zorder=2)
    ax.set_xticks(positions, [f"C{x}" for x in cluster["cluster"]])
    ax.set_ylabel("New-supported fraction")
    ax.set_title("New-supported fraction by expression cluster", loc="left")
    clean_axis(ax, "y")
    fig.subplots_adjust(left=0.11, right=0.98, bottom=0.14, top=0.90)
    save_figure(fig, figures / "cluster_new_supported_fraction", args.dpi)

    gene = pd.read_csv(tables / "cluster_gene_slam.tsv", sep="\t")
    gene = gene[gene["passes_min_informative"].astype(str).str.lower().eq("true")].copy()
    totals = (
        gene.groupby(["gene_id", "gene_name"], observed=True)["informative_molecules"]
        .sum()
        .sort_values(ascending=False)
        .head(args.gene_plot_top)
    )
    selected_ids = set(totals.index.get_level_values("gene_id"))
    gene = gene[gene["gene_id"].isin(selected_ids)]
    if gene.empty:
        return
    gene_order = list(dict.fromkeys(totals.index.get_level_values("gene_id")))
    gene_names = (
        gene.drop_duplicates("gene_id").set_index("gene_id")["gene_name"].fillna("").astype(str).to_dict()
    )
    cluster_order = list(cluster["cluster"])
    x_map = {gene_id: index for index, gene_id in enumerate(gene_order)}
    y_map = {value: index for index, value in enumerate(cluster_order)}
    max_info = max(float(gene["informative_molecules"].max()), 1.0)
    fig, ax = plt.subplots(
        figsize=(max(9.2, 0.43 * len(gene_order) + 2.8), max(5.2, 0.38 * len(cluster_order) + 2.0))
    )
    points = ax.scatter(
        gene["gene_id"].map(x_map),
        gene["cluster"].astype(str).map(y_map),
        s=28 + 210 * np.sqrt(gene["informative_molecules"] / max_info),
        c=gene["new_supported_fraction"],
        cmap=ENRICHMENT_CMAP,
        vmin=0,
        vmax=max(float(gene["new_supported_fraction"].quantile(0.98)), 0.01),
        linewidths=0,
    )
    ax.set_xticks(
        range(len(gene_order)),
        [gene_names.get(gene_id) or gene_id for gene_id in gene_order],
        rotation=50,
        ha="right",
    )
    ax.set_yticks(range(len(cluster_order)), [f"C{x}" for x in cluster_order])
    ax.set_xlabel("Gene")
    ax.set_ylabel("Whole-transcriptome cluster")
    ax.set_title("Gene expression support and SLAM fraction", loc="left")
    ax.invert_yaxis()
    clean_axis(ax)
    colorbar = fig.colorbar(points, ax=ax, fraction=0.035, pad=0.02)
    colorbar.set_label("New-supported fraction")
    fig.subplots_adjust(left=0.09, right=0.92, bottom=0.28, top=0.90)
    save_figure(fig, figures / "cluster_gene_slam_dotplot", args.dpi)


def main() -> None:
    args = parse_args()
    if args.bin_size <= 0 or args.min_informative <= 0 or args.gene_plot_top <= 0:
        raise SystemExit("Bin size and joint-analysis thresholds must be positive")
    apply_style(args.font_family, args.base_size)
    outdir = Path(args.outdir)
    tables, figures, metadata = outdir / "01_tables", outdir / "02_figures", outdir / "05_metadata"
    for directory in (tables, figures, metadata):
        directory.mkdir(parents=True, exist_ok=True)
    clusters = pd.read_csv(args.clusters, sep="\t")
    required = {"bin_id", "cluster"}
    if not required.issubset(clusters.columns):
        raise SystemExit(f"Cluster table is missing: {sorted(required - set(clusters.columns))}")
    if not args.render_only:
        build_tables(args, clusters, tables)
    plot_results(args, tables, figures)


if __name__ == "__main__":
    main()
