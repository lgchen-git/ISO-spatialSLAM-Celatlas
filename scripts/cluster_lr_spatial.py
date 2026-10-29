#!/usr/bin/env python3
"""Cluster spatial bins from Celatlas long-read molecule calls."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from collections import Counter, defaultdict
from pathlib import Path

_MATPLOTLIB_CACHE = Path(os.environ.get("TMPDIR", Path.cwd() / "tmp")) / "matplotlib-cache"
_MATPLOTLIB_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_MATPLOTLIB_CACHE.resolve()))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import PercentFormatter
from scipy import sparse
from scipy.stats import pearsonr

from lr_figure_style import (
    ENRICHMENT_CMAP,
    GENERAL_CMAP,
    SIGNED_CMAP,
    apply_style,
    clean_axis,
    named_palette,
    save_figure,
)

CLASSES = ("new_supported", "no_conversion_detected", "ambiguous")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--molecule-calls")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--sample")
    parser.add_argument("--bin-size", type=int, default=50)
    parser.add_argument("--coordinate-um-per-unit", type=float, default=1.0)
    parser.add_argument("--chip-bounds", type=int, default=9357)
    parser.add_argument("--min-bin-molecules", type=int, default=20)
    parser.add_argument("--min-gene-bins", type=int, default=3)
    parser.add_argument("--hvg", type=int, default=2000)
    parser.add_argument("--pca-components", type=int, default=30)
    parser.add_argument("--neighbors", type=int, default=24)
    parser.add_argument("--regress-log-depth", default="true")
    parser.add_argument("--spatial-graph-weight", type=float, default=0.0)
    parser.add_argument("--min-cluster-bins", type=int, default=1)
    parser.add_argument("--leiden-resolution", type=float, default=0.5)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--marker-top", type=int, default=20)
    parser.add_argument("--marker-display-top", type=int, default=3)
    parser.add_argument("--marker-min-logfoldchange", type=float, default=0.25)
    parser.add_argument("--marker-max-fdr", type=float, default=0.05)
    parser.add_argument("--chunk-rows", type=int, default=250000)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--font-family", default="DejaVu Sans")
    parser.add_argument("--base-size", type=float, default=10.5)
    parser.add_argument("--render-only", action="store_true")
    return parser.parse_args()


def truth(value: str) -> bool:
    return value.lower() in {"1", "true", "yes", "on"}


def build_matrix(path: Path, bin_size: int, chunk_rows: int):
    bin_index: dict[str, int] = {}
    bin_coords: list[tuple[float, float]] = []
    gene_index: dict[str, int] = {}
    gene_names: list[str] = []
    row_parts: list[np.ndarray] = []
    col_parts: list[np.ndarray] = []
    value_parts: list[np.ndarray] = []
    class_counts: dict[str, Counter[str]] = defaultdict(Counter)
    qc = Counter()
    usecols = ["x", "y", "gene_id", "gene_name", "exploratory_class"]
    for chunk_number, chunk in enumerate(
        pd.read_csv(path, sep="\t", compression="infer", usecols=usecols, chunksize=chunk_rows), start=1
    ):
        qc["input_molecules"] += len(chunk)
        valid = (
            pd.to_numeric(chunk["x"], errors="coerce").notna()
            & pd.to_numeric(chunk["y"], errors="coerce").notna()
            & chunk["gene_id"].notna()
            & (chunk["gene_id"].astype(str) != "")
        )
        qc["invalid_coordinate_or_gene"] += int((~valid).sum())
        chunk = chunk.loc[valid].copy()
        chunk["x"] = chunk["x"].astype(int)
        chunk["y"] = chunk["y"].astype(int)
        chunk = chunk[(chunk["x"] >= 0) & (chunk["y"] >= 0)]
        qc["valid_molecules"] += len(chunk)
        chunk["bin_x"] = (chunk["x"] // bin_size) * bin_size + bin_size / 2.0
        chunk["bin_y"] = (chunk["y"] // bin_size) * bin_size + bin_size / 2.0
        chunk["bin_id"] = chunk["bin_x"].astype(int).astype(str) + "_" + chunk["bin_y"].astype(int).astype(str)

        grouped = (
            chunk.groupby(["bin_id", "bin_x", "bin_y", "gene_id", "gene_name"], observed=True, dropna=False)
            .size()
            .reset_index(name="molecules")
        )
        rows = np.empty(len(grouped), dtype=np.int32)
        cols = np.empty(len(grouped), dtype=np.int32)
        for idx, item in enumerate(grouped.itertuples(index=False)):
            if item.bin_id not in bin_index:
                bin_index[item.bin_id] = len(bin_coords)
                bin_coords.append((float(item.bin_x), float(item.bin_y)))
            gene_id = str(item.gene_id)
            if gene_id not in gene_index:
                gene_index[gene_id] = len(gene_names)
                gene_name = "" if pd.isna(item.gene_name) else str(item.gene_name)
                gene_names.append(gene_name)
            rows[idx] = bin_index[item.bin_id]
            cols[idx] = gene_index[gene_id]
        row_parts.append(rows)
        col_parts.append(cols)
        value_parts.append(grouped["molecules"].to_numpy(dtype=np.int32))

        per_class = chunk.groupby(["bin_id", "exploratory_class"], observed=True).size()
        for (bin_id, label), count in per_class.items():
            class_counts[str(bin_id)][str(label)] += int(count)
        print(f"[lr-spatial] chunk={chunk_number} molecules={qc['input_molecules']}", flush=True)

    if not row_parts:
        raise SystemExit("No valid long-read molecules were found")
    matrix = sparse.coo_matrix(
        (np.concatenate(value_parts), (np.concatenate(row_parts), np.concatenate(col_parts))),
        shape=(len(bin_coords), len(gene_names)),
        dtype=np.float32,
    ).tocsr()
    matrix.sum_duplicates()
    gene_ids = [""] * len(gene_index)
    for gene_id, index in gene_index.items():
        gene_ids[index] = gene_id
    bin_ids = [""] * len(bin_index)
    for bin_id, index in bin_index.items():
        bin_ids[index] = bin_id
    return matrix, bin_ids, bin_coords, gene_ids, gene_names, class_counts, qc


def plot_categorical(
    frame: pd.DataFrame,
    x: str,
    y: str,
    output: Path,
    title: str,
    dpi: int,
    invert_y: bool,
    show_labels: bool,
) -> dict[str, str]:
    levels = sorted(frame["cluster"].astype(str).unique(), key=int)
    palette = named_palette(levels)
    fig, ax = plt.subplots(figsize=(7.4, 6.2))
    for level in levels:
        subset = frame[frame["cluster"].astype(str) == level]
        size = 4.4 if invert_y else (1.2 if len(frame) >= 20000 else 1.8)
        marker = "s" if invert_y else "o"
        alpha = 0.86 if invert_y else 0.74
        ax.scatter(subset[x], subset[y], s=size, marker=marker, alpha=alpha, linewidths=0, color=palette[level], label=f"C{level}", rasterized=True)
    if show_labels:
        centers = frame.assign(_cluster=frame["cluster"].astype(str)).groupby("_cluster")[[x, y]].median()
        for level, center in centers.iterrows():
            ax.text(center[x], center[y], f"C{level}", ha="center", va="center", fontsize=7.2, fontweight="bold", color="#222222", bbox={"facecolor": "white", "alpha": 0.82, "edgecolor": "none", "pad": 1.2}, zorder=5)
    ax.set_title(title, loc="left")
    ax.set_xlabel("X (µm; right)" if invert_y else "UMAP1")
    ax.set_ylabel("Y (µm; down)" if invert_y else "UMAP2")
    if invert_y:
        ax.invert_yaxis()
    ax.set_aspect("equal", adjustable="datalim")
    clean_axis(ax)
    ax.legend(
        frameon=False,
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        markerscale=2.2,
        ncol=2,
        columnspacing=0.8,
        handletextpad=0.3,
        borderaxespad=0,
    )
    fig.subplots_adjust(left=0.11, right=0.78, bottom=0.11, top=0.92)
    save_figure(fig, output, dpi)
    return palette


def plot_continuous(
    frame: pd.DataFrame,
    value: str,
    output: Path,
    title: str,
    colorbar_label: str,
    dpi: int,
    cmap,
    percent: bool = False,
) -> None:
    plot = frame[np.isfinite(pd.to_numeric(frame[value], errors="coerce"))].copy()
    plot[value] = pd.to_numeric(plot[value])
    plot = plot.sort_values(value)
    vmin, vmax = np.percentile(plot[value], [2, 98])
    fig, ax = plt.subplots(figsize=(7.4, 6.2))
    points = ax.scatter(
        plot["bin_x"], plot["bin_y"], c=plot[value], cmap=cmap, vmin=vmin, vmax=vmax,
        s=5.2, marker="s", alpha=0.88, linewidths=0, rasterized=True,
    )
    colorbar = fig.colorbar(points, ax=ax, label=colorbar_label)
    if percent:
        colorbar.ax.yaxis.set_major_formatter(PercentFormatter(xmax=1.0))
    ax.set_title(title, loc="left")
    ax.set_xlabel("X (µm; right)")
    ax.set_ylabel("Y (µm; down)")
    ax.invert_yaxis()
    ax.set_aspect("equal", adjustable="datalim")
    clean_axis(ax)
    fig.subplots_adjust(left=0.11, right=0.90, bottom=0.11, top=0.92)
    save_figure(fig, output, dpi)


def marker_table(sc, adata, output: Path, top: int, min_logfc: float, max_fdr: float) -> pd.DataFrame:
    sc.tl.rank_genes_groups(
        adata,
        groupby="cluster",
        method="t-test_overestim_var",
        reference="rest",
        n_genes=min(top, adata.n_vars),
        pts=True,
        use_raw=False,
    )
    rows = []
    for cluster in adata.obs["cluster"].cat.categories:
        frame = sc.get.rank_genes_groups_df(adata, group=cluster)
        for rank, item in enumerate(frame.head(top).itertuples(index=False), start=1):
            significant = bool(float(item.logfoldchanges) >= min_logfc and float(item.pvals_adj) <= max_fdr)
            rows.append(
                {
                    "cluster": str(cluster),
                    "rank": rank,
                    "gene_id": str(item.names),
                    "gene_name": str(adata.var.loc[str(item.names), "gene_name"]),
                    "score": float(item.scores),
                    "logfoldchange": float(item.logfoldchanges),
                    "pvalue_adj": float(item.pvals_adj),
                    "passes_display_filter": significant,
                }
            )
    result = pd.DataFrame(rows)
    result.to_csv(output, sep="\t", index=False)
    return result


def marker_figures(
    adata,
    markers: pd.DataFrame,
    figure_dir: Path,
    display_top: int,
    dpi: int,
) -> None:
    selected = (
        markers[markers["passes_display_filter"]]
        .sort_values(["cluster", "rank"])
        .groupby("cluster", observed=True, group_keys=False)
        .head(display_top)
    )
    if selected.empty:
        return
    gene_ids = list(dict.fromkeys(selected["gene_id"].astype(str)))
    gene_labels = {
        str(row.gene_id): (str(row.gene_name) if str(row.gene_name) not in {"", "nan"} else str(row.gene_id))
        for row in selected.itertuples(index=False)
    }
    levels = [str(value) for value in adata.obs["cluster"].cat.categories]
    expression = adata[:, gene_ids].X
    expression = expression.toarray() if sparse.issparse(expression) else np.asarray(expression)
    means = np.zeros((len(levels), len(gene_ids)))
    fractions = np.zeros_like(means)
    labels = adata.obs["cluster"].astype(str).to_numpy()
    for index, level in enumerate(levels):
        subset = expression[labels == level]
        means[index] = subset.mean(axis=0)
        fractions[index] = (subset > 0).mean(axis=0)
    zscores = (means - means.mean(axis=0)) / np.where(means.std(axis=0) == 0, 1, means.std(axis=0))

    fig, ax = plt.subplots(
        figsize=(max(9.0, 0.40 * len(gene_ids) + 2.8), max(5.0, 0.30 * len(levels) + 1.8))
    )
    x_grid, y_grid = np.meshgrid(np.arange(len(gene_ids)), np.arange(len(levels)))
    points = ax.scatter(
        x_grid.ravel(), y_grid.ravel(),
        s=(20 + 190 * fractions).ravel(),
        c=zscores.ravel(),
        cmap=SIGNED_CMAP,
        vmin=-2,
        vmax=2,
        linewidths=0,
    )
    ax.set_xticks(np.arange(len(gene_ids)), [gene_labels[gene] for gene in gene_ids], rotation=50, ha="right")
    ax.set_yticks(np.arange(len(levels)), [f"C{level}" for level in levels])
    ax.set_title("Data-driven cluster markers", loc="left")
    ax.set_xlabel("")
    ax.set_ylabel("Expression cluster")
    ax.invert_yaxis()
    clean_axis(ax)
    colorbar = fig.colorbar(points, ax=ax, pad=0.015, fraction=0.035)
    colorbar.set_label("Mean-expression z score")
    size_handles = [
        ax.scatter([], [], s=20 + 190 * fraction, color="#6F99AD", linewidths=0, label=f"{fraction:.0%}")
        for fraction in (0.25, 0.50, 0.75)
    ]
    ax.legend(
        handles=size_handles, title="Detected bins", frameon=False, ncol=3,
        loc="upper center", bbox_to_anchor=(0.5, -0.23), columnspacing=1.1, handletextpad=0.3,
    )
    fig.subplots_adjust(left=0.08, right=0.93, bottom=0.32, top=0.92)
    save_figure(fig, figure_dir / "cluster_marker_dotplot", dpi)

    fig, ax = plt.subplots(figsize=(max(8.8, 0.38 * len(gene_ids) + 2.7), max(4.7, 0.28 * len(levels) + 1.6)))
    image = ax.imshow(zscores, aspect="auto", cmap=SIGNED_CMAP, vmin=-2, vmax=2)
    ax.set_xticks(np.arange(len(gene_ids)), [gene_labels[gene] for gene in gene_ids], rotation=50, ha="right")
    ax.set_yticks(np.arange(len(levels)), [f"C{level}" for level in levels])
    ax.set_title("Cluster mean expression of selected markers", loc="left")
    fig.colorbar(image, ax=ax, label="Mean-expression z score", shrink=0.8)
    fig.subplots_adjust(left=0.08, right=0.94, bottom=0.24, top=0.92)
    save_figure(fig, figure_dir / "cluster_marker_heatmap", dpi)


def symmetric_normalize(matrix):
    matrix = matrix.tocsr().astype(np.float64)
    degree = np.asarray(matrix.sum(axis=1)).ravel()
    inverse = np.zeros_like(degree)
    inverse[degree > 0] = 1.0 / np.sqrt(degree[degree > 0])
    scaling = sparse.diags(inverse)
    return (scaling @ matrix @ scaling).tocsr()


def make_spatial_graph(coords: list[tuple[float, float]], step: int):
    lookup = {(int(x), int(y)): index for index, (x, y) in enumerate(coords)}
    rows: list[int] = []
    columns: list[int] = []
    for (x, y), index in lookup.items():
        for neighbor in ((x + step, y), (x, y + step)):
            if neighbor in lookup:
                other = lookup[neighbor]
                rows.extend((index, other))
                columns.extend((other, index))
    values = np.ones(len(rows), dtype=np.float64)
    return sparse.coo_matrix((values, (rows, columns)), shape=(len(coords), len(coords))).tocsr()


def merge_small_clusters(graph, labels: pd.Series, minimum_size: int) -> pd.Series:
    merged = labels.astype(str).to_numpy(copy=True)
    while True:
        counts = pd.Series(merged).value_counts()
        small = counts[counts < minimum_size]
        if small.empty or len(counts) <= 1:
            return pd.Series(merged, index=labels.index, dtype=str)
        source = str(small.index[0])
        mask = merged == source
        connection = np.asarray(graph[mask].sum(axis=0)).ravel()
        scores = {
            target: float(connection[merged == target].sum())
            for target in np.unique(merged)
            if target != source
        }
        target = max(scores, key=scores.get) if max(scores.values(), default=0.0) > 0 else str(counts.drop(source).idxmax())
        merged[mask] = target


def weighted_edge_agreement(graph, labels: pd.Series) -> float:
    values = labels.astype(str).to_numpy()
    upper = sparse.triu(graph, k=1).tocoo()
    if upper.nnz == 0 or float(upper.data.sum()) == 0:
        return math.nan
    return float(upper.data[values[upper.row] == values[upper.col]].sum() / upper.data.sum())


def depth_eta_squared(depth: np.ndarray, labels: pd.Series) -> float:
    values = labels.astype(str).to_numpy()
    overall = float(depth.mean())
    between = sum(
        int((values == label).sum()) * (float(depth[values == label].mean()) - overall) ** 2
        for label in np.unique(values)
    )
    total = float(((depth - overall) ** 2).sum())
    return between / total if total else math.nan


def spatial_coherence(frame: pd.DataFrame, step: int) -> tuple[pd.DataFrame, dict[str, float]]:
    lookup = {(int(row.bin_x), int(row.bin_y)): (str(row.cluster), str(row.bin_id)) for row in frame.itertuples(index=False)}
    edges: list[tuple[str, str, str, str]] = []
    for (x, y), (cluster, bin_id) in lookup.items():
        for neighbor in ((x + step, y), (x, y + step)):
            if neighbor in lookup:
                other_cluster, other_bin = lookup[neighbor]
                edges.append((bin_id, other_bin, cluster, other_cluster))
    incident = Counter()
    same = Counter()
    adjacency: dict[str, set[str]] = defaultdict(set)
    for left, right, left_cluster, right_cluster in edges:
        incident[left_cluster] += 1
        incident[right_cluster] += 1
        if left_cluster == right_cluster:
            same[left_cluster] += 2
            adjacency[left].add(right)
            adjacency[right].add(left)
    rows = []
    for cluster, subset in frame.groupby(frame["cluster"].astype(str), observed=True):
        nodes = set(subset["bin_id"].astype(str))
        unseen = set(nodes)
        sizes = []
        while unseen:
            seed = unseen.pop()
            stack = [seed]
            size = 0
            while stack:
                node = stack.pop()
                size += 1
                neighbors = adjacency.get(node, set()) & unseen
                unseen.difference_update(neighbors)
                stack.extend(neighbors)
            sizes.append(size)
        largest = max(sizes) if sizes else 0
        rows.append(
            {
                "cluster": cluster,
                "bins": len(nodes),
                "same_label_neighbor_fraction": same[cluster] / incident[cluster] if incident[cluster] else math.nan,
                "connected_components": len(sizes),
                "largest_component_bins": largest,
                "largest_component_fraction": largest / len(nodes) if nodes else math.nan,
            }
        )
    result = pd.DataFrame(rows).sort_values("cluster", key=lambda values: values.astype(int))
    observed = sum(left == right for _, _, left, right in edges) / len(edges) if edges else math.nan
    proportions = frame["cluster"].astype(str).value_counts(normalize=True)
    expected = float((proportions * proportions).sum())
    return result, {"neighbor_edges": len(edges), "same_label_edge_fraction": observed, "expected_same_label_fraction": expected}


def qc_figures(frame: pd.DataFrame, coherence: pd.DataFrame, figure_dir: Path, dpi: int) -> None:
    levels = sorted(frame["cluster"].astype(str).unique(), key=int)
    palette = named_palette(levels)
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0))
    data = [np.log1p(frame.loc[frame["cluster"].astype(str) == level, "total_lr_molecules"]) for level in levels]
    box = axes[0].boxplot(data, patch_artist=True, showfliers=False, tick_labels=[f"C{x}" for x in levels])
    for patch, level in zip(box["boxes"], levels):
        patch.set_facecolor(palette[level])
        patch.set_alpha(0.72)
    axes[0].set_ylabel("log1p long-read molecules per bin")
    data = [frame.loc[frame["cluster"].astype(str) == level, "detected_genes"] for level in levels]
    box = axes[1].boxplot(data, patch_artist=True, showfliers=False, tick_labels=[f"C{x}" for x in levels])
    for patch, level in zip(box["boxes"], levels):
        patch.set_facecolor(palette[level])
        patch.set_alpha(0.72)
    axes[1].set_ylabel("Detected genes per bin")
    for axis in axes:
        axis.set_xlabel("Expression cluster")
        clean_axis(axis, "y")
    fig.suptitle("Sequencing-depth distribution by cluster", x=0.06, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save_figure(fig, figure_dir / "cluster_depth_qc", dpi)

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.9))
    positions = np.arange(len(coherence))
    colors = [palette[str(value)] for value in coherence["cluster"]]
    axes[0].bar(positions, coherence["same_label_neighbor_fraction"], color=colors)
    axes[1].bar(positions, coherence["largest_component_fraction"], color=colors)
    for axis, label in zip(axes, ("Same-label neighbor fraction", "Largest connected-component fraction")):
        axis.set_xticks(positions, [f"C{x}" for x in coherence["cluster"]])
        axis.set_ylim(0, 1)
        axis.set_ylabel(label)
        clean_axis(axis, "y")
    fig.suptitle("Spatial continuity of expression clusters", x=0.06, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save_figure(fig, figure_dir / "spatial_coherence", dpi)


def render_existing(outdir: Path, args: argparse.Namespace) -> None:
    import anndata as ad

    table_dir = outdir / "01_tables"
    figure_dir = outdir / "02_figures"
    tool_dir = outdir / "04_tool_output" / "scanpy"
    metadata_dir = outdir / "05_metadata"
    h5ad_path = tool_dir / "long_read_spatial_bin50.h5ad"
    required = [
        h5ad_path,
        table_dir / "long_read_cluster_markers.tsv",
        table_dir / "long_read_spatial_coherence.tsv",
        metadata_dir / "long_read_spatial_clustering.json",
    ]
    missing = [str(path) for path in required if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise SystemExit("Render-only inputs are missing:\n" + "\n".join(missing))
    adata = ad.read_h5ad(h5ad_path)
    frame = adata.obs.reset_index()
    frame["UMAP1"] = adata.obsm["X_umap"][:, 0]
    frame["UMAP2"] = adata.obsm["X_umap"][:, 1]
    counts = adata.layers["counts"]
    frame["detected_genes"] = np.asarray((counts > 0).sum(axis=1)).ravel().astype(int)
    if not (table_dir / "whole_transcriptome_gene_summary.tsv").is_file():
        gene_names = adata.var.get("gene_name", pd.Series(adata.var_names, index=adata.var_names))
        pd.DataFrame(
            {
                "gene_id": adata.var_names.astype(str),
                "gene_name": gene_names.astype(str).to_numpy(),
                "total_molecules": np.asarray(counts.sum(axis=0)).ravel().astype(int),
                "detected_bins": np.asarray((counts > 0).sum(axis=0)).ravel().astype(int),
            }
        ).sort_values(["total_molecules", "gene_id"], ascending=[False, True]).to_csv(
            table_dir / "whole_transcriptome_gene_summary.tsv", sep="\t", index=False
        )
    metadata = json.loads((metadata_dir / "long_read_spatial_clustering.json").read_text(encoding="utf-8"))
    display_frame = frame.copy()
    display_frame["bin_x"] *= float(metadata.get("coordinate_um_per_unit", 1.0))
    display_frame["bin_y"] *= float(metadata.get("coordinate_um_per_unit", 1.0))
    plot_categorical(
        display_frame, "bin_x", "bin_y", figure_dir / "long_read_spatial_clusters",
        "Spatial expression clusters", args.dpi, True, False,
    )
    plot_categorical(
        frame, "UMAP1", "UMAP2", figure_dir / "long_read_umap_clusters",
        "Depth-adjusted expression UMAP", args.dpi, False, False,
    )
    plot_continuous(
        display_frame, "total_lr_molecules", figure_dir / "long_read_molecules",
        "Long-read molecules per spatial bin", "Long-read molecules (2nd–98th percentile)",
        args.dpi, GENERAL_CMAP,
    )
    marker_figures(
        adata,
        pd.read_csv(table_dir / "long_read_cluster_markers.tsv", sep="\t"),
        figure_dir,
        args.marker_display_top,
        args.dpi,
    )
    qc_figures(
        frame,
        pd.read_csv(table_dir / "long_read_spatial_coherence.tsv", sep="\t"),
        figure_dir,
        args.dpi,
    )


def main() -> None:
    args = parse_args()
    apply_style(args.font_family, args.base_size)
    if args.render_only:
        render_existing(Path(args.outdir), args)
        return
    if not args.molecule_calls or not args.sample:
        raise SystemExit("--molecule-calls and --sample are required unless --render-only is used")
    if args.bin_size <= 0 or args.min_bin_molecules <= 0 or args.min_gene_bins <= 0:
        raise SystemExit("Bin size and minimum count thresholds must be positive")
    if not 0 <= args.spatial_graph_weight < 1:
        raise SystemExit("--spatial-graph-weight must be in [0, 1)")
    if args.min_cluster_bins <= 0:
        raise SystemExit("--min-cluster-bins must be positive")
    outdir = Path(args.outdir)
    table_dir = outdir / "01_tables"
    figure_dir = outdir / "02_figures"
    tool_dir = outdir / "04_tool_output" / "scanpy"
    metadata_dir = outdir / "05_metadata"
    for directory in (table_dir, figure_dir, tool_dir, metadata_dir):
        directory.mkdir(parents=True, exist_ok=True)
    numba_cache = outdir / "04_tool_output" / "numba-cache"
    numba_cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("NUMBA_CACHE_DIR", str(numba_cache.resolve()))

    matrix, bin_ids, coords, gene_ids, gene_names, class_counts, qc = build_matrix(
        Path(args.molecule_calls), args.bin_size, args.chunk_rows
    )
    totals = np.asarray(matrix.sum(axis=1)).ravel()
    keep_bins = totals >= args.min_bin_molecules
    matrix = matrix[keep_bins]
    kept_bin_ids = [value for value, keep in zip(bin_ids, keep_bins) if keep]
    kept_coords = [value for value, keep in zip(coords, keep_bins) if keep]
    totals = totals[keep_bins]
    gene_bin_counts = np.asarray((matrix > 0).sum(axis=0)).ravel()
    keep_genes = gene_bin_counts >= args.min_gene_bins
    matrix = matrix[:, keep_genes]
    kept_gene_ids = [value for value, keep in zip(gene_ids, keep_genes) if keep]
    kept_gene_names = [value for value, keep in zip(gene_names, keep_genes) if keep]
    nonzero_bins = np.asarray(matrix.sum(axis=1)).ravel() > 0
    matrix = matrix[nonzero_bins]
    kept_bin_ids = [value for value, keep in zip(kept_bin_ids, nonzero_bins) if keep]
    kept_coords = [value for value, keep in zip(kept_coords, nonzero_bins) if keep]
    totals = totals[nonzero_bins]
    detected_genes = np.asarray((matrix > 0).sum(axis=1)).ravel().astype(int)
    if matrix.shape[0] < 10 or matrix.shape[1] < 50:
        raise SystemExit(f"Insufficient long-read matrix after filtering: bins={matrix.shape[0]} genes={matrix.shape[1]}")

    import anndata as ad
    import scanpy as sc

    obs = pd.DataFrame(index=kept_bin_ids)
    obs.index.name = "bin_id"
    obs["bin_x"] = [value[0] for value in kept_coords]
    obs["bin_y"] = [value[1] for value in kept_coords]
    obs["total_lr_molecules"] = totals.astype(int)
    for label in CLASSES:
        obs[label] = [class_counts[bin_id][label] for bin_id in kept_bin_ids]
    obs["new_supported_fraction"] = np.where(
        obs["total_lr_molecules"] > 0,
        obs["new_supported"] / obs["total_lr_molecules"],
        np.nan,
    )
    var = pd.DataFrame({"gene_name": kept_gene_names}, index=kept_gene_ids)
    var.index.name = "gene_id"
    adata = ad.AnnData(X=matrix, obs=obs, var=var)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=min(args.hvg, adata.n_vars), flavor="seurat")
    hvg_mask = adata.var["highly_variable"].to_numpy()
    cluster_adata = adata[:, hvg_mask].copy()
    sc.pp.scale(cluster_adata, max_value=10)
    n_pcs = max(2, min(args.pca_components, cluster_adata.n_obs - 1, cluster_adata.n_vars - 1))
    sc.tl.pca(cluster_adata, n_comps=n_pcs, svd_solver="arpack", random_state=args.seed)
    raw_pca = np.asarray(cluster_adata.obsm["X_pca"]).copy()
    log_depth = np.log1p(totals.astype(float))
    raw_pc1_depth = pearsonr(raw_pca[:, 0], log_depth)
    use_depth_residual = truth(args.regress_log_depth)
    if use_depth_residual:
        design = np.column_stack([np.ones(len(log_depth)), log_depth])
        pca_for_neighbors = raw_pca - design @ np.linalg.lstsq(design, raw_pca, rcond=None)[0]
        cluster_adata.obsm["X_pca_depth_residual"] = pca_for_neighbors
        neighbor_representation = "X_pca_depth_residual"
    else:
        pca_for_neighbors = raw_pca
        neighbor_representation = "X_pca"
    n_neighbors = max(2, min(args.neighbors, cluster_adata.n_obs - 1))
    sc.pp.neighbors(
        cluster_adata,
        n_neighbors=n_neighbors,
        use_rep=neighbor_representation,
        random_state=args.seed,
    )
    expression_graph = symmetric_normalize(cluster_adata.obsp["connectivities"])
    spatial_graph_raw = make_spatial_graph(kept_coords, args.bin_size)
    spatial_graph = symmetric_normalize(spatial_graph_raw)
    clustering_graph = (
        (1.0 - args.spatial_graph_weight) * expression_graph
        + args.spatial_graph_weight * spatial_graph
    ).tocsr()
    if args.leiden_resolution <= 0:
        raise SystemExit("--leiden-resolution must be positive")
    sc.tl.leiden(
        cluster_adata,
        adjacency=clustering_graph,
        resolution=args.leiden_resolution,
        random_state=args.seed,
        key_added="leiden",
        flavor="leidenalg",
    )
    raw_labels = cluster_adata.obs["leiden"].astype(str)
    selected_labels = merge_small_clusters(clustering_graph, raw_labels, args.min_cluster_bins)
    for legacy in (
        table_dir / "long_read_resolution_scan.tsv",
        figure_dir / "long_read_resolution_scan.png",
        figure_dir / "long_read_resolution_scan.pdf",
    ):
        if legacy.exists():
            legacy.unlink()
    sc.tl.umap(cluster_adata, random_state=args.seed)
    labels = selected_labels
    ordered = labels.value_counts().index.tolist()
    remap = {label: str(index + 1) for index, label in enumerate(ordered)}
    clusters = labels.map(remap)
    if clusters.nunique() > 16:
        raise SystemExit(f"Leiden produced {clusters.nunique()} clusters; lower the configured resolution or provide a larger explicit palette")
    adata.obs["cluster"] = pd.Categorical(clusters, categories=[str(index + 1) for index in range(clusters.nunique())], ordered=True)
    adata.obsm["X_umap"] = cluster_adata.obsm["X_umap"]
    adata.obsm["X_pca"] = raw_pca
    adata.obsm["X_pca_clustering"] = pca_for_neighbors
    adata.obsp["connectivities"] = cluster_adata.obsp["connectivities"]
    adata.obsp["clustering_connectivities"] = clustering_graph
    adata.obsp["distances"] = cluster_adata.obsp["distances"]
    adata.uns["neighbors"] = cluster_adata.uns["neighbors"]
    adata.uns["umap"] = cluster_adata.uns["umap"]

    frame = adata.obs.reset_index()
    frame["UMAP1"] = adata.obsm["X_umap"][:, 0]
    frame["UMAP2"] = adata.obsm["X_umap"][:, 1]
    frame.to_csv(table_dir / "long_read_bin50_clusters.tsv", sep="\t", index=False)
    frame["detected_genes"] = detected_genes
    gene_summary = pd.DataFrame(
        {
            "gene_id": kept_gene_ids,
            "gene_name": kept_gene_names,
            "total_molecules": np.asarray(matrix.sum(axis=0)).ravel().astype(int),
            "detected_bins": np.asarray((matrix > 0).sum(axis=0)).ravel().astype(int),
        }
    ).sort_values(["total_molecules", "gene_id"], ascending=[False, True])
    gene_summary.to_csv(table_dir / "whole_transcriptome_gene_summary.tsv", sep="\t", index=False)
    markers = marker_table(
        sc,
        adata,
        table_dir / "long_read_cluster_markers.tsv",
        args.marker_top,
        args.marker_min_logfoldchange,
        args.marker_max_fdr,
    )
    marker_figures(adata, markers, figure_dir, args.marker_display_top, args.dpi)
    cluster_summary = (
        frame.groupby("cluster", observed=True)
        .agg(
            bins=("bin_id", "size"),
            total_lr_molecules=("total_lr_molecules", "sum"),
            median_lr_molecules=("total_lr_molecules", "median"),
            median_detected_genes=("detected_genes", "median"),
            new_supported=("new_supported", "sum"),
            median_new_supported_fraction=("new_supported_fraction", "median"),
        )
        .reset_index()
    )
    cluster_summary.to_csv(table_dir / "long_read_cluster_summary.tsv", sep="\t", index=False)

    display_frame = frame.copy()
    display_frame["bin_x"] *= args.coordinate_um_per_unit
    display_frame["bin_y"] *= args.coordinate_um_per_unit
    palette = plot_categorical(display_frame, "bin_x", "bin_y", figure_dir / "long_read_spatial_clusters", "Spatial expression clusters", args.dpi, True, False)
    plot_categorical(frame, "UMAP1", "UMAP2", figure_dir / "long_read_umap_clusters", "Depth-adjusted expression UMAP", args.dpi, False, False)
    plot_continuous(
        display_frame, "total_lr_molecules", figure_dir / "long_read_molecules",
        "Long-read molecules per spatial bin", "Long-read molecules (2nd–98th percentile display)",
        args.dpi, GENERAL_CMAP,
    )
    coherence, coherence_overall = spatial_coherence(frame, args.bin_size)
    coherence.to_csv(table_dir / "long_read_spatial_coherence.tsv", sep="\t", index=False)
    qc_figures(frame, coherence, figure_dir, args.dpi)
    adjusted_pc1_depth = pearsonr(pca_for_neighbors[:, 0], log_depth)
    final_depth_eta = depth_eta_squared(log_depth, frame["cluster"].astype(str))
    clustering_qc = {
        "raw_pc1_log_depth_pearson_r": float(raw_pc1_depth.statistic),
        "raw_pc1_log_depth_pearson_pvalue": float(raw_pc1_depth.pvalue),
        "pc1_log_depth_pearson_r": float(adjusted_pc1_depth.statistic),
        "pc1_log_depth_pearson_pvalue": float(adjusted_pc1_depth.pvalue),
        "cluster_log_depth_eta_squared": final_depth_eta,
        "expression_edge_agreement": weighted_edge_agreement(expression_graph, frame["cluster"].astype(str)),
        "spatial_graph_weight": args.spatial_graph_weight,
        "minimum_cluster_bins": args.min_cluster_bins,
        **coherence_overall,
    }
    (metadata_dir / "clustering_qc.json").write_text(json.dumps(clustering_qc, indent=2), encoding="utf-8")

    adata.uns["long_read_spatial_clustering"] = {
        "coordinate_used_for_clustering": bool(args.spatial_graph_weight > 0),
        "coordinate_source": "long-read barcode assignment against the chip whitelist",
        "bin_size_raw_units": args.bin_size,
        "coordinate_um_per_unit": args.coordinate_um_per_unit,
        "bin_size_um": args.bin_size * args.coordinate_um_per_unit,
        "min_bin_molecules": args.min_bin_molecules,
        "min_gene_bins": args.min_gene_bins,
        "hvg": int(hvg_mask.sum()),
        "pca_components": n_pcs,
        "neighbors": n_neighbors,
        "regress_log_depth": use_depth_residual,
        "spatial_graph_weight": args.spatial_graph_weight,
        "minimum_cluster_bins": args.min_cluster_bins,
        "leiden_resolution": args.leiden_resolution,
        "seed": args.seed,
    }
    adata.write_h5ad(tool_dir / "long_read_spatial_bin50.h5ad", compression="gzip")
    coordinate_valid = bool(
        (frame["bin_x"] >= 0).all()
        and (frame["bin_x"] <= args.chip_bounds).all()
        and (frame["bin_y"] >= 0).all()
        and (frame["bin_y"] <= args.chip_bounds).all()
    )
    summary = {
        "sample": args.sample,
        "method": "depth-adjusted Scanpy expression and spatial-graph clustering",
        "coordinate_used_for_clustering": bool(args.spatial_graph_weight > 0),
        "coordinate_source": "long-read barcode assignment against the chip whitelist",
        "canonical_origin": "upper-left",
        "canonical_x_direction": "right",
        "canonical_y_direction": "down",
        "coordinate_bounds_valid": coordinate_valid,
        "bin_size_raw_units": args.bin_size,
        "coordinate_um_per_unit": args.coordinate_um_per_unit,
        "bin_size_um": args.bin_size * args.coordinate_um_per_unit,
        "input_molecules": int(qc["input_molecules"]),
        "valid_coordinate_gene_molecules": int(qc["valid_molecules"]),
        "clustered_molecules": int(frame["total_lr_molecules"].sum()),
        "clustered_bins": int(len(frame)),
        "genes_after_filter": int(adata.n_vars),
        "highly_variable_genes": int(hvg_mask.sum()),
        "clusters": int(frame["cluster"].nunique()),
        "cluster_counts": {str(row.cluster): int(row.bins) for row in cluster_summary.itertuples(index=False)},
        "scanpy_version": sc.__version__,
        "anndata_version": ad.__version__,
        "algorithm": "normalize_total(1e4)/log1p/Seurat-HVG/scale/PCA/depth-residualization/expression-spatial graph/Leiden/minimum-size merge/UMAP",
        "parameters": {
            "min_bin_molecules": args.min_bin_molecules,
            "min_gene_bins": args.min_gene_bins,
            "hvg": args.hvg,
            "pca_components": n_pcs,
            "neighbors": n_neighbors,
            "regress_log_depth": use_depth_residual,
            "spatial_graph_weight": args.spatial_graph_weight,
            "minimum_cluster_bins": args.min_cluster_bins,
            "leiden_resolution": args.leiden_resolution,
            "seed": args.seed,
            "marker_method": "t-test_overestim_var; descriptive ranking only",
        },
        "cluster_color_map": palette,
        "depth_qc": clustering_qc,
        "continuous_display_limits": "2nd-98th percentile; values unchanged",
    }
    (metadata_dir / "long_read_spatial_clustering.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    with (metadata_dir / "long_read_spatial_input.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        source = Path(args.molecule_calls).resolve()
        writer.writerow(["role", "path", "size_bytes", "mtime_ns"])
        writer.writerow(["long_read_molecule_calls", source, source.stat().st_size, source.stat().st_mtime_ns])


if __name__ == "__main__":
    main()
