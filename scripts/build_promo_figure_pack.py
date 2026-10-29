#!/usr/bin/env python3
"""Build a reproducible PPT figure pack from completed LR-Celatlas outputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
import shutil
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle


CLUSTER_COLORS = {
    "1": "#E64B35", "2": "#4DBBD5", "3": "#00A087", "4": "#3C5488",
    "5": "#F39B7F", "6": "#8491B4", "7": "#91D1C2", "8": "#DC0000",
    "9": "#7E6148", "10": "#B09C85", "11": "#6F99AD", "12": "#8A7AA9",
}
CONTINUOUS = LinearSegmentedColormap.from_list(
    "npg_magnitude", ["#E5E7EB", "#4DBBD5", "#3C5488", "#E64B35"]
)
ACTIVITY = LinearSegmentedColormap.from_list(
    "npg_activity", ["#E5E7EB", "#91D1C2", "#00A087", "#DC0000"]
)
LEGACY_GROUPED_STEMS = {
    "02_spatial_qc",
    "03_long_read_length_n50",
    "04_barcode_backtrace",
    "05_mismatch_spectrum_slam",
    "04_spatial_inventory",
    "07_transcript_structure_classes",
    "07_known_novel_transcripts",
    "07_gffcompare_classes",
    "08_total_vs_new_gene_spatial",
    "10_candidate_gene_case",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-dir", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--dpi", default=300, type=int)
    parser.add_argument("--candidate-min-informative", default=1000, type=int)
    parser.add_argument("--bin-size", default=50, type=int)
    parser.add_argument("--scan-candidate-calls", action="store_true", help="Scan the large molecule table for gene-specific per-Bin calls.")
    return parser.parse_args()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_metrics(path: Path) -> dict[str, str]:
    frame = pd.read_csv(path, sep="\t", dtype=str)
    return dict(zip(frame.iloc[:, 0], frame.iloc[:, 1]))


def save_figure(fig, stem: Path, dpi: int) -> tuple[Path, Path]:
    stem.parent.mkdir(parents=True, exist_ok=True)
    png_path = stem.with_suffix(".png")
    pdf_path = stem.with_suffix(".pdf")
    fig.savefig(png_path, dpi=dpi, facecolor="white", bbox_inches="tight")
    try:
        from PIL import Image
        with Image.open(png_path) as image:
            image.convert("RGB").save(pdf_path, resolution=120.0)
    except ImportError:
        fig.savefig(pdf_path, dpi=min(dpi, 120), facecolor="white", bbox_inches="tight")
    plt.close(fig)
    return png_path, pdf_path


def style_axes(ax) -> None:
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=8)


def spatial_axes(ax, bins: pd.DataFrame, title: str) -> None:
    ax.set_title(title, loc="left", fontsize=11, fontweight="bold")
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("y (µm)")
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xlim(-100, max(9457, float(bins["bin_x"].max()) + 100))
    ax.set_ylim(max(9457, float(bins["bin_y"].max()) + 100), -100)
    style_axes(ax)


def plot_spatial_values(ax, bins: pd.DataFrame, column: str, title: str, cmap, log=False) -> None:
    values = pd.to_numeric(bins[column], errors="coerce").fillna(0)
    display = (values + 1).map(math.log10) if log else values
    low, high = display.quantile(0.02), display.quantile(0.98)
    if high <= low:
        low, high = float(display.min()), float(display.max() or 1)
    sc = ax.scatter(
        bins["bin_x"], bins["bin_y"], c=display, s=6, cmap=cmap,
        vmin=low, vmax=high, linewidths=0, rasterized=True,
    )
    spatial_axes(ax, bins, title)
    cb = ax.figure.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cb.ax.tick_params(labelsize=7)
    cb.set_label("log10(value+1)" if log else column, fontsize=8)


def robust_display_limits(values: pd.Series, lower: float = 0.02, upper: float = 0.98) -> tuple[float, float]:
    numeric = pd.to_numeric(values, errors="coerce")
    numeric = numeric[numeric.notna() & numeric.map(math.isfinite)]
    if numeric.empty:
        return 0.0, 1.0
    low, high = float(numeric.quantile(lower)), float(numeric.quantile(upper))
    if high <= low:
        low, high = float(numeric.min()), float(numeric.max())
    if high <= low:
        high = low + 1.0
    return low, high


def spatial_bin_marker_size(bin_size: float) -> float:
    """Keep a square marker close to the displayed physical Bin diameter."""
    return max(2.0, min(8.0, 4.0 * (float(bin_size) / 50.0) ** 2))


def remove_legacy_grouped_figures(figdir: Path) -> None:
    """Keep the delivery directory single-panel after a format migration."""
    for stem in LEGACY_GROUPED_STEMS:
        for suffix in (".png", ".pdf"):
            (figdir / f"{stem}{suffix}").unlink(missing_ok=True)


def bin_id(x: int, y: int, bin_size: int) -> str:
    bx = (x // bin_size) * bin_size + bin_size // 2
    by = (y // bin_size) * bin_size + bin_size // 2
    return f"{bx}_{by}"


def load_bin_umi_counts(path: Path, bin_size: int) -> pd.DataFrame:
    """Count the UMI-level molecule-call rows in each physical Bin."""
    if not path.is_file():
        return pd.DataFrame(columns=["bin_id", "umi_count"])
    counts: defaultdict[str, int] = defaultdict(int)
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            try:
                counts[bin_id(int(float(row["x"])), int(float(row["y"])), bin_size)] += 1
            except (KeyError, TypeError, ValueError):
                continue
    return pd.DataFrame({"bin_id": list(counts), "umi_count": list(counts.values())})


def select_candidates(gene_table: pd.DataFrame, minimum: int) -> list[dict[str, str]]:
    if gene_table.empty:
        return []
    frame = gene_table.copy()
    frame["informative_molecules"] = pd.to_numeric(frame["informative_molecules"], errors="coerce").fillna(0)
    frame["new_supported_fraction"] = pd.to_numeric(frame["new_supported_fraction"], errors="coerce").fillna(0)
    frame = frame[frame["informative_molecules"] >= minimum]
    frame = frame[~frame["gene_name"].astype(str).str.match(r"^(mt-|Rpl|Rps|Gm)", case=False, na=False)]
    frame = frame.sort_values(
        ["new_supported_fraction", "informative_molecules", "gene_name"],
        ascending=[False, False, True],
    )
    return frame.head(20).to_dict("records")


def scan_candidate_calls(path: Path, candidates: list[dict[str, str]], bin_size: int) -> dict[str, dict[str, list[int]]]:
    wanted = {str(row["gene_id"]): row for row in candidates}
    profiles: dict[str, dict[str, list[int]]] = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    if not wanted or not path.is_file():
        return profiles
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            gene_id = row.get("gene_id", "")
            if gene_id not in wanted:
                continue
            try:
                x, y = int(float(row["x"])), int(float(row["y"]))
            except (KeyError, TypeError, ValueError):
                continue
            key = bin_id(x, y, bin_size)
            values = profiles[gene_id][key]
            values[0] += 1
            if row.get("exploratory_class", row.get("final_call", "")) == "new_supported":
                values[1] += 1
    return profiles


def choose_candidate(candidates: list[dict[str, str]], profiles: dict[str, dict[str, list[int]]]) -> dict[str, str]:
    if not candidates:
        return {"gene_id": "", "gene_name": "Unavailable"}
    scored = []
    for row in candidates:
        profile = profiles.get(str(row["gene_id"]), {})
        fractions = [new / total for total, new in profile.values() if total >= 3]
        spread = (max(fractions) - min(fractions)) if fractions else -1
        scored.append((spread, len(profile), int(row.get("informative_molecules", 0)), str(row["gene_name"]), row))
    return max(scored, key=lambda item: item[:4])[4]


def load_gene_spatial_profile(h5ad_path: Path, gene_id: str) -> pd.DataFrame:
    if not h5ad_path.is_file() or not gene_id:
        return pd.DataFrame()
    try:
        import anndata
    except ImportError:
        return pd.DataFrame()
    try:
        adata = anndata.read_h5ad(h5ad_path, backed="r")
        if gene_id not in adata.var_names:
            adata.file.close()
            return pd.DataFrame()
        values = adata[:, gene_id].X
        if hasattr(values, "toarray"):
            values = values.toarray()
        values = values.reshape(-1)
        frame = adata.obs[["bin_x", "bin_y"]].copy()
        frame["total"] = values
        adata.file.close()
        return frame.reset_index(names="bin_id")
    except Exception:
        return pd.DataFrame()


def parse_gtf(gtf: Path, gene_id: str) -> list[tuple[str, int, int]]:
    if not gtf.is_file() or not gene_id:
        return []
    transcripts: dict[str, list[tuple[int, int]]] = defaultdict(list)
    pattern = re.compile(r'gene_id "([^"]+)";.*?transcript_id "([^"]+)"')
    with gtf.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith("#") or "\texon\t" not in line:
                continue
            match = pattern.search(line)
            if not match or match.group(1) != gene_id:
                continue
            fields = line.rstrip("\n").split("\t")
            try:
                transcripts[match.group(2)].append((int(fields[3]), int(fields[4])))
            except (IndexError, ValueError):
                continue
    if not transcripts:
        return []
    transcript, exons = max(transcripts.items(), key=lambda item: (len(item[1]), item[0]))
    return [(transcript, start, end) for start, end in sorted(exons)]


def parse_gff_attributes(value: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in value.split(";"):
        if "=" not in item:
            continue
        key, raw = item.split("=", 1)
        attributes[key] = raw
    return attributes


def parse_isoform_classes(tmap: Path) -> dict[str, str]:
    if not tmap.is_file():
        return {}
    try:
        frame = pd.read_csv(tmap, sep="\t")
    except (OSError, pd.errors.ParserError):
        return {}
    if "qry_id" not in frame or "class_code" not in frame:
        return {}
    return {
        str(row.qry_id).replace("transcript:", ""): str(row.class_code)
        for row in frame[["qry_id", "class_code"]].itertuples(index=False)
    }


def parse_isoform_tracks(gff3: Path, gene_id: str, tmap: Path, limit: int = 6) -> list[dict[str, object]]:
    """Read observed FLAMES transcript/exon tracks for one gene."""
    if not gff3.is_file() or not gene_id:
        return []
    classes = parse_isoform_classes(tmap)
    tracks: dict[str, dict[str, object]] = {}
    try:
        handle = gff3.open(encoding="utf-8", errors="ignore")
    except OSError:
        return []
    with handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in {"transcript", "exon"}:
                continue
            attrs = parse_gff_attributes(fields[8])
            if fields[2] == "transcript":
                parent_gene = attrs.get("Parent", "").replace("gene:", "")
                if attrs.get("gene_id", parent_gene) != gene_id and parent_gene != gene_id:
                    continue
                transcript_id = attrs.get("transcript_id", attrs.get("ID", "")).replace("transcript:", "")
                tracks[transcript_id] = {
                    "transcript_id": transcript_id,
                    "support": int(float(attrs.get("support_count", 0) or 0)),
                    "source": attrs.get("source", fields[1]),
                    "class_code": classes.get(transcript_id, "?"),
                    "exons": [],
                }
            elif fields[2] == "exon":
                parent = attrs.get("Parent", "").replace("transcript:", "")
                if parent not in tracks:
                    continue
                try:
                    tracks[parent]["exons"].append((int(fields[3]), int(fields[4])))
                except (IndexError, ValueError):
                    continue
    selected = [row for row in tracks.values() if row["exons"]]
    selected.sort(key=lambda row: (int(row["support"]), len(row["exons"]), str(row["transcript_id"])), reverse=True)
    return selected[:limit]


def write_manifest(outdir: Path, rows: list[dict[str, str]]) -> None:
    fields = ["figure_id", "level", "filename_png", "filename_pdf", "ppt_use", "source_paths", "interpretation", "limitation"]
    with (outdir / "01_tables" / "figure_manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def figure_record(stem: Path, figure_id: str, level: str, use: str, sources: list[Path], interpretation: str, limitation: str) -> dict[str, str]:
    return {
        "figure_id": figure_id,
        "level": level,
        "filename_png": str(stem.with_suffix(".png").name),
        "filename_pdf": str(stem.with_suffix(".pdf").name),
        "ppt_use": use,
        "source_paths": ";".join(str(path.resolve()) for path in sources),
        "interpretation": interpretation,
        "limitation": limitation,
    }


def main() -> None:
    args = parse_args()
    sample = args.sample_dir
    outdir = args.outdir
    figdir = outdir / "02_figures"
    table_dir = outdir / "01_tables"
    metadata_dir = outdir / "05_metadata"
    for directory in (figdir, table_dir, outdir / "03_logs", outdir / "04_tool_output", metadata_dir):
        directory.mkdir(parents=True, exist_ok=True)
    remove_legacy_grouped_figures(figdir)

    s1, s2, s3, s4, s5, s6, s7 = [sample / name for name in (
        "01_read_qc", "02_barcode_cdna", "03_align", "04_full_length_transcriptome",
        "05_slam", "06_whole_transcriptome", "07_slam_joint",
    )]
    qc_path = s1 / "01_tables/qc_metrics.tsv"
    bins_path = s6 / "01_tables/long_read_bin50_clusters.tsv"
    cluster_path = s6 / "01_tables/long_read_cluster_summary.tsv"
    barcode_level_path = s3 / "01_tables/barcode_level_summary.tsv"
    gene_path = s5 / "01_tables/gene_new_supported_summary.tsv"
    calls_path = s5 / "01_tables" / f"{args.sample_id}.slam_molecule_calls.tsv.gz"
    if not calls_path.is_file():
        calls_path = s5 / "04_tool_output/caller" / f"{args.sample_id}.molecule_calls.tsv.gz"
    required = [qc_path, bins_path, cluster_path, gene_path]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise SystemExit("Missing required completed outputs:\n" + "\n".join(missing))

    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.titlesize": 11, "figure.dpi": args.dpi})
    qc = read_metrics(qc_path)
    bins = pd.read_csv(bins_path, sep="\t")
    cluster_summary = pd.read_csv(cluster_path, sep="\t")
    barcode_level = pd.read_csv(barcode_level_path, sep="\t") if barcode_level_path.is_file() else pd.DataFrame()
    gene_table = pd.read_csv(gene_path, sep="\t")
    structure = read_json(s4 / "05_metadata/isoform_class_summary.json")
    isoform_metrics = read_metrics(s4 / "01_tables/isoform_summary.tsv")
    structure_summary = pd.read_csv(s4 / "01_tables/isoform_structure_summary.tsv", sep="\t")
    barcode = read_json(s2 / "05_metadata/barcode_summary.json")
    extract = read_json(s2 / "05_metadata/extract_summary.json")
    alignment = read_json(s3 / "05_metadata/alignment_summary.json")
    slam = read_json(s5 / "05_metadata" / f"{args.sample_id}.slam_summary.json")
    mismatch = read_json(s5 / "05_metadata/mismatch_summary.json")
    spatial = read_json(s6 / "05_metadata/long_read_spatial_clustering.json")
    joint = read_json(s7 / "05_metadata/slam_joint_summary.json")
    bin_umi = load_bin_umi_counts(calls_path, args.bin_size)
    umi_total = int(bin_umi["umi_count"].sum()) if not bin_umi.empty else int(slam.get("molecules_total", 0))
    median_umi = float(bin_umi["umi_count"].median()) if not bin_umi.empty else float("nan")
    candidates = select_candidates(gene_table, args.candidate_min_informative)
    profiles = scan_candidate_calls(calls_path, candidates, args.bin_size) if args.scan_candidate_calls else {}
    candidate = choose_candidate(candidates, profiles)
    candidate_id = str(candidate.get("gene_id", ""))
    candidate_name = str(candidate.get("gene_name", "Unavailable"))
    candidate_profile = profiles.get(candidate_id, {})
    bin_lookup = bins.set_index("bin_id")
    h5ad_path = s6 / "04_tool_output/scanpy/long_read_spatial_bin50.h5ad"
    (outdir / "03_logs" / "build_promo_figures.log").write_text(
        "\n".join([
            f"[INFO] source_root={sample}",
            f"[INFO] output_root={outdir}",
            "[INFO] mode=derived_from_completed_results",
            f"[INFO] candidate_call_scan={args.scan_candidate_calls}",
            f"[INFO] derived_figure_dpi={args.dpi}",
            f"[INFO] candidate_gene={candidate_name}",
        ]) + "\n",
        encoding="utf-8",
    )
    input_samples = pd.read_csv(s1 / "05_metadata/input_samples.tsv", sep="\t", dtype=str)
    gtf_value = input_samples.iloc[0].get("gtf", "") if not input_samples.empty else ""
    gtf_path = Path(str(gtf_value)) if gtf_value else Path()
    structure_exons = parse_gtf(gtf_path, candidate_id)
    isoform_gff3 = s4 / "04_tool_output/flames/isoforms.gff3"
    isoform_tmap = s4 / "04_tool_output/flames/flames_vs_reference.isoforms.gff3.tmap"
    isoform_tracks = parse_isoform_tracks(isoform_gff3, candidate_id, isoform_tmap)

    records: list[dict[str, str]] = []
    # 01: coordinate-only spatial structure map.
    stem = figdir / "01_spatial_structure"
    fig, ax = plt.subplots(figsize=(7.0, 6.4))
    for cluster, group in bins.groupby(bins["cluster"].astype(str), sort=True):
        ax.scatter(group["bin_x"], group["bin_y"], s=5, color=CLUSTER_COLORS.get(cluster, "#777777"), linewidths=0, rasterized=True)
    spatial_axes(ax, bins, f"{args.sample_id} spatial Bin structure (50 µm)")
    handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=CLUSTER_COLORS.get(str(c), "#777"), markersize=6, label=f"Cluster {c}") for c in sorted(bins["cluster"].astype(str).unique(), key=int)]
    ax.legend(handles=handles, title="Spatial cluster", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7, title_fontsize=8, frameon=False)
    png, pdf = save_figure(fig, stem, args.dpi)
    records.append(figure_record(stem, "01_spatial_structure", "A", "空间结构与常规表达底图", [bins_path], "50 µm Bin 坐标和长读长表达簇保持空间连续性。", "未提供 H&E；这不是组织学图像叠加。"))

    # 02: single-panel spatial QC figures.
    stem = figdir / "02_spatial_molecules_per_bin"
    fig, ax = plt.subplots(figsize=(6.5, 5.8))
    plot_spatial_values(ax, bins, "total_lr_molecules", "Long-read molecules per Bin", CONTINUOUS, log=True)
    fig.tight_layout()
    save_figure(fig, stem, args.dpi)
    records.append(figure_record(stem, "02_spatial_molecules_per_bin", "A", "空间分子深度", [bins_path], "展示每个 50 µm Bin 的长读长分子数。", "没有 H&E；空间坐标图不是组织学叠加。"))

    stem = figdir / "02_bin_coverage_distribution"
    fig, ax = plt.subplots(figsize=(6.5, 4.8))
    ax.hist(pd.to_numeric(bins["total_lr_molecules"]), bins=40, color="#3C5488", alpha=0.88)
    ax.set_xlabel("Long-read molecules per Bin")
    ax.set_ylabel("Bins")
    ax.set_title("Bin coverage distribution", loc="left", fontsize=11, fontweight="bold")
    style_axes(ax)
    fig.tight_layout()
    save_figure(fig, stem, args.dpi)
    records.append(figure_record(stem, "02_bin_coverage_distribution", "A", "空间深度分布", [bins_path], "展示空间 Bin 的长读长分子深度分布。", "分布基于 long-read molecule，不是独立 UMI。"))

    stem = figdir / "02_cluster_depth_qc"
    fig, ax = plt.subplots(figsize=(6.5, 4.8))
    med = cluster_summary.sort_values("cluster")
    ax.barh(med["cluster"].astype(str), med["median_lr_molecules"], color=[CLUSTER_COLORS.get(str(x), "#777") for x in med["cluster"]])
    ax.set_xlabel("Median long-read molecules per Bin")
    ax.set_ylabel("Spatial cluster")
    ax.set_title("Cluster depth QC", loc="left", fontsize=11, fontweight="bold")
    style_axes(ax)
    fig.tight_layout()
    save_figure(fig, stem, args.dpi)
    records.append(figure_record(stem, "02_cluster_depth_qc", "A", "空间 cluster 深度质控", [cluster_path], "比较各空间 cluster 的中位长读长分子深度。", "cluster 是长读长表达/空间图上的分析标签。"))

    stem = figdir / "02_spatial_qc_counts"
    qc_counts = {
        "Raw reads": int(qc["n_reads"]),
        "Spatial bins": len(bins),
        "Genes after filter": int(spatial.get("genes_after_filter", 0)),
        "Spatial clusters": int(spatial.get("clusters", 0)),
    }
    fig, ax = plt.subplots(figsize=(6.5, 4.8))
    keys, values = list(qc_counts), list(qc_counts.values())
    ax.barh(keys[::-1], values[::-1], color="#4DBBD5")
    ax.set_xscale("log")
    ax.set_xlabel("Observed count (log scale)")
    ax.set_title("Spatial QC counts", loc="left", fontsize=11, fontweight="bold")
    for index, value in enumerate(values[::-1]):
        ax.text(value * 1.08, index, f"{value:,}", va="center", fontsize=8)
    style_axes(ax)
    fig.tight_layout()
    save_figure(fig, stem, args.dpi)
    records.append(figure_record(stem, "02_spatial_qc_counts", "A", "空间质控计数", [qc_path, bins_path, s6 / "05_metadata/long_read_spatial_clustering.json"], "汇总 raw reads、空间 Bin、过滤后 gene 和空间 cluster 数量。", "当前结果树没有传统短读长 UMI/saturation 指标。"))

    # 03: single-panel read-length and QC figures.
    source_plot = s1 / "02_figures/raw_read_length.png"
    stem = figdir / "03_read_length_distribution"
    if source_plot.is_file():
        shutil.copy2(source_plot, stem.with_suffix(".png"))
        try:
            from PIL import Image
            with Image.open(source_plot) as image:
                image.convert("RGB").save(stem.with_suffix(".pdf"), resolution=120.0)
        except ImportError:
            fig, ax = plt.subplots(figsize=(7.0, 4.0))
            ax.imshow(plt.imread(source_plot))
            ax.axis("off")
            save_figure(fig, stem, args.dpi)
    else:
        fig, ax = plt.subplots(figsize=(7.0, 4.0))
        ax.axis("off")
        ax.text(0.5, 0.5, "Read-length distribution unavailable", ha="center", va="center")
        save_figure(fig, stem, args.dpi)
    records.append(figure_record(stem, "03_read_length_distribution", "A", "长读长 read length 分布", [qc_path, source_plot], "展示原始 read length 分布。", "该图是原始文库 QC，不是标记前后比较。"))

    stem = figdir / "03_long_read_qc_summary"
    fig, ax = plt.subplots(figsize=(6.5, 4.8))
    ax.axis("off")
    ax.set_title("Long-read QC summary", loc="left", fontsize=11, fontweight="bold")
    ax.text(0, 0.84, "\n".join([
        f"Reads: {int(qc['n_reads']):,}", f"Mean length: {float(qc['mean_len']):,.1f} bp",
        f"Median length: {int(float(qc['median_len'])):,} bp", f"N50: {int(float(qc['N50'])):,} bp",
        f"Mean Q: {float(qc['mean_q']):.2f}",
    ]), va="top", fontsize=12, linespacing=1.8)
    fig.tight_layout()
    save_figure(fig, stem, args.dpi)
    records.append(figure_record(stem, "03_long_read_qc_summary", "A", "长读长 QC 指标", [qc_path], "展示 reads、mean/median length、N50 和平均质量值。", "该图是原始文库 QC，不是标记前后比较。"))

    # 04: single-panel spatial barcode and inventory figures.
    stage_values = {
        "Raw reads": int(qc["n_reads"]),
        "Structure candidates": int(read_json(s2 / "05_metadata/structure_summary.json")["candidate_found"]),
        "Kept cDNA": int(extract["counts"]["kept"]),
        "MAPQ≥20": int(alignment.get("mapped", 0)),
        "SLAM molecules": int(slam["molecules_total"]),
        "Joined to bins": int(joint["molecules_joined_to_whole_transcriptome_bins"]),
    }
    keys, values = list(stage_values), list(stage_values.values())
    stem = figdir / "04_barcode_recovery"
    fig, ax = plt.subplots(figsize=(6.5, 4.8))
    ax.barh(keys[::-1], values[::-1], color="#4DBBD5")
    ax.set_xscale("log")
    ax.set_xlabel("Reads / molecules (log scale)")
    ax.set_title("Barcode-to-coordinate recovery", loc="left", fontsize=11, fontweight="bold")
    style_axes(ax)
    fig.tight_layout()
    save_figure(fig, stem, args.dpi)
    records.append(figure_record(stem, "04_barcode_recovery", "A", "barcode 回溯链路", [s2 / "05_metadata/extract_summary.json", s3 / "05_metadata/alignment_summary.json", s5 / "05_metadata" / f"{args.sample_id}.slam_summary.json", s7 / "05_metadata/slam_joint_summary.json"], "展示从 raw read 到空间 Bin 的回溯链路。", "当前可追溯对象是空间 Bin/分子，不是 isoform-level 定位。"))

    observed_barcode_count = len(barcode_level)
    stem = figdir / "04_spatial_barcode_count"
    fig, ax = plt.subplots(figsize=(6.5, 4.8))
    ax.bar(["Observed barcode\nsequences"], [observed_barcode_count], color="#3C5488")
    ax.set_ylabel("Barcode count")
    ax.set_title("Observed spatial barcodes", loc="left", fontsize=11, fontweight="bold")
    ax.text(0, observed_barcode_count, f"{observed_barcode_count:,}", ha="center", va="bottom", fontsize=9)
    style_axes(ax)
    fig.tight_layout()
    save_figure(fig, stem, args.dpi)
    records.append(figure_record(stem, "04_spatial_barcode_count", "A", "空间 barcode 数量", [barcode_level_path], "展示 alignment 阶段观测到的 unique barcode 序列数量。", "这是观测 barcode 序列数，不等同于独立 UMI 数。"))

    stem = figdir / "04_umi_count"
    fig, ax = plt.subplots(figsize=(6.5, 4.8))
    umi_labels = ["Total UMI\nmolecules", "Median UMI\nper Bin"]
    umi_values = [umi_total, median_umi]
    ax.bar(umi_labels, umi_values, color=["#3C5488", "#4DBBD5"])
    ax.set_ylabel("UMI count")
    ax.set_title("UMI count", loc="left", fontsize=11, fontweight="bold")
    for index, value in enumerate(umi_values):
        ax.text(index, value, f"{value:,.0f}", ha="center", va="bottom", fontsize=9)
    style_axes(ax)
    fig.tight_layout()
    save_figure(fig, stem, args.dpi)
    records.append(figure_record(stem, "04_umi_count", "A", "bin50 UMI 统计", [calls_path], "按 50 µm Bin 汇总 UMI-level molecule-call 行数，并给出总量和每 Bin 中位数。", "UMI 统计基于当前 SLAM molecule-call 表的 UMI-level 行；未额外进行短读长 UMI 去重。"))

    stem = figdir / "04_gene_count"
    gene_count = int(spatial.get("genes_after_filter", 0))
    median_genes = float(pd.to_numeric(cluster_summary["median_detected_genes"], errors="coerce").median())
    fig, ax = plt.subplots(figsize=(6.5, 4.8))
    ax.bar(["Genes after\nfilter", "Median genes\nper Bin"], [gene_count, median_genes], color=["#00A087", "#4DBBD5"])
    ax.set_ylabel("Gene count")
    ax.set_title("Spatial gene statistics", loc="left", fontsize=11, fontweight="bold")
    for index, value in enumerate([gene_count, median_genes]):
        ax.text(index, value, f"{value:,.0f}", ha="center", va="bottom", fontsize=9)
    style_axes(ax)
    fig.tight_layout()
    save_figure(fig, stem, args.dpi)
    records.append(figure_record(stem, "04_gene_count", "A", "空间 gene 统计", [s6 / "05_metadata/long_read_spatial_clustering.json", cluster_path], "展示过滤后 gene 数和每个空间 Bin 的中位 detected gene 数。", "该统计基于 long-read molecule 聚合结果。"))

    stem = figdir / "04_isoform_count"
    isoform_counts = {
        "Candidate\ntranscripts": int(float(isoform_metrics.get("query_transcripts", structure.get("transcripts", 0)))),
        "Candidate\nloci": int(float(isoform_metrics.get("query_loci", 0))),
        "Multi-exon\ntranscripts": int(float(isoform_metrics.get("multi_exon_transcripts", 0))),
    }
    fig, ax = plt.subplots(figsize=(6.5, 4.8))
    ax.bar(list(isoform_counts), list(isoform_counts.values()), color="#3C5488")
    ax.set_ylabel("Count")
    ax.set_title("Isoform statistics", loc="left", fontsize=11, fontweight="bold")
    for index, value in enumerate(isoform_counts.values()):
        ax.text(index, value, f"{value:,}", ha="center", va="bottom", fontsize=9)
    style_axes(ax)
    fig.tight_layout()
    save_figure(fig, stem, args.dpi)
    records.append(figure_record(stem, "04_isoform_count", "A", "空间全长 isoform 统计", [s4 / "01_tables/isoform_summary.tsv"], "展示 candidate transcript、candidate loci 和 multi-exon transcript 数量。", "当前没有 read-to-isoform 空间定位统计。"))

    assignment = barcode.get("assign_counts", {})
    assignment_keys = ["Exact", "Rescued", "Ambiguous", "Not in whitelist"]
    assignment_values = [
        int(assignment.get("exact", 0)),
        int(assignment.get("rescued", 0)),
        int(assignment.get("ambiguous_barcode", 0)),
        int(assignment.get("barcode_not_in_whitelist", 0)),
    ]
    stem = figdir / "04_barcode_assignment"
    fig, ax = plt.subplots(figsize=(6.5, 4.8))
    ax.bar(assignment_keys, assignment_values, color=["#00A087", "#4DBBD5", "#F39B7F", "#B0B7C3"])
    ax.set_ylabel("Reads")
    ax.set_title("Spatial barcode assignment", loc="left", fontsize=11, fontweight="bold")
    ax.tick_params(axis="x", rotation=25)
    for index, value in enumerate(assignment_values):
        ax.text(index, value, f"{value:,}", ha="center", va="bottom", fontsize=8)
    style_axes(ax)
    fig.tight_layout()
    save_figure(fig, stem, args.dpi)
    records.append(figure_record(stem, "04_barcode_assignment", "A", "空间 barcode 分配", [s2 / "05_metadata/barcode_summary.json"], "展示 exact、rescued、ambiguous 和 whitelist 外 barcode 的分配结果。", "barcode assignment 不是 isoform-level 空间定位。"))

    med_genes = pd.to_numeric(cluster_summary["median_detected_genes"], errors="coerce")
    stem = figdir / "04_gene_coverage"
    fig, ax = plt.subplots(figsize=(6.5, 4.8))
    ax.bar(cluster_summary["cluster"].astype(str), med_genes, color=[CLUSTER_COLORS.get(str(x), "#777777") for x in cluster_summary["cluster"]])
    ax.set_xlabel("Spatial expression cluster")
    ax.set_ylabel("Median detected genes per Bin")
    ax.set_title("Gene coverage by spatial cluster", loc="left", fontsize=11, fontweight="bold")
    style_axes(ax)
    fig.tight_layout()
    save_figure(fig, stem, args.dpi)
    records.append(figure_record(stem, "04_gene_coverage", "A", "空间 gene coverage", [cluster_path], "展示各空间 cluster 的中位 detected gene 数量。", "UMI/saturation 指标在当前结果树中不可用。"))

    # 05: single-panel substitution and SLAM evidence figures.
    substitution = pd.read_csv(s5 / "01_tables/base_substitution_comparison.tsv", sep="\t")
    substitution = substitution.sort_values("rate_percent", ascending=False)
    colors = ["#DC0000" if x in {"T>C", "A>G"} else "#B0B7C3" for x in substitution["channel"]]
    stem = figdir / "05_mismatch_spectrum"
    fig, ax = plt.subplots(figsize=(6.5, 4.8))
    ax.bar(substitution["channel"], substitution["rate_percent"], color=colors)
    ax.set_ylabel("Mismatch rate (%)")
    ax.set_title("Base substitution spectrum", loc="left", fontsize=11, fontweight="bold")
    ax.tick_params(axis="x", rotation=45)
    style_axes(ax)
    fig.tight_layout()
    save_figure(fig, stem, args.dpi)
    records.append(figure_record(stem, "05_mismatch_spectrum", "A", "错配谱", [s5 / "01_tables/base_substitution_comparison.tsv"], "展示全部碱基替换类型，并突出 T>C/A>G 通道。", "没有未标记对照；信号只能作探索性解释。"))

    stem = figdir / "05_slam_evidence_summary"
    fig, ax = plt.subplots(figsize=(6.5, 4.8))
    ax.axis("off")
    ax.set_title("SLAM evidence summary", loc="left", fontsize=11, fontweight="bold")
    ax.text(0, 0.84, "\n".join([
        f"T>C rate: {float(substitution.loc[substitution.channel.eq('T>C'), 'rate_percent'].iloc[0]):.3f}%" if (substitution.channel == "T>C").any() else "T>C rate: NA",
        f"A>G rate: {float(substitution.loc[substitution.channel.eq('A>G'), 'rate_percent'].iloc[0]):.3f}%" if (substitution.channel == "A>G").any() else "A>G rate: NA",
        f"Background rate: {float(mismatch.get('background_rate', 0)):.6f}",
        f"Reads evaluated: {int(mismatch.get('processed_primary_reads', 0)):,}",
    ]), va="top", fontsize=12, linespacing=1.8)
    fig.tight_layout()
    save_figure(fig, stem, args.dpi)
    records.append(figure_record(stem, "05_slam_evidence_summary", "A", "SLAM 错配证据", [s5 / "01_tables/base_substitution_comparison.tsv", s5 / "05_metadata/mismatch_summary.json"], "汇总 T>C、A>G、背景错配率和评估 reads。", "没有未标记对照；不能把转换率当作已校准 nascent RNA。"))

    # 06: single-library compatibility chain.
    stem = figdir / "06_single_library_compatibility"
    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    labels = ["Raw reads", "Kept cDNA", "MAPQ≥20", "SLAM molecules", "Joined bins"]
    values = [stage_values["Raw reads"], stage_values["Kept cDNA"], stage_values["MAPQ≥20"], stage_values["SLAM molecules"], stage_values["Joined to bins"]]
    fractions = [value / values[0] for value in values]
    ax.plot(labels, fractions, marker="o", linewidth=2.5, color="#3C5488")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Fraction of raw reads")
    ax.set_title("Information retention through processing", loc="left", fontsize=11, fontweight="bold")
    ax.tick_params(axis="x", rotation=25)
    style_axes(ax)
    fig.tight_layout()
    save_figure(fig, stem, args.dpi)
    records.append(figure_record(stem, "06_single_library_compatibility", "A", "说明空间与长读长信息链路保留", [qc_path, s2 / "05_metadata/extract_summary.json", s3 / "05_metadata/alignment_summary.json", s7 / "05_metadata/slam_joint_summary.json"], "单文库链路显示信息从读段、cDNA、比对到空间联合结果的保留情况。", "不是 SLAM 化学前后兼容性实验比较。"))

    # 07: single-panel transcript structure figures.
    class_counts = structure.get("class_counts", {})
    known_classes = {"=": "Exact reference match", "c": "Contained in reference"}
    novel_classes = {
        "j": "Novel junction combination",
        "o": "Other same-strand overlap",
        "k": "Contains reference",
        "m": "Retained intron",
        "n": "Reference-intron overlap",
    }
    stem = figdir / "07_known_transcripts"
    fig, ax = plt.subplots(figsize=(6.5, 4.8))
    known_labels = list(known_classes.values())
    known_values = [class_counts.get(code, 0) for code in known_classes]
    ax.bar(known_labels, known_values, color="#4DBBD5")
    ax.set_ylabel("Candidate transcript count")
    ax.set_title("Known transcript classes", loc="left", fontsize=11, fontweight="bold")
    ax.tick_params(axis="x", rotation=25)
    for index, value in enumerate(known_values):
        ax.text(index, value, f"{value:,}", ha="center", va="bottom", fontsize=8)
    style_axes(ax)
    fig.tight_layout()
    save_figure(fig, stem, args.dpi)
    records.append(figure_record(stem, "07_known_transcripts", "B", "Known transcript 分类", [s4 / "05_metadata/isoform_class_summary.json"], "展示 gffcompare 的 exact reference match 和 contained in reference transcript 类别。", "当前输出不是 FSM/ISM/NIC/NNC 分类。"))

    stem = figdir / "07_novel_transcripts"
    fig, ax = plt.subplots(figsize=(6.5, 4.8))
    novel_labels = list(novel_classes.values())
    novel_values = [class_counts.get(code, 0) for code in novel_classes]
    ax.bar(novel_labels, novel_values, color="#E64B35")
    ax.set_ylabel("Candidate transcript count")
    ax.set_title("Novel transcript classes", loc="left", fontsize=11, fontweight="bold")
    ax.tick_params(axis="x", rotation=35)
    for index, value in enumerate(novel_values):
        ax.text(index, value, f"{value:,}", ha="center", va="bottom", fontsize=8)
    style_axes(ax)
    fig.tight_layout()
    save_figure(fig, stem, args.dpi)
    records.append(figure_record(stem, "07_novel_transcripts", "B", "Novel transcript 分类", [s4 / "05_metadata/isoform_class_summary.json"], "展示 gffcompare 的 novel junction、overlap、retained intron 等 transcript 类别。", "当前输出不是 FSM/ISM/NIC/NNC 分类。"))

    structure_summary["transcripts"] = pd.to_numeric(structure_summary["transcripts"], errors="coerce")
    structure_plot = structure_summary.groupby("exon_group", as_index=True)["transcripts"].sum()
    stem = figdir / "07_exon_structure_summary"
    fig, ax = plt.subplots(figsize=(6.5, 4.8))
    ax.bar(structure_plot.index, structure_plot.values, color="#00A087")
    ax.set_ylabel("Transcript count")
    ax.set_title("Exon structure summary", loc="left", fontsize=11, fontweight="bold")
    ax.tick_params(axis="x", rotation=45)
    style_axes(ax)
    fig.tight_layout()
    save_figure(fig, stem, args.dpi)
    records.append(figure_record(stem, "07_exon_structure_summary", "B", "外显子结构统计", [s4 / "01_tables/isoform_structure_summary.tsv"], "展示 multi-exon 和 single-exon transcript 统计。", "该统计不等同于 FSM/ISM/NIC/NNC 分类。"))

    # 08: single-panel total and exploratory new-supported spatial figures.
    total_col = "total_lr_molecules"
    new_col = "new_supported"
    for figure_id, column, title, cmap in (
        ("08_total_long_read_spatial", total_col, "Total long-read molecules", CONTINUOUS),
        ("08_new_supported_spatial", new_col, "New-supported proxy", ACTIVITY),
        ("08_new_supported_fraction_spatial", "new_supported_fraction", "New-supported fraction", ACTIVITY),
    ):
        stem = figdir / figure_id
        fig, ax = plt.subplots(figsize=(6.5, 5.8))
        plot_spatial_values(ax, bins, column, title, cmap, log=(column != "new_supported_fraction"))
        fig.tight_layout()
        save_figure(fig, stem, args.dpi)
        records.append(figure_record(stem, figure_id, "B", title, [bins_path, s6 / "05_metadata/long_read_spatial_clustering.json"], "在 50 µm Bin 上展示长读长总量或 exploratory new-supported 空间分布。", "new-supported 不是有对照校准的 nascent RNA。"))

    # 09: cluster-level new-supported fractions with confidence intervals.
    stem = figdir / "09_domain_new_fraction"
    joint_cluster = pd.read_csv(s7 / "01_tables/cluster_slam_summary.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    x = range(len(joint_cluster))
    y = joint_cluster["new_supported_fraction"]
    lower = y - joint_cluster["wilson95_low"]
    upper = joint_cluster["wilson95_high"] - y
    ax.errorbar(list(x), y, yerr=[lower, upper], fmt="o", color="#DC0000", ecolor="#3C5488", capsize=3)
    ax.set_xticks(list(x), joint_cluster["cluster"].astype(str))
    ax.set_xlabel("Spatial expression cluster")
    ax.set_ylabel("New-supported fraction")
    ax.set_title("New-supported fraction by spatial cluster", loc="left", fontsize=11, fontweight="bold")
    ax.set_ylim(0, max(0.3, float(joint_cluster["wilson95_high"].max()) * 1.15))
    style_axes(ax)
    fig.tight_layout()
    save_figure(fig, stem, args.dpi)
    records.append(figure_record(stem, "09_domain_new_fraction", "B", "把动态 proxy 落到空间 domain", [s7 / "01_tables/cluster_slam_summary.tsv"], "按现有表达簇汇总 new-supported fraction 并保留 Wilson 区间。", "单样本且无未标记对照，不作组间差异或因果解释。"))

    # 10: single-panel candidate gene figures.
    profile_rows = []
    for key, (total, new) in candidate_profile.items():
        if key in bin_lookup.index:
            row = bin_lookup.loc[key]
            profile_rows.append({"bin_id": key, "bin_x": row["bin_x"], "bin_y": row["bin_y"], "total": total, "new": new, "fraction": new / total if total else 0})
    profile = pd.DataFrame(profile_rows)
    if profile.empty:
        profile = load_gene_spatial_profile(h5ad_path, candidate_id)
        if not profile.empty:
            profile["new"] = float("nan")
            profile["fraction"] = float("nan")
    point_size = spatial_bin_marker_size(spatial.get("bin_size_um", args.bin_size))
    map_specs = (
        ("10_candidate_expression", "total", f"{candidate_name} expression", CONTINUOUS, None, "候选 gene 的空间 total expression。", "gene-level total expression 来自 H5AD；没有 gene-specific nascent 调用。"),
        ("10_new_supported_proxy", "new", "New-supported proxy", ACTIVITY, "new_supported", "展示 exploratory new-supported 空间 proxy。", "当前没有候选 gene-specific per-Bin new-supported 调用。"),
        ("10_new_supported_fraction", "fraction", "New-supported fraction", ACTIVITY, "new_supported_fraction", "展示 exploratory new-supported fraction 空间分布。", "当前没有候选 gene-specific per-Bin new-supported 调用。"),
    )
    for figure_id, column, title, cmap, fallback_column, interpretation, limitation in map_specs:
        stem = figdir / figure_id
        fig, ax = plt.subplots(figsize=(6.5, 5.8))
        profile_values = profile[column] if column in profile else pd.Series(dtype=float)
        if not profile_values.dropna().empty:
            low, high = robust_display_limits(profile_values)
            sc = ax.scatter(profile["bin_x"], profile["bin_y"], c=profile[column], s=point_size, marker="s", cmap=cmap, vmin=low, vmax=high, linewidths=0, rasterized=True)
            cb = ax.figure.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
            cb.ax.tick_params(labelsize=7)
            cb.set_label("2nd–98th percentile display", fontsize=8)
        elif fallback_column:
            low, high = robust_display_limits(bins[fallback_column])
            sc = ax.scatter(bins["bin_x"], bins["bin_y"], c=bins[fallback_column], s=point_size, marker="s", cmap=cmap, vmin=low, vmax=high, linewidths=0, rasterized=True)
            cb = ax.figure.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
            cb.ax.tick_params(labelsize=7)
            cb.set_label("2nd–98th percentile display", fontsize=8)
        else:
            ax.text(0.5, 0.5, "Candidate expression unavailable", ha="center", va="center")
        spatial_axes(ax, bins, title)
        fig.tight_layout()
        save_figure(fig, stem, args.dpi)
        records.append(figure_record(stem, figure_id, "C", title, [gene_path, h5ad_path, bins_path], interpretation, limitation))

    stem = figdir / "10_candidate_transcript_structure"
    fig, track = plt.subplots(figsize=(8.5, 5.2))
    track.set_title(f"{candidate_name} transcript structures", loc="left", fontsize=11, fontweight="bold")
    track_colors = {"=": "#3C5488", "c": "#4DBBD5", "j": "#E64B35", "?": "#00A087"}
    all_coordinates = [coordinate for _, start, end in structure_exons for coordinate in (start, end)]
    for row in isoform_tracks:
        for start, end in row["exons"]:
            all_coordinates.extend((start, end))
    if all_coordinates:
        reference_y = len(isoform_tracks) + 1
        if structure_exons:
            starts = [item[1] for item in structure_exons]
            ends = [item[2] for item in structure_exons]
            track.hlines(reference_y, min(starts), max(ends), color="#6B7280", linewidth=1.0)
            for _, start, end in structure_exons:
                track.add_patch(Rectangle((start, reference_y - 0.16), end - start, 0.32, color="#6B7280", alpha=0.85))
        labels = ["Reference"] if structure_exons else []
        for index, row in enumerate(isoform_tracks, start=1):
            y = len(isoform_tracks) + 1 - index
            exons = row["exons"]
            starts = [item[0] for item in exons]
            ends = [item[1] for item in exons]
            color = track_colors.get(str(row["class_code"]), "#00A087")
            track.hlines(y, min(starts), max(ends), color="#9CA3AF", linewidth=0.8)
            for start, end in exons:
                track.add_patch(Rectangle((start, y - 0.16), end - start, 0.32, color=color, alpha=0.9))
            labels.append(f"LR-{index} {row['class_code']} (n={int(row['support']):,})")
        track.set_xlim(min(all_coordinates) - 100, max(all_coordinates) + 100)
        track.set_ylim(0, reference_y + 1)
        track.set_xlabel("Reference genomic coordinate")
        track.set_yticks(list(range(reference_y, 0, -1)), labels)
        track.tick_params(axis="y", labelsize=7)
        track.legend(
            handles=[
                Line2D([0], [0], color="#6B7280", linewidth=5, label="Reference"),
                Line2D([0], [0], color="#E64B35", linewidth=5, label="Novel/other class"),
            ], loc="lower right", fontsize=7, frameon=False,
        )
    else:
        track.axis("off")
        track.text(0.5, 0.5, "Long-read transcript track unavailable", ha="center", va="center")
    style_axes(track)
    fig.tight_layout()
    save_figure(fig, stem, args.dpi)
    records.append(figure_record(stem, "10_candidate_transcript_structure", "C", "候选 gene 长读长 transcript 结构", [gtf_path, isoform_gff3, isoform_tmap], "展示参考 exon 和实际 FLAMES 长读长 transcript 结构。", "当前没有 read-to-isoform 空间定位或 isoform-level nascent 调用。"))

    stem = figdir / "10_candidate_gene_summary"
    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    ax.axis("off")
    ax.set_title(f"{candidate_name} candidate summary", loc="left", fontsize=11, fontweight="bold")
    ax.text(0, 0.84, "\n".join([
        f"Gene ID: {candidate_id or 'NA'}",
        f"Informative molecules: {candidate.get('informative_molecules', 'NA')}",
        f"Observed spatial Bins: {len(profile):,}",
        f"Reference exons: {len(structure_exons):,}",
        f"Observed long-read isoforms in track: {len(isoform_tracks):,}",
        "Gene total: H5AD expression column",
        "Gene new-supported per Bin: unavailable by current outputs",
        "Isoform-level spatial/dynamic assignment: unavailable",
    ]), va="top", fontsize=11, linespacing=1.8)
    fig.tight_layout()
    save_figure(fig, stem, args.dpi)
    records.append(figure_record(stem, "10_candidate_gene_summary", "C", "候选 gene 统计", [gene_path, h5ad_path, gtf_path, isoform_gff3, isoform_tmap], "汇总候选 gene 的 informative molecule、空间 Bin、参考 exon 和长读长 transcript 数量。", "该 gene 仅作为候选案例，不是独立验证后的机制结论。"))

    # Compact traceability tables and explanatory handoff.
    metrics = [
        ("sample_id", args.sample_id, "input_samples.tsv", "sample_id", "observed", "single sample"),
        ("species", "Mus musculus", "input_scientific_scope.json", "species", "observed", ""),
        ("chemistry", "BBV2.4", "input_samples.tsv", "chemistry", "observed", ""),
        ("tissue_type", "unknown", "input_scientific_scope.json", "tissue_type", "missing", "not supplied"),
        ("matched_no4su_control", "false", "input_scientific_scope.json", "matched_no4su_control", "missing", "no control"),
        ("label_duration_known", "false", "input_scientific_scope.json", "label_duration_known", "missing", "no 4sU timing"),
        ("bin_size_um", str(spatial.get("bin_size_um", args.bin_size)), "long_read_spatial_clustering.json", "bin_size_um", "observed", ""),
        ("raw_reads", qc.get("n_reads", "NA"), "qc_metrics.tsv", "n_reads", "observed", ""),
        ("kept_cdna", str(extract["counts"]["kept"]), "extract_summary.json", "counts.kept", "observed", ""),
        ("mapq20_mapped", str(alignment.get("mapped", "NA")), "alignment_summary.json", "mapped", "observed", ""),
        ("candidate_transcripts", str(structure.get("transcripts", "NA")), "isoform_class_summary.json", "transcripts", "observed", ""),
        ("slam_molecules", str(slam.get("molecules_total", "NA")), "slam_summary.json", "molecules_total", "observed", ""),
        ("joined_spatial_molecules", str(joint.get("molecules_joined_to_whole_transcriptome_bins", "NA")), "slam_joint_summary.json", "molecules_joined_to_whole_transcriptome_bins", "observed", ""),
        ("observed_spatial_barcodes", str(len(barcode_level)), "barcode_level_summary.tsv", "rows", "observed", "unique observed barcode sequences"),
        ("spatial_bins", str(len(bins)), "long_read_bin50_clusters.tsv", "rows", "observed", ""),
        ("spatial_genes_after_filter", str(spatial.get("genes_after_filter", "NA")), "long_read_spatial_clustering.json", "genes_after_filter", "observed", ""),
        ("candidate_transcript_loci", str(isoform_metrics.get("query_loci", "NA")), "isoform_summary.tsv", "query_loci", "observed", ""),
        ("umi_count", str(umi_total), f"{args.sample_id}.slam_molecule_calls.tsv.gz", "rows", "observed", "UMI-level molecule-call rows"),
        ("median_umi_per_bin", str(median_umi), f"{args.sample_id}.slam_molecule_calls.tsv.gz", "rows per bin", "observed", "50 µm Bin median"),
        ("median_detected_genes_per_bin", str(pd.to_numeric(cluster_summary["median_detected_genes"], errors="coerce").median()), "long_read_cluster_summary.tsv", "median_detected_genes", "observed", "median across spatial clusters"),
        ("spatial_clusters", str(spatial.get("clusters", "NA")), "long_read_spatial_clustering.json", "clusters", "observed", ""),
        ("isoform_spatial_slam", "unavailable", "workflow/report", "read_to_isoform", "missing", "no read-to-isoform mapping"),
    ]
    with (table_dir / "promo_metrics.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value", "source_file", "source_field", "status", "note"])
        writer.writerows(metrics)
    summary = f"""# {args.sample_id} 宣传图包说明

输出目录：`{outdir}`

本图包基于已完成的 01–07 阶段结果，共 {len(records)} 张独立 PNG/PDF 图，交付目录不含组图。直接复用的阶段图保持源文件分辨率；本次派生图按 {args.dpi} DPI 输出。所有空间图使用 50 µm Bin 中心和原始坐标方向，不做插值或平滑，颜色极值按显示层 2%–98% 分位数处理。

## 可以用于宣传的结论

- 长读长数据具备明确的读长、N50、质量值、比对和空间 barcode 回溯证据。
- 空间 barcode、bin50 UMI、gene 和 isoform 已拆成独立统计图；UMI 直接来自 molecule-call 表的 UMI-level 行。
- 50 µm 空间 Bin 上可观察到长读长表达簇、marker 和空间连续性。
- FLAMES/gffcompare 可解析参考匹配和候选转录本结构类别。
- Known transcript 和 Novel transcript 已分别输出为独立柱状图。
- 候选 gene 图使用 FLAMES GFF3 的实际长读长 transcript/exon 结构，并与参考 exon track 对照。
- T>C/A>G 及全部错配谱已被汇总，并可与背景错配率同时展示。

## 必须保留的限制

- 当前为单样本分析；组织类型、4sU 剂量、标记时间、取材时间、重复和切片信息未提供。
- 没有 matched unlabelled/no-4sU control，因此 `new_supported` 只能称为探索性动态 proxy，不能当作已校准 nascent RNA。
- 没有 H&E，图 01 是坐标型空间结构图，不是组织学叠加图。
- 当前结果没有 read-to-isoform 映射，不能输出可靠的 isoform-level total/nascent 空间图或完整机制结论。
- 候选 gene 图是 gene-level 展示案例，不是独立验证后的机制结论。
"""
    (outdir / "promo_summary.zh-CN.md").write_text(summary, encoding="utf-8")
    write_manifest(outdir, records)
    source_paths = [qc_path, bins_path, cluster_path, gene_path, calls_path]
    with (metadata_dir / "input_manifest.tsv").open("w", encoding="utf-8") as handle:
        handle.write("role\tpath\tsize_bytes\n")
        for path in source_paths:
            if path.is_file():
                handle.write(f"source\t{path.resolve()}\t{path.stat().st_size}\n")
    outputs = sorted(path for path in outdir.rglob("*") if path.is_file() and path.name not in {"promo_validation.json", "stage_complete.json", "output_manifest.tsv"})
    validation = {
        "ok": len(records) == 29 and all((figdir / row["filename_png"]).is_file() and (figdir / row["filename_png"]).stat().st_size > 0 for row in records),
        "sample_id": args.sample_id,
        "figure_count": len(records),
        "candidate_gene": candidate_name,
        "limitations": ["no_h_and_e", "no_matched_unlabelled_control", "no_replicates", "no_read_to_isoform_mapping"],
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    (metadata_dir / "promo_validation.json").write_text(json.dumps(validation, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    outputs = sorted(path for path in outdir.rglob("*") if path.is_file() and path.name != "output_manifest.tsv")
    with (metadata_dir / "output_manifest.tsv").open("w", encoding="utf-8") as handle:
        handle.write("path\tsize_bytes\n")
        for path in outputs:
            handle.write(f"{path.resolve()}\t{path.stat().st_size}\n")
    (metadata_dir / "stage_complete.json").write_text(json.dumps({"stage": "09_promo_figures", "sample_id": args.sample_id, "completion_mode": "derived_from_completed_results", "outputs": [str(path.resolve()) for path in outputs]}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(validation, ensure_ascii=False))


if __name__ == "__main__":
    main()
