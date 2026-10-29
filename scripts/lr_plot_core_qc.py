#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from lr_figure_style import apply_style, clean_axis, save_figure

PALETTE = {
    "raw": "#3C5488",
    "status": "#4DBBD5",
    "accent": "#E64B35",
    "good": "#00A087",
    "warn": "#DC0000",
    "neutral": "#91D1C2",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot core QC and summary figures.")
    parser.add_argument("--section", choices=("read", "barcode", "align"), required=True)
    parser.add_argument("--qc-per-read")
    parser.add_argument("--structure-tsv")
    parser.add_argument("--barcode-tsv")
    parser.add_argument("--structure-summary")
    parser.add_argument("--barcode-summary")
    parser.add_argument("--stats-dir")
    parser.add_argument("--barcode-level-summary")
    parser.add_argument("--alignment-summary")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--dpi", type=int, default=150)
    parser.add_argument("--font-family", default="DejaVu Sans")
    parser.add_argument("--base-size", type=float, default=10.5)
    parser.add_argument(
        "--max-hist-rows",
        type=int,
        default=500000,
        help="Maximum rows retained in memory for histogram-style plots. Counters still use all rows.",
    )
    return parser.parse_args()


def read_gz_tsv(path: str):
    with gzip.open(path, "rt") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def plot_hist(values, title, xlabel, outpath, bins=50):
    if not values:
        return
    fig, ax = plt.subplots(figsize=(7.0, 4.3))
    ax.hist(values, bins=bins, color=PALETTE["raw"], edgecolor="white", alpha=0.9)
    ax.set_title(title, loc="left")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Reads")
    clean_axis(ax, "y")
    fig.tight_layout()
    save_figure(fig, os.path.splitext(outpath)[0], int(plt.rcParams["savefig.dpi"]), pdf=False)
    print(f"[PLOT] {outpath}", flush=True)


def plot_counter(counter, title, outpath, annotate_pct=True, order=None, horizontal=True, color=None):
    if not counter:
        return
    keys = order if order is not None else list(counter.keys())
    vals = [counter.get(k, 0) for k in keys]
    total = sum(vals)
    color = color or PALETTE["status"]
    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    if horizontal:
        bars = ax.barh(keys, vals, color=color, edgecolor="none", alpha=0.92)
        ax.invert_yaxis()
        ax.set_xlabel("Count")
        ax.set_ylabel("")
        if annotate_pct and total:
            xmax = max(vals) if vals else 0
            for bar, value in zip(bars, vals):
                pct = value / total * 100
                ax.text(
                    bar.get_width() + max(xmax * 0.015, 1),
                    bar.get_y() + bar.get_height() / 2,
                    f"{value:,} ({pct:.1f}%)",
                    va="center",
                    ha="left",
                    fontsize=9,
                )
            ax.set_xlim(0, max(vals) * 1.24 if max(vals) else 1)
    else:
        bars = ax.bar(keys, vals, color=color, edgecolor="none", alpha=0.92)
        ax.set_ylabel("Count")
        ax.tick_params(axis="x", rotation=30)
        if annotate_pct and total:
            ymax = max(vals) if vals else 0
            for bar, value in zip(bars, vals):
                pct = value / total * 100
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + max(ymax * 0.015, 1),
                    f"{value:,}\n{pct:.1f}%",
                    ha="center",
                    va="bottom",
                    fontsize=9,
                )
    ax.set_title(title, loc="left")
    clean_axis(ax, "x" if horizontal else "y")
    fig.tight_layout()
    save_figure(fig, os.path.splitext(outpath)[0], int(plt.rcParams["savefig.dpi"]), pdf=False)
    print(f"[PLOT] {outpath}", flush=True)


def append_limited(values, value, limit):
    if limit <= 0 or len(values) < limit:
        values.append(value)


def plot_barcode_status_for_source_counts(counts, source, outpath):
    categories = ["exact", "mismatch_1", "mismatch_2", "ambiguous_barcode", "barcode_not_in_whitelist"]
    values = [counts.get(key, 0) for key in categories]
    total = sum(values)
    if total == 0:
        return
    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    bars = ax.barh(categories, values, color=PALETTE["neutral"], edgecolor="none", alpha=0.92)
    ax.invert_yaxis()
    ax.set_title(f"Barcode assignment by structure source: {source}", loc="left")
    ax.set_xlabel("Count")
    ax.set_ylabel("")
    ymax = max(values)
    for bar, value in zip(bars, values):
        pct = value / total * 100
        ax.text(
            bar.get_width() + max(ymax * 0.015, 1),
            bar.get_y() + bar.get_height() / 2,
            f"{pct:.1f}%",
            ha="left",
            va="center",
            fontsize=9,
        )
    ax.set_xlim(0, max(values) * 1.20 if max(values) else 1)
    clean_axis(ax, "x")
    fig.tight_layout()
    save_figure(fig, os.path.splitext(outpath)[0], int(plt.rcParams["savefig.dpi"]), pdf=False)
    print(f"[PLOT] {outpath}", flush=True)


def main() -> None:
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    apply_style(args.font_family, args.base_size)
    plt.rcParams["figure.dpi"] = args.dpi
    plt.rcParams["savefig.dpi"] = args.dpi

    if args.section == "read":
        if not args.qc_per_read:
            raise SystemExit("--qc-per-read is required for --section read")
        read_lengths = []
        read_q = []
        for row in read_gz_tsv(args.qc_per_read):
            append_limited(read_lengths, int(row["read_len"]), args.max_hist_rows)
            append_limited(read_q, float(row["mean_q"]), args.max_hist_rows)
        plot_hist(read_lengths, "Raw Read Length", "Length", os.path.join(args.outdir, "raw_read_length.png"))
        plot_hist(read_q, "Raw Mean Q", "Mean Q", os.path.join(args.outdir, "raw_mean_q.png"))
        return

    if args.section == "barcode":
        compact_mode = bool(args.structure_summary and args.barcode_summary)
        if compact_mode:
            if not args.structure_tsv:
                raise SystemExit("--structure-tsv is required to sample cDNA lengths")
            with open(args.structure_summary, encoding="utf-8") as handle:
                structure_summary = json.load(handle)
            with open(args.barcode_summary, encoding="utf-8") as handle:
                barcode_summary = json.load(handle)
            structure_counts = {
                key: int(structure_summary.get(key, 0))
                for key in ("candidate_found", "short_cdna", "anchor_not_found")
            }
            source_counts = structure_summary.get("source_counts", {})
            anchor_linker_counts = {
                "anchor_not_found": structure_counts["anchor_not_found"],
                "short_cdna": structure_counts["short_cdna"],
                "anchor_5p_pass": int(source_counts.get("anchor_5p", 0)),
                "linker1_fallback_pass": int(source_counts.get("linker1_fallback", 0)),
            }
            cdna_lengths = []
            for row in read_gz_tsv(args.structure_tsv):
                if row.get("scan_status") == "candidate_found":
                    append_limited(cdna_lengths, int(row["cdna_len"]), args.max_hist_rows)
                    if args.max_hist_rows > 0 and len(cdna_lengths) >= args.max_hist_rows:
                        break
            assign_counts = barcode_summary.get("assign_counts", {})
            barcode_counts = {
                "exact": int(assign_counts.get("exact", 0)),
                "mismatch_1": int(assign_counts.get("mismatch_1", assign_counts.get("rescued", 0))),
                "mismatch_2": int(assign_counts.get("mismatch_2", 0)),
                "ambiguous_barcode": int(assign_counts.get("ambiguous_barcode", 0)),
                "barcode_not_in_whitelist": int(assign_counts.get("barcode_not_in_whitelist", 0)),
            }
            edit_distance_counts = {
                str(key): int(value)
                for key, value in barcode_summary.get("edit_distance_counts", {}).items()
            }
            plot_counter(anchor_linker_counts, "Anchor/Linker Match Status", os.path.join(args.outdir, "anchor_linker_status.png"), order=["anchor_5p_pass", "linker1_fallback_pass", "short_cdna", "anchor_not_found"], color=PALETTE["accent"])
            plot_counter(structure_counts, "Structure Scan Status", os.path.join(args.outdir, "structure_status.png"), order=["candidate_found", "short_cdna", "anchor_not_found"], color=PALETTE["status"])
            plot_hist(cdna_lengths, "Valid cDNA Length (sampled)", "cDNA length", os.path.join(args.outdir, "cdna_length.png"))
            plot_counter(barcode_counts, "Barcode Assignment Status", os.path.join(args.outdir, "barcode_status.png"), order=["exact", "mismatch_1", "mismatch_2", "ambiguous_barcode", "barcode_not_in_whitelist"], color=PALETTE["good"])
            plot_counter(edit_distance_counts, "Barcode Edit Distance", os.path.join(args.outdir, "barcode_edit_distance.png"), order=sorted(edit_distance_counts, key=lambda x: int(x)), color=PALETTE["neutral"])
            return
        if not args.structure_tsv or not args.barcode_tsv:
            raise SystemExit(
                "--section barcode requires either structure/barcode tables or compact summary inputs"
            )
        structure_counts = {}
        structure_source_counts = {}
        cdna_lengths = []
        for row in read_gz_tsv(args.structure_tsv):
            structure_counts[row["scan_status"]] = structure_counts.get(row["scan_status"], 0) + 1
            if row["scan_status"] == "candidate_found":
                source = row.get("structure_source", "anchor_5p")
                structure_source_counts[source] = structure_source_counts.get(source, 0) + 1
                append_limited(cdna_lengths, int(row["cdna_len"]), args.max_hist_rows)
        anchor_linker_counts = {
            "anchor_not_found": structure_counts.get("anchor_not_found", 0),
            "short_cdna": structure_counts.get("short_cdna", 0),
            "anchor_5p_pass": structure_source_counts.get("anchor_5p", 0),
            "linker1_fallback_pass": structure_source_counts.get("linker1_fallback", 0),
        }
        plot_counter(anchor_linker_counts, "Anchor/Linker Match Status", os.path.join(args.outdir, "anchor_linker_status.png"), order=["anchor_5p_pass", "linker1_fallback_pass", "short_cdna", "anchor_not_found"], color=PALETTE["accent"])
        plot_counter(structure_counts, "Structure Scan Status", os.path.join(args.outdir, "structure_status.png"), order=["candidate_found", "short_cdna", "anchor_not_found"], color=PALETTE["status"])
        plot_hist(cdna_lengths, "Valid cDNA Length", "cDNA length", os.path.join(args.outdir, "cdna_length.png"))

        barcode_counts = {"exact": 0, "mismatch_1": 0, "mismatch_2": 0, "ambiguous_barcode": 0, "barcode_not_in_whitelist": 0}
        confidence_counts = {"A": 0, "B": 0, "C": 0, "D": 0}
        edit_distance_counts = {}
        source_barcode_counts = {
            "anchor_5p": {key: 0 for key in barcode_counts},
            "linker1_fallback": {key: 0 for key in barcode_counts},
        }
        for row in read_gz_tsv(args.barcode_tsv):
            assign_status = row["assign_status"]
            confidence = row.get("confidence_class", "reject")
            if confidence in confidence_counts:
                confidence_counts[confidence] += 1
            barcode_key = None
            if assign_status == "exact":
                barcode_key = "exact"
            elif assign_status == "rescued":
                barcode_key = f"mismatch_{row['whitelist_ed']}" if row["whitelist_ed"] in {"1", "2"} else None
            elif assign_status in {"ambiguous_barcode", "barcode_not_in_whitelist"}:
                barcode_key = assign_status
            if barcode_key is not None:
                barcode_counts[barcode_key] += 1
                source = row.get("structure_source", ".")
                if row.get("scan_status") == "candidate_found" and source in source_barcode_counts:
                    source_barcode_counts[source][barcode_key] += 1
            if row["whitelist_ed"] != "-1":
                edit_distance_counts[row["whitelist_ed"]] = edit_distance_counts.get(row["whitelist_ed"], 0) + 1
        plot_counter(barcode_counts, "Barcode Assignment Status", os.path.join(args.outdir, "barcode_status.png"), order=["exact", "mismatch_1", "mismatch_2", "ambiguous_barcode", "barcode_not_in_whitelist"], color=PALETTE["good"])
        plot_counter(confidence_counts, "Barcode Confidence Class", os.path.join(args.outdir, "confidence_class.png"), order=["A", "B", "C", "D"], color=PALETTE["accent"])
        plot_counter(edit_distance_counts, "Barcode Edit Distance", os.path.join(args.outdir, "barcode_edit_distance.png"), order=sorted(edit_distance_counts, key=lambda x: int(x)), color=PALETTE["neutral"])
        plot_barcode_status_for_source_counts(source_barcode_counts["anchor_5p"], "anchor_5p", os.path.join(args.outdir, "barcode_status_anchor_5p.png"))
        plot_barcode_status_for_source_counts(source_barcode_counts["linker1_fallback"], "linker1_fallback", os.path.join(args.outdir, "barcode_status_linker1_fallback.png"))
        return

    if not args.stats_dir and not (args.barcode_level_summary and args.alignment_summary):
        raise SystemExit("--stats-dir or both compact alignment inputs are required for --section align")
    barcode_level = args.barcode_level_summary or os.path.join(args.stats_dir, "barcode_level_summary.tsv")
    if os.path.exists(barcode_level):
        totals, mapped = [], []
        with open(barcode_level, "r", encoding="utf-8") as handle:
            for index, row in enumerate(csv.DictReader(handle, delimiter="\t")):
                if index >= 50:
                    break
                totals.append(int(row["reads_total"]))
                mapped.append(int(row["reads_mapped"]))
        if totals:
            fig, ax = plt.subplots(figsize=(7.4, 4.4))
            xs = list(range(len(totals)))
            ax.plot(xs, totals, label="total", color=PALETTE["raw"], linewidth=2.4, marker="o", markersize=3)
            ax.plot(xs, mapped, label="mapped", color=PALETTE["warn"], linewidth=2.4, marker="o", markersize=3)
            ax.fill_between(xs, totals, color=PALETTE["raw"], alpha=0.08)
            ax.set_title("Top 50 barcodes", loc="left")
            ax.set_xlabel("Barcode rank")
            ax.set_ylabel("Reads")
            clean_axis(ax, "y")
            ax.legend(frameon=False, loc="upper right", ncol=2)
            fig.tight_layout()
            save_figure(fig, os.path.join(args.outdir, "top50_barcodes"), int(plt.rcParams["savefig.dpi"]), pdf=False)
    alignment_json = args.alignment_summary or os.path.join(args.stats_dir, "alignment_summary.json")
    if os.path.exists(alignment_json):
        with open(alignment_json, "r", encoding="utf-8") as handle:
            plot_counter(json.load(handle), "Alignment Summary", os.path.join(args.outdir, "alignment_summary.png"), order=["mapped", "unmapped", "secondary", "supplementary"], color=PALETTE["accent"])


if __name__ == "__main__":
    main()
