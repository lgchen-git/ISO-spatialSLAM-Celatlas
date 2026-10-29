#!/usr/bin/env python3
"""Build reproducible bin-level spatial QC tables from molecule calls."""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--molecule-calls", required=True, type=Path)
    parser.add_argument("--clusters", required=True, type=Path)
    parser.add_argument("--isoform-summary", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--bin-size", required=True, type=int)
    parser.add_argument("--chunk-rows", default=250000, type=int)
    parser.add_argument("--max-rows", default=0, type=int, help="Read at most this many rows for a bounded validation run.")
    return parser.parse_args()


def bin_centers(values: pd.Series, bin_size: int) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    return (np.floor(numeric / bin_size) * bin_size + bin_size // 2).astype("Int64")


def main() -> None:
    args = parse_args()
    outdir = args.outdir
    tables = outdir / "01_tables"
    metadata = outdir / "05_metadata"
    tables.mkdir(parents=True, exist_ok=True)
    metadata.mkdir(parents=True, exist_ok=True)

    clusters = pd.read_csv(args.clusters, sep="\t", dtype={"bin_id": str})
    required_clusters = {"bin_id", "bin_x", "bin_y", "cluster"}
    missing = sorted(required_clusters - set(clusters.columns))
    if missing:
        raise SystemExit(f"Missing cluster columns: {', '.join(missing)}")
    clusters = clusters[["bin_id", "bin_x", "bin_y", "cluster"]].copy()
    clusters["bin_id"] = clusters["bin_id"].astype(str)
    clusters = clusters.drop_duplicates("bin_id")
    bin_rows: Counter[str] = Counter()
    new_supported: Counter[str] = Counter()
    no_conversion: Counter[str] = Counter()
    ambiguous: Counter[str] = Counter()
    bin_codes: dict[str, int] = {}
    umi_pairs: list[np.ndarray] = []
    gene_pairs: list[np.ndarray] = []
    total_gene_set: set[str] = set()
    usecols = ["barcode", "x", "y", "gene_id", "umi"]
    optional = ["exploratory_class", "final_call"]
    header = pd.read_csv(args.molecule_calls, sep="\t", compression="gzip", nrows=0).columns.tolist()
    usecols.extend(column for column in optional if column in header)
    rows_read = 0
    for chunk in pd.read_csv(
        args.molecule_calls,
        sep="\t",
        compression="gzip",
        usecols=usecols,
        dtype={"barcode": str, "gene_id": str, "umi": str},
        chunksize=args.chunk_rows,
    ):
        if args.max_rows > 0:
            remaining = args.max_rows - rows_read
            if remaining <= 0:
                break
            chunk = chunk.head(remaining)
        rows_read += len(chunk)
        chunk["x"] = pd.to_numeric(chunk["x"], errors="coerce")
        chunk["y"] = pd.to_numeric(chunk["y"], errors="coerce")
        chunk = chunk.dropna(subset=["x", "y", "gene_id", "umi", "barcode"])
        if chunk.empty:
            continue
        bx = bin_centers(chunk["x"], args.bin_size)
        by = bin_centers(chunk["y"], args.bin_size)
        chunk["bin_id"] = bx.astype(str) + "_" + by.astype(str)
        for bin_id in chunk["bin_id"].astype(str).unique():
            if bin_id not in bin_codes:
                bin_codes[bin_id] = len(bin_codes)
        chunk["bin_code"] = chunk["bin_id"].map(bin_codes).astype(np.int32)
        classes = chunk.get("exploratory_class", chunk.get("final_call", pd.Series("", index=chunk.index))).astype(str)
        for bin_id, group in chunk.groupby("bin_id", sort=False):
            bin_id = str(bin_id)
            bin_rows[bin_id] += len(group)
            genes = group["gene_id"].astype(str)
            total_gene_set.update(genes)
            umi_frame = group[["bin_code", "barcode", "umi"]].astype(str)
            umi_hash = pd.util.hash_pandas_object(umi_frame, index=False).to_numpy(dtype="uint64")
            umi_pairs.append(np.rec.fromarrays([group["bin_code"].to_numpy(dtype=np.int32), umi_hash], names="bin,hash"))
            gene_frame = group[["bin_code", "gene_id"]].astype(str)
            gene_hash = pd.util.hash_pandas_object(gene_frame, index=False).to_numpy(dtype="uint64")
            gene_pairs.append(np.rec.fromarrays([group["bin_code"].to_numpy(dtype=np.int32), gene_hash], names="bin,hash"))
            group_classes = classes.loc[group.index]
            new_supported[bin_id] += int(group_classes.eq("new_supported").sum())
            no_conversion[bin_id] += int(group_classes.eq("no_conversion_detected").sum())
            ambiguous[bin_id] += int(group_classes.eq("ambiguous").sum())

    def unique_counts_by_bin(parts: list[np.ndarray]) -> np.ndarray:
        if not parts:
            return np.zeros(len(bin_codes), dtype=np.int64)
        values = np.concatenate(parts)
        order = np.lexsort((values["hash"], values["bin"]))
        values = values[order]
        unique = np.ones(len(values), dtype=bool)
        if len(values) > 1:
            unique[1:] = (values["bin"][1:] != values["bin"][:-1]) | (values["hash"][1:] != values["hash"][:-1])
        return np.bincount(values["bin"][unique], minlength=len(bin_codes))

    umi_counts = unique_counts_by_bin(umi_pairs)
    gene_counts = unique_counts_by_bin(gene_pairs)
    code_to_bin = {code: bin_id for bin_id, code in bin_codes.items()}
    clusters["umi_count"] = clusters["bin_id"].map(lambda value: int(umi_counts[bin_codes[value]]) if value in bin_codes else 0)
    clusters["lr_molecule_count"] = clusters["bin_id"].map(lambda value: bin_rows.get(value, 0))
    clusters["detected_gene_count"] = clusters["bin_id"].map(lambda value: int(gene_counts[bin_codes[value]]) if value in bin_codes else 0)
    clusters["new_supported"] = clusters["bin_id"].map(lambda value: new_supported.get(value, 0))
    clusters["no_conversion_detected"] = clusters["bin_id"].map(lambda value: no_conversion.get(value, 0))
    clusters["ambiguous"] = clusters["bin_id"].map(lambda value: ambiguous.get(value, 0))
    denominator = clusters["new_supported"] + clusters["no_conversion_detected"]
    clusters["new_supported_fraction"] = np.where(denominator > 0, clusters["new_supported"] / denominator, np.nan)
    columns = [
        "bin_id", "bin_x", "bin_y", "cluster", "umi_count", "lr_molecule_count",
        "detected_gene_count", "new_supported", "no_conversion_detected", "ambiguous",
        "new_supported_fraction",
    ]
    clusters[columns].to_csv(tables / "spatial_bin_qc.tsv", sep="\t", index=False)

    isoform_metrics = {}
    if args.isoform_summary.is_file():
        isoform_frame = pd.read_csv(args.isoform_summary, sep="\t")
        if {"metric", "value"}.issubset(isoform_frame.columns):
            isoform_metrics = dict(zip(isoform_frame["metric"].astype(str), isoform_frame["value"]))
    summary = {
        "bin_size": args.bin_size,
        "rows_scanned": int(rows_read),
        "valid_bin_count": int(len(clusters)),
        "bins_with_molecules": int((clusters["lr_molecule_count"] > 0).sum()),
        "total_lr_molecules": int(clusters["lr_molecule_count"].sum()),
        "total_umi_count": int(clusters["umi_count"].sum()),
        "median_umi_per_bin": float(clusters["umi_count"].median()),
        "median_lr_molecules_per_bin": float(clusters["lr_molecule_count"].median()),
        "median_detected_genes_per_bin": float(clusters["detected_gene_count"].median()),
        "total_detected_genes": int(len(total_gene_set)),
        "candidate_isoforms": int(float(isoform_metrics.get("query_transcripts", 0))),
        "candidate_loci": int(float(isoform_metrics.get("query_loci", 0))),
        "multi_exon_isoforms": int(float(isoform_metrics.get("multi_exon_transcripts", 0))),
        "umi_definition": "unique hash(barcode, umi) within each bin",
        "isoform_definition": "sample-level FLAMES candidate transcripts; no per-bin read-to-isoform mapping",
    }
    (tables / "spatial_qc_summary.tsv").write_text(
        "metric\tvalue\n" + "\n".join(f"{key}\t{value}" for key, value in summary.items()) + "\n",
        encoding="utf-8",
    )
    (metadata / "spatial_qc_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main()
