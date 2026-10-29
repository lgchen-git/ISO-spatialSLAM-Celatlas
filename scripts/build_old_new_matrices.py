#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

import pandas as pd
from scipy import sparse
from scipy.io import mmwrite


def write_matrix(frame, out_prefix, features, barcodes):
    feature_index = {gene: idx for idx, gene in enumerate(features)}
    barcode_index = {barcode: idx for idx, barcode in enumerate(barcodes)}

    if frame.empty:
        matrix = sparse.coo_matrix((len(features), len(barcodes))).tocsr()
    else:
        rows = frame["gene_id"].map(feature_index).to_numpy()
        cols = frame["barcode"].map(barcode_index).to_numpy()
        data = frame["umi_count"].to_numpy()
        matrix = sparse.coo_matrix((data, (rows, cols)), shape=(len(features), len(barcodes))).tocsr()
    mmwrite(f"{out_prefix}.mtx", matrix)


def write_bin50_bundle(frame, outdir, features, feature_names, bin_size):
    matrix_dir = outdir / "filtered_feature_bc_matrix"
    spatial_dir = matrix_dir / "spatial"
    matrix_dir.mkdir(parents=True, exist_ok=True)
    spatial_dir.mkdir(parents=True, exist_ok=True)

    frame = frame[(frame["x"] >= 0) & (frame["y"] >= 0)].copy()
    if frame.empty:
        bin_rows = pd.DataFrame(columns=["bin_id", "spatial1", "spatial2", "umi_count", "gene_count"])
        counts = pd.DataFrame(columns=["gene_id", "bin_id", "umi_count"])
    else:
        frame["bin_x"] = (frame["x"] // bin_size).astype(int) * bin_size
        frame["bin_y"] = (frame["y"] // bin_size).astype(int) * bin_size
        frame["bin_id"] = frame["bin_x"].astype(str) + "_" + frame["bin_y"].astype(str)
        frame["molecule_id"] = (
            frame["barcode"].astype(str)
            + "|"
            + frame["umi"].astype(str)
            + "|"
            + frame["gene_id"].astype(str)
        )
        counts = (
            frame.groupby(["gene_id", "bin_id"], as_index=False)["molecule_id"]
            .nunique()
            .rename(columns={"molecule_id": "umi_count"})
        )
        bin_rows = (
            frame.groupby(["bin_id", "bin_x", "bin_y"], as_index=False)
            .agg(umi_count=("molecule_id", "nunique"), gene_count=("gene_id", "nunique"))
            .sort_values(["bin_y", "bin_x"])
        )
        bin_rows = bin_rows.rename(columns={"bin_x": "spatial1", "bin_y": "spatial2"})

    bin_ids = bin_rows["bin_id"].astype(str).tolist()
    feature_index = {gene: idx for idx, gene in enumerate(features)}
    barcode_index = {barcode: idx for idx, barcode in enumerate(bin_ids)}
    if counts.empty or not features or not bin_ids:
        matrix = sparse.coo_matrix((len(features), len(bin_ids)), dtype="int32")
    else:
        counts = counts[counts["gene_id"].isin(feature_index) & counts["bin_id"].isin(barcode_index)]
        rows = counts["gene_id"].map(feature_index).to_numpy()
        cols = counts["bin_id"].map(barcode_index).to_numpy()
        data = counts["umi_count"].astype("int32").to_numpy()
        matrix = sparse.coo_matrix((data, (rows, cols)), shape=(len(features), len(bin_ids)), dtype="int32")
    mmwrite(str(matrix_dir / "matrix.mtx"), matrix)

    with (matrix_dir / "features.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        for gene in features:
            writer.writerow([gene, feature_names.get(gene, gene), "Gene Expression"])

    with (matrix_dir / "barcodes.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        for bin_id in bin_ids:
            writer.writerow([bin_id])

    with (matrix_dir / "metadata.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["bin_id", "spatial1", "spatial2", "umi_count", "gene_count"],
            delimiter="\t",
        )
        writer.writeheader()
        for row in bin_rows.to_dict("records"):
            writer.writerow(
                {
                    "bin_id": row["bin_id"],
                    "spatial1": int(row["spatial1"]),
                    "spatial2": int(row["spatial2"]),
                    "umi_count": int(row["umi_count"]),
                    "gene_count": int(row["gene_count"]),
                }
            )

    with (spatial_dir / "tissue_positions_list.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        for idx, row in enumerate(bin_rows.to_dict("records")):
            writer.writerow([row["bin_id"], 1, idx, idx, int(row["spatial2"]), int(row["spatial1"])])

    metrics = {
        "matrix_bins": matrix.shape[1],
        "matrix_genes": matrix.shape[0],
        "matrix_nonzero": int(matrix.nnz),
        "molecules_with_spatial_coordinates": int(len(frame)),
    }
    with (outdir / "matrix_metrics.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        for key, value in metrics.items():
            writer.writerow([key, value])

    return metrics


def main():
    parser = argparse.ArgumentParser(description="Build old/new sparse matrices from molecule calls.")
    parser.add_argument("--molecule-calls", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--bin-size", type=int, default=50)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    molecule_df = pd.read_csv(args.molecule_calls, sep="\t")
    features = sorted(molecule_df["gene_id"].dropna().unique())
    feature_names = (
        molecule_df[["gene_id", "gene_name"]]
        .dropna()
        .drop_duplicates("gene_id")
        .set_index("gene_id")["gene_name"]
        .to_dict()
    )
    barcodes = sorted(molecule_df["barcode"].dropna().unique())

    grouped = (
        molecule_df.groupby(["assignment", "gene_id", "barcode"], as_index=False)
        .size()
        .rename(columns={"size": "umi_count"})
    )
    for assignment in ("new", "old", "ambiguous"):
        sub = grouped[grouped["assignment"] == assignment].copy()
        write_matrix(sub, str(outdir / f"{args.sample}.{assignment}_umi"), features, barcodes)

    bin_root = outdir / f"bin{args.bin_size}"
    bundle_inputs = {
        "total": molecule_df,
        "new": molecule_df[molecule_df["assignment"] == "new"].copy(),
        "old": molecule_df[molecule_df["assignment"] == "old"].copy(),
        "ambiguous": molecule_df[molecule_df["assignment"] == "ambiguous"].copy(),
    }
    summary = []
    for label, frame in bundle_inputs.items():
        metrics = write_bin50_bundle(frame, bin_root / label, features, feature_names, args.bin_size)
        summary.append({"matrix": label, **metrics})
    pd.DataFrame(summary).to_csv(bin_root / "matrix_summary.tsv", sep="\t", index=False)

    feature_meta = (
        molecule_df[["gene_id", "gene_name"]]
        .drop_duplicates()
        .sort_values(["gene_id", "gene_name"])
    )
    feature_meta.to_csv(outdir / "features.tsv", sep="\t", index=False)
    Path(outdir / f"{args.sample}.features.tsv").write_text(
        "\n".join(features) + ("\n" if features else ""),
        encoding="utf-8",
    )

    barcode_meta = (
        molecule_df[["barcode", "x", "y"]]
        .drop_duplicates()
        .sort_values(["barcode"])
    )
    barcode_meta.to_csv(outdir / "barcodes.tsv", sep="\t", index=False)
    Path(outdir / f"{args.sample}.barcodes.tsv").write_text(
        "\n".join(barcodes) + ("\n" if barcodes else ""),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
