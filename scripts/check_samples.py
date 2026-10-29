#!/usr/bin/env python3
"""Validate generic and LR-Celatlas sample-sheet contracts."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import shutil
import sys


GENERIC_REQUIRED = ("sample_id", "group", "input_1", "input_2")
LR_REQUIRED = (
    "sample_id",
    "fastq",
    "chip_id",
    "whitelist_csv",
    "species",
    "chemistry",
    "genome_fa",
    "gtf",
)
LR_PATH_COLUMNS = (
    "fastq",
    "whitelist_csv",
    "genome_fa",
    "gtf",
)
LR_OPTIONAL_METADATA = (
    "tissue_type", "experimental_group", "label_4su_dose", "label_4su_duration",
    "harvest_time", "unlabelled_control", "biological_replicate", "section_id",
    "tissue_image",
)
LR_CONFIG_REQUIRED = (
    "project_name",
    "run_id",
    "samples_tsv",
    "data_root",
    "results_root",
    "tmp_root",
    "run_mode",
    "species",
    "chemistry",
    "adapter_python",
    "analysis_python",
    "spatial_python",
    "minimap2_bin",
    "samtools_bin",
    "flames_python_dir",
    "gffcompare_bin",
    "flames_min_support",
    "flames_min_count_fraction",
    "flames_min_support_fraction",
    "flames_strand_specific",
    "flames_downsample_ratio",
    "min_mapq",
    "min_tc_baseq",
    "mismatch_max_reads",
    "min_t_covered",
    "min_tc_count",
    "slam_fdr",
    "max_site_other_mismatch_rate",
    "max_read_other_mismatch_rate",
    "gene_slam_min_informative_molecules",
    "gene_slam_plot_top",
    "bin_size",
    "coordinate_um_per_unit",
    "chip_bounds",
    "lr_spatial_min_bin_molecules",
    "lr_spatial_min_gene_bins",
    "lr_spatial_hvg",
    "lr_spatial_pca_components",
    "lr_spatial_neighbors",
    "lr_spatial_regress_log_depth",
    "lr_spatial_graph_weight",
    "lr_spatial_min_cluster_bins",
    "lr_spatial_leiden_resolution",
    "lr_spatial_seed",
    "lr_spatial_marker_top",
    "lr_spatial_marker_display_top",
    "lr_spatial_marker_min_logfoldchange",
    "lr_spatial_marker_max_fdr",
    "lr_spatial_chunk_rows",
    "joint_min_informative_molecules",
    "joint_gene_plot_top",
    "plot_dpi",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True, help="Path to TSV sample sheet")
    parser.add_argument("--config", default="", help="Optional top-level YAML config")
    return parser.parse_args()


def read_top_level_config(path: str) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path:
        return values
    config_path = Path(path)
    if not config_path.is_file():
        raise ValueError(f"Config not found: {config_path}")
    for raw in config_path.read_text(encoding="utf-8").splitlines():
        if not raw or raw[:1].isspace() or raw.lstrip().startswith("#") or ":" not in raw:
            continue
        key, value = raw.split(":", 1)
        values[key.strip()] = value.split("#", 1)[0].strip().strip("'\"")
    return values


def fail(message: str) -> int:
    print(message, file=sys.stderr)
    return 1


def executable_exists(value: str) -> bool:
    path = Path(value)
    return (path.is_file() and path.stat().st_size > 0) if path.is_absolute() else shutil.which(value) is not None


def main() -> int:
    args = parse_args()
    sample_path = Path(args.samples)
    if not sample_path.is_file():
        return fail(f"Sample sheet not found: {sample_path}")
    try:
        config = read_top_level_config(args.config)
    except ValueError as exc:
        return fail(str(exc))
    analysis_type = config.get("analysis_type", "")
    expected_species = config.get("species", "")
    expected_chemistry = config.get("chemistry", "")
    if analysis_type == "lr_celatlas_spatial_slam":
        missing_config = [name for name in LR_CONFIG_REQUIRED if not config.get(name, "").strip()]
        if missing_config:
            return fail("Missing required LR-Celatlas config keys: " + ", ".join(missing_config))
        placeholders = [name for name in LR_CONFIG_REQUIRED if "CHANGE_ME" in config.get(name, "")]
        if placeholders:
            return fail("Replace CHANGE_ME values in config keys: " + ", ".join(placeholders))
        if config.get("run_mode") not in {"full", "pilot"}:
            return fail("run_mode must be full or pilot")
        if config.get("lr_spatial_resolution_grid", "").strip():
            return fail("lr_spatial_resolution_grid is deprecated; enter lr_spatial_leiden_resolution directly")
        for name in ("adapter_python", "analysis_python", "spatial_python", "minimap2_bin", "samtools_bin", "gffcompare_bin"):
            if not executable_exists(config[name]):
                return fail(f"Configured executable is unavailable: {name}={config[name]}")
        if not Path(config["flames_python_dir"]).is_dir():
            return fail(f"Configured FLAMES Python directory is unavailable: {config['flames_python_dir']}")
        positive_ints = (
            "flames_min_support", "min_mapq", "min_tc_baseq", "mismatch_max_reads",
            "min_t_covered", "min_tc_count", "gene_slam_min_informative_molecules",
            "gene_slam_plot_top", "bin_size", "chip_bounds", "lr_spatial_min_bin_molecules",
            "lr_spatial_min_gene_bins", "lr_spatial_hvg", "lr_spatial_pca_components",
            "lr_spatial_neighbors", "lr_spatial_min_cluster_bins", "lr_spatial_marker_top",
            "lr_spatial_marker_display_top", "lr_spatial_chunk_rows",
            "joint_min_informative_molecules", "joint_gene_plot_top", "plot_dpi",
        )
        try:
            invalid_ints = [name for name in positive_ints if int(config[name]) <= 0]
        except ValueError as exc:
            return fail(f"Scientific integer parameter is invalid: {exc}")
        if invalid_ints:
            return fail("Scientific integer parameters must be positive: " + ", ".join(invalid_ints))
        unit_intervals = (
            "flames_min_count_fraction", "flames_min_support_fraction", "flames_downsample_ratio",
            "slam_fdr", "max_site_other_mismatch_rate", "max_read_other_mismatch_rate",
            "lr_spatial_graph_weight", "lr_spatial_marker_max_fdr",
        )
        try:
            invalid_unit = [
                name for name in unit_intervals
                if not 0 <= float(config[name]) <= 1
            ]
            positive_floats = (
                "coordinate_um_per_unit", "lr_spatial_leiden_resolution",
                "lr_spatial_marker_min_logfoldchange",
            )
            invalid_positive = [name for name in positive_floats if float(config[name]) <= 0]
            int(config["lr_spatial_seed"])
            strand = int(config["flames_strand_specific"])
        except ValueError as exc:
            return fail(f"Scientific numeric parameter is invalid: {exc}")
        if invalid_unit:
            return fail("Scientific proportions must be in [0, 1]: " + ", ".join(invalid_unit))
        if invalid_positive:
            return fail("Scientific parameters must be positive: " + ", ".join(invalid_positive))
        if strand not in {0, 1, 2}:
            return fail("flames_strand_specific must be 0, 1, or 2")
        if config["lr_spatial_regress_log_depth"].lower() not in {"true", "false"}:
            return fail("lr_spatial_regress_log_depth must be true or false")
        variant = config.get("external_variant_vcf", "")
        if variant and (not Path(variant).is_file() or Path(variant).stat().st_size == 0):
            return fail(f"Configured external_variant_vcf is missing or empty: {variant}")

    with sample_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames or []
        required = list(LR_REQUIRED if analysis_type == "lr_celatlas_spatial_slam" else GENERIC_REQUIRED)
        missing = [name for name in required if name not in header]
        if missing:
            return fail("Missing required columns: " + ", ".join(missing))
        if analysis_type == "lr_celatlas_spatial_slam":
            missing_metadata = [name for name in LR_OPTIONAL_METADATA if name not in header]
            if missing_metadata:
                print("WARNING: optional experiment metadata columns are absent: " + ", ".join(missing_metadata))

        seen_ids: set[str] = set()
        seen_chips: dict[str, str] = {}
        row_count = 0
        for row_number, row in enumerate(reader, start=2):
            row_count += 1
            clean = {key: (value or "").strip() for key, value in row.items()}
            sample_id = clean.get("sample_id", "")
            if not sample_id:
                return fail(f"Row {row_number}: sample_id is empty")
            if sample_id in seen_ids:
                return fail(f"Row {row_number}: duplicate sample_id {sample_id}")
            seen_ids.add(sample_id)
            if analysis_type != "lr_celatlas_spatial_slam":
                for name in GENERIC_REQUIRED[1:]:
                    if not clean.get(name):
                        return fail(f"Row {row_number}: {name} is empty")
                for name in ("input_1", "input_2"):
                    if not Path(clean[name]).exists():
                        return fail(f"Row {row_number}: input path does not exist: {clean[name]}")
                continue
            for name in LR_REQUIRED:
                if not clean.get(name):
                    return fail(f"Row {row_number}: {name} is empty for LR-Celatlas analysis")
            for name in LR_PATH_COLUMNS:
                path = Path(clean[name])
                if not path.is_file() or path.stat().st_size == 0:
                    return fail(f"Row {row_number}: missing or empty {name}: {path}")
            for name in ("mmi", "junc_bed"):
                value = clean.get(name, "")
                if value and (not Path(value).is_file() or Path(value).stat().st_size == 0):
                    return fail(f"Row {row_number}: missing or empty {name}: {value}")
            if expected_species and clean["species"] != expected_species:
                return fail(
                    f"Row {row_number}: species {clean['species']} does not match config {expected_species}"
                )
            if expected_chemistry and clean["chemistry"] != expected_chemistry:
                return fail(
                    f"Row {row_number}: chemistry {clean['chemistry']} does not match config {expected_chemistry}"
                )
            if analysis_type == "lr_celatlas_spatial_slam":
                missing_values = [name for name in LR_OPTIONAL_METADATA if name in header and not clean.get(name)]
                if missing_values:
                    print(f"WARNING: Row {row_number}: optional experiment metadata is empty: {', '.join(missing_values)}")
                tissue_image = clean.get("tissue_image", "")
                if tissue_image and (not Path(tissue_image).is_file() or Path(tissue_image).stat().st_size == 0):
                    return fail(f"Row {row_number}: tissue_image is missing or empty: {tissue_image}")
            chip_id = clean["chip_id"]
            prior = seen_chips.get(chip_id)
            if prior and prior != sample_id:
                return fail(f"Row {row_number}: chip_id {chip_id} is already assigned to {prior}")
            seen_chips[chip_id] = sample_id

    if row_count == 0:
        return fail("Sample sheet has no data rows")
    suffix = f", analysis_type={analysis_type}" if analysis_type else ""
    print(f"Sample sheet OK: {sample_path} ({row_count} samples{suffix})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
