#!/usr/bin/env python3
"""Migrate a legacy sample result tree to the reusable six-stage layout."""

from __future__ import annotations

import argparse
import json
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path


STAGES = (
    "01_read_qc",
    "02_barcode_cdna",
    "03_align",
    "04_slam",
    "05_spatial_clustering",
    "06_report",
)
SUBDIRS = ("01_tables", "02_figures", "03_logs", "04_tool_output", "05_metadata")
ALIGN_FIGURES = {"alignment_summary.png", "top50_barcodes.png"}
READ_FIGURES = {"raw_read_length.png", "raw_mean_q.png"}
BARCODE_FIGURES = {
    "structure_status.png",
    "cdna_length.png",
    "barcode_status.png",
    "barcode_edit_distance.png",
    "confidence_class.png",
}
SOURCE_STAGES = {
    "01_read_qc": ["01_input_qc", "02_read_qc", "04_align/02_figures"],
    "02_barcode_cdna": ["03_barcode_cdna", "04_align/02_figures"],
    "03_align": ["04_align"],
    "04_slam": ["06_slam_exploratory"],
    "05_spatial_clustering": ["07_spatial_integration"],
    "06_report": [],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-root", required=True, type=Path)
    parser.add_argument("--tmp-root", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--archive-label", default="pre_reusable_layout_20260720")
    parser.add_argument("--apply", action="store_true", help="Perform the migration; default is dry-run.")
    parser.add_argument(
        "--refresh-metadata",
        action="store_true",
        help="Refresh manifests for an already migrated sample without moving result files.",
    )
    return parser.parse_args()


def link_file(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        raise FileExistsError(f"Target collision: {target}")
    os.link(source, target)


def link_tree(source: Path, target: Path, exclude_figures: bool = False) -> None:
    if not source.is_dir():
        return
    for path in sorted(source.rglob("*")):
        if not path.is_file():
            continue
        relative = path.relative_to(source)
        if exclude_figures and relative.parts[:1] == ("02_figures",):
            continue
        link_file(path, target / relative)


def link_named_figures(source: Path, target: Path, names: set[str]) -> None:
    figure_root = source / "02_figures"
    for name in sorted(names):
        path = figure_root / name
        if path.is_file():
            link_file(path, target / "02_figures" / name)


def rename_if_present(root: Path, old: str, new: str) -> None:
    source = root / old
    if source.is_file():
        target = root / new
        target.parent.mkdir(parents=True, exist_ok=True)
        source.rename(target)


def write_stage_metadata(stage: Path, source_stages: list[str], published_stage: Path) -> None:
    meta = stage / "05_metadata"
    now = datetime.now(timezone.utc).isoformat()
    files = [
        path for path in sorted(stage.rglob("*"))
        if path.is_file() and path.parent != meta
    ]
    (meta / "migration_provenance.tsv").write_text(
        "key\tvalue\n"
        f"migration_schema\t1\n"
        f"migrated_at_utc\t{now}\n"
        f"source_stages\t{','.join(source_stages)}\n"
        "storage_method\thard_link\n",
        encoding="utf-8",
    )
    with (meta / "output_manifest.tsv").open("w", encoding="utf-8") as handle:
        handle.write("path\tsize_bytes\tmtime_ns\n")
        for path in files:
            stat = path.stat()
            published = (published_stage / path.relative_to(stage)).resolve()
            handle.write(f"{published}\t{stat.st_size}\t{stat.st_mtime_ns}\n")
    (meta / "large_file_manifest.tsv").write_text(
        "path\tsize_bytes\tpurpose\n"
        + "".join(
            f"{(published_stage / path.relative_to(stage)).resolve()}\t{path.stat().st_size}\tretained_stage_output\n"
            for path in files
            if path.stat().st_size >= 100 * 1024 * 1024
        ),
        encoding="utf-8",
    )
    (meta / "stage_complete.json").write_text(
        json.dumps(
            {
                "stage": stage.name,
                "completed_at_utc": now,
                "completion_mode": "migrated",
                "source_stages": source_stages,
                "outputs": [
                    str((published_stage / path.relative_to(stage)).resolve())
                    for path in files
                ],
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def build_staging(source: Path, staging: Path, sample_id: str) -> None:
    for stage_name in STAGES:
        for subdir in SUBDIRS:
            (staging / stage_name / subdir).mkdir(parents=True, exist_ok=True)

    read_stage = staging / "01_read_qc"
    barcode_stage = staging / "02_barcode_cdna"
    align_stage = staging / "03_align"
    slam_stage = staging / "04_slam"
    spatial_stage = staging / "05_spatial_clustering"

    link_tree(source / "02_read_qc", read_stage)
    old_input_meta = source / "01_input_qc" / "05_metadata"
    if old_input_meta.is_dir():
        for path in sorted(old_input_meta.iterdir()):
            if path.is_file() and path.name not in {"stage_complete.json", "output_manifest.tsv"}:
                link_file(path, read_stage / "05_metadata" / f"input_{path.name}")
    input_manifest = source / "01_input_qc" / "01_tables" / "input_manifest.tsv"
    if input_manifest.is_file():
        link_file(input_manifest, read_stage / "05_metadata" / "input_manifest.tsv")

    link_tree(source / "03_barcode_cdna", barcode_stage)
    link_tree(source / "04_align", align_stage, exclude_figures=True)
    link_named_figures(source / "04_align", read_stage, READ_FIGURES)
    link_named_figures(source / "04_align", barcode_stage, BARCODE_FIGURES)
    link_named_figures(source / "04_align", align_stage, ALIGN_FIGURES)
    rename_if_present(align_stage, "05_metadata/pilot_gate.json", "05_metadata/quality_gate.json")

    link_tree(source / "06_slam_exploratory", slam_stage)
    rename_if_present(
        slam_stage,
        f"01_tables/{sample_id}.exploratory_molecule_calls.tsv.gz",
        f"01_tables/{sample_id}.slam_molecule_calls.tsv.gz",
    )
    rename_if_present(
        slam_stage,
        f"05_metadata/{sample_id}.exploratory_summary.json",
        f"05_metadata/{sample_id}.slam_summary.json",
    )
    link_tree(source / "07_spatial_integration", spatial_stage)

    for stage_name in STAGES:
        write_stage_metadata(
            staging / stage_name,
            SOURCE_STAGES[stage_name],
            source / stage_name,
        )


def main() -> int:
    args = parse_args()
    results_root = args.results_root.resolve()
    tmp_root = args.tmp_root.resolve()
    source = results_root / args.sample_id
    staging = results_root / f".{args.sample_id}.reusable_migration"
    archive = tmp_root / "results_legacy" / args.archive_label / args.sample_id
    print(f"source\t{source}")
    print(f"archive\t{archive}")
    print(f"replacement\t{source}")
    print("storage\thard links; unchanged file contents are not duplicated")
    if args.apply and args.refresh_metadata:
        raise SystemExit("ERROR: choose either --apply or --refresh-metadata")
    if args.refresh_metadata:
        missing = [stage for stage in STAGES[:5] if not (source / stage).is_dir()]
        if missing:
            raise SystemExit("ERROR: migrated stages missing: " + ", ".join(missing))
        for stage_name in STAGES[:5]:
            write_stage_metadata(
                source / stage_name,
                SOURCE_STAGES[stage_name],
                source / stage_name,
            )
        print("mode\tmetadata-refreshed")
        return 0
    if not args.apply:
        print("mode\tdry-run")
        return 0
    if not source.is_dir():
        raise SystemExit(f"ERROR: source sample directory not found: {source}")
    if archive.exists() or staging.exists():
        raise SystemExit(f"ERROR: archive or staging path already exists: {archive} / {staging}")
    required = ("01_input_qc", "02_read_qc", "03_barcode_cdna", "04_align", "06_slam_exploratory", "07_spatial_integration")
    missing = [name for name in required if not (source / name).is_dir()]
    if missing:
        raise SystemExit("ERROR: legacy stages missing: " + ", ".join(missing))
    try:
        build_staging(source, staging, args.sample_id)
        archive.parent.mkdir(parents=True, exist_ok=True)
        source.rename(archive)
        staging.rename(source)
    except Exception:
        if staging.exists():
            shutil.rmtree(staging)
        raise
    print("mode\tapplied")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
