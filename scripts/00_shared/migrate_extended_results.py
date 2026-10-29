#!/usr/bin/env python3
"""Migrate a six-stage sample tree to the reusable eight-stage layout."""

from __future__ import annotations

import argparse
import json
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path

from repair_result_metadata import repair_sample


STAGES = (
    "01_read_qc",
    "02_barcode_cdna",
    "03_align",
    "04_full_length_transcriptome",
    "05_slam",
    "06_whole_transcriptome",
    "07_slam_joint",
    "08_report",
)
SUBDIRS = ("01_tables", "02_figures", "03_logs", "04_tool_output", "05_metadata")
STALE_SPATIAL = {
    "01_tables/long_read_resolution_scan.tsv",
    "02_figures/long_read_resolution_scan.png",
    "02_figures/long_read_resolution_scan.pdf",
    "02_figures/new_supported_fraction.png",
    "02_figures/new_supported_fraction.pdf",
}
REBUILT_METADATA = {
    "05_metadata/stage_complete.json",
    "05_metadata/output_manifest.tsv",
    "05_metadata/large_file_manifest.tsv",
    "05_metadata/migration_provenance.tsv",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-root", required=True, type=Path)
    parser.add_argument("--tmp-root", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--transcriptome-source", required=True, type=Path)
    parser.add_argument("--archive-label", default="pre_extended_layout")
    parser.add_argument("--apply", action="store_true")
    return parser.parse_args()


def link_tree(source: Path, target: Path, excluded: set[str] | None = None) -> None:
    excluded = excluded or set()
    if not source.is_dir():
        raise RuntimeError(f"Source directory not found: {source}")
    for path in sorted(source.rglob("*")):
        if not path.is_file():
            continue
        relative = path.relative_to(source)
        if relative.as_posix() in excluded:
            continue
        destination = target / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists():
            raise FileExistsError(destination)
        os.link(path, destination)


def rewrite_json(path: Path, transform) -> None:
    if not path.is_file():
        return
    payload = transform(json.loads(path.read_text(encoding="utf-8")))
    temporary = path.with_name(path.name + ".rewrite")
    temporary.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def strip_parameter_scan(payload: dict) -> dict:
    payload.pop("resolution_grid", None)
    parameters = payload.get("parameters")
    if isinstance(parameters, dict):
        parameters.pop("resolution_grid", None)
        parameters.pop("resolution_selection", None)
    return payload


def write_stage_metadata(stage: Path, source_label: str, published_stage: Path) -> None:
    metadata = stage / "05_metadata"
    now = datetime.now(timezone.utc).isoformat()
    files = [
        path
        for path in sorted(stage.rglob("*"))
        if path.is_file()
        and path.parent != metadata
        and path.name not in {"stage_complete.json", "output_manifest.tsv", "large_file_manifest.tsv"}
    ]
    (metadata / "migration_provenance.tsv").write_text(
        "key\tvalue\n"
        "migration_schema\t2\n"
        f"migrated_at_utc\t{now}\n"
        f"source_stage\t{source_label}\n"
        "storage_method\thard_link\n",
        encoding="utf-8",
    )
    with (metadata / "output_manifest.tsv").open("w", encoding="utf-8") as handle:
        handle.write("path\tsize_bytes\tmtime_ns\n")
        for path in files:
            stat = path.stat()
            published = (published_stage / path.relative_to(stage)).resolve()
            handle.write(f"{published}\t{stat.st_size}\t{stat.st_mtime_ns}\n")
    with (metadata / "large_file_manifest.tsv").open("w", encoding="utf-8") as handle:
        handle.write("path\tsize_bytes\tpurpose\n")
        for path in files:
            if path.stat().st_size >= 100 * 1024 * 1024:
                published = (published_stage / path.relative_to(stage)).resolve()
                handle.write(f"{published}\t{path.stat().st_size}\tretained_stage_output\n")
    (metadata / "stage_complete.json").write_text(
        json.dumps(
            {
                "stage": stage.name,
                "completed_at_utc": now,
                "completion_mode": "migrated",
                "source_stage": source_label,
                "outputs": [str((published_stage / path.relative_to(stage)).resolve()) for path in files],
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def main() -> int:
    args = parse_args()
    results_root = args.results_root.resolve()
    tmp_root = args.tmp_root.resolve()
    source = results_root / args.sample_id
    staging = results_root / f".{args.sample_id}.extended_migration"
    archive = tmp_root / "results_legacy" / args.archive_label / args.sample_id
    transcriptome = args.transcriptome_source.resolve()
    print(f"source\t{source}")
    print(f"transcriptome_source\t{transcriptome}")
    print(f"archive\t{archive}")
    print(f"replacement\t{source}")
    if not args.apply:
        print("mode\tdry-run")
        return 0
    required = ("01_read_qc", "02_barcode_cdna", "03_align", "04_slam", "05_spatial_clustering", "06_report")
    missing = [name for name in required if not (source / name).is_dir()]
    if missing:
        raise SystemExit("Missing six-stage inputs: " + ", ".join(missing))
    if not transcriptome.is_dir():
        raise SystemExit(f"Transcriptome source not found: {transcriptome}")
    if staging.exists() or archive.exists():
        raise SystemExit(f"Migration target already exists: {staging} / {archive}")
    mapping = {
        "01_read_qc": (source / "01_read_qc", "01_read_qc"),
        "02_barcode_cdna": (source / "02_barcode_cdna", "02_barcode_cdna"),
        "03_align": (source / "03_align", "03_align"),
        "04_full_length_transcriptome": (transcriptome, str(transcriptome)),
        "05_slam": (source / "04_slam", "04_slam"),
        "06_whole_transcriptome": (source / "05_spatial_clustering", "05_spatial_clustering"),
    }
    try:
        for stage_name in STAGES:
            for subdir in SUBDIRS:
                (staging / stage_name / subdir).mkdir(parents=True, exist_ok=True)
        for stage_name, (origin, source_label) in mapping.items():
            excluded = set(REBUILT_METADATA)
            if stage_name == "06_whole_transcriptome":
                excluded.update(STALE_SPATIAL)
            link_tree(origin, staging / stage_name, excluded)
            if stage_name == "06_whole_transcriptome":
                rewrite_json(
                    staging / stage_name / "05_metadata/long_read_spatial_clustering.json",
                    strip_parameter_scan,
                )
            write_stage_metadata(staging / stage_name, source_label, source / stage_name)
        for stage_name in ("07_slam_joint", "08_report"):
            write_stage_metadata(staging / stage_name, "new_stage", source / stage_name)
        repair_sample(staging, args.sample_id, source)
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
