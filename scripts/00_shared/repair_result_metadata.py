#!/usr/bin/env python3
"""Remove retired scan outputs and repair internal paths in an eight-stage result."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


STALE_SPATIAL = (
    "06_whole_transcriptome/01_tables/long_read_resolution_scan.tsv",
    "06_whole_transcriptome/02_figures/long_read_resolution_scan.png",
    "06_whole_transcriptome/02_figures/long_read_resolution_scan.pdf",
    "06_whole_transcriptome/02_figures/new_supported_fraction.png",
    "06_whole_transcriptome/02_figures/new_supported_fraction.pdf",
)


def atomic_write(path: Path, text: str) -> None:
    temporary = path.with_name(path.name + ".repair")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def strip_retired_json_fields(path: Path, payload: object, sample_id: str) -> object:
    if not isinstance(payload, dict):
        return payload
    if path.name == "long_read_spatial_clustering.json":
        payload.pop("resolution_grid", None)
        parameters = payload.get("parameters")
        if isinstance(parameters, dict):
            parameters.pop("resolution_grid", None)
            parameters.pop("resolution_selection", None)
    if path.name == f"{sample_id}.slam_summary.json":
        for key in ("interpretation", "label_duration_known", "matched_no4su_control"):
            payload.pop(key, None)
    if path.name == "stage_complete.json":
        outputs = payload.get("outputs")
        if isinstance(outputs, list):
            payload["outputs"] = [
                value
                for value in outputs
                if isinstance(value, str)
                and "threshold_sensitivity" not in value
                and "exploratory_summary" not in value
                and "long_read_resolution_scan" not in value
                and not value.endswith("/new_supported_fraction.png")
                and not value.endswith("/new_supported_fraction.pdf")
            ]
    return payload


def repair_sample(
    sample_root: Path,
    sample_id: str,
    published_root: Path | None = None,
) -> dict[str, int]:
    sample_root = sample_root.resolve()
    published_root = (published_root or sample_root).resolve()
    retired_relative = (
        f"05_slam/01_tables/{sample_id}.threshold_sensitivity.tsv",
        f"05_slam/01_tables/{sample_id}.exploratory_summary.json",
        *STALE_SPATIAL,
    )
    retired_published = {str(published_root / relative) for relative in retired_relative}
    removed = 0
    for relative in retired_relative:
        target = sample_root / relative
        if target.is_file():
            target.unlink()
            removed += 1

    old_exploratory_calls = (
        published_root
        / "06_slam_exploratory"
        / "01_tables"
        / f"{sample_id}.exploratory_molecule_calls.tsv.gz"
    )
    current_calls = (
        published_root / "05_slam" / "01_tables" / f"{sample_id}.slam_molecule_calls.tsv.gz"
    )
    replacements = {
        str(old_exploratory_calls): str(current_calls),
        str(published_root / "06_slam_exploratory"): str(published_root / "05_slam"),
        str(published_root / "04_slam"): str(published_root / "05_slam"),
        str(published_root / "05_spatial_clustering"): str(
            published_root / "06_whole_transcriptome"
        ),
        str(published_root / "06_report"): str(published_root / "08_report"),
    }

    rewritten = 0
    metadata_files = sorted(sample_root.glob("*/05_metadata/*"))
    metadata_files.extend(sorted(sample_root.glob("*/04_tool_output/**/*.json")))
    for path in metadata_files:
        if not path.is_file() or path.suffix not in {".json", ".tsv"}:
            continue
        original = path.read_text(encoding="utf-8")
        updated = original
        for old, new in replacements.items():
            updated = updated.replace(old, new)
        if path.name in {"output_manifest.tsv", "large_file_manifest.tsv"}:
            lines = updated.splitlines()
            kept = [
                line
                for line in lines
                if not any(retired in line for retired in retired_published)
                and "long_read_resolution_scan" not in line
                and not line.split("\t", 1)[0].endswith("/new_supported_fraction.png")
                and not line.split("\t", 1)[0].endswith("/new_supported_fraction.pdf")
            ]
            updated = "\n".join(kept) + ("\n" if original.endswith("\n") else "")
        if path.suffix == ".json":
            payload = json.loads(updated)
            payload = strip_retired_json_fields(path, payload, sample_id)
            updated = json.dumps(payload, ensure_ascii=False, indent=2) + "\n"
        if updated != original:
            atomic_write(path, updated)
            rewritten += 1
    return {"removed_files": removed, "rewritten_metadata": rewritten}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-root", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--published-root", type=Path)
    args = parser.parse_args()
    summary = repair_sample(args.sample_root, args.sample_id, args.published_root)
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
