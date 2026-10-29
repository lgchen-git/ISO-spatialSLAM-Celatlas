#!/usr/bin/env python3
"""Archive selected result stages with an auditable manifest."""

from __future__ import annotations

import argparse
import csv
import shutil
from datetime import datetime, timezone
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-dir", action="append", required=True, help="LABEL=absolute sample result directory")
    parser.add_argument("--stage", action="append", required=True)
    parser.add_argument("--archive-root", required=True)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    archive_root = Path(args.archive_root).resolve()
    if not args.apply:
        raise SystemExit("Refusing to move results without --apply")
    if archive_root.exists():
        raise SystemExit(f"Archive target already exists: {archive_root}")
    stages = list(dict.fromkeys(args.stage))
    sources: list[tuple[str, Path, Path]] = []
    for item in args.sample_dir:
        if "=" not in item:
            raise SystemExit(f"Expected LABEL=PATH: {item}")
        label, value = item.split("=", 1)
        sample_dir = Path(value).resolve()
        for stage in stages:
            source = sample_dir / stage
            if not source.is_dir():
                raise SystemExit(f"Missing source stage: {source}")
            sources.append((label, source, archive_root / label / sample_dir.name / stage))
    archive_root.mkdir(parents=True)
    rows = []
    for label, source, target in sources:
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(source), str(target))
        rows.append((label, source, target, datetime.now(timezone.utc).isoformat()))
    with (archive_root / "archive_manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["profile", "source", "archive", "archived_at_utc"])
        writer.writerows(rows)


if __name__ == "__main__":
    main()
