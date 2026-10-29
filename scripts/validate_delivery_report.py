#!/usr/bin/env python3
"""Validate the self-contained final HTML report and delivery manifest."""

from __future__ import annotations

import argparse
import base64
import csv
import re
from pathlib import Path


REQUIRED_SECTIONS = (
    "1 样本与输入",
    "5 全长转录组",
    "6 SLAM 分析",
    "7 长读长空间转录组",
    "8 长读长空间转录组与 SLAM 联合",
    "9 方法与参数",
)
FORBIDDEN_TEXT = (
    "Analysis scope",
    "Same-chip spatial integration",
    "single_sample_report",
    "短读长结果",
    "分析模式",
    "SLAM 仅作相对比较",
    "样本身份需优先复核",
    "空间聚类已完成",
    "参数评估",
    "long_read_resolution_scan",
    "threshold_sensitivity",
    "全场",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--forbid-text", action="append", default=[])
    parser.add_argument("--require-section", action="append", default=[])
    parser.add_argument("--max-unpaginated-rows", type=int, default=12)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    report = args.report.resolve()
    manifest = args.manifest.resolve()
    errors: list[str] = []
    if not report.is_file() or report.stat().st_size == 0:
        errors.append(f"report missing or empty: {report}")
        text = ""
    else:
        text = report.read_text(encoding="utf-8", errors="replace")

    hrefs = re.findall(r'\bhref="([^"]+)"', text)
    sources = re.findall(r'\bsrc="([^"]+)"', text)
    if any(not href.startswith("#") for href in hrefs):
        errors.append("report contains external or file href values")
    if any(not source.startswith("data:") for source in sources):
        errors.append("report contains non-data src values")
    all_ids = re.findall(r'\bid="([^"]+)"', text)
    ids = set(all_ids)
    if len(ids) != len(all_ids):
        errors.append("report contains duplicate element IDs")
    stale = [href for href in hrefs if href.startswith("#") and href[1:] not in ids]
    if stale:
        errors.append("report contains stale anchors: " + ", ".join(stale[:5]))
    for source in sources:
        if ";base64," not in source:
            errors.append("report contains a non-base64 image")
            continue
        try:
            if not base64.b64decode(source.split(";base64,", 1)[1], validate=True):
                errors.append("report contains an empty embedded image")
        except ValueError:
            errors.append("report contains an invalid embedded image")
    if 'class="sidebar"' not in text or 'class="footer"' not in text:
        errors.append("report navigation or footer is missing")
    for attrs, body in re.findall(r"<table\b([^>]*)>(.*?)</table>", text, flags=re.I | re.S):
        rows = len(re.findall(r"<tr\b", body, flags=re.I)) - 1
        if rows > args.max_unpaginated_rows and "data-page-size=" not in attrs:
            errors.append(f"table has {rows} rows without pagination")
    for section in (*REQUIRED_SECTIONS, *args.require_section):
        if section not in text:
            errors.append(f"required section is missing: {section}")
    forbidden = (*FORBIDDEN_TEXT, *args.forbid_text)
    for term in forbidden:
        if term and term in text:
            errors.append(f"forbidden text found in report: {term}")

    if not manifest.is_file() or manifest.stat().st_size == 0:
        errors.append(f"manifest missing or empty: {manifest}")
    else:
        with manifest.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        matching = []
        for row in rows:
            target = Path(row.get("target", "")).resolve()
            if target == report:
                matching.append(row)
            if row.get("status") in {"copied", "up_to_date"} and (
                not target.is_file() or target.stat().st_size == 0
            ):
                errors.append(f"manifest target missing: {target}")
        if len(matching) != 1:
            errors.append(f"manifest must contain exactly one final report target; found {len(matching)}")
        manifest_text = manifest.read_text(encoding="utf-8", errors="replace")
        for term in forbidden:
            if term and term in manifest_text:
                errors.append(f"forbidden text found in manifest: {term}")
    html_files = [path.resolve() for path in report.parent.glob("*.html") if path.is_file()]
    if html_files != [report]:
        errors.append(
            "report directory must contain exactly one final HTML report: "
            + ", ".join(str(path) for path in html_files)
        )

    print(f"report: {report}")
    print(f"ok: {not errors}")
    for error in errors:
        print(f"ERROR: {error}")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
