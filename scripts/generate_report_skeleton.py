#!/usr/bin/env python3
"""Generate a formal report skeleton for pipeline runs."""

from __future__ import annotations

import argparse
from pathlib import Path
import shutil


PALETTE_MAP = {
    "okabe_ito": "#0072B2, #D55E00, #009E73, #CC79A7, #E69F00, #56B4E9",
    "tableau_10": "#4E79A7, #F28E2B, #E15759, #76B7B2, #59A14F, #EDC948",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate a report skeleton.")
    parser.add_argument("--project-name", required=True)
    parser.add_argument("--analysis-type", required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--report-title", required=True)
    parser.add_argument("--report-style", required=True)
    parser.add_argument("--report-logo", required=True)
    parser.add_argument("--plot-palette", required=True)
    parser.add_argument("--plot-font-family", required=True)
    parser.add_argument("--plot-base-size", required=True)
    parser.add_argument("--plot-dpi", required=True)
    parser.add_argument("--sample-count", required=True)
    parser.add_argument("--config-path", required=True)
    parser.add_argument("--samples-path", required=True)
    parser.add_argument("--summary-path", required=True)
    parser.add_argument("--figures-dir", required=True)
    parser.add_argument("--tables-dir", required=True)
    parser.add_argument("--reports-dir", required=True)
    parser.add_argument("--deliverables-dir", required=True)
    parser.add_argument("--editable-dir", required=True)
    return parser.parse_args()


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def build_report(args: argparse.Namespace, logo_ref: str, palette_text: str) -> str:
    return f"""# {args.report_title}

![Logo]({logo_ref})

## Project Summary

| Field | Value |
| --- | --- |
| Project | {args.project_name} |
| Analysis Type | {args.analysis_type} |
| Run ID | {args.run_id} |
| Sample Count | {args.sample_count} |
| Config | `{args.config_path}` |
| Samples | `{args.samples_path}` |
| Style | {args.report_style} |

## Methods

| Item | Description |
| --- | --- |
| Workflow Entry | `workflow/run.sh` |
| Parameter Source | Top-level keys in the project config |
| Resource Control | Use `threads`, `memory_gb`, and `max_parallel_samples` from config |
| Resume Strategy | Reuse run directories and checkpoint state |
| QC Summary Table | `{args.summary_path}` |

Replace this section with the actual methods used for the project.

## Results

### Key Metrics

Fill this section from tables in `{args.tables_dir}`.

### Main Figures

Add figure panels from `{args.figures_dir}`.

## Conclusions

State the main findings in short, reviewable sentences.

## QC And Delivery Notes

| Item | Path | Note |
| --- | --- | --- |
| Summary Table | `{args.summary_path}` | Start here for key metrics |
| Figures | `{args.figures_dir}` | Use publication-ready exports |
| Editable Sources | `{args.editable_dir}` | Keep plot and report sources here |
| Deliverables | `{args.deliverables_dir}` | Final client-facing files |

## Plot Style

| Item | Value |
| --- | --- |
| Palette | {args.plot_palette} |
| Palette Colors | {palette_text} |
| Font Family | {args.plot_font_family} |
| Base Font Size | {args.plot_base_size} |
| DPI | {args.plot_dpi} |

## Attention Points

- Replace placeholders before delivery.
- Prefer tables for key statistics and sample-level review.
- Keep methods, results, and conclusions specific.
"""


def main() -> int:
    args = parse_args()

    reports_dir = Path(args.reports_dir)
    deliverables_dir = Path(args.deliverables_dir)
    editable_dir = Path(args.editable_dir)
    logo_path = Path(args.report_logo)

    if not logo_path.is_file():
        raise SystemExit(f"Report logo not found: {logo_path}")

    report_name = f"{args.run_id}_report.md"
    editable_path = editable_dir / report_name
    deliverable_path = deliverables_dir / report_name
    reports_path = reports_dir / report_name
    style_path = editable_dir / "report_style.tsv"
    summary_stub_path = editable_dir / "report_sections.tsv"

    palette_text = PALETTE_MAP.get(args.plot_palette, args.plot_palette)
    logo_copy_path = editable_dir / logo_path.name
    shutil.copy2(logo_path, logo_copy_path)

    logo_ref = logo_copy_path.name
    report_text = build_report(args, logo_ref, palette_text)

    write_text(editable_path, report_text)
    write_text(reports_path, report_text)
    write_text(deliverable_path, report_text)
    write_text(
        style_path,
        "field\tvalue\n"
        f"palette\t{args.plot_palette}\n"
        f"palette_colors\t{palette_text}\n"
        f"font_family\t{args.plot_font_family}\n"
        f"base_size\t{args.plot_base_size}\n"
        f"dpi\t{args.plot_dpi}\n",
    )
    write_text(
        summary_stub_path,
        "section\tstatus\tnote\n"
        "project_summary\tdraft\tupdate metadata if needed\n"
        "methods\tdraft\treplace with project methods\n"
        "results\tdraft\tfill from tables and figures\n"
        "conclusions\tdraft\twrite concise takeaways\n"
        "attention_points\tdraft\tlist delivery caveats\n",
    )

    print(editable_path)
    print(deliverable_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
