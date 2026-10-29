#!/usr/bin/env python3
"""Capture stage-relevant software versions without importing heavy packages."""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path


def read_config(path: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        if not raw or raw[:1].isspace() or raw.lstrip().startswith("#") or ":" not in raw:
            continue
        key, value = raw.split(":", 1)
        values[key.strip()] = value.split("#", 1)[0].strip().strip("'\"")
    return values


def package_command(python: str, packages: list[str]) -> list[str]:
    code = (
        "import importlib.metadata as m; "
        f"print('; '.join(p + '=' + m.version(p) for p in {packages!r}))"
    )
    return [python, "-c", code]


def first_line(command: list[str]) -> tuple[str, str]:
    try:
        result = subprocess.run(command, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=30)
        line = next((item.strip() for item in result.stdout.splitlines() if item.strip()), "unknown")
        return line, "ok"
    except (OSError, subprocess.SubprocessError) as exc:
        return str(exc).replace("\t", " ").replace("\n", " "), "unavailable"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--stage", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    config = read_config(args.config)
    analysis_python = config.get("analysis_python", "python3")
    spatial_python = config.get("spatial_python", analysis_python)
    commands: list[tuple[str, list[str]]] = [("python", [analysis_python, "--version"])]
    if args.stage == "01_read_qc":
        commands.append(("NanoPlot", [config.get("nanoplot_bin", "NanoPlot"), "--version"]))
    elif args.stage == "03_align":
        commands.extend([
            ("minimap2", [config.get("minimap2_bin", "minimap2"), "--version"]),
            ("samtools", [config.get("samtools_bin", "samtools"), "--version"]),
        ])
    elif args.stage == "04_full_length_transcriptome":
        commands.extend([
            (
                "FLAMES",
                [
                    config.get("rscript_bin", "Rscript"),
                    "-e",
                    'cat(as.character(packageVersion("FLAMES")), "\\n")',
                ],
            ),
            ("gffcompare", [config.get("gffcompare_bin", "gffcompare"), "--version"]),
        ])
    elif args.stage == "05_slam":
        commands.append(("python_packages", package_command(analysis_python, ["pysam", "pandas", "numpy", "scipy", "matplotlib"])))
    elif args.stage == "06_whole_transcriptome":
        commands.append(("python_packages", package_command(spatial_python, ["scanpy", "anndata", "leidenalg", "pandas", "numpy", "scipy", "matplotlib"])))
    elif args.stage == "07_slam_joint":
        commands.append(("python_packages", package_command(analysis_python, ["pandas", "numpy", "scipy", "matplotlib"])))

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["tool", "version", "status", "command"])
        for tool, command in commands:
            version, status = first_line(command)
            writer.writerow([tool, version, status, " ".join(command)])


if __name__ == "__main__":
    main()
