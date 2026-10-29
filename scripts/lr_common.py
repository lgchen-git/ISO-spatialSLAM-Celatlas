#!/usr/bin/env python3
from __future__ import annotations

import csv
import os
import subprocess
import sys
from pathlib import Path
from shutil import which


def die(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def ensure_dir(path: str | Path) -> None:
    Path(path).mkdir(parents=True, exist_ok=True)


def require_file(path: str | Path) -> Path:
    resolved = Path(path)
    if not resolved.exists():
        die(f"file not found: {resolved}")
    return resolved


def resolve_path(root_dir: Path, raw_path: str) -> Path:
    path = Path(raw_path)
    return path if path.is_absolute() else root_dir / path


def command_path(name: str) -> str | None:
    return which(name)


def resolve_python_command(env_name: str = "longread_scrna_10x_qitan") -> list[str]:
    if f"envs/{env_name}/bin/python" in sys.executable:
        return [sys.executable]
    home = Path.home()
    candidates = [
        home / "miniforge3" / "envs" / env_name / "bin" / "python",
        home / "miniconda3" / "envs" / env_name / "bin" / "python",
        home / "anaconda3" / "envs" / env_name / "bin" / "python",
        Path("/media/epibiotek/epi3/miniforge3/envs") / env_name / "bin" / "python",
        Path("/media/epibiotek/epi3/miniconda3/envs") / env_name / "bin" / "python",
    ]
    for candidate in candidates:
        if candidate.exists():
            return [str(candidate)]
    conda = which("conda")
    if conda:
        return [conda, "run", "-n", env_name, "python"]
    return [sys.executable]


def run_command(command: list[str], log_file: Path | None = None) -> None:
    if command:
        py = resolve_python_command()
        if command[0] in {sys.executable, "python", "python3"} and command[0] != py[0]:
            command = py + command[1:]
    if log_file is not None:
        ensure_dir(log_file.parent)
        with log_file.open("a", encoding="utf-8") as handle:
            handle.write(f"[RUN] {' '.join(command)}\n")
            handle.flush()
            proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            assert proc.stdout is not None
            for line in proc.stdout:
                sys.stdout.write(line)
                handle.write(line)
            code = proc.wait()
            if code != 0:
                raise subprocess.CalledProcessError(code, command)
        return
    subprocess.run(command, check=True)


def load_sample_row(samples_tsv: Path, sample_id: str) -> dict[str, str]:
    with samples_tsv.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["sample_id"] == sample_id:
                return row
    die(f"sample not found in sample sheet: {sample_id}")


def load_env_file(path: Path) -> dict[str, str]:
    config: dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            key, value = line.split("=", 1)
            config[key] = value
    return config
