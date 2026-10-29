#!/usr/bin/env python3
"""Run the staged LR-Celatlas spatial SLAM workflow."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


STAGES = [
    "01_read_qc", "02_barcode_cdna", "03_align",
    "04_full_length_transcriptome", "05_slam", "06_whole_transcriptome",
    "07_slam_joint", "08_report",
]


def args_parser() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--logs-dir", required=True)
    parser.add_argument("--tmp-dir", required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--from", dest="from_stage", choices=STAGES, default=STAGES[0])
    parser.add_argument("--to", dest="to_stage", choices=STAGES, default=STAGES[-1])
    parser.add_argument("--force-stage", action="append", choices=STAGES, default=[])
    parser.add_argument("--report-title", default="", help="Override report_title from config.")
    parser.add_argument("--refresh-visuals", action="store_true")
    parser.add_argument("--promo-figures", action="store_true", help="Build the optional PPT figure pack from completed outputs.")
    return parser.parse_args()


def read_config(path: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        if not raw or raw[:1].isspace() or raw.lstrip().startswith("#") or ":" not in raw:
            continue
        key, value = raw.split(":", 1)
        values[key.strip()] = value.split("#", 1)[0].strip().strip("'\"")
    return values


def read_samples(path: str) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return [{key: (value or "").strip() for key, value in row.items()} for row in csv.DictReader(handle, delimiter="\t")]


def truth(value: str, default: bool = False) -> bool:
    if not value:
        return default
    return value.lower() in {"1", "true", "yes", "on"}


def ready(paths: list[Path]) -> bool:
    return all(path.is_file() and path.stat().st_size > 0 for path in paths)


class Runner:
    def __init__(self, root: Path, config_path: Path, config: dict[str, str], logs: Path, tmp: Path, force: set[str]):
        self.root = root
        self.config_path = config_path
        self.config = config
        self.logs = logs
        self.tmp = tmp
        self.force = force
        self.config_sha256 = hashlib.sha256(config_path.read_bytes()).hexdigest()
        self.current_stage = ""
        self.current_stage_root: Path | None = None
        self.tmp.mkdir(parents=True, exist_ok=True)
        os.environ.update({
            "TMPDIR": str(self.tmp), "TMP": str(self.tmp), "TEMP": str(self.tmp),
            "XDG_CACHE_HOME": str(self.tmp / "xdg-cache"), "MPLCONFIGDIR": str(self.tmp / "mplconfig"),
            "NUMBA_CACHE_DIR": str(self.tmp / "numba-cache"),
            "MPLBACKEND": "Agg", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
        })
        bins = {str(Path(value).parent) for key, value in config.items() if key.endswith("_bin") and value.startswith("/")}
        os.environ["PATH"] = os.pathsep.join(sorted(bins)) + os.pathsep + os.environ.get("PATH", "")

    def run(self, command: list[str], label: str) -> None:
        if self.current_stage_root is not None:
            log = self.current_stage_root / "03_logs" / f"{label}.log"
        else:
            log = self.logs / self.current_stage / f"{label}.log"
        log.parent.mkdir(parents=True, exist_ok=True)
        with log.open("a", encoding="utf-8") as handle:
            handle.write("[RUN] " + " ".join(command) + "\n")
            handle.flush()
            proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=os.environ.copy())
            assert proc.stdout is not None
            for line in proc.stdout:
                sys.stdout.write(line)
                handle.write(line)
            code = proc.wait()
        if code:
            raise subprocess.CalledProcessError(code, command)

    def stage(self, stage: str, outputs: list[Path], action) -> None:
        self.current_stage = stage
        stage_root = next((parent for parent in outputs[0].parents if parent.name == stage), None)
        if stage_root is None:
            raise RuntimeError(f"Cannot locate stage root for {outputs[0]}")
        marker = stage_root / "05_metadata" / "stage_complete.json"
        if stage not in self.force and marker.is_file() and ready(outputs):
            try:
                prior = json_load(marker)
            except (OSError, json.JSONDecodeError):
                prior = {}
            if prior.get("config_sha256") == self.config_sha256:
                print(f"[SKIP] {stage}: complete", flush=True)
                return
            print(f"[RERUN] {stage}: configuration changed", flush=True)
        print(f"[STAGE] {stage}", flush=True)
        self.current_stage_root = stage_root
        metadata = stage_root / "05_metadata"
        metadata.mkdir(parents=True, exist_ok=True)
        self.run(
            [
                sys.executable,
                str(self.root / "scripts" / "capture_git_provenance.py"),
                "--repo-root",
                str(self.root),
                "--output",
                str(metadata / "git_provenance.tsv"),
                "--dirty-policy",
                self.config.get("git_dirty_policy", "warn"),
            ],
            "git_provenance",
        )
        self.run(
            [
                sys.executable,
                str(self.root / "scripts" / "00_shared" / "capture_software_versions.py"),
                "--config",
                str(self.config_path),
                "--stage",
                stage,
                "--output",
                str(metadata / "software_versions.tsv"),
            ],
            "software_versions",
        )
        action()
        if not ready(outputs):
            raise RuntimeError(f"{stage} did not create required outputs: {outputs}")
        with (metadata / "output_manifest.tsv").open("w", encoding="utf-8") as handle:
            handle.write("path\tsize_bytes\tmtime_ns\n")
            for path in outputs:
                stat = path.stat()
                handle.write(f"{path.resolve()}\t{stat.st_size}\t{stat.st_mtime_ns}\n")
        threshold = int(float(self.config.get("large_file_threshold_mb", "100")) * 1024 * 1024)
        with (metadata / "large_file_manifest.tsv").open("w", encoding="utf-8") as handle:
            handle.write("path\tsize_bytes\tpurpose\n")
            roots = [(stage_root, "retained_stage_output"), (self.tmp / stage, "non_deliverable_intermediate")]
            seen: set[Path] = set()
            for scan_root, purpose in roots:
                if not scan_root.is_dir():
                    continue
                for path in sorted(scan_root.rglob("*")):
                    if not path.is_file() or path.stat().st_size < threshold:
                        continue
                    resolved = path.resolve()
                    if resolved in seen:
                        continue
                    seen.add(resolved)
                    handle.write(f"{resolved}\t{path.stat().st_size}\t{purpose}\n")
        marker.write_text(
            json.dumps(
                {
                    "stage": stage,
                    "config_sha256": self.config_sha256,
                    "completed_at_utc": datetime.now(timezone.utc).isoformat(),
                    "outputs": [str(path.resolve()) for path in outputs],
                },
                indent=2,
            ),
            encoding="utf-8",
        )


def copy_file(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def copy_tree(source: Path, target: Path) -> None:
    if not source.is_dir():
        raise RuntimeError(f"Import directory not found: {source}")
    for path in sorted(source.rglob("*")):
        if path.is_file():
            copy_file(path, target / path.relative_to(source))


def write_sample_design(row: dict[str, str], target: Path) -> None:
    fields = (
        "sample_id", "chip_id", "species", "chemistry", "tissue_type",
        "experimental_group", "label_4su_dose", "label_4su_duration",
        "harvest_time", "unlabelled_control", "biological_replicate",
        "section_id", "tissue_image",
    )
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8") as handle:
        handle.write("field\tvalue\n")
        for field in fields:
            handle.write(f"{field}\t{row.get(field, '')}\n")


def report_path(outdir: Path, sample: str) -> Path:
    return outdir / f"{sample}_report.html"


def json_load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def make_dirs(stage: Path) -> None:
    for name in ("01_tables", "02_figures", "03_logs", "04_tool_output", "05_metadata"):
        (stage / name).mkdir(parents=True, exist_ok=True)


def require_ready(paths: list[Path], stage: str) -> None:
    missing = [str(path) for path in paths if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise RuntimeError(f"{stage} requires upstream outputs:\n" + "\n".join(missing))


def record_visual_refresh(
    runner: Runner,
    stage: Path,
    sources: list[Path],
    outputs: list[Path],
    font_family: str,
    base_size: str,
    dpi: str,
) -> None:
    metadata = stage / "05_metadata"
    runner.current_stage = stage.name
    runner.current_stage_root = stage
    runner.run(
        [
            sys.executable,
            str(runner.root / "scripts" / "capture_git_provenance.py"),
            "--repo-root",
            str(runner.root),
            "--output",
            str(metadata / "visual_git_provenance.tsv"),
            "--dirty-policy",
            runner.config.get("git_dirty_policy", "warn"),
        ],
        "visual_git_provenance",
    )
    payload = {
        "stage": stage.name,
        "refreshed_at_utc": datetime.now(timezone.utc).isoformat(),
        "style": {
            "name": "nature_npg",
            "font_family": font_family,
            "base_size": float(base_size),
            "dpi": int(dpi),
        },
        "sources": [
            {
                "path": str(path.resolve()),
                "size_bytes": path.stat().st_size,
                "mtime_ns": path.stat().st_mtime_ns,
            }
            for path in sources
        ],
        "outputs": [str(path.resolve()) for path in outputs],
    }
    (metadata / "visual_refresh.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def refresh_completion_metadata(runner: Runner, stage: Path, outputs: list[Path]) -> None:
    metadata = stage / "05_metadata"
    with (metadata / "output_manifest.tsv").open("w", encoding="utf-8") as handle:
        handle.write("path\tsize_bytes\tmtime_ns\n")
        for path in outputs:
            stat = path.stat()
            handle.write(f"{path.resolve()}\t{stat.st_size}\t{stat.st_mtime_ns}\n")
    (metadata / "stage_complete.json").write_text(
        json.dumps(
            {
                "stage": stage.name,
                "completed_at_utc": datetime.now(timezone.utc).isoformat(),
                "completion_mode": "visual_refresh",
                "config_sha256": runner.config_sha256,
                "outputs": [str(path.resolve()) for path in outputs],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def refresh_visuals(args: argparse.Namespace, config: dict[str, str], row: dict[str, str]) -> None:
    root = Path.cwd().resolve()
    scripts = root / "scripts"
    sample = row["sample_id"]
    sample_root = Path(args.results_dir) / sample
    work = Path(args.tmp_dir) / sample
    runner = Runner(
        root,
        Path(args.config),
        config,
        Path(args.logs_dir) / sample,
        work,
        set(),
    )
    analysis_py = config.get("analysis_python", sys.executable)
    promo_py = config.get("spatial_python", analysis_py)
    spatial_py = config.get("spatial_python", analysis_py)
    font_family = config.get("plot_font_family", "DejaVu Sans")
    base_size = config.get("plot_base_size", "10.5")
    dpi = config.get("plot_dpi", "300")
    style_args = ["--font-family", font_family, "--base-size", base_size, "--dpi", dpi]
    stages = {name: sample_root / name for name in STAGES}
    for stage in stages.values():
        make_dirs(stage)

    s1 = stages["01_read_qc"]
    read_source = s1 / "04_tool_output" / "read_qc" / "per_read_qc.tsv.gz"
    require_ready([read_source], "visual refresh 01_read_qc")
    runner.current_stage_root = s1
    runner.current_stage = s1.name
    runner.run(
        [
            analysis_py,
            str(scripts / "lr_plot_core_qc.py"),
            "--section",
            "read",
            "--qc-per-read",
            str(read_source),
            "--outdir",
            str(s1 / "02_figures"),
            *style_args,
        ],
        "refresh_read_qc",
    )
    read_outputs = [s1 / "02_figures/raw_read_length.png", s1 / "02_figures/raw_mean_q.png"]
    record_visual_refresh(runner, s1, [read_source], read_outputs, font_family, base_size, dpi)

    s2 = stages["02_barcode_cdna"]
    structure_source = s2 / "04_tool_output" / "structure" / "read_structure.tsv.gz"
    structure_summary = s2 / "05_metadata" / "structure_summary.json"
    barcode_summary = s2 / "05_metadata" / "barcode_summary.json"
    require_ready([structure_source, structure_summary, barcode_summary], "visual refresh 02_barcode_cdna")
    runner.current_stage_root = s2
    runner.current_stage = s2.name
    runner.run(
        [
            analysis_py,
            str(scripts / "lr_plot_core_qc.py"),
            "--section",
            "barcode",
            "--structure-tsv",
            str(structure_source),
            "--structure-summary",
            str(structure_summary),
            "--barcode-summary",
            str(barcode_summary),
            "--max-hist-rows",
            config.get("plot_hist_max_rows", "500000"),
            "--outdir",
            str(s2 / "02_figures"),
            *style_args,
        ],
        "refresh_barcode_qc",
    )
    barcode_outputs = [
        s2 / "02_figures/anchor_linker_status.png",
        s2 / "02_figures/structure_status.png",
        s2 / "02_figures/cdna_length.png",
        s2 / "02_figures/barcode_status.png",
        s2 / "02_figures/barcode_edit_distance.png",
    ]
    record_visual_refresh(
        runner,
        s2,
        [structure_source, structure_summary, barcode_summary],
        barcode_outputs,
        font_family,
        base_size,
        dpi,
    )
    s3 = stages["03_align"]
    compact_barcode = s3 / "01_tables" / "barcode_level_summary.tsv"
    candidates = [
        compact_barcode,
        work / "03_align" / "core_stats" / "barcode_level_summary.tsv",
        work / "04_align" / "core_stats" / "barcode_level_summary.tsv",
    ]
    barcode_level = next((path for path in candidates if path.is_file() and path.stat().st_size > 0), None)
    if barcode_level is None:
        raise RuntimeError("visual refresh requires barcode_level_summary.tsv")
    if barcode_level != compact_barcode:
        copy_file(barcode_level, compact_barcode)
    alignment_summary = s3 / "05_metadata" / "alignment_summary.json"
    require_ready([compact_barcode, alignment_summary], "visual refresh 03_align")
    runner.current_stage_root = s3
    runner.current_stage = s3.name
    runner.run(
        [
            analysis_py,
            str(scripts / "lr_plot_core_qc.py"),
            "--section",
            "align",
            "--barcode-level-summary",
            str(compact_barcode),
            "--alignment-summary",
            str(alignment_summary),
            "--outdir",
            str(s3 / "02_figures"),
            *style_args,
        ],
        "refresh_alignment_qc",
    )
    align_outputs = [s3 / "02_figures/alignment_summary.png", s3 / "02_figures/top50_barcodes.png"]
    record_visual_refresh(
        runner,
        s3,
        [compact_barcode, alignment_summary],
        align_outputs,
        font_family,
        base_size,
        dpi,
    )

    s4 = stages["04_full_length_transcriptome"]
    tmap_candidates = sorted((s4 / "04_tool_output/flames").glob("*.tmap"))
    if not tmap_candidates:
        raise RuntimeError("visual refresh requires a FLAMES gffcompare TMAP")
    runner.current_stage_root = s4
    runner.current_stage = s4.name
    runner.run(
        [
            analysis_py,
            str(scripts / "summarize_isoforms.py"),
            "--tmap",
            str(tmap_candidates[0]),
            "--outdir",
            str(s4),
            *style_args,
        ],
        "refresh_isoform_plots",
    )
    transcript_outputs = sorted((s4 / "02_figures").glob("*.png"))
    record_visual_refresh(
        runner, s4, [tmap_candidates[0]], transcript_outputs, font_family, base_size, dpi
    )

    s5 = stages["05_slam"]
    spectrum = s5 / "04_tool_output" / "mismatch_spectrum" / "mismatch_spectrum.tsv"
    mismatch_summary = s5 / "04_tool_output" / "mismatch_spectrum" / "mismatch_summary.json"
    gene_summary = s5 / "01_tables" / "gene_new_supported_summary.tsv"
    class_summary = s5 / "05_metadata" / f"{sample}.slam_summary.json"
    require_ready([spectrum, mismatch_summary, gene_summary, class_summary], "visual refresh 05_slam")
    runner.current_stage_root = s5
    runner.current_stage = s5.name
    runner.run(
        [
            analysis_py,
            str(scripts / "plot_mismatch_spectrum.py"),
            "--spectrum",
            str(spectrum),
            "--summary",
            str(mismatch_summary),
            "--outdir",
            str(s5),
            "--figure-only",
            *style_args,
        ],
        "refresh_mismatch_plot",
    )
    runner.run(
        [
            analysis_py,
            str(scripts / "summarize_gene_slam.py"),
            "--summary-table",
            str(gene_summary),
            "--class-summary",
            str(class_summary),
            "--outdir",
            str(s5),
            "--min-informative",
            config["gene_slam_min_informative_molecules"],
            "--top",
            config["gene_slam_plot_top"],
            *style_args,
        ],
        "refresh_gene_slam_plots",
    )
    slam_outputs = sorted((s5 / "02_figures").glob("*.png"))
    record_visual_refresh(
        runner,
        s5,
        [spectrum, mismatch_summary, gene_summary, class_summary],
        slam_outputs,
        font_family,
        base_size,
        dpi,
    )

    s6 = stages["06_whole_transcriptome"]
    spatial_sources = [
        s6 / "04_tool_output/scanpy/long_read_spatial_bin50.h5ad",
        s6 / "01_tables/long_read_cluster_markers.tsv",
        s6 / "01_tables/long_read_spatial_coherence.tsv",
    ]
    require_ready(spatial_sources, "visual refresh 06_whole_transcriptome")
    runner.current_stage_root = s6
    runner.current_stage = s6.name
    runner.run(
        [
            spatial_py,
            str(scripts / "cluster_lr_spatial.py"),
            "--outdir",
            str(s6),
            "--render-only",
            "--marker-display-top",
            config["lr_spatial_marker_display_top"],
            *style_args,
        ],
        "refresh_spatial_plots",
    )
    spatial_outputs = sorted((s6 / "02_figures").glob("*.png"))
    record_visual_refresh(
        runner,
        s6,
        spatial_sources,
        spatial_outputs,
        font_family,
        base_size,
        dpi,
    )

    s7 = stages["07_slam_joint"]
    joint_sources = [
        s7 / "01_tables/cluster_slam_summary.tsv",
        s7 / "01_tables/cluster_gene_slam.tsv",
        s6 / "01_tables/long_read_bin50_clusters.tsv",
    ]
    require_ready(joint_sources, "visual refresh 07_slam_joint")
    runner.current_stage_root = s7
    runner.current_stage = s7.name
    runner.run(
        [
            analysis_py,
            str(scripts / "analyze_spatial_slam_joint.py"),
            "--clusters",
            str(joint_sources[2]),
            "--sample",
            sample,
            "--outdir",
            str(s7),
            "--bin-size",
            config["bin_size"],
            "--min-informative",
            config["joint_min_informative_molecules"],
            "--gene-plot-top",
            config["joint_gene_plot_top"],
            "--render-only",
            *style_args,
        ],
        "refresh_joint_plots",
    )
    joint_outputs = sorted((s7 / "02_figures").glob("*.png"))
    record_visual_refresh(
        runner, s7, joint_sources, joint_outputs, font_family, base_size, dpi
    )

    s8 = stages["08_report"]
    runner.current_stage_root = s8
    runner.current_stage = s8.name
    runner.run(
        [
            analysis_py,
            str(scripts / "00_shared" / "repair_result_metadata.py"),
            "--sample-root",
            str(sample_root),
            "--sample-id",
            sample,
        ],
        "repair_result_metadata",
    )
    report_title = args.report_title or config.get("report_title", "长读长空间 SLAM 分析报告")
    runner.run(
        [
            analysis_py,
            str(scripts / "build_lr_spatial_slam_report.py"),
            "--sample-id",
            sample,
            "--chip-id",
            row["chip_id"],
            "--species",
            row["species"],
            "--chemistry",
            row["chemistry"],
            "--sample-dir",
            str(sample_root),
            "--report-title",
            report_title,
            "--outdir",
            str(s8),
        ],
        "refresh_report",
    )
    report = report_path(s8, sample)
    delivery = s8 / "05_metadata/delivery_manifest.tsv"
    runner.run(
        [
            sys.executable,
            str(scripts / "validate_delivery_report.py"),
            "--report",
            str(report),
            "--manifest",
            str(delivery),
            "--forbid-text",
            "参数评估",
            "--forbid-text",
            "long_read_resolution_scan",
        ],
        "validate_refreshed_report",
    )
    record_visual_refresh(
        runner,
        s8,
        transcript_outputs + slam_outputs + spatial_outputs + joint_outputs + align_outputs + barcode_outputs + read_outputs,
        [report, delivery],
        font_family,
        base_size,
        dpi,
    )
    refresh_completion_metadata(runner, s8, [report, delivery])


def build_promo_figures(args: argparse.Namespace, config: dict[str, str], row: dict[str, str]) -> None:
    root = Path.cwd().resolve()
    sample = row["sample_id"]
    sample_root = Path(args.results_dir) / sample
    dirname = config.get("promo_figures_dirname", "09_promo_figures")
    outdir = sample_root / dirname
    outdir.mkdir(parents=True, exist_ok=True)
    runner = Runner(root, Path(args.config), config, Path(args.logs_dir) / sample, Path(args.tmp_dir) / sample, set())
    runner.current_stage = dirname
    runner.current_stage_root = outdir
    analysis_py = config.get("analysis_python", sys.executable)
    promo_py = config.get("spatial_python", analysis_py)
    runner.run(
        [
            promo_py,
            str(root / "scripts" / "build_promo_figure_pack.py"),
            "--sample-dir",
            str(sample_root),
            "--sample-id",
            sample,
            "--outdir",
            str(outdir),
            "--dpi",
            config.get("plot_dpi", "300"),
            "--candidate-min-informative",
            config.get("promo_candidate_min_informative", "1000"),
        ],
        "build_promo_figures",
    )
    if truth(config.get("promo_isoform_spatial_enabled", "false")):
        runner.run(
            [
                promo_py,
                str(root / "scripts" / "09_promo_figures" / "build_isoform_spatial_candidates.py"),
                "--sample-dir", str(sample_root),
                "--sample-id", sample,
                "--outdir", str(outdir),
                "--dpi", config.get("plot_dpi", "300"),
                "--min-isoform-molecules", config.get("promo_isoform_min_molecules", "100"),
                "--max-genes-to-map", config.get("promo_isoform_max_genes_to_map", "25"),
                "--min-dominant-difference", config.get("promo_isoform_min_dominant_difference", "0.25"),
                "--min-js-divergence", config.get("promo_isoform_min_js_divergence", "0.15"),
            ],
            "build_isoform_spatial_candidates",
        )


def run_sample(args: argparse.Namespace, config: dict[str, str], row: dict[str, str], selected: list[str]) -> None:
    root = Path.cwd().resolve()
    scripts = root / "scripts"
    sample = row["sample_id"]
    sample_root = Path(args.results_dir) / sample
    work = Path(args.tmp_dir) / sample
    runner = Runner(root, Path(args.config), config, Path(args.logs_dir) / sample, work, set(args.force_stage))
    threads = int(row.get("threads") or args.threads)
    analysis_py = config.get("analysis_python", sys.executable)
    adapter_py = config.get("adapter_python", analysis_py)
    spatial_py = config.get("spatial_python", analysis_py)
    run_mode = config.get("run_mode", "full")
    min_mapq = int(config["min_mapq"])
    style_args = [
        "--font-family",
        config.get("plot_font_family", "DejaVu Sans"),
        "--base-size",
        config.get("plot_base_size", "10.5"),
    ]
    active_fastq = work / "01_read_qc" / f"{sample}.pilot.fastq.gz" if run_mode == "pilot" else Path(row["fastq"])
    for stage_name in STAGES:
        make_dirs(sample_root / stage_name)

    s1 = sample_root / "01_read_qc"
    qc_tool = s1 / "04_tool_output" / "read_qc"
    input_manifest = s1 / "05_metadata" / "input_manifest.tsv"
    def stage1() -> None:
        if run_mode == "pilot":
            runner.run([analysis_py, str(scripts / "subsample_fastq.py"), "--fastq", row["fastq"], "--reads", config.get("pilot_reads", "200000"), "--out", str(active_fastq)], "subsample_fastq")
        with input_manifest.open("w", encoding="utf-8") as handle:
            handle.write("role\tpath\tsize_bytes\tmtime_ns\n")
            for role in ("fastq", "genome_fa", "gtf", "whitelist_csv"):
                path = Path(row[role])
                stat = path.stat()
                handle.write(f"{role}\t{path.resolve()}\t{stat.st_size}\t{stat.st_mtime_ns}\n")
        copy_file(Path(args.config), s1 / "05_metadata" / "run_config.yaml")
        copy_file(Path(args.samples), s1 / "05_metadata" / "samples.tsv")
        write_sample_design(row, s1 / "05_metadata" / "sample_design.tsv")
        selection = Path(config.get("software_selection", "conf/software_selection.lr_celatlas_spatial_slam.yaml"))
        if selection.is_file():
            copy_file(selection, s1 / "05_metadata" / "software_selection.yaml")
        command = [analysis_py, str(scripts / "lr_qc_raw_fastq.py"), "--fastq", str(active_fastq), "--outdir", str(qc_tool), "--threads", str(threads), "--skip-fastqc"]
        if not truth(config.get("enable_nanoplot", "false")):
            command.append("--skip-nanoplot")
        runner.run(command, "read_qc")
        copy_file(qc_tool / "qc_metrics.tsv", s1 / "01_tables" / "qc_metrics.tsv")
        runner.run([analysis_py, str(scripts / "lr_plot_core_qc.py"), "--section", "read", "--qc-per-read", str(qc_tool / "per_read_qc.tsv.gz"), "--outdir", str(s1 / "02_figures"), "--dpi", config.get("plot_dpi", "300"), *style_args], "plot_read_qc")
    if "01_read_qc" in selected:
        outputs = [input_manifest, s1 / "05_metadata" / "sample_design.tsv", s1 / "01_tables" / "qc_metrics.tsv", qc_tool / "per_read_qc.tsv.gz", s1 / "02_figures" / "raw_read_length.png"]
        if run_mode == "pilot":
            outputs.append(active_fastq)
        runner.stage("01_read_qc", outputs, stage1)

    s2 = sample_root / "02_barcode_cdna"
    structure_dir = s2 / "04_tool_output" / "structure"
    barcode_dir = s2 / "04_tool_output" / "barcode"
    cdna_dir = work / "02_barcode_cdna" / "cdna"
    cdna_fastq = cdna_dir / f"{sample}.cdna.fastq.gz"
    per_read = cdna_dir / "per_read.tsv.gz"
    def stage2() -> None:
        require_ready([active_fastq], "02_barcode_cdna")
        if not ready([structure_dir / "read_structure.tsv.gz", structure_dir / "structure_summary.json"]):
            runner.run([adapter_py, str(scripts / "lr_scan_read_structure.py"), "--fastq", str(active_fastq), "--outdir", str(structure_dir), "--sample-id", sample, "--anchor-max-ed-5p", config.get("anchor_max_ed_5p", "5"), "--linker-max-ed", config.get("linker_max_ed", "4"), "--linker-pos-tolerance", config.get("linker_pos_tolerance", "8"), "--umi-len", config.get("umi_len", "10"), "--search-window", config.get("search_window", "300"), "--poly-trim-min-run", config.get("poly_trim_min_run", "10"), "--min-read-len", config.get("min_read_len", "80"), "--workers", config.get("structure_scan_workers", "4"), "--batch-size", config.get("structure_scan_batch_size", "1000")], "scan_structure")
        if not ready([barcode_dir / "barcode_assign.tsv.gz", barcode_dir / "barcode_summary.json"]):
            runner.run([adapter_py, str(scripts / "lr_assign_barcodes.py"), "--structure-tsv", str(structure_dir / "read_structure.tsv.gz"), "--whitelist-csv", row["whitelist_csv"], "--outdir", str(barcode_dir), "--anchor-barcode-max-mismatch", config.get("anchor_barcode_max_mismatch", "1"), "--fallback-barcode-max-mismatch", config.get("fallback_barcode_max_mismatch", "1")], "assign_barcodes")
        runner.run([adapter_py, str(scripts / "lr_extract_cdna_fastq.py"), "--fastq", str(active_fastq), "--structure-tsv", str(structure_dir / "read_structure.tsv.gz"), "--barcode-tsv", str(barcode_dir / "barcode_assign.tsv.gz"), "--sample-id", sample, "--outdir", str(cdna_dir)], "extract_cdna")
        for source, name in ((structure_dir / "structure_summary.json", "structure_summary.json"), (barcode_dir / "barcode_summary.json", "barcode_summary.json"), (cdna_dir / "extract_summary.json", "extract_summary.json")):
            copy_file(source, s2 / "05_metadata" / name)
        runner.run([analysis_py, str(scripts / "lr_plot_core_qc.py"), "--section", "barcode", "--structure-tsv", str(structure_dir / "read_structure.tsv.gz"), "--structure-summary", str(structure_dir / "structure_summary.json"), "--barcode-summary", str(barcode_dir / "barcode_summary.json"), "--max-hist-rows", config.get("plot_hist_max_rows", "500000"), "--outdir", str(s2 / "02_figures"), "--dpi", config.get("plot_dpi", "300"), *style_args], "plot_barcode_qc")
    if "02_barcode_cdna" in selected:
        runner.stage("02_barcode_cdna", [s2 / "05_metadata" / "extract_summary.json", cdna_fastq, per_read, s2 / "02_figures" / "barcode_status.png"], stage2)

    s3 = sample_root / "03_align"
    align_dir = s3 / "04_tool_output" / "minimap2"
    stats_dir = work / "03_align" / "core_stats"
    align_work = work / "03_align" / "minimap2"
    shared_ref = Path(config.get("tmp_root", "tmp")) / "reference" / row["species"]
    mmi = Path(row.get("mmi") or shared_ref / (Path(row["genome_fa"]).name + ".mmi"))
    primary_bam = align_dir / f"{sample}.primary.q{min_mapq}.bam"
    quality_gate = s3 / "05_metadata" / "quality_gate.json"
    def stage3() -> None:
        require_ready([cdna_fastq, per_read, structure_dir / "read_structure.tsv.gz", barcode_dir / "barcode_assign.tsv.gz"], "03_align")
        if not mmi.is_file():
            mmi.parent.mkdir(parents=True, exist_ok=True)
            runner.run([config.get("minimap2_bin", "minimap2"), "-d", str(mmi), row["genome_fa"]], "build_minimap2_index")
        command = [analysis_py, str(scripts / "lr_align_reads.py"), "--fastq", str(cdna_fastq), "--mmi", str(mmi), "--outdir", str(align_dir), "--workdir", str(align_work), "--sample-id", sample, "--threads", str(threads), "--min-mapq", str(min_mapq), "--mm2-preset", config.get("minimap2_preset", "splice"), "--minimap2-bin", config.get("minimap2_bin", "minimap2"), "--samtools-bin", config.get("samtools_bin", "samtools")]
        if row.get("junc_bed"):
            command.extend(["--junc-bed", row["junc_bed"]])
        runner.run(command, "align")
        runner.run([analysis_py, str(scripts / "lr_collect_core_stats.py"), "--structure-tsv", str(structure_dir / "read_structure.tsv.gz"), "--barcode-tsv", str(barcode_dir / "barcode_assign.tsv.gz"), "--per-read-tsv", str(per_read), "--bam", str(primary_bam), "--outdir", str(stats_dir)], "collect_stats")
        copy_file(stats_dir / "barcode_level_summary.tsv", s3 / "01_tables" / "barcode_level_summary.tsv")
        runner.run([analysis_py, str(scripts / "lr_plot_core_qc.py"), "--section", "align", "--stats-dir", str(stats_dir), "--outdir", str(s3 / "02_figures"), "--dpi", config.get("plot_dpi", "300"), *style_args], "plot_alignment_qc")
        runner.run([analysis_py, str(scripts / "compute_lr_quality_gate.py"), "--qc-metrics", str(s1 / "01_tables" / "qc_metrics.tsv"), "--structure-summary", str(structure_dir / "structure_summary.json"), "--extract-summary", str(cdna_dir / "extract_summary.json"), "--primary-flagstat", str(align_dir / f"{sample}.primary.flagstat.txt"), "--per-read", str(per_read), "--chip-bounds", config.get("chip_bounds", "9357"), "--min-structure-rate", config.get("pilot_min_structure_rate", "0.5"), "--min-cdna-rate", config.get("pilot_min_cdna_rate", "0.15"), "--min-mapq20-rate", config.get("pilot_min_mapq20_rate", "0.75"), "--output", str(quality_gate)], "compute_quality_gate")
        copy_file(stats_dir / "alignment_summary.json", s3 / "05_metadata" / "alignment_summary.json")
        if run_mode == "pilot" and json_load(quality_gate)["overall"] != "PASS":
            raise RuntimeError(f"Quality gate failed: {json_load(quality_gate)}")
    if "03_align" in selected:
        runner.stage("03_align", [primary_bam, primary_bam.with_suffix(".bam.bai"), quality_gate, s3 / "01_tables" / "barcode_level_summary.tsv", s3 / "02_figures" / "alignment_summary.png"], stage3)

    s4 = sample_root / "04_full_length_transcriptome"
    flames_dir = s4 / "04_tool_output" / "flames"
    compare_dir = s4 / "04_tool_output" / "gffcompare"
    compare_prefix = compare_dir / "flames_vs_reference"
    discovered = flames_dir / "isoforms.gff3"
    tmap = flames_dir / f"{compare_prefix.name}.{discovered.name}.tmap"
    def stage4() -> None:
        require_ready([primary_bam, primary_bam.with_suffix(".bam.bai")], "04_full_length_transcriptome")
        import_dir = config.get("transcriptome_existing_dir", "").strip()
        if import_dir and not ready([discovered, tmap, compare_prefix.with_suffix(".stats")]):
            source = Path(import_dir)
            copy_tree(source / "04_tool_output" / "flames", flames_dir)
            copy_tree(source / "04_tool_output" / "gffcompare", compare_dir)
        if not ready([discovered]):
            runner.run(
                [
                    analysis_py,
                    str(scripts / "run_flames_isoforms.py"),
                    "--bam",
                    str(primary_bam),
                    "--gtf",
                    row["gtf"],
                    "--genome",
                    row["genome_fa"],
                    "--flames-python-dir",
                    config["flames_python_dir"],
                    "--outdir",
                    str(flames_dir),
                    "--min-support",
                    config["flames_min_support"],
                    "--min-count-fraction",
                    config["flames_min_count_fraction"],
                    "--min-support-fraction",
                    config["flames_min_support_fraction"],
                    "--strand-specific",
                    config["flames_strand_specific"],
                    "--downsample-ratio",
                    config["flames_downsample_ratio"],
                    "--threads",
                    str(threads),
                ],
                "flames_isoform_discovery",
            )
        compare_dir.mkdir(parents=True, exist_ok=True)
        if not ready([compare_prefix.with_suffix(".stats"), tmap]):
            runner.run(
                [
                    config["gffcompare_bin"],
                    "-r",
                    row["gtf"],
                    "-o",
                    str(compare_prefix),
                    str(discovered),
                ],
                "gffcompare",
            )
        runner.run(
            [
                analysis_py,
                str(scripts / "00_shared" / "summarize_gffcompare.py"),
                "--stats",
                str(compare_prefix.with_suffix(".stats")),
                "--output",
                str(s4 / "01_tables" / "isoform_summary.tsv"),
            ],
            "summarize_gffcompare",
        )
        runner.run(
            [
                analysis_py,
                str(scripts / "summarize_isoforms.py"),
                "--tmap",
                str(tmap),
                "--outdir",
                str(s4),
                "--dpi",
                config["plot_dpi"],
                *style_args,
            ],
            "summarize_isoforms",
        )
        (s4 / "05_metadata" / "transcriptome_method.json").write_text(
            json.dumps(
                {
                    "discovery": "FLAMES",
                    "reference_comparison": "gffcompare",
                    "imported_existing_output": bool(import_dir),
                    "parameters_source": "run configuration",
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
    if "04_full_length_transcriptome" in selected:
        runner.stage(
            "04_full_length_transcriptome",
            [
                discovered,
                s4 / "01_tables" / "isoform_summary.tsv",
                s4 / "01_tables" / "isoform_class_summary.tsv",
                s4 / "02_figures" / "isoform_class_distribution.png",
                s4 / "05_metadata" / "transcriptome_method.json",
            ],
            stage4,
        )

    s5 = sample_root / "05_slam"
    mismatch_dir = s5 / "04_tool_output" / "mismatch_spectrum"
    calls_dir = s5 / "04_tool_output" / "caller"
    annotation_dir = work / "05_slam" / "annotation"
    barcode_input = work / "05_slam" / "barcode"
    slam_calls = s5 / "01_tables" / f"{sample}.slam_molecule_calls.tsv.gz"
    slam_summary = s5 / "05_metadata" / f"{sample}.slam_summary.json"
    def stage5() -> None:
        require_ready([primary_bam, primary_bam.with_suffix(".bam.bai"), per_read], "05_slam")
        force = "05_slam" in runner.force
        variant_vcf = config.get("external_variant_vcf", "").strip()
        mismatch_command = [analysis_py, str(scripts / "summarize_bam_mismatches.py"), "--bam", str(primary_bam), "--fasta", row["genome_fa"], "--outdir", str(mismatch_dir), "--min-mapq", str(min_mapq), "--min-baseq", config["min_tc_baseq"], "--max-reads", config["mismatch_max_reads"]]
        if variant_vcf:
            mismatch_command.extend(["--exclude-vcf", variant_vcf])
        if force or not ready([mismatch_dir / "mismatch_spectrum.tsv", mismatch_dir / "mismatch_summary.json"]):
            runner.run(mismatch_command, "mismatch_spectrum")
        runner.run([analysis_py, str(scripts / "plot_mismatch_spectrum.py"), "--spectrum", str(mismatch_dir / "mismatch_spectrum.tsv"), "--summary", str(mismatch_dir / "mismatch_summary.json"), "--outdir", str(s5), "--dpi", config["plot_dpi"], *style_args], "plot_mismatch_spectrum")
        background = str(json_load(mismatch_dir / "mismatch_summary.json")["background_rate"])
        raw_calls = calls_dir / f"{sample}.molecule_calls.tsv.gz"
        if force or not ready([raw_calls]):
            runner.run([analysis_py, str(scripts / "build_exon_intron_annotation.py"), "--gtf", row["gtf"], "--outdir", str(annotation_dir)], "build_annotation")
            runner.run([analysis_py, str(scripts / "make_slam_barcode_table.py"), "--per-read-tsv", str(per_read), "--sample-id", sample, "--outdir", str(barcode_input)], "slam_barcode_table")
            command = [analysis_py, str(scripts / "call_new_old_molecules.py"), "--bam", str(primary_bam), "--fasta", row["genome_fa"], "--annotation-dir", str(annotation_dir), "--barcode-table", str(barcode_input / f"{sample}.slam_barcode_table.tsv.gz"), "--sample", sample, "--outdir", str(calls_dir), "--main-call", "umi_site_binomial", "--min-tc-baseq", config["min_tc_baseq"], "--min-mapq", str(min_mapq), "--min-t-covered", config["min_t_covered"], "--bg-tc-rate", background, "--alpha", config["slam_fdr"], "--max-site-other-mismatch-rate", config["max_site_other_mismatch_rate"], "--max-read-other-mismatch-rate", config["max_read_other_mismatch_rate"]]
            if variant_vcf:
                command.extend(["--exclude-vcf", variant_vcf])
            if truth(config.get("slam_partition_by_contig", "true")):
                command.append("--partition-by-contig")
            runner.run(command, "call_slam")
        if force or not ready([slam_calls, s5 / "01_tables" / f"{sample}.slam_summary.json"]):
            runner.run([analysis_py, str(scripts / "reclassify_slam_calls.py"), "--molecule-calls", str(raw_calls), "--sample", sample, "--outdir", str(s5 / "01_tables"), "--min-t-covered", config["min_t_covered"], "--min-tc-count", config["min_tc_count"], "--fdr", config["slam_fdr"], "--max-other-mismatch-rate", config["max_read_other_mismatch_rate"]], "classify_slam")
        runner.run([analysis_py, str(scripts / "summarize_gene_slam.py"), "--molecule-calls", str(slam_calls), "--outdir", str(s5), "--min-informative", config["gene_slam_min_informative_molecules"], "--top", config["gene_slam_plot_top"], "--chunk-rows", config["lr_spatial_chunk_rows"], "--dpi", config["plot_dpi"], *style_args], "summarize_gene_slam")
        copy_file(s5 / "01_tables" / f"{sample}.slam_summary.json", slam_summary)
        copy_file(mismatch_dir / "mismatch_summary.json", s5 / "05_metadata" / "mismatch_summary.json")
    if "05_slam" in selected:
        runner.stage("05_slam", [slam_calls, slam_summary, s5 / "02_figures" / "base_substitution_comparison.png", s5 / "01_tables" / "gene_new_supported_summary.tsv"], stage5)

    s6 = sample_root / "06_whole_transcriptome"
    def stage6() -> None:
        require_ready([slam_calls], "06_whole_transcriptome")
        runner.run([spatial_py, str(scripts / "cluster_lr_spatial.py"), "--molecule-calls", str(slam_calls), "--sample", sample, "--outdir", str(s6), "--bin-size", config["bin_size"], "--coordinate-um-per-unit", config["coordinate_um_per_unit"], "--chip-bounds", config["chip_bounds"], "--min-bin-molecules", config["lr_spatial_min_bin_molecules"], "--min-gene-bins", config["lr_spatial_min_gene_bins"], "--hvg", config["lr_spatial_hvg"], "--pca-components", config["lr_spatial_pca_components"], "--neighbors", config["lr_spatial_neighbors"], "--regress-log-depth", config["lr_spatial_regress_log_depth"], "--spatial-graph-weight", config["lr_spatial_graph_weight"], "--min-cluster-bins", config["lr_spatial_min_cluster_bins"], "--leiden-resolution", config["lr_spatial_leiden_resolution"], "--seed", config["lr_spatial_seed"], "--marker-top", config["lr_spatial_marker_top"], "--marker-display-top", config["lr_spatial_marker_display_top"], "--marker-min-logfoldchange", config["lr_spatial_marker_min_logfoldchange"], "--marker-max-fdr", config["lr_spatial_marker_max_fdr"], "--chunk-rows", config["lr_spatial_chunk_rows"], "--dpi", config["plot_dpi"], *style_args], "whole_transcriptome_clustering")
        runner.run([analysis_py, str(scripts / "06_whole_transcriptome" / "summarize_spatial_qc.py"), "--molecule-calls", str(slam_calls), "--clusters", str(s6 / "01_tables" / "long_read_bin50_clusters.tsv"), "--isoform-summary", str(s4 / "01_tables" / "isoform_summary.tsv"), "--outdir", str(s6), "--bin-size", config["bin_size"], "--chunk-rows", config["lr_spatial_chunk_rows"]], "summarize_spatial_qc")
    if "06_whole_transcriptome" in selected:
        runner.stage("06_whole_transcriptome", [s6 / "01_tables" / "long_read_bin50_clusters.tsv", s6 / "01_tables" / "whole_transcriptome_gene_summary.tsv", s6 / "01_tables" / "spatial_bin_qc.tsv", s6 / "01_tables" / "spatial_qc_summary.tsv", s6 / "05_metadata" / "spatial_qc_summary.json", s6 / "05_metadata" / "long_read_spatial_clustering.json", s6 / "05_metadata" / "clustering_qc.json", s6 / "02_figures" / "long_read_spatial_clusters.png", s6 / "02_figures" / "cluster_marker_dotplot.png"], stage6)

    s7 = sample_root / "07_slam_joint"
    def stage7() -> None:
        require_ready([slam_calls, s6 / "01_tables" / "long_read_bin50_clusters.tsv"], "07_slam_joint")
        runner.run([analysis_py, str(scripts / "analyze_spatial_slam_joint.py"), "--molecule-calls", str(slam_calls), "--clusters", str(s6 / "01_tables" / "long_read_bin50_clusters.tsv"), "--sample", sample, "--outdir", str(s7), "--bin-size", config["bin_size"], "--min-informative", config["joint_min_informative_molecules"], "--gene-plot-top", config["joint_gene_plot_top"], "--chunk-rows", config["lr_spatial_chunk_rows"], "--dpi", config["plot_dpi"], *style_args], "slam_joint_analysis")
    if "07_slam_joint" in selected:
        runner.stage("07_slam_joint", [s7 / "01_tables" / "cluster_slam_summary.tsv", s7 / "01_tables" / "cluster_gene_slam.tsv", s7 / "02_figures" / "cluster_slam_composition.png", s7 / "02_figures" / "cluster_gene_slam_dotplot.png", s7 / "05_metadata" / "slam_joint_summary.json"], stage7)

    s8 = sample_root / "08_report"
    def stage8() -> None:
        require_ready([s1 / "01_tables" / "qc_metrics.tsv", s2 / "05_metadata" / "extract_summary.json", s3 / "05_metadata" / "alignment_summary.json", s4 / "05_metadata" / "transcriptome_method.json", slam_summary, s6 / "05_metadata" / "long_read_spatial_clustering.json", s7 / "05_metadata" / "slam_joint_summary.json"], "08_report")
        runner.run([analysis_py, str(scripts / "00_shared" / "repair_result_metadata.py"), "--sample-root", str(sample_root), "--sample-id", sample], "repair_result_metadata")
        report_title = args.report_title or config.get("report_title", "长读长空间 SLAM 分析报告")
        runner.run([analysis_py, str(scripts / "build_lr_spatial_slam_report.py"), "--sample-id", sample, "--chip-id", row["chip_id"], "--species", row["species"], "--chemistry", row["chemistry"], "--sample-dir", str(sample_root), "--report-title", report_title, "--outdir", str(s8)], "build_report")
        validate = [sys.executable, str(scripts / "validate_delivery_report.py"), "--report", str(report_path(s8, sample)), "--manifest", str(s8 / "05_metadata" / "delivery_manifest.tsv")]
        for text in ("短读长结果", "分析模式", "SLAM 仅作相对比较", "样本身份需优先复核", "空间聚类已完成", "参数评估", "long_read_resolution_scan"):
            validate.extend(["--forbid-text", text])
        runner.run(validate, "validate_report")
    if "08_report" in selected:
        runner.stage("08_report", [report_path(s8, sample), s8 / "05_metadata" / "delivery_manifest.tsv"], stage8)
    print(f"LR-Celatlas spatial SLAM workflow completed: {sample_root}")


def main() -> None:
    args = args_parser()
    config = read_config(args.config)
    if args.refresh_visuals:
        for row in read_samples(args.samples):
            refresh_visuals(args, config, row)
            if args.promo_figures:
                build_promo_figures(args, config, row)
        return
    start, end = STAGES.index(args.from_stage), STAGES.index(args.to_stage)
    if start > end:
        raise SystemExit("--from must not be after --to")
    selected = STAGES[start:end + 1]
    for row in read_samples(args.samples):
        run_sample(args, config, row, selected)
        if args.promo_figures:
            build_promo_figures(args, config, row)


if __name__ == "__main__":
    main()
