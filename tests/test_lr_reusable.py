#!/usr/bin/env python3
"""Smoke tests for the reusable LR-Celatlas workflow contract."""

import ast
import base64
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_STAGES = (
    "01_read_qc",
    "02_barcode_cdna",
    "03_align",
    "04_full_length_transcriptome",
    "05_slam",
    "06_whole_transcriptome",
    "07_slam_joint",
    "08_report",
)
FORBIDDEN = (
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


def report_python() -> str | None:
    candidates = (
        Path(sys.executable),
    )
    for candidate in candidates:
        if not candidate.is_file():
            continue
        result = subprocess.run(
            [str(candidate), "-c", "import pandas"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
        if result.returncode == 0:
            return str(candidate)
    return None


def plot_python() -> str | None:
    for candidate in (
        Path(sys.executable),
    ):
        if not candidate.is_file():
            continue
        result = subprocess.run(
            [str(candidate), "-c", "import matplotlib, pandas"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
        if result.returncode == 0:
            return str(candidate)
    return None


def run(*command: str) -> None:
    subprocess.run(command, cwd=ROOT, check=True)


def write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def test_config_contract(work: Path) -> None:
    inputs = work / "inputs"
    for name in ("reads.fastq.gz", "whitelist.csv", "genome.fa", "genes.gtf"):
        write(inputs / name, "fixture\n")
    flames_python_dir = inputs / "flames_python"
    flames_python_dir.mkdir()
    samples = work / "samples.tsv"
    write(
        samples,
        "sample_id\tfastq\tchip_id\twhitelist_csv\tspecies\tchemistry\tgenome_fa\tgtf\n"
        f"demo\t{inputs/'reads.fastq.gz'}\tchip_demo\t{inputs/'whitelist.csv'}\t"
        f"Mus_musculus\tBBV2.4\t{inputs/'genome.fa'}\t{inputs/'genes.gtf'}\n",
    )
    config = work / "config.yaml"
    write(
        config,
        "\n".join(
            (
                "project_name: lr_reusable_smoke",
                "analysis_type: lr_celatlas_spatial_slam",
                "run_id: smoke",
                f"samples_tsv: {samples}",
                f"data_root: {work}",
                f"results_root: {work/'results'}",
                f"logs_root: {work/'logs'}",
                f"tmp_root: {work/'tmp'}",
                "sample_first_layout: true",
                "run_mode: full",
                "species: Mus_musculus",
                "chemistry: BBV2.4",
                f"adapter_python: {sys.executable}",
                f"analysis_python: {sys.executable}",
                f"spatial_python: {sys.executable}",
                f"minimap2_bin: {sys.executable}",
                f"samtools_bin: {sys.executable}",
                f"flames_python_dir: {flames_python_dir}",
                f"gffcompare_bin: {sys.executable}",
                "flames_min_support: 5",
                "flames_min_count_fraction: 0.001",
                "flames_min_support_fraction: 0.1",
                "flames_strand_specific: 0",
                "flames_downsample_ratio: 1.0",
                "min_mapq: 20",
                "min_tc_baseq: 20",
                "mismatch_max_reads: 1000",
                "min_t_covered: 10",
                "min_tc_count: 2",
                "slam_fdr: 0.05",
                "max_site_other_mismatch_rate: 0.05",
                "max_read_other_mismatch_rate: 0.1",
                "gene_slam_min_informative_molecules: 10",
                "gene_slam_plot_top: 10",
                "bin_size: 50",
                "coordinate_um_per_unit: 1.0",
                "chip_bounds: 9357",
                "lr_spatial_min_bin_molecules: 10",
                "lr_spatial_min_gene_bins: 2",
                "lr_spatial_hvg: 100",
                "lr_spatial_pca_components: 10",
                "lr_spatial_neighbors: 5",
                "lr_spatial_regress_log_depth: true",
                "lr_spatial_graph_weight: 0.2",
                "lr_spatial_min_cluster_bins: 2",
                "lr_spatial_leiden_resolution: 1.0",
                "lr_spatial_seed: 42",
                "lr_spatial_marker_top: 10",
                "lr_spatial_marker_display_top: 2",
                "lr_spatial_marker_min_logfoldchange: 0.25",
                "lr_spatial_marker_max_fdr: 0.05",
                "lr_spatial_chunk_rows: 1000",
                "joint_min_informative_molecules: 5",
                "joint_gene_plot_top: 10",
                "plot_dpi: 72",
                "threads: 1",
                "memory_gb: 2",
                "",
            )
        ),
    )
    run(sys.executable, "scripts/check_samples.py", "--samples", str(samples), "--config", str(config))
    run("bash", "workflow/run.sh", "--config", str(config), "--dry-run")
    overridden = subprocess.run(
        [
            "bash",
            "workflow/run.sh",
            "--config",
            str(config),
            "--dry-run",
            "--report-title",
            "自定义报告标题",
        ],
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=True,
    )
    if "report_title: 自定义报告标题" not in overridden.stdout:
        raise AssertionError("report title override was not applied to workflow metadata")
    incompatible = subprocess.run(
        [
            "bash",
            "workflow/run.sh",
            "--config",
            str(config),
            "--refresh-visuals",
            "--from",
            "03_align",
        ],
        cwd=ROOT,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if incompatible.returncode == 0 or "cannot be combined" not in incompatible.stderr:
        raise AssertionError("incompatible visual-refresh stage selectors unexpectedly passed")

    module = ast.parse((ROOT / "scripts/run_lr_spatial_slam.py").read_text(encoding="utf-8"))
    stage_assignment = next(
        node for node in module.body
        if isinstance(node, ast.Assign)
        and any(isinstance(target, ast.Name) and target.id == "STAGES" for target in node.targets)
    )
    if tuple(ast.literal_eval(stage_assignment.value)) != EXPECTED_STAGES:
        raise AssertionError("LR driver stage contract does not match the reusable eight-stage layout")

    deprecated_config = work / "deprecated_config.yaml"
    deprecated_config.write_text(
        config.read_text(encoding="utf-8") + "lr_spatial_resolution_grid: 0.8,1.0,1.2\n",
        encoding="utf-8",
    )
    deprecated = subprocess.run(
        [sys.executable, "scripts/check_samples.py", "--samples", str(samples), "--config", str(deprecated_config)],
        cwd=ROOT,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if deprecated.returncode == 0 or "deprecated" not in deprecated.stderr:
        raise AssertionError("deprecated automatic resolution scan parameter unexpectedly passed")

    invalid = work / "invalid.tsv"
    write(invalid, "sample_id\tfastq\nbad\tmissing.fastq\n")
    result = subprocess.run(
        [sys.executable, "scripts/check_samples.py", "--samples", str(invalid), "--config", str(config)],
        cwd=ROOT,
        check=False,
    )
    if result.returncode == 0:
        raise AssertionError("invalid LR sample sheet unexpectedly passed")


def test_report(work: Path) -> None:
    renderer = report_python()
    if renderer is None:
        print("SKIP: no configured Python environment with pandas for report rendering")
        return
    sample = "demo"
    sample_root = work / "report_fixture" / sample
    write(sample_root / "01_read_qc/01_tables/qc_metrics.tsv", "metric\tvalue\nn_reads\t10\nN50\t900\n")
    write(sample_root / "02_barcode_cdna/05_metadata/structure_summary.json", '{"candidate_found": 8}\n')
    write(sample_root / "02_barcode_cdna/05_metadata/barcode_summary.json", '{"assign_counts": {"exact": 6}}\n')
    write(sample_root / "02_barcode_cdna/05_metadata/extract_summary.json", '{"counts": {"kept": 6}}\n')
    write(sample_root / "03_align/05_metadata/quality_gate.json", '{"overall": "PASS", "mapq20_rate": 0.9, "coordinates_valid": true}\n')
    write(sample_root / "03_align/05_metadata/alignment_summary.json", '{"mapped": 6}\n')
    write(
        sample_root / "04_full_length_transcriptome/05_metadata/isoform_class_summary.json",
        '{"transcripts": 3, "class_counts": {"=": 2, "j": 1}}\n',
    )
    write(
        sample_root / "04_full_length_transcriptome/01_tables/isoform_class_summary.tsv",
        "class_code\tdescription\ttranscripts\tfraction\n"
        "=\tExact reference match\t2\t0.6667\n"
        "j\tJunction match, new combination\t1\t0.3333\n",
    )
    write(
        sample_root / "04_full_length_transcriptome/01_tables/isoform_summary.tsv",
        "metric\tvalue\nquery_transcripts\t3\nmulti_exon_transcripts\t2\n",
    )
    write(
        sample_root / f"05_slam/05_metadata/{sample}.slam_summary.json",
        '{"molecules_total": 5, "class_counts": {"new_supported": 2, '
        '"no_conversion_detected": 2, "ambiguous": 1}}\n',
    )
    write(sample_root / "05_slam/05_metadata/gene_slam_summary.json", '{"genes_passing_threshold": 1}\n')
    write(
        sample_root / "05_slam/01_tables/base_substitution_comparison.tsv",
        "channel\tcount\treference_coverage\trate_percent\texpected_slam_channel\nT>C\t2\t100\t2.0\ttrue\n",
    )
    write(
        sample_root / "05_slam/01_tables/gene_new_supported_summary.tsv",
        "gene_id\tgene_name\tnew_supported\tno_conversion_detected\tambiguous\tinformative_molecules\t"
        "new_supported_fraction\twilson95_low\twilson95_high\tpasses_min_informative\n"
        "ENSMUSG1\tGene1\t2\t2\t1\t4\t0.5\t0.15\t0.85\tTrue\n",
    )
    write(
        sample_root / "06_whole_transcriptome/05_metadata/long_read_spatial_clustering.json",
        '{"input_molecules": 5, "clustered_bins": 3, "genes_after_filter": 2, '
        '"clusters": 2, "bin_size_um": 50}\n',
    )
    write(
        sample_root / "06_whole_transcriptome/05_metadata/clustering_qc.json",
        '{"spatial_graph_weight": 0.2}\n',
    )
    write(
        sample_root / "06_whole_transcriptome/01_tables/long_read_cluster_summary.tsv",
        "cluster\tbins\ttotal_lr_molecules\tmedian_lr_molecules\tmedian_detected_genes\n"
        "1\t2\t3\t1.5\t1\n2\t1\t2\t2\t2\n",
    )
    write(
        sample_root / "06_whole_transcriptome/01_tables/long_read_cluster_markers.tsv",
        "cluster\trank\tgene_id\tgene_name\tlogfoldchange\tpvalue_adj\tpasses_display_filter\n"
        "1\t1\tENSMUSG1\tGene1\t1.2\t0.01\tTrue\n"
        "2\t1\tENSMUSG2\tGene2\t1.1\t0.02\tTrue\n",
    )
    write(
        sample_root / "07_slam_joint/05_metadata/slam_joint_summary.json",
        '{"molecules_joined_to_whole_transcriptome_bins": 5, "clusters": 2, '
        '"genes": 2, "min_informative_molecules": 2}\n',
    )
    write(
        sample_root / "07_slam_joint/01_tables/cluster_slam_summary.tsv",
        "cluster\tnew_supported\tno_conversion_detected\tambiguous\tinformative_molecules\t"
        "new_supported_fraction\twilson95_low\twilson95_high\n"
        "1\t1\t1\t1\t2\t0.5\t0.095\t0.905\n"
        "2\t1\t1\t0\t2\t0.5\t0.095\t0.905\n",
    )
    write(
        sample_root / "07_slam_joint/01_tables/cluster_gene_slam.tsv",
        "cluster\tgene_id\tgene_name\tnew_supported\tno_conversion_detected\tambiguous\t"
        "informative_molecules\tnew_supported_fraction\twilson95_low\twilson95_high\t"
        "passes_min_informative\n"
        "1\tENSMUSG1\tGene1\t1\t1\t0\t2\t0.5\t0.095\t0.905\tTrue\n"
        "2\tENSMUSG2\tGene2\t1\t1\t0\t2\t0.5\t0.095\t0.905\tTrue\n",
    )
    png = base64.b64decode(
        "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mNk+A8AAQUBAScY42YAAAAASUVORK5CYII="
    )
    figure_paths = [
        "01_read_qc/02_figures/raw_read_length.png",
        "01_read_qc/02_figures/raw_mean_q.png",
        "02_barcode_cdna/02_figures/structure_status.png",
        "02_barcode_cdna/02_figures/cdna_length.png",
        "02_barcode_cdna/02_figures/barcode_status.png",
        "02_barcode_cdna/02_figures/barcode_edit_distance.png",
        "03_align/02_figures/alignment_summary.png",
        "03_align/02_figures/top50_barcodes.png",
        "04_full_length_transcriptome/02_figures/isoform_class_distribution.png",
        "04_full_length_transcriptome/02_figures/isoform_structure_distribution.png",
        "05_slam/02_figures/base_substitution_comparison.png",
        "05_slam/02_figures/slam_class_composition.png",
        "05_slam/02_figures/gene_new_supported_fraction.png",
        "06_whole_transcriptome/02_figures/long_read_spatial_clusters.png",
        "06_whole_transcriptome/02_figures/long_read_umap_clusters.png",
        "06_whole_transcriptome/02_figures/long_read_molecules.png",
        "06_whole_transcriptome/02_figures/cluster_marker_dotplot.png",
        "06_whole_transcriptome/02_figures/cluster_depth_qc.png",
        "06_whole_transcriptome/02_figures/spatial_coherence.png",
        "07_slam_joint/02_figures/cluster_slam_composition.png",
        "07_slam_joint/02_figures/cluster_new_supported_fraction.png",
        "07_slam_joint/02_figures/cluster_gene_slam_dotplot.png",
    ]
    for relative in figure_paths:
        target = sample_root / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(png)
    outdir = sample_root / "08_report"
    run(
        renderer,
        "scripts/build_lr_spatial_slam_report.py",
        "--sample-id",
        sample,
        "--chip-id",
        "chip_demo",
        "--species",
        "Mus_musculus",
        "--chemistry",
        "BBV2.4",
        "--sample-dir",
        str(sample_root),
        "--report-title",
        "长读长空间 SLAM 分析报告",
        "--outdir",
        str(outdir),
    )
    run(
        sys.executable,
        "scripts/validate_delivery_report.py",
        "--report",
        str(outdir / f"{sample}_report.html"),
        "--manifest",
        str(outdir / "05_metadata/delivery_manifest.tsv"),
    )
    report_path = outdir / f"{sample}_report.html"
    report_text = report_path.read_text(encoding="utf-8")
    if report_text.count('class="figure-card') != len(figure_paths):
        raise AssertionError("report did not embed the complete figure set")
    for marker in ('class="figure-card wide"', 'class="figure-card panorama"', "object-fit:contain"):
        if marker not in report_text:
            raise AssertionError(f"optimized report layout marker is missing: {marker}")
    if f"<title>{sample} · 长读长空间 SLAM 分析报告</title>" not in report_text:
        raise AssertionError("report display title is missing the sample prefix")
    for section in (
        "5 全长转录组",
        "7 长读长空间转录组",
        "8 长读长空间转录组与 SLAM 联合",
    ):
        if section not in report_text:
            raise AssertionError(f"corrected report section is missing: {section}")
    if (outdir / "report.html").exists():
        raise AssertionError("legacy unprefixed report was retained")
    manifest_text = (outdir / "05_metadata/delivery_manifest.tsv").read_text(encoding="utf-8")
    if str(report_path.resolve()) not in manifest_text:
        raise AssertionError("delivery manifest does not identify the sample-prefixed report")
    for text in FORBIDDEN:
        if text in report_text:
            raise AssertionError(f"forbidden text in reusable report: {text}")


def test_no_automatic_parameter_scans(work: Path) -> None:
    renderer = plot_python()
    if renderer is None:
        print("SKIP: no configured Python environment with pandas and numpy")
        return
    cluster_source = (ROOT / "scripts/cluster_lr_spatial.py").read_text(encoding="utf-8")
    if "--resolution-grid" in cluster_source or "resolution_grid" in cluster_source:
        raise AssertionError("spatial clustering still exposes an automatic resolution scan")

    fixture = work / "fixed_parameters"
    raw_calls = fixture / "raw_calls.tsv"
    write(
        raw_calls,
        "barcode\tx\ty\tgene_id\tgene_name\tumi\tbinomial_p\tumi_t_covered\t"
        "umi_tc_count\tumi_other_mismatch_rate\n"
        "BC1\t10\t10\tG1\tGene1\tU1\t0.001\t20\t2\t0.0\n"
        "BC2\t60\t10\tG2\tGene2\tU2\t1.0\t20\t0\t0.0\n",
    )
    outdir = fixture / "classified"
    stale = outdir / "demo.threshold_sensitivity.tsv"
    write(stale, "stale\n")
    run(
        renderer,
        "scripts/reclassify_slam_calls.py",
        "--molecule-calls",
        str(raw_calls),
        "--sample",
        "demo",
        "--outdir",
        str(outdir),
        "--min-t-covered",
        "10",
        "--min-tc-count",
        "2",
        "--fdr",
        "0.05",
        "--max-other-mismatch-rate",
        "0.1",
    )
    if stale.exists():
        raise AssertionError("legacy SLAM threshold sensitivity output was not removed")
    if not (outdir / "demo.slam_molecule_calls.tsv.gz").is_file():
        raise AssertionError("fixed-parameter SLAM classification did not create molecule calls")


def test_full_length_and_joint_fixtures(work: Path) -> None:
    renderer = plot_python()
    if renderer is None:
        print("SKIP: no configured Python environment with matplotlib and pandas")
        return
    fixture = work / "extended_analysis"
    tmap = fixture / "isoforms.tmap"
    write(
        tmap,
        "class_code\tnum_exons\tlen\n"
        "=\t3\t1200\n"
        "=\t2\t900\n"
        "j\t4\t1500\n",
    )
    transcriptome = fixture / "04_full_length_transcriptome"
    run(
        renderer,
        "scripts/summarize_isoforms.py",
        "--tmap",
        str(tmap),
        "--outdir",
        str(transcriptome),
        "--dpi",
        "72",
    )
    isoform_summary = json.loads(
        (transcriptome / "05_metadata/isoform_class_summary.json").read_text(encoding="utf-8")
    )
    if isoform_summary.get("transcripts") != 3 or isoform_summary.get("class_counts", {}).get("j") != 1:
        raise AssertionError("full-length transcriptome fixture produced an incorrect class summary")

    clusters = fixture / "whole_transcriptome_clusters.tsv"
    write(clusters, "bin_id\tcluster\n25_25\t1\n75_25\t2\n")
    molecules = fixture / "slam_molecules.tsv"
    write(
        molecules,
        "x\ty\tgene_id\tgene_name\texploratory_class\n"
        "10\t10\tG1\tGene1\tnew_supported\n"
        "11\t11\tG1\tGene1\tno_conversion_detected\n"
        "60\t10\tG2\tGene2\tnew_supported\n"
        "61\t11\tG2\tGene2\tno_conversion_detected\n"
        "62\t12\tG2\tGene2\tambiguous\n",
    )
    joint = fixture / "07_slam_joint"
    run(
        renderer,
        "scripts/analyze_spatial_slam_joint.py",
        "--molecule-calls",
        str(molecules),
        "--clusters",
        str(clusters),
        "--sample",
        "demo",
        "--outdir",
        str(joint),
        "--bin-size",
        "50",
        "--min-informative",
        "1",
        "--gene-plot-top",
        "5",
        "--chunk-rows",
        "2",
        "--dpi",
        "72",
    )
    joint_summary = json.loads(
        (joint / "05_metadata/slam_joint_summary.json").read_text(encoding="utf-8")
    )
    if joint_summary.get("molecules_joined_to_whole_transcriptome_bins") != 5:
        raise AssertionError("SLAM joint fixture did not retain all whole-transcriptome-bin molecules")
    expected = (
        joint / "01_tables/cluster_slam_summary.tsv",
        joint / "01_tables/cluster_gene_slam.tsv",
        joint / "02_figures/cluster_slam_composition.png",
        joint / "02_figures/cluster_new_supported_fraction.png",
        joint / "02_figures/cluster_gene_slam_dotplot.png",
    )
    if not all(path.is_file() and path.stat().st_size > 0 for path in expected):
        raise AssertionError("SLAM joint fixture did not create the expected tables and figures")


def test_figure_only_rendering(work: Path) -> None:
    renderer = plot_python()
    if renderer is None:
        print("SKIP: no configured Python environment with matplotlib and pandas")
        return
    fixture = work / "figure_only"
    spectrum = fixture / "mismatch_spectrum.tsv"
    channels = (
        "A>C", "A>G", "A>T", "C>A", "C>G", "C>T",
        "G>A", "G>C", "G>T", "T>A", "T>C", "T>G",
    )
    write(
        spectrum,
        "channel\tcount\treference_coverage\trate\n"
        + "".join(f"{channel}\t{index + 1}\t1000\t{(index + 1) / 1000:.4f}\n" for index, channel in enumerate(channels)),
    )
    mismatch_summary = fixture / "mismatch_summary.json"
    write(mismatch_summary, '{"background_rate": 0.001}\n')
    gene_summary = fixture / "gene_new_supported_summary.tsv"
    write(
        gene_summary,
        "gene_id\tgene_name\tnew_supported\tno_conversion_detected\tambiguous\tinformative_molecules\t"
        "new_supported_fraction\twilson95_low\twilson95_high\tpasses_min_informative\n"
        "ENSMUSG1\tGene1\t20\t80\t2\t100\t0.2\t0.13\t0.29\tTrue\n"
        "ENSMUSG2\tGene2\t10\t90\t1\t100\t0.1\t0.06\t0.17\tTrue\n",
    )
    class_summary = fixture / "slam_summary.json"
    write(
        class_summary,
        '{"class_counts": {"new_supported": 30, "no_conversion_detected": 170, "ambiguous": 3}}\n',
    )
    before = {
        path: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in (spectrum, mismatch_summary, gene_summary, class_summary)
    }
    mismatch_out = fixture / "mismatch_out"
    gene_out = fixture / "gene_out"
    run(
        renderer,
        "scripts/plot_mismatch_spectrum.py",
        "--spectrum",
        str(spectrum),
        "--summary",
        str(mismatch_summary),
        "--outdir",
        str(mismatch_out),
        "--figure-only",
        "--dpi",
        "72",
    )
    run(
        renderer,
        "scripts/summarize_gene_slam.py",
        "--summary-table",
        str(gene_summary),
        "--class-summary",
        str(class_summary),
        "--outdir",
        str(gene_out),
        "--min-informative",
        "10",
        "--dpi",
        "72",
    )
    expected = (
        mismatch_out / "02_figures/base_substitution_comparison.png",
        mismatch_out / "02_figures/base_substitution_comparison.pdf",
        gene_out / "02_figures/gene_new_supported_fraction.png",
        gene_out / "02_figures/gene_new_supported_fraction.pdf",
        gene_out / "02_figures/slam_class_composition.png",
        gene_out / "02_figures/slam_class_composition.pdf",
    )
    if not all(path.is_file() and path.stat().st_size > 0 for path in expected):
        raise AssertionError("figure-only mode did not create the expected PNG/PDF figures")
    if (mismatch_out / "01_tables/base_substitution_comparison.tsv").exists():
        raise AssertionError("mismatch figure-only mode unexpectedly wrote an analysis table")
    if (gene_out / "01_tables/gene_new_supported_summary.tsv").exists():
        raise AssertionError("gene figure-only mode unexpectedly wrote an analysis table")
    after = {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in before}
    if after != before:
        raise AssertionError("figure-only rendering modified a scientific source file")


def test_extended_migration(work: Path) -> None:
    results = work / "migration/results"
    scratch = work / "migration/tmp"
    sample = "demo"
    source = results / sample
    for stage in (
        "01_read_qc",
        "02_barcode_cdna",
        "03_align",
        "04_slam",
        "05_spatial_clustering",
        "06_report",
    ):
        for subdir in ("01_tables", "02_figures", "03_logs", "04_tool_output", "05_metadata"):
            (source / stage / subdir).mkdir(parents=True, exist_ok=True)
    write(source / "01_read_qc/01_tables/qc_metrics.tsv", "metric\tvalue\nn_reads\t1\n")
    write(source / "02_barcode_cdna/05_metadata/extract_summary.json", "{}\n")
    write(source / "03_align/05_metadata/quality_gate.json", '{"overall": "PASS"}\n')
    write(source / f"04_slam/01_tables/{sample}.slam_molecule_calls.tsv.gz", "calls\n")
    write(source / f"04_slam/01_tables/{sample}.threshold_sensitivity.tsv", "threshold\n2\n")
    write(source / f"04_slam/01_tables/{sample}.exploratory_summary.json", "{}\n")
    write(
        source / f"04_slam/05_metadata/{sample}.slam_summary.json",
        json.dumps(
            {
                "interpretation": "exploratory evidence classes only",
                "parameters": {
                    "molecule_calls": str(
                        source
                        / "06_slam_exploratory"
                        / "01_tables"
                        / f"{sample}.exploratory_molecule_calls.tsv.gz"
                    )
                },
            }
        )
        + "\n",
    )
    write(
        source / f"04_slam/04_tool_output/caller/{sample}.qc.json",
        json.dumps(
            {
                "parameters": {
                    "outdir": str(
                        source
                        / "06_slam_exploratory"
                        / "04_tool_output"
                        / "caller"
                    )
                }
            }
        )
        + "\n",
    )
    write(
        source / "05_spatial_clustering/05_metadata/long_read_spatial_clustering.json",
        '{"resolution_grid": [0.8, 1.0], "parameters": {"resolution_grid": [0.8, 1.0], '
        '"resolution_selection": "automatic"}}\n',
    )
    write(source / "05_spatial_clustering/01_tables/long_read_resolution_scan.tsv", "resolution\n1.0\n")
    write(source / "05_spatial_clustering/02_figures/long_read_resolution_scan.png", "png\n")
    write(source / "06_report/report.html", "<html>legacy</html>\n")

    transcriptome = work / "migration/transcriptome_source"
    for subdir in ("01_tables", "02_figures", "03_logs", "04_tool_output", "05_metadata"):
        (transcriptome / subdir).mkdir(parents=True, exist_ok=True)
    write(
        transcriptome / "01_tables/isoform_class_summary.tsv",
        "class_code\ttranscripts\n=\t1\n",
    )
    write(
        transcriptome / "05_metadata/isoform_class_summary.json",
        '{"transcripts": 1, "class_counts": {"=": 1}}\n',
    )

    source_file = source / "01_read_qc/01_tables/qc_metrics.tsv"
    source_inode = source_file.stat().st_ino
    run(
        sys.executable,
        "scripts/00_shared/migrate_extended_results.py",
        "--results-root",
        str(results),
        "--tmp-root",
        str(scratch),
        "--sample-id",
        sample,
        "--transcriptome-source",
        str(transcriptome),
        "--archive-label",
        "test_extended",
        "--apply",
    )
    migrated = results / sample
    archive = scratch / "results_legacy/test_extended" / sample
    if (migrated / "01_read_qc/01_tables/qc_metrics.tsv").stat().st_ino != source_inode:
        raise AssertionError("migration did not preserve large files by hard link")
    if not (archive / "06_report/report.html").is_file():
        raise AssertionError("six-stage source tree was not archived")
    if not (migrated / f"05_slam/01_tables/{sample}.slam_molecule_calls.tsv.gz").is_file():
        raise AssertionError("SLAM output was not moved to the eight-stage location")
    for stale in (
        migrated / f"05_slam/01_tables/{sample}.threshold_sensitivity.tsv",
        migrated / f"05_slam/01_tables/{sample}.exploratory_summary.json",
    ):
        if stale.exists():
            raise AssertionError(f"retired SLAM output survived migration: {stale}")
    if not (migrated / "04_full_length_transcriptome/01_tables/isoform_class_summary.tsv").is_file():
        raise AssertionError("full-length transcriptome source was not imported")
    for stale in (
        migrated / "06_whole_transcriptome/01_tables/long_read_resolution_scan.tsv",
        migrated / "06_whole_transcriptome/02_figures/long_read_resolution_scan.png",
    ):
        if stale.exists():
            raise AssertionError(f"automatic parameter scan output survived migration: {stale}")
    spatial_metadata = json.loads(
        (migrated / "06_whole_transcriptome/05_metadata/long_read_spatial_clustering.json").read_text(
            encoding="utf-8"
        )
    )
    if "resolution_grid" in spatial_metadata or "resolution_grid" in spatial_metadata.get("parameters", {}):
        raise AssertionError("automatic parameter scan metadata survived migration")
    slam_metadata = json.loads(
        (migrated / f"05_slam/05_metadata/{sample}.slam_summary.json").read_text(
            encoding="utf-8"
        )
    )
    if "interpretation" in slam_metadata:
        raise AssertionError("retired exploratory wording survived migration")
    molecule_calls = slam_metadata.get("parameters", {}).get("molecule_calls", "")
    expected_calls = str(
        migrated / "05_slam" / "01_tables" / f"{sample}.slam_molecule_calls.tsv.gz"
    )
    if molecule_calls != expected_calls:
        raise AssertionError(f"internal SLAM path was not repaired: {molecule_calls}")
    caller_qc = json.loads(
        (migrated / f"05_slam/04_tool_output/caller/{sample}.qc.json").read_text(
            encoding="utf-8"
        )
    )
    expected_caller = str(migrated / "05_slam/04_tool_output/caller")
    if caller_qc.get("parameters", {}).get("outdir") != expected_caller:
        raise AssertionError("tool metadata retained an obsolete result path")
    manifests = "\n".join(
        (migrated / stage / "05_metadata/output_manifest.tsv").read_text(encoding="utf-8")
        for stage in EXPECTED_STAGES
    )
    if ".extended_migration" in manifests:
        raise AssertionError("migration manifest contains a staging path")
    if tuple(sorted(path.name for path in migrated.iterdir())) != tuple(sorted(EXPECTED_STAGES)):
        raise AssertionError("migrated result tree does not match the eight-stage contract")


def main() -> int:
    (ROOT / "tmp").mkdir(exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="lr-reusable-smoke-", dir=ROOT / "tmp"))
    try:
        test_config_contract(work)
        test_report(work)
        test_no_automatic_parameter_scans(work)
        test_full_length_and_joint_fixtures(work)
        test_figure_only_rendering(work)
        test_extended_migration(work)
    finally:
        shutil.rmtree(work)
    print("Reusable LR-Celatlas smoke tests passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
