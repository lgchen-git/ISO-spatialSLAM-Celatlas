#!/usr/bin/env python3
"""Screen genes and build real bin50 isoform spatial candidate figures.

The current long-read workflow keeps genomic alignments and barcode/UMI metadata
but does not retain a read-to-isoform table.  This helper reconstructs a
conservative assignment from the primary genomic BAM and the FLAMES transcript
models.  Ambiguous assignments are excluded from candidate selection.
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import json
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle


CONTINUOUS = LinearSegmentedColormap.from_list(
    "npg_magnitude", ["#E5E7EB", "#4DBBD5", "#3C5488", "#E64B35"]
)
CLUSTER_COLORS = {
    "1": "#E64B35", "2": "#4DBBD5", "3": "#00A087", "4": "#3C5488",
    "5": "#F39B7F", "6": "#8491B4", "7": "#91D1C2", "8": "#DC0000",
    "9": "#7E6148", "10": "#B09C85", "11": "#6F99AD", "12": "#8A7AA9",
}


@dataclass
class Transcript:
    transcript_id: str
    gene_id: str
    chrom: str
    strand: str
    start: int
    end: int
    support: int
    source: str
    exons: list[tuple[int, int]] = field(default_factory=list)
    class_code: str = "?"

    @property
    def junctions(self) -> list[tuple[int, int]]:
        exons = sorted(self.exons)
        return [(exons[i][1], exons[i + 1][0]) for i in range(len(exons) - 1)]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-dir", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--per-read", type=Path)
    parser.add_argument("--bin-size", type=int, default=50)
    parser.add_argument("--min-gene-molecules", type=int, default=1000)
    parser.add_argument("--min-gene-clusters", type=int, default=3)
    parser.add_argument("--min-isoform-molecules", type=int, default=100)
    parser.add_argument("--min-dominant-difference", type=float, default=0.25)
    parser.add_argument("--min-js-divergence", type=float, default=0.15)
    parser.add_argument("--max-genes-to-map", type=int, default=50)
    parser.add_argument("--junction-tolerance", type=int, default=10)
    parser.add_argument("--min-read-coverage", type=float, default=0.40)
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args()


def gff_attributes(value: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in value.split(";"):
        if "=" not in item:
            continue
        key, raw = item.split("=", 1)
        result[key] = raw
    return result


def parse_tmap(path: Path) -> dict[str, str]:
    if not path.is_file():
        return {}
    frame = pd.read_csv(path, sep="\t", dtype=str)
    if "qry_id" not in frame or "class_code" not in frame:
        return {}
    return {
        str(row.qry_id).replace("transcript:", ""): str(row.class_code)
        for row in frame[["qry_id", "class_code"]].itertuples(index=False)
    }


def parse_transcripts(gff3: Path, tmap: Path) -> tuple[dict[str, Transcript], dict[str, list[str]]]:
    classes = parse_tmap(tmap)
    transcripts: dict[str, Transcript] = {}
    genes: dict[str, list[str]] = defaultdict(list)
    with gff3.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in {"transcript", "exon"}:
                continue
            attrs = gff_attributes(fields[8])
            if fields[2] == "transcript":
                tid = attrs.get("transcript_id", attrs.get("ID", "")).replace("transcript:", "")
                gene = attrs.get("gene_id", attrs.get("Parent", "").replace("gene:", ""))
                if not tid or not gene:
                    continue
                try:
                    support = int(float(attrs.get("support_count", 0) or 0))
                except ValueError:
                    support = 0
                transcripts[tid] = Transcript(
                    transcript_id=tid,
                    gene_id=gene,
                    chrom=fields[0],
                    strand=fields[6],
                    start=int(fields[3]) - 1,
                    end=int(fields[4]),
                    support=support,
                    source=attrs.get("source", fields[1]),
                    class_code=classes.get(tid, "?"),
                )
                genes[gene].append(tid)
            else:
                parent = attrs.get("Parent", "").replace("transcript:", "")
                if parent in transcripts:
                    transcripts[parent].exons.append((int(fields[3]) - 1, int(fields[4])))
    # FLAMES can emit multiple IDs for the same exon chain.  They are not
    # distinguishable from genomic read alignments, so retain only the
    # highest-support representative before spatial assignment.
    filtered_genes: dict[str, list[str]] = {}
    for gene, tids in genes.items():
        representatives: dict[tuple[str, str, tuple[tuple[int, int], ...]], str] = {}
        for tid in tids:
            model = transcripts[tid]
            key = (model.chrom, model.strand, tuple(sorted(model.exons)))
            previous = representatives.get(key)
            if previous is None or (model.support, model.transcript_id) > (transcripts[previous].support, previous):
                representatives[key] = tid
        filtered_genes[gene] = sorted(representatives.values())
    return transcripts, filtered_genes


def locate_per_read(sample: Path, explicit: Path | None) -> Path:
    if explicit and explicit.is_file():
        return explicit
    candidates = [
        sample / "02_barcode_cdna/04_tool_output/cdna/per_read.tsv.gz",
        sample / "02_barcode_cdna/04_tool_output/per_read.tsv.gz",
    ]
    config = sample / "01_read_qc/05_metadata/input_run_config.yaml"
    if config.is_file():
        match = re.search(r"^tmp_root:\s*(\S+)", config.read_text(encoding="utf-8"), re.MULTILINE)
        if match:
            tmp_root = Path(match.group(1))
            candidates.extend([
                tmp_root / "work" / sample.name / "03_barcode_cdna/cdna/per_read.tsv.gz",
                tmp_root / "work" / sample.name / "02_barcode_cdna/cdna/per_read.tsv.gz",
            ])
    for path in candidates:
        if path.is_file():
            return path
    raise SystemExit("Missing per-read barcode/UMI table. Checked:\n" + "\n".join(str(p) for p in candidates))


def read_bin_metadata(h5ad_path: Path) -> tuple[pd.DataFrame, dict[str, str]]:
    import anndata

    adata = anndata.read_h5ad(h5ad_path, backed="r")
    obs = adata.obs[["bin_x", "bin_y", "cluster"]].copy()
    obs.index = obs.index.astype(str)
    clusters = obs["cluster"].astype(str).to_dict()
    return obs, clusters


def load_gene_screen(
    h5ad_path: Path,
    gene_summary_path: Path,
    transcripts: dict[str, Transcript],
    genes: dict[str, list[str]],
    args: argparse.Namespace,
) -> pd.DataFrame:
    import anndata

    gene_summary = pd.read_csv(gene_summary_path, sep="\t", dtype={"gene_id": str})
    gene_summary["total_molecules"] = pd.to_numeric(gene_summary["total_molecules"], errors="coerce").fillna(0)
    gene_summary["detected_bins"] = pd.to_numeric(gene_summary["detected_bins"], errors="coerce").fillna(0)
    multi = {
        gene: [tid for tid in tids if transcripts[tid].support >= 5]
        for gene, tids in genes.items()
    }
    strong = {
        gene: [tid for tid in tids if transcripts[tid].support >= 20]
        for gene, tids in genes.items()
    }
    base = gene_summary[gene_summary["total_molecules"] >= args.min_gene_molecules].copy()
    base = base[base["gene_id"].isin(multi)]
    base = base[base["gene_id"].map(lambda gene: len(multi.get(gene, [])) >= 2)]
    base = base[base["gene_id"].map(lambda gene: len(strong.get(gene, [])) >= 2)]
    if base.empty:
        return pd.DataFrame()
    base = base.sort_values("total_molecules", ascending=False).head(max(args.max_genes_to_map * 4, 100))

    adata = anndata.read_h5ad(h5ad_path, backed="r")
    var_ids = [str(x) for x in adata.var_names]
    var_index = {gene: i for i, gene in enumerate(var_ids)}
    selected = [gene for gene in base["gene_id"] if gene in var_index]
    counts = adata.layers["counts"]
    cluster_values = adata.obs["cluster"].astype(str).to_numpy()
    cluster_levels = sorted(set(cluster_values), key=lambda value: int(value) if value.isdigit() else value)
    var_names = dict(zip(
        [str(value) for value in adata.var_names],
        adata.var["gene_name"].astype(str).tolist(),
    )) if "gene_name" in adata.var else {}
    rows: list[dict[str, object]] = []
    for gene in selected:
        values = counts[:, var_index[gene]]
        if hasattr(values, "toarray"):
            values = values.toarray()
        values = np.asarray(values).reshape(-1)
        cluster_totals = {cluster: float(values[cluster_values == cluster].sum()) for cluster in cluster_levels}
        cluster_bins = {cluster: int((values[cluster_values == cluster] > 0).sum()) for cluster in cluster_levels}
        expressed = [cluster for cluster, value in cluster_bins.items() if value > 0]
        top_cluster = max(cluster_totals, key=cluster_totals.get) if cluster_totals else "NA"
        rows.append({
            "gene_id": gene,
            "gene_name": var_names.get(gene, gene),
            "total_molecules": int(values.sum()),
            "detected_bins": int((values > 0).sum()),
            "expressed_clusters": len(expressed),
            "top_cluster": top_cluster,
            "isoform_models": len(multi.get(gene, [])),
            "strong_isoforms": len(strong.get(gene, [])),
            "cluster_totals": json.dumps(cluster_totals, sort_keys=True),
            "cluster_bins": json.dumps(cluster_bins, sort_keys=True),
            "screen_status": "eligible" if len(expressed) >= args.min_gene_clusters else "too_few_clusters",
        })
    adata.file.close()
    result = pd.DataFrame(rows)
    result = result[result["expressed_clusters"] >= args.min_gene_clusters]
    result["screen_score"] = np.log1p(result["total_molecules"]) * np.log1p(result["detected_bins"])
    result["mapping_score"] = result["screen_score"] * np.log1p(result["strong_isoforms"])
    return result.sort_values("mapping_score", ascending=False).reset_index(drop=True)


def make_interval_index(transcripts: dict[str, Transcript], gene_ids: set[str], width: int = 100_000) -> dict[tuple[str, int], set[str]]:
    index: dict[tuple[str, int], set[str]] = defaultdict(set)
    for transcript in transcripts.values():
        if transcript.gene_id not in gene_ids:
            continue
        for bucket in range(transcript.start // width, transcript.end // width + 1):
            index[(transcript.chrom, bucket)].add(transcript.gene_id)
    return index


def candidate_genes_for_read(
    chrom: str,
    start: int,
    end: int,
    interval_index: dict[tuple[str, int], set[str]],
    gene_bounds: dict[str, tuple[str, int, int]],
    width: int = 100_000,
) -> list[str]:
    found: set[str] = set()
    for bucket in range(start // width, end // width + 1):
        for gene in interval_index.get((chrom, bucket), set()):
            _, gene_start, gene_end = gene_bounds[gene]
            if gene_start < end and start < gene_end:
                found.add(gene)
    return sorted(found)


def read_junctions(read: pysam.AlignedSegment) -> list[tuple[int, int]]:
    if not read.cigartuples:
        return []
    junctions: list[tuple[int, int]] = []
    reference = read.reference_start
    for operation, length in read.cigartuples:
        if operation == 3:
            junctions.append((reference, reference + length))
            reference += length
        elif operation in {0, 2, 7, 8}:
            reference += length
    return junctions


def overlap_length(blocks: list[tuple[int, int]], exons: list[tuple[int, int]]) -> int:
    total = 0
    for block_start, block_end in blocks:
        for exon_start, exon_end in exons:
            total += max(0, min(block_end, exon_end) - max(block_start, exon_start))
    return total


def transcript_score(read: pysam.AlignedSegment, transcript: Transcript, tolerance: int, min_coverage: float, read_junc: list[tuple[int, int]] | None = None) -> tuple[float, float] | None:
    blocks = read.get_blocks()
    aligned_bases = sum(end - start for start, end in blocks)
    if not blocks or aligned_bases == 0:
        return None
    if read.reference_name != transcript.chrom or transcript.start >= read.reference_end or read.reference_start >= transcript.end:
        return None
    coverage = overlap_length(blocks, transcript.exons) / aligned_bases
    if coverage < min_coverage:
        return None
    read_junc = read_junc if read_junc is not None else read_junctions(read)
    iso_junc = transcript.junctions
    matched = 0
    for start, end in read_junc:
        if any(abs(start - candidate_start) <= tolerance and abs(end - candidate_end) <= tolerance for candidate_start, candidate_end in iso_junc):
            matched += 1
        else:
            return None
    if read_junc and matched == 0:
        return None
    terminal_distance = abs(read.reference_start - transcript.start) + abs(read.reference_end - transcript.end)
    terminal_penalty = min(4.0, terminal_distance / 1000.0)
    score = matched * 20.0 + coverage * 10.0 - terminal_penalty
    return score, coverage


def load_read_metadata(per_read: Path, wanted: set[str]) -> dict[str, tuple[str, str, int, int]]:
    metadata: dict[str, tuple[str, str, int, int]] = {}
    usecols = ["read_id", "assign_status", "barcode_corr", "umi_seq", "x", "y"]
    for chunk in pd.read_csv(per_read, sep="\t", usecols=usecols, dtype=str, compression="gzip", chunksize=500_000):
        chunk = chunk[chunk["assign_status"].isin({"exact", "rescued"})]
        chunk = chunk[chunk["read_id"].isin(wanted)]
        chunk = chunk[(chunk["barcode_corr"].notna()) & (chunk["umi_seq"].notna())]
        chunk = chunk[(chunk["barcode_corr"] != ".") & (chunk["umi_seq"] != ".")]
        if chunk.empty:
            continue
        chunk["x"] = pd.to_numeric(chunk["x"], errors="coerce")
        chunk["y"] = pd.to_numeric(chunk["y"], errors="coerce")
        chunk = chunk.dropna(subset=["x", "y"])
        metadata.update({
            str(row.read_id): (str(row.barcode_corr), str(row.umi_seq), int(row.x), int(row.y))
            for row in chunk.itertuples(index=False)
        })
    return metadata


def assign_reads(
    bam_path: Path,
    per_read: Path,
    candidate_gene_ids: set[str],
    transcripts: dict[str, Transcript],
    gene_transcripts: dict[str, list[str]],
    bin_clusters: dict[str, str],
    bin_size: int,
    tolerance: int,
    min_coverage: float,
) -> tuple[list[dict[str, object]], dict[str, int]]:
    gene_bounds: dict[str, tuple[str, int, int]] = {}
    for gene in candidate_gene_ids:
        models = [transcripts[tid] for tid in gene_transcripts.get(gene, [])]
        if not models:
            continue
        gene_bounds[gene] = (models[0].chrom, min(x.start for x in models), max(x.end for x in models))
    interval_index = make_interval_index(transcripts, set(gene_bounds))
    junction_index: dict[str, dict[tuple[int, int], set[str]]] = defaultdict(lambda: defaultdict(set))
    for gene, tids in gene_transcripts.items():
        if gene not in gene_bounds:
            continue
        for tid in tids:
            for junction in transcripts[tid].junctions:
                junction_index[gene][junction].add(tid)
    regions = sorted(set(gene_bounds.values()))
    print(f"[INFO] candidate_gene_regions={len(regions)}", flush=True)
    wanted: set[str] = set()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for chrom, start, end in regions:
            for read in bam.fetch(chrom, start, end):
                if read.is_unmapped or read.reference_name is None:
                    continue
                if candidate_genes_for_read(read.reference_name, read.reference_start, read.reference_end, interval_index, gene_bounds):
                    wanted.add(read.query_name)
    print(f"[INFO] candidate_read_ids={len(wanted)}", flush=True)
    metadata = load_read_metadata(per_read, wanted)
    print(f"[INFO] reads_with_spatial_metadata={len(metadata)}", flush=True)
    molecules: dict[tuple[str, str, str, int, int], dict[str, object]] = {}
    qc = Counter({"bam_reads_over_candidate_genes": 0, "reads_with_spatial_metadata": 0, "assigned_reads": 0, "ambiguous_reads": 0})
    processed_reads: set[str] = set()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for region_index, (chrom, start, end) in enumerate(regions, start=1):
            for read in bam.fetch(chrom, start, end):
                if read.is_unmapped or read.query_name not in metadata or read.reference_name is None or read.query_name in processed_reads:
                    continue
                processed_reads.add(read.query_name)
                gene_ids = candidate_genes_for_read(read.reference_name, read.reference_start, read.reference_end, interval_index, gene_bounds)
                if not gene_ids:
                    continue
                qc["bam_reads_over_candidate_genes"] += 1
                barcode, umi, x, y = metadata[read.query_name]
                qc["reads_with_spatial_metadata"] += 1
                read_junc = read_junctions(read)
                if not read_junc:
                    qc["reads_without_splice_junction"] += 1
                    continue
                scored: list[tuple[float, float, str, str]] = []
                for gene in gene_ids:
                    possible: set[str] | None = None
                    for read_start, read_end in read_junc:
                        matching: set[str] = set()
                        for (junction_start, junction_end), tids in junction_index.get(gene, {}).items():
                            if abs(read_start - junction_start) <= tolerance and abs(read_end - junction_end) <= tolerance:
                                matching.update(tids)
                        possible = matching if possible is None else possible.intersection(matching)
                    for tid in sorted(possible or set()):
                        result = transcript_score(read, transcripts[tid], tolerance, min_coverage, read_junc)
                        if result is not None:
                            score, coverage = result
                            scored.append((score, coverage, gene, tid))
                if not scored:
                    continue
                scored.sort(reverse=True)
                best = scored[0]
                tied = [row for row in scored if best[0] - row[0] < 2.0]
                if len(tied) > 1:
                    qc["ambiguous_reads"] += 1
                    continue
                qc["assigned_reads"] += 1
                gene, tid = best[2], best[3]
                bid = f"{(x // bin_size) * bin_size + bin_size // 2}_{(y // bin_size) * bin_size + bin_size // 2}"
                key = (barcode, gene, umi, x, y)
                row = molecules.setdefault(key, {"barcode": barcode, "gene_id": gene, "umi": umi, "x": x, "y": y, "bin_id": bid, "cluster": bin_clusters.get(bid, "NA"), "transcript_id": tid, "read_count": 0, "score": best[0], "coverage": best[1], "ambiguous": False})
                row["read_count"] = int(row["read_count"]) + 1
                if row["transcript_id"] != tid:
                    row["ambiguous"] = True
            print(f"[INFO] assigned_region={region_index}/{len(regions)}", flush=True)
    rows = [row for row in molecules.values() if not row["ambiguous"]]
    qc["unique_molecules"] = len(molecules)
    qc["unique_unambiguous_molecules"] = len(rows)
    return rows, dict(qc)


def js_divergence(a: Counter[str], b: Counter[str], clusters: list[str]) -> float:
    p = np.array([a.get(cluster, 0) for cluster in clusters], dtype=float)
    q = np.array([b.get(cluster, 0) for cluster in clusters], dtype=float)
    if p.sum() == 0 or q.sum() == 0:
        return 0.0
    p, q = p / p.sum(), q / q.sum()
    m = (p + q) / 2
    p_mask = p > 0
    q_mask = q > 0
    p_kl = float(np.sum(p[p_mask] * np.log2(p[p_mask] / m[p_mask])))
    q_kl = float(np.sum(q[q_mask] * np.log2(q[q_mask] / m[q_mask])))
    return 0.5 * (p_kl + q_kl)


def choose_isoform_pair(rows: list[dict[str, object]], screen: pd.DataFrame, transcripts: dict[str, Transcript], args: argparse.Namespace) -> tuple[pd.DataFrame, dict[str, object] | None]:
    clusters = sorted({str(row["cluster"]) for row in rows if str(row["cluster"]) != "NA"}, key=lambda value: int(value) if value.isdigit() else value)
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    totals: Counter[str] = Counter()
    genes: dict[str, str] = {}
    for row in rows:
        tid = str(row["transcript_id"])
        gene = str(row["gene_id"])
        counts[tid][str(row["cluster"])] += 1
        totals[tid] += 1
        genes[tid] = gene
    output: list[dict[str, object]] = []
    best: dict[str, object] | None = None
    for gene in screen["gene_id"].astype(str):
        tids = [tid for tid in counts if genes[tid] == gene and totals[tid] >= args.min_isoform_molecules]
        for i, a in enumerate(tids):
            for b in tids[i + 1:]:
                ca, cb = counts[a], counts[b]
                da, db = ca.most_common(1)[0], cb.most_common(1)[0]
                separation = js_divergence(ca, cb, clusters)
                dominant_diff = abs(da[1] / totals[a] - db[1] / totals[b])
                if da[0] == db[0] or dominant_diff < args.min_dominant_difference or separation < args.min_js_divergence:
                    continue
                enrichment = math.log2((ca[da[0]] + 1) / (cb[da[0]] + 1))
                score = separation * math.log1p(totals[a] + totals[b]) * (1 + dominant_diff)
                record = {
                    "gene_id": gene,
                    "isoform_a": a,
                    "isoform_b": b,
                    "isoform_a_molecules": totals[a],
                    "isoform_b_molecules": totals[b],
                    "isoform_a_top_cluster": da[0],
                    "isoform_b_top_cluster": db[0],
                    "isoform_a_top_fraction": da[1] / totals[a],
                    "isoform_b_top_fraction": db[1] / totals[b],
                    "dominant_fraction_difference": dominant_diff,
                    "js_divergence": separation,
                    "log2_a_over_b_at_a_cluster": enrichment,
                    "candidate_score": score,
                    "class_a": transcripts[a].class_code,
                    "class_b": transcripts[b].class_code,
                }
                output.append(record)
                if best is None or score > float(best["candidate_score"]):
                    best = record
    return pd.DataFrame(output).sort_values("candidate_score", ascending=False) if output else pd.DataFrame(), best


def event_candidates(gene: str, tids: list[str], transcripts: dict[str, Transcript], molecules: Counter[str]) -> list[dict[str, object]]:
    events: list[dict[str, object]] = []
    usable = [transcripts[tid] for tid in tids if molecules[tid] >= 50]
    for i, a in enumerate(usable):
        for b in usable[i + 1:]:
            ja, jb = set(a.junctions), set(b.junctions)
            if ja != jb and (ja.issubset(jb) or jb.issubset(ja)):
                long, short = (a, b) if len(a.exons) > len(b.exons) else (b, a)
                events.append({"event_type": "exon_skipping", "gene_id": gene, "isoform_a": long.transcript_id, "isoform_b": short.transcript_id, "support_a": molecules[long.transcript_id], "support_b": molecules[short.transcript_id]})
            if a.junctions == b.junctions:
                end_delta = abs(a.end - b.end)
                start_delta = abs(a.start - b.start)
                if end_delta >= 50:
                    events.append({"event_type": "APA", "gene_id": gene, "isoform_a": a.transcript_id, "isoform_b": b.transcript_id, "support_a": molecules[a.transcript_id], "support_b": molecules[b.transcript_id], "coordinate_delta": end_delta})
                if start_delta >= 50:
                    events.append({"event_type": "TSS", "gene_id": gene, "isoform_a": a.transcript_id, "isoform_b": b.transcript_id, "support_a": molecules[a.transcript_id], "support_b": molecules[b.transcript_id], "coordinate_delta": start_delta})
    for transcript in usable:
        if transcript.class_code in {"j", "o", "k", "m", "n"} or transcript.source == "new":
            events.append({"event_type": "novel_transcript", "gene_id": gene, "isoform_a": transcript.transcript_id, "isoform_b": "", "support_a": molecules[transcript.transcript_id], "support_b": 0, "coordinate_delta": 0})
    return events


def spatial_axes(ax, obs: pd.DataFrame, title: str) -> None:
    ax.set_title(title, loc="left", fontsize=11, fontweight="bold")
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("y (µm)")
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xlim(-100, max(9457, float(obs["bin_x"].max()) + 100))
    ax.set_ylim(max(9457, float(obs["bin_y"].max()) + 100), -100)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)


def save_figure(fig, stem: Path, dpi: int) -> None:
    fig.savefig(stem.with_suffix(".png"), dpi=dpi, facecolor="white", bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), dpi=min(dpi, 150), facecolor="white", bbox_inches="tight")
    plt.close(fig)


def write_figure_manifest(figdir: Path, outdir: Path, sample_id: str) -> None:
    rows = []
    for path in sorted(figdir.glob("11_*.png")):
        rows.append({"sample_id": sample_id, "figure_id": path.stem, "png": str(path.resolve()), "pdf": str(path.with_suffix(".pdf").resolve())})
    pd.DataFrame(rows, columns=["sample_id", "figure_id", "png", "pdf"]).to_csv(outdir / "isoform_figure_manifest.tsv", sep="\t", index=False)


def draw_spatial(stem: Path, obs: pd.DataFrame, values: pd.Series, title: str, dpi: int, cmap=CONTINUOUS, shared_limits: tuple[float, float] | None = None) -> None:
    fig, ax = plt.subplots(figsize=(6.5, 5.8))
    values = pd.to_numeric(values, errors="coerce").fillna(0)
    low, high = shared_limits or (float(values.quantile(0.02)), float(values.quantile(0.98)))
    if high <= low:
        low, high = float(values.min()), float(values.max() or 1)
    if high <= low:
        high = low + 1
    sc = ax.scatter(obs["bin_x"], obs["bin_y"], c=values, s=8, marker="s", cmap=cmap, vmin=low, vmax=high, linewidths=0, rasterized=True)
    spatial_axes(ax, obs, title)
    cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Molecule count", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    fig.tight_layout()
    save_figure(fig, stem, dpi)


def draw_structure(stem: Path, gene: str, event: dict[str, object], transcripts: dict[str, Transcript], dpi: int) -> None:
    ids = [str(event["isoform_a"])] + ([str(event["isoform_b"])] if event.get("isoform_b") else [])
    models = [transcripts[tid] for tid in ids if tid in transcripts]
    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    event_label = {"exon_skipping": "Exon skipping", "APA": "APA", "TSS": "TSS", "novel_transcript": "Novel transcript"}.get(str(event["event_type"]), str(event["event_type"]))
    ax.set_title(f"{gene} · {event_label}", loc="left", fontsize=11, fontweight="bold")
    coordinates = [coord for model in models for exon in model.exons for coord in exon]
    if not coordinates:
        ax.axis("off")
        ax.text(0.5, 0.5, "Structure unavailable", ha="center", va="center")
    else:
        left, right = min(coordinates) - 100, max(coordinates) + 100
        for index, model in enumerate(models):
            y = len(models) - index
            exons = sorted(model.exons)
            ax.hlines(y, exons[0][0], exons[-1][1], color="#9CA3AF", linewidth=1)
            for start, end in exons:
                ax.add_patch(Rectangle((start, y - 0.18), end - start, 0.36, color="#3C5488" if index == 0 else "#E64B35", alpha=0.9))
            ax.text(left, y + 0.26, f"{model.transcript_id}  n={event['support_a'] if index == 0 else event.get('support_b', 0):,}", fontsize=8, va="bottom")
        ax.set_xlim(left, right)
        ax.set_ylim(0.3, len(models) + 0.8)
        ax.set_yticks([])
        ax.set_xlabel("Reference genomic coordinate")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(False)
    fig.tight_layout()
    save_figure(fig, stem, dpi)


def main() -> None:
    args = parse_args()
    sample = args.sample_dir
    outdir = args.outdir
    table_dir = outdir / "01_tables"
    figdir = outdir / "02_figures"
    metadata_dir = outdir / "05_metadata"
    for path in (table_dir, figdir, metadata_dir):
        path.mkdir(parents=True, exist_ok=True)

    s3 = sample / "03_align"
    s4 = sample / "04_full_length_transcriptome"
    s5 = sample / "05_slam"
    s6 = sample / "06_whole_transcriptome"
    bam = s3 / "04_tool_output/minimap2" / f"{args.sample_id}.primary.q20.bam"
    gff3 = s4 / "04_tool_output/flames/isoforms.gff3"
    tmap = s4 / "04_tool_output/flames/flames_vs_reference.isoforms.gff3.tmap"
    h5ad = s6 / "04_tool_output/scanpy/long_read_spatial_bin50.h5ad"
    gene_summary = s6 / "01_tables/whole_transcriptome_gene_summary.tsv"
    for path in (bam, gff3, tmap, h5ad, gene_summary):
        if not path.is_file():
            raise SystemExit(f"Missing required input: {path}")
    per_read = locate_per_read(sample, args.per_read)

    transcripts, gene_transcripts = parse_transcripts(gff3, tmap)
    obs, bin_clusters = read_bin_metadata(h5ad)
    screen = load_gene_screen(h5ad, gene_summary, transcripts, gene_transcripts, args)
    screen.to_csv(table_dir / "candidate_gene_screen.tsv", sep="\t", index=False)
    if screen.empty:
        (metadata_dir / "isoform_spatial_qc.json").write_text(json.dumps({"status": "no_gene_candidates"}, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({"status": "no_gene_candidates"}))
        return
    map_genes = set(screen.head(args.max_genes_to_map)["gene_id"].astype(str))
    rows, qc = assign_reads(bam, per_read, map_genes, transcripts, gene_transcripts, bin_clusters, args.bin_size, args.junction_tolerance, args.min_read_coverage)
    pd.DataFrame(rows, columns=["barcode", "gene_id", "umi", "x", "y", "bin_id", "cluster", "transcript_id", "read_count", "score", "coverage", "ambiguous"]).to_csv(table_dir / "isoform_spatial_molecules.tsv.gz", sep="\t", index=False, compression="gzip")
    pair_table, best_pair = choose_isoform_pair(rows, screen, transcripts, args)
    pair_table.to_csv(table_dir / "candidate_isoform_pairs.tsv", sep="\t", index=False)
    if best_pair is None:
        qc["status"] = "no_qualifying_isoform_pair"
        qc.update({"min_isoform_molecules": args.min_isoform_molecules, "min_dominant_difference": args.min_dominant_difference, "min_js_divergence": args.min_js_divergence})
        (metadata_dir / "isoform_spatial_qc.json").write_text(json.dumps(qc, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(qc))
        return

    gene = str(best_pair["gene_id"])
    gene_name = screen.loc[screen["gene_id"].astype(str).eq(gene), "gene_name"].iloc[0]
    selected = {"gene_id": gene, "gene_name": str(gene_name), "interpretation": "exploratory_isoform_spatial_proxy", **best_pair}
    (table_dir / "selected_isoform_candidate.json").write_text(json.dumps(selected, indent=2) + "\n", encoding="utf-8")
    selected_rows = [row for row in rows if row["gene_id"] == gene]
    molecule_counts = Counter(str(row["transcript_id"]) for row in selected_rows)
    event_rows = event_candidates(gene, gene_transcripts.get(gene, []), transcripts, molecule_counts)
    event_table = pd.DataFrame(event_rows)
    event_table.to_csv(table_dir / "candidate_structure_events.tsv", sep="\t", index=False)

    import anndata
    adata = anndata.read_h5ad(h5ad, backed="r")
    var_index = {str(value): index for index, value in enumerate(adata.var_names)}
    gene_values = adata.layers["counts"][:, var_index[gene]]
    if hasattr(gene_values, "toarray"):
        gene_values = gene_values.toarray()
    gene_values = pd.Series(np.asarray(gene_values).reshape(-1), index=obs.index)
    draw_spatial(figdir / "11_candidate_gene_total_spatial", obs, gene_values, f"{gene_name} spatial expression", args.dpi)

    iso_a, iso_b = str(best_pair["isoform_a"]), str(best_pair["isoform_b"])
    iso_counts = {}
    for tid in (iso_a, iso_b):
        count = Counter(str(row["bin_id"]) for row in selected_rows if str(row["transcript_id"]) == tid)
        iso_counts[tid] = pd.Series([count.get(str(index), 0) for index in obs.index], index=obs.index)
    shared_high = max(float(pd.concat(iso_counts.values()).quantile(0.98)), 1.0)
    for label, tid in (("A", iso_a), ("B", iso_b)):
        draw_spatial(figdir / f"11_candidate_isoform_{label}_spatial", obs, iso_counts[tid], f"{gene_name} · Isoform {label} spatial distribution", args.dpi, shared_limits=(0.0, shared_high))
    usage_rows = []
    for cluster in sorted(set(obs["cluster"].astype(str)), key=lambda value: int(value) if value.isdigit() else value):
        usage_rows.append({"cluster": cluster, "Isoform A": sum(1 for row in selected_rows if str(row["transcript_id"]) == iso_a and str(row["cluster"]) == cluster), "Isoform B": sum(1 for row in selected_rows if str(row["transcript_id"]) == iso_b and str(row["cluster"]) == cluster)})
    usage = pd.DataFrame(usage_rows)
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    x = np.arange(len(usage))
    width = 0.38
    ax.bar(x - width / 2, usage["Isoform A"], width, label="Isoform A", color="#3C5488")
    ax.bar(x + width / 2, usage["Isoform B"], width, label="Isoform B", color="#E64B35")
    ax.set_xticks(x, usage["cluster"].astype(str))
    ax.set_xlabel("Spatial cluster")
    ax.set_ylabel("UMI-level molecules")
    ax.set_title(f"{gene_name} · Isoform usage by spatial cluster", loc="left", fontsize=11, fontweight="bold")
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    save_figure(fig, figdir / "11_candidate_isoform_cluster_usage", args.dpi)
    usage.to_csv(table_dir / "selected_isoform_cluster_usage.tsv", sep="\t", index=False)

    for event_type in ("exon_skipping", "APA", "TSS", "novel_transcript"):
        matching = event_table[event_table["event_type"].eq(event_type)] if not event_table.empty else pd.DataFrame()
        if matching.empty:
            continue
        event = matching.sort_values(["support_a", "support_b"], ascending=False).iloc[0].to_dict()
        draw_structure(figdir / f"11_candidate_event_{event_type.lower()}", str(gene_name), event, transcripts, args.dpi)

    qc.update({"status": "selected", "interpretation": "exploratory_isoform_spatial_proxy", "selected_gene_id": gene, "selected_gene_name": str(gene_name), "selected_isoform_a": iso_a, "selected_isoform_b": iso_b, "figure_prefix": "11_"})
    qc.update({"min_isoform_molecules": args.min_isoform_molecules, "min_dominant_difference": args.min_dominant_difference, "min_js_divergence": args.min_js_divergence})
    (metadata_dir / "isoform_spatial_qc.json").write_text(json.dumps(qc, indent=2) + "\n", encoding="utf-8")
    write_figure_manifest(figdir, outdir, args.sample_id)
    print(json.dumps(qc, ensure_ascii=False))


if __name__ == "__main__":
    main()
