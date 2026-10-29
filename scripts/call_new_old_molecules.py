#!/usr/bin/env python3
import argparse
import csv
import gzip
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
import pysam

try:
    from scipy.stats import binom as scipy_binom
except ImportError:  # pragma: no cover
    scipy_binom = None

BIN_SIZE = 20000
BASE_EVENTS = ("match_t_count", "tc_count", "ta_count", "tg_count", "deletion_count")
SITE_COLUMNS = [
    "chrom", "pos", "ref_base", "t_depth", "match_t_count", "tc_count",
    "ta_count", "tg_count", "deletion_count", "tc_rate", "ta_rate",
    "tg_rate", "other_mismatch_rate", "site_call",
]
MOLECULE_COLUMNS = [
    "barcode", "x", "y", "gene_id", "gene_name", "umi", "supporting_reads",
    "tc_events", "t_covered", "feature_class", "assignment", "final_call",
    "call_confidence", "ambiguity_reason", "umi_t_covered", "umi_tc_count",
    "umi_ta_count", "umi_tg_count", "umi_tc_rate", "umi_other_mismatch_rate",
    "high_conf_tc_sites_hit", "noisy_sites_hit", "binomial_p", "bg_tc_rate",
    "bg_rate_source", "read_sum_tc_count", "read_sum_t_covered",
]


def load_annotations(annotation_dir):
    exon_df = pd.read_csv(Path(annotation_dir) / "exons.tsv", sep="\t", dtype={"chrom": str})
    intron_df = pd.read_csv(Path(annotation_dir) / "introns.tsv", sep="\t", dtype={"chrom": str})
    annotations = defaultdict(lambda: defaultdict(list))
    for frame in (exon_df, intron_df):
        for row in frame.itertuples(index=False):
            interval = {
                "start": int(row.start),
                "end": int(row.end),
                "gene_id": row.gene_id,
                "gene_name": row.gene_name,
                "strand": row.strand,
                "feature_class": row.feature_class,
            }
            for bin_id in range(interval["start"] // BIN_SIZE, (interval["end"] - 1) // BIN_SIZE + 1):
                annotations[str(row.chrom)][bin_id].append(interval)
    return annotations


def load_barcode_table(path):
    table = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        header = next(handle)
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            table[fields[0]] = (fields[1], fields[2], int(fields[3]), int(fields[4]))
    return table


def load_excluded_sites(path):
    sites = set()
    if not path:
        return sites
    variants = pysam.VariantFile(path)
    for record in variants:
        sites.add((record.contig, record.pos - 1))
    variants.close()
    return sites


def overlap_fraction(read_blocks, interval):
    covered = 0
    for start, end in read_blocks:
        ov_start = max(start, interval["start"])
        ov_end = min(end, interval["end"])
        if ov_end > ov_start:
            covered += ov_end - ov_start
    return covered


def assign_gene_feature(chrom, read_blocks, annotations):
    overlaps = defaultdict(lambda: {"score": 0, "features": set(), "gene_name": "", "strand": ""})
    candidates = {}
    chrom_bins = annotations.get(str(chrom))
    if not chrom_bins:
        return None
    for start, end in read_blocks:
        for bin_id in range(start // BIN_SIZE, (end - 1) // BIN_SIZE + 1):
            for interval in chrom_bins.get(bin_id, []):
                candidates[id(interval)] = interval
    for interval in candidates.values():
        score = overlap_fraction(read_blocks, interval)
        if score <= 0:
            continue
        entry = overlaps[interval["gene_id"]]
        entry["score"] += score
        entry["features"].add(interval["feature_class"])
        entry["gene_name"] = interval["gene_name"]
        entry["strand"] = interval["strand"]
    if not overlaps:
        return None
    sorted_hits = sorted(overlaps.items(), key=lambda item: item[1]["score"], reverse=True)
    top_gene, top = sorted_hits[0]
    if len(sorted_hits) > 1 and sorted_hits[1][1]["score"] == top["score"]:
        return {"gene_id": None, "feature_class": "ambiguous_gene"}
    feature_class = "exon" if "exon" in top["features"] else "intron"
    if len(top["features"]) > 1:
        feature_class = "ambiguous_feature"
    return {
        "gene_id": top_gene,
        "gene_name": top["gene_name"],
        "strand": top["strand"],
        "feature_class": feature_class,
    }


def block_reference_cache(read, fasta):
    cache = []
    for start, end in read.get_blocks():
        cache.append((start, end, fasta.fetch(read.reference_name, start, end).upper()))
    return cache


def empty_base_counter():
    return Counter({key: 0 for key in BASE_EVENTS})


def is_transcript_t(ref_base, strand):
    return (strand == "+" and ref_base == "T") or (strand == "-" and ref_base == "A")


def classify_t_base(ref_base, query_base, strand):
    ref_base = ref_base.upper()
    query_base = query_base.upper()
    if strand == "+" and ref_base == "T":
        return {"T": "match_t_count", "C": "tc_count", "A": "ta_count", "G": "tg_count"}.get(query_base)
    if strand == "-" and ref_base == "A":
        return {"A": "match_t_count", "G": "tc_count", "T": "ta_count", "C": "tg_count"}.get(query_base)
    return None


def check_base_qualities(read):
    if read.query_qualities is None:
        raise SystemExit(
            f"Read {read.query_name} has no query base qualities. "
            "Formal SLAM calling requires base qualities for --min-tc-baseq filtering."
        )


def ref_base_from_blocks(block_cache, rpos):
    for start, end, ref_seq in block_cache:
        if start <= rpos < end:
            return ref_seq[rpos - start]
    return ""


def iter_t_events_from_cigar_md(read, strand, min_baseq):
    query_sequence = (read.query_sequence or "").upper()
    qualities = read.query_qualities
    ref_sequence = read.get_reference_sequence().upper()
    qpos = 0
    rpos = read.reference_start
    ref_idx = 0
    for op, length in read.cigartuples or []:
        if op in {0, 7, 8}:
            for offset in range(length):
                ref_base = ref_sequence[ref_idx + offset]
                if is_transcript_t(ref_base, strand) and qualities[qpos + offset] >= min_baseq:
                    event = classify_t_base(ref_base, query_sequence[qpos + offset], strand)
                    if event:
                        yield read.reference_name, rpos + offset, event
            qpos += length
            rpos += length
            ref_idx += length
        elif op == 2:
            for offset in range(length):
                ref_base = ref_sequence[ref_idx + offset]
                if is_transcript_t(ref_base, strand):
                    yield read.reference_name, rpos + offset, "deletion_count"
            rpos += length
            ref_idx += length
        elif op == 3:
            rpos += length
        elif op in {1, 4}:
            qpos += length
        elif op in {5, 6}:
            continue


def iter_t_events(read, fasta, strand, min_baseq):
    check_base_qualities(read)
    if read.has_tag("MD"):
        yield from iter_t_events_from_cigar_md(read, strand, min_baseq)
        return
    query_sequence = (read.query_sequence or "").upper()
    qualities = read.query_qualities
    block_cache = block_reference_cache(read, fasta)
    for qpos, rpos in read.get_aligned_pairs(matches_only=False):
        if rpos is None:
            continue
        ref_base = ref_base_from_blocks(block_cache, rpos)
        if not ref_base or not is_transcript_t(ref_base, strand):
            continue
        if qpos is None:
            yield read.reference_name, rpos, "deletion_count"
            continue
        if qualities[qpos] < min_baseq:
            continue
        event = classify_t_base(ref_base, query_sequence[qpos], strand)
        if event:
            yield read.reference_name, rpos, event


def consensus_event(events):
    observed = {key: events[key] for key in ("match_t_count", "tc_count", "ta_count", "tg_count") if events[key] > 0}
    if not observed:
        return "deletion_count"
    ranked = sorted(observed.items(), key=lambda item: item[1], reverse=True)
    if len(ranked) > 1 and ranked[0][1] == ranked[1][1]:
        if ranked[0][0] == "match_t_count":
            return "match_t_count"
        return "tg_count"
    return ranked[0][0]


def site_rates(counts):
    t_depth = counts["match_t_count"] + counts["tc_count"] + counts["ta_count"] + counts["tg_count"]
    return {
        "t_depth": t_depth,
        "tc_rate": counts["tc_count"] / t_depth if t_depth else 0.0,
        "ta_rate": counts["ta_count"] / t_depth if t_depth else 0.0,
        "tg_rate": counts["tg_count"] / t_depth if t_depth else 0.0,
        "other_mismatch_rate": (counts["ta_count"] + counts["tg_count"]) / t_depth if t_depth else 0.0,
    }


def classify_site(counts, args):
    rates = site_rates(counts)
    if rates["t_depth"] < args.min_site_depth:
        return "low_depth_site"
    if rates["other_mismatch_rate"] > args.max_site_other_mismatch_rate:
        return "noisy_site"
    if (
        counts["tc_count"] >= args.min_site_tc_count
        and rates["tc_rate"] >= args.min_site_tc_rate
        and rates["other_mismatch_rate"] <= args.max_site_other_mismatch_rate
    ):
        return "high_conf_tc_site"
    return "background_site"


def binomial_tail(k, n, p):
    if n <= 0:
        return math.nan
    if k <= 0:
        return 1.0
    p = min(max(p, 0.0), 1.0)
    if scipy_binom is not None:
        return float(scipy_binom.sf(k - 1, n, p))
    if p <= 0.0:
        return 0.0
    if p >= 1.0:
        return 1.0
    if k <= n * p:
        return 1.0
    log_pmf = (
        math.lgamma(n + 1)
        - math.lgamma(k + 1)
        - math.lgamma(n - k + 1)
        + k * math.log(p)
        + (n - k) * math.log1p(-p)
    )
    term = math.exp(log_pmf) if log_pmf > -745 else 0.0
    tail = 0.0
    x = k
    while x <= n and term > 0.0:
        tail += term
        if tail >= 1.0:
            return 1.0
        x += 1
        if x > n:
            break
        term *= ((n - x + 1) / x) * (p / (1.0 - p))
    return min(max(tail, 0.0), 1.0)


def format_pvalue(value):
    if math.isnan(value):
        return "NA"
    return f"{value:.6g}"


def choose_background_rate(args, totals):
    if args.bg_tc_rate is not None:
        if args.bg_tc_rate < 0 or args.bg_tc_rate > 1:
            raise SystemExit("--bg-tc-rate must be between 0 and 1")
        return args.bg_tc_rate, "provided_bg_tc_rate"
    t_depth = totals["t_depth"]
    global_ta_rate = totals["ta_count"] / t_depth if t_depth else 0.0
    global_tg_rate = totals["tg_count"] / t_depth if t_depth else 0.0
    return max(global_ta_rate, global_tg_rate, 1e-5), "max_global_ta_tg_rate_not_control"


def add_site_counts(target, events_by_site):
    for events in events_by_site.values():
        target[consensus_event(events)] += 1


def molecule_counts(record, site_calls, noisy_site_mode):
    counts = empty_base_counter()
    high_conf_hits = 0
    noisy_hits = 0
    for site, events in record["site_events"].items():
        event = consensus_event(events)
        site_call = site_calls.get(site, "")
        if event == "deletion_count":
            counts[event] += 1
            continue
        if site_call == "noisy_site":
            noisy_hits += int(event in {"tc_count", "ta_count", "tg_count"})
            if noisy_site_mode == "exclude":
                continue
        counts[event] += 1
        high_conf_hits += int(site_call == "high_conf_tc_site")
    counts["t_covered"] = counts["match_t_count"] + counts["tc_count"] + counts["ta_count"] + counts["tg_count"]
    return counts, high_conf_hits, noisy_hits


def call_for_molecule(record, counts, high_conf_hits, noisy_hits, binom_p, args):
    t_covered = counts["t_covered"]
    other_rate = (counts["ta_count"] + counts["tg_count"]) / t_covered if t_covered else 0.0
    reasons = []
    if record["low_mapq_reads"] and record["low_mapq_reads"] == record["supporting_reads"]:
        reasons.append("low_mapq")
    if other_rate > args.max_read_other_mismatch_rate:
        reasons.append("high_other_mismatch")
    if args.noisy_site_mode == "ambiguous" and noisy_hits > 0:
        reasons.append("noisy_site_mismatch")
    if reasons:
        return "ambiguous", ",".join(reasons), "technical_ambiguous"
    if args.main_call == "fixed_threshold":
        if t_covered < args.min_t_covered and counts["tc_count"] == 0 and args.low_t_zero_tc_as_old:
            return "old", "low_t_zero_tc", "low_conf_old"
        if t_covered < args.min_t_covered:
            return "ambiguous", "low_t_covered", "low_t_ambiguous"
        if counts["tc_count"] >= args.new_threshold:
            return "new", "fixed_threshold", "high_conf_new"
        return "old", "no_tc_threshold", "high_conf_old"
    if t_covered < args.min_t_covered:
        if counts["tc_count"] == 0 and args.low_t_zero_tc_as_old:
            return "old", "low_t_zero_tc", "low_conf_old"
        low_t_new = (
            args.low_t_high_conf_tc_as_new
            and counts["tc_count"] > 0
            and high_conf_hits >= args.min_high_conf_tc_sites_hit
            and not math.isnan(binom_p)
            and binom_p < args.alpha
        )
        if low_t_new:
            return "new", "low_t_high_conf_tc", "low_conf_new"
        return "ambiguous", "low_t_covered", "low_t_ambiguous"
    is_new = binom_p < args.alpha
    if args.site_evidence_mode == "required":
        is_new = is_new and high_conf_hits >= args.min_high_conf_tc_sites_hit
    if is_new:
        return "new", "binomial_tc_signal", "high_conf_new"
    return "old", "not_significant", "high_conf_old"


def molecule_row(key, record, site_calls, bg_rate, bg_rate_source, args):
    barcode, x, y, gene_id, gene_name, umi = key
    counts, high_conf_hits, noisy_hits = molecule_counts(record, site_calls, args.noisy_site_mode)
    tc_events = counts["tc_count"]
    t_covered = counts["t_covered"]
    tc_rate = tc_events / t_covered if t_covered else 0.0
    other_rate = (counts["ta_count"] + counts["tg_count"]) / t_covered if t_covered else 0.0
    binom_p = binomial_tail(tc_events, t_covered, bg_rate) if t_covered > 0 else math.nan
    final_call, ambiguity_reason, call_confidence = call_for_molecule(
        record, counts, high_conf_hits, noisy_hits, binom_p, args
    )
    feature_classes = record["feature_classes"]
    if "exon" in feature_classes and "intron" in feature_classes:
        feature_class = "ambiguous_feature"
    elif "exon" in feature_classes:
        feature_class = "exon"
    elif "intron" in feature_classes:
        feature_class = "intron"
    else:
        feature_class = "ambiguous_feature"
    return {
        "barcode": barcode,
        "x": x,
        "y": y,
        "gene_id": gene_id,
        "gene_name": gene_name,
        "umi": umi,
        "supporting_reads": record["supporting_reads"],
        "tc_events": tc_events,
        "t_covered": t_covered,
        "feature_class": feature_class,
        "assignment": final_call,
        "final_call": final_call,
        "call_confidence": call_confidence,
        "ambiguity_reason": ambiguity_reason,
        "umi_t_covered": t_covered,
        "umi_tc_count": tc_events,
        "umi_ta_count": counts["ta_count"],
        "umi_tg_count": counts["tg_count"],
        "umi_tc_rate": f"{tc_rate:.8f}",
        "umi_other_mismatch_rate": f"{other_rate:.8f}",
        "high_conf_tc_sites_hit": high_conf_hits,
        "noisy_sites_hit": noisy_hits,
        "binomial_p": format_pvalue(binom_p),
        "bg_tc_rate": f"{bg_rate:.8g}",
        "bg_rate_source": bg_rate_source,
        "read_sum_tc_count": record["read_sum"]["tc_count"],
        "read_sum_t_covered": (
            record["read_sum"]["match_t_count"]
            + record["read_sum"]["tc_count"]
            + record["read_sum"]["ta_count"]
            + record["read_sum"]["tg_count"]
        ),
    }


def update_molecule_qc(qc, row):
    qc["molecules_total"] += 1
    qc[f"molecules_{row['assignment']}"] += 1
    qc["total_t_covered"] += int(row["umi_t_covered"])
    qc["total_tc_count"] += int(row["umi_tc_count"])
    qc["total_ta_count"] += int(row["umi_ta_count"])
    qc["total_tg_count"] += int(row["umi_tg_count"])


def finalize_qc(qc, feature_stats, site_call_counts, site_totals, excluded_sites, bg_rate, bg_rate_source, args, outdir):
    total_called = qc["molecules_new"] + qc["molecules_old"] + qc["molecules_ambiguous"]
    high_quality = qc["molecules_new"] + qc["molecules_old"]
    qc["new_fraction"] = qc["molecules_new"] / total_called if total_called else 0.0
    qc["old_fraction"] = qc["molecules_old"] / total_called if total_called else 0.0
    qc["ambiguous_fraction"] = qc["molecules_ambiguous"] / total_called if total_called else 0.0
    qc["new_fraction_excluding_ambiguous"] = qc["molecules_new"] / high_quality if high_quality else 0.0
    qc["global_tc_rate"] = qc["total_tc_count"] / qc["total_t_covered"] if qc["total_t_covered"] else 0.0
    qc["global_ta_rate"] = qc["total_ta_count"] / qc["total_t_covered"] if qc["total_t_covered"] else 0.0
    qc["global_tg_rate"] = qc["total_tg_count"] / qc["total_t_covered"] if qc["total_t_covered"] else 0.0
    qc["global_other_mismatch_rate"] = (
        (qc["total_ta_count"] + qc["total_tg_count"]) / qc["total_t_covered"] if qc["total_t_covered"] else 0.0
    )
    qc["bg_tc_rate"] = bg_rate
    qc["bg_rate_source"] = bg_rate_source
    qc["site_call_counts"] = dict(site_call_counts)
    qc["site_totals"] = dict(site_totals)
    qc["main_call"] = args.main_call
    qc["site_evidence_mode"] = args.site_evidence_mode
    qc["calling_parameters"] = vars(args)
    qc["excluded_snp_sites"] = len(excluded_sites)
    (outdir / f"{args.sample}.qc.json").write_text(json.dumps(qc, indent=2), encoding="utf-8")

    with (outdir / f"{args.sample}.exon_intron_conversion.summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", lineterminator="\n", fieldnames=["sample", "feature_class", "reads", "t_covered", "tc_events", "tc_rate"])
        writer.writeheader()
        for feature_class in ("exon", "intron"):
            item = feature_stats[feature_class]
            writer.writerow({
                "sample": args.sample,
                "feature_class": feature_class,
                "reads": item["reads"],
                "t_covered": item["t_covered"],
                "tc_events": item["tc_events"],
                "tc_rate": item["tc_events"] / item["t_covered"] if item["t_covered"] else 0.0,
            })


def run_partitioned(args, outdir, annotations, barcode_table, excluded_sites, bam, fasta):
    if args.bg_tc_rate is None:
        raise SystemExit("--partition-by-contig requires --bg-tc-rate")
    if not bam.has_index():
        raise SystemExit("--partition-by-contig requires an indexed BAM")
    bg_rate, bg_rate_source = choose_background_rate(args, Counter())
    qc = Counter()
    feature_stats = defaultdict(lambda: {"tc_events": 0, "t_covered": 0, "reads": 0})
    site_call_counts = Counter()
    site_totals = Counter()
    site_path = outdir / f"{args.sample}.slam_t_sites.tsv.gz"
    molecule_path = outdir / f"{args.sample}.molecule_calls.tsv.gz"

    with gzip.open(site_path, "wt") as site_handle, gzip.open(molecule_path, "wt") as molecule_handle:
        site_writer = csv.DictWriter(site_handle, delimiter="\t", lineterminator="\n", fieldnames=SITE_COLUMNS)
        molecule_writer = csv.DictWriter(molecule_handle, delimiter="\t", lineterminator="\n", fieldnames=MOLECULE_COLUMNS)
        site_writer.writeheader()
        molecule_writer.writeheader()
        for contig in bam.references:
            site_counts = defaultdict(empty_base_counter)
            molecule_reads = {}
            for read in bam.fetch(contig):
                qc["total_bam_reads"] += 1
                if args.progress_every and qc["total_bam_reads"] % args.progress_every == 0:
                    print(
                        f"[call_new_old] reads={qc['total_bam_reads']} molecules={qc['molecules_total'] + len(molecule_reads)} contig={contig}",
                        file=sys.stderr,
                        flush=True,
                    )
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    qc["skipped_non_primary"] += 1
                    continue
                meta = barcode_table.get(read.query_name)
                if meta is None:
                    qc["missing_barcode_metadata"] += 1
                    continue
                blocks = read.get_blocks()
                assign = assign_gene_feature(read.reference_name, blocks, annotations)
                if assign is None:
                    qc["unassigned_gene"] += 1
                    continue
                if assign.get("gene_id") is None:
                    qc["ambiguous_gene"] += 1
                    continue
                feature_class = assign["feature_class"]
                barcode, umi, x, y = meta
                key = (barcode, x, y, assign["gene_id"], assign["gene_name"], umi)
                entry = molecule_reads.setdefault(key, {"supporting_reads": 0, "low_mapq_reads": 0, "site_events": defaultdict(empty_base_counter), "read_sum": empty_base_counter(), "feature_classes": set()})
                entry["supporting_reads"] += 1
                entry["feature_classes"].add(feature_class)
                if read.mapping_quality < args.min_mapq:
                    entry["low_mapq_reads"] += 1
                    qc["low_mapq_reads"] += 1
                    continue
                read_sites = defaultdict(empty_base_counter)
                for chrom, rpos, event in iter_t_events(read, fasta, assign["strand"], args.min_tc_baseq):
                    if (chrom, rpos) in excluded_sites:
                        qc["excluded_snp_events"] += 1
                        continue
                    read_sites[(chrom, rpos)][event] += 1
                add_site_counts(entry["read_sum"], read_sites)
                for site, events in read_sites.items():
                    for event, value in events.items():
                        if value:
                            entry["site_events"][site][event] += value
                    site_counts[site][consensus_event(events)] += 1
                if feature_class in {"exon", "intron"}:
                    read_sum = empty_base_counter()
                    add_site_counts(read_sum, read_sites)
                    feature_stats[feature_class]["tc_events"] += read_sum["tc_count"]
                    feature_stats[feature_class]["t_covered"] += read_sum["match_t_count"] + read_sum["tc_count"] + read_sum["ta_count"] + read_sum["tg_count"]
                    feature_stats[feature_class]["reads"] += 1
                else:
                    qc["ambiguous_feature_reads"] += 1

            site_calls = {}
            for (chrom, rpos), counts in sorted(site_counts.items()):
                rates = site_rates(counts)
                site_call = classify_site(counts, args)
                site_calls[(chrom, rpos)] = site_call
                site_call_counts[site_call] += 1
                for key_name in BASE_EVENTS:
                    site_totals[key_name] += counts[key_name]
                site_totals["t_depth"] += rates["t_depth"]
                site_writer.writerow({
                    "chrom": chrom,
                    "pos": rpos + 1,
                    "ref_base": fasta.fetch(chrom, rpos, rpos + 1).upper(),
                    "t_depth": rates["t_depth"],
                    "match_t_count": counts["match_t_count"],
                    "tc_count": counts["tc_count"],
                    "ta_count": counts["ta_count"],
                    "tg_count": counts["tg_count"],
                    "deletion_count": counts["deletion_count"],
                    "tc_rate": f"{rates['tc_rate']:.8f}",
                    "ta_rate": f"{rates['ta_rate']:.8f}",
                    "tg_rate": f"{rates['tg_rate']:.8f}",
                    "other_mismatch_rate": f"{rates['other_mismatch_rate']:.8f}",
                    "site_call": site_call,
                })
            for key in sorted(molecule_reads):
                row = molecule_row(key, molecule_reads[key], site_calls, bg_rate, bg_rate_source, args)
                molecule_writer.writerow(row)
                update_molecule_qc(qc, row)
            print(
                f"[call_new_old] completed_contig={contig} reads={qc['total_bam_reads']} molecules={qc['molecules_total']}",
                file=sys.stderr,
                flush=True,
            )
    qc["partition_strategy"] = "contig_streaming"
    finalize_qc(qc, feature_stats, site_call_counts, site_totals, excluded_sites, bg_rate, bg_rate_source, args, outdir)


def main():
    parser = argparse.ArgumentParser(description="Call old/new molecules and exon/intron conversion.")
    parser.add_argument("--bam", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--annotation-dir", required=True)
    parser.add_argument("--barcode-table", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--new-threshold", type=int, default=1)
    parser.add_argument("--min-tc-baseq", type=int, default=15)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--min-t-covered", type=int, default=10)
    parser.add_argument("--bg-tc-rate", type=float, default=None)
    parser.add_argument("--alpha", type=float, default=0.01)
    parser.add_argument("--main-call", choices=["umi_site_binomial", "fixed_threshold"], default="umi_site_binomial")
    parser.add_argument("--min-site-depth", type=int, default=10)
    parser.add_argument("--min-site-tc-count", type=int, default=3)
    parser.add_argument("--min-site-tc-rate", type=float, default=0.2)
    parser.add_argument("--max-site-other-mismatch-rate", type=float, default=0.05)
    parser.add_argument("--max-read-other-mismatch-rate", type=float, default=0.1)
    parser.add_argument("--min-high-conf-tc-sites-hit", type=int, default=1)
    parser.add_argument("--site-evidence-mode", choices=["supportive", "required"], default="supportive")
    parser.add_argument("--noisy-site-mode", choices=["exclude", "ambiguous"], default="exclude")
    parser.add_argument("--low-t-zero-tc-as-old", action="store_true")
    parser.add_argument("--low-t-high-conf-tc-as-new", action="store_true")
    parser.add_argument("--progress-every", type=int, default=1000000)
    parser.add_argument("--exclude-vcf", default="", help="VCF positions excluded from conversion evidence.")
    parser.add_argument("--partition-by-contig", action="store_true", help="Bound memory by processing and writing one BAM contig at a time.")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    print("[call_new_old] loading annotations", file=sys.stderr, flush=True)
    annotations = load_annotations(args.annotation_dir)
    print("[call_new_old] loading barcode table", file=sys.stderr, flush=True)
    barcode_table = load_barcode_table(args.barcode_table)
    excluded_sites = load_excluded_sites(args.exclude_vcf)
    print(
        f"[call_new_old] loaded barcode_records={len(barcode_table)}",
        file=sys.stderr,
        flush=True,
    )
    bam = pysam.AlignmentFile(args.bam, "rb")
    fasta = pysam.FastaFile(args.fasta)

    if args.partition_by_contig:
        run_partitioned(args, outdir, annotations, barcode_table, excluded_sites, bam, fasta)
        bam.close()
        fasta.close()
        return

    molecule_reads = {}
    site_counts = defaultdict(empty_base_counter)
    qc = Counter()
    feature_stats = defaultdict(lambda: {"tc_events": 0, "t_covered": 0, "reads": 0})

    for read in bam.fetch(until_eof=True):
        qc["total_bam_reads"] += 1
        if args.progress_every and qc["total_bam_reads"] % args.progress_every == 0:
            print(
                f"[call_new_old] reads={qc['total_bam_reads']} molecules={len(molecule_reads)}",
                file=sys.stderr,
                flush=True,
            )
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            qc["skipped_non_primary"] += 1
            continue
        meta = barcode_table.get(read.query_name)
        if meta is None:
            qc["missing_barcode_metadata"] += 1
            continue
        blocks = read.get_blocks()
        assign = assign_gene_feature(read.reference_name, blocks, annotations)
        if assign is None:
            qc["unassigned_gene"] += 1
            continue
        if assign.get("gene_id") is None:
            qc["ambiguous_gene"] += 1
            continue
        feature_class = assign["feature_class"]
        barcode, umi, x, y = meta
        key = (barcode, x, y, assign["gene_id"], assign["gene_name"], umi)
        entry = molecule_reads.setdefault(
            key,
            {
                "supporting_reads": 0,
                "low_mapq_reads": 0,
                "site_events": defaultdict(empty_base_counter),
                "read_sum": empty_base_counter(),
                "feature_classes": set(),
            },
        )
        entry["supporting_reads"] += 1
        entry["feature_classes"].add(feature_class)
        if read.mapping_quality < args.min_mapq:
            entry["low_mapq_reads"] += 1
            qc["low_mapq_reads"] += 1
            continue
        read_sites = defaultdict(empty_base_counter)
        for chrom, rpos, event in iter_t_events(read, fasta, assign["strand"], args.min_tc_baseq):
            if (chrom, rpos) in excluded_sites:
                qc["excluded_snp_events"] += 1
                continue
            read_sites[(chrom, rpos)][event] += 1
        add_site_counts(entry["read_sum"], read_sites)
        for site, events in read_sites.items():
            for event, value in events.items():
                if value:
                    entry["site_events"][site][event] += value
            site_counts[site][consensus_event(events)] += 1
        tc_events = entry["read_sum"]["tc_count"]
        t_covered = (
            entry["read_sum"]["match_t_count"]
            + entry["read_sum"]["tc_count"]
            + entry["read_sum"]["ta_count"]
            + entry["read_sum"]["tg_count"]
        )
        if feature_class in {"exon", "intron"}:
            read_sum = empty_base_counter()
            add_site_counts(read_sum, read_sites)
            feature_stats[feature_class]["tc_events"] += read_sum["tc_count"]
            feature_stats[feature_class]["t_covered"] += (
                read_sum["match_t_count"] + read_sum["tc_count"] + read_sum["ta_count"] + read_sum["tg_count"]
            )
            feature_stats[feature_class]["reads"] += 1
        else:
            qc["ambiguous_feature_reads"] += 1

    site_calls = {}
    site_call_counts = Counter()
    site_totals = Counter()
    site_path = outdir / f"{args.sample}.slam_t_sites.tsv.gz"
    with gzip.open(site_path, "wt") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "chrom",
                "pos",
                "ref_base",
                "t_depth",
                "match_t_count",
                "tc_count",
                "ta_count",
                "tg_count",
                "deletion_count",
                "tc_rate",
                "ta_rate",
                "tg_rate",
                "other_mismatch_rate",
                "site_call",
            ]
        )
        for (chrom, rpos), counts in sorted(site_counts.items()):
            rates = site_rates(counts)
            site_call = classify_site(counts, args)
            site_calls[(chrom, rpos)] = site_call
            site_call_counts[site_call] += 1
            for key_name in BASE_EVENTS:
                site_totals[key_name] += counts[key_name]
            site_totals["t_depth"] += rates["t_depth"]
            writer.writerow(
                [
                    chrom,
                    rpos + 1,
                    fasta.fetch(chrom, rpos, rpos + 1).upper(),
                    rates["t_depth"],
                    counts["match_t_count"],
                    counts["tc_count"],
                    counts["ta_count"],
                    counts["tg_count"],
                    counts["deletion_count"],
                    f"{rates['tc_rate']:.8f}",
                    f"{rates['ta_rate']:.8f}",
                    f"{rates['tg_rate']:.8f}",
                    f"{rates['other_mismatch_rate']:.8f}",
                    site_call,
                ]
            )
    bg_rate, bg_rate_source = choose_background_rate(args, site_totals)

    molecule_rows = []
    for key, record in molecule_reads.items():
        barcode, x, y, gene_id, gene_name, umi = key
        counts, high_conf_hits, noisy_hits = molecule_counts(record, site_calls, args.noisy_site_mode)
        tc_events = counts["tc_count"]
        t_covered = counts["t_covered"]
        tc_rate = tc_events / t_covered if t_covered else 0.0
        other_rate = (counts["ta_count"] + counts["tg_count"]) / t_covered if t_covered else 0.0
        binom_p = binomial_tail(tc_events, t_covered, bg_rate) if t_covered > 0 else math.nan
        final_call, ambiguity_reason, call_confidence = call_for_molecule(
            record,
            counts,
            high_conf_hits,
            noisy_hits,
            binom_p,
            args,
        )
        feature_classes = record["feature_classes"]
        if "exon" in feature_classes and "intron" in feature_classes:
            feature_class = "ambiguous_feature"
        elif "exon" in feature_classes:
            feature_class = "exon"
        elif "intron" in feature_classes:
            feature_class = "intron"
        else:
            feature_class = "ambiguous_feature"
        molecule_rows.append(
            {
                "barcode": barcode,
                "x": x,
                "y": y,
                "gene_id": gene_id,
                "gene_name": gene_name,
                "umi": umi,
                "supporting_reads": record["supporting_reads"],
                "tc_events": tc_events,
                "t_covered": t_covered,
                "feature_class": feature_class,
                "assignment": final_call,
                "final_call": final_call,
                "call_confidence": call_confidence,
                "ambiguity_reason": ambiguity_reason,
                "umi_t_covered": t_covered,
                "umi_tc_count": tc_events,
                "umi_ta_count": counts["ta_count"],
                "umi_tg_count": counts["tg_count"],
                "umi_tc_rate": f"{tc_rate:.8f}",
                "umi_other_mismatch_rate": f"{other_rate:.8f}",
                "high_conf_tc_sites_hit": high_conf_hits,
                "noisy_sites_hit": noisy_hits,
                "binomial_p": format_pvalue(binom_p),
                "bg_tc_rate": f"{bg_rate:.8g}",
                "bg_rate_source": bg_rate_source,
                "read_sum_tc_count": record["read_sum"]["tc_count"],
                "read_sum_t_covered": (
                    record["read_sum"]["match_t_count"]
                    + record["read_sum"]["tc_count"]
                    + record["read_sum"]["ta_count"]
                    + record["read_sum"]["tg_count"]
                ),
            }
        )

    molecule_df = pd.DataFrame(molecule_rows).sort_values(["barcode", "gene_id", "umi"])
    molecule_path = outdir / f"{args.sample}.molecule_calls.tsv.gz"
    molecule_df.to_csv(molecule_path, sep="\t", index=False, compression="gzip")

    exon = feature_stats["exon"]
    intron = feature_stats["intron"]
    rate_rows = [
        {
            "sample": args.sample,
            "feature_class": "exon",
            "reads": exon["reads"],
            "t_covered": exon["t_covered"],
            "tc_events": exon["tc_events"],
            "tc_rate": exon["tc_events"] / exon["t_covered"] if exon["t_covered"] else 0.0,
        },
        {
            "sample": args.sample,
            "feature_class": "intron",
            "reads": intron["reads"],
            "t_covered": intron["t_covered"],
            "tc_events": intron["tc_events"],
            "tc_rate": intron["tc_events"] / intron["t_covered"] if intron["t_covered"] else 0.0,
        },
    ]
    pd.DataFrame(rate_rows).to_csv(outdir / f"{args.sample}.exon_intron_conversion.summary.tsv", sep="\t", index=False)

    qc["molecules_total"] = len(molecule_df)
    qc["molecules_new"] = int((molecule_df["assignment"] == "new").sum()) if not molecule_df.empty else 0
    qc["molecules_old"] = int((molecule_df["assignment"] == "old").sum()) if not molecule_df.empty else 0
    qc["molecules_ambiguous"] = int((molecule_df["assignment"] == "ambiguous").sum()) if not molecule_df.empty else 0
    qc["total_t_covered"] = int(molecule_df["umi_t_covered"].astype(int).sum()) if not molecule_df.empty else 0
    qc["total_tc_count"] = int(molecule_df["umi_tc_count"].astype(int).sum()) if not molecule_df.empty else 0
    qc["total_ta_count"] = int(molecule_df["umi_ta_count"].astype(int).sum()) if not molecule_df.empty else 0
    qc["total_tg_count"] = int(molecule_df["umi_tg_count"].astype(int).sum()) if not molecule_df.empty else 0
    total_called = qc["molecules_new"] + qc["molecules_old"] + qc["molecules_ambiguous"]
    high_quality = qc["molecules_new"] + qc["molecules_old"]
    qc["new_fraction"] = qc["molecules_new"] / total_called if total_called else 0.0
    qc["old_fraction"] = qc["molecules_old"] / total_called if total_called else 0.0
    qc["ambiguous_fraction"] = qc["molecules_ambiguous"] / total_called if total_called else 0.0
    qc["new_fraction_excluding_ambiguous"] = qc["molecules_new"] / high_quality if high_quality else 0.0
    qc["global_tc_rate"] = qc["total_tc_count"] / qc["total_t_covered"] if qc["total_t_covered"] else 0.0
    qc["global_ta_rate"] = qc["total_ta_count"] / qc["total_t_covered"] if qc["total_t_covered"] else 0.0
    qc["global_tg_rate"] = qc["total_tg_count"] / qc["total_t_covered"] if qc["total_t_covered"] else 0.0
    qc["global_other_mismatch_rate"] = (
        (qc["total_ta_count"] + qc["total_tg_count"]) / qc["total_t_covered"] if qc["total_t_covered"] else 0.0
    )
    qc["bg_tc_rate"] = bg_rate
    qc["bg_rate_source"] = bg_rate_source
    qc["site_call_counts"] = dict(site_call_counts)
    qc["site_totals"] = dict(site_totals)
    qc["main_call"] = args.main_call
    qc["site_evidence_mode"] = args.site_evidence_mode
    qc["calling_parameters"] = vars(args)
    qc["excluded_snp_sites"] = len(excluded_sites)
    (outdir / f"{args.sample}.qc.json").write_text(json.dumps(qc, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
