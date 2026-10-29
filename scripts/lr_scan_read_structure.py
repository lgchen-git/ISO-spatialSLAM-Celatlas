#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import json
import multiprocessing as mp
from collections import Counter
from dataclasses import dataclass
from itertools import islice
from typing import Iterator, Optional

try:
    import edlib
except ModuleNotFoundError:  # Use the compiled aligner already shipped with cutadapt.
    edlib = None
    from cutadapt._align import Aligner
    from cutadapt.align import EndSkip

ANCHOR_PATTERNS = (
    ("anchor_5p", "CTACACGACGCTCTTCCGATCT"),
    ("anchor_tso", "AAGCAGTGGTATCAACGCAGAGT"),
)
ANCHOR_PREFILTER_SEEDS = {
    "anchor_tso": ("AAGCAGT", "CAACGCA", "GCAGAGT"),
}
LINKER_1_VARIANTS = (
    "CGACTCACTACACGT",
    "TCGACTCACTACACGT",
    "GTCGACTCACTACACGT",
    "AGTCGACTCACTACACGT",
)
LINKER_2 = "TCGCTGACACGATCG"
RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")
SHORTLIST_SENTINEL = 99
CONFIDENT_ANCHOR_MAX_ED = 1
CONFIDENT_LINKER_TOTAL_ED = 2
SOURCE_RANK = {
    "anchor_5p": 0,
    "anchor_tso": 1,
    "linker1_fallback": 2,
}


@dataclass
class Config:
    anchor_max_ed_5p: int
    linker_max_ed: int
    linker_pos_tolerance: int
    umi_len: int
    search_window: int
    poly_trim_min_run: int
    min_read_len: int


@dataclass
class Candidate:
    source: str
    strand: str
    anchor_start: int
    anchor_end: int
    anchor_ed: int
    linker1_start: int
    linker1_end: int
    linker1_ed: int
    linker2_start: int
    linker2_end: int
    linker2_ed: int
    barcode_raw: str
    umi_seq: str
    cdna_start: int
    cdna_end: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Scan read structure and export barcode/linker coordinates.")
    parser.add_argument("--fastq", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--anchor-max-ed-5p", type=int, default=3)
    parser.add_argument("--linker-max-ed", type=int, default=3)
    parser.add_argument("--linker-pos-tolerance", type=int, default=4)
    parser.add_argument("--umi-len", type=int, default=10)
    parser.add_argument("--search-window", type=int, default=120)
    parser.add_argument("--poly-trim-min-run", type=int, default=10)
    parser.add_argument("--min-read-len", type=int, default=80)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--batch-size", type=int, default=1000)
    return parser.parse_args()


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8")


def iter_fastq(path: str) -> Iterator[tuple[str, str, str]]:
    with open_text(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline().strip()
            handle.readline()
            qual = handle.readline().strip()
            yield header.strip()[1:].split()[0], seq, qual


def reverse_complement(sequence: str) -> str:
    return sequence.translate(RC_TABLE)[::-1]


def align_locations(pattern: str, text: str, max_ed: int):
    if edlib is not None:
        result = edlib.align(pattern, text, mode="HW", task="locations", k=max_ed)
        if result["editDistance"] == -1:
            return []
        return [(start, end, result["editDistance"]) for start, end in result["locations"]]
    if not pattern or not text:
        return []
    aligner = Aligner(
        pattern,
        max_error_rate=max_ed / len(pattern),
        flags=EndSkip.SEMIGLOBAL,
        min_overlap=max(1, len(pattern) - max_ed),
    )
    located = aligner.locate(text)
    if located is None:
        return []
    _ref_start, _ref_stop, query_start, query_stop, _score, errors = located
    if errors > max_ed:
        return []
    return [(query_start, query_stop - 1, errors)]


def trim_terminal_homopolymers(sequence: str, min_run: int) -> tuple[int, int]:
    left, right = 0, len(sequence)
    while left < right:
        changed = False
        for base in ("A", "T"):
            run = 0
            while left + run < right and sequence[left + run] == base:
                run += 1
            if run >= min_run:
                left += run
                changed = True
                break
        if not changed:
            break
    while left < right:
        changed = False
        for base in ("A", "T"):
            run = 0
            while right - run - 1 >= left and sequence[right - run - 1] == base:
                run += 1
            if run >= min_run:
                right -= run
                changed = True
                break
        if not changed:
            break
    return left, right


def leading_t_run(sequence: str, start: int) -> int:
    run = 0
    while start + run < len(sequence) and sequence[start + run] == "T":
        run += 1
    return run


def candidate_sort_key(candidate: Candidate) -> tuple[int, int, int, int, int]:
    return (
        SOURCE_RANK.get(candidate.source, SHORTLIST_SENTINEL),
        candidate.anchor_ed if candidate.anchor_ed >= 0 else SHORTLIST_SENTINEL,
        candidate.linker1_ed + candidate.linker2_ed,
        0 if candidate.strand == "+" else 1,
        -(candidate.cdna_end - candidate.cdna_start),
    )


def has_prefilter_seed(text: str, seeds: tuple[str, ...]) -> bool:
    return any(seed in text for seed in seeds)


def is_confident_anchor_candidate(candidate: Candidate, config: Config) -> bool:
    return (
        candidate.source.startswith("anchor")
        and candidate.anchor_ed >= 0
        and candidate.anchor_ed <= min(config.anchor_max_ed_5p, CONFIDENT_ANCHOR_MAX_ED)
        and (candidate.linker1_ed + candidate.linker2_ed) <= CONFIDENT_LINKER_TOTAL_ED
        and (candidate.cdna_end - candidate.cdna_start) >= config.min_read_len
    )


def shortlist_sort_key(candidate: Candidate) -> tuple[int, int, int]:
    return (
        SOURCE_RANK.get(candidate.source, SHORTLIST_SENTINEL),
        candidate.anchor_ed if candidate.anchor_ed >= 0 else SHORTLIST_SENTINEL,
        candidate.linker1_ed + candidate.linker2_ed,
    )


def find_candidates(sequence: str, prefix: str, config: Config, strand: str) -> list[Candidate]:
    out: list[Candidate] = []
    for anchor_source, anchor_pattern in ANCHOR_PATTERNS:
        prefilter_seeds = ANCHOR_PREFILTER_SEEDS.get(anchor_source)
        if prefilter_seeds is not None and not has_prefilter_seed(prefix, prefilter_seeds):
            continue
        for anchor_start, anchor_end, anchor_ed in align_locations(anchor_pattern, prefix, config.anchor_max_ed_5p):
            b1_start = anchor_end + 1
            b1_end = b1_start + 4
            if b1_end >= len(sequence):
                continue
            b1 = sequence[b1_start:b1_end]
            l1_expected = b1_end
            for linker_1 in LINKER_1_VARIANTS:
                w1_start = max(0, l1_expected - 2)
                w1_end = min(len(sequence), l1_expected + len(linker_1) + 2)
                for rel1_start, rel1_end, l1_ed in align_locations(linker_1, sequence[w1_start:w1_end], config.linker_max_ed):
                    l1_start = w1_start + rel1_start
                    if abs(l1_start - l1_expected) > config.linker_pos_tolerance:
                        continue
                    l1_end = w1_start + rel1_end + 1
                    b2_start = l1_end
                    b2_end = b2_start + 4
                    if b2_end >= len(sequence):
                        continue
                    b2 = sequence[b2_start:b2_end]
                    l2_expected = b2_end
                    w2_start = max(0, l2_expected - 2)
                    w2_end = min(len(sequence), l2_expected + len(LINKER_2) + 2)
                    for rel2_start, rel2_end, l2_ed in align_locations(LINKER_2, sequence[w2_start:w2_end], config.linker_max_ed):
                        l2_start = w2_start + rel2_start
                        if abs(l2_start - l2_expected) > config.linker_pos_tolerance:
                            continue
                        l2_end = w2_start + rel2_end + 1
                        b3_start = l2_end
                        b3_end = b3_start + 4
                        umi_start = b3_end
                        umi_end = umi_start + config.umi_len
                        if umi_end > len(sequence):
                            continue
                        b3 = sequence[b3_start:b3_end]
                        poly_t_run = leading_t_run(sequence, umi_end)
                        cdna_start = umi_end + poly_t_run if poly_t_run >= config.poly_trim_min_run else umi_end
                        left_trim, right_trim = trim_terminal_homopolymers(sequence[cdna_start:], config.poly_trim_min_run)
                        cdna_left = cdna_start + left_trim
                        cdna_right = cdna_start + right_trim
                        out.append(
                            Candidate(
                                source=anchor_source,
                                strand=strand,
                                anchor_start=anchor_start,
                                anchor_end=anchor_end + 1,
                                anchor_ed=anchor_ed,
                                linker1_start=l1_start,
                                linker1_end=l1_end,
                                linker1_ed=l1_ed,
                                linker2_start=l2_start,
                                linker2_end=l2_end,
                                linker2_ed=l2_ed,
                                barcode_raw=b1 + b2 + b3,
                                umi_seq=sequence[umi_start:umi_end],
                                cdna_start=cdna_left,
                                cdna_end=cdna_right,
                            )
                        )
    return out


def find_candidates_from_linker1(sequence: str, prefix: str, config: Config, strand: str) -> list[Candidate]:
    out: list[Candidate] = []
    for linker_1 in LINKER_1_VARIANTS:
        for l1_start, l1_end_raw, l1_ed in align_locations(linker_1, prefix, config.linker_max_ed):
            l1_end = l1_end_raw + 1
            b1_start = l1_start - 4
            b1_end = l1_start
            if b1_start < 0:
                continue
            b1 = sequence[b1_start:b1_end]
            b2_start = l1_end
            b2_end = b2_start + 4
            if b2_end >= len(sequence):
                continue
            b2 = sequence[b2_start:b2_end]
            l2_expected = b2_end
            w2_start = max(0, l2_expected - 2)
            w2_end = min(len(sequence), l2_expected + len(LINKER_2) + 2)
            for rel2_start, rel2_end, l2_ed in align_locations(LINKER_2, sequence[w2_start:w2_end], config.linker_max_ed):
                l2_start = w2_start + rel2_start
                if abs(l2_start - l2_expected) > config.linker_pos_tolerance:
                    continue
                l2_end = w2_start + rel2_end + 1
                b3_start = l2_end
                b3_end = b3_start + 4
                umi_start = b3_end
                umi_end = umi_start + config.umi_len
                if umi_end > len(sequence):
                    continue
                b3 = sequence[b3_start:b3_end]
                poly_t_run = leading_t_run(sequence, umi_end)
                cdna_start = umi_end + poly_t_run if poly_t_run >= config.poly_trim_min_run else umi_end
                left_trim, right_trim = trim_terminal_homopolymers(sequence[cdna_start:], config.poly_trim_min_run)
                cdna_left = cdna_start + left_trim
                cdna_right = cdna_start + right_trim
                out.append(
                    Candidate(
                        source="linker1_fallback",
                        strand=strand,
                        anchor_start=-1,
                        anchor_end=-1,
                        anchor_ed=-1,
                        linker1_start=l1_start,
                        linker1_end=l1_end,
                        linker1_ed=l1_ed,
                        linker2_start=l2_start,
                        linker2_end=l2_end,
                        linker2_ed=l2_ed,
                        barcode_raw=b1 + b2 + b3,
                        umi_seq=sequence[umi_start:umi_end],
                        cdna_start=cdna_left,
                        cdna_end=cdna_right,
                    )
                )
    return out


def best_candidate(seq: str, config: Config) -> tuple[str, Optional[Candidate]]:
    prefix = seq[: config.search_window]
    plus_candidates = find_candidates(seq, prefix, config, "+")
    confident_plus = [candidate for candidate in plus_candidates if is_confident_anchor_candidate(candidate, config)]
    if confident_plus:
        return "candidate_found", min(confident_plus, key=candidate_sort_key)

    candidates = list(plus_candidates)

    rc_seq: str | None = None
    if True:
        rc_seq = reverse_complement(seq)
        rc_prefix = rc_seq[: config.search_window]
        candidates.extend(find_candidates(rc_seq, rc_prefix, config, "-"))
    if not candidates:
        candidates = find_candidates_from_linker1(seq, prefix, config, "+")
        if rc_seq is None:
            rc_seq = reverse_complement(seq)
        rc_prefix = rc_seq[: config.search_window]
        candidates.extend(find_candidates_from_linker1(rc_seq, rc_prefix, config, "-"))
    if not candidates:
        return "anchor_not_found", None

    best_short = min(candidates, key=shortlist_sort_key)
    if (best_short.cdna_end - best_short.cdna_start) < config.min_read_len:
        kept = [candidate for candidate in candidates if (candidate.cdna_end - candidate.cdna_start) >= config.min_read_len]
        if not kept:
            return "short_cdna", best_short
        return "candidate_found", min(kept, key=candidate_sort_key)
    kept = [candidate for candidate in candidates if (candidate.cdna_end - candidate.cdna_start) >= config.min_read_len]
    return "candidate_found", min(kept, key=candidate_sort_key)


def render_row(read_id: str, seq: str, status: str, candidate: Optional[Candidate]) -> str:
    if candidate is None:
        return f"{read_id}\t{len(seq)}\t{status}\t.\t.\t.\t.\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t0\n"
    return (
        f"{read_id}\t{len(seq)}\t{status}\t{candidate.source}\t{candidate.strand}\t{candidate.barcode_raw}\t{candidate.umi_seq}\t"
        f"{candidate.anchor_start}\t{candidate.anchor_end}\t{candidate.anchor_ed}\t"
        f"{candidate.linker1_start}\t{candidate.linker1_end}\t{candidate.linker1_ed}\t"
        f"{candidate.linker2_start}\t{candidate.linker2_end}\t{candidate.linker2_ed}\t"
        f"{candidate.cdna_start}\t{candidate.cdna_end}\t{candidate.cdna_end - candidate.cdna_start}\n"
    )


def process_read(read_id: str, seq: str, config: Config) -> tuple[str, Counter, Counter]:
    summary = Counter()
    strand_counts = Counter()
    summary["total_reads"] += 1
    status, candidate = best_candidate(seq, config)
    if candidate is None:
        summary[status] += 1
        return render_row(read_id, seq, status, candidate), summary, strand_counts
    strand_counts[candidate.strand] += 1
    if status == "short_cdna":
        summary["short_cdna"] += 1
    else:
        summary["candidate_found"] += 1
        summary[f"candidate_found_{candidate.source}"] += 1
    return render_row(read_id, seq, status, candidate), summary, strand_counts


def iter_batches(iterator: Iterator[tuple[str, str, str]], batch_size: int) -> Iterator[list[tuple[str, str]]]:
    while True:
        batch = [(read_id, seq) for read_id, seq, _qual in islice(iterator, batch_size)]
        if not batch:
            return
        yield batch


def process_batch(batch: list[tuple[str, str]], config: Config) -> tuple[list[str], Counter, Counter]:
    rows: list[str] = []
    summary = Counter()
    strand_counts = Counter()
    for read_id, seq in batch:
        row, row_summary, row_strands = process_read(read_id, seq, config)
        rows.append(row)
        summary.update(row_summary)
        strand_counts.update(row_strands)
    return rows, summary, strand_counts


def process_batch_star(args: tuple[list[tuple[str, str]], Config]) -> tuple[list[str], Counter, Counter]:
    return process_batch(*args)


def main() -> None:
    args = parse_args()
    config = Config(
        anchor_max_ed_5p=args.anchor_max_ed_5p,
        linker_max_ed=args.linker_max_ed,
        linker_pos_tolerance=args.linker_pos_tolerance,
        umi_len=args.umi_len,
        search_window=args.search_window,
        poly_trim_min_run=args.poly_trim_min_run,
        min_read_len=args.min_read_len,
    )
    out_tsv = f"{args.outdir}/read_structure.tsv.gz"
    out_json = f"{args.outdir}/structure_summary.json"

    summary = Counter()
    strand_counts = Counter()

    import os

    os.makedirs(args.outdir, exist_ok=True)

    with gzip.open(out_tsv, "wt") as out_handle:
        out_handle.write(
            "read_id\tread_len\tscan_status\tstructure_source\tstrand\tbarcode_raw\tumi_seq\tanchor_start\tanchor_end\tanchor_ed\t"
            "linker1_start\tlinker1_end\tlinker1_ed\tlinker2_start\tlinker2_end\tlinker2_ed\tcdna_start\tcdna_end\tcdna_len\n"
        )
        batches = ((batch, config) for batch in iter_batches(iter_fastq(args.fastq), max(1, args.batch_size)))
        if args.workers > 1:
            with mp.Pool(processes=args.workers) as pool:
                for rows, batch_summary, batch_strands in pool.imap(process_batch_star, batches, chunksize=1):
                    summary.update(batch_summary)
                    strand_counts.update(batch_strands)
                    out_handle.writelines(rows)
        else:
            for rows, batch_summary, batch_strands in map(process_batch_star, batches):
                summary.update(batch_summary)
                strand_counts.update(batch_strands)
                out_handle.writelines(rows)

    report = {
        "sample_id": args.sample_id,
        "total_reads": summary["total_reads"],
        "candidate_found": summary["candidate_found"],
        "anchor_not_found": summary["anchor_not_found"],
        "short_cdna": summary["short_cdna"],
        "strand_counts": dict(strand_counts),
        "source_counts": {
            key.replace("candidate_found_", ""): value
            for key, value in sorted(summary.items())
            if key.startswith("candidate_found_")
        },
    }
    with open(out_json, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)


if __name__ == "__main__":
    main()
