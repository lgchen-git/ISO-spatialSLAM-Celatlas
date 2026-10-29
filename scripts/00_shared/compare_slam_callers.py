#!/usr/bin/env python3
"""Compare legacy and partitioned SLAM caller outputs independent of row order."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path


def normalized_hash(path: Path) -> tuple[int, str]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header = next(handle).rstrip("\n")
        rows = sorted(line.rstrip("\n") for line in handle)
    digest = hashlib.sha256()
    digest.update((header + "\n").encode())
    for row in rows:
        digest.update((row + "\n").encode())
    return len(rows), digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy-dir", required=True)
    parser.add_argument("--partitioned-dir", required=True)
    parser.add_argument("--sample", required=True)
    args = parser.parse_args()
    legacy = Path(args.legacy_dir)
    partitioned = Path(args.partitioned_dir)
    old_qc = json.loads((legacy / f"{args.sample}.qc.json").read_text(encoding="utf-8"))
    new_qc = json.loads((partitioned / f"{args.sample}.qc.json").read_text(encoding="utf-8"))
    old_qc.pop("calling_parameters", None)
    new_qc.pop("calling_parameters", None)
    new_qc.pop("partition_strategy", None)
    if old_qc != new_qc:
        differing = sorted(key for key in set(old_qc) | set(new_qc) if old_qc.get(key) != new_qc.get(key))
        raise SystemExit(f"QC mismatch: {differing}")
    for suffix in ("molecule_calls.tsv.gz", "slam_t_sites.tsv.gz"):
        old_count, old_hash = normalized_hash(legacy / f"{args.sample}.{suffix}")
        new_count, new_hash = normalized_hash(partitioned / f"{args.sample}.{suffix}")
        if (old_count, old_hash) != (new_count, new_hash):
            raise SystemExit(f"Table mismatch for {suffix}: legacy={old_count},{old_hash} partitioned={new_count},{new_hash}")
        print(f"MATCH\t{suffix}\trows={old_count}\tsha256={old_hash}")
    print("MATCH\tqc.json")


if __name__ == "__main__":
    main()
