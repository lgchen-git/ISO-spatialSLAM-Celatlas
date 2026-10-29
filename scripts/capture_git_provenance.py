#!/usr/bin/env python3
"""Capture Git provenance, including the workspace's .git-local fallback."""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlsplit, urlunsplit


RELEASE_TAG = re.compile(r"^pipeline-v\d+\.\d+\.\d+$")


def git_env(repo: Path) -> dict[str, str]:
    env = os.environ.copy()
    fallback = repo / ".git-local"
    if fallback.is_dir() and not (repo / ".git" / "HEAD").exists():
        env["GIT_DIR"] = str(fallback)
        env["GIT_WORK_TREE"] = str(repo)
    return env


def git(repo: Path, *args: str, binary: bool = False) -> str | bytes:
    result = subprocess.run(
        ["git", "-C", str(repo), *args],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=not binary,
        env=git_env(repo),
    )
    if result.returncode != 0:
        error = result.stderr if isinstance(result.stderr, str) else result.stderr.decode(errors="replace")
        raise RuntimeError(error.strip() or f"git {' '.join(args)} failed")
    return result.stdout


def sanitize_remote(value: str) -> str:
    value = value.strip()
    if not value or "://" not in value:
        return value.split("@", 1)[-1]
    parsed = urlsplit(value)
    host = parsed.hostname or ""
    if parsed.port:
        host = f"{host}:{parsed.port}"
    return urlunsplit((parsed.scheme, host, parsed.path, parsed.query, parsed.fragment))


def digest(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--dirty-policy", choices=("allow", "warn", "error"), default="warn")
    args = parser.parse_args()
    repo = args.repo_root.resolve()
    commit = str(git(repo, "rev-parse", "--verify", "HEAD")).strip()
    try:
        branch = str(git(repo, "symbolic-ref", "--quiet", "--short", "HEAD")).strip()
    except RuntimeError:
        branch = "DETACHED"
    status = str(git(repo, "status", "--porcelain=v1", "--untracked-files=all"))
    dirty = bool(status)
    if dirty and args.dirty_policy == "error":
        raise SystemExit("ERROR: Git worktree is dirty")
    if dirty and args.dirty_policy == "warn":
        print("WARNING: Git worktree is dirty; hashes are recorded")
    diff = git(repo, "diff", "--binary", "HEAD", "--", binary=True)
    describe = str(git(repo, "describe", "--always", "--dirty", "--tags")).strip()
    try:
        remote = sanitize_remote(str(git(repo, "remote", "get-url", "origin")))
    except RuntimeError:
        remote = ""
    tags = str(git(repo, "tag", "--points-at", "HEAD")).splitlines()
    release = next((tag for tag in sorted(tags) if RELEASE_TAG.fullmatch(tag)), "")
    rows = (
        ("schema_version", "1"),
        ("captured_at_utc", datetime.now(timezone.utc).isoformat()),
        ("repo_name", repo.name),
        ("commit", commit),
        ("branch", branch or "DETACHED"),
        ("describe", describe),
        ("dirty", str(dirty).lower()),
        ("tracked_diff_sha256", digest(diff if isinstance(diff, bytes) else diff.encode())),
        ("status_sha256", digest(status.encode())),
        ("changed_path_count", str(len(status.splitlines()))),
        ("remote_url", remote),
        ("pipeline_release_tag", release),
        ("git_metadata_mode", "fallback_.git-local" if (repo / ".git-local").is_dir() else "standard"),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        "key\tvalue\n" + "".join(f"{key}\t{value}\n" for key, value in rows), encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
