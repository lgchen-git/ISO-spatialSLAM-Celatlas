#!/usr/bin/env python3
"""Read a top-level scalar value from the template YAML config."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Read a top-level YAML value.")
    parser.add_argument("--config", required=True, help="Path to config YAML")
    parser.add_argument("--key", required=True, help="Top-level key to read")
    parser.add_argument(
        "--required",
        action="store_true",
        help="Exit with an error if the key is missing or empty",
    )
    return parser.parse_args()


def strip_inline_comment(value: str) -> str:
    if "#" not in value:
        return value.strip()
    return value.split("#", 1)[0].strip()


def unquote(value: str) -> str:
    if len(value) >= 2 and value[0] == value[-1] and value[0] in {"'", '"'}:
        return value[1:-1]
    return value


def main() -> int:
    args = parse_args()
    config_path = Path(args.config)

    if not config_path.is_file():
        print(f"Config not found: {config_path}", file=sys.stderr)
        return 1

    value = ""

    with config_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            stripped = line.strip()

            if not stripped or stripped.startswith("#"):
                continue
            if raw_line[:1].isspace():
                continue
            if ":" not in line:
                continue

            name, raw_value = line.split(":", 1)
            if name.strip() != args.key:
                continue

            value = unquote(strip_inline_comment(raw_value))
            break

    if args.required and not value:
        print(
            f"Required config key is missing or empty: {args.key}",
            file=sys.stderr,
        )
        return 1

    print(value)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
