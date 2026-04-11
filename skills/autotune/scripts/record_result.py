#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("metric_file")
    parser.add_argument("--commit", required=True)
    parser.add_argument("--status", choices=("keep", "discard", "crash"), required=True)
    parser.add_argument("--description", required=True)
    parser.add_argument("--metric", action="append", default=[])
    args = parser.parse_args()

    entries = []
    for item in args.metric:
        key, value = item.split("=", 1)
        entries.append((key, float(value)))

    lines = [
        "[[metric]]",
        f'commit = "{args.commit}"',
        f'status = "{args.status}"',
        f'description = "{args.description}"',
    ]
    lines.extend(f"{key} = {value}" for key, value in entries)
    with Path(args.metric_file).open("a", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
