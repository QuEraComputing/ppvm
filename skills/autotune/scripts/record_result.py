#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("metric_file")
    parser.add_argument("--commit", required=True)
    parser.add_argument("--status", choices=("keep", "discard", "crash"), required=True)
    parser.add_argument("--description", required=True)
    parser.add_argument("--metric", action="append", default=[])
    args = parser.parse_args()

    if not args.metric:
        parser.error("at least one --metric is required")

    entries = []
    for item in args.metric:
        try:
            key, value = item.split("=", 1)
        except ValueError:
            parser.error(f"invalid --metric value: {item!r}; expected key=value")
        if not key or not value:
            parser.error(f"invalid --metric value: {item!r}; expected key=value")
        try:
            entries.append((key, float(value)))
        except ValueError:
            parser.error(f"invalid --metric value: {item!r}; expected key=value with a numeric value")

    lines = [
        "[[metric]]",
        f"{json.dumps('commit')} = {json.dumps(args.commit)}",
        f"{json.dumps('status')} = {json.dumps(args.status)}",
        f"{json.dumps('description')} = {json.dumps(args.description)}",
    ]
    lines.extend(f"{json.dumps(key)} = {value}" for key, value in entries)
    with Path(args.metric_file).open("a", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
