#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


RESERVED_METRIC_KEYS = {"commit", "status", "description"}


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
    seen_keys = set()
    for item in args.metric:
        try:
            key, value = item.split("=", 1)
        except ValueError:
            parser.error(f"invalid --metric value: {item!r}; expected key=value")
        if not key or not value:
            parser.error(f"invalid --metric value: {item!r}; expected key=value")
        if key in RESERVED_METRIC_KEYS:
            parser.error(f"invalid --metric value: {item!r}; {key!r} is a reserved metric key")
        if key in seen_keys:
            parser.error(f"invalid --metric value: {item!r}; duplicate metric key {key!r} in a single invocation")
        seen_keys.add(key)
        try:
            numeric_value = float(value)
        except ValueError:
            parser.error(f"invalid --metric value: {item!r}; expected key=value with a numeric value")
        if not math.isfinite(numeric_value):
            parser.error(f"invalid --metric value: {item!r}; expected key=value with a finite numeric value")
        entries.append((key, numeric_value))

    lines = [
        "[[metric]]",
        f"{json.dumps('commit')} = {json.dumps(args.commit)}",
        f"{json.dumps('status')} = {json.dumps(args.status)}",
        f"{json.dumps('description')} = {json.dumps(args.description)}",
    ]
    lines.extend(f"{json.dumps(key)} = {value}" for key, value in entries)
    metric_path = Path(args.metric_file)
    with metric_path.open("a", encoding="utf-8") as fh:
        if metric_path.exists() and metric_path.stat().st_size > 0:
            with metric_path.open("rb") as existing:
                existing.seek(-1, 2)
                if existing.read(1) != b"\n":
                    fh.write("\n")
        fh.write("\n".join(lines) + "\n\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
