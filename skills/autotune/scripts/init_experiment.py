#!/usr/bin/env python3
from __future__ import annotations

import argparse
from datetime import date
from pathlib import Path
import re


def slugify(value: str) -> str:
    return re.sub(r"-{2,}", "-", re.sub(r"[^a-z0-9]+", "-", value.lower())).strip("-")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("name")
    parser.add_argument("--root", default="docs/autotune")
    args = parser.parse_args()

    slug = slugify(args.name)
    if not slug:
        parser.error("task name must contain at least one ASCII letter or digit")

    task = f"{date.today().isoformat()}-{slug}"
    root = Path(args.root) / task
    root.mkdir(parents=True, exist_ok=True)

    metric = root / "metric.toml"
    log = root / "log.md"
    if not metric.exists():
        metric.write_text("", encoding="utf-8")
    if not log.exists():
        log.write_text(f"# Log for {task}\n\n## {date.today().isoformat()}\n", encoding="utf-8")

    print(task)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
