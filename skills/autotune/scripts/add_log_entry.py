#!/usr/bin/env python3
from __future__ import annotations

import argparse
from datetime import date
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("log_file")
    parser.add_argument("entry")
    args = parser.parse_args()

    log_file = Path(args.log_file)
    today = date.today().isoformat()
    prefix = f"## {today}\n"
    content = log_file.read_text(encoding="utf-8") if log_file.exists() else ""
    with log_file.open("a", encoding="utf-8") as fh:
        if prefix not in content:
            if content and not content.endswith("\n"):
                fh.write("\n")
            fh.write(f"\n{prefix}")
        fh.write(f"- {args.entry}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
