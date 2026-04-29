"""regen-stim entry point."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="regen-stim", description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)

    sub.add_parser("codes", help="(Task 7) Generate generated/codes/ via stim gen sweeps")
    sub.add_parser("noise-sweeps", help="(Task 8) Generate generated/noise_sweeps/ per-channel")
    sub.add_parser("dialect", help="(Task 9) Generate generated/dialect/ ppvm-specific")
    sub.add_parser("random", help="(Task 10) Generate generated/random/ random-walk")
    sub.add_parser("unsupported", help="(Task 11) Generate unsupported/ phase-1-rejected")
    p_refresh = sub.add_parser("refresh", help="Re-record one fixture by path")
    p_refresh.add_argument("path", type=Path)
    p_verify = sub.add_parser("verify", help="Re-run cross-check; error on mismatch; write nothing")
    p_verify.add_argument("path", type=Path)
    sub.add_parser("all", help="Run every regen subcommand")

    args = parser.parse_args(argv)

    if args.cmd == "codes":
        from . import codes
        return codes.run()
    if args.cmd == "noise-sweeps":
        from . import noise_sweeps
        return noise_sweeps.run()
    if args.cmd == "dialect":
        from . import dialect
        return dialect.run()
    if args.cmd == "random":
        from . import random_walk
        return random_walk.run()
    if args.cmd == "unsupported":
        from . import unsupported
        return unsupported.run()
    if args.cmd == "refresh":
        return _refresh(args.path)
    if args.cmd == "verify":
        return _verify(args.path)
    if args.cmd == "all":
        rc = 0
        for mod_name in ("codes", "noise_sweeps", "dialect", "random_walk", "unsupported"):
            mod = __import__(f"regen_stim.{mod_name}", fromlist=["run"])
            rc |= mod.run()
        return rc

    parser.print_help()
    return 2


def _refresh(path: Path) -> int:
    """Re-emit one fixture given a path to <name>.stim or <name>.expected.json."""
    import json

    from . import core

    stim_path = path if path.suffix == ".stim" else path.with_suffix("").with_suffix(".stim")
    if stim_path.suffix == ".expected" and stim_path.stem.endswith(".expected"):
        stim_path = stim_path.with_name(stim_path.stem.removesuffix(".expected") + ".stim")
    if not stim_path.exists():
        print(f"refresh: no .stim file at {stim_path}", file=sys.stderr)
        return 1

    json_path = stim_path.with_name(stim_path.stem + ".expected.json")
    if not json_path.exists():
        print(f"refresh: no expected.json at {json_path}", file=sys.stderr)
        return 1

    existing = json.loads(json_path.read_text())
    src = stim_path.read_text()
    paths = core.CorpusPaths.default()
    category = stim_path.parent.relative_to(paths.root).as_posix()
    name = stim_path.stem

    mode = existing.get("mode")
    if mode == "distribution":
        meta = core.FixtureMeta(
            name=name,
            category=category,
            source=src,
            test_num_shots=existing["num_shots"],
            stim_num_shots=existing.get("stim_num_shots", core.DEFAULT_STIM_SHOTS),
            stim_seed=existing.get("stim_seed", 0),
            tolerance_sigma=existing.get("tolerance_sigma_at_regen", core.DEFAULT_TOLERANCE_SIGMA),
        )
        core.write_distribution_fixture(meta, paths)
    elif mode == "deterministic":
        meta = core.FixtureMeta(name=name, category=category, source=src, test_num_shots=1)
        core.write_deterministic_fixture(meta, paths)
    elif mode == "unsupported":
        meta = core.FixtureMeta(
            name=name,
            category=category,
            source=src,
            test_num_shots=0,
            stim_num_shots=existing.get("stim_num_shots", core.DEFAULT_STIM_SHOTS),
            stim_seed=existing.get("stim_seed", 0),
        )
        core.write_unsupported_fixture(
            meta, paths, awaiting_phase2_instruction=existing["awaiting_phase2_instruction"]
        )
    else:
        print(f"refresh: unknown mode {mode!r}", file=sys.stderr)
        return 1
    print(f"refreshed {category}/{name}")
    return 0


def _verify(path: Path) -> int:
    """Re-run cross-check without writing; error on mismatch."""
    import json

    from . import core

    stim_path = path if path.suffix == ".stim" else path.with_suffix("").with_suffix(".stim")
    json_path = stim_path.with_name(stim_path.stem + ".expected.json")
    existing = json.loads(json_path.read_text())
    src = stim_path.read_text()
    if existing["mode"] != "distribution":
        print("verify: only meaningful for distribution mode", file=sys.stderr)
        return 0
    ref = core.run_stim(
        src,
        num_shots=existing["stim_num_shots"],
        seed=existing.get("stim_seed", 0),
    )
    ppvm_run = core.run_ppvm(
        src,
        num_shots=existing["num_shots"],
        seed=existing["ppvm_seed"],
    )
    ok = core.within_tolerance(
        ppvm_run.bit_means,
        ref.bit_means,
        existing["num_shots"],
        existing.get("tolerance_sigma_at_regen", core.DEFAULT_TOLERANCE_SIGMA),
    )
    if not ok:
        print(f"verify: mismatch at {stim_path}", file=sys.stderr)
        print(f"  stim_means: {ref.bit_means}", file=sys.stderr)
        print(f"  ppvm_means: {ppvm_run.bit_means}", file=sys.stderr)
        return 1
    print(f"verify: {stim_path} OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
