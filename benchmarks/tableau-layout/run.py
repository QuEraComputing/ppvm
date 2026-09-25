# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
"""Build, validate and benchmark tableau storage with one reproducible command."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import random
import shutil
import statistics
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from reference import OPERATIONS, apply, arbitrary_fixture, fixture, validate_algebra

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
PPVM_REV = "3befaf1b597033138d6c1bab394b967b09cd57fb"
QC_REV = "1e553b25fdf12c67b55ccedfbc2f48f4dea78fcc"
DEFAULT_SIZES = "8,31,32,33,63,64,65,127,128,129,255,256,257,512,1024,2048"


def command(args, **kwargs):
    print("+ " + " ".join(map(str, args)), file=sys.stderr, flush=True)
    return subprocess.run(list(map(str, args)), cwd=REPO, check=True, **kwargs)


def capture(args):
    return command(args, capture_output=True, text=True).stdout.strip()


def metadata(args, env):
    cpu = platform.processor()
    if platform.system() == "Darwin":
        cpu = capture(["sysctl", "-n", "machdep.cpu.brand_string"])
    return {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "ppvm_base": PPVM_REV,
        "harness_commit": capture(["git", "rev-parse", "HEAD"]),
        "harness_status": capture(["git", "status", "--short"]),
        "quantumclifford_commit": capture(
            ["git", "-C", args.quantumclifford, "rev-parse", "HEAD"]
        ),
        "platform": platform.platform(),
        "cpu": cpu,
        "rustc": capture(["rustc", "-Vv"]),
        "julia": capture(["julia", "--version"]),
        "python": sys.version,
        "arguments": vars(args),
        "environment": {
            k: env.get(k)
            for k in (
                "RUSTFLAGS",
                "JULIA_NUM_THREADS",
                "JULIA_CPU_TARGET",
                "OPENBLAS_NUM_THREADS",
            )
        },
        "source_sha256": {
            str(p.relative_to(REPO)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [
                *HERE.glob("*.py"),
                *HERE.glob("*.jl"),
                REPO / "crates/ppvm-tableau-2/examples/tableau_layout.rs",
                REPO / "crates/ppvm-tableau-2/examples/candidate_storage.rs",
            ]
            if p.exists()
        },
    }


def configurations(worker, public):
    if worker == "native":
        result = [
            ("ppvm-native-kernel", "column", 64, op)
            for op in (*OPERATIONS, "transpose_roundtrip")
        ]
        result += [
            ("ppvm-native-kernel", layout, 64, op)
            for layout in ("row", "column-transpose-batch")
            for op in ("comm", "mul")
        ]
        if public:
            result += [
                ("ppvm-public", "column", 64, op) for op in OPERATIONS if op != "comm"
            ]
        return result
    layouts = (
        ("qubit_bits", "generator_bits")
        if worker == "candidate"
        else ("fastrow", "fastcolumn")
    )
    return [
        (worker, layout, width, op)
        for layout in layouts
        for width in (8, 16, 32, 64, 128)
        for op in OPERATIONS
        if (layout, width, op) != ("fastrow", 128, "mul")
    ]


def run_worker(
    worker, argv, directory, expected, samples, min_ms, output, label, env, public
):
    path = output / f"{worker}-{label}.csv"
    log = output / f"{worker}-{label}.log"
    with path.open("w") as stdout, log.open("w") as stderr:
        command(
            [*argv, directory, samples, min_ms], env=env, stdout=stdout, stderr=stderr
        )
    with path.open() as file:
        rows = list(csv.DictReader(file))
    expected_keys = Counter(
        (impl, layout, str(width), str(n), op, str(sample))
        for impl, layout, width, op in configurations(worker, public)
        for n in expected
        for sample in range(1, samples + 1)
    )
    keys = ("implementation", "layout", "word_bits", "n", "operation", "sample")
    actual_keys = Counter(tuple(row[k] for k in keys) for row in rows)
    if actual_keys != expected_keys:
        raise AssertionError(
            f"{path}: missing={expected_keys - actual_keys}; extra={actual_keys - expected_keys}"
        )
    for row in rows:
        n, op = int(row["n"]), row["operation"]
        if row["checksum"] != expected[n][op]:
            raise AssertionError(
                f"{path}: logical result mismatch: {row}; expected {expected[n][op]}"
            )
        count = (
            2 * n
            if op in ("comm", "mul")
            else 3 * n
            if op == "circuit"
            else 1
            if op == "transpose_roundtrip"
            else n
        )
        if (
            int(row["operations_per_iteration"]) != count
            or min(int(row["iterations"]), float(row["elapsed_ns"])) <= 0
        ):
            raise AssertionError(f"{path}: invalid measurement: {row}")
    unsupported = [
        line for line in log.read_text().splitlines() if line.startswith("UNSUPPORTED")
    ]
    if worker == "quantumclifford":
        for n in expected:
            prefix = f"UNSUPPORTED quantumclifford layout=fastrow word_bits=128 n={n} operation=mul:"
            if sum(line.startswith(prefix) for line in unsupported) != 1:
                raise AssertionError(
                    f"{log}: missing expected UInt128 SIMD limitation at n={n}"
                )
        assert len(unsupported) == len(expected)
    else:
        assert not unsupported
    print(
        f"Validated {len(rows)} measurements from {worker} ({label}).",
        file=sys.stderr,
        flush=True,
    )
    return rows


def summarize(rows, output):
    keys = ("implementation", "layout", "word_bits", "n", "operation")
    groups = defaultdict(list)
    for row in rows:
        groups[tuple(row[k] for k in keys)].append(row)
    with (output / "summary.csv").open("w") as file:
        writer = csv.DictWriter(
            file,
            fieldnames=[
                *keys,
                "ns_per_op",
                "min_ns",
                "max_ns",
                "launch_min_ns",
                "launch_max_ns",
                "samples",
            ],
        )
        writer.writeheader()
        for key, values in sorted(groups.items()):
            timings = [
                float(r["elapsed_ns"])
                / int(r["iterations"])
                / int(r["operations_per_iteration"])
                for r in values
            ]
            launches = defaultdict(list)
            for row, timing in zip(values, timings):
                launches[row["launch"]].append(timing)
            medians = [statistics.median(v) for v in launches.values()]
            writer.writerow(
                dict(zip(keys, key))
                | {
                    "ns_per_op": statistics.median(medians),
                    "min_ns": min(timings),
                    "max_ns": max(timings),
                    "launch_min_ns": min(medians),
                    "launch_max_ns": max(medians),
                    "samples": len(timings),
                }
            )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--quantumclifford",
        type=Path,
        required=True,
        help="Checkout of pinned QuantumClifford.jl commit",
    )
    parser.add_argument("--sizes", default=DEFAULT_SIZES)
    parser.add_argument("--samples", type=int, default=5)
    parser.add_argument("--min-ms", type=float, default=10)
    parser.add_argument("--launches", type=int, default=3)
    parser.add_argument("--out", type=Path, default=REPO / "target/tableau-layout")
    parser.add_argument(
        "--julia-project",
        type=Path,
        help="Existing isolated Julia environment; otherwise set up under output",
    )
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--workers", default="native,candidate,quantumclifford")
    args = parser.parse_args()
    selected = args.workers.split(",")
    if set(selected) - {"native", "candidate", "quantumclifford"} or len(
        selected
    ) != len(set(selected)):
        parser.error("Unknown or duplicate worker")
    sizes = sorted(set(map(int, args.sizes.split(","))))
    if (
        not sizes
        or min(sizes) < 2
        or min(args.samples, args.launches, args.min_ms) <= 0
    ):
        parser.error("sizes must be >=2; samples, launches and min-ms must be positive")
    args.quantumclifford = args.quantumclifford.resolve()
    if capture(["git", "-C", args.quantumclifford, "rev-parse", "HEAD"]) != QC_REV:
        parser.error(f"QuantumClifford checkout must be at {QC_REV}")
    command(
        [
            "git",
            "-C",
            args.quantumclifford,
            "diff",
            "--exit-code",
            QC_REV,
            "--",
            "src",
            "lib/QECCore",
        ]
    )
    command(["git", "diff", "--exit-code", PPVM_REV, "--", "crates/ppvm-tableau-2/src"])
    output = args.out.resolve()
    output.mkdir(parents=True, exist_ok=True)
    fixtures = output / "fixtures"
    if fixtures.exists():
        # Remove only harness-owned fixture inputs so a narrower rerun does
        # not silently keep benchmarking sizes from a previous invocation.
        for old in [*fixtures.glob("n*.txt"), *fixtures.glob("n*.gates")]:
            old.unlink()
    fixtures.mkdir(exist_ok=True)
    validate_algebra()
    expected = {}
    for n in sizes:
        print(f"Preparing and checking reference at n={n}", file=sys.stderr, flush=True)
        rows = fixture(n, fixtures)
        expected[n] = {
            op: apply(rows, n, op) for op in (*OPERATIONS, "transpose_roundtrip")
        }
    (output / "expected.json").write_text(json.dumps(expected, indent=2) + "\n")
    validation = output / "validation"
    validation.mkdir(exist_ok=True)
    validation_expected = {}
    for n in (2, 3, 7, 9, 17, 33, 65, 129, 257):
        rows = arbitrary_fixture(n, validation)
        validation_expected[n] = {
            op: apply(rows, n, op) for op in (*OPERATIONS, "transpose_roundtrip")
        }
    (output / "validation-expected.json").write_text(
        json.dumps(validation_expected, indent=2) + "\n"
    )
    env = os.environ | {"JULIA_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}
    env.setdefault("RUSTFLAGS", "-C target-cpu=native")
    env.setdefault("JULIA_CPU_TARGET", "native")
    julia_project = (args.julia_project or output / "julia-env").resolve()
    if not args.skip_build:
        command(
            [
                "cargo",
                "build",
                "--release",
                "-p",
                "ppvm-tableau-2",
                "--example",
                "tableau_layout",
                "--example",
                "candidate_storage",
            ],
            env=env,
        )
        if "quantumclifford" in args.workers and args.julia_project is None:
            command(
                [
                    "julia",
                    "--startup-file=no",
                    HERE / "setup.jl",
                    args.quantumclifford,
                    julia_project,
                ],
                env=env,
            )
    workers = {
        "native": [REPO / "target/release/examples/tableau_layout"],
        "candidate": [REPO / "target/release/examples/candidate_storage"],
        "quantumclifford": [
            "julia",
            "--startup-file=no",
            "-t1",
            f"--project={julia_project}",
            HERE / "quantumclifford.jl",
        ],
    }
    info = metadata(args, env)
    if "quantumclifford" in selected:
        loaded = capture(
            [
                "julia",
                "--startup-file=no",
                f"--project={julia_project}",
                "-e",
                "using QuantumClifford; println(pkgdir(QuantumClifford)); println(pkgversion(QuantumClifford))",
            ]
        ).splitlines()
        if Path(loaded[0]).resolve() != args.quantumclifford:
            raise RuntimeError(
                f"Julia environment loads {loaded[0]}, expected {args.quantumclifford}"
            )
        info["quantumclifford_loaded"] = {"path": loaded[0], "version": loaded[1]}
        info["julia_manifest_sha256"] = hashlib.sha256(
            (julia_project / "Manifest.toml").read_bytes()
        ).hexdigest()
        for name in ("Project.toml", "Manifest.toml"):
            shutil.copyfile(julia_project / name, output / name)
    info["executable_sha256"] = {
        worker: hashlib.sha256(Path(workers[worker][0]).read_bytes()).hexdigest()
        for worker in selected
        if worker != "quantumclifford"
    }
    info["arguments"] = {
        k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()
    }
    (output / "metadata.json").write_text(json.dumps(info, indent=2) + "\n")
    (output / "unsupported.json").write_text(
        json.dumps(
            [
                {
                    "implementation": "quantumclifford",
                    "layout": "fastrow",
                    "word_bits": 128,
                    "operation": "mul",
                    "n": n,
                    "reason": "SIMD.jl rejects UInt128 vector lanes",
                }
                for n in sizes
                if "quantumclifford" in selected
            ],
            indent=2,
        )
        + "\n"
    )
    for worker in selected:
        run_worker(
            worker,
            workers[worker],
            validation,
            validation_expected,
            1,
            0.01,
            output,
            "validation",
            env,
            public=False,
        )
    results = []
    for launch in range(args.launches):
        order = selected.copy()
        random.Random(launch + 204).shuffle(order)
        for worker in order:
            print(
                f"Launch {launch + 1}/{args.launches}: {worker}",
                file=sys.stderr,
                flush=True,
            )
            rows = run_worker(
                worker,
                workers[worker],
                fixtures,
                expected,
                args.samples,
                args.min_ms,
                output,
                launch,
                env,
                public=True,
            )
            for row in rows:
                row["launch"] = launch
            results.extend(rows)
    with (output / "results.csv").open("w") as file:
        writer = csv.DictWriter(file, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)
    summarize(results, output)
    print(f"Validated {len(results)} measurements. Results: {output / 'summary.csv'}")


if __name__ == "__main__":
    main()
