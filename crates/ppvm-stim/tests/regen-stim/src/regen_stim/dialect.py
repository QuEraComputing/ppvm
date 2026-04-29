"""generated/dialect/: ppvm-specific dialect (I[R_X(...)], S[T], etc.).

Stim cannot simulate these so there's no oracle. Where the outcome is uniquely
determined we use deterministic mode and assert a hand-derived bitstring; where
there is measurement randomness we record ppvm's output to lock down regression
behavior.
"""

from __future__ import annotations

import json

from . import core


# (name, source, mode, asserted_bitstring | None, num_shots | 0)
DIALECT_FIXTURES: list[tuple[str, str, str, list[bool] | None, int]] = [
    # Deterministic — rotations that net out to identity or a Pauli.
    (
        "rx_pi_flips",
        "I[R_X(theta=1.0*pi)] 0\nM 0\n",
        "deterministic", [True], 0,
    ),
    (
        "ry_pi_flips",
        "I[R_Y(theta=1.0*pi)] 0\nM 0\n",
        "deterministic", [True], 0,
    ),
    (
        "rz_pi_no_flip",
        "I[R_Z(theta=1.0*pi)] 0\nM 0\n",
        "deterministic", [False], 0,
    ),
    (
        "rx_2pi_no_flip",
        "I[R_X(theta=2.0*pi)] 0\nM 0\n",
        "deterministic", [False], 0,
    ),
    (
        "rx_then_inverse",
        "I[R_X(theta=0.5*pi)] 0\nI[R_X(theta=1.5*pi)] 0\nM 0\n",
        "deterministic", [False], 0,
    ),
    (
        "u3_x_equivalent",
        "I[U3(theta=1.0*pi, phi=0.0, lambda=0.0)] 0\nM 0\n",
        "deterministic", [True], 0,
    ),
    (
        "t_then_t_dag_identity",
        "X 0\nS[T] 0\nS_DAG[T] 0\nM 0\n",
        "deterministic", [True], 0,
    ),
    (
        "t_eight_times_identity",
        "X 0\n" + "S[T] 0\n" * 8 + "M 0\n",
        "deterministic", [True], 0,
    ),
    # Distribution — rotations into a superposition; ppvm-only oracle.
    (
        "rx_half_pi_random",
        "I[R_X(theta=0.5*pi)] 0\nM 0\n",
        "distribution", None, 256,
    ),
    (
        "ry_half_pi_random",
        "I[R_Y(theta=0.5*pi)] 0\nM 0\n",
        "distribution", None, 256,
    ),
    (
        "u3_h_equivalent_random",
        "I[U3(theta=0.5*pi, phi=0.0, lambda=1.0*pi)] 0\nM 0\n",
        "distribution", None, 256,
    ),
    (
        "h_then_t_random",
        "H 0\nS[T] 0\nM 0\n",
        "distribution", None, 256,
    ),
    (
        "rx_pi_quarter_random",
        "I[R_X(theta=0.25*pi)] 0\nM 0\n",
        "distribution", None, 256,
    ),
    (
        "rx_ry_rz_combo",
        "I[R_X(theta=0.5*pi)] 0\nI[R_Y(theta=0.5*pi)] 0\nI[R_Z(theta=0.5*pi)] 0\nM 0\n",
        "distribution", None, 256,
    ),
    (
        "u3_arbitrary",
        "I[U3(theta=0.5*pi, phi=0.5*pi, lambda=0.5*pi)] 0\nM 0\n",
        "distribution", None, 256,
    ),
]


def run() -> int:
    paths = core.CorpusPaths.default()
    cat_dir = paths.category_dir("generated/dialect")
    failures: list[str] = []
    written = 0
    for name, src, mode, bitstring, num_shots in DIALECT_FIXTURES:
        if mode == "deterministic":
            meta = core.FixtureMeta(
                name=name,
                category="generated/dialect",
                source=src,
                test_num_shots=1,
            )
            try:
                core.write_deterministic_fixture(meta, paths)
                json_path = cat_dir / f"{name}.expected.json"
                got = json.loads(json_path.read_text())["bitstring"]
                if got != bitstring:
                    failures.append(
                        f"{name}: ppvm produced {got}, asserted {bitstring} — "
                        "either ppvm has a bug or the asserted bitstring is wrong"
                    )
                written += 1
            except Exception as e:
                failures.append(f"{name}: {e}")
        else:
            # No Stim oracle. Run ppvm at seed=0 and record means directly,
            # bypassing the seed-search loop in core.write_distribution_fixture.
            try:
                ppvm_run = core.run_ppvm(src, num_shots=num_shots, seed=0)
                payload = {
                    "mode": "distribution",
                    "num_shots": num_shots,
                    "ppvm_seed": 0,
                    "ppvm_bit_means": ppvm_run.bit_means,
                    "note": "no Stim oracle (ppvm dialect); means recorded directly from ppvm",
                }
                cat_dir.mkdir(parents=True, exist_ok=True)
                (cat_dir / f"{name}.stim").write_text(src)
                (cat_dir / f"{name}.expected.json").write_text(
                    json.dumps(payload, indent=2) + "\n"
                )
                written += 1
            except Exception as e:
                failures.append(f"{name}: {e}")

    print(f"regen-stim dialect: wrote {written} fixtures")
    if failures:
        print("regen-stim dialect: failures:")
        for f in failures:
            print(f"  {f}")
        return 1
    return 0
