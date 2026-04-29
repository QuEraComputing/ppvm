"""Shared library for the regen-stim CLI.

Layout: every subcommand emits one or more (stim_source, category, name) triples
to ``write_fixture``, which:

  1. Runs Stim on the source to compute reference per-bit means.
  2. Searches ppvm seeds in [0, max_seed) until ppvm's empirical means at
     ``num_shots`` are within tolerance_sigma * sqrt(p*(1-p)/N) of Stim's.
  3. Writes ``<name>.stim`` and ``<name>.expected.json`` into
     ``<corpus_root>/<category>/``.

Distribution mode is the default. ``write_unsupported_fixture`` and
``write_deterministic_fixture`` handle the two narrower cases.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path

import stim

DEFAULT_STIM_SHOTS = 10_000
DEFAULT_MAX_SEED = 32
DEFAULT_TOLERANCE_SIGMA = 5.0


@dataclass(frozen=True)
class CorpusPaths:
    """Where fixtures live on disk. Defaults to the corpus root next to this tool."""

    root: Path

    @classmethod
    def default(cls) -> "CorpusPaths":
        # regen-stim/src/regen_stim/core.py → ../../../data
        here = Path(__file__).resolve()
        root = here.parents[3] / "data"
        return cls(root=root)

    def category_dir(self, category: str) -> Path:
        return self.root / category


@dataclass
class StimReference:
    """Per-bit means computed by running Stim."""

    bit_means: list[float]
    num_shots: int
    seed: int
    stim_version: str


def run_stim(source: str, num_shots: int = DEFAULT_STIM_SHOTS, seed: int = 0) -> StimReference:
    """Run Stim and return per-bit empirical means."""
    circuit = stim.Circuit(source)
    sampler = circuit.compile_sampler(seed=seed)
    samples = sampler.sample(shots=num_shots)
    n_meas = samples.shape[1] if samples.ndim == 2 else 0
    means = [float(samples[:, i].mean()) if n_meas > 0 else 0.0 for i in range(n_meas)]
    return StimReference(
        bit_means=means,
        num_shots=num_shots,
        seed=seed,
        stim_version=stim.__version__,
    )


@dataclass
class PpvmRun:
    """Per-bit means computed by running ppvm."""

    bit_means: list[float]
    num_shots: int
    seed: int


def run_ppvm(source: str, num_shots: int, seed: int) -> PpvmRun:
    """Run ppvm and return per-bit empirical means.

    Counts qubits by inspecting Stim's parsed circuit (Stim has the parser, ppvm
    expects an n_qubits arg) and uses the same value when constructing tableaux.
    """
    import ppvm

    n_qubits = max(64, _max_qubit_in_source(source) + 1)
    prog = ppvm.StimProgram.parse(source)
    shots = ppvm.sample_stim(prog, n_qubits=n_qubits, num_shots=num_shots, seed=seed)
    if not shots:
        return PpvmRun(bit_means=[], num_shots=num_shots, seed=seed)
    n_meas = len(shots[0])
    means: list[float] = []
    for i in range(n_meas):
        s = 0
        for shot in shots:
            v = shot[i]
            if v is None:
                raise ValueError(
                    f"ppvm returned None (loss) for bit {i}; corpus excludes loss "
                    "(see spec Non-Goals)"
                )
            s += 1 if v else 0
        means.append(s / num_shots)
    return PpvmRun(bit_means=means, num_shots=num_shots, seed=seed)


def _max_qubit_in_source(source: str) -> int:
    """Walk Stim's parsed circuit to find the highest qubit index referenced."""
    circuit = stim.Circuit(source)
    max_q = -1
    for inst in circuit.flattened():
        for t in inst.targets_copy():
            if t.is_qubit_target:
                max_q = max(max_q, t.qubit_value)
    return max_q


def per_bit_sigma(stim_means: list[float], num_shots: int) -> list[float]:
    """Worst-case binomial sigma per bit at the test-time num_shots."""
    out: list[float] = []
    for p in stim_means:
        if 0.0 < p < 1.0:
            out.append(math.sqrt(p * (1.0 - p) / num_shots))
        else:
            # p ∈ {0, 1}: empirical mean is always p in the noiseless limit, but
            # any noise term we can't see in stim_means could raise the variance
            # by one order. Use the worst-case bound 1/N as a guardrail.
            out.append(math.sqrt(1.0 / num_shots))
    return out


def within_tolerance(
    ppvm_means: list[float],
    stim_means: list[float],
    test_num_shots: int,
    tolerance_sigma: float = DEFAULT_TOLERANCE_SIGMA,
) -> bool:
    if len(ppvm_means) != len(stim_means):
        return False
    sigmas = per_bit_sigma(stim_means, test_num_shots)
    for p, s, sig in zip(ppvm_means, stim_means, sigmas, strict=True):
        if abs(p - s) > tolerance_sigma * sig:
            return False
    return True


@dataclass
class FixtureMeta:
    name: str
    category: str
    source: str
    test_num_shots: int
    stim_num_shots: int = DEFAULT_STIM_SHOTS
    stim_seed: int = 0
    max_ppvm_seed: int = DEFAULT_MAX_SEED
    tolerance_sigma: float = DEFAULT_TOLERANCE_SIGMA
    extra_metadata: dict = field(default_factory=dict)


def write_distribution_fixture(meta: FixtureMeta, paths: CorpusPaths) -> Path:
    """Generate a distribution-mode fixture; raises on irreconcilable cross-check."""
    ref = run_stim(meta.source, num_shots=meta.stim_num_shots, seed=meta.stim_seed)
    chosen_seed: int | None = None
    chosen_means: list[float] = []
    for seed in range(meta.max_ppvm_seed):
        ppvm_run = run_ppvm(meta.source, num_shots=meta.test_num_shots, seed=seed)
        if within_tolerance(
            ppvm_run.bit_means, ref.bit_means, meta.test_num_shots, meta.tolerance_sigma
        ):
            chosen_seed = seed
            chosen_means = ppvm_run.bit_means
            break
    if chosen_seed is None:
        raise RuntimeError(
            f"{meta.category}/{meta.name}: no ppvm seed in [0, {meta.max_ppvm_seed}) "
            f"agrees with Stim within {meta.tolerance_sigma} sigma at "
            f"num_shots={meta.test_num_shots}. Stim means: {ref.bit_means}. "
            "This is a real correctness divergence — do not commit."
        )

    payload = {
        "mode": "distribution",
        "num_shots": meta.test_num_shots,
        "ppvm_seed": chosen_seed,
        "ppvm_bit_means": chosen_means,
        "stim_seed": meta.stim_seed,
        "stim_num_shots": meta.stim_num_shots,
        "stim_bit_means": ref.bit_means,
        "tolerance_sigma_at_regen": meta.tolerance_sigma,
        "stim_version": ref.stim_version,
        **meta.extra_metadata,
    }
    return _emit(meta, paths, payload)


def write_deterministic_fixture(meta: FixtureMeta, paths: CorpusPaths) -> Path:
    """One-shot ppvm run; record the bitstring."""
    import ppvm

    n_qubits = max(64, _max_qubit_in_source(meta.source) + 1)
    prog = ppvm.StimProgram.parse(meta.source)
    shots = ppvm.sample_stim(prog, n_qubits=n_qubits, num_shots=1, seed=0)
    bitstring = [bool(b) for b in shots[0]] if shots else []

    payload = {
        "mode": "deterministic",
        "ppvm_seed": 0,
        "bitstring": bitstring,
        **meta.extra_metadata,
    }
    return _emit(meta, paths, payload)


def write_unsupported_fixture(
    meta: FixtureMeta,
    paths: CorpusPaths,
    awaiting_phase2_instruction: str,
) -> Path:
    """Pre-record Stim reference for the day phase-2 lifts the restriction."""
    ref = run_stim(meta.source, num_shots=meta.stim_num_shots, seed=meta.stim_seed)
    payload = {
        "mode": "unsupported",
        "awaiting_phase2_instruction": awaiting_phase2_instruction,
        "stim_seed": meta.stim_seed,
        "stim_num_shots": meta.stim_num_shots,
        "stim_bit_means": ref.bit_means,
        "tolerance_sigma_at_regen": meta.tolerance_sigma,
        "stim_version": ref.stim_version,
        **meta.extra_metadata,
    }
    return _emit(meta, paths, payload)


def _emit(meta: FixtureMeta, paths: CorpusPaths, payload: dict) -> Path:
    cat_dir = paths.category_dir(meta.category)
    cat_dir.mkdir(parents=True, exist_ok=True)
    stim_path = cat_dir / f"{meta.name}.stim"
    json_path = cat_dir / f"{meta.name}.expected.json"
    stim_path.write_text(meta.source)
    json_path.write_text(json.dumps(payload, indent=2) + "\n")
    return json_path
