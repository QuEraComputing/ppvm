// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! End-to-end fixture coverage: parse + resolve + run each `.sst` file in
//! this directory via the public `PPVM` API.

use ppvm_vihaco::composite::PPVM;
use ppvm_vihaco::measurements::MeasurementOutcome;

/// Dump a fixture to a `.ssb` file, load it back, and run it. Exercises the
/// full bytecode round-trip through disk: `dump_file` → `load_bytecode_file`.
fn dump_load_run(sst_path: &str, ssb_name: &str) -> PPVM {
    let out = std::env::temp_dir().join(ssb_name);
    let out = out.to_str().expect("utf-8 temp path");
    ppvm_vihaco::dump_file(sst_path, out).unwrap_or_else(|e| panic!("dump {sst_path}: {e:?}"));

    let mut machine = PPVM::default();
    machine
        .load_bytecode_file(out)
        .unwrap_or_else(|e| panic!("load {out}: {e:?}"));
    machine.run().unwrap_or_else(|e| panic!("run {out}: {e:?}"));

    let _ = std::fs::remove_file(out);
    machine
}

#[test]
fn bell_sst_runs_and_records_two_measurements() {
    let mut machine = PPVM::default();
    machine
        .run_file("tests/bell.sst")
        .unwrap_or_else(|e| panic!("run bell.sst: {e:?}"));
    assert_eq!(machine.measurement_record().len(), 2);
}

#[test]
fn hello_circuit_sst_parses_and_runs() {
    let mut machine = PPVM::default();
    machine
        .run_file("tests/hello_circuit.sst")
        .unwrap_or_else(|e| panic!("run hello_circuit.sst: {e:?}"));
    // hello_circuit.sst applies H + CNOT + RX(0.1); no measurements.
    assert_eq!(machine.measurement_record().len(), 0);
}

#[test]
fn rotxy_sst_runs_and_flips_qubit() {
    // `rotxy.sst` applies R(axis_angle=π/2, θ=π) = RY(π) to q0, deterministically
    // sending |0> → |1>, then measures it. Exercises the `gate r` path end to
    // end: parse → resolve (pop θ, axis_angle, qubit) → execute via `tab.r`.
    let machine =
        ppvm_vihaco::run_file("tests/rotxy.sst").unwrap_or_else(|e| panic!("run rotxy.sst: {e:?}"));
    let record = machine.measurement_record();
    assert_eq!(record.len(), 1, "expected exactly one measurement");
    assert_eq!(
        record[0].as_slice(),
        &[MeasurementOutcome::One],
        "R(π/2, π) = RY(π) must flip q0 to 1"
    );
}

#[test]
fn dumped_rotxy_runs_and_flips_qubit() {
    // Same program, but through the bytecode round-trip: confirms the `R`
    // instruction survives dump → `.ssb` → load → execute.
    let machine = dump_load_run("tests/rotxy.sst", "ppvm_dump_rotxy.ssb");
    let record = machine.measurement_record();
    assert_eq!(record.len(), 1);
    assert_eq!(record[0].as_slice(), &[MeasurementOutcome::One]);
}

#[test]
fn run_file_via_library_helper() {
    let machine =
        ppvm_vihaco::run_file("tests/bell.sst").unwrap_or_else(|e| panic!("run bell.sst: {e:?}"));
    assert_eq!(machine.measurement_record().len(), 2);
}

#[test]
fn function_call_jumps_into_callee_body() {
    // `function_call.sst` has main `call` into `@run_circuit`, which puts q1
    // in |+>, measures it, and `halt`s. Verifies CallPatch resolves the
    // symbolic target and op_call actually transfers control there.
    let machine = ppvm_vihaco::run_file("tests/function_call.sst")
        .unwrap_or_else(|e| panic!("run function_call.sst: {e:?}"));
    let record = machine.measurement_record();
    assert_eq!(record.len(), 1, "expected exactly one measurement");
    assert_eq!(record[0].len(), 1);
    assert!(
        record[0][0] != MeasurementOutcome::Lost,
        "measurement should not be lost"
    );
}

#[test]
fn function_call_returns() {
    let machine = ppvm_vihaco::run_file("tests/function_call_ret.sst")
        .unwrap_or_else(|e| panic!("run function_call.sst: {e:?}"));
    let record = machine.measurement_record();
    assert_eq!(record.len(), 1, "expected exactly one measurement");
    assert_eq!(record[0].len(), 1);
    assert!(
        record[0][0] != MeasurementOutcome::Lost,
        "measurement should not be lost"
    );
}

#[test]
fn branch_on_outcome_deterministic_x_path() {
    // `branch_on_outcome_x.sst` applies X to q0 instead of H, so the outcome
    // is deterministically 1. The cond_br must therefore take the @one path,
    // which flips q1 before measuring it, yielding m1 = 1 as well.
    let machine = ppvm_vihaco::run_file("tests/branch_on_outcome_x.sst")
        .unwrap_or_else(|e| panic!("run branch_on_outcome_x.sst: {e:?}"));
    let record = machine.measurement_record();
    assert_eq!(record.len(), 2, "expected exactly two measurements");
    assert_eq!(
        record[0].as_slice(),
        &[MeasurementOutcome::One],
        "X-prepared q0 must measure 1"
    );
    assert_eq!(
        record[1].as_slice(),
        &[MeasurementOutcome::One],
        "branch must have flipped q1"
    );
}

#[test]
fn branch_on_outcome_statistics_balanced_and_invariant_holds() {
    // `branch_on_outcome.sst` puts q0 in |+>, so its measurement is a fair
    // coin. The branch then flips q1 iff the outcome was 1, making m1 == m0
    // an invariant on every shot. Run many shots and check both properties.
    const SHOTS: usize = 400;
    let mut ones = 0usize;
    for _ in 0..SHOTS {
        let machine = ppvm_vihaco::run_file("tests/branch_on_outcome.sst")
            .unwrap_or_else(|e| panic!("run branch_on_outcome.sst: {e:?}"));
        let record = machine.measurement_record();
        assert_eq!(record.len(), 2);
        let m0 = record[0][0];
        let m1 = record[1][0];
        assert_eq!(m0, m1, "branch must steer q1 to match q0 on every shot");
        assert!(
            m0 != MeasurementOutcome::Lost,
            "measurement should not be lost"
        );
        if m0 == MeasurementOutcome::One {
            ones += 1;
        }
    }
    // Fair coin with SHOTS=400: mean=200, stddev=10. ±6σ window catches a
    // truly broken RNG without flaking on a healthy one.
    let lo = SHOTS / 2 - 60;
    let hi = SHOTS / 2 + 60;
    assert!(
        (lo..=hi).contains(&ones),
        "expected {lo}..={hi} ones out of {SHOTS}, got {ones}"
    );
}

#[test]
fn function_call_branch_on_both_returned_values() {
    // `function_call_branch_both.sst`: helper returns the tri-state outcome
    // (0/1/Lost) via `ret 1`. Main first branches on is_lost, then on the
    // 0/1 outcome, steering q1 to |1> on the lost path and on the
    // kept-outcome=1 path, leaving q1 in |0> only on the kept-outcome=0
    // path. With loss prob 0.5 and a |+> prep:
    //   P(m1 = 1) = P(lost) + P(kept) · P(outcome = 1 | kept)
    //             = 0.5     + 0.5 · 0.5  = 0.75
    //   P(m0 = lost) = 0.5
    const SHOTS: usize = 400;
    let mut q0_lost = 0usize;
    let mut q1_ones = 0usize;
    for _ in 0..SHOTS {
        let machine = ppvm_vihaco::run_file("tests/function_call_branch_both.sst")
            .unwrap_or_else(|e| panic!("run function_call_branch_both.sst: {e:?}"));
        let record = machine.measurement_record();
        assert_eq!(record.len(), 2, "expected exactly two measurements");
        assert_eq!(record[0].len(), 1);
        assert_eq!(record[1].len(), 1);
        if record[0][0] == MeasurementOutcome::Lost {
            q0_lost += 1;
        }
        if record[1][0] == MeasurementOutcome::One {
            q1_ones += 1;
        }
    }
    // P(lost) = 0.5, SHOTS=400: mean=200, stddev=10. ±6σ window.
    assert!(
        (140..=260).contains(&q0_lost),
        "expected ~200 lost shots, got {q0_lost}"
    );
    // P(m1=1) = 0.75, SHOTS=400: mean=300, stddev≈8.66. ±6σ → ~248..352.
    assert!(
        (240..=360).contains(&q1_ones),
        "expected ~300 q1=true shots, got {q1_ones}"
    );
}

// ─── Task 12: PauliSum / LossyPauliSum end-to-end via .sst source ────────

#[test]
fn paulisum_bell_zz_trace_through_sst() {
    // Bell-state ⟨ZZ⟩ via PauliSum. Textbook circuit H(0); CNOT(0,1) is
    // emitted reversed for Heisenberg propagation: `gate cnot; gate h`.
    // Conjugating ZZ by CNOT(0,1) gives Z_1 (= IZ); H on q0 leaves IZ
    // untouched. Tracing against |00> matches IZ (pattern `Z?*`) and
    // returns +1.0 — matching ⟨Φ+|ZZ|Φ+⟩ = 1.
    let machine = ppvm_vihaco::run_file("tests/paulisum_bell_trace.sst")
        .unwrap_or_else(|e| panic!("run paulisum_bell_trace.sst: {e:?}"));
    let trace = machine.trace_record();
    assert_eq!(trace.len(), 1, "expected one trace emission");
    assert!(
        (trace[0] - 1.0).abs() < 1e-12,
        "expected ⟨ZZ⟩ = 1.0, got {}",
        trace[0]
    );
}

#[test]
fn paulisum_multi_term_observable_trace_through_sst() {
    // Phase 3 wiring: a 2-term observable `1.0*ZZ+0.5*XX` parses on the
    // header, seeds both terms, and an end-of-program trace picks up their
    // contributions. With no gates applied, tracing `[XZ]?*` matches both
    // ZZ and XX directly → coef sum = 1.0 + 0.5 = 1.5.
    let machine = ppvm_vihaco::run_file("tests/paulisum_multi_term_trace.sst")
        .unwrap_or_else(|e| panic!("run paulisum_multi_term_trace.sst: {e:?}"));
    let trace = machine.trace_record();
    assert_eq!(trace.len(), 1);
    assert!(
        (trace[0] - 1.5).abs() < 1e-12,
        "expected 1.5, got {}",
        trace[0]
    );
}

#[test]
fn paulisum_trotter_matches_pure_rust_reference() {
    // Two Trotter layers of RXX(0.1) + RZZ(0.05), interleaved with explicit
    // `gate truncate`. The .sst-driven path should agree bit-for-bit with a
    // pure Rust PauliSum running the same gates: `indexmap::ByteFxHashF64`
    // gives deterministic iteration order (Decision 7), so truncation order
    // and float accumulation are stable across both paths.

    let machine = ppvm_vihaco::run_file("tests/paulisum_trotter_truncate.sst")
        .unwrap_or_else(|e| panic!("run paulisum_trotter_truncate.sst: {e:?}"));
    let through_sst = machine.trace_record();
    assert_eq!(through_sst.len(), 1, "expected one trace emission");

    // Pure Rust reference: same N=8 / strategy / gate order as the PauliSum
    // Bits64 bucket in `ppvm_vihaco::component`.
    use ppvm_runtime::config::indexmap::ByteFxHashF64;
    use ppvm_runtime::prelude::*;
    use ppvm_runtime::strategy::{CoefficientThreshold, CombinedStrategy, MaxPauliWeight};
    type RefConfig = ByteFxHashF64<8, CombinedStrategy<CoefficientThreshold, MaxPauliWeight>>;

    let mut state: PauliSum<RefConfig> = PauliSum::builder()
        .n_qubits(2)
        .strategy(CombinedStrategy(
            CoefficientThreshold(1e-6),
            MaxPauliWeight(usize::MAX),
        ))
        .build();
    state += ("ZZ", 1.0);
    state.rxx(0, 1, 0.1);
    state.rzz(0, 1, 0.05);
    state.truncate();
    state.rxx(0, 1, 0.1);
    state.rzz(0, 1, 0.05);
    state.truncate();
    let pat = PauliPattern::parse("Z?*").expect("parse pattern");
    let reference = state.trace(&pat);

    assert_eq!(
        through_sst[0], reference,
        ".sst-driven trace must match pure Rust reference bit-for-bit"
    );
}

#[test]
fn paulisum_measure_returns_unsupported_error() {
    // Per Decision 11, Measure on PauliSum hits the dispatch fallback with a
    // clear "not supported on the PauliSum backend" error.
    let mut machine = PPVM::default();
    machine
        .load_file("tests/paulisum_measure_error.sst")
        .unwrap_or_else(|e| panic!("load paulisum_measure_error.sst: {e:?}"));
    let err = machine.run().unwrap_err();
    let msg = err.to_string();
    assert!(
        msg.contains("not supported on the PauliSum backend"),
        "expected PauliSum-rejection error, got: {msg}"
    );
}

// ─── Task 16: Tableau-side Trace, cross-backend agreement ────────────────

#[test]
fn tableau_bell_zz_trace_through_sst() {
    // Bell-state ⟨ZZ⟩ via the Tableau backend's `Trace` instruction:
    // forward H(0); CNOT(0, 1) leaves |ψ⟩ = |Φ+⟩, and `Z0Z1` matches exactly
    // the Pauli word ZZ, so `tab.trace(&pat) = ⟨Φ+|ZZ|Φ+⟩ = 1`.
    let machine = ppvm_vihaco::run_file("tests/tableau_bell_trace.sst")
        .unwrap_or_else(|e| panic!("run tableau_bell_trace.sst: {e:?}"));
    let trace = machine.trace_record();
    assert_eq!(trace.len(), 1, "expected one trace emission");
    assert!(
        (trace[0] - 1.0).abs() < 1e-12,
        "expected ⟨ZZ⟩ = 1.0, got {}",
        trace[0]
    );
}

fn assert_cross_backend_agreement(tableau_sst: &str, paulisum_sst: &str) {
    let tab =
        ppvm_vihaco::run_file(tableau_sst).unwrap_or_else(|e| panic!("run {tableau_sst}: {e:?}"));
    let ps =
        ppvm_vihaco::run_file(paulisum_sst).unwrap_or_else(|e| panic!("run {paulisum_sst}: {e:?}"));
    let tab_v = tab.trace_record();
    let ps_v = ps.trace_record();
    assert_eq!(tab_v.len(), 1, "{tableau_sst}: expected one trace emission");
    assert_eq!(ps_v.len(), 1, "{paulisum_sst}: expected one trace emission");
    assert!(
        (tab_v[0] - ps_v[0]).abs() < 1e-12,
        "tableau {} vs PauliSum {} (|Δ|={}) for {tableau_sst} ↔ {paulisum_sst}",
        tab_v[0],
        ps_v[0],
        (tab_v[0] - ps_v[0]).abs()
    );
}

#[test]
fn tableau_and_paulisum_agree_on_bell_zz_trace() {
    // ⟨Φ+|ZZ|Φ+⟩ = 1. Tableau forward-evolves |0…0⟩ and matches `Z0Z1`;
    // PauliSum Heisenberg-propagates ZZ backward through the reversed
    // circuit, then sums Z/I-only coefficients.
    assert_cross_backend_agreement(
        "tests/tableau_bell_trace.sst",
        "tests/paulisum_bell_trace.sst",
    );
}

#[test]
fn tableau_and_paulisum_agree_on_ghz_xxx_trace() {
    // ⟨GHZ|XXX|GHZ⟩ = 1. Exercises a 3-qubit Clifford chain (H + two
    // CNOTs) and the non-trivial Heisenberg evolution XXX → ZII through
    // the reversed circuit on the PauliSum side.
    assert_cross_backend_agreement(
        "tests/tableau_ghz_xxx_trace.sst",
        "tests/paulisum_ghz_xxx_trace.sst",
    );
}

#[test]
fn tableau_and_paulisum_agree_on_ry_z_trace() {
    // ⟨RY(θ)|0⟩|Z|RY(θ)|0⟩⟩ = cos(θ). Exercises a non-Clifford rotation:
    // the tableau opens a branch via `tab.ry`, while the PauliSum applies
    // the Heisenberg dual RY(θ)†·Z·RY(θ) = cos(θ)·Z + sin(θ)·X in one step.
    assert_cross_backend_agreement(
        "tests/tableau_ry_z_trace.sst",
        "tests/paulisum_ry_z_trace.sst",
    );
}

// ─── Auto-detect via load_file: route by content, not extension ───────────

#[test]
fn is_bytecode_distinguishes_ssb_from_sst() {
    let sst = std::fs::read("tests/bell.sst").expect("read bell.sst");
    assert!(
        !ppvm_vihaco::bytecode::is_bytecode(&sst),
        ".sst source must not be detected as bytecode"
    );

    let ssb = ppvm_vihaco::bytecode::compile_to_bytes(
        &String::from_utf8(sst).expect("bell.sst is utf-8"),
    )
    .expect("compile bell.sst");
    assert!(
        ppvm_vihaco::bytecode::is_bytecode(&ssb),
        "dumped .ssb must be detected as bytecode"
    );

    // Inputs shorter than the 4-byte magic are never bytecode. Note "PPVM" as
    // text also fails: the magic is a little-endian u32, so its on-disk bytes
    // are "MVPP", not "PPVM".
    assert!(!ppvm_vihaco::bytecode::is_bytecode(b"PPV"));
    assert!(!ppvm_vihaco::bytecode::is_bytecode(b""));
    assert!(!ppvm_vihaco::bytecode::is_bytecode(b"PPVM"));
}

#[test]
fn load_file_auto_detects_bytecode_and_text() {
    // Use the deterministic X-prepared fixture: q0 measures 1, the branch
    // flips q1, so both routes must yield exactly [1], [1]. Any divergence —
    // or a binary file mis-parsed as text — fails loudly.
    let from_text = ppvm_vihaco::run_file("tests/branch_on_outcome_x.sst")
        .unwrap_or_else(|e| panic!("run .sst via load_file: {e:?}"));

    // Dump the same fixture to a `.ssb` and run *that file* through the same
    // run_file entry point. If load_file didn't sniff the magic it would try
    // to parse the binary as text and error.
    let out = std::env::temp_dir().join("ppvm_autodetect_branch_x.ssb");
    let out = out.to_str().expect("utf-8 temp path");
    ppvm_vihaco::dump_file("tests/branch_on_outcome_x.sst", out)
        .unwrap_or_else(|e| panic!("dump: {e:?}"));
    let from_binary = ppvm_vihaco::run_file(out).unwrap_or_else(|e| panic!("run .ssb: {e:?}"));
    let _ = std::fs::remove_file(out);

    for (label, machine) in [("text", &from_text), ("binary", &from_binary)] {
        let record = machine.measurement_record();
        assert_eq!(record.len(), 2, "{label}: expected two measurements");
        assert_eq!(
            record[0].as_slice(),
            &[MeasurementOutcome::One],
            "{label}: X-prepared q0 must measure 1"
        );
        assert_eq!(
            record[1].as_slice(),
            &[MeasurementOutcome::One],
            "{label}: branch must flip q1"
        );
    }
}

// ─── Bytecode round-trip: dump → load → execute each fixture ──────────────

#[test]
fn dumped_bell_records_two_measurements() {
    let machine = dump_load_run("tests/bell.sst", "ppvm_dump_bell.ssb");
    assert_eq!(machine.measurement_record().len(), 2);
}

#[test]
fn dumped_hello_circuit_runs_with_no_measurements() {
    let machine = dump_load_run("tests/hello_circuit.sst", "ppvm_dump_hello_circuit.ssb");
    assert_eq!(machine.measurement_record().len(), 0);
}

#[test]
fn dumped_function_call_executes_callee() {
    let machine = dump_load_run("tests/function_call.sst", "ppvm_dump_function_call.ssb");
    let record = machine.measurement_record();
    assert_eq!(record.len(), 1);
    assert_eq!(record[0].len(), 1);
    assert!(record[0][0] != MeasurementOutcome::Lost);
}

#[test]
fn dumped_function_call_ret_executes() {
    let machine = dump_load_run(
        "tests/function_call_ret.sst",
        "ppvm_dump_function_call_ret.ssb",
    );
    let record = machine.measurement_record();
    assert_eq!(record.len(), 1);
    assert_eq!(record[0].len(), 1);
}

#[test]
fn dumped_branch_on_outcome_x_is_deterministic() {
    // X-prepared q0 measures 1, so the branch flips q1 → both outcomes are 1.
    // Confirms branch targets survive the dump/load round-trip.
    let machine = dump_load_run("tests/branch_on_outcome_x.sst", "ppvm_dump_branch_x.ssb");
    let record = machine.measurement_record();
    assert_eq!(record.len(), 2);
    assert_eq!(record[0].as_slice(), &[MeasurementOutcome::One]);
    assert_eq!(record[1].as_slice(), &[MeasurementOutcome::One]);
}

#[test]
fn dumped_branch_on_outcome_preserves_invariant() {
    // q0 in |+> is a fair coin, but the branch steers q1 to match q0 every
    // shot — that invariant must hold after a round-trip.
    let machine = dump_load_run("tests/branch_on_outcome.sst", "ppvm_dump_branch.ssb");
    let record = machine.measurement_record();
    assert_eq!(record.len(), 2);
    assert_eq!(record[0][0], record[1][0]);
}

#[test]
fn dumped_function_call_branch_both_runs() {
    let machine = dump_load_run(
        "tests/function_call_branch_both.sst",
        "ppvm_dump_function_call_branch_both.ssb",
    );
    let record = machine.measurement_record();
    assert_eq!(record.len(), 2);
    assert_eq!(record[0].len(), 1);
    assert_eq!(record[1].len(), 1);
}
