// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Differential parity: old `stim_parser` vs new `stim_parser_2`.
//! Asserts (1) same accept/reject and (2) byte-identical canonical print,
//! over a hand corpus and proptest-generated programs. This is the
//! no-regression gate before the consumer swap.

use proptest::prelude::*;

// ---------------------------------------------------------------------------
// Parity-checking helpers
// ---------------------------------------------------------------------------

fn assert_vanilla_parity(src: &str) {
    let old = stim_parser::prelude::parse(src);
    let new = stim_parser_2::prelude::parse(src);
    assert_eq!(
        old.is_ok(),
        new.is_ok(),
        "accept/reject mismatch (vanilla) for:\n{src}\nold: {old:?}\nnew: {new:?}"
    );
    if let (Ok(o), Ok(n)) = (old, new) {
        assert_eq!(
            format!("{o}"),
            n.to_stim(),
            "vanilla print mismatch for:\n{src}"
        );
    }
}

fn assert_extended_parity(src: &str) {
    let old = stim_parser::extended::parse_extended(src);
    let new = stim_parser_2::prelude::parse_extended(src);
    assert_eq!(
        old.is_ok(),
        new.is_ok(),
        "accept/reject mismatch (extended) for:\n{src}\nold: {old:?}\nnew: {new:?}"
    );
    if let (Ok(o), Ok(n)) = (old, new) {
        assert_eq!(
            format!("{o}"),
            n.to_stim(),
            "extended print mismatch for:\n{src}"
        );
    }
}

// ---------------------------------------------------------------------------
// Hand corpus
// ---------------------------------------------------------------------------

/// Each entry is `(label, src)`. The corpus covers the full feature surface
/// so both the accept-side print check and the reject-side accept/reject check
/// get exercised.
const CORPUS: &[(&str, &str)] = &[
    // ---- Canonical Stim circuits ----
    ("empty", ""),
    ("bell_pair", "H 0\nCX 0 1\nM 0 1\n"),
    ("ghz4", "H 0\nCX 0 1\nCX 1 2\nCX 2 3\nM 0 1 2 3\n"),
    (
        "rep_code_d3_r3",
        "R 0 1 2 3 4\nREPEAT 3 {\n    CX 0 3\n    CX 1 3\n    CX 1 4\n    CX 2 4\n    MR 3 4\n    DETECTOR rec[-1]\n    DETECTOR rec[-2]\n    TICK\n}\nM 0 1 2\nOBSERVABLE_INCLUDE(0) rec[-1]\n",
    ),
    // ---- REPEAT (flat and nested) ----
    ("repeat_flat", "REPEAT 5 {\n    X 0\n}\nM 0\n"),
    (
        "nested_repeat",
        "REPEAT 2 {\n    REPEAT 3 {\n        H 0\n        M 0\n    }\n}\n",
    ),
    // ---- Noise instructions ----
    ("depolarize1", "DEPOLARIZE1(0.05) 0 1\n"),
    ("depolarize2", "DEPOLARIZE2(0.01) 0 1\n"),
    ("pauli_channel_1", "PAULI_CHANNEL_1(0.01, 0.02, 0.03) 0\n"),
    ("x_error", "X_ERROR(0.1) 0 1 2\n"),
    ("y_error", "Y_ERROR(0.05) 0\n"),
    ("z_error", "Z_ERROR(0.05) 0\n"),
    // ---- Measure variants ----
    ("m_noise", "M(0.001) 0 1 2\n"),
    ("mz", "MZ 0 1\n"),
    ("mz_noise", "MZ(0.01) 0\n"),
    ("mr", "MR 0 1\n"),
    ("mr_noise", "MR(0.001) 1\n"),
    ("mx_unsupported", "MX 0 1\n"),
    ("my_unsupported", "MY 0\n"),
    ("mrx_unsupported", "MRX 0\n"),
    ("mry_unsupported", "MRY 0\n"),
    ("mxx_unsupported", "MXX 0 1\n"),
    ("myy_unsupported", "MYY 0 1\n"),
    ("mzz_unsupported", "MZZ 0 1\n"),
    // ---- MPAD ----
    ("mpad_bare", "MPAD 0 1 0\n"),
    ("mpad_noisy", "MPAD(0.1) 1 1 0\n"),
    // ---- Annotations (including rec[-k] targets) ----
    ("tick", "TICK\n"),
    ("detector_bare", "DETECTOR\n"),
    ("detector_rec", "DETECTOR rec[-1]\n"),
    ("observable_include", "OBSERVABLE_INCLUDE(0) rec[-1]\n"),
    ("qubit_coords", "QUBIT_COORDS(0, 0) 0\n"),
    ("qubit_coords_3d", "QUBIT_COORDS(1.5, 2.5, 0) 3\n"),
    ("shift_coords", "SHIFT_COORDS(0, 1) \n"),
    // ---- rec[-k] feed-forward on two-qubit gates ----
    ("cx_rec_ff", "M 0\nCX rec[-1] 1\n"),
    ("cy_rec_ff", "M 0\nCY rec[-1] 1\n"),
    ("cz_rec_ff", "M 0\nCZ rec[-1] 1\n"),
    ("cnot_rec_ff", "M 0\nCNOT rec[-1] 1\n"),
    // ---- MPP single and multi-factor products ----
    ("mpp_single_x", "MPP X0\n"),
    ("mpp_single_z", "MPP Z3\n"),
    ("mpp_multi_factor", "MPP X0*Y1*Z2\n"),
    ("mpp_two_products", "MPP X0*Y1*Z2 Z3*Z4\n"),
    ("mpp_many", "MPP X0 Z1 X0*Z1\n"),
    // ---- Extended sugar (T, T_DAG) ----
    ("s_t_sugar", "S[T] 0 1\n"),
    ("s_dag_t_sugar", "S_DAG[T] 2\n"),
    ("native_t", "T 0\n"),
    ("native_t_dag", "T_DAG 0\n"),
    // ---- Extended sugar: I[R_X/Y/Z] ----
    ("rot_x", "I[R_X(theta=0.5)] 0\n"),
    ("rot_y", "I[R_Y(theta=1.25)] 1\n"),
    ("rot_z", "I[R_Z(theta=-0.5)] 2\n"),
    // ---- Extended sugar: I[U3] ----
    ("u3", "I[U3(theta=0.5, phi=1.0, lambda=1.5)] 0\n"),
    // ---- Extended sugar: I_ERROR[loss] / I_ERROR[correlated_loss] ----
    ("loss", "I_ERROR[loss](0.01) 0 1\n"),
    (
        "correlated_loss_3args",
        "I_ERROR[correlated_loss](0.1, 0.05, 0.05) 0 1 2 3\n",
    ),
    (
        "correlated_loss_1arg",
        "I_ERROR[correlated_loss](0.5) 0 1\n",
    ),
    // ---- SWAP / ISWAP (parseable but unsupported at execution) ----
    ("swap", "SWAP 0 1\n"),
    ("iswap", "ISWAP 0 1 2 3\n"),
    // ---- Mixed extended + vanilla circuit ----
    (
        "mixed_extended",
        "H 0\nS[T] 1\nCX 0 1\nI_ERROR[loss](0.01) 0\nM 0 1\n",
    ),
    (
        "extended_in_repeat",
        "REPEAT 3 {\n    S[T] 0\n    I[R_Z(theta=0.25)] 0\n    M 0\n}\n",
    ),
    // ---- Intentionally-invalid inputs (reject parity) ----
    // unknown instruction
    ("invalid_unknown_instr", "FROBNICATE 0\n"),
    // bad arg count for depolarize1
    ("invalid_bad_arg_count", "DEPOLARIZE1(0.1, 0.2) 0\n"),
    // odd target count for CX (must be pairs)
    ("invalid_odd_cx_targets", "CX 0 1 2\n"),
    // no targets for H (AtLeastOne)
    ("invalid_h_no_targets", "H\n"),
    // MPAD with zero targets
    ("invalid_mpad_no_targets", "MPAD\n"),
    // MPAD with two args
    ("invalid_mpad_two_args", "MPAD(0.1, 0.2) 0\n"),
    // bad MPAD bit value
    ("invalid_mpad_bad_bit", "MPAD 0 2 1\n"),
    // malformed REPEAT (missing closing brace)
    ("invalid_repeat_unclosed", "REPEAT 2 {\n    H 0\n"),
    // unknown instruction in extended context
    ("invalid_extended_unknown", "FROBNICATE 0\n"),
    // T rec[-1] — rec target disallowed on sugar gate
    ("invalid_t_rec", "M 0\nT rec[-1]\n"),
    // I_ERROR with no tag
    ("invalid_i_error_no_tag", "I_ERROR(0.1) 0\n"),
    // I_ERROR with unknown tag
    ("invalid_i_error_unknown_tag", "I_ERROR[bogus](0.1) 0\n"),
    // I[R_X] missing theta
    ("invalid_rot_missing_theta", "I[R_X] 0\n"),
    // S with unknown tag
    ("invalid_s_bad_tag", "S[X] 0\n"),
    // rec[-0] is out of range (k must be >= 1)
    ("invalid_rec_zero", "M 0\nCX rec[-0] 1\n"),
];

#[test]
fn corpus_vanilla_parity() {
    for &(label, src) in CORPUS {
        let _: () = {
            let old = stim_parser::prelude::parse(src);
            let new = stim_parser_2::prelude::parse(src);
            assert_eq!(
                old.is_ok(),
                new.is_ok(),
                "[{label}] accept/reject mismatch (vanilla)\nold: {old:?}\nnew: {new:?}"
            );
            if let (Ok(o), Ok(n)) = (old, new) {
                assert_eq!(
                    format!("{o}"),
                    n.to_stim(),
                    "[{label}] vanilla print mismatch"
                );
            }
        };
    }
}

#[test]
fn corpus_extended_parity() {
    for &(label, src) in CORPUS {
        let _: () = {
            let old = stim_parser::extended::parse_extended(src);
            let new = stim_parser_2::prelude::parse_extended(src);
            assert_eq!(
                old.is_ok(),
                new.is_ok(),
                "[{label}] accept/reject mismatch (extended)\nold: {old:?}\nnew: {new:?}"
            );
            if let (Ok(o), Ok(n)) = (old, new) {
                assert_eq!(
                    format!("{o}"),
                    n.to_stim(),
                    "[{label}] extended print mismatch"
                );
            }
        };
    }
}

// ---------------------------------------------------------------------------
// Proptest fragment-based generator
// (mirrored from proptest_roundtrip.rs; extended with additional fragments
//  to exercise MPP, rec feed-forward, native T, correlated_loss, rotations,
//  and U3 — features missing from the roundtrip generator)
// ---------------------------------------------------------------------------

fn instruction_fragment() -> impl Strategy<Value = String> {
    prop_oneof![
        // Bare Clifford gates
        Just("H 0\n".to_string()),
        Just("X 0\n".to_string()),
        Just("Y 0\n".to_string()),
        Just("Z 0\n".to_string()),
        Just("S 0\n".to_string()),
        Just("S_DAG 1\n".to_string()),
        Just("I 0\n".to_string()),
        Just("CX 0 1\n".to_string()),
        Just("CZ 1 2\n".to_string()),
        Just("CNOT 0 3\n".to_string()),
        // Reset / measure
        Just("R 0 1\n".to_string()),
        Just("M 0\n".to_string()),
        Just("MZ 0 1 2\n".to_string()),
        Just("MR 0\n".to_string()),
        Just("M(0.001) 0\n".to_string()),
        // Tagged sugar
        Just("S[T] 0\n".to_string()),
        Just("S_DAG[T] 1\n".to_string()),
        Just("I[R_X(theta=0.5)] 0\n".to_string()),
        Just("I[R_Y(theta=1.25)] 1\n".to_string()),
        Just("I[R_Z(theta=-0.5)] 2\n".to_string()),
        Just("I[U3(theta=0.5, phi=1.0, lambda=1.5)] 0\n".to_string()),
        // Native T / T_DAG
        Just("T 0\n".to_string()),
        Just("T_DAG 0\n".to_string()),
        // Noise
        Just("DEPOLARIZE1(0.05) 0\n".to_string()),
        Just("DEPOLARIZE2(0.05) 0 1\n".to_string()),
        Just("PAULI_CHANNEL_1(0.01, 0.02, 0.03) 0\n".to_string()),
        Just("X_ERROR(0.1) 0\n".to_string()),
        Just("I_ERROR[loss](0.01) 0\n".to_string()),
        Just("I_ERROR[correlated_loss](0.1, 0.05, 0.05) 0 1\n".to_string()),
        Just("I_ERROR[correlated_loss](0.5) 0 1 2 3\n".to_string()),
        // MPAD
        Just("MPAD 0 1 0\n".to_string()),
        Just("MPAD(0.1) 1 1 0\n".to_string()),
        // Annotations (with and without rec[-1])
        Just("TICK\n".to_string()),
        Just("DETECTOR\n".to_string()),
        Just("DETECTOR rec[-1]\n".to_string()),
        Just("OBSERVABLE_INCLUDE(0) rec[-1]\n".to_string()),
        Just("QUBIT_COORDS(0, 0) 0\n".to_string()),
        // MPP single and multi-factor products
        Just("MPP X0\n".to_string()),
        Just("MPP Z1\n".to_string()),
        Just("MPP X0*Y1*Z2\n".to_string()),
        Just("MPP X0*Y1*Z2 Z3*Z4\n".to_string()),
        // rec[-k] feed-forward on two-qubit gates (always preceded by M 0)
        Just("M 0\nCX rec[-1] 1\n".to_string()),
        Just("M 0\nCY rec[-1] 1\n".to_string()),
        Just("M 0\nCZ rec[-1] 1\n".to_string()),
        Just("M 0\nCNOT rec[-1] 1\n".to_string()),
        // REPEAT
        Just("REPEAT 3 {\n    H 0\n    M 0\n}\n".to_string()),
        Just("REPEAT 2 {\n    REPEAT 3 {\n        X 0\n    }\n}\n".to_string()),
        // Stylistic noise the printer normalizes away
        Just("# leading\n".to_string()),
        Just("\n".to_string()),
        Just("H 0  # trail\n".to_string()),
    ]
}

fn program_source() -> impl Strategy<Value = String> {
    prop::collection::vec(instruction_fragment(), 0..16).prop_map(|frags| frags.concat())
}

proptest! {
    #[test]
    fn proptest_vanilla_parity(src in program_source()) {
        assert_vanilla_parity(&src);
    }

    #[test]
    fn proptest_extended_parity(src in program_source()) {
        assert_extended_parity(&src);
    }
}
