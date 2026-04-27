use ppvm_stim::{execute, normalize, parse};
use ppvm_runtime::config::indexmap::ByteFxHashF64;
use ppvm_tableau::prelude::*;

type Tab = GeneralizedTableau<ByteFxHashF64<1>, usize>;

fn run(src: &str, n_qubits: usize) -> (Vec<Option<bool>>, Tab) {
    let prog = parse(src).expect("parse");
    let tprog = normalize::to_tableau(&prog).expect("normalize");
    let mut tab: Tab = GeneralizedTableau::new(n_qubits, 1e-10);
    let results = execute(&tprog, &mut tab).expect("execute");
    (results, tab)
}

#[test]
fn x_then_measure_returns_one() {
    let (results, _) = run("X 0\nM 0", 1);
    assert_eq!(results, vec![Some(true)]);
}

#[test]
fn fresh_state_measures_zero() {
    let (results, _) = run("M 0", 1);
    assert_eq!(results, vec![Some(false)]);
}

#[test]
fn h_h_returns_zero() {
    let (results, _) = run("H 0\nH 0\nM 0", 1);
    assert_eq!(results, vec![Some(false)]);
}

#[test]
fn rx_pi_flips_qubit() {
    let (results, _) = run("I[R_X(theta=1.0*pi)] 0\nM 0", 1);
    assert_eq!(results, vec![Some(true)]);
}

#[test]
fn u3_pi_flip_via_y_axis() {
    let (results, _) = run("I[U3(theta=1.0*pi, phi=0.0, lambda=0.0)] 0\nM 0", 1);
    assert_eq!(results, vec![Some(true)]);
}

#[test]
fn t_gate_via_s_t_tag_no_op_on_zero() {
    let (results, _) = run("S[T] 0\nM 0", 1);
    assert_eq!(results, vec![Some(false)]);
}

#[test]
fn cnot_creates_bell_correlation() {
    for _ in 0..32 {
        let (results, _) = run("H 0\nCX 0 1\nM 0 1", 2);
        assert_eq!(results.len(), 2);
        assert_eq!(results[0], results[1]);
    }
}

#[test]
fn cnot_alias_equivalents() {
    let (cx, _) = run("X 0\nCX 0 1\nM 1", 2);
    let (cnot, _) = run("X 0\nCNOT 0 1\nM 1", 2);
    let (zcx, _) = run("X 0\nZCX 0 1\nM 1", 2);
    assert_eq!(cx, vec![Some(true)]);
    assert_eq!(cnot, vec![Some(true)]);
    assert_eq!(zcx, vec![Some(true)]);
}

#[test]
fn mr_resets_qubit_after_measure() {
    let (results, _) = run("X 0\nMR 0\nM 0", 1);
    assert_eq!(results, vec![Some(true), Some(false)]);
}

#[test]
fn loss_channel_with_p1_marks_qubit_lost() {
    let prog = parse("I_ERROR[loss](1.0) 0").unwrap();
    let tprog = normalize::to_tableau(&prog).unwrap();
    let mut tab: Tab = GeneralizedTableau::new(1, 1e-10);
    execute(&tprog, &mut tab).unwrap();
    assert!(tab.is_lost[0]);
}

#[test]
fn repeat_executes_body_n_times() {
    let (results, _) = run("REPEAT 2 { X 0 }\nM 0", 1);
    assert_eq!(results, vec![Some(false)]);
    let (results, _) = run("REPEAT 3 { X 0 }\nM 0", 1);
    assert_eq!(results, vec![Some(true)]);
}

#[test]
fn comments_and_annotations_are_no_ops() {
    let (results, _) = run(
        "# c\nQUBIT_COORDS(0,0) 0\nX 0\nTICK\nM 0",
        1,
    );
    assert_eq!(results, vec![Some(true)]);
}

#[test]
fn measurement_buffer_is_pre_sized() {
    let prog = parse("X 0\nM 0 1 2 3 4").unwrap();
    let tprog = normalize::to_tableau(&prog).unwrap();
    assert_eq!(tprog.expected_measurement_count, 5);
    let mut tab: Tab = GeneralizedTableau::new(5, 1e-10);
    let results = execute(&tprog, &mut tab).unwrap();
    assert_eq!(results.len(), 5);
}
