use ppvm_stim::{ExecError, parse_extended, prepare};

fn err_from_src(src: &str) -> ExecError {
    let prog = parse_extended(src).expect("parse_extended");
    prepare(&prog).expect_err("must reject")
}

#[test]
fn unsupported_swap_rejected() {
    let ExecError::Unsupported { name, line } = err_from_src("SWAP 0 1") else {
        panic!("expected ExecError::Unsupported");
    };
    assert_eq!(name, "SWAP");
    assert_eq!(line, 1);
}

#[test]
fn unsupported_mx_rejected() {
    let e = err_from_src("MX 0");
    assert!(matches!(e, ExecError::Unsupported { .. }));
}

#[test]
fn unsupported_heralded_erase_rejected() {
    let e = err_from_src("HERALDED_ERASE(0.1) 0");
    assert!(matches!(e, ExecError::Unsupported { .. }));
}

#[test]
fn unsupported_swap_inside_repeat_rejected() {
    let ExecError::Unsupported { name, line } = err_from_src("REPEAT 3 {\n    SWAP 0 1\n}\n")
    else {
        panic!("expected ExecError::Unsupported");
    };
    assert_eq!(name, "SWAP");
    assert_eq!(line, 2);
}

#[test]
fn supported_structural_instructions_are_not_rejected_by_prepare() {
    let prog =
        parse_extended("MPAD 0 1\nI_ERROR[correlated_loss](0.5) 0 1\n").expect("parse_extended");
    assert_eq!(prepare(&prog), Ok(()));
}
