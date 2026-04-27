use ppvm_stim::{run_file, run_string, Error, ParseError, NormalizeError};
use ppvm_runtime::config::indexmap::ByteFxHashF64;
use ppvm_tableau::prelude::*;

type Tab = GeneralizedTableau<ByteFxHashF64<1>, usize>;

#[test]
fn run_string_executes_one_shot() {
    let mut tab: Tab = GeneralizedTableau::new(2, 1e-10);
    let results = run_string("X 0\nM 0 1", &mut tab).unwrap();
    assert_eq!(results, vec![Some(true), Some(false)]);
}

#[test]
fn run_file_round_trips_with_run_string() {
    let circuit = "X 0\nH 1\nCX 1 2\nM 0 1 2";
    let path = std::env::temp_dir().join("ppvm_stim_test.stim");
    std::fs::write(&path, circuit).unwrap();

    let mut a: GeneralizedTableau<ByteFxHashF64<1>, usize> = GeneralizedTableau::new_with_seed(3, 1e-10, 7);
    let mut b: GeneralizedTableau<ByteFxHashF64<1>, usize> = GeneralizedTableau::new_with_seed(3, 1e-10, 7);
    let r_str = run_string(circuit, &mut a).unwrap();
    let r_file = run_file(&path, &mut b).unwrap();
    assert_eq!(r_str, r_file);
}

#[test]
fn run_string_propagates_parse_error() {
    let mut tab: Tab = GeneralizedTableau::new(1, 1e-10);
    let err = run_string("FROBNICATE 0", &mut tab).unwrap_err();
    assert!(matches!(err, Error::Parse(ParseError::UnknownInstruction { .. })));
}

#[test]
fn run_string_propagates_normalize_error() {
    let mut tab: Tab = GeneralizedTableau::new(2, 1e-10);
    let err = run_string("SWAP 0 1", &mut tab).unwrap_err();
    assert!(matches!(err, Error::Normalize(NormalizeError::Unsupported { .. })));
}

#[test]
fn run_file_missing_file_yields_io_panic_or_error() {
    // Phase 1 surfaces IO failures as a panic from std::fs::read_to_string,
    // matching today's behavior. The Python tests assert `BaseException`,
    // which catches that. We mirror that here with std::panic::catch_unwind.
    let mut tab: Tab = GeneralizedTableau::new(1, 1e-10);
    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        run_file(std::path::Path::new("/nonexistent/x.stim"), &mut tab)
    }));
    assert!(result.is_err() || result.unwrap().is_err());
}
