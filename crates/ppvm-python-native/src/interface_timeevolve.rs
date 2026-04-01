use ppvm_runtime::prelude::{Config, PauliSum};
use ppvm_timeevolve::lindblad::{JumpOp, LadderDirection, LadderOp, LindbladOp, RateMatrix};
use ppvm_timeevolve::product_state::ProductState;
use ppvm_timeevolve::solve::{SolverCache, SolverConfig, solve_cached};
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyAnyMethods;

// ──────────────────────────────────────────────────────────────────
// Helpers
// ──────────────────────────────────────────────────────────────────

fn parse_rate_matrix(rates: &Bound<PyAny>) -> PyResult<RateMatrix> {
    // Try diagonal (list[float]) first — more specific, checked first.
    if let Ok(v) = rates.extract::<Vec<f64>>() {
        return Ok(RateMatrix::Vector(v));
    }
    if let Ok(m) = rates.extract::<Vec<Vec<f64>>>() {
        return Ok(RateMatrix::Dense(m));
    }
    Err(PyTypeError::new_err(
        "rates must be list[float] (diagonal) or list[list[float]] (dense)",
    ))
}

fn build_jump_ops<T: Config>(ops: &[(usize, String)]) -> PyResult<Vec<JumpOp<T>>> {
    ops.iter()
        .map(|(qubit, dir)| {
            let direction = match dir.as_str() {
                "raise" => LadderDirection::Raise,
                "lower" => LadderDirection::Lower,
                other => {
                    return Err(PyValueError::new_err(format!(
                        "invalid direction '{other}': expected 'raise' or 'lower'"
                    )))
                }
            };
            Ok(JumpOp::Ladder(LadderOp { qubit: *qubit, direction }))
        })
        .collect()
}

fn build_config(
    rtol: f64,
    atol: f64,
    h0: Option<f64>,
    hmin: f64,
    hmax: f64,
) -> SolverConfig {
    SolverConfig { rtol, atol, h0, hmin, hmax }
}

// ──────────────────────────────────────────────────────────────────
// solve_timeevolve_states
// Returns (times, list-of-native-PauliSum-objects) — one clone per save point.
// ──────────────────────────────────────────────────────────────────

#[pyfunction]
#[allow(clippy::too_many_arguments)]
pub fn solve_timeevolve_states(
    py: Python<'_>,
    observable: &Bound<PyAny>,
    lindblad_ops: Vec<(usize, String)>,
    rates: &Bound<PyAny>,
    t_span_start: f64,
    t_span_end: f64,
    save_at: Vec<f64>,
    hamiltonian: Option<&Bound<PyAny>>,
    rtol: f64,
    atol: f64,
    h0: Option<f64>,
    hmin: f64,
    hmax: f64,
) -> PyResult<(Vec<f64>, Vec<Py<PyAny>>)> {
    let rate_matrix = parse_rate_matrix(rates)?;
    let config = build_config(rtol, atol, h0, hmin, hmax);

    // Generate one if-let arm per concrete PauliSum type (N = 0..15).
    macro_rules! try_arm {
        ($N:literal) => {
            paste::paste! {
                if let Ok(s) = observable.cast::<crate::interface::[<PauliSumIndexMapFxHash $N>]>() {
                    let s_ref = s.borrow();
                    let ham_opt: Option<PyRef<crate::interface::[<PauliSumIndexMapFxHash $N>]>> =
                        match hamiltonian {
                            Some(h) => Some(
                                h.cast::<crate::interface::[<PauliSumIndexMapFxHash $N>]>()
                                    .map_err(|_| PyTypeError::new_err(
                                        "hamiltonian and observable must have the same native type"
                                    ))?
                                    .borrow(),
                            ),
                            None => None,
                        };
                    let ops = build_jump_ops::<crate::interface::[<IndexMapFxHash $N>]>(
                        &lindblad_ops,
                    )?;
                    let lindblad = LindbladOp::new(ops, rate_matrix);
                    let mut cache = SolverCache::new(&s_ref.inner);
                    let (times, raw) = solve_cached(
                        ham_opt.as_ref().map(|h| &h.inner),
                        &lindblad,
                        &s_ref.inner,
                        (t_span_start, t_span_end),
                        &save_at,
                        |_, p: &PauliSum<crate::interface::[<IndexMapFxHash $N>]>| p.clone(),
                        config,
                        &mut cache,
                    );
                    let results = raw
                        .into_iter()
                        .map(|inner| -> PyResult<Py<PyAny>> {
                            Ok(Py::new(
                                py,
                                crate::interface::[<PauliSumIndexMapFxHash $N>] { inner },
                            )?
                            .into_any())
                        })
                        .collect::<PyResult<Vec<_>>>()?;
                    return Ok((times, results));
                }
            }
        };
    }

    try_arm!(0);
    try_arm!(1);
    try_arm!(2);
    try_arm!(3);
    try_arm!(4);
    try_arm!(5);
    try_arm!(6);
    try_arm!(7);
    try_arm!(8);
    try_arm!(9);
    try_arm!(10);
    try_arm!(11);
    try_arm!(12);
    try_arm!(13);
    try_arm!(14);
    try_arm!(15);

    Err(PyTypeError::new_err(
        "unsupported observable type: expected PauliSumIndexMapFxHash0 through PauliSumIndexMapFxHash15",
    ))
}

// ──────────────────────────────────────────────────────────────────
// solve_timeevolve_expectation
// Returns (times, list[float]) — ⟨O(t)⟩ = Tr(ρ₀ O(t)) at each save point.
// ProductState is constructed once; no PauliSum cloning per save point.
// ──────────────────────────────────────────────────────────────────

#[pyfunction]
#[allow(clippy::too_many_arguments)]
pub fn solve_timeevolve_expectation(
    observable: &Bound<PyAny>,
    bloch_vectors: Vec<f64>,
    lindblad_ops: Vec<(usize, String)>,
    rates: &Bound<PyAny>,
    t_span_start: f64,
    t_span_end: f64,
    save_at: Vec<f64>,
    hamiltonian: Option<&Bound<PyAny>>,
    rtol: f64,
    atol: f64,
    h0: Option<f64>,
    hmin: f64,
    hmax: f64,
) -> PyResult<(Vec<f64>, Vec<f64>)> {
    let rate_matrix = parse_rate_matrix(rates)?;
    let config = build_config(rtol, atol, h0, hmin, hmax);
    // Construct ProductState once before dispatch; captured by reference in callback.
    let ps = ProductState::from_flat(&bloch_vectors);

    macro_rules! try_arm {
        ($N:literal) => {
            paste::paste! {
                if let Ok(s) = observable.cast::<crate::interface::[<PauliSumIndexMapFxHash $N>]>() {
                    let s_ref = s.borrow();
                    let ham_opt: Option<PyRef<crate::interface::[<PauliSumIndexMapFxHash $N>]>> =
                        match hamiltonian {
                            Some(h) => Some(
                                h.cast::<crate::interface::[<PauliSumIndexMapFxHash $N>]>()
                                    .map_err(|_| PyTypeError::new_err(
                                        "hamiltonian and observable must have the same native type"
                                    ))?
                                    .borrow(),
                            ),
                            None => None,
                        };
                    let ops = build_jump_ops::<crate::interface::[<IndexMapFxHash $N>]>(
                        &lindblad_ops,
                    )?;
                    let lindblad = LindbladOp::new(ops, rate_matrix);
                    let mut cache = SolverCache::new(&s_ref.inner);
                    let (times, results) = solve_cached(
                        ham_opt.as_ref().map(|h| &h.inner),
                        &lindblad,
                        &s_ref.inner,
                        (t_span_start, t_span_end),
                        &save_at,
                        |_, p: &PauliSum<crate::interface::[<IndexMapFxHash $N>]>| ps.expectation(p),
                        config,
                        &mut cache,
                    );
                    return Ok((times, results));
                }
            }
        };
    }

    try_arm!(0);
    try_arm!(1);
    try_arm!(2);
    try_arm!(3);
    try_arm!(4);
    try_arm!(5);
    try_arm!(6);
    try_arm!(7);
    try_arm!(8);
    try_arm!(9);
    try_arm!(10);
    try_arm!(11);
    try_arm!(12);
    try_arm!(13);
    try_arm!(14);
    try_arm!(15);

    Err(PyTypeError::new_err(
        "unsupported observable type: expected PauliSumIndexMapFxHash0 through PauliSumIndexMapFxHash15",
    ))
}

// ──────────────────────────────────────────────────────────────────
// product_state_expectation
// Standalone Tr(ρ₀ O) on a raw PauliSum — no ODE solve.
// ──────────────────────────────────────────────────────────────────

#[pyfunction]
pub fn product_state_expectation(
    observable: &Bound<PyAny>,
    bloch_vectors: Vec<f64>,
) -> PyResult<f64> {
    let ps = ProductState::from_flat(&bloch_vectors);

    macro_rules! try_arm {
        ($N:literal) => {
            paste::paste! {
                if let Ok(s) = observable.cast::<crate::interface::[<PauliSumIndexMapFxHash $N>]>() {
                    let s_ref = s.borrow();
                    return Ok(ps.expectation(&s_ref.inner));
                }
            }
        };
    }

    try_arm!(0);
    try_arm!(1);
    try_arm!(2);
    try_arm!(3);
    try_arm!(4);
    try_arm!(5);
    try_arm!(6);
    try_arm!(7);
    try_arm!(8);
    try_arm!(9);
    try_arm!(10);
    try_arm!(11);
    try_arm!(12);
    try_arm!(13);
    try_arm!(14);
    try_arm!(15);

    Err(PyTypeError::new_err(
        "unsupported observable type: expected PauliSumIndexMapFxHash0 through PauliSumIndexMapFxHash15",
    ))
}
