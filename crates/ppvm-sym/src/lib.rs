//! Symbolic, parametric Pauli propagation.
//!
//! `ppvm-sym` provides a [`Term`] type that represents a polynomial in
//! sines and cosines of symbolic parameters. It is used by ppvm to
//! propagate Pauli operators through *parametric* circuits (e.g.,
//! variational ansätze) without committing to specific angle values
//! until the very end.
//!
//! # Quick example
//!
//! Build the symbolic expression `sin(x0) * cos(x1)`, then evaluate it
//! at a concrete `(x0, x1)`:
//!
//! ```
//! use ppvm_sym::Term;
//!
//! let expr = Term::var(0).sin() * Term::var(1).cos();
//!
//! let v = expr.eval(&[0.5, 1.0]).unwrap();
//! let expected = 0.5_f64.sin() * 1.0_f64.cos();
//! assert!((v - expected).abs() < 1e-12);
//! ```

mod add;
mod coeff;
mod display;
mod eval;
mod mul;
mod term;

pub use term::{Inner, Prod, Sum, Term};
