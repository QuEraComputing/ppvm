// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

#[cfg(all(feature = "legacy", feature = "traits-2"))]
compile_error!("features `legacy` and `traits-2` are mutually exclusive");
#[cfg(not(any(feature = "legacy", feature = "traits-2")))]
compile_error!("enable exactly one Python backend: `legacy` or `traits-2`");

#[cfg(feature = "legacy")]
mod legacy;
#[cfg(feature = "traits-2")]
mod traits_2;

#[cfg(feature = "legacy")]
pub use legacy::*;
#[cfg(feature = "traits-2")]
pub use traits_2::*;

#[cfg(feature = "legacy")]
pub const NAME: &str = "legacy";
#[cfg(feature = "traits-2")]
pub const NAME: &str = "traits-2";

// ── Where the randomness lives ──────────────────────────────────────────────
//
// The `-2` crates inject randomness: every stochastic method takes
// `&mut impl Rng`, so the *caller* owns the stream. That is the right Rust API —
// seeding, forking and replaying are all the caller's to control — but it is the
// wrong Python API, where the high-level user should never have to thread a
// generator through `tab.measure(0)`.
//
// So the `#[pyclass]` wrapper is the owner: it holds one `SmallRng` per
// simulator object, seeded from the constructor's `seed=` argument (or OS
// entropy), and `draw!` threads it into every call that needs it. The legacy
// backend keeps its generator inside the tableau, so there the wrapper carries
// no RNG and `draw!` is a plain forward — which is what lets the `interface*`
// macros stay backend-agnostic.

/// Call a stochastic method on the wrapped simulator.
///
/// Written to read exactly like the direct call it replaces:
/// `draw!(self.inner.measure(addr0))`.
#[cfg(feature = "legacy")]
macro_rules! draw {
    ($this:ident . inner . $method:ident ( $($arg:expr),* $(,)? )) => {
        $this.inner.$method($($arg),*)
    };
}

/// Call a stochastic method on the wrapped simulator, passing the wrapper's RNG.
#[cfg(feature = "traits-2")]
macro_rules! draw {
    ($this:ident . inner . $method:ident ( $($arg:expr),* $(,)? )) => {
        $this.inner.$method($($arg),* , &mut $this.rng)
    };
}

/// Build a `#[pyclass]` wrapper around `$inner`, giving it an RNG seeded from
/// `$seed` (`None` = OS entropy) when the backend needs one.
#[cfg(feature = "legacy")]
macro_rules! wrap {
    ($inner:expr, $seed:expr) => {{
        let _: Option<u64> = $seed;
        Self { inner: $inner }
    }};
}

/// Build a `#[pyclass]` wrapper around `$inner` and the RNG it owns.
#[cfg(feature = "traits-2")]
macro_rules! wrap {
    ($inner:expr, $seed:expr) => {
        Self {
            inner: $inner,
            rng: $crate::backend::make_rng($seed),
        }
    };
}

/// Clone a `#[pyclass]` wrapper, carrying the RNG state across unchanged.
///
/// This is `__copy__` / `__deepcopy__` semantics: the copy replays the same
/// stream as the original from this point on. `fork` builds with [`wrap!`]
/// instead, which reseeds.
#[cfg(feature = "legacy")]
macro_rules! wrap_cloned {
    ($this:ident) => {
        Self {
            inner: $this.inner.clone(),
        }
    };
}

/// Clone a `#[pyclass]` wrapper, carrying the RNG state across unchanged.
#[cfg(feature = "traits-2")]
macro_rules! wrap_cloned {
    ($this:ident) => {
        Self {
            inner: $this.inner.clone(),
            rng: $this.rng.clone(),
        }
    };
}
