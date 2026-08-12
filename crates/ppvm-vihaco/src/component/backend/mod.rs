// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

#[cfg(feature = "legacy")]
mod legacy;
#[cfg(feature = "traits-2")]
mod traits_2;

#[cfg(feature = "legacy")]
pub use legacy::*;
#[cfg(feature = "traits-2")]
pub use traits_2::*;

// ── Where the randomness lives ──────────────────────────────────────────────
//
// The `-2` crates inject randomness: every stochastic method takes
// `&mut impl Rng`, so the caller owns the stream. Vihaco is a frontend, not a
// caller with an opinion — a device runs one trajectory and wants one seed at
// construction — so the *executor* is the owner. Each executor carries a
// `SmallRng` (seeded from `new_with_seed`, or OS entropy) and [`draw!`] threads
// it into every op that needs it.
//
// The legacy backends keep their generator inside the simulator, so there the
// executor holds no RNG and `draw!` is a plain forward. That is what lets the
// instruction-dispatch tables stay backend-agnostic.

/// Apply a stochastic op to `$ex.$field`.
///
/// Takes the whole executor rather than the state field so that the RNG the
/// `-2` arm reaches for never has to be named at the call site.
#[cfg(feature = "legacy")]
macro_rules! draw {
    ($ex:expr, $field:ident, $method:ident($($arg:expr),* $(,)?)) => {
        $ex.$field.$method($($arg),*)
    };
}

/// Apply a stochastic op to `$ex.$field`, passing the executor's RNG.
#[cfg(feature = "traits-2")]
macro_rules! draw {
    ($ex:expr, $field:ident, $method:ident($($arg:expr),* $(,)?)) => {
        $ex.$field.$method($($arg),* , &mut $ex.rng)
    };
}

/// Build an executor, giving it an RNG seeded from `$seed` when the backend
/// needs one. `None` seeds from OS entropy.
#[cfg(feature = "legacy")]
macro_rules! executor {
    ($wrapper:ident { $($field:ident: $value:expr),* $(,)? }, $seed:expr) => {{
        let _: Option<u64> = $seed;
        $wrapper { $($field: $value),* }
    }};
}

/// Build an executor together with the RNG it owns.
#[cfg(feature = "traits-2")]
macro_rules! executor {
    ($wrapper:ident { $($field:ident: $value:expr),* $(,)? }, $seed:expr) => {
        $wrapper {
            $($field: $value,)*
            rng: $crate::component::backend::make_rng($seed),
        }
    };
}

/// Construct the tableau an executor wraps; legacy seeds it directly.
#[cfg(feature = "legacy")]
macro_rules! new_tableau {
    ($new:ident, $seeded:ident, $n:expr, $threshold:expr, $seed:expr) => {
        match $seed {
            Some(s) => Backend::$seeded($n, $threshold, s),
            None => Backend::$new($n, $threshold),
        }
    };
}

/// Construct the tableau an executor wraps; under `-2` the seed belongs to the
/// executor's RNG instead, so it is not consulted here.
#[cfg(feature = "traits-2")]
macro_rules! new_tableau {
    ($new:ident, $seeded:ident, $n:expr, $threshold:expr, $seed:expr) => {
        Backend::$new($n, $threshold)
    };
}

pub(crate) use {draw, executor, new_tableau};
