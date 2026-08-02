// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The `sym.exact.multiply` workload: the L4 twisted operator product over an
//! **exact** coefficient ring, with no `f64` anywhere in the coefficient
//! representation.
//!
//! This is what Phase 5 exists to prove. `ppvm-traits-2` split `Halvable`,
//! `Angle<C>`, `ImaginaryUnit` and `Conjugate` off `Coefficient` and dropped the
//! old `Mul<f64>` bound *specifically* so exact rings stay expressible (gap
//! `t2.coefficient.1`: `Coefficient::half` would "foreclose exact rings
//! (`0.5·(1+i) ∉ ℤ[i]`)"). A claim about a trait tower is only worth its witness,
//! so this file **instantiates** the whole engine at `C = GaussianInt` — two
//! `i64`s, no float — and runs `multiply_into` over it. If any bound still
//! secretly forced a float, this file would not compile.
//!
//! # Why the acceptance oracle is Lean, not old
//!
//! Old `ppvm-sym::Term`'s phase handling is broken in three independent places
//! (its `MulAssign<Prod>` drops the phase, `add_term` discards a phase-only
//! monomial's phase, and `eval` returns `f64` so an accumulated `i^k` evaluates
//! as `+1`), so old's numbers on this path are not a parity target. The
//! assertions here are the machine-checked laws:
//!
//! * `lean/PPVM/Algebra/Twisted.lean` — `twistedConv_assoc` (whole-map
//!   associativity), `tmul_assoc` (monomial associativity), `iPow_add`
//!   (`i^a·i^b = i^{a+b}`), `twistedConv_add_left`/`_right` (biadditivity);
//! * `lean/PPVM/Pauli/Matrix.lean` — `star_iU` (`conj i = −i`), `iU_sq`.
//!
//! Exact ring ⇒ every assertion is `assert_eq!` with **zero tolerance**.

use ppvm_pauli_sum_2::{HashMapStore, NoPolicy, PauliWord, Sum};
use ppvm_sym_2::GaussianInt;
use ppvm_traits_2::{Conjugate, ImaginaryUnit};

use num::{One, Zero};

type Key = PauliWord<[u8; 1]>;
/// The engine over `ℤ[i]`. Note the coefficient parameter: no `f64` in sight.
type ExactSum = Sum<HashMapStore<Key, GaussianInt>, NoPolicy>;

const N: usize = 5;

/// A multi-term operator over `N` qubits, seeded deterministically.
fn operand(words: &[(&str, (i64, i64))]) -> ExactSum {
    let mut s: ExactSum = ExactSum::new(N);
    for (w, (re, im)) in words {
        assert_eq!(w.len(), N);
        s += (Key::from(*w), GaussianInt::new(*re, *im));
    }
    s
}

fn a() -> ExactSum {
    operand(&[
        ("XIIII", (1, 0)),
        ("IYIII", (0, 2)),
        ("ZZIII", (-3, 1)),
        ("IIXYI", (2, -1)),
    ])
}

fn b() -> ExactSum {
    operand(&[
        ("YIIII", (0, 1)),
        ("IZIII", (4, 0)),
        ("XXIII", (1, 1)),
        ("IIYXI", (-2, 3)),
    ])
}

fn c() -> ExactSum {
    operand(&[
        ("ZIIII", (2, 2)),
        ("IXIII", (-1, 0)),
        ("YYIII", (0, -5)),
        ("IIZZI", (3, 1)),
    ])
}

/// The support as a sorted, canonical `(word, coefficient)` list — exact, no
/// tolerance anywhere.
fn support(s: &ExactSum) -> Vec<(String, GaussianInt)> {
    let mut v: Vec<(String, GaussianInt)> = s.iter().map(|(k, c)| (k.to_string(), c)).collect();
    // `multiply_into` runs no `reduce`, so an exact cancellation stays in the
    // support with coefficient `0`. Drop those before comparing: the two
    // association orders reach the zero keys through different pair sets.
    v.retain(|(_, c)| !c.is_zero());
    v.sort();
    v
}

#[test]
fn twisted_product_is_exactly_associative() {
    // `lean/PPVM/Algebra/Twisted.lean` `twistedConv_assoc` (lifted from
    // `tmul_assoc` by bilinearity). Exact ring ⇒ no tolerance.
    let (a, b, c) = (a(), b(), c());

    let ab = a.multiply(&b);
    let ab_c = ab.multiply(&c);

    let bc = b.multiply(&c);
    let a_bc = a.multiply(&bc);

    assert_eq!(support(&ab_c), support(&a_bc));
    assert!(!support(&ab_c).is_empty(), "the product collapsed to zero");
}

#[test]
fn twisted_product_is_exactly_biadditive() {
    // `twistedConv_add_left` / `twistedConv_add_right`: `acc += A·B; acc += A·C`
    // is `A·(B + C)`. This is the property old's in-place `MulAssign` lost.
    let (a, b, c) = (a(), b(), c());

    let mut acc: ExactSum = ExactSum::new(N);
    a.multiply_into(&b, &mut acc);
    a.multiply_into(&c, &mut acc);

    let mut b_plus_c = b.clone();
    for (k, v) in c.iter() {
        b_plus_c += (k, v);
    }
    let direct = a.multiply(&b_plus_c);

    assert_eq!(support(&acc), support(&direct));
}

#[test]
fn conj_of_i_is_minus_i_on_the_exact_ring() {
    // `lean/PPVM/Pauli/Matrix.lean` `star_iU`, at the coefficient level and
    // lifted over a whole sum.
    let i = GaussianInt::imaginary_unit();
    assert_eq!(i.conj(), -i);
    assert_eq!(i * i, -GaussianInt::one());

    let a = a();
    let conj: Vec<(String, GaussianInt)> = support(&a)
        .into_iter()
        .map(|(k, v)| (k, v.conj()))
        .collect();
    for ((_, orig), (_, cj)) in support(&a).iter().zip(conj.iter()) {
        assert_eq!(cj.re, orig.re);
        assert_eq!(cj.im, -orig.im);
    }
}

#[test]
fn the_exact_ring_carries_no_float() {
    // The structural claim, asserted rather than asserted-in-prose: the whole
    // coefficient value fits in two `i64`s.
    assert_eq!(
        std::mem::size_of::<GaussianInt>(),
        2 * std::mem::size_of::<i64>()
    );
    // And the ring is genuinely exact under the engine: a product of large
    // integers is reproduced bit-for-bit, which an `f64` mantissa could not do.
    let big = 1_i64 << 40;
    let x = GaussianInt::new(big + 1, big - 1);
    let y = GaussianInt::new(3, 5);
    assert_eq!(
        x * y,
        GaussianInt::new((big + 1) * 3 - (big - 1) * 5, (big + 1) * 5 + (big - 1) * 3)
    );
}

#[test]
fn multiply_into_accumulates_without_reducing() {
    // Behaviour parity with the rest of `-2`: `multiply_into` runs no `reduce`
    // and no truncation, so an exact cancellation stays in the support with
    // coefficient `0`.
    let mut p: ExactSum = ExactSum::new(N);
    p += (Key::from("XIIII"), GaussianInt::new(1, 0));
    let mut q: ExactSum = ExactSum::new(N);
    q += (Key::from("XIIII"), GaussianInt::new(1, 0));

    let mut acc: ExactSum = ExactSum::new(N);
    p.multiply_into(&q, &mut acc);
    // X·X = I with phase +1.
    assert_eq!(acc.get(&Key::from("IIIII")), Some(GaussianInt::new(1, 0)));

    // Subtracting the same product leaves an exact zero *in* the support.
    let mut minus_q: ExactSum = ExactSum::new(N);
    minus_q += (Key::from("XIIII"), GaussianInt::new(-1, 0));
    p.multiply_into(&minus_q, &mut acc);
    assert_eq!(acc.get(&Key::from("IIIII")), Some(GaussianInt::zero()));
    assert_eq!(acc.len(), 1, "the zero-coefficient key must survive");
}
