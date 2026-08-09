// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Behaviour-parity tests for the gate/noise surface ported from `ppvm-traits`:
//! [`CliffordExtensions`] (and its blanket impl), the batched gate traits,
//! [`Reset`], the stim [`PauliError`] aliases, and the channel family.
//!
//! Two distinct obligations are pinned here:
//!
//! 1. **The blanket `CliffordExtensions`** (`pauli.rs`) expresses each extension
//!    gate as a product of the audited [`Clifford`] generators. That is only
//!    behaviour-preserving if the products reproduce the old crate's conjugation
//!    tables *including signs*, so the sweep below replays them on a ℤ₄-phased
//!    stub whose per-generator phase rules are byte-for-byte the ones
//!    `ppvm-phased-pauli-word-2` ports from the old fused kernel (and which
//!    `ppvm-conformance-2` diffs against the old crate). The expected values are
//!    the tables written in `ppvm-traits/src/traits/clifford.rs`:
//!    the six single-qubit rows and the 16-entry `CY` table.
//! 2. **Every default body** (`zcy`, the `*_many` loops, `reset_x` = `reset` then
//!    `h`, `x_error(p)` = `pauli_error([p,0,0])`, …) fires the same underlying
//!    required method, the same number of times, in the same order — the old
//!    crate's observable contract.

use ppvm_traits_2::{
    AmplitudeDamping, AsymmetricLossChannel, BlanketClifford, Clifford, CliffordBatch,
    CliffordExtensions, CliffordExtensionsBatch, CorrelatedLossChannel, Depolarizing,
    Depolarizing2, LossChannel, PauliError, PauliErrorAll, PhaseTrack, Reset, ResetLossChannel,
    SymplecticColumns, TwoQubitPauliError,
};
use rand::SeedableRng;
use rand::rngs::SmallRng;

// ---------------------------------------------------------------------------
// A two-qubit ℤ₄-phased stub: the smallest honest `SymplecticColumns` +
// `PhaseTrack` opt-in, so the blanket `Clifford`/`CliffordExtensions` impls run
// with real bits *and* real signs.
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
struct Z4Word {
    x: [bool; 2],
    z: [bool; 2],
    /// ℤ₄ phase exponent (`+1 → 0`, `i → 1`, `−1 → 2`, `−i → 3`).
    phase: u8,
}

impl Z4Word {
    /// Parse `"XY"`-style two-site strings with an implicit `+` phase.
    fn new(sites: &str) -> Self {
        let mut w = Z4Word {
            x: [false; 2],
            z: [false; 2],
            phase: 0,
        };
        for (q, c) in sites.chars().enumerate() {
            let (x, z) = match c {
                'I' => (false, false),
                'X' => (true, false),
                'Y' => (true, true),
                'Z' => (false, true),
                other => panic!("bad Pauli {other}"),
            };
            w.x[q] = x;
            w.z[q] = z;
        }
        w
    }

    /// Render as `"+XY"` / `"-ZI"`; the extension gates are real ±1 conjugations,
    /// so a `±i` phase would itself be a failure.
    fn render(&self) -> String {
        let sign = match self.phase {
            0 => "+",
            2 => "-",
            other => panic!("Clifford conjugation produced a non-real phase {other}"),
        };
        let sites: String = (0..2)
            .map(|q| match (self.x[q], self.z[q]) {
                (false, false) => 'I',
                (true, false) => 'X',
                (true, true) => 'Y',
                (false, true) => 'Z',
            })
            .collect();
        format!("{sign}{sites}")
    }

    #[inline]
    fn flip_sign_if(&mut self, cond: bool) {
        if cond {
            self.phase = (self.phase + 2) % 4;
        }
    }
}

impl SymplecticColumns for Z4Word {
    fn n_qubits(&self) -> usize {
        2
    }
    fn swap_xz(&mut self, q: usize) {
        let (x, z) = (self.x[q], self.z[q]);
        self.x[q] = z;
        self.z[q] = x;
    }
    fn xor_z_from_x(&mut self, q: usize) {
        self.z[q] ^= self.x[q];
    }
    fn xor_x_col(&mut self, ctrl: usize, tgt: usize) {
        self.x[tgt] ^= self.x[ctrl];
    }
    fn xor_z_col(&mut self, tgt: usize, ctrl: usize) {
        self.z[ctrl] ^= self.z[tgt];
    }
    fn cz_bits(&mut self, a: usize, b: usize) {
        let (xa, xb) = (self.x[a], self.x[b]);
        self.z[a] ^= xb;
        self.z[b] ^= xa;
    }
}

/// The per-generator signs of `ppvm-phased-pauli-word-2` (ported from the old
/// fused kernel, machine-checked in `lean/PPVM/Pauli/Conjugation.lean`): every
/// delta is computed from the *pre-mutation* bits, which is exactly how the
/// blanket sequences them (phase primitive first, column primitive second).
impl PhaseTrack for Z4Word {
    fn flip_phase_where_xz(&mut self, q: usize) {
        self.flip_sign_if(self.x[q] && self.z[q]);
    }
    fn s_phase(&mut self, q: usize) {
        self.flip_sign_if(self.x[q] && !self.z[q]);
    }
    fn cnot_phase(&mut self, ctrl: usize, tgt: usize) {
        self.flip_sign_if(self.x[ctrl] && self.z[tgt] && (self.x[tgt] == self.z[ctrl]));
    }
    fn cz_phase(&mut self, a: usize, b: usize) {
        self.flip_sign_if(self.x[a] && self.x[b] && (self.z[a] ^ self.z[b]));
    }
    fn x_phase(&mut self, q: usize) {
        self.flip_sign_if(self.z[q]);
    }
    fn y_phase(&mut self, q: usize) {
        self.flip_sign_if(self.x[q] ^ self.z[q]);
    }
    fn z_phase(&mut self, q: usize) {
        self.flip_sign_if(self.x[q]);
    }
}

impl BlanketClifford for Z4Word {}

fn conj(input: &str, gate: impl Fn(&mut Z4Word)) -> String {
    let mut w = Z4Word::new(input);
    gate(&mut w);
    w.render()
}

// ---------------------------------------------------------------------------
// 1. The blanket CliffordExtensions reproduces the old conjugation tables.
// ---------------------------------------------------------------------------

/// The single-qubit table of `ppvm-traits/src/traits/clifford.rs`:
///
/// | Gate | X | Y | Z |
/// |:---:|:---:|:---:|:---:|
/// | `s` | `-Y` | `X` | `Z` |
/// | `s_dag` | `Y` | `-X` | `Z` |
/// | `sqrt_x` | `X` | `-Z` | `Y` |
/// | `sqrt_x_dag` | `X` | `Z` | `-Y` |
/// | `sqrt_y` | `Z` | `Y` | `-X` |
/// | `sqrt_y_dag` | `-Z` | `Y` | `X` |
#[test]
fn blanket_clifford_extensions_match_old_conjugation_table() {
    // `S` itself (the audited generator) fixes the convention the rest inherit.
    for (input, want) in [("II", "+II"), ("XI", "-YI"), ("YI", "+XI"), ("ZI", "+ZI")] {
        assert_eq!(conj(input, |w| w.s(0)), want, "s {input}");
    }
    for (input, want) in [("II", "+II"), ("XI", "+YI"), ("YI", "-XI"), ("ZI", "+ZI")] {
        assert_eq!(conj(input, |w| w.s_dag(0)), want, "s_dag {input}");
    }
    for (input, want) in [("II", "+II"), ("XI", "+XI"), ("YI", "-ZI"), ("ZI", "+YI")] {
        assert_eq!(conj(input, |w| w.sqrt_x(0)), want, "sqrt_x {input}");
    }
    for (input, want) in [("II", "+II"), ("XI", "+XI"), ("YI", "+ZI"), ("ZI", "-YI")] {
        assert_eq!(conj(input, |w| w.sqrt_x_dag(0)), want, "sqrt_x_dag {input}");
    }
    for (input, want) in [("II", "+II"), ("XI", "+ZI"), ("YI", "+YI"), ("ZI", "-XI")] {
        assert_eq!(conj(input, |w| w.sqrt_y(0)), want, "sqrt_y {input}");
    }
    for (input, want) in [("II", "+II"), ("XI", "-ZI"), ("YI", "+YI"), ("ZI", "+XI")] {
        assert_eq!(conj(input, |w| w.sqrt_y_dag(0)), want, "sqrt_y_dag {input}");
    }
}

/// Each extension gate is the inverse conjugation of its dagger, on every input —
/// the property the `S³ = S†` / `H·S·H` products must not silently break.
#[test]
fn blanket_extension_daggers_are_inverse_conjugations() {
    for a in ["I", "X", "Y", "Z"] {
        for b in ["I", "X", "Y", "Z"] {
            let input = format!("{a}{b}");
            let start = Z4Word::new(&input).render();
            assert_eq!(
                conj(&input, |w| {
                    w.s(0);
                    w.s_dag(0);
                }),
                start,
                "s∘s_dag {input}"
            );
            assert_eq!(
                conj(&input, |w| {
                    w.sqrt_x(0);
                    w.sqrt_x_dag(0);
                }),
                start,
                "sqrt_x∘sqrt_x_dag {input}"
            );
            assert_eq!(
                conj(&input, |w| {
                    w.sqrt_y(0);
                    w.sqrt_y_dag(0);
                }),
                start,
                "sqrt_y∘sqrt_y_dag {input}"
            );
            // √X applied twice is X-conjugation; likewise √Y and S² = Z.
            assert_eq!(
                conj(&input, |w| {
                    w.sqrt_x(0);
                    w.sqrt_x(0);
                }),
                conj(&input, |w| w.x(0)),
                "sqrt_x² = x {input}"
            );
            assert_eq!(
                conj(&input, |w| {
                    w.sqrt_y(0);
                    w.sqrt_y(0);
                }),
                conj(&input, |w| w.y(0)),
                "sqrt_y² = y {input}"
            );
            assert_eq!(
                conj(&input, |w| {
                    w.s(0);
                    w.s(0);
                }),
                conj(&input, |w| w.z(0)),
                "s² = z {input}"
            );
        }
    }
}

/// The 16-entry `CY` table of `ppvm-traits/src/traits/clifford.rs` (rows =
/// control, columns = target), signs included.
#[test]
fn blanket_cy_matches_old_two_qubit_table() {
    const TABLE: [(&str, &str); 16] = [
        ("II", "+II"),
        ("IX", "+ZX"),
        ("IY", "+IY"),
        ("IZ", "+ZZ"),
        ("XI", "+XY"),
        ("XX", "-YZ"),
        ("XY", "+XI"),
        ("XZ", "+YX"),
        ("YI", "+YY"),
        ("YX", "+XZ"),
        ("YY", "+YI"),
        ("YZ", "-XX"),
        ("ZI", "+ZI"),
        ("ZX", "+IX"),
        ("ZY", "+ZY"),
        ("ZZ", "+IZ"),
    ];
    for (input, want) in TABLE {
        assert_eq!(conj(input, |w| w.cy(0, 1)), want, "cy {input}");
        // The stim alias must be the same gate.
        assert_eq!(conj(input, |w| w.zcy(0, 1)), want, "zcy {input}");
    }
}

// ---------------------------------------------------------------------------
// 1b. The blanket `Clifford` itself (the generators the extensions are built
//     from) against the Lean conjugation oracle.
//
// The tables above pin the *derived* gates; the base generators the blanket
// emits — `h`, `s`, `x`, `y`, `z`, `cnot`, `cz` — are pinned here directly
// against `lean/PPVM/Pauli/Conjugation.lean`, in the backward Heisenberg
// convention this crate uses (the crate's `s` is Lean's `conjSdag`):
//
//   * `conjH_bits`/`conjH_sign`         — swap `(x,z)`, `−1` iff `x∧z`;
//   * `conjSdag_bits`/`conjSdag_sign`   — `z ⊕= x`, `−1` iff `x∧¬z`;
//   * `conjX`/`conjY`/`conjZ`           — word fixed, `−1` iff `z` / `x⊕z` / `x`;
//   * `conjCNOT_bits`/`conjCNOT_sign`   — `z_c ⊕= z_t`, `x_t ⊕= x_c`,
//                                         `−1` iff `x_c ∧ z_t ∧ (x_t = z_c)`;
//   * `conjCZ_bits`/`conjCZ_sign`       — `z_c ⊕= x_t`, `z_t ⊕= x_c`,
//                                         `−1` iff `x_c ∧ x_t ∧ (z_c ≠ z_t)`.
//
// Two distinct things are checked. The *tables* below are literal transcriptions
// of the Lean generator lemmas (`conjH_X`/`conjH_Y`/`conjH_Z`, `conjSdag_X`/`_Y`/
// `_Z`, `conjCNOT_Xc`/`_Xt`/`_Zc`/`_Zt`/`_YcYt`, `conjCZ_Xc`/`_Xt`/`_Zc`/`_Zt`/
// `_YcXt`), so they are an oracle independent of the stub's own formulas. The
// full sweep then checks the blanket's *composition* on every input: the
// blanket runs the phase primitive **before** the column primitive, and every
// delta above is a function of the pre-mutation bits, so reordering the two
// (or dropping one) changes the answer — `cz` in particular has no other
// coverage in this file, and `cnot` is otherwise reached only through `cy`.
// ---------------------------------------------------------------------------

/// Lean `conjH`: swap the bit columns, `−1` iff `x∧z`.
fn oracle_h(w: Z4Word, q: usize) -> Z4Word {
    let mut o = w;
    if w.x[q] && w.z[q] {
        o.phase = (o.phase + 2) % 4;
    }
    o.x[q] = w.z[q];
    o.z[q] = w.x[q];
    o
}

/// Lean `conjSdag` (the crate's backward `s`): `z ⊕= x`, `−1` iff `x∧¬z`.
fn oracle_s(w: Z4Word, q: usize) -> Z4Word {
    let mut o = w;
    if w.x[q] && !w.z[q] {
        o.phase = (o.phase + 2) % 4;
    }
    o.z[q] = w.z[q] ^ w.x[q];
    o
}

/// Lean `conjX`/`conjY`/`conjZ`: the word is fixed, only the sign moves.
fn oracle_pauli(w: Z4Word, flip: bool) -> Z4Word {
    let mut o = w;
    if flip {
        o.phase = (o.phase + 2) % 4;
    }
    o
}

/// Lean `conjCNOT`: `z_c ⊕= z_t`, `x_t ⊕= x_c`, `−1` iff `x_c ∧ z_t ∧ (x_t = z_c)`.
fn oracle_cnot(w: Z4Word, c: usize, t: usize) -> Z4Word {
    let mut o = w;
    if w.x[c] && w.z[t] && (w.x[t] == w.z[c]) {
        o.phase = (o.phase + 2) % 4;
    }
    o.z[c] = w.z[c] ^ w.z[t];
    o.x[t] = w.x[t] ^ w.x[c];
    o
}

/// Lean `conjCZ`: `z_a ⊕= x_b`, `z_b ⊕= x_a`, `−1` iff `x_a ∧ x_b ∧ (z_a ≠ z_b)`.
fn oracle_cz(w: Z4Word, a: usize, b: usize) -> Z4Word {
    let mut o = w;
    if w.x[a] && w.x[b] && (w.z[a] ^ w.z[b]) {
        o.phase = (o.phase + 2) % 4;
    }
    o.z[a] = w.z[a] ^ w.x[b];
    o.z[b] = w.z[b] ^ w.x[a];
    o
}

/// Every two-site word, at every ℤ₄ starting phase (64 inputs).
fn all_words() -> Vec<Z4Word> {
    let mut out = Vec::new();
    for a in ["I", "X", "Y", "Z"] {
        for b in ["I", "X", "Y", "Z"] {
            for phase in 0..4u8 {
                let mut w = Z4Word::new(&format!("{a}{b}"));
                w.phase = phase;
                out.push(w);
            }
        }
    }
    out
}

#[test]
fn blanket_clifford_generators_match_the_lean_generator_tables() {
    // `conjH_X`, `conjH_Z`, `conjH_Y` (+ the trivial `I`).
    for (input, want) in [("II", "+II"), ("XI", "+ZI"), ("ZI", "+XI"), ("YI", "-YI")] {
        assert_eq!(conj(input, |w| w.h(0)), want, "h {input}");
    }
    // `conjSdag_X`, `conjSdag_Y`, `conjSdag_Z` — the *backward* signs.
    for (input, want) in [("II", "+II"), ("XI", "-YI"), ("YI", "+XI"), ("ZI", "+ZI")] {
        assert_eq!(conj(input, |w| w.s(0)), want, "s {input}");
    }
    // `conjX` / `conjY` / `conjZ`: word fixed, sign iff the two anticommute.
    for (input, want) in [("II", "+II"), ("XI", "+XI"), ("YI", "-YI"), ("ZI", "-ZI")] {
        assert_eq!(conj(input, |w| w.x(0)), want, "x {input}");
    }
    for (input, want) in [("II", "+II"), ("XI", "-XI"), ("YI", "+YI"), ("ZI", "-ZI")] {
        assert_eq!(conj(input, |w| w.y(0)), want, "y {input}");
    }
    for (input, want) in [("II", "+II"), ("XI", "-XI"), ("YI", "-YI"), ("ZI", "+ZI")] {
        assert_eq!(conj(input, |w| w.z(0)), want, "z {input}");
    }
    // `conjCNOT_Xc` / `_Xt` / `_Zc` / `_Zt` and the forced sign `conjCNOT_YcYt`
    // (`Y_c Y_t ↦ −X_c ⊗ Z_t`).
    for (input, want) in [
        ("XI", "+XX"),
        ("IX", "+IX"),
        ("ZI", "+ZI"),
        ("IZ", "+ZZ"),
        ("YY", "-XZ"),
    ] {
        assert_eq!(conj(input, |w| w.cnot(0, 1)), want, "cnot {input}");
    }
    // `conjCZ_Xc` / `_Xt` / `_Zc` / `_Zt` and the forced sign `conjCZ_YcXt`
    // (`Y_c X_t ↦ −X_c ⊗ Y_t`).
    for (input, want) in [
        ("XI", "+XZ"),
        ("IX", "+ZX"),
        ("ZI", "+ZI"),
        ("IZ", "+IZ"),
        ("YX", "-XY"),
    ] {
        assert_eq!(conj(input, |w| w.cz(0, 1)), want, "cz {input}");
    }
}

#[test]
fn blanket_clifford_composes_its_primitives_in_the_order_the_lean_deltas_assume() {
    for w in all_words() {
        for q in 0..2 {
            let mut got = w;
            got.h(q);
            assert_eq!(got, oracle_h(w, q), "h({q}) on {w:?}");

            let mut got = w;
            got.s(q);
            assert_eq!(got, oracle_s(w, q), "s({q}) on {w:?}");

            let mut got = w;
            got.x(q);
            assert_eq!(got, oracle_pauli(w, w.z[q]), "x({q}) on {w:?}");

            let mut got = w;
            got.y(q);
            assert_eq!(got, oracle_pauli(w, w.x[q] ^ w.z[q]), "y({q}) on {w:?}");

            let mut got = w;
            got.z(q);
            assert_eq!(got, oracle_pauli(w, w.x[q]), "z({q}) on {w:?}");
        }
        // Both orientations of the two-qubit gates.
        for (c, t) in [(0usize, 1usize), (1, 0)] {
            let mut got = w;
            got.cnot(c, t);
            assert_eq!(got, oracle_cnot(w, c, t), "cnot({c},{t}) on {w:?}");

            let mut got = w;
            got.cz(c, t);
            assert_eq!(got, oracle_cz(w, c, t), "cz({c},{t}) on {w:?}");
        }
    }
}

#[test]
fn blanket_clifford_conjugation_keeps_the_phase_real() {
    // `conjH_isRealPhase` / `conjSdag_isRealPhase` / `conjCNOT_isRealPhase` /
    // `conjCZ_isRealPhase`: every generator's delta is `if … then 2 else 0`, so a
    // real input phase stays real. This is what makes the `±1` coefficient drain
    // in `ppvm-pauli-sum-2` total.
    /// One named generator of the blanket `Clifford`, applied to the stub.
    type NamedGate = (&'static str, fn(&mut Z4Word));

    let gates: [NamedGate; 7] = [
        ("h", |w| w.h(0)),
        ("s", |w| w.s(1)),
        ("x", |w| w.x(0)),
        ("y", |w| w.y(1)),
        ("z", |w| w.z(0)),
        ("cnot", |w| w.cnot(0, 1)),
        ("cz", |w| w.cz(1, 0)),
    ];
    for mut w in all_words() {
        if w.phase % 2 != 0 {
            continue; // start real (`±1`), as the propagation path always does
        }
        for (name, gate) in gates {
            let before = w;
            gate(&mut w);
            assert_eq!(w.phase % 2, 0, "{name} emitted a ±i phase from {before:?}");
        }
    }
}

// ---------------------------------------------------------------------------
// 2. Default bodies: the batched traits, Reset, and the stim noise aliases.
// ---------------------------------------------------------------------------

/// Records every required-method call so the defaults' *sequence* of side
/// effects — not merely their end state — can be diffed against the old bodies.
#[derive(Default)]
struct GateLog {
    log: Vec<String>,
}

impl Clifford for GateLog {
    fn x(&mut self, q: usize) {
        self.log.push(format!("x({q})"));
    }
    fn y(&mut self, q: usize) {
        self.log.push(format!("y({q})"));
    }
    fn z(&mut self, q: usize) {
        self.log.push(format!("z({q})"));
    }
    fn h(&mut self, q: usize) {
        self.log.push(format!("h({q})"));
    }
    fn s(&mut self, q: usize) {
        self.log.push(format!("s({q})"));
    }
    fn cnot(&mut self, c: usize, t: usize) {
        self.log.push(format!("cnot({c},{t})"));
    }
    fn cz(&mut self, a: usize, b: usize) {
        self.log.push(format!("cz({a},{b})"));
    }
}

impl CliffordExtensions for GateLog {
    fn s_dag(&mut self, q: usize) {
        self.log.push(format!("s_dag({q})"));
    }
    fn sqrt_x(&mut self, q: usize) {
        self.log.push(format!("sqrt_x({q})"));
    }
    fn sqrt_x_dag(&mut self, q: usize) {
        self.log.push(format!("sqrt_x_dag({q})"));
    }
    fn sqrt_y(&mut self, q: usize) {
        self.log.push(format!("sqrt_y({q})"));
    }
    fn sqrt_y_dag(&mut self, q: usize) {
        self.log.push(format!("sqrt_y_dag({q})"));
    }
    fn cy(&mut self, c: usize, t: usize) {
        self.log.push(format!("cy({c},{t})"));
    }
}

// The old crate's convention: a type with no specialization opts in with an
// empty impl and inherits the loop defaults.
impl CliffordBatch for GateLog {}
impl CliffordExtensionsBatch for GateLog {}

impl Reset for GateLog {
    fn reset<R: rand::Rng + ?Sized>(&mut self, q: usize, _rng: &mut R) {
        self.log.push(format!("reset({q})"));
    }
}

#[test]
fn clifford_stim_aliases_route_to_the_same_gate() {
    let mut g = GateLog::default();
    g.cx(0, 1);
    g.zcx(0, 1);
    g.zcz(2, 3);
    g.zcy(4, 5);
    assert_eq!(g.log, vec!["cnot(0,1)", "cnot(0,1)", "cz(2,3)", "cy(4,5)"]);
}

#[test]
fn clifford_batch_defaults_loop_in_order() {
    let mut g = GateLog::default();
    g.x_many(&[0, 2]);
    g.y_many(&[1]);
    g.z_many(&[3, 3]);
    g.h_many(&[0, 1]);
    g.s_many(&[4]);
    g.cnot_many(&[(0, 1), (2, 3)]);
    g.cz_many(&[(1, 2)]);
    assert_eq!(
        g.log,
        vec![
            "x(0)",
            "x(2)",
            "y(1)",
            "z(3)",
            "z(3)",
            "h(0)",
            "h(1)",
            "s(4)",
            "cnot(0,1)",
            "cnot(2,3)",
            "cz(1,2)",
        ]
    );

    // Empty index lists are no-ops.
    g.log.clear();
    g.h_many(&[]);
    g.cnot_many(&[]);
    assert!(g.log.is_empty());
}

#[test]
fn clifford_extensions_batch_defaults_loop_in_order() {
    let mut g = GateLog::default();
    g.s_dag_many(&[0, 1]);
    g.sqrt_x_many(&[2]);
    g.sqrt_x_dag_many(&[2]);
    g.sqrt_y_many(&[3]);
    g.sqrt_y_dag_many(&[3]);
    g.cy_many(&[(0, 1), (1, 0)]);
    assert_eq!(
        g.log,
        vec![
            "s_dag(0)",
            "s_dag(1)",
            "sqrt_x(2)",
            "sqrt_x_dag(2)",
            "sqrt_y(3)",
            "sqrt_y_dag(3)",
            "cy(0,1)",
            "cy(1,0)",
        ]
    );
}

/// `reset_x` is *defined* as `reset` then `h`, and `reset_y` as `reset`, `h`,
/// `s` — the old default bodies, reproduced call-for-call.
#[test]
fn reset_defaults_compose_reset_with_basis_change_cliffords() {
    let mut g = GateLog::default();
    let mut rng = SmallRng::seed_from_u64(0);
    g.reset(0, &mut rng);
    g.reset_z(1, &mut rng);
    g.reset_x(2, &mut rng);
    g.reset_y(3, &mut rng);
    assert_eq!(
        g.log,
        vec![
            "reset(0)", "reset(1)", "reset(2)", "h(2)", "reset(3)", "h(3)", "s(3)",
        ]
    );
}

#[test]
fn reset_batch_defaults_loop_in_order() {
    let mut g = GateLog::default();
    let mut rng = SmallRng::seed_from_u64(0);
    g.reset_many(&[0, 1], &mut rng);
    g.reset_z_many(&[2], &mut rng);
    g.reset_x_many(&[3], &mut rng);
    g.reset_y_many(&[4], &mut rng);
    assert_eq!(
        g.log,
        vec![
            "reset(0)", "reset(1)", "reset(2)", "reset(3)", "h(3)", "reset(4)", "h(4)", "s(4)",
        ]
    );
}

// ---------------------------------------------------------------------------
// 3. The noise surface: the stim `*_ERROR` aliases and the channel family.
// ---------------------------------------------------------------------------

#[derive(Default)]
struct NoiseLog {
    pauli: Vec<(usize, [f64; 3])>,
    other: Vec<String>,
}

impl PauliError<f64> for NoiseLog {
    fn pauli_error<R: rand::Rng + ?Sized>(
        &mut self,
        qubit: usize,
        probabilities: [f64; 3],
        _rng: &mut R,
    ) {
        self.pauli.push((qubit, probabilities));
    }
}

impl PauliErrorAll<f64> for NoiseLog {
    fn pauli_error_all<R: rand::Rng + ?Sized>(&mut self, p: [f64; 3], _rng: &mut R) {
        self.other.push(format!("all{p:?}"));
    }
}

impl TwoQubitPauliError<f64> for NoiseLog {
    fn two_qubit_pauli_error<R: rand::Rng + ?Sized>(
        &mut self,
        q0: usize,
        q1: usize,
        p: [f64; 15],
        _rng: &mut R,
    ) {
        self.other.push(format!("two({q0},{q1},{})", p[0]));
    }
}

impl Depolarizing<f64> for NoiseLog {
    fn depolarize1<R: rand::Rng + ?Sized>(&mut self, qubit: usize, p: f64, _rng: &mut R) {
        self.other.push(format!("dep1({qubit},{p})"));
    }
}

impl Depolarizing2<f64> for NoiseLog {
    fn depolarize2<R: rand::Rng + ?Sized>(&mut self, q0: usize, q1: usize, p: f64, _rng: &mut R) {
        self.other.push(format!("dep2({q0},{q1},{p})"));
    }
}

impl AmplitudeDamping<f64> for NoiseLog {
    fn amplitude_damping(&mut self, qubit: usize, gamma: f64) {
        self.other.push(format!("ad({qubit},{gamma})"));
    }
}

impl LossChannel<f64> for NoiseLog {
    fn loss_channel<R: rand::Rng + ?Sized>(&mut self, qubit: usize, p: f64, _rng: &mut R) {
        self.other.push(format!("loss({qubit},{p})"));
    }
}

impl CorrelatedLossChannel<f64> for NoiseLog {
    fn correlated_loss_channel<R: rand::Rng + ?Sized>(
        &mut self,
        q0: usize,
        q1: usize,
        p: [f64; 3],
        _rng: &mut R,
    ) {
        self.other.push(format!("corr({q0},{q1},{p:?})"));
    }
}

impl ResetLossChannel for NoiseLog {
    fn reset_loss_channel(&mut self, qubit: usize) {
        self.other.push(format!("reset_loss({qubit})"));
    }
}

impl AsymmetricLossChannel<f64> for NoiseLog {
    fn asymmetric_loss_channel<R: rand::Rng + ?Sized>(
        &mut self,
        qubit: usize,
        p0: f64,
        p1: f64,
        _rng: &mut R,
    ) {
        self.other.push(format!("asym({qubit},{p0},{p1})"));
    }
}

/// The stim aliases are *defined* as one-hot Pauli channels — the old default
/// bodies, including which slot each probability lands in.
#[test]
fn pauli_error_stim_aliases_are_one_hot_channels() {
    let mut n = NoiseLog::default();
    let mut rng = SmallRng::seed_from_u64(0);
    n.x_error(0, 0.1, &mut rng);
    n.y_error(1, 0.2, &mut rng);
    n.z_error(2, 0.3, &mut rng);
    assert_eq!(
        n.pauli,
        vec![
            (0, [0.1, 0.0, 0.0]),
            (1, [0.0, 0.2, 0.0]),
            (2, [0.0, 0.0, 0.3]),
        ]
    );
}

#[test]
fn pauli_error_batch_defaults_loop_in_order() {
    let mut n = NoiseLog::default();
    let mut rng = SmallRng::seed_from_u64(0);
    n.pauli_error_many(&[0, 1], [0.1, 0.2, 0.3], &mut rng);
    n.x_error_many(&[2, 3], 0.4, &mut rng);
    n.y_error_many(&[4], 0.5, &mut rng);
    n.z_error_many(&[5], 0.6, &mut rng);
    assert_eq!(
        n.pauli,
        vec![
            (0, [0.1, 0.2, 0.3]),
            (1, [0.1, 0.2, 0.3]),
            (2, [0.4, 0.0, 0.0]),
            (3, [0.4, 0.0, 0.0]),
            (4, [0.0, 0.5, 0.0]),
            (5, [0.0, 0.0, 0.6]),
        ]
    );
}

#[test]
fn channel_family_is_callable_and_batches_in_order() {
    let mut n = NoiseLog::default();
    let mut rng = SmallRng::seed_from_u64(0);
    n.pauli_error_all([0.01, 0.02, 0.03], &mut rng);
    n.two_qubit_pauli_error_many(&[(0, 1), (2, 3)], [0.5; 15], &mut rng);
    n.depolarize1_many(&[0, 1], 0.1, &mut rng);
    n.depolarize2_many(&[(0, 1)], 0.2, &mut rng);
    n.amplitude_damping(3, 0.05);
    n.loss_channel(4, 0.01, &mut rng);
    n.correlated_loss_channel(0, 1, [0.1, 0.2, 0.3], &mut rng);
    n.reset_loss_channel(4);
    n.asymmetric_loss_channel(5, 0.01, 0.02, &mut rng);
    assert_eq!(
        n.other,
        vec![
            "all[0.01, 0.02, 0.03]",
            "two(0,1,0.5)",
            "two(2,3,0.5)",
            "dep1(0,0.1)",
            "dep1(1,0.1)",
            "dep2(0,1,0.2)",
            "ad(3,0.05)",
            "loss(4,0.01)",
            "corr(0,1,[0.1, 0.2, 0.3])",
            "reset_loss(4)",
            "asym(5,0.01,0.02)",
        ]
    );
}
