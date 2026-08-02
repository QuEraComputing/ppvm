// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Lean-oracle property tests for `ppvm-tableau-2`.
//!
//! A tableau row **is** a Pauli word with a `ℤ/4` sign, so the machine-checked
//! per-generator semantics already cover the row action:
//!
//! * `lean/PPVM/Pauli/Symplectic.lean` — the Clifford bit maps `hAct`, `sAct`,
//!   `cnotAct`, `czAct`, `cyAct` are `Sp(2n, 2)` **isometries**
//!   (`*_isometry`), **involutive** (`*Act_involutive`) and **bijective**
//!   (`*_bijective`); the loss-guarded `hActL`/`sActL`/`cnotActL`/`czActL`
//!   preserve `LossInv` and reduce to the unguarded map on present qubits.
//! * `lean/PPVM/Pauli/Conjugation.lean` — the exact `ℤ/4` sign deltas:
//!   `conjH_sign`/`conjS_sign` (`x∧z`), `conjSdag_sign` (`x∧¬z`),
//!   `conjCNOT_sign` (`x_c ∧ z_t ∧ (x_t = z_c)`), `conjCZ_sign`
//!   (`x_c ∧ x_t ∧ (z_c ≠ z_t)`), and `conjX`/`conjY`/`conjZ`.
//! * `lean/PPVM/Tableau/Frame.lean` — `IsSymplecticFrame`,
//!   `isSymplecticFrame_identity`, `isSymplecticFrame_hAct`/`sAct`/`cnotAct`/
//!   `czAct`, `frame_linearIndependent`, `measurement_dichotomy`,
//!   `measure_deterministic_iff_xfree`, and — for the one non-unitary frame
//!   mutation — `isSymplecticFrame_projectFrame` (`projectFrame`,
//!   `rowUpdate_eq_ite`).
//!
//! This file reproduces those statements as executable Rust — exhaustively over
//! the finite single-/two-qubit phased-Pauli groups where Lean uses `decide`,
//! randomized over `n` qubits where it uses a general argument — including the
//! **tableau-level** invariant that the `2n` rows remain a symplectic basis after
//! every gate *and* after every measurement. Both halves are now machine-checked:
//! the unitary half by `isSymplecticFrame_*` (via `IsSymplecticFrame.map`), the
//! measurement projection by `isSymplecticFrame_projectFrame`.

use ppvm_conformance_2::seeded_rng;
use ppvm_conformance_2::tableau::*;
use ppvm_traits_2::Pauli as Pauli2;
use rand::RngExt;

// ===========================================================================
// The Lean reference semantics, transcribed
// ===========================================================================

/// A one-site phased Pauli: `i^phase · X^x Z^z` (`PPVM.PauliPhase.PhasedPauli`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct P1 {
    phase: u8,
    x: bool,
    z: bool,
}

/// `Conjugation.lean::conjH` — `⟨phase + 2·(x∧z), z, x⟩`.
fn conj_h(p: P1) -> P1 {
    P1 {
        phase: (p.phase + if p.x && p.z { 2 } else { 0 }) % 4,
        x: p.z,
        z: p.x,
    }
}

/// `Conjugation.lean::conjS` — `⟨phase + 2·(x∧z), x, x⊕z⟩`.
fn conj_s(p: P1) -> P1 {
    P1 {
        phase: (p.phase + if p.x && p.z { 2 } else { 0 }) % 4,
        x: p.x,
        z: p.x ^ p.z,
    }
}

/// `Conjugation.lean::conjSdag` — `⟨phase + 2·(x∧¬z), x, x⊕z⟩`.
fn conj_s_dag(p: P1) -> P1 {
    P1 {
        phase: (p.phase + if p.x && !p.z { 2 } else { 0 }) % 4,
        x: p.x,
        z: p.x ^ p.z,
    }
}

/// `Conjugation.lean::conjX` — word fixed, phase `+2` iff `z`.
fn conj_x(p: P1) -> P1 {
    P1 {
        phase: (p.phase + if p.z { 2 } else { 0 }) % 4,
        ..p
    }
}

/// `Conjugation.lean::conjY` — word fixed, phase `+2` iff `x⊕z`.
fn conj_y(p: P1) -> P1 {
    P1 {
        phase: (p.phase + if p.x ^ p.z { 2 } else { 0 }) % 4,
        ..p
    }
}

/// `Conjugation.lean::conjZ` — word fixed, phase `+2` iff `x`.
fn conj_z(p: P1) -> P1 {
    P1 {
        phase: (p.phase + if p.x { 2 } else { 0 }) % 4,
        ..p
    }
}

/// A two-site phased Pauli (`PPVM.PauliPhase.TwoPauli`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct P2 {
    phase: u8,
    xc: bool,
    zc: bool,
    xt: bool,
    zt: bool,
}

/// `Conjugation.lean::conjCNOT` — bits `x_t ⊕= x_c`, `z_c ⊕= z_t`; sign
/// `cnotDelta = x_c ∧ z_t ∧ (x_t = z_c)`.
fn conj_cnot(p: P2) -> P2 {
    P2 {
        phase: (p.phase + if p.xc && p.zt && (p.xt == p.zc) { 2 } else { 0 }) % 4,
        xc: p.xc,
        zc: p.zc ^ p.zt,
        xt: p.xt ^ p.xc,
        zt: p.zt,
    }
}

/// `Conjugation.lean::conjCZ` — bits `z_c ⊕= x_t`, `z_t ⊕= x_c`; sign
/// `czDelta = x_c ∧ x_t ∧ (z_c ≠ z_t)`.
fn conj_cz(p: P2) -> P2 {
    P2 {
        phase: (p.phase + if p.xc && p.xt && (p.zc ^ p.zt) { 2 } else { 0 }) % 4,
        xc: p.xc,
        zc: p.zc ^ p.xt,
        xt: p.xt,
        zt: p.zt ^ p.xc,
    }
}

/// Every phased one-site Pauli (16 elements) — the `decide` domain.
fn all_p1() -> Vec<P1> {
    let mut v = Vec::new();
    for phase in 0..4u8 {
        for x in [false, true] {
            for z in [false, true] {
                v.push(P1 { phase, x, z });
            }
        }
    }
    v
}

/// Every phased two-site Pauli (64 elements).
fn all_p2() -> Vec<P2> {
    let mut v = Vec::new();
    for phase in 0..4u8 {
        for xc in [false, true] {
            for zc in [false, true] {
                for xt in [false, true] {
                    for zt in [false, true] {
                        v.push(P2 {
                            phase,
                            xc,
                            zc,
                            xt,
                            zt,
                        });
                    }
                }
            }
        }
    }
    v
}

// ===========================================================================
// Symplectic.lean — bijectivity, involutivity, isometry (exhaustive)
// ===========================================================================

/// `hAct_bijective`, `sAct_bijective`, `cnotAct_bijective`, `czAct_bijective`:
/// each conjugation is injective on the finite phased-Pauli group, hence an
/// automorphism.
#[test]
fn lean_conjugations_are_bijections() {
    for (name, f) in [
        ("conjH", conj_h as fn(P1) -> P1),
        ("conjS", conj_s),
        ("conjSdag", conj_s_dag),
        ("conjX", conj_x),
        ("conjY", conj_y),
        ("conjZ", conj_z),
    ] {
        let mut images: Vec<P1> = all_p1().into_iter().map(f).collect();
        images.sort_by_key(|p| (p.phase, p.x, p.z));
        images.dedup();
        assert_eq!(images.len(), 16, "{name} is not injective on 𝒫₁");
    }

    for (name, f) in [("conjCNOT", conj_cnot as fn(P2) -> P2), ("conjCZ", conj_cz)] {
        let mut images: Vec<P2> = all_p2().into_iter().map(f).collect();
        images.sort_by_key(|p| (p.phase, p.xc, p.zc, p.xt, p.zt));
        images.dedup();
        assert_eq!(images.len(), 64, "{name} is not injective on 𝒫₂");
    }
}

/// `hAct_involutive`, `cnotAct_involutive`, `czAct_involutive` (and
/// `conjH_involutive` / `conjCNOT_involutive` / `conjCZ_involutive`, which are
/// involutive **including** the sign); `sAct_involutive` holds on the bits only
/// — `conjS` has order 4 (`conjS_iterate_four`, `conjS_iterate_two_ne_id`).
#[test]
fn lean_conjugations_are_involutive_where_lean_says_so() {
    for p in all_p1() {
        assert_eq!(conj_h(conj_h(p)), p, "conjH² ≠ id at {p:?}");
        // sAct (bits) is involutive; conjS (with sign) has order 4.
        let ss = conj_s(conj_s(p));
        assert_eq!((ss.x, ss.z), (p.x, p.z), "sAct² ≠ id on bits at {p:?}");
        assert_eq!(conj_s(conj_s(conj_s(conj_s(p)))), p, "conjS⁴ ≠ id at {p:?}");
        assert_eq!(conj_s_dag(conj_s(p)), p, "conjSdag ∘ conjS ≠ id at {p:?}");
        assert_eq!(conj_s(conj_s_dag(p)), p, "conjS ∘ conjSdag ≠ id at {p:?}");
    }
    assert!(
        all_p1().into_iter().any(|p| conj_s(conj_s(p)) != p),
        "conjS² must NOT be the identity (conjS_iterate_two_ne_id)"
    );
    for p in all_p2() {
        assert_eq!(conj_cnot(conj_cnot(p)), p, "conjCNOT² ≠ id at {p:?}");
        assert_eq!(conj_cz(conj_cz(p)), p, "conjCZ² ≠ id at {p:?}");
    }
}

/// `hAct_isometry`, `sAct_isometry`, `cnotAct_isometry`, `czAct_isometry`: the
/// symplectic form `ω(p, q) = x_p·z_q ⊕ z_p·x_q` is preserved.
#[test]
fn lean_conjugations_are_symplectic_isometries() {
    let omega1 = |a: P1, b: P1| (a.x & b.z) ^ (a.z & b.x);
    for f in [conj_h as fn(P1) -> P1, conj_s, conj_s_dag] {
        for a in all_p1() {
            for b in all_p1() {
                assert_eq!(
                    omega1(f(a), f(b)),
                    omega1(a, b),
                    "single-site isometry failed at {a:?}, {b:?}"
                );
            }
        }
    }
    let omega2 = |a: P2, b: P2| ((a.xc & b.zc) ^ (a.zc & b.xc)) ^ ((a.xt & b.zt) ^ (a.zt & b.xt));
    for f in [conj_cnot as fn(P2) -> P2, conj_cz] {
        for a in all_p2() {
            for b in all_p2() {
                assert_eq!(
                    omega2(f(a), f(b)),
                    omega2(a, b),
                    "two-site isometry failed at {a:?}, {b:?}"
                );
            }
        }
    }
}

/// `conjH_X`/`conjH_Z`/`conjH_Y`, `conjS_*`, `conjSdag_*`, `conjCNOT_*`,
/// `conjCZ_*`: the generator tables Lean pins with `decide`.
#[test]
fn lean_generator_tables() {
    let px = P1 {
        phase: 0,
        x: true,
        z: false,
    };
    let pz = P1 {
        phase: 0,
        x: false,
        z: true,
    };
    let py = P1 {
        phase: 0,
        x: true,
        z: true,
    };
    assert_eq!(conj_h(px), pz, "conjH_X");
    assert_eq!(conj_h(pz), px, "conjH_Z");
    assert_eq!(conj_h(py), P1 { phase: 2, ..py }, "conjH_Y");
    assert_eq!(conj_s(px), py, "conjS_X");
    assert_eq!(conj_s(pz), pz, "conjS_Z");
    assert_eq!(conj_s(py), P1 { phase: 2, ..px }, "conjS_Y");
    assert_eq!(conj_s_dag(px), P1 { phase: 2, ..py }, "conjSdag_X");
    assert_eq!(conj_s_dag(py), px, "conjSdag_Y");
    assert_eq!(conj_s_dag(pz), pz, "conjSdag_Z");

    let two = |xc, zc, xt, zt| P2 {
        phase: 0,
        xc,
        zc,
        xt,
        zt,
    };
    // CNOT: X_c → X_cX_t, X_t → X_t, Z_c → Z_c, Z_t → Z_cZ_t.
    assert_eq!(
        conj_cnot(two(true, false, false, false)),
        two(true, false, true, false),
        "conjCNOT_Xc"
    );
    assert_eq!(
        conj_cnot(two(false, false, true, false)),
        two(false, false, true, false),
        "conjCNOT_Xt"
    );
    assert_eq!(
        conj_cnot(two(false, true, false, false)),
        two(false, true, false, false),
        "conjCNOT_Zc"
    );
    assert_eq!(
        conj_cnot(two(false, false, false, true)),
        two(false, true, false, true),
        "conjCNOT_Zt"
    );
    assert_eq!(
        conj_cnot(two(true, true, true, true)),
        P2 {
            phase: 2,
            ..two(true, false, false, true)
        },
        "conjCNOT_YcYt"
    );
    // CZ: X_c → X_cZ_t, X_t → Z_cX_t, Z_c → Z_c, Z_t → Z_t.
    assert_eq!(
        conj_cz(two(true, false, false, false)),
        two(true, false, false, true),
        "conjCZ_Xc"
    );
    assert_eq!(
        conj_cz(two(false, false, true, false)),
        two(false, true, true, false),
        "conjCZ_Xt"
    );
    assert_eq!(
        conj_cz(two(false, true, false, false)),
        two(false, true, false, false),
        "conjCZ_Zc"
    );
    assert_eq!(
        conj_cz(two(false, false, false, true)),
        two(false, false, false, true),
        "conjCZ_Zt"
    );
}

// ===========================================================================
// The impl vs the Lean reference: every row of a randomized frame
// ===========================================================================

/// Read site `q` of a `(x-plane, z-plane, phase)` snapshot. The planes come from
/// the `[usize; 2]` (64-bit-element) configuration, so 64 bits per element.
fn site(row: &RowSnapshot, q: usize) -> (bool, bool) {
    let (w, off) = (q / 64, q % 64);
    ((row.0[w] >> off) & 1 == 1, (row.1[w] >> off) & 1 == 1)
}

/// A randomized frame on `n` qubits (Clifford only — the frame is what the Lean
/// spec talks about).
fn random_frame(n: usize, seed: u64, len: usize) -> NewWide {
    let mut rng = seeded_rng(seed);
    let mut tab: NewWide = Driver::new_seeded(n, 1e-12, seed);
    for _ in 0..len {
        let q = rng.random_range(0..n);
        let mut b = rng.random_range(0..n);
        while b == q {
            b = rng.random_range(0..n);
        }
        match rng.random_range(0..8usize) {
            0 => tab.h(q),
            1 => tab.s(q),
            2 => tab.x(q),
            3 => tab.y(q),
            4 => tab.z(q),
            5 => tab.sqrt_y(q),
            6 => tab.cnot(q, b),
            _ => tab.cz(q, b),
        }
    }
    tab
}

/// Every single-qubit generator of the tableau conjugates each of the `2n` rows
/// **exactly** by the Lean sign rule — bits and `ℤ/4` phase.
#[test]
fn tableau_single_qubit_rows_follow_the_lean_conjugation() {
    let n = 6;
    for seed in 0..10u64 {
        let base = random_frame(n, seed, 30);
        for q in 0..n {
            for (name, gate, reference) in [
                (
                    "h",
                    NewWide::h as fn(&mut NewWide, usize),
                    conj_h as fn(P1) -> P1,
                ),
                ("s", NewWide::s, conj_s),
                ("x", NewWide::x, conj_x),
                ("y", NewWide::y, conj_y),
                ("z", NewWide::z, conj_z),
            ] {
                let before = base.rows();
                let mut t = base.fork(Some(0));
                gate(&mut t, q);
                let after = t.rows();
                for (i, (rb, ra)) in before.iter().zip(after.iter()).enumerate() {
                    let (x, z) = site(rb, q);
                    let want = reference(P1 { phase: rb.2, x, z });
                    let (xa, za) = site(ra, q);
                    assert_eq!(
                        (xa, za, ra.2),
                        (want.x, want.z, want.phase),
                        "seed {seed}: {name}({q}) row {i}: bits/phase diverge from the Lean rule"
                    );
                    // Untouched sites must be bit-identical.
                    for other in (0..n).filter(|&o| o != q) {
                        assert_eq!(
                            site(rb, other),
                            site(ra, other),
                            "seed {seed}: {name}({q}) disturbed site {other} of row {i}"
                        );
                    }
                }
            }
        }
    }
}

/// `cnot` and `cz` conjugate each row exactly by `conjCNOT` / `conjCZ`,
/// including the `ℤ/4` sign deltas `conjCNOT_sign` / `conjCZ_sign`.
#[test]
fn tableau_two_qubit_rows_follow_the_lean_conjugation() {
    let n = 6;
    for seed in 0..10u64 {
        let base = random_frame(n, seed, 30);
        for c in 0..n {
            for t_q in 0..n {
                if c == t_q {
                    continue;
                }
                for (name, gate, reference) in [
                    (
                        "cnot",
                        NewWide::cnot as fn(&mut NewWide, usize, usize),
                        conj_cnot as fn(P2) -> P2,
                    ),
                    ("cz", NewWide::cz, conj_cz),
                ] {
                    let before = base.rows();
                    let mut tab = base.fork(Some(0));
                    gate(&mut tab, c, t_q);
                    let after = tab.rows();
                    for (i, (rb, ra)) in before.iter().zip(after.iter()).enumerate() {
                        let (xc, zc) = site(rb, c);
                        let (xt, zt) = site(rb, t_q);
                        let want = reference(P2 {
                            phase: rb.2,
                            xc,
                            zc,
                            xt,
                            zt,
                        });
                        let (xca, zca) = site(ra, c);
                        let (xta, zta) = site(ra, t_q);
                        assert_eq!(
                            (xca, zca, xta, zta, ra.2),
                            (want.xc, want.zc, want.xt, want.zt, want.phase),
                            "seed {seed}: {name}({c},{t_q}) row {i} diverges from the Lean rule"
                        );
                    }
                }
            }
        }
    }
}

// ===========================================================================
// Frame.lean — the symplectic-frame invariant (the tableau-level statement)
// ===========================================================================

/// `ω(p, q) = Σ_k (x_p[k]·z_q[k] ⊕ z_p[k]·x_q[k])` over all sites.
fn omega(a: &RowSnapshot, b: &RowSnapshot) -> bool {
    let mut acc = 0u32;
    for k in 0..a.0.len() {
        acc += (a.0[k] & b.1[k]).count_ones() + (a.1[k] & b.0[k]).count_ones();
    }
    acc % 2 == 1
}

/// `IsSymplecticFrame`: `ω(dᵢ, dⱼ) = 0`, `ω(sᵢ, sⱼ) = 0`, `ω(dᵢ, sⱼ) = δᵢⱼ`.
#[track_caller]
fn assert_symplectic_frame<D: Driver>(tab: &D, ctx: &str) {
    let rows = tab.rows();
    let n = tab.n_qubits();
    assert_eq!(rows.len(), 2 * n);
    for i in 0..n {
        for j in 0..n {
            assert!(
                !omega(&rows[i], &rows[j]),
                "{ctx}: ω(d{i}, d{j}) ≠ 0 — destabilizers must commute"
            );
            assert!(
                !omega(&rows[n + i], &rows[n + j]),
                "{ctx}: ω(s{i}, s{j}) ≠ 0 — stabilizers must commute"
            );
            assert_eq!(
                omega(&rows[i], &rows[n + j]),
                i == j,
                "{ctx}: ω(d{i}, s{j}) ≠ δ — the destabilizer pairing broke"
            );
        }
    }
}

/// `isSymplecticFrame_identity`: the fresh `|0…0⟩` frame is symplectic.
#[test]
fn identity_frame_is_symplectic() {
    for n in [1usize, 2, 5, 17, 64, 65, 85] {
        let tab: NewWide = Driver::new_seeded(n, 1e-12, 0);
        assert_symplectic_frame(&tab, &format!("identity frame n={n}"));
    }
}

/// `isSymplecticFrame_hAct`/`sAct`/`cnotAct`/`czAct` lifted to the whole frame:
/// **every** gate on the surface (including the derived extensions, the batched
/// sweeps and the fused `cz_block` kernels) leaves the frame symplectic.
#[test]
fn every_gate_preserves_the_symplectic_frame() {
    let n = 70; // > 64: exercises the cross-word bit addressing too
    let mut rng = seeded_rng(4242);
    let mut tab: NewWide = Driver::new_seeded(n, 1e-12, 4242);
    assert_symplectic_frame(&tab, "start");
    for step in 0..300 {
        let q = rng.random_range(0..n);
        let mut b = rng.random_range(0..n);
        while b == q {
            b = rng.random_range(0..n);
        }
        let idx: Vec<usize> = (0..n).step_by(1 + rng.random_range(0..5usize)).collect();
        let pairs: Vec<(usize, usize)> = (0..10).map(|i| (i, i + 20)).collect();
        match rng.random_range(0..18usize) {
            0 => tab.x(q),
            1 => tab.y(q),
            2 => tab.z(q),
            3 => tab.h(q),
            4 => tab.s(q),
            5 => tab.s_dag(q),
            6 => tab.sqrt_x(q),
            7 => tab.sqrt_x_dag(q),
            8 => tab.sqrt_y(q),
            9 => tab.sqrt_y_dag(q),
            10 => tab.cnot(q, b),
            11 => tab.cz(q, b),
            12 => tab.cy(q, b),
            13 => tab.zcx(q, b),
            14 => tab.zcy(q, b),
            15 => tab.h_many(&idx),
            16 => tab.cz_many(&pairs),
            _ => tab.cz_block(0, 20, 10),
        }
        assert_symplectic_frame(&tab, &format!("after step {step}"));
    }
}

/// `Frame.lean::isSymplecticFrame_projectFrame`: the measurement
/// projection (`update_tableau_according_to_outcome`) also preserves the
/// symplectic frame — in both the case-a (random) and case-b (deterministic)
/// branches, on a state with real amplitude branching.
#[test]
fn measurement_preserves_the_symplectic_frame() {
    let n = 12;
    for seed in 0..16u64 {
        let mut rng = seeded_rng(seed);
        let mut tab: NewWide = Driver::new_seeded(n, 1e-12, seed);
        for _ in 0..40 {
            let q = rng.random_range(0..n);
            let mut b = rng.random_range(0..n);
            while b == q {
                b = rng.random_range(0..n);
            }
            match rng.random_range(0..6usize) {
                0 => tab.h(q),
                1 => tab.s(q),
                2 => tab.cnot(q, b),
                3 => tab.cz(q, b),
                4 => tab.t(q),
                _ => tab.sqrt_y(q),
            }
        }
        assert_symplectic_frame(&tab, &format!("seed {seed}: pre-measure"));
        for q in 0..n {
            tab.measure(q);
            assert_symplectic_frame(&tab, &format!("seed {seed}: after measure({q})"));
        }
    }
}

/// The whole 85-qubit MSD workload — the deepest real circuit — keeps the frame
/// symplectic from construction through all 85 measurements.
#[test]
fn msd_workload_preserves_the_symplectic_frame() {
    let mut tab: NewWide = msd_state(Some(1));
    assert_symplectic_frame(&tab, "MSD: after the Clifford+T portion");
    for q in 0..MSD_QUBITS {
        tab.measure(q);
        assert_symplectic_frame(&tab, &format!("MSD: after measure({q})"));
    }
}

/// `frame_linearIndependent`: the `2n` rows are linearly independent over
/// `𝔽₂`. Checked directly by Gaussian elimination on the `2n × 2n` bit matrix.
#[test]
fn frame_rows_are_linearly_independent() {
    for seed in 0..8u64 {
        let n = 12;
        let tab = random_frame(n, seed, 60);
        let rows = tab.rows();
        // Each row → a 2n-bit vector (x-plane ‖ z-plane) as a u128.
        let mut m: Vec<u128> = rows
            .iter()
            .map(|r| {
                let mut v = 0u128;
                for q in 0..n {
                    let (x, z) = site(r, q);
                    if x {
                        v |= 1u128 << q;
                    }
                    if z {
                        v |= 1u128 << (n + q);
                    }
                }
                v
            })
            .collect();
        // Gaussian elimination: rank must be exactly 2n.
        let mut rank = 0usize;
        for bit in 0..(2 * n) {
            let Some(pivot) = (rank..m.len()).find(|&r| (m[r] >> bit) & 1 == 1) else {
                continue;
            };
            m.swap(rank, pivot);
            for r in 0..m.len() {
                if r != rank && (m[r] >> bit) & 1 == 1 {
                    m[r] ^= m[rank];
                }
            }
            rank += 1;
        }
        assert_eq!(rank, 2 * n, "seed {seed}: frame rows are not independent");
    }
}

// ===========================================================================
// Frame.lean — the measurement dichotomy
// ===========================================================================

/// `measurement_dichotomy` / `measure_deterministic_iff_xfree`: the outcome is
/// deterministic **exactly** when no stabilizer anticommutes with `Z_q`. The
/// observable consequences: a deterministic measurement consumes no randomness
/// and is idempotent, and a repeated measurement always reproduces its own
/// outcome.
#[test]
fn measurement_dichotomy_holds() {
    let n = 8;
    for seed in 0..24u64 {
        let mut rng = seeded_rng(seed);
        let mut tab: NewWide = Driver::new_seeded(n, 1e-12, seed);
        for _ in 0..30 {
            let q = rng.random_range(0..n);
            let mut b = rng.random_range(0..n);
            while b == q {
                b = rng.random_range(0..n);
            }
            match rng.random_range(0..4usize) {
                0 => tab.h(q),
                1 => tab.s(q),
                2 => tab.cnot(q, b),
                _ => tab.cz(q, b),
            }
        }
        for q in 0..n {
            // Determine the branch from the frame itself: case b ⟺ no stabilizer
            // has an X-bit at `q`.
            let rows = tab.rows();
            let deterministic = (n..2 * n).all(|i| !site(&rows[i], q).0);

            let mut probe = tab.fork(Some(seed + 1));
            let first = probe.measure(q).unwrap();
            if deterministic {
                // No RNG draw: the stream is untouched, so a differently-seeded
                // fork gives the same answer.
                let mut other = tab.fork(Some(seed + 9999));
                assert_eq!(
                    other.measure(q).unwrap(),
                    first,
                    "seed {seed}: deterministic measure({q}) depended on the RNG"
                );
            }
            // Idempotence: measuring again always reproduces the outcome, and
            // the second measurement is deterministic by construction.
            let second = probe.measure(q).unwrap();
            assert_eq!(first, second, "seed {seed}: measure({q}) is not idempotent");
            let rows2 = probe.rows();
            assert!(
                (n..2 * n).all(|i| !site(&rows2[i], q).0),
                "seed {seed}: after measuring {q}, Z_{q} is not a stabilizer"
            );
        }
    }
}

// ===========================================================================
// Symplectic.lean — the loss-guarded actions
// ===========================================================================

/// `hActL`/`sActL`/`cnotActL`/`czActL` and `*_preserves_loss`: a gate whose
/// target is lost is the **identity** on the frame, and a two-qubit gate is the
/// identity when EITHER endpoint is lost (`cnotActL_lost_target_stays_identity`).
#[test]
fn loss_guarded_actions_are_the_identity_on_lost_qubits() {
    let n = 6;
    let mut tab: NewWide = Driver::new_seeded(n, 1e-12, 3);
    for q in 0..n {
        tab.h(q);
    }
    tab.cnot(0, 1);
    tab.loss_channel(2, 1.0);
    assert!(tab.lost()[2]);

    let before = tab.rows();
    for f in [
        NewWide::h as fn(&mut NewWide, usize),
        NewWide::s,
        NewWide::x,
        NewWide::y,
        NewWide::z,
        NewWide::sqrt_x,
        NewWide::sqrt_y,
    ] {
        let mut t = tab.fork(Some(0));
        f(&mut t, 2);
        assert_eq!(
            t.rows(),
            before,
            "a single-qubit gate fired on a lost qubit"
        );
    }
    for f in [
        NewWide::cnot as fn(&mut NewWide, usize, usize),
        NewWide::cz,
        NewWide::cy,
    ] {
        for &(a, b) in &[(2usize, 3usize), (3, 2)] {
            let mut t = tab.fork(Some(0));
            f(&mut t, a, b);
            assert_eq!(
                t.rows(),
                before,
                "a two-qubit gate fired with a lost endpoint ({a},{b})"
            );
        }
    }

    // ...and the guard never touches a present pair.
    let mut t = tab.fork(Some(0));
    t.cnot(0, 1);
    assert_ne!(t.rows(), before, "the guard blocked a fully-present pair");
}

// ===========================================================================
// Bitstring.lean — the XOR relabel is a bijection
// ===========================================================================

/// The branch relabel `idx ↦ idx ⊕ stab_anticomm_bits` is a bijection, which is
/// why the branch stream can be pushed with `unsafe_insert` (no dedup probe) and
/// why a per-key coalesce could never merge anything on the relabel path.
///
/// Observable consequence, asserted here on the real engine: a `T` gate on a
/// **fresh** qubit exactly doubles the support (no collisions), while a `T` on an
/// already-branched qubit maps the support onto itself (all collisions) — and in
/// both regimes every stored index stays distinct.
#[test]
fn xor_relabel_is_a_bijection() {
    // Exhaustive over the relabel map itself.
    for stab in 0u32..256 {
        let mut images: Vec<u32> = (0u32..256).map(|i| i ^ stab).collect();
        images.sort_unstable();
        images.dedup();
        assert_eq!(images.len(), 256, "XOR by {stab} is not a bijection");
    }

    // On the engine: doubling and merging regimes, with distinct keys throughout.
    for j in [1usize, 3, 6, 9] {
        let mut doubling: NewWide = branch_grow(j);
        assert_eq!(doubling.n_coeffs(), 1 << j);
        doubling.h(j);
        doubling.t(j);
        assert_eq!(
            doubling.n_coeffs(),
            1 << (j + 1),
            "a T on a fresh qubit must double the support"
        );
        let mut keys: Vec<u128> = doubling.coeffs().iter().map(|e| e.0).collect();
        keys.sort_unstable();
        let len = keys.len();
        keys.dedup();
        assert_eq!(keys.len(), len, "duplicate index after the relabel");

        let mut merging: NewWide = branch_grow(j);
        merging.t(0);
        assert_eq!(
            merging.n_coeffs(),
            1 << j,
            "a T on an already-branched qubit must map the support onto itself"
        );
        let mut keys2: Vec<u128> = merging.coeffs().iter().map(|e| e.0).collect();
        keys2.sort_unstable();
        let len2 = keys2.len();
        keys2.dedup();
        assert_eq!(keys2.len(), len2, "duplicate index after the merge");
    }
}

// ===========================================================================
// Rotation.lean — the staged branch is a norm-preserving 2-D rotation
// ===========================================================================

/// `lean/PPVM/Instantiations/Rotation.lean` (`rot_norm_sq`, `rot_rot`): each
/// branching rotation is a norm-preserving, angle-additive rotation on the
/// coefficient pair. Observable: with truncation off, `rz(θ₁); rz(θ₂)` equals
/// `rz(θ₁+θ₂)` and the state norm is unchanged (gates never normalize, so the
/// norm staying 1 IS the norm-preservation statement).
#[test]
fn rotations_are_norm_preserving_and_angle_additive() {
    let mut rng = seeded_rng(77);
    for _ in 0..40 {
        let t1: f64 = rng.random_range(-3.0..3.0);
        let t2: f64 = rng.random_range(-3.0..3.0);
        let mut split: NewWide = Driver::new_seeded(3, 0.0, 0);
        let mut joint: NewWide = Driver::new_seeded(3, 0.0, 0);
        for t in [&mut split, &mut joint] {
            t.h(0);
            t.h(1);
            t.cnot(0, 2);
        }
        split.rz(1, t1);
        split.rz(1, t2);
        joint.rz(1, t1 + t2);

        let a = split.coeffs_sorted();
        let b = joint.coeffs_sorted();
        assert_eq!(a.len(), b.len(), "rot_rot: support size");
        for (x, y) in a.iter().zip(b.iter()) {
            assert_eq!(x.0, y.0);
            assert!(
                (x.1 - y.1).norm() < 1e-12,
                "rot_rot: rz({t1});rz({t2}) ≠ rz({}) at index {}",
                t1 + t2,
                x.0
            );
        }
        let norm: f64 = split.coeffs().iter().map(|(_, c)| c.norm_sqr()).sum();
        assert!(
            (norm - 1.0).abs() < 1e-12,
            "rot_norm_sq: the rotation changed the norm ({norm})"
        );
    }
}

// ===========================================================================
// BranchPhase.lean — the frame identity that discharges `SelfInverse`
// ===========================================================================

/// Drive a `NewWide` through a pseudo-random Clifford+`T` circuit, calling
/// `check` on the live frame every few gates.
fn decomposition_sweep(mut check: impl FnMut(&NewWide, usize, u64)) {
    const N: usize = 12;
    for seed in 0..24u64 {
        let mut tab: NewWide = Driver::new_seeded(N, 1e-12, seed);
        let mut rng = seeded_rng(seed ^ 0xD1CE);
        for step in 0..60 {
            let a = rng.random_range(0..N);
            let b = rng.random_range(0..N);
            match rng.random_range(0..6u32) {
                0 => tab.h(a),
                1 => tab.s(a),
                2 => tab.sqrt_y(a),
                3 if a != b => tab.cnot(a, b),
                4 if a != b => tab.cz(a, b),
                _ => tab.t(a),
            }
            if step % 7 == 0 {
                check(&tab, step, seed);
            }
        }
    }
}

/// `lean/PPVM/Tableau/BranchPhase.lean` (`selfInverse_branchPhase_iff`,
/// `frameOp_involutive_iff`, `frameInvolution_zero_iff`).
///
/// Every case-a/case-b theorem in `lean/PPVM/Tableau/Projection.lean` is stated
/// under `SelfInverse s φ`, with `φ` abstract. `BranchPhase.lean` proves that for
/// the `φ` the crate actually computes (`phase_decomp` +
/// `compute_phase_with_mask_static`) that hypothesis is *equivalent* to the
/// single `ℤ/2` frame identity
///
/// ```text
/// phase_decomp + popcount(destab_anticomm ∧ stab_anticomm)
///              + popcount(stab_anticomm ∧ odd_phase_mask) ≡ 0  (mod 2)
/// ```
///
/// — which `frameOp_involutive_iff` identifies as `M² = I`. This checks the
/// implementation really lands on it, so the oracle is not vacuous here.
#[test]
fn decomposition_satisfies_the_lean_frame_identity() {
    let n = 12usize;
    let mut checked = 0usize;
    let mut any_odd_mask = false;
    decomposition_sweep(|tab, step, seed| {
        let mask = tab.odd_phase_destabilizer_mask();
        any_odd_mask |= mask != 0;
        for q in 0..n {
            for p in [Pauli2::X, Pauli2::Y, Pauli2::Z] {
                let (phase, stab, destab) = tab.compute_decomposition(q, p);
                let parity =
                    (u32::from(phase) + (destab & stab).count_ones() + (stab & mask).count_ones())
                        % 2;
                assert_eq!(
                    parity, 0,
                    "FrameInvolution violated (seed {seed}, step {step}, qubit {q}, {p:?}): \
                     phase {phase}, stab {stab:b}, destab {destab:b}, mask {mask:b}"
                );
                // Case b (`stab_anticomm == 0`) specializes to the crate's
                // `debug_assert!(phase_decomp == 0 || phase_decomp == 2)`
                // ("Measurement result cannot be imaginary!"), which
                // `frameInvolution_zero_iff` turns into a theorem.
                if stab == 0 {
                    assert!(
                        phase == 0 || phase == 2,
                        "case-b phase {phase} is imaginary"
                    );
                }
                checked += 1;
            }
        }
    });
    assert!(
        checked > 5_000,
        "sweep degenerated ({checked} decompositions)"
    );
    // The rows of a valid frame are Hermitian, so their `ℤ/4` phases are even
    // and the odd-phase mask never fires; the mask term is computed because old
    // computes it (behaviour preservation), not because it can contribute. The
    // Lean statements hold for arbitrary masks regardless.
    assert!(
        !any_odd_mask,
        "odd-phase destabilizer mask unexpectedly non-empty"
    );
}

/// `lean/PPVM/Tableau/BranchPhase.lean` (`rot2_order_irrelevant`,
/// `dot_crateWeight_order`) with `lean/PPVM/Tableau/Frame.lean`
/// (`omega_eq_frame_coords`, `omega_disjoint_support`).
///
/// `rotate_2` composes two single-site relabels, `b` before `a`. The whole order
/// dependence of the accumulated `ℤ/4` phase is
/// `⟨destab_a, stab_b⟩ + ⟨destab_b, stab_a⟩`, which is `ω(P_a, P_b)` read in
/// frame coordinates — zero for Paulis on distinct qubits. This checks that on
/// the real frames, so the `b`-before-`a` order is provably phase-neutral (it is
/// kept anyway: what it pins is the float summation order).
#[test]
fn rot2_application_order_is_phase_neutral() {
    let n = 12usize;
    decomposition_sweep(|tab, step, seed| {
        for a in 0..n {
            for b in 0..n {
                if a == b {
                    continue;
                }
                for pa in [Pauli2::X, Pauli2::Y, Pauli2::Z] {
                    for pb in [Pauli2::X, Pauli2::Y, Pauli2::Z] {
                        let (_, stab_a, destab_a) = tab.compute_decomposition(a, pa);
                        let (_, stab_b, destab_b) = tab.compute_decomposition(b, pb);
                        let cross = ((destab_a & stab_b).count_ones()
                            + (destab_b & stab_a).count_ones())
                            % 2;
                        assert_eq!(
                            cross, 0,
                            "ω(P_a, P_b) ≠ 0 on distinct sites \
                             (seed {seed}, step {step}, {a}:{pa:?} vs {b}:{pb:?})"
                        );
                    }
                }
            }
        }
    });
}
