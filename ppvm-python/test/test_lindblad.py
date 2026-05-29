# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
"""Cross-checks for the adjoint Pauli-Lindbladian shim (``ppvm.Lindbladian``).

Each test compares the Rust shim's ``action`` / ``leakage`` / ``generator``
against an independent, dense reference built from the single-qubit Pauli
multiplication table, on a few representative Lindbladians (XY + Z dephasing,
TFIM + X dephasing).
"""

from __future__ import annotations

import numpy as np
import scipy.sparse.linalg as spla

from ppvm import Lindbladian, sigma_minus, sigma_plus

# --- independent reference via the dense Pauli multiplication table ---------
I, X, Y, Z = range(4)  # noqa: E741  (I is standard Pauli notation here)
MUL = {
    (I, I): (1, I),
    (I, X): (1, X),
    (I, Y): (1, Y),
    (I, Z): (1, Z),
    (X, I): (1, X),
    (X, X): (1, I),
    (X, Y): (1j, Z),
    (X, Z): (-1j, Y),
    (Y, I): (1, Y),
    (Y, X): (-1j, Z),
    (Y, Y): (1, I),
    (Y, Z): (1j, X),
    (Z, I): (1, Z),
    (Z, X): (1j, Y),
    (Z, Y): (-1j, X),
    (Z, Z): (1, I),
}
CODE = {I: "I", X: "X", Y: "Y", Z: "Z"}
ICODE = {"I": I, "X": X, "Y": Y, "Z": Z}


def _str_to_tuple(s):
    return tuple(ICODE[ch] for ch in s)


def _tuple_to_str(t):
    return "".join(CODE[ch] for ch in t)


def _mul_pauli(p, q):
    """Return (phase, p·q) for two Pauli strings given as code tuples."""
    phase = 1
    r = list(p)
    for i, (pi, qi) in enumerate(zip(p, q)):
        ph, rr = MUL[(pi, qi)]
        phase *= ph
        r[i] = rr
    return phase, tuple(r)


def _reference_action(p_str, h_terms, jump_terms):
    """L*(p) = i[H, p] + sum_k gamma_k (L_k p L_k - p), term by term."""
    p = _str_to_tuple(p_str)
    out: dict = {}
    for h_str, coeff_h in h_terms:
        h = _str_to_tuple(h_str)
        ph_pp, pp = _mul_pauli(h, p)
        ph_pq, _ = _mul_pauli(p, h)
        # i[H, p] = i (Hp - pH) = i (ph_pp - ph_pq) r ; real for Hermitian H, p.
        coeff = (1j * coeff_h * (ph_pp - ph_pq)).real
        if coeff:
            out[pp] = out.get(pp, 0.0) + coeff
    for j_str, gamma in jump_terms:
        j = _str_to_tuple(j_str)
        ph_pp, _ = _mul_pauli(j, p)
        # Hermitian Pauli L, p: Lp has imaginary phase iff they anti-commute,
        # and then L p L = -p, contributing -2 gamma p.
        if abs(ph_pp.imag) > 0.5:
            out[p] = out.get(p, 0.0) + (-2.0 * gamma)
    return {kk: v for kk, v in out.items() if v}


def _to_str_dict(d):
    return {_tuple_to_str(kk): v for kk, v in d.items()}


def _xy_dephasing(L, alpha, gamma):
    pairs = [
        (a, b, 1.0 / min(b - a, L - b + a) ** alpha) for a in range(L) for b in range(a + 1, L)
    ]
    kac = sum(j for _, _, j in pairs) / L
    pairs = [(a, b, j / kac) for a, b, j in pairs]
    h_terms = []
    for a, b, j in pairs:
        for q in "XY":
            term = ["I"] * L
            term[a] = term[b] = q
            h_terms.append(("".join(term), j))
    jump_terms = []
    for i in range(L):
        term = ["I"] * L
        term[i] = "Z"
        jump_terms.append(("".join(term), gamma))
    return h_terms, jump_terms


def _tfim_xdeph(L, J, h, gamma):
    h_terms = []
    for i in range(L - 1):
        term = ["I"] * L
        term[i] = term[i + 1] = "Z"
        h_terms.append(("".join(term), J))
    for i in range(L):
        term = ["I"] * L
        term[i] = "X"
        h_terms.append(("".join(term), h))
    jump_terms = []
    for i in range(L):
        term = ["I"] * L
        term[i] = "X"
        jump_terms.append(("".join(term), gamma))
    return h_terms, jump_terms


def _random_pauli_str(rng, L):
    chars = ["I"] * L
    positions = rng.choice(L, size=int(rng.integers(1, L + 1)), replace=False)
    for q in positions:
        chars[q] = rng.choice(["X", "Y", "Z"])
    return "".join(chars)


def _assert_action_matches(L_op, h_terms, jump_terms, strings):
    for p in strings:
        got = L_op.action(p)
        want = _to_str_dict(_reference_action(p, h_terms, jump_terms))
        for kk in set(got) | set(want):
            assert abs(got.get(kk, 0.0) - want.get(kk, 0.0)) < 1e-12, (
                f"action mismatch at p={p!r} k={kk!r}: "
                f"shim={got.get(kk, 0.0)} ref={want.get(kk, 0.0)}"
            )


def test_action_xy_dephasing():
    L = 8
    h_terms, jump_terms = _xy_dephasing(L, alpha=1.0, gamma=0.3)
    L_op = Lindbladian(L, h_terms, jump_terms)
    rng = np.random.default_rng(42)
    strings = ["I" * L]
    strings += ["I" * i + "Z" + "I" * (L - i - 1) for i in range(L)]
    strings += [_random_pauli_str(rng, L) for _ in range(20)]
    _assert_action_matches(L_op, h_terms, jump_terms, strings)


def test_action_tfim_xdephasing():
    L = 6
    h_terms, jump_terms = _tfim_xdeph(L, J=0.7, h=0.4, gamma=0.2)
    L_op = Lindbladian(L, h_terms, jump_terms)
    rng = np.random.default_rng(7)
    strings = [_random_pauli_str(rng, L) for _ in range(30)]
    _assert_action_matches(L_op, h_terms, jump_terms, strings)


def test_generator_leakage_and_expm():
    L = 5
    dt = 0.1
    h_terms, jump_terms = _xy_dephasing(L, alpha=1.0, gamma=0.4)
    L_op = Lindbladian(L, h_terms, jump_terms)
    basis = ["I" * i + "Z" + "I" * (L - i - 1) for i in range(L)]
    basis += ["YIYII", "IYIYI", "ZZIII"]
    coeffs = np.array([0.5, -0.3, 0.2, 0.4, -0.1, 0.6, 0.2, 0.1])

    M_shim = L_op.generator(basis).toarray()
    M_ref = np.zeros((len(basis), len(basis)))
    index = {p: i for i, p in enumerate(basis)}
    for col, p in enumerate(basis):
        for r_tuple, v in _reference_action(p, h_terms, jump_terms).items():
            r = _tuple_to_str(r_tuple)
            if r in index:
                M_ref[index[r], col] += v
    assert np.max(np.abs(M_shim - M_ref)) < 1e-12

    shim_leak = L_op.leakage(basis, coeffs)
    ref_leak: dict = {}
    for p, cf in zip(basis, coeffs):
        for r_tuple, v in _reference_action(p, h_terms, jump_terms).items():
            r = _tuple_to_str(r_tuple)
            if r not in index:
                ref_leak[r] = ref_leak.get(r, 0.0) + v * cf
    ref_leak = {kk: v for kk, v in ref_leak.items() if v}
    for kk in set(shim_leak) | set(ref_leak):
        assert abs(shim_leak.get(kk, 0.0) - ref_leak.get(kk, 0.0)) < 1e-12

    c_shim = spla.expm_multiply(dt * L_op.generator(basis), coeffs)
    c_ref = spla.expm_multiply(dt * M_ref, coeffs)
    assert np.allclose(c_shim, c_ref, atol=1e-13)


def test_protected_strings_suppressed():
    L = 4
    h_terms, jump_terms = _xy_dephasing(L, alpha=1.0, gamma=0.0)
    L_op = Lindbladian(L, h_terms, jump_terms)
    basis = ["ZIII"]
    coeffs = np.array([1.0])
    leak = L_op.leakage(basis, coeffs)
    assert leak, "expected some leakage"
    protected_key = next(iter(leak))
    leak2 = L_op.leakage(basis, coeffs, protected=[protected_key])
    assert protected_key not in leak2


# --- non-Hermitian dissipators: dense Liouvillian reference -----------------
#
# The dense reference builds full 2^L × 2^L matrices for H and the jumps and
# evaluates L*(p) = i[H, p] + Σ_k γ_k (L_k† p L_k − ½ {L_k† L_k, p}) directly,
# then expands the result in the Pauli basis. Restricted to L ≤ 3 it costs
# 4^L · 8^L FLOPs (a few ms) but it makes no assumption about the shape of
# L_k, so it doubles as a correctness check for the general dissipator path
# in the Rust shim.

_DENSE_PAULI = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}


def _pauli_mat(s):
    """Dense matrix for Pauli string ``s`` (leftmost char = leftmost factor)."""
    M = np.array([[1.0]], dtype=complex)
    for c in s:
        M = np.kron(M, _DENSE_PAULI[c])
    return M


def _all_strings(L):
    if L == 0:
        return [""]
    sub = _all_strings(L - 1)
    return [c + s for c in "IXYZ" for s in sub]


def _dense_action(H, jumps, p_str, L):
    """L*(p) computed densely, returned as a real Pauli-coefficient dict."""
    p_mat = _pauli_mat(p_str)
    # np.errstate quiets spurious "divide by zero" warnings emitted by
    # macOS Accelerate when complex @ contains exact zeros.
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        out_mat = 1j * (H @ p_mat - p_mat @ H)
        for Lop, gamma in jumps:
            Ld = Lop.conj().T
            out_mat += gamma * (
                Ld @ p_mat @ Lop - 0.5 * (Ld @ Lop @ p_mat + p_mat @ Ld @ Lop)
            )
    d = 2**L
    out = {}
    for q_str in _all_strings(L):
        coef = np.trace(_pauli_mat(q_str) @ out_mat) / d
        assert abs(coef.imag) < 1e-9, f"non-real coefficient for {q_str}: {coef}"
        if abs(coef.real) > 1e-12:
            out[q_str] = coef.real
    return out


def _embed_op(op_1q, site, L):
    """Embed a single-qubit dense operator at ``site`` of an L-qubit register."""
    eye = _DENSE_PAULI["I"]
    M = np.array([[1.0]], dtype=complex)
    for j in range(L):
        M = np.kron(M, op_1q if j == site else eye)
    return M


_SIGMA_MINUS = np.array([[0, 0], [1, 0]], dtype=complex)
_SIGMA_PLUS = np.array([[0, 1], [0, 0]], dtype=complex)


def test_amplitude_damping_action():
    L = 3
    gamma = 0.5
    h_terms = [("XXI", 1.0), ("IXX", 0.7), ("ZIZ", 0.3)]
    jump_terms = [(sigma_minus(i, L), gamma) for i in range(L)]
    L_op = Lindbladian(L, h_terms, jump_terms)

    H = sum(c * _pauli_mat(s) for s, c in h_terms)
    jumps_dense = [(_embed_op(_SIGMA_MINUS, i, L), gamma) for i in range(L)]

    rng = np.random.default_rng(11)
    strings = ["III", "ZII", "IZI", "IIZ", "XYZ", "YYZ"]
    strings += [_random_pauli_str(rng, L) for _ in range(10)]

    for p in strings:
        got = L_op.action(p)
        want = _dense_action(H, jumps_dense, p, L)
        for k in set(got) | set(want):
            diff = abs(got.get(k, 0.0) - want.get(k, 0.0))
            assert diff < 1e-10, (
                f"sigma_minus action mismatch at p={p!r} k={k!r}: "
                f"shim={got.get(k, 0.0)} ref={want.get(k, 0.0)} diff={diff}"
            )


def test_thermal_excitation_damping_action():
    """σ⁺ + σ⁻ jumps together (thermal bath at finite temperature)."""
    L = 2
    h_terms = [("XX", 1.0), ("ZI", 0.2), ("IZ", 0.1)]
    jump_terms = [(sigma_minus(i, L), 0.4) for i in range(L)] + [
        (sigma_plus(i, L), 0.1) for i in range(L)
    ]
    L_op = Lindbladian(L, h_terms, jump_terms)

    H = sum(c * _pauli_mat(s) for s, c in h_terms)
    jumps_dense = [(_embed_op(_SIGMA_MINUS, i, L), 0.4) for i in range(L)] + [
        (_embed_op(_SIGMA_PLUS, i, L), 0.1) for i in range(L)
    ]

    for p in _all_strings(L):
        got = L_op.action(p)
        want = _dense_action(H, jumps_dense, p, L)
        for k in set(got) | set(want):
            diff = abs(got.get(k, 0.0) - want.get(k, 0.0))
            assert diff < 1e-10, (
                f"sigma_plus/sigma_minus action mismatch at p={p!r} k={k!r}: "
                f"shim={got.get(k, 0.0)} ref={want.get(k, 0.0)}"
            )


def test_amplitude_damping_generator_and_leakage():
    L = 3
    gamma = 0.3
    h_terms = [("XXI", 0.5), ("IXX", 0.5), ("ZII", 0.2), ("IZI", 0.2), ("IIZ", 0.2)]
    jump_terms = [(sigma_minus(0, L), gamma), (sigma_minus(2, L), gamma)]
    L_op = Lindbladian(L, h_terms, jump_terms)

    basis = ["III", "ZII", "IZI", "IIZ", "ZZI", "IZZ"]
    coeffs = np.array([0.1, 0.5, -0.3, 0.4, 0.2, -0.1])

    H = sum(c * _pauli_mat(s) for s, c in h_terms)
    jumps_dense = [(_embed_op(_SIGMA_MINUS, i, L), gamma) for i in (0, 2)]

    M_shim = L_op.generator(basis).toarray()
    M_ref = np.zeros((len(basis), len(basis)))
    idx = {p: i for i, p in enumerate(basis)}
    leak_ref: dict = {}
    for col, p in enumerate(basis):
        action_p = _dense_action(H, jumps_dense, p, L)
        for q, v in action_p.items():
            if q in idx:
                M_ref[idx[q], col] += v
            else:
                leak_ref[q] = leak_ref.get(q, 0.0) + v * coeffs[col]
    assert np.max(np.abs(M_shim - M_ref)) < 1e-10

    leak_shim = L_op.leakage(basis, coeffs)
    leak_ref = {k: v for k, v in leak_ref.items() if abs(v) > 1e-14}
    for k in set(leak_shim) | set(leak_ref):
        diff = abs(leak_shim.get(k, 0.0) - leak_ref.get(k, 0.0))
        assert diff < 1e-10, f"leakage mismatch at k={k!r}: diff={diff}"


# --- convergence to closed bilinear evolution (NN XY + Z dephasing, OBC) ----
#
# For nearest-neighbour interactions the Jordan-Wigner fermion bilinears
# B_{mn} = c_m† c_n stay closed under the adjoint Lindbladian, so the spin
# correlator C_j(t) = 2^{-L} Tr[Z_j(t) Z_i] obeys a tractable L²×L² ODE.
# We use open boundary conditions: with OBC the JW transformation maps the
# spin XY chain to a bilinear fermion hopping problem exactly, free of the
# parity-twist that makes PBC bilinears match only up to 1/L finite-size
# corrections.


def _bilinear_nn_xy_z_dephasing_obc(L, J, gamma, times, site0):
    """Closed bilinear evolution of C_j(t) for the NN XY + Z-dephasing model
    with open boundary conditions.

    `F_{mn}(t) = 2^{-L} Tr[B_{mn}(t) Z_i]` evolves as
    `∂_t F_{mn} = i·2J·(F_{m+1,n}+F_{m-1,n}-F_{m,n+1}-F_{m,n-1}) - 4γ(1-δ_{mn})F_{mn}`
    on the L×L lattice; OBC means edge terms (m=0 or m=L-1, etc.) drop.
    Z_j = I - 2 n_j gives C_j = -2 F_{jj}.
    """
    dim = L * L

    def idx(m, n):
        return m * L + n

    gen = np.zeros((dim, dim), dtype=complex)
    for m in range(L):
        for n in range(L):
            row = idx(m, n)
            if m + 1 < L:
                gen[row, idx(m + 1, n)] += 1j * 2 * J
            if m - 1 >= 0:
                gen[row, idx(m - 1, n)] += 1j * 2 * J
            if n + 1 < L:
                gen[row, idx(m, n + 1)] += -1j * 2 * J
            if n - 1 >= 0:
                gen[row, idx(m, n - 1)] += -1j * 2 * J
            if m != n:
                gen[row, row] += -4 * gamma

    evals, evecs = np.linalg.eig(gen)
    evecs_inv = np.linalg.inv(evecs)
    f0 = np.zeros(dim, dtype=complex)
    f0[idx(site0, site0)] = -0.5
    coeffs = evecs_inv @ f0

    corr = np.empty((len(times), L))
    diag = [idx(j, j) for j in range(L)]
    for nt, t in enumerate(times):
        ft = evecs @ (np.exp(evals * t) * coeffs)
        corr[nt] = np.real(-2 * ft[diag])
    return corr


def _adaptive_z_correlator(L_op, L, site0, dt, n_steps, tau_add):
    """Adaptive Heisenberg-picture evolution of Z_{site0} on a growing Pauli
    basis. Returns the matrix `C_j(t_n) = ⟨Z_j, Z_{site0}(t_n)⟩` for every
    site `j` and step `n`. Matrix exponential is exact in `dt` within the
    basis; the only approximation is dropping leakage strings of magnitude
    below `tau_add`."""
    z_strings = ["I" * j + "Z" + "I" * (L - j - 1) for j in range(L)]
    basis = [z_strings[site0]]
    coeffs = np.array([1.0])
    protected = [z_strings[site0]]

    corr = np.zeros((n_steps + 1, L))
    corr[0, site0] = 1.0

    for step in range(n_steps):
        leak = L_op.leakage(basis, coeffs, protected=protected)
        new = [k for k, v in leak.items() if abs(v) > tau_add]
        if new:
            basis = basis + new
            coeffs = np.concatenate([coeffs, np.zeros(len(new))])
        M = L_op.generator(basis)
        coeffs = spla.expm_multiply(dt * M, coeffs)
        index = {s: i for i, s in enumerate(basis)}
        for j in range(L):
            if z_strings[j] in index:
                corr[step + 1, j] = coeffs[index[z_strings[j]]]
    return corr


def _adaptive_z_correlator_pc(L_op, L, site0, dt, n_steps, tau_add):
    """Same as :func:`_adaptive_z_correlator` but with a predictor-corrector
    basis expansion at each step.

    1. Expand basis with first-hop leakage from the current state and save
       the (pre-step) coefficients.
    2. Run a predictor step via expm_multiply.
    3. Compute leakage from the predicted state — these are the second-hop
       strings the predictor flowed into but did not have in basis.
    4. Add the second-hop strings, pad the saved pre-step coefficients with
       zeros for them, and redo the step on the enlarged basis from the
       pre-step state.

    The corrector basis covers the second-hop flow, lifting the local
    truncation per step from O(dt²) to O(dt³).
    """
    z_strings = ["I" * j + "Z" + "I" * (L - j - 1) for j in range(L)]
    basis = [z_strings[site0]]
    coeffs = np.array([1.0])
    protected = [z_strings[site0]]

    corr = np.zeros((n_steps + 1, L))
    corr[0, site0] = 1.0

    for step in range(n_steps):
        # Predictor basis: add first-hop leakage from current state.
        leak = L_op.leakage(basis, coeffs, protected=protected)
        new = [k for k, v in leak.items() if abs(v) > tau_add]
        if new:
            basis = basis + new
            coeffs = np.concatenate([coeffs, np.zeros(len(new))])
        coeffs_pre = coeffs.copy()
        # Predict.
        M = L_op.generator(basis)
        coeffs_predict = spla.expm_multiply(dt * M, coeffs)
        # Corrector basis: second-hop strings from the predicted state.
        leak2 = L_op.leakage(basis, coeffs_predict, protected=protected)
        new2 = [k for k, v in leak2.items() if abs(v) > tau_add]
        if new2:
            basis = basis + new2
            coeffs_pre = np.concatenate([coeffs_pre, np.zeros(len(new2))])
        # Redo the step on the enlarged basis from the pre-step state.
        M = L_op.generator(basis)
        coeffs = spla.expm_multiply(dt * M, coeffs_pre)
        index = {s: i for i, s in enumerate(basis)}
        for j in range(L):
            if z_strings[j] in index:
                corr[step + 1, j] = coeffs[index[z_strings[j]]]
    return corr


def test_adaptive_converges_to_nn_xy_z_dephasing_bilinear():
    """Halving dt drives the adaptive shim toward the closed bilinear solution
    for the NN XY chain with Z dephasing (OBC).

    The single-hop adaptive scheme has local truncation O(dt²) per step
    (first-step leakage only adds first-hop strings; second-hop strings appear
    inside the matrix exponential but aren't yet in the basis), so global
    error is O(T·dt). The test exercises this dt-scaling on a fixed time
    window with a leakage threshold τ tight enough that truncation is well
    below the integrator floor.
    """
    L = 4
    J = 1.0
    gamma = 1.0
    site0 = L // 2
    T = 0.05
    tau_add = 1e-12  # tight enough that integrator (T·dt) dominates

    h_terms = []
    for i in range(L - 1):  # OBC: L-1 bonds, no wraparound
        a, b = i, i + 1
        xs = ["I"] * L
        xs[a] = xs[b] = "X"
        ys = ["I"] * L
        ys[a] = ys[b] = "Y"
        h_terms += [("".join(xs), J), ("".join(ys), J)]
    jump_terms = [("I" * j + "Z" + "I" * (L - j - 1), gamma) for j in range(L)]
    L_op = Lindbladian(L, h_terms, jump_terms)

    errors = []
    final_corr = None
    for dt in (0.01, 0.005, 0.0025):
        n_steps = round(T / dt)
        times = np.arange(n_steps + 1) * dt
        exact = _bilinear_nn_xy_z_dephasing_obc(L, J, gamma, times, site0)
        L_op.clear_cache()
        shim = _adaptive_z_correlator(L_op, L, site0, dt, n_steps, tau_add)
        # Compare only the endpoint to keep the comparison decoupled from
        # exactly which step counts we use.
        errors.append(float(np.max(np.abs(shim[-1] - exact[-1]))))
        if dt == 0.0025:
            final_corr = (shim[-1], exact[-1])

    # First-order dt-scaling: halving dt should roughly halve the error.
    # Allow a factor of 2 slack on the predicted 2× drop.
    assert errors[1] < 0.8 * errors[0], f"dt-halving did not help: {errors}"
    assert errors[2] < 0.8 * errors[1], f"dt-halving did not help: {errors}"

    # At dt = 0.0025 the integrator floor is ~T·dt = 1.25e-4. The shim
    # should sit comfortably under 1e-3 against the bilinear reference.
    assert errors[-1] < 1e-3, (
        f"shim vs bilinear at smallest dt: max abs error = {errors[-1]:.3g}; "
        f"shim={final_corr[0]}, exact={final_corr[1]}"
    )


def test_predictor_corrector_lifts_dt_scaling_to_cubic():
    """The predictor-corrector basis expansion lifts the single-hop scheme's
    local O(dt²) truncation to O(dt³), so the PC error at fixed T scales as
    dt² (one order steeper than the dt-scaling of single-hop). The PC error
    is also strictly smaller at every dt we test.
    """
    L = 4
    J = 1.0
    gamma = 1.0
    site0 = L // 2
    T = 0.05
    tau_add = 1e-12

    h_terms = []
    for i in range(L - 1):
        a, b = i, i + 1
        xs = ["I"] * L
        xs[a] = xs[b] = "X"
        ys = ["I"] * L
        ys[a] = ys[b] = "Y"
        h_terms += [("".join(xs), J), ("".join(ys), J)]
    jump_terms = [("I" * j + "Z" + "I" * (L - j - 1), gamma) for j in range(L)]
    L_op = Lindbladian(L, h_terms, jump_terms)

    err_single = []
    err_pc = []
    for dt in (0.01, 0.005, 0.0025):
        n_steps = round(T / dt)
        times = np.arange(n_steps + 1) * dt
        exact = _bilinear_nn_xy_z_dephasing_obc(L, J, gamma, times, site0)
        L_op.clear_cache()
        single = _adaptive_z_correlator(L_op, L, site0, dt, n_steps, tau_add)
        L_op.clear_cache()
        pc = _adaptive_z_correlator_pc(L_op, L, site0, dt, n_steps, tau_add)
        err_single.append(float(np.max(np.abs(single[-1] - exact[-1]))))
        err_pc.append(float(np.max(np.abs(pc[-1] - exact[-1]))))

    # PC is strictly more accurate than single-hop at every dt — by ~100×
    # in the test regime. Leave the threshold loose enough that small
    # changes in expm_multiply tolerance don't flake.
    for s, p, dt in zip(err_single, err_pc, (0.01, 0.005, 0.0025)):
        assert p < s / 50, (
            f"PC ({p:.3e}) not meaningfully better than single-hop ({s:.3e}) at dt={dt}"
        )

    # PC's dt-scaling is one order steeper. Halving dt should drop the error
    # by ~8× (dt³ vs single-hop's ~4× = dt²). Require >5× per halving with
    # some safety margin.
    from itertools import pairwise

    for prev, curr in pairwise(err_pc):
        assert curr < prev / 5, f"PC dt-halving ratio < 5: errors {err_pc}"

    # And the smallest-dt PC error itself should sit at FP noise of the
    # bilinear reference.
    assert err_pc[-1] < 1e-7, f"PC at smallest dt: error = {err_pc[-1]:.3e}"


def _adaptive_z_correlator_pc_rust(L_op, L, site0, dt, n_steps, tau_add):
    """Same as :func:`_adaptive_z_correlator_pc` but the per-step PC work
    (leakage expansion, predictor expm, second-hop expansion, corrector
    expm) all runs in Rust through :meth:`Lindbladian.pc_step`."""
    z_strings = ["I" * j + "Z" + "I" * (L - j - 1) for j in range(L)]
    basis = [z_strings[site0]]
    coeffs = np.array([1.0])
    protected = [z_strings[site0]]

    corr = np.zeros((n_steps + 1, L))
    corr[0, site0] = 1.0

    for step in range(n_steps):
        basis, coeffs = L_op.pc_step(basis, coeffs, dt, tau_add, protected=protected)
        index = {s: i for i, s in enumerate(basis)}
        for j in range(L):
            if z_strings[j] in index:
                corr[step + 1, j] = coeffs[index[z_strings[j]]]
    return corr


def test_pc_step_rust_matches_scipy_pc():
    """The pure-Rust PC step (Al-Mohy & Higham expm + parallel SpMV) agrees
    with the scipy-based PC implementation at FP precision. This pins the
    Rust matrix exponential against the well-tested scipy reference under
    the exact same basis-expansion schedule."""
    L = 4
    J = 1.0
    gamma = 1.0
    site0 = L // 2
    dt = 0.01
    n_steps = 5
    tau_add = 1e-12

    h_terms = []
    for i in range(L - 1):
        a, b = i, i + 1
        xs = ["I"] * L
        xs[a] = xs[b] = "X"
        ys = ["I"] * L
        ys[a] = ys[b] = "Y"
        h_terms += [("".join(xs), J), ("".join(ys), J)]
    jump_terms = [("I" * j + "Z" + "I" * (L - j - 1), gamma) for j in range(L)]
    L_op = Lindbladian(L, h_terms, jump_terms)

    L_op.clear_cache()
    rust = _adaptive_z_correlator_pc_rust(L_op, L, site0, dt, n_steps, tau_add)
    L_op.clear_cache()
    scipy_pc = _adaptive_z_correlator_pc(L_op, L, site0, dt, n_steps, tau_add)

    # Both schemes share the same algorithm and tolerances; agreement should
    # be near FP precision of the matrix-exponential approximation.
    diff = float(np.max(np.abs(rust - scipy_pc)))
    assert diff < 1e-10, f"Rust PC differs from scipy PC by {diff:.3e}"


def test_pc_step_rust_dt_scaling_is_cubic():
    """End-to-end: the Rust-only PC step matches the bilinear reference with
    the same dt³ scaling as the scipy version, confirming that the Rust
    matrix exponential is not the bottleneck of accuracy."""
    L = 4
    J = 1.0
    gamma = 1.0
    site0 = L // 2
    T = 0.05
    tau_add = 1e-12

    h_terms = []
    for i in range(L - 1):
        a, b = i, i + 1
        xs = ["I"] * L
        xs[a] = xs[b] = "X"
        ys = ["I"] * L
        ys[a] = ys[b] = "Y"
        h_terms += [("".join(xs), J), ("".join(ys), J)]
    jump_terms = [("I" * j + "Z" + "I" * (L - j - 1), gamma) for j in range(L)]
    L_op = Lindbladian(L, h_terms, jump_terms)

    err = []
    for dt in (0.01, 0.005, 0.0025):
        n_steps = round(T / dt)
        times = np.arange(n_steps + 1) * dt
        exact = _bilinear_nn_xy_z_dephasing_obc(L, J, gamma, times, site0)
        L_op.clear_cache()
        rust = _adaptive_z_correlator_pc_rust(L_op, L, site0, dt, n_steps, tau_add)
        err.append(float(np.max(np.abs(rust[-1] - exact[-1]))))

    # Halving dt should drop the error by ≥5× (cubic gives 8×).
    from itertools import pairwise

    for prev, curr in pairwise(err):
        assert curr < prev / 5, f"Rust PC dt-halving ratio < 5: errors {err}"
    assert err[-1] < 1e-7, f"Rust PC tight-dt error too large: {err[-1]:.3e}"


def test_lincomb_single_term_matches_hermitian_fast_path():
    """A length-1 real lincomb should be routed to the Hermitian fast path
    and produce numerically identical results to passing the string directly."""
    L = 4
    h_terms, jump_simple = _xy_dephasing(L, alpha=1.0, gamma=0.3)
    L_simple = Lindbladian(L, h_terms, jump_simple)
    # Same operator, expressed as a length-1 complex lincomb.
    jump_lincomb = [([(s, 1.0 + 0.0j)], g) for s, g in jump_simple]
    L_lincomb = Lindbladian(L, h_terms, jump_lincomb)

    rng = np.random.default_rng(99)
    for _ in range(5):
        p = _random_pauli_str(rng, L)
        got_a = L_simple.action(p)
        got_b = L_lincomb.action(p)
        for k in set(got_a) | set(got_b):
            assert got_a.get(k, 0.0) == got_b.get(k, 0.0), (
                f"lincomb fast path mismatch at p={p!r}: {got_a} vs {got_b}"
            )
