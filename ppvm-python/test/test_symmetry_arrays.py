# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0

"""Tests for the array-form symmetry primitives on the ``(basis_arr, coeffs)``
representation used by ``Lindbladian.pc_step_arr``, as exported from the
``ppvm`` package (thin dtype-coercing wrappers over ``ppvm._core``):

- ``canonicalize_basis_arr`` — plain real merge (sums colliding coefficients)
- ``canonicalize_basis_arr_complex`` — momentum-sector projection (averages
  over the distinct orbit members with the character weight)
- ``check_momentum_sector_arr`` — validation that an input really lies in the
  sector it is about to be projected onto

References are computed here in numpy from the group action, independent of
the Rust merge routines.
"""

import cmath

import numpy as np
import pytest

from ppvm import (
    TranslationGroup,
    canonicalize_basis_arr,
    canonicalize_basis_arr_complex,
    check_momentum_sector_arr,
)

_CODE = {"I": 0, "X": 1, "Z": 2, "Y": 3}
_CHAR = {v: k for k, v in _CODE.items()}


def basis_arr(strings):
    return np.array([[_CODE[c] for c in s] for s in strings], dtype=np.uint8)


def string(row):
    return "".join(_CHAR[int(c)] for c in row)


def to_dict(pair):
    words, coeffs = pair
    return {string(w): c for w, c in zip(words, coeffs, strict=True)}


def rep_of(group, s):
    return string(group.canonicalize(np.array([_CODE[c] for c in s], dtype=np.uint8)))


def momentum(*modes):
    """The wrappers coerce momentum for us; most tests pass a plain tuple.

    See `test_wrappers_coerce_argument_dtypes` for the coercion itself.
    """
    return modes


def z_strings(n):
    return ["I" * j + "Z" + "I" * (n - j - 1) for j in range(n)]


# ── argument coercion (the reason the wrappers exist) ────────────────────────
def test_wrappers_coerce_argument_dtypes():
    """The compiled entry points demand exact dtypes — uint8 basis, float64 /
    complex128 coefficients, int32 momentum. numpy's default integer dtype is
    int64, so an unwrapped ``np.array([0])`` momentum is rejected; the wrappers
    accept plain Python sequences and default-dtype arrays.
    """
    n = 4
    g = TranslationGroup.chain_1d(n)
    words = z_strings(n)
    py_basis = [[_CODE[c] for c in s] for s in words]  # list[list[int]]

    real = to_dict(canonicalize_basis_arr(py_basis, [1.0] * n, g))
    assert real == pytest.approx({rep_of(g, words[0]): float(n)})

    # int64 momentum (numpy default) and a plain list of complex.
    cx = to_dict(canonicalize_basis_arr_complex(py_basis, [1 + 0j] * n, g, np.array([0])))
    assert len(cx) == 1
    assert check_momentum_sector_arr(py_basis, [1 + 0j] * n, g, [0]) is None


# ── canonicalize_basis_arr (real, k=0) ───────────────────────────────────────
def test_canonicalize_basis_arr_sums_collisions():
    n = 4
    g = TranslationGroup.chain_1d(n)
    words = z_strings(n)
    coeffs = np.array([1.0, 2.0, 3.0, 4.0])
    merged = to_dict(canonicalize_basis_arr(basis_arr(words), coeffs, g))
    assert merged == pytest.approx({rep_of(g, words[0]): 10.0})


def test_canonicalize_basis_arr_matches_manual_grouping():
    """Reference: group rows by their rep in numpy and sum."""
    n = 4
    g = TranslationGroup.chain_1d(n)
    rng = np.random.default_rng(3)
    words = [*z_strings(n), "XXII", "IXXI", "IIXX", "XIIX", "ZZZZ"]
    coeffs = rng.normal(size=len(words))

    expected: dict[str, float] = {}
    for w, c in zip(words, coeffs, strict=True):
        expected[rep_of(g, w)] = expected.get(rep_of(g, w), 0.0) + c

    merged = to_dict(canonicalize_basis_arr(basis_arr(words), coeffs, g))
    assert merged.keys() == expected.keys()
    for w in expected:
        assert merged[w] == pytest.approx(expected[w])


def test_canonicalize_basis_arr_validates_shapes():
    g = TranslationGroup.chain_1d(4)
    with pytest.raises(ValueError, match="3 qubits per row but group acts on 4"):
        canonicalize_basis_arr(basis_arr(["ZII"]), np.array([1.0]), g)
    with pytest.raises(ValueError, match="coeffs has length 2 but basis has 1 rows"):
        canonicalize_basis_arr(basis_arr(["ZIII"]), np.array([1.0, 2.0]), g)


# ── canonicalize_basis_arr_complex (momentum sectors) ────────────────────────
def _momentum_seed(n, k):
    """``O_k = Σ_a e^{-2πi k a / n} Z_a`` as ``(basis_arr, coeffs)``."""
    words = z_strings(n)
    coeffs = np.array([cmath.exp(-2j * cmath.pi * k * a / n) for a in range(n)])
    return basis_arr(words), coeffs


@pytest.mark.parametrize("k", [0, 1, 2, 3])
def test_complex_merge_of_momentum_eigenstate_has_unit_rep_coefficient(k):
    """The projection *averages* over the orbit, so a normalized momentum
    eigenstate folds to a rep coefficient of modulus 1."""
    n = 4
    g = TranslationGroup.chain_1d(n)
    words, coeffs = _momentum_seed(n, k)
    merged = to_dict(canonicalize_basis_arr_complex(words, coeffs, g, momentum(k)))
    assert len(merged) == 1
    assert abs(next(iter(merged.values()))) == pytest.approx(1.0)


def test_complex_merge_projects_out_other_sectors():
    """A pure k=1 state has zero component in every other sector."""
    n = 4
    g = TranslationGroup.chain_1d(n)
    words, coeffs = _momentum_seed(n, 1)
    for k_other in [0, 2, 3]:
        merged = to_dict(canonicalize_basis_arr_complex(words, coeffs, g, momentum(k_other)))
        for c in merged.values():
            assert abs(c) < 1e-12, f"k=1 state leaked into sector {k_other}: {c}"


def test_complex_merge_matches_character_average():
    """Reference: (1/|orbit|) Σ_{p in orbit} χ_k(g_p) · c_p, computed here
    by walking the cyclic shifts explicitly."""
    n = 4
    k = 1
    g = TranslationGroup.chain_1d(n)
    rng = np.random.default_rng(11)
    words = z_strings(n)
    coeffs = rng.normal(size=n) + 1j * rng.normal(size=n)

    # Z_a is the shift of Z_0 by `a`, so the character weight is e^{2πika/n}.
    by_word = dict(zip(words, coeffs, strict=True))
    rep = rep_of(g, words[0])
    shift_of_rep = words.index(rep)
    expected = (
        sum(
            cmath.exp(2j * cmath.pi * k * ((a - shift_of_rep) % n) / n) * by_word[words[a]]
            for a in range(n)
        )
        / n
    )

    merged = to_dict(canonicalize_basis_arr_complex(basis_arr(words), coeffs, g, momentum(k)))
    assert merged[rep] == pytest.approx(expected)


def test_complex_merge_validates_shapes():
    g = TranslationGroup.chain_1d(4)
    words, coeffs = _momentum_seed(4, 1)
    with pytest.raises(ValueError, match="momentum has 2 entries but group has 1 generators"):
        canonicalize_basis_arr_complex(words, coeffs, g, momentum(0, 0))
    with pytest.raises(ValueError, match="coeffs has length 2 but basis has 4 rows"):
        canonicalize_basis_arr_complex(words, coeffs[:2], g, momentum(1))
    with pytest.raises(ValueError, match="3 qubits per row but group acts on 4"):
        canonicalize_basis_arr_complex(basis_arr(["ZII"]), np.array([1 + 0j]), g, momentum(0))


# ── check_momentum_sector_arr ────────────────────────────────────────────────
@pytest.mark.parametrize("k", [0, 1, 2, 3])
def test_check_momentum_sector_accepts_eigenstate(k):
    n = 4
    g = TranslationGroup.chain_1d(n)
    words, coeffs = _momentum_seed(n, k)
    assert check_momentum_sector_arr(words, coeffs, g, momentum(k)) is None


def test_check_momentum_sector_rejects_wrong_sector():
    n = 4
    g = TranslationGroup.chain_1d(n)
    words, coeffs = _momentum_seed(n, 1)
    with pytest.raises(ValueError, match="not in target momentum sector"):
        check_momentum_sector_arr(words, coeffs, g, momentum(0))


def test_check_momentum_sector_rejects_incomplete_orbit():
    """Orbit members missing from the basis count as zero, so a lone Z_0 is
    not a momentum eigenstate."""
    g = TranslationGroup.chain_1d(4)
    with pytest.raises(ValueError, match="not in target momentum sector"):
        check_momentum_sector_arr(basis_arr(["ZIII"]), np.array([1 + 0j]), g, momentum(0))


def test_check_momentum_sector_flags_incompatible_stabilizer():
    """``ZIZI`` has a period-2 stabilizer, which cannot carry k=1."""
    g = TranslationGroup.chain_1d(4)
    with pytest.raises(ValueError, match="stabilizer incompatible with momentum sector"):
        check_momentum_sector_arr(
            basis_arr(["ZIZI", "IZIZ"]),
            np.array([1 + 0j, -1 + 0j]),
            g,
            momentum(1),
        )


def test_check_momentum_sector_tolerance_is_configurable():
    n = 4
    g = TranslationGroup.chain_1d(n)
    words, coeffs = _momentum_seed(n, 1)
    perturbed = coeffs.copy()
    perturbed[0] += 1e-6
    with pytest.raises(ValueError, match="not in target momentum sector"):
        check_momentum_sector_arr(words, perturbed, g, momentum(1), 1e-9)
    # Same input passes once the tolerance exceeds the perturbation.
    assert check_momentum_sector_arr(words, perturbed, g, momentum(1), 1e-4) is None


def test_check_momentum_sector_rejects_invalid_tolerance():
    n = 4
    g = TranslationGroup.chain_1d(n)
    words, coeffs = _momentum_seed(n, 0)
    with pytest.raises(ValueError, match="invalid tolerance"):
        check_momentum_sector_arr(words, coeffs, g, momentum(0), -1.0)
