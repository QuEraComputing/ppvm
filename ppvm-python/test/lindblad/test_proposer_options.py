# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
"""Independent candidate-strategy and retention options across Python entry points."""

import numpy as np
import pytest

from ppvm import Lindbladian


def model():
    return Lindbladian(4, [("XXII", 0.7), ("IYYI", 0.8), ("IIXX", 0.9), ("ZIII", 0.3)], [])


def evolve(op, dt=0.025, **options):
    basis, coeffs = ["ZIII"], np.array([1.0])
    config = dict(max_basis=4, admit_basis=12, drop_tol=0.0, protected=basis, **options)
    for _ in range(100):
        basis, coeffs = op.pc_step(basis, coeffs, dt, **config)
        assert "ZIII" in basis and len(basis) <= 4
        assert np.isfinite(coeffs).all()
    return basis, coeffs


@pytest.mark.parametrize("proposer", ["streaming", "shortlist"])
def test_retention_options(proposer):
    op = model()
    ordinary = evolve(op, proposer=proposer)
    forward = evolve(op, proposer=proposer, forward_looking=True)
    assert ordinary[0] != forward[0] or not np.array_equal(ordinary[1], forward[1])
    a = evolve(op, dt=0.1, proposer=proposer)
    b = evolve(op, dt=0.1, proposer=proposer, forward_looking=True)
    assert a[0] == b[0]
    np.testing.assert_array_equal(a[1], b[1])
    if proposer == "streaming":
        implicit = evolve(op)
        assert implicit[0] == ordinary[0]
        np.testing.assert_array_equal(implicit[1], ordinary[1])


@pytest.mark.parametrize("proposer", ["streaming", "shortlist"])
@pytest.mark.parametrize("forward_looking", [False, True])
def test_array_and_timed_options_agree(proposer, forward_looking):
    op = model()
    basis = np.array([[2, 0, 0, 0]], dtype=np.uint8)
    coeffs = np.array([1.0])
    options = dict(drop_tol=0.0, admit_basis=12, proposer=proposer, forward_looking=forward_looking)
    a = op.pc_step_arr(basis, coeffs, 0.025, 4, protected_arr=basis, **options)
    b, timings = op._spec.pc_step_timed(basis, coeffs, 0.025, 4, protected=basis, **options)
    np.testing.assert_array_equal(a[0], b[0])
    np.testing.assert_array_equal(a[1], b[1])
    assert all(value >= 0 for value in timings.values())


@pytest.mark.parametrize(
    "options", [{"proposer": "unknown"}, {"forward_looking": True, "lookahead_dt": 0}]
)
def test_invalid_options(options):
    with pytest.raises(ValueError):
        evolve(model(), **options)
