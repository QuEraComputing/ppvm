import math

import pytest

from ppvm import PauliSum


def test_basics():
    state = PauliSum(n_qubits=2, initial_terms=["ZZ"], coefficients=[1.0])  # ZZ

    state.cnot(0, 1)
    state.h(0)

    assert str(state) == "1.000 * IZ"
    assert state.overlap_with_zero() == 1.0


def test_noise():
    state = PauliSum(n_qubits=2, initial_terms=["IZ"], coefficients=[1.0])  # |00><00|

    error_probs = {"ZZ": 0.1, "XX": 0.2}
    error_probs_list = state.two_qubit_pauli_error_probabilities(error_probs)
    state.two_qubit_pauli_error(0, 1, p=error_probs_list)


def test_large_state():
    n = 200
    weight = 80

    terms = ["".join(["Z" if i == j else "I" for i in range(n)]) for j in range(n)]
    large_state = PauliSum(n_qubits=n, max_pauli_weight=weight, initial_terms=terms)

    for i in reversed(range(1, n)):
        large_state.cnot(i - 1, i)

    large_state.h(0)

    assert large_state.overlap_with_zero() == 0.0


def test_copy():

    state = PauliSum(n_qubits=2, initial_terms=["ZZ"], coefficients=[1.0])  # ZZ

    state.cnot(0, 1)
    state.h(0)

    assert str(state) == "1.000 * IZ"
    assert state.overlap_with_zero() == 1.0

    tmp = state.copy()
    assert tmp == state

    assert len(tmp) == len(state) == 1

    assert tmp.terms == [("IZ", 1.0)]


def test_weights():
    state = PauliSum(n_qubits=2, initial_terms=["IZ"])

    assert state.current_max_weight() == 1

    state.cnot(0, 1)
    state.h(0)

    assert state.current_max_weight() == 2

    state2 = PauliSum(n_qubits=2, initial_terms=["ZX", "IY"])
    weights = state2.weights()

    weights.sort(key=lambda w: w[1])
    assert weights == [("IY", 1), ("ZX", 2)]


def test_overlap():
    state = PauliSum(n_qubits=2, initial_terms=["IZ"])

    assert state.overlap(state) == 1.0

    state2 = PauliSum(n_qubits=2, initial_terms=["IX"])
    assert state.overlap(state2) == 0.0


def test_new():
    s = PauliSum.new(2, "IX")
    assert len(s) == 1
    assert s.terms == [("IX", 1.0)]

    s = PauliSum.new(2, ("IX", 0.25))
    assert len(s) == 1
    assert s.terms == [("IX", 0.25)]

    s = PauliSum.new(2, "X1")
    assert s == PauliSum.new(2, "IX")

    s = PauliSum.new(3, [("Y1", 0.1), "ZIZ"])
    assert len(s) == 2
    assert set(s.terms) == {("IYI", 0.1), ("ZIZ", 1.0)}

    n = 17
    s = PauliSum.new(n, [f"Z{i}" for i in range(n)])
    assert len(s) == n
    assert set(s.terms) == {
        ("".join(["Z" if i == j else "I" for i in range(n)]), 1.0) for j in range(n)
    }

    with pytest.raises(ValueError, match="out of range"):
        PauliSum.new(2, "X2")

    n = 5
    terms = ["Z0", "Z1Z2"]
    s = PauliSum.new(n, terms)


def test_gate_methods():
    """Each gate is verified by a known Heisenberg-picture transformation P → G P G†."""
    PI = math.pi

    def t(state):
        return dict(state.terms)

    # x: X commutes with X (invariant); anticommutes with Z (sign flip)
    s = PauliSum(n_qubits=2, initial_terms=["XI"], coefficients=[1.0])
    s.x(0)
    assert t(s) == {"XI": 1.0}

    s = PauliSum(n_qubits=2, initial_terms=["ZI"], coefficients=[1.0])
    s.x(0)
    assert pytest.approx(t(s)["ZI"]) == -1.0

    # y: Y commutes with Y (invariant)
    s = PauliSum(n_qubits=2, initial_terms=["IY"], coefficients=[1.0])
    s.y(1)
    assert t(s) == {"IY": 1.0}

    # z: Z commutes with Z (invariant)
    s = PauliSum(n_qubits=2, initial_terms=["IZ"], coefficients=[1.0])
    s.z(1)
    assert t(s) == {"IZ": 1.0}

    # h: Z → X, X → Z (self-inverse)
    s = PauliSum(n_qubits=2, initial_terms=["ZI"], coefficients=[1.0])
    s.h(0)
    assert t(s) == {"XI": 1.0}

    s = PauliSum(n_qubits=2, initial_terms=["XI"], coefficients=[1.0])
    s.h(0)
    assert t(s) == {"ZI": 1.0}

    # s: X → -Y
    s = PauliSum(n_qubits=2, initial_terms=["XI"], coefficients=[1.0])
    s.s(0)
    assert t(s) == {"YI": -1.0}

    # s_dag: Y → -X (inverse of S)
    s = PauliSum(n_qubits=2, initial_terms=["YI"], coefficients=[1.0])
    s.s_dag(0)
    assert t(s) == {"XI": -1.0}

    # cnot(ctrl=0, tgt=1): IZ → ZZ; IX → IX (X on target is invariant)
    s = PauliSum(n_qubits=2, initial_terms=["IZ"], coefficients=[1.0])
    s.cnot(0, 1)
    assert t(s) == {"ZZ": 1.0}

    s = PauliSum(n_qubits=2, initial_terms=["IX"], coefficients=[1.0])
    s.cnot(0, 1)
    assert t(s) == {"IX": 1.0}

    # cz: IX → ZX, XI → XZ (symmetric)
    s = PauliSum(n_qubits=2, initial_terms=["IX"], coefficients=[1.0])
    s.cz(0, 1)
    assert t(s) == {"ZX": 1.0}

    s = PauliSum(n_qubits=2, initial_terms=["XI"], coefficients=[1.0])
    s.cz(0, 1)
    assert t(s) == {"XZ": 1.0}

    # ry(−π/2): ZI → XI  [cos·Z − sin·X at θ=−π/2 → 0·Z − (−1)·X = X]
    s = PauliSum(n_qubits=2, initial_terms=["ZI"], coefficients=[1.0])
    s.ry(0, theta=-PI / 2)
    assert pytest.approx(t(s).get("XI", 0.0)) == 1.0

    # rx(π/2): ZI → YI  [cos·Z + sin·Y at θ=π/2 → Y]
    s = PauliSum(n_qubits=2, initial_terms=["ZI"], coefficients=[1.0])
    s.rx(0, theta=PI / 2)
    assert pytest.approx(t(s).get("YI", 0.0)) == 1.0

    # rz(−π/2): XI → YI  [cos·X − sin·Y at θ=−π/2 → 0·X − (−1)·Y = Y]
    s = PauliSum(n_qubits=2, initial_terms=["XI"], coefficients=[1.0])
    s.rz(0, theta=-PI / 2)
    assert pytest.approx(t(s).get("YI", 0.0)) == 1.0

    # r(axis_angle=0, θ=π/2) = rx(π/2): ZI → YI
    s = PauliSum(n_qubits=2, initial_terms=["ZI"], coefficients=[1.0])
    s.r(0, 0.0, PI / 2)
    assert pytest.approx(t(s).get("YI", 0.0)) == 1.0

    # r(axis_angle=π/2, θ=π/2) = ry(π/2): ZI → −XI
    s = PauliSum(n_qubits=2, initial_terms=["ZI"], coefficients=[1.0])
    s.r(0, PI / 2, PI / 2)
    assert pytest.approx(t(s).get("XI", 0.0)) == -1.0

    # rxx(π/2): ZI → YX  [cos·ZI + sin·YX at θ=π/2 → YX]
    s = PauliSum(n_qubits=2, initial_terms=["ZI"], coefficients=[1.0])
    s.rxx(0, 1, theta=PI / 2)
    assert pytest.approx(t(s).get("YX", 0.0)) == 1.0

    # ryy(−π/2): ZI → XY  [cos·ZI − sin·XY at θ=−π/2 → 0·ZI − (−1)·XY = XY]
    s = PauliSum(n_qubits=2, initial_terms=["ZI"], coefficients=[1.0])
    s.ryy(0, 1, theta=-PI / 2)
    assert pytest.approx(t(s).get("XY", 0.0)) == 1.0

    # rzz(−π/2): XI → YZ  [cos·XI − sin·YZ at θ=−π/2 → 0·XI − (−1)·YZ = YZ]
    s = PauliSum(n_qubits=2, initial_terms=["XI"], coefficients=[1.0])
    s.rzz(0, 1, theta=-PI / 2)
    assert pytest.approx(t(s).get("YZ", 0.0)) == 1.0

    # sqrt_x (= HSH): X → X, Y → -Z, Z → Y
    s = PauliSum(n_qubits=2, initial_terms=["XI"], coefficients=[1.0])
    s.sqrt_x(0)
    assert t(s) == {"XI": 1.0}

    s = PauliSum(n_qubits=2, initial_terms=["YI"], coefficients=[1.0])
    s.sqrt_x(0)
    assert pytest.approx(t(s).get("ZI", 0.0)) == -1.0

    s = PauliSum(n_qubits=2, initial_terms=["ZI"], coefficients=[1.0])
    s.sqrt_x(0)
    assert pytest.approx(t(s).get("YI", 0.0)) == 1.0

    # sqrt_x_dag (= HS†H): X → X, Y → Z, Z → -Y
    s = PauliSum(n_qubits=2, initial_terms=["XI"], coefficients=[1.0])
    s.sqrt_x_dag(0)
    assert t(s) == {"XI": 1.0}

    s = PauliSum(n_qubits=2, initial_terms=["YI"], coefficients=[1.0])
    s.sqrt_x_dag(0)
    assert pytest.approx(t(s).get("ZI", 0.0)) == 1.0

    s = PauliSum(n_qubits=2, initial_terms=["ZI"], coefficients=[1.0])
    s.sqrt_x_dag(0)
    assert pytest.approx(t(s).get("YI", 0.0)) == -1.0

    # sqrt_y (= S·sqrt_x·S†): X → Z, Y → Y, Z → -X
    s = PauliSum(n_qubits=2, initial_terms=["XI"], coefficients=[1.0])
    s.sqrt_y(0)
    assert pytest.approx(t(s).get("ZI", 0.0)) == 1.0

    s = PauliSum(n_qubits=2, initial_terms=["YI"], coefficients=[1.0])
    s.sqrt_y(0)
    assert t(s) == {"YI": 1.0}

    s = PauliSum(n_qubits=2, initial_terms=["ZI"], coefficients=[1.0])
    s.sqrt_y(0)
    assert pytest.approx(t(s).get("XI", 0.0)) == -1.0

    # sqrt_y_dag (= S·sqrt_x_dag·S†): X → -Z, Y → Y, Z → X
    s = PauliSum(n_qubits=2, initial_terms=["XI"], coefficients=[1.0])
    s.sqrt_y_dag(0)
    assert pytest.approx(t(s).get("ZI", 0.0)) == -1.0

    s = PauliSum(n_qubits=2, initial_terms=["YI"], coefficients=[1.0])
    s.sqrt_y_dag(0)
    assert t(s) == {"YI": 1.0}

    s = PauliSum(n_qubits=2, initial_terms=["ZI"], coefficients=[1.0])
    s.sqrt_y_dag(0)
    assert pytest.approx(t(s).get("XI", 0.0)) == 1.0


def test_noise_methods():
    """Noise channels act as super-operators E(P) = Σ_k p_k N_k† P N_k + (1−Σp_k)P.

    A Pauli term P scales by +1 for each noise operator it commutes with and −1
    for each it anticommutes with, weighted by the probabilities.  The net
    scaling for axis A is (1 − 2p_B − 2p_C) where {B, C} are the other axes.
    """

    def t(state):
        return dict(state.terms)

    # --- pauli_error ---

    # zero probabilities: no change
    s = PauliSum(n_qubits=2, initial_terms=["ZI"], coefficients=[1.0])
    s.pauli_error(0, p=[0.0, 0.0, 0.0])
    assert t(s) == {"ZI": 1.0}

    # pz = 0.5: Z commutes with Z-error → unchanged; X anticommutes → zeroed
    s = PauliSum(n_qubits=2, initial_terms=["ZI"], coefficients=[1.0])
    s.pauli_error(0, p=[0.0, 0.0, 0.5])
    assert t(s) == {"ZI": 1.0}  # Z → (1 − 2·0 − 2·0)·Z = Z

    s = PauliSum(n_qubits=2, initial_terms=["XI"], coefficients=[1.0])
    s.pauli_error(0, p=[0.0, 0.0, 0.5])
    assert len(s) == 0  # X → (1 − 0 − 2·0.5)·X = 0

    # px = 0.5: X commutes with X-error → unchanged; Z anticommutes → zeroed
    s = PauliSum(n_qubits=2, initial_terms=["XI"], coefficients=[1.0])
    s.pauli_error(0, p=[0.5, 0.0, 0.0])
    assert t(s) == {"XI": 1.0}  # X → (1 − 2·0 − 2·0)·X = X

    s = PauliSum(n_qubits=2, initial_terms=["ZI"], coefficients=[1.0])
    s.pauli_error(0, p=[0.5, 0.0, 0.0])
    assert len(s) == 0  # Z → (1 − 2·0.5 − 0)·Z = 0

    # symmetric px=py=pz=1/4: all non-identity Pauli terms vanish
    s = PauliSum(n_qubits=2, initial_terms=["ZI"], coefficients=[1.0])
    s.pauli_error(0, p=[0.25, 0.25, 0.25])
    assert len(s) == 0  # Z → (1 − 0.5 − 0.5)·Z = 0

    # pauli_error on qubit 1 does not affect a term with support only on qubit 0
    s = PauliSum(n_qubits=2, initial_terms=["ZI"], coefficients=[1.0])
    s.pauli_error(1, p=[0.5, 0.0, 0.0])
    assert t(s) == {"ZI": 1.0}

    # --- two_qubit_pauli_error ---

    # zero probabilities: no change
    s = PauliSum(n_qubits=2, initial_terms=["IZ"], coefficients=[1.0])
    s.two_qubit_pauli_error(0, 1, p=[0.0] * 15)
    assert t(s) == {"IZ": 1.0}

    # p_IX = 1: IX anticommutes with IZ (X anticommutes with Z on qubit 1) → sign flip
    p = [0.0] * 15
    p[0] = 1.0  # IX
    s = PauliSum(n_qubits=2, initial_terms=["IZ"], coefficients=[1.0])
    s.two_qubit_pauli_error(0, 1, p=p)
    assert pytest.approx(t(s).get("IZ", 0.0)) == -1.0

    # p_IX = p_"IY" = p_IZ = 0.25 (depolarize qubit 1): kills IZ, leaves ZI intact
    p = [0.0] * 15
    p[0] = p[1] = p[2] = 0.25  # IX, "IY", IZ

    s = PauliSum(n_qubits=2, initial_terms=["IZ"], coefficients=[1.0])
    s.two_qubit_pauli_error(0, 1, p=p)
    assert len(s) == 0  # IZ → 0

    s = PauliSum(n_qubits=2, initial_terms=["ZI"], coefficients=[1.0])
    s.two_qubit_pauli_error(0, 1, p=p)
    assert t(s) == {"ZI": 1.0}  # ZI commutes with all I⊗{X,Y,Z} → unchanged

    # --- depolarize ---

    # zero probability: no change
    s = PauliSum(n_qubits=2, initial_terms=["ZI"], coefficients=[1.0])
    s.depolarize1(0, p=0.0)
    assert t(s) == {"ZI": 1.0}

    # depolarize1(p) with px=py=pz=p/3: non-identity Paulis scale by 1 − 4p/3
    # (anticommutes with 2 axes, each with prob p/3 → scaling = 1 − 2·p/3 − 2·p/3)
    p = 0.1
    s = PauliSum(n_qubits=2, initial_terms=["ZI"], coefficients=[1.0])
    s.depolarize1(0, p=p)
    assert pytest.approx(t(s).get("ZI", 0.0)) == 1.0 - 4 * p / 3

    # at p=3/4 the channel is completely depolarizing: all Paulis vanish
    s = PauliSum(n_qubits=2, initial_terms=["ZI"], coefficients=[1.0])
    s.depolarize1(0, p=0.75)
    assert len(s) == 0

    # depolarize on qubit 1 does not affect a term with support only on qubit 0
    s = PauliSum(n_qubits=2, initial_terms=["ZI"], coefficients=[1.0])
    s.depolarize1(1, p=0.3)
    assert t(s) == {"ZI": 1.0}

    # --- depolarize2 ---

    # zero probability: no change
    s = PauliSum(n_qubits=2, initial_terms=["ZZ"], coefficients=[1.0])
    s.depolarize2(0, 1, p=0.0)
    assert t(s) == {"ZZ": 1.0}

    # ZZ anticommutes with 8 of 15 non-identity two-qubit Paulis → scales by 1 − 16p/15
    p = 0.1
    s = PauliSum(n_qubits=2, initial_terms=["ZZ"], coefficients=[1.0])
    s.depolarize2(0, 1, p=p)
    assert pytest.approx(t(s).get("ZZ", 0.0)) == 1.0 - 16 * p / 15

    # --- amplitude_damping ---

    # gamma=0: no change
    s = PauliSum(n_qubits=2, initial_terms=["ZI"], coefficients=[1.0])
    s.amplitude_damping(0, 0.0)
    assert t(s) == {"ZI": 1.0}

    # E†[Z] = (1−γ)Z + γI: longitudinal Pauli branches into Z and I
    gamma = 0.2
    s = PauliSum(n_qubits=2, initial_terms=["ZI"], coefficients=[1.0])
    s.amplitude_damping(0, gamma)
    assert pytest.approx(t(s).get("ZI", 0.0)) == 1.0 - gamma
    assert pytest.approx(t(s).get("II", 0.0)) == gamma

    # E†[X] = √(1−γ)·X: transverse Paulis decay
    s = PauliSum(n_qubits=2, initial_terms=["XI"], coefficients=[1.0])
    s.amplitude_damping(0, gamma)
    assert pytest.approx(t(s).get("XI", 0.0)) == (1.0 - gamma) ** 0.5

    # T₂ = 2T₁ physics: transverse decay factor squared equals longitudinal decay
    assert pytest.approx(((1.0 - gamma) ** 0.5) ** 2) == 1.0 - gamma

    # amplitude_damping on qubit 1 does not affect a term with support only on qubit 0
    s = PauliSum(n_qubits=2, initial_terms=["ZI"], coefficients=[1.0])
    s.amplitude_damping(1, gamma)
    assert t(s) == {"ZI": 1.0}
