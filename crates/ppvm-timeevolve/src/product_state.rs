use ppvm_runtime::prelude::{ACMapIter, Config, Pauli, PauliIter, PauliSum};

/// A separable initial state ρ₀ = ⊗ᵢ ρ₀⁽ⁱ⁾ encoded as per-qubit Bloch vectors.
///
/// Used to compute expectation values ⟨O(t)⟩ = Tr(ρ₀ O(t)) for a Heisenberg-picture
/// observable O(t). ρ₀ is never propagated — it is evaluated only at output checkpoints.
pub struct ProductState {
    /// `bloch[i] = [bx, by, bz]`.
    /// Convention: bz = +1 for |0⟩, bz = -1 for |1⟩.
    bloch: Vec<[f64; 3]>,
}

impl ProductState {
    /// All qubits in |0⟩: bz = +1.
    pub fn all_zero(n_qubits: usize) -> Self {
        ProductState { bloch: vec![[0.0, 0.0, 1.0]; n_qubits] }
    }

    /// All qubits in |1⟩: bz = -1.
    pub fn all_one(n_qubits: usize) -> Self {
        ProductState { bloch: vec![[0.0, 0.0, -1.0]; n_qubits] }
    }

    /// Computational basis state from a bit slice.
    /// `bits[i] = 0` → |0⟩ (bz=+1);  `bits[i] = 1` → |1⟩ (bz=-1).
    /// Panics if any element is not 0 or 1.
    pub fn bitstring(bits: &[u8]) -> Self {
        let bloch = bits
            .iter()
            .map(|&b| {
                let bz = match b {
                    0 => 1.0,
                    1 => -1.0,
                    _ => panic!("bitstring: bit value {b} is not 0 or 1"),
                };
                [0.0, 0.0, bz]
            })
            .collect();
        ProductState { bloch }
    }

    /// Arbitrary product state. `vectors[i] = [bx, by, bz]`.
    /// Pure states satisfy |b|² = 1; mixed states |b|² < 1.
    /// Prints a warning via `eprintln!` if any |bᵢ|² > 1 + 1e-9.
    pub fn bloch_vectors(vectors: Vec<[f64; 3]>) -> Self {
        for (i, &[bx, by, bz]) in vectors.iter().enumerate() {
            let norm_sq = bx * bx + by * by + bz * bz;
            if norm_sq > 1.0 + 1e-9 {
                eprintln!(
                    "ProductState::bloch_vectors: qubit {i} has |b|² = {norm_sq:.6} > 1 \
                     (not a valid density matrix). Proceeding anyway."
                );
            }
        }
        ProductState { bloch: vectors }
    }

    /// Returns the number of qubits this state is defined for.
    pub fn n_qubits(&self) -> usize {
        self.bloch.len()
    }

    /// Constructs a `ProductState` from a flat array `[bx₀,by₀,bz₀, bx₁,by₁,bz₁, …]`.
    /// Used by the native Python bridge, which passes Bloch vectors as a flat `Vec<f64>`.
    /// Panics if `flat.len()` is not divisible by 3.
    /// Called by the Python native bridge (`ppvm-python-native`), which passes Bloch
    /// vectors as a flat buffer. Must be `pub` so the bridge crate can access it.
    pub fn from_flat(flat: &[f64]) -> Self {
        assert!(flat.len().is_multiple_of(3), "from_flat: length must be divisible by 3");
        let bloch = flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
        ProductState { bloch }
    }

    /// Computes ⟨O⟩ = Tr(ρ₀ O) = Σ_α c_α · Π_i weight(α_i).
    ///
    /// Runs in O(|O| × n). Called only at save checkpoints — not in the RK hot loop.
    pub fn expectation<T>(&self, observable: &PauliSum<T>) -> f64
    where
        T: Config,
        for<'a> T::Map: ACMapIter<'a, Item = (&'a T::PauliWordType, &'a T::Coeff)>,
        T::Coeff: Into<f64> + Copy,
        T::PauliWordType: PauliIter,
    {
        observable.data().iter().map(|(word, coeff)| {
            let weight: f64 = word.iter().enumerate().map(|(i, pauli)| match pauli {
                Pauli::I => 1.0,
                Pauli::X => self.bloch[i][0],
                Pauli::Y => self.bloch[i][1],
                Pauli::Z => self.bloch[i][2],
                _ => 0.0, // Pauli::L (loss) does not contribute
            }).product();
            (*coeff).into() * weight
        }).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ppvm_runtime::prelude::{PauliSum, config::indexmap::ByteFxHashF64};

    type S = ByteFxHashF64<2>;

    fn sum2(terms: &[(&str, f64)]) -> PauliSum<S> {
        let mut s: PauliSum<S> = PauliSum::builder().n_qubits(2).build();
        for &(w, c) in terms {
            s += (w, c);
        }
        s
    }

    #[test]
    fn test_all_zero_n_qubits() {
        assert_eq!(ProductState::all_zero(3).n_qubits(), 3);
    }

    #[test]
    fn test_bitstring_encoding() {
        let ps = ProductState::bitstring(&[0, 1]);
        assert_eq!(ps.bloch[0], [0.0, 0.0, 1.0]);
        assert_eq!(ps.bloch[1], [0.0, 0.0, -1.0]);
    }

    #[test]
    fn test_from_flat_roundtrip() {
        let vecs = vec![[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]];
        let ps = ProductState::bloch_vectors(vecs.clone());
        let flat: Vec<f64> = ps.bloch.iter().flat_map(|v| v.iter().copied()).collect();
        let ps2 = ProductState::from_flat(&flat);
        assert_eq!(ps2.bloch, vecs);
    }

    #[test]
    #[should_panic(expected = "bit value 2 is not 0 or 1")]
    fn test_bitstring_invalid() {
        ProductState::bitstring(&[2]);
    }

    // ── expectation tests ────────────────────────────────────────────────────

    #[test]
    fn test_all_zero_expectation() {
        // all_zero: bz = +1 for both qubits.
        // ZI → weight = (+1)(1) = +1; IZ → (1)(+1) = +1; ZZ → (+1)(+1) = +1; II → 1.
        // Total = 1 + 1 + 1 + 1 = 4.
        let rho0 = ProductState::all_zero(2);
        let obs = sum2(&[("ZI", 1.0), ("IZ", 1.0), ("ZZ", 1.0), ("II", 1.0)]);
        assert!((rho0.expectation(&obs) - 4.0).abs() < 1e-14);
    }

    #[test]
    fn test_bitstring_10_expectation() {
        // bits=[1,0]: bz = [-1, +1].
        // ZI → (-1)(1) = -1; IZ → (1)(+1) = +1; ZZ → (-1)(+1) = -1; II → 1.
        // Total = -1 + 1 - 1 + 1 = 0.
        let rho0 = ProductState::bitstring(&[1, 0]);
        let obs = sum2(&[("ZI", 1.0), ("IZ", 1.0), ("ZZ", 1.0), ("II", 1.0)]);
        assert!((rho0.expectation(&obs) - 0.0).abs() < 1e-14);
    }

    #[test]
    fn test_bloch_x_plus_expectation() {
        // bx=1, by=bz=0 for both qubits.
        // Only {I,X}^⊗2 strings survive; Y and Z terms contribute 0.
        // XI → (1)(1) = 1; IX → (1)(1) = 1; XX → (1)(1) = 1; II → 1. Total = 4.
        // ZI → (0)(1) = 0 (bz=0).
        let rho0 = ProductState::bloch_vectors(vec![[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);
        let obs = sum2(&[("XI", 1.0), ("IX", 1.0), ("XX", 1.0), ("II", 1.0), ("ZI", 1.0)]);
        assert!((rho0.expectation(&obs) - 4.0).abs() < 1e-14);
    }

    #[test]
    fn test_xy_zero_for_bitstring() {
        // Bitstring states have bx=by=0; X and Y terms contribute 0.
        let rho0 = ProductState::bitstring(&[0, 1]);
        let obs = sum2(&[("XI", 1.0), ("YI", 1.0), ("IX", 1.0), ("IY", 1.0), ("XY", 1.0)]);
        assert!((rho0.expectation(&obs) - 0.0).abs() < 1e-14);
    }
}
