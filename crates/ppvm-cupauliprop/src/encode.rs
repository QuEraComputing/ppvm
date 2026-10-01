// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Host-side conversion into cuPauliProp's data layout (no GPU needed).
//!
//! A cuPauliProp term is `2 * words` packed `u64`s: the X mask words followed
//! by the Z mask words, qubit `q` at bit `q % 64` of word `q / 64`, unused
//! high bits zero.

/// Packed `u64`s per X (or Z) mask; equals `cupaulipropGetNumPackedIntegers`.
pub fn num_words(n_qubits: usize) -> usize {
    n_qubits.div_ceil(64)
}

/// Append the packed X/Z words of a dense Pauli string (`I/X/Y/Z` per qubit).
pub fn push_packed_term(out: &mut Vec<u64>, term: &str, words: usize) {
    let start = out.len();
    out.resize(start + 2 * words, 0);
    let (x, z) = out[start..].split_at_mut(words);
    for (q, c) in term.chars().enumerate() {
        let (xb, zb) = match c {
            'I' => (false, false),
            'X' => (true, false),
            'Y' => (true, true),
            'Z' => (false, true),
            _ => panic!("invalid Pauli character {c:?} in {term:?}"),
        };
        let (w, bit) = (q / 64, 1u64 << (q % 64));
        if xb {
            x[w] |= bit;
        }
        if zb {
            z[w] |= bit;
        }
    }
}

/// ppvm `[p_x, p_y, p_z]` → cuPauliProp `[p_I, p_X, p_Y, p_Z]`.
pub fn pauli_channel_probs(p: [f64; 3]) -> [f64; 4] {
    [1.0 - p[0] - p[1] - p[2], p[0], p[1], p[2]]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn word_counts() {
        assert_eq!(num_words(1), 1);
        assert_eq!(num_words(64), 1);
        assert_eq!(num_words(65), 2);
        assert_eq!(num_words(512), 8);
    }

    #[test]
    fn packs_x_then_z_words() {
        let mut out = Vec::new();
        let mut term = vec!['I'; 70];
        term[0] = 'X';
        term[1] = 'Z';
        term[65] = 'Y';
        push_packed_term(&mut out, &term.iter().collect::<String>(), 2);
        assert_eq!(out, [0b01, 0b10, 0b10, 0b10]);
    }

    #[test]
    fn appends_terms_contiguously() {
        let mut out = Vec::new();
        push_packed_term(&mut out, "XI", 1);
        push_packed_term(&mut out, "IZ", 1);
        assert_eq!(out, [0b01, 0, 0, 0b10]);
    }

    #[test]
    fn channel_probs_prepend_identity() {
        assert_eq!(
            pauli_channel_probs([0.1, 0.2, 0.3]),
            [1.0 - 0.1 - 0.2 - 0.3, 0.1, 0.2, 0.3]
        );
    }
}
