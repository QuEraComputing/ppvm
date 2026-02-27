use crate::{
    config::Config, phase::PhasedPauliWord, sum::PauliSum, traits::Clifford, traits::PauliStorage,
};
use std::hash::BuildHasher;

macro_rules! map_scale {
    ($name:ident, $($index:ident),*) => {
        fn $name(&mut self, $($index: usize),*) {
            self.scale(|k, v| {
                let mut p: PhasedPauliWord<T::Storage, T::BuildHasher, <T as Config>::PauliWordType> = k.clone().into();
                p.$name($($index),*);
                if !p.is_positive() {
                    *v *= -1.0;
                }
            })
        }
    };
}

macro_rules! map_word {
    ($name:ident, $($index:ident),*) => {
        fn $name(&mut self, $($index: usize),*) {
            self.map_add(|k, v| {
                let mut p: PhasedPauliWord<T::Storage, T::BuildHasher, <T as Config>::PauliWordType> = k.clone().into();
                p.$name($($index),*);
                if p.is_positive() {
                    (p.word, v.clone())
                } else {
                    (p.word, -v.clone())
                }
            })
        }
    };
}

// NOTE: impl for PauliWord only; not a blanket, since PhasedPauliWord Clifford also isn't
impl<S, H, T> Clifford for PauliSum<T>
where
    S: PauliStorage,
    H: BuildHasher + Clone + Default,
    T: Config<Storage = S, BuildHasher = H>,
    T::PauliWordType: Clifford,
{
    map_scale!(x, index);
    map_scale!(y, index);
    map_scale!(z, index);
    map_word!(h, index);
    map_word!(s, index);
    map_word!(s_adj, index);
    map_word!(cnot, a, b);
    map_word!(cz, a, b);
}

#[cfg(test)]
mod tests {
    use crate::prelude::*;

    type PS = PauliSum<config::indexmap::ByteFxHashF64<1>>;

    fn ps(term: &str, coeff: f64) -> PS {
        let mut p: PS = PauliSum::builder().n_qubits(1).build();
        p += (term, coeff);
        p
    }

    // sqrt_x = H S H  (backward convention P → G† P G)
    // Z→Y, Y→−Z, X→X (invariant)

    #[test]
    fn sqrt_x_z_to_y() {
        let mut p = ps("Z", 1.0);
        p.sqrt_x(0);
        assert_eq!(p, ps("Y", 1.0));
    }

    #[test]
    fn sqrt_x_y_to_neg_z() {
        let mut p = ps("Y", 1.0);
        p.sqrt_x(0);
        assert_eq!(p, ps("Z", -1.0));
    }

    #[test]
    fn sqrt_x_x_invariant() {
        let mut p = ps("X", 1.0);
        p.sqrt_x(0);
        assert_eq!(p, ps("X", 1.0));
    }

    // sqrt_x_adj = H S† H
    // Z→−Y, Y→Z, X→X (invariant)

    #[test]
    fn sqrt_x_adj_z_to_neg_y() {
        let mut p = ps("Z", 1.0);
        p.sqrt_x_adj(0);
        assert_eq!(p, ps("Y", -1.0));
    }

    #[test]
    fn sqrt_x_adj_y_to_z() {
        let mut p = ps("Y", 1.0);
        p.sqrt_x_adj(0);
        assert_eq!(p, ps("Z", 1.0));
    }

    #[test]
    fn sqrt_x_adj_x_invariant() {
        let mut p = ps("X", 1.0);
        p.sqrt_x_adj(0);
        assert_eq!(p, ps("X", 1.0));
    }

    // sqrt_x and sqrt_x_adj are mutual inverses
    #[test]
    fn sqrt_x_sqrt_x_adj_inverse() {
        let mut p = ps("Z", 1.0);
        p.sqrt_x(0);
        p.sqrt_x_adj(0);
        assert_eq!(p, ps("Z", 1.0));
    }

    // (sqrt_x)² conjugates like X: Z→−Z
    #[test]
    fn sqrt_x_squared_acts_as_x() {
        let mut p = ps("Z", 1.0);
        p.sqrt_x(0);
        p.sqrt_x(0);
        assert_eq!(p, ps("Z", -1.0));
    }

    // sqrt_y = S (H S H) S†  and  sqrt_y_adj = S† (H S† H) S
    // Both have identical Pauli conjugation (differ only by global phase):
    // Z→−X, X→Z, Y→Y (invariant)

    #[test]
    fn sqrt_y_x_to_z() {
        let mut p = ps("X", 1.0);
        p.sqrt_y(0);
        assert_eq!(p, ps("Z", 1.0));
    }

    #[test]
    fn sqrt_y_z_to_neg_x() {
        let mut p = ps("Z", 1.0);
        p.sqrt_y(0);
        assert_eq!(p, ps("X", -1.0));
    }

    #[test]
    fn sqrt_y_y_invariant() {
        let mut p = ps("Y", 1.0);
        p.sqrt_y(0);
        assert_eq!(p, ps("Y", 1.0));
    }

    #[test]
    fn sqrt_y_adj_x_to_z() {
        let mut p = ps("X", 1.0);
        p.sqrt_y_adj(0);
        assert_eq!(p, ps("Z", 1.0));
    }

    #[test]
    fn sqrt_y_adj_z_to_neg_x() {
        let mut p = ps("Z", 1.0);
        p.sqrt_y_adj(0);
        assert_eq!(p, ps("X", -1.0));
    }

    // (sqrt_y)² conjugates like Y: Z→−Z
    #[test]
    fn sqrt_y_squared_acts_as_y() {
        let mut p = ps("Z", 1.0);
        p.sqrt_y(0);
        p.sqrt_y(0);
        assert_eq!(p, ps("Z", -1.0));
    }
}
