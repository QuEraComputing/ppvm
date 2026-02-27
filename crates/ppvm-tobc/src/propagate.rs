use num::Zero;
use ppvm_runtime::prelude::*;

/// dt * H(theta)
pub struct HamiltonianDelta<T: Config> {
    pub xx: Vec<Vec<T::Coeff>>,
    pub zz: Vec<Vec<T::Coeff>>,
    pub yy: Vec<Vec<T::Coeff>>,
    pub x: Vec<T::Coeff>,
    pub y: Vec<T::Coeff>,
    pub z: Vec<T::Coeff>,
}

impl<T: Config> HamiltonianDelta<T> {
    pub fn new(n: usize) -> Self {
        Self {
            xx: vec![vec![T::Coeff::zero(); n]; n],
            yy: vec![vec![T::Coeff::zero(); n]; n],
            zz: vec![vec![T::Coeff::zero(); n]; n],
            x: vec![T::Coeff::zero(); n],
            y: vec![T::Coeff::zero(); n],
            z: vec![T::Coeff::zero(); n],
        }
    }

    pub fn set_xx(&mut self, i: usize, j: usize, value: T::Coeff) {
        self.xx[i][j] = value;
    }

    pub fn set_yy(&mut self, i: usize, j: usize, value: T::Coeff) {
        self.yy[i][j] = value;
    }

    pub fn set_zz(&mut self, i: usize, j: usize, value: T::Coeff) {
        self.zz[i][j] = value;
    }

    pub fn set_x(&mut self, i: usize, value: T::Coeff) {
        self.x[i] = value;
    }

    pub fn set_y(&mut self, i: usize, value: T::Coeff) {
        self.y[i] = value;
    }

    pub fn set_z(&mut self, i: usize, value: T::Coeff) {
        self.z[i] = value;
    }

    pub fn evolve(&self, op: &mut PauliSum<T>, steps: usize)
    where
        T::Coeff: std::fmt::Debug,
    {
        let n = op.n_qubits();
        for _ in 0..steps {
            for i in 0..n {
                for j in i + 1..n {
                    op.rxx(i, j, self.xx[i][j].clone());
                    op.ryy(i, j, self.yy[i][j].clone());
                    op.rzz(i, j, self.zz[i][j].clone());
                }
                // println!("site: {}", i);
            }
            op.truncate();

            for i in 0..n {
                op.rx(i, self.x[i].clone());
                op.ry(i, self.y[i].clone());
                op.rz(i, self.z[i].clone());
                // println!("site: {}", i);
            }
            op.truncate();
        }
    }
}
