use super::GeneralizedTableau;
use super::sparsevec::SparseVector;
use super::traits::{TGate, TableauIndex};
use crate::tableau::{CliffordExtensions, LossyMeasure, Reset};
use crate::traits::{Clifford, Depolarizing2, PauliError, TwoQubitPauliError};
use crate::{config::Config, traits::Depolarizing};
use itertools::Itertools;
use num::Integer;
use num::{
    Complex, One, ToPrimitive, Zero,
    complex::{Complex64, ComplexFloat},
};
use std::fmt::Debug;

pub trait RunStim {
    fn run_stim_string(&mut self, circuit: &str);
    fn parse_line(&mut self, line: &str, line_no: &usize);
    fn parse_instruction(instruction: &str) -> (&str, Option<Vec<&str>>, Option<Vec<f64>>);
    // fn run_stim_file(&mut self, file_path: &str);
}

impl<T, I, C> RunStim for GeneralizedTableau<T, I, C>
where
    T: Config,
    C: SparseVector<Complex<T::Coeff>, I> + std::fmt::Debug,
    T::Coeff: One
        + Zero
        + Clone
        + num::Num
        + ToPrimitive
        + std::fmt::Debug
        + std::ops::Mul<f64>
        + PartialOrd<f64>,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat,
    I: TableauIndex + Debug,
{
    fn run_stim_string(&mut self, circuit: &str) {
        for (i, line) in circuit.lines().enumerate() {
            let trimmed_line = line.trim();
            if trimmed_line.is_empty() {
                continue;
            }

            if trimmed_line.starts_with("#") {
                // comment
                continue;
            }

            self.parse_line(line.trim(), &i);
        }
    }

    fn parse_line(&mut self, line: &str, line_no: &usize) {
        let parts: Vec<&str> = line.split(" ").collect();

        let (instruction, tags, parens_args) = Self::parse_instruction(parts[0].trim());
        let addrs = parts[1..parts.len()]
            .iter()
            .map(|&c| c.parse::<usize>().unwrap());

        // instruction
        match instruction {
            "I" => {}

            "R" => {
                for addr in addrs {
                    self.reset(addr);
                }
            }

            "X" => {
                for addr in addrs {
                    self.x(addr);
                }
            }

            "Y" => {
                for addr in addrs {
                    self.y(addr);
                }
            }

            "Z" => {
                for addr in addrs {
                    self.z(addr);
                }
            }

            "H" => {
                for addr in addrs {
                    self.h(addr);
                }
            }

            "S" => {
                let is_t_gate = match tags {
                    Some(ts) => ts.len() == 1 && ts[0] == "T",
                    None => false,
                };
                if is_t_gate {
                    for addr in addrs {
                        self.t(addr);
                    }
                } else {
                    for addr in addrs {
                        self.s(addr);
                    }
                }
            }

            "S_DAG" => {
                let is_t_gate = match tags {
                    Some(ts) => ts.len() == 1 && ts[0] == "T",
                    None => false,
                };
                if is_t_gate {
                    for addr in addrs {
                        self.t_adj(addr);
                    }
                } else {
                    for addr in addrs {
                        self.s_adj(addr);
                    }
                }
            }

            "SQRT_Z" => {
                for addr in addrs {
                    self.s(addr);
                }
            }

            "SQRT_Z_DAG" => {
                for addr in addrs {
                    self.s_adj(addr);
                }
            }

            "SQRT_X" => {
                for addr in addrs {
                    self.sqrt_x(addr);
                }
            }

            "SQRT_X_DAG" => {
                for addr in addrs {
                    self.sqrt_x_adj(addr);
                }
            }

            "SQRT_Y" => {
                for addr in addrs {
                    self.sqrt_y(addr);
                }
            }

            "SQRT_Y_DAG" => {
                for addr in addrs {
                    self.sqrt_y_adj(addr);
                }
            }

            "CNOT" => {
                for (control, target) in addrs.tuples() {
                    self.cnot(control, target);
                }
            }

            "CZ" => {
                for (control, target) in addrs.tuples() {
                    self.cz(control, target);
                }
            }

            "DEPOLARIZE1" => {
                let ps = parens_args.unwrap();
                debug_assert_eq!(ps.len(), 1);
                let p = ps[0];
                for addr in addrs {
                    self.depolarize(addr, p.into());
                }
            }

            "DEPOLARIZE2" => {
                let ps = parens_args.unwrap();
                debug_assert_eq!(ps.len(), 1);
                let p = ps[0];
                for (control, target) in addrs.tuples() {
                    self.depolarize2(control, target, p.into());
                }
            }

            "PAULI_CHANNEL_1" => {
                let ps = parens_args.unwrap();
                debug_assert_eq!(ps.len(), 3);
                let ps_arr: [T::Coeff; 3] = [ps[0].into(), ps[1].into(), ps[2].into()];
                // FIXME: use reference in noise so we don't need to clone here
                for addr in addrs {
                    self.pauli_error(addr, ps_arr.clone());
                }
            }

            "PAULI_CHANNEL_2" => {
                let ps = parens_args.unwrap();
                debug_assert_eq!(ps.len(), 15);
                let ps_arr: [T::Coeff; 15] = std::array::from_fn(|i| ps[i].into());
                debug_assert!(addrs.clone().collect::<Vec<usize>>().len().is_even());
                for (control, target) in addrs.tuples() {
                    self.two_qubit_pauli_error(control.clone(), target.clone(), ps_arr.clone());
                }
            }

            "M" => {
                for addr in addrs {
                    self.measure(addr);
                }
            }

            _ => {
                panic!(
                    "Unknown circuit instruction {} in line {}",
                    parts[0], line_no
                );
            }
        }
    }

    fn parse_instruction(
        instruction_with_parens: &str,
    ) -> (&str, Option<Vec<&str>>, Option<Vec<f64>>) {
        let has_tags = instruction_with_parens.contains("[");
        let has_parens_args = instruction_with_parens.contains(")");

        if !has_tags && !has_parens_args {
            return (instruction_with_parens, None, None);
        } else if !has_parens_args {
            let parts: Vec<&str> = instruction_with_parens.split("[").collect();
            debug_assert_eq!(parts.len(), 2);
            let instruction = parts[0].trim();
            let tags: Vec<&str> = parts[1]
                .strip_suffix("]")
                .unwrap_or(parts[1])
                .split(",")
                .map(|c| c.trim())
                .collect();
            return (instruction, Some(tags), None);
        } else if !has_tags {
            let parts: Vec<&str> = instruction_with_parens.split("(").collect();
            debug_assert_eq!(parts.len(), 2);
            let instruction = parts[0].trim();
            let parens_args = parts[1].strip_suffix(")").unwrap_or(parts[1]).split(",");
            let ps = parens_args.map(|c| c.parse::<f64>().unwrap());
            return (instruction, None, Some(ps.collect()));
        } else {
            let parts: Vec<&str> = instruction_with_parens.split("[").collect();
            debug_assert_eq!(parts.len(), 2);
            let instruction = parts[0];
            let parts_tags_parens: Vec<&str> = parts[1].split("]").collect();
            debug_assert_eq!(parts_tags_parens.len(), 2);
            let tags = parts_tags_parens[0].split(",").map(|c| c.trim());
            let parens_args = parts_tags_parens[1];
            let ps = parens_args
                .trim()
                .strip_prefix("(")
                .unwrap_or(parens_args)
                .strip_suffix(")")
                .unwrap_or(parens_args)
                .split(",")
                .map(|c| c.parse::<f64>().unwrap());
            return (instruction, Some(tags.collect()), Some(ps.collect()));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::RunStim;
    use crate::config::indexmap::ByteFxHashF64;
    use crate::tableau::GeneralizedTableau;

    const test_program: &str = "
    R 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49
    H 18
    DEPOLARIZE1(0.4) 18
    I[R_Y(theta=0.32*pi)] 39
    PAULI_CHANNEL_1(0.3, 0.2, 0.1) 34
    H 13
    Z_ERROR(0.4) 5
    DEPOLARIZE2(0.5) 26 25
    I[U3(theta=0.34*pi, phi=0.21*pi, lambda=0.46*pi)] 31
    X_ERROR(1) 14
    PAULI_CHANNEL_2(0.1, 0.12, 0.2, 0.1, 0.05, 0.03, 0.02, 0.01, 0.005, 0.003, 0.002, 0.001, 0.0005, 0.0003, 0.0002) 40 30
    DEPOLARIZE1(0.4) 8
    SQRT_Y 36
    I[R_X(theta=0.31*pi)] 1
    H 6
    I[R_Y(theta=0.32*pi)] 25
    Z_ERROR(1) 2
    Y_ERROR(0.4) 1
    Z_ERROR(0.4) 23
    I[R_X(theta=0.31*pi)] 0
    S 38
    Y_ERROR(0.4) 8
    S 11
    I[R_X(theta=0.31*pi)] 23
    I[R_Z(theta=0.33*pi)] 22
    DEPOLARIZE2(0.5) 41 16
    SQRT_X 40
    H 12
    CX 23 25
    X_ERROR(1) 0
    I[R_X(theta=0.31*pi)] 46
    Y 49
    I[R_Y(theta=0.32*pi)] 19
    PAULI_CHANNEL_2(0.1, 0.12, 0.2, 0.1, 0.05, 0.03, 0.02, 0.01, 0.005, 0.003, 0.002, 0.001, 0.0005, 0.0003, 0.0002) 18 40
    SQRT_Y 36
    S 38
    DEPOLARIZE1(0.4) 6
    Y 13
    I[R_Y(theta=0.32*pi)] 20
    X_ERROR(1) 0
    Y 21
    X_ERROR(0.4) 26 18
    PAULI_CHANNEL_2(0.1, 0.12, 0.2, 0.1, 0.05, 0.03, 0.02, 0.01, 0.005, 0.003, 0.002, 0.001, 0.0005, 0.0003, 0.0002) 38 18
    Y_ERROR(0.4) 38
    DEPOLARIZE1(0.4) 28
    I[R_X(theta=0.31*pi)] 37
    PAULI_CHANNEL_1(0.3, 0.2, 0.1) 19
    X_ERROR(0.4) 7
    I[U3(theta=0.34*pi, phi=0.21*pi, lambda=0.46*pi)] 17
    I[R_Z(theta=0.33*pi)] 20
    Y 31
    I[U3(theta=0.34*pi, phi=0.21*pi, lambda=0.46*pi)] 15
    PAULI_CHANNEL_2(0.1, 0.12, 0.2, 0.1, 0.05, 0.03, 0.02, 0.01, 0.005, 0.003, 0.002, 0.001, 0.0005, 0.0003, 0.0002) 46 6
    CX 27 28 29 42
    Y 39
    DEPOLARIZE1(0.4) 40
    Y 18
    PAULI_CHANNEL_1(0.3, 0.2, 0.1) 5
    Z_ERROR(0.4) 35
    Z_ERROR(1) 14
    Y_ERROR(0.4) 35
    PAULI_CHANNEL_1(0.3, 0.2, 0.1) 46
    SQRT_X 27
    Z_ERROR(1) 2
    DEPOLARIZE1(0.4) 14
    S[T] 17
    CX 29 22
    S 14
    H 12
    SQRT_Y 27
    Z_ERROR(0.4) 16
    I[R_X(theta=0.31*pi)] 0
    H 4
    CX 12 21
    I[R_Z(theta=0.33*pi)] 2
    X_ERROR(0.4) 0
    Y_ERROR(0.4) 40
    Z_ERROR(0.4) 11
    SQRT_X 40 15
    PAULI_CHANNEL_2(0.1, 0.12, 0.2, 0.1, 0.05, 0.03, 0.02, 0.01, 0.005, 0.003, 0.002, 0.001, 0.0005, 0.0003, 0.0002) 41 48
    I[R_X(theta=0.31*pi)] 29
    SQRT_Y 6
    X_ERROR(0.4) 26
    I[R_X(theta=0.31*pi)] 49
    X_ERROR(1) 15
    Y_ERROR(0.4) 30
    I[U3(theta=0.34*pi, phi=0.21*pi, lambda=0.46*pi)] 10
    I[R_X(theta=0.31*pi)] 42
    I[R_Y(theta=0.32*pi)] 11
    S 25
    Y_ERROR(0.4) 49
    PAULI_CHANNEL_1(0.3, 0.2, 0.1) 22
    SQRT_Y 21
    I[R_Z(theta=0.33*pi)] 23
    CX 4 41
    H 32
    X_ERROR(1) 15
    X_ERROR(0.4) 28
    M 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49
    ";

    #[test]
    fn test_execution() {
        let n = 50;
        let mut tab: GeneralizedTableau<ByteFxHashF64<11>, u128> =
            GeneralizedTableau::new(n, 1e-10);

        tab.run_stim_string(test_program);
    }
}
