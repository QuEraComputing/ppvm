use super::GeneralizedTableau;
use super::sparsevec::SparseVector;
use super::traits::{TGate, TableauIndex};
use crate::tableau::{CliffordExtensions, LossyMeasure, Reset};
use crate::traits::{
    Clifford, CorrelatedLossChannel, Depolarizing2, LossChannel, PauliError, RotationOne,
    TwoQubitPauliError, U3Gate,
};
use crate::{config::Config, traits::Depolarizing};
use itertools::Itertools;
use num::Integer;
use num::{
    Complex, One, ToPrimitive, Zero,
    complex::{Complex64, ComplexFloat},
};
use std::collections::HashMap;
use std::fmt::Debug;

/// Split `s` by commas that are not inside parentheses.
fn split_commas_shallow(s: &str) -> Vec<&str> {
    let mut depth = 0usize;
    let mut start = 0;
    let mut result = Vec::new();
    for (i, c) in s.char_indices() {
        match c {
            '(' | '[' => depth += 1,
            ')' | ']' => depth -= 1,
            ',' if depth == 0 => {
                result.push(s[start..i].trim());
                start = i + 1;
            }
            _ => {}
        }
    }
    result.push(s[start..].trim());
    result
}

/// Parse an expression of the form `<coeff>*pi` or a plain float.
fn parse_pi_expr(s: &str) -> f64 {
    if let Some(coeff) = s.strip_suffix("*pi") {
        coeff.trim().parse::<f64>().unwrap() * std::f64::consts::PI
    } else {
        s.parse::<f64>().unwrap()
    }
}

pub trait RunStim {
    fn run_stim_string(&mut self, circuit: &str) -> HashMap<usize, Option<bool>>;
    fn parse_line(
        &mut self,
        line: &str,
        line_no: &usize,
        results: &mut HashMap<usize, Option<bool>>,
    );
    fn parse_instruction(
        line: &str,
        line_no: usize,
    ) -> (&str, Option<Vec<&str>>, Option<Vec<f64>>, &str);
    fn run_stim_file(&mut self, file_path: &str) -> HashMap<usize, Option<bool>> {
        let circuit = std::fs::read_to_string(file_path)
            .unwrap_or_else(|e| panic!("failed to read stim file {}: {}", file_path, e));
        self.run_stim_string(&circuit)
    }
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
    fn run_stim_string(&mut self, circuit: &str) -> HashMap<usize, Option<bool>> {
        let mut results = HashMap::new();
        for (i, line) in circuit.lines().enumerate() {
            let trimmed_line = line.trim();
            if trimmed_line.is_empty() {
                continue;
            }

            if trimmed_line.starts_with("#") {
                // comment
                continue;
            }

            self.parse_line(line.trim(), &i, &mut results);
        }
        results
    }

    fn parse_line(
        &mut self,
        line: &str,
        line_no: &usize,
        results: &mut HashMap<usize, Option<bool>>,
    ) {
        let (instruction, tags, parens_args, addr_part) = Self::parse_instruction(line, *line_no);
        let addrs = addr_part
            .split_whitespace()
            .map(|c| c.parse::<usize>().unwrap());

        // instruction
        match instruction {
            "I" => {
                if let Some([tag]) = tags.as_deref() {
                    if let Some(paren_start) = tag.find('(') {
                        let gate = tag[..paren_start].trim();
                        let inner = tag[paren_start + 1..].strip_suffix(')').unwrap();
                        let params: Vec<f64> = inner
                            .split(',')
                            .map(|p| parse_pi_expr(p.split('=').nth(1).unwrap().trim()))
                            .collect();
                        match gate {
                            "R_X" => {
                                for addr in addrs {
                                    self.rx(addr, params[0]);
                                }
                            }
                            "R_Y" => {
                                for addr in addrs {
                                    self.ry(addr, params[0]);
                                }
                            }
                            "R_Z" => {
                                for addr in addrs {
                                    self.rz(addr, params[0]);
                                }
                            }
                            "U3" => {
                                for addr in addrs {
                                    self.u3(
                                        addr,
                                        params[0].into(),
                                        params[1].into(),
                                        params[2].into(),
                                    );
                                }
                            }
                            _ => {}
                        }
                    }
                }
            }

            "I_ERROR" => match tags.as_deref() {
                Some(["loss"]) => {
                    let ps = parens_args.unwrap();
                    debug_assert_eq!(ps.len(), 1);
                    for addr in addrs {
                        self.loss_channel(addr, ps[0].into());
                    }
                }
                Some(["correlated_loss"]) => {
                    let ps = parens_args.unwrap();
                    let ps_arr: [T::Coeff; 3] = if ps.len() == 1 {
                        // simple correlated loss
                        [ps[0].into(), T::Coeff::zero(), T::Coeff::zero()]
                    } else {
                        debug_assert_eq!(ps.len(), 3);
                        [ps[0].into(), ps[1].into(), ps[2].into()]
                    };
                    for (control, target) in addrs.tuples() {
                        self.correlated_loss_channel(control, target, ps_arr.clone());
                    }
                }
                _ => {}
            },

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

            "H" | "H_XZ" => {
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

            "CX" | "CNOT" => {
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

            "X_ERROR" => {
                let ps = parens_args.unwrap();
                debug_assert_eq!(ps.len(), 1);
                let ps_arr: [T::Coeff; 3] = [ps[0].into(), T::Coeff::zero(), T::Coeff::zero()];
                for addr in addrs {
                    self.pauli_error(addr, ps_arr.clone());
                }
            }

            "Y_ERROR" => {
                let ps = parens_args.unwrap();
                debug_assert_eq!(ps.len(), 1);
                let ps_arr: [T::Coeff; 3] = [T::Coeff::zero(), ps[0].into(), T::Coeff::zero()];
                for addr in addrs {
                    self.pauli_error(addr, ps_arr.clone());
                }
            }

            "Z_ERROR" => {
                let ps = parens_args.unwrap();
                debug_assert_eq!(ps.len(), 1);
                let ps_arr: [T::Coeff; 3] = [T::Coeff::zero(), T::Coeff::zero(), ps[0].into()];
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
                    results.insert(addr, self.measure(addr));
                }
            }

            "MR" => {
                for addr in addrs {
                    let outcome = self.measure(addr);
                    if outcome == Some(true) {
                        self.x(addr);
                    }
                    results.insert(addr, outcome);
                }
            }

            // no-ops
            "DETECTOR" | "MPAD" | "OBSERVABLE_INCLUDE" | "QUBIT_COORDS" | "SHIFT_COORDS"
            | "TICK" => {}

            _ => {
                panic!(
                    "Unknown circuit instruction {} in line {}",
                    instruction, line_no
                );
            }
        }
    }

    fn parse_instruction(
        line: &str,
        line_no: usize,
    ) -> (&str, Option<Vec<&str>>, Option<Vec<f64>>, &str) {
        // Find the split between instruction token and addresses (first space at depth 0)
        let mut depth = 0usize;
        let mut split = None;
        for (i, c) in line.char_indices() {
            match c {
                '(' | '[' => depth += 1,
                ')' | ']' => depth -= 1,
                ' ' if depth == 0 => {
                    split = Some(i);
                    break;
                }
                _ => {}
            }
        }
        let (instr_token, addr_part) = match split {
            Some(p) => (line[..p].trim(), &line[p + 1..]),
            None => (line.trim(), ""),
        };

        if let Some(bracket_start) = instr_token.find('[') {
            let instruction = instr_token[..bracket_start].trim();
            let after_open = &instr_token[bracket_start + 1..];
            let bracket_end = after_open
                .find(']')
                .expect(&format!("unclosed [ in line {}", line_no));
            let tags: Vec<&str> = split_commas_shallow(&after_open[..bracket_end]);
            let after_bracket = after_open[bracket_end + 1..].trim();
            let parens_args = after_bracket
                .strip_prefix('(')
                .and_then(|s| s.strip_suffix(')'))
                .map(|inner| {
                    inner
                        .split(',')
                        .map(|c| c.trim().parse().unwrap())
                        .collect()
                });
            (instruction, Some(tags), parens_args, addr_part)
        } else if let Some(paren_start) = instr_token.find('(') {
            let instruction = instr_token[..paren_start].trim();
            let inner = instr_token[paren_start + 1..]
                .strip_suffix(')')
                .expect(&format!("unclosed ( in line {}", line_no));
            let ps = inner
                .split(',')
                .map(|c| c.trim().parse().unwrap())
                .collect();
            (instruction, None, Some(ps), addr_part)
        } else {
            (instr_token, None, None, addr_part)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::RunStim;
    use crate::config::indexmap::ByteFxHashF64;
    use crate::tableau::GeneralizedTableau;
    use crate::tableau::traits::LossyMeasure;

    const TEST_PROGRAM: &str = "
    R 0 1 2 3 4 5 6 7 8 9
    CX 0 7
    SQRT_Y 1
    H 5
    S[T] 3
    DEPOLARIZE2(0.5) 1 7
    I[R_Z(theta=0.33*pi)] 7
    DEPOLARIZE1(0.4) 7
    Z_ERROR(0.4) 5
    SQRT_X 6
    H 7
    S[T] 2
    I[R_X(theta=0.31*pi)] 6
    PAULI_CHANNEL_2(0.1, 0.12, 0.2, 0.1, 0.05, 0.03, 0.02, 0.01, 0.005, 0.003, 0.002, 0.001, 0.0005, 0.0003, 0.0002) 4 3
    Y_ERROR(0.4) 1
    CX 3 6
    DEPOLARIZE2(0.5) 0 4 4 1
    Z_ERROR(0.4) 8
    I[R_Z(theta=0.33*pi)] 8
    H 8
    X_ERROR(1) 3
    X_ERROR(0.4) 3
    Y 5
    DEPOLARIZE2(0.5) 2 3 3 1
    I[R_X(theta=0.31*pi)] 8
    I[U3(theta=0.34*pi, phi=0.21*pi, lambda=0.46*pi)] 5
    Y 1
    X_ERROR(0.4) 5
    X_ERROR(1) 2
    X_ERROR(0.4) 7
    X_ERROR(1) 3
    S[T] 8
    X_ERROR(1) 0
    SQRT_Y 5
    X_ERROR(0.4) 0
    Y 5
    X_ERROR(0.4) 4
    H 6
    S[T] 9
    CX 8 0 7 8
    SQRT_Y 0
    X_ERROR(1) 0
    H 2
    I[R_Z(theta=0.33*pi)] 0
    Y 5
    SQRT_Y 8
    S[T] 5
    PAULI_CHANNEL_1(0.3, 0.2, 0.1) 7
    DEPOLARIZE1(0.4) 7
    PAULI_CHANNEL_1(0.3, 0.2, 0.1) 2
    S[T] 4
    X_ERROR(0.4) 4
    S[T] 0
    I[R_Z(theta=0.33*pi)] 3
    CX 9 5
    H 2
    SQRT_Y 3
    H 8
    DEPOLARIZE2(0.5) 3 7
    Z_ERROR(0.4) 5
    DEPOLARIZE1(0.4) 7
    I[R_Y(theta=0.32*pi)] 0
    S[T] 8
    DEPOLARIZE1(0.4) 8
    Y 4
    Z_ERROR(0.4) 8
    I[R_X(theta=0.31*pi)] 6
    Z_ERROR(0.4) 6
    I[R_X(theta=0.31*pi)] 6
    H 1
    Z_ERROR(0.4) 5
    CX 5 4
    SQRT_X 6
    S[T] 9
    SQRT_X 8
    Y_ERROR(0.4) 8
    CX 6 5
    H 6
    PAULI_CHANNEL_2(0.1, 0.12, 0.2, 0.1, 0.05, 0.03, 0.02, 0.01, 0.005, 0.003, 0.002, 0.001, 0.0005, 0.0003, 0.0002) 2 1
    S 9
    PAULI_CHANNEL_1(0.3, 0.2, 0.1) 9
    X_ERROR(1) 6
    CX 0 9
    PAULI_CHANNEL_1(0.3, 0.2, 0.1) 1
    DEPOLARIZE2(0.5) 3 1
    S[T] 3
    Z_ERROR(0.4) 7
    X_ERROR(1) 6
    X_ERROR(0.4) 8 9
    H 7 1
    I[R_X(theta=0.31*pi)] 2
    S 9
    PAULI_CHANNEL_1(0.3, 0.2, 0.1) 8
    X_ERROR(1) 4 0
    H 4
    M 0 1 2 3 4 5 6 7 8 9
    ";

    #[test]
    fn test_execution() {
        let n = 50;
        let mut tab: GeneralizedTableau<ByteFxHashF64<7>, u128> = GeneralizedTableau::new(n, 1e-10);

        tab.run_stim_string(TEST_PROGRAM);
    }

    #[test]
    fn test_i_rotation_x() {
        // I[R_X(theta=1.0*pi)] should flip |0⟩ → |1⟩
        let mut tab: GeneralizedTableau<ByteFxHashF64<1>, usize> =
            GeneralizedTableau::new(1, 1e-10);
        tab.run_stim_string("I[R_X(theta=1.0*pi)] 0");
        assert!(tab.measure(0).unwrap());
    }

    #[test]
    fn test_i_rotation_y() {
        // I[R_Y(theta=1.0*pi)] should flip |0⟩ → |1⟩
        let mut tab: GeneralizedTableau<ByteFxHashF64<1>, usize> =
            GeneralizedTableau::new(1, 1e-10);
        tab.run_stim_string("I[R_Y(theta=1.0*pi)] 0");
        assert!(tab.measure(0).unwrap());
    }

    #[test]
    fn test_i_rotation_z() {
        // I[R_Z(theta=1.0*pi)] leaves |0⟩ unchanged (Z rotation only adds phase)
        let mut tab: GeneralizedTableau<ByteFxHashF64<1>, usize> =
            GeneralizedTableau::new(1, 1e-10);
        tab.run_stim_string("I[R_Z(theta=1.0*pi)] 0");
        assert!(!tab.measure(0).unwrap());
    }

    #[test]
    fn test_i_u3_flip() {
        // I[U3(theta=1.0*pi, phi=0.0*pi, lambda=0.0*pi)] = RY(π): |0⟩ → |1⟩
        let mut tab: GeneralizedTableau<ByteFxHashF64<1>, usize> =
            GeneralizedTableau::new(1, 1e-10);
        tab.run_stim_string("I[U3(theta=1.0*pi, phi=0.0*pi, lambda=0.0*pi)] 0");
        assert!(tab.measure(0).unwrap());
    }

    #[test]
    fn test_loss() {
        let test_program_loss = "I_ERROR[loss](1.0) 0";
        let mut tab: GeneralizedTableau<ByteFxHashF64<1>, usize> =
            GeneralizedTableau::new(1, 1e-10);

        tab.run_stim_string(test_program_loss);

        assert!(tab.is_lost[0])
    }

    #[test]
    fn test_run_stim_string_measurements() {
        let mut tab: GeneralizedTableau<ByteFxHashF64<1>, usize> =
            GeneralizedTableau::new(2, 1e-10);
        let results = tab.run_stim_string("X 0\nM 0 1");
        assert_eq!(results.get(&0), Some(&Some(true)));
        assert_eq!(results.get(&1), Some(&Some(false)));
    }

    #[test]
    fn test_run_stim_string_double_measurement() {
        // Second measurement of qubit 0 overwrites the first in the map
        let mut tab: GeneralizedTableau<ByteFxHashF64<1>, usize> =
            GeneralizedTableau::new(1, 1e-10);
        let results = tab.run_stim_string("X 0\nM 0\nM 0");
        assert_eq!(results.get(&0), Some(&Some(true)));
    }

    #[test]
    fn test_run_stim_file() {
        let path = std::env::temp_dir().join("ppvm_test_stim.stim");
        std::fs::write(&path, "X 0\nM 0 1").unwrap();

        let mut tab: GeneralizedTableau<ByteFxHashF64<1>, usize> =
            GeneralizedTableau::new(2, 1e-10);
        let results = tab.run_stim_file(path.to_str().unwrap());
        assert_eq!(results.get(&0), Some(&Some(true)));
        assert_eq!(results.get(&1), Some(&Some(false)));
    }
}
