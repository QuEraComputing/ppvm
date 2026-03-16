use std::array;

use ppvm_runtime::{config::fxhash::ByteF64, prelude::*, strategy::CoefficientThreshold};
use ppvm_timeevolve::{CollapseOp, LindbladOp, RateMatrix, SolverConfig, solve::solve};

type S = ByteF64<1, CoefficientThreshold>;

fn main() {
    let n = 5;
    let gamma0 = 1.0;
    let d = 0.1;
    let tmax = 10.0;
    const TSTEPS: usize = 100;
    let dt = tmax / TSTEPS as f64;

    let mut rows: Vec<Vec<f64>> = Vec::with_capacity(n);

    for i in 0i32..(n as i32) {
        let mut row: Vec<f64> = Vec::with_capacity(n);
        for j in 0i32..(n as i32) {
            if i == j {
                row.push(gamma0.clone());
            } else {
                let gamma_ij = gamma0 / (1.0 + d * ((i - j) as f64).abs());
                row.push(gamma_ij);
            }
        }
        rows.push(row);
    }

    println!("Gamma: {:?}", rows);

    let gamma_mat = RateMatrix::Dense(rows);

    let mut c_ops: Vec<CollapseOp<S>> = Vec::with_capacity(n);

    let ppw = |pauli: &str,
               phase: u8|
     -> PhasedPauliWord<
        [u8; 1],
        fxhash::FxBuildHasher,
        PauliWord<[u8; 1], fxhash::FxBuildHasher>,
    > {
        PhasedPauliWord::build_from_word(
            PauliWord::<[u8; 1], fxhash::FxBuildHasher>::from(pauli),
            phase,
        )
    };
    let tmp = vec!['I'; n];
    for i in 0..n {
        let mut c = CollapseOp::<S>::new(n);
        let mut paulichars_x = tmp.clone();
        let mut paulichars_y = tmp.clone();
        paulichars_x[i] = 'X';
        paulichars_y[i] = 'Y';
        let pauli_x: String = paulichars_x.iter().collect();
        let pauli_y: String = paulichars_y.iter().collect();
        c.push(ppw(&pauli_x, 0), 1.0);
        c.push(ppw(&pauli_y, 3), 1.0);
        c_ops.push(c);
    }

    let lindblad_op = LindbladOp::new(c_ops, gamma_mat);

    let strat = CoefficientThreshold(1e-6);
    let mut initial: PauliSum<S> = PauliSum::builder().n_qubits(n).strategy(strat).build();
    for i in 0..n {
        let mut zi = tmp.clone();
        zi[i] = 'Z';
        let z: String = zi.iter().collect();
        initial += (z, 1.0);
    }

    let zero_pattern: PauliPattern = "Z?*".into();
    let fout = |_t: f64, p: &PauliSum<S>| {
        // let mut p_ = p.clone();
        // //
        // for i in 0..n {
        //     p_.x(i);
        // }
        // p_.trace(&zero_pattern)
        p.trace(&zero_pattern)
    };

    let save_at: [f64; TSTEPS] = array::from_fn(|i| dt * i as f64);
    let config = SolverConfig::default();
    let (ts, rs) = solve(
        None,
        &lindblad_op,
        &initial,
        (0.0, tmax),
        &save_at,
        fout,
        config,
    );

    println!("tout: {:?}", ts);
    println!("values: {:?}", rs);
}
