use ppvm_runtime::{config::fxhash::ByteF64, prelude::*, strategy::CoefficientThreshold};
use ppvm_timeevolve::{JumpOp, LadderDirection, LadderOp, LindbladOp, RateMatrix, SolverConfig, solve::solve};

type S = ByteF64<1, CoefficientThreshold>;
const N: usize = 5;

fn main() {
    let ops: Vec<JumpOp<S>> = (0..N)
        .map(|i| JumpOp::Ladder(LadderOp { qubit: i, direction: LadderDirection::Lower }))
        .collect();
    let rates: Vec<Vec<f64>> = (0..N).map(|i|
        (0..N).map(|j| 1.0/(1.0+(i as f64-j as f64).abs())).collect()
    ).collect();
    let lindblad = LindbladOp::new(ops, RateMatrix::Dense(rates));
    let strat = CoefficientThreshold(1e-6);
    let mut initial: PauliSum<S> = PauliSum::builder().n_qubits(N).strategy(strat).build();
    let template = vec!['I'; N];
    for i in 0..N {
        let mut zi = template.clone(); zi[i] = 'Z';
        initial += (zi.into_iter().collect::<String>(), 1.0_f64);
    }
    for t in [0.01, 0.05, 0.1, 0.5, 1.0] {
        let (_, states) = solve(None, &lindblad, &initial, (0.0, t), &[t],
            |_, s| s.data().len(), SolverConfig::default());
        println!("t={t}: |p| = {}", states[0]);
    }
}
