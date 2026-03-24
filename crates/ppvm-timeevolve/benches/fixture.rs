use ppvm_runtime::{config::fxhash::ByteF64, prelude::*, strategy::CoefficientThreshold};
use ppvm_timeevolve::{CollapseOp, JumpOp, LindbladOp, RateMatrix};

pub type S = ByteF64<1, CoefficientThreshold>;

const N: usize = 5;

/// Build the benchmark Lindblad operator:
/// - 5 lowering operators c_i = X_i + iY_i (i = 0..4)
/// - Dense 5×5 rate matrix γ_ij = 1 / (1 + |i − j|)
pub fn build_lindblad() -> LindbladOp<S> {
    let ppw = |pauli: &str, phase: u8|
        -> PhasedPauliWord<[u8; 1], fxhash::FxBuildHasher, PauliWord<[u8; 1], fxhash::FxBuildHasher>>
    {
        PhasedPauliWord::build_from_word(
            PauliWord::<[u8; 1], fxhash::FxBuildHasher>::from(pauli),
            phase,
        )
    };

    let template = vec!['I'; N];
    let mut c_ops: Vec<CollapseOp<S>> = Vec::with_capacity(N);
    for i in 0..N {
        let mut c = CollapseOp::<S>::new(N);
        let mut px = template.clone();
        let mut py = template.clone();
        px[i] = 'X';
        py[i] = 'Y';
        let sx: String = px.into_iter().collect();
        let sy: String = py.into_iter().collect();
        c.push(ppw(&sx, 0), 1.0); // X_i (phase 0)
        c.push(ppw(&sy, 1), 1.0); // iY_i (phase 1 = i)
        c_ops.push(c);
    }

    // Dense 5×5 rate matrix: γ_ij = 1 / (1 + |i − j|)
    let rates: Vec<Vec<f64>> = (0..N)
        .map(|i| {
            (0..N)
                .map(|j| 1.0 / (1.0 + (i as f64 - j as f64).abs()))
                .collect()
        })
        .collect();

    LindbladOp::new(c_ops.into_iter().map(JumpOp::Generic).collect(), RateMatrix::Dense(rates))
}

/// Build the benchmark initial state: P = Σ_i Z_i, threshold 1e-6.
pub fn build_initial() -> PauliSum<S> {
    let strat = CoefficientThreshold(1e-6);
    let mut p: PauliSum<S> = PauliSum::builder().n_qubits(N).strategy(strat).build();
    let template = vec!['I'; N];
    for i in 0..N {
        let mut zi = template.clone();
        zi[i] = 'Z';
        let sz: String = zi.into_iter().collect();
        p += (sz, 1.0_f64);
    }
    p
}
