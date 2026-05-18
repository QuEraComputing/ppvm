use crate::config::Config;

/// The non-Clifford `T` gate and its adjoint.
///
/// `T = diag(1, e^{iπ/4})`. Implemented by the simulator backends; see
/// the example in [`ppvm_tableau`](https://docs.rs/ppvm-tableau).
pub trait TGate<T: Config> {
    /// Apply `T` (`diag(1, e^{iπ/4})`) to qubit `addr0`.
    fn t(&mut self, addr0: usize);
    /// Apply `T†` to qubit `addr0`.
    fn t_adj(&mut self, addr0: usize);
}
