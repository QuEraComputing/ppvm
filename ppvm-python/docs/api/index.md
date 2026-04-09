# API Reference

ppvm exposes the following classes:

## Pauli Propagation

| Class | Description |
|-------|-------------|
| [`PauliSum`][ppvm.paulisum.PauliSum] | Pauli sum for standard (lossless) qubit simulations |
| [`LossyPauliSum`][ppvm.paulisum.LossyPauliSum] | Pauli sum with support for qubit loss channels |

Both classes share the same gate and noise channel interface. `LossyPauliSum`
additionally exposes `loss_channel` and `reset_loss_channel`.

Note that the loss simulation comes with a slight memory overhead to track
the information of which qubit was lost.
See [Simulating loss](../pauli_sum/index.md#simulating-loss) for details.

## Generalized Stabilizer Tableau

| Class | Description |
|-------|-------------|
| [`GeneralizedTableau`][ppvm.generalized_tableau.GeneralizedTableau] | Generalized stabilizer tableau for circuits with Clifford + non-Clifford gates, noise, and measurement |
| [`MeasurementResult`][ppvm.generalized_tableau.MeasurementResult] | Measurement outcome enum (`ZERO`, `ONE`, `LOST`) |
