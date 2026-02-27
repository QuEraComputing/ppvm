# API Reference

ppvm exposes two classes:

| Class | Description |
|-------|-------------|
| [`PauliSum`][ppvm.paulisum.PauliSum] | Pauli sum for standard (lossless) qubit simulations |
| [`LossyPauliSum`][ppvm.paulisum.LossyPauliSum] | Pauli sum with support for qubit loss channels |

Both classes share the same gate and noise channel interface. `LossyPauliSum`
additionally exposes `loss_channel` and `reset_loss_channel`.

Note that the loss simulation comes with a slight memory overhead to track
the information of which qubit was lost.
See [../loss.md] for details.
