# Simulating loss with Pauli Propagation


## Implementation

To simulate loss, `ppvm` offers a [`LossyPauliSum`][ppvm.paulisum.LossyPauliSum] class.
This is a dedicated class, which behaves just like a [`PauliSum`][ppvm.paulisum.PauliSum], but adds additional methods for the loss.

Also, this separation is there since we need to extend the Pauli basis in order to
account for loss (see [below](#background)). This comes at a storage overhead. Specifically, we now need 3 bits in order to represent a character in a Pauli string rather than just 2.

Here is a small example:

```python
from ppvm import LossyPauliSum

ps = LossyPauliSum.new(n_qubits = 1, terms=["Z"])

# Reset at the end of the circuit
ps.reset_loss_channel(0)

# Loss after an X gate
ps.loss_channel(0, 0.1)

# Apply an X gate
ps.x(0)

z_exp = ps.overlap_with_zero()

# This will be -0.8: in 10% of cases we have <Z> = 1 instead of -1.
print(f"<Z>: {z_exp}")
```


## Background

### Updated Pauli basis
We can include loss in qubit simulation by adding a third state,
which we will call the leakage state $|L\rangle$.
In order to include this in Pauli Propagation, we extend the
set of Pauli basis operator by an additional operator $L = |L\rangle \langle L|$,
which is the projector on the leakage state.

The complete set of Pauli operators we use to describe any Pauli String is then

$$\{I, X, Y, Z, L\}.$$

Neglecting coherences between the qubit subspace and the leakage state, this basis
fully describes the three-level system.

For clarity, the corresponding matrix definitions are

$$I = \begin{pmatrix} 1 & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & 0
\end{pmatrix}, ~~ X = \begin{pmatrix} 0 & 1 & 0 \\
1 & 0 & 0 \\
0 & 0 & 0
\end{pmatrix}, ~~Y =\begin{pmatrix} 0 & -i & 0 \\
i & 0 & 0 \\
0 & 0 & 0
\end{pmatrix} , \\
Z = \begin{pmatrix} 1 & 0 & 0 \\
0 & -1 & 0 \\
0 & 0 & 0
\end{pmatrix}, ~~ L = \begin{pmatrix} 0 & 0 & 0 \\
0 & 0 & 0 \\
0 & 0 & 1
\end{pmatrix}.
$$

!!! note
    In the basis above, gates **no longer correspond to the matrix** definitions.
    The reason for this is simple: applying e.g. an $X$ gate to a qubit in $|L\rangle$,
    will leave the qubit invariant. However, it is obvious that multiplying any Pauli
    $P \neq L$ from the basis with $L$ will give 0. This would correspond to a
    zero-amplitude state. The implementation takes that into account accordingly.


### Loss channel

The action of an independent loss channel is to map part of the population from
the qubit subspace into the $|L\rangle$ state.
This corresponds to reducing the trace of the density operator.

The Kraus operators for a loss channel are

$$K_0 = \sqrt{1 - p_L}I + |L\rangle\langle L|, ~~ K_1 = \sqrt{p_L} |L\rangle\langle 0|, ~~ K_2 = \sqrt{p_L} |L\rangle\langle 1|,$$

where $p_L$ is the probability of losing a qubit.

The action of the channel on the Pauli basis is

$$\mathcal{E}(P) = (1 - p_L) P, ~~ P\in\{I, X, Y, Z\},\\
\mathcal{E}(L) = L + p_L I.
$$

Note that this channel is not symmetric under forward and backward propagation.
Since for any channel $\mathcal{E}[\rho] = \sum_i K_i \rho K_i ^\dagger$, we have $\langle A \rangle = \text{tr}(A \mathcal{E}[\rho]) = \text{tr}(\sum_i K_i^\dagger A K_i \rho)$,
the above action on Paulis corresponds to the **adjoint** channel.

The action of the channel is clear from a physical perspective:
whenever a qubit undergoes a loss channel, we lose population out of the qubit subspace.

However, we can also see if there is no loss present in the initial state, we will never map any part of a Pauli string to $L$. This might seem a bit surprising, but is accurate. Under the loss channel described here, we simply lose population out of the qubit subspace. Any expectation value of a qubit-subspace operator will tend towards $0$ with increasing loss.

An intuitive picture is also to consider a pure state which describes a single probabilistic trajectory through the circuit. Before loss, a state vector in the qubit subspace is just

$$|\psi\rangle = \begin{pmatrix} c_0 \\ c_1 \\ 0 \end{pmatrix}. $$

Should the qubit be lost, the resulting state vector is

$$ |\psi\rangle = |L\rangle.$$

Since under the described dynamics there is no way to recover a qubit once it has been lost, any trajectory where the qubit has been lost will contribute a 0 overlap for any qubit subspace operator $Q$,

$$\langle \psi| Q|\psi\rangle = 0.$$

While this is conceptually accurate, the resulting loss does not capture what happens on hardware. If there is no loss detection, then a qubit that has been lost will be falsely counted as being in the $|0\rangle$ state when being measured. If there is loss detection, we could either reject the trajectory via post-selection, which is equivalent to having no loss in the system. Or, we can reset the qubit to the $|0\rangle$ during the circuit.


### Reset channel

In order to accurately model the hardware behavior, we will extend the above description by another channel. A reset channel, which incoherently resets a lost qubit into the $|0\rangle$ state.

The underlying dynamics of this channel are equivalent to having amplitude damping from $|L\rangle$ to the $|0\rangle$ state. The Kraus operators for this channel are

$$ K_0 = I, ~~~ K_1 = |0\rangle\langle L|.$$
The corresponding mappings of the Pauli basis are

$$\mathcal{E}(P) = P, ~ P\in\{X, Y\},$$

$$\mathcal{E}(P) = P + L, ~ P\in \{I, Z\},$$

$$\mathcal{E}(L) = 0.$$
Intuitively, these mappings can be understood as
- Coherences remain unchanged.
- Lost qubits are reset: lost population is removed ($L$ is set to 0) and added to the $|0\rangle$ state, which is equivalent to a positive contribution to both $Z$ or $I$.

This reset channel is sufficient in order to model either actually resetting a lost qubit, or the error that arises when falsely counting a lost qubit as 0. For the latter, we simply apply the loss channel to all qubits at the end of the circuit. This ensures that lost qubits to not partake in the dynamics of the circuit, but are then counted as 0 in a measurement.

Note, that since we are doing Pauli Propagation, we always have to apply the adjoint circuit. Thus, resetting lost qubits at the end of the circuit means that we apply the reset channel at the very beginning when propagating a Pauli string.


### Comment on branching and scaling

Both channels above branch, which without truncation would lead to exponential scaling.

For the loss channel, it is easy to see that coefficient truncation can deal with branching,
since it only branches on $L \to p_L I + L$.
Branches therefore scale with $p_L$, which is usually $\ll 1$.

The reset channel, on the other hand, branches on both $I$ and $Z$ into $L$ without
any leading coefficient.
The probability of having lost the $i$-th qubit at the end of the circuit is just
given by $\langle L_i\rangle$.
Since there are no coherences between $L$ and the qubit subspace, expectation values
factorize, i.e. for any Pauli string, where $P_i \neq L_i$,

$$\langle P_1 L_2 P_3 P_4 L_5 P_6 ...\rangle = \langle L_2\rangle\langle L_5\rangle \langle P_1P_3 P_4 P_6...\rangle.$$

In a slight abuse of notation, we note that $\langle L_i\rangle \propto p_L $,
the contribution of Pauli Strings that have an $L$ on many positions is very small.
Therefore, we can truncate these in order to suppress exponential branching.

In practice, this works similar to large weight truncation and is achieved by a designated
truncation strategy, which you can configure by setting the [`max_loss_weight`] on a [`ppvm.paulisum.LossyPauliSum`].
