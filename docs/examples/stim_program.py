"""
Stim interop example.

Parse a Stim program, execute it against a generalized tableau, and sample
it many times via ``sample_stim``.
"""

from ppvm import GeneralizedTableau, StimProgram, sample_stim

prog = StimProgram.parse(
    """
    H 0
    CX 0 1
    M 0 1
    """
)

# Run once against a fresh tableau.
tab = GeneralizedTableau(n_qubits=2)
tab.run(prog)
print("single-run ok")  # → single-run ok

# Multi-shot sampling. The two measurement outcomes are correlated on every shot.
shots = sample_stim(prog, n_qubits=2, num_shots=16, seed=0)
print(f"sampled {len(shots)} shots")  # → sampled 16 shots
print("first shot:", shots[0])  # → first shot: [<MeasurementResult.ONE: 1>, <MeasurementResult.ONE: 1>]
all_correlated = all(s[0] == s[1] for s in shots)
print("all shots correlated:", all_correlated)  # → all shots correlated: True
