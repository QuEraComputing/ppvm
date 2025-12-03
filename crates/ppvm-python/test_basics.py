from src import pauli_sum

state = pauli_sum(2, terms=["ZZ"], coefficients=(1.0,))

print(state)

state.cnot(0, 1)
state.h(0)

print(state)

print(state.overlap_with_zero())
print(state.trace("Z?*"))


n = 200
weight = 80

terms = ["".join(["Z" if i == j else "I" for i in range(n)]) for j in range(n)]
large_state = pauli_sum(n, max_pauli_weight=weight, terms=terms)


for i in reversed(range(1, n)):
    large_state.cnot(i - 1, i)

large_state.h(0)

print(large_state.overlap_with_zero())
