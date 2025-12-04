import math
from typing import Sequence

import ppvm_python


def pauli_sum(
    n_qubits: int,
    min_abs_coeff: float = 1e-10,
    max_pauli_weight: int | None = None,
    terms: Sequence[str] = (),
    coefficients: Sequence[float] = (),
):
    # number of bytes we need
    N = math.ceil(n_qubits / 8.0)

    # number of bytes we have
    possible_interfaces = range(15)
    N_interface = next(n for n in possible_interfaces if 2**n > N)

    interface = getattr(ppvm_python, f"PauliSumIndexMapFxHash{N_interface}")

    if terms and not coefficients:
        coefficients = (1.0,) * len(terms)

    if max_pauli_weight is None:
        # NOTE: let rust handle the default setting for max_pauli_weight
        return interface(
            n_qubits,
            min_abs_coeff=min_abs_coeff,
            terms=terms,
            coefficients=coefficients,
        )
    else:
        return interface(
            n_qubits,
            min_abs_coeff=min_abs_coeff,
            max_pauli_weight=max_pauli_weight,
            terms=terms,
            coefficients=coefficients,
        )
