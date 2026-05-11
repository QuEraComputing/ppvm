from dataclasses import field, dataclass

from ..generalized_tableau import GeneralizedTableau
from kirin import interp
from typing_extensions import Self


@dataclass
class GeneralizedTableauInterpreter(interp.Interpreter):
    keys = ("generalized_tableau", "main")
    backend: GeneralizedTableau

    rng_seed: int | None = None

    current_qubit_addr: int = field(init=False, default=0)

    def initialize(self) -> Self:
        super().initialize()
        self.backend = GeneralizedTableau(
            n_qubits=self.backend.n_qubits,
            min_abs_coeff=self.backend.min_abs_coeff,
            seed=self.rng_seed,
        )
        self.current_qubit_addr = 0
        return self

    def allocate_qubit(self) -> int:
        addr = self.current_qubit_addr
        self.current_qubit_addr += 1
        return addr
