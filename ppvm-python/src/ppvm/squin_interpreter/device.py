from dataclasses import dataclass, field
from typing import Any, ParamSpec, TypeVar, cast

from bloqade.analysis.address import AddressAnalysis, UnknownQubit, UnknownReg
from bloqade.device import AbstractSimulatorDevice, AbstractSimulatorTask
from kirin import ir

from ..generalized_tableau import GeneralizedTableau
from ._interp import GeneralizedTableauInterpreter

RetType = TypeVar("RetType")
Param = ParamSpec("Param")


@dataclass
class GeneralizedTableauSimulatorTask(AbstractSimulatorTask[Param, RetType, GeneralizedTableau]):
    generalized_tableau_interp: GeneralizedTableauInterpreter

    def run(self) -> RetType:
        _, ret = self.generalized_tableau_interp.run(
            self.kernel,
            *self.args,
            **self.kwargs,
        )
        return cast(RetType, ret)

    @property
    def state(self) -> GeneralizedTableau:
        return self.generalized_tableau_interp.backend


@dataclass
class GeneralizedTableauSimulator(AbstractSimulatorDevice[GeneralizedTableauSimulatorTask]):
    n_qubits: int | None = None
    options: dict[str, Any] = field(default_factory=dict)

    def task(
        self,
        kernel: ir.Method[Param, RetType],
        args: tuple[Any, ...] = (),
        kwargs: dict[str, Any] | None = None,
    ):
        address_analysis = AddressAnalysis(dialects=kernel.dialects)
        frame, _ = address_analysis.run(kernel)

        if self.n_qubits is None and any(
            isinstance(a, (UnknownQubit, UnknownReg)) for a in frame.entries.values()
        ):
            raise ValueError(
                "All addresses must be resolved. Or set n_qubits to a positive integer."
            )

        if self.n_qubits is not None and self.n_qubits <= 0:
            raise ValueError(f"n_qubits must be a positive integer, got {self.n_qubits}.")

        n_qubits = max(self.n_qubits or 0, address_analysis.qubit_count)

        tableau_options = dict(self.options)
        seed = tableau_options.pop("seed", None)
        tab = GeneralizedTableau(n_qubits=n_qubits, seed=seed, **tableau_options)
        interp = GeneralizedTableauInterpreter(kernel.dialects, backend=tab, rng_seed=seed)
        return GeneralizedTableauSimulatorTask(kernel, args, kwargs or {}, interp)
