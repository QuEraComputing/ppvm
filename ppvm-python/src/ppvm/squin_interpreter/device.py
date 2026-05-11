from typing import Any, TypeVar, ParamSpec, cast
from dataclasses import field, dataclass

from ppvm import GeneralizedTableau
from kirin import ir
from bloqade.device import AbstractSimulatorTask, AbstractSimulatorDevice
from bloqade.analysis.address import UnknownReg, UnknownQubit, AddressAnalysis

from ._interp import GeneralizedTableauInterpreter

RetType = TypeVar("RetType")
Param = ParamSpec("Param")
InternalState = TypeVar("InternalState")


@dataclass
class GeneralizedTableauSimulatorTask(
    AbstractSimulatorTask[Param, RetType, GeneralizedTableau]
):
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
class GeneralizedTableauSimulator(
    AbstractSimulatorDevice[GeneralizedTableauSimulatorTask]
):
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

        n_qubits = max(self.n_qubits or 0, address_analysis.qubit_count)

        tab = GeneralizedTableau(n_qubits=n_qubits, **self.options)
        interp = GeneralizedTableauInterpreter(
            kernel.dialects, backend=tab, rng_seed=self.options.get("seed")
        )
        return GeneralizedTableauSimulatorTask(kernel, args, kwargs or {}, interp)
