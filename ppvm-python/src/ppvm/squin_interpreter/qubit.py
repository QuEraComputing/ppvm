from typing import Any
from dataclasses import dataclass

from ppvm import MeasurementResult
from kirin import interp
from bloqade.types import Qubit
from kirin.dialects import ilist
from bloqade.decoders.dialects.annotate.types import MeasurementResultValue

from bloqade import qubit

from ._interp import GeneralizedTableauInterpreter


@dataclass
class GeneralizedTableauQubit(Qubit):
    addr: int


def _measurement_result_conversion(result: MeasurementResult):
    if result == MeasurementResult.ZERO:
        return MeasurementResultValue.Zero
    elif result == MeasurementResult.ONE:
        return MeasurementResultValue.One
    return MeasurementResultValue.Lost


@qubit.dialect.register(key="generalized_tableau")
class QubitMethods(interp.MethodTable):

    @interp.impl(qubit.stmts.New)
    def new(
        self,
        interp_: GeneralizedTableauInterpreter,
        frame: interp.Frame,
        stmt: qubit.stmts.New,
    ):
        addr = interp_.allocate_qubit()
        return (GeneralizedTableauQubit(addr),)

    @interp.impl(qubit.stmts.Measure)
    def measure(
        self,
        interp_: GeneralizedTableauInterpreter,
        frame: interp.Frame,
        stmt: qubit.stmts.Measure,
    ):
        qubits: ilist.IList[GeneralizedTableauQubit, Any] = frame.get(stmt.qubits)
        results = [interp_.backend.measure(qbit.addr) for qbit in qubits]
        results_converted = ilist.IList(
            list(map(_measurement_result_conversion, results))
        )
        return (results_converted,)

    @interp.impl(qubit.stmts.Reset)
    def reset(
        self,
        interp_: GeneralizedTableauInterpreter,
        frame: interp.Frame,
        stmt: qubit.stmts.Reset,
    ):
        qubits: ilist.IList[GeneralizedTableauQubit, Any] = frame.get(stmt.qubits)
        for qbit in qubits:
            interp_.backend.reset(qbit.addr)
