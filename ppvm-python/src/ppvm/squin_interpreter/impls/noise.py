# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0

from typing import Any

from bloqade.squin import noise
from kirin import interp
from kirin.dialects import ilist

from .._interp import GeneralizedTableauInterpreter
from ..qubit import GeneralizedTableauQubit


@noise.dialect.register(key="generalized_tableau")
class NoiseMethods(interp.MethodTable):
    @interp.impl(noise.stmts.Depolarize)
    def depolarize(
        self,
        interp_: GeneralizedTableauInterpreter,
        frame: interp.Frame,
        stmt: noise.stmts.Depolarize,
    ):
        p = frame.get(stmt.p)
        qubits: ilist.IList[GeneralizedTableauQubit, Any] = frame.get(stmt.qubits)
        for qbit in qubits:
            interp_.backend.depolarize(qbit.addr, p)

    @interp.impl(noise.stmts.Depolarize2)
    def depolarize2(
        self,
        interp_: GeneralizedTableauInterpreter,
        frame: interp.Frame,
        stmt: noise.stmts.Depolarize2,
    ):
        p = frame.get(stmt.p)
        controls: ilist.IList[GeneralizedTableauQubit, Any] = frame.get(stmt.controls)
        targets: ilist.IList[GeneralizedTableauQubit, Any] = frame.get(stmt.targets)
        for control, target in zip(controls, targets):
            interp_.backend.depolarize2(control.addr, target.addr, p)

    @interp.impl(noise.stmts.SingleQubitPauliChannel)
    def single_qubit_pauli_channel(
        self,
        interp_: GeneralizedTableauInterpreter,
        frame: interp.Frame,
        stmt: noise.stmts.SingleQubitPauliChannel,
    ):
        px = frame.get(stmt.px)
        py = frame.get(stmt.py)
        pz = frame.get(stmt.pz)
        qubits: ilist.IList[GeneralizedTableauQubit, Any] = frame.get(stmt.qubits)

        ps = [px, py, pz]
        for qbit in qubits:
            interp_.backend.pauli_error(qbit.addr, ps)

    @interp.impl(noise.stmts.TwoQubitPauliChannel)
    def two_qubit_pauli_channel(
        self,
        interp_: GeneralizedTableauInterpreter,
        frame: interp.Frame,
        stmt: noise.stmts.TwoQubitPauliChannel,
    ):
        ps = frame.get(stmt.probabilities)
        controls: ilist.IList[GeneralizedTableauQubit, Any] = frame.get(stmt.controls)
        targets: ilist.IList[GeneralizedTableauQubit, Any] = frame.get(stmt.targets)

        for control, target in zip(controls, targets):
            interp_.backend.two_qubit_pauli_error(control.addr, target.addr, ps)

    @interp.impl(noise.stmts.QubitLoss)
    def qubit_loss(
        self,
        interp_: GeneralizedTableauInterpreter,
        frame: interp.Frame,
        stmt: noise.stmts.QubitLoss,
    ):
        p = frame.get(stmt.p)
        qubits: ilist.IList[GeneralizedTableauQubit, Any] = frame.get(stmt.qubits)
        for qbit in qubits:
            interp_.backend.loss_channel(qbit.addr, p)

    @interp.impl(noise.stmts.CorrelatedQubitLoss)
    def correlated_qubit_loss(
        self,
        interp_: GeneralizedTableauInterpreter,
        frame: interp.Frame,
        stmt: noise.stmts.CorrelatedQubitLoss,
    ):
        p = frame.get(stmt.p)
        ps = [p, 0, 0]
        qubits: list[list[GeneralizedTableauQubit]] = frame.get(stmt.qubits)
        for qubit_group in qubits:
            if len(qubit_group) != 2:
                raise ValueError("Correlated loss only supported for qubit pairs!")
            addr0, addr1 = [qbit.addr for qbit in qubit_group]
            interp_.backend.correlated_loss_channel(addr0, addr1, ps)
