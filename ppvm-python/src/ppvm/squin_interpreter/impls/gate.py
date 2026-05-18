# SPDX-FileCopyrightText: 2026 QuEra Computing Inc.
# SPDX-License-Identifier: Apache-2.0

import math
from typing import Any

from bloqade.squin import gate
from kirin import interp
from kirin.dialects import ilist

from .._interp import GeneralizedTableauInterpreter
from ..qubit import GeneralizedTableauQubit


def _turns_to_radian(turns: float):
    return 2 * turns * math.pi


_SQRT_MAP = {
    "sqrtx": "sqrt_x",
    "sqrty": "sqrt_y",
}


@gate.dialect.register(key="generalized_tableau")
class GateMethods(interp.MethodTable):
    @interp.impl(gate.stmts.X)
    @interp.impl(gate.stmts.Y)
    @interp.impl(gate.stmts.Z)
    @interp.impl(gate.stmts.H)
    def single_qubit_gate(
        self,
        interp_: GeneralizedTableauInterpreter,
        frame: interp.Frame,
        stmt: gate.stmts.X | gate.stmts.Y | gate.stmts.Z | gate.stmts.H,
    ):
        qubits: ilist.IList[GeneralizedTableauQubit, Any] = frame.get(stmt.qubits)
        method = getattr(interp_.backend, stmt.name.lower())
        for qbit in qubits:
            method(qbit.addr)

    @interp.impl(gate.stmts.S)
    @interp.impl(gate.stmts.T)
    @interp.impl(gate.stmts.SqrtX)
    @interp.impl(gate.stmts.SqrtY)
    def single_qubit_nh_gate(
        self,
        interp_: GeneralizedTableauInterpreter,
        frame: interp.Frame,
        stmt: gate.stmts.S | gate.stmts.T | gate.stmts.SqrtX | gate.stmts.SqrtY,
    ):
        method_name = stmt.name.lower()
        method_name = _SQRT_MAP.get(method_name, method_name)

        if stmt.adjoint:
            method_name += "_adj"

        method = getattr(interp_.backend, method_name)
        qubits: ilist.IList[GeneralizedTableauQubit, Any] = frame.get(stmt.qubits)
        for qbit in qubits:
            method(qbit.addr)

    @interp.impl(gate.stmts.Rx)
    @interp.impl(gate.stmts.Ry)
    @interp.impl(gate.stmts.Rz)
    def rot(
        self,
        interp_: GeneralizedTableauInterpreter,
        frame: interp.Frame,
        stmt: gate.stmts.Rx | gate.stmts.Ry | gate.stmts.Rz,
    ):
        method = getattr(interp_.backend, stmt.name.lower())

        qubits: ilist.IList[GeneralizedTableauQubit, Any] = frame.get(stmt.qubits)
        angle = _turns_to_radian(frame.get(stmt.angle))

        for qbit in qubits:
            method(qbit.addr, angle)

    @interp.impl(gate.stmts.CX)
    @interp.impl(gate.stmts.CY)
    @interp.impl(gate.stmts.CZ)
    def control(
        self,
        interp_: GeneralizedTableauInterpreter,
        frame: interp.Frame,
        stmt: gate.stmts.CX | gate.stmts.CY | gate.stmts.CZ,
    ):
        method_name = stmt.name.lower()
        if method_name == "cx":
            method_name = "cnot"

        method = getattr(interp_.backend, method_name)

        controls: ilist.IList[GeneralizedTableauQubit, Any] = frame.get(stmt.controls)
        targets: ilist.IList[GeneralizedTableauQubit, Any] = frame.get(stmt.targets)

        for ctrl, target in zip(controls, targets):
            method(ctrl.addr, target.addr)

    @interp.impl(gate.stmts.U3)
    def u3(
        self,
        interp_: GeneralizedTableauInterpreter,
        frame: interp.Frame,
        stmt: gate.stmts.U3,
    ):
        qubits: ilist.IList[GeneralizedTableauQubit, Any] = frame.get(stmt.qubits)
        theta = _turns_to_radian(frame.get(stmt.theta))
        phi = _turns_to_radian(frame.get(stmt.phi))
        lam = _turns_to_radian(frame.get(stmt.lam))

        for qbit in qubits:
            interp_.backend.u3(qbit.addr, theta, phi, lam)

    @interp.impl(gate.stmts.PhasedXZ)
    def phased_xz(
        self,
        interp_: GeneralizedTableauInterpreter,
        frame: interp.Frame,
        stmt: gate.stmts.PhasedXZ,
    ):
        x_exponent = frame.get(stmt.x_exponent)
        z_exponent = frame.get(stmt.z_exponent)
        axis_phase_exponent = frame.get(stmt.axis_phase_exponent)
        qubits: ilist.IList[GeneralizedTableauQubit, Any] = frame.get(stmt.qubits)

        angle_rz_pre = -axis_phase_exponent * math.pi * 2
        angle_rx = x_exponent * math.pi * 2
        angle_rz_post = (axis_phase_exponent + z_exponent) * math.pi * 2

        for qbit in qubits:
            interp_.backend.rz(qbit.addr, angle_rz_pre)
            interp_.backend.rx(qbit.addr, angle_rx)
            interp_.backend.rz(qbit.addr, angle_rz_post)
