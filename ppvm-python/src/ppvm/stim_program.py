"""StimProgram: parsed + normalized Stim circuit, ready for many shots."""

from __future__ import annotations

import ppvm_python_native


class StimProgram:
    """A parsed and normalized Stim circuit.

    Use ``parse`` to construct from a source string, or ``from_file`` to read
    a ``.stim`` file. The resulting object is a thin handle around a
    pre-normalized program; passing it to ``tab.run(prog)`` or
    ``GeneralizedTableau.sample(prog, ...)`` reuses the normalized form across
    shots.
    """

    __slots__ = ("_inner",)

    def __init__(self, _inner: object) -> None:  # internal
        self._inner = _inner

    @staticmethod
    def parse(src: str) -> "StimProgram":
        return StimProgram(ppvm_python_native.StimProgram.parse(src))

    @staticmethod
    def from_file(path: str) -> "StimProgram":
        return StimProgram(ppvm_python_native.StimProgram.from_file(path))

    def __repr__(self) -> str:
        return repr(self._inner)
