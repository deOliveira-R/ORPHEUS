"""The arms of the ``Evidence`` union, one instance each, and a kind OUTSIDE the union (not collected).

Both exhaustiveness gates of the evidence readers (``test_convergence_evidence_reader_is_exhaustive``
and ``tests/gates/sn/operators/test_gauge_evidence_reader_is_exhaustive``) import ``ARM_FIXTURES`` from
here, so that a fifth arm added to the union is given its fixture once.
"""

from __future__ import annotations

from dataclasses import dataclass

from orpheus.numerics.evidence import Asserted, Measured, NotApplicable, NotYet

__all__ = ["ARM_FIXTURES", "SILENT_ARMS", "NewKind"]


@dataclass(frozen=True)
class NewKind:
    """A value outside the union: no reader of ``Evidence`` may treat it as one of the four arms."""


ARM_FIXTURES = {
    Measured: Measured(0.5),
    Asserted: Asserted(1e-9, "t"),
    NotApplicable: NotApplicable("x"),
    NotYet: NotYet(1, "r"),
}

SILENT_ARMS = [arm for arm in ARM_FIXTURES if arm is not Measured]
