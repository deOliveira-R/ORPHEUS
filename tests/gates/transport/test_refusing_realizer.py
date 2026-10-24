"""The refusal-only realizer base: one unreachable ``realize`` body for every table, and the invariant that makes it unreachable.

A table registers a unit it never builds when it has a reason sharper than the unregistered door's (MC's torus refusal on a cylinder,
``Refused(TOPOLOGY, #547)``; a per-method plumbing message). Such a realizer's ``realize`` can never run, because the table's
``realize`` calls ``require`` first. :class:`~orpheus.transport.boundary_conditions.RefusingRealizer` holds that body ONCE (the S4
integration retired six copies in ``cp``, ``moc`` and ``mc``) and refuses, at class creation, a registered subclass whose status is
not a refusal. Mutation, by reading: drop the check in ``RefusingRealizer.__init_subclass__`` and the first test goes green on the
broken class, which reddens it.
"""
from __future__ import annotations

import pytest

from orpheus.transport.boundary_conditions import (
    FaceSite,
    Realized,
    RealizerTable,
    RefusalAxis,
    Refused,
    RefusingRealizer,
    SpecularChain,
    StateCone,
    WallCondition,
    ZeroNetwork,
)
from orpheus.geometry import CoordSystem

pytestmark = pytest.mark.foundation

_SITE = FaceSite("xmax", CoordSystem.CARTESIAN, 1)


def test_a_registered_refusing_realizer_must_state_a_refusal() -> None:
    table: RealizerTable = RealizerTable("probe", state_cone=StateCone.NONNEGATIVE_FLUX, unregistered_issue=654)
    with pytest.raises(TypeError, match="states a Refused status"):

        class Admitting(RefusingRealizer, table=table, realizes=ZeroNetwork):  # noqa: F841 (the class body is the probe)
            status = Realized(True)


def test_the_table_refuses_before_the_body_and_the_body_raises() -> None:
    table: RealizerTable = RealizerTable("probe", state_cone=StateCone.NONNEGATIVE_FLUX, unregistered_issue=654)

    class Refusing(RefusingRealizer, table=table, realizes=SpecularChain):
        status = Refused(RefusalAxis.PLUMBING, "not built", 654)

    table.seal()
    admission = table.admit(WallCondition.specular(0.5), _SITE)
    assert isinstance(admission, Refused) and admission.issue == 654
    with pytest.raises(AssertionError, match="never realized"):
        Refusing().realize(SpecularChain(), 0.5, _SITE, None)
