r"""The harmonic frame's mints refuse a MOMENT space — by name, on both doors
(CS4c step 6 item 6.2c-ii; the memo's gate P6, hazard H-6).

Until 6.2c-ii a moment space was axes-less, so ``moment_space_on`` refused it
through its *"must be axis-built"* check — a guard that silently lost its
subject the moment the head became axis-built (`[M]` the 6.2c verification
round: an axis-built moment space was ACCEPTED and returned a plausible-looking
product). Both doors now refuse through the same channel the blessed frame
chain reads — the leading axis's generator must be the QUADRATURE — so the
refusal's TYPE and MESSAGE are the ones asserted here, on the shipped moment
space of every ledger carrier.

Layer DAG Phase 4 (2026-11-01): posed on the four monomorphic material meshes and the raw constructors. The leg
"the frame's mint on the angular space equals the problem's moment space" left with the problem: here both sides
would be the same call (``retirement-audit`` D.16); the hub-versus-frame agreement is
``tests/gates/sn/mesh/test_hub_and_frame_agree_on_the_moment_space.py::test_P5_…``'s, and P5 alone is the dropped
leg's catcher: ``[M]`` 2026-11-01, the mutation "the frame builds its head undressed" reds P5 (12) and P7 (12) and none
of this file's rows, which are refusal rows outside that mutation's reach (battery ``frame_moment_space_drops_parseval``,
``scratch/layer_dag/phase4/landing_report.md``).

Foundation mark: software invariant; no physics claim.
"""

from __future__ import annotations

import pytest

from orpheus.transport.frames import HarmonicFrame
from orpheus.transport.phase_spaces import angular_bulk_space
from tests.gates.transport import _phase_space_cases as cases

pytestmark = pytest.mark.foundation

_GEOMETRIES = ("slab", "sphere", "cylinder", "cart2d")


@pytest.mark.parametrize("geometry", _GEOMETRIES, ids=list(_GEOMETRIES))
def test_P6_both_doors_refuse_a_moment_space_naming_the_missing_quadrature(geometry: str) -> None:
    mesh, quadrature = cases.monomorphic()[geometry]
    L = 1
    bulk = angular_bulk_space(mesh, quadrature)
    frame = HarmonicFrame.for_space(bulk, L)
    moment = frame.moment_space_on(bulk)
    assert moment.axes is not None, "the moment space IS axis-built — the old axes-less refusal cannot catch it"
    with pytest.raises(ValueError, match="HarmonicFrame.moment_space_on needs the generating Quadrature"):
        frame.moment_space_on(moment)
    with pytest.raises(ValueError, match="HarmonicFrame.for_space needs the generating Quadrature"):
        HarmonicFrame.for_space(moment, L)
    # the positive leg: the angular space passes the frame's door, and the frame is interned per (space, L)
    assert HarmonicFrame.for_space(bulk, L) is frame
