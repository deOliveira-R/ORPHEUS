"""CP refuses a declaration it would drop (#513), as admission.

CP's slab kernel images a mirror at the left face and reads the outer face of a curved body, so before P1 step 2 a law it did not
read was dropped silently (``[M]`` 2026-09-25, the P1 specification's probes ``cp_slab_left_law.py`` and
``partial_law_readers.py``). Step 2 made a hollow body and a declared inner law expressible on a
:class:`~orpheus.geometry.StructuredGeometry`, so the refusal landed with it (the specification's S3.12). Since wave 2 (slice S4)
the refusal is the admission table's: the ``xmin`` face of a Cartesian site admits a ``MirrorDeck`` only, ``Refused(PLUMBING,
#513)``, and the inner face of a hollow body is refused at ``xmin`` the same way.

Dispositions of the rows of the pre-wave-2 file, each with its reason (``retirement-audit`` D.14 and D.16):

* ``wrong-axis-mirror-white`` (a perfect mirror about y on the slab's x-face): DEAD. A ``MirrorDeck`` has no axis, so the subject (a
  guard asking ``kind == "reflective"`` admitting the wrong motion) cannot be constructed. Its replacement witness is the face-local
  condition law (a condition is invariant under the carrier's colouring symmetry or ``bind`` refuses it:
  ``tests/gates/transport/test_binding.py``).
* ``test_the_slab_kernel_drops_the_left_law`` (the RECORD that CP's ``k`` does not depend on the declared left law, measured by bypassing
  ``_refuse_a_law_cp_drops``): DEAD. The guard is an admission cell of the table, and a bypass of an admission would test the
  monkeypatch and not CP. The measured numbers it recorded (``[M]`` 2026-09-30, two groups, fuel A on [0, 1] cm and moderator B on
  [1, 2.5] cm, 8 cells: ``k = 1.0402206390766764`` with a white right face and ``1.2301755431553738`` with a vacuum one, to the last bit,
  for every left law) are kept here. The vacuum figure is the ``k`` that ``solve_cp`` reports, not the eigenvalue of the CP pencil it
  iterates, which is 0.5016: ``compute_keff`` omits the leakage term (ERR-110, #546). The claim the row witnessed (a left face the
  kernel replaces by a mirror is refused) is ``tests/gates/cp/test_cp_realizer_table.py::test_cp_left_face_is_a_forced_mirror`` and the
  CP rows of the admission pin.

Moved from ``tests/gates/mesh/test_dropped_laws_are_refused.py`` by the layer-DAG re-filing (Phase 1b, 2026-10-26): the
refusal is CP's, so the row lives with CP's gates; the bodies are unchanged. The two mesh builders the three methods share
live in ``tests/gates/transport/_dropped_laws_helpers.py``.
"""
from __future__ import annotations

import pytest

from orpheus.cp.solver import CPMesh
from orpheus.geometry import CoordSystem, MirrorDeck
from orpheus.transport.boundary_conditions import BoundaryRefused, RefusalAxis, WallCondition
from tests.gates.transport._dropped_laws_helpers import _hollow, _slab

pytestmark = pytest.mark.foundation


# ── S3.12: CP reads only the outer face, and a mirror at a slab's left ──────────────────────────────────────────────


class TestCP:
    @pytest.mark.parametrize(
        "left, right",
        [
            (WallCondition.vacuum(), WallCondition.white()), (WallCondition.white(), WallCondition.vacuum()),
            (WallCondition.white(), WallCondition.white()), (WallCondition.specular(0.7), WallCondition.white()),
        ],
        ids=["vacuum-white", "white-vacuum", "white-white", "partial-specular-white"],
    )
    def test_a_slab_whose_left_face_is_not_the_mirror_is_refused(self, left, right):
        """CP's slab kernel images a mirror at the left face, so any other left condition, even one equal to the right one, would be
        replaced. The ``white-white`` row was admitted until 2026-09-30, on the false claim that equal laws are what CP computes
        (``[M]``: CP declared white|white 1.212883, SN reflective|white 1.212884, SN white|white 1.212537); here it is ONE tag bound to
        ``white`` on both faces (rule R1), refused at ``xmin``, the first face.

        The ``partial-specular-white`` row is the specular wall ``specular(0.7)``: it permutes ordinates like the mirror and returns 0.7
        of the outflow, and CP would compute a full mirror in its place; the row discriminates a table asking "does it permute
        ordinates" from one asking for the mirror. The refusal is the table's: ``BoundaryRefused`` whose admission is
        ``Refused(PLUMBING, #513)`` at the face ``xmin`` (each row is refused for its own face and issue, not only for its type)."""
        mesh, conditions = _slab(left, right)
        with pytest.raises(BoundaryRefused) as info:
            CPMesh(mesh, boundary_conditions=conditions)
        assert info.value.face == "xmin"
        assert info.value.admission.axis is RefusalAxis.PLUMBING and info.value.admission.issue == 513

    @pytest.mark.parametrize("right", [WallCondition.white(), WallCondition.vacuum()], ids=["white", "vacuum"])
    def test_a_slab_with_the_mirror_on_its_left_builds(self, right):
        mesh, conditions = _slab(MirrorDeck(), right)
        CPMesh(mesh, boundary_conditions=conditions)

    @pytest.mark.parametrize(
        "coord", [CoordSystem.CYLINDRICAL, CoordSystem.SPHERICAL], ids=lambda c: c.name.lower(),
    )
    def test_a_hollow_body_with_an_inner_condition_is_refused(self, coord):
        """The inner face of a hollow body carries a tag whose condition CP would drop: ``BoundaryRefused``, ``Refused(PLUMBING, #513)``
        at ``xmin`` (today's refusal, re-keyed by face; the same cell is pinned by ``test_cp_realizer_table``)."""
        with pytest.raises(BoundaryRefused) as info:
            CPMesh(_hollow(coord), boundary_conditions={1: WallCondition.white(), 0: WallCondition.white()})
        assert info.value.face == "xmin" and info.value.admission.issue == 513
