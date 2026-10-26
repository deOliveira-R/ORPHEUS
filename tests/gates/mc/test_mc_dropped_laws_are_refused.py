"""MC refuses a declaration it would drop (#654) and a body it would fill (#513).

MC applies the periodic wrap to every face of a slab, so any other face (a wall, a mirror) would be dropped. Since wave 2 (slice
S4) the refusal is the admission table's: a slab admits the translation deck only, everything else ``Refused(PLUMBING, #654)``.
The cavity of a hollow cylinder stays a scope guard (``NotImplementedError``): MC's lookup would fill it with the innermost
material (the specification's S3.14). The geometry's half of the old rows (a translation deck needs its partner, ``DeckError``)
is the geometry's pairing rule, pinned by
``tests/gates/geometry/test_face_roles.py::TestACurvedFaceCarriesNoDeck::test_the_slab_admits_both_decks_and_closes_the_pairing`` (a mesh duplicate was deleted on 2026-10-26).

Dispositions of the rows of the pre-wave-2 file, each with its reason (``retirement-audit`` D.14 and D.16):

* the MC rows "a non-periodic left law is refused" on a slab with a periodic right face: a periodic face needs a periodic partner (the
  geometry refuses a lone ``TranslationDeck`` with ``DeckError``), so the subject splits into the geometry's pairing rule (pinned in geometry's ``test_face_roles.py``) and the table's refusal of every non-translation face (the rows here).

Moved from ``tests/gates/mesh/test_dropped_laws_are_refused.py`` by the layer-DAG re-filing (Phase 1b, 2026-10-26): the
refusal is MC's, so the row lives with MC's gates; the bodies are unchanged. The two mesh builders the three methods share
live in ``tests/gates/transport/_dropped_laws_helpers.py``.
"""
from __future__ import annotations

import pytest

from orpheus.geometry import CoordSystem, MirrorDeck, TranslationDeck
from orpheus.mc.solver import MCMesh
from orpheus.transport.boundary_conditions import BoundaryRefused, RefusalAxis, WallCondition
from tests.gates.transport._dropped_laws_helpers import _hollow, _slab

pytestmark = pytest.mark.foundation


# ── S3.14: MC reads the translation of a unit cell ──────────────────────────────────────────────────────────────────


class TestMC:
    @pytest.mark.parametrize(
        "faces", [(WallCondition.vacuum(), WallCondition.vacuum()), (MirrorDeck(), MirrorDeck())], ids=["vacuum", "mirror"],
    )
    def test_a_slab_whose_faces_are_not_the_translation_is_refused(self, faces):
        """MC applies the periodic wrap to every face: any other face (a wall, a mirror) is dropped. After the carve the table admits
        the translation deck on a slab and states everything else as ``Refused(PLUMBING, #654)`` (the registry held nothing else)."""
        mesh, conditions = _slab(*faces)
        with pytest.raises(BoundaryRefused) as info:
            MCMesh(mesh, pitch=2.0, boundary_conditions=conditions)
        assert info.value.admission.axis is RefusalAxis.PLUMBING and info.value.admission.issue == 654

    def test_a_translation_pair_builds(self):
        mesh, conditions = _slab(TranslationDeck(), TranslationDeck())
        MCMesh(mesh, pitch=2.0, boundary_conditions=conditions)

    def test_a_hollow_cylinder_is_refused(self):
        """Its lookup would fill the cavity with the innermost material (the scope guard, #513, is read before the declaration)."""
        mesh = _hollow(CoordSystem.CYLINDRICAL)
        with pytest.raises(NotImplementedError, match="cavity"):
            MCMesh(mesh, pitch=4.0, boundary_conditions={1: WallCondition.specular(1.0), 0: WallCondition.specular(1.0)})
