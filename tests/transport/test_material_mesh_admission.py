r"""Carrier-admission gates: the ``mesh is None`` two-meanings un-weld (S7 G7.1–G7.3).

The three CS4b S7 repairs, gated. `[M]` re-measured at HEAD 2026-08-24
before the repair (matching the census probes exactly): the SN promotion
of the mesh-less carrier died on a MESSAGELESS ``assert`` under plain
python and on a deep ``AttributeError: 'NoneType' object has no
attribute 'coord'`` under the CANONICAL ``-O`` runner (the assert
stripped); ``.areas`` answered every armless carrier with a message
claiming "2-D meshes" — false on two of its three arms; and the
``mesh is None`` sentinel carried two DIFFERENT meanings (the
infinite-medium 1-cell carrier vs the d≥3 axis-native carrier) that
nothing discriminated.

G7.1 — the promotion refusal is a REAL raise (this module runs under the
canonical ``-O``, where the pre-repair assert did not run at all —
Mode-8: a gate against ``AssertionError`` would have been false).
G7.2 — one message per arm, each pinned by its shortest distinctive
fragment; collapsing the three back to one reddens two rows.
G7.3 — the un-weld LAW: the two ``mesh is None`` states are
distinguishable (``ndim``) and take different arms everywhere it
matters; merging the branches reddens the discrimination rows.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.transport.mesh.axis import AxisMesh
from orpheus.transport.mesh.material_mesh import MaterialMesh
from tests.sn._test_helpers import placeholder_materials

pytestmark = [pytest.mark.foundation]


def _meshless() -> MaterialMesh:
    return MaterialMesh.from_materials(placeholder_materials(ng=2))


def _legacy_1d() -> MaterialMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, 6), mat_ids=np.zeros(5, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"), bc_right=BC("vacuum"),
    )
    return MaterialMesh(mesh, placeholder_materials(ng=2))


def _d3_sn() -> SNMesh:
    axes = tuple(
        AxisMesh(edges=np.linspace(0.0, ext, n + 1))
        for ext, n in zip((1.0, 2.0, 3.0), (2, 3, 2))
    )
    return SNMesh.from_axes(
        axes, Quadrature.level_symmetric(sn_order=4),
        placeholder_materials(ng=2),
    )


class TestG71PromotionRefusal:
    def test_a_real_carrier_promotes(self):
        """Positive leg (vv #11): the legacy 1-D carrier promotes to a
        solvable SN phase space."""
        sn = SNMesh.from_material_mesh(
            _legacy_1d(), Quadrature.gauss_legendre(4),
        )
        if sn.ng != 2 or sn.spatial_shape != (5,):
            pytest.fail("promotion must preserve the carrier's data block")

    def test_the_meshless_carrier_is_refused_with_a_real_raise(self):
        """Negative leg, run under the canonical -O: a typed ValueError
        naming the case — NOT the pre-repair messageless assert (which
        -O stripped, leaving a deep AttributeError in the streaming
        constructor)."""
        with pytest.raises(ValueError, match="infinite-medium 1-cell"):
            SNMesh.from_material_mesh(
                _meshless(), Quadrature.gauss_legendre(4),
            )


class TestG72AreasNamesItsOwnCase:
    def test_meshless_arm_names_the_facelessness(self):
        with pytest.raises(AttributeError, match="no faces at all"):
            _meshless().areas

    def test_axis_native_arm_names_the_carrier(self):
        """The d=3 axis-native SN carrier (mesh is None, ndim=3) — the
        pre-repair message blamed '2-D meshes' here, falsely."""
        with pytest.raises(AttributeError, match="axis-native"):
            _d3_sn().areas

    def test_legacy_2d_arm_keeps_its_true_message(self):
        from orpheus.geometry import Mesh2D

        mesh2d = Mesh2D(
            edges_x=np.linspace(0.0, 1.0, 3),
            edges_y=np.linspace(0.0, 2.0, 4),
            mat_map=np.zeros((2, 3), dtype=int),
            bc_xmin=BC("vacuum"), bc_xmax=BC("vacuum"),
            bc_ymin=BC("vacuum"), bc_ymax=BC("vacuum"),
        )
        carrier = MaterialMesh(mesh2d, placeholder_materials(ng=2))
        with pytest.raises(AttributeError, match="not defined for 2-D"):
            carrier.areas


class TestG73TheSentinelIsNotTheMeaning:
    def test_the_two_mesh_none_states_are_discriminated(self):
        """The un-weld LAW: both carriers spell ``mesh is None``, and
        every consumer arm discriminates them by structure (ndim), not
        by the sentinel — different areas refusals, opposite promotion
        outcomes. Merging any branch pair reddens this."""
        meshless, d3 = _meshless(), _d3_sn()
        if meshless.mesh is not None or d3.mesh is not None:
            pytest.fail("both states must spell mesh is None")
        if (meshless.ndim, d3.ndim) != (1, 3):
            pytest.fail("ndim must discriminate the two states")
        msgs = []
        for carrier in (meshless, d3):
            with pytest.raises(AttributeError) as exc:
                carrier.areas
            msgs.append(str(exc.value))
        if msgs[0] == msgs[1]:
            pytest.fail("the two mesh-None areas arms must differ")
        # Promotion: the meshless state refuses; the d≥3 state IS a
        # promoted SN phase space already (from_axes succeeded above).
        with pytest.raises(ValueError):
            SNMesh.from_material_mesh(
                meshless, Quadrature.gauss_legendre(4),
            )
