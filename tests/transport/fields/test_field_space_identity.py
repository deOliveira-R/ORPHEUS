r"""S2 gates for the field layer's space identity (campaign 1 CS4b).

G2.1–G2.3 of the CS4b verification plan
(``scratch/cs4b_verification_plan.md`` §11, step S2): the F2 doctrine's
per-family discrimination table, the scalar family's deliberate
quadrature-blindness (space half), and the role gate's relocation onto the
CLASS arm.

**The doctrine under test** (F2, RATIFIED): partner identity is space
CONTENT equality. A space is its DOF set + Gram — so two carriers differing
only in a BOUNDARY LAW mint EQUAL spaces (a law changes neither DOFs nor
Gram; laws are operator data), and the scalar bulk is additionally
quadrature-blind (no angular axis ⟹ ``S_N`` order is not its data). Volumes,
``ng``, and (for the angular family) quadrature all discriminate, because
each changes an axis's content.

⚠ Step boundary (§6b): at S2 the PARTNER GUARDS are still mesh-keyed — the
cross-mesh ADDITION permissions these equalities imply (the BC-only and
scalar-S4+S8 legs of the doctrine) become spellable only at S3's re-key, and
their gates land there. What S2 owns is the SPACE-level table and the role
consequence (G2.3), which is same-mesh and live now.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.source_sinks import AngularSourceSink
from tests.sn._test_helpers import placeholder_materials

pytestmark = pytest.mark.foundation

#: The reference configuration (verification plan §4): NON-uniform edges,
#: GL(4), ng = 2, vacuum/vacuum.
_EDGES = np.array([0.0, 0.2, 0.5, 0.9, 1.6, 3.0])
#: One interior edge moved — the volumes discriminator.
_EDGES_MOVED = np.array([0.0, 0.2, 0.5, 1.0, 1.6, 3.0])


def _mesh(
    *,
    edges: np.ndarray = _EDGES,
    n_ord: int = 4,
    ng: int = 2,
    bc_left: str = "vacuum",
) -> SNMesh:
    mesh = Mesh1D(
        edges=edges,
        mat_ids=np.zeros(edges.size - 1, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC(bc_left),
        bc_right=BC("vacuum"),
    )
    return SNMesh(
        mesh, Quadrature.gauss_legendre(n_ord), placeholder_materials(ng=ng)
    )


#: (what differs, mesh-builder kwargs, angular ``==`` expected, scalar
#: ``==`` expected) — the `[M]` table of verification plan §4.
_TABLE = [
    ("twin", {}, True, True),
    ("bc_only", {"bc_left": "reflective"}, True, True),
    ("volumes", {"edges": _EDGES_MOVED}, False, False),
    ("quadrature", {"n_ord": 8}, False, True),
    ("ng", {"ng": 3}, False, False),
]


class TestG21DiscriminationTable:
    """G2.1 — space ``==`` discriminates exactly what an axis carries.

    LAW, both legs (vv #19): the EQUAL rows (twin, BC-only, scalar-quad)
    are the negative-control half — they prove the gate is CONTENT
    equality, not provenance identity (every row's meshes are distinct
    instances, so a provenance gate would read False everywhere).
    """

    @pytest.mark.parametrize(
        "label, kwargs, angular_eq, scalar_eq",
        _TABLE,
        ids=[row[0] for row in _TABLE],
    )
    def test_per_family_space_equality(self, label, kwargs, angular_eq, scalar_eq):
        a, b = _mesh(), _mesh(**kwargs)
        assert (a.angular_bulk_space == b.angular_bulk_space) is angular_eq
        assert (a.bulk_space == b.bulk_space) is scalar_eq


class TestG22ScalarQuadratureBlindness:
    """G2.2 — the scalar family is DELIBERATELY quadrature-blind.

    A recorded PERMISSION, not an accident: the scalar bulk has no angular
    axis, so ``S_N`` order is not part of its DOF set or Gram — φ from an
    S4 solve and φ from an S8 solve on the same grid are elements of the
    SAME space (ratified with the F2 ruling; strictly wider than the
    BC-only permission, and the one a physicist will challenge first —
    the DOF+Gram criterion is the rationale). Its control is G2.1's
    angular quadrature row, which must stay UNEQUAL.

    ⚠ The addition leg (``phi_S4 + phi_S8`` SUCCEEDS) lands at S3 with the
    partner-guard re-key — until then the mesh arm still refuses it.
    """

    def test_scalar_space_ignores_quadrature_order_angular_does_not(self):
        s4, s8 = _mesh(n_ord=4), _mesh(n_ord=8)
        assert s4.bulk_space == s8.bulk_space
        assert s4.angular_bulk_space != s8.angular_bulk_space


class TestG23RoleIsClassIdentity:
    """G2.3 — role stays gated by CLASS after the space stops seeing it.

    ADMISSION, both legs (vv #11). Until CS4b the role spaces were
    non-``==`` (per-leaf ``_SPACE_NAME`` tags), so deleting the class arm
    of ``Field._check_partner`` was MASKED by the space arm. After F1-sub
    the family shares ONE cached space instance, so this gate is genuinely
    NEW coverage: the class arm is the sole role enforcement (battery arm
    M5 reds exactly here and nowhere else).
    """

    def test_cross_role_addition_refuses_while_spaces_are_identical(self):
        sn = _mesh()
        flux = AngularFlux.zeros_on(sn)
        source = AngularSourceSink.zeros_on(sn)
        # The premise that makes this NEW coverage: one space, two roles.
        assert flux.space is source.space
        with pytest.raises(TypeError):
            flux + source  # type: ignore[operator]  # the refusal IS the claim

    def test_same_role_addition_still_succeeds(self):
        sn = _mesh()
        a = AngularFlux.from_mesh(np.ones((sn.quad.N, sn.ng, *sn.spatial_shape)), sn)
        b = AngularFlux.from_mesh(2.0 * np.ones_like(a.values), sn)
        c = a + b
        assert isinstance(c, AngularFlux)
        assert np.all(c.values == 3.0)
