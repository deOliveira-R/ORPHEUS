r"""TransportMethod Protocol — structural conformance of the method-meshes (#290 P7b).

The Protocol's whole design is that conformance is STRUCTURAL:
``orpheus.sn`` and ``orpheus.diffusion`` never import
:class:`~orpheus.transport.method.TransportMethod`; pyright checks
each mesh where it calls the shared
:func:`~orpheus.transport.method.realize_boundary` body.
These are the RUNTIME twins of that static check:

* **positive** — both witnesses (:class:`SNProblem`,
  :class:`DiffusionMesh`) satisfy the runtime-checkable Protocol;
* **negative** — the bare :class:`MaterialMesh` data carrier does NOT
  (it has the data block but none of the method surface: no
  ``boundary_realizers``, no ``binding``, no ``bc``, no
  ``realize_boundary_role``) — the Protocol genuinely discriminates
  data from method behavior (vv anti-pattern #11: a conformance claim
  needs the instance that must fail, not just the ones that must
  pass);
* **the shared refusal through both witnesses** — a declaration the
  method's own table does not admit comes from the ONE generic
  ``realize_boundary`` body as a ``BoundaryRefused`` naming the method;
  the per-condition admissions are pinned in each table's own gates
  (``tests/gates/diffusion/test_diffusion_realizer_table.py``,
  ``tests/gates/sn/architecture/test_sn_realizer_table.py``), so HERE both
  witnesses pin that the same body serves both methods.

The per-method behavior of the realized conditions (𝒜-table actions,
face coverage) stays pinned in each
method's own suites (``tests/gates/sn/operators/test_snmesh_realizer_wiring.py``,
``tests/gates/sn/primitives/test_axis_native_construction.py``,
``tests/gates/diffusion/test_augmented_mesh.py``) — this file owns only the
cross-method Protocol facts.
"""

from __future__ import annotations

import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.diffusion import DiffusionMesh
from orpheus.diffusion import ZeroFluxCondition
from orpheus.geometry import StructuredGeometry
from orpheus.geometry.deck import MirrorDeck
from orpheus.mesh import CellsByCount, Mesh1D, Mesher
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.problem import SNProblem
from orpheus.transport.boundary_conditions import BoundaryCondition, BoundaryRefused, WallCondition
from orpheus.transport.mesh.material_mesh import MaterialMesh
from orpheus.transport.method import TransportMethod

pytestmark = [pytest.mark.foundation]

_MATS = {0: get_mixture("A", "2g")}


def _mesh1d(right: "MirrorDeck | int" = MirrorDeck()) -> Mesh1D:
    geometry = StructuredGeometry.slab(
        (0.0, 10.0), (0,), left=MirrorDeck(), right=right,
    )
    return Mesher(geometry).partition(CellsByCount.uniform_width(4)).mesh


def _sn_mesh(*, right: BoundaryCondition | None = None) -> SNProblem:
    """A slab with a mirror on the left; the right face is a mirror, or tag 0 bound to ``right``."""
    return SNProblem(
        _mesh1d(**({} if right is None else {"right": 0})), Quadrature.gauss_legendre(4), _MATS,
        boundary_conditions={} if right is None else {0: right},
    )


def _diffusion_mesh(*, right: BoundaryCondition | None = None) -> DiffusionMesh:
    return DiffusionMesh(
        _mesh1d(**({} if right is None else {"right": 0})), _MATS,
        boundary_conditions={} if right is None else {0: right},
    )


class TestStructuralConformance:
    @pytest.mark.parametrize(
        "build", [_sn_mesh, _diffusion_mesh],
        ids=["SNProblem", "DiffusionMesh"],
    )
    def test_method_mesh_satisfies_the_protocol(self, build):
        """Every method-mesh IS a TransportMethod — checked at runtime
        against the same member list pyright verifies statically at
        the shared realize call."""
        assert isinstance(build(), TransportMethod)

    def test_bare_material_mesh_is_not_a_transport_method(self):
        """The data carrier alone is NOT a method: the Protocol names
        the method surface (admission table, realized ``bc``, the
        ``realize_boundary_law`` arm), none of which the
        data/behavior axis puts on :class:`MaterialMesh`."""
        assert not isinstance(
            MaterialMesh(_mesh1d(), _MATS), TransportMethod,
        )


class TestSharedRealizeBody:
    """The ONE generic ``realize_boundary`` body serves both witnesses: a declaration the method's own table does not admit
    refuses through ``BoundaryRefused`` naming that method (the table's refusal message), whichever mesh class called it."""

    def test_sn_refuses_what_its_table_does_not_admit_naming_the_method(self):
        """Diffusion's zero flux is outside the non-negative cone SN's state lives in: the table of ``sn`` refuses it on the
        state-cone axis, without importing the diffusion type (the unregistered-unit door)."""
        with pytest.raises(BoundaryRefused, match=r"sn realizes .*; not ZeroFluxCondition \(STATE_CONE") as info:
            _sn_mesh(right=ZeroFluxCondition())
        assert info.value.method == "sn" and info.value.face == "xmax"

    def test_diffusion_refuses_what_its_table_does_not_admit_naming_the_method(self):
        """A boundary source has no channel in diffusion (#656): the table of ``diffusion`` refuses it on the plumbing axis."""
        with pytest.raises(BoundaryRefused, match=r"diffusion realizes .*; not WallCondition \(PLUMBING") as info:
            _diffusion_mesh(right=WallCondition.uniform(1.0))
        assert info.value.method == "diffusion" and info.value.face == "xmax"
