r"""TransportMethod Protocol — structural conformance of the method-meshes (#290 P7b).

The Protocol's whole design is that conformance is STRUCTURAL:
``orpheus.sn`` and ``orpheus.diffusion`` never import
:class:`~orpheus.transport.method.TransportMethod`; pyright checks
each mesh where it calls the shared
:func:`~orpheus.transport.method.resolve_boundary_conditions` body.
These are the RUNTIME twins of that static check:

* **positive** — both witnesses (:class:`SNMesh`,
  :class:`DiffusionMesh`) satisfy the runtime-checkable Protocol;
* **negative** — the bare :class:`MaterialMesh` data carrier does NOT
  (it has the data block but none of the method surface: no
  ``BOUNDARY_OPERATOR_REGISTRY``, no ``bc``, no
  ``realize_boundary_law``) — the Protocol genuinely discriminates
  data from method behavior (vv anti-pattern #11: a conformance claim
  needs the instance that must fail, not just the ones that must
  pass);
* **the shared refusal through the SECOND witness** — the
  unsupported-tag ValueError comes from the ONE generic
  ``_law_from_tag`` body and names the concrete mesh class; the
  diffusion side is pinned in
  ``tests/diffusion/test_augmented_mesh.py``, so HERE the SN side
  pins that the same body serves both methods.

The per-method behavior of the resolved laws (𝒜-table actions,
kind tags, default-reflective, face coverage) stays pinned in each
method's own suites (``tests/sn/operators/test_snmesh_realizer_wiring.py``,
``tests/sn/primitives/test_axis_native_construction.py``,
``tests/diffusion/test_augmented_mesh.py``) — this file owns only the
cross-method Protocol facts.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.diffusion import DiffusionMesh
from orpheus.geometry import BC, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.transport.mesh.material_mesh import MaterialMesh
from orpheus.transport.method import TransportMethod

pytestmark = [pytest.mark.foundation]

_MATS = {0: get_mixture("A", "2g")}


def _mesh1d(
    bc_left: BC | None = None, bc_right: BC | None = None,
) -> Mesh1D:
    return Mesh1D(
        np.linspace(0.0, 10.0, 5), np.zeros(4, dtype=int),
        bc_left=bc_left, bc_right=bc_right,
    )


def _sn_mesh(**bc_kwargs: "BC | None") -> SNMesh:
    return SNMesh(
        _mesh1d(**bc_kwargs), Quadrature.gauss_legendre(4), _MATS,
    )


def _diffusion_mesh(**bc_kwargs: "BC | None") -> DiffusionMesh:
    return DiffusionMesh(_mesh1d(**bc_kwargs), _MATS)


class TestStructuralConformance:
    @pytest.mark.parametrize(
        "build", [_sn_mesh, _diffusion_mesh],
        ids=["SNMesh", "DiffusionMesh"],
    )
    def test_method_mesh_satisfies_the_protocol(self, build):
        """Every method-mesh IS a TransportMethod — checked at runtime
        against the same member list pyright verifies statically at
        the shared resolve call."""
        assert isinstance(build(), TransportMethod)

    def test_bare_material_mesh_is_not_a_transport_method(self):
        """The data carrier alone is NOT a method: the Protocol names
        the method surface (admission table, realized ``bc``, the
        ``realize_boundary_law`` arm), none of which the
        data/behavior axis puts on :class:`MaterialMesh`."""
        assert not isinstance(
            MaterialMesh(_mesh1d(), _MATS), TransportMethod,
        )


class TestSharedResolveBody:
    def test_sn_unsupported_tag_refuses_naming_the_mesh_class(self):
        """The unsupported-tag refusal is the ONE generic body
        (``_law_from_tag``) serving both witnesses: through the SN
        witness it names ``SNMesh`` and SN's own admission table.
        (``albedo`` is admitted by diffusion but NOT by SN — the
        method-specific table drives the shared body.)"""
        with pytest.raises(
            ValueError,
            match=(
                r"SNMesh does not support boundary condition 'albedo' "
                r"on face 'xmax'\. Supported: 'reflective', 'vacuum'\."
            ),
        ):
            _sn_mesh(bc_right=BC("albedo", {"albedo": 0.5}))
