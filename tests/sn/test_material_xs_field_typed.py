r"""foundation — MaterialXSField typed CrossSectionField accessors (#257 S2).

The per-cell σ_t / σ_a / νΣ_f views, wrapped as the typed
:class:`~orpheus.transport.fields.cross_section_field.CrossSectionField` — the
field side of the operator promotion ``C = M[σ_t]`` (#257 §5.7). Pins:

* each typed accessor returns a ``CrossSectionField`` on the mesh, units 1/cm;
* ``.values`` IS the SAME cached ndarray as the raw view (zero-copy — a typed
  lens, NOT a second representation; bit-identical by identity, not equality);
* shape ``== (ng, *spatial)``, bound to the mesh;
* the raw ndarray accessors are UNTOUCHED (the existing consumer path).

``foundation`` — a software invariant on the type surface (the S2 bit-identity
gate), no theory ``:label:`` (so no ``verifies(...)``).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.units import CROSS_SECTION_UNITS
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.transport.mesh.material_xs_field import MaterialXSField
from orpheus.transport.fields.cross_section_field import CrossSectionField

from tests.sn._test_helpers import placeholder_materials

pytestmark = [pytest.mark.foundation]


def _mat_xs(nx: int = 4, ng: int = 2) -> tuple[MaterialXSField, SNMesh]:
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    sn_mesh = SNMesh(mesh, quad, placeholder_materials(ng=ng))
    return MaterialXSField.from_mesh(sn_mesh), sn_mesh


#: (typed accessor name, raw ndarray accessor name) for the three macroscopic
#: cross sections that promote to a multiplication operator.
ACCESSORS = [
    ("total_cross_section_field", "total_cross_section"),
    ("absorption_cross_section_field", "absorption_cross_section"),
    ("fission_production_field", "fission_production"),
]


class TestTypedCrossSectionAccessors:
    @pytest.mark.parametrize("typed, _raw", ACCESSORS)
    def test_returns_cross_section_field_in_inverse_cm(self, typed, _raw) -> None:
        mat_xs, _ = _mat_xs()
        field = getattr(mat_xs, typed)
        assert isinstance(field, CrossSectionField)
        assert field.UNITS == CROSS_SECTION_UNITS

    @pytest.mark.parametrize("typed, raw", ACCESSORS)
    def test_values_are_the_same_cached_array(self, typed, raw) -> None:
        """Zero-copy: the typed field wraps the SAME ndarray as the raw view
        (bit-identical by identity — a typed lens, not a second copy)."""
        mat_xs, _ = _mat_xs()
        field = getattr(mat_xs, typed)
        raw_view = getattr(mat_xs, raw)
        assert field.values is raw_view

    @pytest.mark.parametrize("typed, _raw", ACCESSORS)
    def test_shape_and_mesh_binding(self, typed, _raw) -> None:
        mat_xs, sn_mesh = _mat_xs()
        field = getattr(mat_xs, typed)
        assert field.values.shape == (sn_mesh.ng, *sn_mesh.spatial_shape)
        assert field.mesh is sn_mesh

    def test_raw_views_untouched(self) -> None:
        """The raw ndarray accessors still return a bare ``np.ndarray`` — the
        existing consumer path is unchanged (S2 is pure-additive)."""
        mat_xs, _ = _mat_xs()
        raw = mat_xs.total_cross_section
        assert isinstance(raw, np.ndarray)
        assert not isinstance(raw, CrossSectionField)
