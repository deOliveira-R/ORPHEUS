r"""foundation — intrinsic-property gates for the coefficient field leaves.

Pins the defining algebraic structure of the #257 S1 coefficient leaf —
NOT happy-path usage, but the laws/invariants that DEFINE the concept
(``feedback_test_intrinsic_properties``: every math-bearing type ships a test
of its intrinsic properties):

* :class:`~orpheus.transport.fields.cross_section_field.CrossSectionField` —
  a **cone** in a vector space: closed under ``+`` and ``λ≥0·``, with an
  origin ``Σ=0``; nonnegativity is a tested property, NOT a construction
  invariant (so a signed ``Σ−Σ′`` perturbation stays a coefficient).
* it carries :class:`~orpheus.transport.fields._coefficient_role.CoefficientRole`,
  the COMPLEMENT of
  :class:`~orpheus.transport.fields._flux_role.FluxRole` — it has NO
  ``flux+flux→TypeError`` affine gate.

The fission emission spectrum χ is NOT a field here: its probability-simplex
invariant (``Σ_g χ_g = 1``, ``χ≥0``) is a property of the *source* — the
per-fissile-material :class:`~orpheus.data.macro_xs.mixture.Mixture` ``chi``,
the thing that generates the per-cell broadcast — and is enforced there
(deferred to the end of #257). A downstream ``SpectrumField`` would have no
native behaviour beyond that invariant, so it does not exist (see #257).

``foundation`` — software invariants on the algebra's type surface, no
theory-page ``:label:`` (so no ``verifies(...)``).

References
----------

* ``.claude/plans/issue_257_coefficient_field_promotion.md`` — S1.
* ``.claude/agent-memory/cross-domain-attacker/coefficient_field_promotion_frames.md``
  — Frame 2 (the coefficient cone algebra is NOT the flux torsor).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.units import CROSS_SECTION_UNITS
from orpheus.numerics.vector import Vector
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.transport.fields._coefficient_role import CoefficientRole
from orpheus.transport.fields._flux_role import FluxRole
from orpheus.transport.fields.cross_section_field import CrossSectionField

from tests.sn._test_helpers import placeholder_materials

pytestmark = [pytest.mark.foundation]


def _slab_mesh(nx: int = 4, ng: int = 2) -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _sigma(mesh: SNMesh, fill: float) -> CrossSectionField:
    return CrossSectionField.from_ndarray(
        np.full((mesh.ng, *mesh.spatial_shape), fill), mesh
    )


# ═══════════════════════════════════════════════════════════════════════
# CrossSectionField — the cone in a vector space
# ═══════════════════════════════════════════════════════════════════════


class TestCrossSectionConeAlgebra:
    def test_cone_closed_under_addition(self) -> None:
        """``Σ + Σ′`` is a CrossSectionField (cone closed under +) — NO
        ``flux+flux``-style gate, and the result is the pointwise sum."""
        m = _slab_mesh()
        s1, s2 = _sigma(m, 0.5), _sigma(m, 0.3)
        out = s1 + s2  # must NOT raise (unlike FluxRole)
        assert isinstance(out, CrossSectionField)
        np.testing.assert_array_equal(out.values, s1.values + s2.values)
        assert np.all(out.values >= 0.0)  # cone closed under +

    def test_cone_closed_under_nonneg_scalar(self) -> None:
        """``λ·Σ`` for ``λ≥0`` stays in the cone (closed under nonneg ·)."""
        m = _slab_mesh()
        s = _sigma(m, 0.5)
        out = 2.0 * s
        assert isinstance(out, CrossSectionField)
        np.testing.assert_array_equal(out.values, 2.0 * s.values)
        assert np.all(out.values >= 0.0)

    def test_cone_has_origin(self) -> None:
        """``Σ=0`` is a valid CrossSectionField and the additive identity —
        the coefficient space HAS an origin (it promotes to ``M_0``),
        unlike the flux affine torsor."""
        m = _slab_mesh()
        zero = CrossSectionField.zeros_on(m)  # constructs (no invariant rejects 0)
        assert isinstance(zero, CrossSectionField)
        np.testing.assert_array_equal(zero.values, np.zeros((m.ng, *m.spatial_shape)))
        s = _sigma(m, 0.5)
        np.testing.assert_array_equal((s + zero).values, s.values)

    def test_is_vector_space_not_torsor(self) -> None:
        """Coefficients keep the plain vector-space dunders: ``Σ − Σ′`` returns
        a CrossSectionField (NOT a displacement, as FluxRole would), and may be
        SIGNED — nonnegativity is the physical cone (a property), not a
        construction invariant. This keeps the multiplier-algebra domain a full
        vector space."""
        m = _slab_mesh()
        s1, s2 = _sigma(m, 0.3), _sigma(m, 0.5)
        diff = s1 - s2  # 0.3 - 0.5 = -0.2: a signed coefficient, must NOT raise
        assert isinstance(diff, CrossSectionField)
        assert np.all(diff.values < 0.0)  # signed value is representable

    def test_units_are_inverse_cm(self) -> None:
        assert CrossSectionField.UNITS == CROSS_SECTION_UNITS


# ═══════════════════════════════════════════════════════════════════════
# The role distinction — CoefficientRole is NOT FluxRole
# ═══════════════════════════════════════════════════════════════════════


class TestCoefficientRoleVsFluxRole:
    def test_carries_coefficient_role_not_flux_role(self) -> None:
        assert issubclass(CrossSectionField, CoefficientRole)
        assert not issubclass(CrossSectionField, FluxRole)

    def test_cross_section_addition_has_no_affine_gate(self) -> None:
        """The behavioural counterpart: ``Σ + Σ`` does NOT raise (no torsor
        gate), where a flux state's ``ψ + ψ`` would — the role difference made
        observable in the algebra, not just the type tree."""
        m = _slab_mesh()
        s = _sigma(m, 0.5)
        out = s + s  # CoefficientRole adds no affine gate
        assert isinstance(out, CrossSectionField)
        np.testing.assert_array_equal(out.values, 2.0 * s.values)


# ═══════════════════════════════════════════════════════════════════════
# The leaf flows through the #256 operator-algebra vector type
# ═══════════════════════════════════════════════════════════════════════


class TestCoefficientLeafIsVector:
    def test_cross_section_satisfies_vector(self) -> None:
        """``CrossSectionField`` structurally satisfies the #256 ``Vector``
        protocol (the four vector-space dunders), so it flows through the
        operator algebra's carrier type."""
        assert isinstance(_sigma(_slab_mesh(), 0.5), Vector)


# ═══════════════════════════════════════════════════════════════════════
# Mesh binding inherited from BulkField (the coefficient leaf composed it)
# ═══════════════════════════════════════════════════════════════════════


class TestCoefficientMeshBinding:
    def test_cross_section_cross_mesh_arithmetic_rejected(self) -> None:
        """Two CrossSectionFields on DISTINCT meshes cannot combine even when
        shapes match — the inherited BulkField mesh-identity gate fires."""
        m1, m2 = _slab_mesh(), _slab_mesh()
        s1, s2 = _sigma(m1, 0.5), _sigma(m2, 0.5)
        with pytest.raises(ValueError):
            _ = s1 + s2
