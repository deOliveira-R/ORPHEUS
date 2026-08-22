r"""S1 gates for the carrier's angular-bulk space mint (campaign 1 CS4b).

G1.1–G1.5 of the CS4b verification plan
(``scratch/cs4b_verification_plan.md`` §11, step S1) plus the scheme-side
``moment_axis`` admission pair. The step is provably behaviour-neutral —
nothing consumes :attr:`SNMesh.angular_bulk_space` yet — so every gate here
is either a RECORD of the mint's content, a LAW comparing it against the
SHIPPED dense composite interior (the §6c witness that exists today), or an
ADMISSION with both legs (vv #11).

Conventions gated (CS4b crosswalk B1/B5, ``.claude/plans/cs4b_crosswalk.md``):

* the axis order is ``(angular, energy, spatial)``, matching the bulk tensor
  ``(N, ng, *spatial)``;
* the energy and spatial arms are ``bulk_space``'s axes REUSED VERBATIM
  (object identity — the energy-arm rule is spelled once);
* the Gram of the axis product equals ``full_field_space``'s dense interior
  on both the DD and the LD arm (LD composes the scheme-owned MODAL
  ``moment_axis`` carrying ``moment_mass_diagonal``);
* the derived space NAME is never pinned (R4: CS2's typed axis subclasses
  change the digest — every assertion here is per-axis content or relative
  ``is``/``==``);
* cone predicates: nodal bulk families answer ``True``, the harmonic-moment
  and trace families ``None``, the LD moment-tailed product ``False`` (the
  Q6-ratified routing base for ``cone_violations``; the exact vertex test
  is #400).

Fixture: the verification plan's §2 configuration — a NON-uniform 5-cell
slab (uniform volumes would collapse ``V`` to a scalar and blind half the
metric claims), ``gauss_legendre(4)``, ``ng = 2``, vacuum/vacuum.
"""

from __future__ import annotations

import numpy as np
import numpy.testing as npt
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.axis import BasisKind, EnergyAxis
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.space import FunctionSpace
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux
from orpheus.transport.spatial import LinearDiscontinuous
from orpheus.transport.spatial.diamond import DiamondDifference
from tests.sn._test_helpers import placeholder_materials

pytestmark = pytest.mark.foundation

#: NON-uniform edges — ``V = [0.2, 0.3, 0.4, 0.7, 1.4]``, a genuine vector.
_EDGES = np.array([0.0, 0.2, 0.5, 0.9, 1.6, 3.0])
_NG = 2


def _slab(*, scheme=None, ng: int = _NG) -> SNMesh:
    mesh = Mesh1D(
        edges=_EDGES,
        mat_ids=np.zeros(_EDGES.size - 1, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    kwargs = {} if scheme is None else {"scheme": scheme}
    return SNMesh(
        mesh, Quadrature.gauss_legendre(4), placeholder_materials(ng=ng), **kwargs
    )


class TestG11AxisTuple:
    """G1.1 — the axis tuple IS (angular w_n, energy, spatial V). RECORD."""

    def test_axes_are_angular_energy_spatial_with_the_carrier_measures(self):
        sn = _slab()
        space = sn.angular_bulk_space
        assert space.axes is not None and len(space.axes) == 3
        angular, energy, spatial = space.axes

        assert angular.label == "angular"
        assert angular.shape == (sn.quad.N,)
        assert angular.kind is BasisKind.NODAL
        assert angular.weights is not None
        assert np.array_equal(angular.weights, sn.quad.weights)

        assert isinstance(energy, EnergyAxis)
        assert energy.shape == (sn.ng,)

        assert spatial.label == "spatial"
        assert spatial.shape == sn.spatial_shape
        assert spatial.kind is BasisKind.NODAL
        assert spatial.weights is not None
        assert np.array_equal(spatial.weights, sn.volumes)

        assert space.shape == (sn.quad.N, sn.ng, *sn.spatial_shape)

    def test_energy_and_spatial_arms_are_bulk_space_axes_verbatim(self):
        """The scalar arms are REUSED objects, not respelled twins — the
        energy-arm rule (``EnergyAxis.from_materials``) is spelled exactly
        once, in ``bulk_space`` (Pattern 2)."""
        sn = _slab()
        scalar_axes = sn.bulk_space.axes
        assert scalar_axes is not None
        assert sn.angular_bulk_space.axes is not None
        assert sn.angular_bulk_space.axes[1] is scalar_axes[0]
        assert sn.angular_bulk_space.axes[2] is scalar_axes[1]


class TestG12Cache:
    """G1.2 — the mint is CACHED. LAW (the is/== asymmetry IS the gate)."""

    def test_same_carrier_reads_the_same_instance(self):
        sn = _slab()
        assert sn.angular_bulk_space is sn.angular_bulk_space

    def test_twin_carriers_mint_equal_but_distinct_spaces(self):
        a, b = _slab(), _slab()
        assert a.angular_bulk_space == b.angular_bulk_space
        assert a.angular_bulk_space is not b.angular_bulk_space


class TestG13GramEquivalenceDD:
    """G1.3 — the axis product's Gram equals the dense composite interior's
    (``full_field_space.interior_space``, G_bulk = V·w_n), on the DD arm.
    LAW; the dense interior is the §6c witness that ships TODAY."""

    def test_all_three_metric_faces_agree_with_the_dense_interior(self):
        sn = _slab()
        dense = sn.full_field_space.interior_space
        assert dense is not None
        axis_built = sn.angular_bulk_space
        assert axis_built.shape == dense.shape

        rng = np.random.default_rng(0)
        x = rng.standard_normal(dense.shape)
        y = rng.standard_normal(dense.shape)

        # Scalar faces: bit-equal ([M] verification plan G1.3 — rel diff 0.0).
        assert axis_built.inner_product(x, y) == dense.inner_product(x, y)
        assert axis_built.norm(x) == dense.norm(x)
        # Vector face: ≤ 4 ulp ([M] max abs Δ 2.78e-17 on this fixture).
        npt.assert_array_almost_equal_nulp(
            axis_built.apply_metric(x), dense.apply_metric(x), nulp=4
        )


class TestG14GramEquivalenceLD:
    """G1.4 — the LD arm: the dense interior carries the scheme's moment
    mass on the trailing 2^d axis; the axis form composes the scheme-owned
    MODAL ``moment_axis`` and reproduces the same Gram. LAW.

    [M] R9 measured its draw's inner product bit-identical; that was the
    draw's luck, not a law — the two spellings associate the weight
    products differently, and on THIS fixture's ``rng(0)`` draw the
    near-cancelling bilinear form lands 6 ULP apart (measured 2026-08-22).
    The honest bound is nulp ≤ 64 on the cancellation-conditioned scalar
    face (vv #16: never assert tighter than construction gives); a
    mass-placement error is O(θ·value) ≈ 1e14 ULP, so the gate's
    discrimination is unharmed."""

    def test_all_three_metric_faces_agree_on_the_ld_interior(self):
        sn = _slab(scheme=LinearDiscontinuous())
        dense = sn.full_field_space.interior_space
        assert dense is not None
        base = sn.angular_bulk_space
        assert base.axes is not None
        widened = FunctionSpace.of_axes(
            *base.axes, sn.scheme.moment_axis(sn.ndim)
        )
        assert widened.shape == dense.shape

        rng = np.random.default_rng(0)
        x = rng.standard_normal(dense.shape)
        y = rng.standard_normal(dense.shape)

        # Cancellation-conditioned scalar face: standard-normal x·G·y sums
        # ~80 O(1) signed terms to ~0.2, so association differences amplify
        # in result-relative ULPs ([M] 6 here).
        npt.assert_array_almost_equal_nulp(
            np.array([widened.inner_product(x, y)]),
            np.array([dense.inner_product(x, y)]),
            nulp=64,
        )
        # Positive-term faces: no cancellation, tight.
        npt.assert_array_almost_equal_nulp(
            np.array([widened.norm(x)]), np.array([dense.norm(x)]), nulp=4
        )
        npt.assert_array_almost_equal_nulp(
            widened.apply_metric(x), dense.apply_metric(x), nulp=4
        )


class TestMomentAxisAdmission:
    """The scheme-side mint's ADMISSION pair (vv #11: both legs)."""

    def test_ld_mints_the_modal_mass_axis(self):
        scheme = LinearDiscontinuous()
        axis = scheme.moment_axis(1)
        assert axis.label == "spatial_moment"
        assert axis.shape == (2,)
        assert axis.kind is BasisKind.MODAL
        assert axis.weights is not None
        assert np.array_equal(axis.weights, scheme.moment_mass_diagonal(1))

    def test_slopeless_closure_refuses(self):
        with pytest.raises(ValueError, match="no moment axis"):
            DiamondDifference().moment_axis(1)


class TestG15ConePredicates:
    """G1.5 — the cone predicates, stated. RECORD ([M] verification plan
    Finding 8 + R9)."""

    def test_nodal_bulk_families_answer_true(self):
        sn = _slab()
        assert sn.angular_bulk_space.has_coordinate_cone is True
        assert sn.bulk_space.has_coordinate_cone is True

    def test_harmonic_moment_and_trace_families_answer_none(self):
        sn = _slab()
        moment_space = HarmonicMomentFlux.zeros_for_mesh_and_L(sn, 1).space
        assert moment_space.has_coordinate_cone is None
        assert sn.angular_trace.has_coordinate_cone is None

    def test_the_ld_moment_tail_is_modal_so_the_cone_reads_false(self):
        """The Q6-ratified routing base: a moment-tailed LD bulk space
        answers ``False`` (signed slope coefficients are legal on a
        positive function), which ``Field.cone_violations`` turns into the
        typed refusal. The exact modal test (the vertex theorem) is #400."""
        sn = _slab(scheme=LinearDiscontinuous())
        assert sn.angular_bulk_space.axes is not None
        widened = FunctionSpace.of_axes(
            *sn.angular_bulk_space.axes, sn.scheme.moment_axis(sn.ndim)
        )
        assert widened.has_coordinate_cone is False
