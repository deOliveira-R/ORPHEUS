r"""Foundation suite for the optional spatial-moment field-space factor (#240 D5b-S3-A0).

The field-space factories (``AngularField`` / ``ScalarField`` un-windowed
carriers + ``HarmonicMomentFlux`` windowed carrier) gained an OPTIONAL
``spatial_moments`` parameter that composes a
:class:`~orpheus.numerics.spaces.spatial_moment_space.SpatialMomentSpace`
factor onto the field space — gated "append iff > 1" so the default leaves
the space BYTE-IDENTICAL.

These are ``foundation`` (software-invariant) tests. The load-bearing gate
is **byte-identity at default for ALL schemes** (DD, Step, AND LD): no
production field carries the axis yet (construct-general / select-narrow),
so widening the factories must not change ANY live field shape. The widened
path is the CAPABILITY, exercised here by passing ``spatial_moments`` > 1
explicitly.

Mode-8 / L26: every assertion is a FUNCTION CALL (``np.testing.*`` /
``pytest.fail`` / ``pytest.raises``) — bare ``assert`` is a NO-OP under the
canonical ``-O`` invocation. Structural independence (L11): expected shapes
are hand-built from the mesh's own ``(N/ng, *spatial)`` + an
independently-computed ``per_axis ** ndim`` tail, NOT read off the field.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D, Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.spaces import SpatialMomentSpace, SphericalHarmonicSpace
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.transport.spatial import DiamondDifference, LinearDiscontinuous
from orpheus.transport.fields import HarmonicMomentFlux
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.scalar_flux import ScalarFlux

from tests.sn._test_helpers import placeholder_materials


def _check(cond: bool, msg: str) -> None:
    """Mode-8-safe boolean assertion (a function call, fires under ``-O``)."""
    if not cond:
        pytest.fail(msg)


# ─────────────────────────────────────────────────────────────────────
# Mesh fixtures: a 2-D Cartesian mesh under each scheme + a 1-D LD mesh.
# ─────────────────────────────────────────────────────────────────────


def _mesh_2d(scheme):
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 2.0, 4),   # nx = 3
        edges_y=np.linspace(0.0, 3.0, 5),   # ny = 4
        mat_map=np.zeros((3, 4), dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_xmin=BC("reflective"), bc_xmax=BC("reflective"),
        bc_ymin=BC("reflective"), bc_ymax=BC("reflective"),
    )
    return SNMesh(
        mesh, Quadrature.level_symmetric(4), placeholder_materials(ng=2),
        scheme=scheme,
    )


def _mesh_1d(scheme):
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, 6),     # nx = 5
        mat_ids=np.zeros(5, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"), bc_right=BC("vacuum"),
    )
    return SNMesh(
        mesh, Quadrature.gauss_legendre(4), placeholder_materials(ng=2),
        scheme=scheme,
    )


@pytest.fixture
def dd_2d():
    return _mesh_2d(DiamondDifference())


@pytest.fixture
def ld_2d():
    return _mesh_2d(LinearDiscontinuous())


@pytest.fixture
def ld_1d():
    return _mesh_1d(LinearDiscontinuous())


# ─────────────────────────────────────────────────────────────────────
# (c) BYTE-IDENTITY AT DEFAULT — the negative control (load-bearing).
#
# No production field carries the axis yet, so the default factory output
# must be the EXACT pre-S3 space + shape for EVERY scheme (DD/Step AND LD).
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
@pytest.mark.parametrize("scheme_name", ["dd", "ld"])
def test_angular_flux_default_byte_identical_all_schemes(scheme_name, dd_2d, ld_2d):
    r"""``AngularFlux.zeros_on`` default space == the pre-S3 ``(N, ng, *spatial)``.

    The default (``spatial_moments=1``) appends NO factor regardless of the
    mesh's scheme — DD and LD produce the IDENTICAL space. Pinned against
    the independently-built expected shape from the mesh's own dims.
    """
    mesh = {"dd": dd_2d, "ld": ld_2d}[scheme_name]
    field = AngularFlux.zeros_on(mesh)
    expected = (mesh.quad.N, mesh.ng, *mesh.spatial_shape)
    np.testing.assert_equal(field.space.shape, expected)
    np.testing.assert_equal(field.values.shape, expected)
    # the space is a bare FunctionSpace (no tensor-product factor at default)
    _check(not hasattr(field.space, "factors"), "default space must be bare")


@pytest.mark.foundation
@pytest.mark.parametrize("scheme_name", ["dd", "ld"])
def test_scalar_flux_default_byte_identical_all_schemes(scheme_name, dd_2d, ld_2d):
    r"""``ScalarFlux.zeros_on`` default space == the pre-S3 ``(ng, *spatial)``.

    The :class:`ScalarSourceSink` scattering accumulator shares this
    ``ScalarField`` factory; the default must stay byte-identical for both
    schemes (the negative control for the scattering-source widening).
    """
    mesh = {"dd": dd_2d, "ld": ld_2d}[scheme_name]
    field = ScalarFlux.zeros_on(mesh)
    expected = (mesh.ng, *mesh.spatial_shape)
    np.testing.assert_equal(field.space.shape, expected)
    np.testing.assert_equal(field.values.shape, expected)
    _check(not hasattr(field.space, "factors"), "default space must be bare")


@pytest.mark.foundation
@pytest.mark.parametrize("scheme_name", ["dd", "ld"])
def test_harmonic_moment_flux_default_byte_identical(scheme_name, dd_2d, ld_2d):
    r"""``HarmonicMomentFlux.zeros_for_mesh_and_L`` default == pre-S3 shape.

    The windowed iterate carrier. Default ``spatial_moments=1`` →
    ``(L+1, 2L+1, ng, *spatial)`` with NO trailing spatial-moment axis, AND
    the composition tree carries only the angular ``SphericalHarmonicSpace``
    factor (the spatial factor is absent — ``find_factor`` raises for it).
    """
    mesh = {"dd": dd_2d, "ld": ld_2d}[scheme_name]
    L = 1
    field = HarmonicMomentFlux.zeros_for_mesh_and_L(mesh, L)
    expected = (L + 1, 2 * L + 1, mesh.ng, *mesh.spatial_shape)
    np.testing.assert_equal(field.space.shape, expected)
    np.testing.assert_equal(field.values.shape, expected)
    np.testing.assert_equal(field.spatial_moments, 1)
    # the angular factor is present; the spatial factor is NOT (byte-id tree)
    np.testing.assert_equal(field.space.find_factor(SphericalHarmonicSpace).L, L)
    with pytest.raises(KeyError):
        field.space.find_factor(SpatialMomentSpace)


# ─────────────────────────────────────────────────────────────────────
# (d) WIDENED PATH — explicit spatial_moments > 1 composes the factor.
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
@pytest.mark.parametrize(
    "field_factory",
    [
        pytest.param(lambda m, n: AngularFlux.zeros_on(m, spatial_moments=n),
                     id="angular_flux"),
        pytest.param(lambda m, n: ScalarFlux.zeros_on(m, spatial_moments=n),
                     id="scalar_flux"),
    ],
)
def test_bulk_field_widened_2d_shape(field_factory, ld_2d):
    r"""A widened bulk field gets a trailing ``per_axis ** ndim`` axis (d=2).

    For ``spatial_moments=2`` on a 2-D mesh the trailing axis is
    ``2 ** 2 = 4``; the factor is ``find_factor``-queryable. The expected
    trailing length is recomputed inline (not read off the field).
    """
    mesh = ld_2d
    per_axis = 2
    field = field_factory(mesh, per_axis)
    independent_tail = (per_axis ** mesh.ndim,)
    np.testing.assert_equal(field.space.shape[-1:], independent_tail)
    np.testing.assert_equal(field.values.shape, field.space.shape)
    factor = field.space.find_factor(SpatialMomentSpace)
    np.testing.assert_equal(factor.per_axis, per_axis)
    np.testing.assert_equal(factor.ndim, mesh.ndim)
    np.testing.assert_equal(factor.n_moments, per_axis ** mesh.ndim)


@pytest.mark.foundation
def test_angular_flux_widened_1d_shape(ld_1d):
    r"""A widened 1-D LD bulk field gets a trailing ``per_axis ** 1`` axis.

    The 1-D scan carrier: ``spatial_moments=2`` on a 1-D mesh → trailing
    axis of ``2`` (the ``[bar, slope]`` per-cell moment pair).
    """
    mesh = ld_1d
    field = AngularFlux.zeros_on(mesh, spatial_moments=2)
    np.testing.assert_equal(mesh.ndim, 1)
    np.testing.assert_equal(field.space.shape[-1], 2)
    np.testing.assert_equal(field.values.shape, field.space.shape)
    np.testing.assert_equal(field.space.find_factor(SpatialMomentSpace).per_axis, 2)


@pytest.mark.foundation
def test_harmonic_moment_flux_widened_2d_shape(ld_2d):
    r"""A widened windowed iterate gets a trailing ``per_axis ** ndim`` axis.

    ``spatial_moments=2`` on a 2-D mesh → ``(L+1, 2L+1, ng, *spatial, 4)``;
    BOTH moment factors (angular SH + spatial) are queryable by type — the
    orthogonal-axes invariant on the live carrier.
    """
    mesh = ld_2d
    L = 1
    field = HarmonicMomentFlux.zeros_for_mesh_and_L(mesh, L, spatial_moments=2)
    expected = (L + 1, 2 * L + 1, mesh.ng, *mesh.spatial_shape, 2 ** mesh.ndim)
    np.testing.assert_equal(field.space.shape, expected)
    np.testing.assert_equal(field.values.shape, expected)
    np.testing.assert_equal(field.spatial_moments, 2)
    # both moment factors coexist (orthogonal axes)
    np.testing.assert_equal(field.space.find_factor(SphericalHarmonicSpace).L, L)
    np.testing.assert_equal(field.space.find_factor(SpatialMomentSpace).per_axis, 2)


@pytest.mark.foundation
def test_from_mesh_widened_roundtrip(ld_2d):
    r"""``from_mesh`` accepts a pre-shaped widened buffer and validates it.

    The reconstruction path (the iterate carrier S3-A will use): a buffer
    already carrying the trailing moment axis round-trips through
    ``from_mesh(values, mesh, spatial_moments=2)`` and passes the
    ``BulkField`` shape gate (the validator derives the widened expected
    shape from the composed space — Pattern 4).
    """
    mesh = ld_2d
    shape = (mesh.quad.N, mesh.ng, *mesh.spatial_shape, 2 ** mesh.ndim)
    values = np.arange(np.prod(shape), dtype=np.float64).reshape(shape)
    field = AngularFlux.from_mesh(values, mesh, spatial_moments=2)
    np.testing.assert_equal(field.space.shape, shape)
    np.testing.assert_array_equal(field.values, values)


@pytest.mark.foundation
def test_widened_field_rejects_wrong_shape_buffer(ld_2d):
    r"""A buffer whose trailing axis disagrees with ``spatial_moments`` raises.

    Pattern 4: the ``BulkField`` shape gate cross-checks ``values.shape``
    against the widened expected shape; a mismatched trailing axis is an
    illegal state. Production invariant → real ``raise`` (fires under ``-O``).
    """
    mesh = ld_2d
    # request spatial_moments=2 (tail 4) but feed a tail-2 buffer
    bad_shape = (mesh.quad.N, mesh.ng, *mesh.spatial_shape, 2)
    bad = np.zeros(bad_shape)
    with pytest.raises(ValueError):
        AngularFlux.from_mesh(bad, mesh, spatial_moments=2)
