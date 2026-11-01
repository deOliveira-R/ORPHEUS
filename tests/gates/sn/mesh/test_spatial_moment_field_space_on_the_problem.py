r"""Foundation suite for the optional spatial-moment field-space factor (#240 D5b-S3-A0).

Since CS4b S5 the widening is a SPACE selection, not a factory parameter:
the un-windowed carriers construct on the carrier's cached mints — the
width-1 ``angular_bulk_space`` / ``bulk_space``, or the scheme-widened
``angular_trial_space`` (the ``spatial_moments`` int retired with the sugar
tier; only the S6-pending ``HarmonicMomentFlux`` keyed factory still
threads it). The "append iff > 1" gate lives in the trial mint: slopeless
schemes collapse it to the width-1 instance, byte-identical.

These are ``foundation`` (software-invariant) tests. The load-bearing gate
is **scheme-blindness of the width-1 mints** (DD, Step, AND LD read the
same bulk spaces — construct-general / select-narrow: widening is the
CALLER's property choice, never a scheme auto-read). The widened path is
the CAPABILITY, exercised on the trial mint / the composed scalar space.

Mode-8 / L26: every assertion is a FUNCTION CALL (``np.testing.*`` /
``pytest.fail`` / ``pytest.raises``) — bare ``assert`` is a NO-OP under the
canonical ``-O`` invocation. Structural independence (L11): expected shapes
are hand-built from the mesh's own ``(N/ng, *spatial)`` + an
independently-computed ``per_axis ** ndim`` tail, NOT read off the field.

Layer DAG Phase 4 (2026-11-01): this file holds the two rows of
tests/gates/transport/fields/test_spatial_moment_field_space.py whose subject is the SN problem's moment-field
factory (zeros_for_problem_and_L at width one and two); the source file keeps the transport space rows, posed on
material meshes and the raw constructors.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import CoordSystem, StructuredGeometry
from orpheus.mesh import CellsByCount, Mesh2D, Mesher
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.basis import SPATIAL_MOMENT_AXIS_LABEL
from orpheus.numerics.spaces import head_layout
from orpheus.transport.fields._bases import BulkField
from orpheus.sn.problem import SNProblem
from orpheus.transport.spatial import DiamondDifference, LinearDiscontinuous
from orpheus.transport.fields import HarmonicMomentFlux

from tests.gates.transport._materials import placeholder_materials
from orpheus.geometry.deck import MirrorDeck
from orpheus.transport.boundary_conditions import WallCondition
from tests.gates._exterior import conditions_for

#: Exterior conditions of the tags this file's meshes carry (wave 2: a face's role is a deck or an int tag; the condition a tag
#: stands for is declared at the solve, keyword ``boundary_conditions``).
_CONDITIONS = {0: WallCondition.vacuum()}


# ─────────────────────────────────────────────────────────────────────
# Mesh fixtures: a 2-D Cartesian mesh under each scheme + a 1-D LD mesh.
# ─────────────────────────────────────────────────────────────────────


def _mesh_2d(scheme):
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 2.0, 4),   # nx = 3
        edges_y=np.linspace(0.0, 3.0, 5),   # ny = 4
        mat_map=np.zeros((3, 4), dtype=int),
        coord=CoordSystem.CARTESIAN,
        face_roles={
            "xmin": MirrorDeck(),
            "xmax": MirrorDeck(),
            "ymin": MirrorDeck(),
            "ymax": MirrorDeck(),
        },
    )
    return SNProblem(
        mesh, Quadrature.level_symmetric(4), placeholder_materials(ng=2),
        scheme=scheme, boundary_conditions=conditions_for(mesh, _CONDITIONS),
    )


def _mesh_1d(scheme):
    geometry = StructuredGeometry.slab((0.0, 1.0), (0,), left=0, right=0)
    mesh = Mesher(geometry).partition(CellsByCount.uniform_width(5)).mesh     # nx = 5
    return SNProblem(
        mesh, Quadrature.gauss_legendre(4), placeholder_materials(ng=2),
        scheme=scheme, boundary_conditions=conditions_for(mesh, _CONDITIONS),
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
def test_harmonic_moment_flux_default_byte_identical(scheme_name, dd_2d, ld_2d):
    r"""``HarmonicMomentFlux.zeros_for_problem_and_L`` default == pre-S3 shape.

    The windowed iterate carrier. Default ``spatial_moments=1`` →
    ``(L+1, 2L+1, ng, *spatial)`` with NO trailing spatial-moment axis, AND
    the composition tree carries only the angular harmonic head
    factor (no axis labelled ``spatial_moment`` rides the product — since
    CS4c step 6 item 6.2c-iii the tail is the scheme's own axis, never a
    separate class).
    """
    mesh = {"dd": dd_2d, "ld": ld_2d}[scheme_name]
    L = 1
    field = HarmonicMomentFlux.zeros_for_problem_and_L(mesh, L)
    expected = (L + 1, 2 * L + 1, mesh.ng, *mesh.spatial_shape)
    np.testing.assert_equal(field.space.shape, expected)
    np.testing.assert_equal(field.values.shape, expected)
    np.testing.assert_equal(field.spatial_moments, 1)
    # the angular factor is present; the spatial-moment axis is NOT (byte-id tree)
    np.testing.assert_equal(head_layout(field.space.factors[0]).L, L)
    assert field.space.axes is not None
    np.testing.assert_equal(
        [ax.label for ax in field.space.axes if ax.label == SPATIAL_MOMENT_AXIS_LABEL], [],
    )
    np.testing.assert_equal(BulkField.spatial_moments_per_axis_of(field.space), 1)


# ─────────────────────────────────────────────────────────────────────
# (d) WIDENED PATH — explicit spatial_moments > 1 composes the factor.
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_harmonic_moment_flux_widened_2d_shape(ld_2d):
    r"""A widened windowed iterate gets a trailing ``per_axis ** ndim`` axis.

    ``spatial_moments=2`` on a 2-D mesh → ``(L+1, 2L+1, ng, *spatial, 4)``;
    BOTH moment factors coexist — the angular head queryable by type, the
    spatial tail the scheme's own MODAL axis (mass-weighted, labelled
    ``spatial_moment``; CS4c step 6 item 6.2c-iii) read by label — the
    orthogonal-axes invariant on the live carrier.
    """
    mesh = ld_2d
    L = 1
    field = HarmonicMomentFlux.zeros_for_problem_and_L(mesh, L, spatial_moments=2)
    expected = (L + 1, 2 * L + 1, mesh.ng, *mesh.spatial_shape, 2 ** mesh.ndim)
    np.testing.assert_equal(field.space.shape, expected)
    np.testing.assert_equal(field.values.shape, expected)
    np.testing.assert_equal(field.spatial_moments, 2)
    # both moment factors coexist (orthogonal axes)
    np.testing.assert_equal(head_layout(field.space.factors[0]).L, L)
    np.testing.assert_equal(BulkField.spatial_moments_per_axis_of(field.space), 2)
    assert field.space.axes is not None
    tail = [ax for ax in field.space.axes if ax.label == SPATIAL_MOMENT_AXIS_LABEL]
    np.testing.assert_equal(len(tail), 1)
    np.testing.assert_equal(tail[0], mesh.scheme.moment_axis(mesh.axes))


