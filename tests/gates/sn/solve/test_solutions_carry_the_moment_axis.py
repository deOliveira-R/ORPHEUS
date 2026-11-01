r"""#246 P3 / P4 — an SN solution's angular flux carries the spatial-moment axis exactly when its closure is LD.

Moved from ``tests/gates/transport/spatial/test_moment_axis_predicates.py`` by layer DAG Phase 4 (2026-11-01): the two
rows read the output of ``solve_sn_fixed_source`` (the moment-aware output wrap of the SN entry), an SN claim. The
transport file keeps the scheme-level predicate rows and the construct-general pin P4' (a hand-built width-1 field on
an LD closure's bulk carries no moment axis), which is why the inner walk reads the SCHEME, not the field.

``foundation`` (a typed query over a data structure; no ``verifies``). ``-O``-safe: ``pytest.fail``, never a bare
``assert``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import StructuredGeometry
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn import solve_sn_fixed_source
from orpheus.transport.boundary_conditions import WallCondition
from orpheus.transport.spatial import DiamondDifference, LinearDiscontinuous
from tests.gates._exterior import conditions_for

#: Exterior conditions of the tags this file's meshes carry.
_CONDITIONS = {0: WallCondition.vacuum()}

pytestmark = pytest.mark.foundation


def _require(predicate: bool, msg: str) -> None:
    if not predicate:
        pytest.fail(msg)


def _solve_1g_slab(scheme):
    r"""A minimal 1G vacuum slab solve — returns the result whose
    ``angular_flux`` carries (LD) or lacks (DD) the spatial-moment factor."""
    from orpheus.derivations.continuous.mms.sn import _make_1g_mixture

    nx = 8
    materials = {0: _make_1g_mixture(1.0, 0.5)}
    mesh = Mesher(StructuredGeometry.from_homogeneous(1.0, 0)).partition(
        CellsByCount.uniform_width(nx),
    ).mesh
    quad = Quadrature.gauss_legendre(8)
    Q = np.ones((quad.n_ordinates, 1, nx)) / quad.weights.sum()
    return solve_sn_fixed_source(
        materials, mesh, quad, Q, scheme=scheme,
        inner_solver="krylov",
        max_inner=500, inner_tol=1e-10, boundary_conditions=conditions_for(mesh, _CONDITIONS),
    )


def test_spatial_moment_axis_present_on_moment_aware_producer_output() -> None:
    r"""P3 (positive): a field produced THROUGH the LD moment-aware output wrap
    carries the spatial-moment factor → ``spatial_moments_per_axis > 1``."""
    result = _solve_1g_slab(LinearDiscontinuous())
    bulk = result.angular_flux.interior
    _require(
        bulk.spatial_moments_per_axis > 1,
        "LD producer output should carry the spatial-moment axis",
    )


def test_spatial_moment_axis_absent_on_dd_field() -> None:
    r"""P4 (negative): a DD output field lacks the factor → ``spatial_moments_per_axis == 1``."""
    result = _solve_1g_slab(DiamondDifference())
    bulk = result.angular_flux.interior
    _require(
        bulk.spatial_moments_per_axis == 1,
        "DD output should NOT carry the spatial-moment axis",
    )
