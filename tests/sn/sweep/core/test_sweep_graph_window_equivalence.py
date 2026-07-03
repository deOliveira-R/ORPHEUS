r"""``window ≡ full-field`` equivalence — the Phase 5b storage-B oracle.

The rolling :math:`(d{-}1)`-frontier moving-frontier window walks
(:meth:`SweepDependencyGraph.walk_windowed`, the storage-B PRODUCTION
path / the ``MovingFrontierWindow`` strategy) MUST reproduce the full per-axis
face-field walk (:meth:`~SweepDependencyGraph.walk_full`, the storage-A
path, KEPT as the
verification oracle) **bit-for-bit**.

Both share the storage-free cell kernel
(:meth:`DiamondDifference.cell_kernel_batch` /
:meth:`~DiamondDifference.residual_kernel_batch`), so the cell *math* cannot
drift between them — only the storage walk (full field vs window) differs,
which is exactly what this test pins.

Dimension coverage (the sweep-strategy carve, phase S4)
-----------------------------------------------------
The window walk is generalised from the hardcoded 2-diagonal to the general
``frontier_dim = d − 1`` rolling slab.  This file pins ``window ≡ full`` at
EVERY ``d`` the spine admits:

* **d = 1** — the degenerate ``frontier_dim == 0`` point (the window is not the
  production default at d=1 — the cumprod scan is — but its d=1 *capability* is
  the governing principle "construct general, select narrow", verified here).
* **d = 2** — the production default; the ~3× peak-memory + ~0.77× contiguity
  win.  This is the BIT-IDENTITY anchor: the d=2 read selector degenerates to a
  contiguous slice, so the generalised walk reproduces the legacy 2-diagonal
  window byte-for-byte (a regression here means the carve changed d=2 behaviour).
* **d = 3** — the SYNTHETIC admission (no 3-D quadrature — arbitrary
  ``N_oct`` / streaming / source, the B7 ``test_sweep_graph_nd_admission``
  idiom).  Proves the ``frontier_dim`` construction is genuinely general; the
  d=3 *speed* (the contiguity win at the simplex surface) is a SEPARATE
  profiling question, deliberately OUT of this correctness gate.

``foundation`` level: a software invariant (window ≡ full-field) with no
theory ``:label:`` — it pins the storage refactor, not a physics equation.

``-O`` discipline (vv Mode 8): every assertion is a ``np.testing.*`` /
``pytest.fail`` FUNCTION CALL, never a bare ``assert`` — so it fires under the
canonical ``python -O`` invocation.
"""

from __future__ import annotations

import itertools

import numpy as np
import pytest

from orpheus.sn.spatial.diamond import DiamondDifference
from orpheus.sn.loss_representation.sweep_graph import (
    OctantLabel,
    SweepDependencyGraph,
    _CellResidual,
    _CellSolveAngular,
)

pytestmark = pytest.mark.foundation


# (shape, N_oct, ng) — d ∈ {1, 2, 3}.  d=2 carries the non-square / asymmetric
# meshes that catch an x↔y axis swap; d=3 is synthetic shape admission.
SHAPES = [
    ((6,), 4, 2),                  # d=1: the frontier_dim == 0 point
    ((10,), 3, 3),                 # d=1
    ((8, 8), 4, 2),                # d=2: the production default (bit-identity)
    ((12, 7), 6, 3),               # d=2 non-square (axis-swap moat)
    ((5, 9), 4, 2),                # d=2 non-square
    ((16, 16), 8, 2),              # d=2
    ((3, 2, 3), 3, 2),             # d=3 synthetic admission
    ((4, 3, 2), 4, 2),             # d=3 synthetic (asymmetric)
]


def _octants(d):
    """The ``2^d`` sign signatures."""
    return list(itertools.product((-1, +1), repeat=d))


def _inputs(shape, N_oct, ng, seed):
    """Random per-axis sweep inputs (NO real quadrature — synthetic admission).

    ``inflow[a]`` is the domain inflow face perpendicular to axis ``a``, indexed
    by the OTHER spatial coordinates: shape ``(N_oct, ng, *shape∖a)``.
    """
    rng = np.random.default_rng(seed)
    d = len(shape)
    perp = lambda a: tuple(shape[b] for b in range(d) if b != a)
    return dict(
        sig_t=rng.uniform(0.3, 3.0, size=(ng, *shape)),
        Q=rng.uniform(0.0, 2.0, size=(N_oct, ng, *shape)),
        str_axes=tuple(rng.uniform(0.5, 4.0, size=(N_oct, shape[a])) for a in range(d)),
        inflow=tuple(rng.uniform(0.0, 1.0, size=(N_oct, ng, *perp(a))) for a in range(d)),
        weights=rng.uniform(0.1, 1.0, size=(N_oct,)),
        probe=rng.uniform(-1.0, 1.0, size=(N_oct, ng, *shape)),
    )


def _edge_index(d, a, pos):
    """Index into a per-axis face buffer: face position ``pos`` along axis ``a``,
    full on the others — ``(:, :, …, pos@a, …)``."""
    return (slice(None), slice(None)) + tuple(
        pos if b == a else slice(None) for b in range(d)
    )


def _seed_faces(graph, shape, N_oct, ng, inp):
    """Allocate the d full-field per-axis face buffers and ι_*-seed the domain
    inflow at each axis's incoming edge (face 0 for sign ≥ 0, face ``n_a`` else)."""
    d = len(shape)
    signs = graph.label.signs
    faces = []
    for a in range(d):
        fshape = list(shape)
        fshape[a] += 1
        buf = np.zeros((N_oct, ng, *fshape))
        in_pos = 0 if signs[a] >= 0 else shape[a]
        buf[_edge_index(d, a, in_pos)] = inp["inflow"][a]
        faces.append(buf)
    return tuple(faces)


def _shed_faces(graph, shape, faces):
    """ι*-read each axis's domain outflow off the full-field face buffers (the
    opposite edge: face ``n_a`` for sign ≥ 0, face 0 else)."""
    d = len(shape)
    signs = graph.label.signs
    return tuple(
        faces[a][_edge_index(d, a, shape[a] if signs[a] >= 0 else 0)]
        for a in range(d)
    )


def _full_field_apply(graph, shape, N_oct, ng, inp):
    """Oracle: the production full-field solve walk on dense face buffers."""
    faces = _seed_faces(graph, shape, N_oct, ng, inp)
    ang = np.zeros((N_oct, ng, *shape))
    scal = np.zeros((ng, *shape))
    graph.walk_full(
        level_op=_CellSolveAngular(
            scheme=DiamondDifference(),
            weights_octant=inp["weights"],
            angular_flux_octant=ang,
            scalar_flux_buf=scal,
        ),
        psi_faces_octant=faces,
        Q_octant=inp["Q"],
        sig_t=inp["sig_t"],
        str_axes_octant=inp["str_axes"],
    )
    return ang, scal, _shed_faces(graph, shape, faces)


def _windowed_apply(graph, shape, N_oct, ng, inp):
    ang = np.zeros((N_oct, ng, *shape))
    scal = np.zeros((ng, *shape))
    cap = tuple(
        np.zeros((N_oct, ng, *[shape[b] for b in range(len(shape)) if b != a]))
        for a in range(len(shape))
    )
    graph.walk_windowed(
        level_op=_CellSolveAngular(
            scheme=DiamondDifference(),
            weights_octant=inp["weights"],
            angular_flux_octant=ang,
            scalar_flux_buf=scal,
        ),
        inflow=inp["inflow"],
        Q_octant=inp["Q"],
        sig_t=inp["sig_t"],
        str_axes_octant=inp["str_axes"],
        capture=cap,
    )
    return ang, scal, cap


def _full_field_residual(graph, shape, N_oct, ng, inp):
    faces = _seed_faces(graph, shape, N_oct, ng, inp)
    res = np.zeros((N_oct, ng, *shape))
    graph.walk_full(
        level_op=_CellResidual(
            scheme=DiamondDifference(),
            psi_avg_probe_octant=inp["probe"],
            residual_octant=res,
        ),
        psi_faces_octant=faces,
        Q_octant=inp["Q"],
        sig_t=inp["sig_t"],
        str_axes_octant=inp["str_axes"],
    )
    return res, _shed_faces(graph, shape, faces)


def _windowed_residual(graph, shape, N_oct, ng, inp):
    res = np.zeros((N_oct, ng, *shape))
    cap = tuple(
        np.zeros((N_oct, ng, *[shape[b] for b in range(len(shape)) if b != a]))
        for a in range(len(shape))
    )
    graph.walk_windowed(
        level_op=_CellResidual(
            scheme=DiamondDifference(),
            psi_avg_probe_octant=inp["probe"],
            residual_octant=res,
        ),
        inflow=inp["inflow"],
        Q_octant=inp["Q"],
        sig_t=inp["sig_t"],
        str_axes_octant=inp["str_axes"],
        capture=cap,
    )
    return res, cap


@pytest.mark.parametrize("shape,N_oct,ng", SHAPES)
def test_solve_window_equals_full_field(shape, N_oct, ng):
    """walk_windowed ≡ walk_full (solve direction), bit-for-bit, at every
    octant and every ``d``.

    d=2 is the bit-identity anchor; d=1 / d=3 are the synthetic ``frontier_dim``
    admission (the S4 generalisation correctness gate)."""
    d = len(shape)
    for signs in _octants(d):
        graph = SweepDependencyGraph.from_cartesian(shape, label=OctantLabel(signs))
        inp = _inputs(shape, N_oct, ng, seed=hash((shape, signs)) % (2**32))
        ang_f, scal_f, out_f = _full_field_apply(graph, shape, N_oct, ng, inp)
        ang_w, scal_w, cap_w = _windowed_apply(graph, shape, N_oct, ng, inp)
        np.testing.assert_array_equal(
            ang_w, ang_f, err_msg=f"shape={shape} signs={signs}: angular flux")
        np.testing.assert_array_equal(
            scal_w, scal_f, err_msg=f"shape={shape} signs={signs}: scalar flux")
        for a in range(d):
            np.testing.assert_array_equal(
                cap_w[a], out_f[a],
                err_msg=f"shape={shape} signs={signs}: axis-{a} outflow shed")


@pytest.mark.parametrize("shape,N_oct,ng", SHAPES)
def test_residual_window_equals_full_field(shape, N_oct, ng):
    """walk_windowed ≡ walk_full (apply direction — the matvec twin),
    bit-for-bit, every ``d``."""
    d = len(shape)
    for signs in _octants(d):
        graph = SweepDependencyGraph.from_cartesian(shape, label=OctantLabel(signs))
        inp = _inputs(shape, N_oct, ng, seed=(hash((shape, signs)) % (2**32)) ^ 7)
        res_f, out_f = _full_field_residual(graph, shape, N_oct, ng, inp)
        res_w, cap_w = _windowed_residual(graph, shape, N_oct, ng, inp)
        np.testing.assert_array_equal(
            res_w, res_f, err_msg=f"shape={shape} signs={signs}: operator residual")
        for a in range(d):
            np.testing.assert_array_equal(
                cap_w[a], out_f[a],
                err_msg=f"shape={shape} signs={signs}: axis-{a} outflow shed")


@pytest.mark.parametrize("shape", [(16, 16), (16, 64), (3, 2, 3), (4, 3, 5)])
def test_window_backing_is_the_d_minus_1_slab(shape):
    """The window plan's free box is the ``(d−1)``-slab ``shape[:-1]`` — the
    interior backing is ``O(∏_{a<d−1} n_a)``, INDEPENDENT of the determined
    (last) axis ``n_{d−1}``.  The storage-B peak-memory win, generalised: the
    full field is ``O(∏ n_a)``, the window drops the determined-axis factor."""
    graph = SweepDependencyGraph.from_cartesian(
        shape, label=OctantLabel((1,) * len(shape)))
    plan = graph.window_plan
    np.testing.assert_array_equal(
        plan.free_bbox, shape[:-1],
        err_msg=f"shape={shape}: free box {plan.free_bbox} != shape[:-1]={shape[:-1]}")
    np.testing.assert_array_equal(
        plan.det, len(shape) - 1,
        err_msg=f"shape={shape}: determined axis != last")
