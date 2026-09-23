"""The structurally independent oracle for the 2-D DD sweep snapshot (#175).

Not collected (no ``test_`` prefix).  Three things live here, and nowhere
else:

* :func:`solver_2g_case` — the 6 x 4 two-material (fuel A | moderator B,
  2 groups) Lebedev-17 case that ``tests/gates/sn/operators/test_solver_components.py``
  hands to every test through its ``solver_2g`` fixture, and on which the
  frozen snapshot ``tests/gates/sn/sweep_ref_2g.npy`` is taken;
* :func:`snapshot_source` — the seeded source (legacy global RNG, seed 7)
  the snapshot is taken with;
* :func:`hand_sweep` — a from-scratch per-cell 2-D diamond-difference sweep
  written with PLAIN PYTHON LOOPS, the oracle.

Why the oracle is independent of the production sweep (X4, per axis):

* DERIVATION — it shares no code with the production path
  (``StreamingOperator.pose(problem) + MultiplicationOperator`` solved by
  ``(L + C).solve``, whose cell math is
  :meth:`~orpheus.transport.spatial.diamond.DiamondDifference.cell_kernel_batch`
  walked over an anti-diagonal wavefront).  The traversal order differs
  (per ordinate, row by row, column by column), the recurrence is
  transcribed independently from :eq:`dd-cartesian-2d`, and no helper is
  shared.
* INPUT — both sides read the same quadrature (``mu_x``, ``mu_y``,
  ``weights``), the same cell widths and the same ``sigma_t`` array.  Those
  are data, verified on their own rungs (the Lebedev quadrature tests), not
  the discrete model under test.

The discrete model is the same by construction: DD closure
``psi_out = 2 psi_c - psi_in``, vacuum inflow on every face, and the
per-ordinate source ``Q / W`` with ``W = sum(weights)`` (the producer-side
normalisation, lesson L18).

Regenerating the snapshot.  The snapshot is written ONLY through
:func:`regenerate_sweep_snapshot`, which refuses to write unless the
production sweep agrees with the oracle to the tolerance of
``TestTransportSweep::test_matches_independent_hand_sweep``::

    .venv/bin/python -O -c "from tests.gates.sn._independent_dd_sweep import regenerate_sweep_snapshot; regenerate_sweep_snapshot()"

History.  The oracle was written for #175 (2026-06-12) in the diagnostic
script ``derivations/diagnostics/diag_175_sweep_snapshot_regen.py``
(retired, R19), which no collected test ran; that script went unrunnable
twice as the production API moved under it (#347, and again when
``AngularBoundaryFlux.zeros_on`` retired).  Its standing readings, each
the cross-check on this case:

===========  ======================================  ==============  ==============
date         production leg                          ``max |dpsi|``  ``rel |dphi|``
===========  ======================================  ==============  ==============
2026-06-12   ``transport_sweep`` (operator-free)      ``3.5e-17``     ``9.8e-16``
2026-08-09   ``sweep_once``, ``(L + C).solve``        ``2.776e-17``   ``5.152e-16``
2026-09-22   ``sweep_once``, ``(L + C).solve``        ``8.3e-17``     ``1.39e-15``
===========  ======================================  ==============  ==============
"""

from __future__ import annotations

import numpy as np

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.problem import SNProblem
from orpheus.sn.solver import SNSolver
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.source_sinks import AngularSourceSink

from tests.gates.sn._test_helpers import SN_TESTS_ROOT, sweep_once

#: The seed of the snapshot's source.  Legacy global RNG, so the bytes of
#: ``sweep_ref_2g.npy`` reproduce.
SNAPSHOT_SEED = 7

SNAPSHOT_PATH = SN_TESTS_ROOT / "sweep_ref_2g.npy"


def solver_2g_case() -> tuple[SNSolver, dict, SNProblem, Quadrature]:
    """The 6 x 4 fuel|moderator 2-group Lebedev-17 case.

    Returns ``(solver, materials, problem, quadrature)``: x-columns 0..2
    are fuel (mixture A), 3..5 moderator (mixture B), square cells of side
    0.2 cm.
    """
    fuel = get_mixture("A", "2g")
    mod = get_mixture("B", "2g")
    materials = {2: fuel, 0: mod}

    nx, ny, delta = 6, 4, 0.2
    mat = np.zeros((nx, ny), dtype=int)
    mat[:3, :] = 2
    mat[3:, :] = 0
    mesh = Mesh2D(
        edges_x=np.linspace(0, nx * delta, nx + 1),
        edges_y=np.linspace(0, ny * delta, ny + 1),
        mat_map=mat,
    )
    quad = Quadrature.lebedev(order=17)
    problem = SNProblem(mesh, quad, materials)
    return SNSolver(problem), materials, problem, quad


def snapshot_source(ng: int, spatial_shape: tuple[int, ...]) -> np.ndarray:
    """The isotropic volumetric source ``Q`` (``(ng, nx, ny)``, strictly
    positive) the snapshot is taken with."""
    np.random.seed(SNAPSHOT_SEED)
    return np.random.rand(ng, *spatial_shape) + 0.01


def hand_sweep(
    problem: SNProblem, sig_t: np.ndarray, Q: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Per-cell-loop 2-D DD sweep, vacuum inflow on every face.

    Solves ``(Omega . grad + sigma_t) psi_n = Q / W`` per ordinate with the DD
    closure.  Returns ``(psi, phi)``: the angular flux ``(N, ng, nx, ny)``
    and the scalar flux ``phi = sum_n w_n psi_n``, ``(ng, nx, ny)``.
    """
    quad = problem.quad
    N = quad.N
    ng, nx, ny = Q.shape
    W = float(quad.weights.sum())
    dx, dy = (np.asarray(w) for w in problem.axis_widths)

    qn = Q / W  # per-ordinate source magnitude, identical across n
    psi = np.zeros((N, ng, nx, ny))

    for n in range(N):
        mu = float(quad.mu_x[n])
        eta = float(quad.mu_y[n])
        xs = range(nx) if mu >= 0 else range(nx - 1, -1, -1)
        ys = range(ny) if eta >= 0 else range(ny - 1, -1, -1)
        for g in range(ng):
            # Vacuum: zero inflow at the upstream domain faces.
            inflow_y = np.zeros(nx)  # psi on the upstream y-face, per x-column
            for j in ys:
                inflow_x = 0.0  # psi on the upstream x-face of this row
                for i in xs:
                    cx = 2.0 * abs(mu) / dx[i]
                    cy = 2.0 * abs(eta) / dy[j]
                    psi_c = (qn[g, i, j] + cx * inflow_x + cy * inflow_y[i]) / (
                        sig_t[g, i, j] + cx + cy
                    )
                    psi[n, g, i, j] = psi_c
                    inflow_x = 2.0 * psi_c - inflow_x
                    inflow_y[i] = 2.0 * psi_c - inflow_y[i]

    phi = np.einsum("n,ngxy->gxy", quad.weights, psi)
    return psi, phi


def production_sweep(
    problem: SNProblem, sig_t: np.ndarray, Q: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """The production leg: one ``(L + C).solve`` sweep, vacuum inflow.

    The same helper (:func:`tests.gates.sn._test_helpers.sweep_once`) the
    snapshot test calls, so the snapshot is written from the path it is
    compared against.
    """
    ang, phi = sweep_once(
        AngularSourceSink.from_isotropic(Q, problem),
        sig_t,
        problem,
        AngularBoundaryFlux.zeros(problem.angular_trace),
    )
    return np.asarray(ang), np.asarray(phi)


def sweep_agreement_tolerances(problem: SNProblem) -> tuple[float, float]:
    """``(tol_psi, tol_phi)``: the relative tolerances of production vs oracle.

    Both sides evaluate the same recurrence in a different order, so the
    gap is floating-point re-association only.  A cell's angular flux is
    the end of a dependency chain of at most ``nx + ny`` cell updates, and
    the scalar flux adds a reduction of depth ``N`` over the ordinates, so
    the bounds are ``SAFETY * depth * eps`` with ``SAFETY = 10``:
    ``tol_psi = 10 (nx + ny) eps`` (normalised by ``max |psi|``) and
    ``tol_phi = 10 (nx + ny + N) eps`` (per cell, relative).
    `[M]` 2026-09-22 on :func:`solver_2g_case` (``nx + ny = 10``,
    ``N = 110``): the gap is ``max |dpsi| / max |psi| = 5.3e-16`` and
    ``max rel |dphi| = 1.39e-15``, against bounds of ``2.2e-14`` and
    ``2.7e-13``.
    """
    eps = float(np.finfo(float).eps)
    nx, ny = problem.spatial_shape
    depth = nx + ny
    return 10.0 * depth * eps, 10.0 * (depth + problem.quad.N) * eps


def regenerate_sweep_snapshot() -> None:
    """Write ``sweep_ref_2g.npy`` from the production sweep, only if the
    production sweep agrees with :func:`hand_sweep`.

    Raises :class:`RuntimeError`, writing nothing, when it does not: a
    snapshot regenerated from a production sweep the oracle disagrees with
    would freeze the defect as the reference.
    """
    solver, _, problem, _ = solver_2g_case()
    sig_t = solver.problem.mat_xs.total_cross_section
    Q = snapshot_source(solver.ng, problem.spatial_shape)
    psi_prod, phi_prod = production_sweep(problem, sig_t, Q)
    psi_hand, phi_hand = hand_sweep(problem, sig_t, Q)
    tol_psi, tol_phi = sweep_agreement_tolerances(problem)
    dpsi = float(np.abs(psi_prod - psi_hand).max() / np.abs(psi_hand).max())
    dphi = float((np.abs(phi_prod - phi_hand) / np.abs(phi_hand)).max())
    if dpsi > tol_psi or dphi > tol_phi:
        raise RuntimeError(
            "refusing to write the sweep snapshot: the production sweep "
            f"disagrees with the independent hand sweep (max |dpsi|/max|psi| "
            f"= {dpsi:.3e} vs {tol_psi:.3e}; max rel |dphi| = {dphi:.3e} vs "
            f"{tol_phi:.3e})."
        )
    np.save(SNAPSHOT_PATH, phi_prod)
    print(f"wrote {SNAPSHOT_PATH} (max |dpsi|/max|psi| = {dpsi:.3e}, "
          f"max rel |dphi| = {dphi:.3e})")
