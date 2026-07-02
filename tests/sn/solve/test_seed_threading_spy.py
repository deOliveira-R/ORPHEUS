r"""Mode-11 seed-threading spy — the fast-lane guard on SI's ``initial_guess``
contract (#226 taxonomy step-3 PREREQUISITE; verification spec §14).

**What is pinned (the contract, not the route).** Every source-iteration
inner step must receive the PREVIOUS ITERATE as its seed:
:math:`\psi_k = L^{-1}(q + \sum_i g_i\,\psi_{k-1})` with
``initial_guess = psi_{k-1}``.  The curvilinear Morel–Montry Carlson
coupled-pole closure reads the level's :math:`\psi(\mu=-1)` from this
argument (lesson L21), so a dropped / zeroed / stale seed is a
WRONG-FIXED-POINT bug on curvilinear meshes — not a rate change.

**Why a path-spy is the load-bearing fast net (spec §17 F4).**  A simulated
seed-drop on the het-2G sphere moves the eigenvalue by |Δk| ≈ 3.46e-2 — yet
under the canonical ``-m "not slow"`` run the strong value catcher
(``test_si_krylov_eigenvalue_equivalence_sphere``, GL-8 het-2G, marked
``@catches("M-SEED-DROP")``) is DESELECTED as ``@slow``; the only fast
reddening is a fragile 1G monotone margin, and cylinder's seeded-value gates
are structurally VACUOUS for seed-drops (the per-level α-dome telescopes —
the seed cancels exactly).  A fast sphere VALUE gate does not exist: the
config that activates the drop (GL-S16 / 40 cells) is slow, and a small
sphere either NaNs (LS-S4 / 16-cell fixed-source SI) or nulls the
sensitivity — a sub-floor value gate would be Mode-10 false confidence,
worse than none.  Hence: this spy pins the PATH in seconds (Mode 11); the
``@slow`` sphere equivalence gate stays the value catcher.

**Route-invariance across the step-3 driver rewire.**  The spy wraps
:meth:`~orpheus.sn.operators.streaming.InvertibleOperator.solve` — the one
surface both driver generations route through: pre-step-3
:class:`~orpheus.numerics.iteration.SourceIteration` calls
``resolvent.solve(rhs, initial_guess=...)`` (this method, on the 1-D
sphere's Jacobi path); post-step-3 the driver applies ``(L+C).inverse()``
(the invertible loss ``A``) and
:meth:`~orpheus.sn.operators.sweep_operator.SweepOperator.apply` DELEGATES
into ``inner.solve``.  This gate is GREEN before the rewire and MUST stay
green through it — it pins the contract, not the spelling.

**Copy-robustness.**  The assertion is values-equal
(:func:`numpy.testing.assert_array_equal`), NOT object identity — a
defensive copy of the iterate satisfies the contract; zeros, the rhs, or a
stale (cold) seed do not.

**Teeth (mutation-verified 2026-07-01, throwaway in-process patches —
spec §14.4 + §3a):** on the PRE-rewire driver (patching the then-extant
``SourceIteration._solve_with_seed``): M-SEED-DROP (thread no seed) REDs
the ``seed DROPPED`` guard at inner solve 1; M-SEED-ZERO (thread zeros
for ``psi_prev``) and M-SEED-STALE (thread the FIRST iterate every step)
RED the values-equal chain (at k = 1 / k ≥ 2 respectively).  On the
POST-rewire driver: M-PROBE (sever the ``initial_guess`` thread at
``SweepOperator.apply`` — observationally the driver dropping the kwarg)
REDs the same ``seed DROPPED`` guard, in under a second.  The ``@slow``
value companion M-SEED-VALUE (|Δk| ≈ 3.46e-2) stays on the sphere
equivalence gate, marked ``@catches("M-SEED-DROP")``.
"""
from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, CoordSystem
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.operators.streaming import InvertibleOperator
from orpheus.sn.solver import solve_sn_fixed_source
from tests.sn._test_helpers import curvilinear_two_region_mesh

pytestmark = [
    pytest.mark.foundation,
    pytest.mark.sentinel,
    pytest.mark.cap("solve"),
]

#: SI inner-iteration budget: ≥ 2 solves are needed to OBSERVE threading;
#: 6 keeps the run in fractions of a second (convergence is irrelevant —
#: ``inner_tol=0.0`` makes the loop run exactly this many steps).
_MAX_INNER = 6


def _small_het_2g_sphere():
    r"""~10-cell het-2G vacuum sphere on GL-8 — the smallest config where
    seed-threading is load-bearing (spec §14.2).

    Gauss–Legendre, NOT level-symmetric: LS-S4 on a small fixed-source
    sphere SI-diverges to NaN (the μ = −1 pole start); GL avoids the pole.
    Mixtures A (fuel) / B (moderator) give a genuinely heterogeneous 2G
    scattering ratio, so successive iterates differ and a stale seed is
    distinguishable from the fresh one.
    """
    materials = {2: get_mixture("A", "2g"), 0: get_mixture("B", "2g")}
    mesh = curvilinear_two_region_mesh(
        outers=(0.5, 1.0), mat_ids=(2, 0), n_cells=(5, 5),
        coord=CoordSystem.SPHERICAL, bc=BC.vacuum,
    )
    quad = Quadrature.gauss_legendre(8)
    return materials, mesh, quad


def test_si_threads_previous_iterate_as_seed(monkeypatch):
    r"""Every k-th inner solve's ``initial_guess`` IS the (k−1)-th iterate.

    Records ``(initial_guess, return)`` for every
    ``InvertibleOperator.solve`` call issued by a production fixed-source
    SI run, then walks the chain: call k's seed must equal call (k−1)'s
    return by VALUE.  Fails loudly (``pytest.fail`` — ``-O``-proof) on a
    dropped seed and on a chain that never iterates.
    """
    materials, mesh, quad = _small_het_2g_sphere()

    # (initial_guess, return) per inner solve — TimedFullField composites at
    # runtime; ``Any`` because the spy records whatever the duck-typed driver
    # threads (that looseness IS what the assertions then pin down).
    calls: list[tuple[Any, Any]] = []
    real_solve = InvertibleOperator.solve

    def spy(self, rhs, *, initial_guess=None):
        out = real_solve(self, rhs, initial_guess=initial_guess)
        calls.append((initial_guess, out))
        return out

    monkeypatch.setattr(InvertibleOperator, "solve", spy)

    ng, nx = 2, int(np.asarray(mesh.volumes).size)
    source = np.ones((quad.n_ordinates, ng, nx))
    solve_sn_fixed_source(
        materials, mesh, quad, source,
        inner_solver="source_iteration",
        max_inner=_MAX_INNER, inner_tol=0.0,
    )

    if len(calls) < 2:
        pytest.fail(
            f"SI issued only {len(calls)} inner solve(s) — no seed threading "
            f"to observe (expected {_MAX_INNER}; did the driver stop routing "
            f"through InvertibleOperator.solve?)"
        )
    for k in range(1, len(calls)):
        seed_k, _ = calls[k]
        if seed_k is None:
            pytest.fail(
                f"inner solve {k}: seed DROPPED (initial_guess=None) — the "
                f"driver no longer threads the previous iterate (M-SEED-DROP)"
            )
        prev_return = calls[k - 1][1]
        np.testing.assert_array_equal(
            seed_k.bulk.values, prev_return.bulk.values,
            err_msg=(
                f"inner solve {k}: seed ≠ previous iterate — zeros / rhs / "
                f"stale substituted for psi_(k-1) (M-SEED-ZERO / M-SEED-STALE)"
            ),
        )
