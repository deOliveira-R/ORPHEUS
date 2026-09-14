r"""The subcritical multiplying source — the ``(M, q)`` cell's witness (consumers
campaign step 2, C3b-2; plan §5.1; test-architect delta §B5).

``(A − F)ψ = q`` is ``SourcePosing(hub.pencil.at(1), q)`` — the pencil's member at
the physical σ = 1 handed to the fixed-source Strategy with the production LAGGED
as one more explicit gain.  Well posed iff the medium is SUBCRITICAL
(:math:`\rho(A^{-1}F) = k_{\rm eff} < 1`); the driver certifies it with the hub's
k-solve and REFUSES otherwise (RULED 2026-09-13 fork 3 (a)).  ``[M]`` fixtures
(test-architect §B5, 2-group fuel|moderator slab, GL-8, 4+4 cells): ``L = 2.0``
reflective|vacuum k = 0.435195214 (the anchors' own ``_slab_hub``), ``L = 4.0``
k = 0.907457573 (1/(1−k) = 10.8 — the STRONG discriminator), ``L = 8.0``
reflective|reflective k = 1.374233987 (the REFUSAL leg).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry.mesh import BC, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import (
    SupercriticalSourceProblem, solve_sn, solve_sn_fixed_source, solve_sn_multiplying_source,
)

pytestmark = pytest.mark.l1


def _require(cond: object, msg: str) -> None:
    if not cond:
        raise AssertionError(msg)


def _slab(L: float, bc_right: str = "vacuum"):
    mats = {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")}
    mesh = Mesh1D(edges=np.linspace(0.0, L, 9), mat_ids=np.array([0] * 4 + [1] * 4, dtype=int),
                  bc_left=BC("reflective"), bc_right=BC(bc_right))
    return mats, mesh, Quadrature.gauss_legendre(8)


def _uniform_source(quad, ng, nx):
    return np.ones((quad.N, ng, nx))


@pytest.mark.parametrize("L,k_ref,budget", [(2.0, 0.435195214, None), (4.0, 0.907457573, 6000)],
                         ids=["slab_sub", "slab_near_critical"])
def test_the_multiplying_solve_converges_and_its_certificate_holds(L: float, k_ref: float, budget) -> None:
    """The lagged-fission iteration converges at a rate governed by k: ``[M]`` the
    near-critical slab needs 4300 inner iterations at ``inner_tol = 1e-12``
    (the default budget derived from the tolerance alone is 1961 — a starved
    inner is NOT a refusal, it is a budget; #340), the subcritical one fits the
    default."""
    mats, mesh, quad = _slab(L)
    k = float(solve_sn(mats, mesh, quad).keff)
    _require(abs(k - k_ref) < 1e-6, f"the fixture's k moved: {k!r} vs {k_ref!r}")
    sol = solve_sn_multiplying_source(mats, mesh, quad, _uniform_source(quad, 2, 8), inner_tol=1e-12, max_inner=budget)
    _require(sol.converged(), "the multiplying source solve must converge (k < 1)")
    # The balance ⟨1, (A − F)ψ − q⟩ is enforced by the driver's convergence
    # CERTIFICATE (`_certify_within_group_exit`, posed on the SOLVED equation:
    # the lagged production re-enters the certified rhs) — it RAISES on a
    # defect, so reaching this line IS the closure.  `[M]` the fixed-source
    # path records no `balance_defect` (eigenvalue-only), so an assertion on
    # it here would be inert (the archivist's C3b-2 review).


def test_the_multiplying_solution_exceeds_the_pure_transport_one_cellwise() -> None:
    """Krein–Rutman / Neumann: fission adds neutrons — ψ_mult > ψ_pure componentwise on the bulk."""
    mats, mesh, quad = _slab(4.0)
    q = _uniform_source(quad, 2, 8)
    pure = solve_sn_fixed_source(mats, mesh, quad, q, inner_tol=1e-12)
    mult = solve_sn_multiplying_source(mats, mesh, quad, q, inner_tol=1e-12)
    phi_p = np.asarray(getattr(pure.scalar_flux, "values", pure.scalar_flux), dtype=float)
    phi_m = np.asarray(getattr(mult.scalar_flux, "values", mult.scalar_flux), dtype=float)
    _require(np.all(phi_m > phi_p), "the multiplying solution must exceed the pure-transport one cell-wise")
    _require(phi_m.max() / phi_p.max() > 3.0, "near-critical: 1/(1−k) ≈ 10.8 amplifies the fission term strongly")


def test_a_supercritical_hub_is_REFUSED_at_the_driver() -> None:
    mats, mesh, quad = _slab(8.0, bc_right="reflective")
    with pytest.raises(SupercriticalSourceProblem, match=r"k_eff = 1\.37"):
        solve_sn_multiplying_source(mats, mesh, quad, _uniform_source(quad, 2, 8))


def test_the_composition_is_the_loss_minus_the_production() -> None:
    """``SourcePosing(pencil.at(1), q).operator.apply(x) == loss.apply(x) − production.apply(x)``
    bit-identically on a typed coupled state."""
    from orpheus.sn.mesh.augmented_mesh import SNMesh
    mats, mesh, quad = _slab(2.0)
    hub = SNMesh(mesh, quad, mats)
    rec = hub.system
    rng = np.random.default_rng(0)
    x = rec.space.zeros()
    for member in x.systems:  # a seeded coupled state: every member's interior and boundary
        for part in (member.interior, member.boundary):
            values = np.asarray(part.values)
            values[...] = rng.random(values.shape)
    posing = hub.source_posing(rec.space.zeros())
    lhs = np.asarray(posing.operator.apply(x).to_flat(), dtype=float)
    ref = np.asarray((rec.loss.apply(x) - rec.production.apply(x)).to_flat(), dtype=float)
    _require(np.any(lhs != 0.0), "POSITIVE CONTROL: the seeded state must be non-trivial")
    _require(np.array_equal(lhs, ref), "the (M, q) operator is the loss minus the production, bit-identically")


def test_the_zero_D_multiplying_flux_matches_the_closed_form() -> None:
    """REFERENCE (closed form): on the manufactured subcritical 0-D mixture
    (``SigP``/``SigF`` × 0.4 → k_inf = 0.75 exactly) ``(A − F)⁻¹·1 = [60, 70]``
    strictly positive; ``[M]`` scale 0.6 → k_inf = 1.125 and the sign FLIPS — the
    0-D control of the refusal leg."""
    import dataclasses
    from orpheus.homogeneous.solver import HomogeneousProblem
    from orpheus.numerics.pencil import OperatorPencil
    base = get_mixture("A", "2g")
    for scale, positive in ((0.4, True), (0.6, False)):
        sub = dataclasses.replace(base, SigP=base.SigP * scale, SigF=base.SigF * scale)
        problem = HomogeneousProblem(sub)
        pencil = OperatorPencil(problem.loss, problem.production)
        A = np.asarray(pencil.at(1.0).as_matrix(), dtype=float)
        phi = np.linalg.solve(A, np.ones(A.shape[0]))
        _require(bool(np.all(phi > 0)) is positive, f"scale {scale}: (A − F)⁻¹·1 = {phi} — positivity {positive}")
        if positive:
            _require(np.allclose(phi, [60.0, 70.0], rtol=1e-12), f"the closed form [60, 70]; got {phi}")
