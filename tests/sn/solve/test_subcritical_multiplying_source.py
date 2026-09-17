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

import warnings

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry.mesh import BC, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import (
    SupercriticalSourceProblem, _as_sn_mesh, solve_sn, solve_sn_fixed_source, solve_sn_multiplying_source,
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
    k = float(solve_sn(mats, mesh, quad).outcome.keff)
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


# ═══════════════════════════════════════════════════════════════════════
# ERR-086 — the multiplying entry's SILENCE (step 3, 2026-09-17)
# ═══════════════════════════════════════════════════════════════════════
#
# Until step 3 U2 `solve_sn_multiplying_source` returned the SI arm's Solution
# directly and never reached the hoisted `warn_if_unconverged` /
# `warn_if_gauge_freedom` its four siblings emit from their public entry: a
# truncated solve and a gauge-fixed trace went by in silence (the step-3
# anchors recorded both, beside the sibling's warning on the same hub), and
# the warn-site count gate `len(sites) == 7` pinned the inventory that excluded
# it.  The two rows below are the catchers: each pairs the fifth entry with its
# sibling on ONE hub, so "no warning" cannot be mistaken for "nothing to warn".


def _dilute_fissile(ng: int = 2):
    """Fissile enough to be ADMITTED (k_eff > 0), dilute enough to stay SUBCRITICAL on
    an all-reflective box — a library mixture there is supercritical (k∞ = 1.875 for A)
    and the entry would refuse."""
    from orpheus.derivations.common.xs_library import make_mixture
    sig_t = np.linspace(0.8, 1.6, ng)
    sig_f = 0.05 * np.ones(ng)
    return make_mixture(
        sig_t=sig_t, sig_c=sig_t - sig_f, sig_f=sig_f, nu=2.4 * np.ones(ng),
        chi=np.array([1.0] + [0.0] * (ng - 1)), sig_s=np.zeros((ng, ng)),
    )


def _gauge_singular_box():
    """The entry ledger's all-reflective (3, 4) box — ODD first axis, ≥ 2 reflective
    axis pairs: the loss operator is exactly singular and a uniform isotropic source
    EXCITES the kernel (``[M]`` the ledger gate's parity table)."""
    from orpheus.geometry import Mesh2D
    quadrature = Quadrature.level_symmetric(sn_order=4)
    reflective = BC("reflective")
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 1.0, 4), edges_y=np.linspace(0.0, 2.0, 5),
        mat_map=np.zeros((3, 4), dtype=int),
        bc_xmin=reflective, bc_xmax=reflective, bc_ymin=reflective, bc_ymax=reflective,
    )
    source = np.full((quadrature.weights.size, 2, 3, 4), 1.0 / float(quadrature.weights.sum()))
    return {0: _dilute_fissile()}, mesh, quadrature, source


@pytest.mark.catches("ERR-086")
def test_a_truncated_multiplying_solve_is_AUDIBLE_like_its_sibling() -> None:
    """A truncated solve announces itself ONCE, from the PUBLIC entry, blaming the
    caller — exactly as the sibling entry does on the same hub."""
    from orpheus.numerics.convergence import ConvergenceWarning
    mats, mesh, quad = _slab(4.0)
    q = _uniform_source(quad, 2, 8)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        mult = solve_sn_multiplying_source(mats, mesh, quad, q, inner_tol=1e-12, max_inner=5)
    convergence = [w for w in caught if issubclass(w.category, ConvergenceWarning)]
    _require(len(convergence) == 1, f"the multiplying entry must warn exactly once on a truncated solve; got {len(convergence)}")
    _require(convergence[0].filename == __file__, f"the warning must blame the CALLER, not {convergence[0].filename}")
    _require(not mult.converged(), "non-vacuity: the solve was truncated")
    with warnings.catch_warnings(record=True) as sibling:
        warnings.simplefilter("always")
        solve_sn_fixed_source(mats, mesh, quad, q, inner_tol=1e-12, max_inner=5)
    _require(sum(issubclass(w.category, ConvergenceWarning) for w in sibling) == 1, "the sibling warns once on the same hub — the pairing that makes this row non-vacuous")


@pytest.mark.catches("ERR-086")
def test_a_gauge_singular_multiplying_solve_is_AUDIBLE_like_its_sibling() -> None:
    """On a gauge-singular subcritical fissile box the repair FIRES on the fifth entry
    and, since step 3, is SAID — the same ``GaugeFreedomWarning`` its sibling emits,
    and the same recorded displacement (``[M]`` 6.08e-02 of the trace on both)."""
    from orpheus.numerics.outcome import Measured
    from orpheus.sn.operators.loss_kernel_gauge import GaugeFreedomWarning, gauge_freedom
    mats, mesh, quad, source = _gauge_singular_box()
    hub = _as_sn_mesh(mesh, quad, mats, None)
    _require(gauge_freedom(hub).present, "non-vacuity: the hub must be gauge-singular")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        mult = solve_sn_multiplying_source(
            mats, mesh, quad, source, boundary_condition=None,
            inner_schedule="gauss_seidel", inner_tol=1e-13, max_inner=400_000,
        )
    with warnings.catch_warnings(record=True) as sibling_caught:
        warnings.simplefilter("always")
        pure = solve_sn_fixed_source(
            mats, mesh, quad, source, boundary_condition=None, inner_solver="source_iteration",
            inner_schedule="gauss_seidel", inner_tol=1e-13, max_inner=400_000,
        )
    _require(any(issubclass(w.category, GaugeFreedomWarning) for w in caught), "the multiplying entry must say the trace was gauge-fixed")
    _require(any(issubclass(w.category, GaugeFreedomWarning) for w in sibling_caught), "…as its sibling does on the same hub (the pairing)")
    mg, pg = mult.certificate.gauge, pure.certificate.gauge
    _require(isinstance(mg, Measured) and isinstance(pg, Measured), f"both certificates record the displacement; got {mg!r} / {pg!r}")
    _require(mg.value > 1e-3 and abs(mg.value - pg.value) <= 1e-9, f"the same kernel component on both entries: {mg.value:.6e} vs {pg.value:.6e}")
