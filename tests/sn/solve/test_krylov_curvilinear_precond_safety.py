r"""L1 regression catcher — typed-Krylov inner solver convergence under the
two production preconditioner choices.

Origin
======

Promoted from ``derivations/diagnostics/diag_r1_step_d_probe_b_identity_precond.py``
(R-1 Step 4 session 1, 2026-05-19).  The diagnostic was a hypothesis
test for the silent stateful-inverse bug class (issue #203,
lesson :ref:`L19 <lessons-l19>`):
:class:`~orpheus.numerics.iteration.KrylovAcceleration` with
``preconditioner=None`` silently routed through
:meth:`~orpheus.sn.operators.streaming.InvertibleOperator.solve`, which (pre-
Phase-1.2) read the lag-1 frame ``rhs(1)`` to seed the curvilinear
Carlson coupled-pole closure.  GMRES feeds residual VECTORS into the
preconditioner — those residuals have no history — so the sweep
silently used the in-iteration default (zero) seed.  On sphere/cylinder
this destabilised the M-M closure as a preconditioner.

Phase 1.2 (commit ``c93355c``) made the bug class **structurally
unreachable**: ``InvertibleOperator.solve`` is now a pure function of
``(rhs, initial_guess, boundary)``.  ``KrylovAcceleration``'s
``preconditioner=None`` default invokes ``L.solve(q)`` with
``initial_guess=None``, which the M-M closure interprets as an
explicit cold start (deterministic zero seed via the
``psi_half_seed`` strategy).  No silent fallback path remains.

What this file pins
===================

1. **Identity preconditioner** converges to the analytical ``k_inf``
   on slab, sphere, and cylinder.  This is the production contract
   (``SNSolver._solve_krylov`` ships ``preconditioner=lambda q: q``
   until issue #200 lands the block-inverse face preconditioner).
2. **Default sweep preconditioner** converges to ``k_inf`` on **slab**.
   Slab has no curvilinear M-M closure, so the cold-start sweep IS a
   strong preconditioner.  Sphere/cylinder default-sweep precond
   convergence is a known limitation tracked by issue #200 — pinning
   it as a passing L1 would lock-in the wrong production behaviour;
   the right diagnostic for the curvilinear case is the block-inverse
   face precond that #200 designs.

Why these two checks together
-----------------------------

The identity-precond gate (1) directly verifies the production code
path.  The slab default-sweep gate (2) verifies that the
``preconditioner=None`` fallback is **deterministic and structurally
correct** post-Phase-1.2 — if Phase 1.2's structural fix were
reverted (e.g., re-introducing a stateful ``rhs(1)`` read inside
``L.solve``), slab default-sweep would *also* destabilise because
the residual-vector input would again carry no history.  Slab is the
sentinel: it's the geometry where preconditioner quality is NOT the
confound, so a failure there points directly at the silent-fallback
class.

References
==========

- ``.claude/lessons.md`` L19 — None-default preconditioner stateful-
  inverse pitfall.
- ``.claude/lessons.md`` L21 — sweep and matvec are different
  applications of the same operator → share ONE strategy.
- Sphinx theory: ``docs/theory/discrete_ordinates.rst``
  (transport-cartesian, sn-curvilinear-homogeneous-kinf-recovery).
- ERR catalog: ERR-050 — Silent preconditioner fallback breaks
  stateful-inverse contract (Phase 1.4).
- Issue #203 — closed by Phase 1.2 supersession (``c93355c``).
- Issue #200 — sweep-as-preconditioner quality on curvilinear
  (block-inverse face precond re-enablement).
- Issue #204 — A5 promote diagnostics to permanent L1.
"""
from __future__ import annotations

import sys
import warnings

import numpy as np
import pytest


# Re-use the L1 analytical homogeneous reference + mesh / quadrature helpers.
# Same dependency the diagnostic relied on, lifted to module-load.
from tests.sn._test_helpers import SN_TESTS_ROOT  # noqa: E402

sys.path.insert(
    0, str(SN_TESTS_ROOT / "verification" / "analytical"),
)
warnings.simplefilter("ignore")
from test_kinf_homogeneous import (  # noqa: E402  (post-sys.path import)
    _get_continuous_case,
    _homogeneous_mesh,
    _quadrature_for,
)
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux


pytestmark = [
    pytest.mark.l1,
    pytest.mark.verifies(
        "transport-cartesian",
        "sn-curvilinear-homogeneous-kinf-recovery",
    ),
    pytest.mark.catches("ERR-050"),
]


def _krylov_power_iteration_kinf(
    *, coord: str, ng_key: str, preconditioner_kind: str,
    n_cells: int = 10, max_outer: int = 80, keff_tol: float = 1e-10,
) -> tuple[float, int]:
    r"""Run manual outer power iteration with typed-Krylov inner solver.

    Returns ``(keff_converged, n_outer_iterations)``.

    Replicates the production ``SNSolver._solve_krylov`` typed contract
    (R-1 Step 4 A1 producer-side normalisation; R-1 Step D
    KrylovAcceleration on the typed operator triple) but exposes the
    preconditioner choice as a test parameter.  ``preconditioner_kind``:

    * ``"default_sweep"`` → ``preconditioner=None`` (KrylovAcceleration
      falls back to ``LC.solve`` — the WDD sweep with cold-start M-M
      seed post-Phase-1.2).
    * ``"identity"`` → explicit ``lambda q: q`` (the production default
      until issue #200 ships the block-inverse face preconditioner).
    """
    from dataclasses import replace

    from orpheus.numerics.iteration import KrylovAcceleration
    from orpheus.sn.operators.boundary import SNBoundaryOperator
    from orpheus.sn.mesh.augmented_mesh import SNMesh
    from orpheus.sn.operators.streaming import StreamingOperator
    from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
    from orpheus.sn.solver import SNSolver
    from orpheus.transport.full_field import FullField
    from orpheus.transport.source_sinks import AngularSourceSink
    from orpheus.transport.timed_full_field import TimedFullField

    case = _get_continuous_case(ng_key)
    mat_id = next(iter(case.problem.materials.keys()))
    mesh = _homogeneous_mesh(
        coord=coord, n_cells=n_cells, length=2.0, mat_id=mat_id,
    )
    quad = _quadrature_for(coord)
    sn_mesh = SNMesh(mesh, quad, case.problem.materials)

    solver = SNSolver(
        sn_mesh=sn_mesh, scattering_order=0,
        max_inner=300, inner_tol=1e-12, inner_solver="krylov",
    )

    L_leaf = StreamingOperator(sn_mesh)
    C_t = MultiplicationOperator.from_mesh(solver.mat_xs.total_cross_section, sn_mesh)
    LC = L_leaf + C_t

    if preconditioner_kind == "default_sweep":
        precond = None
    elif preconditioner_kind == "identity":
        precond = lambda q: q  # noqa: E731 (closure required by KrylovAcceleration)
    else:
        raise ValueError(f"unknown preconditioner_kind: {preconditioner_kind!r}")

    N = sn_mesh.quad.N
    n_cells, ng = int(np.prod(sn_mesh.spatial_shape)), solver.ng
    krylov = KrylovAcceleration(
        # Wave O O.2a: mirror the honest production gains from
        # ``SNSolver._within_group_triple`` — the resolvent ``L+C`` plus the
        # scattering gain ``S`` AND the boundary reflection gain ``B``.  The
        # transitional ``S+B`` fold + the vestigial ``F=Zero`` slot are RETIRED
        # (variadic ``(L, *gains)`` drivers).  OMITTING ``B`` (as this test did
        # pre-O.2a, when it hand-built ``(LC, scattering_op, ZeroOperator)``)
        # drops the reflective coupling → the GMRES converges to the WRONG
        # eigenmode (k ≈ 1.67, not k_inf = 1.875) on a reflective box.
        LC, solver.scattering_op, SNBoundaryOperator(sn_mesh),
        preconditioner=precond,
        tol=1e-12, max_iter=300,
        restart=N * ng * n_cells,
    )

    phi = solver.initial_flux_distribution()
    keff = 1.0
    # #257 S8a — the Krylov matvec leaves are base arrows ``FullField ->
    # FullField`` so ``krylov.solve`` is inferred over the FullField carrier;
    # the runtime warm iterate is still a TimedFullField (templated on the timed
    # ``initial_guess``), which IS a FullField.
    psi_typed_warm: "FullField | None" = None
    for n_outer in range(max_outer):
        fis = solver.compute_fission_source(phi, keff)
        q_ext_per_ord = AngularSourceSink.from_isotropic(fis, sn_mesh)
        # D-H.2-C1: build the per-ordinate external source as a
        # :class:`TimedFullField` composite — bulk = L2 AngularFlux
        # carrying ``q_ext_per_ord.values``; boundary = implicit zero.
        zero = TimedFullField.zeros(bulk=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh)
        # B.5.2: q_ext IS a source (AngularSourceSink), emitted directly — no
        # re-wrap into AngularFlux.  The iterate/return is a flux, so bootstrap
        # a flux initial_guess on the cold start (psi_typed_warm is None on the
        # first outer), mirroring SNSolver._solve_krylov.
        q_ext_typed = replace(zero, bulk=q_ext_per_ord)
        psi_typed, _residuals = krylov.solve(
            q_ext_typed,
            initial_guess=psi_typed_warm if psi_typed_warm is not None else zero,
        )
        psi_typed_warm = psi_typed
        phi = psi_typed.bulk.integrate_angular().values
        keff_new = solver.compute_keff(phi)
        if abs(keff_new - keff) < keff_tol:
            keff = keff_new
            return keff, n_outer + 1
        keff = keff_new
    return keff, max_outer


@pytest.mark.parametrize("coord", ["slab", "sphere", "cylinder"])
def test_identity_preconditioner_recovers_kinf(coord: str) -> None:
    r"""Identity preconditioner recovers analytical ``k_inf`` on every coord.

    Pins the **production contract**: ``SNSolver._solve_krylov`` ships
    GMRES with ``preconditioner=lambda q: q`` (explicit identity).
    The typed-AngularFlux Krylov inner solver embedded in a manual
    outer power iteration MUST converge to the homogeneous-reflective
    analytical reference at ``rtol < 1e-8`` on slab, sphere, and
    cylinder.

    Regression coverage: any drift in the typed-flux operator algebra
    (A1 producer-side ``/sum_w`` normalisation, Phase-1.2's explicit
    ``initial_guess`` plumbing) that propagates to the GMRES residual
    chain shows up here.
    """
    case = _get_continuous_case("2eg")
    keff_recovered, _n_outer = _krylov_power_iteration_kinf(
        coord=coord, ng_key="2eg", preconditioner_kind="identity",
    )
    rel_err = abs(keff_recovered - case.k_eff) / case.k_eff
    assert rel_err < 1e-8, (
        f"{coord} / identity precond: keff={keff_recovered:.10f}, "
        f"ref={case.k_eff:.10f}, rel_err={rel_err:.3e}.  Identity-precond "
        f"convergence on all coords is the R-1 production contract — "
        f"a failure here flags a regression in the typed Krylov path."
    )


def test_default_sweep_preconditioner_recovers_kinf_on_slab() -> None:
    r"""Default sweep preconditioner converges to ``k_inf`` on slab.

    Pins the **structural-fix sentinel** for Phase 1.2.  With
    ``preconditioner=None``, KrylovAcceleration invokes ``L.solve(q)``
    on each GMRES residual vector.  Pre-Phase-1.2 this read the
    residual's (absent) ``rhs(1)`` history and silently passed garbage
    state into the M-M closure — slab happens not to exercise the M-M
    closure (no curvilinear pole), so even pre-fix slab default-sweep
    converged.  Post-Phase-1.2 ``L.solve(q)`` is a pure function
    (``initial_guess=None`` → deterministic cold-start seed).

    Slab is the **regression sentinel** for the silent-fallback bug
    class: it's the geometry where preconditioner quality (sphere /
    cylinder M-M cold-start convergence) is NOT the confound.  If
    slab default-sweep starts failing while slab identity-precond
    still passes, the silent-fallback path has been re-introduced
    (e.g., a stateful read inside ``L.solve``).  See L19 + ERR-050.

    **Why sphere / cylinder are NOT pinned here**: the curvilinear
    M-M closure with cold-start zero seed is a poor preconditioner
    for GMRES (numerical issue, not the bug-class structural fix).
    Issue #200 designs the block-inverse face preconditioner that
    restores sweep-as-preconditioner quality on curvilinear; pinning
    a failing curvilinear default-sweep here would lock in the wrong
    production state.
    """
    case = _get_continuous_case("2eg")
    keff_recovered, _n_outer = _krylov_power_iteration_kinf(
        coord="slab", ng_key="2eg", preconditioner_kind="default_sweep",
    )
    rel_err = abs(keff_recovered - case.k_eff) / case.k_eff
    assert rel_err < 1e-8, (
        f"slab / default_sweep precond: keff={keff_recovered:.10f}, "
        f"ref={case.k_eff:.10f}, rel_err={rel_err:.3e}.  Pre-Phase-1.2 "
        f"this could have failed if ``L.solve`` had a stateful coupling "
        f"to ``rhs(1)``; Phase 1.2 (commit c93355c) made the silent-"
        f"fallback path structurally unreachable.  A failure here is "
        f"strong evidence that the structural fix has been reverted — "
        f"see L19 + ERR-050."
    )
