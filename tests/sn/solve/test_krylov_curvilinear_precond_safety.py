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
    from orpheus.numerics.coupled_system import CoupledOperator
    from orpheus.sn.coupled_system import (
        _system_a_member,
        build_within_group_system,
    )
    from orpheus.sn.mesh.augmented_mesh import SNMesh
    from orpheus.sn.solver import (
        SNSolver,
        _coupled_flux_state,
        _coupled_source_state,
        _radial_characteristic_source_from_per_ordinate,
    )
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

    if preconditioner_kind == "default_sweep":
        precond = None
    elif preconditioner_kind == "identity":
        precond = lambda q: q  # noqa: E731 (closure required by KrylovAcceleration)
    else:
        raise ValueError(f"unknown preconditioner_kind: {preconditioner_kind!r}")

    # B.2d: consume THE production record (build_within_group_system) —
    # M + N on the coupled pair for a carrying mesh, the bare
    # ``(L+C, (S, B_a))`` pieces seedless (the pre-record hand-build
    # duplicated exactly this composition; omitting B drops the reflective
    # coupling → the WRONG eigenmode, k ≈ 1.67 not 1.875).
    system = build_within_group_system(
        sn_mesh, solver.mat_xs, scattering_op=solver.scattering_op,
    )
    coupled = isinstance(system.resolvent, CoupledOperator)
    zero = TimedFullField.zeros(
        interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh,
    )
    cold = _coupled_flux_state(zero, sn_mesh) if coupled else zero
    krylov = KrylovAcceleration(
        system.resolvent, *system.gains,
        preconditioner=precond,
        tol=1e-12, max_iter=300,
        restart=int(cold.to_flat().size),
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
        # B.5.2: q_ext IS a source (AngularSourceSink), emitted directly — no
        # re-wrap into AngularFlux.  On a carrying mesh (sphere) the coupled
        # rhs pairs the 2-block source with the q½ fold of the (isotropic)
        # fission source — System B's member, mirroring the production
        # ``SNSolver._solve_krylov`` assembly (B.2d); the reflective B arms
        # are folded by ``KrylovAcceleration`` internally (N is in the gains).
        q_a = replace(zero, interior=q_ext_per_ord)
        q_ext_typed = (
            _coupled_source_state(
                q_a,
                _radial_characteristic_source_from_per_ordinate(
                    q_ext_per_ord.values, sn_mesh,
                ),
                sn_mesh, context="test_krylov_precond_safety",
            )
            if coupled else q_a
        )
        psi_typed, _residuals = krylov.solve(
            q_ext_typed,
            initial_guess=psi_typed_warm if psi_typed_warm is not None else cold,
        )
        psi_typed_warm = psi_typed
        phi = _system_a_member(psi_typed).interior.integrate_angular().values
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


@pytest.mark.parametrize("n_cells", [5, 10, 20])
def test_krylov_restart_covers_augmented_composite(n_cells: int) -> None:
    r"""ERR-053-family regression gate (#282/#280 route (a)): the within-group
    GMRES ``restart`` MUST be sized from the FULL composite ravel
    (bulk ⊕ trace ⊕ ψ½-seed), NOT the bulk alone.

    Route (a) grew the Krylov state to a 3-block ``TimedFullField`` on a
    carrying mesh (R12a).  The pre-fix ``n_dof = N·ng·∏spatial`` formula
    (``solver.py`` eigenvalue + fixed-source Krylov drivers) counts only the
    bulk, so on the sphere it is STRICTLY LESS than the raveled ``to_flat``
    dimension — restarted GMRES then cannot span the augmented Krylov
    subspace and the poorly-conditioned curvilinear-eigenvalue inner STALLS
    (residual plateau, scipy ``info > 0``; measured 868 s vs SI ~1 s, and at
    a realistic outer cap it returns a WRONG keff).  The production solver now
    sizes ``n_dof = initial_guess.to_flat().size`` — this gate pins the
    deficit the fix closes so a revert reddens here (fast) instead of
    stalling the sphere eigenvalue wall.  Distinct from issue #200 (the
    identity preconditioner).  numerics-investigator 2026-07-04.
    """
    from orpheus.derivations.common.xs_library import get_mixture
    from orpheus.geometry import CoordSystem
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.sn.mesh.augmented_mesh import SNMesh
    from orpheus.transport.timed_full_field import TimedFullField
    from tests.sn._test_helpers import curvilinear_two_region_mesh

    mesh = curvilinear_two_region_mesh(
        outers=(0.5, 1.0), mat_ids=(2, 0), n_cells=(n_cells, n_cells),
        coord=CoordSystem.SPHERICAL,
    )
    sn = SNMesh(
        mesh, Quadrature.gauss_legendre(8),
        {2: get_mixture("A", "1g"), 0: get_mixture("B", "1g")},
    )
    bulk_only = sn.quad.N * sn.ng * int(np.prod(sn.spatial_shape))
    # B.2d (F3 migration): the Krylov state on a carrying mesh is the COUPLED
    # pair — the honest two-system ravel (2-block ψ_A + System B's composite;
    # NO dead padding).  The pre-fix bulk-only restart under-sizes it by
    # exactly the trace + BOTH System-B legs; the production sizing (the
    # driver state's ``to_flat``) covers it BY CONSTRUCTION — a revert to the
    # bulk formula re-opens this deficit and re-triggers the stall.
    from orpheus.sn.solver import _coupled_flux_state

    pair = _coupled_flux_state(
        TimedFullField.zeros(
            interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn,
        ),
        sn,
    )
    coupled_dim = int(pair.to_flat().size)
    assert bulk_only < coupled_dim, (
        "the sphere coupled ravel no longer exceeds the bulk DOF count — "
        "if System B was removed this gate is stale; otherwise the premise "
        "of the #282 Krylov-restart fix changed."
    )
    deficit = coupled_dim - bulk_only
    assert deficit == (
        int(sn.angular_trace.layout.total_size)
        + int(np.asarray(pair.systems[1].to_flat()).size)
    ), f"restart deficit {deficit} ≠ trace + System B — the ravel layout drifted"


# ── B.2d d3 — G-d3.3: the honest END-TO-END site proof ─────────────────
#
# The count pin above (d2's F3 migration) decomposes the coupled ravel; the
# gates below close the ERR-053 loop END-TO-END: BOTH production Krylov
# drivers (`solve_sn` eigenvalue + `solve_sn_fixed_source`) must CONSTRUCT
# KrylovAcceleration with ``restart`` == that coupled honest ravel — the
# site-plumbing claim (`initial_guess.to_flat()` reaches ``restart``), which
# the count pin alone cannot see.  A revert to the d1-era bulk-only formula
# reds HERE in seconds instead of re-opening the measured 868 s sphere
# stall.


def _carrying_sphere_case(n_cells: int = 6):
    """Homogeneous FISSILE 2G reflective carrying sphere (the kinf case
    materials — fissile, so BOTH production drivers run on it)."""
    from orpheus.sn.mesh.augmented_mesh import SNMesh

    case = _get_continuous_case("2eg")
    mat_id = next(iter(case.problem.materials.keys()))
    mesh = _homogeneous_mesh(
        coord="sphere", n_cells=n_cells, length=2.0, mat_id=mat_id,
    )
    quad = _quadrature_for("sphere")
    sn = SNMesh(mesh, quad, case.problem.materials)
    return case.problem.materials, mesh, quad, sn


def _expected_coupled_restart(sn) -> tuple[int, int]:
    """(coupled honest restart, bulk-only count) from the SPACES — explicit
    arithmetic independent of the driver's own ``to_flat`` read (the d2
    count pin above proves the decomposition; HERE the sum is the oracle
    for the SITE plumbing)."""
    from orpheus.transport.radial_characteristic_field import (
        RadialCharacteristicField,
    )

    bulk = int(sn.quad.N * sn.ng * int(np.prod(sn.spatial_shape)))
    trace = int(sn.angular_trace.layout.total_size)
    size_b = int(np.asarray(
        RadialCharacteristicField.from_mesh(sn).to_flat(),
    ).size)
    return bulk + trace + size_b, bulk


def _install_restart_spy(monkeypatch) -> list[int]:
    """Wrap ``KrylovAcceleration.__init__`` capturing every constructed
    ``restart`` (the Mode-11 sentinel: the drivers must actually route
    through the construction this gate inspects)."""
    from orpheus.numerics.iteration import KrylovAcceleration

    captured: list[int] = []
    orig_init = KrylovAcceleration.__init__

    def _spy_init(self, *args, **kwargs):
        orig_init(self, *args, **kwargs)
        captured.append(int(self.restart))

    monkeypatch.setattr(KrylovAcceleration, "__init__", _spy_init)
    return captured


@pytest.mark.catches("ERR-053")
@pytest.mark.parametrize("driver", ["eigenvalue", "fixed_source"])
def test_g_d3_3_production_sites_size_restart_from_the_coupled_ravel(
    driver: str, monkeypatch,
) -> None:
    r"""G-d3.3 END-TO-END: the production driver constructs its
    within-group GMRES with ``restart`` == the COUPLED honest ravel
    (bulk ⊕ trace ⊕ BOTH System-B legs) on the carrying sphere.

    Loose outer tolerances — the claim is the CONSTRUCTION plumbing, not
    convergence; every KrylovAcceleration the driver builds is captured
    and compared against the space-derived sum."""
    from orpheus.sn.solver import solve_sn, solve_sn_fixed_source

    materials, mesh, quad, sn = _carrying_sphere_case()
    expected, bulk_only = _expected_coupled_restart(sn)
    if not expected > bulk_only:
        pytest.fail("premise broke: the coupled ravel no longer exceeds "
                    "the bulk count — if System B was removed this gate "
                    "is stale")

    captured = _install_restart_spy(monkeypatch)

    if driver == "eigenvalue":
        solve_sn(
            materials, mesh, quad, inner_solver="krylov",
            max_outer=2, keff_tol=1e-3, flux_tol=1e-3,
            max_inner=60, inner_tol=1e-8,
        )
    else:
        solve_sn_fixed_source(
            materials=materials, mesh=mesh, quadrature=quad,
            external_source=np.ones((quad.N, sn.ng, sn.nx)),
            inner_solver="krylov", max_inner=60, inner_tol=1e-8,
        )

    if not captured:
        pytest.fail(f"[{driver}] Mode-11: the KrylovAcceleration spy never "
                    f"fired — the driver no longer constructs the "
                    f"within-group GMRES this gate inspects")
    if set(captured) != {expected}:
        pytest.fail(
            f"[{driver}] production restart(s) {sorted(set(captured))} ≠ "
            f"the coupled honest ravel {expected} (bulk-only = {bulk_only})"
            f" — ERR-053 re-opened at a driver site"
        )


@pytest.mark.catches("ERR-053")
def test_g_d3_3_site_gate_has_teeth(monkeypatch) -> None:
    r"""The stall-deficit tooth: force the d1-era bulk-only ``n_dof``
    through the ONE production helper (``_within_group_krylov``) — the
    captured restart then lands EXACTLY at the bulk-only count and
    mismatches the coupled ravel, proving the site gate above reds on a
    formula revert (the capture is live, not vacuous)."""
    import orpheus.sn.solver as solver_mod

    materials, mesh, quad, sn = _carrying_sphere_case()
    expected, bulk_only = _expected_coupled_restart(sn)

    orig_wgk = solver_mod._within_group_krylov

    def _bulk_only_wgk(LC, *gains, n_dof, max_iter, tol):
        return orig_wgk(LC, *gains, n_dof=bulk_only, max_iter=max_iter,
                        tol=tol)

    monkeypatch.setattr(solver_mod, "_within_group_krylov", _bulk_only_wgk)
    captured = _install_restart_spy(monkeypatch)

    solver_mod.solve_sn_fixed_source(
        materials=materials, mesh=mesh, quadrature=quad,
        external_source=np.ones((quad.N, sn.ng, sn.nx)),
        inner_solver="krylov", max_inner=40, inner_tol=1e-6,
    )

    if not captured:
        pytest.fail("Mode-11: the spy never fired under the mutation")
    if set(captured) == {expected}:
        pytest.fail("the bulk-only mutation did NOT move the captured "
                    "restart — the site gate has no teeth")
    if set(captured) != {bulk_only}:
        pytest.fail(f"mutated captures {sorted(set(captured))} ≠ bulk-only "
                    f"{bulk_only} — the tooth is mis-wired")
