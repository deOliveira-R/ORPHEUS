"""MOVED from tests/gates/numerics/test_iteration.py (layer DAG Phase 2 refresh, specs_2r §5). Bodies unchanged."""

from __future__ import annotations
import warnings
import numpy as np
import pytest
from orpheus.data.macro_xs.synthetic import get_mixture
from orpheus.geometry import StructuredGeometry
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.green_operator import GreenOperator
from orpheus.numerics.iteration import KEigenvalue
from orpheus.numerics.pencil import OperatorPencil
from orpheus.numerics.posing import K_MAP, EigenPosing
from orpheus.numerics.operator import InverseOperator, LinearOperator
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.geometry.deck import MirrorDeck




@pytest.fixture
def rng() -> np.random.Generator:
    return np.random.default_rng(20260509)




def _sn_composite_triple():
    r"""The HONEST typed-composite SN operator triple + solve_sn reference.

    #276 A4 activation fixture: ``build_within_group_system`` supplies the
    production splitting — ``LC`` (the invertible resolvent whose ``solve``
    is the sweep), the gains ``(S, N2N, B_a)`` summed into KEigenvalue's single
    coupling slot, and the production ``FissionOperator`` — all acting on
    typed :class:`FullField` composites (no scalar adapter shims; contrast
    the legacy-shim gate above, kept as the pre-A4 record).
    """
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.sn.coupled_system import build_within_group_system
    from orpheus.sn.problem import SNProblem
    from orpheus.sn.solver import SNSolver, solve_sn
    from orpheus.transport.fields.angular_flux import AngularFlux
    from orpheus.transport.full_field import FullField

    mix = get_mixture("A", "2g")
    materials = {0: mix}
    mesh = Mesher(StructuredGeometry.from_homogeneous(5.0, MirrorDeck())).partition(
        CellsByCount.uniform_width(10),
    ).mesh
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    ref = solve_sn(
        materials, mesh, quad,
        inner_solver="source_iteration", scattering_order=0,
        max_outer=500, keff_tol=1e-9, flux_tol=1e-8,
        max_inner=500, inner_tol=1e-10,
    )
    if ref.outcome.keff is None:  # explicit narrow — fires under -O (Mode 8)
        pytest.fail("solve_sn returned no eigenvalue — the reference leg "
                    "of the A4 activation fixture is broken.")
    sn = SNProblem(mesh, quad, materials, scattering_order=0)
    solver = SNSolver(sn)
    system = build_within_group_system(
        sn, solver.problem.mat_xs,
    )
    from functools import reduce
    from operator import add

    # ALL lagged gains summed into KEigenvalue's single coupling slot —
    # (S, N2N, B_a) since §14.1; a positional pair-sum here silently
    # dropped B_a when the tuple grew (caught by the A4 smoke's own
    # eig(A†)=eig(A) gate).
    from orpheus.sn.splitting import Splitting, resolve_schedule

    splitting = Splitting.from_schedule(system, resolve_schedule(sn, "jacobi"))
    S_total = reduce(add, splitting.explicit)
    # CS4c step 4: the composite triple's F is the ANGULAR binding — the
    # frame-conjugated FissionOperator on the SAME composite space the
    # loss members ride (the pencil peer; its .H composes the Riesz legs
    # around the reversed full_fission_kernel product). The hub's
    # ``fission.isotropic_energy`` is the scalar ENERGY binding the
    # bare-array k-outer feeds — the same datum's other face.
    from orpheus.transport.operators.fission import FissionOperator

    F_composite = FissionOperator.from_solver_data(
        mat_xs=solver.problem.mat_xs, space=sn.full_field_space,
    )
    guess = FullField(
        interior=AngularFlux(values=np.ones((sn.quad.n_ordinates, sn.ng, *sn.spatial_shape)), space=sn.angular_bulk_space),
        boundary=AngularBoundaryFlux.zeros(sn.angular_trace),
    )
    return float(ref.outcome.keff), splitting.implicit, S_total, F_composite, guess, mix






















@pytest.mark.l1
@pytest.mark.verifies("multigroup")
def test_keigenvalue_matches_solve_sn_2g_slab():
    """L1 gate: KEigenvalue on the SN operator triple matches solve_sn.

    Build a 2-group homogeneous-material 1-D slab, run both:

    * :func:`solve_sn` (the legacy power_iteration-based path), and
    * :class:`KEigenvalue` directly on
      ``(StreamingCollisionOperator (= L + C), ScatteringOperator, FissionOperator)``
      with adapter shims that present scalar-flux shapes consistent
      across operators.

    Assert recovered keff agrees to 1e-9.  Both paths use the SAME
    underlying operators — the only difference is whether the iteration
    primitive is the legacy ``power_iteration(solver)`` (Wave-E pre-
    Round-1) or the new ``KEigenvalue`` (Wave-E Round 1).

    R8 rewire (#259 P1, 2026-07-03 — the estimator-injection seam
    retired): the pre-R8 version injected ``solver.compute_keff`` as
    KEigenvalue's ``keff_estimator``.  Post-R8 the k assertion routes
    through the THEOREM instead: the eigen-FLUX fixed point is
    estimator-independent (unit-production renormalisation cancels the
    k scaling), and at a converged eigenpair every CONSISTENT estimator
    agrees — so evaluating the SN method-layer functional
    ``solver.compute_keff`` at KEigenvalue's converged flux must
    reproduce ``solve_sn``'s reported k.  KEigenvalue's OWN reported k
    is deliberately NOT asserted on this stack: the scalar-level
    ``A_inv_adapter`` cannot advertise an honest ``.apply`` (the true
    within-group loss action is angular), so the hardwired Rayleigh
    denominator ``Σ(Aψ)−Σ(Sψ)`` is off-contract here — a limitation the
    injection seam used to paper over, exposed (by design) at R8.  The
    honest-triple version of the reported-k leg lives with the
    operator-level pins in ``test_estimators_as_functionals.py``.
    """
    # Suppress the eigenvalue.py deprecation warning for this test —
    # we are deliberately exercising the legacy path as a reference.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)

        from orpheus.transport.operators.fission import FissionOperator
        from orpheus.sn.problem import SNProblem
        from orpheus.numerics.quadrature import Quadrature
        from orpheus.transport.operators.scattering import ScatteringOperator
        from orpheus.sn.solver import SNSolver, solve_sn
        from tests.gates.sn._test_helpers import reflect_outflow_into_inflow, sweep_once

    # 2-group homogeneous 1-D slab — the same canonical fixture
    # ``test_solver_components.py`` uses for component checks.
    mix = get_mixture("A", "2g")
    materials = {0: mix}
    mesh = Mesher(StructuredGeometry.from_homogeneous(5.0, MirrorDeck())).partition(
        CellsByCount.uniform_width(10),
    ).mesh
    quad = Quadrature.gauss_legendre(n_ordinates=8)

    # Reference: solve_sn (legacy power_iteration path).
    ref = solve_sn(
        materials, mesh, quad,
        inner_solver="source_iteration",
        scattering_order=0,
        max_outer=500,
        keff_tol=1e-9,
        flux_tol=1e-8,
        max_inner=500,
        inner_tol=1e-10,
    )
    expected_keff = ref.outcome.keff
    if expected_keff is None:
        pytest.fail("solve_sn reference returned no eigenvalue.")

    # Build the SN operator triple from the same precomputed solver
    # data used for the reference run.  ``A_inv_adapter`` (defined
    # below) wraps ``sweep_once`` directly; the
    # :class:`StreamingCollisionOperator` (= ``L + C``) on the SNSolver is
    # unused here.  Solver instance retained to provide
    # ``solver.scattering_op`` / ``solver.problem.fission.isotropic_energy`` / ``solver.mat_xs``.
    problem = SNProblem(mesh, quad, materials, scattering_order=0)
    solver = SNSolver(problem)
    # The canonical S, F operators built directly from solver state.
    S = solver.problem.system.factors.scattering
    F = solver.problem.fission.isotropic_energy

    # ── Adapter shims to keep the iteration primitive scalar-flux-only.
    #
    # The SN operator triple has an internal shape-mismatch for
    # historical reasons (Wave D Round 3 docstring §"Vector layout"):
    #   • F.apply takes scalar phi (nx,ny,ng), returns scalar Q (nx,ny,ng)
    #   • S.apply takes angular psi (N,nx,ny,ng), returns angular Q (N,nx,ny,ng)
    #   • L.solve takes Q+psi_bc+Q_aniso, returns (psi, phi) tuple.
    # Round 2 normalises this; for the L1 gate test we wrap each
    # operator into a thin scalar-in/scalar-out facade.
    # Issue #197 PR-TYPED-2 — typed AngularBoundaryFlux replaces psi_bc: dict.
    boundary_flux = AngularBoundaryFlux.zeros(problem.angular_trace)

    class A_inv_adapter(LinearOperator):
        """Adapter: rhs (ng, nx, ny) → phi via the unified sweep.

        Issue #196 PR-INDEX-5: principled layout throughout.  Returns
        scalar flux (drops angular flux for shape consistency with
        F.apply / S.apply scalar facade).

        #226 step 3: carries the invertibility pair (``is_invertible`` +
        ``inverse()``) so :class:`KEigenvalue` — which now BUILDS the
        inverse and hands it to the inner driver — can lift this leaf's
        ``solve`` through the generic :class:`InverseOperator`.
        """
        # S4-amendment: the base DEMANDS an answer from every subclass; this
        # double is a deliberately-unbound probe, so it DECLARES the unbound
        # state instead of inheriting a silent default (which no longer exists).
        domain = None
        codomain = None

        @property
        def is_invertible(self) -> bool:
            return True

        def inverse(self) -> InverseOperator:
            return InverseOperator(self)

        def apply(self, phi):  # not used by the iteration primitive
            return phi

        def solve(self, rhs):
            # R-1 Step 4 A1: single per-ordinate source carrier.
            # ``rhs`` is bare ndarray (ng, nx, ny) — wrap via the
            # canonical iso → per-ord factory at the adapter boundary.
            from orpheus.transport.fields import AngularSourceSink
            source = AngularSourceSink.from_isotropic(rhs, problem)
            # Wave O (#208) O.4a.2 — the bare ``transport_sweep`` no longer
            # re-applies the reflective BC at entry; drive the −B coupling
            # explicitly (reflect the persisted outflow — ``boundary_flux``
            # is the closure-scoped partner-flux carrier — into the inflow
            # slots) before each sweep — the sweep-tier gates' inter-sweep −B
            # (the drivers deliver it as the ``B`` gain; #448).
            reflect_outflow_into_inflow(boundary_flux, problem)
            _angular, scalar = sweep_once(
                source, solver.problem.mat_xs.total_cross_section, problem,
                boundary_flux,
            )
            return scalar

    class S_scalar_adapter(LinearOperator):
        """Adapter: phi (ng, nx, ny) → P0 scattering source (ng, nx, ny).

        Issue #196 PR-INDEX-5: principled end-to-end.
        """
        # S4-amendment: the base DEMANDS an answer from every subclass; this
        # double is a deliberately-unbound probe, so it DECLARES the unbound
        # state instead of inheriting a silent default (which no longer exists).
        domain = None
        codomain = None

        def apply(self, phi):
            Q = np.zeros_like(phi)
            S.transfer.add_p0_source(Q, phi)
            # §14.1: the (n,2n) verb lives on the solver-held N2N binding.
            solver.problem.system.factors.n2n.isotropic_energy.transfer.add_p0_source(Q, phi)
            return Q

    class F_scalar_adapter(LinearOperator):
        """Adapter: phi (ng, nx, ny) → fission source (ng, nx, ny).

        Issue #196 PR-INDEX-5: principled end-to-end.
        """
        # S4-amendment: the base DEMANDS an answer from every subclass; this
        # double is a deliberately-unbound probe, so it DECLARES the unbound
        # state instead of inheriting a silent default (which no longer exists).
        domain = None
        codomain = None

        def apply(self, phi):
            return F.apply(phi)

    A_adapt = A_inv_adapter()
    S_adapt = S_scalar_adapter()
    F_adapt = F_scalar_adapter()

    ng = solver.ng
    # Issue #196 PR-INDEX-5: principled initial guess.
    initial = np.ones((ng, *problem.spatial_shape))

    ke = KEigenvalue(
        EigenPosing(OperatorPencil(A_adapt - S_adapt, F_adapt), K_MAP), GreenOperator(A_adapt, S_adapt),
        max_outer=500, keff_tol=1e-9, flux_tol=1e-8,
        max_inner=500, inner_tol=1e-10,
    )
    _o = ke.solve(initial_guess=initial)
    _keff_rayleigh, keff_history, phi_ke = _o.keff, _o.keff_history, _o.flux_distribution

    # The theorem-form assertion (R8): the flux fixed point is shared, so
    # the SN method-layer functional evaluated at KEigenvalue's converged
    # flux reproduces solve_sn's reported k.  (All-reflective slab —
    # solver.compute_keff needs no boundary trace; the functional is a
    # scale-invariant ratio, so KEigenvalue's unit-production
    # normalisation drops out.)
    keff_at_ke_flux = solver.compute_keff(phi_ke)
    assert abs(keff_at_ke_flux - expected_keff) < 1e-9, (
        f"SN functional at KEigenvalue's converged flux "
        f"keff={keff_at_ke_flux!r} differs from solve_sn "
        f"keff={expected_keff!r} by {abs(keff_at_ke_flux-expected_keff):.2e}; "
        f"Rayleigh history[-3:]={keff_history[-3:]}"
    )


@pytest.mark.l1
@pytest.mark.verifies("multigroup")
def test_keigenvalue_honest_composite_triple_matches_solve_sn():
    r"""#276 A4 activation CONTROL: the honest typed-composite triple.

    ``KEigenvalue(posing, GreenOperator(LC, S+B_a))`` on :class:`FullField` composites — the
    production operators, NO adapter shims — converges to ``solve_sn``'s
    eigenvalue, asserted on KEigenvalue's OWN hardwired Rayleigh estimator
    (on-contract for the first time: ``A.apply`` is the TRUE composite
    ``L+C`` action, which the legacy scalar shim above could not provide).
    This is the FORWARD control leg the daggered smoke below rests on: any
    carrier-honesty regression in the posing layer (the ``_ravel`` sums,
    the ``_l2_norm`` convergence test, the frozen-alias guess stash) reds
    here before any adjoint claim is evaluated.
    """
    ref_keff, LC, S_total, F, guess, _mix = _sn_composite_triple()
    ke = KEigenvalue(
        EigenPosing(OperatorPencil(LC - S_total, F), K_MAP), GreenOperator(LC, S_total),
        max_outer=500, keff_tol=1e-9, flux_tol=1e-8,
        max_inner=500, inner_tol=1e-10,
    )
    _o = ke.solve(initial_guess=guess)
    keff, history, psi = _o.keff, _o.keff_history, _o.flux_distribution
    np.testing.assert_allclose(
        keff, ref_keff, rtol=0, atol=1e-9,
        err_msg=f"honest-composite KEigenvalue k={keff!r} differs from "
        f"solve_sn k={ref_keff!r}; history[-3:]={history[-3:]}",
    )
    from orpheus.transport.full_field import FullField as _FF

    if not isinstance(psi, _FF):
        pytest.fail(
            f"the honest-composite drive must return a typed composite "
            f"iterate; got {type(psi).__name__}."
        )


@pytest.mark.l1
def test_keigenvalue_daggered_triple_adjoint_smoke():
    r"""#276 A4 activation: ``KEigenvalue(posing†, GreenOperator(LC.H, (S+B_a).H))`` — the
    daggered posing through the UNCHANGED ``power_iteration``.

    The adjoint eigenproblem ``A_loss† ψ* = F† ψ*/k`` posed purely by
    dagger-ing the forward triple's members (the operator algebra IS the
    implementation — zero adjoint-specific loop code).  Gates, per Mode-12
    (k is NEVER sole evidence — ``eig(Aᵀ) = eig(A)`` makes every k-level
    functional blind to the whole adjoint mutation class):

    * ``k_adj == k_fwd == solve_sn`` (the exact algebraic equality);
    * the ∞-medium adjoint SPECTRUM equals the corrected closed form
      (``kinf_and_adjoint_spectrum_homogeneous`` — the dominant eigenvector
      of ``(Aᵀ)⁻¹Fᵀ``; on the flat reflective mesh the discrete daggered
      solve reproduces the 0-D energy shape).  This vector-level leg is
      what caught the reference's original ``eig(Mᵀ)`` factor-order
      degeneracy (≡ ν̂Σf, zero A-physics) on this machinery's first run —
      and it is asserted NOT to be that degenerate vector here.

    The full P1.3/P1.4 batteries (4G, heterogeneous slab, sphere, F†/S†
    mutations) land with the solver entries; this row is the posing-layer
    activation smoke.
    """
    from orpheus.derivations.common.eigenvalue import (
        kinf_and_adjoint_spectrum_homogeneous,
    )

    ref_keff, LC, S_total, F, guess, mix = _sn_composite_triple()
    ke_adj = KEigenvalue(
        EigenPosing(OperatorPencil(LC.H - S_total.H, F.H), K_MAP), GreenOperator(LC.H, S_total.H),
        max_outer=500, keff_tol=1e-9, flux_tol=1e-8,
        max_inner=500, inner_tol=1e-10,
    )
    _o = ke_adj.solve(initial_guess=guess)
    k_adj, history, psi_star = _o.keff, _o.keff_history, _o.flux_distribution
    np.testing.assert_allclose(
        k_adj, ref_keff, rtol=0, atol=1e-9,
        err_msg=f"daggered-triple k_adj={k_adj!r} differs from the forward "
        f"solve_sn k={ref_keff!r} — eig(A†)=eig(A) is violated by the "
        f"posing; history[-3:]={history[-3:]}",
    )

    k_cf, phi_star_cf = kinf_and_adjoint_spectrum_homogeneous(
        np.asarray(mix.SigT),
        np.asarray(mix.SigS[0].todense()),
        np.asarray(mix.SigP),
        np.asarray(mix.chi),
    )
    np.testing.assert_allclose(
        k_adj, k_cf, rtol=0, atol=1e-9,
        err_msg="daggered k does not match the closed-form k∞ anchor.",
    )
    bulk = np.asarray(psi_star.interior.values)  # (N, ng, *spatial)
    spec = bulk.mean(axis=(0, *range(2, bulk.ndim)))
    spec = spec / np.linalg.norm(spec)
    np.testing.assert_allclose(
        spec, phi_star_cf, rtol=1e-8,
        err_msg="the SN daggered solve's ∞-medium energy spectrum does not "
        "match the closed-form adjoint eigenvector of (Aᵀ)⁻¹Fᵀ — the "
        "flux-shape leg (the k-blind Mode-12 catcher).",
    )
    nsf_hat = np.asarray(mix.SigP) / np.linalg.norm(np.asarray(mix.SigP))
    if np.allclose(spec, nsf_hat, rtol=1e-6):
        pytest.fail(
            "the daggered spectrum equals ν̂Σf — the eig(Mᵀ) factor-order "
            "degenerate; either the solve or the reference regressed to "
            "the wrong resolvent ordering."
        )


