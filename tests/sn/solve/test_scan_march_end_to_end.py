r"""S5.2 — ScanMarch end-to-end solver gates: G4 FP-invariance + G6 eigenvalue (#222).

The S5.1 ``ScanMarch`` representation was verified at the SWEEP/MATVEC layer
only (G2 oracle nulp + G3 conditioning); these are the deferred-to-after-S6.5
END-TO-END gates from the test-architect plan
(``scan_march_verification.md`` §4 G4 + the gate-table G6 row): drive the FULL
production solvers (``solve_sn`` / ``solve_sn_fixed_source``) with ScanMarch
forced as the operator's representation and pin the converged fixed point.

The row-march schedule is an iteration-INVARIANT change: it MUST NOT move the
converged :math:`\psi^*`/:math:`k` (only the per-sweep FP-association, which
the outer iteration washes out).  vv-principles **Mode 9** demands the
FP-invariance be verified on configs that BREAK the degenerate coincidence —
NEVER the isotropic-reflective box where a wrong formulation is accidentally
exact:

* **G4.a** — anisotropic (P1) + heterogeneous (fuel|moderator) + VACUUM
  (streaming) on the x-pair: the ``si_2d_p1_aniso_het`` config of the affine
  golden, with the non-flat degenerate-gate guard.
* **G4.b** — REFLECTIVE on BOTH axis pairs + a level-symmetric (O_h, shared
  faces) cubature, so every octant's outflow on every face is another
  octant's inflow and the OUTFLOW SHED ORDER is load-bearing (the ERR-056
  failure class — the scan-march sheds x-outflow as the last scan value per
  line and y-outflow as the last row's ψy, vs the wavefront's per-octant
  frontier shed).  LIMITATION (stated per the gate plan): in pure 2-D
  Cartesian the in-plane octants are the 4 sign-quadrants and the shared-face
  coupling is the per-axis-edge reflective inflow=outflow; a genuinely
  diagonal cubature with off-axis ordinates sharing one face is a 3-D /
  curvilinear stressor.  This test pins the reflective shed order to the
  extent d=2 Cartesian admits it; the full ERR-056 stressor is exercised at
  d=3 (deferred — no 3-D quadrature yet).
* **G6** — the end-to-end eigenvalue gates: closed-form
  :math:`k_\infty = \lambda_{\max}(A^{-1}F)` (homogeneous, ≥2G — no 1G
  eigenvalue evidence per the cardinal rule) and SI ≡ Krylov flux-SHAPE
  agreement on the heterogeneous non-flat 2G config, both with ScanMarch
  forced.  SI≡Krylov alone is twin agreement (necessary, not sufficient);
  the k_inf leg supplies the structural-independence anchor.

FORCING (post-S6.5; INVERTED at the S6.9 Fork-B2 flip):
``loss_representation.default_for`` is patched to substitute
``MovingFrontierWindow`` wherever the selector picks ``ScanMarch`` on a
multi-D mesh.  HISTORY: at S5.2 birth the window was the default and
ScanMarch was forced; the 2026-06-11 flip made ScanMarch the production
default, so the polarity inverted — the DEFAULT leg now runs ScanMarch and
the FORCED leg runs the window, which doubles as the window's end-to-end
coverage now that it is a selectable (non-default) peer.  The FP-invariance
claim is symmetric, so the gate's meaning is unchanged: two schedules, one
converged fixed point.  Every door reads the module attribute at call time —
the operator's ``loss_representation`` cached_property (lazy import) and
(via the operator's instance, S6.5) the G-S resolvent; the third,
operator-free ``transport_sweep`` door retired at the coupled-block
campaign step 6 — so ONE patch forces the whole solve.  The forcing is a
CONTEXT MANAGER, not a test-scoped fixture, because the REFERENCE leg of each
FP-invariance pair must run UNFORCED (a fixture would force both legs and
compare the window to itself).  NON-VACUITY is explicit: the context counts
``MovingFrontierWindow._sweep_interior`` / ``loss_action`` executions and
every test asserts the forced path actually RAN (a silent fall-back to the
default would otherwise make these gates scanmarch≡scanmarch tautologies).

The fixed points are compared at SOLVER tolerance (NOT nulp): the schedules
differ in transient FP-association BY CONSTRUCTION; only the converged values
are schedule-invariant.  This is the explicit "do NOT demand bit-identity
across schedules" boundary (vv-principles §bit-identity-vs-principled).
"""
from __future__ import annotations

from contextlib import contextmanager
from unittest.mock import patch

import numpy as np
import pytest

from orpheus.derivations import get
from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import solve_sn, solve_sn_fixed_source

pytestmark = pytest.mark.l1


@contextmanager
def window_forced():
    """Force MovingFrontierWindow wherever the selector picks the ScanMarch
    default on a multi-D mesh — and COUNT its kernel executions so vacuous
    fall-back cannot pass.

    Yields the hit-counter dict; callers assert the relevant counter is
    non-zero after the forced solve.
    """
    from orpheus.sn import loss_representation as lr

    real_default = lr.default_for
    real_interior = lr.MovingFrontierWindow._sweep_interior
    real_action = lr.MovingFrontierWindow.loss_action
    hits = {"sweep_interior": 0, "loss_action": 0}

    def forced(mesh, spatial_closure, angular_closure):
        rep = real_default(mesh, spatial_closure, angular_closure)
        if isinstance(rep, lr.ScanMarch) and not mesh.is_1d:
            return lr.MovingFrontierWindow(mesh, spatial_closure, angular_closure)
        return rep

    def spy_interior(self, *args, **kwargs):
        hits["sweep_interior"] += 1
        return real_interior(self, *args, **kwargs)

    def spy_action(self, *args, **kwargs):
        hits["loss_action"] += 1
        return real_action(self, *args, **kwargs)

    with patch.object(lr, "default_for", forced), \
         patch.object(lr.MovingFrontierWindow, "_sweep_interior", spy_interior), \
         patch.object(lr.MovingFrontierWindow, "loss_action", spy_action):
        yield hits


def _require_ran(hits: dict, counter: str) -> None:
    """Non-vacuity tripwire: the forced leg must have RUN the window
    kernel, else the FP comparison silently degraded to default≡default."""
    if hits[counter] == 0:
        pytest.fail(
            f"the forced leg never executed MovingFrontierWindow.{counter} — "
            "the forcing fell back to the default; this gate is vacuous."
        )


def _build_2d_aniso_het_vacuum():
    """The ``si_2d_p1_aniso_het`` config (mirrors the affine golden): 8×4,
    fuel|moderator across x, P1, 2G, VACUUM-x / reflective-y, LS-4."""
    fuel = get_mixture("A", "2g")
    mod = get_mixture("B", "2g")
    nx, ny = 8, 4
    mat = np.zeros((nx, ny), dtype=int)
    mat[:4, :] = 2
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 2.0, nx + 1),
        edges_y=np.linspace(0.0, 1.0, ny + 1),
        mat_map=mat,
        bc_xmin=BC("vacuum"), bc_xmax=BC("vacuum"),
        bc_ymin=BC("reflective"), bc_ymax=BC("reflective"),
    )
    quad = Quadrature.level_symmetric(sn_order=4)
    sum_w = float(quad.weights.sum())
    q_ext = np.full((quad.N, 2, nx, ny), 1.0 / sum_w)
    return {2: fuel, 0: mod}, mesh, quad, q_ext


def _build_2d_het_all_reflective():
    """G4.b config: same fuel|moderator 8×4 heterogeneity, REFLECTIVE on
    both axis pairs (every face's outflow re-enters as another octant's
    inflow → the shed order is load-bearing), level-symmetric O_h cubature.

    The volumetric source lives in the FUEL HALF ONLY: under all-reflective
    BCs a uniform source flattens the flux (measured max/min = 1.067, which
    trips the non-flat guard — the degenerate regime Mode 9 forbids); the
    half-domain source drives a genuine x-gradient so the reflected
    boundary coupling carries real information.
    """
    fuel = get_mixture("A", "2g")
    mod = get_mixture("B", "2g")
    nx, ny = 8, 4
    mat = np.zeros((nx, ny), dtype=int)
    mat[:4, :] = 2
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 2.0, nx + 1),
        edges_y=np.linspace(0.0, 1.0, ny + 1),
        mat_map=mat,
        bc_xmin=BC("reflective"), bc_xmax=BC("reflective"),
        bc_ymin=BC("reflective"), bc_ymax=BC("reflective"),
    )
    quad = Quadrature.level_symmetric(sn_order=4)
    sum_w = float(quad.weights.sum())
    q_ext = np.zeros((quad.N, 2, nx, ny))
    q_ext[:, :, :4, :] = 1.0 / sum_w  # the fuel half only
    return {2: fuel, 0: mod}, mesh, quad, q_ext


def _solve_fixed_source(builder):
    mats, mesh, quad, q = builder()
    return solve_sn_fixed_source(
        mats, mesh, quad, q, scattering_order=1,
        inner_solver="source_iteration", max_inner=3000, inner_tol=1e-12,
    )


def _assert_nonflat(phi: np.ndarray) -> None:
    """Degenerate-gate guard: the flux MUST be genuinely non-flat, else the
    redistribution / shed-order terms are out of play and FP-equality is
    vacuous (vv Mode 9 anti-recommendation)."""
    prof = phi[0].mean(axis=1)
    ratio = float(prof.max() / prof.min())
    if ratio <= 1.2:
        pytest.fail(
            f"flux too flat (max/min={ratio:.3f}); redistribution not "
            "exercised — the FP-invariance gate is degenerate."
        )


# ═══════════════════════════════════════════════════════════════════
# G4 — Mode-9 FP-invariance of the row-march schedule (fixed-source)
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.verifies("transport-cartesian-2d", "multigroup")
def test_g4a_fixed_source_fp_invariant_aniso_het_vacuum():
    """[G4.a] ScanMarch converged fixed-source ψ* ≡ window fixed point —
    anisotropic + heterogeneous + vacuum (the degeneracy-breaking config)."""
    sol_default = _solve_fixed_source(_build_2d_aniso_het_vacuum)  # ScanMarch
    phi_sm = np.asarray(sol_default.scalar_flux.values, dtype=np.float64)
    _assert_nonflat(phi_sm)

    with window_forced() as hits:
        sol_window = _solve_fixed_source(_build_2d_aniso_het_vacuum)
    _require_ran(hits, "sweep_interior")
    phi_w = np.asarray(sol_window.scalar_flux.values, dtype=np.float64)

    # Converged values agree to solver tolerance — NOT nulp (the schedules
    # differ at transient FP-association by construction).
    np.testing.assert_allclose(
        phi_sm, phi_w, rtol=1e-6, atol=1e-10,
        err_msg=(
            "the two schedules (ScanMarch default vs forced window) "
            "converged to DIFFERENT fixed-source fluxes — the schedule pair "
            "is NOT iteration-invariant (Mode-9 FP violation)."
        ),
    )


@pytest.mark.verifies("transport-cartesian-2d", "multigroup")
@pytest.mark.catches("ERR-056")
def test_g4b_fixed_source_fp_invariant_reflective_shed_order():
    """[G4.b] ScanMarch converged ψ* ≡ window fixed point — ALL-reflective,
    level-symmetric: the d=2 reflective shed-order pin (ERR-056 class).

    Every face's outflow is reflected back as inflow, so the scan-march's
    different shed mechanics (last-scan-value-per-line x-outflow / last-row
    ψy y-outflow) are load-bearing for the converged boundary coupling.  See
    the module docstring for the honest d=2 LIMITATION statement — the full
    diagonal shared-face ERR-056 stressor needs d=3.
    """
    sol_default = _solve_fixed_source(_build_2d_het_all_reflective)  # ScanMarch
    phi_sm = np.asarray(sol_default.scalar_flux.values, dtype=np.float64)
    _assert_nonflat(phi_sm)

    with window_forced() as hits:
        sol_window = _solve_fixed_source(_build_2d_het_all_reflective)
    _require_ran(hits, "sweep_interior")
    phi_w = np.asarray(sol_window.scalar_flux.values, dtype=np.float64)

    np.testing.assert_allclose(
        phi_sm, phi_w, rtol=1e-6, atol=1e-10,
        err_msg=(
            "the two schedules diverged under all-reflective BCs — an "
            "OUTFLOW SHED ORDER defect in one of them (the ERR-056 class: "
            "observable only on a fully-reflected configuration)."
        ),
    )


# ═══════════════════════════════════════════════════════════════════
# G6 — end-to-end eigenvalue with ScanMarch forced
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.verifies("transport-cartesian-2d", "matrix-eigenvalue", "multigroup")
def test_g6_eigenvalue_hits_kinf_window():
    """[G6] ``solve_sn`` (default SI inner) with the WINDOW forced hits the
    closed-form k_inf — the structural-independence anchor (≥2G; a 1G
    eigenvalue is flux-shape-independent and proves nothing, cardinal rule).

    Post-flip polarity: the ScanMarch default is pinned to k_inf by the
    standing ``test_keff_2d`` suite (which now runs it as the default); THIS
    test keeps the closed-form anchor on the non-default window peer.
    """
    case = get("sn_slab_2eg_1rg")
    mix = next(iter(case.materials.values()))
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 1.0, 3),
        edges_y=np.linspace(0.0, 1.0, 3),
        mat_map=np.zeros((2, 2), dtype=int),
    )
    with window_forced() as hits:
        sol = solve_sn(
            {0: mix}, mesh, Quadrature.level_symmetric(sn_order=4),
            keff_tol=1e-12, flux_tol=1e-10, max_inner=500, inner_tol=1e-10,
        )
    _require_ran(hits, "sweep_interior")
    assert np.isfinite(sol.keff)
    assert abs(sol.keff - case.k_inf) < 1e-8, (
        f"window-forced SI keff={sol.keff:.10f} vs closed-form "
        f"k_inf={case.k_inf:.10f}"
    )


@pytest.mark.slow
@pytest.mark.verifies("transport-cartesian-2d", "matrix-eigenvalue", "multigroup")
def test_g6_si_krylov_heterogeneous_window():
    """[G6] SI ≡ Krylov with the WINDOW forced on the heterogeneous non-flat
    2G config — the wavefront drives BOTH inner methods (SI through
    ``_sweep_interior``, Krylov through the ``loss_action`` matvec) and they
    must converge to the same eigenpair (anchored by the k_inf leg above).
    The ScanMarch default's SI≡Krylov is the standing test_keff_2d gate.
    """
    materials = {2: get_mixture("A", "2g"), 0: get_mixture("B", "2g")}
    nx, ny = 8, 4
    mat = np.zeros((nx, ny), dtype=int)
    mat[:4, :] = 2
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, nx * 0.25, nx + 1),
        edges_y=np.linspace(0.0, ny * 0.25, ny + 1),
        mat_map=mat,
    )
    quad = Quadrature.level_symmetric(sn_order=4)

    with window_forced() as hits:
        sol_si = solve_sn(
            materials, mesh, quad, inner_solver="source_iteration",
            keff_tol=1e-12, flux_tol=1e-10, max_inner=500, inner_tol=1e-10,
        )
    _require_ran(hits, "sweep_interior")

    with window_forced() as hits:
        sol_kry = solve_sn(
            materials, mesh, quad, inner_solver="krylov",
            keff_tol=1e-12, flux_tol=1e-10, max_inner=4000, inner_tol=1e-8,
        )
    _require_ran(hits, "loss_action")

    phi_si = np.asarray(sol_si.scalar_flux.values, dtype=np.float64)
    phi_kry = np.asarray(sol_kry.scalar_flux.values, dtype=np.float64)
    _assert_nonflat(phi_si)

    assert abs(sol_si.keff - sol_kry.keff) < 1e-7, (
        f"window-forced SI keff={sol_si.keff:.10f} vs "
        f"Krylov keff={sol_kry.keff:.10f}"
    )
    phi_si_n = phi_si / phi_si.mean()
    phi_kry_n = phi_kry / phi_kry.mean()
    np.testing.assert_allclose(
        phi_si_n, phi_kry_n, rtol=1e-6, atol=1e-8,
        err_msg="window-forced SI vs Krylov flux SHAPE diverged",
    )
