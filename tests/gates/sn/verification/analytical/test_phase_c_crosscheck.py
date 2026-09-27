r"""Issue #168 Phase C/D/E — Gate Set 4: the curvilinear discrete-ordinates solve against a semi-analytical reference.

* Gate 4.1 — homogeneous-reflective :math:`k_\infty` recovery by the SN
  eigenvalue solver, against the closed form :math:`k_\infty =
  \nu\Sigma_f/\Sigma_a`.
* Gate 4.2 — the SN eigenvalue (Phase D) and flux shape (Phase E) against the
  trajectory-resolvent Variant α Green's-function solvers, called bare (the
  Billiard facade does not route the multi-region variants; GH #190).

The rows, with the claim each makes:

- ``sphere_2g_homogeneous_dd_n20`` and the two ``cyl_1g_homogeneous`` rows:
  on a uniform reflective medium both methods return :math:`k_\infty` exactly
  (the V_α1 / V_α1_cyl identities), so :math:`k` agrees to ``1e-9``. These are
  the edge rows the heterogeneous rows rest on. They are blind to the angular
  closure: a homogeneous medium's angular flux is near-flat, which nulls the
  redistribution the closure feeds (``[M]`` the cylinder rows' flux moves
  1.1e-10 under a deliberate ``tau := 0.7`` mutation, against 8.8e-2 for the
  2-group 3-region cylinder; ``vv-principles`` anti-pattern #3).
- ``sphere_2g_3reg``: the frozen SN snapshot (Gauss-Legendre 8, 40 cells,
  pinned bit-exactly by ``tests/gates/sn/regression/test_dd_regression.py``)
  against the reference at :math:`(n_r, n_\mu) = (36, 96)`, eigenvalue and
  flux shape, with bounds derived below from both methods' measured ladders.
- ``cyl_2g_3reg``: a LIVE SN solve at folded 16x32, 40 cells, against the
  reference at :math:`(n_r, n_{\mu,\rm axial}, n_\varphi) = (24, 16, 32)`. The
  rows that carry the bound are strict ``xfail`` on #516: the reference's
  azimuthal error is larger than a tenth of any bound that would verify the
  SN solve. A RECORD row pins today's readings.

What changed on 2026-09-26 (ERR-090, and two defects of this file):

1. The reference fitted one cubic spline to the emission density across the
   material interfaces; it now fits one per region. Its eigenvalue moved by
   up to 2 % (cylinder 1.20693 to 1.23104; sphere 1.35808 to 1.38374 at the
   old resolutions), and it converges in :math:`n_r`.
2. The Phase D sphere row compared the reference with a hand-typed
   :math:`k = 1.3578153` that had been stale since e30d8d14 (2026-05-12, the
   Phase F snapshot regeneration: the snapshot holds 1.3816447). It agreed at
   2e-4 only because the defective reference happened to read 1.35808. The SN
   eigenvalue is now read from the snapshot file.
3. The Phase E rows evaluated the reference at uniformly spaced cell centres,
   but the snapshot meshes are equal-volume (sphere) and equal-area (cylinder)
   (``RegionMesh``'s default since b5e85c2d). The reference is now
   compared as volume-weighted cell averages over the SN's own cells, read
   from the mesh that produced the SN flux.

The flux-shape metric changed with item 3. Each profile was normalised to its
own maximum per group, which discards the ratio between the groups and
amplifies an error in the maximum's cell. Both profiles are now scaled to unit
total fission production, one gauge for both groups, so the group ratio is
part of the claim; the metric is the largest cell-average difference relative
to the largest reference cell average.

trajectory_resolvent is the **semi-analytical pillar** (``vv-principles``,
the three pillars): chords integrated with scipy/numpy quadrature along
characteristics, sharing no project primitive with the SN sweep above the
trusted-library line. ORPHEUS SN is the production discretisation under test.
"""
from __future__ import annotations

import functools

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_xs
from orpheus.derivations.common.eigenvalue import kinf_homogeneous
from orpheus.derivations.continuous.trajectory_resolvent.greens_function import (
    solve_greens_function_sphere_mg,
    solve_greens_function_sphere_mr,
)
from orpheus.derivations.continuous.trajectory_resolvent.greens_function_cylinder import (
    solve_greens_function_cylinder,
)
from orpheus.geometry import (
    BC,
    CoordSystem,
    Mesh1D,
)
from orpheus.derivations.continuous.trajectory_resolvent.chord_oracle import (
    _regionwise_cubic_spline,
)
from tests.gates.derivations._trajectory_resolvent_ladders import (
    CYLINDER_3REG_SN_STEPS,
    SPHERE_3REG_SN_STEPS,
    sn_residual,
    sphere_3reg_reference_bound,
    tolerance_for,
)
from tests.gates.sn.verification.analytical._certified_agreement import (
    ABA_RADII,
    CYLINDER_3REG_RECORD,
    CYLINDER_3REG_RECORD_BAND,
    CYLINDER_3REG_RECORD_RELATIVE,
    CYLINDER_3REG_REFERENCE_BOUND,
    aba_xs_2g,
    assert_record,
    awaits_cylinder_bound,
    certify_agreement,
    cylinder_3reg_reference,
)


def _make_2g_mixture(
    sigma_t,
    sigma_s_matrix,
    nu_sigma_f,
    chi,
):
    """Build a 2-group Mixture from explicit XS arrays."""
    from orpheus.derivations.common.xs_library import make_mixture
    sigma_t = np.asarray(sigma_t, dtype=float)
    sig_s = np.asarray(sigma_s_matrix, dtype=float)
    nu_sig_f = np.asarray(nu_sigma_f, dtype=float)
    chi = np.asarray(chi, dtype=float)
    # Derive sig_c + nu + sig_f from sig_t + nu_sig_f + scatter.
    sig_a = sigma_t - sig_s.sum(axis=1)
    nu = np.ones_like(nu_sig_f)
    sig_f = nu_sig_f.copy()
    sig_c = sig_a - sig_f
    return make_mixture(
        sig_t=sigma_t,
        sig_c=sig_c,
        sig_f=sig_f,
        nu=nu,
        chi=chi,
        sig_s=sig_s,
    )


# ═══════════════════════════════════════════════════════════════════════
# Gate 4.1 — k_∞ recovery (homogeneous reflective sphere)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("sn-curvilinear-homogeneous-kinf-recovery")
def test_sn_spherical_homogeneous_kinf_recovery_2g():
    r"""Gate 4.1 — 2G homogeneous reflective sphere recovers analytical k_∞.

    On a homogeneous reflective domain the eigenvalue is shape-
    independent: ``k_∞ = (νΣ_f^T φ) / (Σ_a^T φ)`` with the dominant
    eigenvector of the (Σ_a^{-1} νΣ_f^T) matrix. The closed-form
    reference uses :func:`kinf_homogeneous`.

    Phase C: per plan §5 Gate 4.1, target rtol ≤ 5e-4. Eigenvalues
    are essentially exact on uniform reflective curvilinear (the
    only mode is the homogeneous-reflective flat eigenmode).
    """
    from orpheus.sn import solve_sn
    from orpheus.numerics.quadrature import Quadrature

    # Simple 2G test material (no upscatter, no fission spectrum split).
    sigma_t = [0.5, 1.0]
    sig_s = [[0.3, 0.05], [0.0, 0.7]]  # (g_from, g_to)
    nu_sig_f = [0.4, 0.6]
    chi = [1.0, 0.0]
    mat = _make_2g_mixture(sigma_t, sig_s, nu_sig_f, chi)

    k_analytical = kinf_homogeneous(
        np.asarray(sigma_t),
        np.asarray(sig_s),
        np.asarray(nu_sig_f),
        np.asarray(chi),
    )

    nx = 20
    mesh = Mesh1D(
        edges=np.linspace(0.0, 2.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("reflective"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    result = solve_sn(
        materials={0: mat},
        mesh=mesh,
        quadrature=quad,
        max_outer=200, keff_tol=1e-9, flux_tol=1e-8,
        max_inner=200, inner_tol=1e-10,
    )
    keff_sn = result.outcome.keff
    rel = abs(keff_sn - k_analytical) / k_analytical
    print(f"k_analytical={k_analytical:.10f}, k_sn={keff_sn:.10f}, rel={rel:.2e}")
    assert rel < 5e-4, (
        f"Phase C k_∞ recovery target rtol<5e-4 violated: "
        f"k_analytical={k_analytical:.8f}, k_sn={keff_sn:.8f}, "
        f"rel={rel:.2e}"
    )




# ═══════════════════════════════════════════════════════════════════════
# Gate 4.2 — trajectory_resolvent cross-check (BARE function calls)
# ═══════════════════════════════════════════════════════════════════════

_THIS = "tests/gates/sn/verification/analytical/test_phase_c_crosscheck.py"
_REGIONWISE = (
    "tests/gates/derivations/test_trajectory_resolvent_regionwise_source.py"
    "::test_mr_oracle_first_leg_matches_the_line_integral"
)


def _snapshot(snapshot_id: str):
    r"""The frozen SN regression snapshot ``snapshot_id`` (``scalar_flux`` stored ``(ng, nx)``, ``keff``).

    The SN side of a snapshot row is READ from the file that
    ``test_dd_regression`` pins, never typed here: a hand-typed copy of the
    sphere eigenvalue went stale for four months (see the module docstring).
    """
    from tests.gates.sn._test_helpers import SN_TESTS_ROOT
    path = SN_TESTS_ROOT / "regression" / "snapshots" / f"{snapshot_id}.npz"
    if not path.exists():
        pytest.skip(f"snapshot {snapshot_id!r} not present at {path}")
    return np.load(path)


def _run_sphere_2g_homogeneous_closed() -> float:
    """Bare ``solve_greens_function_sphere_mg``, uniform A, R = 2 cm."""
    A2 = get_xs("A", "2g")
    res = solve_greens_function_sphere_mg(
        R=2.0,
        sigma_t=A2["sig_t"],
        sigma_s=A2["sig_s"],
        nu_sigma_f=A2["nu"] * A2["sig_f"],
        chi=A2["chi"],
        alpha=1.0,                                    # closed sphere
        n_r=16, n_mu=16, n_traj_quad=32,              # V_α1 exact at α=1
        max_iter=20, tol=1e-10,
    )
    return float(res.k_eff)


def _run_cyl_1g_homogeneous_closed() -> float:
    """Bare ``solve_greens_function_cylinder``, uniform A, R = 2 cm."""
    A1 = get_xs("A", "1g")
    res = solve_greens_function_cylinder(
        R=2.0,
        sigma_t=float(A1["sig_t"][0]),
        sigma_s=float(A1["sig_s"][0, 0]),
        nu_sigma_f=float(A1["nu"][0] * A1["sig_f"][0]),
        alpha=1.0,                                    # closed cylinder
        n_r=16, n_mu_axial=12, n_phi_az=24, n_traj_quad=32,
        max_iter=20, tol=1e-10,
    )
    return float(res.k_eff)


# (snapshot_id, runner, rtol, rationale): the homogeneous edge rows. On a
# uniform reflective medium k = k_inf = νΣ_f/Σ_a exactly for both methods
# (V_α1 / V_α1_cyl, derive_T00_equals_P_ss_sphere in the SymPy origins), so
# rtol 1e-9 is the two iterations' floor with headroom.
_GATE_4_2_CASES: tuple[tuple[str, object, float, str], ...] = (
    (
        "sphere_2g_homogeneous_dd_n20",
        _run_sphere_2g_homogeneous_closed,
        1e-9,
        "V_α1 algebraic identity — k=k_∞ exact",
    ),
    (
        "cyl_1g_homogeneous_folded_4x8_dd_n20",
        _run_cyl_1g_homogeneous_closed,
        1e-9,
        "V_α1_cyl algebraic identity — k=k_∞ exact",
    ),
    (
        "cyl_1g_homogeneous_folded_2x4_dd_n20",
        _run_cyl_1g_homogeneous_closed,
        1e-9,
        "V_α1_cyl algebraic identity — k=k_∞ exact (folded 2x4 split)",
    ),
)


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("sn-curvilinear-trajectory-resolvent-crosscheck")
@pytest.mark.parametrize(
    "snapshot_id, runner, rtol, rationale",
    _GATE_4_2_CASES,
    ids=[case[0] for case in _GATE_4_2_CASES],
)
def test_phase_d_trajectory_resolvent_crosscheck(
    snapshot_id, runner, rtol, rationale,
) -> None:
    r"""Gate 4.2, the edge rows: SN snapshot k against the trajectory resolvent on a homogeneous medium.

    Both sides must return :math:`k_\infty` of the uniform medium; the
    reference does so by the V_α1 / V_α1_cyl identities
    (:mod:`orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function`).
    These rows are 1-group or homogeneous and so blind to every spatial and
    angular operator (``vv-principles`` anti-pattern #3); they are the
    foundation the heterogeneous rows below rest on, not evidence about them.
    """
    expected_keff = float(_snapshot(snapshot_id)["keff"])
    k_ref = runner()
    rel = abs(k_ref - expected_keff) / expected_keff
    print(
        f"{snapshot_id}: k_sn={expected_keff:.10f}  "
        f"k_trajres={k_ref:.10f}  rel={rel:.2e}  target={rtol:.0e}  "
        f"({rationale})"
    )
    assert rel < rtol, (
        f"Gate 4.2 cross-check for {snapshot_id!r} exceeded tolerance: "
        f"k_sn_snapshot={expected_keff:.8f}, "
        f"k_trajectory_resolvent={k_ref:.8f}, rel={rel:.2e}, "
        f"target rtol={rtol:.0e}. Rationale: {rationale}"
    )


# ════════════════════════════════════════════════════════════════════════
# The heterogeneous rows — eigenvalue and flux shape
# ════════════════════════════════════════════════════════════════════════
#
# The observable of the shape rows is the scalar flux as CELL AVERAGES over
# the SN solve's own cells (the mesh object that produced the SN flux), with
# both profiles scaled to unit total fission production. The reference is a
# nodal field on composite per-region Gauss-Legendre nodes; its cell average
# is the volume-weighted integral of a per-region cubic spline through those
# nodes (the flux is continuous across an interface but its derivative is
# not, so the interpolant is per region, as the reference's own emission


def _reference_cell_averages(result, mesh: Mesh1D) -> np.ndarray:
    r"""The reference's scalar flux averaged over each SN cell, ``(G, N)``.

    :math:`\bar\phi_{g,i} = \int_{r_i}^{r_{i+1}} \phi_g(r)\, r^p\,dr \big/
    \int_{r_i}^{r_{i+1}} r^p\,dr` with :math:`p = 2` (sphere) or 1 (cylinder),
    by 16-point Gauss-Legendre per cell on the reference's regionwise cubic
    spline (the one the reference reads its own emission density through).
    Each SN cell lies in one region (the region meshes are built per region).
    """
    p = 2 if mesh.coord is CoordSystem.SPHERICAL else 1
    x, w = np.polynomial.legendre.leggauss(16)
    edges = np.asarray(mesh.edges, dtype=float)
    region_of_cell = np.searchsorted(ABA_RADII, 0.5 * (edges[1:] + edges[:-1]))
    phi = np.asarray(result.phi_g, dtype=float)
    averages = np.empty((phi.shape[0], len(edges) - 1))
    for g in range(phi.shape[0]):
        pieces = _regionwise_cubic_spline(result.r_nodes, phi[g], result.region_at_node, len(ABA_RADII))
        for i, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
            r = 0.5 * (b - a) * (x + 1.0) + a
            weight = 0.5 * (b - a) * w * r**p
            averages[g, i] = np.sum(weight * pieces[region_of_cell[i]](r)) / np.sum(weight)
    return averages


def _fission_gauged(phi: np.ndarray, mesh: Mesh1D) -> np.ndarray:
    r"""``phi`` ``(G, N)`` scaled so that :math:`\sum_{g,i} \nu\Sigma_{f,g,i}\,\phi_{g,i} V_i = 1`."""
    centres = 0.5 * (np.asarray(mesh.edges)[1:] + np.asarray(mesh.edges)[:-1])
    nu_sigma_f = aba_xs_2g()[2][np.searchsorted(ABA_RADII, centres)].T  # (G, N)
    production = float(np.sum(nu_sigma_f * phi * np.asarray(mesh.volumes)))
    return phi / production


def _shape_gap(phi_sn: np.ndarray, reference_result, mesh: Mesh1D) -> np.ndarray:
    r"""Per group: :math:`\max_i |\hat\phi^{\rm SN}_{g,i} - \hat\phi^{\rm ref}_{g,i}| / \max_{g,i} \hat\phi^{\rm ref}_{g,i}`, fission-gauged cell averages."""
    sn = _fission_gauged(np.asarray(phi_sn, dtype=float), mesh)
    ref = _fission_gauged(_reference_cell_averages(reference_result, mesh), mesh)
    return np.max(np.abs(sn - ref), axis=1) / float(np.max(ref))


@functools.cache
def _sphere_3reg_reference():
    """The sphere reference at (n_r, n_mu) = (36, 96), solved once per session."""
    sigma_t, sigma_s, nu_sigma_f, chi = aba_xs_2g()
    return solve_greens_function_sphere_mr(
        radii=ABA_RADII, sigma_t=sigma_t, sigma_s=sigma_s,
        nu_sigma_f=nu_sigma_f, chi=chi, alpha=1.0,
        n_r=36, n_mu=96, n_traj_quad=64,
        max_iter=2000, tol=1e-9, initial_k=1.38,
    )


@functools.cache
def _cylinder_3reg_sn_16x32():
    r"""A live SN solve of the snapshot's cylinder problem at folded 16x32: ``(k, phi (G, N), mesh)``.

    The problem is ``cyl_2g_3reg_folded_4x8_dd_n40``'s (same materials, same
    40-cell equal-area mesh) with the angular grid refined from 4x8 to 16x32
    and the inner budget raised so the inner iteration converges (the
    snapshot's own 300 exits best-effort). The 4x8 snapshot stays pinned
    bit-exactly by ``test_dd_regression``, whose τ sensitivity (8.8e-2 in the
    flux under ``tau := 0.7``) makes it the cylinder angular-closure
    catcher; at 4x8 the SN solve's own angular error in this metric is about
    0.13, larger than any bound this comparison can certify, so the
    cross-check is posed where SN is angularly resolved.
    """
    from orpheus.numerics.quadrature import Quadrature
    from tests.gates.sn.regression import _generate_snapshots as generator
    config = {
        **generator._cylinder_3region("2g", 40, "folded_4x8"),
        "quadrature": Quadrature.folded_product(n_mu=16, n_phi=32),
        "max_inner": 2000,
    }
    result = generator.run_case(config)
    return (
        float(result.outcome.keff),
        np.asarray(result.scalar_flux.values, dtype=float),
        config["mesh"],
    )


@functools.cache
def _sphere_3reg_sn_gl32():
    r"""A live SN solve of the snapshot's sphere problem at Gauss-Legendre 32: ``(k, phi (G, N), mesh)``.

    ``sphere_2g_3reg_dd_n40``'s problem (same materials, same 40-cell
    equal-volume mesh) with 32 ordinates instead of 8 and an inner budget
    the inner iteration converges within. The 8-ordinate snapshot stays
    pinned bit-exactly by ``test_dd_regression``; its own angular error in
    the shape metric (1.3e-2) would consume a bound this comparison can
    otherwise certify ten times tighter.
    """
    from orpheus.numerics.quadrature import Quadrature
    from tests.gates.sn.regression import _generate_snapshots as generator
    config = {
        **generator._sphere_3region("2g", 40),
        "quadrature": Quadrature.gauss_legendre(n_ordinates=32),
        "max_inner": 2000,
    }
    result = generator.run_case(config)
    return (
        float(result.outcome.keff),
        np.asarray(result.scalar_flux.values, dtype=float),
        config["mesh"],
    )




# ── the sphere: certified ─────────────────────────────────────────────────
#
# Bounds and tolerances are COMPUTED from the ladders in
# tests/gates/derivations/_trajectory_resolvent_ladders.py (which also holds
# the command that re-measures each): the reference's bound at (36, 96) is its
# radial step extrapolated at the measured second order plus its alternating
# mu step; the SN residual is the sum of its mesh and angular steps; the
# tolerance is tolerance_for(SN residual, reference bound). [M] 2026-09-26:
# reference bound 3.2e-4 (k) and 1.4e-3 (shape); SN residual 1.5e-5 (k) and
# 7.3e-3 (shape); tolerances 4e-3 and 2e-2.
_SPHERE_REFERENCE_BOUND = sphere_3reg_reference_bound()
_SPHERE_TOLERANCE = {
    observable: tolerance_for(sn_residual(SPHERE_3REG_SN_STEPS[observable]), _SPHERE_REFERENCE_BOUND[observable])
    for observable in ("k", "shape")
}

_SPHERE_SUPPORTS = (
    f"{_REGIONWISE}[sphere]",
    "tests/gates/derivations/test_peierls_greens_function_mr.py::test_mr_sphere_k_converges_in_n_r",
    f"{_THIS}::test_phase_d_trajectory_resolvent_crosscheck[sphere_2g_homogeneous_dd_n20]",
    f"{_THIS}::test_sn_spherical_homogeneous_kinf_recovery_2g",
)


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("sn-curvilinear-trajectory-resolvent-crosscheck")
@pytest.mark.rests_on(*_SPHERE_SUPPORTS)
def test_sphere_3reg_k_against_trajectory_resolvent() -> None:
    r"""The heterogeneous closed sphere's eigenvalue: SN at GL32 against the reference at (36, 96).

    Fuel A | moderator B | fuel A at 0.5, 1.5, 2.0 cm, 2 groups, reflective
    at r = R. The bound, derived above (:data:`_SPHERE_TOLERANCE`), is the
    reference's: 4e-3 against the SN solve's own residual of 1.5e-5.
    ``[M]`` reading 6.5e-5; the one-spline reference (ERR-090) reads 7.9e-3.

    Blind to the SN angular-closure defect class: ``tau := 0.7`` reads 8.7e-4
    here, inside the reference-limited bound; that defect's catcher is
    ``tests/gates/sn/regression/test_dd_regression.py``. The SN-defect
    witness of this row is a wrong boundary law (the reflective face realised
    as vacuum reads 2.6e1).

    Until 2026-09-26 this row compared the reference at (24, 24) with a
    hand-typed k four months stale (1.3578153, against the snapshot's
    1.3816447) and read 2e-4 under a 2e-2 bound: the stale number and the
    one-spline reference had drifted to the same value.
    """
    k_sn, _, _ = _sphere_3reg_sn_gl32()
    k_ref = float(_sphere_3reg_reference().k_eff)
    reading = abs(k_ref - k_sn) / k_sn
    print(f"sphere 3-region: k_sn={k_sn:.8f} k_ref={k_ref:.8f} rel={reading:.3e}")
    certify_agreement("sphere k", reading, _SPHERE_TOLERANCE["k"], _SPHERE_REFERENCE_BOUND["k"]).require()


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("sn-curvilinear-trajectory-resolvent-crosscheck")
@pytest.mark.rests_on(*_SPHERE_SUPPORTS)
def test_sphere_3reg_flux_shape_against_trajectory_resolvent() -> None:
    r"""The heterogeneous closed sphere's flux shape, fission-gauged cell averages over the SN mesh.

    Activates the spatial and angular redistribution in both groups and the
    group ratio (one gauge for both groups); the material interfaces are
    inside the comparison. Bound 2e-2 (derived above), set by the
    reference's 1.4e-3 and the SN solve's own 7.3e-3. ``[M]`` reading
    4.4e-3 (fast) / 1.7e-3 (thermal). Blind to ``tau := 0.7`` (1.0e-2) like
    the eigenvalue row; red under the vacuum-for-reflective law (6.1) and a
    1 % perturbation of the reference's moderator density (7.2e-2), green
    under the one-spline reference (7.6e-3).
    """
    _, phi_sn, mesh = _sphere_3reg_sn_gl32()
    gap = _shape_gap(phi_sn, _sphere_3reg_reference(), mesh)
    print(f"sphere 3-region shape gap per group: {gap}")
    certify_agreement(
        "sphere flux shape", float(gap.max()), _SPHERE_TOLERANCE["shape"], _SPHERE_REFERENCE_BOUND["shape"],
    ).require()


# ── the cylinder: not yet certifiable (#516) ─────────────────────────────
#
# The tolerances the cylinder rows are held to once the reference is
# certified: tolerance_for(SN residual, None), the reference assumed at the
# floor. [M] SN residual at folded 16x32, 40 cells: 3.0e-5 (k), 5.0e-3
# (shape); tolerances 8e-5 and 2e-2. The reference carries no bound
# (_certified_agreement.CYLINDER_3REG_REFERENCE_BOUND), so the floor fails.
_CYLINDER_TOLERANCE = {
    observable: tolerance_for(sn_residual(CYLINDER_3REG_SN_STEPS[observable]), None)
    for observable in ("k", "shape")
}

_CYLINDER_SUPPORTS = (
    f"{_REGIONWISE}[cylinder]",
    f"{_THIS}::test_phase_d_trajectory_resolvent_crosscheck[cyl_1g_homogeneous_folded_4x8_dd_n20]",
    "tests/gates/derivations/test_peierls_greens_function_cylinder_mr.py::test_mr_K3_uniform_reduces_to_mg_2g",
    "tests/gates/derivations/test_peierls_greens_function_cylinder_mr_xverif.py::test_mr_single_region_vacuum_matches_wm72",
)


def _cylinder_readings() -> dict[str, float]:
    """Today's cylinder readings, keyed as :data:`CYLINDER_3REG_RECORD` keys them."""
    k_sn, phi_sn, mesh = _cylinder_3reg_sn_16x32()
    reference = cylinder_3reg_reference()
    gap = _shape_gap(phi_sn, reference, mesh)
    return {
        "k_ref": float(reference.k_eff),
        "phase_c_k_sn": k_sn,
        "phase_c_k_gap": abs(float(reference.k_eff) - k_sn) / k_sn,
        "phase_c_shape_fast": float(gap[0]),
        "phase_c_shape_thermal": float(gap[1]),
    }


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("sn-curvilinear-trajectory-resolvent-crosscheck")
@pytest.mark.rests_on(*_CYLINDER_SUPPORTS)
@awaits_cylinder_bound
def test_cylinder_3reg_k_against_trajectory_resolvent() -> None:
    r"""The heterogeneous closed cylinder's eigenvalue: live SN at folded 16x32 against the reference at (24, 16, 32).

    Held to 8e-5 (derived above), which needs a reference bound of at most
    8e-6. The cylinder reference carries none (its ladder is not monotone in
    the azimuthal order), so the floor fails first: the expected failure.
    The row XPASSes when a bound derived from a converging ladder replaces
    ``CYLINDER_3REG_REFERENCE_BOUND["k"]`` (after #516), not when a repair
    alone lands; the RECORD row below keeps the reading live meanwhile.
    """
    readings = _cylinder_readings()
    print(f"cylinder 3-region k: {readings}")
    certify_agreement(
        "cylinder k", readings["phase_c_k_gap"], _CYLINDER_TOLERANCE["k"], CYLINDER_3REG_REFERENCE_BOUND["k"],
    ).require()


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("sn-curvilinear-trajectory-resolvent-crosscheck")
@pytest.mark.rests_on(*_CYLINDER_SUPPORTS)
@awaits_cylinder_bound
def test_cylinder_3reg_flux_shape_against_trajectory_resolvent() -> None:
    r"""The heterogeneous closed cylinder's flux shape, fission-gauged cell averages over the SN mesh.

    Re-posed on 2026-09-26 from the 4x8 snapshot onto a live 16x32 solve: at
    4x8 SN's own angular error in this metric is 6.1e-2, which no comparison
    can separate from a defect; the snapshot keeps its job in
    ``test_dd_regression``. Held to 2e-2 (derived above); the expected
    failure is the eigenvalue row's.
    """
    readings = _cylinder_readings()
    print(f"cylinder 3-region shape: {readings}")
    certify_agreement(
        "cylinder flux shape", max(readings["phase_c_shape_fast"], readings["phase_c_shape_thermal"]),
        _CYLINDER_TOLERANCE["shape"], CYLINDER_3REG_REFERENCE_BOUND["shape"],
    ).require()


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.rests_on(*_CYLINDER_SUPPORTS)
def test_cylinder_3reg_crosscheck_record() -> None:
    r"""RECORD: today's cylinder readings, green until either side moves.

    Not verification. Red under ``tau := 0.7`` (the SN k moves 7.8e-4), the
    vacuum-for-reflective law, the one-spline reference and a 1 % reference
    perturbation.
    """
    readings = _cylinder_readings()
    print(f"cylinder 3-region readings: {readings}")
    assert_record(
        readings, {name: CYLINDER_3REG_RECORD[name] for name in readings},
        CYLINDER_3REG_RECORD_BAND, relative=CYLINDER_3REG_RECORD_RELATIVE,
    )
