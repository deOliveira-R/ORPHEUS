r"""Plan-2 follow-on A3 — Multi-group Variant α extension test gate.

Closes the 1G Cardinal-Rule-6 gap from the V&V audit. The Plan 2
Part B prototype was verified only at 1G, where k = νΣ_f/Σ_a is
flux-shape independent — angular, scattering-matrix, and group-coupling
errors are invisible. A3 Option A extends to ≥2G with arbitrary
scattering matrix and tests:

1. **G = 1 special case bit-equivalence** — MG solver with G=1 inputs
   reproduces the existing 1G solver to machine precision. Sanity
   check that the MG refactor didn't break the 1G path.
2. **Closed sphere k_eff = k_inf at 2G with downscatter** — MG analog
   of V_α1.numerical. Closed sphere has no leakage; the multi-group
   eigenmode is the infinite-medium eigenmode with k = k_inf
   (transfer-matrix dominant eigenvalue) and group flux spectrum
   matching the dominant right eigenvector.
3. **Closed sphere k_eff = k_inf at 2G with upscatter** —
   exercises the more general MG iteration (without the lower-
   triangular shortcut). Same identity holds.
4. **Vacuum 2G sphere k_eff < k_inf** — leakage reduces effective
   multiplication. Spatial mode non-trivial in both groups.
5. **2G group-flux spectrum** — for closed sphere, the per-group
   scalar flux ratio matches the analytical downscatter relation
   :math:`\phi_2/\phi_1 = \Sigma_{s,1\to 2}/\Sigma_{a,2}`.
6. **Reciprocity sanity** — k_eff at G=1 with a synthetic
   "single-group equivalent" XS matches the 2G result when the
   1G XS is correctly homogenised. (Skipped for now; relies on
   exact homogenisation theory.)

Predecessors:

- :mod:`.test_peierls_greens_function_solver` (V_α1.numerical 1G)
- :mod:`.test_peierls_greens_function_vacuum` (vacuum BC 1G)
- :mod:`.test_peierls_greens_function_xverif_ps1982` (PS-1982 1G)

Reference for closed sphere k_inf:
:func:`orpheus.derivations.common.eigenvalue.kinf_homogeneous` and
:func:`kinf_and_spectrum_homogeneous` (transfer-matrix dominant
eigenvalue + spectrum).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.eigenvalue import (
    kinf_and_spectrum_homogeneous,
)
from orpheus.derivations.common.xs_library import get_xs
from orpheus.derivations.continuous.peierls.greens_function import (
    solve_greens_function_sphere,
    solve_greens_function_sphere_mg,
)


# ════════════════════════════════════════════════════════════════════════
# Fixtures: 1G fuel-A-like + 2G test cases
# ════════════════════════════════════════════════════════════════════════


@pytest.fixture(scope="module")
def fuelA_1G():
    """Same fuel-A-like 1G XS used by 1G prototype tests."""
    return {
        "sigma_t": np.array([0.5]),
        "sigma_s": np.array([[0.38]]),
        "nu_sigma_f": np.array([0.025]),
        "chi": np.array([1.0]),
    }


@pytest.fixture(scope="module")
def two_group_downscatter():
    """2G downscatter-only test case.

    Group 1 (fast):    σ_t = 0.6,   σ_{s,1→1} = 0.25, σ_{s,1→2} = 0.15,
                       νσ_{f,1} = 0.05
    Group 2 (thermal): σ_t = 1.0,   σ_{s,2→2} = 0.7,  νσ_{f,2} = 0.15
    χ = (1, 0): all-fast emission.

    Removal cross-sections (Σ_R,g = Σ_t,g − Σ_{s,g→g}):
        Σ_R,1 = 0.6 − 0.25 = 0.35 (out-scatter from g=1: σ_{s,1→2})
        Σ_R,2 = 1.0 − 0.70 = 0.30 (just absorption since no upscatter)

    Closed-form 2G k_inf with χ = (1, 0), no upscatter:
        k_inf = (νΣ_f,1 · Σ_R,2 + νΣ_f,2 · Σ_{s,1→2})
                / (Σ_R,1 · Σ_R,2)
              = (0.05·0.30 + 0.15·0.15) / (0.35·0.30)
              = 0.0375 / 0.105
              ≈ 0.357

    Subcritical, well-conditioned, with non-trivial group coupling.
    """
    sigma_t = np.array([0.6, 1.0])
    sigma_s = np.array([
        [0.25, 0.15],   # group 1: in-scatter 0.25, downscatter to 2 = 0.15
        [0.0,  0.7],    # group 2: no upscatter, in-scatter 0.7
    ])
    nu_sigma_f = np.array([0.05, 0.15])
    chi = np.array([1.0, 0.0])
    return {
        "sigma_t": sigma_t,
        "sigma_s": sigma_s,
        "nu_sigma_f": nu_sigma_f,
        "chi": chi,
    }


@pytest.fixture(scope="module")
def fuelA_4G():
    """Fuel-A 4G XS from the project library — production-grade test
    fixture with realistic χ spectrum.

    Groups (fast → thermal): σ_t = (0.4, 0.61, 0.85, 1.05),
    full lower-triangular σ_s (4×4 downscatter matrix), χ = (0.6,
    0.35, 0.05, 0) → 60 % of fission neutrons in group 0, 35 % in
    group 1, 5 % in group 2, none in group 3 (thermal).

    Realistic light-water-reactor-like spectrum with non-trivial
    fission distribution. Closed sphere k_inf is computed via
    transfer-matrix.
    """
    xs = get_xs("A", "4g")
    return {
        "sigma_t": xs["sig_t"],
        "sigma_s": xs["sig_s"],
        "nu_sigma_f": xs["nu"] * xs["sig_f"],
        "chi": xs["chi"],
    }


@pytest.fixture(scope="module")
def two_group_with_upscatter():
    """2G case with upscatter — exercises the general MG iteration.

    Same as downscatter case but with a small upscatter term σ_{2→1}.
    k_inf must be recomputed via transfer-matrix.
    """
    sigma_t = np.array([0.6, 1.0])
    sigma_s = np.array([
        [0.25, 0.15],
        [0.05, 0.65],   # upscatter 2→1 at 0.05
    ])
    nu_sigma_f = np.array([0.05, 0.15])
    chi = np.array([1.0, 0.0])
    return {
        "sigma_t": sigma_t,
        "sigma_s": sigma_s,
        "nu_sigma_f": nu_sigma_f,
        "chi": chi,
    }


# ════════════════════════════════════════════════════════════════════════
# Tests
# ════════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_a3_mg_g1_matches_1g_solver(fuelA_1G):
    """A3.1 — MG solver with G=1 inputs reproduces 1G solver to
    machine precision.

    The MG solver is implemented in terms of `_apply_operator_with_
    source_profile`; the 1G solver wraps the same underlying machinery.
    This test pins that the MG → 1G reduction is bit-identical.
    """
    fix = fuelA_1G

    # 1G via existing solver
    res_1g = solve_greens_function_sphere(
        R=5.0,
        sigma_t=float(fix["sigma_t"][0]),
        sigma_s=float(fix["sigma_s"][0, 0]),
        nu_sigma_f=float(fix["nu_sigma_f"][0]),
        alpha=1.0,
        n_r=12, n_mu=12, n_traj_quad=24, max_iter=20, tol=1e-12,
    )

    # Same XS via MG solver with G=1
    res_mg = solve_greens_function_sphere_mg(
        R=5.0,
        sigma_t=fix["sigma_t"],
        sigma_s=fix["sigma_s"],
        nu_sigma_f=fix["nu_sigma_f"],
        chi=fix["chi"],
        alpha=1.0,
        n_r=12, n_mu=12, n_traj_quad=24, max_iter=20, tol=1e-12,
    )

    np.testing.assert_allclose(
        res_mg.k_eff, res_1g.k_eff, rtol=1e-13,
        err_msg=(
            f"A3.1: MG (G=1) k_eff = {res_mg.k_eff} differs from 1G "
            f"k_eff = {res_1g.k_eff} beyond machine precision"
        ),
    )


@pytest.mark.foundation
def test_a3_mg_2g_closed_sphere_equals_k_inf_downscatter(
    two_group_downscatter,
):
    """A3.2 — 2G closed sphere k_eff = k_inf (transfer-matrix dominant
    eigenvalue), with downscatter only.

    Multi-group analog of V_α1.numerical: closed sphere has no
    leakage, eigenmode is the infinite-medium spectrum, k_eff =
    k_inf. The transfer-matrix k_inf is the dominant eigenvalue of
    A^{-1}·F where A = diag(σ_t) - σ_s.T and F = outer(χ, νσ_f).

    Cardinal-Rule-6 closure: this test exercises non-trivial
    multi-group coupling (downscatter, fission spectrum) which 1G
    cannot detect.
    """
    fix = two_group_downscatter

    k_inf_ref, phi_spectrum_ref = kinf_and_spectrum_homogeneous(
        fix["sigma_t"], fix["sigma_s"], fix["nu_sigma_f"], fix["chi"],
    )

    res = solve_greens_function_sphere_mg(
        R=5.0, **{k: fix[k] for k in ("sigma_t", "sigma_s", "nu_sigma_f", "chi")},
        alpha=1.0,
        n_r=12, n_mu=12, n_traj_quad=24, max_iter=500, tol=1e-12,
    )

    assert res.converged, (
        f"A3.2: MG iteration did not converge in 500 iter; "
        f"k_eff = {res.k_eff:.10f}"
    )
    np.testing.assert_allclose(
        res.k_eff, k_inf_ref, rtol=1e-9,
        err_msg=(
            f"A3.2: Variant α MG closed sphere k_eff = {res.k_eff} "
            f"differs from k_inf = {k_inf_ref} (transfer matrix)"
        ),
    )

    # Group flux spectrum check — closed sphere has uniform per-group
    # flux equal to the dominant transfer-matrix eigenvector
    # (sign-normalised, ℓ²-normalised).
    phi_per_group = res.phi_g.mean(axis=1)  # (G,) — average over r
    phi_per_group_normed = phi_per_group / np.linalg.norm(phi_per_group)
    np.testing.assert_allclose(
        phi_per_group_normed, phi_spectrum_ref, rtol=1e-6,
        err_msg=(
            f"A3.2 spectrum: variant α gives normalised φ = "
            f"{phi_per_group_normed}, transfer matrix gives "
            f"{phi_spectrum_ref}"
        ),
    )


@pytest.mark.foundation
def test_a3_mg_2g_closed_sphere_equals_k_inf_upscatter(
    two_group_with_upscatter,
):
    """A3.3 — same as A3.2 but with upscatter present.

    Tests the general MG iteration path (Σ_s is not lower-triangular).
    The transfer-matrix solve is the reference regardless of the
    scatter pattern.
    """
    fix = two_group_with_upscatter

    k_inf_ref, _ = kinf_and_spectrum_homogeneous(
        fix["sigma_t"], fix["sigma_s"], fix["nu_sigma_f"], fix["chi"],
    )

    res = solve_greens_function_sphere_mg(
        R=5.0, **{k: fix[k] for k in ("sigma_t", "sigma_s", "nu_sigma_f", "chi")},
        alpha=1.0,
        n_r=12, n_mu=12, n_traj_quad=24, max_iter=500, tol=1e-12,
    )

    assert res.converged
    np.testing.assert_allclose(
        res.k_eff, k_inf_ref, rtol=1e-9,
        err_msg=(
            f"A3.3 (upscatter): k_eff = {res.k_eff} ≠ k_inf = "
            f"{k_inf_ref}"
        ),
    )


@pytest.mark.foundation
def test_a3_mg_2g_vacuum_below_k_inf(two_group_downscatter):
    """A3.4 — vacuum 2G sphere k_eff < k_inf (leakage), spatial modes
    non-trivial in both groups.

    Stress-tests the per-group operator on a non-degenerate
    eigenmode. Closed sphere has uniform φ_g; vacuum BC drains both
    groups, with a non-trivial radial profile per group.
    """
    fix = two_group_downscatter

    k_inf, _ = kinf_and_spectrum_homogeneous(
        fix["sigma_t"], fix["sigma_s"], fix["nu_sigma_f"], fix["chi"],
    )

    res = solve_greens_function_sphere_mg(
        R=5.0, **{k: fix[k] for k in ("sigma_t", "sigma_s", "nu_sigma_f", "chi")},
        alpha=0.0,
        n_r=24, n_mu=24, n_traj_quad=48, max_iter=500, tol=1e-9,
    )
    assert res.converged

    assert res.k_eff < k_inf, (
        f"A3.4: vacuum k_eff = {res.k_eff} should be < k_inf = "
        f"{k_inf}"
    )
    # Non-trivial leakage: less than 0.95 · k_inf for τ_R ≈ 3.
    assert res.k_eff < 0.95 * k_inf, (
        f"A3.4: vacuum should have substantial leakage; got "
        f"ratio = {res.k_eff/k_inf:.4f}"
    )

    # Per-group spatial modes are non-trivial.
    for g in range(2):
        rel_spread = res.phi_g[g].std() / res.phi_g[g].mean()
        assert rel_spread > 0.3, (
            f"A3.4: vacuum group-{g} φ should be spatially non-trivial; "
            f"got std/mean = {rel_spread:.3f}"
        )

    # Per-group flux at center > flux at surface
    for g in range(2):
        assert res.phi_g[g, 0] > res.phi_g[g, -1], (
            f"A3.4: vacuum group-{g} φ should peak at centre; got "
            f"φ(centre) = {res.phi_g[g, 0]:.4e}, "
            f"φ(surface) = {res.phi_g[g, -1]:.4e}"
        )


@pytest.mark.foundation
def test_a3_mg_2g_downscatter_flux_ratio_closed(two_group_downscatter):
    r"""A3.5 — closed-sphere 2G downscatter-only flux ratio satisfies
    the analytical balance.

    For closed sphere with downscatter only, χ = (1, 0):

    .. math::

       \Sigma_{a,2}\,\phi_2 \;=\; \Sigma_{s,1\to 2}\,\phi_1
       \quad\Longrightarrow\quad
       \frac{\phi_2}{\phi_1} \;=\; \frac{\Sigma_{s,1\to 2}}
                                          {\Sigma_{a,2}}.

    For the fixture: :math:`\phi_2/\phi_1 = 0.15 / 0.30 = 0.5`.

    This is a strict algebraic identity that the MG iteration must
    reproduce — independent of the radial geometry (no leakage).
    """
    fix = two_group_downscatter

    res = solve_greens_function_sphere_mg(
        R=5.0, **{k: fix[k] for k in ("sigma_t", "sigma_s", "nu_sigma_f", "chi")},
        alpha=1.0,
        n_r=12, n_mu=12, n_traj_quad=24, max_iter=500, tol=1e-12,
    )
    assert res.converged

    expected_ratio = fix["sigma_s"][0, 1] / (
        fix["sigma_t"][1] - fix["sigma_s"][1, 1]
    )
    actual_ratio = res.phi_g[1].mean() / res.phi_g[0].mean()
    np.testing.assert_allclose(
        actual_ratio, expected_ratio, rtol=1e-9,
        err_msg=(
            f"A3.5: φ_2/φ_1 = {actual_ratio:.6f} ≠ analytical "
            f"Σ_{{1→2}}/Σ_a,2 = {expected_ratio:.6f}"
        ),
    )


# ════════════════════════════════════════════════════════════════════════
# A3 Option B — 4G with non-trivial χ spectrum (production-grade)
# ════════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_a3_mg_4g_closed_sphere_equals_k_inf(fuelA_4G):
    r"""A3.B1 — 4G fuel-A closed sphere k_eff = k_inf with realistic
    χ = (0.6, 0.35, 0.05, 0.0) fission spectrum.

    Production-grade test: 4 groups, full downscatter matrix, non-
    trivial χ that distributes fission across multiple groups. The
    transfer-matrix k_inf is the dominant eigenvalue of A^{-1}·F where
    F = outer(χ, νΣ_f) has rank > 1 (unlike the (1, 0) special case
    where F is rank 1).

    Variant α MG must reproduce both the k_eff and the dominant
    flux spectrum φ_g (the right eigenvector of A^{-1}·F).
    """
    fix = fuelA_4G

    k_inf_ref, phi_spectrum_ref = kinf_and_spectrum_homogeneous(
        fix["sigma_t"], fix["sigma_s"], fix["nu_sigma_f"], fix["chi"],
    )

    res = solve_greens_function_sphere_mg(
        R=5.0, **{k: fix[k] for k in ("sigma_t", "sigma_s", "nu_sigma_f", "chi")},
        alpha=1.0,
        n_r=12, n_mu=12, n_traj_quad=24, max_iter=1000, tol=1e-11,
    )

    assert res.converged, (
        f"A3.B1: 4G iteration did not converge in 1000 iter; "
        f"k_eff = {res.k_eff:.10f}"
    )
    np.testing.assert_allclose(
        res.k_eff, k_inf_ref, rtol=1e-8,
        err_msg=(
            f"A3.B1: 4G k_eff = {res.k_eff} differs from k_inf = "
            f"{k_inf_ref}"
        ),
    )

    # Spectrum check — every group's normalised flux matches the
    # transfer-matrix dominant eigenvector.
    phi_per_group = res.phi_g.mean(axis=1)
    phi_per_group_normed = phi_per_group / np.linalg.norm(phi_per_group)
    np.testing.assert_allclose(
        phi_per_group_normed, phi_spectrum_ref, rtol=1e-5,
        err_msg=(
            f"A3.B1 spectrum: variant α gives normalised φ = "
            f"{phi_per_group_normed}, transfer matrix gives "
            f"{phi_spectrum_ref}"
        ),
    )


@pytest.mark.foundation
def test_a3_mg_4g_vacuum_subcritical(fuelA_4G):
    r"""A3.B2 — 4G fuel-A vacuum sphere has k_eff < k_inf, all groups
    spatially non-trivial.

    Stress-tests the per-group operator across 4 groups simultaneously
    with realistic χ. Convergence is slower than 2G (more groups
    coupled) but the spatial mode is well-defined and the k-eigenvalue
    drops below k_inf due to leakage.
    """
    fix = fuelA_4G

    k_inf, _ = kinf_and_spectrum_homogeneous(
        fix["sigma_t"], fix["sigma_s"], fix["nu_sigma_f"], fix["chi"],
    )

    res = solve_greens_function_sphere_mg(
        R=10.0, **{k: fix[k] for k in ("sigma_t", "sigma_s", "nu_sigma_f", "chi")},
        alpha=0.0,
        n_r=24, n_mu=24, n_traj_quad=48, max_iter=1000, tol=1e-7,
    )
    assert res.converged, (
        f"A3.B2: 4G vacuum iteration did not converge in 1000 iter"
    )
    assert res.k_eff < k_inf, (
        f"A3.B2: vacuum 4G k_eff = {res.k_eff} should be < k_inf = "
        f"{k_inf}"
    )

    # All four group fluxes spatially non-trivial.
    for g in range(4):
        rel_spread = res.phi_g[g].std() / res.phi_g[g].mean()
        assert rel_spread > 0.1, (
            f"A3.B2: vacuum group-{g} φ should be non-trivial; got "
            f"std/mean = {rel_spread:.3f}"
        )
