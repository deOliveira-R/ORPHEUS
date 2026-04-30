"""L1 tests for the rank-1 white-BC closure of the Peierls cylinder solver.

C7 of the Phase-4.2 campaign. Exercises the ``boundary='white'``
path of
:func:`~orpheus.derivations.peierls_geometry.solve_peierls_1g`
(with ``geometry=_pg.CYLINDER_1D``) and the two helpers
:func:`~orpheus.derivations.peierls_geometry.compute_P_esc` and
:func:`~orpheus.derivations.peierls_geometry.compute_G_bc`.

.. important::

   The rank-1 white-BC closure is the CP-flat-source-equivalent
   (Mark / isotropic re-entry) closure at the pointwise-Nyström
   level. It is an **approximation** that degrades as the cell
   becomes thinner — the Wigner-Seitz exact identity
   ``k_eff(white) = k_inf`` holds only asymptotically. See the
   caveat block in
   :func:`~orpheus.derivations.peierls_geometry.build_white_bc_correction`
   for the error table. Tests here gate on the **thick-limit**
   behaviour, where the closure is quantitatively accurate.

Structural checks validated here:

1. ``P_esc`` behaves correctly: decreasing with :math:`r_i \\to R`
   (neutrons born near the boundary escape with high probability)
   and bounded by :math:`[0, 1]`.
2. ``G_bc`` is positive and non-zero for every interior node.
3. White-BC :math:`k_{\\rm eff}` exceeds vacuum-BC :math:`k_{\\rm eff}`
   at every finite :math:`R` (less leakage ⇒ higher eigenvalue).
4. White-BC :math:`k_{\\rm eff}` converges to
   :math:`k_\\infty = \\nu\\Sigma_f / \\Sigma_a` in the thick-cylinder
   limit, at a rate consistent with the rank-1 closure error.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations import peierls_geometry as _pg
from orpheus.derivations.peierls_cylinder import GEOMETRY
from orpheus.derivations.peierls_geometry import (
    compute_G_bc,
    compute_P_esc,
)


_SIG_T = np.array([1.0])
_SIG_S = np.array([0.5])
_NU_SIG_F = np.array([0.75])
_K_INF = _NU_SIG_F[0] / (_SIG_T[0] - _SIG_S[0])  # = 1.5


@pytest.mark.l0
@pytest.mark.verifies("peierls-equation")
class TestPescProperties:
    """``P_esc(r_i)`` — uncollided escape probability from interior."""

    def test_bounded_in_unit_interval(self):
        r_nodes = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
        P_esc = compute_P_esc(
            GEOMETRY, r_nodes, np.array([1.0]), np.array([1.0]),
            n_angular=16, dps=20,
        )
        assert np.all(P_esc >= 0.0)
        assert np.all(P_esc <= 1.0)

    def test_increases_toward_boundary(self):
        """For uniform Σ_t, points closer to the surface escape more
        easily (shorter mean chord to surface)."""
        r_nodes = np.linspace(0.05, 0.95, 7)
        P_esc = compute_P_esc(
            GEOMETRY, r_nodes, np.array([1.0]), np.array([1.0]),
            n_angular=16, dps=20,
        )
        diffs = np.diff(P_esc)
        assert np.all(diffs > -1e-6), (
            f"P_esc not monotone in r: diffs = {diffs}"
        )


@pytest.mark.l0
@pytest.mark.verifies("peierls-equation")
class TestGbcProperties:
    """``G_bc(r_i)`` — flux at :math:`r_i` from unit uniform surface current."""

    def test_positive_everywhere(self):
        r_nodes = np.array([0.1, 0.5, 0.95])
        G_bc = compute_G_bc(
            GEOMETRY, r_nodes, np.array([1.0]), np.array([1.0]),
            n_surf_quad=16, dps=20,
        )
        assert np.all(G_bc > 0.0)

    def test_increases_toward_boundary(self):
        """Uniform surface source contributes more to points nearer
        the boundary (unattenuated path length is shorter)."""
        r_nodes = np.linspace(0.05, 0.95, 7)
        G_bc = compute_G_bc(
            GEOMETRY, r_nodes, np.array([1.0]), np.array([1.0]),
            n_surf_quad=16, dps=20,
        )
        diffs = np.diff(G_bc)
        assert np.all(diffs > -1e-6), (
            f"G_bc not monotone in r: diffs = {diffs}"
        )


@pytest.mark.l1
@pytest.mark.verifies("peierls-equation", "one-group-kinf")
class TestWhiteBCEigenvalue:
    """White-BC eigenvalue approaches k_inf in the thick limit,
    and exceeds vacuum-BC everywhere."""

    def test_k_eff_thick_approaches_kinf(self):
        """At R = 10 MFP, the rank-1 white-BC closure gives k_eff
        within 2 % of k_inf = 1.5."""
        sol = _pg.solve_peierls_1g(
            _pg.CYLINDER_1D,
            radii=np.array([10.0]),
            sig_t=_SIG_T, sig_s=_SIG_S, nu_sig_f=_NU_SIG_F,
            boundary="white",
            n_panels_per_region=3, p_order=5,
            n_angular=20, n_rho=20, n_surf_quad=20, dps=20,
        )
        err = abs(sol.k_eff - _K_INF) / _K_INF
        assert err < 2e-2, (
            f"Thick-limit k_eff(white, R=10) = {sol.k_eff:.6f}, "
            f"err vs k_inf = {err:.3e} (>2%)"
        )

    def test_white_exceeds_vacuum(self):
        """White BC suppresses leakage ⇒ k_eff(white) > k_eff(vacuum)
        at every finite R."""
        for R in (2.0, 5.0, 10.0):
            sol_vac = _pg.solve_peierls_1g(
                _pg.CYLINDER_1D,
                radii=np.array([R]),
                sig_t=_SIG_T, sig_s=_SIG_S, nu_sig_f=_NU_SIG_F,
                boundary="vacuum",
                n_panels_per_region=3, p_order=5,
                n_angular=18, n_rho=18, dps=20,
            )
            sol_wht = _pg.solve_peierls_1g(
                _pg.CYLINDER_1D,
                radii=np.array([R]),
                sig_t=_SIG_T, sig_s=_SIG_S, nu_sig_f=_NU_SIG_F,
                boundary="white",
                n_panels_per_region=3, p_order=5,
                n_angular=18, n_rho=18, n_surf_quad=18, dps=20,
            )
            assert sol_wht.k_eff > sol_vac.k_eff, (
                f"At R={R}: white k_eff = {sol_wht.k_eff:.6f} should "
                f"exceed vacuum k_eff = {sol_vac.k_eff:.6f}"
            )

    @pytest.mark.slow
    def test_k_eff_monotone_approach_to_kinf(self):
        """k_eff(white) → k_inf monotonically from below as R increases."""
        keffs = []
        for R in (3.0, 5.0, 10.0, 20.0):
            sol = _pg.solve_peierls_1g(
                _pg.CYLINDER_1D,
                radii=np.array([R]),
                sig_t=_SIG_T, sig_s=_SIG_S, nu_sig_f=_NU_SIG_F,
                boundary="white",
                n_panels_per_region=3, p_order=5,
                n_angular=16, n_rho=16, n_surf_quad=16, dps=20,
            )
            keffs.append(sol.k_eff)
        # All below k_inf
        assert all(k < _K_INF for k in keffs), keffs
        # Monotone increasing
        assert all(keffs[i+1] >= keffs[i] - 1e-5 for i in range(len(keffs)-1)), keffs
        # Last point within 2 % of k_inf
        assert abs(keffs[-1] - _K_INF) / _K_INF < 2e-2


# ──────────────────────────────────────────────────────────────────────
# NEGATIVE-regression gate: scalar Hébert (1-P_ss)⁻¹ alone is INSUFFICIENT
# on cylinder Class B (Issue #132 WONTFIX)
# ──────────────────────────────────────────────────────────────────────
# Issue #132 investigated whether the Hébert (1 - P_ss^cyl)⁻¹ factor —
# which closes the sphere Class B Mark gap to <1.5 % on 1G/1R, 2G/1R,
# 2G/2R — could be applied to cylinder. Result: it does NOT close the
# cylinder Mark gap (the bare Mark closure has a -22 % to -82 % error
# on Class B 1G-2G/1R-2R cases; the patched closure inherits the same
# error because the cylinder Mark calibration limit is upstream of the
# geometric-series factor — most likely the ``compute_G_bc`` cylinder
# branch's 2-D-projected-cosine limitation, which needs the 3-D
# Knyazev correction (Issue #112 Phase C / Hypothesis 3).
#
# This is a NEGATIVE-regression gate: it pins that the Mark-only path
# on cylinder Class B has a documented error in the WONTFIX range
# (-10 % to -90 %). Any future closure improvement that closes the gap
# below 10 % is a structural fix and must update this gate.
#
# The pss diagnostic that supplies the multi-region P_ss^cyl
# integration lives in
# ``derivations/diagnostics/diag_cylinder_hebert_pss.py`` (an
# investigation primitive, not a diag scratchpad — kept).
#
# Promoted from
# ``derivations/diagnostics/diag_cylinder_hebert_keff.py`` (2026-04-30
# triage).


@pytest.mark.foundation
class TestHebertCylinderInsufficient:
    """Cylinder Mark closure on Class B has a documented -10 % to -90 %
    error that the Hébert (1 - P_ss)⁻¹ factor alone does NOT close.

    The test is a NEGATIVE-regression gate: it asserts the gap on each
    Class B fixture is in the published WONTFIX range. If a future
    closure improvement reduces the gap below the lower bound, this
    test SHOULD start failing — that is the signal the structural fix
    has been integrated and the gate must be updated.
    """

    @pytest.mark.parametrize(
        "ng_key,n_regions,err_min_pct,err_max_pct",
        [
            # Empirical numbers from 2026-04-30 triage of the keff diag:
            # 1G/1R: -21.85 %, 1G/2R: -22.52 %, 2G/1R: -82.06 %, 2G/2R: -76.65 %.
            # Tolerance: ±5 % around each pinned value to absorb quadrature
            # / numerical jitter while still failing on a ≥10 % structural
            # change in the error magnitude.
            ("1g", 1, 16.0, 27.0),
            ("1g", 2, 17.0, 28.0),
            ("2g", 1, 75.0, 88.0),
            ("2g", 2, 70.0, 83.0),
        ],
    )
    def test_mark_only_cylinder_class_B_gap_in_wontfix_range(
        self, ng_key, n_regions, err_min_pct, err_max_pct,
    ):
        """rank-1 Mark cylinder Class B gap stays in the documented
        WONTFIX range (~-10 % to -90 %).

        Bug it would catch: a regression that *worsens* the Mark
        cylinder error (>err_max_pct%) — likely a sign error or
        quadrature drift in ``compute_G_bc`` cylinder branch — OR a
        silent structural fix that *closes* the gap below err_min_pct%
        without updating this gate.
        """
        from orpheus.derivations import cp_cylinder
        from orpheus.derivations._xs_library import LAYOUTS, get_xs
        from orpheus.derivations.peierls_geometry import (
            CYLINDER_1D, solve_peierls_mg,
        )

        layout = LAYOUTS[n_regions]
        xs_list = [get_xs(region, ng_key) for region in layout]
        sig_t = np.vstack([xs["sig_t"] for xs in xs_list])
        sig_s = np.stack([xs["sig_s"] for xs in xs_list], axis=0)
        nu_sig_f = np.vstack([xs["nu"] * xs["sig_f"] for xs in xs_list])
        chi = np.vstack([xs["chi"] for xs in xs_list])
        radii = np.array(cp_cylinder._RADII[n_regions])
        kinf = cp_cylinder._build_case(ng_key, n_regions).k_inf

        sol = solve_peierls_mg(
            CYLINDER_1D, radii=radii, sig_t=sig_t, sig_s=sig_s,
            nu_sig_f=nu_sig_f, chi=chi,
            boundary="white_rank1_mark", n_bc_modes=1,
            n_panels_per_region=2, p_order=3,
            n_angular=24, n_rho=24, n_surf_quad=24, dps=15,
        )
        err_pct = abs(sol.k_eff - kinf) / kinf * 100.0
        assert err_min_pct <= err_pct <= err_max_pct, (
            f"Cylinder Mark Class-B {ng_key}/{n_regions}r error "
            f"{err_pct:.3f} % outside WONTFIX range "
            f"[{err_min_pct:.1f}, {err_max_pct:.1f}] %. "
            f"k_eff={sol.k_eff:.6f}, k_inf={kinf:.6f}. "
            f"If err_pct < {err_min_pct:.1f} %, a structural closure "
            f"improvement has landed — update this gate (Issue #132 / "
            f"#112 Phase C)."
        )
