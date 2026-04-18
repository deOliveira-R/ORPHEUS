"""L1 white-BC rank-1 closure tests for the spherical Peierls solver.

The rank-1 Schur closure treats the outgoing/incoming partial
currents as uniform-isotropic (Mark's closure), which is exact only
in the flat-source (region-averaged) limit. Applied at the pointwise
Nyström level, it produces an error in :math:`k_{\\rm eff}` that
shrinks with cell optical thickness.

This test file pins the observed sphere behaviour so a regression
surfaces as a test failure rather than a silent degradation:

=====  ===============  =============
R/MFP  k_eff(white-sph) err vs k_inf
=====  ===============  =============
1.0    ~1.10            ~27 %
2.0    ~1.39            ~7 %
5.0    ~1.49            ~0.7 %
10.0   ~1.50            ~0.3 %
=====  ===============  =============

These are parallel to the cylinder values (1.19 / 1.40 / 1.48 /
1.49) documented in the Phase-4.2 campaign — rank-1 Mark closure
suffers the same inverse-cell-size divergence in both geometries,
and Issue #103 (N1) is the architectural fix.

The rigorous floor here is *monotonic improvement with R*; the
tabulated percentages are pinned loosely (generous tolerances) so
the test does not fail under incidental quadrature drift."""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.peierls_sphere import solve_peierls_sphere_1g


# ═══════════════════════════════════════════════════════════════════════
# Rank-1 closure error scan
# ═══════════════════════════════════════════════════════════════════════

@pytest.mark.l1
@pytest.mark.verifies("peierls-unified")
class TestWhiteBCRank1ErrorScan:
    """k_eff(white, R) monotonically approaches k_inf as R grows."""

    @staticmethod
    def _k_inf():
        # Homogeneous, pure scatter + small fission
        sig_t, sig_s, nu_sig_f = 1.0, 0.5, 0.75
        return sig_t, sig_s, nu_sig_f, nu_sig_f / (sig_t - sig_s)

    @classmethod
    def _solve(cls, R: float):
        sig_t, sig_s, nu_sig_f, _ = cls._k_inf()
        return solve_peierls_sphere_1g(
            np.array([R]),
            np.array([sig_t]), np.array([sig_s]), np.array([nu_sig_f]),
            boundary="white",
            n_panels_per_region=2, p_order=5,
            n_theta=20, n_rho=20, n_phi=32, dps=25,
        )

    def test_k_eff_monotone_in_R(self):
        """|k_eff(R) − k_inf| decreases from R = 2 MFP to R = 10 MFP
        (the regime where the rank-1 closure error dominates).

        Beyond R ≈ 10 MFP the rank-1 error is already below ~1 % and
        quadrature / iteration noise take over; the **bounded-bias**
        part of the test only checks that the quadrature-noise regime
        sits near k_inf (loose tolerance)."""
        _, _, _, k_inf = self._k_inf()
        R_values = [2.0, 5.0, 10.0]
        errors = [abs(self._solve(R).k_eff - k_inf) for R in R_values]
        diffs = np.diff(errors)
        assert np.all(diffs <= 1e-3), (
            f"Error in k_eff(white, R) not monotone-decreasing from "
            f"R=2 MFP to R=10 MFP: R={R_values}, errors={errors}"
        )
        # Asymptote floor: beyond R=10 MFP the error is already small;
        # refining further keeps it under 1 %.
        err_20 = abs(self._solve(20.0).k_eff - k_inf)
        assert err_20 < 1e-2, (
            f"R=20 MFP rank-1 error = {err_20*100:.2f} % exceeds 1 %; "
            f"quadrature regime expected"
        )

    def test_thin_sphere_rank1_error_bounded(self):
        """R = 1 MFP: rank-1 error is bounded above at ~30 %
        (matches the cylinder's ~21 % at the same size within the
        geometry-dependent factor)."""
        _, _, _, k_inf = self._k_inf()
        sol = self._solve(1.0)
        err = abs(sol.k_eff - k_inf) / k_inf
        assert err < 0.35, (
            f"R=1 MFP rank-1 error = {err*100:.1f} % exceeds 35 %; "
            f"expected ~27 % (k_eff ≈ 1.10 vs k_inf = 1.5)"
        )
        # And significantly away from the regression-test floor of 1 %
        # — a passing test at R=1 would indicate the rank-1 error has
        # mysteriously vanished, which would be suspicious.
        assert err > 0.05, (
            f"R=1 MFP rank-1 error suspiciously small: {err*100:.1f} % "
            f"(expected ~27 % from Phase-4.3 pinning)"
        )

    def test_medium_sphere_rank1_error_bounded(self):
        """R = 5 MFP: rank-1 error below ~3 % (pinning the ~0.7 %
        observation with generous headroom)."""
        _, _, _, k_inf = self._k_inf()
        sol = self._solve(5.0)
        err = abs(sol.k_eff - k_inf) / k_inf
        assert err < 3e-2, (
            f"R=5 MFP rank-1 error = {err*100:.2f} % exceeds 3 %; "
            f"expected ~0.7 %"
        )

    def test_thick_sphere_rank1_near_k_inf(self):
        """R = 10 MFP: rank-1 error below 1.5 %."""
        _, _, _, k_inf = self._k_inf()
        sol = self._solve(10.0)
        err = abs(sol.k_eff - k_inf) / k_inf
        assert err < 1.5e-2, (
            f"R=10 MFP rank-1 error = {err*100:.2f} % exceeds 1.5 %; "
            f"expected ~0.3 %"
        )
