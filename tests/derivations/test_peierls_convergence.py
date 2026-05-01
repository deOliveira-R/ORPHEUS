"""L0 self-convergence of the Peierls Nyström solver.

Before trusting the Peierls reference for verifying CP, we must verify
that the Nyström discretisation itself converges under quadrature
refinement. These tests check:

1. The Nyström eigenvalue converges to a fixed value as panels increase.
2. The Nyström eigenvalue agrees with the existing CP eigenvalue
   (which uses a completely different method — E₃ flat-source matrices).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations import reference_values
from orpheus.derivations.common.xs_library import LAYOUTS, get_xs
from orpheus.derivations.continuous.peierls.slab import solve_peierls_eigenvalue


@pytest.mark.l0
@pytest.mark.verifies("peierls-equation")
class TestPeierlsSelfConvergence:
    """The Peierls Nyström solver converges under quadrature refinement."""

    def test_nystrom_eigenvalue_converges_1g(self):
        """1G k_eff from successive refinements converges monotonically."""
        xs_list = [get_xs("A", "1g")]
        keffs = []
        for n_panels in [4, 8, 16]:
            sol = solve_peierls_eigenvalue(
                [xs["sig_t"] for xs in xs_list],
                [xs["sig_s"] for xs in xs_list],
                [xs["nu"] * xs["sig_f"] for xs in xs_list],
                [xs["chi"] for xs in xs_list],
                thicknesses=[0.5],
                n_panels_per_region=n_panels,
                p_order=4,
                precision_digits=20,
                boundary="white",
            )
            keffs.append(sol.k_eff)

        d1 = abs(keffs[1] - keffs[0])
        d2 = abs(keffs[2] - keffs[1])
        assert d2 < d1, (
            f"Not converging: d1={d1:.3e}, d2={d2:.3e}, "
            f"keffs={[f'{k:.8f}' for k in keffs]}"
        )

    @pytest.mark.parametrize("case_name", [
        "cp_slab_1eg_1rg",
        "cp_slab_1eg_2rg",
        "cp_slab_2eg_1rg",
    ])
    def test_nystrom_matches_cp_eigenvalue(self, case_name: str):
        """Peierls k_eff agrees with the CP eigenvalue to < 1%.

        The Peierls solver uses E₁ via Nyström quadrature; the CP
        solver uses E₃ flat-source matrices. Agreement confirms both
        methods solve the same integral transport equation.

        The 2G 2-region case is excluded from this parametric test
        because it is too slow for routine CI; it is covered by the
        dedicated ``test_nystrom_matches_cp_2g_2r`` below.
        """
        case = reference_values.get(case_name)
        ng_key = f"{case.n_groups}g"
        layout = LAYOUTS[case.n_regions]
        xs_list = [get_xs(r, ng_key) for r in layout]
        thicknesses = case.geom_params["thicknesses"]

        sol = solve_peierls_eigenvalue(
            [xs["sig_t"] for xs in xs_list],
            [xs["sig_s"] for xs in xs_list],
            [xs["nu"] * xs["sig_f"] for xs in xs_list],
            [xs["chi"] for xs in xs_list],
            thicknesses=thicknesses,
            n_panels_per_region=8,
            p_order=6,
            precision_digits=20,
            boundary="white",
        )

        rel_err = abs(sol.k_eff - case.k_inf) / case.k_inf
        assert rel_err < 0.01, (
            f"Peierls k={sol.k_eff:.8f} vs CP k={case.k_inf:.8f}, "
            f"relative error {rel_err:.3e} > 1%"
        )

    @pytest.mark.slow
    def test_nystrom_matches_cp_2g_2r(self):
        """2G 2-region Peierls k_eff agrees with CP to < 2%.

        This test uses a coarser quadrature than the 1G tests because
        the 2G solve is significantly slower (O(minutes) with mpmath).
        """
        case = reference_values.get("cp_slab_2eg_2rg")
        xs_list = [get_xs(r, "2g") for r in LAYOUTS[2]]

        sol = solve_peierls_eigenvalue(
            [xs["sig_t"] for xs in xs_list],
            [xs["sig_s"] for xs in xs_list],
            [xs["nu"] * xs["sig_f"] for xs in xs_list],
            [xs["chi"] for xs in xs_list],
            thicknesses=[0.5, 0.5],
            n_panels_per_region=4,
            p_order=4,
            precision_digits=20,
            boundary="white",
        )

        rel_err = abs(sol.k_eff - case.k_inf) / case.k_inf
        assert rel_err < 0.02, (
            f"Peierls k={sol.k_eff:.8f} vs CP k={case.k_inf:.8f}, "
            f"relative error {rel_err:.3e} > 2%"
        )
