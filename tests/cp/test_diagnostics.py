"""Verify CP solver diagnostics: residual tracking, Gauss-Seidel mode, inner iterations."""

import numpy as np
import pytest

from orpheus.derivations import get
from orpheus.geometry import CoordSystem, Mesh1D
from orpheus.cp.solver import solve_cp, CPParams

# Diagnostic *term* tests (residuals, inner-iter counts) are L0.
# NOT file-level — the parametrized eigenvalue tests below inherit L1
# from their VerificationCase via the conftest hook, and a file-level
# l0 mark would incorrectly override that inheritance.


# ═══════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════

def _make_mesh(case):
    """Build a Mesh1D from a VerificationCase."""
    gp = case.geom_params
    if case.geometry == "slab":
        thicknesses = np.array(gp["thicknesses"])
        edges = np.concatenate([[0.0], np.cumsum(thicknesses)])
        coord = CoordSystem.CARTESIAN
    elif case.geometry == "cyl1D":
        radii = np.array(gp["radii"])
        edges = np.concatenate([[0.0], radii])
        coord = CoordSystem.CYLINDRICAL
    elif case.geometry == "sph1D":
        radii = np.array(gp["radii"])
        edges = np.concatenate([[0.0], radii])
        coord = CoordSystem.SPHERICAL
    else:
        raise ValueError(f"Unknown geometry: {case.geometry}")
    return Mesh1D(edges=edges, mat_ids=np.array(gp["mat_ids"]), coord=coord)


def _tolerance_for(case):
    """Slab gets 1e-6 tolerance, radial geometries get 1e-5."""
    return 1e-6 if case.geometry == "slab" else 1e-5


# ═══════════════════════════════════════════════════════════════════════
# Phase 2: Jacobi diagnostics
# ═══════════════════════════════════════════════════════════════════════

@pytest.mark.l0
class TestJacobiDiagnostics:
    """Jacobi: residual_history is populated and there is NO inner level."""

    @pytest.fixture()
    def result(self):
        case = get("cp_slab_2eg_2rg")
        mesh = _make_mesh(case)
        params = CPParams(keff_tol=1e-7, flux_tol=1e-6, solver_mode="jacobi")
        return solve_cp(case.materials, mesh, params), case

    def test_residual_history_populated(self, result):
        res, _ = result
        assert len(res.residual_history) > 0

    def test_final_residual_small(self, result):
        res, _ = result
        assert res.residual_history[-1] < 1e-3

    def test_jacobi_records_no_inner_level(self, result):
        """#340 N4: successor to the retired ``n_inner is None`` assertion.

        The absence is DELIBERATE, and the record states it structurally
        rather than by a sentinel: Jacobi lags the scattering source, so no
        inner level exists to converge and an empty ``children`` is the
        honest tree.
        """
        res, _ = result
        assert res.record.children == ()
        assert res.record.fully_converged is res.record.converged

    def test_keff_matches_analytical(self, result):
        res, case = result
        err = abs(res.keff - case.k_inf)
        assert err < 1e-6, (
            f"Jacobi keff={res.keff:.10f} analytical={case.k_inf:.10f} err={err:.2e}"
        )


# ═══════════════════════════════════════════════════════════════════════
# Phase 3a: Gauss-Seidel single-case diagnostics
# ═══════════════════════════════════════════════════════════════════════

@pytest.mark.l0
class TestGaussSeidelDiagnostics:
    """Verify GS mode produces correct keff and diagnostics."""

    @pytest.fixture()
    def result(self):
        case = get("cp_slab_2eg_2rg")
        mesh = _make_mesh(case)
        params = CPParams(
            keff_tol=1e-7, flux_tol=1e-6, solver_mode="gauss_seidel",
        )
        return solve_cp(case.materials, mesh, params), case

    def test_keff_matches_analytical(self, result):
        res, case = result
        err = abs(res.keff - case.k_inf)
        assert err < 1e-6, (
            f"GS keff={res.keff:.10f} analytical={case.k_inf:.10f} err={err:.2e}"
        )

    def test_gs_records_one_inner_per_outer_per_group(self, result):
        """#340 N4: successor to the retired ``n_inner.shape`` assertion."""
        res, _ = result
        record = res.record
        ng = 2  # 2 energy groups
        assert record.n_iterations > 0
        assert len(record.children) == ng * record.n_iterations
        assert {
            child.label.rsplit("g=", 1)[1].rstrip(")")
            for child in record.children
        } == {str(g) for g in range(ng)}

    def test_every_inner_ran_AND_measured_something(self, result):
        """Counts positive, AND a residual behind each one.

        The retired array could only say ``>= 1``.  A count-only child cannot
        claim convergence at all (the primitive refuses it), so asserting the
        trajectory exists is asserting the thing that makes a count mean
        anything.
        """
        res, _ = result
        for child in res.record.children:
            assert child.n_iterations >= 1
            assert child.criteria[0].trajectory

    def test_residual_history_populated(self, result):
        res, _ = result
        assert len(res.residual_history) > 0
        assert res.residual_history[-1] < 1e-3


# ═══════════════════════════════════════════════════════════════════════
# Phase 3b: Gauss-Seidel — all 27 CP cases
# ═══════════════════════════════════════════════════════════════════════

_ALL_CP_SLAB = [f"cp_slab_{g}eg_{r}rg" for g in (1, 2, 4) for r in (1, 2, 4)]
_ALL_CP_CYL = [f"cp_cyl1D_{g}eg_{r}rg" for g in (1, 2, 4) for r in (1, 2, 4)]
_ALL_CP_SPH = [f"cp_sph1D_{g}eg_{r}rg" for g in (1, 2, 4) for r in (1, 2, 4)]


@pytest.mark.parametrize("case_name", _ALL_CP_SLAB + _ALL_CP_CYL + _ALL_CP_SPH)
def test_gauss_seidel_eigenvalue(case_name):
    """GS mode must produce the same eigenvalue as Jacobi for all 27 CP cases."""
    case = get(case_name)
    mesh = _make_mesh(case)
    tol = _tolerance_for(case)
    params = CPParams(
        keff_tol=1e-7, flux_tol=1e-6, solver_mode="gauss_seidel",
    )
    result = solve_cp(case.materials, mesh, params)
    err = abs(result.keff - case.k_inf)
    assert err < tol, (
        f"{case_name} GS: solver={result.keff:.10f} "
        f"analytical={case.k_inf:.10f} err={err:.2e}"
    )


# ═══════════════════════════════════════════════════════════════════════
# Phase 4: Thermal groups need more inner iterations than fast groups
# ═══════════════════════════════════════════════════════════════════════

@pytest.mark.l1
def test_thermal_group_needs_more_inner_iterations():
    """Thermal groups (high index) should need >= inner iterations than fast (low index).

    The 4-group problem has self-scatter in thermal groups that requires
    multiple inner iterations.  The fast group (g=0) should converge
    in fewer (or equal) inner iterations than the thermal group (g=3).
    """
    case = get("cp_slab_4eg_4rg")
    mesh = _make_mesh(case)
    params = CPParams(
        keff_tol=1e-7, flux_tol=1e-6, solver_mode="gauss_seidel",
    )
    result = solve_cp(case.materials, mesh, params)

    # #340 N4: group by each child's own LABEL — the retired ``n_inner``
    # array needed an agreed axis order; the record carries the group itself.
    def mean_count(g: int) -> float:
        counts = [
            child.n_iterations
            for child in result.record.children
            if child.label.endswith(f"g={g})")
        ]
        assert counts, f"no inner records for group {g}"
        return sum(counts) / len(counts)

    fast, thermal = mean_count(0), mean_count(3)

    # Thermal group (g=3) should need >= inner iters than fast group (g=0)
    assert thermal >= fast, (
        f"Expected thermal group to need more inner iterations: "
        f"fast(g=0)={fast:.1f} thermal(g=3)={thermal:.1f}"
    )
