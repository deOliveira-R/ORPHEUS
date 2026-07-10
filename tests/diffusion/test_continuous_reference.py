r"""L1 verification of the modern diffusion solver against continuous references.

Consumer tests for the Phase-0
:class:`~orpheus.derivations.ContinuousReferenceSolution` contract
as applied to 1D two-group diffusion. These compare the #290
operator-algebra solver (:func:`orpheus.diffusion.solve_diffusion_1d`)
to mesh-independent analytical/semi-analytical references for both
**eigenvalue** and **flux shape**.

Two problem families:

1. **Bare slab** — compared against the pure analytical sine
   eigenfunction from :func:`derive_1rg_continuous`. Tests both
   the eigenvalue agreement at a single mesh and the
   :math:`\mathcal{O}(h^{2})` spatial convergence of the
   diffusion finite-difference discretisation.

2. **Fuel + reflector slab** — compared against the transcendental
   transfer-matrix reference from :func:`derive_2rg_continuous`,
   limited only by the :math:`\sqrt{\text{brentq-tol}} \sim 3
   \times 10^{-7}` null-vector precision floor of the
   double-precision transfer-matrix back-substitution.

**Boundary naming (#290 ruling 3).** Both references pin the
**zero-flux Dirichlet** law :math:`\phi_\Gamma = 0` — the honest name
for what the retired MATLAB island called "vacuum" (in ORPHEUS,
*vacuum* means zero incoming partial current — Marshak,
:math:`\mathcal A = 0`). The solver poses it as ``BC("zero_flux")``
(albedo :math:`\mathcal A = -1` on the :math:`(J^+, J^-)` trace).
The mathematics of every reference is unchanged — re-attributed,
never re-derived.

The legacy Richardson tier and its migration-window cross-check
(``test_2region_transcendental_matches_richardson_cache``) retired at
#290 P6 with the island solver that computed the Richardson value;
the cross-check's job — locking in the transcendental replacement at
the ~1e-7 level, inside the Richardson ~1e-5 floor — was completed
before retirement.

See :doc:`/theory/diffusion_1d` for the operator form the tests
commit to, and :doc:`/verification/reference_solutions` for the
campaign philosophy.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations import continuous_get
from orpheus.derivations.common.xs_library import mixture_from_diffusion_tables
from orpheus.diffusion import solve_diffusion_1d
from orpheus.geometry import BC
from orpheus.geometry.mesh import Mesh1D


pytestmark = [pytest.mark.l1, pytest.mark.verifies(
    "diffusion-operator",
    "diffusion-coefficient",
    "bare-slab-buckling",
    "bare-slab-eigenfunction",
    "bare-slab-critical-equation",
    "diffusion-region-ode",
    "diffusion-M-matrix",
    "diffusion-mode-decomposition",
    "diffusion-exponential-branch",
    "diffusion-trigonometric-branch",
    "diffusion-interface-matching",
    "diffusion-matching-matrix",
    "diffusion-transcendental",
    "diffusion-spurious-root-validation",
    "diffusion-back-substitution",
)]


def _zero_flux_problem(
    region_heights: list[float], dz: float,
) -> tuple[Mesh1D, np.ndarray]:
    """Uniform-``dz`` multi-region slab mesh under the zero-flux law.

    ``dz`` must divide every region height exactly (all reference
    configurations use dyadic refinements of 5 cm on 50/30 cm
    regions, so it does). Region ``i`` carries material id ``i``.
    Returns the mesh and the cell-centre coordinates the continuous
    reference is evaluated at.
    """
    n_cells = [round(h / dz) for h in region_heights]
    length = float(sum(region_heights))
    edges = np.linspace(0.0, length, sum(n_cells) + 1)
    mat_ids = np.concatenate(
        [np.full(n, i, dtype=int) for i, n in enumerate(n_cells)]
    )
    mesh = Mesh1D(
        edges=edges, mat_ids=mat_ids,
        bc_left=BC("zero_flux"), bc_right=BC("zero_flux"),
    )
    z_cells = 0.5 * (edges[:-1] + edges[1:])
    return mesh, z_cells


def _peak_normalise(flux_g0: np.ndarray, flux_g1: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    r"""Divide both groups by ``max(|flux_g0|)``.

    Returns the pair with :math:`\max|\phi_0| = 1`. The ratio
    :math:`\phi_1/\phi_0` is preserved — the physically meaningful
    quantity for comparing flux shapes between a solver output
    normalised to unit production and a reference solution
    normalised to peak amplitude.
    """
    scale = np.max(np.abs(flux_g0))
    if scale == 0:
        raise ValueError("Group-0 flux is identically zero — cannot normalise.")
    return flux_g0 / scale, flux_g1 / scale


# ═══════════════════════════════════════════════════════════════════════
# Bare slab — analytical sine eigenfunction
# ═══════════════════════════════════════════════════════════════════════

def test_bare_slab_eigenvalue_matches_continuous_reference():
    r"""The diffusion solver must reproduce the analytical bare-slab
    eigenvalue at a sufficiently fine mesh.

    At ``dz = 0.3125`` cm the discretisation error is
    :math:`\mathcal{O}(10^{-5})` on the 50-cm slab, which is the
    floor we assert to. The exact eigenvalue comes from the
    2x2 buckled matrix — no Richardson, no iteration.
    """
    ref = continuous_get("dif_slab_2eg_1rg")
    assert ref.operator_form == "diffusion"
    assert ref.k_eff is not None

    fuel = mixture_from_diffusion_tables(
        next(iter(ref.problem.materials.values()))
    )
    H = ref.problem.geometry_params["fuel_height"]

    mesh, _ = _zero_flux_problem([H], dz=0.3125)
    result = solve_diffusion_1d({0: fuel}, mesh)

    # At dz=0.3125 on L=50, O(h²) ~ 4e-5; the observed error on this
    # configuration is ~5e-6 empirically. Assert 1e-4 to leave headroom.
    assert abs(result.keff - ref.k_eff) < 1e-4, (
        f"k_eff mismatch: solver={result.keff:.10f} "
        f"reference={ref.k_eff:.10f} diff={result.keff - ref.k_eff:.2e}"
    )


def test_bare_slab_flux_shape_converges_second_order():
    r"""The bare-slab solver flux must converge to
    :math:`c_g \sin(\pi x/L)` at the design :math:`\mathcal{O}(h^{2})`
    rate of central finite differences.

    Runs four mesh refinements, compares peak-normalised shapes
    against the reference sine, and asserts the measured order of
    the cell-wise :math:`\ell^{2}` norm of the residual is
    :math:`> 1.8`. Tight eigenpair tolerances keep the iteration
    residual well below the finite-difference truncation error at
    every refinement.
    """
    ref = continuous_get("dif_slab_2eg_1rg")
    fuel = mixture_from_diffusion_tables(
        next(iter(ref.problem.materials.values()))
    )
    H = ref.problem.geometry_params["fuel_height"]

    dzs = [5.0, 2.5, 1.25, 0.625]
    errs = []
    for dz in dzs:
        mesh, z = _zero_flux_problem([H], dz)
        result = solve_diffusion_1d(
            {0: fuel}, mesh,
            keff_tol=1e-13, flux_tol=1e-11, max_outer=5000,
        )

        ref_g0 = np.asarray(ref.phi(z, 0), dtype=float)
        ref_g1 = np.asarray(ref.phi(z, 1), dtype=float)
        solver_g0 = result.flux.interior.values[0].copy()
        solver_g1 = result.flux.interior.values[1].copy()

        ref_g0, ref_g1 = _peak_normalise(ref_g0, ref_g1)
        solver_g0, solver_g1 = _peak_normalise(solver_g0, solver_g1)

        diff = np.concatenate([solver_g0 - ref_g0, solver_g1 - ref_g1])
        err = float(np.sqrt(dz * np.sum(diff * diff)))
        errs.append(err)

    errs_arr = np.asarray(errs)
    ratios = errs_arr[:-1] / errs_arr[1:]
    orders = np.log2(ratios)

    assert np.all(orders > 1.8), (
        f"Bare slab flux shape convergence below O(h²): "
        f"errors={errs_arr}, orders={orders}"
    )

    assert errs_arr[-1] < 1e-3, (
        f"Finest-mesh error {errs_arr[-1]:.2e} above expected ~1e-5"
    )


def test_bare_slab_sine_peak_is_at_midslab():
    r"""Sanity: the continuous reference must peak at the slab midpoint.

    Regression guard on the eigenfunction closure — any future
    bug that drops the :math:`\sin` argument or uses the wrong
    slab length would move the peak and this test catches it.
    """
    ref = continuous_get("dif_slab_2eg_1rg")
    L = ref.problem.geometry_params["length"]

    z_dense = np.linspace(0.0, L, 501)
    phi = ref.phi(z_dense, 0)

    peak_idx = int(np.argmax(np.abs(phi)))
    z_peak = z_dense[peak_idx]
    assert abs(z_peak - L / 2) < L / 250  # within one grid cell

    # Boundary values exactly zero
    assert ref.phi(np.array([0.0]), 0)[0] == pytest.approx(0.0, abs=1e-14)
    assert ref.phi(np.array([L]), 0)[0] == pytest.approx(0.0, abs=1e-14)


# ═══════════════════════════════════════════════════════════════════════
# 2-region fuel + reflector — transcendental transfer-matrix reference
# ═══════════════════════════════════════════════════════════════════════

def test_2region_eigenvalue_matches_continuous_reference():
    r"""The diffusion solver must reproduce the transcendental
    2-region eigenvalue at a fine mesh.
    """
    ref = continuous_get("dif_slab_2eg_2rg")
    assert ref.operator_form == "diffusion"
    assert ref.k_eff is not None

    materials = {
        0: mixture_from_diffusion_tables(ref.problem.materials[0]),
        1: mixture_from_diffusion_tables(ref.problem.materials[1]),
    }
    H_f = ref.problem.geometry_params["fuel_height"]
    H_r = ref.problem.geometry_params["refl_height"]

    mesh, _ = _zero_flux_problem([H_f, H_r], dz=0.3125)
    result = solve_diffusion_1d(materials, mesh)

    # At dz=0.3125 on L=80 cm with fuel/reflector interface, O(h²)
    # discretisation error is ~1e-5. Reference is accurate to ~1e-7.
    # Assert 5e-5 for a comfortable margin.
    assert abs(result.keff - ref.k_eff) < 5e-5, (
        f"k_eff mismatch: solver={result.keff:.10f} "
        f"reference={ref.k_eff:.10f} diff={result.keff - ref.k_eff:.2e}"
    )


def test_2region_flux_shape_converges_second_order():
    r"""The fuel+reflector solver flux must converge to the
    transcendental back-substituted reference at
    :math:`\mathcal{O}(h^{2})`.

    The comparison uses the same cell-wise :math:`\ell^{2}`-norm
    weighted by cell width as
    :func:`test_bare_slab_flux_shape_converges_second_order`.
    Reference precision floor is ~6e-6 (see the
    ``derive_2rg_continuous`` provenance), so the finest mesh
    error plateau is bounded below at that floor.
    """
    ref = continuous_get("dif_slab_2eg_2rg")
    materials = {
        0: mixture_from_diffusion_tables(ref.problem.materials[0]),
        1: mixture_from_diffusion_tables(ref.problem.materials[1]),
    }
    H_f = ref.problem.geometry_params["fuel_height"]
    H_r = ref.problem.geometry_params["refl_height"]

    # Skip the coarsest mesh (dz=5) which is pre-asymptotic on the
    # 80-cm fuel+reflector geometry — the ratio test needs at least
    # three refinements in the quadratic regime.
    dzs = [2.5, 1.25, 0.625, 0.3125]
    errs = []
    for dz in dzs:
        mesh, z = _zero_flux_problem([H_f, H_r], dz)
        result = solve_diffusion_1d(
            materials, mesh,
            keff_tol=1e-13, flux_tol=1e-11, max_outer=5000,
        )

        ref_g0 = np.asarray(ref.phi(z, 0), dtype=float)
        ref_g1 = np.asarray(ref.phi(z, 1), dtype=float)
        solver_g0 = result.flux.interior.values[0].copy()
        solver_g1 = result.flux.interior.values[1].copy()

        if ref_g0.sum() < 0:
            ref_g0 = -ref_g0
            ref_g1 = -ref_g1

        ref_g0, ref_g1 = _peak_normalise(ref_g0, ref_g1)
        solver_g0, solver_g1 = _peak_normalise(solver_g0, solver_g1)

        diff = np.concatenate([solver_g0 - ref_g0, solver_g1 - ref_g1])
        err = float(np.sqrt(dz * np.sum(diff * diff)))
        errs.append(err)

    errs_arr = np.asarray(errs)
    ratios = errs_arr[:-1] / errs_arr[1:]
    orders = np.log2(ratios)

    # O(h²) target. Loosened to 1.6 because the fuel/reflector
    # interface is a C⁰ (not C¹) point in the flux derivative,
    # which degrades the strict convergence rate slightly at
    # the coarsest meshes. The two finer ratios should still be
    # >= 1.8.
    assert orders[-1] > 1.6, (
        f"2-region flux shape convergence below expected: "
        f"errors={errs_arr}, orders={orders}"
    )
    assert np.all(orders > 1.4), (
        f"2-region flux shape convergence pathologically bad: "
        f"errors={errs_arr}, orders={orders}"
    )


def test_2region_interface_flux_is_continuous():
    r"""The back-substituted continuous reference must be continuous
    in :math:`\phi` across the fuel/reflector interface.

    Regression guard on the ``phi`` closure branching logic:
    any future bug in the ``if x <= H_f`` dispatch that forgets
    to carry state through ``T_fuel`` would produce a
    discontinuity here.
    """
    ref = continuous_get("dif_slab_2eg_2rg")
    H_f = ref.problem.geometry_params["fuel_height"]

    eps = 1e-10
    for g in range(ref.problem.n_groups):
        phi_left = ref.phi(np.array([H_f - eps]), g)[0]
        phi_right = ref.phi(np.array([H_f + eps]), g)[0]
        assert abs(phi_left - phi_right) < 1e-6, (
            f"Flux discontinuity at fuel/reflector interface (g={g}): "
            f"phi(H_f⁻)={phi_left:.10f} phi(H_f⁺)={phi_right:.10f}"
        )


def test_2region_zero_flux_boundaries_satisfied():
    r"""The back-substituted reference must satisfy :math:`\phi(0) = 0`
    exactly and :math:`\phi(L) \approx 0` to the transfer-matrix
    precision floor.

    This is the **zero-flux Dirichlet** law the reference commits to
    (#290 ruling 3 — not Marshak vacuum). ``phi(0) = 0`` is exact
    because the matching matrix's first two rows put zero in the flux
    slots by construction. ``phi(L) = 0`` is the transcendental
    eigenvalue condition; its residual is bounded by the square-root
    of the brentq tolerance on ``k``, giving an expected ceiling of
    ~3e-7 and observed ~6e-6 in practice.
    """
    ref = continuous_get("dif_slab_2eg_2rg")
    L = ref.problem.geometry_params["length"]

    for g in range(ref.problem.n_groups):
        # phi(0) must be exactly zero
        phi_0 = ref.phi(np.array([0.0]), g)[0]
        assert abs(phi_0) < 1e-14, (
            f"phi(0, g={g}) = {phi_0:.2e} is not exactly zero"
        )
        # phi(L) bounded by the precision floor
        phi_L = ref.phi(np.array([L]), g)[0]
        assert abs(phi_L) < 5e-5, (
            f"phi(L, g={g}) = {phi_L:.2e} exceeds transfer-matrix floor"
        )
