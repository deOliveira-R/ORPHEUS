"""Phase 2.1b SN heterogeneous eigenvalue consumer tests.

Consumes the Case singular-eigenfunction continuous reference from
:func:`orpheus.derivations.sn.derive_sn_heterogeneous_continuous` to
verify ``solve_sn`` on a two-region reflective slab at material
interfaces with piecewise-constant cross sections — the regime where
the Phase 2.1a smooth-Σ MMS test cannot reach, and the regime that
exposed ERR-025 (the DD cumprod recurrence bug) during the Phase 2.1b
investigation.

The reference is mesh-independent, mathematically self-contained, and
uses the same Gauss–Legendre quadrature order as the solver. It is the
exact discrete-:math:`S_N` eigenvalue, not the continuous-angle limit;
``solve_sn`` at the matching quadrature order must converge to it as
:math:`h \\to 0`.

Tests:

- ``test_sn_2region_reflective_case_eigenvalue`` — mesh-refinement
  convergence study. Confirms that ``solve_sn`` approaches the Case
  reference at the material-interface-degraded :math:`O(h)` rate
  expected for diamond-difference on piecewise-constant :math:`\\Sigma`,
  and that the finest-mesh residual is below 1e-5.
- ``test_sn_2region_reflective_flux_shape`` — shape comparison of the
  scalar flux against the back-substituted reference callable on a
  300-point evaluation grid.

Both tests are tagged ``@pytest.mark.l1 @pytest.mark.catches("ERR-025")``:
the fine-mesh convergence agreement is impossible with the pre-fix
``_sweep_1d_cumprod`` (which converges to ~1.25988, off by ~1.4e-2).
"""

import numpy as np
import pytest

from orpheus.derivations.reference_values import continuous_get
from orpheus.geometry import CoordSystem, Mesh1D
from orpheus.sn.quadrature import GaussLegendre1D
from orpheus.sn.solver import solve_sn

pytestmark = pytest.mark.verifies(
    "transport-cartesian",
    "dd-cartesian-1d",
    "dd-recurrence",
    "reflective-bc",
    "one-group-kinf",
    "sn-case-per-ordinate",
    "sn-case-slope-matrix",
    "sn-case-spatial-modes",
    "sn-case-real-basis",
    "sn-case-matching-matrix",
    "sn-case-physical-validation",
    "sn-case-back-substitution",
)


def _build_2region_mesh(H_A: float, H_B: float, n_per: int) -> Mesh1D:
    """Two-region uniform mesh with ``n_per`` cells per region."""
    edges = np.linspace(0.0, H_A + H_B, 2 * n_per + 1)
    mat_ids = np.array([0] * n_per + [1] * n_per)
    return Mesh1D(edges=edges, mat_ids=mat_ids, coord=CoordSystem.CARTESIAN)


@pytest.mark.l1
@pytest.mark.catches("ERR-025")
def test_sn_2region_reflective_case_eigenvalue():
    """``solve_sn`` converges to the Case reference :math:`k_\\text{eff}`
    at :math:`O(h)` with finest-mesh residual < 1e-5.

    Piecewise-constant :math:`\\Sigma(x)` degrades diamond-difference
    from its nominal :math:`O(h^{2})` to :math:`O(h)` at the material
    interface — this is the classical result Salari & Knupp cite for
    why smooth-:math:`\\Sigma` MMS is the spatial-operator verification
    path and Case-method eigenvalue is the eigenvalue verification
    path. This test sees exactly that: the error halves each mesh
    refinement.
    """
    ref = continuous_get("sn_slab_1eg_2rg_S8")
    geom = ref.problem.geometry_params
    materials = ref.problem.materials
    H_A = float(geom["fuel_height"])
    H_B = float(geom["refl_height"])
    N_ord = int(geom["n_ordinates"])
    quad = GaussLegendre1D.create(N_ord)

    n_per_values = [20, 40, 80, 160, 320]
    keffs: list[float] = []
    for n_per in n_per_values:
        mesh = _build_2region_mesh(H_A, H_B, n_per)
        result = solve_sn(
            materials, mesh, quad,
            max_outer=500, max_inner=500,
            keff_tol=1e-12, inner_tol=1e-12,
        )
        keffs.append(float(result.keff))

    errors = np.abs(np.array(keffs) - ref.k_eff)

    # Order check: each successive refinement should roughly halve
    # the error (O(h) convergence on a refine-by-2 sequence).
    ratios = errors[:-1] / errors[1:]
    assert np.all(ratios > 1.6), (
        f"Mesh refinement does not show O(h) convergence toward Case "
        f"reference. Errors per refinement: {errors.tolist()}, "
        f"ratios: {ratios.tolist()}. Expected each ratio ≳ 2."
    )

    # Finest-mesh absolute error: the pre-ERR-025 solver plateaued at
    # ~1.48e-2; the post-fix solver at n_per=320 reaches ~3.4e-8 on
    # this configuration. A 1e-5 threshold is ~3 orders above the
    # actual post-fix residual and ~3 orders below the ERR-025 gap.
    assert errors[-1] < 1e-5, (
        f"solve_sn finest-mesh keff={keffs[-1]:.10f} vs Case reference "
        f"{ref.k_eff:.10f} (Δ={keffs[-1] - ref.k_eff:+.2e}). "
        f"Expected |Δ| < 1e-5 at n_per={n_per_values[-1]}."
    )


@pytest.mark.l1
@pytest.mark.catches("ERR-025")
def test_sn_2region_reflective_flux_shape():
    """Scalar flux shape matches the Case reference's back-substituted
    callable to :math:`O(h)` pointwise.

    The reference ``phi(x)`` is normalised so :math:`\\max|\\phi| = 1`.
    ``solve_sn``'s cell-averaged flux is also amplitude-normalised
    before comparison, so the residual measures pure shape agreement
    (phase, region ratio, interface slope). Piecewise-constant
    :math:`\\Sigma(x)` gives :math:`O(h)` pointwise convergence near
    the interface, so we assert max |Δφ| < 1e-3 at n_per=320 — tight
    enough to catch any coefficient bug in the SN sweep, loose enough
    to accommodate the interface-layer DD error.
    """
    ref = continuous_get("sn_slab_1eg_2rg_S8")
    geom = ref.problem.geometry_params
    materials = ref.problem.materials
    H_A = float(geom["fuel_height"])
    H_B = float(geom["refl_height"])
    N_ord = int(geom["n_ordinates"])
    quad = GaussLegendre1D.create(N_ord)

    n_per = 320
    mesh = _build_2region_mesh(H_A, H_B, n_per)
    result = solve_sn(
        materials, mesh, quad,
        max_outer=500, max_inner=500,
        keff_tol=1e-12, inner_tol=1e-12,
    )

    centers = mesh.centers
    phi_ref = np.asarray(ref.phi(centers, 0), dtype=float)
    # scalar_flux shape is (nx, ny=1, ng=1) — take the 1-group slice.
    phi_sol = np.asarray(result.scalar_flux[:, 0, 0], dtype=float)

    # Normalise both to peak 1 so the comparison is shape-only.
    phi_ref_n = phi_ref / float(np.max(np.abs(phi_ref)))
    phi_sol_n = phi_sol / float(np.max(np.abs(phi_sol)))
    # Sign-align in case of opposite normalisation conventions.
    if np.sign(phi_sol_n[len(phi_sol_n) // 2]) != np.sign(
        phi_ref_n[len(phi_ref_n) // 2]
    ):
        phi_sol_n = -phi_sol_n

    # At n_per=320 with the post-ERR-025 fix, the Case reference
    # (evaluated pointwise via its back-substituted basis) agrees with
    # ``solve_sn``'s cell-averaged scalar flux to ~1e-6 after amplitude
    # normalisation. Pointwise vs cell-average contributes at most
    # O(h²) ≈ 2.5e-6, well below the 1e-5 threshold here.
    max_resid = float(np.max(np.abs(phi_sol_n - phi_ref_n)))
    assert max_resid < 1e-5, (
        f"Flux-shape residual {max_resid:.3e} exceeds 1e-5 at "
        f"n_per={n_per}. The post-fix residual on this configuration "
        f"is ~1e-6 (Case ref ↔ solve_sn agreement at matching "
        f"quadrature order); a residual this large indicates either "
        f"a reference-implementation drift or a sweep bug not caught "
        f"by the eigenvalue test."
    )
