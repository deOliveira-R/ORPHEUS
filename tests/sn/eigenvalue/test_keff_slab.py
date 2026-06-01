"""Slab SN k-eigenvalue + spatial/angular convergence verification.

Split from the legacy ``tests/sn/test_cartesian.py`` (SN taxonomy
reorg): the eigenvalue claims (homogeneous-exact k_inf, heterogeneous
absolute k_eff against external references, spatial O(h²) and angular
spectral convergence). The per-cell DD recurrence L0 gate moved to
``tests/sn/sweep/slab/test_dd_recurrence.py``.
"""

import numpy as np
import pytest

from orpheus.derivations import get
from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import (
    BC,
    Mesh1D,
    Region,
    RegionMesh,
    StructuredGeometry,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import solve_sn

# Equation-coverage list preserved verbatim from the legacy
# test_cartesian module so no verifies(...) edge is lost in the split.
pytestmark = pytest.mark.verifies(
    "transport-cartesian",
    "dd-cartesian-1d",
    "dd-solve",
    "dd-recurrence",
    "multigroup",
    "reflective-bc",
    "one-group-kinf",
    "matrix-eigenvalue",
    "mg-balance",
)


def _homogeneous_slab_mesh(n_cells: int, total_width: float, mat_id: int = 0) -> Mesh1D:
    """Single-region Cartesian mesh helper (reflective BCs — eigenvalue convention)."""
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(Region(mat_id=mat_id, outer_thickness_cm=total_width),),
        bcs=(BC.reflective, BC.reflective),
    )
    return Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=n_cells),))


def _slab_fuel_moderator_mesh(
    n_fuel: int, n_mod: int, t_fuel: float, t_mod: float,
) -> Mesh1D:
    """Two-region fuel-moderator slab mesh (reflective BCs; 2 = fuel, 0 = moderator)."""
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(
            Region(mat_id=2, outer_thickness_cm=t_fuel),
            Region(mat_id=0, outer_thickness_cm=t_mod),
        ),
        bcs=(BC.reflective, BC.reflective),
    )
    return Mesh1D.from_geometry(geom, region_meshes=(
        RegionMesh(n_cells=n_fuel),
        RegionMesh(n_cells=n_mod),
    ))


def _convergence_order(values, spacings, reference):
    """Observed convergence order between successive refinements."""
    orders = []
    for i in range(1, len(values)):
        err_prev = abs(values[i - 1] - reference)
        err_curr = abs(values[i] - reference)
        if err_prev > 0 and err_curr > 0:
            orders.append(
                np.log(err_prev / err_curr)
                / np.log(spacings[i - 1] / spacings[i])
            )
    return orders


# ─── Homogeneous infinite medium (SN with reflective BCs) ────────────

@pytest.mark.parametrize("case_name", [
    "sn_slab_1eg_1rg",
    # Sentinel for the ``eigenvalue`` capability node: the cheapest
    # NON-DEGENERATE solver eigenvalue (2G; 1G k=νΣf/Σa is flux-shape
    # independent per vv-principles 1-group degeneracy). A flip here
    # localizes to the k-eigenvalue / source-iteration outer loop.
    # See .claude/plans/sn_sentinel_harness.md.
    pytest.param("sn_slab_2eg_1rg", marks=pytest.mark.sentinel),
    "sn_slab_4eg_1rg",
])
def test_homogeneous_exact(case_name):
    """SN 1D with reflective BCs on a homogeneous slab must match
    the analytical infinite-medium eigenvalue."""
    case = get(case_name)
    mix = next(iter(case.materials.values()))
    materials = {0: mix}
    mesh = _homogeneous_slab_mesh(20, 2.0, mat_id=0)
    quad = Quadrature.gauss_legendre(8)
    result = solve_sn(materials, mesh, quad,
                      max_inner=500, inner_tol=1e-10)

    assert abs(result.keff - case.k_inf) < 1e-8, (
        f"keff={result.keff:.8f} vs analytical={case.k_inf:.8f}"
    )


# ─── Spatial convergence O(h²) ───────────────────────────────────────

@pytest.mark.l1
def test_spatial_convergence():
    """Diamond-difference scheme must show O(h²) spatial convergence."""
    fuel = get_mixture("A", "1g")
    mod = get_mixture("B", "1g")
    materials = {2: fuel, 0: mod}
    t_fuel, t_mod = 0.5, 0.5

    keffs = []
    dxs = []
    for n_per in [5, 10, 20, 40]:
        mesh = _slab_fuel_moderator_mesh(
            n_fuel=n_per, n_mod=n_per, t_fuel=t_fuel, t_mod=t_mod,
        )
        quad = Quadrature.gauss_legendre(16)
        result = solve_sn(
            materials, mesh, quad,
            max_outer=300, max_inner=500, inner_tol=1e-10,
        )
        keffs.append(result.keff)
        dxs.append(t_fuel / n_per)

    # Richardson extrapolation reference
    k_ref = keffs[-1] + (keffs[-1] - keffs[-2]) / 3.0
    orders = _convergence_order(keffs, dxs, k_ref)

    assert orders[-1] > 1.7, (
        f"Expected O(h²) convergence, got order {orders[-1]:.2f}"
    )


# ─── Heterogeneous absolute eigenvalue regression (ERR-025) ──────────

@pytest.mark.l1
@pytest.mark.catches("ERR-025")
def test_heterogeneous_absolute_keff():
    """2-region A+B reflective slab must match external references.

    Regression for the DD face-flux recurrence bug (ERR-025): the
    legacy ``_sweep_1d_cumprod`` used the wrong recurrence coefficients

        a = 2μ / (2μ + Δx·Σ_t)          (wrong)
        s = 0.5·Δx·Q / (2μ + Δx·Σ_t)    (wrong)

    instead of those derived in
    ``orpheus.derivations.discrete.sn.balance.derive_cumprod_recurrence``

        a = (2μ − Δx·Σ_t) / (2μ + Δx·Σ_t)
        b = 2·Δx·(Q/W) / (2μ + Δx·Σ_t)

    Both wrongs are off by a factor of two in opposite directions and
    cancelled for eigenvalue problems with a single material (a scale
    on φ leaves k = ν·Σ_f·φ / Σ_a·φ invariant). At material interfaces
    the cancellation broke because the factor depended on Σ_t(x),
    shifting k_eff by ~1.4e-2.

    Two independent references — a Case singular-eigenfunction
    expansion and the ORPHEUS CP slab solver (E₃ kernel) — agree on
    k ≈ 1.27461 for this config, against which the fixed DD
    recurrence now matches to 5e-5 at n_per=320.
    """
    from orpheus.derivations.reference_values import continuous_get
    from orpheus.geometry import Mesh1D, CoordSystem

    ref = continuous_get("sn_slab_1eg_2rg_S8")
    geom = ref.problem.geometry_params
    materials = ref.problem.materials
    H_A = float(geom["fuel_height"])
    H_B = float(geom["refl_height"])
    N_ord = int(geom["n_ordinates"])

    n_per = 320
    edges = np.linspace(0.0, H_A + H_B, 2 * n_per + 1)
    mat_ids = np.array([0] * n_per + [1] * n_per)
    mesh = Mesh1D(edges=edges, mat_ids=mat_ids, coord=CoordSystem.CARTESIAN)
    quad = Quadrature.gauss_legendre(N_ord)

    result = solve_sn(
        materials, mesh, quad,
        max_outer=500, max_inner=500,
        keff_tol=1e-11, inner_tol=1e-11,
    )

    # 5e-4 tolerance is loose enough to accommodate the O(h) spatial
    # truncation at material interfaces with DD on a piecewise-constant
    # Σ(x), and tight enough to catch the 1.4e-2 ERR-025 gap.
    assert abs(result.keff - ref.k_eff) < 5e-4, (
        f"keff={result.keff:.10f} vs Case reference={ref.k_eff:.10f} "
        f"(Δ={result.keff - ref.k_eff:+.2e})"
    )


# ─── Angular spectral convergence ────────────────────────────────────

@pytest.mark.l1
def test_angular_convergence():
    """Gauss-Legendre quadrature must show spectral convergence in angle."""
    fuel = get_mixture("A", "1g")
    mod = get_mixture("B", "1g")
    materials = {2: fuel, 0: mod}

    keffs = []
    n_ords = [4, 8, 16, 32]
    for N in n_ords:
        mesh = _slab_fuel_moderator_mesh(
            n_fuel=40, n_mod=40, t_fuel=0.5, t_mod=0.5,
        )
        quad = Quadrature.gauss_legendre(N)
        result = solve_sn(
            materials, mesh, quad,
            max_outer=300, max_inner=500, inner_tol=1e-10,
        )
        keffs.append(result.keff)

    k_ref = keffs[-1]
    orders = _convergence_order(keffs, [1 / N for N in n_ords], k_ref)
    assert len(orders) >= 2
    assert max(orders[:-1]) > 1.5, (
        f"Expected spectral convergence, got orders {orders}"
    )
