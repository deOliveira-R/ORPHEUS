"""Verify the 1D SN solver: homogeneous exact, spatial O(h²), angular spectral."""

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


def _homogeneous_slab_mesh(n_cells: int, total_width: float, mat_id: int = 0) -> Mesh1D:
    """Single-region Cartesian mesh helper (replaces legacy ``homogeneous_1d``).

    Reflective BCs on both ends — the eigenvalue / lattice convention.
    """
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(Region(mat_id=mat_id, outer_thickness_cm=total_width),),
        bcs=(BC.reflective, BC.reflective),
    )
    return Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=n_cells),))


def _slab_fuel_moderator_mesh(
    n_fuel: int, n_mod: int, t_fuel: float, t_mod: float,
) -> Mesh1D:
    """Two-region fuel-moderator slab mesh (replaces legacy ``slab_fuel_moderator``).

    Reflective BCs on both ends. Material IDs: 2 = fuel, 0 = moderator.
    """
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
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import solve_sn

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


# ─── Homogeneous infinite medium (SN with reflective BCs) ────────────

@pytest.mark.parametrize("case_name", [
    "sn_slab_1eg_1rg",
    "sn_slab_2eg_1rg",
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


# ─── Heterogeneous: replaced by MMS continuous reference ─────────────
#
# Phase 2.1a of the verification campaign removed the legacy
# ``test_heterogeneous_convergence`` test that consumed
# ``sn_slab_Neg_Nrg`` (N > 1) references from
# ``orpheus.derivations.sn._derive_sn_heterogeneous``. Those
# references were Richardson-extrapolated from the SN solver
# itself (T3 circular self-verification) and have been deleted.
#
# The new heterogeneous SN spatial-operator verification lives in
# ``tests/sn/test_mms_heterogeneous.py`` and consumes the
# ``sn_mms_slab_2g_hetero`` Phase-0 ContinuousReferenceSolution
# from ``orpheus.derivations.sn_mms`` — the Method of Manufactured
# Solutions with smooth cross sections. See the heterogeneous MMS
# section of ``docs/theory/discrete_ordinates.rst`` for why.
#
# The eigenvalue-heterogeneous verification that the deleted test
# was nominally covering (but did not actually verify, because it
# compared the solver to its own extrapolant) will be restored in
# Phase 2.1b by a Case singular-eigenfunction reference.


# ─── Spatial convergence O(h²) ───────────────────────────────────────

def _convergence_order(values, spacings, reference):
    """Compute observed convergence order between successive refinements."""
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


# ─── L0 term verification of the DD cumprod recurrence (ERR-025) ────

@pytest.mark.l0
@pytest.mark.catches("ERR-025")
def test_dd_per_cell_recurrence_matches_symbolic_derivation():
    """Term-level verification that ``DiamondDifference.update``'s
    per-cell recurrence matches the symbolic derivation in
    :func:`orpheus.derivations.discrete.sn.balance.derive_cumprod_recurrence`.

    Issue #196 Step 2.5: ``_sweep_1d_cumprod`` was retired in favour
    of a unified fold over DAG visits that delegates per-cell
    algebra to :class:`DiamondDifference`.  This test directly
    invokes ``DiamondDifference.update`` on a synthetic per-cell
    input and checks the ``ψ_cell = ½(ψ_in + a·ψ_in + b·Q/W)``
    identity that gated the legacy ERR-025 coefficient drift.  The
    DD strategy's per-cell algebra IS the recurrence, so this is the
    most direct gate possible — no sweep / BC / boundary scaffolding
    in the way.
    """
    import sympy as sp
    from orpheus.derivations.discrete.sn.balance import derive_cumprod_recurrence
    from orpheus.geometry import CoordSystem, Mesh1D
    from orpheus.geometry.reduced_operator import slab_streaming
    from orpheus.sn.spatial.cell_update import CellVisit, UpstreamState
    from orpheus.sn.spatial.diamond import DiamondDifference

    # Symbolic coefficients, captured silently.
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        a_sym, b_sym = derive_cumprod_recurrence()

    mu_sym, dx_sym, Sig_t_sym, S_sym = sp.symbols(
        "mu dx Sigma_t S", positive=True
    )

    ng = 1
    sig_t_val = 1.5
    dx_val = 0.7
    Q_val = 3.0

    quad = Quadrature.gauss_legendre(4)
    edges = np.array([0.0, dx_val])
    mesh = Mesh1D(
        edges=edges,
        mat_ids=np.zeros(1, dtype=int),
        coord=CoordSystem.CARTESIAN,
    )
    op = slab_streaming(mesh, quad)
    W = quad.weights.sum()
    n_half = quad.N // 2
    mu_pos = np.abs(quad.mu_x[n_half:])

    # Synthetic ψ_in per positive ordinate.
    psi_in_per_ordinate = [0.4, 0.9]
    strat = DiamondDifference()
    total_xs = np.array([sig_t_val])

    for n in range(n_half):
        mu_val = float(mu_pos[n])
        direction_idx = n_half + n
        st = op.streaming_terms(cell_idx=0, direction_idx=direction_idx)

        # The contract source is Q · V · weight_norm = Q · dx / W.
        source = np.array([Q_val * dx_val / W])
        psi_in = np.array([psi_in_per_ordinate[n]])
        upstream = UpstreamState(
            spatial_upstream=psi_in, angular_upstream=None,
        )
        visit = CellVisit(
            cell_idx=0, streaming_terms=st, face_area_downstream=1.0,
        )
        result = strat.update(visit, total_xs, source, upstream)
        cell_avg_code = float(result.cell_average_flux[0])

        # Symbolic-derived reference.
        a_num = float(a_sym.subs({mu_sym: mu_val, dx_sym: dx_val,
                                  Sig_t_sym: sig_t_val}))
        b_num = float(b_sym.subs({mu_sym: mu_val, dx_sym: dx_val,
                                  Sig_t_sym: sig_t_val,
                                  S_sym: Q_val / W}))
        psi_out_expected = a_num * psi_in[0] + b_num
        cell_avg_expected = 0.5 * (psi_in[0] + psi_out_expected)

        # Issue #196 Step 2.5: DD's unified body computes
        # ``(source + numer_upstream)/denom`` (algebraically
        # identical to the cumprod's ``a·ψ_in + 2q/denom; ½(ψ_in +
        # ψ_out)`` but ULP-different).  Re-baseline rtol=1e-13.
        assert abs(cell_avg_code - cell_avg_expected) < 1e-12, (
            f"Ordinate n={n} (μ={mu_val:.6f}): "
            f"DD gave {cell_avg_code:.10e}, "
            f"derivation gives {cell_avg_expected:.10e}, "
            f"Δ={cell_avg_code - cell_avg_expected:+.2e}. "
            "DiamondDifference does not match "
            "sn_balance.derive_cumprod_recurrence."
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
    # Σ(x), and tight enough to catch the 1.4e-2 ERR-025 gap. Phase 2.1b
    # proper (the convergence study in
    # ``tests/sn/test_heterogeneous_transport.py``) asserts the actual
    # asymptotic agreement to ~1e-5.
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
