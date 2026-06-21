r"""C5.5 (#225) — d=3 end-to-end value gates through the production entries.

The 3-axis SNMesh runs on the d-generic ``FullFieldWavefront`` oracle
spine from day one. These gates are the acceptance criteria for the
admission, each driven through the production ``solve_sn`` /
``solve_sn_fixed_source`` entry (C5-G17 — the axes tuple IS the 3-D
entry surface):

* **k_inf 3-D ≡ 2-D ≡ 1-D** (the headline, L0 eigenvalue → closed
  form): k_inf = λ_max(A⁻¹F) with A = diag(Σ_t) − Σ_s₀ᵀ is
  flux-shape-independent — the same value at every dimension on a
  homogeneous reflective box, computed by a reference that never
  touches the streaming operator (structurally independent). ≥2 groups
  ONLY (1G is degenerate — k = νΣ_f/Σ_a regardless of flux shape).
* **Infinite-medium flux shape** (L1 → closed form): pure absorber
  per-ordinate ψ_{n,g} = Q_g/(W·Σ_{t,g}) — DD is flat-flux exact and
  c=0 converges immediately, so the bound is solver-tolerance-tight;
  PLUS the scattering companion φ = (diag(Σ_t) − Σ_s₀ᵀ)⁻¹Q (group
  coupling active — a Mode-6 convention-drift catcher the pure
  absorber is blind to; convergence-limited tolerance).
* **Mode-9 G-S ≡ Jacobi FP-invariance** (vv-principles Mode 9): the
  boundary-G-S splitting must converge to the Jacobi fixed point on a
  config that BREAKS the degenerate coincidences — vacuum/reflective
  MIXED faces, nx≠ny≠nz, heterogeneous 2G, DIAGONAL cubature
  (level-symmetric: shared faces across octants — the ERR-056
  discipline). An all-reflective isotropic box would pass with a wrong
  splitting.

Mode-2 asymmetry is everywhere: distinct per-axis cell counts AND
extents, per-group distinct sources, axis-asymmetric BCs — a mu_y↔mu_z
swap or a transposed mat_map is a detectable difference, not a
symmetric coincidence.

Assertions are ``np.testing`` only (Mode-8: bare asserts are inert
under the canonical ``python -O`` invocation).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations import get
from orpheus.derivations.common.xs_library import get_mixture, make_mixture
from orpheus.geometry import BC, CoordSystem, Mesh1D, Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.axis import AxisMesh
from orpheus.sn.solver import (
    _apply_default_bcs,
    _maybe_window,
    _select_si_resolvent,
    _GaussSeidelResolvent,
    solve_sn,
    solve_sn_fixed_source,
)
from orpheus.sn.geometry import SNMesh


def _d3_axes(extents=(1.0, 2.0, 3.0), cells=(3, 4, 5), bcs=None):
    """Distinct counts AND extents per axis (Mode-2 asymmetry)."""
    bcs = bcs or [(None, None)] * 3
    return tuple(
        AxisMesh(
            edges=np.linspace(0.0, ext, n + 1), bc_low=lo, bc_high=hi,
        )
        for ext, n, (lo, hi) in zip(extents, cells, bcs)
    )


# ─── (b)(iii) THE HEADLINE — k_inf 3-D ≡ 2-D ≡ 1-D ───────────────────────


@pytest.mark.l1
@pytest.mark.verifies("matrix-eigenvalue", "multigroup")
@pytest.mark.parametrize("ng_key", ["2g", "4g"])
def test_kinf_3d_equals_2d_equals_1d_homogeneous_reflective(ng_key) -> None:
    """k_inf on a homogeneous all-reflective box at d=1, 2, 3.

    Closed-form reference: ``case.k_inf`` (matrix eigenvalue, never
    touches the sweep). All-reflective with no explicit BCs — the
    SNMesh reflective default — on every surface. Cell counts are
    deliberately small AND asymmetric (k_inf is mesh-independent, so
    the asymmetry is free Mode-2 insurance).
    """
    case = get(f"sn_slab_{ng_key[0]}eg_1rg")
    mix = next(iter(case.materials.values()))
    materials = {0: mix}
    quad3 = Quadrature.level_symmetric(sn_order=4)

    sol3 = solve_sn(
        materials,
        _d3_axes(extents=(0.5, 0.75, 0.5), cells=(2, 3, 2)),
        quad3, keff_tol=1e-10, inner_tol=1e-11,
    )
    sol2 = solve_sn(
        materials,
        Mesh2D(
            edges_x=np.linspace(0.0, 0.5, 3),
            edges_y=np.linspace(0.0, 0.75, 4),
            mat_map=np.zeros((2, 3), dtype=int),
        ),
        quad3, keff_tol=1e-10, inner_tol=1e-11,
    )
    sol1 = solve_sn(
        materials,
        Mesh1D(
            edges=np.linspace(0.0, 0.5, 3),
            mat_ids=np.zeros(2, dtype=int),
            coord=CoordSystem.CARTESIAN,
        ),
        Quadrature.gauss_legendre(n_ordinates=8),
        keff_tol=1e-10, inner_tol=1e-11,
    )

    for label, sol in (("3-D", sol3), ("2-D", sol2), ("1-D", sol1)):
        np.testing.assert_allclose(
            sol.keff, case.k_inf, atol=1e-8, rtol=0,
            err_msg=f"{ng_key} {label}: keff vs closed-form k_inf",
        )
    np.testing.assert_allclose(sol3.keff, sol2.keff, atol=1e-8, rtol=0)
    np.testing.assert_allclose(sol2.keff, sol1.keff, atol=1e-8, rtol=0)


# ─── (b)(i) infinite-medium flux shape — closed form ─────────────────────


@pytest.mark.l1
@pytest.mark.verifies("transport-cartesian")
def test_d3_pure_absorber_per_ordinate_psi_exact() -> None:
    """Pure absorber, all-reflective: ψ_{n,g} = Q_g/(W·Σ_{t,g}) per ordinate.

    The sharpest Mode-1/3/4 probe: DD is exact for flat flux and c=0
    needs no iteration, so EVERY ordinate must carry the closed-form
    value to solver-tolerance precision. Per-group distinct Q and
    distinct Σ_t make a group swap observable.
    """
    quad = Quadrature.level_symmetric(sn_order=4)
    mix = make_mixture(
        sig_t=np.array([0.8, 1.6]),
        sig_c=np.array([0.8, 1.6]),       # pure absorber: Σ_c = Σ_t
        sig_f=np.array([0.0, 0.0]),
        nu=np.array([0.0, 0.0]),
        chi=np.zeros(2),  # non-fissile ⇒ null spectrum (S10a __post_init__ guard)
        sig_s=np.zeros((2, 2)),
    )
    Q_g = np.array([1.0, 0.5])
    W = float(np.sum(quad.weights))
    q = np.broadcast_to(
        (Q_g / W)[None, :, None, None, None], (quad.N, 2, 3, 4, 5),
    ).copy()

    sol = solve_sn_fixed_source(
        {0: mix}, _d3_axes(), quad,
        external_source=q, boundary_condition="reflective",
        inner_tol=1e-13,
    )
    psi = np.asarray(sol.angular_flux.bulk.values)   # (N, ng, 3, 4, 5)
    sig_t = np.asarray(mix.SigT)
    for g in range(2):
        np.testing.assert_allclose(
            psi[:, g], Q_g[g] / (W * sig_t[g]), rtol=1e-10,
            err_msg=f"per-ordinate flat-flux identity, group {g}",
        )
    phi = np.asarray(sol.scalar_flux.values)
    for g in range(2):
        np.testing.assert_allclose(phi[g], Q_g[g] / sig_t[g], rtol=1e-10)


@pytest.mark.l1
@pytest.mark.verifies("transport-cartesian", "multigroup")
def test_d3_scattering_infinite_medium_matches_multigroup_balance() -> None:
    """Scattering medium, all-reflective: φ = (diag(Σ_t) − Σ_s₀ᵀ)⁻¹ Q.

    The group-coupling companion (Mode-6 convention-drift catcher —
    Σ_s vs Σ_sᵀ is observable because mixture C's scattering matrix is
    asymmetric). Tolerance is SI-convergence-limited, not
    discretization-limited (DD stays flat-flux exact): the iterate-
    delta criterion under-estimates the true error by ~1/(1−c).
    """
    quad = Quadrature.level_symmetric(sn_order=4)
    mix = get_mixture("C", "2g")
    Q_g = np.array([1.0, 0.5])
    W = float(np.sum(quad.weights))
    q = np.broadcast_to(
        (Q_g / W)[None, :, None, None, None], (quad.N, 2, 3, 4, 5),
    ).copy()

    sol = solve_sn_fixed_source(
        {0: mix}, _d3_axes(), quad,
        external_source=q, boundary_condition="reflective",
        inner_tol=1e-13,
    )
    phi = np.asarray(sol.scalar_flux.values)
    sig_s0 = np.asarray(mix.SigS[0].todense())
    A = np.diag(np.asarray(mix.SigT)) - sig_s0.T
    phi_exact = np.linalg.solve(A, Q_g)
    for g in range(2):
        np.testing.assert_allclose(
            phi[g], phi_exact[g], rtol=1e-7,
            err_msg=f"multigroup infinite-medium balance, group {g}",
        )


# ─── (b)(ii) Mode-9 — G-S ≡ Jacobi on the degenerate-breaking box ────────


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.catches("ERR-056")
@pytest.mark.verifies("transport-cartesian", "multigroup")
def test_d3_gauss_seidel_jacobi_fixed_point_invariance() -> None:
    """Boundary-G-S and Jacobi converge to the SAME d=3 fixed point.

    The config breaks every degenerate coincidence (vv Mode 9):
    vacuum/reflective MIXED faces (x-refl / y-vac / z-refl —
    axis-asymmetric, so a wrong reflection partner shifts the answer),
    nx≠ny≠nz, heterogeneous 2G fuel|moderator split across x (non-flat
    flux activates the wavefront redistribution), and a DIAGONAL
    cubature (level-symmetric — faces shared across octant groups, the
    ERR-056 reflect-order discipline). A splitting bug that is exact on
    the reflective isotropic box is WRONG here.
    """
    quad = Quadrature.level_symmetric(sn_order=4)
    refl, vac = BC("reflective"), BC("vacuum")
    axes = _d3_axes(
        extents=(2.0, 1.5, 1.0), cells=(5, 3, 4),
        bcs=[(refl, refl), (vac, vac), (refl, refl)],
    )
    mat_map = np.zeros((5, 3, 4), dtype=int)
    mat_map[:2] = 2
    materials = {2: get_mixture("A", "2g"), 0: get_mixture("B", "2g")}

    def _keff_flux(schedule: str):
        sol = solve_sn(
            materials, axes, quad,
            keff_tol=1e-10, inner_tol=1e-11, inner_schedule=schedule,
            mat_map=mat_map,
        )
        return sol.keff, np.asarray(sol.scalar_flux.values)

    k_j, phi_j = _keff_flux("jacobi")
    k_gs, phi_gs = _keff_flux("gauss_seidel")

    # Non-degeneracy guard: the flux must be genuinely non-flat.
    profile = phi_j[0]
    if float(profile.max() / profile.min()) <= 1.2:
        pytest.fail("degenerate (near-flat) flux — Mode-9 config broken")

    np.testing.assert_allclose(k_gs, k_j, atol=1e-8, rtol=0)
    np.testing.assert_allclose(
        phi_gs / phi_gs.mean(), phi_j / phi_j.mean(), rtol=1e-6, atol=1e-8,
    )


# ─── C5-G15/G16 real-mesh twins + C5-G18 entry-surface semantics ─────────


@pytest.mark.foundation
def test_d3_real_mesh_window_passthrough_and_gs_admissible() -> None:
    """Real-mesh twins of the C5.4 duck-typed gate pins.

    d=3: NO moment windowing (the in-sweep moment kernel is 2-D-only)
    and boundary-G-S IS admissible (the schedule + scheduled sweep are
    d-generic; FP-invariance proven by the Mode-9 gate above).
    """
    quad = Quadrature.level_symmetric(sn_order=4)
    mats = {0: get_mixture("A", "2g")}
    sn = SNMesh.from_axes(_d3_axes(), quad, mats)

    base = object()
    scattering_stub = type(
        "S", (), {"quadrature": quad, "scattering_order": 0},
    )()
    wrapped, windowed = _maybe_window(base, scattering_stub, sn)
    np.testing.assert_equal(windowed, False)
    np.testing.assert_equal(wrapped is base, True)

    from types import SimpleNamespace
    lc_stub = SimpleNamespace(sn_mesh=sn)
    resolvent, gains = _select_si_resolvent(
        lc_stub, "S_gain", "B_gain", sn, "gauss_seidel",
    )
    np.testing.assert_equal(isinstance(resolvent, _GaussSeidelResolvent), True)
    np.testing.assert_equal(gains, ("S_gain",))


@pytest.mark.foundation
def test_axes_default_bc_semantics() -> None:
    """C5-G18: the entry default fills axes ONLY when no BC is declared."""
    vac = BC("vacuum")
    bare = _d3_axes()
    filled = _apply_default_bcs(bare, "vacuum")
    for ax in filled:
        np.testing.assert_equal(
            [b.kind for b in ax.bc.values()], ["vacuum", "vacuum"],
        )
    partially = _d3_axes(bcs=[(vac, None), (None, None), (None, None)])
    np.testing.assert_equal(
        _apply_default_bcs(partially, "reflective") is partially, True,
    )
