r"""L1 MMS — Linear-Discontinuous (LD) SN on a 1-D slab, via the affine scan.

LD is the first non-Diamond-Difference spatial scheme (#158).  Since Increment B
it is **affine-scannable** (the Schur complement eliminates the slope, leaving
the single-upstream recurrence ``ψ_out = a·ψ_in + b``), so a 1-D slab LD mesh
rides the fast DAG-free ``CumprodScan`` via the coefficient model — it supplies
``(a, inverse_denom, w)`` through ``affine_scan_coefficients`` and the generic
``affine_closure`` ops do the rest.  The polymorphic ``FullFieldWavefront`` DAG
oracle (Increment A, via ``cell_kernel_batch``) remains the verification
reference — the two-paths gate pins ``CumprodScan``-LD ≡ ``FullFieldWavefront``-LD.
This file pins:

* **Routing** — a 1-D slab LD mesh dispatches to ``CumprodScan`` by
  ``default_for`` (and a DD mesh still picks ``CumprodScan`` too — both affine,
  no selection regression).
* **Spatial order** — LD converges :math:`\mathcal{O}(h^2)` on the manufactured
  slab problem, run end-to-end through
  ``solve_sn_fixed_source(..., cell_update=LinearDiscontinuous())``.
* **Two-paths (#158 Inc B)** — ``CumprodScan``-LD ≡ ``FullFieldWavefront``-LD
  (principled-equivalent at nULP; the ×V scan vs ÷V kernel conventions agree).
* **Matvec ≡ sweep** — the Krylov apply (LD's group-2 ``residual_kernel_batch``
  via ``_apply_walk``) matches the SI sweep.

The per-cell exactness-on-linears (the strong structurally-independent
correctness oracle) lives in ``tests/sn/spatial/test_linear_discontinuous.py``.

Related: ``orpheus.sn.spatial.linear_discontinuous`` (the occupant);
``orpheus.sn.spatial.affine_closure`` (the generic coefficient-model ops);
``orpheus.sn.loss_representation.{CumprodScan, FullFieldWavefront}``;
``.claude/plans/mellow-swinging-breeze.md`` (Increment B).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.mms.sn import build_1d_slab_mms_case
from orpheus.sn import solve_sn_fixed_source
from orpheus.sn.geometry import SNMesh
from orpheus.sn.loss_representation import (
    CumprodScan,
    FullFieldWavefront,
    default_for,
)
from orpheus.sn.spatial import LinearDiscontinuous


def _l2_error(phi_num: np.ndarray, phi_ref: np.ndarray, widths: np.ndarray) -> float:
    diff = phi_num - phi_ref
    return float(np.sqrt(np.sum(widths * diff * diff)))


@pytest.mark.foundation
def test_ld_slab_mesh_routes_to_cumprod_scan() -> None:
    """A 1-D slab LD mesh dispatches to the fast affine scan (#158 Inc B); DD
    stays on the scan too (both affine — no selection regression)."""
    case = build_1d_slab_mms_case()
    mesh = case.build_mesh(16)
    ld_mesh = SNMesh(mesh, case.quadrature, case.materials,
                     cell_update=LinearDiscontinuous())
    dd_mesh = SNMesh(mesh, case.quadrature, case.materials)
    if not isinstance(default_for(ld_mesh), CumprodScan):
        pytest.fail(
            "LD slab mesh routed to "
            f"{type(default_for(ld_mesh)).__name__}, expected CumprodScan "
            "(LD is affine-scannable since Increment B)"
        )
    if not isinstance(default_for(dd_mesh), CumprodScan):
        pytest.fail(
            "DD slab mesh no longer routes to CumprodScan "
            f"({type(default_for(dd_mesh)).__name__}) — selection regression"
        )


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("ld-cartesian-1d", "ld-slab", "transport-cartesian")
def test_sn_1d_slab_ld_mms_converges_second_order() -> None:
    r"""LD SN on a 1-D slab shows measured :math:`\mathcal{O}(h^2)` end-to-end.

    Drives the manufactured slab problem through
    ``solve_sn_fixed_source(cell_update=LinearDiscontinuous())`` — which routes
    to ``FullFieldWavefront`` calling LD's ``cell_kernel_batch`` — at four mesh
    refinements.  Asserts the order converges to 2 (the finest pair > 1.95; the
    pre-asymptotic coarse pair > 1.85) AND the converged value (the magnitude
    band against the manufactured ``phi_exact`` — rate alone is necessary, never
    sufficient, vv-principles §5).
    """
    case = build_1d_slab_mms_case()
    n_cells = [20, 40, 80, 160]
    errors = []
    for nc in n_cells:
        mesh = case.build_mesh(nc)
        Q = case.external_source(mesh)
        result = solve_sn_fixed_source(
            case.materials, mesh, case.quadrature, Q,
            boundary_condition="vacuum", max_inner=500, inner_tol=1e-13,
            cell_update=LinearDiscontinuous(),
        )
        phi_num = result.scalar_flux.values[0, :]
        phi_ref = case.phi_exact(mesh.centers)
        errors.append(_l2_error(phi_num, phi_ref, mesh.widths))

    errors = np.asarray(errors)
    orders = np.log2(errors[:-1] / errors[1:])
    # O(h²) converging to 2: finest pair clears 1.95, all pairs clear 1.85
    # (the smooth sin ansatz is mildly pre-asymptotic on the coarse pair; the
    #  strong exactness-on-linears oracle is the per-cell gate-1 companion).
    assert orders[-1] > 1.95, f"LD finest order {orders[-1]:.3f} ≯ 1.95 ({orders})"
    assert np.all(orders > 1.85), f"LD order dropped below 1.85: {orders}"
    assert 1e-9 < errors[-1] < 1e-2, (
        f"LD finest-mesh error {errors[-1]:.3e} outside the value band — "
        "rate may be O(h²) to the WRONG limit (vv-principles §5)"
    )


@pytest.mark.l1
@pytest.mark.verifies("ld-cartesian-1d")
def test_sn_1d_slab_ld_mms_krylov_matches_si() -> None:
    r"""LD forward matvec ≡ SI sweep (the L14 matvec≡sweep leg).

    Since Increment B LD routes to ``CumprodScan``: the Krylov inner exercises
    LD's ``residual_kernel_batch`` via the scan walk's ``_apply_walk`` cartesian
    arm (the apply twin of the SI sweep's group-3 ``affine_scan_coefficients`` +
    generic ``cell_average``); both must converge to the SAME manufactured flux.
    """
    case = build_1d_slab_mms_case()
    mesh = case.build_mesh(80)
    Q = case.external_source(mesh)
    si = solve_sn_fixed_source(
        case.materials, mesh, case.quadrature, Q, boundary_condition="vacuum",
        max_inner=500, inner_tol=1e-13, cell_update=LinearDiscontinuous(),
    )
    kry = solve_sn_fixed_source(
        case.materials, mesh, case.quadrature, Q, boundary_condition="vacuum",
        max_inner=500, inner_tol=1e-13, inner_solver="krylov",
        cell_update=LinearDiscontinuous(),
    )
    np.testing.assert_allclose(
        kry.scalar_flux.values[0], si.scalar_flux.values[0],
        rtol=1e-9, atol=1e-11,
    )


@pytest.mark.l1
@pytest.mark.verifies("ld-cartesian-1d")
def test_ld_two_paths_scan_equals_dag_oracle() -> None:
    r"""#158 Inc B headline — ``CumprodScan``-LD ≡ ``FullFieldWavefront``-LD.

    The fast affine scan (``CumprodScan``, group-3 ``affine_scan_coefficients``
    in the ×V convention) and the polymorphic DAG oracle
    (``FullFieldWavefront``, group-2 ``cell_kernel_batch`` in the ÷V ``g=|μ|/Δ``
    convention) solve the SAME LD discrete system via two different schedules
    AND two different scalings.  They are **principled-equivalent** (NOT
    byte-identical: ``S_scan = V·S_kernel`` + distinct reduction trees), so the
    two-paths cross-check pins LD's affine coefficients against the trusted
    Increment-A kernel — the system-level analogue of the unit group2≡group3
    gate, and the slope-row sign-trap catcher at the schedule level.

    STRESSING config (Mode-9 — a flat/homogeneous two-paths gate would be
    FP-coincidentally exact): heterogeneous per-cell-per-group ``σ_t`` (≥2G) +
    a non-flat random per-ordinate source.  One sweep each, on a fresh boundary
    (the sweep mutates it in place).
    """
    from orpheus.geometry import BC, CoordSystem, Mesh1D
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.transport.fields.boundary_flux import BoundaryFlux
    from orpheus.derivations.continuous.mms.sn import _make_2g_asymmetric_mixture

    nx = 24
    materials = {0: _make_2g_asymmetric_mixture(
        np.array([1.0, 1.5]), np.array([[0.3, 0.2], [0.0, 0.6]]),
    )}
    mesh = Mesh1D(
        edges=np.linspace(0.0, 5.0, nx + 1), mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN, bc_left=BC("vacuum"), bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(8)
    ld_mesh = SNMesh(mesh, quad, materials, cell_update=LinearDiscontinuous())
    N, ng = quad.N, ld_mesh.ng

    rng = np.random.default_rng(20260614)
    sig_t = rng.uniform(0.3, 3.0, size=(ng, nx))          # het, ≥2G
    Q = rng.standard_normal((N, ng, nx))                  # non-flat per-ordinate

    _, phi_scan = CumprodScan(ld_mesh).sweep(
        Q, sig_t, BoundaryFlux.zeros_on(ld_mesh),
    )
    _, phi_dag = FullFieldWavefront(ld_mesh).sweep(
        Q, sig_t, BoundaryFlux.zeros_on(ld_mesh),
    )
    # Principled-equivalent (×V scan vs ÷V kernel): tight nULP-scale band, well
    # below any algorithmic-difference signature (a sign-trap would be O(1)).
    np.testing.assert_allclose(phi_scan, phi_dag, rtol=1e-12, atol=1e-13)


@pytest.mark.foundation
def test_ld_curvilinear_scan_rejected() -> None:
    """Slab-only guard (#158 Inc B): a 1-D *curvilinear* LD mesh would match
    ``CumprodScan.supports`` (``is_1d and is_affine_scannable``), but the
    curvilinear LD scan closure is unpublished — ``affine_scan_coefficients``
    must raise (fail-fast at the ``CollisionCache`` build in
    ``SNSolver.__init__``) rather than silently run DD-shaped curvature math."""
    from orpheus.geometry import BC, CoordSystem, Mesh1D
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.derivations.continuous.mms.sn import _make_1g_mixture

    nx = 8
    materials = {0: _make_1g_mixture(1.0, 0.5)}
    sphere = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1), mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL, bc_left=BC("reflective"), bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(8)
    Q = np.ones((quad.N, 1, nx)) / quad.weights.sum()
    with pytest.raises(NotImplementedError, match="slab/Cartesian"):
        solve_sn_fixed_source(
            materials, sphere, quad, Q, boundary_condition="vacuum",
            cell_update=LinearDiscontinuous(),
        )


@pytest.mark.l1
@pytest.mark.xfail(
    strict=False,
    reason="#158 Increment C: flat-source LD (Q̂=0) loses the diffusion limit; "
    "needs the canonical slope source Σ_s·φ̂ in the iterate. DD is "
    "interior-diffusion-consistent (LMM-1987); flat-source LD is not yet — "
    "this tripwire flips to xpass when Increment C threads the moment source.",
)
def test_ld_thick_diffusive_limit_xfail() -> None:
    r"""Forward tripwire: thick diffusive slab — flat-source LD ≠ DD interior.

    On an optically thick, scattering-dominated mesh, DD reproduces the interior
    diffusion solution (LMM-1987 Table I); the flat-source LD shipped in
    Increment A does NOT (the slope source — incl. the scattering slope — is
    Increment C).  This gate ASPIRES to LD ≈ DD; it currently FAILS (xfail) and
    will flip to xpass when the canonical source lands.  Documents the
    flat-source limitation as a GATE, not only a docstring caution.
    """
    from orpheus.derivations.continuous.mms.sn import _make_1g_mixture
    from orpheus.geometry import BC, CoordSystem, Mesh1D
    from orpheus.numerics.quadrature import Quadrature

    nx, length, sigma_t, c = 4, 1.0, 40.0, 0.99   # σ_t·h = 10/cell (thick), c→1
    materials = {0: _make_1g_mixture(sigma_t, c * sigma_t)}
    mesh = Mesh1D(
        edges=np.linspace(0.0, length, nx + 1), mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN, bc_left=BC("vacuum"), bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(8)
    Q = np.ones((quad.N, 1, nx)) / quad.weights.sum()   # uniform iso source
    dd = solve_sn_fixed_source(
        materials, mesh, quad, Q, boundary_condition="vacuum",
        inner_solver="krylov", max_inner=2000, inner_tol=1e-10,
    )
    ld = solve_sn_fixed_source(
        materials, mesh, quad, Q, boundary_condition="vacuum",
        inner_solver="krylov", max_inner=2000, inner_tol=1e-10,
        cell_update=LinearDiscontinuous(),
    )
    mid = nx // 2
    dd_mid = float(dd.scalar_flux.values[0, mid])
    ld_mid = float(ld.scalar_flux.values[0, mid])
    rel = abs(ld_mid - dd_mid) / abs(dd_mid)
    assert rel < 0.05, (
        f"flat-source LD midpoint {ld_mid:.4f} vs DD {dd_mid:.4f} "
        f"(rel {rel:.1%}) — the diffusion limit needs the slope source (Inc. C)"
    )
