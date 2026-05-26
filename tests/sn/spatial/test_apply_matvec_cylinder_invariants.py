r"""L0 apply-matvec flat-flux invariant — cylinder (Phase G Step 2 cylinder fix).

Issue #196 Phase G Step 2 cylinder fix (2026-05-13).  Promoted from the
numerics-investigator's diagnostic
``derivations/diagnostics/diag_phase_g_step2_cyl_residual_pytest.py``.

Test premise
------------

On flat angular flux ``ψ_flat = const`` everywhere on a homogeneous
reflective cylinder, the curvilinear apply-matvec must satisfy

.. math::

   L \cdot \psi_{\text{flat}} = \Sigma_t \cdot \psi_{\text{flat}}

per ordinate, per cell.  Streaming + angular redistribution telescope
to zero on a flat field (Bailey's α-dome invariant + Pomraning
structural-singularity isotropy preservation).  The remaining term is
the collision operator ``Σ_t · ψ``.

Why this catches ERR-048 manifestation #3
-----------------------------------------

The pre-fix ``CarlsonInwardSweep`` ``Q_bar = 0.5 · Σ_t · φ_0`` hardcoded
the sphere quadrature convention ``Σw = 2`` (the literal ``0.5`` is
``1/Σw_sphere``).  For cylinder ProductQuadrature where the per-level
weight sum is ``2π``, the literal is wrong by factor ``2π ≈ 6.28``,
producing a divergent residual that **grows** with mesh refinement
(Signature 1 fingerprint from ``numerical-bug-signatures`` SKILL.md).

The fix ``Q_bar = Σ_t · φ_0 / weights.sum()`` is bit-identical for
sphere and correct for any 3D quadrature.

Pre-fix this test produces residual ~6 (cylinder n_mu=2, n_phi=2,
n_cells=40) GROWING to ~12 at n_cells=160.  Post-fix the residual is
at FP-noise (~ O(nx · ULP)) at every mesh.

References
----------

* Hébert §3.9.4 Eqs. (3.432)-(3.435) — Carlson coupled-pole seed.
* Pomraning 1989 NSE 102:317-336 — structural-singularity isotropy.
* ERR-048 manifestation #3 entry in
  ``.claude/skills/vv-principles/error_catalog.md``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.sn.geometry import SNMesh
from orpheus.sn.operator import (
    build_equation_map_cylindrical,
    transport_operator_matvec_unified,
)
from orpheus.numerics.quadrature import Quadrature
from tests.sn._test_helpers import legacy_proxy_matvec, placeholder_materials


@pytest.mark.l0
@pytest.mark.verifies("hebert-3-432")
@pytest.mark.catches("ERR-048")
@pytest.mark.parametrize("n_cells", [10, 20, 40])
@pytest.mark.parametrize("n_mu", [2, 4])
@pytest.mark.parametrize("n_phi", [2, 4])
def test_cylinder_apply_matvec_preserves_flat_psi(
    n_cells: int, n_mu: int, n_phi: int,
) -> None:
    r"""``L · ψ_flat = Σ_t · ψ_flat`` per ordinate per cell.

    The cylindrical apply-matvec must satisfy the flat-flux invariant
    at machine precision: streaming and angular redistribution telescope
    to zero, leaving the collision term ``Σ_t · ψ``.  Pre-fix the
    Carlson seed ``Q_bar`` carried a sphere-only ``0.5 = 1/Σw_sphere``
    normalisation that breaks the invariant on 3D quadratures.  Post-fix
    ``1/weights.sum()`` is geometry-general; this test pins it.
    """
    fuel = get_mixture("B", "1g")
    R = 2.0
    geom = StructuredGeometry(
        geometry="CYL",
        regions=(Region(mat_id=0, outer_thickness_cm=R),),
        bcs=(BC.reflective,),
    )
    mesh = Mesh1D.from_geometry(
        geom, region_meshes=(RegionMesh(n_cells=n_cells),),
    )
    quad = Quadrature.product(n_mu=n_mu, n_phi=n_phi)
    sn_mesh = SNMesh(mesh, quad, placeholder_materials())
    nx = n_cells
    ng = 1
    sig_t = np.full((ng, nx, 1), 2.0)  # (ng, nx, ny) — PR-INDEX-3
    psi_flat = 10.0 / quad.weights.sum()

    # PR-TYPED-6c Step 7: route through ``transport_operator_matvec_unified``
    # — the legacy ``transport_operator_matvec_cylindrical`` retired.
    psi_view = np.full((quad.N, ng, nx, 1), psi_flat)
    m_cell = legacy_proxy_matvec(psi_view, sn_mesh, sig_t)
    eq_map = build_equation_map_cylindrical(nx, quad, ng)
    lhs = m_cell[
        eq_map.ordinate, :, eq_map.ix, eq_map.iy,
    ].T.ravel(order="F")

    expected = 2.0 * psi_flat  # = Σ_t · ψ_flat
    residual = float(np.max(np.abs(lhs - expected)))
    assert residual < 1e-10, (
        f"Cylinder apply-matvec residual on flat ψ: {residual:.6e}, "
        f"expected < 1e-10.  Configuration: n_cells={n_cells}, n_mu={n_mu}, "
        f"n_phi={n_phi}, Σw={quad.weights.sum():.6f}, ψ_flat={psi_flat:.6f}. "
        f"Pre-fix the residual grows with mesh refinement (Signature 1 "
        f"fingerprint, ERR-048 manifestation #3)."
    )


@pytest.mark.l0
@pytest.mark.verifies("hebert-3-432")
@pytest.mark.catches("ERR-048")
@pytest.mark.parametrize("n_cells", [10, 20, 40])
@pytest.mark.parametrize("n_mu", [2, 4])
@pytest.mark.parametrize("n_phi", [2, 4])
def test_cylinder_three_way_standoff(
    n_cells: int, n_mu: int, n_phi: int,
) -> None:
    r"""Krylov ≡ analytical ≡ SI on homogeneous reflective cylinder.

    The 4-leg correctness standard from the user-stated correctness
    framing: Krylov vs analytical, SI vs analytical, Krylov vs SI, and
    refinement convergence.  This test pins the first three at machine
    precision on every cell-quadrature combination.

    Analytical: ``ψ_n = Q / (Σ_t · (1−c) · Σw)`` from the homogeneous
    reflective fixed-source mass balance + Pomraning isotropy.  For
    mixture B 1G (Σ_t=2, c=0.95) with Q=1: ``ψ_n = 1 / (0.1 · Σw)``.
    """
    from orpheus.sn.solver import solve_sn_fixed_source

    fuel = get_mixture("B", "1g")
    R = 2.0
    geom = StructuredGeometry(
        geometry="CYL",
        regions=(Region(mat_id=0, outer_thickness_cm=R),),
        bcs=(BC.reflective,),
    )
    mesh = Mesh1D.from_geometry(
        geom, region_meshes=(RegionMesh(n_cells=n_cells),),
    )
    quad = Quadrature.product(n_mu=n_mu, n_phi=n_phi)
    N = quad.N
    nx = mesh.N
    ng = 1
    # R-1 Step 4 A1 — ``external_source`` is per-ordinate density.
    # Iso scalar magnitude 1 ⇒ per-ord density ``1/sum_w``.
    sum_w = float(quad.weights.sum())
    Q = np.full((N, ng, nx, 1), 1.0 / sum_w)

    psi_ref = 1.0 / (0.1 * sum_w)

    res_si = solve_sn_fixed_source(
        materials={0: fuel}, mesh=mesh, quadrature=quad, external_source=Q,
        boundary_condition="reflective", inner_solver="source_iteration",
    )
    res_k = solve_sn_fixed_source(
        materials={0: fuel}, mesh=mesh, quadrature=quad, external_source=Q,
        boundary_condition="reflective", inner_solver="krylov",
    )

    psi_si = res_si.angular_flux.values[:, :, 0, :]
    psi_k = res_k.angular_flux.values[:, :, 0, :]

    err_si_vs_ref = float(np.max(np.abs(psi_si - psi_ref)))
    err_k_vs_ref = float(np.max(np.abs(psi_k - psi_ref)))
    err_si_vs_k = float(np.max(np.abs(psi_si - psi_k)))

    assert err_si_vs_ref < 1e-9, (
        f"SI ≠ analytical: max|ψ_SI − ψ_ref| = {err_si_vs_ref:.6e}. "
        f"Configuration: n_cells={n_cells}, n_mu={n_mu}, n_phi={n_phi}."
    )
    assert err_k_vs_ref < 1e-9, (
        f"Krylov ≠ analytical: max|ψ_K − ψ_ref| = {err_k_vs_ref:.6e}. "
        f"Configuration: n_cells={n_cells}, n_mu={n_mu}, n_phi={n_phi}."
    )
    assert err_si_vs_k < 1e-9, (
        f"SI ≠ Krylov: max|ψ_SI − ψ_K| = {err_si_vs_k:.6e}. "
        f"Configuration: n_cells={n_cells}, n_mu={n_mu}, n_phi={n_phi}."
    )
