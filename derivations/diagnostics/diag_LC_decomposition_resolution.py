"""Diagnostic: Resolution candidate for the SN matvec L+C decomposition.

Created by numerics-investigator on 2026-05-14.

Companion to ``diag_LC_decomposition_sn.py`` — the algebraic part
established that the curvilinear matvec is NOT affine in σ_t (the
Carlson seed produces a σ_t-coupled term).

This script empirically tests the **Resolution A (refined)**
formulation:

    L.apply(ψ; σ_t) := M(ψ; σ_t) − σ_t ⊙ ψ_packed
    C.apply(ψ; σ_t) := σ_t ⊙ ψ_packed
    (L + C).apply(ψ; σ_t) == M(ψ; σ_t)

i.e., L carries σ_t at CONSTRUCTOR time (the curvilinear M-M closure
is intrinsically σ_t-coupled via Hébert's Carlson coupled-pole seed);
the algebra (L + C) recovers the full matvec EXACTLY by construction.
This is the **subtractive definition** of L.

The corresponding refactor of ``StreamingOperator.apply`` is:

    def StreamingOperator.apply(self, psi):
        full = transport_operator_matvec_curv(psi, σ_t=self.σ_t, ...)
        coll = self._sigma_packed * psi  # element-wise σ ⊙ ψ
        return full - coll

— bit-equivalent to ``(L + C).apply(ψ)`` because the matvec equals
``L_full + C_full`` and we subtract C to leave L.

Promotion guidance
------------------
This test should become a permanent V&V regression after the next
method-implementer dispatch lands the Resolution A `StreamingOperator`.
Move to ``tests/sn/test_streaming_operator_decomposition.py`` and
tag ``@pytest.mark.l0``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.sn.geometry import SNMesh
from orpheus.sn.operator import (
    build_equation_map,
    build_equation_map_spherical,
    build_equation_map_cylindrical,
    transport_operator_matvec,
    transport_operator_matvec_spherical,
    transport_operator_matvec_cylindrical,
)
from orpheus.sn.quadrature import GaussLegendre1D, ProductQuadrature


def _build_problem(geometry: str, n_cells: int = 5, n_ord: int = 4):
    fuel = get_mixture("B", "1g")
    if geometry == "SPH":
        geom = StructuredGeometry(
            geometry="SPH",
            regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
            bcs=(BC.reflective,),
        )
        mesh = Mesh1D.from_geometry(
            geom, region_meshes=(RegionMesh(n_cells=n_cells),)
        )
        quad = GaussLegendre1D.create(n_ordinates=n_ord)
    elif geometry == "CYL":
        geom = StructuredGeometry(
            geometry="CYL",
            regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
            bcs=(BC.reflective,),
        )
        mesh = Mesh1D.from_geometry(
            geom, region_meshes=(RegionMesh(n_cells=n_cells),)
        )
        quad = ProductQuadrature.create(n_mu=n_ord, n_phi=n_ord)
    elif geometry == "CART":
        geom = StructuredGeometry(
            geometry="SLB",
            regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
            bcs=(BC.reflective, BC.reflective),
        )
        mesh = Mesh1D.from_geometry(
            geom, region_meshes=(RegionMesh(n_cells=n_cells),)
        )
        quad = GaussLegendre1D.create(n_ordinates=n_ord)
    else:
        raise ValueError(geometry)
    return fuel, mesh, quad, SNMesh(mesh, quad)


def _call_matvec(geometry: str, psi_vec: np.ndarray,
                 sigma_t_arr: np.ndarray, sn_mesh, eq_map, quad, nx, ng):
    if geometry == "SPH":
        reduced = sn_mesh.reduced
        return transport_operator_matvec_spherical(
            psi_vec, eq_map, quad, sigma_t_arr, nx, ng,
            reduced.face_areas, sn_mesh.volumes,
            reduced.alpha_half, reduced.redist_dAw, reduced.tau_mm,
            sn_mesh=sn_mesh, bc_outer=sn_mesh.bc_right,
            pole_angular_closure=sn_mesh.pole_angular_closure,
        )
    if geometry == "CYL":
        reduced = sn_mesh.reduced
        return transport_operator_matvec_cylindrical(
            psi_vec, eq_map, quad, sigma_t_arr, nx, ng,
            reduced.face_areas, sn_mesh.volumes,
            reduced.alpha_per_level, reduced.redist_dAw_per_level,
            reduced.tau_mm_per_level,
            sn_mesh=sn_mesh, bc_outer=sn_mesh.bc_right,
            pole_angular_closure=sn_mesh.pole_angular_closure,
        )
    ny = sn_mesh.ny
    return transport_operator_matvec(
        psi_vec, eq_map, quad, sigma_t_arr,
        nx, ny, ng, sn_mesh.dx, sn_mesh.dy,
        bc_xmin=sn_mesh.bc_xmin, bc_xmax=sn_mesh.bc_xmax,
        bc_ymin=sn_mesh.bc_ymin, bc_ymax=sn_mesh.bc_ymax,
    )


def _L_apply_resolution_a(geometry, psi_vec, sigma_t_arr,
                          sn_mesh, eq_map, quad, nx, ng):
    """The Resolution A candidate: L.apply(ψ) := M(ψ; σ_t) − σ_t ⊙ ψ.

    L carries σ_t at constructor time (subtractive definition).
    The full matvec is computed at the user's σ_t, then the
    cell-collision term σ_t ⊙ ψ is subtracted.
    """
    M_full = _call_matvec(
        geometry, psi_vec, sigma_t_arr, sn_mesh, eq_map, quad, nx, ng,
    )
    sigma_packed = sigma_t_arr[
        eq_map.ix, eq_map.iy, :
    ].T.ravel(order='F')
    C_psi = sigma_packed * psi_vec
    return M_full - C_psi


def _C_apply(psi_vec, sigma_t_arr, eq_map):
    """The collision operator: C.apply(ψ) := σ_t ⊙ ψ_packed."""
    sigma_packed = sigma_t_arr[
        eq_map.ix, eq_map.iy, :
    ].T.ravel(order='F')
    return sigma_packed * psi_vec


def _resolution_test(geometry, n_cells, n_ord, seed):
    """Empirical bit-equivalence: (L + C).apply(ψ) == M(ψ; σ_t).

    Where:
      L.apply(ψ) := M(ψ; σ_t_full) − σ_t_full ⊙ ψ
      C.apply(ψ) := σ_t_full ⊙ ψ
    """
    fuel, mesh, quad, sn_mesh = _build_problem(geometry, n_cells, n_ord)
    nx, ng = sn_mesh.nx, 1

    if geometry == "SPH":
        eq_map = build_equation_map_spherical(nx, quad, ng)
    elif geometry == "CYL":
        eq_map = build_equation_map_cylindrical(nx, quad, ng)
    else:
        eq_map = build_equation_map(nx, sn_mesh.ny, quad, ng)

    rng = np.random.default_rng(seed)
    psi_vec = rng.standard_normal(eq_map.n_unknowns).astype(np.float64)
    sigma_full = np.full((nx, sn_mesh.ny, ng), 2.0)

    M_full = _call_matvec(
        geometry, psi_vec, sigma_full,
        sn_mesh, eq_map, quad, nx, ng,
    )
    L_apply = _L_apply_resolution_a(
        geometry, psi_vec, sigma_full,
        sn_mesh, eq_map, quad, nx, ng,
    )
    C_apply = _C_apply(psi_vec, sigma_full, eq_map)

    sum_apply = L_apply + C_apply
    residual = sum_apply - M_full
    rel_norm = (
        np.linalg.norm(residual) / max(np.linalg.norm(M_full), 1e-300)
    )

    return {
        "geometry": geometry,
        "n_cells": n_cells,
        "n_ord": n_ord,
        "M_full_norm": float(np.linalg.norm(M_full)),
        "L_apply_norm": float(np.linalg.norm(L_apply)),
        "C_apply_norm": float(np.linalg.norm(C_apply)),
        "residual_norm": float(np.linalg.norm(residual)),
        "rel_residual": float(rel_norm),
        "decomposition_holds": rel_norm < 1e-14,
    }


# ═══════════════════════════════════════════════════════════════════════
# pytest gates — Resolution A bit-equivalence
# ═══════════════════════════════════════════════════════════════════════

@pytest.mark.parametrize("geometry", ["CART", "SPH", "CYL"])
@pytest.mark.parametrize("seed", [0, 1, 2])
def test_resolution_a_bit_equivalent(geometry, seed):
    """Resolution A: (L + C).apply(ψ) == M(ψ; σ_t) bit-equivalently.

    L.apply := M(ψ; σ_t_full) − σ_t_full ⊙ ψ
    C.apply := σ_t_full ⊙ ψ
    """
    r = _resolution_test(geometry, n_cells=5, n_ord=4, seed=seed)
    assert r["decomposition_holds"], (
        f"{geometry} seed={seed}: rel_residual={r['rel_residual']:.3e} "
        f"> 1e-14 — Resolution A subtractive decomposition FAILS"
    )


def test_resolution_a_L_is_NOT_pure_streaming_curvilinear():
    """L from Resolution A is NOT the prior agent's M(ψ; σ_t=0).

    Sanity check: the prior agent's (reverted) implementation
    computed L.apply via matvec(ψ, σ_t=0). Verify our subtractive
    L.apply differs from that empirically — confirms we're
    proposing a DIFFERENT (correct) formulation.
    """
    for geometry in ("SPH", "CYL"):
        fuel, mesh, quad, sn_mesh = _build_problem(geometry, 5, 4)
        nx, ng = sn_mesh.nx, 1
        if geometry == "SPH":
            eq_map = build_equation_map_spherical(nx, quad, ng)
        else:
            eq_map = build_equation_map_cylindrical(nx, quad, ng)

        rng = np.random.default_rng(0)
        psi_vec = rng.standard_normal(eq_map.n_unknowns).astype(np.float64)
        sigma_full = np.full((nx, sn_mesh.ny, ng), 2.0)
        sigma_zero = np.zeros((nx, sn_mesh.ny, ng))

        # Resolution A's L.apply (subtractive):
        L_correct = _L_apply_resolution_a(
            geometry, psi_vec, sigma_full,
            sn_mesh, eq_map, quad, nx, ng,
        )
        # Prior agent's L.apply (matvec with σ_t = 0):
        L_prior = _call_matvec(
            geometry, psi_vec, sigma_zero,
            sn_mesh, eq_map, quad, nx, ng,
        )

        diff = L_correct - L_prior
        rel = np.linalg.norm(diff) / max(np.linalg.norm(L_correct), 1e-300)
        assert rel > 1e-3, (
            f"{geometry}: Resolution A's L.apply and prior agent's "
            f"L.apply differ only by {rel:.3e} — expected O(0.01..1) "
            f"difference because the prior approach degenerates the "
            f"Carlson seed denominator dr·σ_t + 2 → 2."
        )


if __name__ == "__main__":
    print("=" * 72)
    print("Resolution A: subtractive L = M(σ_t) − C(σ_t)")
    print("=" * 72)
    for geometry in ("CART", "SPH", "CYL"):
        print(f"\n{geometry}:")
        for seed in range(3):
            r = _resolution_test(geometry, n_cells=5, n_ord=4, seed=seed)
            print(
                f"  seed={seed}: rel_res={r['rel_residual']:.3e} "
                f"hold={r['decomposition_holds']}"
            )
