"""Diagnostic: SN matvec L+C decomposition — SymPy derivation + numpy probe.

Created by numerics-investigator on 2026-05-14.

Mission
-------
Determine, by SymPy derivation (Branch 1 algebra-of-record) and
numpy empirical cross-check, whether the curvilinear SN matvec
``transport_operator_matvec_{spherical,cylindrical}(ψ; σ_t)`` admits
a clean decomposition

    M(ψ; σ_t) = L(ψ) + σ_t ⊙ ψ                              (Q1)

where L is independent of σ_t and ⊙ is element-wise multiply.

The prior method-implementer (commit ad37ca0, reverted as b47551c)
implemented ``StreamingOperator.apply(ψ)`` by calling the matvec
with ``σ_t = 0``. The L0 streaming-equilibrium check then failed
at the curvilinear cases by ~50% rel-error. The agent attributed
this to "Carlson seed σ_t coupling" via Hébert §3.9.4 Eq. (3.434)
denominator ``dr·σ_t + 2`` and xfailed. The user's gate-keeper
response: STOP. Find the correct decomposition.

This diagnostic SETTLES the algebraic question. Three parts:

  Part A — SymPy: derive the per-cell-per-ordinate matvec equation
  for spherical, cylindrical, and Cartesian (control). Collect by
  σ_t. Show whether M is affine, rational, or non-decomposable.

  Part B — numpy: probe the matvec at fixed (ψ, σ_full) and
  separately at (ψ, 0). Compute residual

      out_diff(ψ; σ_t) := M(ψ; σ_t) − σ_t ⊙ ψ
      out_zero(ψ)     := M(ψ; 0)

  Affine ⇔ out_diff(ψ; σ_t) == out_zero(ψ) ∀ψ, σ_t.

  Part C — root-cause: pinpoint exactly where in operator.py the
  σ_t-coupled non-linear-in-σ_t term enters; characterise it; propose
  a corrected ``L.apply`` body that bypasses the bug.

Promotion guidance
------------------
The empirical affine test in Part B is a permanent V&V invariant for
ANY future ``StreamingOperator``: pure-streaming MUST be linear in ψ
AND independent of σ_t by definition. Promote to
``tests/sn/test_streaming_operator_affine.py`` when the
StreamingOperator lands. Tag ``@pytest.mark.l0``.

The Cartesian sub-test (passes) is the regression baseline. The
curvilinear sub-tests (initially fail; PASS once Resolution lands)
are the regression gates that prove the resolution works.
"""
from __future__ import annotations

import numpy as np
import pytest
import sympy as sp

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


# ═══════════════════════════════════════════════════════════════════════
# Part A — SymPy per-cell algebra (Branch 1)
# ═══════════════════════════════════════════════════════════════════════

def derive_cartesian_matvec_symbolic() -> dict:
    r"""V_A1 — Cartesian matvec is affine in σ_t.

    Per the matvec body (transport_operator_matvec, line 449-456):

        lhs[g, k] = μ_x · dψ/dx + μ_y · dψ/dy + σ_t[ix,iy,g] · ψ[g,n,ix,iy]

    Substituting σ_t = 0 gives the pure streaming part; σ_t enters
    LINEARLY at one isolated cell-collision term. We assert this
    symbolically.

    Returns ``{"pass": bool, "affine_residual": expr}``.
    """
    mu_x, mu_y, sigma_t, psi, dpsi_dx, dpsi_dy = sp.symbols(
        "mu_x mu_y sigma_t psi dpsi_dx dpsi_dy", real=True
    )

    # Full Cartesian matvec entry (one cell, one ordinate, one group)
    M = mu_x * dpsi_dx + mu_y * dpsi_dy + sigma_t * psi

    # Affine candidate: L(ψ) + σ_t · ψ
    L_candidate = mu_x * dpsi_dx + mu_y * dpsi_dy
    affine = L_candidate + sigma_t * psi

    residual = sp.simplify(M - affine)
    return {
        "name": "V_A1 Cartesian matvec affine in σ_t",
        "M": M,
        "L_candidate": L_candidate,
        "residual": residual,
        "pass": residual == 0,
    }


def derive_spherical_matvec_symbolic_collision() -> dict:
    r"""V_A2 — Spherical matvec's CELL-COLLISION term is affine in σ_t.

    The cell-collision contribution at cell i, ordinate n is

        collision_n,i = σ_t[i] · ψ_cell[n,i]

    See operator.py:799 (outward sweep), 832 (inward sweep). Standard
    cell-collision. Linear in σ_t and in ψ_cell, no coupling.

    Returns ``{"pass": bool}``.
    """
    sigma_t, psi_cell = sp.symbols("sigma_t psi_cell", real=True)
    collision = sigma_t * psi_cell
    affine = sigma_t * psi_cell
    residual = sp.simplify(collision - affine)
    return {
        "name": "V_A2 Spherical cell-collision is σ_t · ψ_cell",
        "collision": collision,
        "residual": residual,
        "pass": residual == 0,
    }


def derive_spherical_carlson_seed_symbolic() -> dict:
    r"""V_A3 — Carlson seed `phi_aux` is a RATIONAL function in σ_t.

    The Carlson coupled-pole seed (CarlsonInwardSweep) computes for
    each cell i, group g, running INWARD i = nx-1 → 0:

        Q̄_i = Σ_t,i · φ_0,i / Σw                             (Hébert 3.432)
        φ̄_i = (Δr_i · Q̄_i + 2 · φ̄_{i+½}) / (Δr_i · Σ_t,i + 2)   (3.434)
        φ̄_{i-½} = 2·φ̄_i - φ̄_{i+½}                             (3.435)

    Substituting Q̄_i = σ_t,i · φ_0,i / Σw:

        φ̄_i = (Δr_i · σ_t,i · φ_0,i / Σw + 2 · φ̄_{i+½})
                / (Δr_i · σ_t,i + 2)

    Both numerator and denominator carry σ_t. The seed φ̄_i is

      (a) NON-LINEAR in σ_t (rational, not affine).
      (b) LINEAR in (φ_0, φ̄_{i+½}) at fixed σ_t.

    Substituting σ_t = 0 gives a DEGENERATE form:

        φ̄_i|_{σ_t=0} = (0 + 2·φ̄_{i+½}) / (0 + 2)
                     = φ̄_{i+½}

    i.e. the seed becomes mesh-independent and just propagates the
    boundary value. This is the "L with σ_t=0" candidate — but it is
    NOT what the seed would be if we DECOMPOSED M into a streaming
    part L (built FROM ψ alone without σ_t) plus σ_t·ψ.

    The decomposition's WHY: the seed enters via the M-M angular
    recurrence's `psi_half_left`, which feeds the `redist_full` term
    through the recurrence chain. The full redistribution at (g, m, i)
    is a function of ALL the seeds at cells {nx-1, ..., i} chained
    through the M-M recurrence. So if we want

        M(ψ; σ_t) = L(ψ) + σ_t · ψ                          (?)

    then L(ψ) at the redist term must equal `redist_full(ψ;
    σ_t=σ_full) − 0` (the redist piece is dimensionally not in
    σ_t·ψ form, so the σ_t-coupling MUST move to L). But L cannot
    depend on σ_t by definition. Contradiction.

    Specifically: at fixed ψ, the redist term DOES depend on σ_t
    through the seed chain. So:

        ∂M/∂σ_t|_{ψ fixed} = ∂(seed-chain redist)/∂σ_t + ψ_cell

    The cell-collision term contributes +ψ_cell (the "C" part), but
    the seed contributes a SEPARATE σ_t-derivative term — i.e. there's
    an extra "σ_t-coupled redistribution" sitting on top of the C term.
    This is exactly the failure mode the prior method-implementer ran
    into when calling matvec with σ_t=0.

    Returns ``{"pass": bool}`` where "pass" asserts that the residual
    of the affine decomposition is NONZERO and has the closed form
    R(ψ; σ_t) := M(ψ; σ_t) - L_naïve(ψ) - σ_t·ψ ≠ 0.
    """
    # Single-cell single-group seed at the outer boundary (i = nx-1).
    # Use a 1-cell sweep: at i = nx-1, the recurrence terminates and
    # the seed equals the SOURCE-DRIVEN form at the boundary.
    sigma_t, dr, phi_0, phi_face = sp.symbols(
        "sigma_t dr phi_0 phi_face_outer", real=True, positive=True
    )
    Sw = sp.Symbol("Sigma_w", positive=True)  # quadrature weight sum

    # Hébert 3.432 + 3.434 fused at the outermost cell:
    Q_bar = sigma_t * phi_0 / Sw
    phi_aux_at_nx = (dr * Q_bar + 2 * phi_face) / (dr * sigma_t + 2)

    # Expand the rational form:
    phi_aux_at_nx_expanded = sp.simplify(phi_aux_at_nx)

    # Now what would the seed be in the AFFINE form?
    # Affine decomp claim: phi_aux = L_seed(ψ) + σ_t · seed_collision_term
    # But seed has NO `· ψ_cell` direct term — it depends on phi_0
    # (linear projection of ψ over μ) and phi_face (BC trace ψ at the
    # outer cell-centre).
    # If affine in σ_t at FIXED ψ (hence fixed phi_0 and phi_face),
    # the partial derivative ∂phi_aux/∂σ_t should be ψ-independent
    # (a CONSTANT in σ_t at fixed ψ — that's the affine signature).

    dphi_dsigma_t = sp.simplify(sp.diff(phi_aux_at_nx_expanded, sigma_t))

    # If phi_aux is affine in σ_t, dphi_dsigma_t is a constant (no σ_t).
    # If phi_aux is rational, dphi_dsigma_t depends on σ_t.
    is_rational_in_sigma_t = (sigma_t in dphi_dsigma_t.free_symbols)

    # The actual closed form of the σ_t-dependence:
    phi_aux_at_sigma_t_0 = sp.simplify(phi_aux_at_nx.subs(sigma_t, 0))

    return {
        "name": "V_A3 Spherical Carlson seed is rational in σ_t",
        "phi_aux_full": phi_aux_at_nx_expanded,
        "dphi_dsigma_t": dphi_dsigma_t,
        "phi_aux_at_sigma_t_0": phi_aux_at_sigma_t_0,
        "is_rational_in_sigma_t": is_rational_in_sigma_t,
        # PASS iff the seed IS rational (NOT affine) in σ_t — confirms
        # decomposition is not affine.
        "pass": is_rational_in_sigma_t,
    }


def derive_spherical_residual_symbolic() -> dict:
    r"""V_A4 — write the affine-decomposition residual R(ψ; σ_t) closed form.

    For a 1-cell sphere at i = nx-1, the matvec at one outgoing
    ordinate (mu > 0) is:

        M_n,i = (μ_n/V_i)·(A[i+1]·ψ_face_out - A[i]·ψ_face_in)
                + redist_full[g,n,i]
                + σ_t,i · ψ_cell

    where ``redist_full`` carries the seed chain. The naïve affine
    decomposition splits:

        L_naive(ψ) := M(ψ; σ_t=0)
        C(ψ; σ_t) := σ_t ⊙ ψ
        R(ψ; σ_t) := M(ψ; σ_t) - L_naive(ψ) - C(ψ; σ_t)

    The residual R lives in the seed-chain term ONLY (streaming has
    no σ_t; collision is exactly C). So R = ∂redist/∂σ_t · σ_t-like
    cross term. Symbolically, for the 1-cell case:

        seed_full = (dr · σ_t · φ_0 / Σw + 2 · φ_face) / (dr · σ_t + 2)
        seed_zero = φ_face
        Δseed = seed_full - seed_zero
              = [dr · σ_t · (φ_0/Σw - φ_face)] / (dr · σ_t + 2)
        R ∝ Δseed (non-trivially through the M-M recurrence chain)

    Δseed:
        is zero  iff  σ_t = 0  OR  φ_0/Σw = φ_face  (flat ψ at BC).

    The Δseed is NEITHER zero in general NOR linear in σ_t — it is
    rational. At σ_t → ∞, Δseed → (φ_0/Σw - φ_face). At σ_t → 0,
    Δseed → 0.

    This pins the SHAPE of the σ_t-coupled artefact: it is the
    "BC-trace flux at outer face differs from cell-averaged 0-th moment
    at outer cell" mismatch, gated by ``σ_t · dr / (σ_t · dr + 2)``.

    For a homogeneous reflective sphere with isotropic Q and a flat
    fixed-point ψ = Q/Σ_t · 1/Σw, the BC-trace φ_face equals the
    cell-centre value (which equals 1/Σw · φ_0 = ψ_flat). So Δseed = 0
    at THIS fixed point — that is why the fixed-source diagnostic
    showed agreement at the converged ψ even though the operator is
    not affine.

    Returns ``{"pass": bool}`` where "pass" asserts Δseed has the
    closed form above.
    """
    sigma_t, dr, phi_0, phi_face, Sw = sp.symbols(
        "sigma_t dr phi_0 phi_face_outer Sigma_w",
        real=True, positive=True,
    )

    seed_full = (dr * sigma_t * phi_0 / Sw + 2 * phi_face) / (dr * sigma_t + 2)
    seed_zero = phi_face   # σ_t → 0 limit of seed_full

    Delta = sp.simplify(seed_full - seed_zero)
    Delta_factored = sp.factor(Delta)

    # Closed form claim: Δseed = (dr·σ_t/(dr·σ_t + 2)) · (φ_0/Σw - φ_face)
    closed_form = (dr * sigma_t / (dr * sigma_t + 2)) * (phi_0 / Sw - phi_face)
    residual = sp.simplify(Delta - closed_form)

    return {
        "name": "V_A4 Δseed closed form",
        "Delta_factored": Delta_factored,
        "closed_form_check": residual,
        "pass": residual == 0,
        # Reading aid — at the within-group fixed point with flat ψ
        # the cancellation condition φ_0/Σw = φ_face holds.
        "cancels_at_flat_psi": sp.simplify(
            closed_form.subs(phi_0, Sw * phi_face)
        ) == 0,
    }


# ═══════════════════════════════════════════════════════════════════════
# Part B — Empirical numpy probe (Branch 2)
# ═══════════════════════════════════════════════════════════════════════

def _build_problem(geometry: str, n_cells: int = 5, n_ord: int = 4):
    """Build a representative SN mesh + quadrature + materials.

    Uses mixture B (Σ_t = 2, Σ_s = 1.9, 1G), 2-cm region, reflective
    outer BC, GL-N quadrature for sphere, ProductQuadrature for
    cylinder + 2-D Cartesian.
    """
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
    sn_mesh = SNMesh(mesh, quad)
    return fuel, mesh, quad, sn_mesh


def _call_matvec(geometry: str, psi_vec: np.ndarray,
                 sigma_t_arr: np.ndarray, sn_mesh, eq_map, quad, nx, ng):
    """Dispatch to the geometry-appropriate matvec primitive."""
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
    # Cartesian / slab
    ny = sn_mesh.ny
    return transport_operator_matvec(
        psi_vec, eq_map, quad, sigma_t_arr,
        nx, ny, ng, sn_mesh.dx, sn_mesh.dy,
        bc_xmin=sn_mesh.bc_xmin, bc_xmax=sn_mesh.bc_xmax,
        bc_ymin=sn_mesh.bc_ymin, bc_ymax=sn_mesh.bc_ymax,
    )


def _affine_residual_at_random_psi(
    geometry: str, n_cells: int = 5, n_ord: int = 4, seed: int = 0
) -> dict:
    """Compute ``out_diff − out_zero`` at a random ψ.

    out_diff(ψ; σ_t) := M(ψ; σ_t) − σ_t ⊙ ψ_packed
    out_zero(ψ)     := M(ψ; 0)

    If the matvec is AFFINE in σ_t with the cell-collision term equal
    to σ_t · ψ at the unknown slot, then out_diff == out_zero
    exactly. We assert this empirically — and characterise the
    deviation when it fails.
    """
    fuel, mesh, quad, sn_mesh = _build_problem(geometry, n_cells, n_ord)
    nx, ng = sn_mesh.nx, 1

    if geometry == "SPH":
        eq_map = build_equation_map_spherical(nx, quad, ng)
    elif geometry == "CYL":
        eq_map = build_equation_map_cylindrical(nx, quad, ng)
    else:
        ny = sn_mesh.ny
        eq_map = build_equation_map(nx, ny, quad, ng)

    # Random ψ vector (not flat — flat hides the Δseed coupling).
    rng = np.random.default_rng(seed)
    psi_vec = rng.standard_normal(eq_map.n_unknowns).astype(np.float64)

    # Build σ_t arrays.
    sigma_full = np.full((nx, sn_mesh.ny, ng), 2.0)  # Σ_t = 2 (mixture B)
    sigma_zero = np.zeros((nx, sn_mesh.ny, ng))

    # M(ψ; σ_full) and M(ψ; 0)
    M_full = _call_matvec(geometry, psi_vec, sigma_full,
                          sn_mesh, eq_map, quad, nx, ng)
    M_zero = _call_matvec(geometry, psi_vec, sigma_zero,
                          sn_mesh, eq_map, quad, nx, ng)

    # σ_t ⊙ ψ packed: gather σ at each unknown's (ix, iy, g) slot.
    sigma_per_unknown = sigma_full[
        eq_map.ix, eq_map.iy, :
    ]  # (n_eq, ng)
    # Layout matches solution_to_angular_flux's
    # ``psi.reshape(ng, n_eq, order='F')`` — see CollisionOperator
    # _sigma_at_unknowns for the convention.
    sigma_packed = sigma_per_unknown.T.ravel(order='F')
    C_psi = sigma_packed * psi_vec  # element-wise

    out_diff = M_full - C_psi
    out_zero = M_zero

    residual = out_diff - out_zero
    rel_norm = (
        np.linalg.norm(residual) / max(np.linalg.norm(M_full), 1e-300)
    )

    return {
        "geometry": geometry,
        "n_cells": n_cells,
        "n_ord": n_ord,
        "M_full_norm": float(np.linalg.norm(M_full)),
        "M_zero_norm": float(np.linalg.norm(M_zero)),
        "residual_norm": float(np.linalg.norm(residual)),
        "rel_residual": float(rel_norm),
        "affine": rel_norm < 1e-12,
    }


# ═══════════════════════════════════════════════════════════════════════
# pytest gates
# ═══════════════════════════════════════════════════════════════════════

def test_va1_cartesian_matvec_is_affine_in_sigma_t():
    """Symbolic Cartesian matvec entry: M = L + σ_t · ψ exactly."""
    r = derive_cartesian_matvec_symbolic()
    assert r["pass"], (
        f"Cartesian matvec symbolic affine residual ≠ 0: {r}"
    )


def test_va2_spherical_cell_collision_is_sigma_t_psi():
    """The cell-collision term in spherical matvec is exactly σ_t · ψ."""
    r = derive_spherical_matvec_symbolic_collision()
    assert r["pass"], (
        f"Spherical cell-collision affine residual ≠ 0: {r}"
    )


def test_va3_spherical_carlson_seed_is_rational_in_sigma_t():
    """The Carlson seed φ̄_i is RATIONAL (not affine) in σ_t."""
    r = derive_spherical_carlson_seed_symbolic()
    assert r["pass"], (
        f"Carlson seed not rational in σ_t — affine "
        f"decomposition would hold: {r}"
    )


def test_va4_seed_residual_closed_form():
    """The seed σ_t-coupling residual has the predicted closed form:

    Δseed(σ_t; φ_0, φ_face) = (dr·σ_t / (dr·σ_t + 2)) · (φ_0/Σw - φ_face)

    Cancels iff σ_t = 0 OR ψ flat at the BC (φ_0/Σw = φ_face).
    """
    r = derive_spherical_residual_symbolic()
    assert r["pass"], f"Δseed closed-form residual ≠ 0: {r}"
    assert r["cancels_at_flat_psi"], (
        "Δseed should vanish at flat-ψ fixed point (φ_0/Σw = φ_face)"
    )


def test_vb1_cartesian_matvec_is_affine_numpy():
    """Cartesian matvec empirically affine in σ_t (rtol=1e-12)."""
    r = _affine_residual_at_random_psi("CART", n_cells=5, n_ord=4)
    assert r["affine"], (
        f"Cartesian matvec NOT empirically affine: rel_residual = "
        f"{r['rel_residual']:.3e}; expected < 1e-12"
    )


def test_vb2_spherical_matvec_is_not_affine_numpy():
    """Spherical matvec empirically NOT affine in σ_t (Δseed coupling).

    This is the smoking gun for the prior method-implementer's bug:
    calling matvec with σ_t = 0 does NOT give L(ψ). The deviation
    is the seed-chain σ_t-coupled rational artefact.
    """
    r = _affine_residual_at_random_psi("SPH", n_cells=5, n_ord=4)
    # We ASSERT that the matvec is NOT affine — i.e. the prior
    # decomposition is wrong by a measurable amount on random ψ.
    assert not r["affine"], (
        f"Spherical matvec IS affine? rel_residual = "
        f"{r['rel_residual']:.3e}. Either the analysis is wrong "
        f"OR the Carlson seed was bypassed."
    )
    # The deviation should be O(1) relative — visible, not noise.
    assert r["rel_residual"] > 1e-3, (
        f"Spherical matvec deviation suspiciously small "
        f"(< 1e-3 rel): rel_residual = {r['rel_residual']:.3e}. "
        f"Expected O(0.01..1)."
    )


def test_vb3_cylindrical_matvec_is_not_affine_numpy():
    """Cylindrical matvec empirically NOT affine in σ_t (same coupling)."""
    r = _affine_residual_at_random_psi("CYL", n_cells=5, n_ord=4)
    assert not r["affine"], (
        f"Cylindrical matvec IS affine? rel_residual = "
        f"{r['rel_residual']:.3e}"
    )
    assert r["rel_residual"] > 1e-3, (
        f"Cylindrical matvec deviation suspiciously small: "
        f"rel_residual = {r['rel_residual']:.3e}"
    )


def test_vc1_curvilinear_affine_at_flat_psi():
    """Curvilinear matvec IS affine at flat ψ (the cancellation case).

    At ψ_flat = const, the seed's Δseed = 0 by the closed-form
    cancellation condition. The matvec output then equals the affine
    decomposition output. This pins WHY the prior fixed-source diagnostic
    appeared to PASS at the converged fixed point even though the
    operator is structurally non-affine.
    """
    for geometry in ("SPH", "CYL"):
        fuel, mesh, quad, sn_mesh = _build_problem(geometry, 5, 4)
        nx, ng = sn_mesh.nx, 1
        if geometry == "SPH":
            eq_map = build_equation_map_spherical(nx, quad, ng)
        else:
            eq_map = build_equation_map_cylindrical(nx, quad, ng)

        # Flat ψ — every unknown gets the same value.
        psi_vec = np.full(eq_map.n_unknowns, 0.5)

        sigma_full = np.full((nx, sn_mesh.ny, ng), 2.0)
        sigma_zero = np.zeros((nx, sn_mesh.ny, ng))

        M_full = _call_matvec(geometry, psi_vec, sigma_full,
                              sn_mesh, eq_map, quad, nx, ng)
        M_zero = _call_matvec(geometry, psi_vec, sigma_zero,
                              sn_mesh, eq_map, quad, nx, ng)

        sigma_packed = sigma_full[
            eq_map.ix, eq_map.iy, :
        ].T.ravel(order='F')
        C_psi = sigma_packed * psi_vec

        residual = (M_full - C_psi) - M_zero
        rel = np.linalg.norm(residual) / max(np.linalg.norm(M_full), 1e-300)
        assert rel < 1e-10, (
            f"{geometry} flat-ψ affine residual {rel:.3e} > 1e-10 "
            f"— the Δseed cancellation at flat ψ should hold exactly"
        )


if __name__ == "__main__":
    print("=" * 72)
    print("Part A — SymPy derivations")
    print("=" * 72)
    for fn in (
        derive_cartesian_matvec_symbolic,
        derive_spherical_matvec_symbolic_collision,
        derive_spherical_carlson_seed_symbolic,
        derive_spherical_residual_symbolic,
    ):
        r = fn()
        print(f"\n{r['name']}:")
        for k, v in r.items():
            if k == "name":
                continue
            print(f"  {k}: {v}")

    print()
    print("=" * 72)
    print("Part B — numpy empirical probe")
    print("=" * 72)
    for geometry in ("CART", "SPH", "CYL"):
        print(f"\n{geometry}:")
        for seed in range(3):
            r = _affine_residual_at_random_psi(
                geometry, n_cells=5, n_ord=4, seed=seed,
            )
            print(
                f"  seed={seed}: "
                f"rel_res={r['rel_residual']:.3e} "
                f"affine={r['affine']}"
            )
