r"""Phase-3C-2 annulus (hollow cylinder) **Variant α** Green's function
reference (homogeneous shell :math:`r \in [R_{\rm in}, R_{\rm out}]`,
**independent** specular reflectivities :math:`\alpha_{\rm in} \in
[0, 1]` on the inner cavity surface and :math:`\alpha_{\rm out} \in
[0, 1]` on the outer surface) — research-grade prototype.

Standalone module, parallel to :mod:`.greens_function_cylinder` (solid
cylinder) and :mod:`.greens_function_hollow_sphere` (hollow sphere).
This is the **last 2-BC topology** in the Variant α plan; once it
ships, the 6-geometry × 2-topology family (sphere, cylinder, slab,
slab-asym, hollow sphere, annulus) is complete on the unified
rank-1 / rank-2 framework via :mod:`.variant_alpha_core`.

Architecture
------------

The annulus is the cylindrical analog of the hollow sphere with the
**same closure algebra** (rank-2 BIE block resolvent on through-rays
+ rank-1 outer-only resolvent — both byte-equal-shared with hollow
sphere via :mod:`.variant_alpha_core`) and a different chord meaning:

- **Phase-space**: cylinder uses :math:`(r, \mu_{\rm axial},
  \varphi_{\rm az})` instead of sphere's :math:`(r, \mu)`. Conserved
  quantities under specular reflection at radial-normal surfaces:
  :math:`\mu_{\rm axial}` (axial cosine) and :math:`b = r\,|\sin\varphi_{
  \rm az}|` (in-plane impact parameter).
- **3D chord scaling** (Issue #129 angle-resolved discipline): the
  2D in-plane chord is lifted to 3D arclength by the axial-correction
  factor :math:`1/\sqrt{1 - \mu_{\rm axial}^2}`. The cylinder Bickley-
  Naylor :math:`\mathrm{Ki}_n` form pre-integrates axial and produces
  ~22 % planar-limit mismatch — Variant α stays angle-resolved to
  avoid this.

Phase-space partition by impact parameter
------------------------------------------

For each phase-space point :math:`(r, \mu_{\rm axial}, \varphi_{\rm
az})` with :math:`r \in [R_{\rm in}, R_{\rm out}]`, the **conserved
in-plane impact parameter** is :math:`b = r\,|\sin\varphi_{\rm az}|`.
Two cases:

- :math:`b > R_{\rm in}` — outer-only ray. The 2D in-plane projection
  ray does not intersect the inner-radius circle; the particle bounces
  between two points on the OUTER cylinder only. **Topologically
  rank-1**, structurally identical to a solid-cylinder ray at the
  same outer radius and impact parameter (with axial correction).
- :math:`b \le R_{\rm in}` — through-ray. The 2D in-plane ray crosses
  the inner cavity circle. Under interpretation (A) — the inner
  surface is a specular reflector with reflectivity :math:`\alpha_{\rm
  in}` — the particle bounces alternately inner ↔ outer.
  **Topologically rank-2**.

Trajectory case analysis (first-leg backward)
---------------------------------------------

For each :math:`(r, \mu_{\rm axial}, \varphi_{\rm az})`:

- **Outer-only** (:math:`b > R_{\rm in}`): backward hits OUTER first
  at 2D arclength :math:`L_{\rm 2D, first} = r\cos\varphi_{\rm az} +
  \sqrt{R_{\rm out}^2 - b^2}`. Surface = outer.
- **Through-ray with** :math:`\cos\varphi_{\rm az} > 0`: backward goes
  inward; hits INNER first at :math:`L_{\rm 2D, first} =
  r\cos\varphi_{\rm az} - \sqrt{R_{\rm in}^2 - b^2}` (non-negative
  since :math:`r \ge R_{\rm in}`). Surface = inner.
- **Through-ray with** :math:`\cos\varphi_{\rm az} \le 0`: backward
  goes outward; hits OUTER first at :math:`L_{\rm 2D, first} =
  r\cos\varphi_{\rm az} + \sqrt{R_{\rm out}^2 - b^2}` (non-negative).
  Surface = outer.

3D first-leg arclength: :math:`L_{\rm first}^{\rm 3D} = L_{\rm 2D,
first} / \sqrt{1 - \mu_{\rm axial}^2}`.

Outer-only branch (rank-1)
--------------------------

Rank-1 closure on bounce-period chord:

.. math::

   L_{\rm period}^{\rm 3D}(b, \mu_{\rm axial})
        = \frac{2\sqrt{R_{\rm out}^2 - b^2}}{\sqrt{1 - \mu_{\rm axial}^2}}.

Identical to solid cylinder. Position along antipodal chord at 2D
arclength :math:`s_{\rm 2D, p}`: :math:`r_{\rm chord}^2 = b^2 +
(s_{\rm 2D, p} - \sqrt{R_{\rm out}^2 - b^2})^2`.

Through-ray branch (rank-2)
---------------------------

Rank-2 closure with single-transit shell-traversal optical depth
(Eq. :eq:`peierls-greens-annulus-tau-step`):

.. math::

   \tau_{\rm step}^{\rm annulus}(b, \mu_{\rm axial})
        = \Sigma_t \cdot
          \frac{\sqrt{R_{\rm out}^2 - b^2} - \sqrt{R_{\rm in}^2 - b^2}}
               {\sqrt{1 - \mu_{\rm axial}^2}}.

Two single-transit B integrals:

- :math:`B_{\rm out}` (inner→outer shell chord): position along
  trailing shell segment of the outer-cylinder antipodal chord,
  :math:`r^2(u) = b^2 + (\sqrt{R_{\rm in}^2 - b^2} + u)^2` where
  :math:`u \in [0, L_{\rm step}^{\rm 2D}]` is the 2D arclength.
- :math:`B_{\rm in}` (outer→inner shell chord): position along
  leading shell segment (reversed), :math:`r^2(u) = b^2 +
  (\sqrt{R_{\rm out}^2 - b^2} - u)^2`.

Closure mapping (slab-asym → annulus): :math:`\alpha_L \mapsto \alpha_{
\rm in}`, :math:`\alpha_R \mapsto \alpha_{\rm out}`, :math:`B_{LR}
\mapsto B_{\rm out}`, :math:`B_{RL} \mapsto B_{\rm in}`,
:math:`\psi_L^+ \mapsto \psi_{\rm in}^{\rm out}`, :math:`\psi_R^-
\mapsto \psi_{\rm out}^{\rm in}`. Surface arrival → ``surface='left'``
for inner, ``'right'`` for outer.

Total angular flux:
:math:`\psi(r, \mu_{\rm axial}, \varphi_{\rm az}) = F + e^{-\Sigma_t
L_{\rm first}^{\rm 3D}}\,\psi_{\rm surface}`.

Scalar flux:
:math:`\phi(r) = \int_{-1}^{1}\!\mathrm d\mu_{\rm axial}
\int_0^{2\pi}\!\mathrm d\varphi_{\rm az}\,\psi`.

Assumptions for Phase 3C-2 prototype
-------------------------------------

- Homogeneous shell (single :math:`\Sigma_t, \Sigma_s, \nu\Sigma_f`);
  multi-region deferred (the 6-geometry family completes here on
  homogeneous; multi-region for hollow cylinder is a follow-on).
- Isotropic scattering.
- Independent reflective specular BCs :math:`\alpha_{\rm in},
  \alpha_{\rm out} \in [0, 1]`.

References
----------

- Sanchez, R. (1986). *Transp. Theor. Stat. Phys.* 14.
- :file:`/.claude/plans/peierls-greens-cylinder-and-2bc.md` — Phase
  3C-2 annulus plan.
- :mod:`orpheus.derivations.continuous.peierls_greens_function.origins.specular.greens_function_annulus`
  — V_α1_annulus / V_α2_annulus / V_α2_annulus.aux / V_α3_annulus
  SymPy verifications.
- :mod:`.greens_function_hollow_sphere` — Phase-3C-1 hollow sphere
  reference solver (the rank-2 + impact-parameter template lifted
  here to cylinder cross-section).
- :mod:`.greens_function_cylinder` — Phase-1 solid cylinder reference
  solver (the cylinder 3D angular framework + rank-1 outer-only
  template).
- Issue #129 — cylinder Bickley-Naylor planar-limit physics, motivating
  the angle-resolved (NOT axial-pre-integrated) discipline.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import CubicSpline

from orpheus.derivations.continuous.peierls_greens_function.variant_alpha_core import (
    apply_variant_alpha_closure,
    apply_variant_alpha_closure_rank2,
)


@dataclass(frozen=True)
class AnnulusGreensResult:
    """Result of Variant α power iteration on homogeneous annulus
    with INDEPENDENT specular reflectivities
    :math:`\\alpha_{\\rm in}, \\alpha_{\\rm out}`."""

    k_eff: float
    psi: np.ndarray  # (n_r, n_mu_axial, n_phi_az)
    phi: np.ndarray  # (n_r,)
    r_nodes: np.ndarray
    mu_axial_nodes: np.ndarray
    phi_az_nodes: np.ndarray
    iterations: int
    converged: bool


@dataclass(frozen=True)
class AnnulusGreensMGResult:
    """Result of multi-group Variant α power iteration on homogeneous
    annulus with INDEPENDENT specular BCs."""

    k_eff: float
    psi_g: np.ndarray  # (G, n_r, n_mu_axial, n_phi_az)
    phi_g: np.ndarray  # (G, n_r)
    r_nodes: np.ndarray
    mu_axial_nodes: np.ndarray
    phi_az_nodes: np.ndarray
    iterations: int
    converged: bool


# ═══════════════════════════════════════════════════════════════════════
# Per-group Variant α operator (annulus, rank-2 + impact-parameter
# partition with cylinder 3D angular phase-space)
# ═══════════════════════════════════════════════════════════════════════


def _apply_operator_annulus(
    source_profile: np.ndarray,
    r_nodes: np.ndarray,
    mu_axial_nodes: np.ndarray,
    phi_az_nodes: np.ndarray,
    R_in: float,
    R_out: float,
    sigma_t: float,
    alpha_in: float,
    alpha_out: float,
    *,
    n_traj_quad: int,
) -> np.ndarray:
    r"""Per-group annulus Variant α operator with rank-2 + rank-1
    closure on the cylinder 3D angular phase-space.

    For each :math:`(r_i, \mu_{\rm axial}, \varphi_{\rm az})`:

    1. Compute :math:`b = r_i\,|\sin\varphi_{\rm az}|`. Partition:
       outer-only (:math:`b > R_{\rm in}`) or through-ray (:math:`b
       \le R_{\rm in}`).

    2. Outer-only: identical to solid-cylinder algebra. Rank-1 closure
       with :math:`\alpha = \alpha_{\rm out}` and 3D-lifted bounce-
       period chord :math:`L_{\rm period}^{\rm 3D} = 2\sqrt{R_{\rm
       out}^2 - b^2} / \sqrt{1 - \mu_{\rm axial}^2}`.

    3. Through-ray: rank-2 closure with single-transit shell-traversal
       :math:`\tau_{\rm step} = \Sigma_t\,(\sqrt{R_{\rm out}^2 - b^2} -
       \sqrt{R_{\rm in}^2 - b^2}) / \sqrt{1 - \mu_{\rm axial}^2}`.
       First-leg backward goes to inner surface for :math:`\cos\varphi_{
       \rm az} > 0` and outer surface for :math:`\cos\varphi_{\rm az}
       \le 0`.

    Parameters
    ----------
    source_profile : (n_r,) ndarray
        :math:`q_g(r_i)/(4\pi)` — already-divided isotropic source
        per steradian per unit volume. Cubic-spline-interpolated for
        evaluation along trajectory points.
    r_nodes : (n_r,) ndarray
        Radial nodes on :math:`[R_{\rm in}, R_{\rm out}]`.
    mu_axial_nodes : (n_mu,) ndarray
        Axial-cosine nodes on :math:`[-1, 1]`.
    phi_az_nodes : (n_phi,) ndarray
        Azimuthal nodes on :math:`[0, 2\pi)`.
    R_in, R_out : float
        Inner cavity radius and outer cylinder radius.
    sigma_t : float
        Per-group total cross section in the shell.
    alpha_in, alpha_out : float
        Per-surface reflectivities in :math:`[0, 1]`.
    n_traj_quad : int
        Trajectory + bounce-chord quadrature order (in 2D arclength).

    Returns
    -------
    (n_r, n_mu, n_phi) ndarray
    """
    source_interp = CubicSpline(r_nodes, source_profile, extrapolate=True)

    s_quad_raw, w_quad_raw = np.polynomial.legendre.leggauss(n_traj_quad)
    s_unit = 0.5 * (s_quad_raw + 1.0)
    w_unit = 0.5 * w_quad_raw

    n_r = len(r_nodes)
    n_mu = len(mu_axial_nodes)
    n_phi = len(phi_az_nodes)
    psi_new = np.zeros((n_r, n_mu, n_phi))

    R_out_sq = R_out * R_out
    R_in_sq = R_in * R_in

    for i in range(n_r):
        r = r_nodes[i]
        for q_idx in range(n_mu):
            mu_axial = mu_axial_nodes[q_idx]
            # In-plane velocity fraction: sqrt(1 - μ_axial²).
            s_in_plane = np.sqrt(max(1.0 - mu_axial * mu_axial, 1e-300))
            inv_s_in_plane = 1.0 / s_in_plane

            for p_idx in range(n_phi):
                phi_az = phi_az_nodes[p_idx]
                cos_phi = np.cos(phi_az)
                sin_phi = np.sin(phi_az)

                # In-plane impact parameter b = r·|sin(φ)|.
                b = r * abs(sin_phi)
                b_sq = b * b

                # Outer-cylinder discriminant (always real).
                disc_out = R_out_sq - b_sq
                sqrt_disc_out = np.sqrt(max(disc_out, 0.0))

                if b > R_in:
                    # ───────── Outer-only branch (rank-1) ─────────
                    # First-leg 2D backward chord to outer cylinder.
                    L_2D_first = r * cos_phi + sqrt_disc_out
                    L_first_3D = L_2D_first * inv_s_in_plane

                    # First-leg integral parametrised by 2D arclength.
                    s_pts_2D = s_unit * L_2D_first
                    r_traj_sq = (
                        r * r - 2.0 * r * cos_phi * s_pts_2D
                        + s_pts_2D * s_pts_2D
                    )
                    r_traj = np.sqrt(np.clip(r_traj_sq, R_in_sq, R_out_sq))
                    tau_3D = sigma_t * s_pts_2D * inv_s_in_plane
                    integrand_F = source_interp(r_traj) * np.exp(-tau_3D)
                    F = L_first_3D * np.sum(w_unit * integrand_F)

                    if alpha_out == 0.0:
                        psi_new[i, q_idx, p_idx] = F
                        continue

                    # Bounce-period 2D in-plane chord (full antipodal).
                    L_2D_period = 2.0 * sqrt_disc_out
                    if L_2D_period <= 0.0:
                        psi_new[i, q_idx, p_idx] = F
                        continue
                    L_period_3D = L_2D_period * inv_s_in_plane

                    s_pts_2D_p = s_unit * L_2D_period
                    shifted = s_pts_2D_p - 0.5 * L_2D_period
                    r_chord_sq = b_sq + shifted * shifted
                    r_chord = np.sqrt(
                        np.clip(r_chord_sq, R_in_sq, R_out_sq)
                    )
                    tau_3D_p = sigma_t * s_pts_2D_p * inv_s_in_plane
                    integrand_B = source_interp(r_chord) * np.exp(-tau_3D_p)
                    B = L_period_3D * np.sum(w_unit * integrand_B)

                    psi_new[i, q_idx, p_idx] = apply_variant_alpha_closure(
                        F=F, B=B,
                        tau_first_leg=sigma_t * L_first_3D,
                        tau_period=sigma_t * L_period_3D,
                        alpha=alpha_out,
                    )
                    continue

                # ───────── Through-ray branch (rank-2) ─────────
                disc_in = R_in_sq - b_sq
                sqrt_disc_in = np.sqrt(max(disc_in, 0.0))

                # Single-transit 2D shell-traversal length:
                #   L_step_2D = sqrt(R_out² - b²) - sqrt(R_in² - b²).
                # 3D-lifted: L_step_3D = L_step_2D / s_in_plane.
                L_step_2D = sqrt_disc_out - sqrt_disc_in
                L_step_3D = L_step_2D * inv_s_in_plane

                # First-leg backward depending on sign of cos(φ).
                if cos_phi > 0.0:
                    # Backward goes INWARD; first arrival is INNER surface.
                    # 2D arclength: L = r·cos(φ) - sqrt(R_in² - b²).
                    L_2D_first = r * cos_phi - sqrt_disc_in
                    surface = "inner"
                else:
                    # cos_phi <= 0: backward goes OUTWARD;
                    # first arrival is OUTER surface.
                    # 2D arclength: L = r·cos(φ) + sqrt(R_out² - b²).
                    L_2D_first = r * cos_phi + sqrt_disc_out
                    surface = "outer"

                # Numerical guard for tangent rays.
                L_2D_first = max(L_2D_first, 0.0)
                L_first_3D = L_2D_first * inv_s_in_plane

                # First-leg trajectory points along the 2D backward
                # chord; r²(s) = r² - 2·r·s·cos(φ) + s².
                s_pts_2D = s_unit * L_2D_first
                r_traj_sq = (
                    r * r - 2.0 * r * cos_phi * s_pts_2D
                    + s_pts_2D * s_pts_2D
                )
                r_traj = np.sqrt(np.clip(r_traj_sq, R_in_sq, R_out_sq))
                tau_3D = sigma_t * s_pts_2D * inv_s_in_plane
                integrand_F = source_interp(r_traj) * np.exp(-tau_3D)
                F = L_first_3D * np.sum(w_unit * integrand_F)

                # Vacuum-vacuum branch — both surfaces absorbing.
                if alpha_in == 0.0 and alpha_out == 0.0:
                    psi_new[i, q_idx, p_idx] = F
                    continue

                # Single-transit B integrals.
                #
                # B_out: shell chord from INNER → OUTER at conserved b.
                # 2D parametrisation: u ∈ [0, L_step_2D], position along
                # trailing shell segment of the outer-cylinder antipodal
                # chord:  r²(u) = b² + (sqrt_disc_in + u)².
                # 3D arclength weight: integrand uses 3D attenuation
                # exp(-Σ_t · u / s_ip) and total is L_step_3D · Σ w · ig.
                #
                # B_in: shell chord from OUTER → INNER (reversal).
                # 2D parametrisation: r²(u) = b² + (sqrt_disc_out - u)².
                s_pts_2D_step = s_unit * L_step_2D

                # B_out integrand: r²(u) = b² + (sqrt_disc_in + u)².
                r_chord_out_sq = b_sq + (sqrt_disc_in + s_pts_2D_step) ** 2
                r_chord_out = np.sqrt(
                    np.clip(r_chord_out_sq, R_in_sq, R_out_sq)
                )
                tau_3D_step = sigma_t * s_pts_2D_step * inv_s_in_plane
                integrand_B_out = (
                    source_interp(r_chord_out) * np.exp(-tau_3D_step)
                )
                B_out = L_step_3D * np.sum(w_unit * integrand_B_out)

                # B_in integrand: r²(u) = b² + (sqrt_disc_out - u)².
                r_chord_in_sq = b_sq + (sqrt_disc_out - s_pts_2D_step) ** 2
                r_chord_in = np.sqrt(
                    np.clip(r_chord_in_sq, R_in_sq, R_out_sq)
                )
                integrand_B_in = (
                    source_interp(r_chord_in) * np.exp(-tau_3D_step)
                )
                B_in = L_step_3D * np.sum(w_unit * integrand_B_in)

                # Rank-2 closure. Mapping (slab-asym → annulus):
                #   α_L → α_in, α_R → α_out
                #   ψ_L^+ → ψ_in^out (outgoing from inner)
                #   ψ_R^- → ψ_out^in (outgoing from outer)
                #   B_LR  → B_out  (inner-toward-outer shell chord)
                #   B_RL  → B_in   (outer-toward-inner shell chord)
                #   surface='left' for inner first, 'right' for outer first.
                psi_new[i, q_idx, p_idx] = apply_variant_alpha_closure_rank2(
                    F=F, B_RL=B_in, B_LR=B_out,
                    tau_first_leg=sigma_t * L_first_3D,
                    tau_single_transit=sigma_t * L_step_3D,
                    alpha_left=alpha_in,
                    alpha_right=alpha_out,
                    surface=("left" if surface == "inner" else "right"),
                )

    return psi_new


# ═══════════════════════════════════════════════════════════════════════
# Scalar-flux extractor — same shape as cylinder
# ═══════════════════════════════════════════════════════════════════════


def _scalar_flux_from_psi(
    psi: np.ndarray,
    mu_axial_weights: np.ndarray,
    phi_az_weights: np.ndarray,
) -> np.ndarray:
    r"""Reduce angular flux to scalar flux on the radial grid.

    .. math::

        \phi(r) = \int_{-1}^{1}\!\mathrm d\mu_{\rm axial}
                    \int_0^{2\pi}\!\mathrm d\varphi_{\rm az}\,
                    \psi(r, \mu_{\rm axial}, \varphi_{\rm az}).
    """
    return np.einsum(
        '...rmp,m,p->...r', psi, mu_axial_weights, phi_az_weights,
    )


# ═══════════════════════════════════════════════════════════════════════
# 1G homogeneous solver (asymmetric BC)
# ═══════════════════════════════════════════════════════════════════════


def solve_greens_function_annulus(
    R_in: float,
    R_out: float,
    sigma_t: float,
    sigma_s: float,
    nu_sigma_f: float,
    *,
    alpha_in: float = 1.0,
    alpha_out: float = 1.0,
    n_r: int = 16,
    n_mu_axial: int = 12,
    n_phi_az: int = 24,
    n_traj_quad: int = 32,
    max_iter: int = 200,
    tol: float = 1e-10,
    initial_psi: np.ndarray | None = None,
) -> AnnulusGreensResult:
    r"""Power iteration on the annulus Variant α operator (homogeneous
    shell, isotropic scattering, rank-2 + outer-only closure).

    Solves the k-eigenvalue problem on a homogeneous infinite annulus
    (hollow cylinder) of inner radius :math:`R_{\rm in}` and outer
    radius :math:`R_{\rm out}` with INDEPENDENT specular reflectivities
    on the inner and outer surfaces. Phase-space is partitioned by
    impact parameter :math:`b = r\,|\sin\varphi_{\rm az}|` into
    outer-only (:math:`b > R_{\rm in}`, rank-1 closure) and through-ray
    (:math:`b \le R_{\rm in}`, rank-2 closure) subsets.

    Parameters
    ----------
    R_in : float
        Inner cavity radius (cm). Must be > 0 and < R_out.
    R_out : float
        Outer cylinder radius (cm).
    sigma_t, sigma_s, nu_sigma_f : float
        Group cross sections (1G). All in the shell material.
    alpha_in : float, default 1.0
        Inner-surface specular reflectivity in :math:`[0, 1]`.
    alpha_out : float, default 1.0
        Outer-surface specular reflectivity in :math:`[0, 1]`.
    n_r : int
        Radial Gauss-Legendre quadrature order on
        :math:`(R_{\rm in}, R_{\rm out})`.
    n_mu_axial : int
        Axial-cosine Gauss-Legendre order on :math:`[-1, 1]`.
    n_phi_az : int
        Azimuthal Gauss-Legendre order on :math:`[0, 2\pi)`.
    n_traj_quad : int
        Trajectory + bounce-chord quadrature order.
    max_iter, tol : int, float
        Power-iteration parameters.
    initial_psi : (n_r, n_mu_axial, n_phi_az) ndarray, optional

    Returns
    -------
    :class:`AnnulusGreensResult`
    """
    if not (0.0 <= alpha_in <= 1.0):
        raise ValueError(f"alpha_in = {alpha_in} must lie in [0, 1]")
    if not (0.0 <= alpha_out <= 1.0):
        raise ValueError(f"alpha_out = {alpha_out} must lie in [0, 1]")
    if not (0.0 < R_in < R_out):
        raise ValueError(
            f"Need 0 < R_in < R_out; got R_in={R_in}, R_out={R_out}"
        )

    # Radial grid: GL on [R_in, R_out] (open quadrature avoids the
    # surfaces at endpoints).
    r_quad_pts, r_quad_wts = np.polynomial.legendre.leggauss(n_r)
    r_nodes = R_in + (R_out - R_in) * 0.5 * (r_quad_pts + 1.0)
    r_weights = (R_out - R_in) * 0.5 * r_quad_wts

    # Axial-cosine grid: GL on [-1, 1].
    mu_quad_pts, mu_quad_wts = np.polynomial.legendre.leggauss(n_mu_axial)
    mu_axial_nodes = mu_quad_pts
    mu_axial_weights = mu_quad_wts

    # Azimuthal grid: GL on [0, 2π) — open quadrature avoids the
    # tangential ±π/2 grazing rays.
    phi_quad_pts, phi_quad_wts = np.polynomial.legendre.leggauss(n_phi_az)
    phi_az_nodes = np.pi * (phi_quad_pts + 1.0)
    phi_az_weights = np.pi * phi_quad_wts

    if initial_psi is not None:
        psi = np.asarray(initial_psi, dtype=float).copy()
        if psi.shape != (n_r, n_mu_axial, n_phi_az):
            raise ValueError(
                f"initial_psi shape must be ({n_r}, {n_mu_axial}, "
                f"{n_phi_az}); got {psi.shape}"
            )
    else:
        psi = np.ones((n_r, n_mu_axial, n_phi_az))

    sigma_a = sigma_t - sigma_s
    if sigma_a <= 0:
        raise ValueError(
            f"sigma_a = sigma_t - sigma_s = {sigma_a} ≤ 0; non-absorbing "
            "medium not supported."
        )
    k_eff = nu_sigma_f / sigma_a

    # Volume-integrated fission rate per unit cylinder length:
    # ∫_{R_in}^{R_out} νΣ_f φ(r) · 2π r dr.
    def fission_rate(phi_r: np.ndarray) -> float:
        return float(
            nu_sigma_f * 2.0 * np.pi
            * np.sum(phi_r * r_nodes * r_weights)
        )

    iterations = 0
    converged = False
    inv_4pi = 1.0 / (4.0 * np.pi)

    for it in range(max_iter):
        iterations = it + 1

        phi_r = _scalar_flux_from_psi(psi, mu_axial_weights, phi_az_weights)
        source_profile = inv_4pi * (sigma_s + nu_sigma_f / k_eff) * phi_r

        psi_new = _apply_operator_annulus(
            source_profile, r_nodes, mu_axial_nodes, phi_az_nodes,
            R_in, R_out, sigma_t, alpha_in, alpha_out,
            n_traj_quad=n_traj_quad,
        )

        phi_new = _scalar_flux_from_psi(
            psi_new, mu_axial_weights, phi_az_weights,
        )

        Frate_old = fission_rate(phi_r)
        Frate_new = fission_rate(phi_new)
        if Frate_old < 1e-30:
            raise RuntimeError(
                f"Fission rate vanished at iter {iterations}; "
                "non-multiplying medium."
            )
        k_new = k_eff * Frate_new / Frate_old
        psi_normed = psi_new / Frate_new

        rel_dk = abs(k_new - k_eff) / max(abs(k_eff), 1e-30)
        psi = psi_normed
        k_eff = k_new

        if rel_dk < tol:
            converged = True
            break

    phi_r = _scalar_flux_from_psi(psi, mu_axial_weights, phi_az_weights)

    return AnnulusGreensResult(
        k_eff=float(k_eff),
        psi=psi,
        phi=phi_r,
        r_nodes=r_nodes,
        mu_axial_nodes=mu_axial_nodes,
        phi_az_nodes=phi_az_nodes,
        iterations=iterations,
        converged=converged,
    )


# ═══════════════════════════════════════════════════════════════════════
# Multi-group homogeneous solver (asymmetric BC)
# ═══════════════════════════════════════════════════════════════════════


def solve_greens_function_annulus_mg(
    R_in: float,
    R_out: float,
    sigma_t: np.ndarray,        # (G,)
    sigma_s: np.ndarray,        # (G, G), [g_from, g_to]
    nu_sigma_f: np.ndarray,     # (G,)
    chi: np.ndarray | None = None,  # (G,)
    *,
    alpha_in: float = 1.0,
    alpha_out: float = 1.0,
    n_r: int = 16,
    n_mu_axial: int = 12,
    n_phi_az: int = 24,
    n_traj_quad: int = 32,
    max_iter: int = 300,
    tol: float = 1e-9,
    initial_psi: np.ndarray | None = None,
    initial_k: float | None = None,
) -> AnnulusGreensMGResult:
    r"""Multi-group annulus Variant α power iteration (homogeneous
    shell, isotropic scattering, rank-2 + outer-only closure).

    Multi-group analog of :func:`solve_greens_function_annulus`.
    Convention: ``sigma_s[g_from, g_to]``. At
    :math:`\alpha_{\rm in} = \alpha_{\rm out} = 1` reduces exactly to
    :func:`orpheus.derivations.common.eigenvalue.kinf_homogeneous`.

    Parameters
    ----------
    R_in, R_out : float
    sigma_t : (G,) ndarray
    sigma_s : (G, G) ndarray
    nu_sigma_f : (G,) ndarray
    chi : (G,) ndarray, optional
        Fission spectrum. Default: all-fast emission.
    alpha_in, alpha_out : float
    n_r, n_mu_axial, n_phi_az, n_traj_quad : int
    max_iter, tol : int, float
    initial_psi : (G, n_r, n_mu_axial, n_phi_az) ndarray, optional
    initial_k : float, optional

    Returns
    -------
    :class:`AnnulusGreensMGResult`
    """
    sigma_t = np.asarray(sigma_t, dtype=float)
    sigma_s = np.asarray(sigma_s, dtype=float)
    nu_sigma_f = np.asarray(nu_sigma_f, dtype=float)
    G = len(sigma_t)
    if sigma_s.shape != (G, G):
        raise ValueError(
            f"sigma_s shape must be ({G}, {G}); got {sigma_s.shape}"
        )
    if len(nu_sigma_f) != G:
        raise ValueError(
            f"nu_sigma_f length must be {G}; got {len(nu_sigma_f)}"
        )
    if chi is None:
        chi = np.zeros(G)
        chi[0] = 1.0
    else:
        chi = np.asarray(chi, dtype=float)
        if len(chi) != G:
            raise ValueError(f"chi length must be {G}; got {len(chi)}")
    if not (0.0 <= alpha_in <= 1.0):
        raise ValueError(f"alpha_in = {alpha_in} must lie in [0, 1]")
    if not (0.0 <= alpha_out <= 1.0):
        raise ValueError(f"alpha_out = {alpha_out} must lie in [0, 1]")
    if not (0.0 < R_in < R_out):
        raise ValueError(
            f"Need 0 < R_in < R_out; got R_in={R_in}, R_out={R_out}"
        )

    r_quad_pts, r_quad_wts = np.polynomial.legendre.leggauss(n_r)
    r_nodes = R_in + (R_out - R_in) * 0.5 * (r_quad_pts + 1.0)
    r_weights = (R_out - R_in) * 0.5 * r_quad_wts

    mu_quad_pts, mu_quad_wts = np.polynomial.legendre.leggauss(n_mu_axial)
    mu_axial_nodes = mu_quad_pts
    mu_axial_weights = mu_quad_wts

    phi_quad_pts, phi_quad_wts = np.polynomial.legendre.leggauss(n_phi_az)
    phi_az_nodes = np.pi * (phi_quad_pts + 1.0)
    phi_az_weights = np.pi * phi_quad_wts

    if initial_psi is not None:
        psi = np.asarray(initial_psi, dtype=float).copy()
        if psi.shape != (G, n_r, n_mu_axial, n_phi_az):
            raise ValueError(
                f"initial_psi shape must be ({G}, {n_r}, {n_mu_axial}, "
                f"{n_phi_az}); got {psi.shape}"
            )
    else:
        psi = np.ones((G, n_r, n_mu_axial, n_phi_az))

    if initial_k is not None:
        k_eff = float(initial_k)
    else:
        A = np.diag(sigma_t) - sigma_s.T
        F_mat = np.outer(chi, nu_sigma_f)
        eig = np.linalg.eigvals(np.linalg.solve(A, F_mat))
        k_eff = float(np.real(eig).max())
        if k_eff <= 0:
            raise ValueError(f"k_inf estimate non-positive ({k_eff}).")

    def total_fission_rate(phi_g_r: np.ndarray) -> float:
        F_r = np.einsum('g,gr->r', nu_sigma_f, phi_g_r)
        return float(2.0 * np.pi * np.sum(F_r * r_nodes * r_weights))

    iterations = 0
    converged = False
    inv_4pi = 1.0 / (4.0 * np.pi)

    for it in range(max_iter):
        iterations = it + 1

        phi_g = np.einsum(
            'grmp,m,p->gr', psi, mu_axial_weights, phi_az_weights,
        )  # (G, n_r)

        F_r = np.einsum('g,gr->r', nu_sigma_f, phi_g)
        scatter_source = np.einsum('sg,sr->gr', sigma_s, phi_g)
        fission_source = (chi[:, None] / k_eff) * F_r[None, :]
        source_profile_g = inv_4pi * (scatter_source + fission_source)

        psi_new = np.zeros_like(psi)
        for g in range(G):
            psi_new[g] = _apply_operator_annulus(
                source_profile_g[g], r_nodes, mu_axial_nodes,
                phi_az_nodes, R_in, R_out, float(sigma_t[g]),
                alpha_in, alpha_out, n_traj_quad=n_traj_quad,
            )

        phi_g_new = np.einsum(
            'grmp,m,p->gr', psi_new, mu_axial_weights, phi_az_weights,
        )
        Frate_old = total_fission_rate(phi_g)
        Frate_new = total_fission_rate(phi_g_new)
        if Frate_old < 1e-30:
            raise RuntimeError(
                f"Fission rate vanished at iter {iterations}."
            )
        k_new = k_eff * Frate_new / Frate_old
        psi_normed = psi_new / Frate_new

        rel_dk = abs(k_new - k_eff) / max(abs(k_eff), 1e-30)
        psi = psi_normed
        k_eff = k_new

        if rel_dk < tol:
            converged = True
            break

    phi_g = np.einsum(
        'grmp,m,p->gr', psi, mu_axial_weights, phi_az_weights,
    )

    return AnnulusGreensMGResult(
        k_eff=float(k_eff),
        psi_g=psi,
        phi_g=phi_g,
        r_nodes=r_nodes,
        mu_axial_nodes=mu_axial_nodes,
        phi_az_nodes=phi_az_nodes,
        iterations=iterations,
        converged=converged,
    )
