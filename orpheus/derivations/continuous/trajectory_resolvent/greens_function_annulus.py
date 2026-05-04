r"""Phase-3C-2 annulus (hollow cylinder) **Variant α** Green's function
reference (homogeneous shell :math:`r \in [R_{\rm in}, R_{\rm out}]`,
**independent** specular reflectivities :math:`\alpha_{\rm in} \in
[0, 1]` on the inner cavity surface and :math:`\alpha_{\rm out} \in
[0, 1]` on the outer surface) — research-grade prototype.

Standalone module, parallel to :mod:`.greens_function_cylinder` (solid
cylinder) and :mod:`.greens_function_hollow_sphere` (hollow sphere).
This is the **last two-surface orbit-space class instance** in the
Variant α plan; once it ships, the 6-geometry × 2-orbit-space-class
family (sphere, cylinder, slab, slab-asym, hollow sphere, annulus)
is complete on the unified rank-1 / rank-2 framework via
:mod:`.variant_alpha_core`. See Sphinx
§\ ``orbit-space-m-g-classification`` for the structural signature
that determines closure rank from M/G endpoint count.

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
  between two points on the OUTER cylinder only. **Closure rank-1**
  (orbit-space class: one-surface-compact), structurally identical
  to a solid-cylinder ray at the same outer radius and impact
  parameter (with axial correction).
- :math:`b \le R_{\rm in}` — through-ray. The 2D in-plane ray crosses
  the inner cavity circle. Under interpretation (A) — the inner
  surface is a specular reflector with reflectivity :math:`\alpha_{\rm
  in}` — the particle bounces alternately inner ↔ outer.
  **Closure rank-2** (orbit-space class: two-surface — M/G interval
  with two distinct BC endpoints under :math:`\mathbb{R} \times
  SO(2)` action with inner-radius cut).

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
- :mod:`orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function_annulus`
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

from orpheus.derivations.continuous.trajectory_resolvent.chord_oracle import (
    AnnulusChordOracle,
)
from orpheus.derivations.continuous.trajectory_resolvent.power_iteration import (
    power_iterate_variant_alpha,
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

    Thin facade over :class:`AnnulusChordOracle.apply_operator` (the
    R3 ChordOracle Protocol). The :math:`b`-partition + 3D arclength
    lift + dual rank-1 (outer-only) / rank-2 (through-ray) closure
    routing live in the oracle; this function preserves the legacy
    call signature.

    Bit-equal with the pre-R3 inlined body.
    """
    oracle = AnnulusChordOracle(
        r_nodes=r_nodes,
        mu_axial_nodes=mu_axial_nodes,
        phi_az_nodes=phi_az_nodes,
        R_in=R_in, R_out=R_out,
        alpha_in=alpha_in, alpha_out=alpha_out,
    )
    return oracle.apply_operator(
        source_profile, sigma_t=sigma_t, n_traj_quad=n_traj_quad,
    )


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

    inv_4pi = 1.0 / (4.0 * np.pi)

    def _step(psi_iter, k_iter):
        phi_r = _scalar_flux_from_psi(
            psi_iter, mu_axial_weights, phi_az_weights,
        )
        source_profile = inv_4pi * (sigma_s + nu_sigma_f / k_iter) * phi_r

        psi_new = _apply_operator_annulus(
            source_profile, r_nodes, mu_axial_nodes, phi_az_nodes,
            R_in, R_out, sigma_t, alpha_in, alpha_out,
            n_traj_quad=n_traj_quad,
        )

        phi_new = _scalar_flux_from_psi(
            psi_new, mu_axial_weights, phi_az_weights,
        )

        return psi_new, fission_rate(phi_r), fission_rate(phi_new)

    pi_result = power_iterate_variant_alpha(
        _step, psi, initial_k=k_eff, max_iter=max_iter, tol=tol,
    )
    psi = pi_result.psi
    k_eff = pi_result.k_eff
    iterations = pi_result.iterations
    converged = pi_result.converged

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

    inv_4pi = 1.0 / (4.0 * np.pi)

    def _step(psi_iter, k_iter):
        phi_g = np.einsum(
            'grmp,m,p->gr', psi_iter, mu_axial_weights, phi_az_weights,
        )  # (G, n_r)

        F_r = np.einsum('g,gr->r', nu_sigma_f, phi_g)
        scatter_source = np.einsum('sg,sr->gr', sigma_s, phi_g)
        fission_source = (chi[:, None] / k_iter) * F_r[None, :]
        source_profile_g = inv_4pi * (scatter_source + fission_source)

        psi_new = np.zeros_like(psi_iter)
        for g in range(G):
            psi_new[g] = _apply_operator_annulus(
                source_profile_g[g], r_nodes, mu_axial_nodes,
                phi_az_nodes, R_in, R_out, float(sigma_t[g]),
                alpha_in, alpha_out, n_traj_quad=n_traj_quad,
            )

        phi_g_new = np.einsum(
            'grmp,m,p->gr', psi_new, mu_axial_weights, phi_az_weights,
        )
        return psi_new, total_fission_rate(phi_g), total_fission_rate(phi_g_new)

    pi_result = power_iterate_variant_alpha(
        _step, psi, initial_k=k_eff, max_iter=max_iter, tol=tol,
    )
    psi = pi_result.psi
    k_eff = pi_result.k_eff
    iterations = pi_result.iterations
    converged = pi_result.converged

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
