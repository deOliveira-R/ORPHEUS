r"""Phase-3C-1 hollow sphere **Variant α** Green's function reference
(homogeneous shell :math:`r \in [R_{\rm in}, R_{\rm out}]`, **independent**
specular reflectivities :math:`\alpha_{\rm in} \in [0, 1]` on the inner
cavity surface and :math:`\alpha_{\rm out} \in [0, 1]` on the outer
surface) — research-grade prototype.

Standalone module, parallel to :mod:`.greens_function` (solid sphere)
and :mod:`.greens_function_slab_asymmetric` (asymmetric slab). Implements
the **rank-2 boundary-to-boundary scattering resolvent** on the through-
ray phase-space subset (impact parameter :math:`b \le R_{\rm in}`) and
the **rank-1 outer-only resolvent** on the outer-only subset
(:math:`b > R_{\rm in}`).

Phase-space partition by impact parameter
------------------------------------------

For each interior point :math:`(r, \mu)` with :math:`r \in [R_{\rm in},
R_{\rm out}]` and :math:`\mu \in [-1, 1]`, the **impact parameter**
:math:`b = r\sqrt{1 - \mu^2}` is conserved along straight-line travel.
Two cases:

- :math:`b > R_{\rm in}` — outer-only rays. The trajectory line does not
  intersect the inner sphere; the particle bounces between two points
  on the OUTER surface only. **Topologically rank-1**, structurally
  identical to a solid-sphere ray at the same outer radius.
- :math:`b \le R_{\rm in}` — through-rays. The trajectory line crosses
  the inner cavity. Under interpretation (A) — see plan §"Hollow
  sphere geometry" — the inner surface is a specular reflector with
  reflectivity :math:`\alpha_{\rm in}`, so the particle bounces
  alternately inner ↔ outer. **Topologically rank-2**, structurally
  identical to the asymmetric slab with curvilinear chord algebra.

Interpretation (A): :math:`\alpha_{\rm in} = 0` means "particle is
lost to the cavity" (perfect inner absorber); :math:`\alpha_{\rm in}
= 1` means "perfect inner reflector". This makes the hollow sphere a
clean topological analog of the asymmetric slab.

Architecture
------------

Phase-space: :math:`(r, \mu)` with :math:`r \in [R_{\rm in}, R_{\rm
out}]`, signed :math:`\mu \in [-1, 1]`. For each interior point
:math:`(r_i, \mu_q)`:

1. **Compute** :math:`b = r_i\sqrt{1 - \mu_q^2}` and partition the
   trajectory into outer-only (:math:`b > R_{\rm in}`) or through-ray
   (:math:`b \le R_{\rm in}`).

2. **First-leg trajectory integral** along the backward chord to the
   first surface arrival:

   - Outer-only: :math:`L_{\rm first} = r\mu + \sqrt{R_{\rm out}^2 -
     b^2}`, surface = outer (rank-1 sphere formula).
   - Through-ray with :math:`\mu > 0`: backward goes inward; first
     arrival is INNER at :math:`L_{\rm first} = r\mu - \sqrt{R_{\rm
     in}^2 - b^2}`, surface = inner.
   - Through-ray with :math:`\mu < 0`: backward goes outward; first
     arrival is OUTER at :math:`L_{\rm first} = r\mu + \sqrt{R_{\rm
     out}^2 - b^2}`, surface = outer.

3. **Outer-only branch** (:math:`b > R_{\rm in}`): rank-1 closure on
   bounce-period chord :math:`L_p = 2\sqrt{R_{\rm out}^2 - b^2}`
   identical to solid sphere.

4. **Through-ray branch** (:math:`b \le R_{\rm in}`): rank-2 closure
   with single-transit (shell-traversal) optical depth
   :math:`\tau_{\rm step} = \Sigma_t\,(\sqrt{R_{\rm out}^2 - b^2}
   - \sqrt{R_{\rm in}^2 - b^2})`. The two single-transit B integrals
   :math:`B_{\rm in}` (outer→inner shell chord) and :math:`B_{\rm out}`
   (inner→outer shell chord) feed the rank-2 resolvent.

5. **Total angular flux**:
   :math:`\psi(r, \mu) = F + e^{-\Sigma_t L_{\rm first}}\,\psi_{\rm
   surface}` with the appropriate :math:`\psi_{\rm surface}` component.

6. **Scalar flux**: :math:`\phi(r) = 2\pi \int_{-1}^{1} \psi(r, \mu)\,
   \mathrm d\mu`.

Assumptions for Phase 3C-1 prototype
-------------------------------------

- Homogeneous shell (single :math:`\Sigma_t, \Sigma_s, \nu\Sigma_f`);
  multi-region deferred (Phase 3C-2 annulus is the next priority).
- Isotropic scattering.
- Independent reflective specular BCs :math:`\alpha_{\rm in} \in
  [0, 1]` and :math:`\alpha_{\rm out} \in [0, 1]`.

References
----------

- Sanchez, R. (1986). *Transp. Theor. Stat. Phys.* 14.
- :file:`/.claude/plans/peierls-greens-cylinder-and-2bc.md` — Phase
  3C-1 hollow sphere plan.
- :mod:`orpheus.derivations.continuous.peierls_greens_function.origins.specular.greens_function_hollow_sphere`
  — V_α1_hollow_sph/V_α2_hollow_sph/V_α3_hollow_sph SymPy verifications.
- :mod:`.greens_function_slab_asymmetric` — Phase-3B asymmetric slab
  reference solver (the rank-2 template lifted here to curvilinear
  2-surface).
- :mod:`.greens_function` — solid-sphere reference solver (the
  rank-1 template the outer-only branch reuses).
- :file:`.claude/agent-memory/cross-domain-attacker/variant_alpha_2surface_bie_frame.md`
  — frame match: rank-2 BIE block resolvent at 2-surface topologies.
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
class HollowSphereGreensResult:
    """Result of Variant α power iteration on homogeneous hollow sphere
    with INDEPENDENT specular reflectivities
    :math:`\\alpha_{\\rm in}, \\alpha_{\\rm out}`."""

    k_eff: float
    psi: np.ndarray  # (n_r, n_mu) angular flux on the grid
    phi: np.ndarray  # (n_r,) scalar flux on the radial grid
    r_nodes: np.ndarray
    mu_nodes: np.ndarray
    iterations: int
    converged: bool


@dataclass(frozen=True)
class HollowSphereGreensMGResult:
    """Result of multi-group Variant α power iteration on homogeneous
    hollow sphere with INDEPENDENT specular BCs."""

    k_eff: float
    psi_g: np.ndarray   # (G, n_r, n_mu)
    phi_g: np.ndarray   # (G, n_r)
    r_nodes: np.ndarray
    mu_nodes: np.ndarray
    iterations: int
    converged: bool


# ═══════════════════════════════════════════════════════════════════════
# Per-group Variant α operator (hollow sphere, rank-2 + impact-parameter
# partition)
# ═══════════════════════════════════════════════════════════════════════


def _apply_operator_hollow_sphere(
    source_profile: np.ndarray,
    r_nodes: np.ndarray,
    mu_nodes: np.ndarray,
    R_in: float,
    R_out: float,
    sigma_t: float,
    alpha_in: float,
    alpha_out: float,
    *,
    n_traj_quad: int,
) -> np.ndarray:
    r"""Per-group hollow-sphere Variant α operator with rank-2 + rank-1
    closure.

    For each :math:`(r_i, \mu_q)` with :math:`r_i \in [R_{\rm in},
    R_{\rm out}]`:

    1. Compute :math:`b = r_i\sqrt{1-\mu_q^2}`. Partition: outer-only
       (:math:`b > R_{\rm in}`) or through-ray (:math:`b \le R_{\rm
       in}`).

    2. Outer-only: identical to solid-sphere algebra. Rank-1 closure
       with :math:`\alpha = \alpha_{\rm out}`.

    3. Through-ray: rank-2 closure with single-transit shell-traversal
       :math:`\tau_{\rm step}`. First-leg backward goes to inner
       surface for :math:`\mu > 0` and outer surface for :math:`\mu <
       0`.

    Parameters
    ----------
    source_profile : (n_r,) ndarray
        :math:`q_g(r_i)/(4\pi)` — already-divided isotropic source
        per steradian per unit volume. Cubic-spline-interpolated for
        evaluation along trajectory points.
    r_nodes : (n_r,) ndarray
        Radial nodes on :math:`[R_{\rm in}, R_{\rm out}]`.
    mu_nodes : (n_mu,) ndarray
        Direction-cosine nodes on :math:`[-1, 1]`.
    R_in, R_out : float
        Inner cavity radius and outer surface radius.
    sigma_t : float
        Per-group total cross section in the shell.
    alpha_in, alpha_out : float
        Per-surface reflectivities in :math:`[0, 1]`.
    n_traj_quad : int
        Trajectory + bounce-chord quadrature order.

    Returns
    -------
    (n_r, n_mu) ndarray
    """
    # Cubic-spline interpolant on the shell. Source defined on r ∈
    # [R_in, R_out]; clip evaluations to keep within domain.
    source_interp = CubicSpline(r_nodes, source_profile, extrapolate=True)

    s_quad_raw, w_quad_raw = np.polynomial.legendre.leggauss(n_traj_quad)
    s_unit = 0.5 * (s_quad_raw + 1.0)
    w_unit = 0.5 * w_quad_raw

    n_r = len(r_nodes)
    n_mu = len(mu_nodes)
    psi_new = np.zeros((n_r, n_mu))

    R_out_sq = R_out * R_out
    R_in_sq = R_in * R_in

    for i in range(n_r):
        r = r_nodes[i]
        for q_idx in range(n_mu):
            mu = mu_nodes[q_idx]

            # Impact parameter b = r·sqrt(1 - μ²).
            b_sq = r * r * (1.0 - mu * mu)
            b = np.sqrt(max(b_sq, 0.0))

            # Outer-surface chord parameters (always defined).
            disc_out = R_out_sq - b_sq
            sqrt_disc_out = np.sqrt(max(disc_out, 0.0))

            if b > R_in:
                # ───────── Outer-only branch (rank-1) ─────────
                # First-leg backward: same as solid sphere.
                L_first = r * mu + sqrt_disc_out
                # Bounce-period chord across outer sphere.
                L_p = 2.0 * sqrt_disc_out

                # First-leg trajectory points.
                s_pts = s_unit * L_first
                r_traj_sq = r * r - 2.0 * r * mu * s_pts + s_pts * s_pts
                r_traj = np.sqrt(np.clip(r_traj_sq, R_in_sq, R_out_sq))
                integrand_F = (
                    source_interp(r_traj) * np.exp(-sigma_t * s_pts)
                )
                F = L_first * np.sum(w_unit * integrand_F)

                if alpha_out == 0.0:
                    psi_new[i, q_idx] = F
                    continue

                # Bounce-period chord (antipodal at conserved b).
                # Position: r_chord(s)² = b² + (s - L_p/2)².
                s_pts_p = s_unit * L_p
                r_chord_sq = b_sq + (s_pts_p - 0.5 * L_p) ** 2
                r_chord = np.sqrt(np.clip(r_chord_sq, R_in_sq, R_out_sq))
                integrand_B = (
                    source_interp(r_chord) * np.exp(-sigma_t * s_pts_p)
                )
                B = L_p * np.sum(w_unit * integrand_B)

                psi_new[i, q_idx] = apply_variant_alpha_closure(
                    F=F, B=B,
                    tau_first_leg=sigma_t * L_first,
                    tau_period=sigma_t * L_p,
                    alpha=alpha_out,
                )
                continue

            # ───────── Through-ray branch (rank-2) ─────────
            disc_in = R_in_sq - b_sq
            sqrt_disc_in = np.sqrt(max(disc_in, 0.0))

            # Single-transit shell-traversal length:
            #   L_step = sqrt(R_out² - b²) - sqrt(R_in² - b²).
            # τ_step = Σ_t · L_step.
            L_step = sqrt_disc_out - sqrt_disc_in

            # First-leg backward depending on sign of μ.
            if mu > 0:
                # Backward goes INWARD; first arrival is INNER surface.
                # s_arrival = r·μ - sqrt(R_in² - b²).
                L_first = r * mu - sqrt_disc_in
                # Surface label: outgoing-from-INNER (ψ_in^out).
                surface = "inner"
            else:
                # mu < 0 (or zero, treated as outward): backward goes
                # OUTWARD; first arrival is OUTER surface.
                # s_arrival = r·μ + sqrt(R_out² - b²).
                L_first = r * mu + sqrt_disc_out
                surface = "outer"

            # Numerical guard for tangent rays (μ exactly such that
            # b → R_in from below, sqrt_disc_in → 0). L_first remains
            # well-defined; just ensure non-negative.
            L_first = max(L_first, 0.0)

            # First-leg trajectory points along the backward chord.
            s_pts = s_unit * L_first
            r_traj_sq = r * r - 2.0 * r * mu * s_pts + s_pts * s_pts
            r_traj = np.sqrt(np.clip(r_traj_sq, R_in_sq, R_out_sq))
            integrand_F = (
                source_interp(r_traj) * np.exp(-sigma_t * s_pts)
            )
            F = L_first * np.sum(w_unit * integrand_F)

            # Vacuum-vacuum branch — both surfaces absorbing.
            if alpha_in == 0.0 and alpha_out == 0.0:
                psi_new[i, q_idx] = F
                continue

            # Single-transit B integrals.
            #
            # B_out: shell chord from INNER surface to OUTER surface
            # at conserved b. Length L_step. The chord starts at the
            # inner-surface intersection and ends at the outer-surface
            # intersection.
            #
            # Position along inner→outer chord: parametrise by
            # arclength s ∈ [0, L_step]. The chord is the half of the
            # full antipodal chord on the OUTER sphere that lies
            # between r = R_in (start) and r = R_out (end).
            #
            # Using the antipodal-chord parametrisation of the outer
            # sphere where r²(s_p) = b² + (s_p - L_p_outer/2)² with
            # L_p_outer = 2·sqrt(R_out² - b²), the inner surface is
            # crossed at s_p = L_p_outer/2 ± sqrt(R_in² - b²) =
            # sqrt_disc_out ± sqrt_disc_in. The two crossings define
            # the cavity transit; the TWO shell segments are
            # (0, L_step) and (L_step + 2·sqrt_disc_in, L_p_outer).
            #
            # For B_out (inner → outer) we use the TRAILING shell
            # segment: s_p ∈ [L_step + 2·sqrt_disc_in, L_p_outer],
            # which has length L_step. Reparametrise as u ∈ [0, L_step]:
            #   s_p(u) = L_step + 2·sqrt_disc_in + u
            #   r²(u) = b² + (s_p(u) - L_p_outer/2)²
            #         = b² + (sqrt_disc_in + u)²    [shell traversal
            #                                        going outward]
            # Source-line integral with chord-arclength s = u:
            #   B_out(b; q) = ∫_0^{L_step} q(r(u)) · e^{-Σ_t u} du
            #
            # For B_in (outer → inner), by reversal symmetry the
            # source-line integral has the same r(u) shape (outward
            # from inner means r increases — same chord, traversed in
            # opposite direction), but parametrised so the source
            # contribution AT the inner surface comes from points that
            # are FARTHER along: r²(u) = b² + (sqrt_disc_out - u)².
            s_pts_step = s_unit * L_step

            # B_out integrand: r²(u) = b² + (sqrt_disc_in + u)².
            r_chord_out_sq = b_sq + (sqrt_disc_in + s_pts_step) ** 2
            r_chord_out = np.sqrt(
                np.clip(r_chord_out_sq, R_in_sq, R_out_sq)
            )
            integrand_B_out = (
                source_interp(r_chord_out)
                * np.exp(-sigma_t * s_pts_step)
            )
            B_out = L_step * np.sum(w_unit * integrand_B_out)

            # B_in integrand: r²(u) = b² + (sqrt_disc_out - u)².
            r_chord_in_sq = b_sq + (sqrt_disc_out - s_pts_step) ** 2
            r_chord_in = np.sqrt(
                np.clip(r_chord_in_sq, R_in_sq, R_out_sq)
            )
            integrand_B_in = (
                source_interp(r_chord_in)
                * np.exp(-sigma_t * s_pts_step)
            )
            B_in = L_step * np.sum(w_unit * integrand_B_in)

            # Rank-2 closure. The slab-asymmetric primitive expects
            # surface ∈ {'left', 'right'} — map:
            #   inner ↔ left  (outgoing from inner toward outer = ψ_L^+)
            #   outer ↔ right (outgoing from outer toward inner = ψ_R^-)
            #
            # The slab closure equations:
            #   ψ_L^+ = α_L · B_LR + α_L · e^{-τ} · α_R · B_RL ...
            #   ψ_R^- = α_R · e^{-τ} · α_L · B_LR + α_R · B_RL ...
            # Mapping:
            #   α_L → α_in, α_R → α_out
            #   ψ_L^+ → ψ_in^out (outgoing from inner)
            #   ψ_R^- → ψ_out^in (outgoing from outer)
            #   B_LR (chord from x=0 going +x = "from left wall")
            #        ↔ B_out (chord from inner going outward — emitted
            #          at inner, integrated outward toward outer);
            #          this feeds ψ_in^out (the "left" surface).
            #   B_RL (chord from x=L going -x = "from right wall")
            #        ↔ B_in (chord from outer going inward — emitted
            #          at outer, integrated inward toward inner);
            #          this feeds ψ_out^in (the "right" surface).
            psi_new[i, q_idx] = apply_variant_alpha_closure_rank2(
                F=F, B_RL=B_in, B_LR=B_out,
                tau_first_leg=sigma_t * L_first,
                tau_single_transit=sigma_t * L_step,
                alpha_left=alpha_in,
                alpha_right=alpha_out,
                surface=("left" if surface == "inner" else "right"),
            )

    return psi_new


# ═══════════════════════════════════════════════════════════════════════
# Scalar-flux extractor — shared with sphere
# ═══════════════════════════════════════════════════════════════════════


def _scalar_flux_from_psi(
    psi: np.ndarray,
    mu_weights: np.ndarray,
) -> np.ndarray:
    r"""Reduce angular flux to scalar flux on the radial grid.

    .. math::

        \phi(r) = 2\pi \int_{-1}^{1} \psi(r, \mu)\,\mathrm d\mu.
    """
    return 2.0 * np.pi * np.einsum('...rm,m->...r', psi, mu_weights)


# ═══════════════════════════════════════════════════════════════════════
# 1G homogeneous solver (asymmetric BC)
# ═══════════════════════════════════════════════════════════════════════


def solve_greens_function_hollow_sphere(
    R_in: float,
    R_out: float,
    sigma_t: float,
    sigma_s: float,
    nu_sigma_f: float,
    *,
    alpha_in: float = 1.0,
    alpha_out: float = 1.0,
    n_r: int = 24,
    n_mu: int = 24,
    n_traj_quad: int = 64,
    max_iter: int = 200,
    tol: float = 1e-10,
    initial_psi: np.ndarray | None = None,
) -> HollowSphereGreensResult:
    r"""Power iteration on the hollow-sphere Variant α operator
    (homogeneous shell, isotropic scattering, rank-2 + outer-only
    closure).

    Solves the k-eigenvalue problem on a homogeneous shell with
    INDEPENDENT specular reflectivities on the inner and outer
    surfaces. Phase-space is partitioned by impact parameter
    :math:`b = r\sqrt{1 - \mu^2}` into outer-only (:math:`b > R_{\rm
    in}`, rank-1 closure) and through-ray (:math:`b \le R_{\rm in}`,
    rank-2 closure) subsets.

    Parameters
    ----------
    R_in : float
        Inner cavity radius (cm). Must be > 0 and < R_out.
    R_out : float
        Outer surface radius (cm).
    sigma_t, sigma_s, nu_sigma_f : float
        Group cross sections (1G). All in the shell material.
    alpha_in : float, default 1.0
        Inner-surface specular reflectivity in :math:`[0, 1]`.
    alpha_out : float, default 1.0
        Outer-surface specular reflectivity in :math:`[0, 1]`.
    n_r, n_mu, n_traj_quad : int
    max_iter, tol : int, float
    initial_psi : (n_r, n_mu) ndarray, optional

    Returns
    -------
    :class:`HollowSphereGreensResult`
    """
    if not (0.0 <= alpha_in <= 1.0):
        raise ValueError(f"alpha_in = {alpha_in} must lie in [0, 1]")
    if not (0.0 <= alpha_out <= 1.0):
        raise ValueError(f"alpha_out = {alpha_out} must lie in [0, 1]")
    if not (0.0 < R_in < R_out):
        raise ValueError(
            f"Need 0 < R_in < R_out; got R_in={R_in}, R_out={R_out}"
        )

    # Radial grid: GL on [R_in, R_out] (open quadrature avoids
    # endpoints which are surfaces).
    r_quad_pts, r_quad_wts = np.polynomial.legendre.leggauss(n_r)
    r_nodes = R_in + (R_out - R_in) * 0.5 * (r_quad_pts + 1.0)
    r_weights = (R_out - R_in) * 0.5 * r_quad_wts

    # Direction-cosine grid: GL on [-1, 1].
    mu_quad_pts, mu_quad_wts = np.polynomial.legendre.leggauss(n_mu)
    mu_nodes = mu_quad_pts
    mu_weights = mu_quad_wts

    if initial_psi is not None:
        psi = np.asarray(initial_psi, dtype=float).copy()
        if psi.shape != (n_r, n_mu):
            raise ValueError(
                f"initial_psi shape must be ({n_r}, {n_mu}); "
                f"got {psi.shape}"
            )
    else:
        psi = np.ones((n_r, n_mu))

    sigma_a = sigma_t - sigma_s
    if sigma_a <= 0:
        raise ValueError(
            f"sigma_a = sigma_t - sigma_s = {sigma_a} ≤ 0; non-absorbing "
            "medium not supported."
        )
    k_eff = nu_sigma_f / sigma_a

    # Volume-integrated fission rate for Rayleigh quotient.
    # Shell volume integral: ∫_{R_in}^{R_out} f(r) · 4π r² dr.
    def fission_rate(phi_r: np.ndarray) -> float:
        return float(
            nu_sigma_f * 4.0 * np.pi
            * np.sum(phi_r * r_nodes ** 2 * r_weights)
        )

    iterations = 0
    converged = False
    inv_4pi = 1.0 / (4.0 * np.pi)

    for it in range(max_iter):
        iterations = it + 1

        phi_r = _scalar_flux_from_psi(psi, mu_weights)
        source_profile = inv_4pi * (sigma_s + nu_sigma_f / k_eff) * phi_r

        psi_new = _apply_operator_hollow_sphere(
            source_profile, r_nodes, mu_nodes, R_in, R_out,
            sigma_t, alpha_in, alpha_out,
            n_traj_quad=n_traj_quad,
        )

        phi_new = _scalar_flux_from_psi(psi_new, mu_weights)

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

    phi_r = _scalar_flux_from_psi(psi, mu_weights)

    return HollowSphereGreensResult(
        k_eff=float(k_eff),
        psi=psi,
        phi=phi_r,
        r_nodes=r_nodes,
        mu_nodes=mu_nodes,
        iterations=iterations,
        converged=converged,
    )


# ═══════════════════════════════════════════════════════════════════════
# Multi-group homogeneous solver (asymmetric BC)
# ═══════════════════════════════════════════════════════════════════════


def solve_greens_function_hollow_sphere_mg(
    R_in: float,
    R_out: float,
    sigma_t: np.ndarray,        # (G,)
    sigma_s: np.ndarray,        # (G, G), [g_from, g_to]
    nu_sigma_f: np.ndarray,     # (G,)
    chi: np.ndarray | None = None,  # (G,)
    *,
    alpha_in: float = 1.0,
    alpha_out: float = 1.0,
    n_r: int = 24,
    n_mu: int = 24,
    n_traj_quad: int = 64,
    max_iter: int = 300,
    tol: float = 1e-9,
    initial_psi: np.ndarray | None = None,
    initial_k: float | None = None,
) -> HollowSphereGreensMGResult:
    r"""Multi-group hollow-sphere Variant α power iteration (homogeneous
    shell, isotropic scattering, rank-2 + outer-only closure).

    Multi-group analog of :func:`solve_greens_function_hollow_sphere`.
    Convention: ``sigma_s[g_from, g_to]``. At :math:`\alpha_{\rm in} =
    \alpha_{\rm out} = 1` reduces exactly to
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
    n_r, n_mu, n_traj_quad : int
    max_iter, tol : int, float
    initial_psi : (G, n_r, n_mu) ndarray, optional
    initial_k : float, optional

    Returns
    -------
    :class:`HollowSphereGreensMGResult`
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

    mu_quad_pts, mu_quad_wts = np.polynomial.legendre.leggauss(n_mu)
    mu_nodes = mu_quad_pts
    mu_weights = mu_quad_wts

    if initial_psi is not None:
        psi = np.asarray(initial_psi, dtype=float).copy()
        if psi.shape != (G, n_r, n_mu):
            raise ValueError(
                f"initial_psi shape must be ({G}, {n_r}, {n_mu}); "
                f"got {psi.shape}"
            )
    else:
        psi = np.ones((G, n_r, n_mu))

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
        return float(
            4.0 * np.pi * np.sum(F_r * r_nodes ** 2 * r_weights)
        )

    iterations = 0
    converged = False
    inv_4pi = 1.0 / (4.0 * np.pi)

    for it in range(max_iter):
        iterations = it + 1

        phi_g = 2.0 * np.pi * np.einsum(
            'grm,m->gr', psi, mu_weights,
        )

        F_r = np.einsum('g,gr->r', nu_sigma_f, phi_g)
        scatter_source = np.einsum('sg,sr->gr', sigma_s, phi_g)
        fission_source = (chi[:, None] / k_eff) * F_r[None, :]
        source_profile_g = inv_4pi * (scatter_source + fission_source)

        psi_new = np.zeros_like(psi)
        for g in range(G):
            psi_new[g] = _apply_operator_hollow_sphere(
                source_profile_g[g], r_nodes, mu_nodes, R_in, R_out,
                float(sigma_t[g]), alpha_in, alpha_out,
                n_traj_quad=n_traj_quad,
            )

        phi_g_new = 2.0 * np.pi * np.einsum(
            'grm,m->gr', psi_new, mu_weights,
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

    phi_g = 2.0 * np.pi * np.einsum('grm,m->gr', psi, mu_weights)

    return HollowSphereGreensMGResult(
        k_eff=float(k_eff),
        psi_g=psi,
        phi_g=phi_g,
        r_nodes=r_nodes,
        mu_nodes=mu_nodes,
        iterations=iterations,
        converged=converged,
    )
