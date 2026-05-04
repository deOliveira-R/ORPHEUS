r"""Chord oracle — geometry-specific chord arithmetic + Variant α operator.

This module institutionalises the **chord oracle abstraction** named by
the cross-domain-attacker memo
``.claude/agent-memory/cross-domain-attacker/trajectory_resolvent_foreign_frames.md``
as the *fiber-bundle base atlas* of the Variant α family. A *chord
oracle* is the per-geometry primitive that maps the 1D abstract orbit
space M/G into ray-traversal data on the original manifold M:

1. **First-leg backward chord** — :math:`L_{\rm first}` and the
   parametrised position-along-chord :math:`r(s)` (or :math:`x(s)`).
2. **Bounce-period chord (rank-1)** OR **single-transit shell chord
   (rank-2)** — the optical depth :math:`\tau` of one closure period
   plus the position-along-chord parametrisation needed to compute the
   :math:`B` source-line integral.
3. **Closure routing** — for two-surface orbit-spaces the :math:`b`-
   partition selects rank-1 vs rank-2 closure, and rank-2 needs a
   *surface-arrival-first* label (`'left'` / `'right'` /
   `'inner'` / `'outer'`).
4. **3D angular Jacobian** — for cylinder / annulus, the 2D in-plane
   chord lifts to 3D via :math:`s_{\rm 3D} = s_{\rm 2D} / \sqrt{1 -
   \mu_{\rm axial}^{2}}`.

Each geometry exposes a concrete :class:`ChordOracle` that hides these
geometry-specific details behind a single
:meth:`~ChordOracle.apply_operator` method that produces the new
angular flux array :math:`\psi^{(n+1)}` from a source profile + total
cross section. The closure ranks (rank-1 / rank-2) and the multi-bounce
geometric series are applied through the SHARED
:mod:`.variant_alpha_core` primitives that are invariant across
geometries (Phase-2 unification verdict).

R3 hindsight refactor (algebra-of-record discipline)
----------------------------------------------------

This module is the third hindsight refactor (R3) of the Variant α
family — see ``.claude/plans/trajectory_resolvent_hindsight_refactor.md``
for the full plan. R1 collapsed the 11 byte-identical power-iteration
outer loops into :func:`~.power_iteration.power_iterate_variant_alpha`;
R2 unified the 12 result dataclasses behind the
:class:`.billiard.Billiard` math-rich facade and the SHARED
:class:`~orpheus.derivations.common.solution_types.CriticalSolution`
type; R3 (this module) extracts the chord-arithmetic primitives.

Bit-equality is preserved exactly. The body of every oracle's
:meth:`apply_operator` is **byte-equivalent** to the per-geometry
``_apply_operator_*`` function it replaces — only the location and
encapsulation change. Every FP operation runs in the same order;
IEEE-754 reproducibility guarantees identical numerical output. The 93
trajectory_resolvent tests + 84 cross-method tests must pass with
identical k_eff and ψ values pre- and post-refactor (verified at the
exact-bit level via ``float.hex(k_eff)`` for the canonical 12
configurations).

Why the abstraction lives here (not in :mod:`derivations.common` yet)
--------------------------------------------------------------------

Per :file:`.claude/agent-memory/method-implementer/feedback_unify_after_two_instances.md`
discipline (validated across the family's own Phase-1 / Phase-2 history),
abstractions promote to ``common/`` only after at least TWO consumers
have validated them. Today only the Variant α family consumes the
oracle Protocol; the second consumer is :mod:`.peierls_nystrom`
(production CP variant) and the third is MoC. Promotion is gated on
those second / third consumers landing — see plan §R4.

Cross-pollination map (deferred — for the promotion):

- :mod:`orpheus.derivations.common.kernels.chord_half_lengths` — the
  proto-ChordOracle for concentric annuli / shells. Specialises the
  conserved-impact-parameter chord :math:`L = 2\sqrt{R^2 - b^2}` that
  this module's :class:`SphereChordOracle` and :class:`CylinderChordOracle`
  use. Will be absorbed when the Protocol moves to ``common/``.
- :mod:`orpheus.derivations.continuous.peierls_nystrom` — production
  flat-source CP chord arithmetic. Similar primitives (impact-parameter
  segmentation, chord positions) but expressed differently; the oracle
  abstraction is the natural unifier.

Cross-domain frame (the BaseAtlas in the fiber-bundle picture)
--------------------------------------------------------------

In the fiber-bundle frame validated by
``variant_alpha_family_hindsight.md``, the Variant α family is one
*total space* with three coordinated layers:

1. **BaseAtlas** (this module) — the per-chart parametrisation of the
   spatial manifold M; produces chord lengths + chord positions.
2. **AngularFiber** — the residual symmetry group's parametrisation:
   :math:`\mu` for sphere / hollow_sphere; :math:`(\mu_{\rm axial},
   \varphi_{\rm az})` for cylinder / annulus; signed :math:`\mu` for
   slab. Lives in each geometry's solver entry point (mesh
   construction).
3. **ClosurePrimitives** — :mod:`.variant_alpha_core` carries the
   rank-1 / rank-2 resolvent + closure invariants.

References
----------

- Sanchez, R. (1986). *Transp. Theor. Stat. Phys.* **14**, 69–93. The
  Variant α formalism (Eqs. A3-A5 specular leg).
- Cross-domain-attacker memo:
  ``.claude/agent-memory/cross-domain-attacker/trajectory_resolvent_foreign_frames.md``
  (the BaseAtlas / fiber-bundle naming).
- Method-implementer hindsight memo:
  ``.claude/agent-memory/method-implementer/trajectory_resolvent_hindsight_review.md``
  (the R3 specification).
- :mod:`.variant_alpha_core` — the shared rank-1 / rank-2 closure
  primitives that every oracle composes with.
- :func:`orpheus.derivations.common.kernels.chord_half_lengths` — the
  proto-ChordOracle for concentric annuli / shells.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import numpy as np
from scipy.interpolate import CubicSpline

from orpheus.derivations.continuous.trajectory_resolvent.variant_alpha_core import (
    apply_variant_alpha_closure,
    apply_variant_alpha_closure_rank2,
)


__all__ = [
    "ChordOracle",
    "SphereChordOracle",
    "MultiRegionSphereChordOracle",
    "CylinderChordOracle",
    "SlabAsymmetricChordOracle",
    "HollowSphereChordOracle",
    "AnnulusChordOracle",
]


# ─────────────────────────────────────────────────────────────────────
# Protocol — what every chord oracle exposes
# ─────────────────────────────────────────────────────────────────────


@runtime_checkable
class ChordOracle(Protocol):
    r"""The fiber-bundle base atlas — geometry-specific chord arithmetic.

    A *chord oracle* owns three pieces of geometry-specific knowledge
    that, taken together, fully specify how the Variant α operator is
    applied to a source profile on this geometry:

    1. **Chord lengths** along the backward characteristic and the
       bounce period (or single-transit) at each phase-space point.
    2. **Chord position parametrisations** :math:`r(s)` (or
       :math:`x(s)`) for evaluating the source-interpolant along the
       chord.
    3. **Closure routing** — which rank (1 or 2) of the
       :mod:`.variant_alpha_core` closure applies, and (for rank-2)
       which surface is the first arrival.

    The oracle exposes a single high-level entry point
    :meth:`apply_operator` that produces the new angular flux array
    :math:`\psi^{(n+1)}` from a source profile and total cross section.
    Geometry-specific details (3D Jacobian for cylinder/annulus,
    :math:`b`-partition for hollow geometries, sign-of-:math:`\mu`
    branching for slab) are encapsulated.

    Concrete oracles are constructed once per solver call (with the
    geometry's grid + alpha payload baked in) and then queried by the
    power-iteration driver via :meth:`apply_operator` at every
    iteration.

    Notes
    -----
    The Protocol is :func:`runtime_checkable` so call sites can use
    ``isinstance(oracle, ChordOracle)`` for defensive routing.
    Concrete oracles are frozen dataclasses (immutable, hashable) so
    they can be safely cached across power-iteration calls.
    """

    def apply_operator(
        self,
        source_profile: np.ndarray,
        sigma_t: float,
        *,
        n_traj_quad: int,
    ) -> np.ndarray:
        r"""Apply the Variant α operator to the source profile.

        Evaluates :math:`\psi^{(n+1)}_g = K_{\alpha}\,q_g/(4\pi)` for
        the geometry-specific Variant α operator
        :math:`K_{\alpha}`. The closure rank, the surface routing, and
        the 3D Jacobian (if any) are handled internally.

        Parameters
        ----------
        source_profile : (n_spatial,) ndarray
            The per-spatial-node isotropic source per steradian per unit
            volume :math:`q_g(r_i)/(4\pi)`. Cubic-spline-interpolated
            for evaluation along the chord.
        sigma_t : float
            Per-group total cross section.
        n_traj_quad : int
            Gauss-Legendre order for the trajectory + bounce-period
            (or single-transit) integrals.

        Returns
        -------
        ndarray
            Updated angular flux for this group, shape geometry-
            specific (``(n_r, n_mu)`` for sphere/hollow_sphere,
            ``(n_r, n_mu_axial, n_phi_az)`` for cylinder/annulus,
            ``(n_x, n_mu)`` for slab).
        """
        ...


# ─────────────────────────────────────────────────────────────────────
# Sphere oracle — solid sphere, rank-1 closure
# ─────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class SphereChordOracle:
    r"""Chord oracle for the homogeneous solid sphere.

    Phase-space: :math:`(r, \mu)` with :math:`r \in (0, R)` and
    :math:`\mu \in [-1, 1]`. Geometry primitives:

    - First-leg backward chord :math:`L_{\rm back} = r\mu + \sqrt{R^2
      - r^2(1-\mu^2)}`.
    - Bounce-period chord :math:`L_p = 2 R \mu_{\rm surf}` with
      :math:`\mu_{\rm surf} = \sqrt{R^2 - r^2(1-\mu^2)}/R` (the
      conserved cosine at the surface).
    - Rank-1 closure with reflectivity :math:`\alpha`.

    Bit-equality contract: the body of :meth:`apply_operator` is a
    verbatim relocation of the pre-R3
    ``_apply_operator_with_source_profile`` in
    :mod:`.greens_function`. Every FP operation runs in the same order.

    Attributes
    ----------
    r_nodes : (n_r,) ndarray
        Radial Gauss-Legendre nodes on :math:`(0, R)`.
    mu_nodes : (n_mu,) ndarray
        Direction-cosine Gauss-Legendre nodes on :math:`[-1, 1]`.
    R : float
        Outer radius.
    alpha : float
        Surface specular reflectivity in :math:`[0, 1]`.
    """

    r_nodes: np.ndarray
    mu_nodes: np.ndarray
    R: float
    alpha: float

    def apply_operator(
        self,
        source_profile: np.ndarray,
        sigma_t: float,
        *,
        n_traj_quad: int,
    ) -> np.ndarray:
        r"""Sphere Variant α operator. See
        :meth:`ChordOracle.apply_operator`.

        Bit-equal with the pre-R3
        :func:`._apply_operator_with_source_profile` body.
        """
        r_nodes = self.r_nodes
        mu_nodes = self.mu_nodes
        R = self.R
        alpha = self.alpha

        source_interp = CubicSpline(r_nodes, source_profile, extrapolate=True)

        s_quad_raw, w_quad_raw = np.polynomial.legendre.leggauss(n_traj_quad)
        s_unit = 0.5 * (s_quad_raw + 1.0)
        w_unit = 0.5 * w_quad_raw

        n_r = len(r_nodes)
        n_mu = len(mu_nodes)
        psi_new = np.zeros((n_r, n_mu))

        for i in range(n_r):
            r = r_nodes[i]
            for q_idx in range(n_mu):
                mu = mu_nodes[q_idx]

                disc = R * R - r * r * (1.0 - mu * mu)
                sqrt_disc = np.sqrt(max(disc, 0.0))
                L_back = r * mu + sqrt_disc
                mu_surf = sqrt_disc / R
                L_p = 2.0 * R * mu_surf

                s_pts = s_unit * L_back
                r_traj_sq = r * r - 2.0 * r * mu * s_pts + s_pts * s_pts
                r_traj = np.sqrt(np.clip(r_traj_sq, 0.0, R * R))
                integrand_F = source_interp(r_traj) * np.exp(-sigma_t * s_pts)
                F = L_back * np.sum(w_unit * integrand_F)

                if alpha == 0.0:
                    psi_new[i, q_idx] = F
                    continue

                s_pts_p = s_unit * L_p
                h_sq = R * R * (1.0 - mu_surf * mu_surf)
                r_chord_sq = h_sq + (s_pts_p - 0.5 * L_p) ** 2
                r_chord = np.sqrt(np.clip(r_chord_sq, 0.0, R * R))
                integrand_B = source_interp(r_chord) * np.exp(-sigma_t * s_pts_p)
                B = L_p * np.sum(w_unit * integrand_B)

                psi_new[i, q_idx] = apply_variant_alpha_closure(
                    F=F, B=B,
                    tau_first_leg=sigma_t * L_back,
                    tau_period=sigma_t * L_p,
                    alpha=alpha,
                )

        return psi_new


# ─────────────────────────────────────────────────────────────────────
# Multi-region sphere oracle — piecewise-:math:`\Sigma_t` decomposition
# ─────────────────────────────────────────────────────────────────────


def _region_at_radius_oracle(r: float, radii: np.ndarray) -> int:
    """Return region index containing radius :math:`r`.

    Helper used by :class:`MultiRegionSphereChordOracle`. Convention:
    region :math:`k` occupies :math:`(r_{k-1}, r_k]` for :math:`k \\ge 1`
    and :math:`(0, r_0]` for :math:`k = 0`.
    """
    idx = int(np.searchsorted(radii, r, side="left"))
    return min(idx, len(radii) - 1)


def _trajectory_segments_oracle(
    r_start: float, mu: float, R: float, radii: np.ndarray,
) -> tuple[list[tuple[float, float, int]], float]:
    """Backward-trajectory segments crossing region boundaries.

    Helper used by :class:`MultiRegionSphereChordOracle`. Bit-equal
    with the pre-R3 :func:`_trajectory_segments` in :mod:`.greens_function`.
    """
    L_back = r_start * mu + np.sqrt(
        R * R - r_start * r_start * (1.0 - mu * mu)
    )
    crossings = {0.0, L_back}
    for R_k in radii[:-1]:
        disc = R_k * R_k - r_start * r_start * (1.0 - mu * mu)
        if disc <= 0:
            continue
        sqrt_disc = np.sqrt(disc)
        for s in (r_start * mu + sqrt_disc, r_start * mu - sqrt_disc):
            if 1e-12 < s < L_back - 1e-12:
                crossings.add(s)
    crossings = sorted(crossings)

    segments: list[tuple[float, float, int]] = []
    for i in range(len(crossings) - 1):
        s_a, s_b = crossings[i], crossings[i + 1]
        s_mid = 0.5 * (s_a + s_b)
        r_mid_sq = (
            r_start * r_start
            - 2.0 * r_start * s_mid * mu
            + s_mid * s_mid
        )
        r_mid = np.sqrt(max(0.0, r_mid_sq))
        region_idx = _region_at_radius_oracle(r_mid, radii)
        segments.append((s_a, s_b, region_idx))
    return segments, L_back


def _chord_segments_oracle(
    h: float, R: float, radii: np.ndarray,
) -> tuple[list[tuple[float, float, int]], float]:
    """Antipodal-chord segments crossing region boundaries.

    Helper used by :class:`MultiRegionSphereChordOracle`. Bit-equal
    with the pre-R3 :func:`_chord_segments` in :mod:`.greens_function`.
    """
    if h >= R:
        return [], 0.0
    L_p = 2.0 * np.sqrt(R * R - h * h)
    crossings = {0.0, L_p}
    for R_k in radii[:-1]:
        if R_k <= h:
            continue
        sqrt_disc = np.sqrt(R_k * R_k - h * h)
        for s in (L_p / 2.0 - sqrt_disc, L_p / 2.0 + sqrt_disc):
            if 1e-12 < s < L_p - 1e-12:
                crossings.add(s)
    crossings = sorted(crossings)
    segments: list[tuple[float, float, int]] = []
    for i in range(len(crossings) - 1):
        s_a, s_b = crossings[i], crossings[i + 1]
        s_mid = 0.5 * (s_a + s_b)
        r_mid_sq = h * h + (s_mid - L_p / 2.0) ** 2
        r_mid = np.sqrt(max(0.0, r_mid_sq))
        region_idx = _region_at_radius_oracle(r_mid, radii)
        segments.append((s_a, s_b, region_idx))
    return segments, L_p


@dataclass(frozen=True)
class MultiRegionSphereChordOracle:
    r"""Chord oracle for a piecewise-homogeneous sphere with concentric
    regions.

    The MR oracle is a *decorator* over the basic sphere oracle: each
    chord (both first-leg and bounce-period) is decomposed into segments
    by region boundaries, with per-region :math:`\Sigma_{t,k}` applied
    to each segment's contribution to the optical-depth + source-line
    integrals. The shared :func:`apply_variant_alpha_closure` then sees
    the cumulative optical depths as if the medium were homogeneous.

    Bit-equality contract: the body of :meth:`apply_operator` is a
    verbatim relocation of the pre-R3 :func:`._apply_operator_mr` in
    :mod:`.greens_function`.

    Attributes
    ----------
    r_nodes : (n_r,) ndarray
        Radial Gauss-Legendre nodes on :math:`(0, R)`.
    mu_nodes : (n_mu,) ndarray
        Direction-cosine Gauss-Legendre nodes on :math:`[-1, 1]`.
    R : float
        Outer radius (= ``radii[-1]``).
    radii : (n_regions,) ndarray
        Region outer radii, ascending; ``radii[-1] = R``.
    sigma_t_per_region : (n_regions,) ndarray
        Per-region :math:`\Sigma_{t,g}` for THIS group. The signature
        differs from :class:`SphereChordOracle` (which carries
        ``sigma_t`` as a ``apply_operator`` argument) because the MR
        case must close over the per-region cross-section vector.
    alpha : float
        Surface specular reflectivity at the outermost surface
        :math:`r = R`.
    """

    r_nodes: np.ndarray
    mu_nodes: np.ndarray
    R: float
    radii: np.ndarray
    sigma_t_per_region: np.ndarray
    alpha: float

    def apply_operator(
        self,
        source_profile: np.ndarray,
        sigma_t: float,
        *,
        n_traj_quad: int,
    ) -> np.ndarray:
        r"""MR-sphere Variant α operator. See
        :meth:`ChordOracle.apply_operator`.

        Note that the ``sigma_t`` argument is **ignored** for MR — the
        per-region :math:`\Sigma_{t,k}` carried by the oracle's
        :attr:`sigma_t_per_region` is used instead. The argument is
        present for Protocol conformance.

        Bit-equal with the pre-R3 :func:`._apply_operator_mr` body.
        """
        r_nodes = self.r_nodes
        mu_nodes = self.mu_nodes
        R = self.R
        radii = self.radii
        sigma_t_per_region = self.sigma_t_per_region
        alpha = self.alpha

        source_interp = CubicSpline(r_nodes, source_profile, extrapolate=True)
        s_quad_raw, w_quad_raw = np.polynomial.legendre.leggauss(n_traj_quad)
        s_unit = 0.5 * (s_quad_raw + 1.0)
        w_unit = 0.5 * w_quad_raw

        n_r = len(r_nodes)
        n_mu = len(mu_nodes)
        psi_new = np.zeros((n_r, n_mu))

        for i in range(n_r):
            r = r_nodes[i]
            for q_idx in range(n_mu):
                mu = mu_nodes[q_idx]

                traj_segs, L_back = _trajectory_segments_oracle(
                    r, mu, R, radii,
                )
                F = 0.0
                tau_back = 0.0
                for s_a, s_b, region_idx in traj_segs:
                    sigma_t_k = float(sigma_t_per_region[region_idx])
                    seg_len = s_b - s_a
                    s_pts = s_a + s_unit * seg_len
                    tau_at_pts = tau_back + sigma_t_k * (s_pts - s_a)
                    r_traj_sq = (
                        r * r - 2.0 * r * s_pts * mu + s_pts * s_pts
                    )
                    r_traj = np.sqrt(np.clip(r_traj_sq, 0.0, R * R))
                    source_at_pts = source_interp(r_traj)
                    integrand = source_at_pts * np.exp(-tau_at_pts)
                    F += seg_len * np.sum(w_unit * integrand)
                    tau_back += sigma_t_k * seg_len

                tau_first_leg = tau_back

                if alpha == 0.0:
                    psi_new[i, q_idx] = F
                    continue

                disc = R * R - r * r * (1.0 - mu * mu)
                sqrt_disc = np.sqrt(max(disc, 0.0))
                mu_surf = sqrt_disc / R
                h = R * np.sqrt(max(0.0, 1.0 - mu_surf * mu_surf))
                chord_segs, L_p = _chord_segments_oracle(h, R, radii)

                B = 0.0
                tau_p_partial = 0.0
                for s_a, s_b, region_idx in chord_segs:
                    sigma_t_k = float(sigma_t_per_region[region_idx])
                    seg_len = s_b - s_a
                    s_pts = s_a + s_unit * seg_len
                    tau_at_pts = tau_p_partial + sigma_t_k * (s_pts - s_a)
                    r_chord_sq = h * h + (s_pts - L_p / 2.0) ** 2
                    r_chord = np.sqrt(np.clip(r_chord_sq, 0.0, R * R))
                    source_at_pts = source_interp(r_chord)
                    integrand = source_at_pts * np.exp(-tau_at_pts)
                    B += seg_len * np.sum(w_unit * integrand)
                    tau_p_partial += sigma_t_k * seg_len

                tau_p = tau_p_partial

                psi_new[i, q_idx] = apply_variant_alpha_closure(
                    F=F, B=B,
                    tau_first_leg=tau_first_leg,
                    tau_period=tau_p,
                    alpha=alpha,
                )

        return psi_new


# ─────────────────────────────────────────────────────────────────────
# Cylinder oracle — 3D angular phase-space, rank-1 closure
# ─────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class CylinderChordOracle:
    r"""Chord oracle for the homogeneous solid cylinder.

    Phase-space: :math:`(r, \mu_{\rm axial}, \varphi_{\rm az})` with
    :math:`r \in (0, R)`, :math:`\mu_{\rm axial} \in [-1, 1]`, and
    :math:`\varphi_{\rm az} \in [0, 2\pi)`. Geometry primitives:

    - In-plane backward chord :math:`L_{\rm 2D, first} = r\cos\varphi
      + \sqrt{R^2 - r^2 \sin^2\varphi}`.
    - Conserved impact parameter :math:`b = r\,|\sin\varphi_{\rm az}|`.
    - Bounce-period chord :math:`L_{\rm 2D, period} = 2\sqrt{R^2 - b^2}`.
    - 3D arclength lift :math:`s_{\rm 3D} = s_{\rm 2D}/\sqrt{1 -
      \mu_{\rm axial}^{2}}`.
    - Rank-1 closure with reflectivity :math:`\alpha`.

    Bit-equality contract: the body of :meth:`apply_operator` is a
    verbatim relocation of the pre-R3 :func:`._apply_operator_cylinder`
    in :mod:`.greens_function_cylinder`.

    Attributes
    ----------
    r_nodes : (n_r,) ndarray
    mu_axial_nodes : (n_mu,) ndarray
    phi_az_nodes : (n_phi,) ndarray
    R : float
    alpha : float
    """

    r_nodes: np.ndarray
    mu_axial_nodes: np.ndarray
    phi_az_nodes: np.ndarray
    R: float
    alpha: float

    def apply_operator(
        self,
        source_profile: np.ndarray,
        sigma_t: float,
        *,
        n_traj_quad: int,
    ) -> np.ndarray:
        r"""Cylinder Variant α operator. See
        :meth:`ChordOracle.apply_operator`.

        Bit-equal with the pre-R3 :func:`._apply_operator_cylinder` body.
        """
        r_nodes = self.r_nodes
        mu_axial_nodes = self.mu_axial_nodes
        phi_az_nodes = self.phi_az_nodes
        R = self.R
        alpha = self.alpha

        source_interp = CubicSpline(r_nodes, source_profile, extrapolate=True)

        s_quad_raw, w_quad_raw = np.polynomial.legendre.leggauss(n_traj_quad)
        s_unit = 0.5 * (s_quad_raw + 1.0)
        w_unit = 0.5 * w_quad_raw

        n_r = len(r_nodes)
        n_mu = len(mu_axial_nodes)
        n_phi = len(phi_az_nodes)
        psi_new = np.zeros((n_r, n_mu, n_phi))

        for i in range(n_r):
            r = r_nodes[i]
            for q_idx in range(n_mu):
                mu_axial = mu_axial_nodes[q_idx]
                s_in_plane = np.sqrt(max(1.0 - mu_axial * mu_axial, 1e-300))
                inv_s_in_plane = 1.0 / s_in_plane

                for p_idx in range(n_phi):
                    phi_az = phi_az_nodes[p_idx]

                    sin_phi = np.sin(phi_az)
                    cos_phi = np.cos(phi_az)
                    disc = R * R - r * r * sin_phi * sin_phi
                    L_2D_first = r * cos_phi + np.sqrt(max(disc, 0.0))
                    L_first_3D = L_2D_first * inv_s_in_plane

                    s_pts_2D = s_unit * L_2D_first
                    r_traj_sq = (
                        r * r - 2.0 * r * cos_phi * s_pts_2D
                        + s_pts_2D * s_pts_2D
                    )
                    r_traj = np.sqrt(np.clip(r_traj_sq, 0.0, R * R))
                    tau_3D = sigma_t * s_pts_2D * inv_s_in_plane
                    integrand_F = source_interp(r_traj) * np.exp(-tau_3D)
                    F = L_first_3D * np.sum(w_unit * integrand_F)

                    if alpha == 0.0:
                        psi_new[i, q_idx, p_idx] = F
                        continue

                    b = r * abs(np.sin(phi_az))
                    if b >= R:
                        L_2D_period = 0.0
                    else:
                        L_2D_period = 2.0 * np.sqrt(R * R - b * b)

                    if L_2D_period <= 0.0:
                        psi_new[i, q_idx, p_idx] = F
                        continue

                    L_period_3D = L_2D_period * inv_s_in_plane
                    s_pts_2D_p = s_unit * L_2D_period
                    shifted = s_pts_2D_p - 0.5 * L_2D_period
                    r_chord_sq = b * b + shifted * shifted
                    r_chord = np.sqrt(np.clip(r_chord_sq, 0.0, R * R))
                    tau_3D_p = sigma_t * s_pts_2D_p * inv_s_in_plane
                    integrand_B = source_interp(r_chord) * np.exp(-tau_3D_p)
                    B = L_period_3D * np.sum(w_unit * integrand_B)

                    psi_new[i, q_idx, p_idx] = apply_variant_alpha_closure(
                        F=F, B=B,
                        tau_first_leg=sigma_t * L_first_3D,
                        tau_period=sigma_t * L_period_3D,
                        alpha=alpha,
                    )

        return psi_new


# ─────────────────────────────────────────────────────────────────────
# Asymmetric slab oracle — 1D, rank-2 closure
# ─────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class SlabAsymmetricChordOracle:
    r"""Chord oracle for the homogeneous 1D slab with independent
    specular reflectivities :math:`\alpha_L, \alpha_R \in [0, 1]`.

    Phase-space: :math:`(x, \mu)` with :math:`x \in (0, L)` and signed
    :math:`\mu \in [-1, 1]`. Geometry primitives:

    - First-leg backward chord depends on sign of :math:`\mu`:
      :math:`L_{\rm first} = x/\mu` for :math:`\mu > 0` (came from
      :math:`x=0`); :math:`L_{\rm first} = (L-x)/|\mu|` for
      :math:`\mu < 0` (came from :math:`x=L`).
    - Single-transit chord :math:`L_{\rm transit} = L/|\mu|`.
    - Two single-transit B integrals: :math:`B_{LR}` (chord from
      :math:`x=0`) and :math:`B_{RL}` (chord from :math:`x=L`).
    - Rank-2 closure with surface label `'left'` (:math:`\mu > 0`)
      or `'right'` (:math:`\mu < 0`).

    The slab-symmetric solver delegates to this oracle with
    :math:`\alpha_L = \alpha_R = \alpha`; the rank-2 algebra at the
    symmetric corner is structurally equivalent to (but not bit-equal
    to) the rank-1 closure with :math:`\alpha^2 e^{-2\tau}` per-period
    reflection product.

    Bit-equality contract: verbatim relocation of pre-R3
    :func:`._apply_operator_slab_asymmetric` in
    :mod:`.greens_function_slab_asymmetric`.

    Attributes
    ----------
    x_nodes : (n_x,) ndarray
    mu_nodes : (n_mu,) ndarray
    L : float
    alpha_left : float
    alpha_right : float
    """

    x_nodes: np.ndarray
    mu_nodes: np.ndarray
    L: float
    alpha_left: float
    alpha_right: float

    def apply_operator(
        self,
        source_profile: np.ndarray,
        sigma_t: float,
        *,
        n_traj_quad: int,
    ) -> np.ndarray:
        r"""Slab-asymmetric Variant α operator. See
        :meth:`ChordOracle.apply_operator`.

        Bit-equal with the pre-R3
        :func:`._apply_operator_slab_asymmetric` body.
        """
        x_nodes = self.x_nodes
        mu_nodes = self.mu_nodes
        L = self.L
        alpha_left = self.alpha_left
        alpha_right = self.alpha_right

        source_interp = CubicSpline(x_nodes, source_profile, extrapolate=True)

        s_quad_raw, w_quad_raw = np.polynomial.legendre.leggauss(n_traj_quad)
        s_unit = 0.5 * (s_quad_raw + 1.0)
        w_unit = 0.5 * w_quad_raw

        n_x = len(x_nodes)
        n_mu = len(mu_nodes)
        psi_new = np.zeros((n_x, n_mu))

        for i in range(n_x):
            x = x_nodes[i]
            for q_idx in range(n_mu):
                mu = mu_nodes[q_idx]
                abs_mu = abs(mu)

                if mu > 0:
                    L_first = x / mu
                    s_pts_first = s_unit * L_first
                    x_traj = x - mu * s_pts_first
                    surface = "left"
                else:
                    L_first = (L - x) / abs_mu
                    s_pts_first = s_unit * L_first
                    x_traj = x + abs_mu * s_pts_first
                    surface = "right"

                x_traj = np.clip(x_traj, 0.0, L)

                L_transit = L / abs_mu
                s_pts_transit = s_unit * L_transit

                x_chord_LR = np.clip(abs_mu * s_pts_transit, 0.0, L)
                x_chord_RL = np.clip(L - abs_mu * s_pts_transit, 0.0, L)

                integrand_F = (
                    source_interp(x_traj) * np.exp(-sigma_t * s_pts_first)
                )
                F = L_first * np.sum(w_unit * integrand_F)

                if alpha_left == 0.0 and alpha_right == 0.0:
                    psi_new[i, q_idx] = F
                    continue

                integrand_B_LR = (
                    source_interp(x_chord_LR) * np.exp(-sigma_t * s_pts_transit)
                )
                B_LR = L_transit * np.sum(w_unit * integrand_B_LR)

                integrand_B_RL = (
                    source_interp(x_chord_RL) * np.exp(-sigma_t * s_pts_transit)
                )
                B_RL = L_transit * np.sum(w_unit * integrand_B_RL)

                psi_new[i, q_idx] = apply_variant_alpha_closure_rank2(
                    F=F, B_RL=B_RL, B_LR=B_LR,
                    tau_first_leg=sigma_t * L_first,
                    tau_single_transit=sigma_t * L_transit,
                    alpha_left=alpha_left,
                    alpha_right=alpha_right,
                    surface=surface,
                )

        return psi_new


# ─────────────────────────────────────────────────────────────────────
# Hollow sphere oracle — b-partition (rank-1 outer-only / rank-2 through-ray)
# ─────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class HollowSphereChordOracle:
    r"""Chord oracle for the homogeneous hollow sphere with independent
    specular reflectivities :math:`\alpha_{\rm in}, \alpha_{\rm out}`.

    Phase-space: :math:`(r, \mu)` with :math:`r \in [R_{\rm in}, R_{\rm
    out}]`. The conserved impact parameter :math:`b = r\sqrt{1-\mu^2}`
    partitions phase space:

    - **Outer-only branch** (:math:`b > R_{\rm in}`): trajectory does
      NOT cross the inner cavity. Algebra identical to a solid sphere
      of radius :math:`R_{\rm out}` — rank-1 closure with
      :math:`\alpha_{\rm out}`.
    - **Through-ray branch** (:math:`b \le R_{\rm in}`): trajectory
      crosses the inner cavity. Single-transit shell-traversal length
      :math:`L_{\rm step} = \sqrt{R_{\rm out}^2 - b^2} - \sqrt{R_{\rm
      in}^2 - b^2}`. Rank-2 closure with surface routing: first
      arrival is INNER if :math:`\mu > 0` (backward goes inward),
      OUTER otherwise.

    Bit-equality contract: verbatim relocation of pre-R3
    :func:`._apply_operator_hollow_sphere` in
    :mod:`.greens_function_hollow_sphere`.

    Attributes
    ----------
    r_nodes : (n_r,) ndarray
    mu_nodes : (n_mu,) ndarray
    R_in : float
    R_out : float
    alpha_in : float
    alpha_out : float
    """

    r_nodes: np.ndarray
    mu_nodes: np.ndarray
    R_in: float
    R_out: float
    alpha_in: float
    alpha_out: float

    def apply_operator(
        self,
        source_profile: np.ndarray,
        sigma_t: float,
        *,
        n_traj_quad: int,
    ) -> np.ndarray:
        r"""Hollow sphere Variant α operator. See
        :meth:`ChordOracle.apply_operator`.

        Bit-equal with the pre-R3
        :func:`._apply_operator_hollow_sphere` body.
        """
        r_nodes = self.r_nodes
        mu_nodes = self.mu_nodes
        R_in = self.R_in
        R_out = self.R_out
        alpha_in = self.alpha_in
        alpha_out = self.alpha_out

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

                b_sq = r * r * (1.0 - mu * mu)
                b = np.sqrt(max(b_sq, 0.0))

                disc_out = R_out_sq - b_sq
                sqrt_disc_out = np.sqrt(max(disc_out, 0.0))

                if b > R_in:
                    L_first = r * mu + sqrt_disc_out
                    L_p = 2.0 * sqrt_disc_out

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

                disc_in = R_in_sq - b_sq
                sqrt_disc_in = np.sqrt(max(disc_in, 0.0))

                L_step = sqrt_disc_out - sqrt_disc_in

                if mu > 0:
                    L_first = r * mu - sqrt_disc_in
                    surface = "inner"
                else:
                    L_first = r * mu + sqrt_disc_out
                    surface = "outer"

                L_first = max(L_first, 0.0)

                s_pts = s_unit * L_first
                r_traj_sq = r * r - 2.0 * r * mu * s_pts + s_pts * s_pts
                r_traj = np.sqrt(np.clip(r_traj_sq, R_in_sq, R_out_sq))
                integrand_F = (
                    source_interp(r_traj) * np.exp(-sigma_t * s_pts)
                )
                F = L_first * np.sum(w_unit * integrand_F)

                if alpha_in == 0.0 and alpha_out == 0.0:
                    psi_new[i, q_idx] = F
                    continue

                s_pts_step = s_unit * L_step

                r_chord_out_sq = b_sq + (sqrt_disc_in + s_pts_step) ** 2
                r_chord_out = np.sqrt(
                    np.clip(r_chord_out_sq, R_in_sq, R_out_sq)
                )
                integrand_B_out = (
                    source_interp(r_chord_out)
                    * np.exp(-sigma_t * s_pts_step)
                )
                B_out = L_step * np.sum(w_unit * integrand_B_out)

                r_chord_in_sq = b_sq + (sqrt_disc_out - s_pts_step) ** 2
                r_chord_in = np.sqrt(
                    np.clip(r_chord_in_sq, R_in_sq, R_out_sq)
                )
                integrand_B_in = (
                    source_interp(r_chord_in)
                    * np.exp(-sigma_t * s_pts_step)
                )
                B_in = L_step * np.sum(w_unit * integrand_B_in)

                psi_new[i, q_idx] = apply_variant_alpha_closure_rank2(
                    F=F, B_RL=B_in, B_LR=B_out,
                    tau_first_leg=sigma_t * L_first,
                    tau_single_transit=sigma_t * L_step,
                    alpha_left=alpha_in,
                    alpha_right=alpha_out,
                    surface=("left" if surface == "inner" else "right"),
                )

        return psi_new


# ─────────────────────────────────────────────────────────────────────
# Annulus oracle — cylinder + b-partition + 3D angular phase-space
# ─────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class AnnulusChordOracle:
    r"""Chord oracle for the homogeneous annulus (hollow cylinder shell)
    with independent specular reflectivities :math:`\alpha_{\rm in},
    \alpha_{\rm out}`.

    Phase-space: :math:`(r, \mu_{\rm axial}, \varphi_{\rm az})` with
    :math:`r \in [R_{\rm in}, R_{\rm out}]`. Combines:

    - The cylinder's 3D angular phase-space (:math:`\mu_{\rm axial} \in
      [-1,1]` and :math:`\varphi_{\rm az} \in [0, 2\pi)`), with the
      same 3D arclength lift :math:`s_{\rm 3D} = s_{\rm 2D}/\sqrt{1 -
      \mu_{\rm axial}^2}`.
    - The hollow sphere's :math:`b`-partition routing (:math:`b > R_{\rm
      in}` → outer-only rank-1; :math:`b \le R_{\rm in}` → through-ray
      rank-2). Here :math:`b = r\,|\sin\varphi_{\rm az}|` is the in-
      plane impact parameter (cylinder convention).

    First-arrival surface for through-ray branch is INNER if
    :math:`\cos\varphi_{\rm az} > 0` (backward goes inward), OUTER
    otherwise — the cylinder analog of hollow_sphere's :math:`\mu>0`
    branch.

    Bit-equality contract: verbatim relocation of pre-R3
    :func:`._apply_operator_annulus` in :mod:`.greens_function_annulus`.

    Attributes
    ----------
    r_nodes : (n_r,) ndarray
    mu_axial_nodes : (n_mu,) ndarray
    phi_az_nodes : (n_phi,) ndarray
    R_in : float
    R_out : float
    alpha_in : float
    alpha_out : float
    """

    r_nodes: np.ndarray
    mu_axial_nodes: np.ndarray
    phi_az_nodes: np.ndarray
    R_in: float
    R_out: float
    alpha_in: float
    alpha_out: float

    def apply_operator(
        self,
        source_profile: np.ndarray,
        sigma_t: float,
        *,
        n_traj_quad: int,
    ) -> np.ndarray:
        r"""Annulus Variant α operator. See
        :meth:`ChordOracle.apply_operator`.

        Bit-equal with the pre-R3 :func:`._apply_operator_annulus` body.
        """
        r_nodes = self.r_nodes
        mu_axial_nodes = self.mu_axial_nodes
        phi_az_nodes = self.phi_az_nodes
        R_in = self.R_in
        R_out = self.R_out
        alpha_in = self.alpha_in
        alpha_out = self.alpha_out

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
                s_in_plane = np.sqrt(max(1.0 - mu_axial * mu_axial, 1e-300))
                inv_s_in_plane = 1.0 / s_in_plane

                for p_idx in range(n_phi):
                    phi_az = phi_az_nodes[p_idx]
                    cos_phi = np.cos(phi_az)
                    sin_phi = np.sin(phi_az)

                    b = r * abs(sin_phi)
                    b_sq = b * b

                    disc_out = R_out_sq - b_sq
                    sqrt_disc_out = np.sqrt(max(disc_out, 0.0))

                    if b > R_in:
                        L_2D_first = r * cos_phi + sqrt_disc_out
                        L_first_3D = L_2D_first * inv_s_in_plane

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

                    disc_in = R_in_sq - b_sq
                    sqrt_disc_in = np.sqrt(max(disc_in, 0.0))

                    L_step_2D = sqrt_disc_out - sqrt_disc_in
                    L_step_3D = L_step_2D * inv_s_in_plane

                    if cos_phi > 0.0:
                        L_2D_first = r * cos_phi - sqrt_disc_in
                        surface = "inner"
                    else:
                        L_2D_first = r * cos_phi + sqrt_disc_out
                        surface = "outer"

                    L_2D_first = max(L_2D_first, 0.0)
                    L_first_3D = L_2D_first * inv_s_in_plane

                    s_pts_2D = s_unit * L_2D_first
                    r_traj_sq = (
                        r * r - 2.0 * r * cos_phi * s_pts_2D
                        + s_pts_2D * s_pts_2D
                    )
                    r_traj = np.sqrt(np.clip(r_traj_sq, R_in_sq, R_out_sq))
                    tau_3D = sigma_t * s_pts_2D * inv_s_in_plane
                    integrand_F = source_interp(r_traj) * np.exp(-tau_3D)
                    F = L_first_3D * np.sum(w_unit * integrand_F)

                    if alpha_in == 0.0 and alpha_out == 0.0:
                        psi_new[i, q_idx, p_idx] = F
                        continue

                    s_pts_2D_step = s_unit * L_step_2D

                    r_chord_out_sq = b_sq + (sqrt_disc_in + s_pts_2D_step) ** 2
                    r_chord_out = np.sqrt(
                        np.clip(r_chord_out_sq, R_in_sq, R_out_sq)
                    )
                    tau_3D_step = sigma_t * s_pts_2D_step * inv_s_in_plane
                    integrand_B_out = (
                        source_interp(r_chord_out) * np.exp(-tau_3D_step)
                    )
                    B_out = L_step_3D * np.sum(w_unit * integrand_B_out)

                    r_chord_in_sq = b_sq + (sqrt_disc_out - s_pts_2D_step) ** 2
                    r_chord_in = np.sqrt(
                        np.clip(r_chord_in_sq, R_in_sq, R_out_sq)
                    )
                    integrand_B_in = (
                        source_interp(r_chord_in) * np.exp(-tau_3D_step)
                    )
                    B_in = L_step_3D * np.sum(w_unit * integrand_B_in)

                    psi_new[i, q_idx, p_idx] = apply_variant_alpha_closure_rank2(
                        F=F, B_RL=B_in, B_LR=B_out,
                        tau_first_leg=sigma_t * L_first_3D,
                        tau_single_transit=sigma_t * L_step_3D,
                        alpha_left=alpha_in,
                        alpha_right=alpha_out,
                        surface=("left" if surface == "inner" else "right"),
                    )

        return psi_new
