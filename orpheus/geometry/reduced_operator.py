r"""Reduced streaming operator — connection coefficients on the coordinate chart.

This module hosts the **connection-coefficient primitive** shared by SN,
MoC, and CP curvilinear sweeps.  In differential-geometric language,
the spherical :math:`(1-\mu^2)/r\,\partial_\mu` and cylindrical
:math:`-(1/r)\,\partial_\varphi(\xi\cdot)` redistribution terms are
**the same connection-coefficient operator** on SO(3) viewed in two
coordinate charts (polar-on-sphere vs. azimuthal-on-cylinder).  Each
solver — SN, MoC, CP — needs the same numerical data: chord lengths,
face areas, the :math:`\Delta A/w` geometry factor, the Bailey 2009
:math:`\alpha` dome recursion, and the Morel--Montry angular closure
weights :math:`\tau_{mm}`.

Per Cardinal Rule 2 (architecture), this primitive **MUST NOT** be
duplicated across solvers.  The historical home was
:class:`orpheus.sn.geometry.SNMesh._setup_spherical` and
:class:`~orpheus.sn.geometry.SNMesh._setup_cylindrical`; this module is
the geometry-layer parallel implementation that those consumers will
migrate to in Wave G of the SN reshape campaign (Issue 10).

The output is **bit-identical** to the legacy SNMesh setup methods
(verified by hash-equality assertions in
``tests/geometry/test_reduced_operator.py``).  Bit-identity is the
load-bearing contract — both call sites compute connection coefficients
the same way; this module simply lifts the math to a place where MoC
and CP can share it.

Mathematical content
====================

The Bailey 2009 dome recursion (sphere):

.. math::
   :label: bailey-dome-recursion

   \alpha_{n+\tfrac12} = \alpha_{n-\tfrac12} - w_n\,\mu_n,
   \qquad \alpha_{\tfrac12} = 0.

For Gauss--Legendre quadrature with :math:`\mu` sorted ascending in
:math:`[-1, 1]`, the recursion produces a non-negative dome that closes
back to zero at the upper boundary by GL antisymmetry.

The cylindrical analog is per-:math:`\mu`-level: each level :math:`p`
carries its own :math:`(M+1)`-tuple of :math:`\alpha` with
:math:`\alpha^{(p)}_{m+1/2} = \alpha^{(p)}_{m-1/2} - w_m\,\eta_m`,
where :math:`\eta` is the radial direction cosine.

The Morel--Montry angular closure
(Lewis & Miller 1984, §4.5; Bailey et al. 2009 Eq. 74):

.. math::
   :label: morel-montry-clamp

   \tau_m = \mathrm{clip}\!\left(
       \frac{\mu_m - \mu_{m-1/2}}{\mu_{m+1/2} - \mu_{m-1/2}},
       \;\tfrac12,\; 1
   \right),

with the :math:`[1/2, 1]` clamp keeping the M-M weighting positive.

References
==========

* Bailey, T. S., Adams, M. L., Yang, B., & Zika, M. R. (2009).
  *A piecewise linear finite element discretization of the diffusion
  equation for arbitrary polyhedral grids*. JCP 227, 3738–3757.
  Equations 50 (sphere/cylinder dome recursion) and 74 (Morel--Montry
  closure).
* Lewis, E. E., & Miller, W. F. (1984). *Computational Methods of
  Neutron Transport*.  §4.5 — angular redistribution closure;
  the :math:`\tau_{mm}` clamp keeps the closure stable.

See also
========

* :class:`orpheus.sn.geometry.SNMesh` — the legacy in-line implementation.
  Issue 10 (Wave G) of the SN reshape campaign refactors :class:`SNMesh`
  to consume this primitive.  Issues 11/12 (Wave H) make
  ``SNStreamingOperator.apply`` consume it.  MoC and CP campaigns
  (post-Wave-1) reuse this primitive with their own consumption
  patterns.
* :doc:`/theory/structured_geometry` — the architecture page;
  see "Connection coefficients (reduced streaming operator)".
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal, Protocol

import numpy as np

from .coord import CoordSystem
from .mesh import Mesh1D


# ═══════════════════════════════════════════════════════════════════════
# Angular-measure adapter (Issue 4 contract — local, structural)
# ═══════════════════════════════════════════════════════════════════════

class AngularMeasure(Protocol):
    """Minimal contract this module needs from an angular quadrature.

    The full angular-measure / discrete-measure type is owned by
    Wave A (:mod:`orpheus.sn.discrete_measure`); here we only depend on
    the attributes the connection-coefficient math reads.  This avoids
    importing ``AngularQuadrature`` from :mod:`orpheus.sn.quadrature`
    so the geometry layer remains free of SN-specific imports.
    """

    mu_x: np.ndarray
    """(N,) primary direction cosine — :math:`\\mu` for sphere,
    :math:`\\eta` (radial) for cylindrical product quadrature."""

    weights: np.ndarray
    """(N,) quadrature weights."""

    N: int
    """Number of ordinates."""


# ═══════════════════════════════════════════════════════════════════════
# StreamingTerms — geometry-dependent shape per (cell, direction)
# ═══════════════════════════════════════════════════════════════════════

@dataclass(frozen=True, slots=True)
class StreamingTerms:
    """Per-(cell, direction) numerical inputs for a sweep cell update.

    The set of populated fields is **geometry-dependent**:

    * **Slab**: ``chord_length``, ``mu``, ``volume``, ``abs_mu``.
      No curvature math.
    * **Sphere**: ``chord_length``, ``face_area_in``, ``face_area_out``,
      ``delta_A_over_w``, ``alpha_in``, ``alpha_out``, ``tau_mm``,
      ``volume``, ``abs_mu``.
    * **Cylinder**: same shape as sphere, with the per-level
      :math:`\\alpha` and :math:`\\tau_{mm}` indexed via
      ``mu_level_idx`` at the call site.

    Fields not used in a given geometry retain their default of
    ``None`` (object-typed) so a downstream ``StreamingTerms`` can be
    inspected with ``is None`` checks without raising.

    The trailing ``volume`` and ``abs_mu`` fields are populated by
    **all three factories** so that a downstream
    :class:`~orpheus.sn.spatial.cell_update.CellUpdate` strategy
    receives a self-contained per-cell, per-direction packet and need
    not reach back into ``SNMesh`` or the ``AngularQuadrature``.  The
    ``alpha_in is None`` test discriminates slab from curvilinear
    geometry inside cell-update strategies.
    """

    chord_length: float
    """Cell radial width (slab/sphere/cylinder all use ``mesh.widths[i]``)."""

    mu: float | None = None
    """Direction cosine in the primary axis (slab uses this directly).

    Signed; slab DD's flow direction reads off the sign of this field.
    Cell-update strategies that need only the magnitude use
    ``abs_mu`` instead.
    """

    face_area_in: float | None = None
    """:math:`A_{i-1/2}` — area of the upstream radial face."""

    face_area_out: float | None = None
    """:math:`A_{i+1/2}` — area of the downstream radial face."""

    delta_A_over_w: float | None = None
    """:math:`\\Delta A_i / w_n` — the geometry-redistribution factor."""

    alpha_in: float | None = None
    """:math:`\\alpha_{n-1/2}` — incoming half-angle redistribution.

    ``None`` for slab; populated for sphere/cylinder.  Cell-update
    strategies use ``alpha_in is None`` as the slab vs curvilinear
    discriminator.
    """

    alpha_out: float | None = None
    """:math:`\\alpha_{n+1/2}` — outgoing half-angle redistribution."""

    tau_mm: float | None = None
    """Morel--Montry angular closure weight on this ordinate."""

    volume: float | None = None
    """Cell volume :math:`V_i`.

    Populated by all three factories from
    :attr:`~orpheus.geometry.mesh.Mesh1D.volumes`: slab uses
    ``mesh.volumes[i]`` (which equals ``widths[i]`` for unit
    cross-section in 1-D Cartesian); sphere uses
    :math:`\\tfrac{4}{3}\\pi(r_{i+1}^3 - r_i^3)`; cylinder uses
    :math:`\\pi(r_{i+1}^2 - r_i^2)` (per unit axial length).  Carried
    on the streaming-terms packet so that the
    :class:`~orpheus.sn.spatial.cell_update.CellUpdate` cell-update
    contract receives :math:`V_i` directly without needing access to
    the underlying ``SNMesh``.
    """

    abs_mu: float | None = None
    """Absolute primary direction cosine.

    :math:`|\\mu|` for slab and sphere (axial direction cosine);
    :math:`|\\eta|` for cylindrical 1-D radial sweeps (the radial
    direction cosine).  All three factories compute this as
    ``abs(quadrature.mu_x[direction_idx])`` — by ORPHEUS convention,
    ``mu_x`` carries :math:`\\mu` for sphere and :math:`\\eta` for
    cylinder, and both are real-valued, so the absolute value is
    always well-defined.

    The cylindrical pure-azimuthal degenerate case (``abs_mu <
    1e-15``) is signalled by cell-update strategies via
    :attr:`~orpheus.sn.spatial.cell_update.CellResult.outgoing_spatial_flux`
    set to ``None`` (no radial face flow on the cell).
    """


# ═══════════════════════════════════════════════════════════════════════
# ReducedStreamingOperator — the primitive
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class ReducedStreamingOperator:
    """Connection-coefficient data for a curvilinear (or slab) sweep.

    Carries all geometry- and quadrature-dependent precomputed arrays
    that an SN, MoC, or CP sweep needs to march through the angular
    redistribution.  In Cardinal Rule 2 framing, this primitive holds
    the data the **same connection-coefficient operator** demands, in
    whichever coordinate chart the consumer's geometry lives in.

    Three factory functions construct concrete instances:

    * :func:`slab_streaming` — Cartesian 1-D, no curvature.
    * :func:`spherical_streaming` — 1-D spherical, with
      :math:`\\alpha_{n+1/2}` dome recursion (Bailey 2009 Eq. 50).
    * :func:`cylindrical_streaming` — 1-D cylindrical, with
      per-:math:`\\mu`-level :math:`\\alpha` and :math:`\\tau_{mm}`
      structures.

    Output shape contract (bit-identical to ``SNMesh._setup_*``)::

        slab:        face_areas == None,
                     alpha_half == None,
                     redist_dAw == None,
                     tau_mm == None.
        sphere:      face_areas (nx+1,),
                     delta_A    (nx,),
                     alpha_half (N+1,),
                     redist_dAw (nx, N),
                     tau_mm     (N,).
        cylinder:    face_areas (nx+1,),
                     delta_A    (nx,),
                     alpha_per_level     [(M+1,)] · n_levels,
                     redist_dAw_per_level [(nx, M)] · n_levels,
                     tau_mm_per_level     [(M,)] · n_levels.

    Attributes
    ----------
    coord :
        Coordinate system tag (CARTESIAN / SPHERICAL / CYLINDRICAL).
    mesh :
        The :class:`~orpheus.geometry.mesh.Mesh1D` this operator was
        built from.  Held by reference; not copied.
    requires_upstream_angular_state :
        ``True`` for sphere/cylinder (the sweep needs ψ at the
        upstream half-angle to apply the redistribution); ``False`` for
        slab.
    angular_marching_axis :
        ``"mu"`` for sphere/cylinder (the sweep marches over the μ
        ordinate index per radial cell); ``None`` for slab.
    """

    coord: CoordSystem
    mesh: Mesh1D
    requires_upstream_angular_state: bool
    angular_marching_axis: Literal["mu"] | None

    # Common (curvilinear)
    face_areas: np.ndarray | None = None
    delta_A: np.ndarray | None = None

    # Spherical
    alpha_half: np.ndarray | None = None
    redist_dAw: np.ndarray | None = None
    tau_mm: np.ndarray | None = None

    # Cylindrical (per-level)
    alpha_per_level: list[np.ndarray] | None = None
    redist_dAw_per_level: list[np.ndarray] | None = None
    tau_mm_per_level: list[np.ndarray] | None = None

    # Quadrature reference (kept for streaming_terms() extraction)
    _quadrature: AngularMeasure | None = field(default=None, repr=False)

    # ── Per-direction extraction ───────────────────────────────────

    def streaming_terms(
        self,
        cell_idx: int,
        direction_idx: int,
        mu_level_idx: int | None = None,
    ) -> StreamingTerms:
        """Pack the per-cell, per-direction streaming inputs.

        For cylindrical geometry, ``mu_level_idx`` selects which
        :math:`\\mu`-level the ``direction_idx`` belongs to; the
        ``direction_idx`` is then the within-level azimuthal index
        :math:`m \\in [0, M)`.  For slab and sphere, ``mu_level_idx``
        is ignored.
        """
        chord = float(self.mesh.widths[cell_idx])
        # Common per-(cell, direction) volume + |μ| (or |η|).
        # mesh.volumes returns shape (N,) for 1-D meshes; abs() of
        # the primary direction cosine is well-defined on the real line.
        assert self._quadrature is not None
        volume = float(self.mesh.volumes[cell_idx])
        abs_mu = float(abs(self._quadrature.mu_x[direction_idx]))

        if self.coord is CoordSystem.CARTESIAN:
            # Slab — minimal content; mu carried for the consumer.
            return StreamingTerms(
                chord_length=chord,
                mu=float(self._quadrature.mu_x[direction_idx]),
                volume=volume,
                abs_mu=abs_mu,
            )

        if self.coord is CoordSystem.SPHERICAL:
            assert self.face_areas is not None
            assert self.alpha_half is not None
            assert self.redist_dAw is not None
            assert self.tau_mm is not None
            return StreamingTerms(
                chord_length=chord,
                face_area_in=float(self.face_areas[cell_idx]),
                face_area_out=float(self.face_areas[cell_idx + 1]),
                delta_A_over_w=float(
                    self.redist_dAw[cell_idx, direction_idx]
                ),
                alpha_in=float(self.alpha_half[direction_idx]),
                alpha_out=float(self.alpha_half[direction_idx + 1]),
                tau_mm=float(self.tau_mm[direction_idx]),
                volume=volume,
                abs_mu=abs_mu,
            )

        if self.coord is CoordSystem.CYLINDRICAL:
            if mu_level_idx is None:
                raise ValueError(
                    "cylindrical streaming_terms() requires mu_level_idx "
                    "(which μ-level the direction_idx belongs to)."
                )
            assert self.face_areas is not None
            assert self.alpha_per_level is not None
            assert self.redist_dAw_per_level is not None
            assert self.tau_mm_per_level is not None
            alpha_lv = self.alpha_per_level[mu_level_idx]
            dAw_lv = self.redist_dAw_per_level[mu_level_idx]
            tau_lv = self.tau_mm_per_level[mu_level_idx]
            return StreamingTerms(
                chord_length=chord,
                face_area_in=float(self.face_areas[cell_idx]),
                face_area_out=float(self.face_areas[cell_idx + 1]),
                delta_A_over_w=float(dAw_lv[cell_idx, direction_idx]),
                alpha_in=float(alpha_lv[direction_idx]),
                alpha_out=float(alpha_lv[direction_idx + 1]),
                tau_mm=float(tau_lv[direction_idx]),
                volume=volume,
                abs_mu=abs_mu,
            )

        raise ValueError(  # pragma: no cover — exhaustive match above
            f"Unknown coord system: {self.coord!r}"
        )


# ═══════════════════════════════════════════════════════════════════════
# Factory: slab
# ═══════════════════════════════════════════════════════════════════════

def slab_streaming(
    mesh: Mesh1D,
    angular_measure: AngularMeasure,
) -> ReducedStreamingOperator:
    """Build the slab :class:`ReducedStreamingOperator`.

    Slab geometry has no curvature — the connection coefficients
    vanish.  All ``alpha_*``, ``redist_dAw``, ``tau_mm`` arrays remain
    ``None``; the operator advertises
    ``requires_upstream_angular_state = False`` and
    ``angular_marching_axis = None``.

    Parameters
    ----------
    mesh :
        Slab :class:`~orpheus.geometry.mesh.Mesh1D`
        (``coord == CoordSystem.CARTESIAN``).
    angular_measure :
        Any :class:`AngularMeasure`-shaped object; only ``mu_x`` and
        ``weights`` are read at extraction time.
    """
    if mesh.coord is not CoordSystem.CARTESIAN:
        raise ValueError(
            f"slab_streaming requires CARTESIAN mesh, got {mesh.coord!r}"
        )
    return ReducedStreamingOperator(
        coord=CoordSystem.CARTESIAN,
        mesh=mesh,
        requires_upstream_angular_state=False,
        angular_marching_axis=None,
        _quadrature=angular_measure,
    )


# ═══════════════════════════════════════════════════════════════════════
# Factory: spherical
# ═══════════════════════════════════════════════════════════════════════

def spherical_streaming(
    mesh: Mesh1D,
    angular_measure: AngularMeasure,
) -> ReducedStreamingOperator:
    r"""Build the spherical :class:`ReducedStreamingOperator`.

    Implements Bailey et al. (2009) Eq. 50 (dome recursion) +
    Eq. 74 (Morel--Montry closure) verbatim, producing arrays
    bit-identical to :class:`SNMesh._setup_spherical`.

    The :math:`\Delta A/w` geometry factor (Cardinal Rule 2 — the
    connection-coefficient data, common to SN/MoC/CP) is precomputed
    here once.

    Parameters
    ----------
    mesh :
        Spherical :class:`~orpheus.geometry.mesh.Mesh1D`
        (``coord == CoordSystem.SPHERICAL``).
    angular_measure :
        Quadrature with ``mu_x``, ``weights``, ``N``.
    """
    if mesh.coord is not CoordSystem.SPHERICAL:
        raise ValueError(
            f"spherical_streaming requires SPHERICAL mesh, got {mesh.coord!r}"
        )

    mu = angular_measure.mu_x
    w = angular_measure.weights
    N = angular_measure.N

    # Cell face areas: A_{i+1/2} = 4πr² at each edge — sourced from the
    # mesh, which routes through coord.compute_surfaces_1d().
    face_areas = mesh.surfaces  # (nx+1,)

    # Cell face-area differences: ΔA_i = A_{i+1/2} − A_{i-1/2}
    delta_A = face_areas[1:] - face_areas[:-1]

    # Bailey et al. (2009) Eq. 50 dome recursion
    # α_{n+1/2} = α_{n-1/2} − w_n μ_n
    # For GL quadrature, this gives a non-negative dome closing to 0
    # at α_{N+1/2} by GL antisymmetry.
    alpha = np.zeros(N + 1)
    for n in range(N):
        alpha[n + 1] = alpha[n] - w[n] * mu[n]

    # Verify GL antisymmetry: α_{N+1/2} ≈ 0
    assert abs(alpha[N]) < 1e-12, (
        f"GL antisymmetry violated: α_{{N+1/2}} = {alpha[N]:.2e}"
    )

    # ΔA_i / w_n — the geometry redistribution factor (nx, N).
    redist_dAw = delta_A[:, None] / w[None, :]

    # Morel–Montry closure (Bailey Eq. 74) with [1/2, 1] clamp
    # keeping the M-M weighting positive (Lewis & Miller §4.5).
    # μ_{1/2} = -1, μ_{N+1/2} = +1.
    mu_edge = np.zeros(N + 1)
    mu_edge[0] = -1.0
    for n in range(N):
        mu_edge[n + 1] = mu_edge[n] + w[n]
    tau_mm = np.empty(N)
    for n in range(N):
        dmu = mu_edge[n + 1] - mu_edge[n]
        tau_raw = (mu[n] - mu_edge[n]) / dmu if abs(dmu) > 1e-15 else 0.5
        tau_mm[n] = max(0.5, min(1.0, tau_raw))

    return ReducedStreamingOperator(
        coord=CoordSystem.SPHERICAL,
        mesh=mesh,
        requires_upstream_angular_state=True,
        angular_marching_axis="mu",
        face_areas=face_areas,
        delta_A=delta_A,
        alpha_half=alpha,
        redist_dAw=redist_dAw,
        tau_mm=tau_mm,
        _quadrature=angular_measure,
    )


# ═══════════════════════════════════════════════════════════════════════
# Factory: cylindrical
# ═══════════════════════════════════════════════════════════════════════

def cylindrical_streaming(
    mesh: Mesh1D,
    angular_measure: AngularMeasure,
) -> ReducedStreamingOperator:
    r"""Build the cylindrical :class:`ReducedStreamingOperator`.

    Cylindrical SN sweeps march in azimuth on each :math:`\mu`-level
    independently; the connection coefficients are therefore stored
    per-level.

    Requires ``angular_measure`` to have a ``level_indices`` attribute
    (a list of index arrays, one per :math:`\mu`-level).
    :class:`~orpheus.sn.quadrature.LevelSymmetricSN` and
    :class:`~orpheus.sn.quadrature.ProductQuadrature` provide this.
    The quadrature must also expose ``mu_z`` (axial direction cosine,
    used to derive :math:`\sin\theta` for the level's azimuthal extent).

    Output is bit-identical to
    :class:`SNMesh._setup_cylindrical`.

    Parameters
    ----------
    mesh :
        Cylindrical :class:`~orpheus.geometry.mesh.Mesh1D`
        (``coord == CoordSystem.CYLINDRICAL``).
    angular_measure :
        Quadrature with ``mu_x`` (radial η), ``mu_z`` (axial),
        ``weights``, ``level_indices``.
    """
    if mesh.coord is not CoordSystem.CYLINDRICAL:
        raise ValueError(
            f"cylindrical_streaming requires CYLINDRICAL mesh, "
            f"got {mesh.coord!r}"
        )
    if not hasattr(angular_measure, "level_indices"):
        raise ValueError(
            "Cylindrical streaming requires a quadrature with level "
            "structure (LevelSymmetricSN or ProductQuadrature), "
            f"got {type(angular_measure).__name__}"
        )

    face_areas = mesh.surfaces  # (nx+1,)
    delta_A = face_areas[1:] - face_areas[:-1]

    # Per-level azimuthal redistribution coefficients
    # Bailey et al. (2009) Eq. 50: α_{m+1/2} = α_{m-1/2} − w_m · η_m
    # Ordinates are ordered by increasing η within each level.
    alpha_per_level: list[np.ndarray] = []
    for level_idx in angular_measure.level_indices:
        eta = angular_measure.mu_x[level_idx]
        w = angular_measure.weights[level_idx]
        M = len(level_idx)
        alpha = np.zeros(M + 1)
        for m in range(M):
            alpha[m + 1] = alpha[m] - w[m] * eta[m]
        alpha_per_level.append(alpha)

    # Per-level ΔA_i / w_m — same factor structure as sphere, but the
    # weights are the level's azimuthal weights.
    redist_dAw_per_level: list[np.ndarray] = []
    for level_idx in angular_measure.level_indices:
        w_level = angular_measure.weights[level_idx]
        redist_dAw_per_level.append(
            delta_A[:, None] / w_level[None, :]  # (nx, M)
        )

    # Morel–Montry closure (Bailey Eq. 74) per μ-level.
    # For cylindrical, ordinates are η-sorted but weights come from
    # φ-space (not η-space), so the weight-sum edge approach is wrong.
    # Instead, cell edges are at midpoints of consecutive η values
    # with endpoints at ±sin θ.  This gives a proper η-partition.
    #
    # (See SNMesh._setup_cylindrical for the same construction;
    # GitHub Issue #3 tracks the φ-based edge refinement.)
    mu_z = angular_measure.mu_z  # type: ignore[attr-defined]
    tau_mm_per_level: list[np.ndarray] = []
    for level_idx in angular_measure.level_indices:
        eta = angular_measure.mu_x[level_idx]
        M = len(level_idx)
        sin_theta = np.sqrt(1.0 - mu_z[level_idx[0]] ** 2)
        eta_edge = np.zeros(M + 1)
        eta_edge[0] = -sin_theta
        for m in range(M - 1):
            eta_edge[m + 1] = 0.5 * (eta[m] + eta[m + 1])
        eta_edge[M] = sin_theta
        tau = np.empty(M)
        for m in range(M):
            deta = eta_edge[m + 1] - eta_edge[m]
            tau_raw = (
                (eta[m] - eta_edge[m]) / deta if abs(deta) > 1e-15 else 0.5
            )
            tau[m] = max(0.5, min(1.0, tau_raw))
        tau_mm_per_level.append(tau)

    return ReducedStreamingOperator(
        coord=CoordSystem.CYLINDRICAL,
        mesh=mesh,
        requires_upstream_angular_state=True,
        angular_marching_axis="mu",
        face_areas=face_areas,
        delta_A=delta_A,
        alpha_per_level=alpha_per_level,
        redist_dAw_per_level=redist_dAw_per_level,
        tau_mm_per_level=tau_mm_per_level,
        _quadrature=angular_measure,
    )


__all__ = [
    "AngularMeasure",
    "ReducedStreamingOperator",
    "StreamingTerms",
    "cylindrical_streaming",
    "slab_streaming",
    "spherical_streaming",
]
