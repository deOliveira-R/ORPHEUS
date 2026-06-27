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
:class:`orpheus.sn.mesh.augmented_mesh.SNMesh._setup_spherical` and
:class:`~orpheus.sn.mesh.augmented_mesh.SNMesh._setup_cylindrical`; this module is
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

The :math:`\alpha`-dome recursion (sphere) — Hébert (2009)
§3.9.4 Eqs. 3.423-3.424, in the ORPHEUS factor-of-2-absorbed
normalization:

.. math::
   :label: bailey-dome-recursion

   \alpha_{n+\tfrac12} = \alpha_{n-\tfrac12} - w_n\,\mu_n,
   \qquad \alpha_{\tfrac12} = 0.

(Hébert's Eq. 3.424 reads
:math:`\alpha^{H}_{n+1/2} = \alpha^{H}_{n-1/2} - 2\mathcal{W}_n\,\mu_n`
with the corresponding redistribution divisor :math:`\Delta S_i /
(2\,\mathcal{W}_n)` in Eq. 3.428; the ORPHEUS arrays carry
:math:`\alpha^{O} = \alpha^{H}/2` so the factor of 2 is absorbed
into the recurrence and the redistribution divisor reads
:math:`\Delta A_i / w_n`.  Both forms are mathematically equivalent.)

For Gauss--Legendre quadrature with :math:`\mu` sorted ascending in
:math:`[-1, 1]`, the recursion produces a non-negative dome that closes
back to zero at the upper boundary by GL antisymmetry.

The cylindrical analog is per-:math:`\mu`-level: each level :math:`p`
carries its own :math:`(M+1)`-tuple of :math:`\alpha` with
:math:`\alpha^{(p)}_{m+1/2} = \alpha^{(p)}_{m-1/2} - w_m\,\eta_m`,
where :math:`\eta` is the radial direction cosine.

The Morel--Montry angular closure weight (Bailey-Morel-Chang 2010
Eq. 43) is the *fractional position* of the ordinate :math:`\mu_m` in
its half-angle interval :math:`[\mu_{m-1/2}, \mu_{m+1/2}]`:

.. math::

   \tau_m = \frac{\mu_m - \mu_{m-1/2}}{\mu_{m+1/2} - \mu_{m-1/2}}.

This raw weight is the UNIQUE choice exact for a flux linear in
:math:`\mu` (Bailey-Morel-Chang 2010 Eq. 43), with admissible range
:math:`\tau \in [0, 1]`.  The canonical definition-of-record, the
:math:`[1/2, 1]`-clamp history, the W1 mis-citation finding, and the
geometry split are on the theory page — see
:eq:`morel-montry-clamp` (the equation-of-record) in
:doc:`/theory/structured_geometry` and the comprehensive vindication
at :ref:`sn-curvilinear-aniso-norm-reconciliation` in
:doc:`/theory/discrete_ordinates`.

**Geometry split (W1, 2026-06-13; see** :func:`spherical_streaming` /
:func:`cylindrical_streaming` **for the implementations).**

* **SPHERE** uses the raw :math:`\tau_m` above directly — unclamped.
  The former :math:`[1/2, 1]` clamp was an over-conservative
  positivity floor (mis-cited to Lewis & Miller §4.5) that re-floored
  the anisotropic solution while being 100 % spurious on physical
  fields (unclamped sphere SI is stable + strictly positive on
  thick-absorber / near-vacuum / :math:`c = 0.999` / S64 stress
  configs).  On Gauss--Legendre quadrature :math:`\tau_{m}^{\rm raw}
  \in [0.39, 0.61]` (never 0), so the unclamped weight is always
  admissible.
* **CYLINDER** retains the clamp
  :math:`\tau_m = \mathrm{clip}(\tau_m^{\rm raw}, \tfrac12, 1)`,
  because product / level-symmetric quadratures place the most-inward
  azimuthal ordinate exactly on :math:`\eta = -\sin\theta`, giving
  :math:`\tau_m^{\rm raw} = 0` exactly — a structural :math:`\div 0`
  block in the unclamped recurrence.  Removing the clamp here needs a
  2-D :math:`(\eta, \varphi)` angular closure (out of scope), tracked
  by Issue #229.

References
==========

* Hébert, A. (2009). *Applied Reactor Physics*.  Ch. 3 §3.9.4 (pp.
  141-144), Eqs. 3.418-3.439.  **Primary source** for the curvilinear
  S\ :sub:`N` discretization (cell balance + DD difference relations
  + Carlson starting-direction).  Local copy:
  ``scratch/literature/Hebert(2009)Chapter3.pdf``.
* Bailey, T. S., Morel, J. E., & Chang, J. H. (2010).  *Asymptotic
  Diffusion-Limit Accuracy of Sn Angular Differencing Schemes*.
  NSE 165(2):149-169 (LLNL preprint LLNL-JRNL-420356; OA at
  https://www.osti.gov/servlets/purl/1020346).  **Eq. 43** gives the
  Morel--Montry weight :math:`\tau_m` above as the unique weight exact
  for a flux linear in :math:`\mu`, admissible range :math:`\tau \in
  [0, 1]` — the W1 source for dropping the spherical clamp.  The paper
  also frames the :math:`\beta`-contamination diffusion-limit analysis
  via formal-:math:`\varepsilon` expansion.

  *Citation correction (Issue #168 Phase B)*: this module's
  pre-Phase-B docstring cited "Bailey, T. S., Adams, M. L., Yang, B.,
  & Zika, M. R. (2009).  *A piecewise linear finite element
  discretization of the diffusion equation for arbitrary polyhedral
  grids*. JCP 227, 3738-3757", which is the **wrong** Bailey paper —
  Bailey-Adams-Yang-Zika is a piecewise-linear FE diffusion paper
  unrelated to curvilinear S\ :sub:`N` :math:`\alpha`-recursion.  The
  intended reference is Bailey-Morel-Chang 2010 (above).
* Lewis, E. E., & Miller, W. F. (1984). *Computational Methods of
  Neutron Transport*.  §4.5 — angular redistribution closure.  (NB:
  the historical :math:`[1/2, 1]` clamp was MIS-cited to this §4.5 —
  L&M §4.5 does not prescribe it; W1 traced the exact-on-linear weight
  to Bailey-Morel-Chang 2010 Eq. 43.)

See also
========

* :class:`orpheus.sn.mesh.augmented_mesh.SNMesh` — consumes this primitive
  via the curvilinear connection-coefficients path.
  :class:`~orpheus.sn.operators.streaming.StreamingOperator` consumes it through
  :func:`~orpheus.sn.operators.streaming.transport_operator_matvec_unified`.
  MoC and CP campaigns (post-Wave-1) reuse this primitive with their
  own consumption patterns.
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

    # Read-only properties (not plain attributes): this module only READS
    # the angular measure, and concrete implementers (``Quadrature``) expose
    # these as ``@property``. A mutable attribute in the Protocol is invariant
    # and rejects a read-only ``@property`` ("property is not assignable to
    # ndarray"); a read-only Protocol member accepts both a property and a
    # plain attribute, matching the actual usage (covariant read contract).

    @property
    def mu_x(self) -> np.ndarray:
        """(N,) primary direction cosine — :math:`\\mu` for sphere,
        :math:`\\eta` (radial) for cylindrical product quadrature."""
        ...

    @property
    def weights(self) -> np.ndarray:
        """(N,) quadrature weights."""
        ...

    @property
    def N(self) -> int:
        """Number of ordinates."""
        ...

    @property
    def eta(self) -> np.ndarray:
        r"""(N,) cylindrical-frame radial direction cosine
        :math:`\eta = \Omega \cdot \hat{r}` (read by the curvilinear
        connection-coefficient recursion only)."""
        ...

    @property
    def mu_z(self) -> np.ndarray:
        """(N,) axial direction cosine (curvilinear / level-structure side)."""
        ...

    @property
    def level_indices(self) -> list[np.ndarray]:
        """Per-:math:`\\mu`-level ordinate-index partition: ``[arange(N)]`` for
        a slab quadrature (single level), one index array per level for
        cylindrical-compatible cubatures."""
        ...


# ═══════════════════════════════════════════════════════════════════════
# StreamingTerms — geometry-dependent shape per (cell, direction)
# ═══════════════════════════════════════════════════════════════════════

@dataclass(frozen=True, slots=True)
class StreamingTerms:
    """Per-(cell, direction) **purely geometric** inputs for a sweep cell update.

    All curvature fields are populated for **every** geometry
    (Issue #196 Phase G Step 2.5):

    * **Slab**: neutral-curvature values — ``face_area_inner =
      face_area_outer = 1.0``, ``delta_A_over_w = 0.0``,
      ``alpha_in = alpha_out = 0.0``, ``tau_mm = 1.0`` (synthetic
      neutral element of the M-M angular closure — slab has no
      half-angles; the closure would be the identity if ever
      evaluated, but ``upstream_state.angular_upstream`` is ``None``
      for slab so it never is).  Plus the always-populated
      ``chord_length``, ``mu``, ``volume``, ``abs_mu``.
    * **Sphere / cylinder**: physically-populated curvature fields
      from the dome recursion + M-M clamp.

    Cylindrical ``mu`` / ``abs_mu`` are read from the global ordinate
    ``mu_x[level_indices[mu_level_idx][direction_idx]]`` because
    cylindrical ``direction_idx`` is the within-level azimuthal index
    :math:`m \\in [0, M)`, not the global ordinate.

    Before Step 2.5, slab left the curvature fields as ``None`` and
    cell-update strategies branched on ``alpha_in is None`` to
    discriminate slab from curvilinear.  Step 2.5 retires that
    branch — slab carries neutral curvature so the unified
    cell-balance helper consumes the same packet regardless of
    geometry.  Geometry is data, not control-flow (Cardinal Rule 2
    + Pattern 4: make illegal states unrepresentable).

    Geometric, not direction-resolved
    =================================

    Per Cardinal Rule 2 (architecture), this primitive carries
    **purely geometric** labels — labels that are independent of the
    sweep's marching direction.  The two face-area fields are named
    by their geometric position relative to :math:`r=0`:

    * :attr:`face_area_inner` is :math:`A_{i-1/2}` — the face closer
      to the centre of the geometry (smaller :math:`r`).
    * :attr:`face_area_outer` is :math:`A_{i+1/2}` — the face farther
      from the centre (larger :math:`r`).

    These labels do **not** depend on which way the sweep is
    marching.  For an outward sweep (centre → boundary,
    :math:`\\mu > 0`) the inner face is the upstream / incoming face;
    for an inward sweep (boundary → centre, :math:`\\mu < 0`) the
    outer face is the upstream / incoming face.  The
    sweep-direction resolution — *which* of the two faces is
    "downstream" for a given visit — lives in the SN module
    (:class:`~orpheus.sn.spatial.scheme.CellVisit` packs the
    geometric :class:`StreamingTerms` together with the
    sweep-resolved :attr:`face_area_downstream`).  This module is
    geometry-only and is reusable by future MoC / CP / diffusion,
    none of which has SN's sweep-direction concept.

    The trailing ``volume`` and ``abs_mu`` fields are populated by
    **all three factories** so that a downstream
    :class:`~orpheus.sn.spatial.scheme.DiscretizationScheme` strategy
    receives a self-contained per-cell, per-direction packet and need
    not reach back into ``SNMesh`` or the ``AngularQuadrature``.  Every
    curvature field is populated for **every** geometry (Step 2.5; slab
    carries neutral values), so cell-update strategies consume the same
    packet regardless of chart — geometry is **data, not control-flow**.
    There is no ``alpha_in is None`` slab-vs-curvilinear branch; the
    type makes the all-populated invariant unrepresentable-otherwise
    (every field is a required ``float``, Pattern 4).
    """

    chord_length: float
    """Cell radial width (slab/sphere/cylinder all use ``mesh.widths[i]``)."""

    mu: float
    """Signed primary direction cosine for this ordinate.

    :math:`\\mu` for slab and sphere (axial); :math:`\\eta` for
    cylindrical 1-D radial sweeps (the radial direction cosine,
    with the global ordinate index resolved through
    :attr:`AngularMeasure.level_indices`).  Signed.

    Slab DD's flow direction reads off the sign of this field.
    Cell-update strategies that need only the magnitude use
    ``abs_mu`` instead.
    """

    face_area_inner: float
    """:math:`A_{i-1/2}` — area of the **inner** radial face
    (closer to :math:`r=0`).

    Geometric label, independent of sweep direction.  See class
    docstring "Geometric, not direction-resolved".
    """

    face_area_outer: float
    """:math:`A_{i+1/2}` — area of the **outer** radial face
    (farther from :math:`r=0`).

    Geometric label, independent of sweep direction.  See class
    docstring "Geometric, not direction-resolved".
    """

    delta_A_over_w: float
    """:math:`\\Delta A_i / w_n` — the geometry-redistribution factor."""

    alpha_in: float
    """:math:`\\alpha_{n-1/2}` — incoming half-angle redistribution.

    Step 2.5 (Issue #196 Phase G): always populated.  Slab carries
    ``0.0`` (neutral); sphere/cylinder carry the dome value.  The
    historical ``alpha_in is None`` slab-vs-curvilinear branch is
    retired in favour of geometry-as-data via the neutral values.
    """

    alpha_out: float
    """:math:`\\alpha_{n+1/2}` — outgoing half-angle redistribution."""

    tau_mm: float
    """Morel--Montry angular closure weight on this ordinate."""

    mu_start: float
    """Starting-direction angular edge of this ordinate's M-M level.

    The direction the half-angle thread's seed flux lives at:
    sphere ``-1.0``; cylinder :math:`-\\sqrt{1-\\xi_p^2}` (the level's
    most-inward azimuthal edge); slab ``-1.0`` (neutral — the identity
    closure never reads it).  Constant within a level.  Consumed by
    :class:`~orpheus.sn.spatial.psi_half_angle_seed.AngularEdgeExtrapolation`
    via :class:`~orpheus.sn.spatial.psi_half_angle_seed.CarlsonSweepContext`.
    """

    volume: float
    """Cell volume :math:`V_i`.

    Populated by all three factories from
    :attr:`~orpheus.geometry.mesh.Mesh1D.volumes`: slab uses
    ``mesh.volumes[i]`` (which equals ``widths[i]`` for unit
    cross-section in 1-D Cartesian); sphere uses
    :math:`\\tfrac{4}{3}\\pi(r_{i+1}^3 - r_i^3)`; cylinder uses
    :math:`\\pi(r_{i+1}^2 - r_i^2)` (per unit axial length).  Carried
    on the streaming-terms packet so that the
    :class:`~orpheus.sn.spatial.scheme.DiscretizationScheme` cell-update
    contract receives :math:`V_i` directly without needing access to
    the underlying ``SNMesh``.
    """

    abs_mu: float
    """Absolute primary direction cosine.

    :math:`|\\mu|` for slab and sphere (axial direction cosine);
    :math:`|\\eta|` for cylindrical 1-D radial sweeps (the radial
    direction cosine).  Slab and sphere factories compute this as
    ``abs(quadrature.mu_x[direction_idx])`` — for those geometries
    ``direction_idx`` is the global ordinate index.  The cylindrical
    factory resolves the global ordinate through
    :attr:`AngularMeasure.level_indices` because cylindrical
    ``direction_idx`` is the within-level azimuthal index
    :math:`m \\in [0, M)`, not the global ordinate index.  By
    ORPHEUS convention, ``mu_x`` carries :math:`\\mu` for sphere and
    :math:`\\eta` for cylinder, and both are real-valued, so the
    absolute value is always well-defined.

    The cylindrical pure-azimuthal degenerate case (``abs_mu <
    1e-15``) is signalled by cell-update strategies via
    :attr:`~orpheus.sn.spatial.scheme.CellResult.outgoing_spatial_flux`
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
    mu_start: float | None = None
    """Starting-direction angular edge :math:`\\mu_{1/2}` of the (single)
    M-M level — the direction the half-angle thread's seed flux lives
    at.  Sphere: ``-1.0`` (the Hébert §3.9.4 starting direction).
    Defined HERE, at the same construction site as the α-dome and τ
    (single source of truth for the angular cell partition)."""

    # Cylindrical (per-level)
    alpha_per_level: list[np.ndarray] | None = None
    redist_dAw_per_level: list[np.ndarray] | None = None
    tau_mm_per_level: list[np.ndarray] | None = None
    mu_start_per_level: list[float] | None = None
    """Per-level starting-direction angular edge — cylinder:
    :math:`\\eta_{1/2} = -\\sin\\theta_p = -\\sqrt{1-\\xi_p^2}` (the
    most-inward azimuthal edge of level *p*)."""

    # Quadrature reference (kept for streaming_terms() extraction)
    _quadrature: AngularMeasure | None = field(default=None, repr=False)

    # ── Per-direction extraction ───────────────────────────────────

    def streaming_terms(
        self,
        cell_idx: int,
        direction_idx: int,
        mu_level_idx: int | None = None,
    ) -> StreamingTerms:
        """Pack the per-cell, per-direction **purely geometric** streaming inputs.

        For slab and sphere, ``direction_idx`` is the global ordinate
        index and the signed primary direction cosine is read directly
        as ``quadrature.mu_x[direction_idx]``.

        For cylindrical geometry, ``mu_level_idx`` selects which
        :math:`\\mu`-level the ``direction_idx`` belongs to; the
        ``direction_idx`` is then the within-level azimuthal index
        :math:`m \\in [0, M)`.  The global ordinate index is resolved
        as ``quadrature.level_indices[mu_level_idx][direction_idx]``;
        the signed :math:`\\eta` and absolute :math:`|\\eta|` are read
        from that global index.

        The returned :class:`StreamingTerms` carries **purely
        geometric** labels — :attr:`face_area_inner` /
        :attr:`face_area_outer` are the inner / outer radial faces of
        the cell relative to :math:`r = 0`.  Sweep-direction
        resolution (which face is downstream for a given visit) lives
        in the SN module — see
        :class:`orpheus.sn.spatial.scheme.CellVisit`.
        """
        chord = float(self.mesh.widths[cell_idx])
        # Common per-(cell, direction) volume.
        # mesh.volumes returns shape (N,) for 1-D meshes.
        assert self._quadrature is not None
        volume = float(self.mesh.volumes[cell_idx])

        if self.coord is CoordSystem.CARTESIAN:
            # Slab — neutral curvature values populate the curvilinear
            # fields so cell_balance_terms_unified can consume the
            # packet without geometry dispatch (Issue #196 Phase G
            # Step 2.5).  Slab carries:
            #   face_area_inner = face_area_outer = 1.0  (so A_total =
            #       A_inner + A_outer = 2, and 2*|μ|*A_down = 2|μ|·1
            #       reproduces the slab denominator's "2|μ|" term);
            #   delta_A_over_w = 0.0  (no curvature redistribution);
            #   alpha_in = alpha_out = 0.0  (no half-angle dome);
            #   tau_mm = 1.0  (closure is the identity for slab —
            #       the synthetic neutral element of M-M; documented
            #       in dd_polymorphism.md §5.3).
            # ``volume == chord`` is already true for slab (unit
            # cross-section in 1-D Cartesian).
            mu_n = float(self._quadrature.mu_x[direction_idx])
            return StreamingTerms(
                chord_length=chord,
                mu=mu_n,
                face_area_inner=1.0,
                face_area_outer=1.0,
                delta_A_over_w=0.0,
                alpha_in=0.0,
                alpha_out=0.0,
                tau_mm=1.0,
                mu_start=-1.0,
                volume=volume,
                abs_mu=abs(mu_n),
            )

        if self.coord is CoordSystem.SPHERICAL:
            assert self.face_areas is not None
            assert self.alpha_half is not None
            assert self.redist_dAw is not None
            assert self.tau_mm is not None
            assert self.mu_start is not None
            # Sphere: ``direction_idx`` IS the global ordinate index.
            mu_n = float(self._quadrature.mu_x[direction_idx])
            return StreamingTerms(
                chord_length=chord,
                mu=mu_n,
                face_area_inner=float(self.face_areas[cell_idx]),
                face_area_outer=float(self.face_areas[cell_idx + 1]),
                delta_A_over_w=float(
                    self.redist_dAw[cell_idx, direction_idx]
                ),
                alpha_in=float(self.alpha_half[direction_idx]),
                alpha_out=float(self.alpha_half[direction_idx + 1]),
                tau_mm=float(self.tau_mm[direction_idx]),
                mu_start=float(self.mu_start),
                volume=volume,
                abs_mu=abs(mu_n),
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
            assert self.mu_start_per_level is not None
            # Cylinder: ``direction_idx`` is the within-level azimuthal
            # index; the global ordinate is read through
            # ``level_indices``.  ``mu_x[global_n]`` carries η (the
            # radial direction cosine).
            level_indices = self._quadrature.level_indices
            global_n = int(level_indices[mu_level_idx][direction_idx])
            eta_n = float(self._quadrature.eta[global_n])
            alpha_lv = self.alpha_per_level[mu_level_idx]
            dAw_lv = self.redist_dAw_per_level[mu_level_idx]
            tau_lv = self.tau_mm_per_level[mu_level_idx]
            return StreamingTerms(
                chord_length=chord,
                mu=eta_n,
                face_area_inner=float(self.face_areas[cell_idx]),
                face_area_outer=float(self.face_areas[cell_idx + 1]),
                delta_A_over_w=float(dAw_lv[cell_idx, direction_idx]),
                alpha_in=float(alpha_lv[direction_idx]),
                alpha_out=float(alpha_lv[direction_idx + 1]),
                tau_mm=float(tau_lv[direction_idx]),
                mu_start=float(self.mu_start_per_level[mu_level_idx]),
                volume=volume,
                abs_mu=abs(eta_n),
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

    Implements Hébert (2009) §3.9.4 Eqs. 3.423-3.424 (α-dome
    recursion, in the ORPHEUS factor-of-2-absorbed normalization) +
    Bailey-Morel-Chang (2010) Eq. 5 (Morel--Montry :math:`\tau`
    clamp), producing arrays bit-identical to
    :class:`SNMesh._setup_spherical`.

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
    # mesh, which routes through coord.compute_areas_1d().
    face_areas = mesh.areas  # (nx+1,)

    # Cell face-area differences: ΔA_i = A_{i+1/2} − A_{i-1/2}
    delta_A = face_areas[1:] - face_areas[:-1]

    # Hébert (2009) §3.9.4 Eqs. 3.423-3.424 α-dome recursion (ORPHEUS
    # factor-of-2-absorbed normalization):
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

    # Morel–Montry / weighted-diamond angular weight: the fractional
    # position of μ_n within its half-angle interval [μ_{n-1/2}, μ_{n+1/2}].
    # This is Bailey–Morel–Chang (2010) Eq. 43 — the UNIQUE weight that is
    # exact for a flux linear in μ.  μ_{1/2} = -1, μ_{N+1/2} = +1.
    #
    # NO [1/2, 1] clamp (W1, branch fix/curvilinear-aniso-pole-and-clamp):
    # Bailey allows τ ∈ [0, 1] and recommends *exactly this* unclamped
    # weight; the former clamp was an over-conservative positivity floor
    # mis-cited to Lewis & Miller §4.5.  It re-floored the anisotropic
    # solution (it threads a linear-in-μ field only to O(clamp gap), not
    # to machine precision) WITHOUT preventing any real negativity —
    # verified 100% spurious on every physical converged sphere solve
    # (thick absorber, near-vacuum, c=0.999, S64; all strictly positive).
    # The sphere dome is non-singular (dmu = w_n > 0 on GL, τ_raw ∈ ~[0.39,
    # 0.61]) so the ½ degeneracy fallback never fires.  The CYLINDER keeps
    # its clamp — there τ_raw = 0 exactly (structural ÷0 block; see
    # ``cylindrical_streaming`` and #229).
    mu_edge = np.zeros(N + 1)
    mu_edge[0] = -1.0
    for n in range(N):
        mu_edge[n + 1] = mu_edge[n] + w[n]
    tau_mm = np.empty(N)
    for n in range(N):
        dmu = mu_edge[n + 1] - mu_edge[n]
        tau_mm[n] = (mu[n] - mu_edge[n]) / dmu if abs(dmu) > 1e-15 else 0.5

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
        mu_start=float(mu_edge[0]),
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
    if getattr(angular_measure, "level_structure", None) is None:
        raise ValueError(
            "Cylindrical streaming requires a quadrature with level "
            "structure (Quadrature.level_symmetric or .product); "
            "slab quadratures (Quadrature.gauss_legendre) and pure "
            "sphere cubatures (Quadrature.lebedev) carry no "
            "LevelStructure side-channel."
        )

    face_areas = mesh.areas  # (nx+1,)
    delta_A = face_areas[1:] - face_areas[:-1]

    # Per-level azimuthal redistribution coefficients
    # Hébert §3.9.4 (cylindrical analog): α_{m+1/2} = α_{m-1/2} − w_m · η_m
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

    # Morel–Montry closure (Bailey-Morel-Chang 2010 Eq. 43) per μ-level.
    # For cylindrical, ordinates are η-sorted but weights come from
    # φ-space (not η-space), so the weight-sum edge approach is wrong.
    # Instead, cell edges are at midpoints of consecutive η values
    # with endpoints at ±sin θ.  This gives a proper η-partition.
    #
    # The [1/2, 1] clamp is RETAINED here (unlike the sphere, where W1
    # removed it — see ``spherical_streaming``).  The clamp is STRUCTURAL
    # for the cylinder: the most-inward azimuthal ordinate sits exactly on
    # the level boundary (``eta[0] == eta_edge[0] == -sin θ``), so its raw
    # weight is ``τ_raw = (eta[0] − eta_edge[0])/deta = 0`` bit-exactly —
    # the recurrence ``ψ_out = (ψ_avg − (1−τ)ψ_in)/τ`` would divide by zero
    # unclamped.  Removing the clamp needs a 2-D (η, φ) angular closure
    # (the product/LS quadratures carry duplicate azimuthal η that a 1-D
    # η-thread cannot represent); that is out of scope and tracked by #229.
    #
    # (See SNMesh._setup_cylindrical for the same construction;
    # GitHub Issue #3 tracks the φ-based edge refinement.)
    mu_z = angular_measure.mu_z
    tau_mm_per_level: list[np.ndarray] = []
    mu_start_per_level: list[float] = []
    for level_idx in angular_measure.level_indices:
        eta = angular_measure.mu_x[level_idx]
        M = len(level_idx)
        sin_theta = np.sqrt(1.0 - mu_z[level_idx[0]] ** 2)
        mu_start_per_level.append(float(-sin_theta))
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
        mu_start_per_level=mu_start_per_level,
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
