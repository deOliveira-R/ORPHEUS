r"""Reduced streaming operator — connection coefficients on the coordinate chart.

This module hosts the **connection-coefficient primitive** that a
curvilinear SN sweep marches through.  In differential-geometric
language, the spherical :math:`(1-\mu^2)/r\,\partial_\mu` and cylindrical
:math:`-(1/r)\,\partial_\varphi(\xi\cdot)` redistribution terms are
**the same connection-coefficient operator** on SO(3) viewed in two
coordinate charts (polar-on-sphere vs. azimuthal-on-cylinder).  The
numerical data the sweep needs is chord lengths, face areas, the
:math:`\Delta A/w` geometry factor, and the :math:`\alpha` dome
recursion (Hébert Eqs. 3.423-3.424, after Lathrop &
Carlson 1966 — **not** "Bailey 2009", the wrong-paper citation this
module's References corrected at Issue #168 Phase B).  The Morel--Montry
angular closure weight :math:`\tau` is a *consumer* of these, not a
product of this module — see the ⚠ note at the end of "Mathematical
content".

⛔ **This module read "shared by SN, MoC, and CP curvilinear sweeps"
and "Each solver — SN, MoC, CP — needs the same numerical data" until
2026-08-27.  Both were false, and structurally so.**  `[M]` no file
under ``orpheus/moc/``, ``orpheus/cp/`` or ``orpheus/mc/`` names this
primitive under **any** of eight independent spellings — three symbol
spellings, three concept spellings, an attribute-access spelling and
the prose paraphrase — while ``orpheus/sn/`` and ``orpheus/geometry/``
name every one of them.  That is not an unfinished migration.

⚠ The exact patterns, the root, the occurrence semantics and the
positive control are published **as a runnable predicate rather than a
table of counts**, at ``docs/theory/foundations/structured_geometry.rst``
``§connection-coefficient-census``.  Do not copy counts back here: a
control's value carries no part of the argument (only its non-zeroness
does), and it moves under any edit to those packages — including this
very correction, which names the module and so raises several of them.
An earlier revision of this block did print them and was wrong three
ways at once (pre-edit values, a column set that did not match the
spelling list printed beside it, and an unanchored case-insensitive
``redistribut`` column that absorbed every ``AngularRedistribution``).

An angular redistribution term exists only where an angular unknown is
indexed in a *local rotating frame* and its angular derivative is
*collocated*: MoC fails the middle condition (Ω is fixed in the global
frame, so Ω·∇ = d/ds is chart-free and the curvature moves into track
segmentation) and CP fails the first (angle is integrated out into the
kernel before anything is discretised).  The curvilinear
counter-examples are decisive, because "not migrated yet" predicts
there is nothing to migrate: MoC ray-traces concentric annuli on a
cylindrical mesh (``moc/geometry.py`` ``MOCMesh``), CP ships a real
sphere (``cp/solver.py`` ``_setup_spherical``) and MC a real cylinder
(``mc/solver.py`` ``MCMesh``) — all three with no α anywhere.  Full
argument and falsification conditions:
``docs/theory/foundations/structured_geometry.rst``,
"Who needs a connection coefficient — and who does not".

Per Cardinal Rule 2 (architecture), this primitive **MUST NOT** be
duplicated across solvers.  The historical home was a pair of in-line
setup methods on :class:`~orpheus.sn.mesh.augmented_mesh.SNMesh`
(``_setup_spherical`` / ``_setup_cylindrical``); the migration has
LANDED — those methods are gone and ``SNMesh.__init__`` now calls
:func:`spherical_streaming` / :func:`cylindrical_streaming` here,
storing the result as its ``reduced`` attribute.

Bit-identity with those retired methods was the migration's
load-bearing contract, established at the Wave-B carve while both
implementations were live.

.. warning::

   That contract is **no longer independently gated**.  The
   factory-binding tests in ``tests/geometry/test_reduced_operator.py``
   compare this module's factory output against ``sn_mesh.reduced.*``
   (and against the deprecated ``SNMesh.face_areas`` / ``delta_A``
   read-throughs, which forward to the same object) — i.e. against
   the value this module itself produced, through the mesh
   constructor.  They pin the WIRING (SNMesh really does bind to
   these factories) and would redden if that wiring broke; they can
   no longer detect a change in the connection-coefficient math,
   because there is no second implementation left to disagree with.
   Measured 2026-08-03: garbaging every array these factories emit
   leaves **all 47** tests in that file green.

   ⛔ **The per-field table below is a 2026-08-03 SNAPSHOT and is no
   longer a description of this class.**  The 2026-08-26 un-weld
   retired ``alpha_half`` and ``redist_dAw`` / ``redist_dAw_per_level``
   outright, and re-homed ``alpha_per_level`` onto
   :class:`AngularRedistribution`; ``delta_A`` is now read through
   :attr:`ReducedStreamingOperator.redistribution_pairing`.  So three of
   the five rows name a field that no longer exists here.

   It is kept, tensed, rather than deleted or quietly refreshed: the
   mutation that produced it was real, and the *reasoning* — which
   catcher is structurally independent, and which is a blind snapshot —
   is what a future audit needs.  ⚠ What it must **not** be read as is
   current coverage: the per-field mutation has **not** been re-run
   since the un-weld, so whether each pin survived its field's
   re-homing is **unmeasured**.  Re-measure before citing any row as
   evidence (`coding-standards`: a retirement can silently promote or
   demote a gate's claim class without touching one line of the test).

   As measured on 2026-08-03, the surviving numerical pins on the math
   were, per field — every one of them **structurally independent**
   (a closed form), with the SN curvilinear regression snapshots
   (``tests/sn/regression/test_dd_regression.py``) corroborating but
   nowhere the sole evidence:

   * ``delta_A`` — ``tests/sn/primitives/test_quadrature.py::
     TestL0TermVerification::test_delta_A_magnitude`` (closed form
     :math:`4\pi\,\Delta(r^2)` / :math:`2\pi\,\Delta r`), sole catcher;
     the snapshots are blind, correctly — ``delta_A`` has no
     production consumer.
   * ``alpha_half`` — the L0 per-ordinate flat-flux identity
     ``test_per_ordinate_flat_flux_consistency[SPHERICAL]``
     (``catches("ERR-006", "ERR-007")``) + the sphere snapshots.
   * ``alpha_per_level`` — ``tests/sn/sweep/curvilinear/
     test_alpha_closed_form.py`` (Dirichlet-kernel closed form) + the
     cylindrical flat-flux arm + the cylinder snapshots.
   * ``redist_dAw`` / ``redist_dAw_per_level`` — ``tests/sn/sweep/
     curvilinear/test_streaming_equilibrium_curvilinear.py``, the L0
     closed-form :math:`\varphi = Q/(\Sigma_t(1-c))` gate (15 and 12
     of its 27 cases redden respectively) + both snapshot families.
     The flat-flux identity does NOT cover these — it recomputes
     ``dA / w`` instead of reading the production array.
   * ``face_areas`` — ``tests/geometry/test_geometry.py`` pins the
     producer :func:`~orpheus.geometry.coord.compute_areas_1d` against
     the closed form; the snapshots pin the forwarding.

   ``tests/sn/sweep/curvilinear/test_tau_producer_equivalence.py`` is
   **not** among them, despite an earlier revision of this warning
   naming it: Issue #236 Step C moved :math:`\tau` to the angular
   closure, which derives it from :math:`(\mu, w)` alone, so that gate
   passes untouched (5 passed, 0.03 s) under fully-garbaged factories.

Mathematical content
====================

The :math:`\alpha`-dome recursion (sphere) — Hébert (2009)
§3.9.4 Eqs. 3.423-3.424, in the ORPHEUS factor-of-2-absorbed
normalization:

.. math::
   :label: alpha-dome-recursion

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
Eq. 43 = Lathrop 2000 Eq. 23) is the *barycentric coordinate* of the
ordinate :math:`\mu_m` between the two edges of its own angular cell
:math:`[\mu_{m-1/2}, \mu_{m+1/2}]`:

.. math::

   \tau_m = \frac{\mu_m - \mu_{m-1/2}}{\mu_{m+1/2} - \mu_{m-1/2}}.

It is the UNIQUE choice exact for a flux affine in the radial cosine,
with admissible range :math:`\tau \in [0, 1]`.  The formula carries **no
geometry**: the geometry lives entirely in the *angular cell partition*
it reads.  The canonical definition-of-record and the partition are on
the theory page — see :eq:`morel-montry-closure` (the closure) and
:eq:`angular-cell-partition` (the partition) in
:doc:`/theory/foundations/structured_geometry`, the retired-absorber
account at :ref:`sn-tau-absorber-retirement`, and the comprehensive W1
vindication at :ref:`sn-curvilinear-aniso-norm-reconciliation` in
:doc:`/theory/verification/sn`.

**No clamp, in either arm** (sphere unclamped at W1, 2026-06-13;
cylinder at Q5.6.4, 2026-08-11).  On Gauss--Legendre the sphere's
:math:`\tau_m \in [0.39, 0.61]`; on the shipped folded cylinder arc
:math:`\tau_m \in [\tfrac14, \tfrac34]`.  The former
:math:`[\tfrac12, 1]` clamp was an over-conservative positivity floor
(mis-cited to Lewis & Miller §4.5) on the sphere, and on the cylinder a
compensation for a *wrong cell partition* — the retired chord-midpoint
edges, which disagreed with :math:`\alpha`'s half-angle edges by
~17.5 % in :math:`\omega`-width.  ``morel_montry_tau_raw_per_level``
retired with it: there is one :math:`\tau`.

⚠ This module no longer produces :math:`\tau` at all — #236 Step C
excised the geometry-side producer.  :math:`\tau` is single-sourced in
:func:`~orpheus.sn.angular.closure.morel_montry_tau_per_level`,
reading
:func:`~orpheus.sn.angular.closure.angular_cell_edges_per_level`.
What stays here is genuinely geometric: the :math:`\alpha`-dome, the
redistribution factor, the face areas, and the starting-direction edges.

References
==========

* Hébert, A. (2009). *Applied Reactor Physics*.  Ch. 3 — **§3.9.3
  (cylinder, printed pp. 137-141)** and **§3.9.4 (sphere, printed pp.
  141-144)**; the whole Eq. 3.418-3.439 range is *spherical*.  The
  authority for the curvilinear S\ :sub:`N` **cell balance, DD difference
  relations, sweep ordering and Carlson starting direction** — i.e. for
  everything this module still produces.  He defines **no**
  :math:`\tau` anywhere in chapter 3, in either geometry, so he is
  **not** an authority for the weighted angular closure above; that
  mis-citation produced a wrong cylinder :math:`\tau` and is recorded at
  ``docs/theory/methods/sn/curvilinear_one_group.rst
  §sn-citation-corrections``.  Local copy:
  ``scratch/literature/Hebert(2009)Chapter3.pdf``.
* Morel, J. E., & Montry, G. R. (1984).  *Analysis and Elimination of
  the Discrete-Ordinates Flux Dip*.  Transport Theory and Statistical
  Physics 13(5):615-633, doi:10.1080/00411458408211661.  **PRIMARY**
  source for the weighted angular closure :math:`\tau` above.  Local
  copy: ``scratch/literature/Morel-Montry(1984)Analysis and elimination
  of the discrete-ordinates flux dip.pdf``.
* Bailey, T. S., Morel, J. E., & Chang, J. H. (2010).  *The Asymptotic
  Diffusion-Limit Accuracy of Sn Angular Differencing Schemes*.
  NSE 165(2):149-169 (LLNL preprint LLNL-JRNL-420356; OA at
  https://www.osti.gov/servlets/purl/1020346).  **Eqs. (42)/(43)** are
  the form implemented: the Morel--Montry weight :math:`\tau_m` above as
  the unique weight exact for a flux affine in the radial cosine,
  admissible range :math:`\tau \in [0, 1]` — the W1 source for dropping
  the spherical clamp.  Their **Eq. (41)** is the first-order
  diffusion-limit condition :math:`\beta = 0`, and forcing it to zero is
  what DETERMINES the weights; the paper frames that
  :math:`\beta`-contamination analysis via formal-:math:`\varepsilon`
  expansion, and its Eq. 53 names :math:`\tau = \tfrac12` "the diamond
  scheme".
* Reed, W. H., & Lathrop, K. D. (1970).  *Truncation Error Analysis of
  Finite Difference Approximations to the Transport Equation*.  NSE
  41(2):237-248, doi:10.13182/NSE70-A20710.  Their **Eq. (13c)** IS BMC
  Eq. (43), forty years earlier; their **Eqs. (15)/(16)** give the
  sharpest accuracy criterion on :math:`\tau` — second order iff the
  ordinate is the :math:`\mu`-midpoint of its own cell to
  :math:`O(w^2)`, i.e. iff :math:`\tau = \tfrac12 + O(w)`.  Unlike
  :math:`\beta` that criterion is POINTWISE.  Local copy:
  ``scratch/literature/07-Truncation Error Analysis of Finite Difference
  Approximations to the Transport Equation.pdf``.

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
  to Bailey-Morel-Chang 2010 Eq. 43.  ⭐ The interval's real origin is
  Grant, I. P. (1968), *J. Comp. Phys.* 2(4):381-402,
  doi:10.1016/0021-9991(68)90044-2, where :math:`[\tfrac12, 1]` bounds
  the **SPATIAL** weighted-diamond parameter; Reed & Lathrop's footnote
  8 says so and adds that "Grant does not determine angular weights".
  ⚠ Grant 1968 is not in the local library and has not been read.)

See also
========

* :class:`orpheus.sn.mesh.augmented_mesh.SNMesh` — consumes this primitive
  via the curvilinear connection-coefficients path (as ``SNMesh.reduced``).
  :class:`~orpheus.sn.operators.streaming.StreamingOperator` consumes it
  through the loss-representation walk (its
  :meth:`~orpheus.sn.operators.streaming.StreamingOperator.apply` reads
  ``self.mesh.coord`` inside the walk).  (The per-geometry
  ``transport_operator_matvec_*`` matvecs and their unified successor,
  which consumed it in the Depth-B era, were deleted in the typed-field
  (#197) / walk-unification (#280 campaigns) refactors.)
  ⛔ This bullet read "MoC and CP campaigns (post-Wave-1) reuse this
  primitive with their own consumption patterns" until 2026-08-27.  It
  is retracted as NOT APPLICABLE rather than pending — neither method
  forms an angular redistribution term at all.  See the ⛔ block in
  "Mathematical content" above and, for the full argument,
  ``structured_geometry.rst`` "Who needs a connection coefficient — and
  who does not".
* :doc:`/theory/foundations/structured_geometry` — the architecture page;
  see "Connection coefficients (reduced streaming operator)".
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property
from typing import Protocol

import numpy as np

from .coord import CoordSystem
from .mesh import Mesh1D


# ═══════════════════════════════════════════════════════════════════════
# Angular-measure adapter (Issue 4 contract — local, structural)
# ═══════════════════════════════════════════════════════════════════════

class AngularMeasure(Protocol):
    """Minimal contract this module needs from an angular quadrature.

    The concrete implementer is
    :class:`~orpheus.numerics.quadrature.Quadrature`, which wraps a
    :class:`~orpheus.numerics.measure.DiscreteMeasure`; here we only
    depend on the attributes the connection-coefficient math reads, so
    the geometry layer needs no import from the quadrature package at
    all.  This structural Protocol outlived the type it was written
    against: the ``AngularQuadrature`` Protocol and its four
    per-family adapter classes at ``orpheus.sn.quadrature`` were
    retired into classmethod factories on the single ``Quadrature``
    type (R-1 Phase A detour-C), and the read contract below did not
    have to move.
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
      face_area_outer = 1.0``, ``delta_A_over_w = 0.0``.  The
      Morel–Montry angular weight is NOT carried here (Issue #236
      Step C — it is closure-owned); the neutral slab closure
      (IdentityAngularClosure) supplies τ = 1, α = 0 and stamps the
      derived c on the CellVisit, but ``upstream_state.angular_upstream``
      is ``None`` for slab so the M-M contribution never engages.  Plus
      the always-populated ``chord_length``, ``mu``, ``volume``,
      ``abs_mu``.
    * **Sphere / cylinder**: physically-populated curvature fields
      from the dome recursion (the M-M angular weight is closure-owned,
      not carried here — Issue #236 Step C).

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
    (:class:`~orpheus.transport.spatial.scheme.CellVisit` packs the
    geometric :class:`StreamingTerms` together with the
    sweep-resolved :attr:`face_area_downstream`).  This module is
    geometry-only precisely because the sweep-direction concept is
    SN's alone — MoC, CP, diffusion and MC have no sweep.  (⛔ That
    separation is the true claim; the clause "and is reusable by future
    MoC / CP / diffusion" stood here until 2026-08-27 and is withdrawn.
    A geometry-only module is not thereby a shared one — see the ⛔
    block in the module docstring.)

    The trailing ``volume`` and ``abs_mu`` fields are populated by
    **all three factories** so that a downstream
    :class:`~orpheus.transport.spatial.scheme.DiscretizationScheme` strategy
    receives a self-contained per-cell, per-direction packet and need
    not reach back into ``SNMesh`` or the ``Quadrature``.  Every
    surviving curvature field is populated for **every** geometry
    (Step 2.5; slab carries neutral values) as a required ``float``
    (Pattern 4), so cell-update strategies consume the same packet
    regardless of chart — geometry is **data, not control-flow**.  Slab
    vs curvilinear is discriminated downstream by
    ``upstream_state.angular_upstream is None`` (slab has no half-angle
    state), NOT by any field on this geometry packet — Issue #236 Step C
    retired the former M-M ``alpha_in`` / ``alpha_out`` / ``tau_mm``
    fields (the angular closure now owns that data).
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

    # Issue #236 Step C: the Morel–Montry ``alpha_in`` / ``alpha_out`` /
    # ``tau_mm`` are NO LONGER carried on the per-cell geometry packet.  The
    # angular weight :math:`\\tau` is an angular-scheme property owned by the
    # MorelMontryAngularSweep angular closure; the derived constants
    # ``c_in`` / ``c_out`` and ``tau`` are stamped on
    # :class:`~orpheus.transport.spatial.scheme.CellVisit`.  The α-dome itself
    # survives on :class:`ReducedStreamingOperator` (``alpha_half`` /
    # ``alpha_per_level``).

    # ``mu_start`` retired 2026-08-26.  It was the middle link of a
    # three-link dead chain: AngularRedistribution.mu_start_per_level ->
    # StreamingTerms.mu_start -> StreamingCoefficientCache.mu_start -> nothing.
    # [M] the terminal had ZERO readers of any kind (no attribute access,
    # no getattr by name), so this field's only production consumer was
    # the WRITE into it -- while its own docstring claimed
    # ``MorelMontryAngularSweep`` consumed it.  The closure does consume a
    # starting direction; it reads the OWNER
    # (``AngularRedistribution.mu_start_per_level``), which is why the
    # claim read as true.  The chain was ERR-058's threading; the un-weld
    # gave the datum one owner and the thread became dead weight.

    volume: float
    """Cell volume :math:`V_i`.

    Populated by all three factories from
    :attr:`~orpheus.geometry.mesh.Mesh1D.volumes`: slab uses
    ``mesh.volumes[i]`` (which equals ``widths[i]`` for unit
    cross-section in 1-D Cartesian); sphere uses
    :math:`\\tfrac{4}{3}\\pi(r_{i+1}^3 - r_i^3)`; cylinder uses
    :math:`\\pi(r_{i+1}^2 - r_i^2)` (per unit axial length).  Carried
    on the streaming-terms packet so that the
    :class:`~orpheus.transport.spatial.scheme.DiscretizationScheme` cell-update
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
    :attr:`~orpheus.transport.spatial.scheme.CellResult.outgoing_spatial_flux`
    set to ``None`` (no radial face flow on the cell).
    """


# ═══════════════════════════════════════════════════════════════════════
# ReducedStreamingOperator — the primitive
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class ReducedStreamingOperator:
    """Connection-coefficient data for a curvilinear (or slab) sweep.

    Carries all geometry- and quadrature-dependent precomputed arrays
    that a curvilinear SN sweep needs to march through the angular
    redistribution.  In Cardinal Rule 2 framing, this primitive holds
    the data the **same connection-coefficient operator** demands, in
    whichever coordinate chart the consumer's geometry lives in — the
    point being that ONE object serves the sphere and the cylinder, not
    that several solver families consume it (they do not; see the module
    docstring's ⛔ block).

    Three factory functions construct concrete instances:

    * :func:`slab_streaming` — Cartesian 1-D, no curvature.
    * :func:`spherical_streaming` — 1-D spherical, with
      :math:`\\alpha_{n+1/2}` dome recursion (Hébert Eqs. 3.423-3.424,
      after Lathrop & Carlson 1966 — ⛔ NOT "Bailey 2009 Eq. 50", which
      this line claimed until 2026-08-26 while the module header 490
      lines above already condemned exactly that citation).
    * :func:`cylindrical_streaming` — 1-D cylindrical, with
      per-:math:`\\mu`-level :math:`\\alpha` and :math:`\\tau_{mm}`
      structures.

    Output shape contract.  The Morel–Montry τ is owned by the angular
    closure (Issue #236 Step C) and the α-dome / starting direction by
    :class:`AngularRedistribution` (2026-08-26, the un-weld arc's Phase
    B), so what remains here is the SPATIAL chart data::

        slab:        face_areas == None,
                     delta_A    == None.
        sphere:      face_areas (nx+1,),
                     delta_A    (nx,).
        cylinder:    face_areas (nx+1,),
                     delta_A    (nx,).

    ⭐ ``redist_dAw`` retired with them: it was the *fused product*
    ``ΔA ⊗ 1/w`` of a geometric factor and a quadrature factor, cached on
    the geometry object and read by two consumers that each wanted a
    different one of the two.  Both now form what they need from
    :attr:`delta_A` and the measure's weights (Pattern 2 — the stored
    product was a second spelling of a quantity neither side owned).

    Attributes
    ----------
    mesh :
        The :class:`~orpheus.geometry.mesh.Mesh1D` this operator was
        built from.  Held by reference; not copied.  **It is also where
        the chart lives** — see the note below.

    Note
    ----
    A ``coord: CoordSystem`` field was retired here on 2026-08-27 (P4.1a).
    It was a copy of ``mesh.coord``, and the copy was redundant *by
    construction* rather than merely in practice: each of the three
    factories below refuses a mesh whose chart differs from the literal
    it then stored (``slab_streaming`` raises unless
    ``mesh.coord is CARTESIAN``, and so on for the other two), so no
    reachable state could ever have made ``op.coord`` and
    ``op.mesh.coord`` disagree.

    Retiring it is what lets **L2 stop reaching through this bundle**:
    ``transport/radial_characteristic_field.py`` read
    ``mesh.reduced.coord`` for the single fact
    ``coord is CYLINDRICAL``, which the ``SNMesh`` it already held
    answers directly.  That was ``transport``'s only read of
    ``.reduced``; with it gone, this module has no L2 consumer left to
    break when it dissolves.

    Two flags -- ``requires_upstream_angular_state`` and
    ``angular_marching_axis`` -- were retired here on 2026-08-26.  Both
    were exactly ``coord is not CoordSystem.CARTESIAN``, and both had
    **zero production readers**: the concept is spelled twice already,
    by ``upstream_state.angular_upstream is None`` (the gate the DD and
    LD cell bodies actually branch on) and by ``SNMesh.is_cartesian``.
    Their 12 test assertions each sat one line below an assertion on
    ``coord`` that already pinned the same fact.
    """

    mesh: Mesh1D
    angular: AngularRedistribution
    """The angular factor — the dome and the starting direction, produced
    once by :func:`angular_redistribution` and shared with the angular
    closure.  Non-optional on every chart: Cartesian carries the NEUTRAL
    element (zero dome, diameter-ray start), which is what let the
    per-coordinate ``Optional`` union die."""

    # Common (curvilinear)
    face_areas: np.ndarray | None = None
    delta_A: np.ndarray | None = None

    # Issue #236 Step C retired ``tau_mm`` / ``tau_mm_per_level``: the
    # Morel–Montry angular weight is owned by the angular closure
    # (``morel_montry_tau_per_level``), not by the streaming geometry.
    # The α-dome went the same way at the 2026-08-26 un-weld (Phase B) —
    # it is the ANGULAR factor, and now lives on :attr:`angular` beside
    # the starting direction, with ONE producer.

    # ``_quadrature`` -- a SECOND reference to the same measure the angular
    # factor already holds (``_quadrature is angular.quadrature`` was True) --
    # was retired 2026-08-26.  ``streaming_terms`` reads
    # ``self.angular.quadrature``, which is NON-optional, so the ``| None``
    # union and its narrowing ``assert`` went with it (Pattern 4: the twin
    # existed only because the angular factor had no owner).

    # ── Per-direction extraction ───────────────────────────────────

    @cached_property
    def redistribution_pairing(self) -> np.ndarray:
        r"""The SPATIAL factor of the redistribution operator, ``(nx, n_mom, n_thread)``.

        The other half of the tensor product
        :math:`\mathcal{R} = R_{\rm spatial} \otimes A_{\rm angular}` —
        what the angular closure is handed so that neither side owns the
        product (:class:`AngularRedistribution` is the angular half).

        **The two axes are the two bases being paired**, and they are
        independent:

        * ``n_mom`` — how many spatial moments the SCHEME carries
          (diamond difference: 1; a linear-discontinuous cell: 2);
        * ``n_thread`` — how much of the spatial representation the
          ANGULAR device propagates through its half-angle recurrence.

        so the general entry is the one-measure-down pairing
        :math:`R_{kj} = \int b_k^{\rm scheme}\, b_j^{\rm thread}\,
        r\,\mathrm{d}r` — rectangular whenever the two differ.  Both
        published families are realizable in it: closing the angular
        index per spatial moment (Adams--Martin 1992 App. A) gives a
        square Gram, while closing it on the cell AVERAGE only (ONETRAN,
        Hill 1975 Eq. 32) gives a rank-1 column.  ORPHEUS must pick one
        explicitly; the axes exist so the pick is expressible.

        Today every shipped scheme is single-moment, so this is
        ``(nx, 1, 1)`` and its single entry is :math:`\Delta A_i` —
        exactly what the retired ``redist_dAw`` cache held, before the
        per-ordinate :math:`1/w_n` was fused onto it.

        Cartesian returns the neutral element (zeros): a slab's face
        areas are constant, so :math:`\Delta A \equiv 0` is the physical
        value, not a placeholder.

        **Cached for IDENTITY, not for speed.**  ``delta_A[:, None, None]``
        is a view, so recomputing it costs nothing — but a plain
        ``@property`` mints a fresh object per access (`[M]` ``g1 is not
        g2`` was True), and an object with no stable identity cannot
        participate in the two-strata cache the sweep is built on
        (``sn/sweep/cache.py``: *geometry × quadrature* surviving
        :math:`\Sigma_t` rebinds).  Since the pairing is a mesh-bound S1
        quantity, one value per operator is the honest lifetime.

        ⚠ Safe because this class is de-facto immutable — `[M]` nothing in
        ``orpheus/`` or ``tests/`` assigns to a constructed instance's
        fields.  It is ``@dataclass`` rather than ``@dataclass(frozen=True)``
        (unlike its sibling :class:`AngularRedistribution`), so that
        immutability is currently an audited property rather than a
        constructed one.  Freezing is the right end state and is deferred
        to the re-home, where the class definition is being touched anyway
        — it is not free: ``frozen=True`` synthesises ``__hash__``, and
        these fields hold ndarrays, so hashing would raise.  That wants a
        decision (``eq=False``? per-field ``compare=False``?), not a
        drive-by.
        """
        nx = int(self.mesh.widths.size)
        if self.delta_A is None:
            return np.zeros((nx, 1, 1))
        return np.asarray(self.delta_A, dtype=float)[:, None, None]

    def _weight_of(self, global_ordinate: int) -> float:
        r"""The measure's weight :math:`w_n` for one GLOBAL ordinate index.

        The other factor of the retired ``redist_dAw`` cache.  Callers form
        :math:`\Delta A_i / w_n` from :attr:`delta_A` and this, rather than
        reading a stored product that fused a geometric with a quadrature
        quantity (Pattern 2 — neither side owned the fusion).
        """
        return float(np.asarray(self.angular.quadrature.weights)[global_ordinate])

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
        :class:`orpheus.transport.spatial.scheme.CellVisit`.
        """
        chord = float(self.mesh.widths[cell_idx])
        # Common per-(cell, direction) volume.
        # mesh.volumes returns shape (N,) for 1-D meshes.
        volume = float(self.mesh.volumes[cell_idx])

        if self.mesh.coord is CoordSystem.CARTESIAN:
            # Slab — neutral curvature values populate the curvilinear
            # fields so cell_balance_terms_unified can consume the
            # packet without geometry dispatch (Issue #196 Phase G
            # Step 2.5).  Slab carries:
            #   face_area_inner = face_area_outer = 1.0  (so A_total =
            #       A_inner + A_outer = 2, and 2*|μ|*A_down = 2|μ|·1
            #       reproduces the slab denominator's "2|μ|" term);
            #   delta_A_over_w = 0.0  (no curvature redistribution).
            # The Morel–Montry α / τ are NO LONGER packed here (Issue #236
            # Step C): the neutral slab closure (IdentityAngularClosure)
            # supplies τ = 1, α = 0 and stamps the derived c on CellVisit.
            # ``volume == chord`` is already true for slab (unit
            # cross-section in 1-D Cartesian).
            mu_n = float(self.angular.quadrature.mu_x[direction_idx])
            return StreamingTerms(
                chord_length=chord,
                mu=mu_n,
                face_area_inner=1.0,
                face_area_outer=1.0,
                delta_A_over_w=0.0,
                volume=volume,
                abs_mu=abs(mu_n),
            )

        if self.mesh.coord is CoordSystem.SPHERICAL:
            assert self.face_areas is not None
            assert self.delta_A is not None
            # Sphere: ``direction_idx`` IS the global ordinate index.
            # The Morel–Montry α / τ are NO LONGER packed here (Issue #236
            # Step C): the angular closure owns τ and the derived c
            # (stamped on CellVisit); this packet carries geometry only.
            mu_n = float(self.angular.quadrature.mu_x[direction_idx])
            return StreamingTerms(
                chord_length=chord,
                mu=mu_n,
                face_area_inner=float(self.face_areas[cell_idx]),
                face_area_outer=float(self.face_areas[cell_idx + 1]),
                delta_A_over_w=float(
                    self.delta_A[cell_idx]
                    / self._weight_of(direction_idx)
                ),
                volume=volume,
                abs_mu=abs(mu_n),
            )

        if self.mesh.coord is CoordSystem.CYLINDRICAL:
            if mu_level_idx is None:
                raise ValueError(
                    "cylindrical streaming_terms() requires mu_level_idx "
                    "(which μ-level the direction_idx belongs to)."
                )
            assert self.face_areas is not None
            assert self.delta_A is not None
            # Cylinder: ``direction_idx`` is the within-level azimuthal
            # index; the global ordinate is read through
            # ``level_indices``.  ``mu_x[global_n]`` carries η (the
            # radial direction cosine).  The Morel–Montry α / τ are NO
            # LONGER packed here (Issue #236 Step C): the angular closure
            # owns τ and the derived c, stamped on CellVisit; this packet
            # carries geometry only.  (There is no clamp on either arm —
            # the cylinder [1/2, 1] absorber retired at Q5.6.4.)
            level_indices = self.angular.quadrature.level_indices
            global_n = int(level_indices[mu_level_idx][direction_idx])
            eta_n = float(self.angular.quadrature.eta[global_n])
            return StreamingTerms(
                chord_length=chord,
                mu=eta_n,
                face_area_inner=float(self.face_areas[cell_idx]),
                face_area_outer=float(self.face_areas[cell_idx + 1]),
                delta_A_over_w=float(
                    self.delta_A[cell_idx] / self._weight_of(global_n)
                ),
                volume=volume,
                abs_mu=abs(eta_n),
            )

        raise ValueError(  # pragma: no cover — exhaustive match above
            f"Unknown coord system: {self.mesh.coord!r}"
        )


# ═══════════════════════════════════════════════════════════════════════
# The α-dome — ONE recursion, and its admission contract
# ═══════════════════════════════════════════════════════════════════════
#
# Both curvilinear factories below need the Lathrop--Carlson angular
# redistribution coefficients, and until 2026-08-12 each spelled the
# recursion out in its own body (sphere and cylinder, identical arithmetic
# — the sphere IS the cylinder's single-level case).  That twin path is
# why the closure contract was enforced on ONE arm only, and why the one
# check that did exist was a bare ``assert`` (see below).  Cardinal Rule 2:
# one concept, one body.
#
# The RECURSION and the CONTRACT are deliberately two functions, not one.
# ``orpheus.derivations.discrete.sn.angular_differencing`` explores
# hypothetical and deliberately-non-closing measures (its P0/P4 predicate
# ladder) and must be able to run the recursion on them; a guard welded
# into the recursion would make that analysis unspellable.  The contract
# belongs at the production ADMISSION point — the factories — where a
# non-closing measure is genuinely illegal.


def alpha_dome(
    cosines: np.ndarray,
    weights: np.ndarray,
) -> np.ndarray:
    r"""The Lathrop--Carlson angular redistribution recursion.

    .. math::
       \alpha_{m+1/2} \;=\; \alpha_{m-1/2} \;-\; w_m\,\mu_m ,
       \qquad \alpha_{1/2} = 0

    in the ORPHEUS factor-of-2-absorbed normalization (Hébert 2009
    §3.9.4 Eqs. 3.423-3.424 for the sphere; his cylinder is §3.9.3, where
    the :math:`\eta_{p,q\pm1/2}` construction is credited to Alcouffe &
    O'Dell 1986 — not local, not read).  The recursion itself traces to
    Lathrop & Carlson, *J. Comp. Phys.* **1**:173 (1966) — likewise not
    local, not read; cited by Reed & Lathrop 1970 (their ref. 7) as "a
    requirement commonly invoked to define the α coefficients".

    ⭐ **The recursion is strictly ONE-SIDED**: :math:`\alpha_{1/2} = 0` is
    an axiom of the construction, but the closure at the far end,
    :math:`\alpha_{M+1/2} = 0`, is **not** — it is a *consequence* of the
    measure's antisymmetry, i.e. a genuine contract on every admitted
    quadrature, and nothing here enforces it.  Callers that require a
    physically admissible dome must run
    :func:`_assert_alpha_dome_closes` on the result.

    Parameters
    ----------
    cosines :
        Shape ``(M,)`` — the marching direction cosine per ordinate,
        **ordered by increasing cosine**: :math:`\mu` on the sphere, the
        radial :math:`\eta` within one :math:`\mu`-level on the cylinder.
    weights :
        Shape ``(M,)`` — the matching quadrature weights (the level's
        azimuthal share on the cylinder).

    Returns
    -------
    np.ndarray
        Shape ``(M+1,)``, the dome on the angular EDGES, with
        ``alpha[0] == 0.0`` exactly.
    """
    cosines = np.asarray(cosines, dtype=float)
    weights = np.asarray(weights, dtype=float)
    alpha = np.zeros(cosines.size + 1)
    for m in range(cosines.size):
        alpha[m + 1] = alpha[m] - weights[m] * cosines[m]
    return alpha


#: Closure tolerance for :func:`_assert_alpha_dome_closes`.  `[M]` 2026-08-12
#: the worst residual over the shipped rules is ``2.78e-16``
#: (``folded_product(4, 32)``); Gauss--Legendre reads ``5.6e-17`` (N=4)
#: … ``8.2e-17`` (N=64), i.e. ≈ 1 ULP of the dome peak and not drifting
#: with N.  ``1e-12`` therefore clears every shipped rule by ~4 orders
#: while still refusing any real antisymmetry violation.
_ALPHA_CLOSURE_ATOL = 1e-12


def _assert_alpha_dome_closes(
    alpha: np.ndarray,
    *,
    coord: CoordSystem,
    level: int | None = None,
) -> None:
    r"""The admission contract :math:`\alpha_{M+1/2} = 0`.

    A dome that does not close means the measure's marching cosines are
    not antisymmetric about the level's centre — a mis-ordered, truncated,
    duplicated or simply inadmissible rule.  The consequence is not a
    small error, and it is not local: the closure's coefficient
    precompute reads :math:`\alpha_{M+1/2}` as the LAST entry of its
    ``alpha_out`` slice, so it lands in **both** cell-balance
    coefficients of the final ordinate —
    ``c_out[M-1] = α_{M+1/2}/τ`` (a denominator term) and
    ``c_in[M-1] = (1-τ)/τ·α_{M+1/2} + α_{M-1/2}`` (an
    upstream-numerator term); see
    :class:`~orpheus.sn.angular.closure.MorelMontryAngularSweep`'s
    ``_c_in_per_level`` / ``_c_out_per_level`` precompute.  A closing dome
    makes that denominator term vanish — angular redistribution stops at
    the level's top edge, which is the whole point of the edge.  A
    non-zero value instead redistributes flux past it, into nothing.
    Refuse it here.

    ⛔ **This guard replaced a bare ``assert`` on 2026-08-12, and the
    ``assert`` was INERT where it mattered.** The sphere arm carried
    ``assert abs(alpha[N]) < 1e-12``; the cylinder arm carried nothing at
    all.  The project's canonical runner is ``python -O -m pytest``
    (``.claude/rules/vv-testing.md``), and ``-O`` sets ``__debug__ =
    False`` and strips every ``assert`` statement — so on the one arm that
    had a check, that check did not run in the canonical test suite.
    `[M]` measured on the verbatim recursion: a measure closing at
    ``alpha[N] = +0.2000`` is REFUSED under plain ``python`` and
    **ACCEPTED** under ``python -O``.  A contract that a compiler flag can
    remove is not a contract; this raises.

    Sibling of
    :func:`~orpheus.sn.angular.closure._assert_tau_within_unit_interval`
    — the same shape (an admission predicate on the angular scheme that
    refuses upstream rather than absorbing downstream), for the same
    reason.
    """
    residual = float(abs(alpha[-1]))
    if not residual < _ALPHA_CLOSURE_ATOL:
        where = "" if level is None else f" on level {level}"
        raise ValueError(
            f"the alpha dome does not close{where}: "
            f"alpha[M+1/2] = {float(alpha[-1])!r}, which exceeds the "
            f"{_ALPHA_CLOSURE_ATOL:g} admission tolerance for a "
            f"{coord.name} rule. The recursion "
            f"alpha_{{m+1/2}} = alpha_{{m-1/2}} - w_m*mu_m starts at "
            f"alpha_1/2 = 0 and can only return there if the marching "
            f"cosines are antisymmetric about the level's centre, so a "
            f"non-zero endpoint certifies an inadmissible measure — "
            f"mis-ordered, truncated, or duplicated ordinates upstream. "
            f"alpha[M+1/2] multiplies the last ordinate's outflow "
            f"half-face, so absorbing it here would leak flux out of the "
            f"level's top edge. Fix the quadrature upstream."
        )


# ═══════════════════════════════════════════════════════════════════════
# The angular factor — the redistribution structure of a measure in a chart
# ═══════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class AngularRedistribution:
    r"""The angular measure's redistribution structure on a coordinate chart.

    **The object exists because the redistribution operator FACTORS.**  The
    curvilinear angular-redistribution term of the transport operator is a
    tensor product

    .. math::

       \mathcal{R} \;=\; R_{\rm spatial} \;\otimes\;
                          A_{\rm angular}(\tau,\ \alpha,\ w)

    — a spatial pairing against an angular operator — and this type is the
    **member-independent half of the angular factor**: the :math:`\alpha`
    dome and the starting direction, which every angular-closure member
    shares.  What is NOT here is the member's own choice: :math:`\tau` and
    the derived :math:`c_{\rm in}` / :math:`c_{\rm out}` belong to the
    closure (Morel--Montry's barycentric weight, plain diamond's
    :math:`\tfrac12`, the neutral :math:`\tau \equiv 1`), and a shared
    object holding them would forbid a second member by construction.

    (Derived 2026-08-26, ``scratch/tau_under_ld_dip_analysis.md`` §3 and
    §7: the factorization is why the angular closure's weight cannot
    acquire spatial content, and the member-independent /
    member-dependent split is that memo's §7.1/§7.2 table.)

    **One owner, two consumers.**  Until this type existed, :math:`\alpha`
    and :math:`\mu_{\rm start}` were stored as six ``Optional`` fields on
    :class:`ReducedStreamingOperator` — a per-coordinate union discharged
    by ``assert`` — and were read both by the angular closure and by
    :meth:`ReducedStreamingOperator.streaming_terms`.  They are produced
    HERE, once, by :func:`angular_redistribution`.

    **Cartesian is the NEUTRAL element, not a special case.**  A slab has
    no curvature, so its dome is identically zero and its starting
    direction is the diameter ray; spelling those values rather than
    ``None`` is what lets the per-coordinate union die (Pattern 4 — the
    "no redistribution" state is not separately representable).

    **Sphere is the one-level case.**  Every field is
    per-:math:`\mu`-level; the sphere carries exactly one level.  The
    closure's own constructor already normalized to this shape internally,
    so the normal form is moved here rather than invented.

    Attributes
    ----------
    coord :
        The chart the redistribution was built on.
    alpha_per_level :
        The Lathrop--Carlson dome :math:`\alpha_{m\pm 1/2}` per level,
        each of shape ``(M_p + 1,)`` and closing at
        :math:`\alpha_{M+1/2} = 0` (:func:`_assert_alpha_dome_closes`).
    mu_start_per_level :
        The starting-direction angular edge per level — the direction the
        half-angle thread's seed flux lives at.  Sphere: ``-1.0`` (the
        Hébert §3.9.4 diameter ray).  Cylinder:
        :math:`-\sin\theta_p = -\sqrt{1-\xi_p^2}`.
    quadrature :
        The measure this was built from.  Held so a consumer that needs
        the weights, the level partition or the cosines does not have to
        be handed the quadrature separately alongside its redistribution.
    """

    coord: CoordSystem
    alpha_per_level: tuple[np.ndarray, ...]
    mu_start_per_level: tuple[float, ...]
    quadrature: AngularMeasure = field(repr=False)

    @property
    def n_levels(self) -> int:
        r"""Number of :math:`\mu`-levels (sphere and slab: ``1``)."""
        return len(self.alpha_per_level)


def angular_redistribution(
    quadrature: AngularMeasure,
    coord: CoordSystem,
) -> AngularRedistribution:
    r"""Build the :class:`AngularRedistribution` for a measure on a chart.

    THE single producer of :math:`\alpha` and :math:`\mu_{\rm start}`.
    Both are functions of ``(quadrature, coord)`` alone — no cell, no
    mesh, no material enters — which is measured, not assumed: every
    value this returns is **bit-identical** to what the three streaming
    factories stored before the carve (`[M]` 2026-08-26, both geometries,
    all levels, ``max|delta| = 0.000e+00``).

    Per level, :math:`\alpha` is :func:`alpha_dome` over that level's
    cosines and weights, and the admission contract
    :math:`\alpha_{M+1/2} = 0` is checked here — at the one site that
    produces it, rather than on whichever arm remembered to.
    """
    if coord is CoordSystem.CARTESIAN:
        # The neutral element: no curvature ⇒ no redistribution.  The dome
        # is zero (so the closure's c_in = c_out = 0 fall out rather than
        # being special-cased) and the starting direction is the diameter.
        n = int(np.asarray(quadrature.mu_x).size)
        return AngularRedistribution(
            coord=coord,
            alpha_per_level=(np.zeros(n + 1),),
            mu_start_per_level=(-1.0,),
            quadrature=quadrature,
        )

    mu_x = np.asarray(quadrature.mu_x)
    weights = np.asarray(quadrature.weights)

    if coord is CoordSystem.SPHERICAL:
        # One level: every ordinate marches on the single polar dome, and
        # the seed rides the diameter ray (Hébert §3.9.4).
        alpha = alpha_dome(mu_x, weights)
        _assert_alpha_dome_closes(alpha, coord=coord)
        return AngularRedistribution(
            coord=coord,
            alpha_per_level=(alpha,),
            mu_start_per_level=(-1.0,),
            quadrature=quadrature,
        )

    if coord is CoordSystem.CYLINDRICAL:
        # One dome per μ-level; the level's seed rides its most-inward
        # azimuthal edge, η_{1/2} = -sin θ_p = -sqrt(1 - ξ_p²).
        alphas: list[np.ndarray] = []
        starts: list[float] = []
        xi = np.asarray(quadrature.mu_z)
        for p, level_idx in enumerate(quadrature.level_indices):
            level_idx = np.asarray(level_idx)
            alpha = alpha_dome(mu_x[level_idx], weights[level_idx])
            _assert_alpha_dome_closes(alpha, coord=coord, level=p)
            alphas.append(alpha)
            starts.append(-float(np.sqrt(1.0 - float(xi[level_idx[0]]) ** 2)))
        return AngularRedistribution(
            coord=coord,
            alpha_per_level=tuple(alphas),
            mu_start_per_level=tuple(starts),
            quadrature=quadrature,
        )

    raise ValueError(  # pragma: no cover — exhaustive match above
        f"Unknown coord system: {coord!r}"
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
    vanish.  The SPATIAL chart is genuinely absent (``face_areas`` and
    ``delta_A`` stay ``None``), but the ANGULAR factor is present and
    carries the NEUTRAL element: a dome that is identically zero and a
    starting direction on the diameter ray.  That is what lets the
    angular closure's ``c_in = c_out = 0`` fall out of the general body
    instead of being special-cased.

    (Until 2026-08-26 this paragraph said "all ``alpha_*`` and
    ``redist_dAw`` arrays remain ``None``" and quoted two flags.  Both
    halves went with the un-weld: the ``Optional`` angular fields became
    the neutral element, and the flags were exactly
    ``coord is not CARTESIAN`` with no production reader.)

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
        mesh=mesh,
        angular=angular_redistribution(angular_measure, CoordSystem.CARTESIAN),
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
    recursion, in the ORPHEUS factor-of-2-absorbed normalization),
    producing arrays bit-identical to the retired
    ``SNMesh._setup_spherical`` it replaced.  It does **not** produce the
    Morel--Montry angular weight :math:`\tau` — that moved to the angular
    closure at Issue #236 Step C (see the body comment below).

    The :math:`\Delta A/w` geometry factor (Cardinal Rule 2 — chart
    data, so it belongs to the geometry layer rather than to a
    solver-side mesh class) is precomputed here once.  ⛔ This line read
    "common to SN/MoC/CP" until 2026-08-27; see the module docstring's
    ⛔ block.

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

    # Cell face areas: A_{i+1/2} = 4πr² at each edge — sourced from the
    # mesh, which routes through coord.compute_areas_1d().
    face_areas = mesh.areas  # (nx+1,)

    # Cell face-area differences: ΔA_i = A_{i+1/2} − A_{i-1/2}
    delta_A = face_areas[1:] - face_areas[:-1]

    # The α-dome and the starting direction are NOT produced here — they
    # are the ANGULAR factor, built once by :func:`angular_redistribution`
    # (which also runs the α_{M+1/2} = 0 admission contract, at the one
    # site that produces the dome).

    # Morel–Montry angular weight τ is NO LONGER produced here (Issue #236
    # Step C): τ is an angular-scheme property (a function of (μ, w) only),
    # owned by the MorelMontryAngularSweep angular closure
    # (``morel_montry_tau_per_level``) and stamped on each CellVisit.  This
    # factory keeps the GEOMETRY data only (face areas, α-dome, redistribution
    # factor).
    #
    # The Hébert §3.9.4 starting direction μ_{1/2} = -1.0 of the (single)
    # M-M level — the direction the half-angle thread's seed flux lives at —
    # is carried by the ANGULAR factor as
    # ``angular.mu_start_per_level[0]``, not by this packet.

    return ReducedStreamingOperator(
        mesh=mesh,
        face_areas=face_areas,
        delta_A=delta_A,
        angular=angular_redistribution(angular_measure, CoordSystem.SPHERICAL),
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

    Requires ``angular_measure`` to carry per-:math:`\mu`-level
    structure: the guard below rejects a quadrature whose
    ``level_structure`` is ``None``, and the math then reads
    ``level_indices`` (a list of index arrays, one per level).
    ⛔ **This paragraph used to name** :meth:`Quadrature.level_symmetric`
    **and** :meth:`Quadrature.product` **as "the two cylindrical-compatible
    factories" — corrected 2026-08-26: both are REFUSED** by
    :func:`~orpheus.sn.angular.closure.assert_carrying_quadrature`,
    which this factory's only caller runs twelve lines after calling it (an
    ABS_MU_Z level carries both hemispheres, so η is degenerate on the
    most-inward node and the seed slot is a rank duplicate).  The rule that
    IS admitted is :meth:`Quadrature.folded_product
    <orpheus.numerics.quadrature.Quadrature.folded_product>`, the σ_y
    quotient — which appeared in this module only inside an error string.
    The slab (``gauss_legendre``) and pure-sphere (``lebedev``) factories
    carry no ``LevelStructure`` side-channel. These are named classmethod
    factories on the single ``Quadrature`` type — the per-family
    adapter classes this docstring used to name were retired in R-1
    Phase A detour-C. The quadrature must also expose ``mu_z`` (axial
    direction cosine, used to derive :math:`\sin\theta` for the
    level's azimuthal extent).

    Output is bit-identical to the retired
    ``SNMesh._setup_cylindrical`` it replaced.

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
            "structure (e.g. Quadrature.folded_product); "
            "slab quadratures (Quadrature.gauss_legendre) and pure "
            "sphere cubatures (Quadrature.lebedev) carry no "
            "LevelStructure side-channel."
        )

    face_areas = mesh.areas  # (nx+1,)
    delta_A = face_areas[1:] - face_areas[:-1]

    # Per-level azimuthal redistribution coefficients — the SAME
    # :func:`alpha_dome` recursion the sphere runs, once per μ-level, on
    # the level's radial cosine η and its azimuthal weight share.
    # Ordinates are ordered by increasing η within each level.
    #
    # ⭐ Each level is held to the closure contract INDIVIDUALLY (the arm
    # carried no check at all until 2026-08-12 — see
    # :func:`_assert_alpha_dome_closes`).  A per-level raise is what makes
    # the failure LOCATABLE: a whole-measure check would only say
    # "somewhere", and per-level costs nothing.
    #
    # ⛔ The justification here used to be "a level-symmetric rule can
    # lose antisymmetry on one level alone".  That is unreachable prose:
    # a level-symmetric rule cannot get this far, because an ABS_MU_Z
    # level carries both hemispheres, so η is 4-fold degenerate and
    # `assert_carrying_quadrature` refuses on `eta_level[0] ==
    # eta_level[1]` (`[M]` every order S2..S12, both coord systems;
    # #336 tracks the refuse-or-reduce design for the sphere arm).  The
    # per-level form is still right — locatability is reason enough —
    # but it is not held up by that case.
    # Per-level starting-direction angular edge η_{1/2} = -sin θ_p — the
    # most-inward azimuthal edge of level p, the direction the half-angle
    # thread's seed flux lives at.  The Morel–Montry angular weight τ is NO
    # LONGER produced here (Issue #236 Step C): τ is an angular-scheme
    # property owned by the MorelMontryAngularSweep angular closure
    # (``morel_montry_tau_per_level``) and stamped on each CellVisit.  This
    # factory keeps the GEOMETRY data only.  (The cylinder [½, 1] absorber
    # that used to live alongside that producer retired at Q5.6.4 — there
    # is one τ, unclamped, on both arms.)
    return ReducedStreamingOperator(
        mesh=mesh,
        face_areas=face_areas,
        delta_A=delta_A,
        angular=angular_redistribution(angular_measure, CoordSystem.CYLINDRICAL),
    )


__all__ = [
    "AngularMeasure",
    "ReducedStreamingOperator",
    "StreamingTerms",
    "cylindrical_streaming",
    "slab_streaming",
    "spherical_streaming",
]
