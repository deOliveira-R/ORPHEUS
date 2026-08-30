r"""The chart connection — the coefficient primitive a curvilinear sweep marches.

:class:`ReducedStreamingOperator` is the per-``(mesh, quadrature)`` producer of
the numbers a curvilinear S\ :sub:`N` sweep needs at each cell and direction:
chord length, the two face areas, the :math:`\Delta A/w` redistribution factor,
and the radial cosine.  In differential-geometric language the spherical
:math:`(1-\mu^2)/r\,\partial_\mu` and cylindrical
:math:`-(1/r)\,\partial_\varphi(\xi\cdot)` terms are **the same connection
coefficient** on SO(3) read in two charts (polar-on-sphere vs.
azimuthal-on-cylinder), which is why one object serves both.

**Why this lives in ``sn/`` and not in ``geometry/``** — the un-weld arc's P4.4,
2026-08-28.  It lived in ``orpheus/geometry/reduced_operator.py`` until then, and
the measurement that moved it is the **layer test**:

    *A datum belongs to the layer that can define it without naming a method.
    Everything else is posing, and posing belongs to the method head that
    poses it.*

`[M]` exactly one datum in the old module survived that test — ``face_areas`` —
and it was a verbatim copy of :attr:`orpheus.geometry.mesh.Mesh1D.areas`, already
single-sourced in :func:`orpheus.geometry.coord.compute_areas_1d`.  Everything
else is posing: ``delta_A`` is ``np.diff(mesh.areas)`` with **zero** non-SN
consumers, and :class:`~orpheus.transport.spatial.scheme.StreamingTerms` carries
``mu``, ``abs_mu`` and :math:`\Delta A` divided by a **quadrature weight** — a
per-direction quantity over a quadrature weight is not geometry.  `[M]` the old
module was also an island inside its own package: every genuine geometry
primitive (``CoordSystem``, ``Mesh1D``, ``StructuredGeometry``, ``RigidMotion``)
had 1-4 intra-``geometry/`` consumers and this one had **0**.

✅ **The un-weld completed in three ordered moves** (order forced by the
layer contract): the connection operator and its three factories came here
first (P4.4); then :class:`~orpheus.sn.angular.redistribution.AngularMeasure`
and the five :math:`\alpha` symbols
(:func:`~orpheus.sn.angular.redistribution.alpha_dome`, its closure contract,
:class:`~orpheus.sn.angular.redistribution.AngularRedistribution`,
:func:`~orpheus.sn.angular.redistribution.angular_redistribution`) to
``sn/angular/redistribution.py`` (P4.2 — they could not move first: `[M]` the
three factories below CALL ``angular_redistribution`` at runtime, so moving
:math:`\alpha` while the factories sat in ``geometry/`` would have created
``geometry -> sn``, a declared ``FORBIDDEN_EDGES`` violation *and* a hard
circular import); finally
:class:`~orpheus.transport.spatial.scheme.StreamingTerms` to
:mod:`orpheus.transport.spatial.scheme` (P4.3, 2026-08-28), beside the
scheme contract that consumes it — and ``orpheus/geometry/reduced_operator.py``
was deleted.

The imports below run ``sn -> geometry`` and ``sn -> transport``, both legal
(``geometry`` is ``INPUT_PACKAGES``; ``sn -> transport`` is the established
direction; see ``tests/test_layer_imports.py``).

⚠ ``ReducedStreamingOperator`` is a **provisional name** — the object is not a
streaming operator and never was; it produces the coefficients one marches.
Renaming it is P3c, deliberately sequenced after P4 so the name describes what
remains once ``streaming_terms`` has moved out of it.

Provenance of the primitive
===========================
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

See also
========

* :class:`orpheus.sn.mesh.augmented_mesh.SNMesh` — the **only** constructor and
  holder: ``SNMesh.__init__`` calls :func:`slab_streaming` /
  :func:`spherical_streaming` / :func:`cylindrical_streaming` and stores the
  result as its ``reduced`` attribute.  That is why this module sits in
  ``sn/mesh/``.
* :class:`~orpheus.sn.operators.streaming.StreamingOperator` — a *different*
  object despite the similar name: the operator-algebra leaf :math:`L` in
  :math:`A_{\rm wg} = L + C - S`.  It consumes this primitive through the
  loss-representation walk.
* :mod:`orpheus.sn.angular.redistribution` — the angular-measure adapter and
  the :math:`\alpha` dome with its closure contract, carrying the full
  :math:`\alpha`/:math:`\tau` mathematical content and literature references.
* :class:`~orpheus.transport.spatial.scheme.StreamingTerms` — the packet this
  operator's :meth:`~ReducedStreamingOperator.streaming_terms` produces, now
  beside the scheme contract that consumes it.
* :doc:`/theory/foundations/structured_geometry` — the architecture page; see
  "Connection coefficients (reduced streaming operator)".
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property

import numpy as np

from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.mesh import Mesh1D
from orpheus.transport.spatial.scheme import StreamingTerms
from orpheus.sn.angular.redistribution import (
    AngularMeasure,
    AngularRedistribution,
    angular_redistribution,
)


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

        every chart:  face_areas (nx+1,)   == mesh.areas
                      delta_A    (nx,)     == diff(mesh.areas)

    ⛔ Until 2026-08-27 (P4.1b) this block read ``slab: face_areas ==
    None, delta_A == None`` and the two were stored fields.  They are
    derived now, and the slab is not a special case: its ``face_areas``
    is ``ones(nx+1)`` (``compute_areas_1d`` returns a real unit
    cross-section on CARTESIAN) and its ``delta_A`` is ``zeros(nx)`` —
    no area change, which is exactly what "no curvature" means.

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
    **zero production readers**: the concept was spelled twice already —
    at the time by ``upstream_state.angular_upstream is None`` (the gate
    the DD and LD cell bodies then branched on; itself retired at P4.9a
    for assembled angular arguments) and by ``SNMesh.is_cartesian``.
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

    # ── The spatial chart — DERIVED, not stored (P4.1b, 2026-08-27) ──
    #
    # ``face_areas`` and ``delta_A`` were ``np.ndarray | None`` fields,
    # populated by the two curvilinear factories and left ``None`` by the
    # slab.  They are functions of the mesh, so they are read from it.
    #
    # ⭐ Storing them was a Pattern-2 triplicate waiting to happen: the
    # sphere and cylinder factories each spelled ``face_areas = mesh.areas;
    # delta_A = diff(face_areas)``, and "populate the slab too" would have
    # written the same two lines a third time.  Deriving them writes it
    # once and no factory computes anything.
    #
    # ⭐ The slab's ``None`` was never "no value" — it was the value nobody
    # asked for.  ``compute_areas_1d`` returns a real unit cross-section on
    # CARTESIAN, so ``face_areas`` is ``ones(nx+1)`` and ``delta_A`` is
    # ``zeros(nx)``: a slab has no area change, which IS what "no
    # curvature" means.  Spelling it is a STRICTLY STRONGER claim than
    # ``None`` (the same upgrade P1 made on the ANGULAR half — the zero
    # dome), and it is what collapses ``streaming_terms``' three arms.

    @property
    def face_areas(self) -> np.ndarray:
        r"""Face areas :math:`A_{i\pm1/2}`, ``(nx+1,)`` — the mesh's own.

        One attribute hop: :attr:`Mesh1D.areas` is computed eagerly in
        ``__post_init__`` (the mesh is frozen), so this is the same array
        object on every read, not a recompute and not a copy.
        """
        return self.mesh.areas

    @cached_property
    def delta_A(self) -> np.ndarray:
        r"""The connection integral :math:`\Delta A_i = \int \nabla\cdot\hat{e}_r \, dV`,
        ``(nx,)`` — realized as the face-area difference.

        ⚠ **Cached, and it must be.** :meth:`streaming_terms` is called per
        ``(cell, direction)``, so a plain ``@property`` would recompute an
        ``nx``-element :func:`numpy.diff` inside the sweep's hot loop.

        ⭐ This is the first step of the ruled dissolution of ``delta_A``
        into the redistribution pairing :math:`R` (2026-08-27): ``ΔA`` is
        R's rank-1 realization, not a second object, so R's producer will
        derive the connection integral from ``mesh.areas`` exactly as this
        does.  Deriving it here means what later moves is a derivation
        rather than an array.
        """
        return np.diff(self.mesh.areas)

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
        # ⛔ An ``if self.delta_A is None: return np.zeros((nx, 1, 1))``
        # branch stood here until P4.1b.  It was a value wearing a
        # conditional (Pattern 4): the slab's ``delta_A`` IS ``zeros(nx)``,
        # so the general body returns the same array the branch built.
        return np.asarray(self.delta_A, dtype=float)[:, None, None]

    def streaming_terms(
        self,
        cell_idx: int,
        direction_idx: int,
        mu_level_idx: int | None = None,
    ) -> StreamingTerms:
        """Pack the per-cell, per-direction streaming inputs (the scheme's evaluation point).

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
        # Common per-(cell, direction) volume.
        # mesh.volumes returns shape (N,) for 1-D meshes.
        volume = float(self.mesh.volumes[cell_idx])

        # ── ONE body, every chart (P4.1b, 2026-08-27) ────────────────
        #
        # ⛔ This was three arms.  The CARTESIAN one differed from the
        # SPHERICAL one only in three hardcoded literals —
        # ``face_area_inner = face_area_outer = 1.0`` and
        # ``delta_A_over_w = 0.0`` — which were hand-transcriptions of what
        # the spherical body computes: on a slab ``mesh.areas`` IS
        # ``ones(nx+1)`` and its difference IS ``zeros(nx)``.  Measured over
        # the 5x8 (cell x ordinate) grid with the spatial chart derived
        # rather than stored: 40 of 40 packets bit-identical, with a
        # perturbed ``delta_A`` detected as the positive control.  The slab
        # IS the sphere's zero-curvature case, exactly as the sphere is the
        # cylinder's single-level case (see the α-dome note below, whose
        # twin path retired 2026-08-12 on the same argument).
        #
        # What survives the collapse is NOT chart-dispatched arithmetic —
        # it is what ``direction_idx`` MEANS, which differs by caller:
        # the GLOBAL ordinate on slab and sphere, the WITHIN-LEVEL
        # azimuthal index on the cylinder.  One parameter, two contracts.
        #
        # ⚠ It cannot be unified by always going through ``level_indices``.
        # That reduces to ``direction_idx`` only when the level list is
        # ``arange(N)``, which holds by construction for a quadrature with
        # no ``LevelStructure`` — but NOT for every single-level rule:
        # measured over the 40-rule shipped registry, 9 are single-level
        # and ``level_symmetric(2)`` (N = 8) carries a PERMUTED list
        # ``[2,3,0,1,6,7,4,5]``.  Routing a slab or sphere through it would
        # silently re-index, and nothing in the factories forbids posing one
        # on that rule.
        if self.mesh.coord is CoordSystem.CYLINDRICAL:
            if mu_level_idx is None:
                raise ValueError(
                    "cylindrical streaming_terms() requires mu_level_idx "
                    "(which μ-level the direction_idx belongs to)."
                )
            level_indices = self.angular.quadrature.level_indices
            ordinate = int(level_indices[mu_level_idx][direction_idx])
        else:
            ordinate = direction_idx

        # The direction cosine along the SWEEP AXIS — the radius on a
        # curvilinear chart, x on a slab.  ``mu_x`` and ``eta`` are the same
        # accessor (both ``axis_cosines(0)``); ``mu_x``'s own docstring says
        # "the column index, not the name, is the actual semantic".  Neither
        # spelling is chart-neutral, which is a naming item for the
        # ``face_area_*`` family pass; the local name carries the meaning.
        radial_cosine = float(self.angular.quadrature.mu_x[ordinate])

        # The Morel–Montry α / τ are NOT packed here (Issue #236 Step C):
        # the angular closure owns τ and the derived c, and stamps them on
        # CellVisit.  This packet carries geometry only.  There is no clamp
        # on any chart — the cylinder [1/2, 1] absorber retired at Q5.6.4.
        # ΔA/w is NOT packed (P4.7): the fusion's owner is the angular
        # closure (``_dAw_per_level``, P4.9a); the scan cache interns the
        # strategy-side copy from ``delta_A`` and the measure's weights at
        # build — both ΔA/w formers read the weights directly, so no
        # per-ordinate weight accessor lives here (``_weight_of`` existed
        # until 2026-08-29 and was retired with zero callers).
        return StreamingTerms(
            face_area_inner=float(self.face_areas[cell_idx]),
            face_area_outer=float(self.face_areas[cell_idx + 1]),
            volume=volume,
            abs_mu=abs(radial_cosine),
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
    vanish.  The SPATIAL chart is the NEUTRAL one (``face_areas`` is the
    unit cross-section ``ones(nx+1)`` and ``delta_A`` its zero
    difference, both derived through the one generic body — P4.1a/b;
    ⛔ "stay ``None``" stood here until P4.7 and was false since that
    collapse), and the ANGULAR factor is present and
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

    # Cell face-area differences: ΔA_i = A_{i+1/2} − A_{i-1/2}

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
        angular=angular_redistribution(angular_measure, CoordSystem.CYLINDRICAL),
    )


__all__ = [
    "ReducedStreamingOperator",
    "cylindrical_streaming",
    "slab_streaming",
    "spherical_streaming",
]
