.. _wavefront-cochain:

==============================
The Interior Face-Flux Cochain
==============================

.. contents:: Contents
   :local:
   :depth: 2


.. Machine header — the ``nexus-meta`` schema for this page (PROVISIONAL).
.. This page was extracted verbatim from ``operator_algebra.rst`` (#231
.. Phase 3); content is unchanged. The schema is provisional pending a
.. full re-audit of the split corpus.

.. dropdown:: Machine header — ``nexus-meta`` schema (PROVISIONAL)
   :color: muted

   .. code-block:: yaml

      module: transport
      concept: wavefront_cochain
      role: "the interior face-flux cochain C1_int — the sweep-internal cochain, its biproduct C1 = C1_int (+) C1_boundary and trace algebra, and why the typed WavefrontFlux carrier retired (the concept survives in two native realizations)"
      depends_on: [operator_algebra]
      related: [operator_adjoint, boundary_conditions]
      status: "extracted from operator_algebra.rst; content verbatim, provisional header"


This page develops the **interior face-flux cochain**
:math:`C^1_{\rm int}` — the sweep-internal record of the angular flux on
interior cell faces, distinct from the *iterate state* of the field
algebra. It is a companion to the operator algebra developed in
:doc:`/theory/foundations/operator_algebra`: where that page types the
operators and :doc:`/theory/foundations/field_algebra` types the flux
**states** they act on, this page types the sweep-internal cochain that
carries flux across faces during a transport sweep — its biproduct
decomposition :math:`C^1 = C^1_{\rm int} \oplus C^1_\partial`, its trace
algebra (:math:`\iota_*` / :math:`\iota^*`), and why the typed
``WavefrontFlux`` carrier retired at S6.4(f) while the cochain
*mathematics* survives in its two native realizations.

.. _wavefront-flux-cochain:

The interior face-flux cochain — :math:`C^1_{\rm int}`
======================================================

.. note:: **Succession (S6.4(f), Issue #222, 2026-06-10) — the typed
   carrier retired; the concept survives in its two native
   realizations.**

   The interior face cochain :math:`C^1_{\rm int}` was, from #205 Phase 5
   through S6.4, carried by a dedicated named field ``WavefrontFlux``
   (on the space ``InteriorFaceSpace``). **That type is retired** — the
   modules ``orpheus/transport/fields/wavefront_flux.py`` and
   ``orpheus/numerics/spaces/interior_face_space.py``, and the 25
   foundation tests in ``tests/transport/fields/test_wavefront_flux.py``,
   are deleted. The cochain **mathematics** below — the biproduct
   :math:`C^1 = C^1_{\rm int} \oplus C^1_\partial`, the trace algebra
   :math:`\iota_*` / :math:`\iota^*`, the flux-only single-role
   rationale, the storage × role × locus grid — **remains valid theory**.
   Only the Python carrier is gone, and the historical derivation is
   kept below in past tense because it explains *why* the cochain frame
   is the right one.

   **Why it retired.** The S6.4 walk re-layering (see
   :ref:`sweep-dispatch-relayering` in :doc:`/theory/methods/sn/index`)
   dissolved the type's two load-bearing verbs — the whole-trace
   :math:`\iota_*` ``seed`` (read the octant inflow into the
   domain-edge slots) and the whole-trace :math:`\iota^*` ``absorb``
   (capture the domain-edge outflow) — into the shared ``_OctantWalk``
   frame, where the per-octant inflow read, the per-octant outflow
   shed, and the **single** :ref:`O.4b <bc-extraction>` boundary block
   live once for every walk. The type's whole-:math:`N` **mesh-bound**
   storage then had no place in the **per-octant** walk: an octant
   transient cannot be a persistent mesh-bound field.

   **Where the concept lives now** — two native realizations, each
   truer to a different facet of what the type held:

   * **The values AT the moving wavefront** —
     ``orpheus.sn.loss_representation.sweep_graph._MovingFrontier``, the rolling
     :math:`(d{-}1)`-frontier (per-level ``seed`` /
     ``shed``). Arguably the *truer* realization: the retired type
     held the wavefront's complete **history**; the frontier holds the
     **front itself**, ping-ponged by parity across two active levels.
   * **The full cochain history (the oracle's fuller view)** —
     ``FullFieldWavefront._octant_face_cochain`` raw per-axis buffers
     in ``orpheus.sn.loss_representation`` (the in-edge :math:`\iota_*`
     seed; ``_edge_outflow`` is the :math:`\iota^*` extraction),
     retained for verification cross-checks of the production
     ``_MovingFrontier`` window.

   The type was *a useful concept while it lived* — it named the
   :math:`\iota_*` / :math:`\iota^*` trace operators the raw-numpy
   sweep applied by hand, and it is the substrate on which the
   biproduct decomposition was first articulated. The succession keeps
   that articulation; only the now-redundant carrier is dropped.

Wave O step #205 Phase 5 (`Issue #208
<https://github.com/deOliveira-R/ORPHEUS/issues/208>`_, commits
``478723d`` mint / ``992b0c0`` 2-D sweep / ``0e3e16c`` 2-D matvec,
2026-06-04) typed the SN wavefront sweep's **interior** cell-face
angular fluxes — historically the raw ephemeral numpy arrays ``psi_x``
``(N, ng, nx{+}1, ny)`` / ``psi_y`` ``(N, ng, nx, ny{+}1)`` — as a
named field ``WavefrontFlux``.
This was the **face-locus** sibling of the boundary-block typing of
:ref:`bc-extraction` (cell + trace) and the operator-output typing of
:ref:`bc-extraction-operator-output-typing`: where those typed the
*cell* state and the *operator outputs*, this typed the *interior
face* state that the wavefront sweep propagated between them. It killed
the ``coding-elegance`` Pattern-3 anti-pattern (a flux-bearing tensor
with no type identity) and **named the trace operator the sweep
applied by hand** — the algebra that the S6.4 walk frame later
absorbed.


The native frame — discrete exterior calculus / cochains
--------------------------------------------------------

The native mathematical structure of the SN face fluxes (validated by
the cross-domain-attacker, frame memo
``field_role_typing_faceflux_frames.md``) is **discrete exterior
calculus**. The per-ordinate angular flux crossing faces is a primal
**1-cochain** — an assignment of a value to each oriented face:

.. (vv-status rationale) Structural framing of the SN face fluxes as a
   primal 1-cochain. This is a definitional / representational identity
   (the named-field typing — now carried by _MovingFrontier /
   _octant_face_cochain after the WavefrontFlux carrier retired at
   S6.4(f)), not a solver claim; the verifiable content is the
   bit-identity of the typed walk against the raw psi_x/psi_y walk,
   pinned by the octant-equivalence and Gate-K suites below (the 25
   foundation tests of the retired test_wavefront_flux.py went with the
   type; the window ≡ full oracle now pins the cochain walk).
.. vv-status: wavefront-cochain-primal documented

.. math::
   :label: wavefront-cochain-primal

   \psi^{(1)}_\Omega \;\in\; C^1
   \qquad (\text{a value per oriented face, per ordinate } \Omega).

The cell-average :math:`\bar\psi` is a **0-cochain**
:math:`\bar\psi^{(0)} \in C^0`, and the diamond-difference closure
:math:`\psi^{\rm out} = 2\bar\psi - \psi^{\rm in}` is the discrete
statement that :math:`\bar\psi` is the *midpoint* of the face-pair
bounding the cell — i.e. the **averaging map**
:math:`C^1 \to C^0` co-restricted to the cell (the lowest-order
Whitney 1-form → 0-form reduction, the trapezoidal/midpoint Hodge star
on a 1-D edge). The :ref:`balance-cartesian-2d` DD already realises
this; the cochain frame simply names the spaces.


The biproduct :math:`C^1 = C^1_{\rm int} \oplus C^1_\partial`
-------------------------------------------------------------

The full face cochain splits as a **biproduct** into the interior and
the boundary 1-chains:

.. (vv-status rationale) The face-cochain biproduct decomposition. A
   structural identity mirroring the cell+trace direct sum of
   :eq:`bc-extraction-direct-sum-state` one locus down; the verifiable
   content is the round-trip + projection laws below, pinned by
   test_absorption_is_identity / test_pi_int_after_injection_is_zero_2d.
.. vv-status: wavefront-cochain-biproduct documented

.. math::
   :label: wavefront-cochain-biproduct

   C^1 \;=\; C^1_{\rm int} \;\oplus\; C^1_\partial,
   \qquad
   \begin{cases}
     C^1_{\rm int} & (\text{interior faces}), \\
     C^1_\partial  & (\text{domain-edge faces, } \texttt{AngularBoundaryFlux}),
   \end{cases}

coupled by the injection / projection at the domain edges. This is the
**same biproduct** as the cell+trace direct-sum transport state
:eq:`bc-extraction-direct-sum-state`
(:math:`V = V_{\rm bulk} \oplus V_{\rm boundary}`), realised **one
locus down** at the face level. The two loci nest:

.. list-table:: The two biproduct loci (cell+trace ‖ face)
   :header-rows: 1
   :widths: 22 30 26 22

   * - Locus
     - Interior summand
     - Boundary summand
     - Coupling
   * - **cell + trace**
       (:eq:`bc-extraction-direct-sum-state`)
     - :math:`V_{\rm bulk}`
       (:class:`~orpheus.transport.fields.angular_flux.AngularFlux`)
     - :math:`V_{\rm inflow} \oplus V_{\rm outflow}`
       (:class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`)
     - the streaming sweep + sibling :math:`-B`
   * - **face**
       (:eq:`wavefront-cochain-biproduct`)
     - :math:`C^1_{\rm int}`
       (``_MovingFrontier`` front / ``_octant_face_cochain`` history)
     - :math:`C^1_\partial`
       (:class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`)
     - the trace operators :math:`\iota_*` / :math:`\iota^*`

:class:`AngularBoundaryFlux` is the boundary summand at **both** loci — it is
literally the domain-edge faces, which are exactly the
boundary trace of the cell+trace state. The interior summand differs:
the cell biproduct carries the cell-centre flux, the face biproduct
carries the interior *face* flux. The boundary persists (it carries
the SI / Krylov iterate across calls); the interior is **ephemeral**
(rebuilt each sweep — which is *why*, at S6.4(f), no mesh-bound type
survived for it), so the two summands have different lifetimes —
:eq:`wavefront-cochain-biproduct` is therefore a lifetime split as
well as a spatial split.


The trace operators :math:`\iota_*` (seed) / :math:`\iota^*` (absorb)
---------------------------------------------------------------------

The boundary trace is the **pullback** :math:`\iota^*` of the interior
cochain under the inclusion :math:`\iota \colon \partial\Omega
\hookrightarrow \Omega`. In the discrete setting :math:`\iota^*` is
"select the domain-edge faces"; its adjoint injection :math:`\iota_*`
is "write the boundary trace into the domain-edge faces". These are
exactly what the pre-typing sweep did by hand —
``psi_x[:, :, 0, :] = boundary.face_view("xmin")`` (the seed) and the
write-back (the absorb). The retired ``WavefrontFlux`` named them as
``seed`` / ``absorb`` methods; after S6.4(f) the same verbs live on the
two native realizations — ``_MovingFrontier.seed`` / ``.shed`` (the
per-level front) and ``FullFieldWavefront._octant_face_cochain`` (the
in-edge :math:`\iota_*` seed) / ``_edge_outflow`` (the :math:`\iota^*`
extraction):

.. list-table:: The trace operators (now on the realizations)
   :header-rows: 1
   :widths: 14 22 28 36

   * - Symbol
     - Realization verb
     - Direction
     - Role in the biproduct
   * - :math:`\iota_*`
     - ``_MovingFrontier.seed`` /
       ``FullFieldWavefront._octant_face_cochain`` in-edge
     - :math:`C^1_\partial.\text{inflow} \to C^1_{\rm int}`
       domain-edge slots
     - the **injection** :math:`\iota_\partial` of the biproduct
   * - :math:`\iota^*`
     - ``_MovingFrontier.shed`` /
       ``FullFieldWavefront._edge_outflow``
     - :math:`C^1_{\rm int}` domain-edge slots
       :math:`\to C^1_\partial.\text{outflow}`
     - the **projection** :math:`\pi_\partial` of the biproduct

Historically the retired type routed both through a single
single-source-of-truth ``_edge_slot`` face-to-edge map, so the
injection and projection could not desync, and exposed a third read
accessor ``edge_view`` that returned the domain-edge slot as a
zero-copy view *without* copying into a :class:`AngularBoundaryFlux` — the
2-D matvec used it to difference ``edge_view(face) − given`` when
emitting its active-trace boundary residual (so there was no hardcoded
``psi_x[:, :, 0, :]`` literal at the call site). After S6.4(f) the same
single-map discipline is the shared ``_OctantWalk`` frame's: the
per-octant inflow read and outflow shed are the ONE boundary block
every walk uses, so the injection/projection cannot desync there
either.

The two biproduct laws follow, and are **provable** rather than
coincidental:

.. (vv-status rationale) The two biproduct identities — absorption =
   identity (project-after-inject) and projection-annihilates-the-
   strictly-interior (the off-diagonal block is zero). Structural laws
   of the biproduct; were pinned by test_absorption_is_identity (slab /
   sphere / 2-D box) and test_pi_int_after_injection_is_zero_2d in the
   retired test_wavefront_flux.py. After S6.4(f) the same two laws are
   the content of the window ≡ full-field oracle (the seed/shed
   round-trip equals the full-cochain seed/extract, bit-identically).
.. vv-status: wavefront-cochain-biproduct-laws documented

.. math::
   :label: wavefront-cochain-biproduct-laws

   \iota^* \circ \iota_* \;=\; \mathrm{id}
   \quad (\text{on the boundary chain}),
   \qquad
   \pi_{\rm int} \circ \iota_\partial \;=\; 0.

The first — :math:`\iota^* \circ \iota_* = \mathrm{id}` — IS the
"absorption = identity" fact (``seed`` then ``absorb``, with no
wavefront walk between, returns the boundary trace unchanged). It was
an *observation* under the raw-numpy seed/absorb, then a **named
biproduct law** under the typed ``WavefrontFlux``; after S6.4(f) it is
the seed/shed round-trip the ``_OctantWalk`` frame performs once per
octant. The second — :math:`\pi_{\rm int} \circ
\iota_\partial = 0` — is the biproduct's off-diagonal-zero condition:
injecting the boundary leaves the **strictly**-interior faces
(positions :math:`1 \ldots n{-}1` along each axis) untouched at zero
(in the full-cochain realization, ``_octant_face_cochain`` zero-inits
every interior slot and seeds only the in-edge, so the strictly-interior
faces are untouched by construction).


Why the interior cochain is flux-only (single role)
---------------------------------------------------

The interior cochain carries the **flux** role only — unlike the
boundary trace, it has no source / residual leaves (the retired
``WavefrontFlux`` was, accordingly, flux-only; the surviving
``_MovingFrontier`` / ``_octant_face_cochain`` buffers are likewise
plain flux arrays). The reason is
structural, and it explains the role grid of
:ref:`bc-extraction-operator-output-typing` from the cochain side. Per
the second native frame (sparse triangular factorisation), the interior
face flux is the **off-diagonal coupling of the per-octant
lower-triangular streaming factor** :math:`L_{\rm oct}`: the entry
coupling a cell to its upwind neighbour is the flux on the shared
interior face. That factor is **re-formed each sweep** (forward
substitution down the upwind DAG), so the interior face flux is a
*transient of the factorisation*, not a persistent state with a defect
of its own. The role grid (flux / source / residual) is a **0-cochain
(cell) concept**: a residual arises from a balance of cell-level
reaction rates against a cell-level source. Only the **boundary
1-chain** inherits the role grid — and only because the boundary
consistency residual

.. math::

   r_\Gamma \;=\; \psi.\text{inflow}
                  \;-\; B\,\psi.\text{outflow}
                  \;-\; q.\text{inflow}

(:eq:`bc-extraction-two-residuals`) lives on the boundary faces. The
strictly-interior faces carry no such balance, so the interior cochain
is flux-only by construction — ``illegal-states-unrepresentable``: there
was no ``InteriorFaceResidual`` to mistype an interior face as (and
under the surviving plain-array realizations there is no role grid to
violate at all).


The storage × role × locus grid (Issue #205), extended with the face locus
---------------------------------------------------------------------------

`Issue #205
<https://github.com/deOliveira-R/ORPHEUS/issues/205>`_ frames the
cross-method field architecture as a **storage × role × locus**
classification: a typed transport field is pinned by *where* it lives
(the **locus** — cell, face, boundary), *what role* it plays (flux /
source / residual), and *how it is stored* (persistent vs ephemeral).
Phase 5 adds the **face** locus row:

.. list-table:: Storage × role × locus — the typed field grid (#205)
   :header-rows: 1
   :widths: 16 18 22 22 22

   * - Locus
     - Storage
     - Flux role
     - Source/sink role
     - Residual role
   * - **cell**
       (0-cochain :math:`C^0`)
     - persistent (the iterate)
     - :class:`~orpheus.transport.fields.angular_flux.AngularFlux`
       / :class:`~orpheus.transport.fields.scalar_flux.ScalarFlux`
     - :class:`~orpheus.transport.source_sinks.angular_source_sink.AngularSourceSink`
     - :class:`~orpheus.transport.residuals.angular_residual.AngularResidual`
       (minted, O.2 consumer)
   * - **face — interior**
       (1-cochain :math:`C^1_{\rm int}`)
     - **ephemeral** (rebuilt each sweep)
     - ``_MovingFrontier`` front /
       ``_octant_face_cochain`` history
       (**#205 Phase 5** ``WavefrontFlux``, retired S6.4(f))
     - — (off-diagonal of :math:`L_{\rm oct}`, no source role)
     - — (no cell-balance defect on interior faces)
   * - **face — boundary**
       (1-cochain :math:`C^1_\partial`)
     - persistent (carries the iterate trace)
     - :class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`
     - :class:`~orpheus.transport.source_sinks.angular_boundary_source_sink.AngularBoundarySourceSink`
       (B.5.2)
     - :class:`~orpheus.transport.residuals.angular_boundary_residual.AngularBoundaryResidual`
       (minted, O.2 consumer)

The grid makes the flux-only-single-role rationale (above) visible: the
interior-face row has only the flux cell populated, because the
off-diagonal-of-:math:`L_{\rm oct}` reading forbids a source role and
the no-cell-balance reading forbids a residual role. The boundary-face
row is fully populated only because the boundary 1-chain hosts the BC
consistency residual. The face locus mirrors the cell locus on storage
(interior-face ephemeral ‖ boundary-face persistent is the face-level
analogue of the bulk iterate persisting while operator outputs are
transient).


Field + views, NOT per-face objects
------------------------------------

The retired ``WavefrontFlux`` stored a **single flat backing buffer**
(``space.layout.total_size``); the per-axis face fields were zero-copy
reshape views (a ``face(axis)`` accessor).
The cross-domain-attacker **rejected** a per-face Python object on three
independent grounds — a rejection that **still binds** the surviving
realizations (``_MovingFrontier`` and ``_octant_face_cochain`` are both
plain dense arrays indexed by the batch, never per-face objects):

1. **Vectorisation (load-bearing).** The unit of operation is the
   ``(N_oct, ng, n_diag)`` wavefront batch; a per-face object would
   reintroduce the per-cell / per-face Python fold that caused a
   10–20× slab regression in earlier work. The hot-path walk indexes
   the ``face(axis)`` views with byte-identical fancy-indexing to the
   legacy raw ``psi_x`` / ``psi_y``.
2. **The cochain frame is storage-granularity-indifferent.** A cochain
   is an assignment of values to faces; whether stored as one flat
   buffer or per-face objects is an implementation choice the frame
   does not constrain — so the vectorisation argument settles it for
   the flat buffer.
3. **Biproduct consistency.** The
   :eq:`wavefront-cochain-biproduct-laws` projections are clean
   slice-views only if both summands live on the same flat-buffer
   substrate.

This mirrored :class:`AngularBoundaryFlux`, which still uses
``flat-buffer + FaceLayout + face_view``. The interior space was the
retired ``InteriorFaceSpace`` — the **layout-on-space** (A.5) sibling
of :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace` minus the
``omega_dot_n`` directional selectors (the interior cochain has no
inflow/outflow partition — it is flux-only). It was axis-parametric: its
``interior_layout`` builder took the axis count as a parameter (one
face-normal field per active axis). The dimension-genericity that this
parametricity foreshadowed is now carried *directly* by the
realizations: ``_MovingFrontier`` is built from a mesh-time
``_FrontierPlan`` whose ``(d{-}1)``-frontier slab is constructed for any
``d`` (a point at d=1, a line at d=2, a surface at d=3), and
``FullFieldWavefront._octant_face_cochain`` allocates one ``n_a + 1``
buffer per axis for any ``ndim`` — so a future 3-D Cartesian wavefront
sweep is a *parameter*, not a new code path
(``feedback_unify_after_two_instances``: 1-D + 2-D are two working
instances, 3-D is the validating third).


.. _wavefront-flux-honest-scope:

Honest scope — a representation win, NOT a speed/rate/parallelism win
---------------------------------------------------------------------

.. warning::

   The interior-cochain **typing** was a **representation / elegance win
   only**. It did **NOT**:

   * change the asymptotic cost of the sweep,
   * recover the source-iteration convergence rate,
   * enable parallelism on its own.

   The cross-domain-attacker stated this as an *absence*, not a hedge:
   the wavefront DAG already gives the optimal sweep schedule, and the
   typing relocated no arithmetic. The seed/absorb stays an **inherent
   cheap** :math:`O(\text{boundary faces})` copy — negligible against
   the :math:`O(\text{volume})` sweep — at the **persistent-boundary /
   ephemeral-interior lifetime split**. True zero-copy between the two
   summands is *precluded* (and unnecessary) because
   :class:`AngularBoundaryFlux` persists across SI iterations while the
   interior cochain is rebuilt each sweep — the very lifetime split that
   later left it with **no** mesh-bound type after S6.4(f).

   The asymptotic-cost / peak-memory wins are sought elsewhere: the
   **persistent** SI iterate is shrunk by the orthogonal
   :ref:`angular windowing <sn-angular-windowing>` of Phase 5a (the
   iterate lives in moment space, :math:`N \to (L{+}1)(2L{+}1)`), and the
   **per-sweep transient** interior cochain
   itself — the dominant peak-memory cost noted here — was the target of
   Phase 5b's storage-B rolling moving-frontier window, which never
   materializes the whole interior cochain at once and is the 3-D
   enabler. That window IS the ``_MovingFrontier`` that the interior
   cochain now lives in: the production realization *is* the
   memory-win, with the full-cochain ``_octant_face_cochain`` kept only
   as the verification oracle.

The payoff *while the type lived* was the **type**: a named field,
typed :math:`\iota_*` / :math:`\iota^*`, code that reads like the
cochain math, and illegal-states-unrepresentable (the flux-only
constraint of the interior locus enforced by there being no
interior-face residual leaf). It was also framed as the **clean
substrate** the SI Gauss–Seidel recovery would land on: the
``(octant × face)`` reflective-graph G-S schedule as an explicit
composition (``sweep octant → ι* absorb → −B reflect → ι_* seed next
octant``) rather than an implicit buffer-timing dance. The S6.4 walk
re-layering delivered exactly that explicit composition — but as the
shared ``_OctantWalk`` frame's per-octant ``read inflow → shed outflow
→ −B boundary block`` sequence, **without** a standalone interior-face
type. The composition survived; the carrier did not. That recovery
remains where the actual convergence-rate win is sought (a separate,
research-tagged step).


Numerical evidence — type-only ⟹ bit-identical
-----------------------------------------------

Like the B.5.2 retype, Phase 5 wraps the **same** buffers: the
``face(axis)`` views share memory with the flat backing, and the
wavefront walk indexes them with byte-identical fancy-indexing, so
``.values`` are unchanged — only the type and the named
:math:`\iota_*` / :math:`\iota^*` differ. The change is verified by
**bit-identity (``np.array_equal``) by inheritance** from the
already-verified raw-numpy path:

.. list-table:: Phase 5 bit-identity gates
   :header-rows: 1
   :widths: 30 38 32

   * - Gate
     - Pins
     - Evidence
   * - **Phase 0 de-risk**
       (``diag_phase0_wavefront_derisk.py``, diagnostic)
     - typed :math:`\iota_*` / :math:`\iota^*` round-trip + full bare
       sweep on a flat-buffer backing
     - :math:`\max|\Delta| = 0.0` on angular / scalar / boundary,
       all reflective + vacuum cases; ``shares_memory`` zero-copy;
       transposed-seed negative control breaks (non-vacuous)
   * - **2-D octant equivalence**
       (``test_2d_octant_sweep_equivalence.py``)
     - the typed 2-D sweep reproduces the legacy octant snapshots
     - bit-identical octant snapshots
   * - **Gate-K** :math:`k_\infty = 1.875`
       (``test_2d_l2_matvec_correctness``)
     - the 2-D matvec recovers the closed-form homogeneous eigenvalue
       (:math:`k_\infty = \nu\Sigma_f / \Sigma_a`, closed-form pillar)
     - bit-identical to the pre-typing value
   * - **matvec ≡ sweep** + **BC extraction**
       (``test_bc_extraction_2d`` / ``_matvec``)
     - the typed matvec and sweep are ONE discretisation; the bare-2-D
       BC block matches
     - 126 passed
   * - **foundation suite**
       (``tests/transport/fields/test_wavefront_flux.py`` — RETIRED
       at S6.4(f) with the type)
     - units / class identity / field+views / the two biproduct laws /
       the round-trip pin + L11 negative control / axis-parametricity
       (1-D / 2-D / 3-D one path)
     - was all ``@pytest.mark.foundation`` (software invariants, no
       theory ``:label:``); the biproduct-law content these pinned now
       lives in the ``window ≡ full-field`` oracle
   * - **L16 perf**
       (``diag_l16_wavefront_microbench.py``, diagnostic)
     - no per-cell / per-face Python crept into the hot path
     - typed / raw median wall-clock ratio :math:`\approx 1.00`
       (within the +5 % gate)

The closed-form Gate-K eigenvalue is the only *correctness* (not
merely *equivalence*) anchor here: it is verified against
:math:`k_\infty = \nu\Sigma_f / \Sigma_a` (the closed-form pillar — MMS
does **not** prove eigenvalues). The remaining gates are bit-identity
by inheritance (``vv-principles`` § "Bit-identity vs
principled-equivalence", criterion: the implementation is unchanged —
the wrapper wraps the same buffers, so the values are byte-identical,
not merely close). The **A2D-1 source-hash** deliberate-edit tripwire
was refreshed (``f683f229…`` → ``12697ab3…``, behaviour-neutral — the
hash pins that any future edit to the 2-D matvec source is intentional).


The 1-D sweep is a scan, not a wavefront — deferred to ``nd_foundation``
------------------------------------------------------------------------

.. note::

   The retired ``WavefrontFlux`` typed the **2-D** wavefront sweep +
   matvec only. The **1-D** sweep is a parallel-prefix **scan** (a
   different fold): its interior fluxes are transient chain-ordered scan
   output (``(nx, K, ng)``, no persistent interior buffer at all), the
   structural antithesis of the cochain ``(N, ng, nx{+}1)``. Forcing the
   type into the 1-D scan would have been a wrong-fit (it would risk the
   L16 ``cumprod`` efficiency and multiply concepts), so the 1-D
   *implementation* kept its scan fold. The one shared seam (the
   boundary-trace exchange + DD-closure averaging) was deferred to a
   future ``nd_foundation`` session that would re-express **both** folds
   as one :math:`d`-generic walk, under the hard constraint that the
   collapse must not regress the :math:`d=1` scan's parallel-prefix
   efficiency. **That collapse has since landed** (the ``CumprodScan``
   d=1 sibling and the ``_DAGWavefront`` family share the per-octant
   DAG; see :ref:`sweep-dispatch-relayering` in
   :doc:`/theory/methods/sn/index`), which is part of *why* the standalone
   ``WavefrontFlux`` carrier became redundant at S6.4(f): the
   :math:`d`-generic walk supplies the front (``_MovingFrontier``) and
   the oracle history (``_octant_face_cochain``) for every :math:`d`. A
   future 3-D *Cartesian* sweep IS a wavefront and rides the same
   ``_MovingFrontier`` window directly (the :math:`(d{-}1)`-frontier
   slab is built for any :math:`d`).
