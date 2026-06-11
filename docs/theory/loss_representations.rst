.. _loss-representations:

==================================================================
Selectable Representations of the S\ :sub:`N` Loss Operator
==================================================================

This is the **capstone architecture page** for the within-group
S\ :sub:`N` loss operator :math:`(L+C)` — the streaming-plus-collision
operator whose inverse *is* the transport sweep. It documents the
final state of issue #222 (the *sweep-strategy carve* + the S6
operator/representation re-layering, complete 2026-06-11 at the
Fork-B2 default flip): one lower-triangular operator, several
*algorithms* that realise it, a single source of truth for which
algorithm a mesh gets, and the L21 theorem that makes "matvec is the
same operator as the sweep" a type fact rather than a coincidence.

It deliberately does **not** re-derive the cell-level mathematics.
Those live in :doc:`discrete_ordinates`:

* the WDD cell balance and closure — :eq:`dd-cartesian-1d`,
  :eq:`dd-cartesian-2d`, :eq:`dd-2d-balance-form`;
* the 1-D cumprod recurrence — :ref:`sweep-cumprod`;
* the 2-D anti-diagonal wavefront — :ref:`sweep-wavefront`;
* the three-layer **walk / level-op / kernel** stack — the
  :ref:`sweep-dispatch-relayering` section (S6.4(e)).

…and the operator algebra lives in :doc:`operator_algebra`:

* the typed forward/solve/adjoint actions —
  :eq:`operator-apply`, :eq:`operator-solve`,
  :eq:`operator-apply-transpose`;
* the interior face-flux cochain :math:`C^1_{\rm int}` and its
  succession after the typed ``WavefrontFlux`` retirement —
  :ref:`wavefront-flux-cochain`.

This page is the layer *above* both: the **representation** layer that
selects and unifies them.

.. contents:: On this page
   :local:
   :depth: 2


Key Facts
=========

.. admonition:: Key Facts — the loss-representation architecture
   :class: important

   * **One operator, several schedules.** :math:`(L+C)` is
     **lower-triangular** under the upwind (sweep) cell ordering. Its
     two actions are one object: :math:`SOLVE = (L+C)^{-1}q` is forward
     substitution (the *transport sweep*); :math:`APPLY = (L+C)\psi`
     is the row action (the *Krylov matvec*). A
     :class:`~orpheus.sn.loss_representation.LossRepresentation` is a
     **schedule** for traversing that triangular structure, not a
     different operator.

   * **The four representations**
     (:mod:`orpheus.sn.loss_representation`), each a stateless frozen
     ``@dataclass`` carrying only the mesh:

     - :class:`~orpheus.sn.loss_representation.CumprodScan` — the 1-D
       Blelloch parallel-prefix scan (slab + sphere + cylinder via one
       body); the **1-D production default**.
     - :class:`~orpheus.sn.loss_representation.ScanMarch` —
       :math:`\mathrm{scan}(x)\circ\mathrm{march}(y)`; the **multi-D
       Cartesian production default since the S6.9 Fork-B2 flip**.
     - :class:`~orpheus.sn.loss_representation.MovingFrontierWindow` —
       the anti-diagonal wavefront carrying a rolling
       :math:`(d{-}1)`-frontier; a **selectable peer** (the d=2 default
       through S6.9).
     - :class:`~orpheus.sn.loss_representation.FullFieldWavefront` —
       the same DAG schedule retaining the **whole** interior cochain;
       the **verification oracle**, explicit-select only.

   * **Selection is a single source of truth.**
     :meth:`~orpheus.sn.loss_representation.LossRepresentation.supports`
     returns :class:`~orpheus.sn.loss_representation.Compatibility`
     ``(ok, reason)``; :func:`~orpheus.sn.loss_representation.default_for`
     returns the first match in the ordered registry
     :data:`~orpheus.sn.loss_representation.LOSS_REPRESENTATIONS`.
     Illegal ``(representation, mesh)`` pairings are unrepresentable —
     ``__post_init__`` re-checks ``supports`` and raises
     :class:`~orpheus.sn.loss_representation.IncompatibleRepresentation`.

   * **L21 — one walk, one instance.** A single d-generic
     ``_OctantWalk._interior_walk`` frame serves sweep AND matvec for
     every multi-D representation, forked only by a **kernel object** +
     **emit policy** (never a boolean). The operator holds ONE
     representation instance consumed by ``apply``, ``solve``, and the
     Gauss–Seidel resolvent. Both are *type facts*, pinned by spy tests.

   * **The Fork-B2 decision (the WHY of the default).** ScanMarch is
     **0.55–0.84×** the window's per-sweep time at **identical** peak
     memory; the window's memory claim only ever held against the
     full-field oracle (1.3–1.4× both). The window was **kept**, not
     retired — a genuinely different schedule over the same operator is
     the point of selectability.

   * **Bit-identity vs principled-equivalence.** Different *schedules*
     are NOT bit-comparable (FP-association differs by construction) →
     ``nulp`` / solver-tol Mode-9 gates. Different *storage policies* of
     the same schedule ARE bit-identical → ``array_equal`` oracles. The
     cell kernel's explicit left fold ``((Σ_t + s_0) + s_1)`` is the FP
     reduction tree of record (sha256-source-pinned).

   * **Gotcha — the operator subtracts** :math:`C` **once.**
     ``loss_action`` MUST return the **full** :math:`(L+C)\psi`, NOT
     :math:`L\psi`. The operator applies the only glue
     :math:`L = (L+C) - C` (Resolution A). A leaf returning
     :math:`L\psi` would double-count the collision diagonal.


.. _loss-rep-native-frame:

The native mathematical frame: a lower-triangular operator
==========================================================

The within-group discrete-ordinates balance, for one ordinate
:math:`\Omega_n` and one cell, is the WDD relation derived in
:doc:`discrete_ordinates` (:eq:`dd-cartesian-2d`). Collect every cell
and every ordinate into the within-group operator

.. math::
   :label: loss-rep-LpC

   (L+C)\,\psi \;=\; q,
   \qquad
   L \;=\; \Omega\cdot\nabla\big|_{\rm WDD},
   \quad
   C \;=\; \sigma_t\,\odot,

where :math:`L` is the discretised streaming operator (the upwind
WDD difference relations) and :math:`C` is the collision diagonal
(multiply-by-:math:`\sigma_t`). The decisive structural fact:

.. admonition:: :math:`(L+C)` is lower-triangular under the sweep ordering
   :class: note

   Order the unknowns cell-by-cell in **upwind (inflow-to-outflow)
   order** for a given ordinate. Then cell :math:`i`'s balance depends
   only on faces that are *already known* — the domain inflow, or the
   outflow of strictly-upstream cells. In that ordering the matrix of
   :math:`(L+C)` is **lower-triangular** (block-lower-triangular over
   groups). A lower-triangular system has two canonical operations,
   and they are the same operator viewed two ways:

   .. list-table::
      :header-rows: 1
      :widths: 14 30 30 26

      * - Operation
        - Linear-algebra name
        - S\ :sub:`N` name
        - Code entry
      * - :math:`SOLVE`
        - forward substitution :math:`(L+C)^{-1}q`
        - the **transport sweep**
        - ``LossRepresentation.sweep``
      * - :math:`APPLY`
        - the row action :math:`(L+C)\psi`
        - the **Krylov matvec**
        - ``LossRepresentation.loss_action``

This is the **L21 invariant** referenced throughout the SN code: *the
sweep and the matvec are two actions of the SAME operator* — they are
not independent code paths that happen to agree, and they must never be
allowed to drift. Forward substitution and the row action share the
triangular factor; they differ only in which is the unknown.

Resolution A — the operator's only glue
---------------------------------------

The representation returns the **full** within-group loss action
:math:`(L+C)\psi`. The operator
(:meth:`~orpheus.sn.operator.StreamingOperator.apply`,
:eq:`operator-apply`) then applies the *only* remaining algebra glue,

.. math::
   :label: loss-rep-resolution-a

   L\,\psi \;=\; (L+C)\,\psi \;-\; \sigma_t \odot \psi.\mathrm{bulk}
   \qquad\text{(Resolution A: } L = (L+C) - C\text{)},

subtracting the collision diagonal :math:`C = \sigma_t\,\odot`
**exactly once**. Before the S6.3 re-layering this :math:`-C`
subtraction was duplicated five times across five private ``_apply_*``
bodies; collapsing it to one site is a single-source-of-truth win, but
it imposes the **load-bearing contract** on every representation:

.. warning::

   ``loss_action`` MUST return :math:`(L+C)\psi`, **not** :math:`L\psi`.
   A leaf that returned the bare streaming action :math:`L\psi` would
   make the operator subtract :math:`C` a *second* time — a
   double-counted collision diagonal, a silent sign-and-magnitude error
   (vv-principles failure Mode 3, *missing/duplicated factor*). The
   convention is pinned by ``test_loss_action_convention.py``: the
   non-tautological anchor checks that for a flat reflective field
   :math:`L\psi_{\rm flat} = 0`, so :math:`(L+C)\psi = \sigma_t\psi` —
   proving the action is the *full* :math:`(L+C)` loss, not bare
   :math:`L`, and cross-checks the :math:`-C` glue against an
   independent :class:`~orpheus.sn.operator.CollisionOperator`.

The two actions bottom out, for every Cartesian representation, in **one
pure kernel pair** on
:class:`~orpheus.sn.spatial.diamond.DiamondDifference`:

* :meth:`~orpheus.sn.spatial.diamond.DiamondDifference.cell_kernel_batch`
  (solve) — :math:`\bar\psi = \mathrm{numer}/\mathrm{denom}`;
* :meth:`~orpheus.sn.spatial.diamond.DiamondDifference.residual_kernel_batch`
  (apply) — :math:`r = \mathrm{denom}\cdot\bar\psi - \mathrm{numer}`,

with the explicit **left fold**
:math:`((\Sigma_t + s_0) + s_1) + \cdots` as the IEEE-754 reduction
tree of record (see :ref:`loss-rep-bit-vs-principled`).


.. _loss-rep-four:

The four representations — four schedules of one operator
=========================================================

A :class:`~orpheus.sn.loss_representation.LossRepresentation` is a
stateless frozen dataclass (its only field is the
:class:`~orpheus.sn.geometry.SNMesh` it was selected for). Each is a
distinct **schedule** over the same lower-triangular
:math:`(L+C)` — a different topological linearisation of the identical
cell dependencies. They are *algorithms*, not operators.

.. list-table:: The four loss representations
   :header-rows: 1
   :widths: 22 16 22 20 20

   * - Representation
     - ``supports``
     - Schedule
     - Storage
     - Role
   * - :class:`~orpheus.sn.loss_representation.CumprodScan`
     - 1-D, any geometry
     - Blelloch prefix scan
     - chain recurrence
     - 1-D **production default**
   * - :class:`~orpheus.sn.loss_representation.ScanMarch`
     - 1-D OR Cartesian
     - ``scan(x) ∘ march(y)``
     - :math:`(d{-}1)`-slab per line
     - multi-D **production default**
   * - :class:`~orpheus.sn.loss_representation.MovingFrontierWindow`
     - Cartesian, d = 2
     - anti-diagonal wavefront
     - rolling :math:`(d{-}1)`-frontier
     - selectable **peer**
   * - :class:`~orpheus.sn.loss_representation.FullFieldWavefront`
     - Cartesian, any d
     - anti-diagonal wavefront
     - full interior cochain
     - verification **oracle**


CumprodScan — the 1-D chain
---------------------------

A 1-D mesh is a total order, so the sweep is a **chain prefix scan**.
:class:`~orpheus.sn.loss_representation.CumprodScan` is intrinsically
1-D — a prefix scan needs a total order, and there is no "2-D chain
scan" (this is *legitimate* d-specificity by the algorithm's nature,
not a narrow crutch). The geometry difference (slab vs sphere vs
cylinder) is absorbed by the two-stratum sweep cache, so all three
geometries share **one body**
(:func:`~orpheus.sn.loss_representation._sweep_1d_unified` →
:func:`~orpheus.sn.spatial.scan.ordinate_scan`); the curvilinear
Morel–Montry angular redistribution folds into the scan's affine
source. The recurrence and its closed-form cumprod/cumsum solution are
derived at :ref:`sweep-cumprod`; the pair-monoid algebra that justifies
the closed form is documented on
:func:`~orpheus.sn.spatial.scan.ordinate_scan`.

CumprodScan is conditioning-robust by construction: the closed-form
backend handles the pole reset (ERR-054) and denormal underflow
(ERR-057) that a naive cumprod would lose to gradual underflow.


ScanMarch — the row-march that reuses the scan
----------------------------------------------

:class:`~orpheus.sn.loss_representation.ScanMarch` reframes the d-D DD
sweep as forward substitution along the sweep axis — *the same
first-order linear scan* — **marched over the transverse axes**:

.. math::
   :label: loss-rep-scanmarch

   \mathrm{ScanMarch}
   \;=\;
   \begin{cases}
     \mathrm{scan}(x) & d = 1 \quad(\text{the } s_y = 0 \text{ degeneration})\\[2pt]
     \mathrm{scan}(x)\circ\mathrm{march}(y) & d = 2 \\[2pt]
     \mathrm{scan}(x)\circ\mathrm{march}(y,z) & d = 3.
   \end{cases}

Within each transverse row the diamond-difference x-face recurrence is
the **same Blelloch scan** that
:class:`~orpheus.sn.loss_representation.CumprodScan` uses
(:func:`~orpheus.sn.spatial.scan._scanmarch_row`); the transverse
coupling rides the affine source. The solve coefficients are

.. math::
   :label: loss-rep-scanmarch-solve

   \alpha \;=\; \frac{2 s_x}{D} - 1,
   \qquad
   \beta \;=\; \frac{2\,(Q + s_y\,\psi_{y,\rm in})}{D},
   \qquad
   D \;=\; \Sigma_t + s_x + s_y,

where :math:`s_a = 2|\mu_a|/\Delta a` is the per-axis streaming
coefficient. This **unifies** the 1-D
:class:`~orpheus.sn.loss_representation.CumprodScan` (its degenerate
:math:`s_y = 0` case) and the 2-D row-march in one primitive, and is
DAG-free (no graph is built).

The matvec twin reconstructs the interior x-faces from the *known*
probe :math:`\bar\psi`: because :math:`\bar\psi` is known, the WDD
closure :math:`\psi^{\rm out}_x = 2\bar\psi - \psi^{\rm in}_x` is itself
a first-order recurrence — a **pure-reflection scan** with
:math:`\alpha = -1`, :math:`\beta = 2\bar\psi`
(:func:`~orpheus.sn.spatial.scan._x_scan_faces`) — and the per-cell
residual is

.. math::
   :label: loss-rep-scanmarch-apply

   r_{i,j} \;=\;
   (\Sigma_t + s_x + s_y)\,\bar\psi_{i,j}
   \;-\; s_x\,\psi^{\rm in}_{x,i,j}
   \;-\; s_y\,\psi^{\rm in}_{y,i,j}
   \;\equiv\; (L+C)\bar\psi
   \quad(\text{at zero source}),

from which :eq:`loss-rep-resolution-a` subtracts :math:`\Sigma_t\bar\psi`
to give :math:`L\bar\psi`. ScanMarch additionally inherits the
conditioning robustness of ``ordinate_scan`` (ERR-054 / ERR-057 handled
per line for free) and is the natural home for the flux-independent
``a_attenuation`` two-stratum cache the wavefront lacks (the #206
follow-on).


MovingFrontierWindow — the rolling-frontier wavefront
-----------------------------------------------------

:class:`~orpheus.sn.loss_representation.MovingFrontierWindow` is the
anti-diagonal (level-scheduled) wavefront over the per-octant
:class:`~orpheus.sn.sweep_graph.SweepDependencyGraph` derived at
:ref:`sweep-wavefront`. It carries only a **rolling**
:math:`(d{-}1)`-frontier of interior face fluxes (a 2-diagonal at d=2),
advanced anti-hyperplane by anti-hyperplane via
:meth:`~orpheus.sn.sweep_graph.SweepDependencyGraph.walk_windowed`. The
frontier is the moving realisation of the interior face cochain
:math:`C^1_{\rm int}` — its theory and the post-``WavefrontFlux``
succession are at :ref:`wavefront-flux-cochain`. Its historical claim
to fame was a **~30 % peak-memory win over the full-field oracle**; the
Fork-B2 measurement (:ref:`loss-rep-fork-b2`) showed that the *same*
memory profile is shared by the row-march, which is why the window is
now a selectable peer rather than the default.


FullFieldWavefront — the verification oracle
--------------------------------------------

:class:`~orpheus.sn.loss_representation.FullFieldWavefront` walks the
**same** anti-diagonal DAG schedule as the window, but retains the
**whole** interior face cochain — the *fuller view*. It seeds only the
octant's domain in-edge slot (``_octant_face_cochain``); by the upwind
invariant every other slot is written before any read, so the
zero-initialised buffer is byte-identical to the historical whole-trace
:math:`\iota_*` seed. It is slower and more memory-hungry by design —
its purpose is to be the reference the d-specific production
optimisations are cross-checked against:

* the 1-D :class:`~orpheus.sn.loss_representation.CumprodScan`
  (principled-equivalence at nulp);
* the 2-D :class:`~orpheus.sn.loss_representation.MovingFrontierWindow`
  (``window ≡ full`` bit-identity);
* the multi-D :class:`~orpheus.sn.loss_representation.ScanMarch`
  (principled-equivalence at nulp).

Keeping a *fuller view of the concept* as a pinned verification oracle
is the deliberate **aggressive-retirement exception**: the production
paths could not be cross-checked structurally without it. It is the one
genuinely d-generic body (``supports`` is any-d Cartesian), so it is
also the admission spine for synthetic d=3 correctness tests before any
3-D quadrature exists.


.. _loss-rep-selection:

Selection: one predicate, three consumers
==========================================

Applicability is a **declared, queryable capability** — "make illegal
states unrepresentable" applied to method selection. Each representation
answers one classmethod:

.. code-block:: python

   class CumprodScan(_LossRepresentation):
       @classmethod
       def supports(cls, mesh):
           return Compatibility(mesh.is_1d, "requires a 1-D mesh")

   class ScanMarch(_LossRepresentation):
       @classmethod
       def supports(cls, mesh):
           return Compatibility(
               mesh.is_1d or mesh.is_cartesian,
               "requires a 1-D mesh (any geometry) or Cartesian geometry",
           )

   class _DAGWavefront(_LossRepresentation):       # MovingFrontierWindow's base
       @classmethod
       def supports(cls, mesh):
           return Compatibility(
               mesh.is_cartesian and mesh.ndim == 2,
               "requires Cartesian geometry, d = 2",
           )

   class FullFieldWavefront(_DAGWavefront):
       @classmethod
       def supports(cls, mesh):
           return Compatibility(mesh.is_cartesian, "requires Cartesian geometry")

The compatibility signal is the *genuine* criterion — the coordinate
system (:attr:`~orpheus.sn.geometry.SNMesh.is_cartesian`, i.e.
``curvature is None``) and the dimensionality
(:attr:`~orpheus.sn.geometry.SNMesh.ndim`) — **not** the
``sweep_graphs is None`` substrate proxy that the pre-carve code keyed
on. :class:`~orpheus.sn.loss_representation.Compatibility` is an
``(ok, reason)`` pair; the ``reason`` lets a teaching frontend gray-out
a method *and explain why* ("Moving-frontier window — requires Cartesian
geometry, d = 2"), which is pedagogically load-bearing — ORPHEUS teaches
reactor physics.

**One predicate, three consumers (single source of truth):**

#. **Frontend** —
   ``[R for R in LOSS_REPRESENTATIONS if R.supports(mesh).ok]`` lists
   the applicable methods. A cylinder (non-Cartesian) → only
   ``CumprodScan`` and ``ScanMarch``; the dropdown shows exactly those.

#. **Factory default** —
   :func:`~orpheus.sn.loss_representation.default_for` returns the
   **first** entry in the ordered registry
   :data:`~orpheus.sn.loss_representation.LOSS_REPRESENTATIONS` whose
   ``supports`` admits the mesh, falling back to the oracle so it is
   never stuck:

   .. list-table:: ``default_for`` outcomes (first-supports-match)
      :header-rows: 1
      :widths: 22 40 22

      * - mesh
        - applicable (registry order)
        - ``default_for``
      * - Cart-1D
        - ``{CumprodScan, ScanMarch, FullField}``
        - ``CumprodScan``
      * - Cart-2D
        - ``{ScanMarch, MovingFrontierWindow, FullField}``
        - ``ScanMarch``
      * - Cart-3D
        - ``{ScanMarch, FullField}`` (window is d=2 only)
        - ``ScanMarch``
      * - Cyl/Sph-1D
        - ``{CumprodScan, ScanMarch}``
        - ``CumprodScan``

   The **registry order is the policy**:

   .. code-block:: python

      LOSS_REPRESENTATIONS = (
          CumprodScan,            # 1-D production default
          ScanMarch,              # multi-D Cartesian production default
          MovingFrontierWindow,   # selectable peer
          FullFieldWavefront,     # never-stuck oracle fallback
      )

   At d=1 ``CumprodScan`` wins (registered first; ``ScanMarch`` would
   also apply but degenerates to the same scan with a march shell). At
   d≥2 ``ScanMarch`` wins — the **Fork-B2 flip** (2026-06-11) that moved
   ``ScanMarch`` ahead of ``MovingFrontierWindow`` in the tuple. The day
   a d=3 window lands, widening its ``supports`` is the *only* change
   needed for Cart-3D's available set to grow — one predicate, no caller
   touched.

#. **Construction guard** — ``_LossRepresentation.__post_init__``
   re-runs ``supports(mesh)`` and raises
   :class:`~orpheus.sn.loss_representation.IncompatibleRepresentation`
   on a false verdict, so even a bypassed UI cannot build an illegal
   pairing. Combined with the frozen-dataclass immutability, the
   ``(representation, mesh)`` pairing is *correct by construction*.

That ``supports`` predicate **is** the ``is_1d`` / ``curvature``
dispatch that the pre-carve code scattered across ``transport_sweep``
plus five operator gates — now declared **once** per representation.
"Add 3-D window support" becomes "extend one representation, widen one
predicate," not a hunt through every call site.


.. _loss-rep-one-walk-one-instance:

The one-walk and one-instance theorems
======================================

The structural payoff of the S6 re-layering is two theorems, each a
*type fact* enforced by construction and pinned by a discriminating spy
test (not a tautology — each test FAILED at the pre-carve HEAD and
flipped to PASS in its landing commit).

One walk (S6.4)
---------------

ONE d-generic frame —
:meth:`orpheus.sn.loss_representation._OctantWalk._interior_walk` —
serves **both** the sweep and the matvec for every multi-D
representation. For each octant it projects the quadrature octant to its
in-plane signs, branches the pure-z degenerate octants, derives the
per-axis in/out domain faces, reads the octant's inflow, runs the
interior traversal, and sheds the outflow. The two directions fork
**only** at:

* the **cell kernel** — the per-octant interior traversal the calling
  representation supplies (the window's frontier walk, the scan-march's
  row-march, the oracle's full cochain), in its solve
  (:meth:`~orpheus.sn.spatial.diamond.DiamondDifference.cell_kernel_batch`)
  or apply
  (:meth:`~orpheus.sn.spatial.diamond.DiamondDifference.residual_kernel_batch`)
  direction; and
* the **emit policy** — what the direction accumulates (the sweep's
  angular/moment output via ``_SweepEmit``; the matvec's
  :math:`(L+C)\psi` bulk plus the O.4b boundary defect).

.. admonition:: The direction is an OBJECT, never a boolean
   :class: warning

   The solve/apply fork is carried by **kernel and emit objects**
   (``_CellSolve`` / ``_CellResidual``,
   :class:`~orpheus.sn.loss_representation._SweepEmit`), never by an
   ``is_solve`` flag threaded down the walk. A boolean direction flag is
   the coding-elegance Smell #3 anti-pattern (special-case dispatch
   masquerading as abstraction). ``test_one_octant_walk.py`` carries an
   **AST tripwire** that parses ``_OctantWalk``'s source and fails if any
   ``is_solve`` / ``is_apply`` / ``is_matvec`` name appears — so the
   carve cannot silently degrade back into the flag pattern.

The discriminating test
``test_sweep_and_loss_action_hit_one_octant_walk`` spies the call-time
``self`` of ``_interior_walk`` and asserts that both ``L.apply`` and
``A.solve`` on a 2-D Cartesian mesh exercise it. The three-layer stack
beneath the walk (storage walk / level operation / pure kernel pair) is
documented at :ref:`sweep-dispatch-relayering`; the graph layer
(:class:`~orpheus.sn.sweep_graph.SweepDependencyGraph.for_shape`,
per-shape ``lru_cache`` of immutable ``MappingProxyType`` octant→DAG
maps) is family-owned, so the mesh stays pure geometry.

One instance (S6.5)
-------------------

The operator holds **one** representation instance —
:attr:`StreamingOperator.loss_representation
<orpheus.sn.operator.StreamingOperator.loss_representation>` (a
``cached_property`` = ``default_for(mesh)``) — consumed by:

* :meth:`~orpheus.sn.operator.StreamingOperator.apply` (the matvec
  :math:`(L+C)\psi`);
* :meth:`~orpheus.sn.operator.InvertibleOperator.solve` (the forward
  substitution :math:`(L+C)^{-1}q`), via the delegating
  :attr:`InvertibleOperator.loss_representation
  <orpheus.sn.operator.InvertibleOperator.loss_representation>`
  property; and
* the boundary Gauss–Seidel resolvent.

Because representations are stateless frozen dataclasses, "one instance"
is a **structural type-fact goal** — the L21 invariant becomes
construction-enforced — not a performance fix (the per-shape DAG cache
already de-duplicated the heavy state). The discriminating test
``test_apply_and_solve_share_one_representation_instance`` spies the
call-time ``self`` of *both* doors and asserts object identity; it
**failed pre-S6.5** (the doors each called ``default_for`` independently
— two distinct frozen-dataclass ids per solve), and flips PASS once both
consume the operator's one instance.

.. note::

   **Deliberate scope boundary.** The module-level
   :func:`~orpheus.sn.loss_representation.transport_sweep` REMAINS the
   operator-free functional entry — it still selects via ``default_for``
   because its one production caller (the ``solve_sn``
   post-convergence reconstruction) has no operator in scope. The
   one-instance theorem is about the *operator's* doors; an
   operator-free functional caller legitimately mints its own.


.. _loss-rep-fork-b2:

The S6.9 measurement and the Fork-B2 decision
=============================================

The carve built ``ScanMarch`` as an **opt-in** representation (Fork B1)
with the default unchanged, precisely so the default could be decided on
*measured* numbers rather than a plausibility argument — the governing
principle: *construct each strategy as general as its algorithm
naturally allows; select narrow; specialize only on measured internal
cost.* The S6.9 benchmark
(``derivations/diagnostics/diag_s69_scanmarch_vs_window_bench.py``,
median over repeats, ``python -O``; full table in #222 comment
4683241855) measured ``ScanMarch`` / ``MovingFrontierWindow`` ratios:

.. list-table:: ScanMarch ÷ MovingFrontierWindow (lower is faster / leaner; median, ``-O``)
   :header-rows: 1
   :widths: 30 16 16 16

   * - config
     - sweep
     - matvec
     - peak mem
   * - 24×24 LS4 2g
     - 0.61
     - 0.55
     - 0.98
   * - 48×48 LS4 2g
     - 0.57
     - 0.59
     - 0.99
   * - 96×96 LS4 2g
     - 0.58
     - 0.60
     - 0.99
   * - 48×48 LS8 2g
     - 0.69
     - 0.67
     - 0.99
   * - 48×48 LS16 2g
     - 0.84
     - 0.78
     - 0.99
   * - 48×48 LS8 4g
     - 0.71
     - 0.72
     - 0.99

End-to-end fixed-source (48×48 LS8 2G heterogeneous): **10.5 s vs
12.8 s = 0.82×**.

Three findings shaped the decision:

#. **No memory edge at d=2.** The rolling frontier and the row-march are
   *both* :math:`O(\text{row})` working set at d=2 (``peak ≈ 0.98–0.99``).
   The window's memory advantage only ever held against the
   **full-field oracle** (which is ~1.3–1.4× both). The window's reason
   for being the default — peak-memory — does **not** distinguish it from
   the row-march.

#. **The ScanMarch advantage narrows with angular order** (LS4 ~0.58 →
   LS16 ~0.84). The scan's per-ordinate closed-form work scales with the
   ordinate count :math:`N`, while the level-batched wavefront amortises
   across ordinates per anti-hyperplane — so the gap closes as :math:`N`
   grows. The win is real across the measured range but is not
   asymptotically unbounded.

#. **End-to-end dilution to 0.82×.** The per-sweep kernel win is
   amortised by the scattering / moment-projection / reflect overhead
   that surrounds the within-group solve, so the wall-clock win on a full
   fixed-source solve is smaller than the bare-sweep ratio.

The decision (Fork-B2, USER, 2026-06-11):

.. admonition:: Fork-B2 decision (verbatim)
   :class: important

   *"There is no need to retire the moving frontier window method. We can
   have multiple methods (this is the whole point of having them
   selectable), as long as they are different methods, with slightly
   different principles, and the code is proper. But you can flip the
   default to ScanMarch."*

So the flip is the **one-line registry reorder** (move ``ScanMarch``
ahead of ``MovingFrontierWindow`` in
:data:`~orpheus.sn.loss_representation.LOSS_REPRESENTATIONS`);
``default_for`` is unchanged (still first-supports-match); the 1-D
default is untouched (``CumprodScan`` stays first). The window is
**kept** as a genuinely-different schedule (anti-diagonal wavefront vs
row-march) over the same lower-triangular operator — its end-to-end
coverage rides the forced-window leg of
``test_scan_march_end_to_end.py`` plus the explicit ``window ≡ full``
oracles. This is the **architecture-over-implementation** discipline
(Cardinal Rule 2): selectability is the whole point, and two proper
methods with different principles are an asset, not redundancy to be
pruned.


.. _loss-rep-bit-vs-principled:

Verification architecture: bit-identity vs principled-equivalence
=================================================================

The carve's verification turns on the vv-principles distinction between
**bit-identity** (an implementation property) and
**principled-equivalence** (a math property), applied at two different
granularities.

Across *schedules*: principled-equivalence
------------------------------------------

Different schedules (CumprodScan vs ScanMarch vs the wavefront) are
**NOT bit-comparable**: the row-march and the anti-diagonal reconstruct
the same cell dependencies *in a different order*, and IEEE-754 addition
is non-associative, so the converged values differ at FP-association.
Demanding ``array_equal`` across schedules would be the *wrong* gate. The
cross-schedule gates are therefore:

* **nulp / absolute-tolerance oracle** —
  ``test_scan_march_equivalence.py`` pins ``ScanMarch`` against the
  unconditionally-stable
  :class:`~orpheus.sn.loss_representation.FullFieldWavefront` oracle
  (G2.c), on ≥2G heterogeneous **anisotropic** configs with **non-square**
  meshes (the x↔y swap moat — failure Mode 2). The tolerance is absolute
  (``rtol=1e-11``, ``atol=1e-12``), not nulp, because near-zero boundary
  shed elements amplify a ~1e-15 absolute difference into a spurious
  ~16000-ULP reading.

* **solver-tol Mode-9 FP-invariance** —
  ``test_scan_march_end_to_end.py`` drives the full production solvers
  with the schedules swapped and asserts the **converged fixed point** is
  schedule-invariant to solver tolerance (a row-march MUST NOT move
  :math:`\psi^*` or :math:`k` — only the per-sweep FP-association, which
  the outer iteration washes out). Per vv-principles failure Mode 9, the
  FP-invariance is verified on configs that *break* the degenerate
  coincidence, never the isotropic-reflective box:

  - **G4.a** — P1-anisotropic + heterogeneous (fuel|moderator) + vacuum
    streaming, with the non-flat degenerate-gate guard;
  - **G4.b** — all-reflective + a level-symmetric cubature, so the
    outflow shed order is load-bearing (the ERR-056 shared-face failure
    class); the honest **d=2 limitation** is stated in the test — the
    full diagonal-cubature shared-face stressor is a d=3 case, deferred
    until a 3-D quadrature exists;
  - **G6** — the closed-form
    :math:`\kinf = \lambda_{\max}(A^{-1}F)` (homogeneous, ≥2G — no
    1-group eigenvalue evidence, per the cardinal-rule degeneracy) plus
    SI ≡ Krylov flux-shape agreement on the heterogeneous non-flat 2G
    config.

  The forcing is a **context manager** (a fixture would force the
  *reference* leg too and compare the window to itself), with explicit
  non-vacuity counters asserting the forced path actually ran. Since the
  Fork-B2 flip the polarity is inverted: the *default* leg runs
  ``ScanMarch`` and the *forced* leg runs ``MovingFrontierWindow`` — so
  the window peer gets its end-to-end coverage from the same gate, and
  the FP-invariance claim (being symmetric) is unchanged in meaning.

.. note::

   **MMS reaches the flux-shape layer, not the eigenvalue.** The
   ScanMarch verification matches the vv-principles claim taxonomy: the
   oracle and FP-invariance gates establish *convergence-order* and
   *flux-shape* invariance under the schedule swap; the **eigenvalue**
   anchor is the *closed-form* :math:`\kinf` leg (a structural-independence
   ground), not the SI≡Krylov twin agreement (which is necessary but not
   sufficient — two ORPHEUS paths agreeing is cross-implementation
   agreement, not correctness).

Within a *schedule*: bit-identity
---------------------------------

Two storage *policies* of the **same** schedule must agree to the byte —
they run the identical cell math in the identical level order, differing
only in how much of the cochain they retain. These ARE ``array_equal``
oracles:

* **window ≡ full** — ``test_2d_full_field_oracle.py`` (end-to-end
  sweep + matvec) and ``test_sweep_graph_window_equivalence.py`` (the
  graph-level walk, d=1 / d=2 / synthetic d=3) pin
  :class:`~orpheus.sn.loss_representation.MovingFrontierWindow` exactly
  against :class:`~orpheus.sn.loss_representation.FullFieldWavefront`.

* **the kernel pair is the FP reduction tree of record** —
  ``test_cell_kernel_batch.py::TestKernelSourceOfRecord`` freezes a
  ``sha256`` of the source of
  :meth:`~orpheus.sn.spatial.diamond.DiamondDifference.cell_kernel_batch`
  and
  :meth:`~orpheus.sn.spatial.diamond.DiamondDifference.residual_kernel_batch`.
  Their explicit left fold ``((Σ_t + s_0) + s_1) + …`` is
  bit-identity-load-bearing: an algebraically-equivalent rearrangement
  (a ``sum()`` instead of the fold) passes every value-tolerance test yet
  silently invalidates the 1-ULP regression contract. This is the **one
  legitimate source-hash pin** in the SN stack — the kernels *are* the
  reduction tree every byte-identity anchor inherits from.

* **the affine converged-bytes golden** —
  ``test_affine_carve_bit_identity.py`` freezes a ``sha256`` of the
  production default's converged ``angular_flux`` / ``scalar_flux``
  bytes. At the Fork-B2 flip the four 2-D hashes were **regenerated in
  the flip commit** (a schedule change shifts the converged bytes at
  FP-association level — principled-equivalent, not a numerics change),
  with a history block naming the output-identity evidence (the G4
  Mode-9 gates + the G2.c nulp oracle); the 1-D slab hashes were
  byte-unchanged (the flip's blast-radius pin). The discipline is
  **regenerate-in-commit, never pin stale**.

Structural spies
----------------

The L21 theorems are pinned by the one-walk
(``test_one_octant_walk.py``) and one-instance
(``test_one_representation_instance.py``) call-time-``self`` spies, both
``foundation``-tagged software-structure invariants (no theory
``:label:``), both ``-O``-safe (``pytest.fail``, never bare ``assert``,
so they fire under the canonical ``python -O`` invocation — vv-principles
failure Mode 8).


.. _loss-rep-history:

History and rationale: what was tried, measured, decided
========================================================

The carve arc (the WHY of the final shape)
------------------------------------------

The carve replaced a *scattered, procedural* dispatch — the same 1-D vs
multi-D decision spelled three different ways (``transport_sweep``
branching on ``reduced is not None``; the matvec branching in five
operator gates on ``not is_1d``; the oracle reachable only through
hand-built test adapters) — with one polymorphic abstraction. The
phases, each independently bit-identical-gated:

.. list-table::
   :header-rows: 1
   :widths: 10 90

   * - Phase
     - What it did
   * - **S0**
     - ``test-architect`` verification plan (the proactive dispatch for
       an operator-algebra carve crossing a subsystem boundary).
   * - **S1**
     - The skeleton: protocol + ``_DAGWavefront`` base + the four
       leaves as **thin wrappers** over the existing sweep code +
       ``is_cartesian`` + ``supports`` / ``default_for`` / registry;
       ``transport_sweep`` rewired. Bit-identical.
   * - **S2**
     - The matvec side: ``strategy.residual`` collapses the **5 matvec
       gates** to one delegating call. Bit-identical.
   * - **S3**
     - Solve-vs-solve equivalence **retires the hand adapters**; the
       full-field spine becomes the genuine d-generic oracle.
   * - **S4**
     - The window generalised to ``frontier_dim = d−1`` (a point at d=1,
       a line at d=2, a surface at d=3); d=2 stays bit-identical.
   * - **S5.1**
     - ``ScanMarch`` built oracle-pinned as **opt-in** (Fork B1) — the
       default unchanged until measured.
   * - **S6**
     - The re-layering: S6.2 rename ``SweepStrategy → LossRepresentation``
       / ``residual → loss_action``; S6.3 the walk moved **off** the
       operator (``operator.py`` became pure algebra, the ``−C`` glue
       collapsed 5×→1×); S6.4 **one walk** + family-owned DAG +
       ``diamond.py`` pure kernel pair (``sweep.py`` dissolved,
       ``WavefrontFlux`` retired with the cochain succession); S6.5
       **one instance**.
   * - **S6.9**
     - Measure + the Fork-B2 default flip (this page's
       :ref:`loss-rep-fork-b2`).

Rejected alternatives (and why)
-------------------------------

* **An enum threaded into** ``transport_sweep``. Rejected: it adds a
  *second* branch axis to a function that already branches on
  dimensionality — cyclomatic complexity, not abstraction. The
  polymorphic-dispatch representation declares the choice once per type.

* **A boolean** ``is_solve`` **flag in the shared walk.** Rejected
  (coding-elegance Smell #3): the direction is carried by kernel and
  emit **objects**; an AST tripwire enforces the absence of the flag.

* **Retiring the window.** Rejected by the Fork-B2 decision: a
  genuinely-different schedule over the same operator is the *point* of
  selectability, and the window is correct and proper on its own.

The ``WavefrontFlux`` succession
--------------------------------

The S6.4 carve **retired** the typed ``WavefrontFlux`` field and its
``InteriorFaceSpace``, but the *concept* — the interior 1-cochain
:math:`C^1_{\rm int}` (with :math:`C^1 = C^1_{\rm int}\oplus C^1_\partial`
remaining valid theory) — survives in two realisations: the rolling
``_MovingFrontier`` front (the window) and the
``_octant_face_cochain`` history (the oracle). The cochain theory anchor
:ref:`wavefront-flux-cochain` is **kept**; its derivation is preserved as
the theory of both realisations.

The deferred extension point
----------------------------

The designed-but-deferred fifth representation is an **ExplicitMatrix**:
the sparse-assembled :math:`(L+C)`, whose ``sweep`` is
``scipy.sparse.linalg.spsolve_triangular``. It is the proof that the
representation abstraction is genuinely a *set of schedules over one
operator* — an assembled lower-triangular matrix is the most literal
realisation of the native frame. It slots into the registry with one
``supports`` predicate when a use case (e.g. a teaching demonstration of
the triangular structure, or a direct-solve cross-check) motivates it.


Literature
==========

The cell-level mathematics these representations *schedule* is sourced
from the standard discrete-ordinates references, anchored in the
:class:`~orpheus.sn.spatial.diamond.DiamondDifference` docstring and
:doc:`discrete_ordinates`:

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Reference
     - Used for
   * - Lewis & Miller (1984), *Computational Methods of Neutron
       Transport*, §4.5 / §5.3 / §6.4
     - The Morel–Montry angular closure (§4.5), the Diamond / weighted-DD
       / Step / LD closures and the negative-flux failure mode (§5.3),
       and the canonical MMS ansatz set (§6.4).
   * - Hébert (2009), *Applied Reactor Physics*, Ch. 3 §3.9.4
     - The curvilinear S\ :sub:`N` cell-balance and DD difference
       relations.
   * - Blelloch (1990), *Prefix Sums and Their Applications*,
       CMU-CS-90-190 §1.5
     - The first-order linear-recurrence scan
       (:func:`~orpheus.sn.spatial.scan.ordinate_scan`) that
       ``CumprodScan`` and ``ScanMarch`` both reuse — the pair-monoid
       closed form.
   * - Adams & Larsen (2002), *Fast iterative methods for
       discrete-ordinates particle transport calculations*, §III
     - The within-group iteration framing in which the sweep
       :math:`(L+C)^{-1}` is the transport operator and the matvec
       :math:`(L+C)` its twin.


See also
========

* :doc:`discrete_ordinates` — the cell-level WDD algebra, the 1-D
  cumprod recurrence (:ref:`sweep-cumprod`), the 2-D anti-diagonal
  wavefront (:ref:`sweep-wavefront`), and the three-layer walk /
  level-op / kernel stack (:ref:`sweep-dispatch-relayering`).
* :doc:`operator_algebra` — the Wave-O typed operator algebra
  (:eq:`operator-apply`, :eq:`operator-solve`,
  :eq:`operator-apply-transpose`) and the interior face-flux cochain
  (:ref:`wavefront-flux-cochain`).
* :mod:`orpheus.sn.loss_representation` — the module that implements
  every representation, the selection layer, and the orchestration
  entry.
* ``.claude/plans/sn_sweep_strategy.md`` — the authoritative locked
  design (decisions, verification strategy, phases S0–S6).
