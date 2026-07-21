.. _sn-cartesian-multid:

Cartesian multi-D: space enters the walk
========================================

This chapter broadens exactly one axis of the slab chapters:
**space**. Streaming becomes a true directional gradient — two (in 3-D,
three) sink terms instead of one — and with it the sweep stops being a
chain and becomes a **wavefront over a causal dependency graph**. What
does *not* change is the algebra: Cartesian geometry keeps a neutron's
direction constant in flight, so there is still no angular coupling
outside the sources, the group axis of :doc:`slab_multigroup` rides
along untouched, and the within-group operator keeps its honest shape
:math:`A = L + C - S - B` with :math:`L+C` invertible in one pass.
The chain of the book repeats on the new axis:

1. **the invariant** — sinks = sources, now on a rectangular cell with
   faces on every axis → *pose* the 2-D balance;
2. **the operator** — :math:`L+C` remains **lower-triangular in sweep
   order**, but "sweep order" is now a *partial* order: each octant
   induces a causal DAG whose anti-diagonal levels are mutually
   independent;
3. **the matrix picture** — the loss inverse factors as a **direct sum
   over octants** (:eq:`streaming-inverse-direct-sum`), each block a
   forward substitution over its DAG;
4. **the strategy-encoding operators** — the wavefront walk realized as
   frozen primitives (graph, storage walks, level operations, kernel
   pair), and the representation layer that selects among schedules
   (:doc:`loss_representation`).

The angular quadrature now genuinely uses both direction cosines
(:doc:`angular_quadrature`); the generic discretization machinery is
still :doc:`/theory/foundations/discretization`, cross-linked never
re-derived.

.. admonition:: Key Facts
   :class: tip

   * The 2-D transport equation is :eq:`transport-cartesian-2d` — two
     streaming sinks, and **still no angular coupling between
     ordinates**: space broadens the walk, not the angle algebra.
   * The 2-D DD balance :eq:`dd-cartesian-2d` applies the diamond
     closure on **both** axes simultaneously; the streaming
     coefficients are :math:`s_x = 2|\mu_x|/\Delta x`,
     :math:`s_y = 2|\mu_y|/\Delta y`.
   * Cells on an anti-diagonal :math:`i + j = k` share no faces —
     they are solved **simultaneously**. The sweep is a per-octant
     batched forward substitution over a precomputed causal DAG, and
     the loss inverse is the **direct sum over octants**
     :eq:`streaming-inverse-direct-sum`.
   * The walk factors into three layers — **storage walk × level
     operation × kernel pair** — so a closure supplies only its
     storage-free cell algebra and inherits both storage policies
     (full cochain / rolling frontier) and both directions
     (solve / apply) for free.
   * The multi-D LD closure is the **tensor-product bilinear (UBLD)**
     :math:`\{1, x, y, xy\}` cell system :eq:`ld-ubld-cell-system` —
     the :math:`xy` cross moment is diffusion-limit-load-bearing
     (simplex-P1 fails the thick limit on quadrilaterals), and sweep
     and matvec share one :math:`d`-generic kernel through the octant
     moment-frame involution :eq:`ld-ubld-octant-moment-frame-signs`
     (ERR-061).
   * Boundary conditions apply **once per octant per axis**, never per
     ordinate — the L7 trap (ERR-003): per-ordinate BC application is
     redundant in cost and order-sensitive in correctness.
   * The per-cell operation order is **bit-identity-load-bearing**:
     different *schedules* of the same operator are compared by
     principled-equivalence gates; different *storage policies* of one
     schedule are bit-identical (``array_equal``).


The posing: a second streaming sink
===================================

The invariant is unchanged — **sinks = sources** on every region of
phase space (:doc:`slab_one_group`, The Posing). Space enters through
the streaming term alone: on a 2-D Cartesian phase space the beam at
:math:`(x, y, \hat\Omega_n)` leaks through faces on both axes.

In two Cartesian dimensions the angular flux depends on two direction
cosines :math:`\mu_x` and :math:`\mu_y`:

.. math::
   :label: transport-cartesian-2d

   \mu_x \frac{\partial \psi}{\partial x}
   + \mu_y \frac{\partial \psi}{\partial y}
   + \Sigt{} \, \psi
   = \frac{Q}{W}

There is no angular coupling between ordinates --- each direction is
solved independently.  The two streaming terms are the only difference
from the 1D case.


Nothing else moved. Collision is the same multiplication operator,
scattering and fission the same group-coupling operators of
:doc:`slab_multigroup` — the sources see a longer spatial index, not a
new structure. This is the payoff of the Cartesian structural fact:
**the spatial axis broadens the walk, and only the walk.**


.. _balance-cartesian-2d:

The discrete balance on a rectangular cell
==========================================

The cell balance now carries face terms on every axis — the
one-equation-*three*-unknowns shape (cell average + two downstream
faces) that the closure must reduce per axis
(:doc:`/theory/foundations/discretization` §3):

Integrating :eq:`transport-cartesian-2d` over a rectangular cell
:math:`\Delta x_i \times \Delta y_j`:

.. math::

   \mu_{x,n}\bigl[\psi_{i+\frac12,j} - \psi_{i-\frac12,j}\bigr] \Delta y_j
   + \mu_{y,n}\bigl[\psi_{i,j+\frac12} - \psi_{i,j-\frac12}\bigr] \Delta x_i
   + \Sigt{} \Delta x_i \Delta y_j\, \psi_{n,i,j}
   = S_{i,j}\, \Delta x_i \Delta y_j

Dividing through by :math:`\Delta x_i \Delta y_j` and applying
diamond-difference closures in **both** directions simultaneously:

.. math::

   \psi_{n,i} &= \tfrac{1}{2}(\psi^x_{\rm in} + \psi^x_{\rm out})
   \qquad\text{(x-closure)} \\
   \psi_{n,i} &= \tfrac{1}{2}(\psi^y_{\rm in} + \psi^y_{\rm out})
   \qquad\text{(y-closure)}

yields the 2D DD equation:

.. math::
   :label: dd-cartesian-2d

   \psi_{n,i,j}
   = \frac{S_{i,j}
     + s_x\, \psi^x_{\rm in}
     + s_y\, \psi^y_{\rm in}}
     {\Sigt{} + s_x + s_y}

where the streaming coefficients are:

.. math::

   s_x = \frac{2|\mu_{x,n}|}{\Delta x_i}, \qquad
   s_y = \frac{2|\mu_{y,n}|}{\Delta y_j}

Both outgoing face fluxes are then updated from the DD closure:

.. math::

   \psi^x_{\rm out} = 2\psi_{n,i,j} - \psi^x_{\rm in}, \qquad
   \psi^y_{\rm out} = 2\psi_{n,i,j} - \psi^y_{\rm in}

These are precomputed by :class:`SNMesh` as ``streaming(0)[n, i]`` and
``streaming(1)[n, j]``, so the inner loop in
:func:`_sweep_jacobi` reduces to a single vectorised division per
diagonal.


.. _sweep-wavefront:

Cartesian 2D: Anti-Diagonal Wavefront Sweep
=============================================

In 2D, the DD equation :eq:`dd-cartesian-2d` creates a data dependency:
cell :math:`(i, j)` requires incoming face fluxes from its upwind
neighbours in both :math:`x` and :math:`y`.  Cells along an
**anti-diagonal** :math:`i + j = k` are mutually independent because
they share no incoming faces, so they can be solved simultaneously.

The wavefront sweep is implemented as a **per-octant batched
forward-substitution** over a precomputed causal cell DAG (Wave 2 of
the SN performance plan, closing Issue #4).  This subsection states
the algebraic framing; the primitives that realise it
(:class:`~orpheus.sn.loss_representation.sweep_graph.OctantLabel`,
:class:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph` and its two
storage walks, the level-operation pair ``_CellSolve`` /
``_CellResidual``, and the discretization's kernel pair
:meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.cell_kernel_batch`
/ :meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.residual_kernel_batch`)
are documented in detail at
:ref:`sweep-octant-dependency-graph` immediately below.

The §15A.2 sum-of-tensor-products framing
-----------------------------------------

Following Grand Report v3 §15A.2 (lines 2137–2171), the loss inverse
(the sweep) :math:`(L+C)^{-1}` on the 2-D Cartesian SN field
space decomposes as a **direct sum over angular octants** — the block
structure is streaming-induced, since each octant sweeps in a fixed
direction (in this section's equations :math:`A` abbreviates the loss
composite :math:`L+C`, the invertible sub-composite of the chapter's
:math:`A = L+C-S-B`):

.. math::
   :label: streaming-inverse-direct-sum

   A^{-1} \;=\; \bigoplus_{\sigma \in \mathcal{O}} A^{-1}_{\sigma},
   \qquad
   \mathcal{O} \;=\; \{\sigma = (\mathrm{sgn}\,\mu_x,\,
                                  \mathrm{sgn}\,\mu_y) :
                       \sigma \neq (0,0)\}
                  \,\cup\, \{(0,0)\},

.. vv-status: streaming-inverse-direct-sum documented

acting on the octant-restricted tensor space :math:`(N_\sigma,\,n_x,\,n_y,\,n_g)`.
The direction-cosine partition (Eq. :eq:`octant-sign-predicate`) is
the predicate the
:class:`~orpheus.numerics.quadrature.Quadrature` class exposes as
its cached :attr:`~orpheus.numerics.quadrature.Quadrature.octants`
property — a tuple of
:class:`~orpheus.numerics.measure.DiscreteMeasurePartition`
entries realised by
:meth:`~orpheus.numerics.measure.DiscreteMeasure.partition_by`
(see :ref:`tensorial-framing` and the
:doc:`/theory/foundations/discrete_measures` consumer table).

For each non-degenerate octant :math:`\sigma`, the action of
:math:`A^{-1}_\sigma` is a **forward substitution along a per-octant
causal cell DAG** — the topological order is structural (anti-diagonal
sweep on the Cartesian grid), the per-level cell update is one
vectorised einsum.  The pure-:math:`z` degenerate octant
:math:`\sigma = (0,0)` (ordinates with :math:`\mu_x = \mu_y = 0`,
which appear in 3-D angular cubatures projected to the in-plane
2-D problem) has no spatial streaming and reduces to a per-cell
balance :math:`\psi = Q / \Sigma_t` — the wavefront sweep handles
it via a short-circuit and skips the dependency graph.

**The four quadrant sweeps.**  Each non-degenerate octant
:math:`\sigma = (\mathrm{sgn}\,\mu_x, \mathrm{sgn}\,\mu_y) \in \{-1,+1\}^2`
determines a sweep direction:

.. list-table::
   :header-rows: 1
   :widths: 20 20 30 30

   * - :math:`\mu_x`
     - :math:`\mu_y`
     - *x*-direction
     - *y*-direction
   * - :math:`+`
     - :math:`+`
     - left :math:`\to` right
     - bottom :math:`\to` top
   * - :math:`-`
     - :math:`+`
     - right :math:`\to` left
     - bottom :math:`\to` top
   * - :math:`+`
     - :math:`-`
     - left :math:`\to` right
     - top :math:`\to` bottom
   * - :math:`-`
     - :math:`-`
     - right :math:`\to` left
     - top :math:`\to` bottom

For each octant, the sweep visits topological levels
(anti-diagonals) :math:`k = 0, 1, \ldots, n_x + n_y - 2`.  On level
:math:`k`, the cells :math:`(i, j)` satisfying :math:`i + j = k`
(in the per-octant traversal index space) are gathered into a numpy
batch and solved with a single vectorised evaluation of
:eq:`dd-cartesian-2d` — vectorised across the **ordinate axis**
(:math:`N_\sigma` — every ordinate in the octant), the
**anti-diagonal axis** (:math:`n_{\rm diag}` — number of cells on
this level), and the **group axis** (:math:`n_g`) simultaneously.

**Vectorisation within each level.**  Each level contains up to
:math:`\min(n_x, n_y)` cells.  The incoming face fluxes
``psi_in_x`` and ``psi_in_y`` are gathered by advanced indexing;
the DD equation is evaluated as one numpy operation; and the
outgoing face fluxes are scattered back into the persistent face-
flux buffers.  There is **no Python-level cell loop within a level**
and **no Python-level ordinate loop within an octant** — both axes
are internal to the einsum.

**Reflective BCs in 2D.**  At each boundary face, the incoming flux
for ordinate :math:`n` is set to the outgoing flux of its reflected
partner.  For the left/right boundaries (*x*-reflection), the partner
is ``ref_x[n]`` (negating :math:`\mu_x`); for the top/bottom boundaries
(*y*-reflection), the partner is ``ref_y[n]`` (negating :math:`\mu_y`).
The reflection indices are precomputed by the quadrature's
:meth:`reflection_index` method.  Crucially, the BC apply happens
**once per octant per axis** (not once per ordinate per axis) —
see :ref:`sweep-octant-dependency-graph-l7-trap` for the rationale.

Implemented in :func:`~orpheus.sn.loss_representation._sweep_jacobi`, which
is a thin orchestrator over the
:class:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph` primitives
described next.


.. _sweep-octant-dependency-graph:

Cartesian 2D: Octant Dependency Graph (Wave 2)
===============================================

This section documents the **§15A.2 "upwind trace complex / causal
transport DAG / direction sweep ordering" primitive** as it lives in
:mod:`orpheus.sn.loss_representation.sweep_graph` after Wave 2 of the SN performance plan
(branch ``feature/sn-octant-sweep-graph``, closes Issue #4).  The
shipped architecture replaces the legacy per-ordinate ``for n in
range(N)`` loop in :func:`~orpheus.sn.loss_representation._sweep_jacobi` with
a per-octant batched dispatch, lifting the per-call ``_diag_cache``
build to mesh-time work, and isolating the per-cell DD algebra in the
discretization's pure kernel pair
(:meth:`~orpheus.transport.spatial.diamond.DiamondDifference.cell_kernel_batch`
/ :meth:`~orpheus.transport.spatial.diamond.DiamondDifference.residual_kernel_batch`)
that LD / EC / Step closures can override later.

.. note::

   **Architecture history — the dispatch surface re-layered twice.**
   Wave 2 (the original closure of Issue #4) routed the sweep through
   a per-level *packet* (the ``SweepCellSlice`` dataclass) consumed by
   four direction×storage methods — ``update_batch`` / ``residual_batch``
   on the strategy (full-field) plus their ``apply_windowed`` /
   ``residual_windowed`` siblings on the graph.  S6.4(e) **collapsed
   that surface**: the four walk methods became TWO storage walks
   (:meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_full`,
   :meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_windowed`)
   each parameterised by a level-operation OBJECT (``_CellSolve`` for
   the solve direction, ``_CellResidual`` for the apply direction —
   direction is never a boolean flag); the per-level ``SweepCellSlice``
   packet was retired (it existed only to feed the now-deleted storage
   adapters); and the strategy's ``update_batch`` / ``residual_batch``
   were replaced by the **storage-free kernel pair**
   ``cell_kernel_batch`` / ``residual_kernel_batch`` (pure cell
   algebra — no gather/scatter).  The historical names ``update_batch``
   / ``residual_batch`` / ``SweepCellSlice`` appear below only as
   *history*; the current contract is the kernel pair + the level
   operations.  See :ref:`sweep-dispatch-relayering` for the WHY.

The primitives
--------------

The architecture is a small set of frozen, individually unit-tested
primitives plus a mesh-time precompute step.

.. list-table::
   :header-rows: 1
   :widths: 28 16 56

   * - Primitive
     - Lives in
     - Role
   * - :class:`~orpheus.sn.loss_representation.sweep_graph.OctantLabel`
     - :mod:`orpheus.sn.loss_representation.sweep_graph`
     - Frozen + slotted dataclass carrying one direction sign per
       spatial axis (``signs[axis] ∈ {-1, 0, +1}``) — a single type
       labels a 1-D (``(±1,)``), 2-D (``(±1, ±1)``), or 3-D octant.
       Hashable; used as the key in the per-shape graph family
       :meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.for_shape`
       (owned by the ``_DAGWavefront`` representation family since
       S6.4(c) — historically a mesh attribute).  An all-zero
       signature denotes the pure-:math:`z` degenerate octant — no
       graph is built for it
       (:attr:`~orpheus.sn.loss_representation.sweep_graph.OctantLabel.streams` is
       ``False``).  The 3-D ``sign_z`` is dropped by the 2-D Cartesian
       orchestration: the in-plane sweep is invariant under the
       out-of-plane axis, so multiple ordinates with the same in-plane
       ``signs`` but different ``sign_z`` share a single graph instance.
   * - :class:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph` (+ its
       two storage walks)
     - :mod:`orpheus.sn.loss_representation.sweep_graph`
     - Frozen dataclass holding the per-octant topological levels
       (anti-diagonals) and the per-axis face-index offsets.  Built
       once per ``(shape, octant)`` pair in the
       :meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.for_shape`
       cache (S6.4(c); historically at mesh construction); reused
       across every source iteration / Krylov matvec / outer
       iteration.  Exposes TWO storage walks (S6.4(e)):
       :meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_full`
       carries the COMPLETE per-axis interior face cochain (the
       verification-oracle policy);
       :meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_windowed`
       advances a rolling :math:`(d{-}1)`-frontier window (the
       production policy, ``O(N·n_g·∏ n_a)`` shrunk to
       ``O(N·n_g·∏_{a<d−1} n_a)`` backing).  The walk owns the level
       loop, the storage, and the per-level operand extraction; it
       dispatches the cell algebra to a level operation (next two rows).
   * - The level-operation pair ``_CellSolve`` / ``_CellResidual``
     - :mod:`orpheus.sn.loss_representation.sweep_graph`
     - The **direction fork, as OBJECTS** (S6.4(e); direction is never
       a boolean flag).  Exactly ONE is constructed per octant walk; the
       storage walk calls ``level_op.cell(...)`` per topological level.
       ``_CellSolve`` runs the solve direction — calls the strategy's
       ``cell_kernel_batch`` then performs the Phase-5c angular-XOR-
       moment per-level emit (write the angular flux + accumulate the
       scalar flux, OR accumulate the harmonic-moment tensor, never
       both).  ``_CellResidual`` runs the apply direction — calls
       ``residual_kernel_batch`` then writes the per-level residual.
       The per-level *emit* expressions and their order are
       bit-identity-load-bearing — relocated verbatim from the four
       retired walk methods.
   * - The kernel pair
       :meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.cell_kernel_batch`
       / :meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.residual_kernel_batch`
     - :mod:`orpheus.transport.spatial.scheme`
     - The **storage-free extension point** on the
       :class:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase` ABC
       (S6.4(e); historically the ``SweepCellSlice``-packeted
       ``update_batch`` / ``residual_batch``).  Each takes the per-axis
       incoming face fluxes + streaming coefficients + the level's cross
       section and source and returns ``(psi_avg, psi_out)`` (solve) or
       ``(residual, psi_out)`` (apply) — PURE cell algebra, no
       gather/scatter (that is the walk's job).  Default raises
       :exc:`NotImplementedError` — additive capability, not a contract
       change.  :class:`~orpheus.transport.spatial.diamond.DiamondDifference`
       overrides the pair; LD / EC / Step closures override it later to
       join the batched wavefront walks (their per-cell
       :meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.update`
       stays the canonical reference contract).

Per-shape precompute pattern (family-owned since S6.4(c))
---------------------------------------------------------

The dependency graph is a **derived object** — the
``(shape × octant)`` joint property.  It depends only on cell topology
and the octant sign convention; it does **not** depend on fluxes,
sources, BCs, quadrature, cross sections, or iteration state.  So the
graph build is paid once **per spatial shape** in the cached accessor
:meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.for_shape`, owned
by the DAG-consuming ``_DAGWavefront`` representation family:

.. code-block:: python

   class _DAGWavefront(_LossRepresentation):
       @property
       def sweep_graphs(self) -> dict[OctantLabel, SweepDependencyGraph]:
           # cached per shape; same-shape meshes share byte-identical graphs
           return SweepDependencyGraph.for_shape(self.mesh.spatial_shape)

**Ownership history** (two relocations, each a pure refactor):

#. *Wave 2 / C2.4* lifted the per-call ``_diag_cache`` build that
   previously lived inside the 2-D wavefront sweep (rebuilt once per
   sweep call) to **mesh-construction** time — a measurable but
   second-order saving on the 421-group benchmark; the structurally
   important effect was making the graphs named, inspectable state.
#. *S6.4(c)* moved ownership **off the mesh onto the representation
   family**: the mesh is pure geometry, and only the two DAG-walking
   representations (the window + the full-field oracle) ever mention
   the substrate.  This retired the curvilinear
   ``mesh.sweep_graphs = None`` slot — an illegal state (a mesh
   carrying a "no DAG here" marker for a structure it never owned) —
   and replaced mesh-lifetime caching with per-SHAPE caching, so
   same-shape meshes share one graph family (the graphs carry no
   mesh-identity information).  DAG-free representations
   (``CumprodScan``, ``ScanMarch``) and curvilinear meshes simply
   never touch the accessor; curvilinear sweeps walk the cell graph
   differently (per-ordinate march; see
   :meth:`~orpheus.sn.mesh.augmented_mesh.SNMesh.dag_walk`).

The closed-form precompute lives in
:meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.from_cartesian`
and never appears in the sweep loop.  This is structural, not
hand-rolled — the "library version" (a generic topological-sort over
an explicit DAG) would be over-engineering for a regular pattern that
collapses to ~5 lines of ``arange`` + anti-hyperplane extraction.  The
builder is dimension-generic (``d = len(shape) ∈ {1, 2, 3}``): a d=1
chain, the d=2 anti-diagonal, a d=3 anti-hyperplane.

The §15A.2 invariant set
-------------------------

The Grand Report v3 §15A.2 (lines 2165–2171) prescribes a fixed set
of L0 invariants every ``SweepDependencyGraph`` instance must satisfy.
These are pinned by ``tests/sn/test_sweep_graph.py`` (63 L0 tests):

* **Upwind orientation** — for each octant
  :math:`\sigma = (\mathrm{sgn}\,\mu_x, \mathrm{sgn}\,\mu_y)`, the
  ``face_in_x`` and ``face_out_x`` offsets satisfy
  ``face_in_x + face_out_x == 1`` and
  ``face_in_x = 0`` iff :math:`\mathrm{sgn}\,\mu_x \ge 0` (and
  analogously on :math:`y`).  Asserted by
  ``test_face_pairing_consistent`` and ``test_upwind_orientation``.
* **Topological sort** — every level's cells depend only on cells in
  strictly earlier levels (under the per-octant orientation).  No
  intra-level dependencies; no back-edges.  Asserted by
  ``test_topologically_sorted``.
* **Cell coverage** — every cell :math:`(i, j) \in [0, n_x) \times
  [0, n_y)` appears in **exactly one** level.  Disjoint union over
  the topological levels reconstructs the full grid.  Asserted by
  ``test_cell_coverage``.
* **Face-pairing consistency** — the incoming-face index of cell
  :math:`(i, j)` on level :math:`k` matches the outgoing-face index
  of its upwind neighbour on level :math:`k - 1` (under the per-
  octant orientation).  Asserted by
  ``test_face_pairing_consistent``.

These four invariants are the **load-bearing correctness floor** of
the wavefront sweep.  Any future closure (LD, EC, Step) plugged in
via the kernel pair consumes the same invariants — they describe
the topology, not the algebra, so the strategy contract is orthogonal
to the graph correctness.

.. _sweep-dispatch-relayering:

The dispatch boundary: walk (scheduler) vs cell update (closure)
-------------------------------------------------------------------

A central architectural decision is the **separation between the
scheduler and the closure**.  Three layers stack from storage outward
to algebra (S6.4(e)):

#. **The storage walk** —
   :meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_full` (full
   cochain) or
   :meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_windowed`
   (rolling frontier).  Owns the topological-level loop and the
   per-axis face gather/scatter (full cochain) or the frontier
   seed/incoming/emit/shed cochain trace algebra (window).  Storage is
   the walk's concern — the SAME two walks serve every closure and
   both directions.

#. **The level operation** — ``_CellSolve`` or ``_CellResidual``,
   constructed once per octant walk and called as ``level_op.cell(...)``
   per level.  Owns the direction fork (solve vs apply) and the
   per-level *emit* (angular/moment write, or residual write).
   Direction is an OBJECT, never a boolean flag passed down the walk.

#. **The kernel pair** —
   :meth:`~orpheus.transport.spatial.diamond.DiamondDifference.cell_kernel_batch`
   /
   :meth:`~orpheus.transport.spatial.diamond.DiamondDifference.residual_kernel_batch`.
   Owns the **pure cell algebra** and nothing else — no gather, no
   scatter, no storage. This is the ONLY direction-aware math left in
   the SN spatial stack.

**Why this layering (the WHY behind S6.4(e)).**  Wave 2 carried the
storage concern *inside* the strategy: the DD ``update_batch`` /
``residual_batch`` methods gathered the cell's face inputs from the
``SweepCellSlice`` packet, ran the algebra, and scattered the outgoing
faces back — a four-method direction×storage product
(``update_batch`` / ``residual_batch`` full-field +
``apply_windowed`` / ``residual_windowed`` windowed).  That entangled
two orthogonal concerns: a NEW closure (LD / EC / Step) would have had
to re-implement the gather/scatter (storage) plumbing four times just
to supply its cell math.  S6.4(e) lifts storage to the walk layer
**once, above every strategy**, so a closure supplies ONLY its
storage-free kernel pair (pure algebra over the per-axis incoming face
fluxes) and inherits both storage policies (full + window) and both
directions (solve + apply) for free.  The ``SweepCellSlice`` packet —
which existed only to feed the retired storage adapters — is gone with
them.  This is the Cardinal-Rule-2 "build primitives, not products"
discipline: the four-method product collapses to a 2 (walks) × 1
(level-op pair, direction-by-object) × 1 (kernel pair) factoring where
each factor varies independently.

This means: **DD is the only shipping closure today**, but Step / LD
/ EC override the kernel pair later without touching the walk driver
or the level operations.  The Wave C-extension rollout (Issues #157 /
#158) ships the per-cell :meth:`update` method first as the canonical
reference contract; the batched kernel pair is the parallel
level-vectorised capability for closures whose per-cell algebra
vectorises across an ``(N_oct, n_diag, ng)`` slice without per-cell
branching.

The DD ``cell_kernel_batch`` reproduces the legacy 2-D wavefront DD
math **bit-identically** (operation order matters; see
:class:`~orpheus.transport.spatial.diamond.DiamondDifference` docstring on
bit-identity).  The math is the **balance form** of WDD on a 2-D
Cartesian cell:

.. math::
   :label: dd-2d-balance-form

   \overline{\psi}_{i,j}
   \;=\; \frac{Q_{i,j}
               + s_{x,i}\,\psi^{\rm in}_{x,i,j}
               + s_{y,j}\,\psi^{\rm in}_{y,i,j}}
              {\Sigma_{t,i,j} + s_{x,i} + s_{y,j}},
   \qquad
   s_{x,i} = \frac{2|\mu_x|}{\Delta x_i},
   \quad s_{y,j} = \frac{2|\mu_y|}{\Delta y_j},

.. vv-status: dd-2d-balance-form documented

with the spatial closure
:math:`\psi^{\rm out}_x = 2\overline{\psi} - \psi^{\rm in}_x`
(and analogously on :math:`y`).  The operation order is fixed at:

.. code-block:: text

   denom    = sig_t + sx + sy
   psi_avg  = (Q + sx * psi_in_x + sy * psi_in_y) / denom
   psi_out  = 2 * psi_avg - psi_in

Algebraically equivalent rearrangements (e.g., reordering
``sig_t + sx + sy`` to ``sx + sy + sig_t``) break the 1-ULP regression
contract even though the math is identical.  This is the canonical
**bit-identity vs principled-equivalence** instance from the
``vv-principles`` skill — the regression contract is bit-identity
gated on the existing snapshots; deviations are admissible only when
backed by a structurally-independent reference (e.g., :math:`k_\infty`
analytical limit on a homogeneous reflective problem).

.. _sweep-octant-dependency-graph-l7-trap:

The L7-trap fix: BC apply once per octant
------------------------------------------

Wave 2 closes a class of bugs that the test-architect dispatch
identified as the **L7 trap** — the design pattern where a sweep
driver re-applies a boundary operator at each ordinate iteration.
The legacy ``_sweep_jacobi`` had this shape:

.. code-block:: python

   # legacy — the L7 trap
   for n in range(N):
       psi_x = sn_mesh.bc_xmin.apply(psi_x, quad)[n]   # per-ordinate apply
       psi_y = sn_mesh.bc_ymin.apply(psi_y, quad)[n]   # per-ordinate apply
       # ... walk cells, sweep, etc. ...

Each ``bc.apply`` call sees the FULL ``(N, ny, ng)`` face buffer (so
reflective partners can read across rows) and returns an updated full
buffer.  Calling this :math:`N` times per sweep is wrong on two
counts:

1. **Cost** — :math:`N` redundant invocations of the same boundary
   operator on the same buffer.  For a 2-D ``LS-N`` quadrature with
   :math:`N \sim 30`–:math:`80`, this is the dominant per-sweep cost
   on small meshes.
2. **Correctness** — when reflective BCs interact with mid-sweep
   reflective-buffer state, **the order of BC apply vs ordinate
   sweep matters**.  The legacy code's behaviour is sensitive to
   ordinate iteration order; reorderings that algebraically should
   be no-ops (e.g., octant batching, parallel ordinate evaluation)
   silently change the converged solution.

The Wave-2 form applies each boundary operator **once per octant per
axis** — :math:`O(\text{octants}) = 4` calls, not :math:`O(N)`:

.. code-block:: python

   # Wave 2 — L7-trap closed by construction
   for octant in quad.octants:                    # 4 iterations, structural
       sx, sy = octant.label
       ...
       # Apply BC once for this octant on each axis
       if sx_eff >= 0:
           full_face_x = sn_mesh.bc_xmin.apply(psi_x[:, 0, :, :], quad)
           psi_x[oct_idx, 0, :, :] = full_face_x[oct_idx]
       else:
           full_face_x = sn_mesh.bc_xmax.apply(psi_x[:, nx, :, :], quad)
           psi_x[oct_idx, nx, :, :] = full_face_x[oct_idx]
       # ... analogously on y ...
       sweep_graph.walk_windowed(level_op=_CellSolve(...), ...)  # all N_oct batched

The architectural argument: the boundary operator's *semantics* are
"map outgoing partner-octant fluxes to incoming this-octant fluxes".
That mapping is per-octant by construction — applying it once per
ordinate within an octant is redundant; applying it once per octant
is structurally correct.

.. note::

   The ``sn_mesh.bc_xmin.apply(..., quad)`` spellings in the two
   code blocks above are **historical** (the Wave-2 era 2-arg
   ``apply`` on a per-attribute BC surface). Both spellings are
   retired: the 2-arg ``apply`` by Issue #186 (the law is now a pure
   descriptor; the realizer produces a strict 1-arg operator — see
   :ref:`bc-trace-law-descriptor-model`), and the per-attribute
   ``bc_<face>`` surface by C4 / #220 in favour of the
   face-name-keyed :attr:`SNMesh.bc` dict
   (``sn_mesh.bc["xmin"].apply(psi)`` — see
   :ref:`bc-face-name-carve`). The blocks are preserved verbatim
   because they document the *L7-trap structure* the Wave-2 carve
   closed, which is independent of the storage spelling.

The L7-trap detector test
``tests/sn/test_2d_octant_sweep_equivalence.py::case-3`` is the
load-bearing regression gate — a TESTS-FIRST harness (case 3 with
mixed reflective + vacuum BCs, 2G heterogeneous, ``n_sweeps=2``)
designed to fail if any future refactor reintroduces the per-ordinate
BC apply pattern.  The case-3 design uses the post-sweep ``psi_x`` /
``psi_y`` buffer state as bit-identity oracles (rather than the
converged scalar flux), because the L7-trap is invisible in
single-iteration tests: the FIRST iteration's reflective-buffer state
is zero, so per-ordinate vs once-per-octant give the same answer; the
trap surfaces only on the SECOND iteration when the first iteration's
outgoing-face writes feed the second iteration's BC apply.  The case
also explicitly tags ``@pytest.mark.catches("ERR-003")``: ERR-003 is
the catalogued instance where reflective-BC ordering coupled with
ordinate batching produced a converged-but-wrong solution.

Bit-identity to the legacy implementation
-----------------------------------------

For LS-family quadratures (``LevelSymmetricSN``,
``ProductQuadrature``) whose ordinate ordering is octant-grouped in
lexicographic order, the Wave-2 implementation is **bit-identical**
to the legacy per-ordinate loop on every regression snapshot — the
existing
``tests/sn/regression/snapshots/2d_1g_LS4_dd_15x15.npz``,
``test_apply_2d_cartesian_bit_identical_to_legacy``, and
``test_unified_sweep_dispatch`` snapshots all pass with
``np.array_equal``.  The argument has three parts:

1. **BC apply equivalence.**  The boundary operator for octant
   :math:`\sigma` reads partner-octant rows of the persistent
   ``psi_x`` / ``psi_y`` buffer.  For LS, the lex order of
   ``quad.octants`` matches the legacy n-order at the
   partner-state granularity, so the same iteration's value is
   observed at the same point.  Per-octant BC apply produces the
   same ``psi_x`` / ``psi_y`` octant-row contents as :math:`N_\sigma`
   copies of the legacy per-ordinate apply.
2. **Per-cell sweep equivalence.**  Within an octant, per-ordinate
   cell sweeps are independent — different rows of ``psi_x`` /
   ``psi_y``, different rows of ``angular_flux``.  Batching is
   therefore bit-identical to any per-ordinate sequencing of the
   same set, modulo the per-cell DD operation order which is
   pinned (see :ref:`sweep-octant-dependency-graph` dispatch
   boundary above).
3. **Lebedev (octant ordering not lex).**  For Lebedev quadrature
   (case 6 in the C2.5 harness) the converged scalar flux matches
   the legacy code, but the iter-to-iter values differ on the
   inner iteration (different traversal order ⇒ different
   Gauss-Seidel updates).  Case 6 uses **vacuum BCs**, where the
   partner-state semantics don't matter; this is a deliberate
   choice in the harness — Lebedev with reflective BCs would
   require redesign of the bit-identity gate, but the converged
   answer is still verified at the snapshot level via
   ``test_unified_sweep_dispatch``.

This taxonomy follows the ``vv-principles`` skill's
**bit-identity vs principled-equivalence** discipline: bit-identity
where structurally trivial (LS-family, octant-grouped lex order);
principled equivalence (closed-form L1 anchor + MMS regression suite)
where bit-identity would require more work than the engineering value
returns.

.. _sweep-octant-architecture-cardinal-rule-2:

Architectural framing (Cardinal Rule 2)
---------------------------------------

Per the project memory note ``project_moc_structure.md`` and the
:ref:`cell-update-strategies` discussion,
:class:`SweepDependencyGraph` is **SN-specific by design**.  MoC
will define its own analog (per-ray traversal) — different DAG
shape, different mathematical structure (fiber bundles + solution
sheaves over characteristic curves rather than a topological sort
over a cell graph).  There is **no shared SweepGraph Protocol**
because there is no shared mathematical structure.  Cardinal Rule 2
(architecture) prefers **late unification** ("unify after two
instances" — see ``feedback_unify_after_two_instances`` in agent
memory) to premature abstraction; the sweep DAG lives in
:mod:`orpheus.sn` and stays there until a second mathematically-
similar consumer arrives, which by current understanding is **never**
for MoC and only conjectural for any other deterministic transport
solver.

By contrast, the **angular octant partition** primitive
(:meth:`~orpheus.numerics.measure.DiscreteMeasure.partition_by`) is
genuinely shared infrastructure — see the cross-method consumer
table at :doc:`/theory/foundations/discrete_measures` (octant partition consumed by SN
2-D, MoC track-bundle direction grouping, MC boundary-current
hemisphere scoring, future SN boundary realiser).  The split is
**measure-level primitives are shared, sweep-level orchestration
is SN-specific**.

Performance
-----------

The Wave-2 plan target for Issue #4 closure was 3–10× speedup on
the 421-group benchmark (the canonical ``test_profile_421g``
smoking-gun probe).  The shipped speedups:

.. list-table::
   :header-rows: 1
   :widths: 35 20 15 30

   * - Configuration
     - Speedup
     - Target?
     - Comment
   * - 421-group LS4 ``31×31``
     - 1.7×
     - Below
     - numpy-dispatch-overhead-dominated regime
   * - 2-group LS4 ``31×31``
     - 2.78×
     - At lower bound
     - per-octant batching wins more for small ``ng``
   * - 1-group LS4 ``15×15`` (regression snapshot)
     - bit-identical
     - n/a
     - regression contract preserved

The headline 421-group speedup is below the 3-10× target.  The
honest analysis: the Wave-2 implementation eliminates the
:math:`N`-fold ordinate loop overhead but the per-octant per-level
kernel calls still number :math:`O(\text{levels} \times
\text{octants}) \approx (n_x + n_y - 1) \times 4 \approx 88` per
sweep on a ``31 × 31`` mesh, each carrying its own numpy dispatch
cost.  At 421 groups, the per-call work scales linearly so the
ratio of useful work to dispatch overhead remains modest.  The
**follow-up direction** noted at Wave 2 was to carry full-:math:`N`
buffers plus an ``octant_indices`` field so the kernel calls become
level-only (~ 60 calls / sweep) rather than ``levels × octants``
(~ 240 calls / sweep), eliminating the per-octant copy round-trip.
The subsequent Phase 5 / S6.4 work took a different route to the
same end: the rolling-frontier window
(:meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_windowed`)
holds the interior cochain on a contiguous :math:`(d{-}1)`-frontier
slab, turning the per-level gather into a basic-slice zero-copy view
(a measured ``~0.77×`` contiguity speedup AND a ``~3×`` peak-memory
win at d=2) — see :ref:`wavefront-flux-cochain`.

Closing the smoking gun by construction is itself a load-bearing
result: the legacy ``for n in range(N)`` is gone, the metric
(angular-flux tensor; see :ref:`theory-sn-index-convention` for the
canonical ``(N, ng, nx, ny)`` storage) now knows its iterative
structure, and any future numpy-dispatch-cost reduction benefits
all closures uniformly through the strategy contract.

Verification
------------

The Wave-2 verification chain (per the ``algebra-of-record`` skill
discipline):

* **L0 unit tests** — on the primitives:

  - ``tests/sn/test_octants_property.py`` (60 tests across 8
    quadrature factories) — disjoint union, weight conservation,
    sign-signature correctness, pure-axis ordinates labelled
    ``sign=0``.
  - ``tests/sn/test_cell_kernel_batch.py`` (S6.4(e) successor of
    ``test_cell_update_batch.py``) — term-level L0 on the storage-free
    kernel pair (``cell_kernel_batch`` / ``residual_kernel_batch``):
    bit-identity against per-cell :meth:`update` on a
    single-cell-per-batch reduction; standalone tests against
    analytical DD recurrence on a 1×3 strip; 4-octant bit-identity vs
    the per-ordinate Python loop; plus a ``sha256`` source-of-record
    pin on the two kernel bodies (the explicit left-fold order is
    bit-identity-load-bearing).
  - ``tests/sn/test_sweep_graph.py`` — the §15A.2 invariant set above;
    anti-diagonal cell coverage; topo-order acyclicity per octant sign;
    BC face conventions; and the ``walk_full`` / ``walk_windowed`` ×
    level-operation walks (with ``window ≡ full`` bit-identity oracles).
  - ``tests/sn/primitives/test_dag_ownership.py`` (S6.4(c) successor of ``test_snmesh_sweep_graphs.py``) — graph
    contents agree with hand-derived schedule on a 3×3 mesh; dict
    keys equal ``quad.octants`` labels; cache invalidates when mesh
    changes.

* **L1 closed-form anchor + L7-trap detector** — the C2.5 TESTS-
  FIRST harness ``tests/sn/test_2d_octant_sweep_equivalence.py``
  (7/7 pass), tagged ``@pytest.mark.l1`` and
  ``@pytest.mark.catches("ERR-003")``.  Includes:

  - **case 3 (L7 trap)** — mixed BC + 2G heterogeneous +
    ``n_sweeps=2``, the primary L7-trap detector.
  - **case 7 (closed-form)** — 1G homogeneous reflective with
    :math:`k_\infty = \nu\Sigma_f / \Sigma_a`, the structural-
    independence anchor.
  - cases 1–6 covering BC mixes, ordinate batching corners, and
    Lebedev (vacuum-BC variant).

* **L2 regression** — existing ``tests/sn/verification/mms/test_mms_2d.py``,
  ``test_discrete_ordinates_2d.py``, ``test_streaming_operator.py``,
  ``test_streaming_operator_decomposition.py``,
  ``test_unified_sweep_dispatch.py``, ``tests/sn/regression/``: 56/56
  pass, 6 slow-marked skipped.

The verification chain is the canonical
**L0 (primitive invariants) → L1 (closed-form anchor + bug catcher)
→ L2 (integration regression)** ladder from the ``vv-principles``
skill.

References and pointers
-----------------------

* Grand Report v3 §15A.2 (lines 2137–2171) — the "upwind trace
  complex / causal transport DAG / direction sweep ordering"
  primitive description with the ``assert_*`` invariant set this
  module's tests pin.  Plan file at
  ``.claude/plans/neutron_transport_grand_report_v3.md``.
* Wave 2 plan at ``.claude/plans/transient-giggling-cake.md`` (C2.1
  through C2.8) — the architectural primitives plan, sequencing,
  verification-first harness design.
* Wave 0 :meth:`~orpheus.numerics.measure.DiscreteMeasure.partition_by`
  primitive — the measure-level partition that the SN ``octants``
  property delegates to.  See :doc:`/theory/foundations/discrete_measures`.
* Wave 1 :math:`R \circ \Lambda \circ M` Galerkin scattering
  composition — the parallel "metric knows its iterative structure"
  refactor for the scattering source build.  See
  :ref:`sn-scattering-fission-operators`.
* :class:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase` — the
  strategy ABC carrying the per-cell :meth:`update` reference contract
  and the storage-free batched kernel pair
  :meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.cell_kernel_batch`
  / :meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.residual_kernel_batch`
  (S6.4(e); was the ``SweepCellSlice``-packeted ``update_batch`` /
  ``residual_batch``).
* :class:`~orpheus.transport.spatial.diamond.DiamondDifference` — the only
  shipping closure that overrides the kernel pair; the reference for
  the bit-identity contract (pure cell algebra — the ONLY
  direction-aware math in the SN spatial stack since S6.4(e) lifted
  storage to the walk layer).
* :mod:`orpheus.sn.loss_representation.sweep_graph` — the two storage walks
  (:meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_full`,
  :meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_windowed`)
  and the ``_CellSolve`` / ``_CellResidual`` level operations.
* C2.5 TESTS-FIRST harness:
  ``tests/sn/test_2d_octant_sweep_equivalence.py``.


Choosing the schedule: the representation layer
===============================================

Everything above is ONE schedule — the DAG wavefront — for one
operator. The full selection story is the capstone architecture page
:doc:`loss_representation` (the representation layer of
:math:`(L+C)`); this chapter states only what the multi-D walk adds to
it:

* :math:`(L+C)` is one lower-triangular object with two actions
  (``solve`` = forward substitution, ``apply`` = the row action); a
  :class:`~orpheus.sn.loss_representation.LossRepresentation` is a
  **schedule** for traversing that triangular structure, not a
  different operator.
* Four schedules ship:
  :class:`~orpheus.sn.loss_representation.CumprodScan` (the 1-D
  parallel-prefix scan of :doc:`slab_one_group`),
  :class:`~orpheus.sn.loss_representation.ScanMarch` (scan the
  contiguous axis, march the others — the **multi-D Cartesian
  production default**; the Fork-B2 rationale lives on the capstone),
  :class:`~orpheus.sn.loss_representation.MovingFrontierWindow` (the
  anti-diagonal wavefront above with a rolling
  :math:`(d{-}1)`-frontier — a selectable peer), and
  :class:`~orpheus.sn.loss_representation.FullFieldWavefront` (the
  same DAG schedule retaining the whole interior cochain — the
  verification oracle, explicit-select only).
* Selection is a single source of truth:
  :func:`~orpheus.sn.loss_representation.default_for` returns the
  first compatible entry of the ordered registry, and an illegal
  ``(representation, mesh)`` pairing is unrepresentable — the
  constructor re-checks
  :meth:`~orpheus.sn.loss_representation.LossRepresentation.supports`
  and raises.
* Whatever the schedule, **one d-generic walk frame serves sweep AND
  matvec** (the L21 invariant), forked only by a kernel object and an
  emit policy — never a boolean. The dependency graph and the kernel
  pair documented above are the shared substrate every multi-D
  schedule rides.

The historical Wave-D ``transport_sweep`` consolidation that first
unified the sweep paths — since retired in favour of this
representation polymorphism — is preserved as origin history in
:doc:`index` (the superseded "Unified sweep dispatch" section).


.. _ld-ubld-multidim:

Multi-dimensional LD: the tensor-product bilinear (UBLD) cell system
=====================================================================

The multi-dimensional analog of Linear Discontinuous on a **Cartesian**
cell is **NOT** the simplex-P1 :math:`\{1, x, y\}` object
(:math:`1+d` moments).  Adams (2001) proved simplex-LD *fails* the thick
diffusion limit on quadrilaterals, while the **bilinear / trilinear
DG-P1** (UBLD) — basis :math:`\{1, x, y, xy\}` (:math:`2^d` moments) —
*passes*.  The :math:`xy` cross moment is diffusion-limit-load-bearing.

The :math:`d`-generic per-cell Galerkin system is assembled as Kronecker
products of the verified 1-D LD factor operators (the streaming
:math:`\Omega\cdot\nabla = \sum_a \mu_a \partial_a` is a sum over axes;
the tensor-product basis separates):

.. math::
   :label: ld-ubld-cell-system

   A\,\vec\psi = \vec R, \qquad
   A = G + F_{\rm out} + \Sigma_t M, \qquad
   \vec R = M\,\vec S + F_{\rm in}\,\psi_{\rm in}^{\rm traces},

a :math:`2^d \times 2^d` dense non-symmetric solve, with
:math:`M = M_1 \otimes \cdots \otimes M_d` (mass),
:math:`G = \sum_a \mu_a\,(M_1 \otimes \cdots \otimes G_{1d} \otimes
\cdots \otimes M_d)` (streaming: gradient on the active axis, mass on the
transverse axes), and :math:`F_{\rm out}` likewise from the per-axis
downstream-face trace.

.. math::
   :label: ld-ubld-d1-reduction

   A\big|_{d=1} =
   \begin{bmatrix} \Sigma_t h + |\mu| & |\mu| \\
                   -|\mu| & \Sigma_t \theta h + |\mu| \end{bmatrix}

The :math:`d=1` reduction (Kronecker-with-one-factor identity) recovers
the production slab 2×2 :eq:`dd-cartesian-1d`-sibling exactly; the
:math:`xy` coupling falls out of the algebra for :math:`d \ge 2`.

.. math::
   :label: ld-ubld-exact-on-bilinear

   \psi(x,y) = a + bx + cy + dxy
   \;\Longrightarrow\;
   \vec\psi_{\rm solved} = \vec\psi_{\rm exact-projections}

The UBLD is **exact on any bilinear flux** (the multi-D analog of the
1-D "exact on linear-in-x" oracle), the :math:`xy` cross moment
exercised — the structurally-independent correctness gate for the
:math:`d \ge 2` closure.

The Branch-1 algebra-of-record (the UBLD weak form)
------------------------------------------------------

The canonical symbolic reference for the :math:`d`-generic UBLD system is the
SymPy module :mod:`orpheus.derivations.discrete.sn.ld_ubld` (the
algebra-of-record, State 1A closed-form): the Kronecker assembler
``assemble_ubld`` plus five ``derive_*`` verification functions, each proven by
``sympy.simplify(diff) == 0``.  It is the discrete-SN sibling of
``orpheus.derivations.discrete.sn.balance`` — a *symbolic discretization the
production solver must satisfy*, NOT a continuous reference.

The per-cell system descends from the Galerkin weak form of the within-group
streaming–collision operator (Maginot, Ragusa & Morel 2016, "Non-negative
Methods for Bilinear Discontinuous Differencing of the :math:`S_N` Equations on
Quadrilaterals", NSE 185(1):17–42, Eqs. 1–12).  Multiplying the transport
equation :math:`\Omega\cdot\nabla\psi + \Sigma_t\psi = S` by each basis function
:math:`B_i` and integrating over the cell :math:`K`, then integrating the
streaming term by parts (MRM-2016 Eq. 6),

.. math::
   :label: ld-ubld-weak-form

   \underbrace{(\Omega\cdot)\!\oint_{\partial K} \hat n\,B_i\,\psi\,d\ell}_{\text{surface (upwind)}}
   \;-\; \int_K \psi\,(\Omega\cdot\nabla B_i)\,dV
   \;+\; \Sigma_t\!\int_K B_i\,\psi\,dV
   \;=\; \int_K B_i\,S\,dV,

gives three operators per cell — the **mass** :math:`M_{ij} = \int_K B_i B_j`
(the collision term), the **gradient/stiffness** :math:`G_{ij} = \int_K B_i\,(\Omega\cdot\nabla B_j)`
(the volumetric streaming term, coupling all :math:`2^d` moments), and the
**surface** matrix split per face into an OUTFLOW block (:math:`\Omega\cdot\hat n > 0`,
implicit, the cell's own unknowns) and an INFLOW block (:math:`\Omega\cdot\hat n < 0`,
**upwind**: the incoming face value is the upstream neighbour's outflow trace,
moved to the RHS).  Assembling gives exactly the dense system
:eq:`ld-ubld-cell-system`, :math:`A = G + F_{\rm out} + \Sigma_t M`,
:math:`\vec R = M\vec S + F_{\rm in}\,\psi_{\rm in}^{\rm traces}`.

Why bilinear, not simplex-P1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The naïve multi-D analog of 1-D LD — "cell average plus one slope per axis",
the simplex-P1 basis :math:`\{1, x, y\}` of :math:`1+d` moments — is the
**wrong object on a Cartesian cell**.  Adams (2001, "Discontinuous Finite
Element Transport Solutions in Thick Diffusive Problems", NSE 137(3):298–333)
proved that simplex-LD *fails* the thick diffusion limit on quadrilaterals,
while the **bilinear / trilinear DG-P1** (UBLD) — the tensor-product basis
:math:`\{1, x, y, xy\}` of :math:`2^d` moments — *passes* it.  The reason is the
:math:`xy` **cross moment**: it is exactly what the simplex basis lacks, and it
is the term the leading-order asymptotic diffusion balance needs (Börgers,
Larsen & Adams 1992, "The asymptotic diffusion limit of a linear discontinuous
discretization of a two-dimensional linear transport equation", JCP
98(2):285–300, give the 2-D rectangular analysis explicitly).  The simplex-P1
basis *does* preserve the limit on a genuine simplex (triangle/tetra) mesh
(Wareing, McGhee, Morel & Pautz 2001, NSE 138(3):256–268) — but that is a
different cell topology, not a quadrilateral.  ORPHEUS builds Cartesian cells,
so the :math:`2^d` tensor-product object is the diffusion-limit-consistent
choice; the choice is **load-bearing**, not a convenience (see
:ref:`ld-ubld-scattering-moment-lift` for the companion half of the same
asymptotic argument).

The Kronecker single-source build
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The three matrices are NOT hand-transcribed entry-by-entry (that would be a
:math:`4\times4` / :math:`8\times8` transcription waiting for a sign error).
They are assembled as **Kronecker products of the verified 1-D LD factor
operators** in the Legendre moment basis :math:`\{1, P_1\}` on width :math:`h`:

.. math::
   :label: ld-ubld-kronecker-factors

   M_{1d} = \mathrm{diag}(h,\ \theta h),
   \qquad
   G_{1d} = |\mu|\begin{bmatrix} 0 & 0 \\ -2 & 0 \end{bmatrix},
   \qquad
   F_{\rm out}^{1d} = |\mu|\begin{bmatrix} 1 & 1 \\ 1 & 1 \end{bmatrix},

with the streaming :math:`\Omega\cdot\nabla = \sum_a \mu_a\,\partial_a` a sum
over axes (the tensor-product basis separates), so

.. math::
   :label: ld-ubld-kronecker-assembly

   M = M_1 \otimes \cdots \otimes M_d,
   \qquad
   G = \sum_a \mu_a\,(M_1 \otimes \cdots \otimes G_{1d}^{(a)} \otimes \cdots \otimes M_d),

i.e. the gradient acts on the active axis and **the mass on every transverse
axis** (the volume-integral factorization — this is the load-bearing build
choice).  The :math:`F_{\rm out}` surface matrix is assembled likewise; the
inflow is a :math:`B(-1) = [1, -1]` test-weighting on the active axis (mass on
transverse axes) times :math:`|\mu_{\rm axis}|`.  The diagonal mass weights are
then the Kronecker product of the per-axis diagonals — a power of
:math:`\theta = \tfrac13` equal to the **number of active (slope) axes** of each
moment:

.. math::
   :label: ld-ubld-mass-weights

   M_{ii} = \theta^{|i|},
   \qquad
   |i| = \#\{a : o_a = 1\}
   \;\Longrightarrow\;
   \begin{cases}
     1        & \bar\psi \quad (\text{no slope axis}) \\
     \theta   & \hat\psi_x,\ \hat\psi_y \quad (\text{one slope axis}) \\
     \theta^2 & \hat\psi_{xy} \quad (\text{two slope axes})
   \end{cases}

so the 2-D weights are :math:`(1, \theta, \theta, \theta^2)` and the 3-D
:math:`xyz` cross moment carries :math:`\theta^3`.  These weights re-appear in
the matvec mass-normalization (:eq:`ld-ubld-unified-moment-residual`) — they are
the SAME diagonal Legendre mass.  The :math:`d=1` case is a Kronecker product
with a single factor (an identity), so it reduces EXACTLY to the production
slab 2×2 :eq:`ld-ubld-d1-reduction`; the :math:`xy` coupling *emerges* from the
algebra for :math:`d \ge 2` — no entry is hand-written.

The two oracles
~~~~~~~~~~~~~~~

The module proves the construction with two structurally distinct oracles
(both ``sympy.simplify(diff) == 0``):

* **Oracle (i) — the :math:`d=1` reduction.**
  ``derive_d1_reduction_to_production`` shows the assembled :math:`d=1` system
  equals the production
  :mod:`orpheus.transport.spatial.linear_discontinuous` 2×2 entry-for-entry, with the
  Schur complement :math:`S` and the effective slope denominator
  :math:`D_2' = \Sigma_t h\theta + |\mu|` recovered as the production closed
  forms.  Two further reductions
  (``derive_d1_kernel_view_equals`` / ``derive_d1_scan_view_equals``) prove the
  same :math:`d=1` reduces to BOTH the ÷V DAG kernel ``_kernel_terms`` and the
  ×V scan ``affine_scan_coefficients`` views — the "single-source the math"
  proof that Branch 2's three production views are the SAME algebra
  (:eq:`ld-ubld-rule-of-three-collapse`).

* **Oracle (ii) — exact-on-bilinear at :math:`d=2`.**
  ``derive_d2_exact_on_bilinear`` feeds an exactly-bilinear flux
  :math:`\psi = a + bx + cy + dxy` through the DG-exact upstream face moments and
  the projected source moments, and shows the 4 solved moments equal the exact
  projections (:eq:`ld-ubld-exact-on-bilinear`).  The :math:`xy` cross moment is
  genuinely exercised (:math:`d \ne 0` symbolically) — this is the multi-D
  analog of the 1-D "exact on linear-in-x" oracle and the structurally
  independent correctness gate for the :math:`d \ge 2` closure.

The foundation gate is :mod:`tests.transport.spatial.test_ld_ubld_symbolic` (6
``@pytest.mark.foundation`` tests, one per ``derive_*`` plus an anchor to the
live production ``LinearDiscontinuous.update``); the literature contract is
recorded in
``.claude/agent-memory/literature-researcher/multi_d_ld_closure.md`` (MRM-2016
Eqs. 1–12; the Adams-2001 thick-diffusion verdict; BLA-1992); the closeout is
``.claude/agent-memory/method-implementer/issue_240_d5b_s1_ld_ubld_branch1_closeout.md``.

.. admonition:: ERR-060 — the oracle that earned its keep
   :class: tip

   The first draft of ``assemble_inflow_axis`` dropped the :math:`|\mu_{\rm axis}|`
   streaming factor on the inflow RHS (failure Mode 3, a missing factor).  The
   bug was INVISIBLE to all three :math:`d=1` oracles — the :math:`d=1` RHS is
   built inline, never routed through the multi-axis inflow assembler — and was
   caught by Oracle (ii), the :math:`d=2` exact-on-bilinear gate, which is the
   first consumer of ``assemble_inflow_axis``.  Mutation-verified
   :math:`-O`-safe: re-dropping the factor turns the :math:`d=2` test red while
   the :math:`d=1` tests stay green (proving they are blind to the bug class).
   This is the algebra-of-record discipline working as designed — the bug never
   reached production.  (The ERR-060 marker belongs on the *exact-on-bilinear*
   gates, NOT on the cell-matrix ``A == A`` pin, which checks
   ``assemble_ubld``'s A/M/G/:math:`F_{\rm out}` and is structurally blind to the
   dropped inflow factor — see ``error_catalog.md`` ERR-060.)

.. _ld-ubld-branch2-primitive:

The Branch-2 production primitive + the single-sourced d=1 fast path
------------------------------------------------------------------------

The numpy production counterpart of the symbolic algebra-of-record above
is :mod:`orpheus.transport.spatial._ubld`, in two layers that share ONE source
of truth.  Layer 1 is the :math:`d`-generic dense primitive
(``assemble_ubld`` / ``per_cell_solve``): a batched-over-cells Kronecker
build of the :math:`2^d \times 2^d` system :eq:`ld-ubld-cell-system`,
solved with a batched :func:`numpy.linalg.solve`.  It is the CANONICAL
:math:`d`-generic source for both :math:`d=1` (today) and :math:`d \ge 2`
(S2 wires the bilinear cell-batch kernel onto it); in production
:math:`d=1` does **not** route through this dense solve (that would be
the per-cell-solve performance regression).

Layer 2 is the shared :math:`d=1` closed form ``d1_closed_form`` — the
analytic Schur complement of the primitive's :math:`d=1` 2×2, VECTORIZED
over the cell / ordinate / group stack (no dense solve), so the
production :math:`d=1` path stays on the fast path.  The entire closure
rides two SCALE-FREE invariants:

.. math::
   :label: ld-ubld-scale-free-invariants

   k = \frac{g/\theta}{g/\theta + \Sigma_t}, \qquad
   w = \frac{1}{1 + k}, \qquad g = \frac{|\mu|\,A_{\rm down}}{V},

with ``w`` the cell-average blend weight
(:math:`\bar\psi = (1-w)\psi_{\rm in} + w\,\psi_{\rm out}`).  Every
production view's coefficients are an algebraic function of
:math:`(g, \Sigma_t, k, w)` times a power of the cell volume :math:`V`
(the ×V vs ÷V choice applied at the call site).

.. math::
   :label: ld-ubld-rule-of-three-collapse

   \texttt{\_schur\_terms}\;(\times V), \quad
   \texttt{\_kernel\_terms}\;(\div V), \quad
   \texttt{affine\_scan\_coefficients}\;(\times V\ \text{scan})
   \;\longleftarrow\; \texttt{d1\_closed\_form}

The three production 1-D views in
:mod:`orpheus.transport.spatial.linear_discontinuous` — the ×V per-cell Schur
(``_schur_terms``), the ÷V DAG kernel (``_kernel_terms``), and the ×V
scan (``affine_scan_coefficients``) — now ALL derive their coefficients
from ``d1_closed_form``, applying their ×V / ÷V / ×V-scan scaling at the
call site.  The LD 2×2 algebra (the Rule-of-Three) collapses to ONE
place, proven ``==`` the dense primitive's :math:`d=1` reduction
(symbolically by the Branch-1 oracles, numerically end-to-end by the
Branch-2 gate).

The numpy production counterpart descends from the SAME specialized SymPy
ancestor as the Branch-1 algebra-of-record above; only the evaluation strategy
differs (Branch 1 closes the algebra symbolically, Branch 2 evaluates it on
arrays).  The discipline is **construct general, select narrow, specialize only
on measured cost**:

* **Construct general — the dense primitive.**  ``assemble_ubld`` /
  ``per_cell_solve`` build and solve the :math:`2^d \times 2^d` system
  :eq:`ld-ubld-cell-system` for every :math:`d`, batched over the cell /
  ordinate / group stack with a single :func:`numpy.linalg.solve`.  This is the
  canonical :math:`d`-generic source — :math:`d=1` (today), :math:`d \ge 2`
  (S2 wires the bilinear cell batch onto it), :math:`d = 3` (trilinear).

* **Select narrow — the :math:`d=1` closed form.**  ``d1_closed_form`` /
  :class:`~orpheus.transport.spatial._ubld.D1ClosedForm` is the analytic Schur
  complement of the primitive's :math:`d=1` 2×2, VECTORIZED over the stack with
  no dense solve.  Both scale-free invariants of :eq:`ld-ubld-scale-free-invariants`
  drive it.

* **Specialize on measured cost.**  In production :math:`d=1` does **not**
  route through the dense solve — that would be the per-cell-solve performance
  regression (the L16 constraint).  The closed form keeps the production
  :math:`d=1` sweep on the vectorized fast path
  (:class:`~orpheus.sn.loss_representation.CumprodScan` rides the ×V scan view's
  :math:`(a, \mathrm{inverse\_denom}, w)`; the DAG kernel rides the ÷V arrays).

The Rule-of-Three collapse
~~~~~~~~~~~~~~~~~~~~~~~~~~

Before the carve, the LD 2×2 algebra was transcribed in three production views
that had drifted into three independent copies.  All three now derive their
coefficients from the single ``d1_closed_form`` source
(:eq:`ld-ubld-rule-of-three-collapse`), applying only their scaling at the call
site — the Cardinal-Rule-2 / `coding-elegance` Pattern-2 single-source collapse:

.. list-table:: The three production 1-D views — one algebra, three scalings
   :header-rows: 1
   :widths: 30 16 54

   * - Production view (in :mod:`orpheus.transport.spatial.linear_discontinuous`)
     - Scaling
     - Consumer
   * - ``_schur_terms``
     - :math:`\times V`
     - the per-cell Schur (the matvec / ``update`` / ``residual`` path)
   * - ``_kernel_terms``
     - :math:`\div V`
     - the scale-free DAG wavefront kernel (the :math:`d \ge 2` arm rides the
       ÷V system, :eq:`ld-ubld-divv-scale-free-kernel`)
   * - ``affine_scan_coefficients``
     - :math:`\times V` scan
     - the Blelloch parallel-prefix scan (the production :math:`d=1` sweep)

The ×V / ÷V / ×V-scan choice is the volume scaling applied to the same
coefficients: dividing the Galerkin balance by the cell volume :math:`V` leaves
a scale-free system in the per-axis streaming :math:`g_a = |\mu_a| A_{\rm down}/V`
and :math:`\Sigma_t` alone (the form the :math:`d \ge 2` kernel consumes — fed
unit widths and :math:`\mathrm{mus} = (g_0, \ldots)`, it reduces EXACTLY to
``d1_closed_form``); multiplying restores the volume-weighted per-cell Schur
(:math:`D_2' = \theta V\,d_2`, :math:`S_{\times V} = V\cdot\mathrm{eff\_denom}`).
Each view is proven ``==`` the dense primitive's :math:`d=1` reduction —
symbolically by the Branch-1 oracles above (``derive_d1_kernel_view_equals`` /
``derive_d1_scan_view_equals``), numerically end-to-end by the Branch-2 gate.

.. note::

   **A principled ~1-ULP re-baseline, not a bit-identity break.**  Routing the
   three LD views through the shared helper changes the floating-point
   *reduction tree* relative to the legacy inline associations: the helper
   computes :math:`(g, \Sigma_t, k, w)` once and forms each coefficient as an
   algebraic function of them, whereas the old inline code interleaved the
   multiplies and adds differently.  In exact arithmetic the values are
   identical; in IEEE-754 they differ at ~1 ULP because addition is not
   associative.  This satisfies all three `vv-principles` criteria for
   accepting a non-bit-exact change: every intermediate
   (:math:`g`, :math:`k`, :math:`w`) is a *named, inspectable* physics
   quantity; the value is verified against the structurally-independent
   Branch-1 symbolic oracle; and the drift is FP-non-associativity bounded by
   the reduction depth.  The LD gates carry ``rtol = 1e-12`` (far above the
   ULP-scale drift); the DD-only strict gate remains the **bit-identical
   negative control** (DD never reaches the LD helper — its :math:`w=\tfrac12`
   reconstruction is the exact power-of-two doubling that commutes with
   round-to-nearest).

The verification is :mod:`tests.transport.spatial.test_ld_ubld_primitive` (10 tests):
the numpy primitive :math:`==` the SymPy oracle at :math:`d=1` (matrices +
moments) and exact-on-bilinear at :math:`d=2`; the shared closed form
:math:`==` the dense :math:`d=1` solve in all three views; and the LIVE
production scheme (``update`` / ``cell_kernel_batch`` /
``affine_scan_coefficients``) :math:`==` the dense primitive (the link proof).
The closeout is
``.claude/agent-memory/method-implementer/issue_240_d5b_s1_ld_ubld_branch2_closeout.md``.
The :math:`d \ge 2` hand-off (the bilinear cell-batch kernel + the
:math:`2^{d-1}`-moment face cochain wiring onto this primitive) is S2, the next
subsection.

.. _ld-ubld-d2-wavefront-wiring:

Wiring the d≥2 UBLD kernel onto the DAG wavefront (S2)
--------------------------------------------------------

Sub-step **D5b-S2** closes the :math:`d = 1`-only kernel raise so
Linear-Discontinuous runs in :math:`d \ge 2` on the DAG wavefront,
consuming the verified dense primitive.  Three contract widenings, all
GATED on a single scheme trait so Diamond Difference / Step stay
byte-identical:

.. math::
   :label: ld-ubld-n-spatial-moments

   \text{per-cell} = (\text{per\_axis})^{d}, \qquad
   \text{per-face} = (\text{per\_axis})^{d-1}, \qquad
   \text{per\_axis} =
   \begin{cases} 1 & \text{DD / Step} \\ 2 & \text{LD} \end{cases}

The class-level trait
:attr:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.spatial_basis_per_axis`
(the 1-D moment-basis size) indexes the whole contract via the
tensor-product structure: the per-cell unknown is
:math:`(\text{per\_axis})^d` (LD-2D: 4) and each downstream face carries
:math:`(\text{per\_axis})^{d-1}` transverse moments (LD-2D: 2).  The
boolean ``per_axis > 1`` gates the multi-moment face-cochain trailing
axis and the moment-reducing emit; DD/Step at ``per_axis == 1`` keep the
rank-:math:`r` scalar face and rank-3 ``psi_avg`` EXACTLY.

.. math::
   :label: ld-ubld-divv-scale-free-kernel

   A_{\div V}\,\vec\psi = M_{\div V}\,\vec S + \sum_a F_{\rm in}^{(a)},
   \qquad
   \psi_{\rm out}^{(a)}[t] = \psi[o_a{=}0,\,t] + \psi[o_a{=}1,\,t]

The :math:`d \ge 2` arm rides the **scale-free ÷V** form of the dense
system: dividing the Galerkin balance by the cell volume leaves a system
depending only on the per-axis ÷V streaming :math:`g_a = |\mu_a|/\Delta_a`
(the ``s_axes`` the kernel already receives) and :math:`\Sigma_t` — so the
dense assembler is fed unit widths and ``mus = (g_0, \ldots)``, reducing
EXACTLY to ``d1_closed_form`` at :math:`d=1`.  The :math:`d` downstream
faces are the trace of the tensor-Legendre solution at the downstream node
(:math:`P_0(+1) = P_1(+1) = 1` sums the :math:`o_a{=}0` and
:math:`o_a{=}1` blocks), in the :math:`2^{d-1}` transverse-Kronecker order
the next cell's upwind inflow consumes (out-of-cell == in-of-next-cell —
the closure consistency the matvec twin verifies).

The scale-free ÷V system fed to the dense primitive
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :math:`d \ge 2` arm rides the **scale-free ÷V** form of the dense system
:eq:`ld-ubld-divv-scale-free-kernel`.  Dividing the Galerkin balance by the
cell volume leaves a system depending only on the per-axis ÷V streaming
:math:`g_a = |\mu_a|/\Delta_a` (the ``s_axes`` the kernel already receives) and
:math:`\Sigma_t`.  In code this means the dense assembler ``_ubld_system`` /
``per_cell_solve`` is fed **unit widths** (``hs = [1, …]``) and
``mus = (g_0, …, g_{d-1})`` — so at :math:`d = 1` it reduces EXACTLY to
``d1_closed_form`` (the ÷V view of the Rule-of-Three above), and at
:math:`d \ge 2` it is the same dense object the Branch-1 oracle proves
exact-on-bilinear.  The kernel dispatch lives in
:mod:`orpheus.transport.spatial.linear_discontinuous` (``cell_kernel_batch`` /
``residual_kernel_batch``): the :math:`d=1` closed-form fast path vs the
:math:`d \ge 2` dense ``_ubld_system`` / ``per_cell_solve``.

The downstream faces are the trace of the tensor-Legendre cell solution at the
downstream node: since :math:`P_0(+1) = P_1(+1) = 1`, the outgoing face on axis
:math:`a` sums the :math:`o_a = 0` and :math:`o_a = 1` blocks of the cell moment
vector, producing a :math:`2^{d-1}`-moment face object (average + transverse
slopes) in the transverse-Kronecker order the next cell's upwind inflow
consumes.  *Out-of-cell == in-of-next-cell* is the closure consistency the
matvec twin verifies.

The moment-ordering crosswalk
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The cell moment vector is the tensor (Kronecker) product of the per-axis 1-D
Legendre basis, ordered **x-outer / y-inner** so the all-:math:`P_0` cell
average is always slot 0 (the same convention the :ref:`spatial-moment-space`
factor surfaces, :eq:`spatial-moment-kronecker-order`).  The Kronecker layout in
2-D is :math:`[\bar\psi,\ \hat\psi_y,\ \hat\psi_x,\ \hat\psi_{xy}]` (indexing
:math:`[o_x, o_y]` with :math:`o_x` outer); each downstream face carries its
:math:`2^{d-1}` transverse moments in the matching per-axis order.  The
crosswalk between the cell-moment order, the per-face transverse order, and the
downstream-node trace reconstruction is the design record's load-bearing detail
(``.claude/plans/issue_240_d5b_s2_crosswalk.md``; recovery anchor
``.claude/plans/issue_240_phase2_step_d5b_ubld.md``).

The DD bit-identity backward-compat invariant
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All three contract widenings — the dense kernel arm, the multi-moment
face-cochain trailing axis (:mod:`orpheus.sn.loss_representation.sweep_graph` ``_MovingFrontier``;
the ``_CellSolve`` / ``_CellResidual`` moment-reducing emit), and the window
zero-pad (:mod:`orpheus.sn.loss_representation`
``FullFieldWavefront._octant_face_cochain``, the ``_inflow_to_moments`` pad) —
are GATED on the single scheme trait
:attr:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.spatial_basis_per_axis`
of :eq:`ld-ubld-n-spatial-moments`.  At ``per_axis == 1`` (DD / Step) the tail
is the EMPTY tuple (:eq:`spatial-moment-append-policy`), so the scalar face and
the rank-3 ``psi_avg`` emit are kept byte-identical — DD backward-compatibility
falls out of the ``face_moment_tail`` formula, NOT an ``if scheme is DD`` branch.
This is the negative control: the DD/Step gate stays at its exact pre-S2 golden.

S2 scope boundaries
~~~~~~~~~~~~~~~~~~~

S2 wires the :math:`d \ge 2` UBLD kernel so 2-D LD *runs* on the DAG wavefront,
but it is deliberately scoped to the **average-moment iterate** only.  Three
things remain owed to the later sub-steps, and naming them here is the honest
boundary:

* The full ``loss_action`` / Krylov 2-D LD needs the spatial-moment iterate
  :math:`\hat\phi` to travel between sweeps so the scattering slope source
  :math:`\Sigma_s\hat\phi` couples the slopes globally — **S3**
  (:ref:`ld-ubld-unified-moment-matvec`).  Without it the S2 closure is
  :math:`O(h^2)` but diffusion-limit-inconsistent (the flat-source signature).

* The non-vanishing domain-inflow moment trace (the ``AngularBoundaryFlux`` /
  ``mesh.angular_trace`` widening to a :math:`2^{d-1}` transverse face moment) — **S4**
  (and its honest-scope caveat, :ref:`ld-cartesian-2d`).

* The strengthened vv Mode-7 stress-ansatz MMS and the thick-diffusive
  tripwire — **S4** and **S3** respectively.

The verification is the kernel round-trip + matvec-twin face reconstruction
(:mod:`tests.transport.spatial.test_linear_discontinuous` ``TestLDKernel``), the
end-to-end two-paths FFW :math:`\equiv` MFW, the DD :math:`\ne` LD routing-flip,
and the :math:`O(h^2)` convergence smoke
(:mod:`tests.sn.verification.mms.test_mms_ld_2d`), plus the :math:`d=2`
numpy↔symbolic entry-wise ``A == A`` cell-assembly pin and the
``test_d2_exact_on_bilinear`` ERR-060 catcher
(:mod:`tests.transport.spatial.test_ld_ubld_primitive`).

.. _two-moment-axes:

Two kinds of "moment": angular vs spatial
-------------------------------------------

.. warning::

   The word **moment** denotes two ORTHOGONAL things in this solver, and
   the collision is the single most common source of confusion when reading
   the multi-dimensional Linear-Discontinuous (LD) code.  An **angular
   moment** reduces the *direction* dependence; a **spatial moment** resolves
   the *within-cell position* dependence.  They are independent tensor
   factors of the flux, NOT two names for the same axis.

The discrete-ordinates flux is a function of three independent kinds of
variable: direction :math:`\Omega`, position :math:`\vec r` (which the mesh
splits into a cell index plus a *within-cell* coordinate), and energy
group :math:`g`.  Each admits its own moment expansion, and the LD scheme
in :math:`d \ge 2` carries two of them simultaneously.  Distinguishing them
is the prerequisite for reading the next two subsections.

**Angular moment** :math:`\phi_\ell^m` — *how the flux varies with
direction.*  Projecting the per-ordinate angular flux
:math:`\psi(\Omega_n)` onto the real spherical harmonics
:math:`\{Y_\ell^m\}` collapses the :math:`N` discrete directions into the
:math:`(\ell, m)` harmonic coefficients,

.. math::
   :label: two-moment-angular

   \phi_\ell^m(\vec r, g)
   \;=\;
   \sum_{n=1}^{N} w_n\, Y_\ell^m(\Omega_n)\, \psi_n(\vec r, g),

the typed home of which is
:class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`
on the space
:math:`\mathrm{SphericalHarmonicSpace}(L) \otimes \mathrm{CellGroupSpace}`.
The angular moment is a **replacement representation** of the angular flux:
a calculation holds EITHER the per-ordinate field
:math:`\psi` of shape :math:`(N, ng, *\text{spatial})` OR its harmonic
moments :math:`\phi_\ell^m` of shape
:math:`(L{+}1, 2L{+}1, ng, *\text{spatial})`, bridged by the
spherical-harmonic :class:`~orpheus.numerics.frame.GalerkinFrame`'s two faces —
its ``analysis`` face (:math:`M`, the
:math:`\psi \to \phi` reduction :eq:`two-moment-angular`) and its
``reconstruction`` face
(:math:`R`, the :math:`\phi \to \psi` lift).  You never carry both;
windowing the 2-D Cartesian iterate as :math:`\phi_\ell^m` instead of
:math:`\psi` is the harmonic-moment-projection memory win
(the :math:`N \to (L{+}1)(2L{+}1)` collapse, :eq:`harmonic-moment-projection`).
The :math:`\ell = 0` moment IS the scalar flux exactly.

**Spatial moment** :math:`\hat\psi` — *how the flux varies in space inside
one cell.*  A finite-volume / Diamond-Difference closure carries a single
number per cell (the cell average :math:`\bar\psi`).  The
Linear-Discontinuous closure additionally resolves the **sub-cell slope**:
on a Cartesian cell it expands :math:`\psi` in the tensor product of a 1-D
Legendre basis :math:`\{1, P_1\}` per axis,

.. math::
   :label: two-moment-spatial

   \psi(x, y)\big|_{\rm cell}
   \;=\;
   \bar\psi\,
   + \hat\psi_x\, P_1(\xi_x)
   + \hat\psi_y\, P_1(\xi_y)
   + \hat\psi_{xy}\, P_1(\xi_x) P_1(\xi_y),
   \qquad \xi_a \in [-1, 1],

the four within-cell coefficients of the bilinear (UBLD) basis
:math:`\{1, x, y, xy\}` of :eq:`ld-ubld-cell-system`.  Unlike the angular
moment, the spatial moment is an **additional axis** that rides on whatever
angular representation is in play — it does NOT replace anything.  Its typed
home is the :ref:`spatial-moment-space`
(:class:`~orpheus.numerics.spaces.spatial_moment_space.SpatialMomentSpace`),
a tensor factor of length :math:`(\text{per\_axis})^d` that composes onto a
field's space alongside the cell/group/angular factors.

The two notions are summarised in the contrast table:

.. list-table:: The two "moment" axes — orthogonal tensor factors
   :header-rows: 1
   :widths: 18 26 28 28

   * - Property
     - Angular moment :math:`\phi_\ell^m`
     - Spatial moment :math:`\hat\psi`
     - Shared
   * - Resolves
     - direction :math:`\Omega`
     - within-cell position :math:`x`
     - —
   * - Basis
     - real spherical harmonics :math:`\{Y_\ell^m\}`
     - tensor-Legendre :math:`\{1, P_1\}` per axis
     - both orthogonal polynomial families
   * - Truncation knob
     - :math:`L` (max Legendre order)
     - ``per_axis`` (1 = DD/Step, 2 = LD)
     - —
   * - Count
     - :math:`(L{+}1)(2L{+}1)`
     - :math:`(\text{per\_axis})^d`
     - —
   * - Typed space
     - :class:`SphericalHarmonicSpace`
     - :class:`SpatialMomentSpace`
     - both :class:`FunctionSpace` factors
   * - Role on the flux
     - **replacement** (hold :math:`\psi` OR :math:`\phi_\ell^m`)
     - **additional** (rides on either)
     - both tensor factors of one field
   * - Set by
     - the angular Pℓ order requested
     - the spatial discretization scheme
     - —

Because they are independent factors, a fully-resolved LD-Pℓ angular flux
carries BOTH indices simultaneously — an angular index :math:`(\ell, m)`
and a spatial-moment index :math:`p`:

.. math::
   :label: two-moment-tensor-product

   \phi_\ell^{m, p}(\vec r_{\rm cell}, g),
   \qquad
   (\ell, m) \in \text{angular harmonics}, \quad
   p \in \{\bar{\,}, \hat x, \hat y, \widehat{xy}\}
   \ \text{(spatial moments)}.

The carrier space is the tensor product of the two moment spaces with the
cell/group space,

.. math::
   :label: two-moment-carrier-space

   \mathrm{SphericalHarmonicSpace}(L)
   \;\otimes\;
   \mathrm{CellGroupSpace}(ng, *\text{spatial})
   \;\otimes\;
   \mathrm{SpatialMomentSpace}(\text{per\_axis}, d),

so the stored ndarray gains a trailing :math:`(\text{per\_axis})^d`
spatial-moment axis after the :math:`(\ell, m, g, *\text{spatial})` prefix.
The orthogonality is what makes the architecture clean: the scattering
operator :math:`\Sigma_s` couples energy groups and (for anisotropic
scattering) angular moments, but is a **spectator** to the spatial-moment
axis — it scatters every spatial moment independently
(:eq:`ld-ubld-scattering-moment-lift`, next subsection).  Conversely the
spatial discretization (the sweep / cell solve) acts on the spatial moments
but is a spectator to the angular index.  Two operators, two axes, no
cross-talk except through the physics each is responsible for.

.. note::

   **Why an LD-P3 calculation needs both.**  Anisotropic scattering up to
   :math:`P_3` is an *angular*-resolution choice — it carries
   :math:`\phi_\ell^m` for :math:`\ell \le 3`.  The Linear-Discontinuous
   spatial closure is a *spatial*-resolution choice — it carries
   :math:`\hat\psi` for the within-cell slope.  An LD-P3 transport
   calculation makes both choices at once and so carries the full
   :math:`\phi_\ell^{m, p}` object of :eq:`two-moment-tensor-product`.
   Collapsing either axis to its average (P0 angular, or DD spatial)
   degrades a *different* physical fidelity: the angular collapse loses the
   flux's directional anisotropy; the spatial collapse loses the
   diffusion-limit accuracy that the :math:`xy` cross-moment provides
   (:eq:`ld-ubld-cell-system`, the load-bearing moment).

.. _ld-ubld-scattering-moment-lift:

The Σ_s ⊗ I spatial-moment scattering lift (S3-A, partial)
-------------------------------------------------------------

Sub-step **D5b-S3** completes the *physics* of the multi-dimensional UBLD
Linear-Discontinuous scheme.  Now that the two moment axes are clearly
distinguished (:ref:`two-moment-axes`), the completion is statable in one
sentence: the scattering source must scatter EVERY spatial moment, not just
the cell average.  Where S2 ships an O(h²) but diffusion-limit-INCONSISTENT
closure (it scatters only the spatial-AVERAGE moment — the slope rows of the
source are zero), S3 threads the canonical slope source so the converged
operator becomes the diffusion-limit-CONSISTENT one.

The load-bearing bridge is the scattering operator's
:math:`\Sigma_s \otimes I_{\rm spatial}` lift: :math:`\Sigma_s` carries no
spatial-moment index (it is an energy-group :math:`\to` energy-group matrix
per Legendre order), so it is applied to EVERY spatial moment of the scalar
flux INDEPENDENTLY,

.. math::
   :label: ld-ubld-scattering-moment-lift

   \bigl(\Sigma_s \otimes I_{\rm spatial}\bigr)\,
   (\bar\phi,\ \hat\phi_x,\ \hat\phi_y,\ \hat\phi_{xy})
   \;=\;
   (\Sigma_s\,\bar\phi,\ \Sigma_s\,\hat\phi_x,\
    \Sigma_s\,\hat\phi_y,\ \Sigma_s\,\hat\phi_{xy}),

so the spatial-moment axis is a SPECTATOR to the energy-group matmul,
exactly as the cell axis is.  In code this is the per-material group
contraction with a trailing ``...`` spectator
(``"fg,fc...->gc..."``): at the single-moment closures (Diamond
Difference / Step, ``per_axis == 1``) the trailing axis is ABSENT and the
``...`` matches nothing, so the lift is BYTE-IDENTICAL to the pre-S3
scattering (the negative-control bit-identity, verified rank-2-exact).
At ``per_axis == 1`` :math:`S_{\rm full} \equiv S_{\rm flat}`; only an LD
multi-moment closure activates the slope rows.

.. admonition:: Status — S3-A is PARTIAL
   :class: caution

   The :math:`\Sigma_s \otimes I_{\rm spatial}` lift documented here is the
   LANDED half of S3-A.  The :math:`\hat\phi` spatial-moment **iterate
   carrier** that FILLS the slope rows the lift now accepts is OWED (it was
   blocked on the typed-field space widening — the
   :ref:`spatial-moment-space` subsection, the prerequisite that was minted
   next).  The lift therefore scatters a slope source that, in the production
   path, is still zero (no field carries :math:`\hat\phi` yet); the converged
   fixed point does not change UNTIL the iterate carrier lands.  This page
   marks what is wired (the lift) versus what is owed (the iterate, the
   cell-emit accumulation, the source seams) so a future reader knows the
   S3-A wiring is mid-flight, not complete.  (**Since completed** — the
   iterate carrier, the cell-emit accumulation, and both source seams
   landed with the unified moment matvec and the :math:`d{=}1` moment
   scan: :ref:`ld-ubld-unified-moment-matvec` and
   :ref:`ld-ubld-moment-scan` below.  This admonition is preserved as
   the campaign-time boundary record.)

Physics-completion, not an iteration-only change
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The distinction matters for verification, so it is worth stating precisely.
Most changes to the iteration machinery — a Gauss-Seidel splitting, a σ\
:sub:`r`-removal, a synthetic accelerator (DSA), a preconditioner — MUST NOT
change the converged fixed point; they change only the *rate* at which the
iteration reaches it.  The correctness gate for such a change is
**FP-invariance**: the accelerated solve and the plain solve converge to the
same flux (`vv-principles` Mode 9).

The :math:`\Sigma_s \otimes I_{\rm spatial}` slope source is **not** that
kind of change.  S2 and S3 solve DIFFERENT operators:

.. math::
   :label: ld-ubld-s2-s3-operators

   \text{S2:}\quad (L + C - S_{\rm flat})\,\psi = Q_{\rm ext},
   \qquad
   S_{\rm flat} = \Sigma_s \otimes e_0 e_0^{\mathsf T},
   \\[4pt]
   \text{S3:}\quad (L + C - S_{\rm full})\,\psi = Q_{\rm ext},
   \qquad
   S_{\rm full} = \Sigma_s \otimes I_{\rm spatial}.

:math:`S_{\rm flat}` (the rank-1 projector :math:`e_0 e_0^{\mathsf T}` onto
the cell-average moment) scatters ONLY the spatial average — the slope rows
:math:`\hat\phi_x, \hat\phi_y, \hat\phi_{xy}` of the scattering source are
identically zero.  :math:`S_{\rm full}` (the identity on the spatial-moment
axis) scatters all of them.  The two operators have DIFFERENT spectra, hence
DIFFERENT fixed points.  The converged flux CHANGES — and that is the POINT:
the thick-diffusion-limit tripwire (the ``test_ld_thick_diffusive_limit``
xfail) flips xfail :math:`\to` pass *because* the limit becomes correct, not
because the iteration was accelerated.  S3 is therefore **NOT** verified
against the S2 fixed point; verifying it that way would be the Mode-9
mis-application (asserting FP-invariance of a change that legitimately moves
the FP).  The genuine Mode-9 invariant for S3 is the within-group analog:
source-iteration with a lagged moment iterate :math:`\equiv` direct / Krylov
solve of the **same** :math:`(L + C - S_{\rm full})` operator.

Why the slope rows are diffusion-limit-load-bearing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the optically-thick, scattering-dominated (:math:`c \to 1`,
:math:`\Sigma_t h \gg 1`) limit, the transport solution must collapse to the
diffusion solution.  Adams (2001) and Larsen, Morel & Miller (1987) showed
that a spatial discretization passes this asymptotic limit only if its
discrete diffusion limit is a valid (consistent and stable) diffusion
discretization.  For the bilinear UBLD cell the diffusion limit couples the
slope moments :math:`\hat\phi` — the leading-order asymptotic balance is a
relation between the cell-average and the slopes, not the average alone
(Border, Lewis & Adams 1992 give the 2-D asymptotic analysis explicitly).
If the *scattering source* feeds only the cell average (S2's
:math:`S_{\rm flat}`), the slope rows of the within-cell balance see no
scattering re-supply, the discrete diffusion limit is the WRONG diffusion
operator, and the thick-limit error does not vanish under refinement — the
xfail tripwire stays red.  Threading :math:`\Sigma_s \hat\phi` into the slope
rows (:math:`S_{\rm full}`) restores the correct discrete diffusion limit.
This is why the completion is *physics* (the converged answer becomes right
in a regime where it was wrong), not iteration bookkeeping.

.. note::

   This is the same asymptotic reasoning that selects the bilinear
   :math:`\{1, x, y, xy\}` basis over the simplex :math:`\{1, x, y\}` in the
   first place (:eq:`ld-ubld-cell-system` and the parent section): Adams
   (2001) proved the simplex-LD discrete diffusion limit is invalid on
   quadrilaterals.  The :math:`xy` cross-moment carries the limit; the
   :math:`\Sigma_s \otimes I_{\rm spatial}` lift makes sure that cross-moment
   (and the axis slopes) actually receive scattering.  The basis choice and
   the scattering lift are two halves of the *same* diffusion-limit argument.

The producer-side spectator-broadcast (Pattern 7)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The lift is implemented as a one-character change to the einsum subscripts of
the three scattering producers in
:class:`~orpheus.transport.mesh.material_xs_field.MaterialXSField`:

.. list-table:: The :math:`\Sigma_s \otimes I_{\rm spatial}` lift in code
   :header-rows: 1
   :widths: 34 30 36

   * - Producer
     - Subscript (pre-S3 :math:`\to` S3)
     - What it scatters
   * - :meth:`~orpheus.transport.mesh.material_xs_field.MaterialXSField.apply_p0_in_scatter`
     - ``"fg,fc->gc"`` :math:`\to` ``"fg,fc...->gc..."``
     - the P0 in-scatter :math:`\Sigma_{s,0}^{\mathsf T}\phi`
   * - :meth:`~orpheus.transport.mesh.material_xs_field.MaterialXSField.apply_n2n`
     - ``"fg,fc->gc"`` :math:`\to` ``"fg,fc...->gc..."``
     - the :math:`(n,2n)` source :math:`2\Sigma_{2n}^{\mathsf T}\phi`
   * - :meth:`~orpheus.transport.mesh.material_xs_field.MaterialXSField.apply_legendre_scattering_moments`
     - ``"mfc,fg->mgc"`` :math:`\to` ``"mfc...,fg->mgc..."``
     - the per-:math:`\ell` block-diagonal :math:`\Lambda\phi`

The trailing ``...`` is the **spectator broadcast**: it matches the
spatial-moment axis (if present) and contracts nothing over it — exactly the
:math:`\otimes I_{\rm spatial}` of :eq:`ld-ubld-scattering-moment-lift`.  This
is the `coding-elegance` Pattern-7 producer-side lift: the convention
(scatter each spatial moment independently) is normalised at the producer, so
no consumer special-cases the axis.  Two properties follow by construction:

* **Byte-identical at the single-moment closures.**  When ``phi`` is rank-2
  (``(ng, n_cells)``, the DD/Step shape) the trailing axis is ABSENT, ``...``
  matches nothing, and ``"fg,fc...->gc..."`` is the SAME contraction as
  ``"fg,fc->gc"`` — verified rank-2-exact
  (``np.array_equal`` of the two einsums when no trailing axis is present).
  No re-baseline of the DD/Step path; this is the negative-control
  bit-identity.

* **The projection pair needed no change.**  The spherical-harmonic
  :class:`~orpheus.numerics.frame.GalerkinFrame`'s ``analysis`` and
  ``reconstruction`` faces already carry ``...`` for their trailing
  axes, so :math:`M` and :math:`R` are spatial-moment-agnostic out of
  the box.  The angular reduction
  :eq:`two-moment-angular` and its inverse ride a spatial-moment axis as a
  spectator, which is the architectural payoff of the orthogonal-factor
  framing — the two moment axes never need to know about each other.

.. warning::

   The crosswalk and the original brief ASSUMED ``apply_p0_in_scatter``
   already broadcast over a trailing axis.  It did NOT: the bare ``"fg,fc->gc"``
   hard-codes the cell axis as a single index ``c``, so a rank-3
   ``phi (ng, n_cells, 2^d)`` RAISES (``operand has more dimensions than
   subscripts``).  The fix is the explicit ``...`` spectator, not a reshape.
   Documented so a future reader does not re-derive the (false) assumption.

What is still owed (the iterate carrier and the source seams)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The lift accepts a slope source, but nothing in the production path yet
PRODUCES one.  Filling the slope rows requires (all S3-A proper, owed):

* The :math:`\hat\phi` **iterate carrier** — the between-sweep flux must
  carry the spatial-moment axis.  This was the make-or-break design decision:
  the typed-field spaces validate ``shape == (ng, *spatial)`` with no slot for
  a trailing :math:`(\text{per\_axis})^d` axis, so a slope-carrying field was
  an *illegal state* (Pattern 4 firing correctly).  The resolution — minting
  the first-class :ref:`spatial-moment-space` factor — is the subject of the
  next subsection.

* The **cell-emit moment accumulation** — the wavefront cell solve already
  computes a :math:`(\text{per\_axis})^d`-moment ``psi_avg``, but the
  between-sweep emit currently drops to slot 0 (the cell average); it must
  accumulate the full moment vector.

* The **two source seams** — the :math:`d \ge 2` wavefront genuine
  :math:`(2^d, ng)` moment source through the dense ``_ubld_system``, and the
  :math:`d = 1` scan slope source threaded via
  :meth:`~orpheus.transport.spatial._ubld.D1ClosedForm.kernel_rhs` and the Schur
  ``schur_xV`` term.

The verification chain for the completed S3 is the thick-diffusion-limit
VALUE anchor (the continuous diffusion solution at :math:`\varepsilon \to 0`,
structurally independent of the LD kernel — Adams 2001 / Border-Lewis-Adams
1992 / Larsen-Morel-Miller 1987), the convergence-order MMS smoke (the slope
source exercised), and the genuine Mode-9 SI :math:`\equiv` Krylov on
:math:`(L + C - S_{\rm full})`.  The design record is
``.claude/plans/issue_240_d5b_s3_crosswalk.md``; the verification spec is
``.claude/agent-memory/test-architect/d5b_s3_inc_c_moment_iterate_verification.md``;
the lift's landed-half closeout is
``.claude/agent-memory/method-implementer/issue_240_d5b_s3_a_inc_c_closeout.md``.

All three owed items have SINCE LANDED in the subsections that follow:
the :math:`\hat\phi` iterate carrier and the cell-emit moment
accumulation with the unified moment matvec
(:ref:`ld-ubld-unified-moment-matvec`), and the two source seams with
that matvec's :math:`d \ge 2` dense path and the :math:`d{=}1` moment
scan (:ref:`ld-ubld-moment-scan`).  The owed-list above is preserved as
the campaign-time boundary record.

.. _spatial-moment-space:

The SpatialMomentSpace: a first-class within-cell DG moment carrier (S3-A0)
-------------------------------------------------------------------------------

The typed-field-space half of S3-A (the half the scattering-lift TODO above
flagged as a hard prerequisite).  The within-cell tensor-Legendre DG moment
axis — how :math:`\psi` varies in space WITHIN a cell — is minted as a
first-class function space,
:class:`~orpheus.numerics.spaces.spatial_moment_space.SpatialMomentSpace`,
the **spatial** sibling of the **angular**
:class:`~orpheus.numerics.spaces.spherical_harmonic_space.SphericalHarmonicSpace`.
The two "moment" notions are ORTHOGONAL axes (angular harmonics over
direction :math:`\Omega` vs spatial Legendre over within-cell position
:math:`x`); naming each as its own typed factor keeps the distinction
type-visible and dispels the collision.

.. math::
   :label: spatial-moment-space-size

   \dim(\text{SpatialMomentSpace}) \;=\; (\text{per\_axis})^{d},
   \qquad
   \text{per\_axis} =
   \begin{cases}
     1 & \text{DD / Step (cell-average } \{1\}\text{)} \\
     2 & \text{LD (linear } \{1, P_1\}\text{)}
   \end{cases}

The factor composes via the tensor product ``*`` into the bulk-field spaces
EXACTLY as the angular factor does, and is recovered by type via
``space.find_factor(SpatialMomentSpace).per_axis`` (#207).  The
field-space factories
(:meth:`~orpheus.transport.fields._bases.AngularField.from_mesh`,
:meth:`~orpheus.transport.fields._bases.ScalarField.from_mesh`,
:meth:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux.from_mesh_and_L`)
gained an OPTIONAL ``spatial_moments`` parameter (default ``1``) that
appends the factor **iff the within-cell count exceeds 1** — the
"append iff > 1" gate single-sourced from
:func:`~orpheus.numerics.spaces.spatial_moment_space.spatial_moment_tail`
(the cell analogue that delegates to
:func:`orpheus.numerics.moment_layout.face_moment_tail`, so the cell-moment tail
and the per-face cochain tail can never disagree).  At the default the
field space is BYTE-IDENTICAL to its pre-S3 shape for EVERY scheme (DD,
Step, AND LD): this step builds the CAPABILITY only (construct-general /
select-narrow), and no production field selects the axis yet.

Why a first-class typed factor, not a bare int axis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The slope-carrying flux could, in principle, be stored as a plain ndarray
with a trailing :math:`(\text{per\_axis})^d` axis and an integer remembered
somewhere for its width.  That was rejected (the user's design choice
"option b") in favour of a first-class typed
:class:`~orpheus.numerics.spaces.spatial_moment_space.SpatialMomentSpace`
factor, for the `coding-elegance` Pattern-4 reason
(*make illegal states unrepresentable*).

The typed-field layer validates ``values.shape == space.shape`` at
construction (:meth:`Field.__post_init__`).  Before this step the SN field
spaces were rigidly ``(ng, *spatial)`` / ``(N, ng, *spatial)`` /
``(L+1, 2L+1, ng, *spatial)`` — there was **no slot** for a trailing
spatial-moment axis, so a :math:`\hat\phi`-carrying field FAILED the gate.
A trailing slope axis was, literally, an illegal state (the gate was firing
*correctly* — see the scattering-lift status admonition above, which flagged
this as the make-or-break prerequisite).  Two ways to make the slope-carrying
field legal:

.. list-table:: Bare-int axis vs first-class typed factor
   :header-rows: 1
   :widths: 22 39 39

   * - Aspect
     - Bare ``int`` trailing axis
     - Typed :class:`SpatialMomentSpace` factor
   * - Field validity
     - widen the shape gate to accept *any* trailing axis (loses the
       illegal-state guard)
     - the space DECLARES the axis; the gate stays exact — a slope field
       is now a *legal, declared* shape
   * - Querying the width
     - thread an ``int`` parameter through every call site, or re-derive
       it from a raw ``.shape[-1]``
     - ``space.find_factor(SpatialMomentSpace).per_axis`` — query by
       TYPE, position-independent (#207)
   * - Self-description
     - the axis is anonymous; reading code cannot tell a spatial-moment
       axis from any other trailing axis
     - the factor's type IS its documentation; the
       angular/spatial collision is dispelled at the type level
   * - Precedent
     - none — a one-off convention
     - the EXACT mold of the angular
       :class:`SphericalHarmonicSpace` factor (one architecture, two axes)

The typed factor is the same pattern the angular moment already uses: the
harmonic factor is a :class:`SphericalHarmonicSpace` whose ``L`` is recovered
by ``space.find_factor(SphericalHarmonicSpace).L``, NOT a bare integer
threaded through the API.  Minting the spatial sibling keeps the two axes
*symmetric* — the orthogonal-factor framing of :ref:`two-moment-axes` is then
literally how the carrier space is built (:eq:`two-moment-carrier-space`).

.. note::

   Closing #207 as a side effect.  The
   ``space.find_factor(SphericalHarmonicSpace).L`` query was DOCUMENTED in the
   :class:`HarmonicMomentFlux` docstrings (issue #207) but had never been
   IMPLEMENTED — three docstrings referenced a method that did not exist.  The
   spatial-moment work needed exactly this composition-tree query, so
   :meth:`~orpheus.numerics.space.TensorProductSpace.find_factor` was minted
   now: it returns the first tensor factor that ``isinstance(factor, T)`` and
   raises :exc:`KeyError` if absent (a structural assertion — the caller
   believes the composed space carries the factor — not a silent ``None``,
   Pattern 4).  Both moment factors (angular and spatial) are now queryable
   by type, and the latent broken claim in the docstrings is made true.

The Kronecker moment ordering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The within-cell basis is the tensor (Kronecker) product of the per-axis 1-D
Legendre basis, ordered **x-outer / y-inner** (matching the UBLD assembler
:func:`orpheus.transport.spatial._ubld.assemble_ubld`).  The convention is fixed so
that the all-:math:`P_0` cell average is ALWAYS slot 0:

.. math::
   :label: spatial-moment-kronecker-order

   d{=}1:&\quad [\,\bar\psi,\ \hat\psi_x\,]
   \\[2pt]
   d{=}2:&\quad [\,\bar\psi,\ \hat\psi_y,\ \hat\psi_x,\ \hat\psi_{xy}\,]
   \\[2pt]
   d{=}3:&\quad [\,\bar\psi,\ \hat\psi_z,\ \hat\psi_y,\ \hat\psi_{yz},\
                  \hat\psi_x,\ \hat\psi_{xz},\ \hat\psi_{xy},\ \hat\psi_{xyz}\,]

The slot-0 (cell-average) convention is single-sourced from
:data:`orpheus.numerics.moment_layout.AVERAGE_MOMENT` (the constant every moment
consumer reduces on) — the :class:`SpatialMomentSpace` surfaces it via
:attr:`~orpheus.numerics.spaces.spatial_moment_space.SpatialMomentSpace.average_moment_index`
rather than re-spelling the literal ``0`` (so a layout change happens in ONE
place, not at the scattered ``[..., 0]`` call sites).  Slot 0 is the link
between the two scales: the cell-average moment :math:`\bar\psi` is what the
DD/Step closure carries in full, and it is the moment the
:math:`\Sigma_s \otimes I_{\rm spatial}` lift scatters at every closure; the
remaining slots are the LD-only slope rows the lift activates.

The "append iff > 1" byte-identity policy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The backward-compatibility invariant (#240 D5b) is that DD/Step field shapes
must be UNCHANGED.  This is enforced by a single policy, single-sourced so the
cell-moment tail and the per-face cochain tail can never drift apart:

.. math::
   :label: spatial-moment-append-policy

   \texttt{tail}(n) =
   \begin{cases}
     ()    & n = 1 \quad\text{(DD/Step — NO length-1 axis appended)} \\
     (n,)  & n > 1 \quad\text{(LD — a genuine trailing moment axis)}
   \end{cases}

The critical detail is that ``n == 1`` returns the EMPTY tuple, NOT
``(1,)`` — a length-1 axis is NOT appended.  Appending ``(1,)`` would
broadcast-equal the old shape numerically but would change ``ndarray.shape``
and ``ndim``, breaking every byte-identity gate and every consumer that reads
``.ndim``.  The empty-tuple branch keeps the DD/Step field space *literally
identical* to its pre-S3 self.  :func:`orpheus.numerics.moment_layout.face_moment_tail`
owns the policy; the cell analogue
:func:`~orpheus.numerics.spaces.spatial_moment_space.spatial_moment_tail`
delegates to it (Pattern 7 — normalise the convention at one site), and
:meth:`BulkField._compose_spatial_moments`
(:mod:`orpheus.transport.fields._bases`) returns the space UNCHANGED when the
tail is ``()``.

Construct-general, select-narrow — what this step does and does NOT do
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This step builds the **capability** to carry the spatial-moment axis and
nothing more.  The discipline is deliberate (`coding-elegance` — construct
general, select narrow, specialize only on measured need):

* **The axis exists.**  The :class:`SpatialMomentSpace` is minted, composes
  into every bulk-field space, and is ``find_factor``-queryable.

* **No production field selects it.**  The ``spatial_moments`` factory
  parameter defaults to ``1`` at EVERY call site and is NOT auto-read from
  ``mesh.scheme.spatial_basis_per_axis``.  So DD, Step, AND LD field shapes
  are unchanged this step — even LD does not yet carry the slope axis.

* **Why default-OFF even for LD.**  Auto-reading the scheme would silently
  widen LD field shapes BEFORE the consumers that FILL the axis exist — the
  iterate carrier, the cell-emit accumulation, the source seams (all S3-A
  proper, owed; see the scattering-lift subsection).  A widened axis that no
  producer fills is precisely the illegal state Pattern 4 forbids; turning the
  capability on before its producers exist would re-introduce it.  The gate
  has teeth on exactly this mistake: making
  :meth:`BulkField._compose_spatial_moments` auto-read the scheme turns the LD
  byte-identity foundation tests RED (mutation-verified).

The S3-A iterate / cell-emit / source seams that thread the scheme's
``spatial_basis_per_axis`` here (selecting the axis for LD) are the NEXT
sub-step.  When they land, the only change at the factory call sites is
passing ``spatial_moments=scheme.spatial_basis_per_axis``; the validator
already accepts the widened space, and the scattering lift already scatters
its slopes.

Verification (foundation-level)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The space and the factory widening are software invariants, so they are
verified at the **foundation** level (data-structure / factory-output
invariants, not an L0/L1/L2 solver claim — they carry no eigenvalue or flux
assertion), in two test modules:

* :mod:`tests.numerics.test_spatial_moment_space` — the space layer: the
  :math:`(\text{per\_axis})^d` size law, the
  :meth:`~orpheus.numerics.space.TensorProductSpace.find_factor` round-trip
  (and the :exc:`KeyError`-when-absent assertion), the composition shape, the
  ``per_axis == 1`` no-widening case, and
  :attr:`~orpheus.numerics.spaces.spatial_moment_space.SpatialMomentSpace.average_moment_index`
  :math:`==` :data:`~orpheus.numerics.moment_layout.AVERAGE_MOMENT`.

* ``tests/numerics/test_spatial_moment_field_space.py`` — the factory
  widening: the **byte-identity-at-default negative control** for DD AND LD on
  all three carriers (:class:`AngularField`, :class:`ScalarField`,
  :class:`HarmonicMomentFlux`), the widened :math:`d{=}1` / :math:`d{=}2`
  shapes, the both-moment-factors-coexist case, and the wrong-shape rejection.
  The mutation check — auto-reading the scheme turns the LD byte-identity
  cases red — is what proves the construct-general gate has teeth.

The design record (the angular-vs-spatial distinction, the FP resolution) is
``.claude/plans/issue_240_d5b_s3_crosswalk.md``; the closeout is
``.claude/agent-memory/method-implementer/issue_240_d5b_s3_a0_spatial_moment_space_closeout.md``.

.. _ld-ubld-unified-moment-matvec:

The unified moment matvec: a forward apply is intrinsically moment-valued (S3)
----------------------------------------------------------------------------------

Sub-step **D5b-S3** completes the apply direction: applying the per-cell
:math:`2^d \times 2^d` UBLD operator (:eq:`ld-ubld-cell-system`) to a moment
vector is intrinsically moment-valued, so the matvec carries the full
:math:`(\bar\psi, \hat\psi)` moment vector in every dimension.  The
architectural payoff is a branch removal; the *physics* payoff is the recovery
of the thick-cell diffusion limit, which hinges on a single
frame-consistency identity (ERR-061) that the rest of this subsection derives.
The source files are :mod:`orpheus.transport.spatial.linear_discontinuous`
(``cell_kernel_batch`` / ``residual_kernel_batch`` — now ONE :math:`d`-generic
moment path), :mod:`orpheus.sn.loss_representation.sweep_graph` (``_CellSolve`` / ``_CellResidual``
— the ``len(s_axes) > 1`` moment gate retired), and
:mod:`orpheus.sn.loss_representation` (the ``_spatial_moment_tail`` buffer
widening); the closeout is
``.claude/agent-memory/method-implementer/issue_240_d5b_s3_unified_matvec_closeout.md``.

The earlier increments (A/B) made the :math:`d{=}1` LD **matvec** Schur-reduce
to a scalar residual: the slope :math:`\hat\psi` was eliminated, leaving a
scalar cell-average unknown.  That was a *flat-source artifact* — with
:math:`\hat Q = 0` the slope had no global coupling, so the Krylov unknown could
stay scalar.  Increment C makes the scattering slope source
:math:`\Sigma_s\hat\phi` couple the slope GLOBALLY (the diffusion-limit-
consistent operator :eq:`ld-ubld-scattering-moment-lift`), so the slope becomes a
genuine global degree of freedom in **every** dimension.

.. math::
   :label: ld-ubld-unified-moment-residual

   (L+C)\,\vec\psi
   \;=\;
   M^{-1}\bigl(A\,\vec\psi - F_{\rm in}\bigr),
   \qquad
   A = (L+C-S)\ \Longleftrightarrow\
   (L+C)\,\vec\psi - S\,\vec\psi = \vec q_{\rm ext}

A matvec is a forward APPLY: applying the per-cell
:math:`2^d \times 2^d` UBLD operator to the moment vector is intrinsically
moment-valued, so ``cell_kernel_batch`` and ``residual_kernel_batch`` collapse to
ONE d-generic dense path for every :math:`d` (the former :math:`d{=}1`
Schur-reduced scalar arm — and the :math:`d \ge 2` raise — are both retired).
The :math:`M^{-1}` factor in :eq:`ld-ubld-unified-moment-residual` is the
matvec/sweep moment-source consistency: the UBLD RHS folds the cell source
mass-weighted (:math:`R = M\vec S`, the test-function projection), but the
operator algebra :math:`A = (L+C) - S` subtracts :math:`S\vec\psi` RAW at the
``OperatorSum`` level, so the residual is divided by the diagonal Legendre mass
to put :math:`(L+C)\vec\psi` in raw per-moment units (the slope rows would
otherwise disagree by :math:`M_{ii} = \theta^{|i|}`).

The architectural headline is **branch removal** (Cardinal Rule 2): the
``len(s_axes) > 1`` moment gate at the cell-solve / cell-residual emit is GONE
(replaced by the pure scheme trait ``spatial_basis_per_axis > 1``), the
:math:`d{=}1` scalar kernel twin is retired, and there are ZERO
``isinstance(scheme, ...)`` branches — dispatch stays via the scheme PROTOCOL +
geometry-keyed ``supports()``.

.. _ld-ubld-sweep-global-frame:

The sweep-frame / global-frame involution (ERR-061 — the diffusion-limit root cause)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is the load-bearing result of the whole multi-dimensional LD campaign, and
it is the kind of "what failed and why" Cardinal Rule 3 demands be archived in
full.  Threading the moment matvec made the operator *internally consistent*
(matvec :math:`\equiv` sweep round-trip to :math:`10^{-16}`, source-iteration
:math:`\equiv` Krylov on the SAME operator) — and yet the converged flux was
**wrong**: on a thick scattering slab (:math:`\Sigma_t = 40`, :math:`c = 0.99`,
:math:`\Sigma_t h = 10`/cell, vacuum) at :math:`n_x = 4` the LD scalar flux was
1.47 against the diffusion solution 2.31 (relative error 36 %), and the error
did not grow under refinement — it *shrank* as the cells thinned (the classic
flat-source-LD signature, persisting THROUGH the slope-source thread).

The cause is a frame-consistency error between two individually-correct
components.  The per-cell LD kernel produces and consumes the :math:`2^d` moment
vector in the per-ordinate **sweep frame**: each axis :math:`a` is oriented so
the *downstream* face is at the local coordinate :math:`+1`.  For an ordinate
sweeping in the NEGATIVE global direction on axis :math:`a`
(:math:`\mathrm{octant\_sign}_a = -1`) the sweep coordinate is the *reverse* of
the global coordinate, so the slope (:math:`P_1`) moment on that axis is
sign-FLIPPED relative to the global-:math:`x` slope.  But the iterate
:math:`\hat\phi` and its scattering source :math:`\Sigma_s\hat\phi` live in the
**global frame** — the angular reduction sums slopes across ordinates of BOTH
sweep directions,

.. math::
   :label: ld-ubld-slope-angular-reduction

   \hat\phi(\vec r, g) \;=\; \sum_{n=1}^{N} w_n\,\hat\psi_n(\vec r, g).

The producer (``_CellSolve.cell`` emit) stored the raw sweep-frame slope; the
consumer (``integrate_angular`` / the scattering apply) summed it as if it were
global-frame.  So the backward ordinates' opposite-signed slopes partially
CANCELLED the forward ones: at a cell with a positive global-:math:`x` gradient
the forward ordinate had :math:`\hat\psi_n = +0.048` but the backward had
:math:`-0.028` — opposite signs, the smoking gun.  The summed
:math:`\hat\phi` was :math:`\sim 6\times` too small to satisfy the LM-1989
discrete diffusion continuity (Larsen & Morel 1989, JCP 83(1):212–236, Eq. 4.9b,
:math:`\bar\phi_j + \hat\phi_j = \bar\phi_{j+1} - \hat\phi_{j+1}`), the slope
source was under-driven, and the discrete diffusion limit was the wrong
diffusion operator.

The fix is a single-sourced :math:`2^d` moment-frame **involution**,
:func:`~orpheus.transport.spatial._ubld.octant_moment_frame_signs`,

.. math::
   :label: ld-ubld-octant-moment-frame-signs

   \mathrm{sign}[o_0, \ldots, o_{d-1}]
   \;=\;
   \prod_{a=0}^{d-1} (\mathrm{octant\_sign}_a)^{\,o_a},
   \qquad o_a \in \{0, 1\},

indexed in the tensor-Legendre Kronecker layout (:math:`o_a` = the :math:`P_0` /
:math:`P_1` selector on axis :math:`a`).  The **average** moment (all
:math:`o_a = 0`) is sign-invariant (the empty product is :math:`1`); a per-axis
**slope** flips once if that axis sweeps backward; the 2-D **cross** moment
:math:`\hat\psi_{xy}` flips when an ODD number of its active axes reverse.  The
map is its own inverse, so the SAME sign vector converts global :math:`\to`
sweep on the source/probe INPUT and sweep :math:`\to` global on the
moment/residual OUTPUT.  It is applied through the shared ``_reframe`` helper at
the cell ops; the OUTGOING FACE (``psi_out``) stays sweep-frame — it propagates
along the wavefront and never crosses into the global-frame iterate (so it is
left untouched).  DD/Step (``spatial_basis_per_axis == 1``) get ``None`` (the
sign-invariant average-only moment), so they pass through ``_reframe``
untouched and stay byte-identical (the negative control); a flat scalar source
(matvec zero / flat external — only the average moment) is frame-invariant and
skipped by the ``arr.shape[-1] != frame_signs.shape[0]`` guard, so it is never
broadcast into a spurious moment axis.

After the fix the diffusion limit is recovered on BOTH the matvec (Krylov) and
the sweep (source-iteration) paths:

.. list-table:: Thick-slab LD vs DD relative error, before/after the frame fix
   :header-rows: 1
   :widths: 14 22 22 22

   * - Mesh
     - Cell optical depth
     - Before (sweep-frame slope)
     - After (global-frame slope)
   * - 1-D, :math:`n_x = 4`
     - :math:`\Sigma_t h = 10`
     - 38.9 %
     - 4.1 %
   * - 1-D, :math:`n_x = 16`
     - :math:`\Sigma_t h = 2.5`
     - 7.9 %
     - 0.2 %
   * - 1-D, :math:`n_x = 64`
     - :math:`\Sigma_t h = 0.6`
     - 0.9 %
     - 0.0 %
   * - 2-D, :math:`n = 4/8/16`
     - thick :math:`\to` thin
     - 8.4 %
     - 1.7 % :math:`\to` 0.4 %

.. warning::

   **The matvec-self-consistency gate is necessary but NEVER sufficient for a
   moment-iterate fold.**  Every component here was individually correct (the
   2×2 matched LM-1989 Eq. 4.3, the dense UBLD matched the analytic 2×2, the
   scattering produced :math:`\Sigma_s\hat\phi` at full strength, the matvec
   round-trip vanished to :math:`10^{-16}`, and source-iteration :math:`\equiv`
   Krylov on the SAME operator).  The bug was the frame consistency *between*
   two correct components — a wrong fixed point that the round-trip and the
   SI :math:`\equiv` Krylov gates are structurally BLIND to: they prove the
   operator is internally consistent, NOT that its fixed point is the
   physically-correct one (`vv-principles` §5 — "O(h²) to the wrong limit is
   still O(h²)").  The decisive evidence was a structurally-independent
   from-scratch LD-SN solver (a direct LM-1989 2×2 + source iteration, no
   ORPHEUS kernel) that reproduced ORPHEUS's WRONG value bit-for-bit when it
   summed sweep-frame slopes and RECOVERED the diffusion limit when it stored
   global-frame slopes — pinning the root cause independent of ORPHEUS's code.
   The lesson: gate the converged VALUE against a structurally-independent
   reference (the continuous diffusion solution + the independent from-scratch
   kernel), never the round-trip.  This is failure Mode 1 (sign flip) +
   Mode 6 (convention drift) — see ``error_catalog.md`` ERR-061.

The thick-cell diffusion tripwire is
``tests/sn/verification/mms/test_mms_ld_slab.py::test_ld_thick_diffusive_limit``
(1G) and ``::test_ld_thick_diffusive_limit_2g`` (2G heterogeneous, a
group-coupled slope source — Mode 6), both ``@pytest.mark.l1
@pytest.mark.catches("ERR-061")`` and both Mode-8-safe
(``np.testing.assert_array_less`` fires under ``-O``).  The slope-frame
fingerprint is pinned by
``derivations/diagnostics/diag_240_d5b_s3_probe_11_root_cause.py`` (forward and
backward ordinate slopes must share sign in the global frame), and the
structurally-independent confirmation by
``diag_240_d5b_s3_probe_08_independent_ld.py``.

.. _ld-ubld-pure-z-collision-twin:

The pure-z collision-only twin — sweep :math:`\equiv` matvec single source
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::
   :label: ld-ubld-pure-z-collision

   \text{(solve)}\quad \bar\psi = \frac{Q}{\sigma_t}
   \qquad\Longleftrightarrow\qquad
   \text{(matvec)}\quad (L+C)\,\bar\psi = \sigma_t\,\bar\psi

The pure-z degenerate ordinates (:math:`\mu_x = \mu_y = 0`, the :math:`\pm z`
poles of a Lebedev or product cubature in a 2-D Cartesian sweep) have no
in-plane streaming, so the cell is **collision-only**: the loss couples to
:math:`\sigma_t` alone, and the sweep balance :math:`\bar\psi = Q/\sigma_t` and
its matvec twin :math:`(L+C)\bar\psi = \sigma_t\bar\psi` are two applications of
the SAME operator (:eq:`ld-ubld-pure-z-collision`).  This is the L21 twin-path
relationship — sweep and matvec are the same physics evaluated in opposite
directions — and it is exactly the kind of paired closure that drifts when a new
axis lands on one side and is forgotten on the other.

At a multi-moment closure (LD) the source :math:`Q` and the probe
:math:`\bar\psi` carry the trailing :math:`2^d` spatial-moment axis that
:math:`\sigma_t` of shape :math:`(ng, *\text{spatial})` lacks; each moment
scales by the SAME scalar (:math:`1/\sigma_t` on the solve, :math:`\sigma_t` on
the apply), so :math:`\sigma_t` must gain a length-1 trailing axis to broadcast.
This reshape is single-sourced through
:func:`~orpheus.sn.loss_representation._moment_broadcast_sigma`
(:math:`\sigma \mapsto \sigma[\ldots, \text{None}]` iff the moment-valued
operand out-ranks :math:`\sigma`), called by BOTH the sweep ``pure_z`` arm
(:math:`Q\,/\,\texttt{\_moment\_broadcast\_sigma}(\sigma_t, Q)`) and the matvec
``pure_z`` arm
(:math:`\texttt{\_moment\_broadcast\_sigma}(\sigma, \bar\psi)\cdot\bar\psi`), so
the twin cannot diverge on the moment-axis convention.  DD/Step (no moment axis)
:math:`\to` :math:`\sigma_t` unchanged, byte-identical.

.. admonition:: ERR-062 — the matvec twin forgot the guard the sweep already had
   :class: warning

   Before this fix the sweep arm HAD the moment-broadcast guard but the matvec
   arm wrote the bare ``sigma * probe[oct_idx]``, so a moment-valued probe
   broadcast-FAILED.  The consequence:
   ``solve_sn_fixed_source(scheme=LinearDiscontinuous(), inner_solver="krylov")``
   on ANY 2-D Cartesian LD mesh whose quadrature carries pure-z ordinates raised
   ``ValueError`` at the first Krylov matvec.  The bug hid through the whole
   D5b-S3 development because every committed 2-D LD test used
   ``level_symmetric`` — which has NO pure-z ordinates — while the production
   MMS uses a Lebedev quadrature that does.  This is the canonical L21
   twin-path asymmetry recurring a THIRD time ("the matvec needs a committed
   gate, not a round-trip"): the round-trip and FFW :math:`\equiv` MFW gates
   ran on ``level_symmetric`` and never exercised the pure-z arm at all.  The
   gate is
   ``tests/sn/verification/mms/test_mms_ld_2d.py::test_ld_2d_krylov_equals_si_pure_z_quadrature``
   (``@pytest.mark.foundation @pytest.mark.catches("ERR-062")``), on a Mode-9
   degeneracy-break config: a pure-z-bearing Lebedev order-5 quadrature
   (:math:`N = 14`, genuine :math:`\mu_y` + the 2 :math:`\pm z` poles),
   heterogeneous 2-material map, 2-group asymmetric XS with non-zero self-scatter,
   NON-SQUARE :math:`5\times4`, vacuum edges.  Mutation-verified: re-introducing
   the bare ``sigma * probe[oct_idx]`` makes the gate FAIL with the exact
   ``ValueError``; with the fix Krylov :math:`\equiv` SI to :math:`\sim10^{-11}`
   (the same :math:`(L+C-S_{\rm full})` fixed point).  See ``error_catalog.md``
   ERR-062.

The source is :mod:`orpheus.sn.loss_representation` (``_moment_broadcast_sigma``,
the ``loss_action`` matvec ``pure_z`` arm and its source-iteration sweep twin);
the closeout is
``.claude/agent-memory/method-implementer/issue_240_d5b_s3_purez_gate_closeout.md``.

.. _ld-ubld-moment-scan:

The :math:`d{=}1` moment SCAN (the production sweep) — D5b-S3 OWED-2
------------------------------------------------------------------------

The unified moment *matvec* above is the APPLY direction; the production
:math:`d{=}1` LD SWEEP (source iteration) rides the fast Blelloch parallel-prefix
scan (:class:`~orpheus.sn.loss_representation.CumprodScan`), NOT the dense
per-cell solve (L16).  Sub-step **D5b-S3 OWED-2** threads the spatial-moment
iterate :math:`\hat\phi` through that scan so the SI path recovers the SAME
diffusion-limit-consistent operator the matvec does.

.. math::
   :label: ld-ubld-moment-scan-source

   b \;=\; \underbrace{\bar S\,\frac{\mathrm{inv}}{w}}_{\text{flat (average) emission}}
       \;+\; \underbrace{\frac{\theta\,\hat S}{D_2'}
              \;-\; \frac{\theta\,|\mu| A_{\rm down}\,\hat S}{D_2'}\,
                    \frac{\mathrm{inv}}{w}}_{\text{slope source}\ \Sigma_s\hat\phi}

The scan propagates the scalar downstream FACE
:math:`\psi_{\rm out} = a\,\psi_{\rm in} + b` along the cell chain with the
**slope-augmented** affine source :eq:`ld-ubld-moment-scan-source`: the flat
(cell-average) emission :math:`\bar S\,\mathrm{inv}/w` plus the slope-source
contribution :math:`\theta\hat S/D_2' - (\theta|\mu|A_{\rm down}\hat S/D_2')\,\mathrm{inv}/w`
that carries :math:`\Sigma_s\hat\phi` into the recurrence.  Then it reconstructs
the per-cell :math:`(\bar\psi, \hat\psi)` moments from the chained upstream face.
The slope-row :math:`\hat S` algebra is single-sourced through
:meth:`~orpheus.transport.spatial._ubld.D1ClosedForm._slope_fold`, shared by the
per-cell matvec Schur (:meth:`~orpheus.transport.spatial._ubld.D1ClosedForm.schur_xV`)
AND the scan
(:meth:`~orpheus.transport.spatial._ubld.D1ClosedForm.scan_slope_face_source` for the
face-chain term, :meth:`~orpheus.transport.spatial._ubld.D1ClosedForm.scan_reconstruct`
for the per-cell moments).

Why the face/cell split is necessary (the load-bearing math)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A moment-carrying parallel-prefix scan is **NOT** a drop-in widening of the
scalar scan.  For *flat-source* LD the cell average is the convex blend of the
two faces, :math:`\bar\psi = (1-w)\psi_{\rm in} + w\,\psi_{\rm out}`, so the
scalar scan can reconstruct :math:`\bar\psi` directly from the chained faces.
With a *slope* source, that closure **decouples**: :math:`\bar\psi` and
:math:`\psi_{\rm out}` no longer satisfy the convex blend, because the slope
source enters the cell balance without entering the face propagation the same
way.  The scan therefore splits the work in two:

#. the FACE chain :math:`\psi_{\rm out} = a\,\psi_{\rm in} + b` propagates with
   the slope-augmented :math:`b` (:eq:`ld-ubld-moment-scan-source`), so the next
   cell's :math:`\psi_{\rm in}` is the correct dense :math:`\bar\psi + \hat\psi`;

#. the CELL moments :math:`(\bar\psi, \hat\psi)` are reconstructed per cell from
   the chained :math:`\psi_{\rm in}` via the per-cell Schur
   (``scan_reconstruct``), **not** via ``cell_average``.

Conflating the two (using ``cell_average`` for :math:`\bar\psi` on the moment
scan) gives the WRONG cell average while the face chain still looks right — the
silent trap this split avoids.  The reconstruction was verified against a
from-scratch dense :math:`d=1` chain (face / :math:`\bar\psi` / :math:`\hat\psi`
all match to :math:`10^{-12}`) and against the live DAG (scan :math:`\equiv` DAG
to :math:`10^{-16}` on a 2G-heterogeneous non-flat config).

The same global-frame involution as the matvec
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Like the matvec, the scan applies the SAME
:func:`~orpheus.transport.spatial._ubld.octant_moment_frame_signs` involution
(:eq:`ld-ubld-octant-moment-frame-signs`) through the shared ``_reframe``
helper: the source moments are mapped global :math:`\to` sweep on INPUT and the
reconstructed :math:`(\bar\psi, \hat\psi)` sweep :math:`\to` global on OUTPUT, so
the angular reduction :math:`\hat\phi = \sum_n w_n \hat\psi_n`
(:eq:`ld-ubld-slope-angular-reduction`) is frame-consistent and the diffusion
limit is recovered on the source-iteration path too (the sweep-side analog of
the ERR-061 matvec fix).  The backward sweep flips the slope so forward and
backward ordinates reinforce rather than cancel.  The scalar OUTGOING FACE stays
sweep-frame — it propagates along the chain and never crosses into the global
iterate (the :math:`d=1` face cochain is :math:`2^{d-1} = 1`, scalar, so it is
not reframed).  The scan is a **consumer** of the matvec's machinery, never a
twin: the same ``_slope_fold`` powers the per-cell Schur and the scan, and the
same involution powers the DAG and the scan.  DD/Step (``per_axis == 1``) get no
moment axis, so the scan runs the existing flat slab body verbatim and stays
byte-identical (the negative control).

The SymPy module / live scheme is :mod:`orpheus.transport.spatial._ubld`
(``D1ClosedForm``),
:meth:`orpheus.transport.spatial.linear_discontinuous.LinearDiscontinuous.moment_scan_closure`,
and :meth:`orpheus.sn.loss_representation._OneDimScanWalk._run` (the slab
joint-batch moment branch); the gates are
``tests/sn/verification/mms/test_mms_ld_slab.py::test_ld_two_paths_scan_equals_dag_oracle``
(scan :math:`\equiv` DAG) and ``::test_ld_thick_diffusive_limit`` (the diffusion
limit on the SI path, the same ERR-061 catcher the matvec uses); the closeout is
``.claude/agent-memory/method-implementer/issue_240_d5b_s3_owed2_scan_closeout.md``.


What broadens next
==================

* **Curvature** (spherical/cylindrical): the angular cell balance
  activates — direction is no longer constant in flight, so the
  angular axis acquires its own cell balance, redistribution
  coefficients, and starting-direction state. The walk becomes a
  sequential per-ordinate march; the closure machinery of
  :doc:`/theory/foundations/discretization` §5 is applied on the
  angular axis. That is Part B of this book.
