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


What broadens next
==================

* **Curvature** (spherical/cylindrical): the angular cell balance
  activates — direction is no longer constant in flight, so the
  angular axis acquires its own cell balance, redistribution
  coefficients, and starting-direction state. The walk becomes a
  sequential per-ordinate march; the closure machinery of
  :doc:`/theory/foundations/discretization` §5 is applied on the
  angular axis. That is Part B of this book.
