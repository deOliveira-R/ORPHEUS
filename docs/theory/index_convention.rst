.. _theory-sn-index-convention:

==========================================================
SN Index Convention --- ``(N, ng, nx, ny)``
==========================================================

.. contents:: Contents
   :local:
   :depth: 3


Key Facts
=========

**Read this before reading or writing SN array code.**

This page is the canonical statement of the storage layout for every
SN solver array.  It supersedes every per-file docstring that
described a different layout; if a per-file docstring disagrees with
this page, this page is correct.

- **Angular flux** :math:`\psi`: ``(N, ng, nx, ny)`` --- ordinate
  index first, energy second, spatial last.  1-D problems use
  ``ny = 1`` (the trailing axis is preserved as a singleton, NOT
  squeezed).
- **Scalar flux** :math:`\phi` and **cross sections** (:math:`\Sigma_t`,
  :math:`\Sigma_a`, :math:`\Sigma_p`, :math:`\nu\Sigma_f`, :math:`\chi`):
  ``(ng, nx, ny)`` --- energy first, spatial last.  Same trailing
  ``ny = 1`` rule for 1-D.
- **External source** :math:`q`: ``(ng, nx, ny)`` (isotropic) or
  ``(N, ng, nx, ny)`` (anisotropic / per-ordinate).
- **Cell-flattening invariant**: the principled storage round-trips
  with the legacy one under transpose:
  ``xs.sig_t.T.reshape(ng, nx, ny)[g, i, j] ==
  xs.sig_t.reshape(nx, ny, ng)[i, j, g]``.  Asserted in ``__debug__``
  at :class:`~orpheus.sn.solver.SNSolver` construction.
- **The four-operator algebra** :math:`(L + C - S - F/k)\psi = q`
  consumes and returns :math:`\psi` shaped as ``(N, ng, nx, ny)`` at
  every leaf (:class:`~orpheus.sn.operator.StreamingOperator`,
  :class:`~orpheus.sn.operator.CollisionOperator`,
  :class:`~orpheus.sn.scattering.ScatteringOperator`,
  :class:`~orpheus.sn.fission.FissionOperator`).
- **One internal exception is deliberately preserved**: the
  FD-matvec packed-vector helper
  :meth:`~orpheus.sn.operator.SNStreamingOperator.solution_to_angular_flux`
  returns ``(ng, N, nx, ny)`` internally.  The Krylov decode path
  absorbs the leading-axis-swap via two zero-copy
  :func:`numpy.transpose` adapters at
  :func:`~orpheus.sn.solver.solve_sn_fixed_source`.  PR-INDEX-7 retires
  this exception; see :ref:`sn-index-convention-future-work`.

.. admonition:: Authoritative origin

   The derivation of this layout is encoded in
   :ref:`sn-index-convention-derivation`.  The full migration
   narrative --- six commits between ``c21c2ef`` and ``3356cec`` that
   flipped the codebase from ``(N, nx, ny, ng)`` to
   ``(N, ng, nx, ny)`` --- is recorded in
   :ref:`sn-index-convention-history`.  Numerical evidence
   (11/11 regression snapshots bit-identical at ``rtol=1e-12`` via the
   load-bearing Step-1 bit-identity-via-transpose gate) is in
   :ref:`sn-index-convention-numerical-evidence`.


Overview
========

The SN solver discretises four indices:

- :math:`n` --- ordinate (angular direction), :math:`n = 0, \ldots, N-1`.
- :math:`g` --- energy group, :math:`g = 0, \ldots, n_g - 1`.
- :math:`i, j` --- spatial cell, :math:`i = 0, \ldots, n_x - 1`,
  :math:`j = 0, \ldots, n_y - 1` (with :math:`n_y = 1` for 1-D).

A storage decision picks one ordering of these axes for every flux,
source, and cross-section array.  The decision is consequential: it
affects how every operator-leaf ``apply`` body indexes, how the
sweep's hot loop traverses memory, and how a future JAX or GPU port
maps batched dimensions to the device grid.

Before Issue #196 the codebase carried ``(N, nx, ny, ng)`` for
:math:`\psi` (energy trailing) and ``(nx, ny, ng)`` for everything
scalar.  That layout was historical: it descended from a 2-D-Cartesian
prototype in which group iteration was the outermost Python loop, so
``g`` trailing was the most natural place to put it.  The principled
order --- derived in :ref:`sn-index-convention-derivation` from the
block-diagonality of the within-group system plus the
streaming-dependency-only-on-cell-axes --- is the opposite:
:math:`g` belongs *second*, immediately after :math:`n`, with the
must-iterate cell axes last.

The migration from the historical layout to the principled one was
done in six PRs over a single branch (``refactor/sn-operator-algebra``)
in May 2026.  Every intermediate commit kept
``tests/sn/regression/`` green at ``rtol=1e-12``.  The final commit
(``3356cec``) regenerated the regression snapshots under the
principled layout *only after* a bit-identity-via-transpose gate
demonstrated that every one of the 11 snapshots agreed with the
legacy layout to ``rtol=1e-12, atol=1e-13``.


.. _sn-index-convention-derivation:

Derivation --- why ``(N, ng, nx, ny)`` is principled
====================================================

The choice of axis order is dictated by the structure of the
within-group transport system.  For a single energy group :math:`g`,
the discrete ordinates equation is

.. math::
   :label: sn-within-group-system

   (L_g + C_g)\,\psi_g \;=\; S_g\,\phi + \frac{1}{k}\,\chi_g\,F\,\phi
   + q_g\,,

where :math:`L_g` is the streaming operator on group :math:`g`,
:math:`C_g = \Sigma_{t,g}\,\mathbb{I}` is the diagonal collision
operator, :math:`S_g` accumulates the in-scatter contribution, and
the right-hand side carries the fission source and external source.
The *within-group* system is the per-:math:`g` problem when the
scattering source is held fixed at the current outer iterate ---
exactly what each inner source iteration solves.

Two observations make the storage decision:

1. **No cross-group coupling within a sweep.** The within-group
   problem :eq:`sn-within-group-system` is solved for each :math:`g`
   independently --- the operator :math:`L_g + C_g` is
   *block-diagonal in g*.  Once :math:`\phi^{(k)}` is fixed for the
   inner iteration, the energy group is the obvious
   joint-batch axis: every group-independent quantity (the
   streaming operator, the collision diagonal, the per-ordinate
   intermediate arrays) can be processed in lockstep across all
   groups using a single numpy expression.  This is the same
   structural reason Galerkin spectral methods batch over the
   variational basis.

2. **No cross-ordinate coupling for the within-group P\ :sub:`0`
   problem.** When scattering is :math:`P_0` (isotropic), the
   right-hand side is a scalar flux :math:`\phi(r)` that does not
   depend on the outgoing ordinate.  Different ordinates therefore
   solve *the same source* through different streaming directions,
   and an outer Krylov batch can compute their residuals in parallel.
   The block structure is therefore **block-diagonal in both g and
   n** at the within-group P\ :sub:`0` level.  Curvilinear angular
   redistribution and :math:`P_\ell` anisotropic scattering reduce
   this independence but do not destroy it: the redistribution is
   one tight ordinate band (the M--M angular thread per
   :math:`\mu`-level), not the full :math:`N \times N` dense
   coupling.

3. **The spatial axes are the only must-iterate axes.** The
   streaming term :math:`\mu_x \partial_x \psi + \mu_y \partial_y
   \psi` connects every cell to its upwind neighbour.  Sweep order is
   dictated by the DAG of cell dependencies; the sweep
   *fundamentally cannot* be parallelised across cells along the
   streaming direction.  The cell axes are therefore the innermost
   axes --- the ones we *want* to traverse sequentially.

These three observations give the priority ordering for storage axes:

.. list-table:: Axis-priority table (principled storage)
   :header-rows: 1
   :widths: 12 32 56

   * - Index
     - Within-group coupling
     - Storage role
   * - :math:`n` (ordinate)
     - None for within-group P\ :sub:`0`; tight band for
       curvilinear / P\ :sub:`\ell`
     - **Outermost** --- sweep iterates over chains
       outside the per-ordinate body; Krylov batches
       across :math:`n`
   * - :math:`g` (group)
     - None for within-group
     - **Second** --- block-diagonal axis, joint-batched in
       every per-ordinate kernel
   * - :math:`i, j` (cells)
     - Streaming dependency chain
     - **Innermost** --- the only axes that *must* iterate
       sequentially

The principle generalises:

   **In a tensor-product discretisation, axes with no cross-coupling
   for the within-group system belong before the axes that carry a
   sequential dependency.**

[LewisMiller1984]_ §4.5 ("Source Iteration") gives the same
block-diagonal structure as the textbook proof that the within-group
problem decouples; [AdamsLarsen2002]_ §III confirms the same picture
for the SAILOR preconditioner family.  The block structure is the
mathematical reason every modern transport code (PARTISN, Denovo,
JAGUAR, OpenMOC) carries the same ``(angular, energy, spatial)``
priority --- the storage layout is dictated by the operator algebra,
not by a historical implementation choice.

Algorithmic consequence
-----------------------

Under the principled layout, the per-sweep hot path

.. code-block:: python

   # For SLAB: ordinates within a chain are not coupled (no M-M thread).
   # Joint-batch over (chain_size, ng, nx) — one scan per chain.
   psi_face_chain = ordinate_scan(
       a_atten_chain.T,    # (nx, K, ng) — scan axis leads (Blelloch)
       b_chain.T,          # (nx, K, ng)
       psi_in,
   )

becomes a single
:func:`~orpheus.sn.spatial.scan.ordinate_scan` call per chain (two
chains per slab problem), rather than ``N/2`` per-ordinate calls.
For ``N = 16``, this saves 14 Python invocations per sweep ---
roughly 28 % mean speedup on the
``slab_2g_3reg_dd_n40`` regression case (PR-INDEX-1 benchmark; see
:ref:`sn-index-convention-numerical-evidence`).

The same principle would let curvilinear sweeps joint-batch over
groups too, but the M--M angular thread is sequential *across*
ordinates within a :math:`\mu`-level, so the curvilinear sweep keeps
the per-ordinate scan with ``(ng, nx)`` joint batching only.  A
parallel-prefix reformulation of the M--M recurrence (research-level
algorithm work) could lift this restriction; see
:ref:`sn-index-convention-future-work`.


Cross-section convention
========================

Cross sections follow the same priority as the scalar flux:
:math:`g` first, then spatial.  Per-cell cross sections are stored as
``(ng, nx, ny)`` numpy arrays on
:class:`~orpheus.sn.solver.SNSolver`:

.. code-block:: python

   class SNSolver:
       sig_t: np.ndarray   # (ng, nx, ny) total
       sig_a: np.ndarray   # (ng, nx, ny) absorption
       sig_p: np.ndarray   # (ng, nx, ny) production (νΣ_f)
       chi:   np.ndarray   # (ng, nx, ny) fission spectrum

The producer :func:`~orpheus.data.macro_xs.assemble_cell_xs`
emits the flat ``(N_cells, ng)`` shape (CP also consumes that flat
shape --- the producer is *unchanged* by the SN migration).  The
``.T.reshape(ng, nx, ny)`` bridge lives at exactly one site,
:meth:`SNSolver.__init__`, and the
:ref:`cell-flattening invariant <sn-cell-flattening-invariant>`
asserts the round-trip in ``__debug__`` builds.

.. _sn-cell-flattening-invariant:

Cell-flattening invariant
-------------------------

The principled storage must agree with the legacy storage under a
pure-transpose round-trip.  The check is

.. math::
   :label: sn-cell-flatten-roundtrip

   \texttt{sig\_t}_{\text{principled}}[g, i, j]
   \;=\;
   \texttt{sig\_t}_{\text{legacy}}[i, j, g]
   \qquad \forall (g, i, j)\,,

implemented at :meth:`SNSolver.__init__`:

.. code-block:: python

   xs = assemble_cell_xs(materials, sn_mesh.mat_map)
   self.sig_t = xs.sig_t.T.reshape(self.ng, nx, ny)
   if __debug__:
       _sig_t_old = xs.sig_t.reshape(nx, ny, self.ng)
       assert np.array_equal(
           _sig_t_old, self.sig_t.transpose(1, 2, 0)
       ), "PR-INDEX-3 cell-flattening invariant broke"

The invariant is load-bearing: it detects accidental mat-ids ravel
order changes (Fortran vs C order) that would silently corrupt the
spatial-to-group mapping.  An assertion failure here would surface as
a clean test failure rather than a flux distribution that looks
plausible but is wrong by a permutation of cells.


.. _sn-index-convention-history:

History --- the six-PR migration
================================

The migration unfolded as a six-commit chain on the
``refactor/sn-operator-algebra`` branch between 2026-05-14 and
2026-05-15.  Each PR kept ``tests/sn/regression/`` green at
``rtol=1e-12`` by inserting temporary bridge transposes at the
boundary between flipped (principled) and unflipped (legacy) layers.
The bridges were named ``BRIDGE_*_to_principled`` /
``BRIDGE_*_to_legacy`` so a grep-tag retired them as the migration
progressed.

The proposal that was wrong
---------------------------

The initial typed-field contract memo
(``.claude/agent-memory/explorer/typed_field_contracts_for_phase_g.md``,
committed at ``9d74184``) proposed ``(N, nx, ny, ng)`` as the
canonical storage --- energy trailing.  The memo's argument was
operational: numpy ``block_op @ flux.values`` works with the last
axis as the contraction axis, and many group-block matrices act on
:math:`g` as the inner axis.  That argument is locally correct but
**inverts the coupling priority**: it puts the must-iterate cell axes
``(nx, ny)`` *before* the block-diagonal group axis ``ng``, which
forces every per-:math:`g` numpy operation to reach across the cell
axes to find its group block.  For a strided memory traversal this
wastes cache; for a future GPU port it would map the
block-diagonal-but-stride-:math:`n_x \cdot n_y` axis to the wrong
grid dimension.

The wrong proposal was caught **before implementation** by re-reading
the derivation table in this page's §1 against the memo's §1.1.  The
discovery point is documented in the migration plan (`§1.1 "Why we
paused the typed-field contract plan to do this first"
<../../.claude/plans/principled_index_migration.md>`_) and is the
canonical example of:

   **When a refactor calls for new types AND a layout change, do the
   layout change first on bare arrays.  Types ossify the layout;
   flipping bare arrays is mechanical, flipping a layout that's
   encoded in twelve dataclass ``__init__``\ s + their dunder
   consumers is not.**

This is a concrete instance of the *defer abstraction until you have
evidence* principle from the project's ``coding-elegance`` skill
(Pattern 6): the team had **one** concrete instance of "the layout
we want" (the four-operator algebra acceptance criterion), but the
layout itself was still wrong.  Build the layout first; build the
abstraction on top.

The six PRs
-----------

.. list-table:: Migration commit chain
   :header-rows: 1
   :widths: 14 12 74

   * - PR
     - Commit
     - Scope
   * - PR-INDEX-1
     - ``e09b9f8``
     - ``_run_1d_sweep`` internal layout flipped to principled
       ``(N, ng, nx, ny)``; slab joint-batch ``ordinate_scan`` over
       ``(chain_size, ng, nx)`` --- 2 scan calls per sweep replacing
       ``N/2``.  Public ``transport_sweep`` signature unchanged
       (entry / exit transposes carry the legacy boundary).  11/11
       regression bit-identical at ``rtol=1e-12``; 26/26 L0
       streaming-equilibrium curvilinear; 312/312 spatial.
       ~28 % mean speedup on the slab benchmark.
   * - PR-INDEX-2
     - ``6cfdfd4``
     - :class:`~orpheus.sn.spatial.sweep_cache.CollisionCache` field
       layout flipped to ``(N, ng, nx)`` natively (``a_attenuation``,
       ``inverse_denom``, etc.); cumprod axis updated 1→2; slab
       ``np.swapaxes`` and curvilinear ``.T`` bridges at cache-read
       sites retired.  ``GeometryCoefficients`` untouched (no group
       axis).  New transient bridge at the
       ``CollisionCache.from_geometry`` callers in
       :class:`SNSolver` (PR-INDEX-3 removes).  Mean
       0.149 ms/sweep on the slab benchmark (down from 0.21 ms at
       PR-INDEX-1); benchmark variance tightened from 2× to 1.06×.
   * - PR-INDEX-3
     - ``313f510``
     - :class:`SNSolver` cross-section storage flipped:
       ``sig_t / sig_a / sig_p / chi`` from ``(nx, ny, ng)`` to
       ``(ng, nx, ny)`` via
       ``xs.<field>.T.reshape(ng, nx, ny)`` at ``__init__``.
       Producer :func:`~orpheus.data.macro_xs.assemble_cell_xs`
       **unchanged** (CP no-regression guaranteed by construction).
       PR-INDEX-2 transient bridges removed; new transients added at
       :meth:`FissionOperator.apply` legacy return contract and
       2-D wavefront ``transport_sweep`` body (PR-INDEX-4 removes).
       ``np.einsum`` rewrites at every reduction site (named-
       intermediate Pattern 3).  ``__debug__`` cell-flattening
       invariant assertion added at ``__init__``.
   * - PR-INDEX-4
     - ``fa41767``
     - Operator-leaf ``apply`` PUBLIC contracts flipped to principled:
       :meth:`FissionOperator.apply` returns ``(ng, nx, ny)``;
       :meth:`ScatteringOperator.apply` returns ``(N, ng, nx, ny)``;
       :class:`LegendreMomentScattering` consumes/returns the
       principled moment layout
       ``(L+1, 2L+1, ng, nx, ny)``;
       :meth:`DiamondDifference.update_batch` consumes
       ``(ng, ...)`` slices;
       ``_sweep_2d_wavefront`` body principled.  PR-INDEX-3 bridges
       at ``fission.py:175`` and ``sweep.py:127`` retired.  Fourteen
       new ``BRIDGE_*`` named intermediates at 11 :class:`SNSolver`
       consumption sites + 3 ``sweep.py`` entry/exit points
       (PR-INDEX-5 removes).  ``EquationMap`` packed-vector
       traversal **deferred** to PR-INDEX-7
       (see :ref:`sn-index-convention-future-work`).
   * - PR-INDEX-5
     - ``3356cec``
     - Public API flip:
       :class:`SNFixedSourceResult` / :class:`SNResult` storage flipped
       to ``(N, ng, nx, ny)`` / ``(ng, nx, ny)``;
       :func:`~orpheus.sn.sweep.transport_sweep` PUBLIC contract
       principled;
       :func:`~orpheus.sn.solver.solve_sn` /
       :func:`~orpheus.sn.solver.solve_sn_fixed_source` return shapes
       flipped; external-source contract ``(N, ng, nx, ny)``.
       Eleven regression snapshots regenerated under principled
       layout via the load-bearing
       :ref:`Step-1 bit-identity-via-transpose gate
       <sn-index-convention-step1-gate>` (ALL 11 cases PASS at
       ``rtol=1e-12``, max abs diff 1.75 × 10\ :sup:`-14`).  ALL
       fourteen PR-INDEX-4 ``BRIDGE_*`` named intermediates retired.
   * - PR-INDEX-6
     - **this PR**
     - Documentation deliverable: this page; cross-references from
       :ref:`theory-discrete-ordinates` and
       :ref:`operator-algebra`; sweep audit of legacy-shape
       mentions in code docstrings + comments + tests.  No
       production code semantics touched.

What stayed deliberately legacy: the FD-matvec internal contract
----------------------------------------------------------------

The FD-matvec internal helpers
:func:`~orpheus.sn.operator.solution_to_angular_flux` and
:func:`~orpheus.sn.operator.transport_operator_matvec_cylindrical`
(plus Cartesian / spherical analogues) carry a separate internal
packed-vector convention: the flat layout
:math:`\texttt{flux}[g + n_g \cdot k_{\text{eq}}]` where
:math:`k_{\text{eq}}` enumerates cells in the order
``for iy: for ix: for n:`` so :math:`n` is the next-fastest axis
after :math:`g`.  The unpacked helper returns
``fi.shape == (ng, N, nx, ny)`` --- a *third* layout, distinct from
both the legacy ``(N, nx, ny, ng)`` and the principled
``(N, ng, nx, ny)``.

Flipping the FD-matvec internal layout would touch 200+ lines of code
across :mod:`orpheus.sn.operator` (30+ ``fi[:, n, i, j]`` indexing
sites) plus the two ``np.transpose`` axis-swap adapters at
``solver.py:1361`` and ``solver.py:1408``.  The migration plan
defers this to **PR-INDEX-7** for two reasons:

1. The packed vector is **not user-facing**.  Callers see the public
   ``angular_flux`` and ``scalar_flux`` arrays (now principled) ---
   the packed vector lives between the Krylov solver and the matvec
   primitive.  Its convention is an implementation detail of the
   FD-matvec path.
2. Bit-identity at the public boundary is preserved by two
   *zero-copy* :func:`numpy.transpose` adapter calls at the Krylov
   decode sites (``solver.py:1361`` and ``solver.py:1408``).  A
   transpose of leading axes is a stride-only view, NOT a memory
   copy --- the runtime cost is one numpy header allocation per
   GMRES iteration, well below measurement noise.

See :ref:`sn-index-convention-future-work` for the PR-INDEX-7 scope.


.. _sn-index-convention-step1-gate:

The load-bearing Step-1 bit-identity gate
=========================================

The migration's most consequential single step was the
**bit-identity-via-transpose verification** that ran at PR-INDEX-5
*before* the regression snapshots were regenerated.  Without this
gate, the principled-layout snapshots would have been written to disk
without independent verification that they corresponded to the same
flux distribution as the legacy snapshots --- and any subsequent
regression test would have been measuring agreement with a
potentially-wrong reference.

The gate is the following Python:

.. code-block:: python

   for case in CASES:
       snap_file = SNAPSHOT_DIR / f'{case.name}.npz'
       old = np.load(snap_file)
       old_sf = np.asarray(old['scalar_flux'], dtype=np.float64)
       cfg = case.builder()
       result = run_case(cfg)
       new_sf = np.asarray(result.scalar_flux, dtype=np.float64)
       # OLD layout (nx, ny, ng); NEW layout (ng, nx, ny); transpose-check:
       new_sf_legacy = new_sf.transpose(1, 2, 0)
       np.testing.assert_allclose(
           old_sf, new_sf_legacy, rtol=1e-12, atol=1e-13, equal_nan=True,
       )

Every one of the 11 regression cases passed; the maximum absolute
difference observed was **1.75 × 10\ :sup:`-14`**, which is the
FP-non-associativity ULP scale predicted for layout flips by the
project's ``vv-principles`` skill (the
"bit-identity vs principled-equivalence" boundary).  Eigenvalue
agreement was at machine precision (max ``keff`` delta
6.66 × 10\ :sup:`-16` --- one ULP for ``keff ≈ 1``).

Only after every case passed did the migration proceed to step 2
(snapshot regeneration via
``tests.sn.regression._generate_snapshots``).  This sequence ---
**verify first, then regenerate** --- is the canonical pattern for
any future layout flip and is enshrined in the migration plan's
risk register.


.. _sn-index-convention-numerical-evidence:

Numerical evidence
==================

The migration's correctness rests on the following gates, each run
post-merge of the relevant PR.  All numbers below are verbatim from
the PR closeout memos
(``.claude/agent-memory/method-implementer/issue_196_pr_index_*_closeout.md``).

Regression snapshots (rtol=1e-12)
---------------------------------

The 11 ``tests/sn/regression/`` snapshots cover:

- Slab 2-group homogeneous DD (``slab_2g_homogeneous_dd_n20``).
- Slab 2-group 3-region DD (``slab_2g_3reg_dd_n40``).
- Sphere 2-group homogeneous DD (``sphere_2g_homogeneous_dd_n20``).
- Sphere 2-group 3-region DD (``sphere_2g_3reg_dd_n40``).
- Cylinder 1-group homogeneous LS\ :sub:`4` DD (``cyl_1g_homogeneous_LS4_dd_n20``).
- Cylinder 1-group homogeneous product-quadrature DD (``cyl_1g_homogeneous_product_dd_n20``).
- Cylinder 2-group 3-region LS\ :sub:`4` DD (``cyl_2g_3reg_LS4_dd_n40``).
- Slab 2-group P\ :sub:`1` anisotropic DD (``slab_2g_p1_aniso_dd_n20``).
- Sphere 2-group P\ :sub:`1` anisotropic DD (``sphere_2g_p1_aniso_dd_n20``).
- 2-D Cartesian 1-group LS\ :sub:`4` DD 15×15 (``2d_1g_LS4_dd_15x15``).
- Slab fixed-source DD (``slab_fixed_source_dd_n20``).

Step-1 transpose-check residuals across all 11 cases (PR-INDEX-5):

.. list-table:: Bit-identity-via-transpose residuals
   :header-rows: 1
   :widths: 56 22 22

   * - Case
     - max ``rtol``
     - ``keff`` delta
   * - ``slab_2g_homogeneous_dd_n20``
     - ≤ 1e-12
     - 0.00 × 10\ :sup:`+00`
   * - ``slab_2g_3reg_dd_n40``
     - ≤ 1e-12
     - 4.44 × 10\ :sup:`-16`
   * - ``sphere_2g_homogeneous_dd_n20``
     - ≤ 1e-12
     - 0.00 × 10\ :sup:`+00`
   * - ``sphere_2g_3reg_dd_n40``
     - ≤ 1e-12
     - 2.22 × 10\ :sup:`-16`
   * - ``cyl_1g_homogeneous_LS4_dd_n20``
     - ≤ 1e-12
     - 4.44 × 10\ :sup:`-16`
   * - ``cyl_1g_homogeneous_product_dd_n20``
     - ≤ 1e-12
     - 2.22 × 10\ :sup:`-16`
   * - ``cyl_2g_3reg_LS4_dd_n40``
     - ≤ 1e-12
     - 4.44 × 10\ :sup:`-16`
   * - ``slab_2g_p1_aniso_dd_n20``
     - NaN-bit-identity
     - 0.00 × 10\ :sup:`+00`
   * - ``sphere_2g_p1_aniso_dd_n20``
     - NaN-bit-identity
     - 0.00 × 10\ :sup:`+00`
   * - ``2d_1g_LS4_dd_15x15``
     - ≤ 1e-12
     - 6.66 × 10\ :sup:`-16`
   * - ``slab_fixed_source_dd_n20`` (no ``keff``)
     - ≤ 1e-12
     - n/a

Max absolute difference across the 11 cases:
**1.75 × 10\ :sup:`-14`** --- FP-non-associativity ULP scale.

Post-regeneration, the regression suite passes 11/11 at
``rtol=1e-12`` against the new snapshots (``248.63 s`` wall-clock for
the full regression run on the development host).

L0 streaming-equilibrium curvilinear
------------------------------------

The L0 curvilinear gate at
``tests/sn/spatial/test_streaming_equilibrium_curvilinear.py``
asserts the streaming-equilibrium identity
:math:`\phi = q / \Sigma_t` to machine precision under refinement.
It is the strongest L0 test for the sphere and cylinder sweep, and
the canonical detector for the historical curvilinear bugs
(ERR-004, ERR-025, ERR-026, weight-normalisation slips).

Pre-migration (``c21c2ef`` baseline): 26 passed in 1044 s.
Post-PR-INDEX-4 (``fa41767``): 26 passed in 1044 s.
Post-PR-INDEX-5 (``3356cec``): the Step-1 bit-identity gate on the
six curvilinear regression snapshots (4 sphere + 2 cylinder) is the
strong proxy that the curvilinear sweep math is unchanged --- those
snapshots exercise the same per-cell algebra at the same granularity.

Performance benchmark
---------------------

Slab sweep benchmark on the ``slab_2g_3reg_dd_n40`` configuration
(N=16, ng=2, nx=160):

.. list-table:: Wall-clock per sweep
   :header-rows: 1
   :widths: 36 32 32

   * - Step
     - Mean / sweep
     - Variance (max/min)
   * - Pre-migration baseline (``c21c2ef``)
     - ~0.21 ms
     - 2×
   * - PR-INDEX-1 (slab joint-batch)
     - ~0.21 ms (variance ↓)
     - --
   * - PR-INDEX-2 (cache layout flip)
     - **0.149 ms**
     - **1.06×**
   * - PR-INDEX-5 (public API flip)
     - 0.149 ms
     - 1.06×

The ~28 % mean speedup at PR-INDEX-1 came from the joint-batch
reduction (2 ``ordinate_scan`` calls per sweep replacing N/2).  The
variance tightening at PR-INDEX-2 (from 2× to 1.06×) came from the
cache layout flip eliminating per-cell strided reads.  PR-INDEX-3,
PR-INDEX-4, and PR-INDEX-5 are layout view changes (zero-copy
transposes) --- no measurable wall-clock change.

2-D wavefront equivalence
-------------------------

The 2-D Cartesian octant-equivalence suite at
``tests/sn/test_2d_octant_sweep_equivalence.py`` exercises six
bit-identity cases plus one closed-form L1 anchor.  Post-PR-INDEX-5,
all 7 pass.  The six bit-identity cases agree at ``nulp=64``
(~1.4 × 10\ :sup:`-14`), which is principled-equivalence per
``vv-principles`` --- the principled layout produces the same value
as the legacy layout up to FP-non-associativity at the ULP regime.


Layout-by-array reference table
===============================

Every array a future maintainer encounters in the SN codebase
matches one of these shapes.  The reference table consolidates what
lives in scattered docstrings.

.. list-table:: SN array shapes (post-PR-INDEX-5)
   :header-rows: 1
   :widths: 32 28 40

   * - Array
     - Shape
     - Defined at
   * - :class:`SNSolver`.\ ``sig_t``, ``sig_a``, ``sig_p``, ``chi``
     - ``(ng, nx, ny)``
     - :meth:`SNSolver.__init__`
   * - :class:`SNSolver`.\ ``scalar_flux``
     - ``(ng, nx, ny)``
     - :class:`SNResult`
   * - :class:`SNSolver`.\ ``angular_flux``
     - ``(N, ng, nx, ny)``
     - :class:`SNResult`
   * - :func:`transport_sweep` input ``Q``
     - ``(ng, nx, ny)``
     - :func:`~orpheus.sn.sweep.transport_sweep`
   * - :func:`transport_sweep` input ``sig_t``
     - ``(ng, nx, ny)``
     - :func:`~orpheus.sn.sweep.transport_sweep`
   * - :func:`transport_sweep` input ``Q_aniso``
     - ``(N, ng, nx, ny)``
     - :func:`~orpheus.sn.sweep.transport_sweep`
   * - :func:`transport_sweep` return ``angular_flux``
     - ``(N, ng, nx, ny)``
     - :func:`~orpheus.sn.sweep.transport_sweep`
   * - :func:`transport_sweep` return ``scalar_flux``
     - ``(ng, nx, ny)``
     - :func:`~orpheus.sn.sweep.transport_sweep`
   * - :func:`solve_sn_fixed_source` input ``external_source``
     - ``(N, ng, nx, ny)``
     - :func:`~orpheus.sn.solver.solve_sn_fixed_source`
   * - :class:`CollisionCache` fields (``a_attenuation``,
       ``inverse_denom``, …)
     - ``(N, ng, nx)``
     - :class:`~orpheus.sn.spatial.sweep_cache.CollisionCache`
       (1-D; collapses ``ny=1``)
   * - :class:`ScatteringOperator`.\ ``apply`` in/out
     - ``(N, ng, nx, ny)``
     - :meth:`~orpheus.sn.scattering.ScatteringOperator.apply`
   * - :class:`FissionOperator`.\ ``apply`` in/out
     - ``(ng, nx, ny)``
     - :meth:`~orpheus.sn.fission.FissionOperator.apply`
   * - :class:`StreamingOperator`.\ ``apply`` in/out (Resolution A)
     - ``(N, ng, nx, ny)``
     - :class:`~orpheus.sn.operator.StreamingOperator`
   * - :class:`CollisionOperator`.\ ``apply`` in/out (Resolution A)
     - ``(N, ng, nx, ny)``
     - :class:`~orpheus.sn.operator.CollisionOperator`
   * - ``LegendreMomentScattering`` moment field
     - ``(L+1, 2L+1, ng, nx, ny)``
     - :mod:`orpheus.sn.scattering`
   * - FD-matvec internal ``fi`` (deferred to PR-INDEX-7)
     - ``(ng, N, nx, ny)``
     - :func:`~orpheus.sn.operator.solution_to_angular_flux`

Two arrays do **not** follow the priority order:

- :class:`CollisionCache` fields drop ``ny`` (they are 1-D-only and
  consume the cell axis as a single innermost contraction).
- The FD-matvec packed-vector internal helper carries ``(ng, N,
  nx, ny)`` --- ``ng`` first, then ``N``, then spatial.  This is the
  PR-INDEX-7 deferred work.

Three array shapes that look like exceptions but are not:

- ``ordinate_scan`` consumes ``(nx, K, ng)`` --- the scan axis
  (cell, here ``nx``) leads.  This is a **primitive contract**
  required by Blelloch's parallel-prefix algorithm (the scan axis
  *must* be the outermost iteration).  The principled-storage
  ``(ng, nx)`` slice is transposed to ``(nx, ng)`` at the call
  site (one ``.T`` per ordinate-chain).
- 1-D internal arrays drop the trailing ``ny`` (slice
  ``[:, :, 0]`` at the boundary).  This is layout-consistent: the
  ``ny=1`` singleton is preserved at the public-API surface, but
  internal 1-D primitives work on ``(ng, nx)`` for clarity.
- The 1-D scratch buffer ``Q_p = Q[:, :, 0]`` in
  :func:`_sweep_1d_unified` is a zero-copy view of the public
  ``(ng, nx, ny)`` source.  No layout decision is made here ---
  it's a slice of the principled storage.


Gotchas and subtleties
======================

ny=1 singleton --- do NOT squeeze
---------------------------------

1-D problems are stored with the trailing ``ny = 1`` axis preserved:
:math:`\psi.\text{shape} == (N, n_g, n_x, 1)`, NOT ``(N, ng, nx)``.

The reason is uniform broadcasting: most SN operations broadcast
across the trailing axis (``ng \cdot V[:, None]``, ``Σ_t \cdot \phi``,
etc.).  With ``ny`` preserved as a singleton, every 1-D operation
uses the same numpy expression as its 2-D counterpart.  Squeezing
``ny`` would force per-dimension branching in every consumer.

**The single exception** is the per-cell cross-section slice in
:func:`_sweep_1d_unified`: ``sig_t_1d = sig_t[:, :, 0]`` strips the
trailing axis for the cache primitive's ``(ng, nx)`` contract.
This is a localised slice, not a layout decision; the result is
re-broadcast to ``ny=1`` at the sweep's public-API exit.

SigS scattering convention --- still ``[g_from, g_to]``
-------------------------------------------------------

The :ref:`scattering-matrix-convention` is unchanged by the migration.
:attr:`Mixture.SigS` matrices are stored as
``SigS[l][g_from, g_to]``; the in-scatter source uses the transpose:
``Q_scatter = SigS^T @ phi``.

The layout migration affects the **storage** of the resulting flux
arrays, not the **convention** of the cross-section matrices.

Per-material vs per-cell cross sections
---------------------------------------

:attr:`Mixture.SigT`, :attr:`Mixture.NuSigF`, etc. are stored per-mixture
as shape ``(ng,)`` (group-only).  The per-cell arrays on
:class:`SNSolver` (``sig_t``, ``sig_a``, ``sig_p``, ``chi``) are
shape ``(ng, nx, ny)`` --- group first, then spatial.  The bridge is
:func:`~orpheus.data.macro_xs.assemble_cell_xs`, which lifts the
per-mixture group-only array to the per-cell flat shape
``(N_cells, ng)``; :meth:`SNSolver.__init__` then transposes and
reshapes to the principled per-cell ``(ng, nx, ny)``.  CP consumes
the producer's flat shape directly --- CP is unaffected by the SN
migration.

Test fixture construction order
-------------------------------

Test fixtures construct sources, cross sections, and fluxes in the
principled order directly:

.. code-block:: python

   Q = np.zeros((ng, nx, ny))                  # principled
   sig_t = np.full((ng, nx, ny), sigma_t)      # principled
   psi = np.empty((N, ng, nx, ny))             # principled
   external = rng.standard_normal((N, ng, nx, ny))

No test should construct in legacy order and then transpose.  The two
remaining transposes in
:file:`tests/sn/test_2d_octant_sweep_equivalence.py` (cases 4--5) are
documented adapters that build sources via a broadcast against
``np.array([...])[None, None, :]`` (per-group profile times spatial
profile) and then transpose to principled --- this is a readability
choice, not a layout slip.


.. _sn-index-convention-future-work:

Future work
===========

PR-INDEX-7 --- EquationMap packed-vector traversal flip
-------------------------------------------------------

The FD-matvec internal helpers
(:func:`~orpheus.sn.operator.solution_to_angular_flux` and the
``transport_operator_matvec_*`` family) carry a
Fortran-flatten layout for the packed vector
:math:`\texttt{solution}[g + n_g \cdot k_{\text{eq}}]` with
:math:`k_{\text{eq}}` enumerating cells via
``for iy: for ix: for n:``.  The unpacked result is
``fi.shape == (ng, N, nx, ny)``.

PR-INDEX-7 will flip this to the principled
``(N, ng, nx, ny)`` traversal:

- Reverse the ``EquationMap`` enumeration order to put :math:`n`
  outermost.
- Update the Fortran-reshape pair at every
  ``solution.reshape(ng, n_eq, order='F')`` /
  ``lhs.ravel(order='F')`` site
  (``operator.py:292``, ``operator.py:460``, etc.).
- Update the 30+ ``fi[:, n, i, j]`` indexing sites in
  ``transport_operator_matvec_*``.
- Retire the two ``np.transpose(fi, (1, 0, 2, 3))`` axis-swap
  adapters at ``solver.py:1361`` and ``solver.py:1408``.

Estimated size: ~200 LoC, all internal, no public API surface
affected.  PR-INDEX-7 is **not** a blocker for the typed-field
contract resume; the two efforts are independent.

Typed-field contract resume
---------------------------

After PR-INDEX-6 lands, the typed-field contract plan (see the memo
at
``.claude/agent-memory/explorer/typed_field_contracts_for_phase_g.md``
and the public API reference at :doc:`../api/discrete_ordinates`)
resumes, with the layout corrected.  The dataclasses become:

.. code-block:: python

   @dataclass(frozen=True, slots=True)
   class AngularFlux:
       values: np.ndarray   # (N, ng, nx, ny) — principled
       sn_mesh: "SNMesh"    # by-reference

   @dataclass(frozen=True, slots=True)
   class ScalarFlux:
       values: np.ndarray   # (ng, nx, ny) — principled
       sn_mesh: "SNMesh"

The dataclasses land on the principled foundation laid by
PR-INDEX-5.  Every operator-leaf's ``apply`` signature becomes
``apply(psi: AngularFlux) -> AngularFlux``, with the four-operator
algebra :math:`(L + C - S - F/k).\texttt{apply}(\psi)` distributing
through :class:`~orpheus.numerics.operator.OperatorSum` unchanged.
This closes the Issue #197 Wave 1 partial as documented in the memo's
§6.

Joint-batch ordinate_scan for curvilinear
-----------------------------------------

The M--M angular thread in curvilinear sweeps is sequential across
ordinates within a :math:`\mu`-level, which forces the curvilinear
sweep to keep a per-ordinate scan.  A research-level reformulation
of the M--M recurrence as a parallel-prefix scan over ordinates
would unlock joint-batching for curvilinear at the same scale
PR-INDEX-1 unlocked for slab.  The estimated win is 3--10× sweep
speedup on cylindrical pin-cell problems.  See the migration plan's
§7 deferred-work register.

JAX / GPU port
--------------

Under the principled layout, the ``(N, ng)`` leading batch maps
cleanly to a GPU grid dimension or to ``jax.vmap(scan, axes=(0,
0))``.  The cell axes :math:`(n_x, n_y)` map to the block dimension.
The migration plan's §7 lists this as the natural follow-up to the
typed-field contract resume.


Cross-references
================

- :ref:`theory-discrete-ordinates` --- SN method theory page; the
  Key Facts header references this convention.
- :ref:`operator-algebra` --- four-operator algebra; every leaf's
  ``apply`` consumes / returns arrays in the convention defined
  here.
- :ref:`scattering-matrix-convention` --- the cross-section matrix
  convention, unchanged by the migration.
- :ref:`synthetic-xs-library` --- the verification cross sections
  used by the regression snapshots.
- ``.claude/plans/principled_index_migration.md`` --- the migration
  plan with the per-PR scope and the deferred-work register.
- ``.claude/agent-memory/method-implementer/issue_196_pr_index_*_closeout.md``
  --- the per-PR closeout memos with verbatim test paste-back.


References
==========

The layout derivation is grounded in the standard SN textbook
treatment of the within-group source iteration.

- [LewisMiller1984]_ §4.5 ("Source Iteration") --- block-diagonality
  of the within-group system.
- [AdamsLarsen2002]_ §III --- the SAILOR / Larsen-Adams
  preconditioned-Krylov framework that motivates the
  ``(N, ng)`` joint-batch storage.
- [Bailey2009]_ Eq. 50 (curvilinear :math:`\alpha` recursion) and
  Eq. 74 (Morel--Montry weights) --- the curvilinear M--M angular
  thread that obstructs joint-batch over ordinates and motivates
  the principled :math:`n` leading layout.
