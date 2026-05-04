.. _theory-reference-solvers:

==========================================================================
Reference Solvers — Pillar-2 continuous truth set
==========================================================================

.. contents:: Contents
   :local:
   :depth: 2


Key Facts
=========

**Read this page first when working on any continuous reference
solver.** Reference solvers are NOT production solvers. They live in
:mod:`orpheus.derivations.continuous` and exist to supply analytical
or semi-analytical truth values that the verification suite (and
therefore the production solvers) trust.

- **Pillar classification** (per the
  ``vv-principles`` skill — closed-form / MMS / semi-analytical):
  every reference is exactly one of these three kinds. Mixing the
  three breaks the verification chain. See
  :ref:`reference-solver-pillar-table` below for the per-solver
  classification.
- **Structural independence** is the load-bearing property:
  a reference solver's L1 cross-check against a production solver
  is meaningful **only** if the two paths exercise different
  integrands or different identities. Two derivations sharing an
  upstream identity will agree to machine precision while both
  being wrong (ERR-032 in the project's error catalog).
- **Folder names follow method-canonical naming**, not author names
  (refactor commit ``d7fa25b``). The folder ``trajectory_resolvent/``
  documents what the method *does* (trajectory tracking + resolvent
  closure), not who introduced it. Author names are reserved for
  *case registries* like ``sood_registry/``.
- **Empty reserve folders carry intent**: the seven folders
  ``spectral_resolvent/``, ``pn_method/``, ``spn_method/``,
  ``escape_probability/``, ``bn_method/``, ``galerkin_sn_hybrid/``,
  ``spectral_collocation/`` exist as placeholder names with READMEs
  that pin the canonical references. Stub theory pages exist with
  ``:label:`` anchors so future code lifts cleanly into the existing
  cross-reference graph.


.. _orbit-space-m-g-classification:

Orbit-space M/G classification — the structural signature
==========================================================

Several reference-solver families partition by **topology, not by
shape** — most prominently :ref:`theory-trajectory-resolvent` and
:ref:`theory-peierls-nystrom`. The precise mathematical signature
that organises this partition is the **orbit-space M/G
classification**.

Given a manifold :math:`M` (the physical 3-D problem domain) and a
symmetry group :math:`G` acting on :math:`M` by isometries, the
**orbit space** :math:`M/G` is the quotient — equivalence classes
of points related by a group element. ORPHEUS reference solvers
exploit the symmetry by reducing the transport problem from
:math:`M` to its 1-D orbit space.

.. list-table:: Orbit-space M/G reductions for ORPHEUS reference solvers
   :header-rows: 1
   :widths: 18 25 18 25 14

   * - Geometry (M)
     - Symmetry group G
     - M/G dimension
     - M/G shape
     - Endpoints
   * - Slab (infinite in :math:`y, z`)
     - :math:`\mathbb{R}^2` translation
     - 1
     - Interval :math:`[-L/2, L/2]`
     - 2
   * - Slab asymmetric
     - :math:`\mathbb{R}^2` translation (BCs distinguish endpoints)
     - 1
     - Interval :math:`[0, L]`
     - 2 (distinct BCs)
   * - Sphere (solid)
     - :math:`SO(3)`
     - 1
     - Interval :math:`[0, R]`
     - 1 (origin + outer)
   * - Hollow sphere
     - :math:`SO(3)` (with inner cut)
     - 1
     - Interval :math:`[R_{\rm in}, R_{\rm out}]`
     - 2
   * - Cylinder (solid, infinite axis)
     - :math:`\mathbb{R}` translation × :math:`SO(2)`
     - 1
     - Interval :math:`[0, R]`
     - 1 (axis + outer)
   * - Annulus (hollow cylinder, infinite axis)
     - :math:`\mathbb{R}` translation × :math:`SO(2)` (with inner cut)
     - 1
     - Interval :math:`[R_{\rm in}, R_{\rm out}]`
     - 2

The **number of distinct boundary endpoints** of M/G is the structural
invariant that determines the boundary-closure rank:

- **One-surface compact** (sphere, cylinder solid) — the M/G interval
  has one outer endpoint (the inner endpoint at :math:`r = 0` is a
  *coordinate* singularity, not a *physical* surface; nothing reflects
  off the origin). Boundary-closure rank is **1**: a single scalar
  resolvent :math:`T = (1 - \alpha\,e^{-\tau_{\rm period}})^{-1}` closes
  the multi-bounce sum.
- **Two-surface** (slab asymmetric, hollow sphere, annulus) — the
  M/G interval has two physical endpoints, each carrying its own BC.
  Boundary-closure rank is **2**: the resolvent
  :math:`T = (I - S)^{-1}` is a :math:`2 \times 2` block on the
  surface state vector.

This is **why** :mod:`orpheus.derivations.continuous.trajectory_resolvent.variant_alpha_core`
treats the slab asymmetric, hollow sphere, and annulus as a single
**rank-2 family** despite their entirely different curvatures:
their orbit-space M/G classification is identical (1-D interval with
two distinct BC endpoints). Only the chord algebra and the
G-equivariant lift back to :math:`M` (the higher-dim Jacobian) differ
between the three.

Where this page (and downstream pages) say "topology" or
"topological class", read **orbit-space M/G classification**. The
``CurvilinearGeometry.topology`` property's two values
(``"one_surface_compact"``, ``"two_surface"``) name the two
endpoint counts that ORPHEUS currently supports. Future extensions
to more endpoints (rank-:math:`N` for :math:`N \ge 3`, e.g.
axially-capped cylinder) would correspond to richer orbit-space
boundary structure (1-D interval with three or more endpoints
under :math:`SO(2) \times` reflection-axis-truncation, etc.).

The orbit-space framing also clarifies what the rank-N falsifications
documented at :ref:`peierls-rank-n-per-face-closeout` and
:ref:`peierls-rank-n-class-b-mr-mg-falsification` are saying: the
rank-:math:`N \ge 2` Marshak / DP\ :sub:`N` closures attempted on
**one-surface-compact** orbit spaces (Class B: solid cylinder /
sphere) failed because there is only one physical endpoint to host
the higher-rank moments, not because the higher-order Legendre
expansion is structurally wrong on multi-endpoint geometries
(Class A: slab, hollow cyl/sph, annulus, where multi-mode per-face
closures are fine).

.. _reference-solvers-three-meanings:

The three meanings of "Green's function"
========================================

This taxonomy is **load-bearing for the entire reference-solver
section** and prevents the most common conceptual mistake when
reading the literature: three different mathematical objects all get
called "the Green's function" in transport papers, and they belong
to three structurally independent families. Mixing them corrupts
the V&V chain.

.. list-table::
   :header-rows: 1
   :widths: 18 28 28 26

   * - Meaning
     - What it constructs
     - How it constructs it
     - Folder
   * - **(α) Trajectory resolvent**
     - Scalar Green's kernel
       :math:`G(\rho \to \rho')` for the medium.
     - Trace characteristic rays, integrate optical depth
       :math:`\tau` along each ray, close multi-bounce trajectories
       via the resolvent :math:`T = (I - S)^{-1}`. Sanchez 2002 family
       — Variant α specialisation in ORPHEUS.
     - :ref:`theory-trajectory-resolvent`
       (``trajectory_resolvent/``)
   * - **(β) Spectral resolvent**
     - Same scalar Green's kernel
       :math:`G(\rho \to \rho')`.
     - Closed-form spectral μ-integration of the within-medium
       angular Green's function. Sanchez 1986 Eq. A1 / PS-1982
       Eq. 21 closed-form. Combinations of :math:`E_n`, :math:`Ki_n`,
       X-function residues.
     - ``spectral_resolvent/`` (reserved; stub at
       :ref:`theory-spectral-resolvent`)
   * - **(γ) Singular-eigenfunction angular Green's**
     - Angular Green's function
       :math:`G(\tau, \tau'; \mu, \mu')` directly.
     - Case ν-spectrum + X-function half-range completeness;
       reconstruct :math:`\psi(r,\mu)` via singular-eigenfunction
       expansion (KLL 1974 Fredholm iteration for slab + sphere).
     - :ref:`theory-singular-eigenfunction`
       (``singular_eigenfunction/``)

**Triple-match correctness ladder.** When all three constructions
exist for the same problem, the canonical L1 verification claim is
agreement at all three points: (α) ≈ (β) ≈ (γ). The three paths
exercise structurally distinct integrands (ray-traced vs spectral-μ
vs ν-spectrum), so triple agreement is L1-grade evidence per the
structural-independence rule in the project's ``vv-principles``
skill. The verification matrix in :doc:`/verification/matrix`
indicates which problems currently realise the triple match.

.. note::

   In ORPHEUS today: (α) is in production for slab / cylinder /
   sphere / annulus / hollow sphere via ``trajectory_resolvent/``.
   (γ) is in production for slab / sphere / cylinder via
   ``singular_eigenfunction/`` (criticality only) and
   :mod:`orpheus.derivations.continuous.fn_method` (interior flux
   reconstruction via KLL 1974). (β) is the **gap** — its closed-form
   spectral kernel is currently obtained indirectly through
   ``trajectory_resolvent/`` rather than via the direct PS-1982
   Eq. (21) / Sanchez 1986 Eq. (A6) evaluator. See the literature
   memo at ``.claude/scratch/sanchez_chandrasekhar_gap.md`` for the
   full implementation roadmap.


.. _reference-solver-pillar-table:

Pillar classification
=====================

Every reference solver belongs to exactly one of the three
verification pillars described in the ``vv-principles`` skill. This
table is the canonical assignment.

.. list-table::
   :header-rows: 1
   :widths: 25 18 22 35

   * - Reference solver
     - Pillar
     - Operator form
     - Notes
   * - :ref:`theory-peierls-nystrom` (Peierls Nyström)
     - Semi-analytical
     - Integral
       (:math:`\phi(\rho) = (c/2)\int E_1(|\rho-\rho'|)\phi(\rho')\,d\rho'`)
     - SymPy-derived kernel + Atkinson product Nyström quadrature
       (ERR-036). Tanh substitution for log-singular kernel
       (ERR-037).
   * - :ref:`theory-trajectory-resolvent` (Variant α MoC)
     - Semi-analytical
     - Integral (Peierls form via trajectory tracking)
     - Bouncing characteristics + multi-bounce resolvent
       :math:`T = (I-S)^{-1}` (Sanchez 1986 Eq. A4 / PS-1982 Eq. 14).
   * - :ref:`theory-fn-method` (F_N method)
     - Semi-analytical
     - Differential transport (boundary-collocated)
     - Siewert-Benoist 1979 / Grandjean-Siewert 1979 / Siewert-Thomas
       1986. Slab + sphere; cylinder via :ref:`theory-singular-eigenfunction`.
   * - :ref:`theory-singular-eigenfunction` (Case 1960 family)
     - Closed-form (criticality) /
       semi-analytical (interior flux)
     - Differential transport (singular eigenfunctions)
     - Atalay 1997 anisotropic slab + sphere via parity flip.
       Westfall-Metcalf 1973 cylinder (isotropic). KLL 1974 Fredholm
       iteration.
   * - :ref:`theory-galerkin-spectral` (Galerkin spectral)
     - Closed-form (matrix eigenvalue)
     - Integral (Carlvik) via Legendre-Galerkin projection
     - Dahl-Sjöstrand 1979 with Carlvik 1968 recurrences as named
       primitives.
   * - :ref:`theory-sood-registry` (Sood-family case registry)
     - n/a (registry, not a method)
     - n/a
     - LA-13511 / LA-UR-03-1987 / Atalay 1997 / KLL 1974 truth
       values. Decouples cases from methods.
   * - :ref:`theory-spectral-resolvent` (reserved)
     - Semi-analytical
     - Integral (Peierls) via spectral μ-integration
     - Stub. Sanchez 1986 Eq. (A6) / PS-1982 Eq. (21) closed-form.
   * - :ref:`theory-pn-method` (reserved)
     - Closed-form (matrix eigenvalue) / Semi-analytical
     - Differential transport (P_N expansion)
     - Stub. Garcia 2017 / 2019 / 2021 multi-region sphere.
   * - :ref:`theory-spn-method` (reserved)
     - Semi-analytical
     - Asymptotic-reduced P_N
     - Stub. Makine 2018; Brantley-Larsen 2000.
   * - :ref:`theory-escape-probability` (reserved)
     - Semi-analytical
     - Integral (escape probability)
     - Stub. Carlvik 1965 / 1967; Benoist 1981.
   * - :ref:`theory-bn-method` (reserved)
     - Semi-analytical
     - Differential transport (boundary collocation)
     - Stub (lower priority). Brockmann 1981.
   * - :ref:`theory-galerkin-sn-hybrid` (reserved)
     - Closed-form / Semi-analytical
     - Differential transport (Galerkin angular projection)
     - Stub (lower priority, research-grade). Morel 1989.
   * - :ref:`theory-spectral-collocation` (reserved)
     - Semi-analytical
     - Spectral collocation
     - Stub (no concrete reference yet).


Pages
=====

Production reference solvers
----------------------------

.. toctree::
   :maxdepth: 1

   peierls
   peierls_nystrom
   trajectory_resolvent
   fn_method
   singular_eigenfunction
   galerkin_spectral
   sood_registry

Reserved reference solvers (stubs)
-----------------------------------

These pages are deliberate placeholders. The Python folders exist
under :mod:`orpheus.derivations.continuous` with READMEs pinning the
canonical references. The stubs reserve the ``:label:`` anchors so
that downstream documentation can already link to the equation
positions while implementation is queued.

.. toctree::
   :maxdepth: 1

   spectral_resolvent
   pn_method
   spn_method
   escape_probability
   bn_method
   galerkin_sn_hybrid
   spectral_collocation


Cross-references
================

- :ref:`theory-transport-methods` — production (discrete) solvers
  whose verification consumes the references catalogued here.
- :doc:`/verification/index` — auto-generated V&V matrix, equation
  coverage, and L0 error catalog.
- :doc:`/verification/reference_solutions` — vocabulary discipline,
  operator-form taxonomy, kernel primitives.
- :ref:`theory-verification` — verification suite overview and
  case grid.
