.. _sn-solver-operator-algebra-coordinator:

SNSolver as an operator-algebra coordinator
============================================

This chapter is where the book's machinery converges: the operators
the preceding chapters built — the swept :math:`(L+C)`, the scattering
:math:`S`, the boundary gain :math:`B`, the fission :math:`F` —
coordinated into the production eigenvalue solve, and, from the
converged flux, the frame projections (homogenisation, condensation)
that hand a coarse problem back to the same solver.

At construction time :class:`~orpheus.sn.solver.SNSolver` caches exactly
**two** operators — the two that are cross-section read-through, and
therefore survive a rebind untouched:

* :attr:`SNSolver.scattering_op` —
  :class:`~orpheus.transport.operators.scattering.ScatteringOperator`
  carrying the P0 + (n,2n) + Pℓ Galerkin reconstruction (Wave D
  Issue 13).
* :attr:`SNSolver.fission_op` —
  :class:`~orpheus.transport.operators.fission.FissionOperator`
  carrying the rank-1-in-energy fission emission (Wave D Issue 13).

The loss composite :math:`L+C` is deliberately **not** cached on the
solver.  Its one spelling is
:func:`~orpheus.sn.coupled_system.build_streaming_collision`, reached
through :func:`~orpheus.sn.coupled_system.build_within_group_system`,
which builds the composite it actually inverts — so a solver-held second
copy would be a twin free to drift from the operand the sweep uses.

.. note::

   Before 2026-07-28 the solver also exposed an ``(L, S, F)`` "operator
   triple" (``SNSolver.L`` / ``.S`` / ``.F``).  It was **retired**: the
   attributes had no production reader, and ``SNSolver.L`` was a
   misnomer — it held the *composite* :math:`L+C`, whereas :math:`L`
   throughout this book is the :math:`\sigma`-free streaming leaf
   (:ref:`the affine collision split <operator-algebra>`).  Consumers
   needing the composite call
   :func:`~orpheus.sn.coupled_system.build_streaming_collision`
   directly.

Each of the cached operators is a
:class:`~orpheus.numerics.operator.LinearOperator`
in the Wave A operator-algebra sense: predicate-typed, composable
under :class:`~orpheus.numerics.operator.OperatorSum` and
:class:`~orpheus.numerics.operator.OperatorProduct`, and protocol-
conforming so the iteration primitives in
:mod:`orpheus.numerics.iteration` consume them without SN-specific
plumbing.  The within-group inner solve is built once from a single
source of truth — the :func:`~orpheus.sn.coupled_system.build_within_group_system`
builder assembles the :class:`~orpheus.sn.coupled_system.WithinGroupSystem`
record, the honest within-group decomposition :math:`(L+C,\ S,\ B)` as a
named splitting :math:`A = M - N`: the invertible resolvent
:math:`M = (L+C)` plus its two lagged coupling gains :math:`N = (S,\ B_a)`
(the bulk scattering :math:`S` and the trace boundary reflection :math:`B`;
zero within-group fission), handed to the **variadic** driver
:math:`\text{Driver}(L_{\rm resolvent},\,*\text{gains})` (Wave O step
O.2a — the transitional :math:`S + B` fold is retired; see
:ref:`bc-extraction-variadic-driver` in :doc:`/theory/foundations/boundary_conditions`).
:func:`_within_group_krylov` wraps the matching
:class:`~orpheus.numerics.iteration.KrylovAcceleration` — and the
decomposition is shared verbatim across the eigenvalue source-iteration
inner (:meth:`SNSolver._solve_source_iteration`), the eigenvalue Krylov
inner (:meth:`SNSolver._solve_krylov`), and both fixed-source paths.

.. admonition:: Key Facts
   :class: tip

   * The within-group system is built ONCE, from a single source of
     truth (:func:`~orpheus.sn.coupled_system.build_within_group_system`):
     the named splitting :math:`A = M - N` — the invertible resolvent
     :math:`M = (L+C)` plus the lagged gains :math:`N = (S,\ B)`.
     Fission is never inside the swept operator; it enters as the
     :math:`1/k`-scaled outer source.
   * Both inner paths — source iteration and Krylov — consume that SAME
     decomposition over the SAME one-walk discretization (matvec ≡
     :term:`sweep`, #206 Phase C): the same fixed point, different rate and
     memory.
   * :meth:`~orpheus.sn.solver.SNSolver.compute_keff` reports **fission
     production over net removal** (:eq:`sn-keff-update`) — the
     eigenvalue of the map the inner solve actually poses (#291, #259).
     Leakage reads the typed boundary trace; reflective faces are a
     structural zero, so lattice anchors hold to the ULP.
   * Homogenisation and condensation are the frame page's
     Petrov-Galerkin consumers; the SN layer only orchestrates:
     :meth:`Solution.homogenize <orpheus.sn.solution.Solution.homogenize>`
     returns a mesh-coupled ``MaterialMesh``;
     :meth:`Solution.condense <orpheus.sn.solution.Solution.condense>`
     returns a portable ``dict[int, Mixture]`` — the
     condense/homogenize asymmetry law
     (:ref:`sn-condense-homogenize-asymmetry`).


The within-group inner solve consumes the primitives directly
-------------------------------------------------------------

:class:`SNSolver`'s within-group inner solve **is** the
:class:`~orpheus.numerics.iteration.SourceIteration` /
:class:`~orpheus.numerics.iteration.KrylovAcceleration` primitive — not
a verbatim replica of its loop.
:meth:`SNSolver._solve_source_iteration` constructs a
:class:`SourceIteration` from the :func:`~orpheus.sn.coupled_system.build_within_group_system` SSoT and
runs it; :meth:`SNSolver._solve_krylov` constructs a
:class:`KrylovAcceleration` from :func:`_within_group_krylov` and runs
that.  The Layer-3 resolvent of the SN row in the
:ref:`eigenvalue-posing` architecture is exactly these primitive
instances.

The primitive is **type-agnostic and angular-capable**: it operates on
the typed :class:`~orpheus.transport.timed_full_field.TimedFullField`
composite, which carries the full :term:`angular flux` on its bulk.  Pℓ
anisotropic scattering therefore rides the angular bulk with no special
plumbing — :meth:`ScatteringOperator.apply` on the timeless
:class:`~orpheus.transport.full_field.FullField` operator carrier (the
driver's :class:`~orpheus.transport.timed_full_field.TimedFullField` iterate
reaches it via MRO) reads the angular moments off the composite and builds
the anisotropic source via :meth:`ScatteringOperator.build_aniso_source`,
all inside the primitive's normal RHS path.  There is **no scalar-flux
limitation** and **no pending "Approach A" cleanup**: the earlier
framing — that :class:`SourceIteration` carried only :term:`scalar flux` and SN
had to replicate the loop verbatim until the angular state could be
threaded through — was a property of an interim scalar-only carrier
that the typed composite retired.  The
``.claude/skills/algebra-of-record`` "Branch 2 implements the same
operator algebra" discipline is satisfied: SN is the discretized
Branch-2 consumer of the shared primitive, not a parallel loop.

The (L + C − S − B)·ψ = (1/k)·F·ψ framing at the solver level
-------------------------------------------------------------

Beyond driving the within-group inner solve, the :math:`(L+C,\ S,\ F)`
framing organises the solver's outer API surface:

* :meth:`SNSolver.compute_fission_source` returns
  :math:`F\,\phi/k` — a thin delegator to :meth:`F.apply` with the
  :math:`1/k` outer-loop scaling applied at the solver level.
* :meth:`SNSolver.solve_fixed_source` solves
  :math:`(L+C-S-B)\,\psi = q_{\rm ext}` (with :math:`q_{\rm ext}` the
  fission source built by ``compute_fission_source``).  Two paths:

  * ``inner_solver="source_iteration"`` — sweep-driven fixed-point
    iteration; the resolvent :math:`(L+C)^{-1}` is realised by the
    one-walk WDD sweep.
  * ``inner_solver="krylov"`` — GMRES on the honest within-group
    matvec, with the sweep resolvent as preconditioner — the same
    one-walk discretization either way (matvec ≡ sweep, #206 Phase C).

* :meth:`SNSolver.compute_keff` returns **fission production over net
  removal**, :eq:`sn-keff-update` — the volume-weighted method-layer
  functional :math:`R_{\nu\Sigma_f}(\phi) / (R_{\Sigma_a}(\phi) + L -
  E_{2n}(\phi))`, derived in :ref:`sn-keff-estimator` below.  The
  SN-specific volume weighting lives here (in the typed
  :class:`~orpheus.transport.reaction_rate_functional.IntegratedReactionRate`
  fields) — one honest realization of the same discipline the
  operator-form :meth:`KEigenvalue.compute_keff` spells with the
  measure absorbed into the operators' action.  (Pre-#291 this method
  returned the leakage-blind :math:`\sum F\phi V / \sum \Sigma_a\phi V`
  ratio; see :ref:`sn-keff-estimator` for why that was a
  non-eigenvalue on any vacuum-bounded problem.)

The solver-level :math:`1/k` scaling (in
:meth:`~SNSolver.compute_fission_source`) and the volume-weighted
eigenvalue estimate (in :meth:`~SNSolver.compute_keff`) are exactly the
points where SN's specifics live; the rest of the solver is
operator-algebra coordination over the canonical
:func:`~orpheus.numerics.eigenvalue.power_iteration` boundary.  These
two K-specific hooks are also precisely *why* the Layer-4 loop is not
yet literally K/α-agnostic — relocating the eigenvalue scaling into the
algorithm is the first step of the α-wave (see the honest-scope caveat
in :ref:`eigenvalue-posing`).

The eigenvalue :math:`\keff` is determined by **power iteration**: an
outer loop updates :math:`k` from the net-removal balance
:eq:`sn-keff-update` (fission production over absorption + leakage −
:math:`(n,2n)` emission), with an inner loop that solves the
within-group scattering problem.

.. _sn-keff-estimator:

The reported eigenvalue: fission production over net removal
------------------------------------------------------------

:meth:`~orpheus.sn.solver.SNSolver.compute_keff` reports the eigenvalue
of the problem the inner solve **actually poses**.  This is the SN
symptom (#291) and the MoC/CP/homogeneous root (#259) of a single
discipline: *the reported* :math:`k` *must be the eigenvalue of the
fixed-source map every method scales only fission by* :math:`1/k`
*through* — scattering and the :math:`(n,2n)` emission are plain gains
assembled **inside** :meth:`~orpheus.sn.solver.SNSolver.solve_fixed_source`,
never scaled by :math:`1/k`.  An estimator that disagrees with its own
posed problem converges cleanly and silently to a **non-eigenvalue
ratio**.

.. math::
   :label: sn-keff-update

   k \;=\; \frac{R_{\nu\Sigma_f}(\phi)}
                {R_{\Sigma_a}(\phi) \;+\; L \;-\; E_{2n}(\phi)}

.. (V&V scope note) Governing/definitional identity: the reported k
   IS the eigenvalue of the posed fixed-source map, not a solver
   eigenvalue-correctness claim against an external analytical reference
   (that rests on the multi-group heterogeneous L1/L2 references
   elsewhere on this page). The label is wired to the cross-engine
   consistency gate tests/sn/eigenvalue/test_keff_estimator_gate.py
   (reported k == the converged fixed-point map ratio k* = P(Mφ*)/P(φ*),
   map-ratio ground-truth noise ≤ 2e-11) with in-file mutation teeth.

The three terms are typed volume-integrated reaction-rate functionals
and one boundary functional:

* **Numerator** :math:`R_{\nu\Sigma_f}(\phi) = \int_V \nu\Sigma_f\,\phi\,dV`
  — the fission production, the typed
  :class:`~orpheus.transport.reaction_rate_functional.IntegratedReactionRate`
  over :math:`\nu\Sigma_f` (the :math:`\phi^\dagger\!=\!1` degenerate of
  the homogenization Petrov–Galerkin bilinear).  The :math:`(n,2n)`
  emission is **not** production here — contrast
  :meth:`~orpheus.sn.solver.SNSolver.compute_production_rate`, the
  ERR-052 renormalisation scale anchor, which keeps *total* physical
  production (fission **plus** the :math:`(n,2n)` emission).  The role
  split is the load-bearing #259 correction: the same physical
  :math:`(n,2n)` neutrons are a **scale** quantity for the outer
  renormalisation but a **removal-reduction** in the eigenvalue balance.
* **Absorption** :math:`R_{\Sigma_a}(\phi) = \int_V \Sigma_a\,\phi\,dV`,
  with :math:`\Sigma_a = \Sigma_f + \Sigma_c + \Sigma_L +
  \sum_{g'}\Sigma_{2,g\to g'}` — i.e. ``absorption_xs`` counts the
  :math:`(n,2n)` **collision once** (the neutron is removed from its
  incident group by the collision).  See
  :attr:`~orpheus.data.macro_xs.mixture.Mixture.absorption_xs`.
* **Leakage** :math:`L` — the net vacuum-boundary outflow (below).  On a
  reflective (lattice) problem it is a **structural zero**.
* **Emission** :math:`E_{2n}(\phi) = \int_V \sum_{g,g'} 2\,\Sigma_{2,g'\to
  g}\,\phi_{g'}\,dV` — the :math:`(n,2n)` **emission** (two neutrons out
  per collision; the factor 2).  A gain, so it **reduces** net removal.

The net :math:`(n,2n)` effect on removal is therefore
:math:`\underbrace{\sum_{g'}\Sigma_{2,g\to g'}}_{\text{collision, in }\Sigma_a}
- \underbrace{2\Sigma_2}_{E_{2n}} = -\Sigma_2` — **one extra neutron
gained** per collision, exactly the physics of a neutron-doubling
reaction.

**Balance identity (divergence-telescoping).**  The angle- and
group-summed discrete cell balance for cell :math:`i` in the posed
eigenproblem is

.. math::
   :label: sn-keff-cell-balance

   \underbrace{\sum_{f\in\partial i}\!\bigl(\textstyle\sum_g J_g\bigr)\,\Delta A_f}
              _{\text{net face flow}}
   \;+\; \Sigma_{t,i}\,\phi_i\,V_i
   \;=\; \frac{1}{k}\,R_{\nu\Sigma_f,i}
        \;+\; \Sigma_{s,i}\,\phi_i\,V_i
        \;+\; E_{2n,i}

.. (vv-status rationale) Derivation step (the divergence-telescoping cell
   balance). Its terminal result sn-keff-update is verified by the k* map-ratio
   gate (tests/sn/eigenvalue/test_keff_estimator_gate.py); definitional.
.. vv-status: sn-keff-cell-balance documented

(streaming + total collision on the left; the isotropic fission source
scaled by :math:`1/k`, plus the *unscaled* scatter and :math:`(n,2n)`
gains, on the right).  Summing over **all** cells, every interior face
is shared by two cells with opposite outward normals and equal current
(continuity), so the interior face-flow terms **telescope to zero** —
only the domain-boundary faces survive, and their sum is the net
leakage :math:`L`.  With :math:`\Sigma_t - \Sigma_s = \Sigma_a` this
collapses to

.. math::

   \frac{R_{\nu\Sigma_f}(\phi)}{k} \;=\; R_{\Sigma_a}(\phi) \;+\; L
                                        \;-\; E_{2n}(\phi),

which is :eq:`sn-keff-update` rearranged.  This is the same discrete
divergence discipline the diffusion page states as
:math:`\mathbf 1^{\mathsf T}(C-S)=\Sigma_a` with interior-face
telescoping (see :ref:`diffusion-leakage-boundary-leaves`); SN and
diffusion report the *same* balance-law eigenvalue, differing only in
how the streaming operator is discretised.

The leakage functional
~~~~~~~~~~~~~~~~~~~~~~~~

.. math::
   :label: sn-leakage-functional

   L \;=\; \sum_{f\,\in\,\text{vacuum}} \oint_{f} dA\,
           \sum_g J_g(\mathbf{r}_f)\,,
   \qquad
   J_g \;=\; \sum_m (\Omega_m\cdot\hat n_f)\, w_m\, \psi_{m,g}

is the face-area integral of the boundary trace's **net outward
current**.  The angular-to-scalar reduction :math:`J_g` is
:meth:`AngularBoundaryFlux.net_current
<orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux.net_current>`
— the single source of the :math:`\Omega\cdot\hat n\,w` contraction, the
angular sibling of the scalar trace's :math:`J = J^+ - J^-`
(:meth:`ScalarBoundaryFlux.net_current
<orpheus.transport.fields.scalar_boundary_flux.ScalarBoundaryFlux.net_current>`).
It is spelled through the trace space's own atoms — the signed
projection table
:attr:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.omega_dot_n`
and the :math:`|\Omega\cdot\hat n|\odot w` partial-current metric (using
the identity :math:`\operatorname{sign}(\Omega\cdot\hat n)\cdot
|\Omega\cdot\hat n|\,w = \Omega\cdot\hat n\,w`) — so no consumer
re-derives the cosine weighting.  Tangential :term:`ordinates <ordinate>` carry zero weight
and drop out.

The face measure :math:`dA` is supplied by
:meth:`SNSolver._face_area_of`, matching the cell
``volume_measure`` exactly so the balance identity closes:

.. list-table:: Boundary-face measure by geometry
   :header-rows: 1
   :widths: 30 30 40

   * - Geometry
     - Face measure :math:`\Delta A`
     - Source
   * - 1-D slab
     - :math:`1` (per unit cross-section)
     - :attr:`MaterialMesh.areas <orpheus.transport.mesh.material_mesh.MaterialMesh.areas>`
   * - 1-D cylinder
     - :math:`2\pi R` (per unit height)
     - ``MaterialMesh.areas``
   * - 1-D sphere
     - :math:`4\pi R^2`
     - ``MaterialMesh.areas``
   * - 2-D Cartesian
     - transverse edge-cell widths (unit depth)
     - ``mesh.axes`` transverse extent
   * - 3-D Cartesian
     - :math:`\Delta A_{\mathbf c} = \prod_{j\ne a}\Delta_j[c_j]`
       (transverse-area outer product)
     - ``mesh.axes`` transverse extents

The :math:`d \ge 2` Cartesian arms are ONE generic body: the outer
product of the *other* axes' edge widths in **ascending axis order** —
the same codimension-1 enumeration as
:func:`~orpheus.transport.mesh.axis.face_shape`, so the measure array
broadcasts cell-for-cell against the ``(ng, *face_spatial)`` net
current, and the 2-D width vector is just the single-transverse-axis
degenerate (bit-identical to the pre-3-D spelling).

The 3-D arm originally shipped as a **typed refusal**
(``NotImplementedError``): guessing the transverse product's cell
ordering would silently mis-weight the leakage sum, and Cardinal Rule 1
forbids returning a wrong-but-plausible number.  The wire landed
2026-07-13 when the first 3-D vacuum eigenvalue consumer arrived (the
d=3 Mode-9 G-S≡Jacobi gate), with the ordering pinned twice in
``tests/sn/eigenvalue/test_keff_estimator_gate.py``: an **object-level
pin** (face measure ≡ the boundary layer's ``volumes / Δ_axis``, the
mesh's own ascending-axis enumeration — vv Mode-12 discipline: pin the
object, not only the k functional) and the **k* map-ratio gate** on a
Mode-2 asymmetric all-vacuum box, whose teeth are proven by permanent
in-process mutants — a reversed transverse enumeration moves the
reported k by a measured **13.9 %** against the estimator-independent
:math:`k^*` (clean agreement :math:`6\times10^{-10}`), and a transposed
enumeration crash-REDs on the broadcast.  A trace carrying a
``#251`` transverse face-moment tail is refused loudly at the
consumption site (the face integral must consume ONLY the
transverse-average moment — higher Legendre face moments integrate to
zero over each face cell — and that slot-0 read has no consumer yet).

Reflective faces are a structural zero
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The reflective law equates a face's inflow to its reflected outflow
**exactly**, so the net current there vanishes *by construction*.
:meth:`~orpheus.sn.solver.SNSolver._boundary_leakage_rate` therefore
**skips** reflective faces — it never accumulates them, rather than
accumulating a value that ought to be zero but carries
:math:`\pm`-cancelling angular-sum floating-point noise.

This is a deliberate design choice with a bit-level payoff.  On an
all-reflective (lattice) problem :math:`L` is a structural ``0.0``, and
on a :math:`\Sigma_2`-free mixture :math:`E_{2n}` is exactly ``0.0`` (the
per-material :math:`(n,2n)` loop adds nothing), so
:eq:`sn-keff-update` reduces **bit-identically** to the historical
lattice functional ``production / absorption``.  Every pre-existing
reflective eigenvalue anchor is preserved to the last ULP — the
unification adds terms that vanish structurally, not numerically, on the
cases it must not perturb.

The scale bridge: trace of the last inner solve
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The leakage term reads the **typed** boundary trace of the last inner
solve (``self._psi_typed.boundary``), whereas the numerator/denominator
reaction rates consume the bare-array flux :math:`\phi` the estimator is
handed.  These two representations can be at **different scales**:
:func:`~orpheus.numerics.eigenvalue.power_iteration` renormalises
:math:`\phi` to unit production rate *between* the inner solve and the
:math:`k`-update (ERR-052), so the stored trace belongs to the
**un-renormalised** last iterate while the estimator's :math:`\phi` is
its renormalised multiple.

Leakage is degree-1 homogeneous in :math:`\psi`, so the fix is a single
rescale by the fission-production ratio of the two fluxes
(``self._phi_of_trace``, stored alongside the trace at **both**
inner-path returns) — exactly ``1.0`` when the caller passes the
returned flux itself.  The **contract** is therefore: the flux handed to
:meth:`~orpheus.sn.solver.SNSolver.compute_keff` must be a scalar
multiple of the last inner solve's flux (true for ``power_iteration`` and
for every manual solve-then-estimate loop).

If a vacuum face exists but **no** inner solve has stored a trace,
:meth:`~orpheus.sn.solver.SNSolver._boundary_leakage_rate` raises
``RuntimeError`` — the leakage cannot be answered honestly, and silently
returning it as zero would *reproduce the #291 omission*.  Fail loud;
never return a non-eigenvalue.

The R7 :math:`(n,2n)` convention fork
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The historical spelling put the :math:`(n,2n)` emission in the
**numerator** as production,

.. math::
   :label: sn-keff-old-n2n

   k_{\text{old}} \;=\; \frac{R_{\nu\Sigma_f} + E_{2n}}{R_{\Sigma_a}},

.. (vv-status rationale) Definitional/historical contrast: the superseded
   (n,2n)-in-numerator spelling. No current code implements it; its bias is
   characterised in the pre-fix table (#291 commit d1daaac).
.. vv-status: sn-keff-old-n2n documented

which is a **non-eigenvalue** of the posed map whenever
:math:`\Sigma_2 \neq 0` *and* :math:`k \neq 1`.  The reason is exactly
the posing asymmetry: the inner solve scales **only** fission by
:math:`1/k`; the :math:`(n,2n)` emission is an *unscaled* gain in the
sweep source.  So the eigenvalue of that map is
:math:`k^\star = R_{\nu\Sigma_f}/(R_{\Sigma_a} - E_{2n})` (reflective,
:math:`L=0`), and putting the *unscaled* emission in the numerator does
not recover it.  Writing :math:`f = R_{\nu\Sigma_f}`,
:math:`a = R_{\Sigma_a}`, :math:`e = E_{2n} = 2s_2` and substituting
:math:`f = k^\star (a - e)`:

.. math::
   :label: sn-keff-old-bias

   k_{\text{old}}
   \;=\; \frac{k^\star (a - e) + e}{a}
   \;=\; k^\star \;+\; \frac{2 s_2\,(1 - k^\star)}{a}.

.. (vv-status rationale) Mathematical identity (the derived bias of the
   superseded estimator, k_old - k* = 2 s_2 (1 - k*) / a). Historical
   characterisation; no current implementing code.
.. vv-status: sn-keff-old-bias documented

The two agree only when :math:`s_2 = 0` (no :math:`(n,2n)`) or
:math:`k^\star = 1` (critical).  For a supercritical
:math:`k^\star > 1` the correction is negative
(:math:`k_{\text{old}} < k^\star`); for subcritical, positive.  The MoC
and CP pages carry the same fork
(:eq:`moc-keff-update`, :eq:`cp-keff-update`); CP was the one member
already spelled on net removal.

What was tried and found
~~~~~~~~~~~~~~~~~~~~~~~~~~

The #291 bias was characterised pre-fix (commit ``d1daaac``, Gauss–
Legendre :math:`n=8`, map-ratio ground truth noise :math:`\le 2\times
10^{-11}`) across the five gate configurations:

.. list-table:: Pre-fix reported :math:`k` vs the posed-problem eigenvalue :math:`k^\star`
   :header-rows: 1
   :widths: 40 18 18 24

   * - Configuration
     - Pre-fix reported
     - Posed :math:`k^\star`
     - Bias
   * - homog. 2G vacuum slab (width 8)
     - 1.83767525
     - 0.98163269
     - :math:`+87.2\%` (:math:`L/A = 0.872`)
   * - het. vacuum sphere P\ :sub:`0`
     - 0.86484694
     - 0.70601977
     - :math:`+22.5\%`
   * - het. vacuum sphere P\ :sub:`1`
     - 0.85080423
     - 0.67876772
     - :math:`+25.3\%`
   * - reflective control (:math:`\Sigma_2=0`)
     - 1.87500000
     - 1.87500000
     - :math:`\equiv` (bias :math:`1.2\times 10^{-10}`)
   * - reflective :math:`\Sigma_2\neq 0`
     - 1.92857143
     - 2.61278195
     - :math:`-26.2\%` (R7 defect)

The two failure classes are visible in one table: the vacuum rows are
the **leakage omission** (#291) — the reported :math:`k` overshoots by
the leakage-to-absorption ratio :math:`L/A`; the last row is the **R7
:math:`(n,2n)` convention** — zero leakage, yet a
:math:`-26.2\%` error because the emission was mis-posed as production.
The reflective-control row is exactly the bit-identity guarantee above.
The exact check on the R7 row is
:math:`0.78/(0.5185 - 0.2200) = 2.61278`, and
:math:`(0.78 + 0.2200)/0.5185 = 1.92857` reproduces the old value —
matching :math:`k_{\text{old}} = k^\star + 2s_2(1-k^\star)/a` term for
term.

Post-fix, reported :math:`k` and the map-ratio :math:`k^\star` agree to
:math:`\le 6\times 10^{-10}` on all five.  The P\ :sub:`0`\ –P\ :sub:`1`
sphere gap :math:`\Delta` roughly **doubled** (:math:`1.404\times
10^{-2} \to 2.725\times 10^{-2}`) but stays inside the diagnostic
:math:`(10^{-3}, 5\times 10^{-2})` band — the P\ :sub:`1` anisotropic
correction is now measured against the correct eigenvalue on both
solves.

The V&V decision was a **principled re-baseline** (per ``vv-principles``
bit-identity-vs-principled-equivalence): the old reported :math:`k` was a
*different functional* from the posed problem's eigenvalue, so the new
value is not a regression to be tolerance-matched but a correction to be
verified against a structurally-independent reference (the fixed-point
map ratio).

Verification
~~~~~~~~~~~~

The permanent gate is
``tests/sn/eigenvalue/test_keff_estimator_gate.py``: it asserts the
reported :math:`k` equals the converged fixed-point map ratio
:math:`k^\star = P(M\phi^\star)/P(\phi^\star)` across the four physics
regimes — {vacuum slab, vacuum sphere (pinning the :math:`4\pi R^2`
face-area convention), reflective bitwise-degenerate, reflective
:math:`\Sigma_2\neq 0`} — with **in-file mutation teeth**: a
leakage-drop mutation reds the vacuum legs while staying bitwise-green
on reflective; a leakage sign-flip crash-reds through the scale-bridge
guard; and the old :math:`(n,2n)`-in-numerator convention reds the
:math:`\Sigma_2\neq 0` leg.

This is a **consistency** gate: the map ratio is the structurally-
independent ground truth for "does the estimator return the eigenvalue
of its own posed map", and is blind by construction to *which*
eigenvalue that is.  The SN solver's eigenvalue **correctness** — that
the posed map's eigenvalue is the *physically right* :math:`k` — rests
on the multi-group heterogeneous L1/L2 references in
:doc:`/theory/verification/sn`, not on this gate.

Two Inner Solvers
-----------------

**Source iteration (sweep-based):**

- Operator: :math:`(L+C)^{-1}` (the one-walk WDD sweep)
- Iterate: the typed field composite (angular bulk + boundary trace)
- Fixed-point: :math:`\psi^{(k+1)} = (L+C)^{-1}(S\,\psi^{(k)} +
  B\,\psi^{(k)} + q_{\rm ext})`
- Convergence rate: spectral radius of :math:`(L+C)^{-1}(S+B)` — the
  :term:`scattering ratio` :math:`c` (:doc:`slab_one_group`)
- Cost per iteration: one transport sweep
- Works for all geometries

**Krylov (direct operator):**

- Operator: the honest :math:`(L+C-S-B)` applied matrix-free (its
  :math:`(L+C)` piece via :meth:`StreamingCollisionOperator.apply` — the same
  one-walk discretization the sweep realises; L21 matvec ≡ sweep)
- Iterate: the same typed composite; GMRES additionally stores its
  Krylov basis (``restart`` × the composite's ``n_dof`` — the ERR-053
  sizing family)
- System: :math:`(L+C-S-B)\,\psi = q_{\rm ext}` — scattering and the
  boundary gain live in the operator, not the lagged source;
  :math:`q_{\rm ext}` is the :math:`1/k`-scaled fission source
- Convergence: GMRES with sweep preconditioner, typically ~100
  Krylov iterations at ``tol=1e-4`` (always converges)
- Available for all geometries (Cartesian, spherical, cylindrical)

Wave E Round 2 (Issue #164) replaced the legacy BiCGSTAB FD-operator
path with this Krylov path.  See the Krylov alternative in
:doc:`slab_one_group` for the full discussion.

The two paths share the **one** loss-representation discretization
(matvec ≡ sweep, #206 Phase C), so they converge to the same fixed
point; the Wave-D-era design in which they carried different spatial
closures — and disagreed on coarse-mesh :math:`\keff` — is recorded
in the two-closure history (:ref:`loss-rep-history`).

.. _sn-convergence-contract:

The convergence contract — a best-effort answer says so
--------------------------------------------------------

Every entry above can stop for two structurally different reasons: it
**converged**, or its **budget ran out** and the returned field is a
best-effort iterate, mid-descent.  Both come back as the same type from
the same call, so the distinction has to be carried explicitly or it is
lost.

**Where the fact comes from.**  The loop that stops knows why it stopped,
and since 2026-08-08 it says so rather than discarding it:

* the **outer** (eigenvalue) fact is
  :attr:`~orpheus.numerics.eigenvalue.PowerIterationOutcome.converged`,
  recorded by :func:`~orpheus.numerics.eigenvalue.power_iteration` at the
  ``break``;
* the **inner** (fixed-source) fact is one shared predicate over the
  residual history — strictly below tolerance, the same test the
  within-group certificate uses to decide whether a claim exists at all.

Both surface on :attr:`Solution.history.converged
<orpheus.sn.solution.IterationHistory.converged>`, which is a **required**
field: it has no default, so a producer cannot claim convergence by
omission.

.. warning::

   **A value gate that does not assert convergence is asserting an
   arbitrary iterate.**  This is not hypothetical.  `[M]`
   ``test_d3_pure_absorber_per_ordinate_psi_exact`` asserted a closed-form
   identity to ``rtol=1e-10`` on an all-reflective 3-D box that needs
   **1631** sweeps, against a default ``max_inner`` of **1000**.  It read
   the 999th iterate, never read the flag the solver had honestly set to
   ``False``, and passed for months because the truncated error happened to
   land inside the tolerance — until a *correct* quadrature change (#337)
   moved it out.  The one-line defence is to assert
   ``sol.history.converged`` **before** reading any value.

   The diagnostic tell is worth memorising, because it points the wrong
   way: the error was **bit-identical at every** ``inner_tol`` **from 1e-9
   to 1e-15**, which reads as a discretization floor.  It is the opposite —
   the running residual never fell below even the loosest tolerance, so
   every run hit the same cap and returned the same bytes.
   Tolerance-insensitivity means the tolerance never *bound*; read the
   iteration count against the budget before concluding anything about the
   discretization.

**Loudness.**  A truncated exit emits
:class:`~orpheus.numerics.convergence.ConvergenceWarning`, naming the
budget that ran out, the tolerance missed, and how far the last iterate
was — the distance between "one more sweep" and "diverging" wants opposite
responses.  It is a warning rather than an exception by the ERR-053
precedent (legitimate callers harvest the residual history of a
deliberately-truncated solve), and it escalates to a hard failure with::

    python -O -m pytest -W error::ConvergenceWarning

**Budget sizing.**  ``max_inner`` is a fixed default, and the required
budget is not: `[M]` an all-reflective box needs ~an order of magnitude
more sweeps per added dimension (d=1 **32**, d=2 **258**, d=3 **1631**),
and the cost scales as :math:`\Sigma_t \cdot n_{\rm inner} \approx`
constant.  One vacuum face collapses the d=3 figure to 208 — the expensive
corner is specifically zero-leakage, weakly-absorbing, and 3-D.  Deriving
the default rather than hardcoding it is tracked on `Issue #340
<https://github.com/deOliveira-R/ORPHEUS/issues/340>`_.

Gates: ``tests/sn/solve/test_convergence_contract.py`` — each honesty claim
is a PAIR (a converging configuration and a deliberately-starved one),
because asserting ``converged is True`` on a solve that converges is
satisfied by the very hardcoded ``True`` the contract forbids; it is the
starved leg that has teeth.  `[M]` a 6-mutation battery, positive control
first, reddens as designed — including re-introducing #342 verbatim.


.. _sn-consuming-the-frame:

Consuming the frame in SN
=========================

Spatial homogenisation and energy condensation are **discrete-frame
projections** — the Petrov-Galerkin coefficient extraction
:math:`G^{-1}M` of a flux- (or spectrum-) weighted frame. All of that
theory — rate preservation, the source-group / sink-sum matrix rules, the
metric-fold-vs-bilinear adjoint argument, fractional-overlap re-binning,
the condense/homogenize asymmetry law, and the verification gates — is
the frame page's headline **Petrov-Galerkin** consumer; see
:ref:`sn-spatial-homogenization` and :ref:`sn-energy-condensation`
(:doc:`/theory/foundations/frame`). This section keeps only the **SN-layer orchestration**:
how the SN :class:`~orpheus.sn.solution.Solution` drives that machinery
from a converged flux.

Homogenisation: the solve → homogenize → re-solve loop
------------------------------------------------------

:meth:`Solution.homogenize <orpheus.sn.solution.Solution.homogenize>`
takes a coarse mesh (:class:`~orpheus.geometry.mesh.Mesh1D` or
:class:`~orpheus.geometry.mesh.Mesh2D`) and returns a
:class:`~orpheus.transport.mesh.material_mesh.MaterialMesh` — the coarse
geometry already carrying one freshly-homogenised effective
:class:`~orpheus.data.macro_xs.mixture.Mixture` per coarse cell. The SN
:class:`~orpheus.sn.solution.Solution` owns the converged flux, so the SN
layer is what builds the flux-weighted **test** basis the frame consumes;
the frame itself, and the rate-preservation theory that *forces* the flux
weighting (rather than a plain volume average), live in
:ref:`sn-spatial-homogenization` (:doc:`/theory/foundations/frame`). The returned
``MaterialMesh`` is re-promoted to a solvable phase space by
:meth:`SNMesh.from_material_mesh
<orpheus.sn.mesh.augmented_mesh.SNMesh.from_material_mesh>`, closing the
**solve → homogenize → re-solve** loop. The return type is
**mesh-coupled** (geometry and materials born together) — the space half
of the condense/homogenize asymmetry law
(:ref:`sn-condense-homogenize-asymmetry`, :doc:`/theory/foundations/frame`).

Condensation: per-material representative spectra
-------------------------------------------------

:meth:`Solution.condense <orpheus.sn.solution.Solution.condense>` is the
SN-layer orchestration of energy condensation. It condenses **each
material with its own representative spectrum** — the flux·volume-weighted
flux over the cells where the material appears:

.. math::
   :label: energy-condensation-representative-spectrum

   \varphi^{(m)}_g \;=\;
   \sum_{i:\,\mathrm{mat}(i)=m} V_i\,\phi_{i,g},

.. (vv-status rationale) Representational identity: the per-material
   representative spectrum used as the condense test weight — the
   flux·volume-weighted flux over the material's cells (mirrors how
   ``homogenize`` derives its flux weight). A definition consumed by
   :meth:`Mixture.condense`; the end-to-end rate preservation it feeds is
   the L1 gate, not a separate claim.
.. vv-status: energy-condensation-representative-spectrum documented

used as the test weight in :meth:`Mixture.condense
<orpheus.data.macro_xs.mixture.Mixture.condense>` — the data-layer
collapse verb, whose spectrum-weighted-collapse theory is
:ref:`sn-energy-condensation` (:doc:`/theory/foundations/frame`) — mirroring how
:meth:`Solution.homogenize` derives its flux weight from the same solved
flux. The result is a **portable** ``dict[int, Mixture]`` keyed by
material id — few-group cross sections carrying the coarse ``eg``, not
bound to any mesh (the **mesh-decoupled** half of the asymmetry law,
:ref:`sn-condense-homogenize-asymmetry`, :doc:`/theory/foundations/frame`). A material with
no flux in a fine group contributes zero weight there; the condense
frame's Moore–Penrose Gram handles any empty coarse group.
