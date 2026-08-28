.. _sn-development-history:

Development history
===================

S\ :sub:`N` is ORPHEUS's **architectural prototype**: the typed-field
algebra and the composable loss composite :math:`A = L + C - S - B`
(posed :math:`A\,\psi = \tfrac{1}{k}\,F\,\psi`)
were developed here *first*, and are the standard the other solvers (CP,
MoC, MC, diffusion) inherit as they are built — which is why this page
carries far more development history than any other theory chapter.

This is a reverse-chronological (latest first) changelog of the major
**architectural** milestones. Iteration-rate work, gate counts, and
intermediate replans are deliberately omitted — see the GitHub issues
and the per-phase plan files for that granularity. Each entry names the
architectural change, its issue, and the commit/merge where the work
lives.  Every entry below is **merged to main** — merge-status is a git
question, never a frozen note (the 2026-07-24 repair: 22 entries still
said *"in development"* for campaigns that had merged weeks earlier;
``.claude/rules/process-discipline.md`` "trust git for merge-status" is
the standing rule this page now follows — a new entry lands with its
merge hash or not at all).

.. list-table::
   :header-rows: 1
   :widths: 8 52 12 28

   * - When
     - Architectural milestone
     - Issue
     - Where
   * - 2026-08-24
     - **A field is an element of a SPACE — the mesh binding retires from
       the field layer, construction goes space-primary, and the space's
       own reductions become frame-induced operators** (campaign 1, phase
       CS4b, steps S1–S7).  Five consequences the S\ :sub:`N` solver sees
       directly.  **(1) The carrier mints, the leaves read** (S1/S2):
       :attr:`SNMesh.full_field_space
       <orpheus.sn.mesh.augmented_mesh.SNMesh.full_field_space>`'s
       interior IS the carrier's axis-built angular mint rather than a
       parallel hand-spelled shape, and every field leaf sources its space
       from a cached carrier mint instead of naming one itself.
       **(2) Partnering keys on space CONTENT, never on mesh IDENTITY**
       (S3, the *F2 doctrine*): the partner gate's third tier — object
       identity of the two meshes — retires, so twin carriers and
       BC-only-differing carriers legitimately mix, while a moved cell
       edge, a different group structure or a different quadrature refuse
       exactly as before; trace and ray spaces gain content-digest names,
       which is what makes ``(name, shape)`` equality BE content equality
       (:ref:`cone-fiber-discipline`).  **(3) The binding itself retires**
       (S4/S5): :class:`~orpheus.numerics.field.Field` now carries
       ``values`` and ``space`` and nothing else, and the per-family
       ``_phase_space_shape`` hook collapses into ``Field``'s own
       ``values.shape == space.shape`` check — it was a twin of the
       space's own content, kept alive only by the binding.  The
       mesh-keyed leaf **sugar tier is DELETED** (``zeros_on`` /
       ``from_mesh`` / ``from_ndarray`` are gone from every field leaf),
       the composite allocators take ``space=``, and the carrier grows the
       named replacement for the ``spatial_moments=`` integer:
       :attr:`SNMesh.angular_trial_space
       <orpheus.sn.mesh.augmented_mesh.SNMesh.angular_trial_space>`, the
       scheme-widened angular mint, ``is``-shared at every width.  Two
       survivors are NOT sugar and stay by design —
       :meth:`MaterialXSField.from_mesh
       <orpheus.transport.mesh.material_xs_field.MaterialXSField.from_mesh>`
       (assembly tier) and :meth:`MultiplicationOperator.from_mesh
       <orpheus.transport.operators.multiplication_operator.MultiplicationOperator.from_mesh>`
       (operator tier).  **(4) The angular reduction and the isotropic
       source projection have ONE realization each** (S6): the space mints
       the axis **collapse pair** — :meth:`FunctionSpace.retraction
       <orpheus.numerics.space.FunctionSpace.retraction>` and
       :meth:`~orpheus.numerics.space.FunctionSpace.section` — and both
       are the *induced output* of a single-region indicator frame built
       and discarded at one site (the stage-2 generator discipline), so
       the section's normalising divisor is READ OFF that frame's
       :math:`1 \times 1` Parseval metric instead of being chosen by hand.
       ``integrate_angular`` and ``from_isotropic`` re-key onto the pair,
       and the per-face packing loop re-homes to its layout's own
       :meth:`FaceLayout.pack
       <orpheus.numerics.face_layout.FaceLayout.pack>` (native place).
       Full account: :ref:`spaces-collapse-pair` on
       :doc:`/theory/foundations/spaces`; the field-layer narrative is on
       :doc:`/theory/foundations/field_algebra`.  **(5) The mesh-less
       carrier's two meanings un-weld** (S7): promoting the
       infinite-medium 1-cell carrier to an S\ :sub:`N` phase space now
       raises a typed ``ValueError`` naming the reason — pre-repair a
       messageless bare ``assert`` that the canonical ``python -O``
       runner **strips**, leaving a deep ``AttributeError`` in the
       streaming constructor (``vv-principles`` failure mode 8, live in
       production) — and :attr:`MaterialMesh.areas
       <orpheus.transport.mesh.material_mesh.MaterialMesh.areas>` names
       its own three arms instead of blaming a 2-D mesh for all three.
       ⭐ **Two honest-tier corrections landed with the docs audit, and
       the method behind each is the reusable part.** (a) The
       frame-induced divisor was pinned ``np.array_equal`` against the
       retired ``weights.sum()`` spelling over a ladder of eight
       fixtures — a ladder that *skipped* ``gauss_legendre(8)``, the one
       shipped member where the two spellings diverge.  For a finite
       SHIPPED family a ladder is a
       **sample**, and the member you skip is where the counterexample
       lives; the divergence (1 ULP on the divisor, :math:`\leq` 4 nulp
       on the reconstructed kernel) is now pinned as a *bound* by its own
       falsifiable gate row, under the standing
       principled-over-bit-identical criterion.  (b) The pair's
       round-trip rows advertised "bit-exact" from a single random draw:
       `[M]` on that fixture ``np.array_equal`` fails on **844 of 2000**
       seeds (float re-association of :math:`\sum_n w_n / \sum w`), so
       they were re-pinned at ``nulp=1``.  **"Bit-exact" is a property of
       the DRAW** until a sweep makes it a property of the fixture — and
       the shipped S\ :sub:`N` carrier *is* exact (200 of 200), which is
       why the production-facing rows keep ``array_equal`` honestly.
     - #399
     - ``4069155b`` … ``a82d31e4`` (the carrier mints and the
       space-content re-key), ``554ff10b`` / ``1333135e`` (the mesh
       binding retires), ``b00bf2d7`` … ``2690a434`` (space-primary
       construction), ``048144db`` / ``19b85775`` (the collapse pair and
       its frame induction), ``ffb8f286`` (space-derived truncation —
       closes #399), ``78925753`` / ``53e7d207`` (the re-keys and the
       packer), ``1f8e0323`` / ``2e054bfc`` (the S7 un-welds and the
       typed rate co-vector), ``6734bf15`` (the ULP-honesty corrections)
       — merged @ ``55bb47b9``
   * - 2026-08-23
     - **The frame owns the metric its coefficient side carries, and its
       faces are the BOUND operators it mints** (campaign 1, phase CS4b,
       steps F-0 and F-1).  **F-0 — the metric truth, and ERR-039's third
       chapter.**  Every earlier step of the frame campaign moved
       *operators*; none had asked whether the metric the coefficient
       codomain carries is the right one, and it was not.  The frame
       exposed the basis's space unchanged, so the codomain carried the
       continuum Gram :math:`g_C` where the covariant moments that
       analysis emits need its **inverse**.  The theorem is exact and
       unconditional — :math:`\varphi = M\psi = Gc` identically — so the
       Parseval metric is the **inverse of the frame's DISCRETE Gram**, a
       property of the *(basis* :math:`\otimes` *measure)* pair rather
       than a constant of the basis.  ``FrameBase`` gained
       :attr:`~orpheus.numerics.frame.FrameBase.discrete_gram`, a
       MEASURED
       :attr:`~orpheus.numerics.frame.FrameBase.discrete_gram_structure`
       verdict (deliberately distinct from the basis's DECLARED one), and
       a :attr:`~orpheus.numerics.frame.FrameBase.basis_space` that
       dresses with :math:`G^{-1}` on a diagonal frame and **refuses** on
       a dense one — slab Gauss–Legendre is the standing witness that no
       diagonal candidate can satisfy Parseval there.  Nothing about the
       design was wrong; *what was stored* was.  ⭐ **Three shields
       explain why no gate could see it, and the third is the one to
       read.**  The defining adjoint identity held at the round-off floor
       because ``.H`` is BUILT FROM the stored metric — true for every SPD
       metric, so the reading carries zero information about which one is
       installed; composed chains are immune because interior metrics
       cancel; and no end-of-chain adjoint consumer existed yet.  That
       third shield is **latency, not safety** — which is precisely why
       the metric had to be right *before* the S6 adjoint gates landed.
       Full derivation, the declared-vs-measured Gram, the dense refusal
       and the family-wide residual table: :ref:`frame-parseval-metric`.
       **F-1 — the mint.**  With the metric right, the remaining asymmetry
       was ownership.  A frame is not an operator, it is an operator
       **factory**, and it is *shared*: S\ :sub:`N` scattering, the
       windowed accumulation, DSA's :math:`\ell = 1` row and the
       loss-kernel gauge all read one frame.  So
       :class:`~orpheus.transport.frames.harmonic_frame.HarmonicFrame`
       reverts to the two-argument ``(basis, measure)`` factory and MINTS
       the transport-level faces already bound to their two full field
       spaces (:meth:`~orpheus.transport.frames.harmonic_frame.HarmonicFrame.flux_analysis_on`,
       :meth:`~orpheus.transport.frames.harmonic_frame.HarmonicFrame.source_reconstruction_on`).
       The pre-F-1 "not yet typed" ``codomain is None`` debt on
       :class:`~orpheus.sn.operators.windowing.BulkAnalysisOperator` dies
       with the mint, and the composition guards now check that end.
       ⚠ **The one blast-radius miss is worth recording**, because no
       symbol grep could have found it: a duck-typed test surrogate
       spelled the raised mint contract as a keyword argument, so only
       running the suite exposed it.  A call-site set that is *complete by
       symbol grep* is not complete for a **contract** change — its
       consumers' test run is part of the enumeration, not a post-hoc
       check on it.
     - —
     - ``0317373d`` (F-0, the metric truth), ``3dfea889`` (F-1, the mint),
       ``4aa7f951`` (the surrogate contract) — merged @ ``55bb47b9``
   * - 2026-08-22
     - **An operator is not an operator without its two spaces** (campaign
       1, phase CS4b, the S4 amendment).  The root
       :class:`~orpheus.numerics.operator.LinearOperator` returned ``None``
       from ``domain`` and ``codomain`` by default, so a leaf could inherit
       *"I have no spaces"* by saying nothing — and an
       adjointable-but-unbound leaf's ``.H`` then silently returned the
       **Euclidean transpose wearing the Hilbert adjoint's name** — the
       hazard the monomorphic-leaves suite had catalogued, and the same
       conflation family as ERR-076.  The
       amendment makes the base *demand*: both properties are
       ``@property`` ``@abstractmethod``, every class must answer — with
       spaces, by derivation, by the pointwise law, or by an explicit
       documented override naming its owning campaign — and the generic
       adjoint wrapper now **refuses** an unbound, non-metric-free inner
       rather than degrading.  Four structures fall out.  **(a) The
       pointwise family.**
       :class:`~orpheus.numerics.operator.PointwiseOperator` is the
       space-polymorphic stratum of the multiplier algebra — the identity
       (:math:`\times 1`), the endomorphic zero (:math:`\times 0`) and the
       diagonal (:math:`\times f`): endomorphic at the operand's space, no
       stored pair BY LAW, and ``.H = self``, because a real multiplier
       commutes with every diagonal metric.  Being *bound* and being
       *metric-free* are **orthogonal** properties — the bound
       :class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`
       declares itself metric-free while keeping its spaces — and
       :attr:`~orpheus.numerics.operator.LinearOperator.is_metric_free_adjoint`
       is derived recursively through Sum / Product / Scaled /
       TensorProduct, so the new refusal never fires on a composite that
       genuinely needs no metric.  **(b) The zero splits in two:** the
       natural endomorphic
       :class:`~orpheus.numerics.operator.ZeroOperator` (a stateless
       :math:`\times 0` echo, a pointwise member) and the born-bound
       :class:`~orpheus.numerics.operator.ZeroMorphism`, which mints its
       zeros from its own declared pair — which is what retires the
       ``codomain_zero`` / ``transpose_zero`` closures both production
       consumers had been passing.  **(c) The fission pencil un-welds.**
       Its :math:`(B,B)` hooked zero existed only because a coupled
       operator refuses an all-``None`` column, and the honest reading is
       that the *missing operator was the restriction*: :math:`F` is now
       posed as a stack composed with
       :class:`~orpheus.numerics.coupled_system.SystemRestrictionOperator`,
       a born-bound member projection whose adjoint is extension-by-zero
       through the space's own materialization seam.  In execution this
       forced that seam's **cotangent twin** — `[M]` a flux-classed ray
       zero minted into the daggered chain was refused by the cross-class
       gate, so the member CLASS is load-bearing, not incidental.
       **(d) The frame binds** — the amendment's first step bound
       :class:`~orpheus.transport.frames.harmonic_frame.HarmonicFrame`
       itself to its two field spaces at construction.
       ⛔ **Superseded the next day by F-1** (the row above): the binding
       belongs on the *faces* a frame mints, not on the shared factory.
       The amendment's *demand* stands unchanged — it is what made the
       misplacement visible in the first place.
     - —
     - ``33950d81`` (the sweep-straggler gates), ``6e04a749`` (the frame
       binding — superseded, see F-1), ``6fc247fb`` (the restriction and
       the pencil un-weld), ``aa508d3f`` (the pointwise family, the zero
       split, and the ``.H`` refusal) — merged @ ``55bb47b9``
   * - 2026-08-21
     - **Scattering, fission and the isotropic pair are KERNELS; the
       operators carrying them are born BOUND** (campaign 1, phase CS4a,
       with its clear-context adversarial review).  A cross-section
       *datum* and the *operator* that acts with it had been one object.
       A new :mod:`orpheus.transport.kernels` mints
       :class:`~orpheus.transport.kernels.ScatteringKernel`,
       :class:`~orpheus.transport.kernels.N2NKernel` and
       :class:`~orpheus.transport.kernels.FissionKernel` — frozen,
       read-only views over ONE
       :class:`~orpheus.data.macro_xs.mixture.Mixture`, with the emission
       spectrum applied at the single place it belongs
       (:func:`~orpheus.data.emission_spectrum.enforce_emission_spectrum`)
       and a truncation that **refuses** beyond the carried order rather
       than silently padding.  The energy arm gets ONE rule,
       :meth:`EnergyAxis.from_materials
       <orpheus.numerics.axis.EnergyAxis.from_materials>`, which
       :attr:`MaterialMesh.bulk_space
       <orpheus.transport.mesh.material_mesh.MaterialMesh.bulk_space>`
       reads.  On the operator side the pose is minted FROM the mixture,
       the integrated rates read **the space's** measure rather than a
       hand-carried volume vector, one shared energy-extent conformity
       guard is wired into all four constructors, and
       :class:`~orpheus.transport.operators.fission.FissionOperator`'s
       space becomes MANDATORY — the in-tree precedent every later phase
       of the campaign was measured against.  That none of this moved a value
       was **measured, not assumed**: the byte gate
       ``tests/homogeneous/test_byte_stability.py`` holds the homogeneous
       solve bit-exactly across the rewiring on **8 of 8** rows —
       exhaustive over the producing mixtures the tree ships — and the
       frozen reaction-rate references moved 0 ULP.  ⭐ **The review round's own findings are the
       durable part.**  A per-site witness table for the hoisted
       conformity guard found that only **1 of 4** call sites had a
       witness at all — and the witness-less site was the one whose
       operand expression differed from its siblings', which is now a
       ``vv-principles`` #17 sharpening (*a guard hoisted to one shared
       home has as many arms as it has CALL SITES*).  A physics-semantics
       question that had been recorded merely as an *observation* was
       decided: **condensed cross sections are INTENSIVE**,
       :math:`\bar\sigma_x = \langle \Sigma_x, \varphi\rangle /
       \langle 1, \varphi \rangle` — the shipped spelling had scaled with
       the point weight.  And thirteen first-phase attacks were
       **WITHDRAWN with structural reasons**, three of them refuted by
       their own author's probe; a refuted candidate is first-class
       output.
     - —
     - ``069e2caa`` / ``15bbf935`` / ``9f1d4190`` / ``49b29391`` (the
       kernels, the mixture-minted pose, F's mandatory space),
       ``c7f9fa8d`` / ``d61e097b`` (the review round's production repairs
       and gate strengthening) — merged @ ``55bb47b9``
   * - 2026-08-20
     - **The space layer gains AXES, and the homogeneous solver poses its
       problem on a real Energy space** (campaign 1, phase CS1) — the
       first factor of what will become the S\ :sub:`N` composite.
       :mod:`orpheus.numerics.axis` mints
       :class:`~orpheus.numerics.axis.Axis` (frozen, structural identity
       per subclass, canonical measure storage) and
       :class:`~orpheus.numerics.axis.EnergyAxis`, and
       :meth:`FunctionSpace.of_axes
       <orpheus.numerics.space.FunctionSpace.of_axes>` composes a space as
       the ordered product of its axes, with a per-axis metric path (no
       densification) and a deterministic, injective derived name — the
       **identity bridge** that makes *"metric differences imply space
       differences"* true today rather than aspirationally.  The
       S\ :sub:`N`-facing consequence is the operator slot: :math:`S` and
       :math:`F` stop carrying a ``full_field_space`` and carry a
       ``space`` of the family type, so an operator can be posed on the
       composite OR on a bulk factor without a second spelling — the slot the
       kernel-binding phase then tightens to MANDATORY on :math:`F`.  The
       homogeneous solver poses :math:`A = C - K_{\rm iso}` and
       :math:`F` on that real space, retiring both production
       ``basis_shape=(ng, 1)`` spellings and turning the ``OperatorSum``
       space guard from *skipped* into *validating*.  Full account, with
       the counting-measure theorem that explains why none of it moved a
       value: :ref:`spaces-the-axis` and
       :ref:`spaces-counting-measure-theorem` on
       :doc:`/theory/foundations/spaces`.
     - —
     - ``1afff47b`` / ``f4876354`` (the axis and ``of_axes``),
       ``e8769897`` / ``24a991ba`` (the operator space slot),
       ``6bd782ab`` (the homogeneous pose), ``6da1b23c`` (the cone
       consult), ``37122fd6`` (the corpus seed) — merged @ ``55bb47b9``
   * - 2026-08-15
     - **The reflective boundary trace has a canonical value, and every solve
       exit returns it.**  :math:`A = L + C - S - B` is **exactly singular** on
       any :math:`d \geq 2` Cartesian **diamond-difference** mesh with
       :math:`\geq 2` reflective axis pairs — which is the default
       :math:`k_\infty` lattice, so this was the ordinary case, not a corner.
       `[M]` :math:`\dim\ker A = 12` at d=2 LS4 ng=2 (**mesh-independent**:
       :math:`n_g N/4`) and :math:`138` at d=3 (3,4,5).  **Mechanism:** with
       :math:`\psi_c \equiv 0` the DD closure
       :math:`\psi_{\rm out} = 2\psi_c - \psi_{\rm in}` degenerates to the
       involution :math:`\psi_{\rm out} = -\psi_{\rm in}`, so the face-to-face
       transmission :math:`\Sigma = (2/D)\mathbf{1}w^{\mathsf T} - I` carries
       eigenvalue :math:`-1` with multiplicity :math:`d-1` on
       :math:`\{v : w^{\mathsf T}v = 0\}` — an undamped **face sawtooth**
       :math:`\psi_a(k, i_\perp) = (-1)^k \varphi_a(i_\perp)` invisible to
       :math:`\Sigma_t V \psi_c`, and around a closed reflective loop the
       :math:`-1`\ s compose to :math:`+1`.  `[M]`
       :class:`~orpheus.transport.spatial.linear_discontinuous.LinearDiscontinuous`
       on the IDENTICAL box has :math:`\dim\ker A = 0` — that substitution is
       what proves the mechanism rather than arguing it.  **What a user saw:**
       the returned boundary trace was a function of the cold start, up to
       **27.3 %** apart between starts differing only inside :math:`\ker A`,
       while the bulk was bit-stable at :math:`7\mathrm{e}{-16}` and BOTH
       convergence certificates read converged — a spurious current
       *tangential* to a mirror face, which is physically impossible.
       ⭐ **The repair is CANONICAL, not conventional, and that is a theorem:**
       every kernel mode carries a non-trivial sign character on every axis, so
       every **mirror-even** functional annihilates :math:`\ker A` exactly;
       :math:`\psi_{\rm exact}` is mirror-even, hence :math:`G`-orthogonal to
       the kernel (`[M]` :math:`1.27\mathrm{e}{-15}`), hence the
       **minimum-**\ :math:`\lVert\cdot\rVert_G` **member IS the physical
       answer**.  That same theorem explains why nothing could have caught it:
       every summed trace functional is blind by symmetry.
       ⟹ :class:`~orpheus.sn.operators.loss_kernel_gauge.LossKernelGauge`, a
       :math:`G`-orthogonal projector built from a **closed form** (the ANOVA
       pair generators — no SVD of :math:`A`), as a **direct sum over
       (ordinate orbit** :math:`\times` **group) blocks** with disjoint
       supports, each block :math:`G`-orthonormalised by one fused
       :math:`\sqrt{G}`-weighted SVD that does rank reduction and
       orthonormalisation together — which is what earns
       ``GramStructure.DIAGONAL``.  Applicability is **DERIVED, never
       tabulated** (the ask-don't-tabulate ruling): the closure is asked
       whether it leaves a zero-mean face mode undamped
       (:meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.face_transmission_spectrum`)
       and the mesh is asked how many axis pairs close
       (:attr:`~orpheus.sn.mesh.augmented_mesh.SNMesh.reflective_axes`), so a
       future scheme answers for itself with no edit.  The gauge fires at
       **all four public entries**, records its magnitude on
       :class:`~orpheus.sn.solution.IterationHistory`, and emits a
       ``GaugeFreedomWarning`` that names the **root** fix (switch to a damping
       closure, or break a reflective axis pair) rather than only reporting the
       projection.  Scope is the underdetermined remainder :math:`R`; the
       tangential component :math:`T` (where :math:`G \equiv 0`) has no
       minimum-norm representative and is deliberately untouched.
       ⛔ **Two premises of the campaign's own plan were refuted by
       measurement** and the corrections are the durable part: "excited iff the
       first axis has an ODD cell count" was drawn from 11 meshes that all
       carried a *uniform isotropic source* — `[M]` :math:`\dim\ker A` is the
       same at **every** parity, and an anisotropic source excites an even
       (4,4) mesh at :math:`1.756\mathrm{e}{-02}`, so **kernel-freedom is never
       inferable from a mesh property**; and the Krylov arm's returned boundary
       was documented as a "face residual" when `[M]` it is a flux trace,
       measured three independent ways.
     - #344
     - ``5def63b0`` (the derived predicate +
       :class:`~orpheus.numerics.operator.InverseMetricOperator`),
       ``f934ff57`` (the closed-form kernel + the gauge), ``b51bc802`` (fired
       at every exit, warned, recorded), ``1a2be025`` (the characterization
       promoted into CI)
   * - 2026-08-11
     - **The azimuthal cell partition is taken in ω, not in η — and the
       :math:`[\tfrac12, 1]` absorber retired with the defect it was
       hiding** (Q5.6.4, the quadrature machinery campaign).  The
       campaign's own acceptance test — *"retire the absorber and the
       azimuthal floor falls"* — was **REFUTED by measurement**: the
       naive retirement makes the anisotropic cylinder MMS floor
       :math:`1.8`--:math:`3.4\times` WORSE at every rung.  Chasing that
       found the real defect, one level below.  :math:`\alpha` and
       :math:`\tau` both reference **ONE object** — the boundary between
       azimuthal cell :math:`m` and :math:`m+1` — and each derived it
       independently, in disagreement: :math:`\alpha` at the real
       half-angle :math:`\omega_{m-1/2}`, :math:`\tau` at the CHORD
       midpoint :math:`(\eta_m + \eta_{m+1})/2`.  Because
       :math:`\cos` is nonlinear, every interior chord edge is
       :math:`\cos(\Delta\omega/2) \times` the arc edge while the two
       ENDPOINTS stay pinned at :math:`\mp\sin\theta` unscaled — so the
       outer cells stretch to absorb the shrink.  `[M]` the :math:`\eta`
       error vanishes as :math:`\Delta\omega \to 0` but the implied
       :math:`\omega`-width spread does NOT: it converges to
       :math:`\approx 17.45\,\%` (18.71 / 17.59 / 17.48 / 17.46 % at
       :math:`n_\varphi = 8/16/32/64`) against a quadrature whose own
       cells are bit-exactly equal.  **That** :math:`O(1)`
       **inconsistency is what the absorber was compensating for**, which
       is why removing it alone regressed.  The fix is one word in the
       right vocabulary: the azimuthal march is a march in
       :math:`\omega`, arc by arc, so the cell boundary is the
       **midpoint in** :math:`\omega`.  With a partition chosen,
       predicate **P2** (BMC Eq. 43 = Lathrop Eq. 23 — :math:`\tau` is
       the barycentric coordinate of the ordinate between its own cell's
       edges) *determines* :math:`\tau`; closed form `[M]` verified to
       :math:`1.67\mathrm{e}{-16}`:
       :math:`\tau_m = \tfrac12 + \tfrac12\cot\omega_m\tan(\Delta\omega/4)`.
       ⟹ **one partition producer**
       (:func:`~orpheus.sn.angular.closure.angular_cell_edges_per_level`,
       the only place a cell boundary is defined) and **one τ producer
       with a geometry-FREE body** — ``morel_montry_tau_raw_per_level``
       retired, because the raw/clamped distinction it named no longer
       exists.  Verified solve-free by the **ν-closure** diagnostic (does
       the march implied BY :math:`\tau` land on the level's own
       endpoint?): ``1.000000000000`` for the derived :math:`\tau`, where
       the clamp overshoots by 1.6 % and :math:`\tau \equiv \tfrac12` by
       16.5 % — i.e. neither corresponds to any partition of the level.
       The SPHERE is untouched (cumulative-WEIGHT edges, BMC Eq. 12
       verbatim, literature-confirmed); the sphere's convention provably
       cannot be transplanted (an arc cell's :math:`\eta`-measure
       :math:`\propto \sin\omega` while a trapezoid weight is constant, so
       accumulating weights in :math:`\eta` violates P3 and worsens with
       refinement — `[M]` 0/4 → 4/8 → 12/16 → 28/32 ordinates outside
       their own cell, NaN from :math:`n_\varphi \ge 16`).  ⭐ Cylinder-P3
       became a **theorem** (on a monotone arc the ω-midpoint edges
       bracket their own node, so :math:`\tau \in (0,1)` is forced), which
       reduces it to the fold criterion :math:`\Sigma = \emptyset`; a
       non-monotone (full-circle) level is now refused *by the partition
       producer*, naming the double cover.  ⚠ **Honest cost, ratified not
       hidden**: at :math:`n_x = 320` the principled :math:`\tau` is
       BETTER at :math:`n_\varphi = 8` (3.128e-3 vs 3.511e-3) and
       :math:`\sim 1.8`--:math:`2\times` WORSE at 16/32/64.  Principled
       :math:`\ne` more accurate — and the L2 norm measures truncation
       order, which is exactly what :math:`\tau \equiv \tfrac12`
       optimises and exactly what is blind to the diffusion limit
       :math:`\tau` exists to fix.  ⚠ **α and τ must NOT be "unified"**:
       they share the partition but impose different conditions (α the
       first moment, τ the zeroth); forcing α onto the geometric
       tangential cosine drives :math:`\delta \to 0`, i.e. the angular
       *diamond* scheme.  Companion commit gave quadrature/closure
       analysis a home
       (:mod:`orpheus.derivations.discrete.sn.angular_differencing` — the
       P0--P4 predicate ladder, the τ/β nomenclature, and a written record
       of which diagnostics are BLIND on which rules) and retired
       ``contamination.py``, whose cylindrical arm had become
       present-tense wrong (it built the retired η-midpoint edges, so
       `[M]` its τ disagreed with production by up to 6.8e-2).  Full
       treatment: :ref:`sn-tau-absorber-retirement` and
       :eq:`angular-cell-partition` in
       :doc:`/theory/foundations/structured_geometry`.
     - #229 (record) · #327 (naming)
     - ``3dda18ca`` · ``d5067c4d``
   * - 2026-08-08
     - **The cylindrical admission flip — SNMesh(CYLINDRICAL) admits
       exactly the carrying quadrature rules** (Q5.6.3, the quadrature
       machinery campaign).  The R12a march-start predicate was
       promoted from a classifier to the **admission law**:
       ``SNMesh`` construction calls
       :func:`~orpheus.sn.angular.closure.assert_carrying_quadrature`
       (offender positions via
       :func:`~orpheus.sn.angular.closure.non_carrying_levels`),
       refusing any cylinder rule with a non-carrying μ-level and
       naming the facts true on the first offender plus the remedy
       (:meth:`Quadrature.folded_product
       <orpheus.numerics.quadrature.Quadrature.folded_product>`).
       Admission reads **structure, never provenance** — a hand-built
       σ_y-quotient with the right arrays admits; a tagged quotient of
       a node-aligned parent refuses (pinned by the admission module's
       pincer pair, plus a foreign hand-built GL-in-φ arc rule that
       admits and cold-solves as a true inverse at 4.5e-16).  With
       every admitted level carrying, the **#280 2.5b direct-seed fold
       was retired** with its whole family (builders, transpose twin,
       every ``is_seed_ord`` branch): its precondition — a
       non-carrying admitted cylinder — became unrepresentable, so
       every admitted geometry rides route (a)'s forward substitution
       with genuine independent seeds.  Guarded by the 16-row
       admission module
       (``tests/sn/mesh/test_cylindrical_quadrature_admission.py``)
       and a 10/10 mutation battery (positive control first).
       Follow-ons filed: #338 (Gauss–Lobatto sphere-side admission
       interaction), #339 (the LS double-fold capability).
     - #280 · #326
     - ``1689faf4`` · ``1f220c41``
   * - 2026-07-28
     - **Naming honesty in the operator algebra — three symbols stopped
       lying about what they hold** (Tier 1 of the operator-realization
       plan review). (a) ``SNSolver``'s ``(L, S, F)`` "operator triple"
       was **retired**, not renamed: it had no production reader, and
       ``self.L`` held the *composite* :math:`L+C` while :math:`L`
       everywhere else is the :math:`\sigma`-free streaming leaf — three
       meanings of ``L`` inside one ``__init__``. Every solve builds its
       own composite through
       :func:`~orpheus.sn.coupled_system.build_within_group_system`, so a
       solver-held copy was a twin free to drift from the operand the
       sweep inverts. (b) ``InvertibleOperator`` →
       :class:`~orpheus.sn.operators.streaming.StreamingCollisionOperator`:
       the class wraps exactly a ``StreamingOperator`` plus a
       ``MultiplicationOperator`` — it IS :math:`L+C` — and its old name
       stated its *capability* (invertibility is advertised via
       ``is_invertible``/``inverse()``), leaving the sole factory
       ``build_streaming_collision`` out of step with the type it returns.
       (c) :class:`~orpheus.sn.coupled_system.WithinGroupSystem`'s
       ``resolvent``/``gains`` → ``implicit_operator``/``explicit_gains``:
       the field held :math:`M`, the un-inverted forward operator, where a
       resolvent is inverse-like. The new names state the role in the
       splitting (:math:`M` solved implicitly, :math:`N` evaluated
       explicitly) and free ``resolvent`` for its two honest uses — the
       corpus :math:`K_{\rm pm}=A_{\rm loss}^{-1}M` and the future
       ``A.resolvent(z)`` factory. The numerics drivers keep ``A`` /
       ``*gains`` deliberately: that is the documented solver ↔ numerics
       bridge of :doc:`/theory/conventions/notation` row 8, which now
       states both ends instead of naming a dead field.
     - —
     - ``fc5c41bc`` · ``8367346f`` · ``e2c8b32e``
   * - 2026-07-24
     - **The adjoint matvec completed the scheme × representation grid —
       the #280 residue retired** (adjoint completion campaign, C1–C5).
       The transpose kernel became a scheme-REGISTERED pair
       (``residual_kernel_batch_transpose``, the exact VJP mirror of the
       forward batch kernel; ``has_transpose_kernel`` derived from
       registration — declared-True-with-no-kernel unrepresentable), DD's
       hand-coded reverse relocated into it bit-identically (C1) and LD
       generated from the UBLD algebra-of-record as :math:`A^{\mathsf T}
       M^{-1}` with the SymPy transpose-oracle + the ``θ``-mass bulk
       metric ruling (C2, the LD-slab adjoint).  The multi-D reverse
       landed as the **mirror-octant realization** — the reverse of
       octant :math:`o` is the UNCHANGED ``walk_full``/``walk_windowed``
       over the mirror graph :math:`\mathrm{graph}(-\text{signs})`, whose
       levels/face-roles/frontier ARE the forward's reversed (zero
       forward edits; C3 oracle, C4 windowed production bit-identical to
       it) plus the ScanMarch row-march reverse (the forward's
       sign-reading spellings fed the mirror label; the kernel's x-out
       cotangent structurally zero; ``_x_scan_faces_transpose`` +
       ``reflect_scan_coefficients_transpose``), verified d-generic
       (d=3 spine row) and scheme-complete (C5: the LD-2D moment-tailed
       cochain through the same frame — dense-``Mᵀ``, the
       assembled-``Mᵀ`` keystone, reverse ``window ≡ full`` at
       ``n_face_moments = 2``, the exact moment-drop and cross-moment
       frame-sign parity teeth).  ``.H`` now constructs on the whole
       registered grid; the deferral ledger holds only the G-S
       schedule-reverse (R7, no consumer).  See
       :ref:`loss-rep-orientation-two-frames`.
     - #310 (+#311)
     - C1–C4 merged @ ``59830618``; C5 = this entry's own merge
   * - 2026-07-12
     - **The walk's fused ψ½ joint channel was retired — every sweep
       surface IS the ray-decoupled** :math:`(A,A)` **block, presence
       structural** (coupled-block campaign, step 6). With the four blocks
       posed and the solve routed through the grid, the walk's fused
       joint-leg channel became dead production code and was deleted
       (net −803 lines): the two forward/transpose ψ½ kwarg pairs across
       the six representation signatures, the three-function presence-guard
       family, the walk's in-solve System-B engines, and the operator-free
       ``transport_sweep`` wrapper. **(6a, ``03e275e8``)** the eigenvalue
       finalize re-routes through the SAME
       :func:`~orpheus.sn.coupled_system.build_within_group_system`
       ``.resolvent`` every driver consumes. **(6b, ``015dcc73``)** the
       estate cut: every walk surface is now the ray-decoupled
       :math:`(A,A)` diagonal block on *every* mesh — the matvec
       substitutes a zero seed into the Morel–Montry thread, and the
       transpose **discards** the thread cotangent (a fixed input's
       cotangent propagates nowhere; the
       :class:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicSeeding`
       grid block carries the coupling). :attr:`SweepOperator.is_adjointable
       <orpheus.sn.operators.sweep_operator.SweepOperator.is_adjointable>`
       dropped its carrying-mesh third factor (R-6.2). The mesh stays the
       single authority on presence — nothing *checks* against it anymore,
       the type system carries the biconditional. Sphere re-baselines are
       sphere-only FP-grain (the fused → block-sum association class
       ~1e-15); slab/cylinder byte-identical through a full re-save. See
       :ref:`coupled-block-operator` in :doc:`/theory/foundations/coupled_block_operator`.
     - #280
     - ``015dcc73`` (merged @ ``3f0b8c74``)
   * - 2026-07-12
     - **The four BC-layer invariants (ERR-041/042/045/047) went
       production-live at the realize seam** (coupled-block campaign).
       ``BoundaryTraceLaw.assert_realizable`` became the certification
       template method — ``SNBoundaryRealizer.realize`` fires it ONCE at
       entry, so every ``SNMesh`` construction certifies its BC laws:
       the reflection table must be measure-preserving under the
       :math:`w\,\lvert\mu_{\rm axis}\rvert` trace metric (ERR-042,
       independent of involution), every non-tangential ordinate's
       partner must map inflow→outflow with opposite sign on the law's
       axis (ERR-045), the vacuum arm cross-checks claimed
       ``inflow_indices`` against the orientation the face name alone
       implies (ERR-041), and a nonzero boundary source without an
       inflow mask is uncertifiable while a masked one is delivered
       exactly on :math:`\Gamma_-` (ERR-047). Three measured-independent
       invariants (a mutant passing involution AND measure reds only the
       inflow→outflow check); all four catchers mutation-verified;
       ``tests._harness.audit`` ERR coverage 64/68 → **68/68**.
     - —
     - ``51c22396`` (merged @ ``3f0b8c74``)
   * - 2026-07-12
     - **The within-group coupled solve landed — block-triangular
       substitution + a materialised-LU EXTRACT, with a ρ-honest stop and
       a lag-death certificate** (coupled-block campaign, step 5). The
       numerics :class:`~orpheus.numerics.coupled_system.CoupledOperator`
       gained the structure-keyed DIRECT solve: **(5a, ``6732778a``)** the
       triangular block substitution (ONE body, four orientation × transpose
       combos) and the materialise/LU EXTRACT via
       :class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`;
       **(5b, ``899ee06a``)** the resolvent :math:`M` re-poses as the honest
       upper-triangular grid ``[[LC, Seeding], [None, march]]``, whose
       ``solve`` is the block back-substitution (System B's march first,
       then the bulk sweep on the ray-decoupled :math:`q_A -
       \text{Seeding}\,\psi_B`) and whose ``inverse()`` is the
       :class:`~orpheus.numerics.coupled_system.CoupledSubstitutionOperator`;
       **(5c, ``c98a23d8``)** :class:`~orpheus.numerics.iteration.SourceIteration`
       stops on the free-identity equation residual
       :math:`r_n = \mathrm{rhs}_{n-1} - \mathrm{rhs}_n = A\psi_n - q`, and
       every full-angular arm carries the driver-level
       :class:`~orpheus.sn.solver.ConvergenceCertificateError` — one honest
       :func:`~orpheus.sn.solver.evaluate_residual` per claimed exit closes
       the exact-:math:`M` hole the free identity assumes (the #282
       lag-death class, measured cold-residual defect ~5e5); **(5d,
       ``88e226ff``)** the fused ``CoupledInvertibleOperator`` bridge was
       DELETED (every joint-``M`` fixture rides the one grid helper); **(5e,
       ``0e03c304``)** the :math:`A.H.\text{inverse}() \equiv
       A.\text{inverse}().H` swap-law arm went live on the grid. The EXTRACT
       is principled-equivalent (:math:`5.5\text{e-}16` vs the dense-LU
       oracle), not bit-identical. See :ref:`coupled-block-operator` in
       :doc:`/theory/foundations/coupled_block_operator`.
     - #280 #282
     - ``0e03c304`` (merged @ ``3f0b8c74``)
   * - 2026-07-11
     - **The ψ½ ray solve was un-woven from the walk — System B is marched
       split-native and its ray solve routes through the named ``A_BB``
       resolvent** (coupled-block campaign, Phase C 4e). Three moves
       complete the Cardinal-Rule-2 single source for the curvilinear ray
       orchestration. **(e1, ``63702e7``)** the fused :math:`(L+C)` 1-D walk
       marches System B's ``interior ⊕ boundary`` composite **natively**,
       retiring the historical **unified** ψ½ leaf family (cells ⊕ corner
       interleaved on one ``FaceField[(level, sign, part)]`` buffer), its
       unified space, and the ``from_unified`` / ``to_unified`` bridge (the
       demotion licence ``to_unified ∘ from_unified == id`` discharged by
       the walk-slot rewrite landing bit-identically). **(e1b,
       ``ea7f919c``)** the freed ``RadialCharacteristicField`` name was
       reminted onto System B's composite —
       :class:`~orpheus.transport.radial_characteristic_field.RadialCharacteristicField`
       (``Composite[interior, boundary]``), the exact
       :class:`~orpheus.transport.full_field.FullField` mirror (ONE
       System-B carrier name). **(e2, ``98fe2e36``)** the walk's two welded
       ray-solve orchestrations were **deleted**: the :math:`(L+C)` walk
       now routes System B through the named resolvent
       :meth:`RadialCharacteristicOperator.solve <orpheus.sn.operators.radial_characteristic.RadialCharacteristicOperator.solve>`
       / :meth:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicOperator.solve_transpose`
       (``A_BB``, built inside the walk over its own :math:`\sigma_t`), so
       the two-leg march **orchestration** lives in ONE place and the DD
       **engine**
       (:func:`~orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source`)
       in one other — the walk's ``carlson_inward_sweep_*`` references went
       **8 → 0**. Under the **H1-narrow** ruling only the ``A_BB`` *solve*
       legs were extracted; the ``A_AB`` seed → bulk coupling stays fused in
       the within-group :math:`M` (the transpose augments the seed cotangent
       with the Morel–Montry thread cotangent — the fused
       :math:`A_{AB}^{\mathsf T}` feed — before one ``solve_transpose``
       call). Landed **bit-identical**: the frozen walk baselines hold at
       ``nulp = 1`` (zero drift) and :math:`k` matches to 15 digits across
       sphere/slab × SI/Krylov (full tree 6340/0). See
       :ref:`sn-direct-seed-solve`.
     - #280
     - ``98fe2e36`` (merged @ ``3f0b8c74``)
   * - 2026-07-11
     - **The ψ½ ray was evicted from the SN composite carrier —
       ``FullField`` is pure 2-block and System B is its own composite**
       (coupled-block campaign, Phase B.2d d2, "the atomic eviction"). The
       transitional optional-third-block ψ½ passenger on
       :class:`~orpheus.transport.full_field.FullField` — with a mesh-keyed
       *mixed-presence law* and runtime presence pins — is **retired**:
       ``FullField`` is now the pure 2-block
       ``Composite[BulkField, BoundaryField]``, and the ψ½ ray is
       **System B**, its own 2-block
       :class:`~orpheus.transport.radial_characteristic_field.RadialCharacteristicField`
       coupled to System A through the within-group 2×2 grid as a
       :class:`~orpheus.numerics.coupled_system.CoupledField`. A live-ray
       ``ψ_A`` is now **unrepresentable** (the type system is the guard),
       so the B.2c dead-slot double-count hazard dissolved structurally and
       the coupled flat dimension is the honest two-system sum. The six
       loss-representation walk signatures re-typed to **explicit ψ½ leaf
       kwargs** (forward ``radial_characteristic_flux`` /
       ``radial_characteristic_source``; transpose ``seed_cot`` /
       ``seed_cot_out``); the unified ``RadialCharacteristicResidual`` leaf
       split into a typed interior ⊕ boundary residual pair; and the
       converged ray is returned as
       :attr:`Solution.radial_characteristic <orpheus.sn.solution.Solution.radial_characteristic>`,
       System B's own member with a presence biconditional. See
       :ref:`sn-direct-seed-solve`.
     - #280
     - ``e5d1acf`` (merged @ ``3f0b8c74``)
   * - 2026-07-11
     - **The within-group solve became block-native — one
       ``WithinGroupSystem`` record carrying the named** :math:`A = M - N`
       **splitting** (coupled-block campaign, Phase B.2d d1). The former
       ``_within_group_triple`` / ``_lagged_gains`` construction pair
       retired into a single builder,
       :func:`~orpheus.sn.coupled_system.build_within_group_system`, which
       returns the frozen
       :class:`~orpheus.sn.coupled_system.WithinGroupSystem` record — the
       loss grid together with its Hackbusch splitting
       :math:`A = M - N` (``resolvent`` = :math:`M`, the sweepable part
       inverted each step; ``gains`` = :math:`N`, the lagged couplings —
       both fields renamed at the 2026-07-28 row). This entry read
       *"Hackbusch* **regular** *splitting"* until 2026-08-09 (#341); the
       word is **struck**, not tombstoned, because it never described
       anything that landed — it asserted a property the splitting does
       not have (:ref:`sn-boundary-gs-not-regular`). On
       a carrying sphere :math:`M` is the honest upper-triangular grid
       ``[[LC, Seeding], [None, march]]``
       (:class:`~orpheus.numerics.coupled_system.CoupledOperator` — since
       step 5 its ``solve`` is the numerics block back-substitution and
       its ``inverse()`` the
       :class:`~orpheus.numerics.coupled_system.CoupledSubstitutionOperator`;
       the fused ``CoupledInvertibleOperator`` bridge deleted at 5d) and
       :math:`N` the coupled gain grid ``[[S+B_a, ∅], [Emission, B_b]]``
       (the ``(A,B)`` slot structurally zero — seeding lives in :math:`M`);
       on a seedless mesh it degrades to the bare ``((L+C), (S, B_a))`` the
       Gauss-Seidel / windowing paths consume zero-touch. The four
       within-group solve sites consume the one record. See
       :ref:`bc-extraction-variadic-driver` in :doc:`/theory/foundations/boundary_conditions`.
     - #280
     - ``c0f23f6`` (merged @ ``3f0b8c74``)
   * - 2026-07-05
     - **The adjoint-inverse swap law wired the reverse-scan into the
       operator algebra** (#280 Phase 2.5c). With the reverse-scan in
       hand, the missing adjoint was surfaced through the landed #226
       inverse-as-operator taxonomy *without new solve machinery*:
       :class:`~orpheus.sn.operators.sweep_operator.SweepOperator`
       (:math:`(L+C)^{-1}`) gained
       :meth:`apply_transpose <orpheus.sn.operators.sweep_operator.SweepOperator.apply_transpose>`
       delegating to the inner's
       :meth:`~orpheus.sn.operators.streaming.StreamingCollisionOperator.solve_transpose`
       (the 2.5b reverse-scan), and
       :meth:`_AdjointOperator.inverse() <orpheus.numerics.operator._AdjointOperator.inverse>`
       returns ``inner.inverse().H`` — making the **swap law**
       :math:`(A^{*})^{-1}=(A^{-1})^{*}`
       (``A.H.inverse() ≡ A.inverse().H``,
       :eq:`loss-rep-adjoint-inverse-swap`) an **object identity** of
       the algebra, not a numerical coincidence. The metric
       adjoint-solve then falls out of the existing
       :meth:`_AdjointOperator.apply <orpheus.numerics.operator._AdjointOperator.apply>`
       **for free** — no ``_AdjointOperator.solve``, no metric code in
       the sweep (the A3 "Deliverable 3" dissolved). Companion: the
       vestigial ``initial_guess`` seed-threading retired — direct exact
       inverses accept-and-drop it, the genuine warm start lives at the
       iteration layer. See
       :ref:`loss-rep-inverse-adjoint-swap` and
       :ref:`loss-rep-initial-guess-warm-start`.
     - #280 #226
     - ``8cf5215`` (merged @ ``3f0b8c74``)
   * - 2026-07-05
     - **The adjoint inner solve** :math:`(L+C)^{-\mathsf T}` **landed —
       the empty cell of the 2×2 filled** (#280 Phase 2.5b). The four
       faces ``{forward, transpose} × {solve, apply}`` are a 2×2 whose
       transpose-solve cell never existed;
       :meth:`~orpheus.sn.loss_representation.LossRepresentation.sweep_transpose`
       fills it as the coherent **reverse-scan** of the forward
       sweep-scan — built on
       :func:`~orpheus.sn.sweep.scan.ordinate_scan_transpose` (the
       affine scan's own adjoint, one source of truth — *not* a
       reverse-loop bolted onto the apply path) plus the Hébert §3.9.4
       seed-march adjoint
       :func:`~orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_transpose`.
       1-D DD, all geometries (slab / sphere / cylinder); the keystone
       catcher is the dense :math:`(L+C)^{-\mathsf T}` oracle built from
       the **forward** ``apply`` alone (structurally independent of the
       reverse-walk code). It carries the #284 **source-subspace**
       faithfulness — matching :math:`M_{\rm solve}^{\mathsf T}` on
       every source-carried slot, deviating only on the provably-zero
       outflow column. LD-slab and multi-D adjoint stay typed
       deferrals. This is the substrate for the adjoint campaign A3
       (#276): A4's daggered eigenvalue is the first consumer of
       ``sweep_transpose``. See :ref:`loss-rep-two-frames`.
     - #280 #276
     - ``f1ddeb6`` (merged @ ``3f0b8c74``)
   * - 2026-07-05
     - **The product-cylinder cold solve became a single-pass direct
       inverse — the direct-seed fold** (#280 Phase 2.5b). For a **product**
       quadrature the starting direction coincides with the first-swept
       ordinate (:math:`\tau_{{\rm raw},0}=0`, the #229 clamp fact), so the
       Morel–Montry seed is a **live self-coupling** on the :math:`m_0`
       diagonal (:math:`c_{\rm in}[m_0]\ne 0`) and the cold ``(L+C).solve``
       was seed-**lagged** (cold error :math:`\approx 0.57`). The fold folds
       :math:`\kappa=(\Delta A/w)\,c_{\rm in}` into the diagonal
       (:math:`c_{\rm out}\to c_{\rm out}-c_{\rm in}`), so the cold solve is
       now a single-pass direct inverse for **every** geometry (cold error
       :math:`\to 4.4\times10^{-16}`; scattering-fixed-point
       baseline-neutral). This also **corrected** the long-standing
       "cylinder :math:`\alpha`-dome telescoping absorbs the wrong seed"
       mis-attribution — a level-symmetric-only artefact (the *dead*
       first-ordinate weight :math:`c_{\rm in}[m_0]=0` at raw
       :math:`\tau=1`), **false for a product quadrature**. Companion to the
       sphere route-(a) fix (#282/2.5d) directly below; see
       :ref:`sn-phase-d-gate-1-1-empirical`.
     - #280
     - ``ba202a1`` (merged @ ``3f0b8c74``)
   * - 2026-07-04
     - **The curvilinear starting-direction ψ½ seed became first-class
       typed state — the spherical solve is now a single-pass direct
       inverse** (Issue #282 route (a), #280 Phase 2.5d). The lagged
       Morel–Montry half-angle pole seed (a two-point extrapolation of the
       *previous* source-iteration iterate) was a **walk-order back
       edge** that made the spherical within-group SOLVE a *non-direct*
       inverse (cold residual :math:`5.18\times10^5`). Route (a) promotes
       :math:`\psi_{1/2}` to a **typed** ψ½ role quadruple (flux / source /
       displacement / residual) carrying the SPD :math:`V_{\rm cell}` state
       metric, and the sweep marches it **directly** from the true q½
       source (the
       Hébert §3.9.4 recurrence
       :func:`~orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source`
       + the **full** Legendre fold
       :func:`~orpheus.numerics.spaces.radial_characteristic_space.fold_moments_to_radial_characteristic`,
       :eq:`sn-direct-seed-anisotropic-source`). The augmented :math:`(L+C)` is
       block-lower-triangular in the seed-first walk order
       (:eq:`sn-direct-seed-block-triangular`), so the back edge is dead by
       construction (cold residual :math:`\to 2.5\times10^{-16}`; seed
       insensitivity :math:`4.57\times10^{-2}\to 0` bitwise). The whole
       ``PsiHalfAngleSeed`` strategy zoo retired (851 → 161 lines). The
       eigenvalue re-pose is principled by an angular-order N-sweep (a
       seed is a closure — h→0 is the wrong test). See
       :ref:`sn-direct-seed-solve`.
     - #282 #280
     - ``a29ab2d`` (merged @ ``3f0b8c74``)
   * - 2026-07-04
     - **The assembly mode landed — the sweep's per-cell closure algebra
       emitted a third way, as a sparse matrix** (stencil-assembly
       campaign, Phase 2b). Beside *solve* (sweep) and *apply* (matvec),
       the same per-ordinate ``L(+C)`` coefficients are emitted as
       ``(row, col, value)`` by a **closure-generic symbolic walk** of the
       :class:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph`
       (:func:`~orpheus.sn.loss_representation.assembly.assemble_ordinate_blocks`;
       DD + LD). The emitter owns **no** stencil spelling — it extracts
       every coefficient by unit probes of the production kernel
       (:meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.residual_kernel_batch`,
       exact by ``is_linear``), so solve/apply/assemble share ONE source;
       the LD multi-moment block is conjugated sweep→global by the
       :func:`~orpheus.transport.spatial._ubld.octant_moment_frame_signs`
       involution (:eq:`ld-ubld-octant-moment-frame-signs`). The assembled
       block is lower-triangular in walk order, so LAPACK
       ``solve_triangular`` reproduces the production sweep at
       :math:`\sim 6\times10^{-16}` — **#284 discharged object-level** (the
       sweep IS forward substitution on the source subspace). Curvilinear
       assembly stays OUT (**#282 characterized** as a walk-order back
       edge — the lagged pole seed reads later-ordinate columns). Numerics
       carrier + three-layer ``assemble()`` surface + composer laws land
       in :doc:`/theory/foundations/operator_algebra`; every diffusion loss leaf emits and the
       resolvent runs assembled (bit-identical). See
       :ref:`loss-rep-three-modes`.
     - #272 #284
     - merged @ ``b058083e``
   * - 2026-07-03
     - **The unified k-estimator law: the reported :math:`k` IS the
       eigenvalue of the fixed-source map every method scales only
       fission by :math:`1/k` through (#259 root / #291 SN symptom).**
       SN's :meth:`~orpheus.sn.solver.SNSolver.compute_keff` gained the
       #291 boundary-leakage term :math:`L` (the
       :meth:`AngularBoundaryFlux.net_current
       <orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux.net_current>`
       :math:`\sum_m(\Omega_m\cdot\hat n)w_m\psi_m` face-integral, a
       structural zero on reflective faces so lattice anchors stay
       bit-identical) and the R7 :math:`(n,2n)` convention (emission on
       the removal side, ``absorption_xs`` counts the collision once) —
       see :ref:`sn-keff-estimator` (:eq:`sn-keff-update`).  MoC flipped
       to the same net-removal spelling (:eq:`moc-keff-update`;
       principled re-baseline 1.125 → 1.25 on the L0 case); CP was
       already the operator-consistent member (:eq:`cp-keff-update`).
       The ``KEigenvalue`` estimator-**injection** seam (``keff_estimator``
       / ``production_estimator`` kwargs, ``_default_*`` functions)
       retired (R8): the estimators are hardwired methods, dead by design
       because every estimator consistent with the posed problem agrees
       at a converged eigenpair.
     - #259 #291
     - merged @ ``a4952c3e``
   * - 2026-07-03
     - **The ``TransportMethod`` Protocol minted over the method-meshes;
       BC resolution unified into ONE shared body; the rank-N walker
       moved to geometry; the Wave-5 realizer registry dissolved
       (#290 P7b)** — the two recorded witnesses (the homogenization
       method-layer note on
       :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh` and
       the boundary-realizer seam) are discharged by the structural
       :class:`~orpheus.transport.method.TransportMethod` Protocol over
       ``SNMesh`` / ``DiffusionMesh`` (instance surface only; neither
       mesh imports it — conformance is checked where each calls the
       shared body). The twin ``_resolve_bcs`` loops collapse into
       :func:`~orpheus.transport.method.resolve_boundary_conditions`
       (face loop + reflective default + tag → law parse, method-generic
       over each mesh's ``BOUNDARY_OPERATOR_REGISTRY``); each mesh keeps
       only its ``realize_boundary_law`` arm, and ``SNMesh`` builds its
       unified trace in the construction body (the trace-inside-BC-
       resolution wart retired). ``realize_recursively`` moved to its
       method-blind home :mod:`orpheus.geometry.boundary` (leaf realizer
       REQUIRED); the ``_BoundBoundaryOperator`` kind tag now reads the
       law's own registry key. ``BoundaryRealizerRegistry`` + the
       MoC/MC/CP ``NotImplementedError`` stub realizers deleted — no
       consumer ever resolved a realizer by name (you hold the
       method-mesh → you have its realizer), and the import-side-effect
       registration hazard class dies with the indirection. See
       :ref:`bc-dual-registry` and the walker-placement retrospective at
       :ref:`bc-realize-recursively`.
     - #290
     - ``feature/diffusion-integration`` (P7b)
   * - 2026-07-03
     - **The Morel–Montry unbound legacy mode retired — every
       pole-angular closure is mesh-bound by construction; the sweep
       output mode split into types (pyright burn-down C5)** — the
       ``MorelMontryAngularSweep(sn_mesh=None)`` test-compatibility mode
       (and the ``| None`` widenings it forced on the whole
       :class:`~orpheus.sn.angular.closure.PoleAngularClosureBase`
       state contract, plus four runtime "unbound" guards) is deleted:
       ``sn_mesh`` is REQUIRED, the family's ``cls(sn_mesh)`` construction
       contract is total, and the pure-algebra recurrence surface moved
       back to module level
       (:func:`~orpheus.sn.angular.closure.compute_psi_half_per_level`
       — hand-built-coefficient verification needs no instance).  The
       ``SNMesh`` closure override became a **class** parameter
       (``pole_angular_closure: type[PoleAngularClosureBase]``) — an
       instance could only ever be unbound or foreign-bound, since the
       mesh it must bind to does not exist yet (Pattern 4).  The matvec's
       dead ``is None`` closure fallback and four stale ``type: ignore``
       comments retired with it.  Companion split: the solve-direction
       output DI ``_SweepEmit`` became a closed type family
       (``_SweepEmitAngular`` / ``_SweepEmitMoment`` with REQUIRED
       buffers, polymorphic ``pure_z``), the sweep chain's return became
       the honest mode-keyed pair ``(angular, scalar) | (moments,
       None)``, and the 1-D scan walk binds each geometry arm's face
       views inside its own arm (no cross-arm Optionals).  Landed the SN
       tree at ZERO pyright errors.
     - #226
     - ``refactor/pyright-burndown`` (C5)
   * - in dev
       (2026-07-02)
     - **The (spatial ⊗ angular) product — closure-owned :math:`\tau`,
       the pairing-validity surface, and the space–angle separability
       capstone (#236 Phases 1–3); the ``PoleAngularClosure`` Protocol
       retired to a single ABC (#248 / #249)** — the S\ :sub:`N` sweep's
       two curvilinear discretisation axes (the spatial cell update and
       the angular closure) are made a genuine tensor product in the
       type system. **Phase 1** adds the pairing-validity surface — the
       ``diffusion_limit_consistent`` / ``supports_curvilinear``
       predicates gate a (cell-update, angular-closure) pair on *both*
       single-axis diffusion-limit conditions (LMM-1987 spatial ×
       BMC-2010 angular) — and makes the curvilinear
       linear-discontinuous selection honest. **Phase 2** moves the
       Morel–Montry weight :math:`\tau` off the geometry-owned
       ``StreamingTerms`` onto the angular closure that owns it
       (:attr:`~orpheus.sn.angular.closure.PoleAngularClosureBase.tau_per_ordinate`);
       :math:`\tau` and the derived redistribution constants
       :math:`c_{\rm in}` / :math:`c_{\rm out}` travel to the stateless
       diamond scheme as ``CellVisit``
       *data* (never a closure dependency), stamped at the single
       production site ``SNMesh._make_cell_visit`` (both the three
       stamped fields and the stamping method were retired at P4.9a,
       2026-08-28 — :ref:`sn-p49a-closure-owns-the-march`);
       **Step C** deletes the parallel geometry-side :math:`\tau`
       producer, leaving ``StreamingTerms`` with no closure field (the
       stronger "purely geometric" reading was refuted 2026-08-28 — the
       packet keeps its direction cosines and :math:`\Delta A/w`, and
       moved to ``transport/spatial/scheme.py`` at P4.3; see
       :ref:`sn-tau-c-on-cellvisit-live`). **Phase 3** characterises the
       error surface — Cartesian **separates** (:math:`E \approx E_{\rm
       space} + E_{\rm angle}`), curvilinear **gates** (:math:`E \approx
       \max(E_{\rm space}, E_{\rm angle})`) — pinned by the L1 MMS gate
       :eq:`sn-space-angle-separability` (see
       :ref:`sn-space-angle-separability-section`). **#248 / #249**
       retire the orphaned ``PoleAngularClosure`` Protocol and its dead
       ``__call__`` bundle, hoisting the three strategy methods to the
       sole
       :class:`~orpheus.sn.angular.closure.PoleAngularClosureBase`
       ABC (the L2 single-source contract; see
       :ref:`sn-pole-angular-closure-protocol`).
     - #236 / #248 / #249
     - merged @ ``607b548``
       ``feature/sn-spatial-angular-restage``
       (``82cd210``, ``cdb3cd1``, ``9b93db7``, ``4f9e9b3``,
       ``81422f4``)
   * - in dev
       (2026-07-02)
     - **Operator-inverse taxonomy step 6 tail (carve P5) —
       composition-algebra return types pinned statically; the
       ``inverter`` posing narrative refreshed; #280 redesigned onto
       ``A.H.inverse()`` (#226)** — the composition surfaces now carry
       **precise static return types**, closing the taxonomy's
       type-legibility charter:
       :meth:`ScaledOperator.inverse() <orpheus.numerics.operator.ScaledOperator.inverse>`
       is parametrised on the **swapped carriers**
       ``ScaledOperator[Codomain, Domain]`` (an inverse maps the
       forward's codomain back to its domain), joining the sum / scaled
       / product / ``.H``-swap annotations. A **pyright-only pin bank**
       ``_composition_algebra_return_type_static_pins``
       (``tests/sn/operators/test_operators_apply_typed.py``)
       ``assert_type``-pins what *every* composition surface returns —
       sums, scaled, products, the ``.H`` adjoint carrier-swap, and the
       **algebra-closed vs wrap-delegate** inverse kinds (permutation /
       identity / scaled / tensor invert into their own forward type;
       product / diagonal return the wrap-delegate
       :class:`~orpheus.numerics.operator.InverseOperator`) — with
       **M-SCALED-BARE** ``reportAssertTypeFailure`` teeth:
       un-parametrising the ``ScaledOperator`` swap to a bare
       ``LinearOperator`` reddens the matching pin. The ``inverter``
       posing narrative on this page was **refreshed to the
       pre-inversion model** — the drivers no longer take an
       ``inverter`` callable
       (:class:`~orpheus.numerics.iteration.SourceIteration` applies a
       pre-inverted ``A_inv``;
       :class:`~orpheus.numerics.iteration.KrylovAcceleration`'s
       ``inverter`` → ``preconditioner``), the surviving "caller
       controls :math:`A^{-1}`" principle recast as an
       inverse-operator-family **type choice**
       (:ref:`choosing-inverse-realisation`). And **GitHub #280** (the
       deferred adjoint-inverse) was **redesigned** onto the now-live
       ``A.H.inverse()`` spelling — the swap law
       :math:`(A^{H})^{-1} = (A^{-1})^{H}` with
       ``SweepOperator.apply_transpose`` (the reverse-scan
       solve-transpose) as its concrete deliverable.
     - #226 / #280
     - merged @ ``1729647``
       ``refactor/inverse-as-operator``
   * - in dev
       (2026-07-02)
     - **Operator-inverse taxonomy step 5b — homogeneous spells the full
       resolvent in the operator algebra; MatrixInverseOperator's first
       production consumer (#226)** —
       :func:`~orpheus.homogeneous.solver.solve_homogeneous_infinite`
       re-spelled onto ``K = MatrixInverseOperator(loss, basis_shape=(ng,1))
       @ production``, the eigenpair extracted by the NEW shared
       Perron–Frobenius primitive
       :func:`~orpheus.numerics.eigenvalue.dominant_eigenpair` — the ONE
       home of the argmax-real selection, complex-dominant rejection, and
       ``φ.sum() ≥ 0`` sign convention
       (:func:`~orpheus.numerics.eigenvalue.direct_eigenvalue` now
       delegates and keeps zero production consumers as the
       ``(A, F)``-posed sibling engine; RQI's inline sign-flip folded onto
       the same ``_sign_normalised``). The explicit
       ``MatrixInverseOperator(loss)`` construction is the
       **strategy-choice-as-type** seam realized: the structure-keyed
       ``loss.inverse()`` would return the ITERATIVE Green splitting — the
       exactly-solvable 0-D problem earns the dense direct realization.
       ``_assemble_loss_matrix`` → ``_assemble_loss_operator`` (returns the
       UN-materialized ``C − K_iso``; MatrixInverseOperator's eager-LU ctor
       densifies). Re-baseline **principled-equivalence** (batched ``gesv``
       → held ``lu_factor`` + per-column ``lu_solve``; measured
       bit-identical on this host, gated rtol=1e-12). V&V finding: a
       factor swap ``F·A⁻¹`` and a resolvent transpose are **spectrally
       invisible** (similar matrices ⇒ :math:`|\Delta k| = 0` exactly) —
       every k-level gate is blind to them; the committed catcher is the
       object-level ``K.as_matrix() ≡ solve(A, F)`` intrinsic gate. See
       :ref:`spectral-invisibility` (:doc:`/theory/foundations/infinite_medium`).
     - #226
     - merged @ ``1729647``
       ``refactor/inverse-as-operator``
   * - in dev
       (2026-07-02)
     - **Operator-inverse taxonomy step 6 — the capability frozenset retired,
       both axes (#226, carve P4, W1+W2)** — the stringly-typed
       ``capabilities: frozenset[str]`` advertisement (``CAP_APPLY`` /
       ``CAP_SOLVE`` / ``CAP_APPLY_TRANSPOSE`` + ``MissingCapability`` + both
       ``_has`` copies) is **retired from every operator** — leaves, composers,
       aggregators, shims. Each axis (inverse, adjoint) now carries the
       **three-layer structural surface** (user-locked "Design C + B"): a
       runtime **predicate**
       (:attr:`~orpheus.numerics.operator.LinearOperator.is_invertible` /
       :attr:`~orpheus.numerics.operator.LinearOperator.is_adjointable`,
       value- and structure-aware, recursive on composites); an
       **operator-returning method** (``inverse()`` per-class /
       :attr:`~orpheus.numerics.operator.LinearOperator.H` on the base); and a
       **realization verb** (``solve`` / ``apply_transpose``) present only
       where a native realization exists. **Design C** resolves the
       false dichotomy between structural and value-dependent
       non-invertibility: a *structural* leaf (``ZeroOperator``, masks, source
       dyads, ``L`` / ``S`` / ``F`` / ``B``) declares **no** ``inverse()`` —
       misuse is a *static* pyright error — while a *value-dependent* operator
       (zero-coefficient
       :class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`,
       non-invertible-headed
       :class:`~orpheus.numerics.operator.OperatorSum`) declares it and raises
       :class:`~orpheus.numerics.operator.NotInvertible` **eagerly**. The
       runtime→static bridge is a pair of PEP-647 ``TypeGuard`` free functions
       :func:`~orpheus.numerics.operator.invertible` /
       :func:`~orpheus.numerics.operator.adjointable` (``TypeGuard`` NOT
       ``TypeIs`` — the value-dependent predicate makes only the one-directional
       promise); the runtime check IS the static permission, so all four
       ``cast(SupportsInverse, …)`` sites and all ten ``solve`` /
       ``apply_transpose`` ``# type: ignore``\ s **deleted**. **Design B**
       prunes ``solve`` to native realizations (deleted on Sum / Identity /
       Permutation / Scaled / TensorProduct / the boundary shim; kept on
       Diagonal / Multiplication / the sweep composites / the mixin / Product —
       whose ``solve`` **re-routes** through each factor's ``.inverse().apply``,
       bit-identical per kind and total over the solve-retired kinds).
       ``MissingCapability`` split into two ``TypeError`` successors —
       :class:`~orpheus.numerics.operator.NotInvertible` (inverse axis) /
       :class:`~orpheus.numerics.operator.MissingAdjoint` (adjoint axis). The
       **one behavior change** (spec §38): ``.H`` on a non-adjointable operator
       raises :class:`~orpheus.numerics.operator.MissingAdjoint` **eagerly at
       construction**, not lazily at the first ``.apply``. Two real production
       bugs the net caught and fixed root-cause: ``_seeded_inverse`` crashing on
       algebra-closed preconditioner heads (fixed with a two-kinds dispatch),
       and :class:`~orpheus.numerics.operator.RankOneOperator` losing its only
       adjointability advertisement (restored via an ``is_adjointable``
       override). Verification: **keystone v2** (``tests/_harness/predicates.py``
       rewritten to the permanent two-axis / three-valued contract; the
       coexistence scaffold ``assert_capability_faithful`` **deleted**, its
       licensing job done) is the standing faithfulness gate, plus the §40
       ``OperatorProduct.solve`` re-route gates (Mode-11 sentinel + dense
       anchor + five-row factor-kind matrix vs a pre-carve baseline) and the
       full 127-read / 36-file migration with a ZERO-hit completeness re-grep.
       Full tier ``-O`` serial **3853 / 0**; pyright ratchet re-baselined DOWN
       **148 → 145**. See :ref:`capability-set-semantics` (:doc:`/theory/foundations/operator_algebra`).
     - #226
     - merged @ ``1729647``
       ``f4919b1``
       ``refactor/inverse-as-operator``
   * - in dev
       (2026-07-02)
     - **Operator-inverse taxonomy step 5 — the materialising functor and
       the dense direct inverse; #285 closed for products (#226)** — adds the
       *materialising* family. :meth:`LinearOperator.as_matrix <orpheus.numerics.operator.LinearOperator.as_matrix>`
       is promoted to a **universal base method** — the functor **out** of
       the operator category (:math:`\mathrm{Op}\to\mathrm{Mat}`, the
       taxonomy §2 fourth arrow, NOT an endofunctor), realized as the
       apply-to-basis loop lifted from the homogeneous solver's retired
       ``_as_dense`` (C-order columns, rectangular-honest, dense
       :class:`numpy.ndarray` return). Its size gate
       :class:`~orpheus.numerics.operator.MatrixTooLarge` is a **RuntimeError**
       — a *resource* effect on a **total** functor (§17 A2: every operator
       HAS a matrix), hence **no** ``is_materializable`` predicate, and it is
       class-distinct from the ``ValueError`` an un-derivable basis raises
       (§27.C). :class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`
       is the **4th** :class:`~orpheus.numerics.operator.InverseWrapMixin`
       sibling — the inverse of a **structureless small** operator (eager
       ``lu_factor(inner.as_matrix())`` + per-``apply`` back-solve), direct
       construction only (no ``.inverse()`` routes to it yet). Its **name is
       earned** at the **precision grain** (spec §27.A, which **supersedes**
       the §13 M-row parenthetical now that ``as_matrix`` is universal): an
       iterative Green *also* satisfies :math:`[G][A]\approx I`, but only to
       driver-tol — M-materialise + M-direct hold at **machine·cond**, with an
       explicit Green contrast proving the invariant *distinguishing*. The
       constructor reads **values, not structure** — it does **not** consult
       ``inner.is_invertible``; the witness is :math:`(-S_{\rm ao})+D`, which
       :class:`~orpheus.numerics.green_operator.GreenOperator` **refuses** at
       construction (leading term non-invertible) yet materializes to an
       invertible matrix that this class inverts (the §3 strategy-override
       seam as explicit construction, ndarray analog of the ``FullField``
       :math:`(-S)+(L+C)`, out of ``as_matrix``'s ndarray scope). **#285
       closed for products**:
       :meth:`OperatorProduct.inverse <orpheus.numerics.operator.OperatorProduct.inverse>`
       now returns :class:`~orpheus.numerics.operator.InverseOperator`
       (bit-identical action, canonical seeded ``apply``, object-identity
       involution) — the **two-kinds** taxonomy (wrap-delegate family vs
       algebra-closed Permutation/Identity/Scaled inverses that stay
       unwrapped); the mixin ``_ForwardT`` bound relaxed
       ``_InvertibleForward`` → ``_WrappedForward``. Homogeneous ``_as_dense``
       **retired** → both call sites on ``as_matrix(basis_shape=(ng,1))``
       (byte-identical :math:`k_\infty` / flux; the landed SymPy pins
       untouched); ``dense_per_material`` re-documented as the storage-side
       oracle (zero production consumers). Gates:
       ``tests/numerics/test_matrix_inverse_operator.py`` + extensions;
       **14 mutations verified**, pyright ratchet exactly 148. See
       :ref:`matrix-inverse-operator` (:doc:`/theory/foundations/operator_inverse_family`).
     - #226 / #285
     - merged @ ``1729647``
       ``refactor/inverse-as-operator``
   * - in dev
       (2026-07-02)
     - **Operator-inverse taxonomy step 4 — the Green operator (the first
       *iterative* inverse) + the wrap-delegate mixin (#226)** — the generic
       sum inverse.
       :meth:`OperatorSum.inverse <orpheus.numerics.operator.OperatorSum.inverse>`
       now returns a
       :class:`~orpheus.numerics.green_operator.GreenOperator` — the
       :math:`A`-preconditioned splitting :math:`(A-B)^{-1}=\sum_k
       (A^{-1}B)^k A^{-1}` (the multiple-scattering Neumann series) wrapping
       :class:`~orpheus.numerics.iteration.SourceIteration` (§11.2 — it
       re-implements no iteration math; the left-spine head becomes the
       preconditioner through its *own* structure-keyed ``.inverse()``, the
       remaining terms ride as negated gains). It is the **first** family
       member with **no legacy ``.solve`` to inherit** (the sum was not
       invertible before step 4), so its correctness rests on
       **structural-independence anchors only** (dense-LU + the G-Neumann
       expansion) — no bit-identity twin, and (an iterative sum inverse
       being neither an eigenvalue solver nor source-driven) **no
       eigenvalue and no MMS claim**. The name is **earned** by G-Neumann
       (the splitting a generic :math:`A^{-1}` cannot satisfy) + G-reciprocity
       (the **Euclidean** transposed-operand Green, no ``.H`` — the metric
       adjoint-inverse is the separate #280) + G-kernel (folded into the
       :math:`\delta_j` anchor). The ``OperatorSum`` contract changed:
       :attr:`is_invertible <orpheus.numerics.operator.OperatorSum.is_invertible>`
       → ``self.a.is_invertible`` ("leading-term-preconditionable at this
       spelling", spec §18.B) with lockstep ``CAP_SOLVE`` — flipping two
       frozen ``is_invertible is False`` pins to ``True`` by design. The
       **ordering ruling** landed with all four edges gated: ``L+C`` →
       :class:`~orpheus.sn.operators.sweep_operator.SweepOperator` (the MRO
       shadow); ``(L+C)−S`` → Green, converges; ``C+L`` → Green constructs
       then raises
       :class:`~orpheus.numerics.green_operator.ConvergenceFailure` **loudly**
       (never a silent wrong answer); ``(−S)+A`` → refused at construction
       naming the canonical order. The wrap-delegate back-half extracted into
       :class:`~orpheus.numerics.operator.InverseWrapMixin` at the **third**
       sibling (``_SolveBackedLeaf`` → ``_InvertibleForward``), whose abstract
       ``apply(x, /, *, initial_guess=None)`` **resolves #285 STRUCTURAL**
       (pyright rejects a kwarg-less override, ``ABCMeta`` blocks a missing
       one; residue = the composed ``OperatorProduct.inverse()`` at step 5).
       ``ConvergenceFailure`` reads the **TRUE** relative residual (not the
       driver's ρ-blind increment, Signature 9) and **drives** a refinement
       loop to meet it — a check-only design would false-raise for every
       :math:`\rho>1/2`. Gates:
       ``tests/numerics/test_green_operator.py`` +
       ``tests/sn/operators/test_green_operator_sn.py`` (het-2G vacuum slab
       with a trace-consistent manufactured anchor resolving the #284 source
       subspace); **14 mutations verified**. See :ref:`green-operator`
       (:doc:`/theory/foundations/operator_inverse_family`).
     - #226 / #285
     - merged @ ``1729647``
       ``refactor/inverse-as-operator``
   * - in dev
       (2026-07-01)
     - **Operator-inverse taxonomy step 3 — the solver builds the inverse,
       the driver applies it; the resolvent concept fully dissolved
       (#226)** — completes the steps-1–3 arc.
       :class:`~orpheus.numerics.iteration.SourceIteration` now consumes the
       inverse-application **operator** directly (first parameter ``A_inv``,
       the new static contract
       :class:`~orpheus.numerics.iteration.SupportsSeededApply`
       ``apply(rhs, *, initial_guess=None)``): the solver layer builds it
       once (``base_resolvent.inverse()``, or the windowed product
       ``P @ A.inverse()``) and the loop step is the unconditional
       ``A_inv.apply(rhs, initial_guess=psi_prev)``. The former
       ``inspect.signature`` seed-probe (``_solve_accepts_seed`` /
       ``_solve_with_seed``) and the constructor ``CAP_SOLVE`` gate are
       **deleted** — the invertibility obligation moved to the inverse
       *builder* (``.inverse()`` on a non-invertible leaf raises), so the
       driver checks ``CAP_APPLY`` only and an apply-only step operator (the
       coisometry-factored windowed product) is accepted by design. The
       transitional ``_MomentWindowedResolvent`` adapter dissolved;
       :func:`_maybe_window <orpheus.sn.solver._maybe_window>` is now the
       product factory. :class:`~orpheus.numerics.iteration.KEigenvalue`
       guards ``A.is_invertible`` at construction and builds the inner step
       via the single-source
       :func:`~orpheus.numerics.iteration.seeded_inverse` (spelled
       ``_seeded_inverse`` at the time; promoted to the public name by
       #276 A4);
       :class:`~orpheus.numerics.iteration.KrylovAcceleration` keeps the
       forward ``A`` and rewires its default preconditioner from a
       ``CAP_SOLVE`` probe to ``seeded_inverse(A).apply``. Whether
       seeded-apply becomes a structural mixin or stays per-leaf convention
       is `#285 <https://github.com/deOliveira-R/ORPHEUS/issues/285>`_
       (folded into steps 4–5). Gates:
       ``tests/sn/solve/test_seed_threading_spy.py`` (Mode-11 path spy,
       route-invariant across the rewire; M-SEED-DROP/ZERO/STALE + M-PROBE
       teeth) and
       ``test_2d_windowed_product_over_gauss_seidel_M_equals_post_projection``
       (the windowed×G-S corner). See :ref:`inverse-application-driver`
       (:doc:`/theory/foundations/operator_inverse_family`).
     - #226
     - merged @ ``1729647``
       ``refactor/inverse-as-operator``
   * - in dev
       (2026-07-01)
     - **Operator-inverse taxonomy step 2 — the G-S resolvent reified,
       the windowed path retyped (#226)** — the duck-typed
       ``_GaussSeidelResolvent`` (which paired
       ``apply``\ :math:`=(L+C)\psi` with
       ``solve``\ :math:`=(L+C-B_{\rm lower})^{-1}` — inverses of
       *different* operators, round-trip defect O(1) :math:`=2.667`)
       dissolves into the honest **matrix splitting**
       :math:`(L+C-B)=M-B_{\rm upper}` (this entry said *regular* matrix
       splitting until 2026-08-09; struck per #341 —
       :ref:`sn-boundary-gs-not-regular`).
       :math:`M=(L+C)-B_{\rm lower}` is
       reified as
       :class:`~orpheus.sn.operators.scheduled_invertible.ScheduledInvertibleOperator`
       (via :meth:`StreamingCollisionOperator.__sub__
       <orpheus.sn.operators.streaming.StreamingCollisionOperator.__sub__>`), with
       :math:`B_{\rm lower}` the schedule-split half
       (:meth:`SNBoundaryOperator.split
       <orpheus.sn.operators.boundary.SNBoundaryOperator.split>` /
       :meth:`SweepSchedule.lower_inflow_rows
       <orpheus.sn.loss_representation.sweep_schedule.SweepSchedule.lower_inflow_rows>`
       → :class:`~orpheus.sn.operators.boundary.SNMaskedBoundaryOperator`)
       and :math:`B_{\rm upper}` riding the SI driver as an ordinary gain
       — ``M.inverse()`` now round-trips ``M.apply`` at
       :math:`5.2\times10^{-16}` (bulk) / :math:`4.4\times10^{-16}`
       (trace).  In parallel the ``solve_moments`` output-mode method (a
       codomain change wearing a config) retires into the typed windowed
       composition ``P @ A.inverse()``
       (:class:`~orpheus.sn.operators.windowing.WindowedSweep` =
       :class:`~orpheus.sn.operators.windowing.BulkAnalysisOperator` block
       coisometry ``@`` the forward's inverse), whose fused ``apply`` ≡ the
       deforested oracle at :math:`1.8\times10^{-16}`.  Gates:
       ``tests/sn/solve/test_gauss_seidel_reification.py`` (W2 round-trip /
       split-exactness / FP-invariance + M-SPLIT-DIR / M-SPLIT-PART
       mutations) and ``test_2d_windowed_product_equals_post_projection``.
       See :ref:`si-gauss-seidel-reification` and
       :ref:`windowing-retyped`.
     - #226
     - merged @ ``1729647``
       ``refactor/inverse-as-operator``
   * - in dev
       (2026-06-28)
     - **SN scattering adjoint** :math:`S^{T}` **landed (#118 closed)** —
       :meth:`ScatteringOperator.apply_transpose
       <orpheus.transport.operators.scattering.ScatteringOperator.apply_transpose>`
       is now LIVE as :math:`(1/W)\,\mathrm{full\_scatter\_kernel}^{T}`, the
       harmonic-frame conjugation
       :math:`R\circ(\Lambda_{\ell\ge0}+N_{2n})\circ M` whose transpose
       :math:`M^{T}\circ(\Lambda+N_{2n})^{T}\circ R^{T}` falls out of
       :meth:`OperatorProduct.apply_transpose
       <orpheus.numerics.operator.OperatorProduct.apply_transpose>` for free;
       the operator advertises ``CAP_APPLY_TRANSPOSE`` and the
       "no ``apply_transpose``" confession is retired.  Because the forward
       keeps the scalar fast-path, the frame-form :math:`S^{T}` is
       Euclidean-reciprocity-pinned
       (:math:`\langle S\psi,\chi\rangle=\langle\psi,S^{T}\chi\rangle`)
       against the *independent* forward — a genuine cross-check.  This is
       the discrete adjoint the :math:`\psi^{*}` chain (adjoint-weighted
       homogenisation, perturbation theory) builds on; its forward companion
       (P2) routes the SN isotropic source through the model-shared
       :math:`K_\mathrm{iso}` operators (0-ULP bit-identical).  See
       :ref:`sn-scattering-adjoint`.
     - #118
     - merged @ ``15185e5``
       ``feature/sn-adjoint-transport``
   * - 2026-06-27
     - **Energy condensation landed (the energy-axis transpose of
       homogenization)** — :meth:`Solution.condense
       <orpheus.sn.solution.Solution.condense>` collapses fine-group
       cross sections onto a coarse :class:`~orpheus.data.energy_grid.EnergyGrid`,
       spectrum-weighted, returning **portable** few-group XS
       (``dict[int, Mixture]``, mesh-DECOUPLED — the asymmetry-law
       counterpart of ``homogenize``'s mesh-COUPLED ``MaterialMesh``).
       Reuses the *unchanged* Petrov-Galerkin frame
       (:meth:`frame.project <orpheus.numerics.frame.FrameBase.project>`,
       flux = test weight, counting measure). The production case
       (421-group → WIMS-69/172) is **non-nested**, so the partition is a
       fractional, partition-of-unity overlap table
       (:class:`~orpheus.numerics.basis.OverlapBasis`) with a selectable
       within-group flux model
       (:class:`~orpheus.data.energy_grid.WithinGroupSpectrum`, 1/E
       first) — the flux-weighted average is the rank-0 moment of
       Generalized Energy Condensation (:cite:`Rahnema2008`). Downsampling-only
       (the upscaling guard refuses a finer target). See
       :ref:`sn-energy-condensation`.
     - #274
     - ``e9a6a49`` → ``68ceb9a``
   * - 2026-06-26
     - **Rank-N boundary walker co-located** — ``realize_recursively``
       (the descriptor-tree → operator-tree walker) merged into
       :mod:`orpheus.sn.boundary.realizer` next to the
       :class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer` it
       dispatches to; the near-twin ``boundary_realize`` module retired.
       It stays honestly SN-specific; the method-agnostic generalization
       (registry-resolved leaf + ``MethodSpace`` Protocol, walker moves to
       ``geometry/boundary/``) is deferred to the second functional
       realizer — the **same** trigger that mints the ``TransportMethod``
       Protocol behind :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh`
       (two witnesses to one missing type). See
       :ref:`bc-rank-n-algebra`.
     - —
     - ``b0e5ba1`` → ``932e769``
   * - 2026-06-24
     - **Homogenization re-framed as Petrov-Galerkin** (unified
       Frame-projection campaign, P3) — the spatial-homogenization theory
       and production code are corrected from the forward-only
       ":math:`L^2(\phi V)`-Galerkin" reading to the honest
       **Petrov-Galerkin** framing: the flux is a *test-weighting the
       solution emits* (carried on an explicit
       :class:`~orpheus.numerics.basis.weighted_indicator_basis.WeightedIndicatorBasis`
       test side), **never folded into the measure**, which carries only
       the geometric volume + fixed :math:`L^2` metric. Homogenization
       becomes the coefficient extraction :math:`G^{-1}M`
       (:meth:`FrameBase.project <orpheus.numerics.frame.FrameBase.project>`)
       of a :class:`~orpheus.numerics.frame.PetrovGalerkinFrame`; the
       forward (:math:`\varphi^*=\varphi`) case is the *Galerkin
       degenerate* of the eigenvalue-consistent adjoint-weighted
       (:math:`\varphi^*\ne\varphi`) homogenization (P6). The whole XS
       field projects as one object
       (:meth:`MaterialXSField.project_through
       <orpheus.transport.mesh.material_xs_field.MaterialXSField.project_through>`)
       through *two* frames — reaction-rate-preserving :math:`\Sigma` and
       emission-rate-preserving :math:`\chi`. The 1-D guard is dropped
       (now n-D). See :ref:`sn-homogenization-petrov-galerkin-frame`.
     - #268
     - ``7e8ca2a`` → ``932e769``
   * - 2026-06
     - **Mesh + materials data/behavior split (``MaterialMesh``) + L2
       promotion** — the method-agnostic *mesh + materials* carrier
       :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh` is
       minted as the missing middle type between a bare geometry
       :class:`~orpheus.geometry.mesh.Mesh1D` (material *ids*, no XS) and
       a method phase space, and :class:`SNMesh` becomes
       ``SNMesh(MaterialMesh)`` — *data* (axes + materials + ``mat_map`` +
       volumes + ``ng``) on the base, *behavior* (quadrature + sweep
       stencil + boundary trace + closures) on the SN subclass.
       ``axis`` and ``material_xs_field`` are promoted out of
       ``orpheus.sn`` to L2 :mod:`orpheus.transport.mesh` (both are
       method-agnostic: pure coordinate geometry / mesh+materials reads).
       Cross-section **spatial homogenization** lands as a domain
       operation :meth:`Solution.homogenize
       <orpheus.sn.solution.Solution.homogenize>` (flux·volume-weighted,
       reaction-rate-preserving → :class:`MaterialMesh`), re-promoted to a
       solvable phase space by :meth:`SNMesh.from_material_mesh`. See
       :ref:`sn-spatial-homogenization`.
     - #267
     - ``5bcb1ce`` → ``932e769``
   * - in dev
     - **Spatial ⊗ angular product** — τ becomes closure-owned
       (``morel_montry_tau_per_level``); the geometry-τ producers and
       the ``StreamingTerms`` τ/α fields are deleted, and a separability
       MMS gate pins Cartesian-separates / curvilinear-couples. The
       space-angle tensor structure made explicit.
     - #236
     - merged @ ``607b548``
       ``feature/sn-spatial-angular-product``
   * - 2026-06
     - **Field-typed operator algebra + data XS layer** — cross sections
       become ``CoefficientField`` leaves; operators become *promotions*
       (``C = M[\sigma_t]``); the carrier is the timeless ``FullField``;
       the §5.6 Operator / Kernel / Functional taxonomy and the
       ``@overload`` "Pattern M" apply-fibration. Closes the typed-field
       campaign and pushes the data layer down to validated invariants
       (χ-simplex, production-weighted χ\ :sub:`mix`, XS-balance).
     - #257
     - ``505e1b7`` → ``f62b895``
   * - 2026-06
     - **Linear-Discontinuous on the DAG + principled polymorphism** —
       LD lands as a polymorphic discretization protocol (the
       coefficient model), the ``matvec_via_kernel`` favoritism is
       reverted, σ is single-sourced from the ``StreamingCollisionOperator``'s
       ``C``, and the unified all-d LD moment matvec + diffusion-limit
       closure ships.
     - #240 / #158
     - ``fde76ac`` → ``e74eafb``
   * - 2026-06
     - **Sweep / matvec re-layering — one walk** — the sweep strategy
       is carved into a first-class ``LossRepresentation``; the
       ScanMarch row-march twin makes *matvec ≡ sweep* literally one
       ``_OctantWalk`` traversal rather than two parallel paths.
     - #222
     - ``8913229`` → ``1b4b0c0``
   * - 2026-06
     - **3-D Cartesian admission** — axis-native ``from_axes`` admits
       *d = 3* end-to-end with no ``Mesh3D``; the constructor data-flow
       inverts so axes are primary, and windowing / G-S gates key on
       genuine dimensionality rather than the reduced-proxy.
     - #225
     - ``1da1e2f``
   * - 2026-06
     - **Foundation-cleanup cluster** — break the numerics→SN import
       cycle (moment-layout policy homed in numerics), retire the
       ``M_spatial`` / ``M_angular_redist`` operator-leaf split (the
       fused ``loss_action`` is the only path), and key the moment-frame
       involution on intent rather than coincidental shape.
     - #243 / #245 / #246 / #238
     - ``66b1bbd`` → ``dd4f542``
   * - 2026-06
     - **LD boundary slope source + transverse face-slope trace** — the
       fixed-source solver consumes a moment-resolved external slope
       source (Leg A) and the boundary trace carries the transverse
       face-slope moment (Leg B); both close the LM-1989 trap as
       structural-teeth vv Mode-10 cases.
     - #247 / #251
     - ``d9396a2`` / ``e5f2b1c``
   * - 2026-06
     - **Wave-O affine-typed operator algebra** — ``FluxDisplacement``
       gave the affine difference space :math:`V`; ``flux + flux``
       became a ``TypeError`` while ``flux − flux`` and the typed
       residual ``from_balance`` stayed legal, so :math:`(L+C-S-F)\psi = q`
       typed coherently and the residual was a typed defect.
       ⛔ The **affine** half was **OVERTURNED on 2026-08-19** — flux lives
       in the positive cone of a vector space, ``flux + flux`` is legal,
       and the displacement family retired
       (:ref:`cone-the-overturned-affine-design`). The **typed-residual**
       half stands unchanged.
     - #208 / #201
     - ``8c2f355`` / ``04e2859``
   * - 2026-06
     - **G-adjoint, reciprocity, and operator role-typing** — the
       analytic reverse-sweep ``StreamingOperator.apply_transpose``; the
       Hilbert-adjoint metric is owned by the ``FunctionSpace``
       (``apply_metric``); composite block-roles are *derived* from the
       operands, retiring the ``StreamingCollisionOperator`` FULL stamp.
     - #20
     - ``0efd233`` → ``7ccc14a``
   * - 2026-06
     - **Angular + face windowing** — the held SI iterate becomes
       ``HarmonicMomentFlux`` moments (angular window), the interior
       face cochain becomes a rolling ``_MovingFrontier`` (face window),
       and moments are accumulated in-sweep — eliminating the
       full-angular per-sweep transient (a measured ~3× peak-memory win).
       Full-field sweep/matvec oracles retained for the equivalence gate.
     - #205 / #218
     - ``b97d4f9`` → ``c7be111``
   * - 2026-06
     - **Honest variadic L+C−S−B driver** — the transitional ``S+B``
       fold is retired and boundary reflection ``B`` becomes a
       first-class coupling gain, so the within-group drivers generalize
       over a variadic operator list instead of a hard-coded fold.
     - #18
     - ``8563f4b`` → ``83a4ae6``
   * - 2026-06
     - **Source-iteration boundary Gauss-Seidel** — a polymorphic
       ``SweepSchedule`` (Jacobi / octant-group G-S) plus the
       ``_GaussSeidelResolvent`` and ``inner_schedule`` selector; a
       modest reflective-SI accelerator (the dominant scattering c-mode
       stays Krylov / DSA territory). Surfaced ERR-056 (diagonal-cubature
       shared-face fan-in).
     - #19
     - ``514ae21`` → ``a39905a``
   * - 2026-06
     - **1-D forward + transpose matvec carve** — the 1-D SN matvec is
       relocated off the operator into ``_OneDimScanWalk``, making
       *matvec ≡ sweep* a single-sourced code fact for the 1-D path.
     - #206
     - ``eaafbe1`` / ``7300f3e``
   * - 2026-06
     - **Curvilinear closure-seed fix (ERR-058) and eigenvalue
       verification** — coupled-pole spatial + half-angle angular
       closure seeds restore :math:`O(h^2)` for curvilinear sweeps, and
       the SI ≡ Krylov eigenvalue equivalence is pinned, closing the
       last ERR-026 manifestation.
     - #195 / #196
     - ``3b088ee`` → ``8609282``
   * - 2026-06
     - **2-D matvec through the DD closure** — the 2-D forward matvec is
       routed through the diamond-difference cell closure and the legacy
       FD stencil is retired, so the operator and the sweep share one
       discretization.
     - #208 (O.4b)
     - ``2288ea4``
   * - 2026-05
     - **Four-operator typed apply + typed-field foundation** — the
       ``F → S → C → L`` apply overload on ``AngularFlux``, founded on
       the typed ``AngularFlux`` / ``ScalarFlux`` / ``AngularBoundaryFlux``
       fields (``psi_bc`` retired). The birth of the typed-field
       architecture this whole history builds on.
     - #197
     - ``d8ddb03`` → ``eeac45f``
   * - 2026-05
     - **Depth-B field factoring** — ``ScalarFlux`` migrates to a pure
       ``Field`` subtype and the bare-``ndarray`` operator arms are
       retired (D-I / D-J), so typed-field dispatch becomes the *only*
       apply path through the operators.
     - #197 (Depth B)
     - ``c97897d`` → ``4a53737``

.. note::

   The four ``Phase N`` milestones (#18, #19, #20, #205) carry internal
   phase labels in their commit subjects rather than GitHub-issue
   trailers; the issue mapping above is from the SN development-sequence
   campaign record.  (#236 merged @ ``607b548``, 2026-07-02 — the former
   "pending merge" note here was a frozen-status claim git had already
   outdated.)


