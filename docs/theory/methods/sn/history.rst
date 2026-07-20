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
lives. Entries marked *(in development)* live on an unmerged feature
branch and have no landed hash yet.

.. list-table::
   :header-rows: 1
   :widths: 8 52 12 28

   * - When
     - Architectural milestone
     - Issue
     - Where
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
     - ``refactor/sn-walk-unification`` *(in development, 015dcc73)*
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
     - ``refactor/sn-walk-unification`` *(in development, 51c22396)*
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
     - ``refactor/sn-walk-unification`` *(in development, 0e03c304)*
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
       :ref:`sn-282-direct-starting-direction-solve`.
     - #280
     - ``refactor/sn-walk-unification`` *(in development, 98fe2e36)*
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
       :ref:`sn-282-direct-starting-direction-solve`.
     - #280
     - ``refactor/sn-walk-unification`` *(in development, e5d1acf)*
   * - 2026-07-11
     - **The within-group solve became block-native — one
       ``WithinGroupSystem`` record carrying the named** :math:`A = M - N`
       **splitting** (coupled-block campaign, Phase B.2d d1). The former
       ``_within_group_triple`` / ``_lagged_gains`` construction pair
       retired into a single builder,
       :func:`~orpheus.sn.coupled_system.build_within_group_system`, which
       returns the frozen
       :class:`~orpheus.sn.coupled_system.WithinGroupSystem` record — the
       loss grid together with its Hackbusch regular splitting
       :math:`A = M - N` (``resolvent`` = :math:`M`, the sweepable part
       inverted each step; ``gains`` = :math:`N`, the lagged couplings). On
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
     - ``refactor/sn-walk-unification`` *(in development, c0f23f6)*
   * - 2026-07-05
     - **The adjoint-inverse swap law wired the reverse-scan into the
       operator algebra** (#280 Phase 2.5c). With the reverse-scan in
       hand, the missing adjoint was surfaced through the landed #226
       inverse-as-operator taxonomy *without new solve machinery*:
       :class:`~orpheus.sn.operators.sweep_operator.SweepOperator`
       (:math:`(L+C)^{-1}`) gained
       :meth:`apply_transpose <orpheus.sn.operators.sweep_operator.SweepOperator.apply_transpose>`
       delegating to the inner's
       :meth:`~orpheus.sn.operators.streaming.InvertibleOperator.solve_transpose`
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
     - ``refactor/sn-walk-unification`` *(in development, 8cf5215)*
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
     - ``refactor/sn-walk-unification`` *(in development, f1ddeb6)*
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
     - ``refactor/sn-walk-unification`` *(in development, ba202a1)*
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
       :eq:`sn-282-anisotropic-source`). The augmented :math:`(L+C)` is
       block-lower-triangular in the seed-first walk order
       (:eq:`sn-282-block-triangular`), so the back edge is dead by
       construction (cold residual :math:`\to 2.5\times10^{-16}`; seed
       insensitivity :math:`4.57\times10^{-2}\to 0` bitwise). The whole
       ``PsiHalfAngleSeed`` strategy zoo retired (851 → 161 lines). The
       eigenvalue re-pose is principled by an angular-order N-sweep (a
       seed is a closure — h→0 is the wrong test). See
       :ref:`sn-282-direct-starting-direction-solve`.
     - #282 #280
     - ``refactor/sn-walk-unification`` *(in development, a29ab2d)*
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
     - ``refactor/spatial-promotion-assembly`` *(in development)*
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
     - ``refactor/k-estimator-unification`` *(in development)*
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
       :class:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase`
       state contract, plus four runtime "unbound" guards) is deleted:
       ``sn_mesh`` is REQUIRED, the family's ``cls(sn_mesh)`` construction
       contract is total, and the pure-algebra recurrence surface moved
       back to module level
       (:func:`~orpheus.sn.sweep.pole_angular_closure.compute_psi_half_per_level`
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
       (:attr:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.tau_per_ordinate`);
       :math:`\tau` and the derived redistribution constants
       :math:`c_{\rm in}` / :math:`c_{\rm out}` travel to the stateless
       diamond scheme as :class:`~orpheus.transport.spatial.scheme.CellVisit`
       *data* (never a closure dependency), stamped at the single
       production site
       :meth:`~orpheus.sn.mesh.augmented_mesh.SNMesh._make_cell_visit`;
       **Step C** deletes the parallel geometry-side :math:`\tau`
       producer so ``StreamingTerms`` is now purely geometric (see
       :ref:`sn-tau-c-on-cellvisit-live`). **Phase 3** characterises the
       error surface — Cartesian **separates** (:math:`E \approx E_{\rm
       space} + E_{\rm angle}`), curvilinear **gates** (:math:`E \approx
       \max(E_{\rm space}, E_{\rm angle})`) — pinned by the L1 MMS gate
       :eq:`sn-space-angle-separability` (see
       :ref:`sn-space-angle-separability-section`). **#248 / #249**
       retire the orphaned ``PoleAngularClosure`` Protocol and its dead
       ``__call__`` bundle, hoisting the three strategy methods to the
       sole
       :class:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase`
       ABC (the L2 single-source contract; see
       :ref:`sn-pole-angular-closure-protocol`).
     - #236 / #248 / #249
     - *(in development)*
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
     - *(in development)*
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
     - *(in development)*
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
     - *(in development)*
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
     - *(in development)*
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
     - *(in development)*
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
       :func:`_seeded_inverse <orpheus.numerics.iteration._seeded_inverse>`;
       :class:`~orpheus.numerics.iteration.KrylovAcceleration` keeps the
       forward ``A`` and rewires its default preconditioner from a
       ``CAP_SOLVE`` probe to ``_seeded_inverse(A).apply``. Whether
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
     - *(in development)*
       ``refactor/inverse-as-operator``
   * - in dev
       (2026-07-01)
     - **Operator-inverse taxonomy step 2 — the G-S resolvent reified,
       the windowed path retyped (#226)** — the duck-typed
       ``_GaussSeidelResolvent`` (which paired
       ``apply``\ :math:`=(L+C)\psi` with
       ``solve``\ :math:`=(L+C-B_{\rm lower})^{-1}` — inverses of
       *different* operators, round-trip defect O(1) :math:`=2.667`)
       dissolves into the honest **regular matrix splitting**
       :math:`(L+C-B)=M-B_{\rm upper}`.  :math:`M=(L+C)-B_{\rm lower}` is
       reified as
       :class:`~orpheus.sn.operators.scheduled_invertible.ScheduledInvertibleOperator`
       (via :meth:`InvertibleOperator.__sub__
       <orpheus.sn.operators.streaming.InvertibleOperator.__sub__>`), with
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
     - *(in development)*
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
     - *(in development)*
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
       Generalized Energy Condensation ([Rahnema2008]_). Downsampling-only
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
     - *(in development)*
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
       reverted, σ is single-sourced from the ``InvertibleOperator``'s
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
       gives the affine difference space :math:`V`; ``flux + flux``
       becomes a ``TypeError`` while ``flux − flux`` and the typed
       residual ``from_balance`` stay legal, so :math:`(L+C-S-F)\psi = q`
       types coherently and the residual is a typed defect.
     - #208 / #201
     - ``8c2f355`` / ``04e2859``
   * - 2026-06
     - **G-adjoint, reciprocity, and operator role-typing** — the
       analytic reverse-sweep ``StreamingOperator.apply_transpose``; the
       Hilbert-adjoint metric is owned by the ``FunctionSpace``
       (``apply_metric``); composite block-roles are *derived* from the
       operands, retiring the ``InvertibleOperator`` FULL stamp.
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
   campaign record. Issue #236 is the only entry not yet on ``main`` — it
   lives on ``feature/sn-spatial-angular-product`` pending merge.


