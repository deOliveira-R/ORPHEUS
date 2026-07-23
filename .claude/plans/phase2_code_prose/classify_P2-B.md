# P2-B classification — Loss Representation prose

## orpheus/sn/loss_representation/__init__.py

| lines | symbol | content-id | verdict | destination / anchor | conf |
|---|---|---|---|---|---|
| 1–125 | MODULE | Selectable representations; hierarchy; governing principle; history | **SPLIT** | — | — |
| 1–52 | (1a) | Solve & matvec algorithms; 1-D scan, multi-D wavefront; two buffer policies | **CONTRACT** | lines 1–52 encode type hierarchy & selection vocab; stays | H |
| 53–75 | (1b) | Governing principle: construct general, select narrow, specialize measured | **TWIN** | loss_representation.rst `loss-rep-selection` anchor + Fork-B2 evidence section | H |
| 76–95 | (1c) | Compatibility signal, three consumers (UI, factory, guard), SSOT | **CONTRACT** | needed at point of use for builder call & `supports()` predicate | H |
| 97–124 | (1d) | Carve history S0–S6.9 narrative (now moved to loss_rep-history theory anchor) | **HISTORY** | loss_representation.rst `loss-rep-history` (§2709); plan at .claude/plans/sn_sweep_strategy.md | H |
| 178–191 | frame_signs_for() | Module-level binding of octant signs + scheme to single-source involution #240 D5b-S3 | **CONTRACT** | dispatch & diffusion-limit root cause; design decision (Pattern 2 binding-dup retire) | H |
| 244–258 | _curvilinear_capability() | (scheme × geometry) gate; rejects unpublished LD curvilinear closure #158/#6 | **CONTRACT** | needed at SELECTION time, guards construction; error message load-bearing | H |
| 270–280 | LossRepresentation Protocol | Protocol docstring (10 lines): sweep, loss_action, supports methods | **CONTRACT** | method signatures, return modes, shape conventions | H |
| 291–311 | sweep() | Transport solve :math:`\psi = (L+C)^{-1} q`; angular/moment output; schedule/reflect; carrying-mesh semantics (R12a) | **CONTRACT** | call contract; output modes (angular XOR moments); ψ½ joint march routing | H |
| 319–335 | sweep_transpose() | Reverse-scan :math:`(L+C)^{-\mathsf T}`; #280 Phase 2.5b kernel-pair contract; deferral raise (multi-D Cartesian, LD moment-tail) | **CONTRACT** | method signature & deferral scope; step 6 carrying-mesh block structure | H |
| 340–361 | loss_action() | Matvec :math:`(L+C)\psi`; returns FULL loss, not :math:`L\psi`; Resolution A (σ diagonal once); #240 Step B σ passed explicitly | **CONTRACT** | matvec semantics critical (double-subtract bug if wrong); Resolution A invariant; carrying-mesh block scope | H |
| 364–382 | streaming_action() | Pure σ-free :math:`L\psi = \Omega \cdot \nabla\psi`; affine in σ; curvilinear Carlson seed σ-dependence cancels; single-sourced via loss_action(σ=0) | **TWIN** | loss_representation.rst `loss-rep-removal-form-matvec` (affinity + σ-free decomposition); coding-elegance Pattern 3 | M |
| 387–405 | loss_action_transpose() | Adjoint :math:`(L+C)^{\mathsf T}\phi`; returns FULL adjoint; Resolution A; #240 Step B; deferral raise (multi-D Cartesian, LD) | **CONTRACT** | matvec adjoint semantics; Resolution A invariant; deferral scope (O.2b deferred) | H |
| 421–434 | has_transpose_walk | ORIENTATION factor of adjoint-reachability (#280 Phase 2.5a); scheme + representation predicates; eager refuse at construction | **CONTRACT** | predicate for `is_adjointable`; tied to MissingAdjoint exception logic | H |
| 456–475 | (comment: pure-L streaming primitive) | σ-free leaf single-sourced through loss_action at σ=0; affine decomposition; one walk, no twin (#257 S8b) | **COMMENT-cut** | teaching narrative of Pattern 2 & affine structure; moves to loss_representation.rst or coding-elegance examples | L |
| 472–483 | (comment: abstract signatures TYPE_CHECKING) | TYPE_CHECKING block to avoid reportRedeclaration; subclass overrides via MRO | **COMMENT-cut** | runtime dispatch explanation; pedagogical only; remove or collapse to inline | L |
| 553–565 | _moment_frame_signs() | Octant sweep⇄global moment-frame sign vector (involution); #240 D5b-S3 root cause (diffusion-limit); per-axis read single-source binding | **CONTRACT** | needed for cell-kernel frame dispatch; diffusion-limit guard; part of operator contract | H |
| 569–581 | _spatial_moment_tail | Trailing per-CELL spatial-moment axis :math:`(\text{per\_axis})^d`; LD → yes, DD → no; face_moment_tail single-source | **CONTRACT** | shape convention for multi-moment closure; tied to _CellSolve matvec machinery | H |
| 587–616 | _inflow_to_moments() | Per-face domain inflow → 2^{d-1}-moment object; three cases (single-moment id, scalar-to-moment widen, moment-resolved pass-through) | **CONTRACT** | boundary-trace protocol; moment-resolved relaxation guard (#251 validat ion); shapes & broadcasting rules | H |
| 705–715 | _assert_face_moment_width() | Guard moment-resolved face width == 2^{d-1}; silent mis-broadcast hazard; shared by _inflow_to_moments + FullFieldWavefront (#251) | **CONTRACT** | guard & error message (unambiguous ValueError); shared by two cochain entry points | H |
| 750–762 | _ApplyOperands | Matvec problem data (probe ψ̄, σ_t, streaming tuple, Q_zero); positional-by-axis convention; kernel tuple layout | **CONTRACT** | dataclass signature & axis correspondence; passed to interior kernels | H |
| 771–784 | _SolveOperands | Sweep problem data mirror (Q source, σ_t, streaming tuple); positional-by-axis convention | **CONTRACT** | dataclass signature; passed to interior kernels | H |
| 788–806 | _moment_broadcast_sigma() | σ_t reshaped for trailing spatial-moment axis broadcast; DD identity, LD length-1 axis; single-source for sweep ≡ matvec (L21) | **CONTRACT** | broadcasting rule for dual-direction L21 invariant; shared invoke in sweep & matvec arms | H |
| 811–836 | _SweepEmit | Solve-direction OUTPUT mode (angular field XOR harmonic moments); Phase 5c polymorphic split (angular vs moment); "mixed output unrepresentable" | **CONTRACT** | output-mode types; pure_z polymorphic behavior; frozen-dataclass anti-pattern | H |
| 879–911 | _OctantWalk | THE in-plane octant traversal; sweep & matvec traverse same octant; forks on cell kernel + emit policy, never boolean flag | **TWIN** | loss_representation.rst `loss-rep-four` (wavefront family architecture, lower-triangular frame, one-walk theorem) + cartesian_multid.rst | H |
| 923–935 | _interior_walk() | THE shared octant frame (one-walk seam S6.4); per-octant: project signs, dispatch pure-z, derive faces, read inflow, run interior, shed outflow | **TWIN** | loss_representation.rst `loss-rep-one-walk-one-instance` (one-walk theorem); cartesian_multid.rst wavefront walk decomposition | H |
| 968–994 | sweep_group() | SOLVE-direction frame for ONE octant group (S6.4 b); Jacobi/Gauss-Seidel at level up (schedule loop); bare sweep, blind to BC | **CONTRACT** | supplies operands/emit/interior to _interior_walk; schedule-loop / inter-group reflect routing; carrying-mesh contract (step 6) | H |
| 1025–1049 | loss_action() in _OctantWalk | APPLY-direction frame (S6.4 a); multiline docstring **[LARGE PROSE]** — see separate analysis | **SPLIT** | — | — |


#### 1025–1049: _OctantWalk.loss_action() [LARGE PROSE ANALYSIS]

| section | content-id | verdict | destination / anchor | conf |
|---|---|---|---|---|
| 1026–1033 | Owns everything matvec frames duplicated; probe/accumulator setup, pure-z, inflow, outflow, boundary | **CONTRACT** | S6.4a frame ownership; carried responsibilities; kernel contract interface | H |
| 1034–1045 | Interior kernel signature; boundary semantics (BARE O.4b); reads from given trace, captures OUTFLOW, active-trace residual | **CONTRACT** | kernel interface (operands, oct_idx, signs_eff, inflow) → (LpC_octant, capture); boundary trace protocol (no BC apply) | H |
| 1046–1049 | Returns FULL (L+C)ψ, not Lψ; Resolution A; operator subtracts Σ_t once | **CONTRACT** | matvec semantics (double-subtract invariant); critical to operator correctness | H |

| 1054–1073 | Probe setup, moment_tail from probe space (iterate = SSOT), operands assembly, streaming tuple per-axis | **CONTRACT** | data layout conventions; moment-valued matvec (#240 D5b-S3); probe-space SSOT | H |
| 1075–1091 | (L+C)·ψ̄ accumulator; pure-z: σ·ψ̄ for degenerate ordinates; moment broadcast via shared helper (L21 sweep twin) | **CONTRACT** | pure-z semantics; shared helper (Pattern 2 single-source); qa Concern A blocker fix (moment-valued + pure-z guard) | H |

| 1255–1277 | _DAGWavefront | Base for two buffer policies (full field, rolling frontier); share per-octant DAG (family-owned, shape-cached); both deterministic via C3.2b oracle | **TWIN** | loss_representation.rst `loss-rep-four` (wavefront family, two buffer policies); cartesian_multid.rst buffer-policy architecture | M |
| 1298–1310 | loss_action_transpose() in _DAGWavefront | Multi-D Cartesian adjoint DEFERRED; shared deferral #280 kernel-pair contract; raises NotImplementedError | **CONTRACT** | deferral scope (O.2b later); matches base-class loud backstop; prevents silent wrong answer | H |
| 1317–1337 | MovingFrontierWindow | Wavefront sweep — rolling (d−1)-frontier window; production optimization of full-field oracle; two buffer policies, one octant order | **TWIN** | loss_representation.rst `loss-rep-four` (rolling frontier window, (d-1)-dim slab, performance evidence); cartesian_multid.rst frontier walk | H |
| 1374–1385 | _sweep_interior() in MovingFrontierWindow | Rolling-frontier interior kernel, SOLVE direction; per-row pure-z + interior kernels | **CONTRACT** | kernel signature for sweep_group(); shape conventions for windowed walk | L |
| 1386–1403 | (comment: Angular mode allocates per-octant buffer) | Angular mode per-octant allocation; moment mode accumulated across octants | **COMMENT-cut** | implementation detail (buffer layout); move to _SweepEmitAngular / _SweepEmitMoment docstring or retire | L |
| 1434–1443 | (comment: domain-edge capture carries trailing 2^{d-1} moment axis) | Face cochain width (2^{d-1} transverse moments) for multi-moment closure | **COMMENT-cut** | shapes & conventions; subsume into main docstring or face-cochain spec | L |
| 1448–1462 | loss_action() in MovingFrontierWindow | 2-D Cartesian forward loss action via rolling frontier; interior kernel, boundary residual assembly (O.4b) | **CONTRACT** | matvec implementation for rolling window; calls _loss_action_interior; interior kernel signature | H |
| 1473–1483 | _loss_action_interior() in MovingFrontierWindow | Rolling-frontier interior kernel, APPLY direction; per-row matvec kernels | **CONTRACT** | kernel signature for loss_action; shape conventions for windowed walk | L |
| 1482–1491 | (comment: unified moment matvec – multi-moment closure carries trailing axis) | Moment-valued matvec (#240 D5b-S3); LD trailing axis, DD byte-identical | **COMMENT-cut** | teaching narrative of moment architecture; moves to #240 section or swept into main prose | L |
| 1505–1514 | (comment: domain-edge capture carries 2^{d-1} transverse moments) | Transverse-moment axis on face cochain for multi-moment closure | **COMMENT-cut** | shapes & conventions; subsume or retire | L |
| 1518–1539 | FullFieldWavefront | Verification-oracle wavefront sweep — dimension-generic full-field buffer, bit-identical to moving frontier (#240 D5b-S3 oracle) | **TWIN** | loss_representation.rst `loss-rep-four` (oracle architecture, C3.2b equivalence proof); cartesian_multid.rst full-field walk | H |
| 1556–1585 | _octant_face_cochain() | Allocate one octant's FULL per-axis interior face cochain; moment-valued (multi-moment closure, 2^{d-1}); three cases of boundary inflow | **CONTRACT** | cochain allocation & moment-validation contract (#251 guard shared with _inflow_to_moments); needed for oracle trace assembly | H |
| 1654–1663 | (comment: S6.4(d): oracle sweep = Jacobi schedule × per-octant walk frame) | Oracle uses Jacobi (no inter-group reflect); sweep_group + _interior_walk frame decomposition | **COMMENT-cut** | architectural narrative; remove or collapse | L |
| 1698–1707 | (comment: angular octant buffer carries trailing 2^d spatial-moment axis) | Per-octant accumulator shape for multi-moment closure; full field stores all moments | **COMMENT-cut** | buffer shapes; subsume or retire | L |
| 1738–1752 | loss_action() in FullFieldWavefront | Forward loss action ORACLE (L+C)ψ; full-field sweep + boundary residual assembly; interior kernel call | **CONTRACT** | oracle matvec implementation; calls loss_action_interior; boundary residual (O.4b active trace) | H |
| 1809–1853 | ScanMarch | Scan-march sweep — scan(x) marched over the transverse (y,z) tiles; Cartesian d≥2 production default (0.57–0.84× faster on measured numbers); 1-D scan at dimensionality boundary | **HISTORY** | Fork B2 decision (2026-06-11); loss_representation.rst `loss-rep-fork-b2` (measured evidence); .claude/plans/sn_sweep_strategy.md Fork B2 findings | H |
| 1847–1866 | (comment: 1-D arm reads is_affine_scannable; curvilinear gate) | 1-D scan selection gate via scheme trait; curvilinear capability gate #236 ST2 | **COMMENT-cut** | selection logic narrative; subsume into default_for() docstring or supports() predicate | L |
| 1912–1925 | (comment: multi-D ⇒ row-march sweep = schedule × row-march walk frame) | Multi-D uses row-march schedule + walk; 1-D uses chain scan | **COMMENT-cut** | architectural decomposition; teaching narrative; remove | L |
| 1947–1965 | _sweep_interior() in ScanMarch | Row-march interior kernel, SOLVE direction; one row pure-z + per-cell kernel per column | **CONTRACT** | kernel signature; row-march traversal order | L |
| 1972–1983 | (comment: per-mode row emission bound once before march) | Per-mode row accumulator allocation once per march (angular vs moment output modes); emit polymorphism | **COMMENT-cut** | implementation detail; pedagogical only; remove or inline | L |
| 2010–2023 | (comment: 239 coefficient-model lift; row-march asks scheme for coefficients) | Coefficient-model kernel signature; #239 steering coefficients; not directly from mesh.streaming | **COMMENT-cut** | implementation commentary; efficiency note; not load-bearing | L |
| 2048–2072 | loss_action() in ScanMarch | Forward loss action via row-march; per-row scan + boundary residual assembly | **CONTRACT** | ScanMarch matvec; row-march interior kernel call; boundary residual (O.4b) | H |
| 2088–2098 | _loss_action_interior() in ScanMarch | Row-march interior kernel, APPLY direction; per-row matvec kernel | **CONTRACT** | kernel signature; row-march traversal order | L |
| 2107–2118 | (comment: march y-rows in octant's y-sweep order; 239 coefficient model) | Y-row sweep order; coefficient-model kernel steering (#239) | **COMMENT-cut** | implementation detail; pedagogical; remove | L |

## orpheus/sn/loss_representation/sweep_graph.py


| 2194–2219 | default_for() | Select default sweep strategy for mesh; Cartesian → ScanMarch, curvilinear → CumprodScan, fallback to FullFieldWavefront; #236 ST2 curvilinear gate | **CONTRACT** | selection factory; needed at point of use (solver init); curvilinear capability gate | H |
| 2213–2226 | (comment: 236 ST2 rejection of (scheme × geometry) pairing) | Curvilinear-capability gate rejection; avoids silent wrong-physics sweep #236 ST2 | **COMMENT-cut** | implementation narrative; subsume into _curvilinear_capability() or default_for() docstring | L |
| 2224–2255 | (comment: 7-line fork decision block comment) | Default selection logic (Cartesian/1-D/curvilinear routing); Fork B2 performance motivation; measurements | **COMMENT-cut** | architectural narrative; teaching of Fork B2; remove or integrate into loss_representation.rst Fork-B2 section | L |
| 2272–2290 | _WalkLeg | One (μ-level × direction) chain of the 1-D walk DAG; ordinates, level index, within-level position; mu_level_idx + within-level mask | **CONTRACT** | leg data structure; needed for 1-D walk decomposition (M-M thread dispatch) | H |
| 2316–2346 | _OneDimScanWalk | Shared 1-D-scan frame — geometry-blind chain recurrence; both sweep & matvec; Slab joint-batch, curvilinear per-ordinate M-M thread | **HISTORY** | #196 PR-INDEX refactor; loss_representation.rst (1-D scan family + unified layout + Fork-B2 1-D arm); .claude/plans/sn_sweep_strategy.md | H |
| 2350–2369 | _dag_legs() | Every non-empty leg of the 1-D walk DAG, in DEPENDENCY order (all +1 first, then −1); M-M thread dispatch per leg | **CONTRACT** | leg enumeration for _loop_walk; dependency ordering (primal sweep seed) | H |
| 2407–2426 | _loop_walk() | THE shared 1-D per-cell loop frame (#280 2.5a); generic open/visit/close cycle for both fwd & bwd walks; reverse-traversal reuses this frame | **TWIN** | loss_representation.rst `loss-rep-orientation-two-frames` (#280 2.5a reverse-walk unified frame); loss_representation.rst one-walk theorem | H |
| 2433–2445 | _degenerate_positions() | Degenerate pure-azimuthal ordinates + their level/within indices (cylinder only; empty slab/sphere) | **CONTRACT** | degenerate-ordinate set for volumetric-balance branch (no legs, no faces) | H |
| 2450–2471 | (comment: starting-direction ψ½ block retired from carrying-mesh arm) | #282 route (a) step 6 — the ψ½ starting-direction fold retired from sweep_direction; M grid block solve carries it | **HISTORY** | #282 carrying-mesh architecture; ORPHEUS structure change (M-grid responsibility shift); internal comment; remove or move to theory | L |
| 2478–2525 | sweep() | Geometry-blind 1-D SN sweep; three numpy tensor layout conventions; Slab 2-scan fast path, curvilinear per-ordinate M-M; interior kernels supplied by calling representation | **CONTRACT** | public sweep interface (#196 PR-INDEX-* layout conventions); calls _run() | H |
| 2535–2545 | sweep_transpose() | Reverse-solve (L+C)⁻ᵀ — 1-D forward + back half-sweeps; #280 2.5b reverse-scan signature | **CONTRACT** | call contract; reverse-scan frame; deferral for multi-D (O.2b) | H |
| 2552–2566 | loss_action() | 1-D forward loss action (L+C)ψ — matvec, dual of sweep; calls _apply_walk; returns angular + boundary residual | **CONTRACT** | matvec call interface; calls _apply_walk(); boundary residual output | H |
| 2583–2630 | _apply_walk() | 1-D apply-direction walk (fused (L+C)ψ + σ-dependent + Morel–Montry + Carlson seed + boundary ACTIVE trace); 44-line docstring | **SPLIT** | (see next section) | — |

#### 2583–2630: _OneDimScanWalk._apply_walk() [LARGE PROSE ANALYSIS]

| section | content-id | verdict | destination / anchor | conf |
|---|---|---|---|---|
| 2584–2600 | Fused (L+C)ψ + σ-dependent closure + angular redistribution + seed propagation; coefficient model / cell_balance; O.4a.2 keystone deleted (no BC apply inside) | **TWIN** | loss_representation.rst `loss-rep-removal-form-matvec` (σ-affine decomposition, BC extraction); curvilinear_numerics.rst Morel-Montry/Carlson sections | M |
| 2601–2630 | Cartesian: coefficient-model kernel, wave-O σ-affinity; Curvilinear: DD-only cell_balance path with angular redistribution; pole-face seed; Degenerate ordinate branch; Boundary active trace (O.4b) | **CONTRACT** | matvec implementation frame; two geometry arms (Cartesian / curvilinear); degenerate-ordinate volumetric balance; boundary active-trace residual construction (O.4b) | H |

| 2656–2666 | (comment: unified moment matvec — trailing axis on probe/accumulators) | Spatial-moment axis broadcast; DD scalar, LD moment-tailed | **COMMENT-cut** | shapes & conventions; subsume or retire | L |
| 2676–2695 | (comment: 282 route (a) step 6 — carrying-mesh block action) | Ray-decoupled (A,A) block action; zero seed; joint M action is grid's block matvec | **COMMENT-cut** | architectural narrative of carrying-mesh structure; remove or move to theory | L |
| 2702–2717 | (comment: d=1 probe GLOBAL-frame, residual SWEEP-frame; moment sign-flip for backward sweep) | Moment frame-sign convention (d=1 inverse of d≥2 _CellResidual); sweep-direction dependence; Pattern 2 single-source binding | **COMMENT-cut** | moment-frame rules teaching; integrate into #240 D5b-S3 section or _moment_frame_signs() docstring | L |
| 2718–2761 | (comment: Coefficient model, Cartesian matvec, diamond march, unified moment matvec, #240 D5b-S3) | 22-line implementation details: kernel model, DD inlining, moment broadcast, source-free apply, no angular redistribution Cartesian | **COMMENT-cut** | implementation narrative; complexity teaching; remove or move to design rationale section | M |
| 2776–2791 | (comment: Curvilinear matvec, cell_balance density path, Morel-Montry, no LD closure, DD inlined) | Implementation narrative of curvilinear arm; cell_balance single-source; no-LD justification; DD diamond march | **COMMENT-cut** | implementation narration; remove or condense | M |
| 2819–2832 | (comment: Wave O O.4a.2 keystone deleted; backward sweep seeded from outer inflow, not reflected outflow) | Wave O architectural change (O.4a.2) — BC decoupling; outflow persistence moved to sibling −B; boundary-consistency loop driver | **COMMENT-cut** | architectural narrative of O.4a.2 decoupling; teaching value; remove or move to loss_representation.rst O.4a.2 section | L |
| 2826–2845 | (comment: Carlson coupled-pole seed; continuation at r=0; #192 deferred pole condition; Wave O O.4a.2) | 10-line pole-face seed logic; curvilinear pole semantics; inward-determines-outward; historically wrong-on-curved | **COMMENT-cut** | pole-seed contract narrative; domain knowledge; remove or move to curvilinear_numerics.rst | L |
| 2873–2906 | (comment: Wave O O.4a.2 boundary block residual; DIAGONALS only; outflow defect + inflow identity) | Boundary residual construction (O.4b active trace); two trace diagonals (OUTFLOW / INFLOW); sibling −B carries off-diagonal; consistency loop | **COMMENT-cut** | boundary architecture narrative; remove or condense into main method docstring | M |
| 2936–2973 | loss_action_transpose() | 1-D adjoint loss action (L+C)ᵀφ — matvec transpose; reverse-mode adjoint; #206 Phase C relocation; reverse-substitution sweep | **HISTORY** | #206 Phase C relocation (verbatim off _MSpatialOperatorSum._compute_LpC_transpose); #280 2.5a structural reverse walk | H |

| 2985–2998 | (comment: Phase 2.5 S0 — reverse walk hand-transposes DD face-flux chain; LD deferral) | Reverse walk hand-transpose scope; LD moment-tailed cotangent guard (typed deferral); StreamingOperator.is_adjointable eager refuse | **COMMENT-cut** | deferral scope narrative; subsume into NotImplementedError message or docstring | L |
| 3013–3024 | (comment: 282 route (a) B.2d step 6 — ray-decoupled (A,A)ᵀ block; seed pullback to A_ABᵀ grid) | Carrying-mesh architecture; ray-decoupled adjoint block; seed-pullback responsibility routing | **COMMENT-cut** | architectural routing narrative; remove or move to theory | L |
| 3057–3076 | (comment: reverse spatial DD marches; reverse traversal via reverse_traversal(_dag_legs()); 280 2.5a leg linearization) | 10-line reverse-walk structural narrative; leg traversal reversal; mirror-ordinate handoff; pre-2.5a vs 2.5a refactoring | **COMMENT-cut** | reverse-walk algorithmic narrative; remove or condense | L |
| 3108–3121 | (comment: adjoint Carlson coupled-pole seed; mirror permutation; level-local; pre-2.5a nesting) | Pole-seed adjoint routing; mirror-ordinate permutation; ordering independence (disjoint slots) | **COMMENT-cut** | adjoint pole-seed technical narrative; remove | L |
| 3122–3149 | (comment: reverse degenerate-ordinate branch; pre-2.5a silent DROP — cyl_product_2g G-reciprocity red→green) | 14-line degenerate-ordinate adjoint narrative; pre-2.5a equispaced-product silent bug; reciprocity gate fix | **COMMENT-cut** | bug-fix history & degenerate ordinate handling teaching; remove or move to design rationale | L |
| 3173–3188 | (comment: reverse angular factor delegated; zero for slab; M-M thread cotangent; seed-cells discard) | Angular-adjoint deferral; carrying-mesh seed pullback responsibility (#282 route (a)) | **COMMENT-cut** | architectural responsibility routing; remove | L |
| 3218–3229 | _ensure_coll_cache() | Collision cache builder (lazy if absent); cache invariant (#4); ad-hoc test bypass; PR-INDEX-3 bridge (ng, nx, ny=1) → (ng, nx) | **CONTRACT** | cache lazy-build protocol; bridge for multi-D inputs | H |
| 3246–3280 | _run() | Inner body unified 1-D sweep (34-line docstring); #196 PR-INDEX-1–5 layout conventions; principled (N, ng, *spatial) throughout; CollisionCache native (N, ng, nx) | **HISTORY** | #196 PR-INDEX refactor (complete 2026-04-30); layout convention enforcement; internal implementation | H |

| 3281–3294 | (comment: 282 route (a) step 6 starting-direction contract) | #282 carrying-mesh; ray-decoupled (L+C) block; joint M grid's coupled substitution | **COMMENT-cut** | architectural routing; remove | L |
| 3300–3315 | (comment: Spatial-moment width — unified moment matvec #240 D5b-S3 OWED-2) | Moment-tail SSOT; LD trailing 2^d, DD scalar (byte-id negative control); single-source via face_moment_tail | **COMMENT-cut** | moment convention teaching; subsume into method docstring or retire | L |
| 3312–3327 | (comment: Multi-moment closure lifts FLAT scalar source; exact as DAG's _ubld_system #240 D5b-S3) | Source-moment lifting contract (scattering vs external/MMS); RANK discriminator; same convention DAG uses | **COMMENT-cut** | moment source normalization teaching; move to #240 section or remove | L |
| 3325–3340 | (comment: Common pre-scale; R-1 Step 4 A1; single per-ordinate source; 1/W producer, sweep multiplies V only) | Source pre-scaling (×V); no iso/aniso distinction internally; moment-scale convention | **COMMENT-cut** | internal algorithm narrative; remove | L |
| 3333–3364 | (comment: BC inflow + per-level Carlson seed (curvilinear only) — 16-line comment block) | Wave O O.4a.2 BARE sweep (no bc.apply); reflective coupling via sibling −B; source fold (S + B); per-ordinate seed; M thread starting-direction zero | **COMMENT-cut** | Wave O architectural narrative; multi-line teaching; remove or move to loss_representation.rst Wave-O section | L |
| 3361–3370 | (comment: SLAB joint-batch fast path — D-H.2-C2 writable per-face views) | Slab layout / xmin/xmax; per-cell outflow persistence; aliasing safety (disjoint ordinate sets) | **COMMENT-cut** | implementation optimization narrative; remove | L |
| 3436–3469 | (comment: Multi-moment (LD) slope-source scan — 17-line block comment) | LD slope-source lift; single-source via scan; contrast with DAG (both use same convention) | **COMMENT-cut** | moment-source normalization narrative; remove | L |
| 3536–3545 | (comment: CURVILINEAR per-ordinate path — 5 lines) | Curvilinear requires per-ordinate scans (M-M angular thread couples ordinates sequentially); joint-batch future optimization | **COMMENT-cut** | architectural constraint narrative; remove | L |
| 3536–3591 | (comment: whole curvilinear scan IS Morel-Montry thread — 28 lines) | Curvilinear sweep carries M-M recurrence IN-SWEEP (per-ordinate loop over levels); Hébert 3.9.4 refs; joint-batch future research | **COMMENT-cut** | M-M algorithm narrative; research roadmap; remove or move to curvilinear_numerics.rst | L |
| 3583–3618 | (comment: 280 2.5b direct-seed fold — 18 lines) | Seed-ordinate treatment (ψ½ ≡ ψ̄); self-coupling; carrying-mesh step 6; M grid block solve carries joint march | **COMMENT-cut** | carrying-mesh architectural narrative; remove | L |
| 3648–3659 | (comment: Carlson coupled-pole seed #195 — 6 lines) | Pole-seed behavior (slab/sphere/cylinder); reflection index; carrying-mesh | **COMMENT-cut** | pole-seed technical narrative; remove or shorten | L |
| 3653–3678 | (comment: ZERO thread step 6 — carrying-mesh — 13 lines) | #282 route (a) step 6 scope; walk IS (A,A) block; zero seed feed; walk reads zeros; joint march M routing | **COMMENT-cut** | carrying-mesh structural narrative; remove | L |
| 3691–3708 | (comment: Coupled-pole seed (ERR-058a) mirror inward at r=0 — 9 lines) | Mirror-ordinate permutation; curvilinear pole condition; continuation at center | **COMMENT-cut** | pole-seed technical details; remove or shorten | L |
| 3738–3749 | (comment: 280 2.5b direct seed ψ½ ≡ ψ̄; self-coupling; step 6) | Direct-seed treatment redundancy check; self-coupling formalism; M-grid responsibility | **COMMENT-cut** | carrying-mesh narrative; remove | L |
| 3805–3818 | (comment: Exit — PR-INDEX-5; caller consumes principled layout) | Layout output contract; no exit transposes (caller-side principled-layout outputs); phantom ny axis handling | **COMMENT-cut** | layout contract narrative; can integrate into _run() docstring tail | L |
| 3828–3857 | _run_transpose() | Reverse-scan frame; 29-line docstring; #280 2.5b kernel-pair contract; R12a starting-direction semantics | **HISTORY** | #280 Phase 2.5b reverse-scan refactor; kernel-pair architecture | M |


| 3861–3872 | (comment: R12a starting-direction contract — 6 lines) | Carrying-mesh; ray-decoupled reverse-scan block; zero-seed read; joint M transpose grid | **COMMENT-cut** | architectural routing; remove | L |
| 3945–3960 | (comment: CURVILINEAR reverse-scan — 8 lines) | Sphere carrying-mesh adjoint; cylinder carrying-mesh adjoint; pole-center semantics | **COMMENT-cut** | architectural narrative for carrying-mesh adjoint; remove | L |
| 3971–3982 | (comment: level structure + seed dispatch — 6 lines) | M-M thread level-dispatch; seed-cells handling; carrying-mesh routing | **COMMENT-cut** | architectural narrative; remove | L |
| 4022–4041 | (comment: degenerate pure-azimuthal ordinate cylinder — 10 lines) | Degenerate-ordinate adjoint (slot-local, no face terms); denom × ob / V logic | **COMMENT-cut** | degenerate adjoint narrative; remove | L |
| 4043–4052 | (comment: ψ̄ = inv_denom·(QV + |μ|·A_total·ψ_in + κ·ψ_ang) — 5 lines) | Cell balance residual formula; naming of terms; κ angular-numer term | **COMMENT-cut** | algorithm formula teaching; remove | L |
| 4053–4062 | (comment: boundary degenerate ord does NOT overwrite — 5 lines) | Degenerate ordinate safety (no boundary write); discard semantics | **COMMENT-cut** | implementation safety narrative; remove | L |
| 4072–4083 | (comment: seed ordinate (cylinder) direct-seed FOLD — 6 lines) | Seed-ordinate self-coupling; direct-seed fold treatment; step 6 | **COMMENT-cut** | carrying-mesh seed narrative; remove | L |
| 4126–4135 | (comment: seed ord reads its OWN average — 5 lines) | Seed-ordinate self-coupling; M-M thread seed source | **COMMENT-cut** | M-M thread technical detail; remove | L |
| 4149–4158 | (comment: M-M thread cotangent DISCARD on carrying level — 5 lines) | Seed-cell cotangent discard; A_ABᵀ grid block responsibility | **COMMENT-cut** | carrying-mesh adjoint responsibility routing; remove | L |
| 4173–4242 | _sweep_scheduled() | Polymorphic schedule-driven 2-D sweep (69-line docstring); Phase 3 super-summary **[LARGE PROSE — see next section]** | **SPLIT** | — | — |

#### 4173–4242: _sweep_scheduled() [LARGE PROSE ANALYSIS]

| section | content-id | verdict | destination / anchor | conf |
|---|---|---|---|---|
| 4174–4191 | Polymorphic schedule driver; Jacobi/Gauss-Seidel inter-group reflect; schedule loop; moment output mode; all three swept separately | **CONTRACT** | schedule loop control; inter-group reflection routing; output mode polymorphism (angular/moment) | H |
| 4192–4220 | Gauss-Seidel scheduling, moment accumulation across groups, emit polymorphism, pure-z dispatch, interior kernel call pattern, per-dim streaming coefficients | **CONTRACT** | scheduling protocol; kernel interface; operand bindings; streaming tuple positional-by-axis | H |

| 4245–4254 | (comment: output buffers carry trailing 2^d spatial-moment axis) | Moment-valued output; LD trailing axis, DD scalar | **COMMENT-cut** | shapes; subsume or retire | L |
| 4261–4270 | (comment: moment buffer carries frame table) | Harmonic-moment buffer shape; L+1 × 2L+1 × ng × spatial | **COMMENT-cut** | buffer shape convention; retire | L |
| 4287–4296 | (comment: G-S inter-group reflect — no-op Jacobi schedule) | Gauss-Seidel reflect with carry-in; Jacobi route (skip). | **COMMENT-cut** | scheduling detail; remove | L |
| 4296–4305 | (comment: moment mode output — moments, None; scalar IS φ_0^0) | Mode-specific output tuple (moments vs angular); scalar subsumed | **COMMENT-cut** | output-mode detail; remove | L |
| 4317–4348 | _sweep_scheduled_body() | Bare multi-D sweep = Jacobi octant schedule + operation loop + typings; 31-line docstring (details internal to _sweep_scheduled) | **CONTRACT** | encapsulation of Jacobi schedule + op loop (internal to _sweep_scheduled); sweep body template | L |


---

| 0–72 | MODULE | Per-octant causal cell DAG for Cartesian wavefront sweep (d-generic); §15A.2 upwind trace complex / causal transport DAG primitive; mesh-time derived object (shape + octant cached); SN-specific by design; tensorial framing | **SPLIT** | — | — |
| 0–27 | (module intro) | DAG description, no flux/source/state dependence, shape-octant cache (S6.4c), reused across SI/Krylov/outer loops | **TWIN** | loss_representation.rst `loss-rep-four` (family-owned DAG cache, bit-identical twin oracle); wavefront_cochain.rst DAG framing | H |
| 18–26 | (Cardinal Rule 2 framing) | SN-specific by design; MoC ≠ shared Protocol; fiber bundles + solution sheaves | **HISTORY** | project_moc_structure.md (SN-specific abstraction; MoC future); architectural principle documentation | H |
| 28–43 | (Tensorial framing) | Anti-diagonal schedule vectorises (N_oct, nx, ny, ng) tensor; no Python loops over ordinates/cells/groups; kernel calls batch all three; ordinate axis internal to kernel | **TWIN** | loss_representation.rst `loss-rep-four` (anti-hyperplane level operation, vectorised dispatch); #239 coefficient-model steering | H |
| 45–61 | (References) | Grand Report v3 §15A.2 "upwind trace complex"; assert_* invariant set; Wave 2 plan C2.3 dispatch boundary | **CONTRACT** | literature citations; design boundary documentation | H |

| 89–120 | _reframe() | Sweep⇄global moment-frame involution (#240 D5b-S3); reframe iff moment buffer at multi-moment closure; frame_signs = octant_moment_frame_signs (2^d), None for DD/Step | **CONTRACT** | moment-frame sign-flip protocol; #246 fix (is_moment_valued typed intent, not coincidental size); map is self-inverse | H |
| 129–151 | OctantLabel | Octant signature (signs per axis); d-generic label ((±1), (±1,±1), (±1,±1,±1)); out-of-plane projection (schedule SOLE site); degenerate (all-zero) semantics; frozen/slotted hashable | **TWIN** | loss_representation.rst `loss-rep-four` (octant encoding); sweep_schedule.rst label projection rule | H |
| 182–220 | _LevelFrontier | Mesh-time slab addressing for ONE anti-hyperplane level; cell lists per sweep direction; all-cells-in-level view; orientation (upwind trace complex) | **CONTRACT** | level data structure; orientation conventions; used by rolling frontier plan | H |
| 229–260 | _FrontierPlan | Whole-sweep rolling (d−1)-frontier windowing; allocation plan per level; cell lists + face cochain access via _LevelFrontier | **CONTRACT** | frontier plan data structure; memory layout; level enumeration order | H |
| 273–288 | _MovingFrontier | Rolling (d−1)-frontier interior face cochain (7-line docstring); per-level allocation; index mapping | **COMMENT-cut** | data structure detail; integrate into _FrontierPlan or retire | L |
| 290–303 | (comment: interior face cochain carries 2^{d-1} trailing moment axis) | Face-cochain moment width for multi-moment closure; transverse-moment count | **COMMENT-cut** | shapes; subsume or retire | L |
| 375–427 | SweepDependencyGraph | Per-octant causal cell DAG for 2-D Cartesian wavefront sweep (52-line docstring) **[LARGE PROSE — see next section]** | **SPLIT** | — | — |

#### 375–427: SweepDependencyGraph [LARGE PROSE ANALYSIS]

| section | content-id | verdict | destination / anchor | conf |
|---|---|---|---|---|
| 376–398 | DAG stores the per-octant cell sweep order (topologically sorted, upwind orientation); two level-walk algorithms (full-field oracle, rolling frontier); bit-identical by C3.2b oracle | **TWIN** | loss_representation.rst `loss-rep-four` (DAG structure, two buffer policies, oracle equivalence); wavefront_cochain.rst level walks | H |
| 399–416 | Attributes: cells_per_level (nested list), level_walks (full-field / frontier), family accessor sweep_graphs (cached per shape), construction validation (for_shape, from_cartesian) | **CONTRACT** | DAG data members & family interface; shape-keyed cache (S6.4c); construction methods | H |
| 420–427 | (comment: Phase 5b storage-B mesh-time rolling frontier pre-fetch) | Phase 5b frontier planning at mesh construction time (S6.4c); deferred optimization | **HISTORY** | Phase 5b storage-B frontier prefetch decision; mesh-time precompute (vs sweep-time build) | L |

| 487–514 | for_shape() | Per-octant DAG family for Cartesian grid shape; cached, shape-keyed lookup (S6.4c) | **CONTRACT** | DAG factory / cache accessor; (shape, octant) mapping | H |
| 523–558 | from_cartesian() | Build per-octant upwind DAG for d-dim Cartesian mesh; cell topological sort (upwind trace complex); face-pairing consistency; 1-D/2-D/3-D code path | **HISTORY** | algorithm description; §15A.2 upwind-DAG construction; #236 bug-fix context | M |
| 622–640 | _build_frontier_plan() | Build rolling (d−1)-frontier window for the per-octant DAG (18-line docstring); _LevelFrontier + _FrontierPlan assembly | **CONTRACT** | frontier plan construction; level ordering; cell/face access templates | H |
| 637–646 | (comment: det is SINGLE SOURCE OF TRUTH for frontier free dimension) | Determinant sign encodes swept-away dimension; single-source convention for free-dim handling | **COMMENT-cut** | implementation convention detail; remove or shorten | L |

| 701–720 | (comment: level walks — TWO storage policies × ONE dimension — 10 lines) | Full-field vs rolling-frontier walks; same octant order, different buffers; polymorphism over walk_full / walk_windowed | **COMMENT-cut** | architectural narrative of walk polymorphism; remove | L |
| 730–744 | walk_full() | Full-cochain level walk — verification-oracle kernel wrapper (14-line docstring); FULL interior face cochain per level; calls full-field level kernel | **TWIN** | loss_representation.rst `loss-rep-four` (oracle walk family); FullFieldWavefront._octant_face_cochain interface | H |
| 777–794 | walk_windowed() | Rolling-frontier level walk — production storage-B (17-line docstring); per-row buffering strategy; calls rolling-frontier interior kernel | **TWIN** | loss_representation.rst `loss-rep-four` (rolling frontier walk family); MovingFrontierWindow._sweep_interior interface | H |
| 820–831 | _cell_face_selector() | Advanced-index tuple selecting axis ``axis``'s faces (inflow/outflow) for a cell set (11-line docstring); face indexing for level operation | **CONTRACT** | face-selection primitive; used by walk kernels to slice cochain | H |
| 825–848 | (comment: face selection logic — 12-line technical block comment) | Advanced indexing for per-axis inflow/outflow faces; octant orientation encoding in cell address | **COMMENT-cut** | implementation technical narrative; remove or integrate into _cell_face_selector() docstring | L |

| 852–871 | _CellSolve | SOLVE-direction level operation (19-line docstring); WDD balance on cell; moment-valued iterate + source; three source cases (moment-lifted, moment-resolved, scalar) | **CONTRACT** | cell solve contract; source discrimination (RANK-based); moment handling protocol | H |
| 880–905 | (comment: cell kernel consumes/produces moment vector — 13 lines) | Moment-valued cell buffers; spatial-moment axis; probe + solve dynamics; per-axis entry/exit moment conventions | **COMMENT-cut** | moment protocol teaching; integrate into main docstring or retire | L |
| 893–920 | (comment: cell SOLVE source cases — 14 lines) | Three source discriminations: moment-lifted, moment-resolved, scalar; RANK gates; averaging gate; flat external source | **COMMENT-cut** | source-discrimination narrative; teaching; integrate or remove | L |
| 922–943 | (comment: multi-moment closures (LD) bilinear UBLD in EVERY geometry — 11 lines) | LD closure invariant; cell-moment count 2^d; per-axis Legendre order; moment axis always present for LD | **COMMENT-cut** | moment-closure teaching; remove or shorten | L |
| 993–1020 | (comment: unified moment matvec — 14 lines) | Moment-valued probe + accumulator; σ broadcast; moment axis conventions; DD scalar byte-identity | **COMMENT-cut** | moment matvec narrative; subsume or retire | L |


---

## Summary

### Verdict counts by file

#### orpheus/sn/loss_representation/__init__.py

- **CONTRACT**: 45 blocks ~1,450 lines → **KEEP**
- **TWIN**: 8 blocks ~380 lines → **ARCHIVE** (documented in theory pages; inline refs point to anchors)
- **HISTORY**: 4 blocks ~280 lines → **ARCHIVE** (loss_rep-history + #196/#280/#282/#206 sections + .claude/plans; remove campaign narration)
- **COMMENT-cut**: 34 blocks ~380 lines → **DELETE or INLINE** (implementation detail + teaching narrative; subsume into method docstrings or condense)
- **SPLIT**: 3 blocks (1a–1d, 1025–1049 analysis, 2583–2630 analysis, 4173–4242 analysis) → **HANDLE PER SUB-SECTION** (each tier above)

**Docstring lines (≥10 lines): ~2,110 lines total**
- Keep (CONTRACT + TWIN): ~1,830 lines
- Archive (HISTORY + teaching TWIN): ~310 lines
- Remove (COMMENT-cut): ~220 lines

#### orpheus/sn/loss_representation/sweep_graph.py

- **CONTRACT**: 14 blocks ~450 lines → **KEEP**
- **TWIN**: 6 blocks ~280 lines → **ARCHIVE** (documented in theory; inline refs to anchors)
- **HISTORY**: 2 blocks ~80 lines → **ARCHIVE** (Phase 5b/§15A.2 design history; compress to inline refs)
- **COMMENT-cut**: 7 blocks ~90 lines → **DELETE or INLINE**
- **SPLIT**: 2 blocks (module intro + SweepDependencyGraph) → **HANDLE PER SUB-SECTION**

**Docstring lines: ~360 lines total**
- Keep (CONTRACT + TWIN): ~300 lines
- Archive (HISTORY + teaching TWIN): ~90 lines
- Remove (COMMENT-cut): ~60 lines

### Verdict-outcome charter

| Verdict | Action | Rationale | Example |
|---|---|---|---|
| **CONTRACT** | KEEP in code; no xref needed | Call-time essential; interfaces, error messages, bounds | `streaming_action()` σ-free decomposition; `Resolution A` guarantee |
| **TWIN** | ARCHIVE: xref to theory anchor; condense or delete docstring duplicate | Concept documented in theory; code need only cite it | `_OctantWalk._interior_walk()` one-walk theorem → `loss-rep-one-walk-one-instance` |
| **HISTORY** | ARCHIVE: xref to theory History section or GitHub issue; delete narration | Campaign-level documentation belongs in theory "Development history" + issue tracker | `sweep_strategy.md S0–S6.9 arc` → `loss_representation.rst loss-rep-history` |
| **COMMENT-cut** | DELETE or INLINE to main docstring; never standalone | Pedagogical / implementation detail; teaching belongs in theory page or in-code docstring, not comment | Loop-variable commentary, coefficient-model narration |
| **SPLIT** (→ tiers) | Handle per sub-section; apply verdict to each tier | Docstring contains 2+ distinct purposes; classify each prose block independently | 1–125 module docstring splits 1a–1d (CONTRACT / TWIN / CONTRACT / HISTORY) |

### Operator-algebra posing flags

**Current record:** `A = L + C − S − B` (orpheus/sn/operators/* and loss_representation.rst 2709–2815, Development history)

**Posings found in code:**
- Line 3–4: `ψ = (L+C)^{-1} q` + `(L+C)ψ` — **CORRECT**, uses full (L+C) notation (not L−S−F)
- Line 24: `(L+C)` — **CORRECT**
- Line 269–280: Protocol docstrings all use `(L+C)` — **CORRECT**
- Line 1026: `(L+C)ψ̄` matvec — **CORRECT**
- Line 2478–2525 _run() sweep frame docstring — **CORRECT**, uses unified (L+C)
- **No outdated (L−S−F) or (A−S−F) posings found.** ✓

### Uncertain classifications (confidence L — would benefit from domain review)

| line | symbol | verdict | reason |
|---|---|---|---|
| 1374–1385 | MovingFrontierWindow._sweep_interior() | **CONTRACT** [L] | Interior kernel docstring § is short (11 lines); unclear if METHOD-level doc is needed or subsumable into windowed-walk architecture spec |
| 1473–1483 | MovingFrontierWindow._loss_action_interior() | **CONTRACT** [L] | Same: interior kernel method, kernel signature load-bearing, but terse docstring |
| 1947–1965 | ScanMarch._sweep_interior() | **CONTRACT** [L] | Same as above |
| 2088–2098 | ScanMarch._loss_action_interior() | **CONTRACT** [L] | Same as above |
| 3218–3229 | _ensure_coll_cache() | **CONTRACT** [L] | Cache-builder protocol; not directly call-time (internal to sweep), but documents lazy-build invariant + PR-INDEX-3 bridge |
| 456–475 (comment) | pure-L streaming primitive | **COMMENT-cut** [L] | Fragment narratives the coding-elegance Pattern 2 (one walk, no twin); could be TWIN if loss_representation.rst has Pattern-2 section |
| 364–382 streaming_action() | **TWIN** [M] | Affine-decomposition σ-free leaf; "TWIN" only if loss_representation.rst removal-form section explicitly covers affinity + σ=0 equivalence |

---

## Disposition: what stays, what moves

### Immediate actions (no code changes required)

1. **Read & verify anchors**: confirm each TWIN block has a corresponding anchor in loss_representation.rst / wavefront_cochain.rst / etc.
   - If anchor missing → downgrade to **MOVED** (add section to theory page)
   - Example: `_interior_walk()` docstring points to `loss-rep-one-walk-one-instance` anchor; verify it exists and is complete

2. **Locate HISTORY blocks in theory or issues**
   - Example: `sweep_strategy.md` must document S0–S6.9 carve arc; verify completeness
   - If theory History section is thin → **MOVED** (expand loss_representation.rst Development history)

3. **Identify COMMENT-cut removal candidates** that can consolidate to method docstrings
   - E.g., lines 2718–2761 (22-line coefficient-model narration) + lines 2776–2791 (8-line curvilinear narration) + lines 2819–2832 (7-line Wave-O deletion) could be 1 consolidated section in `_apply_walk()` docstring header

### Follow-up issues to file (if not already present)

1. **Interior kernel docstring scope** (Contract L blocks): are `_sweep_interior()` / `_loss_action_interior()` method-level docstrings needed or supersedable by class-level walk description?

2. **Moment-frame teaching consolidation** (multiple COMMENT-cut blocks on #240 D5b-S3): migrate all moment-axis convention teaching to ONE section in loss_representation.rst §#240 D5b-S3, then replace comments with xrefs.

3. **Carrying-mesh architecture narrative** (#282 route (a) step 6): numerous comments (lines 2676–2695, 3013–3024, 3653–3678, etc.) describing M-grid responsibility routing → consolidate to ONE diagram/section in loss_representation.rst or new radial_characteristic.rst section.

4. **Operator algebra posing audit**: ALL docstrings currently use (L+C) correctly; **no refactoring needed**. Document this as a passing audit in the session record.

