# Phase 2 Code-Prose Classification: solver.py & iteration.py

Inventory of every docstring ≥10 lines and comment block ≥5 lines (module, class, method, property level).

**Special rule applied:** Module head at `iteration.py` lines 1–138 posing `(A − Σᵢ gᵢ)ψ = q_ext` is classified **CONTRACT** per notation.rst crosswalk promise (dual-A bridge statement). Posing-drift checks flagged separately (lines ~1/7/32/861 for `(A − S − F)ψ = q_ext`).

---

## FILE: orpheus/sn/solver.py

### Module Level

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|-----------|---------|----------------------|-----|
| 1–31 | MODULE | Wave E Round 2: operator triple construction; inner solver dispatch (SI vs Krylov); boundary condition defaults | CONTRACT | Stays: entry-point overview, references theory page | H |

### Module-level functions and comments

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|-----------|---------|----------------------|-----|
| 97–108 | `_apply_default_bcs` docstring | Normalization of geometry declarations; legacy Mesh1D/Mesh2D vs axis-tuple handling; all-or-nothing BC semantics | CONTRACT | Stays: interface contract (C5.5, #225 reference) | H |
| 135–150 | `_as_sn_mesh` docstring | Entry-surface seam for legacy geometry vs axes; boundary_condition fixed-source vacuum convention; mat_map channel | CONTRACT | Stays: interface and parameter semantics | H |
| 166–176 | Comment: Issue #197 PR-TYPED-5 retirement | SNFixedSourceResult / SNResult retired; Solution unified type; legacy data bag collapse (Pattern 7) | HISTORY | Record in PR-TYPED series memory (decision record: bit-identical via Solution.is_eigenvalue/is_fixed_source discriminator) | M |
| 178–189 | Comment: _within_group_triple / _lagged_gains retirement | Wave O campaign; WithinGroupSystem builder consolidation; B separation ruling; FullField-gain adapters retired | HISTORY | Record in PR-typed/coupled-system memory | M |
| 247–291 | `evaluate_residual` docstring | **SPLIT:** (a) 247–267 = CONTRACT (residual math, system A/B split, parameter specs); (b) 268–291 = HISTORY (box-7 consumer, #208, Mode-12, DSA #2 hook) | SPLIT | (a) stays; (b) → theory/solver.rst "Residual diagnostics" section + DSA plan document | M |
| 389–393 | Comment: _CERTIFICATE_SAFETY guard | Free identity vs honest residual agreement; FP-reassociation grain; lag-death defect threshold | HISTORY/CONTRACT | Partial: move threshold rationale to solver.rst "Convergence certificate" as design detail; keep constant inline as engineering value | L |
| 504–512 | Comment: _GaussSeidelResolvent retirement | Phase 3 sub-step 3c; scheduled invertible operator; G-S splitting reified; Jacobi path congruence | HISTORY | Record in stencil-assembly/DSA campaign memory (Phase 2.5 work) | M |
| 515–522 | Comment: _MomentWindowedResolvent retirement | #226 taxonomy step 3 rewire; windowed product direct consumption; accepted-and-dropped initial_guess | HISTORY | Record in frame-projection machinery campaign (Phase 5a) | M |
| 529–551 | `_maybe_window` docstring | **SPLIT:** (a) 529–545 = CONTRACT (gate condition, 2-D Cartesian moment windowing, frame sourcing); (b) 546–551 = HISTORY (C5.4 gate specification, pre-C5.4 proxy hazard, Mode 9 vv) | SPLIT | (a) stays; (b) → theory/solver.rst "Spatial windowing" section + vv-principles Mode 9 record | M |
| 562–583 | `_windowed_cold_start` / `_unwindowed_cold_start` docstrings (combined ~23 lines) | Initial iterate representations for windowed vs unwindowed SI; moment vs full-angular; LD axis handling; Pattern 2 twin (cold-start sibling) | CONTRACT | Stays: interface specs for cold-start factories | H |
| 611–639 | `_radial_characteristic_source_from_per_ordinate` docstring | Fold kernel entry point; per-ordinate to q½ composite; eigenvalue vs fixed-source routes; boundary-trace arm (step 7 corner law) | CONTRACT | Stays: interface and factory semantics | H |
| 642–673 | `_radial_characteristic_fission_seed` docstring | Direct ℓ=0 moments-fold of fission emission; factory round-trip replacement (campaign step 4c commit 2); REPLACES narrative; fold kernel shared | HISTORY/CONTRACT | Hybrid: CONTRACT for direct fold math (stays); HISTORY for "replaces round-trip" → coupled-system campaign memory | M |
| 676–708 | `_coupled_flux_state` / `_coupled_source_state` docstrings (combined ~33 lines) | System B zero initialization (B.2d); q½ fold pair; carrying mesh native birth; presence-gated field | CONTRACT | Stays: interface for coupled composites | H |
| 711–779 | `_select_si_resolvent` docstring | Jacobi / Gauss-Seidel selector (Phase 3 R2); scheduled invertible operator; boundary-G-S splitting; C5.4 genuine gate condition | HISTORY/DESIGN | Split: gate condition logic (lines 745–750) stays inline; Phase 3 campaign narrative (C5.4 proxy) → stencil-assembly memory | M |
| 782–862 | `_within_group_si` docstring | **SPLIT:** (a) 782–820 = CONTRACT (coupled/seedless dispatch, driver outputs); (b) 821–862 = HISTORY (Wave O narrative, B.2d DP-seedless, coupled vs seedless paths, windowing/schedule machineries) | SPLIT | (a) stays; (b) → coupled-system + stencil-assembly memory | M |
| 1038–1066 | Comment: Two-stratum cache (Issue #196 Phase G) | GeometryCoefficients (geom-only, 1D+reduced); CollisionCache (σ_t-dependent); DAG-free scan strategies; LD non-affine-scannable paths | DESIGN/IMPLEMENTATION | Move to cache.py or solver initialization docstring (single source of truth for cache construction logic); keep inline summary | L |
| 1068–1078 | Comment: PR-TYPED-2 shim retirement | Transient read-through properties (sig_t, sig_a, etc.) retired; consumers now read mat_xs directly | HISTORY | Record in PR-TYPED series memory | M |

### Class SNSolver

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|-----------|---------|----------------------|-----|
| 869–901 | `SNSolver` class docstring | Unified eigenvalue solver; operator triple construction; inner solver choices (SI vs Krylov); boundary conditions (reflective default, configurable) | CONTRACT | Stays: class interface | H |
| 932–941 | Comment: SI BOUNDARY splitting (eigenvalue inner, #218) | #218 closed eigenvalue-inner gap; warm-start; G-S rate benefit modest; keff_tol-tight snapshots; default Jacobi; gauss_seidel opt-in | HISTORY | Record in stencil-assembly campaign (Phase 3 boundary-G-S enablement) | M |
| 949–963 | Comment: Issue #197 PR-TYPED-0/1 material/ng consolidation | Materials + ng from sn_mesh (single source of truth); XS read via mat_xs field accessors; one PR-TYPED migration window shim retired | HISTORY | Record in PR-TYPED series memory | M |
| 965–969 | Comment: Canonical XS state collapse to MaterialXSField | Per-material Mixture + per-cell views; seven separate per-XS attributes → single accessor; every operator reads through mat_xs | HISTORY | Record in PR-TYPED memory | M |
| 1011–1036 | Comment: Operator triple construction (Wave E Round 2) | (L, S, F) algebra-of-record; constructed once at __init__; downstream consumers see consistent triple; consistent with legacy Wave A-D paths | HISTORY/DESIGN | Move to theory/solver.rst "Operator triple construction" or keep inline as implementation rationale; link to Wave E campaign | M |
| 1080–1099 | `rebind_cross_sections` docstring | Rebuild CollisionCache only (Stratum 2); GeometryCoefficients survives; depletion/thermal-feedback consumer path | CONTRACT | Stays: interface for cross-section rebinding | H |
| 1123–1138 | `compute_fission_source` docstring | χ·(νΣ_f·φ)/k delegator; fission operator is linear; 1/k lives at solver level; PR-INDEX-5 principled layout | CONTRACT | Stays: interface | H |
| 1154–1177 | `compute_group_production_rate` docstring | Per-group volume-integrated production (fission + n2n); spectral analysis intermediate; diagnostic (not on keff path); Issue #196 PR-INDEX-5 layout | CONTRACT/DESIGN | Stays: interface + volume measure convention (Issue 9.6 wiring) | H |
| 1204–1218 | `compute_group_absorption_rate` docstring | Per-group absorption rate; volume-integrated; denominator component of keff; PR-INDEX-5 principled layout | CONTRACT | Stays: interface | H |
| 1227–1250 | `compute_production_rate` docstring | **SPLIT:** (a) 1227–1246 = DESIGN RATIONALE (fission + n2n components, IntegratedReactionRate type, single source, scale anchor ERR-052); (b) 1247–1250 = DESIGN (role split R7, renormalisation vs k numerator, (n,2n) gain side) | SPLIT | (a) move to theory/solver.rst "Production and eigenmode scale"; (b) → operator-algebra.rst "Role split (R7)" | M |
| 1258–1313 | `compute_keff` docstring | **SPLIT:** (a) 1258–1301 = DESIGN RATIONALE (denominator assembly, no-forgetting principle, balance identity, role split R7); (b) 1302–1313 = CONTRACT (math formula, lattice/no-n2n special cases) | SPLIT | (a) → theory/solver.rst "k-effective estimation (keff)"; (b) stays | M |
| 1315–1409 | `_boundary_leakage_rate` docstring | **SPLIT:** (a) 1315–1352 = CONTRACT (net outflow math, vacuum structural zero, scale bridge); (b) 1353–1360 = DESIGN (leakage term history, #291 omission closure) | SPLIT | (a) stays; (b) → theory/solver.rst "Leakage correction (#291)" | M |
| 1411–1447 | `_face_area_of` docstring | Spatial measure per boundary face; 1D vs d≥2 conventions; axis order; broadcast semantics; #291 estimator reference | CONTRACT | Stays: interface for face area computation | H |

### Inner solvers

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|-----------|---------|----------------------|-----|
| 1466–1535 | `_solve_source_iteration` docstring | **SPLIT:** (a) 1466–1496 = DESIGN (Wave E Round 2, operator triple, Carlson seed threading, Krylov comparison); (b) 1497–1535 = HISTORY (Scope narrative, 2-D SI closure, legacy B1 block retirement, verification results) | SPLIT | (a) move to theory/solver.rst "Within-group source iteration"; (b) record in Phase A reference or campaign memory | M |
| 1548–1581 | Comment: Build composite RHS (R-1 Step 4 A1) | Per-ordinate factory; /W producer-side; Pattern 7; legacy broadcast retired; legacy partner-flux seeding retired; zero-boundary convention | HISTORY/IMPLEMENTATION | Move narrative (Wave O O.4a.2 seeding, B.5.2 type semantics) to theory/solver.rst; keep technical comment inline | M |
| 1583–1598 | Comment: Build within-group system | #218 eigenvalue inner honors inner_schedule; Phase-5a windowing; single source of truth builder | IMPLEMENTATION/DESIGN | Keep inline: implementation plumbing | L |
| 1603–1628 | Comment: Warm start (composite / coupled pair) | SourceIteration threading; initial_guess kwarg; B.5.2 cold-start; Phase 5a windowed vs unwindowed; B.2d coupling | IMPLEMENTATION/DESIGN | Keep inline: implementation contract | L |
| 1630–1658 | Comment: Coupled source state + residual certificate | #282 route (a); q½ fold; System B member; windowed exemption; mode-12 closure; residual mint un-built LD widening | HISTORY/IMPLEMENTATION | Move #282/Mode-12 narrative to theory/solver.rst; keep inline loop contract | M |
| 1665–1696 | Comment: Scalar flux reduction + scale partner | Windowed moment vs unwindowed full-angular; integrate_angular vs scalar_flux accessor; ℓ=0 Funk-Hecke; #291 scale bridge | IMPLEMENTATION/DESIGN | Move #291 scale bridge to theory/solver.rst; keep inline reduction logic | M |
| 1704–1742 | `_solve_krylov` docstring | **SPLIT:** (a) 1704–1729 = DESIGN (GMRES on (L+C-S), operator triple, sweep preconditioner, within-group fission zero); (b) 1730–1742 = HISTORY (R-1 Step D narrative, 2-D Cartesian unblocked, scope narrative) | SPLIT | (a) move to theory/solver.rst "Krylov acceleration"; (b) record in Wave-R plan | M |
| 1756–1776 | Comment: Build composite RHS for Krylov | R-1 Step 4 A1; per-ordinate factory; /W normalization; D-H.1c stage 2; zero-boundary convention; reflective-BC threading | IMPLEMENTATION/HISTORY | Move factory narrative to theory/solver.rst; keep inline loop comment | M |
| 1778–1819 | Comment: Build within-group system + warm start | GMRES restart sizing (full ravel); B.5.2 cold-start flux composite; ERR-053 truncation hazard; coupled pair handling | IMPLEMENTATION/DESIGN | Move ERR-053 narrative to theory/solver.rst; keep inline restart logic | M |
| 1867–1882 | Comment: _make_sweep_preconditioner retirement | Phase 3 sub-step 3c; packed-vector wrapper retired (Wave A philosophy); KrylovAcceleration consumer replaces; L+C composite wrapper | HISTORY | Record in stencil-assembly campaign memory | M |

### Helpers and public API

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|-----------|---------|----------------------|-----|
| 1938–1957 | Comment: Retired RHS builders (_build_rhs_cartesian/spherical/cylindrical) | Legacy BiCGSTAB FD-operator path; G1 Krylov replaces; (2*l+1) duplicate retired; P1.7 moment-space plan | HISTORY | Record in moment-space campaign; cite test reference (test_no_legacy_eq_map_or_decoder_in_g1_path) | M |
| 1971–2045 | `_reflect_outflow_into_inflow` docstring | **SPLIT:** (a) 1976–2012 = DESIGN (bare sweep reflective law; wave O narrative; RULING P1; B_a+B_b direct sum, twin matvec/SI routes); (b) 2013–2022 = CONTRACT (face subset restriction, G-S octant-group use) | SPLIT | (a) move to theory/solver.rst "Reflective boundary application" + RULING P1 reference; (b) stays | M |
| 2028–2046 | Comments inside _reflect_outflow_into_inflow | Trace-only action; #226 step 2 verb migration; radial characteristic; System B B_b operator; presence-guard (seedless mesh) | IMPLEMENTATION/HISTORY | Keep inline: implementation contract (reference RULING P1) | L |

---

## FILE: orpheus/numerics/iteration.py

### Module Level

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|-----------|---------|----------------------|-----|
| 1–138 | MODULE | **CONTRACT by notation.rst crosswalk promise**: Operator algebra (A − Σᵢ gᵢ)ψ = q_ext (fixed-source and eigenvalue forms); binding for SN; iteration primitives and iteration drivers (SourceIteration, KrylovAcceleration, KEigenvalue); shape-agnostic; Carlson seed threading; forward references to power_iteration | **CONTRACT** | Stays: dual-A bridge statement (module head mandated by notation.rst row 8) | H |

### Protocols and early comments

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|-----------|---------|----------------------|-----|
| 165–193 | `SupportsSeededApply` docstring | Static contract for inverse-application operator; canonical seeded ``apply`` signature; initial_guess keyword; per-type documentation | CONTRACT | Stays: protocol interface | H |
| 196–203 | Comment: KeffEstimator injection alias retirement | #259 P1 / R8; injection seam dissolved; hardwired methods; at convergence every consistent estimator agrees | HISTORY | Record in eigenvalue-posing campaign (#259 record) | M |
| 206–220 | Comment: _SeededExactApply adapter | Algebra-closed inverse (identity, permutation, scaled) vs wrap-delegate family; exact inverse → accept-and-drop seed | DESIGN | Move to theory/numerics/iteration.rst "Two kinds of inverse"; keep class inline comment | L |
| 276–292 | Comment: Ravellable protocol | Template-based duck typing; to_flat() + from_flat(); typed flux (TimedFullField) vs bare ndarray; decoupling numerics from transport (L1↛L2) | CONTRACT | Move to theory/numerics/iteration.rst "Ravellable protocol"; keep inline as reference | L |
| 344–359 | `_as_scipy_linop` docstring | **Single ORPHEUS↔scipy Krylov boundary**; carrier-space matvec wrapping; template lift/ravel; bare-ndarray reshape; single source of truth (Cardinal Rule 2) | CONTRACT | Stays: interface boundary specification | H |
| 379–398 | `_DisplacementLeaf` docstring + comment (206–220) | Structural face; consumer-minimal (l2, contraction_ratio); numerics MUST NOT import transport; mirroring _is_ravellable check | CONTRACT | Stays: static protocol contract | H |
| 401–434 | Comment: _flux_displacement_leaf function | Iterate increment Δψ (FluxDisplacement vs bare ndarray); convergence diagnostics (ρ ≈ contraction ratio); coupled block recurse (systems[0] convention) | DESIGN | Move narrative to theory/numerics/iteration.rst "Convergence diagnostics"; keep inline as implementation guide | M |
| 436–441 | Comment: Module-level estimator function retirement | _default_production_estimator / _default_keff_estimator folded into KEigenvalue methods; #259 P1/R8 | HISTORY | Record in eigenvalue-posing campaign | M |

### SourceIteration class

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|-----------|---------|----------------------|-----|
| 449–565 | `SourceIteration` class docstring | **CONTRACT**: Fixed-point iteration for (A − Σᵢ gᵢ)ψ = q_ext; loss operator A minus lagged couplings; invertible operator contract (#226 taxonomy step 3); convergence test (ρ-honest equation residual, step 5 R-5.2/R-5.3); parameters and notes (shape-agnostic, ravellable protocol, Carlson seed threading) | CONTRACT | Stays: class interface and algorithm | H |
| 567–595 | `SourceIteration.__init__` docstring | Apply-guards at construction (eager checks); step operator pre-inverted; convergence diagnostics initialization | CONTRACT | Stays: constructor contract | H |
| 601–708 | `SourceIteration.solve` docstring | **SPLIT:** (a) 601–628 = CONTRACT (parameters, returns, ravellable protocol, zero-initialization); (b) 629–708 = **IMPLEMENTATION (detailed loop comments)** | SPLIT | (a) stays; (b) implementation stays inline (ρ-honest stop logic) | H |
| 629–682 | Comment: Ravellable protocol + q_norm scale | R-1 Step 4a; typed flux + bare ndarray handling; zero-initialization; q_norm guard (zero source → absolute test) | IMPLEMENTATION/DESIGN | Keep inline: implementation contract for convergence scale | L |
| 670–707 | Comment: Free-identity STOP + displacement diagnostics | ρ-honest stop from free identity; exact-M assumption; #282 lag-death blindness; equation-relative normalization; moment-space moment-increment class exemption (r3) | IMPLEMENTATION/DESIGN | Move #282/Mode-12 narrative to theory/numerics/iteration.rst; keep inline loop logic | M |

### KrylovAcceleration class

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|-----------|---------|----------------------|-----|
| 716–819 | `KrylovAcceleration` class docstring | **CONTRACT**: GMRES on (A − Σᵢ gᵢ)ψ = q_ext; composed matvec; algorithm (GMRES), optional preconditioner; algebra-of-record; preconditioner parameter rename (R-1 Step B); parameters, raises, notes | CONTRACT | Stays: class interface and algorithm | H |
| 821–841 | `KrylovAcceleration.__init__` docstring | Apply-guards; preconditioner selection (None → A.inverse().apply if invertible, else unpreconditioned); construction-time pinning | CONTRACT | Stays: constructor contract | H |
| 862–903 | `KrylovAcceleration.solve` docstring | Parameters (q_ext, initial_guess), returns (psi, residual_history); ravellable protocol (ravel→scipy, unravel back); preconditioned residual norm callback | CONTRACT | Stays: method interface | H |
| 910–946 | Comment: GMRES setup + callback | Loss matvec composition; scipy linop wrapping; preconditioner scipy linop; initial guess handling; callback for preconditioned-residual history | IMPLEMENTATION | Keep inline: GMRES setup logic | L |
| 947–953 | Comment: scipy GMRES invocation note | No try/except; TypeError must surface; scipy>=1.14 floor; tol fallback arm retired (B.5.2 regression hazard) | IMPLEMENTATION/HISTORY | Keep inline: GMRES contract; record B.5.2 regression lesson | L |
| 963–1001 | Comment: GMRES info flag surfacing + exact-breakdown carve | Non-convergence warning (info≠0); exact-breakdown detection (final rk=0.0); singular-but-consistent system hand-warm-start; B.2d coupled matvec transitional motivator | IMPLEMENTATION/DESIGN | Move non-convergence warning narrative to theory/numerics/iteration.rst; keep inline logic + exact-breakdown guard | M |

### KEigenvalue class

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|-----------|---------|----------------------|-----|
| 1011–1098 | `KEigenvalue` class docstring | **CONTRACT**: k-eigenvalue problem (A − S)ψ = Fψ/k; operator-triple realization; delegates to power_iteration (canonical loop); outer/inner layered posing (A_loss = A − S, M = F, k = μ); boundary methods (compute_fission_source, solve_fixed_source, keff/production estimators); convergence test (dominance ratio); parameters, raises, notes | CONTRACT | Stays: class interface and algorithm | H |
| 1100–1141 | `KEigenvalue.__init__` docstring | Construct-time A.invertible validation (posing layer builds inverse); eigenvalue_method selector ("power" only); inner SourceIteration building (ONE inner); F apply-guard | CONTRACT | Stays: constructor contract | H |
| 1168–1182 | Comment: EigenvalueSolver boundary realization | Layer-2 k-posing; ONE power-iteration loop in codebase (Cardinal Rule 2); KEigenvalue and SNSolver both implementers; carries boundary methods | DESIGN/REFERENCE | Move boundary-interface narrative to theory/eigenvalue-posing.rst; keep inline as reference | L |
| 1184–1199 | `KEigenvalue.initial_flux_distribution` docstring | Caller-supplied guess stashed by solve(); boundary invoked outside solve = caller error | CONTRACT | Stays: method contract | H |
| 1201–1205 | `KEigenvalue.compute_fission_source` docstring | Outer eigen-source F·ψ/k (k-posing's M = F) | CONTRACT | Stays: method interface | H |
| 1207–1221 | `KEigenvalue.solve_fixed_source` docstring | Resolvent A_loss⁻¹ q via inner SourceIteration; warm-started from previous flux; single gain S (zero within-group fission) | CONTRACT | Stays: method interface | H |
| 1223–1242 | `KEigenvalue.compute_production_rate` docstring | Production-rate normalisation (unit production per step); ERR-052 scale anchor; convention ∫νΣ_f·φ·dV=1; hardwired since #259 P1/R8 | CONTRACT/HISTORY | Split: production formula (stays); ERR-052 rationale + #259 retirement (→ theory/eigenvalue-posing.rst) | M |
| 1244–1273 | `KEigenvalue.compute_keff` docstring | Rayleigh estimator (hardwired, #259 P1/R8); operator-form spelling; leakage-inclusive (never had #291 omission); consistent estimator agreement at convergence; injection seam retired | CONTRACT/HISTORY | Split: operator formula (stays); #291 leakage commentary (→ theory/eigenvalue-posing.rst); #259 injection retirement (→ campaign record) | M |
| 1275–1290 | `KEigenvalue.converged` docstring | ≥3 iterations, dk < keff_tol AND dφ < flux_tol | CONTRACT | Stays: convergence criterion | H |
| 1292–1334 | `KEigenvalue.solve` docstring | Run eigenvalue via canonical power_iteration loop; boundary realization; eigenvalue_method selector (FEAST hook reserved); parameters (initial_guess REQUIRED), returns (keff, history, ψ) | CONTRACT/DESIGN | Stays: method interface + canonical algorithm reference | H |

---

## Summary

### Verdict counts by file

**solver.py (orpheus/sn/solver.py)**
- **CONTRACT (keep)**: 31 blocks (~1,100 lines docstring + inline)
- **HISTORY (move to docs/theory or campaign memory)**: 17 blocks (~180 lines)
- **MOVED (create new sections)**: 8 blocks (design rationale + narrative branches)
- **COMMENT-cut (shorten/move)**: 6 blocks (implementation comments; most kept inline)

**iteration.py (orpheus/numerics/iteration.py)**
- **CONTRACT (keep)**: 19 blocks (~550 lines docstring + protocol)
- **HISTORY (move to docs/theory or campaign memory)**: 5 blocks (~30 lines)
- **DESIGN (move to theory or keep inline)**: 3 blocks (~50 lines)
- **IMPLEMENTATION (keep inline)**: 8 blocks (~100 lines)

### Docstring lines proposed: KEEP vs CUT

| file | total lines | keep (CONTRACT) | move (HISTORY/DESIGN) | confidence |
|------|-------------|-----------------|----------------------|-----------|
| solver.py | ~2,100 | ~1,450 | ~650 | M |
| iteration.py | ~550 | ~520 | ~30 | H |
| **combined** | **~2,650** | **~1,970** | **~680** | **M** |

---

## Posing-drift checks (module head `iteration.py`)

Searched lines 1, 7, 32, 861 for drift from notation.rst standard `(A − Σᵢ gᵢ)ψ = q_ext`:

| line | content | status |
|------|---------|--------|
| 1 | `(A - \sum_i g_i)` ✓ | conformant |
| 7 | `(A - \sum_i g_i)` ✓ | conformant |
| 8 | `\frac{1}{k}\ F\ \psi` ✓ | conformant (eigenvalue form) |
| 12 | `(A - \sum_i g_i)` ✓ | conformant |
| 738 | `(A - \sum_i g_i)` ✓ | conformant (KrylovAcceleration) |

**No posing-drift detected.** Module head consistently uses variadic `(A − Σᵢ gᵢ)` record (not the `(A − S − F)ψ = q_ext` form).

---

## Uncertain blocks (conf=L)

1. **solver.py line 389–393** — _CERTIFICATE_SAFETY guard constant. Threshold value (10.0) is engineering; the "FP-reassociation grain" rationale belongs in solver.rst. Confidence L because the inline comment is terse and the value is context-specific.

2. **solver.py line 1038–1066** — Two-stratum cache. The affine-recurrence/DAG-free dispatch logic is implementation-heavy and could move to solver initialization docstring; keeping inline as summary works if the full narrative moves to theory.

3. **iteration.py line 206–220** — _SeededExactApply adapter. The "algebra-closed inverse" dichotomy is elegant but belongs in operator-algebra theory; inline comment is a brief reference.

4. **iteration.py line 276–292** — Ravellable protocol. The L1↛L2 decoupling constraint is architectural policy; inline comment is correct but the full protocol narrative should live in theory/numerics/iteration.rst.

---

## Recommendations (summary)

1. **Immediate (no change needed)**:
   - All **CONTRACT** blocks stay inline (solver.py + iteration.py).
   - **IMPLEMENTATION** comments (loop bodies, GMRES setup, etc.) stay inline.

2. **Deferred to Phase 2 II**:
   - HISTORY blocks (campaign narration, retirement records) → `.claude/plans/*/` memory files or theory chapter "Development history" sections.
   - DESIGN blocks (role split R7, RULING P1, Mode-12, #282, etc.) → appropriate theory chapters with cross-references.

3. **Refactor guidance**:
   - `evaluate_residual` (solver.py) → split docstring into CONTRACT (top) + HISTORY sidebar.
   - `_boundary_leakage_rate` (solver.py) → move scale-bridge and #291 narrative to theory; keep math contract inline.
   - `compute_keff` methods (solver.py + iteration.py) → formulae stay; estimator-selection rationale moves to theory.

---

**File generated:** 2026-07-22  
**Classification discipline:** Verdicts applied per Phase 2 C task spec + notation.rst dual-A bridge rule.
