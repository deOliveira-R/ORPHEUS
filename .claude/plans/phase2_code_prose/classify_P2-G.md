# Phase 2 Code-Prose Classification: streaming.py, augmented_mesh.py, boundary.py

Inventory of docstrings ≥10 lines and comment blocks ≥5 lines across three SN operator files. Each row proposes a verdict: **CONTRACT** (essential for calling/modifying), **TWIN** (duplicates theory book content), **MOVED** (teaching not in book; name target chapter), **HISTORY** (campaign narration), **COMMENT-cut** (narrative comment, not constraint).

---

## File 1: orpheus/sn/operators/streaming.py (1162 lines)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 1–40 | `module` | Four-operator algebra leaves: L, C, (L+C), carrier contract | CONTRACT | — operator algebra | H |
| 42–60 | `module` | Wave E/D-H/D-I/D-J retired helpers; EquationMap codec family | HISTORY | git log `a614610` / commit history | H |
| 61–80 | `module` | Symmetric-closure invariant: DD upwind, curvilinear avg, τ-weighted interp | TWIN | docs/theory/loss_representation.rst / symmetric closure note | H |
| 81–120 | `module` | Boundary-condition handling: bare sweep, psi.inflow defect, reflective coupling, -B sibling | CONTRACT | docs/theory/boundary_conditions.rst + :ref:`bc-extraction` ¶ O.4a.2 | H |
| 181–211 | `_require_typed_composite()` | Shared matvec input contract: timeless composite + mesh identity, D-I.3d typed guard | CONTRACT | — call guard | H |
| 232–243 | `StreamingOperator` class | Pure L operator: streaming + angular redistribution, σ-free structure | CONTRACT | — class contract | H |
| 244–269 | `StreamingOperator` class | Pure σ-free apply via loss_representation.streaming_action, affinity in σ, (L+C) recovery | CONTRACT | — apply semantics | H |
| 271–284 | `StreamingOperator` class | Why L carries NO σ (Pattern 4): Carlson pole seed σ-dependence cancels into C | CONTRACT | — design rationale | H |
| 286–303 | `StreamingOperator` class | Capability set: pure L not invertible, (L+C) invertible, apply_transpose available Wave O | CONTRACT | — capability surface | H |
| 304–312 | `StreamingOperator` class | Parameters: sn_mesh sole parameter, no σ | CONTRACT | — constructor contract | H |
| 313–322 | `StreamingOperator` class | Depth B D-I/D-J deprecations: packed-vector codec retired, typed composite only | HISTORY | git log D-I/D-J waves | H |
| 340–358 | `is_adjointable` property | Two-factor honest factorization: kernel (DD yes/LD no) + orientation (1-D yes/multi-D DAG/wavefront no) | COMMENT-cut | — can extract to theory "adjoint feasibility" section | M |
| 361–400 | `domain` property | Composite V_bulk ⊕ V_trace (Wave O O.2b), full operator, block-diagonal metric for G-adjoint | CONTRACT | — domain algebra | H |
| 407–450 | `apply()` | Pure σ-free L·ψ via streaming_action, σ-freedom byte-stable probe, no cell-collision term | CONTRACT | — apply contract | H |
| 454–484 | `apply_transpose()` | Hilbert transpose L^T·φ via streaming_action_transpose, σ-free, Euclidean (not metric-conjugated) | CONTRACT | — adjoint apply contract | H |
| 489–508 | `loss_representation` property | THE loss-operator representation (S6.5, L21), ONE instance shared by matvec + sweep | CONTRACT | — polymorphic dispatch anchor | H |
| 521–537 | `__add__()` | Compose L + X → InvertibleOperator iff X is collision multiplier; dispatch on L only (one-way) | CONTRACT | — algebra overload rule | H |
| 554–584 | `InvertibleOperator` class | (L+C)^-1 ≈ WDD sweep algebraic identity, MRO type-as-structure, owns full action algebra | CONTRACT | — algebraic identity + typing discipline | H |
| 586–610 | `InvertibleOperator` class | Construction: two paths (L + C dispatch or explicit), structurally identical, call-site readability | CONTRACT | — constructor paths | H |
| 612–622 | `InvertibleOperator` class | Capability set: is_invertible=True, apply_transpose via OperatorSum closure law OVERRIDDEN to composite's own | CONTRACT | — capability surface | H |
| 623–650 | `InvertibleOperator` class | solve() API: timeless rhs carrier (W-C P4.5), history-bearing TimedFullField inherited, curvilinear ψ½ direct seed | CONTRACT | — solve contract | H |
| 652–671 | `InvertibleOperator` class | Parameters & Notes: streaming + diagonal inputs, mesh-identity validation, σ > 0 guard | CONTRACT | — construction contract | H |
| 677–715 | `__init__()` (InvertibleOperator) | Type/instance guards (StreamingOperator, MultiplicationOperator), mesh-identity invariant (kernel geometry), σ > 0 validation | CONTRACT | — guard contract | H |
| 769–796 | `apply()` (InvertibleOperator) | Matvec (L+C)·ψ = M(σ)·ψ OWNS composite action; #240 Step B single-sources σ from diagonal, not leaf sum | HISTORY + CONTRACT | #240 Step B / commit history | H |
| 802–821 | `apply_transpose()` (InvertibleOperator) | Adjoint (L+C)^T·φ = M(σ)^T·φ OVERRIDES OperatorSum, Euclidean (Wave O O.2b metric via _AdjointOperator) | CONTRACT | — adjoint action | H |
| 865–884 | `inverse()` | Returns SweepOperator inverse; BIT-IDENTICAL to solve(); coextensive with ScheduledInvertibleOperator (defer-until-≥2 trigger) | HISTORY + CONTRACT | #226 step 2 / Wave 2–3 coexistence | M |
| 896–950 | `solve()` | Invert (L+C)·ψ = rhs via WDD sweep, ψ½ direct seed (#282 route a 2.5d), accept-and-drop initial_guess | CONTRACT | — direct-inverse contract | H |
| 970–1002 | `_solve_timed_full_field()` | Composite TimedFullField body (D-H.1c stage 1), runs WDD on loss_representation (S6.5 L21), L2 boundary plumbing | CONTRACT + HISTORY | D-H / S6.5 wave | H |
| 1016–1022 | `_solve_timed_full_field()` comment | W-C (P4.5): timeless FullField contract, legacy AngularFlux retired D-H.2-C3, single guard site both solve/solve_moments | HISTORY | D-H / P4.5 wave | H |
| 1037–1043 | `_solve_timed_full_field()` comment | L2 boundary buffer (D-H.2-C2), mutable write-through, sweep mutates in place, result re-wrapped | COMMENT-cut | — implementation detail | L |
| 1044–1051 | `_solve_timed_full_field()` comment | Wave O O.4a.2 BARE SWEEP: inflow seed = rhs.boundary (GIVEN), no bc.apply, ψ½ direct (#282 route a 2.5d) | HISTORY | Wave O O.4a.2 | H |
| 1062–1072 | `_solve_timed_full_field()` code block | Sweep on operator's ONE representation (S6.5 #222), SAME instance as matvec, L21 type fact | HISTORY | S6.5 / #222 / commit history | H |
| 1114–1138 | `solve_transpose()` | Invert (L+C)^T·x = b via REVERSE-SCAN (#280 Phase 2.5b), analytic adjoint, (L+C)^T = L^T + C diagonal | CONTRACT + HISTORY | #280 Phase 2.5 / Wave O | H |
| 954–960 | code comment | solve_moments (Phase 5c) retired #226 step 2, moment-emitting via windowed product P @ A.inverse() | HISTORY | #226 step 2 / Phase 5c | H |

**Summary for streaming.py:**
- **CONTRACT:** 21 blocks keeping docstrings (operator algebra, guards, apply/solve/adjoint contracts, capability surface)
- **HISTORY:** 11 blocks moving to campaign archives or commit record (D-I/D-J waves, O.4a.2, #240/#280 step references, Phase 5c retirement)
- **TWIN:** 1 block (symmetric-closure invariant moves to loss_representation theory if not present)
- **COMMENT-cut:** 2 blocks (narrative comments on implementation detail, factor decomposition)
- **Docstring lines keep:** ~850; proposed cut/moved: ~80

---

## File 2: orpheus/sn/mesh/augmented_mesh.py (1606 lines)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 1–15 | `module` | Augmented geometry for SN, axis-primary (C5.1), precomputes streaming stencil, 3 coordinate systems | CONTRACT | — module scope | H |
| 76–82 | TYPE_CHECKING comment | NOTE (B.5.A): transport-field TYPE_CHECKING imports retired, mesh supplies shape only (not types) | HISTORY | B.5.A wave / commit history | H |
| 97–175 | `SNMesh` class | Axis-primary (C5.1), precomputes streaming stencil, legacy Mesh1D/2D adapter, mesh-identity invariant | CONTRACT | — class contract | H |
| 208–217 | `__init__()` comment | Legacy inbound surface (C5.1 axis-primary #225), convert Mesh1D/2D to axes ONCE at boundary | HISTORY | C5.1 / #225 inversion | H |
| 241–258 | `_init_core()` comment | Method-agnostic DATA block → MaterialMesh._init_data (ng-consistency, validation at construction) | HISTORY + COMMENT-cut | C5.1 design narrative | M |
| 261–273 | `_init_core()` comment | SN method layer (BEHAVIOR atop MaterialMesh data): quadrature, scheme, pole closure | COMMENT-cut | — design structure narrative | L |
| 274–279 | `_init_core()` comment | Cell-update strategy (Wave D Round 2 #161), DD default, LD/EC/Step deferred, user pass at construction | HISTORY | Wave D #161 | H |
| 280–330 | `_init_core()` comment | Angular-redistribution closure (Issue #168 Phase D default flip), MM canonical (Hébert §3.9.4), Carlson seed (#282 2.5d) | HISTORY | #168 Phase D / #282 route (a) | H |
| 339–361 | `_init_core()` comment | Stencil dispatch by coordinate system, reduced operator (shared MoC/CP), slab variant for completeness | COMMENT-cut | — design narrative | M |
| 389–405 | `_init_core()` comment | Boundary trace + realized laws (C5.3 unconditional, #290 P7b built HERE phase-space substrate, not inside BC resolve) | HISTORY | C5.3 / #290 P7b | H |
| 407–422 | `_init_core()` comment | Resolve BC declarations via ONE shared TransportMethod (#290 P7b), face inventory = BC inventory (C4 #220) | HISTORY | C4 #220 / #290 P7b | H |
| 429–443 | `_init_core()` comment | Pole-angular closure binding (PR-TYPED-6.5 Phase 2.9), class not instance, bind at construction time | HISTORY | PR-TYPED-6.5 Phase 2.9 | H |
| 445–452 | comment | Boundary condition resolution: see resolve_boundary_conditions body (#290 P7b) | HISTORY | #290 P7b | H |
| 459–487 | `realize_boundary_law()` | Realize one typed boundary law on face (TransportMethod hook), build SNMethodSpace, wrap in _BoundBoundaryOperator | CONTRACT | — SN arm of BC realization | H |
| 509–521 | `is_1d` property | True iff ndim == 1 (genuine dimensionality, not ny==1 phantom), single source for 1-D vs multi-D dispatch | CONTRACT | — sweep dispatch criterion | H |
| 524–538 | `is_cartesian` property | True iff curvature is None (Cartesian), orthogonal to is_1d, both axes required for sweep strategy selection | CONTRACT | — coordinate-system criterion | H |
| 541–574 | `streaming()` | Per-axis RAW streaming g = \|μ_axis\|·A_down/V, scheme-agnostic (closure factor owned by DD/LD), Cartesian-only | CONTRACT | — down-face coefficient | H |
| 584–596 | `face_labels` property | Canonical boundary-face inventory (FaceLabel per axis), iteration order = canonical for face flattening (C3 packing) | CONTRACT | — boundary face enumeration | H |
| 600–605 | `face_shape()` | Spatial shape of boundary face (codimension-1 hyperplane), axis-agnostic | CONTRACT | — face geometry accessor | H |
| 609–620 | `face_outflow_ordinates()` | Ordinate indices whose Ω·n is OUTWARD at face (±10^-15 tangent threshold), producer for per-face outflow mask (C3/C4/C5) | CONTRACT | — ordinate filtering | H |
| 647–689 | `from_axes()` classmethod | Axis-native surface (C5.1, axis-primary inversion), canonical endpoints (min/max/outer), synthesize legacy adapter d≤2 | CONTRACT + HISTORY | C5.1 #225 / axis-primary inversion | H |
| 722–751 | `from_material_mesh()` classmethod | Data/behavior join (homogenized MaterialMesh → solvable SN phase space), axes/materials pass-through, re-derive data block | CONTRACT | — promotion path | H |
| 767–783 | `angular_trace` property | Unified boundary AngularTraceSpace (geometry-blind since C5.3 #225), single source for Ω·n per face, inflow/outflow selectors | CONTRACT + HISTORY | C5.3 / #225 | H |
| 786–794 | `radial_characteristic_levels` constant comment | FP-noise guard for R12a strict-interior test, τ_raw ∈ (0,1) seed presence predicate | HISTORY | R12a / #282 route (a) | H |
| 798–821 | `radial_characteristic_levels` property | μ-level indices consuming INDEPENDENT starting-direction state (R12a), τ_raw threshold, sphere-GL instance | HISTORY + CONTRACT | R12a / #282 route (a) | H |
| 834–847 | `_radial_characteristic_for_levels_args()` | Shared args for ψ½ spaces, R12a gate, (ng, nx, cell_volumes), single-sources across split spaces | CONTRACT | — space builder input | H |
| 861–872 | `radial_characteristic_interior_space` property | ψ½ INTERIOR (cells) space (System B's interior block, R12a), SPD metric G_sd = V_cell, cached | CONTRACT | — space definition | H |
| 888–899 | `radial_characteristic_boundary_space` property | ψ½ BOUNDARY (corner) space (System B's boundary block, R12a), G = V(r=R) corner gauge, cached | CONTRACT | — space definition | H |
| 915–936 | `radial_characteristic_field_space` property | System B's interior ⊕ boundary composite space, FullFieldSpace instance (post-eviction one class), identity signal | CONTRACT + HISTORY | B.2d eviction / system-B typing | H |
| 950–976 | `full_field_space` property | Composite V_bulk ⊕ V_trace (Wave O O.2b), G-adjoint metric (V·w_n bulk, \|Ω·n\|·w trace), radial_characteristic member separate (B.2d) | CONTRACT + HISTORY | Wave O / B.2d | H |
| 1004–1026 | `boundary_face_layout` property | Flat FaceLayout of boundary faces (per-geometry name/shape/slot), axis-count generic, no hand-list per geometry | CONTRACT | — face layout descriptor | H |
| 1027–1045 | `boundary_face_layout` property "Spatial-moment tail" | Multi-moment closure (LD bilinear UBLD), trailing 2^(d-1) transverse-moment axis per face (#251 Leg B #247), DD byte-identical | CONTRACT + HISTORY | #251 / #247 | H |
| 1047–1062 | `boundary_face_layout` property notes | Layout contains ONLY boundary face slots (interior cache cells on SweepScratch post-D-G) | CONTRACT | — layout scope | H |
| 1089–1160 | `dag_walk()` | Per-ordinate cell DAG walk (Issue #196 Phase G Step 2.6 Q3), canonical 1-D sweep primitive, XOR ordinate_idx/direction_sign | CONTRACT | — sweep iteration protocol | H |
| 1219–1240 | `dag_walk_cell_indices()` | Lightweight dag_walk twin (yields cell indices only, ~14% matvec time save on slab), cell order matches dag_walk | HISTORY + CONTRACT | profiling result (PR-TYPED-6c) | M |
| 1282–1295 | `_require_mu_level()` | Narrow mu_level_idx to int, cylindrical requires it, slab/sphere pass None, single source contract, Pattern 2 | CONTRACT | — guard contract | H |
| 1303–1311 | `_representative_ordinate()` | Pick non-degenerate ordinate (degenerate \|η\| < 10^-15 excluded), cell ordering depends only on direction_sign + level | CONTRACT | — ordinate selection | H |
| 1352–1373 | `_make_cell_visit()` | Single production site stamping M-M τ/c_in/c_out onto visit (Issue #236 Phase 2 B2/B3), sources from pole_angular_closure, Pattern 2 | HISTORY + CONTRACT | #236 Phase 2 | H |
| 1389–1395 | `_iter_cartesian_visits()` | Yield slab visits (1-D Cartesian) forward/backward sweep direction (Issue #196 Phase G Step 2.5, face_area=1.0 neutral) | HISTORY + CONTRACT | #196 Phase G | M |
| 1419–1424 | `_iter_spherical_visits()` | Yield spherical visits (outward 0→nx-1 downstream outer; inward nx-1→0 downstream inner) | CONTRACT | — sweep order | H |
| 1453–1467 | `_iter_cylindrical_visits()` | Cylindrical visits (within-level azimuthal m index), ordinate_idx within level, degenerate \|η\| < 10^-15 forward regardless | CONTRACT | — sweep order | H |
| 1475–1481 | `_iter_cylindrical_visits()` comment | Pure-azimuthal degenerate (Issue #196 Phase G Step 2.5), face_area=0.0 signals "no spatial flow", replaced None | HISTORY | #196 Phase G | H |
| 1526–1549 | `_setup_cartesian()` | Precompute raw per-axis streaming g = \|μ_a\|/Δa, scheme-agnostic (closure factor by scheme), built over range(ndim) | CONTRACT | — stencil setup | H |
| 1562–1574 | backward-compat comment | Properties route to self.reduced, Wave D Round 2 #12 migration path, DeprecationWarning (not -O error), 6 production sites | COMMENT-cut | — deprecation strategy narrative | L |
| 1578–1584 | `face_areas` property | [Deprecated] Use self.reduced.face_areas | HISTORY | deprecation → removed | H |
| 1590–1597 | `delta_A` property | [Deprecated] Use self.reduced.delta_A | HISTORY | deprecation → removed | H |

**Summary for augmented_mesh.py:**
- **CONTRACT:** 28 blocks (mesh contract, accessors, space definitions, sweep primitives, face enumeration)
- **HISTORY:** 18 blocks (C5.1 inversion, #168/#225/#290/#282/#236/#196 waves, PR-TYPED-6.5, B.2d, deprecation)
- **COMMENT-cut:** 5 blocks (design structure, stencil dispatch, backward-compat strategy)
- **Docstring lines keep:** ~1,000; proposed cut/moved: ~200

---

## File 3: orpheus/sn/operators/boundary.py (848 lines)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 1–49 | `module` | Realized boundary law B (Wave O #208), whole-trace BOUNDARY-block operator, 2×2 block structure, per-face laws | CONTRACT | — operator algebra | H |
| 11–34 | `module` "Block structure" | Block matrix visual (C/S/F BULK, L_full FULL, B BOUNDARY), per-system boundary direct sum (RULING P1) | CONTRACT | — block algebra | H |
| 92–100 | `_RULED_CORNER_KINDS` comment | Outer-face law kinds with RULED ψ½ corner action (RULING P1, R12a), vacuum/reflective only, white/albedo/periodic deferred | HISTORY | R12a / 2.5d plan-of-record | H |
| 104–114 | `_zero_bulk_source()` | Zero-bulk A_ss carrier, sized from mesh (full-angular OR harmonic-moment via scheme.spatial_basis_per_axis), single source (Pattern 2) | CONTRACT | — carrier builder | H |
| 124–174 | `SNBoundaryOperator` class | B_a (System A trace boundary), block-diagonal per-face laws, zero bulk + reflected trace, on composite full_field_space | CONTRACT | — class contract | H |
| 136–142 | `SNBoundaryOperator` class | On seed-carrying: B = B_a + B_b (System A + System B boundaries), block-diagonal (RULING P1), off-diagonal = face physics | HISTORY + CONTRACT | RULING P1 | H |
| 184–191 | `_face_laws` property | Map boundary face → per-face realized law from sn_mesh.bc, keys by trace.layout.faces (C4 shared derivation), single source | CONTRACT | — law accessor | H |
| 200–207 | `is_adjointable` property | B = ⊕ per-face laws; composite adjoint exists iff EVERY face law is (white NOT adjointable) | COMMENT-cut | — composition predicate logic | L |
| 211–216 | `domain` property | Composite carrier (NOT bare trace), full FullField for OperatorSum guard, block-diagonal G-adjoint metric | COMMENT-cut | — domain scope reasoning | L |
| 228–263 | `_reflect_trace()` | Core A_ss trace action, apply per-face law, project onto inflow rows (consistency row), single source of truth (Pattern 2) | HISTORY + CONTRACT | #257 S8a / O.4a.2 / Phase C insight | H |
| 267–288 | `_reflect_trace()` comment | Phase 3 G-S schedule (subset faces/rows), block-diagonal exactness (no cross-face), row-granular restriction per face | COMMENT-cut | — schedule semantics narrative | M |
| 333–344 | `_apply_faces()` | Lift trace-only reflect onto full FullField with zero bulk (B_a System A trace block), ray-corner separate B_b (RULING P1, B.2d eviction) | HISTORY + CONTRACT | #257 S8a / RULING P1 / B.2d | H |
| 381–400 | `reflect_into_inflow()` | Trace-only forward reflection B·ψ.outflow→inflow row (AngularBoundarySourceSink), used by direct SI/reconstruction sweep | CONTRACT | — trace-only reflection | H |
| 409–425 | `reflect_inflow_inplace()` | MUTATING reflect_into_inflow, phase-3 G-S signature, both SI loop + reconstruction sweep call this (#226 step 2, RULING P1) | HISTORY + CONTRACT | #226 step 2 / RULING P1 | H |
| 443–456 | `split()` | Split B = B_lower + B_upper under schedule (regular matrix splitting G-S, #226 §17 W2), exact partition inflow rows | HISTORY + CONTRACT | #226 §17 W2 | H |
| 475–481 | `apply_transpose()` | Euclidean Bᵀ·ψ reachable iff every face law is adjointable (white deferred to B.H under \|Ω·n\|·w metric, Wave O step O.2) | CONTRACT + HISTORY | Wave O step O.2 | H |
| 486–540 | `RadialCharacteristicBoundaryOperator` class | B_b (System B ψ½ ray-corner boundary), System B's own carrier, off-quadrature μ = ±1 ray corner action (RULING P1, one law two carriers) | CONTRACT + HISTORY | RULING P1 / B.2b re-type | H |
| 512–519 | RadialCharacteristicBoundaryOperator capabilities | apply always, apply_transpose iff outer-ray law is adjointable (vacuum/reflective yes; white/albedo/periodic deferred) | CONTRACT | — capability surface | H |
| 521–531 | RadialCharacteristicBoundaryOperator adjoint metric | Euclidean apply_transpose = G_sd Hilbert adjoint (symmetric corner gauge g₊=g₋=V(R), RULING P2) | HISTORY + CONTRACT | RULING P2 / B.2b | H |
| 584–610 | `_reflect_corner()` | A_ss CORNER action (R13, 2.5d), specular (μ=+1 mirror = μ=-1), vacuum (zero), white/albedo/periodic deferred (no ruled off-quadrature action) | HISTORY + CONTRACT | R13 / 2.5d plan-of-record | H |
| 616–626 | `_reflect_corner()` comment | Seed-carrying = 1-D curvilinear (exactly one outer face xmax), per-KIND law dispatch (not operator composition tree) | COMMENT-cut | — mesh geometry narrative | L |
| 641–648 | `_apply_faces()` (RadialCharacteristicBoundaryOperator) | B_b on System B's own carrier (ray-corner action + zero-source interior), B.2b re-type no bulk/trace padding (Pattern 4) | HISTORY + CONTRACT | B.2b re-type / Pattern 4 | H |
| 661–672 | `_apply_faces()` comment | Role parse (boundary MEMBER = flux corner, #289 F2 discipline), source-role member is caller error | COMMENT-cut | — role-parsing narrative | L |
| 691–696 | `apply_transpose()` | Euclidean Bᵀ·ψ½ mirror-image corner swap, reachable iff adjointable, Euclidean = G_sd Hilbert (RULING P2 symmetric gauge) | CONTRACT + HISTORY | RULING P2 / symmetric gauge | H |
| 703–713 | `reflect_corner_inplace()` | MUTATING reflect_corner, System B's corner analogue of SNBoundaryOperator.reflect_inflow_inplace (#282 route (a)), final sweep + direct SI call BOTH | HISTORY + CONTRACT | #282 route (a) / RULING P1 | H |
| 723–727 | code comment | B.2b transient _RayBoundaryFullFieldGain adapter RETIRED at B.2d (driver CoupledField, B_b native at (B,B) grid slot) | HISTORY | B.2d / adapter retirement | H |
| 731–750 | `SNMaskedBoundaryOperator` class | One half B_lower/B_upper (schedule split #226 §17 W2), row-restricted projection per face inflow ordinates, NOT invertible | CONTRACT + HISTORY | #226 §17 W2 | H |
| 761–767 | `__init__()` (SNMaskedBoundaryOperator) docstring | Three attributes: inner (whole-trace law), rows (per-face inflow ordinates), schedule (octant-order) | CONTRACT | — constructor contract | H |
| 788–817 | `reflect_rows_inplace()` | ADDITIVE per-face-per-row inflow overwrite (reified forward substitution #226 §17 W2), M·z exact on inflow rows only | HISTORY + CONTRACT | #226 §17 W2 | H |
| 843–847 | `BoundarySplit` NamedTuple | Named pair (lower, upper) from SNBoundaryOperator.split | CONTRACT | — named return type | H |

**Summary for boundary.py:**
- **CONTRACT:** 17 blocks (operator block structure, apply/split/reflect contracts, space definitions, reflected action)
- **HISTORY:** 14 blocks (Wave O #208, #226/#257/#282/#290 steps, RULING P1/P2, B.2d/B.2b eviction, 2.5d plan-of-record)
- **COMMENT-cut:** 4 blocks (property reasoning, schedule/role narrative)
- **Docstring lines keep:** ~550; proposed cut/moved: ~120

---

## Posing Flags

The **operator-algebra posing** record is `A = L + C − S − B` (honest within-group algebra). Dated poses detected:

| file | lines | operator family | dating | context | confidence |
|------|-------|-----------------|--------|---------|------------|
| streaming.py | 1–4 | L, C, (L+C), -B | Phase G | module docstring | H |
| streaming.py | 554–584 | (L+C)^-1 ≈ WDD | Wave 2–3 | InvertibleOperator class | H |
| augmented_mesh.py | — | — | — | (operator algebra in boundary/streaming files, not here) | — |
| boundary.py | 1–49 | B whole-trace block | Wave O | module docstring | H |
| boundary.py | 486–540 | B_a + B_b direct sum | B.2d | RadialCharacteristicBoundaryOperator class | H |

All poses are **correctly stated**. No ERR-NNN gotchas found in prose blocks.

---

## Uncertainty Notes (conf L)

1. **streaming.py line 340–358** (`is_adjointable` comment): narrative explanation of two-factor honest predicate; marginal whether kernel-vs-orientation decomposition belongs in docstring vs. theory page on adjoint feasibility. Keeping as CONTRACT is conservative; could move if the decomposition is theoretically motivated. Confidence L.

2. **augmented_mesh.py lines 241–258** (`_init_core` comment): design narrative (data/behavior split) is illustrative but not essential to understanding construction. Recommending COMMENT-cut if space is urgent; confidence M (mid-confidence that this is purely narrative).

3. **augmented_mesh.py lines 1219–1240** (`dag_walk_cell_indices` docstring): flagging the ~14% time-save profiling result as HISTORY context rather than CONTRACT. The method contract (lightweight variant, yields cell indices only) is essential; the profiling context is design narrative. Splitting: CONTRACT core + HISTORY cite. Confidence M.

4. **boundary.py lines 267–288** (`_reflect_trace` comment): Phase 3 Gauss-Seidel schedule semantics (face/row subsetting, block-diagonal exactness); marginal whether this belongs in code or in a theory/algorithm page on schedule-split reflection. Recommending COMMENT-cut if algorithm docs exist elsewhere. Confidence M.

---

## Summary Counts

### By Verdict

| Verdict | Count | Lines | Keep? |
|---------|-------|-------|-------|
| CONTRACT | 66 | ~2,400 | ✓ All |
| HISTORY | 43 | ~250 | → Campaign archive / commit log |
| COMMENT-cut | 11 | ~80 | → Remove (narrative only) |
| TWIN | 1 | ~20 | → Loss_representation.rst if not present |
| MOVED | 0 | — | — |

### By File

| File | CONTRACT | HISTORY | COMMENT-cut | TWIN | Total Lines | Keep % |
|------|----------|---------|-------------|------|-------------|--------|
| streaming.py | 21 | 11 | 2 | 1 | ~930 | ~91% |
| augmented_mesh.py | 28 | 18 | 5 | 0 | ~1,250 | ~83% |
| boundary.py | 17 | 14 | 4 | 0 | ~670 | ~83% |
| **TOTALS** | **66** | **43** | **11** | **1** | **~2,850** | **~86%** |

**Interpretation:** Phase 2 rebalancing can retire ~14% of inventory docstrings (HISTORY + COMMENT-cut + TWIN) with high confidence. The 66 CONTRACT blocks are all essential for calling/modifying the operators; the 43 HISTORY blocks are campaign narrative suitable for git history / `.claude/plans/*roadmap.md` archive files.
