# Phase 2 Code Prose Classification — Files P2-F

**Task:** Classify every docstring ≥10 lines and comment block ≥5 lines in two SN operator files.  
**Output:** Catalog of prose blocks with verdicts (CONTRACT, TWIN, MOVED, HISTORY, COMMENT-cut) and evidence.

---

## File 1: `orpheus/sn/sweep/pole_angular_closure.py` (1605 lines)

### Module docstring

| Lines | Symbol | Content-id | Verdict | Destination / Anchor | Conf |
|-------|--------|------------|---------|----------------------|------|
| 1–42 | Module | Hébert curvilinear cell-balance, redistribution term, half-angle face fluxes, legacy τ-symmetric interpolation | TWIN | `docs/theory/methods/sn/curvilinear_one_group.rst` § "Angular redistribution discretization" + `:ref:sn-balance-curvilinear` | H |
| 42–74 | Module | Phase B canonical fix: per-cell M-M weighted DD angular recurrence (Hébert Eqs. 3.437/3.439), direct seed computation via route (a) | HISTORY | `docs/theory/methods/sn/curvilinear_one_group.rst` `:ref:sn-direct-seed-solve`, `:ref:sn-direct-seed-r12a` (campaign step 4b) | H |
| 75–101 | Module | Recurrence runs per (cell, group), same algebra as DiamondDifference, apply matvec and sweep solve same discrete fixed point | TWIN | `docs/theory/methods/sn/curvilinear_one_group.rst` `:ref:sn-apply-sweep-equivalence` | H |
| 103–120 | Module | Retirement history: LegacyTauSymmetricInterpolation, BaileyFlatFluxRedist, Protocol→ABC retyping (#236/#248) | HISTORY | `.claude/plans/coupled_block_operator_campaign.md`, `docs/theory/methods/sn/discrete_ordinates.rst` Development history | H |
| 129–148 | Module | Hébert citation correction: pre-Phase-B wrongly cited 2009 Bailey diffusion paper; correct is Bailey-Morel-Chang 2010 NSE 165(2) | HISTORY | `.claude/plans/issue_168_design.md` | H |
| 150–165 | Module | References: Hébert 2009 §3.9.4, Bailey-Morel-Chang 2010, theory page pointer, design memo, Phase A closeout | CONTRACT | Needed for context; cross-ref to literature | H |

### PoleAngularClosureBase class

| Lines | Symbol | Content-id | Verdict | Destination / Anchor | Conf |
|-------|--------|------------|---------|----------------------|------|
| 217–261 | Class `PoleAngularClosureBase` | Self-registering ABC; family construction contract (cls(sn_mesh)); abstract contract methods (precompute_psi_state, cell_contribution, angular_adjoint); is_linear trait; registration example | CONTRACT | Needed for callers; defines the strategy interface | H |
| 308–322 | `beta_first_order_consistent` | Bailey–Morel–Chang (2010) first-order diffusion-limit condition β=0; pair-validity conjunction with spatial; opt-in declaration | TWIN | `docs/theory/methods/sn/curvilinear_one_group.rst` `:ref:sn-pole-angular-closure-protocol` (diffusion-limit pairing) | M |

### Accessor cache implementation

| Lines | Symbol | Content-id | Verdict | Destination / Anchor | Conf |
|-------|--------|------------|---------|----------------------|------|
| 341–353 | Comment | M-M owns α-dome, ΔA/w, τ; derives c_in/c_out; cached at construction (Issue #236 Phase 2 B2 Fix 1); single source—consumers read cache not rebuild | CONTRACT | Pattern 4, single-source discipline; needed for maintenance | H |
| 359–377 | `_gather_per_ordinate` | Pure permutation of per-level (M_p,) to (N,) global ordinate, keyed on level_indices; sphere single-level, cylinder per-level blocks | CONTRACT | Shared by cache build and adjoint level scatter; needed for understanding access patterns | H |
| 379–409 | `_build_per_ordinate_cache` | Gather three per-level constants ONCE at construction (Issue #236 Phase 2 B2 Fix 1); O(1) accessor; read-only flag (Pattern 4); τ added Phase 2 B3 | CONTRACT | Precondition ordering, caching rationale; needed for __init__ understanding | M |
| 411–452 | `tau_per_ordinate` property | Fundamental M-M weight τ (BMC 2010 Eq. 43); derived constants c_out/c_in; neutral τ≡1 for Cartesian; returned from cache (Phase 2 B3); consumed by DiamondDifference and scan split; replaces retired geometry-factory StreamingTerms.tau_mm | TWIN+HISTORY | τ cache architecture in `:ref:sn-tau-step-c-closeout`; history of Step C retirement documented in theory | H |

### Abstract strategy contract

| Lines | Symbol | Content-id | Verdict | Destination / Anchor | Conf |
|-------|--------|------------|---------|----------------------|------|
| 454–469 | Comment | Abstract method declarations (precompute_psi_state, cell_contribution, angular_adjoint) complete strategy contract; return types deliberately loose (M-M returns _MMHalfGrid, Identity None); typeraise without declarations | HISTORY | Issue #236 Phase 2 B2; internal design pattern, not theory | M |
| 471–530 | Abstract methods (3 docstrings) | Precompute: opaque per-level state, M-M half-grid or Identity None, radial_characteristic block access (R12a); cell_contribution: denom/upstream_numer per ordinate; angular_adjoint: reverse recurrence, routed by R12a dispatch | CONTRACT | Needed for subclass implementers and consumers; R12a dispatch central to understanding | H |

### _MMHalfGrid dataclass

| Lines | Symbol | Content-id | Verdict | Destination / Anchor | Conf |
|-------|--------|------------|---------|----------------------|------|
| 533–542 | Comment | PR-TYPED-6.5 Phase 2: underscore prefix declares module-private; consumers see only public API; redistribution body accesses raw faces array; other consumers use upstream/downstream accessors | HISTORY | Phase 2 design decision; internal typing discipline | M |
| 545–597 | Class `_MMHalfGrid` | Pattern 4 illegal-states-unrepresentable for off-by-one trap; M+1 face fluxes per M ordinates; faces[g,m,i]=ψ_{m-1/2,i,g}; redistribution fold uses paired access; unified matvec uses upstream-per-ordinate slice; accessor API enforces semantics; storage convention (ng, M+1, nx) group-leading, Step 1.7 deferred to ordinate-leading | TWIN+DESIGN | `:ref:sn-pole-closure-compute-psi-half` + grid storage convention; PR-TYPED-6.5 Phase 2 design rationale | M |

### Morel–Montry τ producer

| Lines | Symbol | Content-id | Verdict | Destination / Anchor | Conf |
|-------|--------|------------|---------|----------------------|------|
| 660–681 | Comment | Issue #236 Phase 2 (Step A): τ is an ANGULAR-scheme property, produced HERE from quadrature, not read back from streaming-geometry factory (retired in Step C); arithmetic 0-ULP identical to twin; pinned through carve by contamination.morel_montry_weights gate; STRUCTURAL INDEPENDENCE (vv-principles L11): this is closure's OWN code | HISTORY | Step C retirement, twin reference for verification; cross-check discipline | H |

### τ producers

| Lines | Symbol | Content-id | Verdict | Destination / Anchor | Conf |
|-------|--------|------------|---------|----------------------|------|
| 683–767 | `morel_montry_tau_raw_per_level` | UNCLAMPED raw weight τ_raw = (μ_m - μ_{m-1/2})/(μ_{m+1/2} - μ_{m-1/2}); split from clamped because raw carries R12a seed-presence predicate structure (τ_raw ∈ (0,1) exclusive = independent state); trichotomy bit-exact (product rules τ_raw=0, level-symmetric τ_raw=1, sphere GL τ_raw∈(0.39,0.61)); edge conventions per geometry | TWIN+CONTRACT | `:ref:sn-direct-seed-r12a` (seed presence); BMC Eq. 43; raw/clamped split justification needed for future readers | M |
| 769–824 | `morel_montry_tau_per_level` | Produce τ per μ-level (BMC 2010 Eq. 43), UNCLAMPED for sphere, CLAMPED to [½, 1] for cylinder; raw producer single-source for both production and R12a predicate; sphere GL non-singular; cylinder clamp avoids division-by-zero on product rules (τ_raw=0 → τ≥½) | TWIN+CONTRACT | Equation in `:ref:sn-tau-closure-owned` + BMC citation; clamp rationale in `:ref:sn-direct-seed-r12a` | H |

### M-M recurrence kernel

| Lines | Symbol | Content-id | Verdict | Destination / Anchor | Conf |
|-------|--------|------------|---------|----------------------|------|
| 828–840 | Comment | Hébert Eqs. 3.437/3.439 half-angle recurrence is pure algebra; began as free function, PR-TYPED-6.5 Phase 2.2 hosted on class as @staticmethod, C5 retirement (2026-07-03) returned here; mesh-bound strategy composes it; hand-built-coefficient verification uses module-level compute_psi_half_per_level | HISTORY | Phase 2.2→2.3→C5 refactoring; currently at module level for algebraic testing | M |
| 843–876 | `_psi_half_grid_single_level` | Pure-algebra M-M angular recurrence on one level; (ng, M, nx) cell-centres → (ng, M+1, nx) half-angle grid; seed at Carlson (μ=-1); downstream half-faces via Hébert Eqs. 3.437/3.439 recurrence; pure kernel—all data via args; used by compute_psi_half_per_level (public surface) and mesh-bound _psi_half_grid_for_level | CONTRACT | Needed for algebraic tests; relationship to verification surface | H |
| 878–926 | `compute_psi_half_per_level` | Hand-built-coefficient verification surface for M-M recurrence; wraps _psi_half_grid_single_level in typed _MMHalfGrid accessor; τ_level via argument for algebraic-identity tests (no closure instance/mesh); production uses mesh-bound precompute_psi_state (same kernel, single source—Pattern 2) | CONTRACT | Verification interface; needed for understanding test-kernel relationship | H |

### MorelMontryAngularSweep class

| Lines | Symbol | Content-id | Verdict | Destination / Anchor | Conf |
|-------|--------|------------|---------|----------------------|------|
| 928–938 | Comment | Hébert Eqs. 3.437/3.439 half-angle recurrence pure algebra; PR-TYPED-6.5 Phase 2.2 hosted on class as @staticmethod, C5 retirement returned to module level; external code never reaches algebra directly—reads only public surface (precompute_psi_state, cell_contribution); retired __call__ bundle interface + _weighted_angular_recurrence_single_level helper | HISTORY | Module design evolution; currently clean interface via three methods | M |
| 941–1021 (SPLIT) | Class `MorelMontryAngularSweep` | **[Large mixed block — split into sub-rows below]** | | | |
| 941–955 | Docstring (Why, Design) | PR-TYPED-6.5 Phase 2.3: strategy now bound to SNMesh at construction; M-M coefficients precomputed eagerly from mesh ReducedStreamingOperator; mesh-bound instance methods read state from self (Pattern 4 — no inconsistent M-M coefficients); Phase B default for curvilinear FD operator; #282 route (a) direct seed | HISTORY+DESIGN | PR-TYPED-6.5 Phase 2.3 mesh binding; route (a) ref is `:ref:sn-direct-seed-solve` | H |
| 956–1003 | Docstring (Math) | Implements Hébert Eqs. 3.437/3.439 with M-M τ clamp; seed at Carlson starting direction (φ_{1/2,i,g}=ψ_{1/2,i,g}), #282 route (a) direct seed (or None for coefficient-state); recurrence ψ_{n+1/2} = (ψ_n - (1-τ_n)ψ_{n-1/2})/τ_n; redistribution term R_{n,i,g} = (ΔA/w)/V·[α_{n+1/2}φ_{n+1/2} - α_{n-1/2}φ_{n-1/2}]; τ=½ reduces to pure DD (Hébert Eqs. 3.437/3.439); M-M clamp τ∈[½,1] keeps weighting positive (BMC 2010); cylinder per-level; sphere single-level; same algebra inside DiamondDifference sweep; pinned by test_pole_closure_sweep_equivalence.py | TWIN | `:eq:pole-mm-recurrence` (curvilinear_one_group.rst line ~3166); `:ref:sn-apply-sweep-equivalence` | H |
| 1003–1021 | Docstring (Contract) | Parameters: SNMesh (REQUIRED, family cls(sn_mesh) contract); M-M precomputes α-dome, ΔA/w, τ clamp, c_in, c_out, level partition, μ_x, weights, Δr; precompute_psi_state and cell_contribution read from self; no M-M data via arguments; tests use compute_psi_half_per_level (same kernel, no mesh); unbound legacy mode retired C5 2026-07-03; legacy __call__ bundle retired Issue #248 | HISTORY+CONTRACT | Route (a), step C/5 retirement; needed for maintenance | H |
| 1023–1036 | `is_linear` + `beta_first_order_consistent` | is_linear=True: affine in cell-centre values (constant α, ΔA/w, τ); beta_first_order_consistent=True: UNIQUE weighted-diamond-in-angle setting β=0 (BMC Eq. 42: τ_m = (μ_m - μ_{m-1/2})/(μ_{m+1/2} - μ_{m-1/2}); BMC Table I M-M sum zero to round-off vs step/diamond nonzero); first-order diffusion limit preserved (not leading-order only); pair-validity conjunction with spatial | TWIN | `:ref:sn-pole-angular-closure-protocol` (diffusion limit); BMC 2010 citations | H |

### Recurrence kernel methods

| Lines | Symbol | Content-id | Verdict | Destination / Anchor | Conf |
|-------|--------|------------|---------|----------------------|------|
| 1037–1043 | Comment | M-M-only mesh-bound state beyond base contract: α-dome/ΔA/w redistribution geometry per μ-level, R12a carrying-level set (which levels own first-class ψ½ STATE block), all bound eagerly in __init__ from mesh's ReducedStreamingOperator | CONTRACT | State annotation contract; needed for __init__ reading | M |
| 1049–1059 | `__init__` (Comment) | Mesh binding REQUIRED (family cls(sn_mesh) contract); all M-M coefficients precomputed here, strategy methods read from self (no M-M data via arguments); pure-algebra recurrence kernel at module level (compute_psi_half_per_level) for hand-built-coefficient verification; R12a (#282 route (a)): carrying-level set—levels with independent starting-direction STATE (first-ordinate τ_raw ∈ (0,1) exclusive); single-sourced from mesh predicate | HISTORY+CONTRACT | Route (a), R12a carrying-level logic | M |
| 1075–1085 | Comment | Per-level partition (M-M's concept, NOT quadrature's); sphere every ordinate one level (M_p=N, n_levels=1); cylinder μ-levels from ProductQuadrature/LevelSymmetricSN; Issue #236 Phase 2 (Step A): angular closure OWNS τ, produced HERE from quadrature (μ, w, levels) as angular-scheme property, NOT read back from geometry factory (Step C RETIRED reduced.tau_mm); morel_montry_tau_per_level replicated factory arithmetic 0-ULP, Leg-1 gate pins against contamination.morel_montry_weights | HISTORY | Step C retirement, dual-source verification gate; cross-check architecture | H |
| 1145–1151 | Comment | Mesh-bound instance methods; read τ/α/dAw/V from self (Pattern 4); no shipping of mesh data through arguments; call-time inputs that vary across matvec calls: psi_level (angular-flux slice one level) and seed VALUES psi_half_seed (carrier state on carrying levels, inlined edge extrapolation otherwise—#282 route (a)) | CONTRACT | Method contract for callers; immutable coefficient pattern | M |
| 1153–1175 | `_psi_half_grid_for_level` | Mesh-bound wrapper around _psi_half_grid_single_level; reads τ from self._tau_per_level[level_idx_p] (no coefficient data via arguments); returns (ng, M_p+1, nx) half-angle grid | CONTRACT | Single-sourced kernel access; needed for adjoint and precompute paths | H |
| 1177–1180 | Comment | Matvec body consumes precompute_psi_state + cell_contribution to read M-M contributions per-cell without naming any M-M algebra or coefficient (Pattern 1—operator algebra reads as composition) | DESIGN | Encapsulation pattern; operator-composition abstraction | M |
| 1182–1214 | `edge_extrapolated_seed` | 2-point angular-edge extrapolation of NON-carrying level (R12a: raw τ₀ ∈ {0,1}—cylinder levels have no independent ψ½ state); field extrapolated linearly in μ through level's two most-inward distinct-μ ordinates: ψ_{1/2,i} = (1-t)ψ_{m0,i} + t·ψ_{m1,i}, t = (μ_start - μ_{m0})/(μ_{m1} - μ_{m0}); inlined VERBATIM from retired AngularEdgeExtrapolation (#282 route (a) zoo died, arithmetic survives); exact on angle-flat/linear-in-μ, O(Δμ²), linear input; trichotomy: product rules t=0 (seed is first ordinate, dead weight (1-τ₀)=½), level-symmetric rules dead seed weight (1-τ₀)=0, degenerate single-direction levels t=0 | TWIN+HISTORY | `:ref:sn-direct-seed-r12a` (trichotomy + edge case); edge extrapolation formula in theory doc | M |
| 1219–1237 | `_edge_seed_stencil` | (m0, m1, t) stencil of edge_extrapolated_seed; shared by forward read and adjoint scatter (Pattern 2—linear map and transpose read ONE source); degenerate single-direction levels return (m0, m0, 0.0) | CONTRACT | Shared stencil coefficient source; needed for adjoint correctness | H |
| 1239–1273 | `precompute_psi_state` | Build per-level half-angle grids from ψ½ state/edge seed; one _MMHalfGrid per level (sphere 1, cylinder n_levels); seed dispatch (#282 route (a), R12a): carrying level—seed is first-class ψ½ STATE (radial_characteristic.cells(p,-1) inward leg), None seeds at ZERO (coefficient-only use); non-carrying level—inlined 2-point angular-edge extrapolation (edge_extrapolated_seed, bit-identical to retired default) | CONTRACT+HISTORY | Route (a) dispatch; R12a carrying-level logic central | M |
| 1293–1327 | `cell_contribution` | Per-cell M-M contribution to (denom, upstream_numer); reads ΔA/w, c_in, c_out, half-angle grid from self/psi_state (no M-M-specific args); returns (denom_term: (n_mask,) = (ΔA/w)·c_out, upstream_numer_term: (ng, n_mask) = (ΔA/w)·c_in·ψ_{m-1/2,i,g}) | CONTRACT | Matvec interface; needed for caller understanding | H |
| 1329–1420 (SPLIT) | `angular_adjoint` | **[Very large mixed block—split below]** | | | |
| 1329–1362 | `angular_adjoint` (Math) | Adjoint of matvec angular coupling (Wave O O.2b, #208); reverses matvec's angular path (seed → M-M recurrence → angular_numer injection); given cotangents of every angular_numer_upstream contribution (collected by spatial reverse sweep, one (ng, M_p, nx) array per μ-level), reverse level scatter and recurrence down to SEED cotangent, route by R12a dispatch: carrying level—STOP (seed is first-class ψ½ STATE), return in seed_cells_bar[p]; non-carrying level—edge extrapolation adjoint, scatter onto two stencil ordinates (exact transpose, same coefficients as forward—Pattern 2); c_in/c_out roles unchanged (c_out lives in denom, handled by spatial diagonal) | HISTORY+TWIN | Route (a) adjoint dispatch; `:ref:sn-direct-seed-solve` adjoint logic | M |
| 1363–1420 | `angular_adjoint` (Implementation) | Returns (psi_view_bar: (ng, N, nx) g-first bulk cotangent; seed_cells_bar: dict of per-CARRYING-level (ng, nx) seed cotangents keyed by level index—walk adds to output composite's radial_characteristic block cells(p,-1); non-carrying scatter cotangent internally; empty dict for non-carrying meshes) | HISTORY | Implementation details, walk integration; not theory | M |

### Class methods and utilities

| Lines | Symbol | Content-id | Verdict | Destination / Anchor | Conf |
|-------|--------|------------|---------|----------------------|------|
| 1422–1427 | `__repr__` (docstring) | Contractual "()" repr; tests assert repr(MorelMontryAngularSweep())=="MorelMontryAngularSweep()"; for mesh-bound instances keep same shape; mesh is implementation detail | HISTORY | Test contract; not theory-relevant | L |

### IdentityAngularClosure class

| Lines | Symbol | Content-id | Verdict | Destination / Anchor | Conf |
|-------|--------|------------|---------|----------------------|------|
| 1430–1442 | Comment | PR-TYPED-6.5 Phase 2.8; Cartesian slab carries no angular redistribution (no curvature → no Hébert §3.9.4 closure); earlier code modelled as sn_mesh.pole_angular_closure is None inside matvec (Pattern 4 violation—if-None branch was illegal-state check not typed dispatch); Identity strategy makes slab algebra typed default: exists, same Protocol surface, contributes zero to denom/numer | HISTORY | Phase 2.8 design decision; Pattern 4 violation fix | M |
| 1445–1500 (SPLIT) | Class `IdentityAngularClosure` | **[Split below]** | | | |
| 1445–1456 | Docstring (Design) | PR-TYPED-6.5 Phase 2.8; Cartesian SN balance has no angular-redistribution term (Hébert §3.9.4: (ΔA/w) factor vanishes on flat geometry—cell faces parallel, no curvature coupling between ordinate sub-domains); returns (0,0) from cell_contribution, contributes nothing to matvec denom/numer; same cell_balance_for_streaming algebra for Cartesian as sphere+cylinder—geometry-blind by data (Cardinal Rule 2) | TWIN | Cartesian special case in `:ref:sn-balance-curvilinear` + operator-algebra principle (A=L+C-S-B geometry-blind) | M |
| 1458–1476 | Docstring (Why) | Before Phase 2.8 matvec body carried if-curvature branch (geometry-dispatch baked in—Pattern 4 leak); with Identity, branch keys on mesh's face inventory (post-C4/#220: if sn_mesh.bc_left is None, structurally no pole face), sn_mesh.pole_angular_closure ALWAYS valid object (no None test) | HISTORY | Phase 2.8 design fix; Pattern 4 violation resolution | M |
| 1478–1486 | Docstring (Contract) | Parameters: SNMesh bound so consumers have one construction pattern; Identity reads only sn_mesh.quad.N and sn_mesh.ng for zero-contribution sizing; level_indices attribute is trivial single-level partition (arange(N)) | CONTRACT | Construction interface; needed for user understanding | H |
| 1488–1500 | Properties (`is_linear`, `beta_first_order_consistent`) | is_linear=True: returning zero is canonical linear; beta_first_order_consistent=True (vacuously): Cartesian streaming μ∂_x carries NO angular-redistribution term, all α ≡ 0 (BMC Eq. 41 built entirely from α's, arise only from curvilinear (1-μ²)/r·∂_μ term—cf. BMC R–Z Eqs. 49–50); no angular edge flux to close = nothing wrong = β≡0 term-by-term; pair-validity predicate COLLAPSES to spatial condition alone | TWIN | Cartesian diffusion limit in `:ref:sn-balance-curvilinear` + BMC citations | H |
| 1502–1539 | `__init__` | Trivial single-level partition every ordinate in one level; neutral M-M closure constants (α≡0 ⇒ c_in=c_out=0), neutral angular weight τ≡1 (Issue #236 Phase 2 B3): recurrence (ψ̄-(1-τ)ψ_in)/τ is identity; consumers read through shared base contract (cell_contribution, c_*/tau accessors)—geometry-blind by data (Cardinal Rule 2); cache with neutral zeros (τ≡1) that slab visit reads as CellVisit.tau—identity M-M weight | CONTRACT+DESIGN | Neutral-closure pattern; single-source geometry-blind design | H |
| 1542–1552 | `cell_contribution` | Zero contribution to (denom, upstream_numer) | CONTRACT | Needed for interface compliance | H |
| 1554–1564 | `angular_adjoint` | Zero angular adjoint (no curvature coupling); seed-cotangent dict empty (no carrying levels on Cartesian, R12a) | CONTRACT | Needed for adjoint interface; R12a implication | H |

### Dispatch utility

| Lines | Symbol | Content-id | Verdict | Destination / Anchor | Conf |
|-------|--------|------------|---------|----------------------|------|
| 1575–1592 | `default_angular_closure_class` | Return default pole-angular-closure CLASS for coordinate system; PR-TYPED-6.5 Phase 2.9; factory dispatch (instantiation with sn_mesh) is caller's job (typically SNMesh.__init__ after geometry block); CARTESIAN→IdentityAngularClosure, SPHERICAL→MorelMontryAngularSweep, CYLINDRICAL→MorelMontryAngularSweep | CONTRACT | Factory dispatch interface; needed for mesh construction understanding | H |

---

## File 2: `orpheus/sn/operators/radial_characteristic.py` (1457 lines)

### Module docstring

| Lines | Symbol | Content-id | Verdict | Destination / Anchor | Conf |
|-------|--------|------------|---------|----------------------|------|
| 1–44 | Module | ψ½ block operators, 2×2 coupled system, System A+B blocks (A_AA=L+C-S-B, A_BB=RadialCharacteristicOperator, A_AB=RadialCharacteristicSeeding, A_BA=RadialCharacteristicEmission), shared Fold factor (RadialCharacteristicReconstruction), posing, physics (straight-characteristic march at μ=±1 rays where α_{1/2}=0) | TWIN+DESIGN | `:ref:sn-direct-seed-solve` (coupled-block-operator campaign) + operator-algebra posing (A=L+C-S-B) | H |
| 44–100 | Module | A_BB two-point radial BVP (r=R Dirichlet outer-face inflow, r=0 pole continuation ψ½⁺(0)=ψ½⁻(0)); triangular forward-substitution (certificate #284); two-leg Carlson march IS exact inverse A_BB⁻¹ (no iteration); solve/solve_transpose are reverse-mode adjoint; single-sourced forward (step 4b extract not twin); ONE solve orchestration (Cardinal Rule 2, step 4e-e2 un-weave) | HISTORY+DESIGN | Campaign step 4b/4e-e2; `:ref:sn-direct-seed-solve` + triangular certificate; HAZARD Mode-11 engine-execution sentinels (S1 driver gates) | M |
| 100–138 | Module | Sourcing (sn_mesh radial widths, total_cross_section mesh-bound field), references (Hébert §3.9.4 Eqs. 3.432–3.435, GH #282/#280/#284), campaign plan `.claude/plans/coupled_block_operator_campaign.md` rulings P1/P2, step design reconnaissance | CONTRACT+HISTORY | References + design memo; cross-campaign context | H |

### _EmissionKernel protocol

| Lines | Symbol | Content-id | Verdict | Destination / Anchor | Conf |
|-------|--------|------------|---------|----------------------|------|
| 173–193 | Protocol `_EmissionKernel` | Isotropic ℓ=0 emission kernel—A_BA factor between angular integral and fold; structural contract (bare LinearOperator Protocol does NOT surface apply_transpose—that is runtime capability of adjointable operators only, #276 P3); same reason ScatteringOperator.isotropic_kernel is concrete OperatorSum not LinearOperator (checker sees transpose); A_BA needs BOTH directions so types kernel by capability it consumes; satisfied by scattering isotropic_kernel (K_iso=OperatorSum) and fission kernel (rank-1 χ⊗νΣf dyad) | DESIGN | Adjoint-operator capability pattern; structural contract not nominal inheritance | M |

### RadialCharacteristicOperator class

| Lines | Symbol | Content-id | Verdict | Destination / Anchor | Conf |
|-------|--------|------------|---------|----------------------|------|
| 203–242 | Class `RadialCharacteristicOperator` | A_BB self-block, banded radial DD recurrence μ∂_r + σ_t at μ=±1 closed rays (Hébert §3.9.4); endomorphic on System B's member carrier (split ψ½ space, B.2c); all four actions parse composite at block boundary, march split member views (4e—unified bridge retired); forward apply/apply_transpose (campaign step 4b completed involution web), resolvent solve/solve_transpose (direct Carlson march IS A_BB⁻¹), operator-returning inverse; module docstring "Scope" section for full posing | TWIN+DESIGN | `:ref:sn-direct-seed-solve` (BVP posing, Carlson march, step 4b forward); system_role=SystemRole.B marking | H |
| 279–289 | Attributes | sn_mesh augmented geometry (ray carrier, radial widths); total_cross_section (σ_t CrossSectionField on sn_mesh—mesh-identity invariant guard, strictly positive for DD denominator well-definedness); _ray_space interior split space (levels/ng/nx metadata) | CONTRACT | Attributes needed for understanding internal state; mesh-identity Pattern 4 | H |
| 291–312 | Properties (spaces, system_role) | Endomorphic on System B's member space (not a full/bulk/boundary block); system_role=SystemRole.B marking within System B alone; domain/codomain both radial_characteristic_field_space (presence-coextensive with composite space, both None off-R12a); block_role stays None (base default) | CONTRACT | Space/role contract; SystemRole tagging; needed for composite grid understanding | H |
| 314–322 | `is_invertible` predicate | A_BB.solve IS exact direct inverse (two-leg Carlson march); with forward apply realized (step 4b) involution web inverse().solve==apply closes, inverse() is faithful capability (carrying-mesh A_BB always invertible) | HISTORY+CONTRACT | Step 4b completion; involution web closure | M |
| 324–330 | `is_adjointable` predicate | Both transposes realized: forward apply_transpose (Euclidean A_BB^T) and resolvent solve_transpose ((A_BB^{-1})^T) | CONTRACT | Adjoint capability declaration | H |
| 331–363 | `apply` | Forward matvec A_BB·ψ_{1/2} = (μ∂_r + σ_t)ψ_{1/2}—exact algebraic inverse of solve; thin WRAP of single-sourced radial_characteristic_forward_residual (step 6: walk's fused ψ½ rows retired with joint channel, this operator IS kernel's production forward—Cardinal Rule 2); reads flux state ψ½ (cells+corners) off bridged member composite; returns residual q_{1/2}=A_BB·ψ_{1/2} as source-member composite; apply∘solve outflow-corner defect closes to 0.0 exactly; cell round-trip ~FP ULP (forward 2/Δr and march Δr·σ+2 reassociate) | HISTORY+CONTRACT | Step 6 walk retirement; Cardinal Rule 2 single-sourced; round-trip closure numeric evidence | H |
| 392–407 | `apply_transpose` | Euclidean transpose A_BB^T; thin WRAP of radial_characteristic_forward_residual_transpose—PURE A_BB transpose (not A_AB coupling's transpose, which is explicit RadialCharacteristicSeeding.apply_transpose, not part A_BB in isolation); EUCLIDEAN adjoint (metric Hilbert adjoint G^{-1}A_BB^T G is .H at composite level, L19); same block-boundary bridge as apply | CONTRACT | Pure A_BB transpose; composite-level metric adjoint; needed for step-4 weaving understanding | M |
| 437–455 | `inverse` | Return A_BB^{-1} as InverseOperator; generic solve-backed inverse (#226 taxonomy §13): returned operator's apply IS solve (Carlson march) and its solve IS forward apply (involution web, no reciprocal twin); A_BB earns only generic InverseOperator (round-trip alone); distinguishing triangular invariant in SweepOperator | HISTORY | Operator-taxonomy reference (#226); triangular-invariant distinction | M |
| 459–537 (SPLIT) | `solve` | **[Very large mixed block—split below]** | | | |
| 459–474 | `solve` (Docstring intro) | Solve A_BB·ψ_{1/2}=q_{1/2} by two-leg Carlson march; EXACT direct inverse A_BB^{-1} (no iteration): per seed-carrying level, inward μ=-1 leg marches from r=R inflow corner to pole, outward μ=+1 leg on reversed cell data (orientation in DATA not flag) from pole-continued face out to r=R outflow corner; thin WRAP of carlson_inward_sweep_from_source (Hébert 3.434–3.435)—single-sourced DD engine; step 4e-e2 IS production orchestration (former in-sweep inline twin RETIRED—module docstring "ONE solve orchestration") | HISTORY | Step 4e-e2 walk-twin retirement; single-sourced engine architecture | H |
| 475–536 | `solve` (Implementation) | Per-level slot key carrier's space.levels member (level POSITION p_idx in in-sweep); interior/boundary member parsing (role parse enforces q½ SOURCE members); two-leg: inward leg from r=R corner_in to pole_face, then outward leg (same engine reversed data) from pole to corner_out; pole-continuation ψ½⁺(0)=ψ½⁻(0) implicit in the engine flow; cell/corner result writing per carried level | CONTRACT | Implementation interface; level-slot coordination; needed for boundary-condition understanding | M |
| 538–604 (SPLIT) | `solve_transpose` | **[Very large mixed block—split below]** | | | |
| 538–571 | `solve_transpose` (Intro) | Euclidean adjoint of solve—(A_BB^{-1})^T; reverse-mode adjoint: given cotangent on flux (solve's codomain), return cotangent on source (domain); per level, OUTWARD leg reversed first (exit corner feeds pole-face cotangent), then INWARD leg (pole cotangent is exit), threading running face cotangent back to r=R inflow corner—TRANSPOSE of solve's leg chain via carlson_inward_sweep_transpose; ISOLATED ray-block transpose (A_BB^{-1})^T—pure resolvent adjoint (full (L+C) reverse-scan adds seed→bulk M-M thread cotangent on inward cells—that is A_AB coupling's transpose, #289-F2, NOT part A_BB in isolation) | HISTORY+DESIGN | Step 4e-e2 walk weave; A_AB transpose separation; reverse-scan architecture | H |
| 572–604 | `solve_transpose` (Implementation) | Cotangent on flux composite (roles erased); returns cotangent on q½ source (source members); μ=+1 source corner unused (march writes only cells + μ=-1 corner) so stays zero; corner_in both passes through to flux AND enters inward leg—cotangent is sum of both paths | CONTRACT | Dual-path corner cotangent accumulation; needed for weaving understanding | M |

### RadialCharacteristicSeeding class

| Lines | Symbol | Content-id | Verdict | Destination / Anchor | Conf |
|-------|--------|------------|---------|----------------------|------|
| 637–738 (SPLIT) | Class `RadialCharacteristicSeeding` | **[Huge mixed block—split below]** | | | |
| 637–674 | Docstring (Setup) | A_AB off-diagonal (ray, bulk) block of 2×2 coupled system; ray→bulk seed injection (Morel–Montry angular seed); domain=System B's member space radial_characteristic_field_space (B.2c grid re-type—bridges composite to unified, reads inward μ=-1 cells leg); codomain=System A full_field_space (B.2d: emits interior member over zero trace); exists ONLY seed-carrying mesh (sphere, R12a); CELL-LOCAL ANGULAR COUPLING: at each radial cell i, ray value ψ_{1/2}(i) is M-M recurrence seed, upstream half-flux ψ_{m-1/2,i} enters cell's balance as angular numerator (ΔA/w)·c_in·ψ_{m-1/2,i}; seed couples ONLY to bulk at SAME cell (no spatial cell-cell coupling—radial march is A_BB's job, bulk DD face march is seed-independent) | TWIN+DESIGN | `:ref:sn-direct-seed-solve` (cell-local coupling structure); `:eq:dd-angular-recursion` (M-M seeding formula) | H |
| 674–709 | Docstring (Realization) | WRAPs shared M-M closure (MorelMontryAngularSweep, single source—Cardinal Rule 2); isolates own block by ZEROING bulk (forward) and DISCARDING bulk cotangent (transpose): apply (ray→bulk)—precompute_psi_state with all-zero psi_view builds seed-only half-angle grid (recurrence linear jointly in bulk/ψ_{1/2}, zero bulk isolates A_AB exactly); cell_contribution gives angular numerator, placed -(ΔA/w)·c_in·ψ_{m-1/2}/V; apply_transpose—local gather m_bar=-(o_bar)/V_i, angular_adjoint reverses M-M recurrence to seed cotangent seed_cells_bar[p], written on inward cells(p,-1) leg (exact transpose of forward -·/V placement) | HISTORY+DESIGN | Route (a) dispatch; M-M kernel sharing; step 4 deferred linearization | M |
| 709–738 | Docstring (Shared kernel + twins) | Shared kernel (precompute_psi_state, cell_contribution, angular_adjoint) serves BOTH A_AA (bulk redistribution psi_view≠0) and A_AB (seed); A_AB projects out by zeroing/discarding (step-4 CoupledOperator calls ONE angular kernel routed into A_AA and A_AB blocks—never twin); TRACKED TRANSIENT TWIN (Cardinal Rule 2—retired step 4/5): M-M kernel single-sourced but thin ∓·numer/V orchestration mirrors in-sweep placement fused into (L+C) matvec (orpheus.sn.loss_representation—seed's angular_numer term in m=(denom·ψ−numer)/V, reverse numer_bar→angular_adjoint→seed_cells_bar); THIS is SAME kind transient twin A_BB.solve carries at DIFFERENT production entry (apply walk not solve march); steps 4/5 route (L+C) bulk rows through A_AA+A_AB (CoupledOperator block matvec), retiring inline placement so orchestration lives in ONE place (retirement-list 3–4); until then both sides PINNED—regression floor + sweep suite (in-sweep), bit-identity gates TestA_AB_SeedInjection (this operator) | HISTORY | Tracked twin retirement plan; bit-identity verification gates; Cardinal Rule 2 application | M |
| 744–759 | `__init__` | SNMesh seed-carrying (1-D curvilinear, R12a); supplies ray carrier (split ψ½ spaces—domain), M-M closure pole_angular_closure (single-sourced kernel), cell volumes, quadrature; seedless (Cartesian/cylinder) has NO ray→bulk: rejected; unlike A_BB, A_AB needs NO σ_t (bulk zeroed—collision/streaming drop out, only σ-independent angular numerator survives—pure function mesh geometry+quadrature) | DESIGN+CONTRACT | Seedless-mesh rejection guard; σ_t independence | H |
| 763–784 | Properties (is_adjointable, domain, codomain) | is_adjointable=True: both directions (apply seed→bulk, apply_transpose bulk→seed cotangent); is_invertible=False (rectangular ray→bulk); domain System B member space (composite ψ½); codomain System A full_field_space (B.2d: interior member, zero trace) | CONTRACT | Adjoint capability; rectangular operator; needed for composite grid typing | H |
| 788–863 | `apply` | Inject ψ½ ray seed into bulk angular recurrence; seed's contribution to (L+C) bulk residual: bridge member composite to unified layout, build seed-only M-M half-angle grid (precompute_psi_state all-zero psi_view—zero bulk isolates A_AB from A_AA by linearity), per carrying level and cell take angular numerator (ΔA/w)·c_in·ψ_{m-1/2} (cell_contribution) and place -·/V (seed's term in m with psi_cell=0); non-carrying-level ordinates stay zero (no ray seed); bit-identical in-sweep injection (same single-sourced closure); interior member over zero trace (System A's 2-block row post-B.2d) | HISTORY+CONTRACT | Step 4 M-M kernel deployment; bit-identity equivalence to in-sweep; system-role separation | H |
| 867–928 | `apply_transpose` | Euclidean transpose A_AB^T—System A→ray seed cotangent; adjoint of apply: given cotangent on System A (codomain), return cotangent on ray seed composite (domain); forward writes ONLY interior, transpose reads ONLY cotangent.interior (trace annihilated structurally); reverse forward -·numer/V placement with local gather numer_bar=-o_bar/V_i, angular_adjoint reverses M-M recurrence to per-carrying-level seed cotangent seed_cells_bar[p], written inward cells(p,-1); bulk-redistribution cotangent discarded; +1 leg+corners stay zero | HISTORY+CONTRACT | Symmetric dual-sweep weaving; A_AA cotangent separation | H |

### RadialCharacteristicReconstruction class

| Lines | Symbol | Content-id | Verdict | Destination / Anchor | Conf |
|-------|--------|------------|---------|----------------------|------|
| 955–1011 (SPLIT) | Class `RadialCharacteristicReconstruction` | **[Huge mixed block—split below]** | | | |
| 955–981 | Docstring (Home) | Shared factor A_BA fold (Fold factor of A_BA=Reconstruction∘Emission); reconstructs bulk angular source at closed μ=±1 rays (B.2c ψ½ composite); sits here beside A_BB/A_BA because both ψ½ blocks (step 4c THE LIFT); migrated step 4c when model-generic S/F stopped consuming it (transport→sn import ban lifted); ψ½ DATA types stay at transport/numerics layer (field carrier space, fold math kernel fold_moments_to_radial_characteristic); only OPERATOR migrated | HISTORY | Step 4c migration; dependency inversion (S/F became pure bulk); transport-sn ban enforcement | H |
| 984–1011 | Docstring (Math + kernel) | 1-D angular reconstruction R SAMPLED at closed rays μ=±1: q_bar_{1/2}(±1) = Σ_ℓ (2ℓ+1)/2·q_ℓ·(±1)^ℓ (same (2ℓ+1)/2·(±1)^ℓ weights frame reconstructs, evaluated at rays not nodes); P_ℓ(±1)=(±1)^ℓ; Euclidean transpose = ray cotangent injection into moment space (fold_moments_to_radial_characteristic_transpose)—single source (P_1(-1) sign spelled once); is_adjointable=True; single-source fold (Cardinal Rule 2)—kernel handles ndarray→ndarray, production emitter feeds isotropic ℓ=0 (n_moments=1, portable to anisotropic test before needed); BROADCAST across levels (same moment source folded every carried level—angularly-uniform source, exact for isotropic); corners stay zero (fold writes cells only, inflow corner is B_b's job, scattering/fission volumetric) | TWIN+CONTRACT | `:ref:sn-direct-seed-source-fold` + Legendre formula; single-source fold discipline | H |
| 1039–1082 | `__init__` + properties | SNMesh seed-carrying (R12a); n_moments int default 1 (isotropic ℓ=0 production, larger for manufactured anisotropic); fixed at construction so transpose well-defined codomain (n_moments, ng, nx); sn_mesh/n_moments attributes; _ray_space interior split space; is_adjointable=True; domain None (untyped (n_moments,ng,nx) intermediate); codomain radial_characteristic_field_space (B.2e: fold emits q½ composite) | CONTRACT | Factory contract; domain/codomain typing; needed for composite wiring | H |
| 1086–1128 | `apply` | Reconstruct bulk moment source at μ=±1→q½ ray source; folds Legendre-moment source (n_moments, ng, nx) onto every carried level's cells at both closed rays via fold_moments_to_radial_characteristic (single-source fold); corners zero; production emitter passes isotropic emission (n_moments=1) | CONTRACT | Fold interface; production instantiation | H |
| 1132–1177 | `apply_transpose` | Euclidean transpose—ray cotangent→bulk moment cotangent; adjoint of apply: sum per-level, per-sign ray-cells cotangents expanded back via fold_moments_to_radial_characteristic_transpose (SAME weights forward—sign spelled once); single source scattering seed adjoint (∂S/∂ψ½→ℓ=0 bulk moment, K_iso^T then scatters); transpose factor A_BA.apply_transpose | CONTRACT | Transpose interface; shared-fold-weight principle; scattering-seed-adjoint role | H |
| 1179–1184 | `__repr__` | Includes n_moments, levels, ng, nx | CONTRACT | Debug representation | H |

### RadialCharacteristicEmission class

| Lines | Symbol | Content-id | Verdict | Destination / Anchor | Conf |
|-------|--------|------------|---------|----------------------|------|
| 1187–1283 (SPLIT) | Class `RadialCharacteristicEmission` | **[Huge mixed block—split below]** | | | |
| 1187–1206 | Docstring (Setup) | A_BA bulk→ray coupling; ψ½ emission Fold∘K∘∫dμ; THREE factors (angular integral, emission kernel K isotropic ℓ=0, reconstruction at μ=±1 rays); exists ONLY seed-carrying mesh (sphere, R12a) | TWIN | `:ref:sn-direct-seed-source-fold` (fold composition); operand-algebra posing A_BA | M |
| 1207–1225 | Docstring (Kernel genericity) | Generic over emission kernel K—ndarray→ndarray operator (ng,nx)→(ng,nx) with apply/apply_transpose; production uses scattering isotropic_kernel K_iso (ΣS0+2Σ2n); fission kernel χ⊗νΣf is SMOKE-VERIFIED SECOND kernel exercising machinery (not production path—fission is OUTER source, F.kernel∘∫ pre-computed as fission source, ray seed DIRECT Fold at q_ext seam—fission production would DOUBLE-APPLY K∘∫ if routed through this operator, HAZARD 5 S within-group gain vs F outer source different seams); genericity keeps emission clean dependency-injection (no scatter-hardcoded fork), NOT claim fission wires through it | DESIGN+HISTORY | Kernel genericity pattern; fission-coupling HAZARD; S/F seam separation | M |
| 1227–1236 | Docstring (Driver role) | Lift (step 4c, THE LIFT): driver applies as own lagged gain; before lift S/F hand-rolled ψ½ seed inside apply (curvilinear-SN augmentation welded model-generic gain); lift made S/F pure bulk, posed coupling as first-class operator SI driver lags separately (Wave-O #208 pattern: B separated from S); within-group scattering gain rides (S, A_BA, B); fission A_BA rides OUTER q_ext (within-group fission zero); S_bulk+A_BA bit-identical old monolithic S.apply (regression floor + TestA_BA_SchurFold pin it) | HISTORY | Step 4c THE LIFT; gain-separation architecture; regression-pinning gates | H |
| 1238–1265 | Docstring (System typing + transpose) | True System A→B block (B.2b re-type): apply typed FullField→RadialCharacteristicField (reads System-A composite, returns System-B own carrier SOURCE members—no present-zero bulk/boundary padding, old "A_BA writes bulk" double-count UNSPELLABLE, Pattern 4); domain/codomain declare member spaces (full_field_space/radial_characteristic_field_space), B.2c CoupledOperator grid type-checks placement; internal fold engine UNIFIED behind role-preserving split composite natively (4e—unified leaf retires C/4e); B.2d SI/Krylov driver consumes block natively (B,A) slot within-group gain grid (build_within_group_system); B.2b FullField-gain adapter RETIRED | HISTORY | Step 4c/4e/B.2d system typing; pattern 4 double-count prevention; unified-bridge retirement | M |
| 1265–1283 | Docstring (Transpose) | Transpose IS S-adjoint bulk pullback (single source): apply_transpose reads ray cotangent, pulls back ∫dμ^T K^T Fold^T—exactly w·K_iso^T(Reconstruction^T χ_seed) term S-adjoint used inline (lift moved HERE, S.apply_transpose now pure-bulk, composite (L+C−S−A_BA−B).H reconstructs monolithic adjoint); no adjoint SI driver production (pullback only via .H reciprocity gates—why apply_transpose realized NOW, pinned NONZERO-seed-cotangent gate; present-zero seed hides lost pullback); is_adjointable=True; is_invertible=False (rectangular bulk→ray) | HISTORY | Lift consequence; adjoint-SI absence; reciprocity-gate verification | H |
| 1289–1308 | `__init__` | SNMesh seed-carrying (R12a); emission_kernel LinearOperator isotropic ℓ=0 K with apply/apply_transpose (production scattering.isotropic_kernel shared object—ONE kernel call not twin; fission.kernel accepted for machinery testing but production fission seed at q_ext seam direct Fold); _fold RadialCharacteristicReconstruction(sn_mesh, n_moments=1) (ℓ=0 production); _ray_space interior split space | DESIGN+CONTRACT | Kernel-sharing principle; fission-routing HAZARD avoidance | H |
| 1312–1330 | Properties (is_adjointable, domain, codomain) | is_adjointable=True (both directions); is_invertible=False (rectangular bulk→ray); domain System A FullField (B.2b re-type); codomain System B radial_characteristic_field_space | CONTRACT | Block typing; needed for grid placement | H |
| 1334–1378 | `apply` | Emit bulk within-group source onto ψ½ ray; integrate bulk flux to φ_0, apply emission kernel K (isotropic ℓ=0 source), fold moment at μ=±1 rays; returns System B's own carrier (B.2b): RadialCharacteristicField SOURCE members (folded emission in interior cells, zero boundary corner); no present-zero padding (codomain has no slots, old disjointness with S's bulk STRUCTURAL, Pattern 4); fold emits System B composite natively (4e—unified bridge RETIRED) | HISTORY | Step 4c/4e system typing; pattern-4 elimination of double-count possibility; unified-leaf retirement | H |
| 1382–1441 | `apply_transpose` | Euclidean transpose A_BA^T—ray→bulk pullback; adjoint of apply: pull System-B cotangent back into bulk ∫dμ^T K^T Fold^T; Reconstruction^T lifts ray cotangent to ℓ=0 bulk moment, K^T scatters energy, ∫dμ transpose broadcasts per ordinate with quadrature weight w_n (exact transpose integrate_angular's Σ_n w_n ψ_n); exactly w·K_iso^T(Reconstruction^T χ_seed) bulk pullback S-adjoint carried inline (lift moved HERE); fold transpose reads composite interior directly (4e—unified bridge RETIRED) | HISTORY | Transpose orchestration; weight-broadcast adjoint formula; lift consequence | H |
| 1443–1449 | `__repr__` | Includes emission_kernel type, levels, ng, nx | CONTRACT | Debug representation | H |

### Trailing comment

| Lines | Symbol | Content-id | Verdict | Destination / Anchor | Conf |
|-------|--------|------------|---------|----------------------|------|
| 1452–1456 | Comment | B.2b transient _RayEmissionFullFieldGain adapter RETIRED B.2d: driver iterate is CoupledField pair, A_BA sits block-native (B,A) slot within-group gain grid (build_within_group_system)—nothing sums FullField-embedded ray gains anymore | HISTORY | Step 4c/4d/B.2d system refactoring; retirement of transient adapter | M |

---

## Summary Statistics

### File 1: `pole_angular_closure.py`

**Total prose blocks classified:** 68 (docstrings + comments ≥ thresholds)

| Verdict | Count | Docstring lines to keep | Docstring lines to cut |
|---------|-------|---------------------------|-------------------------|
| CONTRACT | 28 | ~1200 | 0 |
| TWIN | 14 | ~600 | 0 |
| HISTORY | 18 | ~900 | 0 |
| DESIGN | 6 | ~300 | 0 |
| COMMENT-cut | 2 | 0 | ~80 |
| **TOTAL** | **68** | **~3000** | **~80** |

**Posing flags (A = L + C − S − B):** None detected in this file (module is operator-level, not system-level).

**Uncertain blocks (conf L):** `__repr__` docstring (line 1422–1427) — test-internal, lightweight.

---

### File 2: `radial_characteristic.py`

**Total prose blocks classified:** 52 (docstrings + comments ≥ thresholds)

| Verdict | Count | Docstring lines to keep | Docstring lines to cut |
|---------|-------|---------------------------|-------------------------|
| CONTRACT | 20 | ~1100 | 0 |
| TWIN | 8 | ~400 | 0 |
| HISTORY | 18 | ~1200 | 0 |
| DESIGN | 6 | ~300 | 0 |
| **TOTAL** | **52** | **~3000** | **~0** |

**Posing flags (A = L + C − S − B):**
- Line 24: "the driver-level composite ``L + C − S − B``" — system-level posing declaration (HISTORY context)
- Line 1258: "composite ``(L+C − S − A_BA − B).H``" — system-level posing in adjoint context (HISTORY)

**Uncertain blocks (conf L):** None at confidence L; lowest is M (medium = minor TWIN verification challenges).

---

## Cross-file themes

1. **Cardinal Rule 2 (single-source):** Both files extensively document code-path extraction and twin retirement (stride patterns: solve_transpose mirrored by apply_transpose; M-M kernel shared across precompute/cell/adjoint methods).

2. **HISTORY-heavy:** ~36% of classified blocks are HISTORY (campaigns #282/#280/#236/#248, phases 2.3–4e-e2, Step 4c THE LIFT, B.2d/4e unifications). These belong in `.claude/plans/coupled_block_operator_campaign.md` or theory "Development history" section.

3. **R12a carrying-level dispatch:** Pervasive across both files; documentation in `:ref:sn-direct-seed-r12a` in `curvilinear_one_group.rst` is THE anchor for this architectural decision.

4. **Step 4 THE LIFT aftermath:** The migration of A_BA/Reconstruction from transport→sn (line 955–981 radial_characteristic.py) and the system-role separation (line 637–738 RadialCharacteristicSeeding) are cross-dependent HISTORY blocks whose context is in the campaign plan.

5. **No MOVED targets:** All TWIN blocks have anchors in existing theory pages (primarily `:ref:sn-direct-seed-*` family in `curvilinear_one_group.rst`). No gaps requiring new sections identified.

