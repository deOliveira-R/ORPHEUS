# Phase 2 Code-Prose Classification (P2-E)

Inventory of docstrings ≥10 lines and comments ≥5 lines in spatial discretization schemes.

Landing theory chapters for TWIN adjudication:
- `docs/theory/foundations/discretization.rst` (cell balance, closure theory)
- `docs/theory/methods/sn/slab_one_group.rst` (DD closure in slab context)
- `docs/theory/methods/sn/cartesian_multid.rst` (2-D LD / UBLD)
- `docs/theory/foundations/path_integral.rst` §6 (Padé/step/DD table — positivity)

---

## File 1: `orpheus/transport/spatial/scheme.py` (1,356 lines)

### Module docstring (lines 1–164)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 1–27 | module (Why abstraction) | Round 1 Wave C strategy abstraction; per-cell closure algebra lifted out | HISTORY | Record in #157/#158 issue trace; the architectural motivation is LIVE in the codebase comment; no doc needed | M |
| 37–72 | module (What each dataclass holds) | UpstreamState/CellResult field docs (spatial_upstream, angular_upstream, cell_average_flux, outgoing_spatial_flux, outgoing_angular_state) | TWIN | `:class:`~orpheus.transport.spatial.scheme.CellVisit` in Sphinx is the definitive description; code docstring mirrors Sphinx | M |
| 74–110 | module (Geometry-as-data §) | Issue #236 Phase G Step 2.5 collapsed runtime branches; neutral curvature for slab; StreamingTerms carries the data | TWIN | Anchor: `discretization.rst` §"dimension-agnostic"; geometry-blind algebra is the core contribution | H |
| 112–138 | module (Where consumers will call) | Wave D sweep refactor; per-visit packet from SNMesh.dag_walk; strategy dispatching | HISTORY | The Wave D refactor is in-flight work; this is the prospective consumer design; link to #159 plan | M |
| 140–154 | module (References) | Lewis & Miller §5.3, Bailey et al. 2009, link to theory/methods/sn/index | CONTRACT | These are the canonical citations for closure theory; keep as module-level references | H |

### CellVisit class docstring (lines 187–306)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 187–214 | CellVisit class | One visit packet during sweep; sweep-direction-resolved; per-cell guidance | TWIN | Anchor: schema in `discretization.rst` line 34 (CellVisit type); field-by-field contract is all in docstring | H |
| 216–291 | CellVisit Attributes | Detailed field descriptions: cell_idx, streaming_terms, face_area_downstream, c_in, c_out, tau | CONTRACT | These are critical call-site docs; kept inline because they specify DATA CONTRACTS for each field | H |
| 293–303 | CellVisit Notes | Sweep-direction pre-resolution; spatial_upstream + face_area_downstream give unified view | TWIN | Anchor: `discretization.rst` §"Pose by the invariant" (the invariant enforcement on data) | H |

### UpstreamState class (lines 322–338)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 322–338 | UpstreamState class | Input state: spatial_upstream (face flux) + angular_upstream (half-flux for curvilinear) | CONTRACT | Call-site input shapes and conventions; needed to invoke update() | H |

### CellResult class (lines 351–377)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 351–377 | CellResult class | Output state: cell_average_flux, outgoing_spatial_flux (None for cylindrical degenerate), outgoing_angular_state (None for slab) | CONTRACT | Call-site output shapes and conventions; needed to interpret update() results | H |

### DiscretizationScheme Protocol (lines 386–510)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 386–510 | DiscretizationScheme Protocol | Seven class-level traits (is_linear, is_positivity_preserving, is_affine_scannable, transverse_coupling_is_facewise, spatial_basis_per_axis, diffusion_limit_consistent, supports_curvilinear) + method signatures | TWIN | Anchor: schema in `discretization.rst` line 34 (key_types: DiscretizationScheme); trait descriptions match closure spectrum | H |

### update() method (lines 520–566)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 520–566 | DiscretizationScheme.update | Protocol method; compute per-cell average + downstream states | CONTRACT | Defines what a strategy must implement; call signature + contract for any concrete strategy | H |

### residual() method (lines 568–643)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 568–643 | DiscretizationScheme.residual | Apply-direction operator residual; round-trip invariant; use case (matvec form); FULL docstring | TWIN | Anchor: `slab_one_group.rst` §"Reciprocity invariant" (line 243+); the adjoint/apply-direction structure | H |

### DiscretizationSchemeBase class (lines 663–831)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 663–675 | DiscretizationSchemeBase intro | Round 1 shipped Protocol; Round 2 shipped DiamondDifference extraction; self-registering via RegistryMixin | HISTORY | #157/#158 Wave C milestones; this is progress narration, not load-bearing doc | L |
| 702–730 | Generic advection-reaction interface | Σ-stateless design; reaction_xs parameter instead of tied cross-section; scheme instance carries NO state | MOVED | Teaching content: "diffusion-readiness contract" — move to discretization.rst § "Advection-reaction reading" (extend the CFD subsection) | M |
| 742–831 | Reconstruction ops + Bit-identity section | Generic affine algebra (a, inverse_denom, w); source_emission/cell_average/outgoing_face_from_average staticmethods; why forms are universal; byte-identity vs principled-equivalence | TWIN | Anchor: `slab_one_group.rst` § "Generic affine outflow reconstruction" (lines 393–441); the bit-identity discussion mirrors code comments | H |

### Attribute docstrings (lines 835–979)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 837–853 | is_affine_scannable ClassVar | Opt-out default False; affine-scannable scheme supplies affine_scan_coefficients; generic reconstruction staticmethods | TWIN | Anchor: `slab_one_group.rst` § "Cartesian 1D: cumprod recurrence" (the affine recurrence ψ_out = a·ψ_in + b) | H |
| 855–909 | transverse_coupling_is_facewise ClassVar | 1-D vs cross-axis trait; facewise (0th order) vs slope-wise (1st order); contrast with is_affine_scannable; DD/Step True, LD False | MOVED | Teaching content belongs in discretization.rst § "Multi-dimensional extensions" (or new subsection on tensor-product structure) | M |
| 911–950 | spatial_basis_per_axis ClassVar | Per-axis moment count; DD/Step = 1, LD = 2; per-cell unknown = per_axis^d; per-face transverse = per_axis^{d-1}; gate on multi-moment wiring | MOVED | Move to discretization.rst or new `_ubld` theory page; the tensor-product Kronecker structure is mathematical, not just a code detail | M |
| 952–967 | diffusion_limit_consistent ClassVar | Thick-diffusion limit consistency; DD (LMM-1987 Eq 4.24) vs Step (no limit); PAIR validity with angular closure; ⚠ spatial vs angular clarity | TWIN | Anchor: `discretization.rst` § "The closure schemes" (lines 418–450) discusses diffusion limits implicitly; explicit statement needed in theory | H |
| 969–979 | supports_curvilinear ClassVar | Opt-in False default; curvilinear closure vs slab-only; rejection at strategy selection, not mid-sweep | CONTRACT | Selection logic; call-site guards depend on this value | H |
| 981–999 | has_transpose_kernel ClassVar | Opt-in False; reverse-mode VJP; DD True (hand-transpose), LD False (deferral to #280) | MOVED | Belongs in a new section on adjoint/reverse-mode operations (when that theory page lands) | M |
| 1002–1019 | is_multi_moment property | True iff per_axis > 1; LD (DG-P1) → True; DD/Step → False; source-of-truth for moment-buffer checks | TWIN | Anchor: `spatial_basis_per_axis` property definition above; SAME as per_axis > 1 | H |

---

## File 2: `orpheus/transport/spatial/diamond.py` (785 lines)

### Module docstring (lines 1–76)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 1–25 | diamond.py intro | Issue #196 Phase G Step 2.5; geometry-polymorphic via data (not branches); ONE body for slab/sphere/cylinder | HISTORY | The architectural milestone (Issue #196 Phase G Step 2.5) is the campaign record; link to it | M |
| 26–51 | The unified body | Three structural observations: cell-balance algebra (one formula, neutral curvature for slab), spatial closure same formula, angular closure same formula | TWIN | Anchor: `discretization.rst` § "Dimension-agnostic by construction" (lines 404–413); the unification principle | H |
| 52–76 | References | Hébert 2009 §3.9.4 (curvilinear SN), Lewis & Miller §4.5/§5.3 (DD/weighted-DD/Step/LD), Bailey–Morel–Chang 2010 (asymptotic diffusion limit + M-M weight clamp) | TWIN | Module-level citations; these anchor the theory pages | H |

### DiamondDifference class attributes (lines 104–193)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 104–119 | DiamondDifference class intro | Single body, geometry-polymorphic via data; no internal geometry dispatch | TWIN | Anchor: `discretization.rst` § "The closure schemes" (DD row, line 440) | M |
| 121–125 | is_linear | DD linear in source + upstream_state; Lewis & Miller §5.3 | TWIN | Anchor: `discretization.rst` closure table (line 441: "central / Keller box") | H |
| 127–134 | is_positivity_preserving | DD NOT positivity-preserving; can produce negative flux in thin/large-source cells | TWIN | Anchor: `discretization.rst` § "The closure schemes" (line 444: "can go negative") + Lewis & Miller §5.3 counter-example | H |
| 136–150 | is_affine_scannable | DD affine recurrence ψ_out = a·ψ_in + b; CumprodScan/ScanMarch via affine_scan_coefficients + generic base staticmethods; w = ½ | TWIN | Anchor: `slab_one_group.rst` § "Cartesian 1D: cumprod recurrence" (Eq. 340–353) | H |
| 152–164 | transverse_coupling_is_facewise | DD transverse coupling facewise (0th order); tensor-product separable; ScanMarch admitted; contrast with LD bilinear | TWIN | Anchor: `discretization.rst` § on tensor-product structure (if exists; else MOVED) | M |
| 166–179 | diffusion_limit_consistent | DD thick-diffusion limit IS consistent; LMM-1987 Eq. (4.24); SPATIAL axis only, not angular β-failure | TWIN | Anchor: `discretization.rst` closure table + LMM-1987 citation | H |
| 181–185 | supports_curvilinear | DD curvilinear closure exists; True; Morel-Montry angular redistribution | CONTRACT | Selection logic guard | H |
| 187–193 | has_transpose_kernel | DD 1-D reverse walk hand-transposes diamond chain; Wave O / O.2b; reciprocity pinned by test_g_adjoint_reciprocity | HISTORY | Wave O milestone; campaign record | L |

### update() method (lines 195–262)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 202–206 | update() docstring | One body, no geometry dispatch; see module docstring; geometry is data | CONTRACT | Call contract; minimal because module docstring covers the structure | M |
| 207–220 | Comment block (cell-balance solve) | Issue #236 Phase 2 B3: M-M constants c_in/c_out angular-closure-owned; cell_balance_terms unified helper | HISTORY | #236 Phase 2 B3 milestone; the architectural ownership rule | L |
| 223–241 | Comment block (spatial closure) | Outputs None for cylindrical degenerate (face_area_downstream == 0.0); slab/non-degenerate share formula ψ^s_out = 2·ψ_avg − ψ^s_in | MOVED | Teaching: explain WHY face_area_downstream == 0.0 is the degenerate signal (not numerical threshold); add to discretization.rst § "Geometry-as-data" | M |
| 238–251 | Comment block (angular closure) | Issue #236 Phase 2 B3: τ angular-closure-owned, sourced off visit; M-M formula; bit-identical to old geometry-side τ_mm (Leg-1 gate) | HISTORY | #236 Phase 2 B3 provenance claim; the ownership rule | L |

### residual() method (lines 266–358)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 274–299 | residual() docstring | Apply-direction companion; at converged cell_avg, residual = 0 to FP; round-trip contract; Issue #197 PR-TYPED-6 delegates to cell_balance_for_streaming | TWIN | Anchor: `slab_one_group.rst` § "Reciprocity invariant" (line 243+) + the "Reciprocity gating" testing narrative (line 271+) | H |
| 300–329 | Comment block (M-M pattern-2 migration) | Issue #236 Phase 2 B2: constants NO LONGER rebuilt inline; come from visit (via closure's accessors); Pattern 2 — single algebra source | HISTORY | #236 Phase 2 B2 refactoring milestone; the ownership migration | L |

### Shared DD Cartesian helpers (lines 360–410)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 360–368 | Comment block (2 helpers, Pattern 2) | Streaming-diagonal fold + w=½ reflection; ONE place diamond 2 enters Cartesian; Cardinal Rule 2 / Pattern 2 | COMMENT-cut | Narrates the refactoring; constraint is "2 = 1/w_DD folds once" (keep as comment near implementation) | H |
| 375–404 | _cartesian_streaming_diagonal docstring | DD's ×V streaming diagonal S = Σ_t + Σ_a 2 g_a + per-axis couplings; generic in wave-speed/reaction-rate; model-agnostic advection-reaction interface | MOVED | Teaching: diffusion readiness + CFD perspective belong in discretization.rst § "Advection-reaction reading" | M |
| 412–428 | _reflection_coeffs docstring | Apply-direction reflection coefficients (α, β) from the affine form ψ_out = β + α·ψ_in; shared arithmetic, byte-identical at w=½ | CONTRACT | Method signature + output contract | H |

### cell_kernel_batch (lines 432–500)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 432–489 | cell_kernel_batch docstring | Pure batched WDD cell update; dimension-generic; SINGLE source of DD cell math (Cardinal Rule 2 / Pattern 2); operation-order discipline for bit-identity | TWIN | Anchor: `slab_one_group.rst` § "Cartesian 1D" + the closure paragraph; the bit-identity pin (vv-principles Bit-identity vs principled-equivalence) is already cited | H |

### residual_kernel_batch (lines 502–543)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 502–543 | residual_kernel_batch docstring | Pure batched DD operator residual; apply-direction companion of cell_kernel_batch; same axis-convention + operation-order discipline | TWIN | Anchor: Reciprocity section + matvec discussion in `slab_one_group.rst` | H |

### affine_scan_coefficients (lines 547–678)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 547–678 | affine_scan_coefficients docstring | DD Σ_t-epoch scan coefficients (a_attenuation, inverse_denom, face_blend_weight = ½); Lewis & Miller §5.3; Hébert §3.9.4; curvilinear formulas; bit-identity discipline; why SEPARATE diagonal from Cartesian (#242) | TWIN | Anchor: `slab_one_group.rst` § "Cartesian 1D: cumprod recurrence" (lines 336–354) for the recurrence coefficients; the curvilinear detail mirrors Hébert | H |

### cartesian_scan_coefficients (lines 680–736)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 680–736 | cartesian_scan_coefficients docstring | 2-D/3-D row-march scan coefficients; DD applies 2 HERE; per-axis breakdown; operation-order discipline | MOVED | The 2-D scan-march extension belongs in cartesian_multid.rst or a new dedicated ScanMarch theory page | M |

### reflect_scan_coefficients (lines 738–753)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 738–753 | reflect_scan_coefficients docstring | DD apply-direction reflection scan coefficients (α = −1, β = 2ψ̄); delegates to _reflection_coeffs; byte-identical to legacy | CONTRACT | Method signature + output contract; the legacy equivalence is a verification gate | M |

### Comment block (storage adapters retired) (lines 755–781)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 755–781 | Storage adapter retirement | S6.4(e): update_batch/residual_batch moved to walk layer; now pure cell algebra in THREE groups (per-cell pair, batched kernel, scan coefficients) | HISTORY | Wave S6 milestone; the architectural layering (discretization vs walk storage) | L |

---

## File 3: `orpheus/transport/spatial/linear_discontinuous.py` (811 lines)

### Module docstring (lines 1–169)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 1–31 | LD module intro + diffusion-limit note | First higher-order O(h²) strategy; FULL LD (slope SOURCE moment Σ_s·φ̂ threaded) is diffusion-limit-consistent; Increment C #240 D5b-S3; test_ld_thick_diffusive_limit pinned; ERR-061 frame fix | HISTORY | #240 D5b-S3 milestone (slope source threading); Issue #247 (external slope source gap); the test pinning | L |
| 32–76 | Why LD carries two moments | LM-1989 linear in-cell function; two moments (average + slope); upwind discontinuous closure; 2×2 system + Schur complement; SymPy verification | TWIN | Anchor: `linear_discontinuous.rst` (theory page, if it exists); LM-1989 Eqs. (4.1a-c) + (4.3c) are the canonical references | H |
| 77–104 | Schur-complement scalar contract | Slope eliminated locally; 1-D affine recurrence ψ_out = a·ψ_in + b; fixed-source first; threading slope through source-iteration global iterate (deferred) | TWIN | Anchor: `linear_discontinuous.rst` theory page (if exists); the Schur form is the link to the affine scan coefficients | H |
| 106–144 | Scope + Traits | Slab/Cartesian only (curvilinear unpublished); is_linear/is_positivity_preserving/is_affine_scannable/spatial_basis_per_axis traits | TWIN | Anchor: `discretization.rst` closure table (LD row, line 445–449) | H |
| 146–158 | References | Larsen & Morel 1989 JCP (slab-LD + diffusion limits), LMM-1987 (LD has 4 diffusion limits, Step has 0), Lewis & Miller §5.3 (LD/positivity/negatives) | TWIN | Module citations; these anchor the theory pages | H |

### Docstring note block (lines 12–30, inside module docstring)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 12–30 | Module note: Diffusion-limit status | Scattering slope Σ_s·φ̂ IS threaded (Increment C); full LD recovers thick-diffusion limit; test_ld_thick_diffusive_limit pinned; EXTERNAL slope source still zeroed (ERR-061 / #247) | HISTORY | Increment C milestone + Issue #247 gap; the scope of "full LD" vs current implementation | L |

### _require_slab (lines 200–214)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 200–214 | _require_slab docstring | Guard: slab/Cartesian only; angular_upstream presence is geometry gate (same as DD) | CONTRACT | Enforcement logic; call-site guard | H |

### _ld_source_moments (lines 217–237)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 217–237 | _ld_source_moments docstring | Split (2, ng) moment source into (average, slope); DD single-moment (ng,) case; source weight-normalisation convention | CONTRACT | Input validation + call-site contract | H |

### _LDCellTerms dataclass (lines 241–276)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 241–276 | _LDCellTerms docstring | Per-cell Schur intermediates; slope elimination; ONE site of slope-row sign convention (LM-1989 memo §1.4/§6 correctness trap) | CONTRACT | Data structure for Schur algebra; docstring specifies the single sign-convention site | H |

### LinearDiscontinuous class (lines 279–370)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 280–291 | LinearDiscontinuous intro | Two spatial moments per cell; upwind-discontinuous closure; slope eliminated locally by Schur; fits scalar contract | TWIN | Anchor: `discretization.rst` closure table (LD row) + Schur form reference | H |
| 293–307 | is_affine_scannable attribute | LD admits single-upstream affine recurrence via Schur complement; supplies scan coefficients; CumprodScan admitted; Slab/Cartesian only | TWIN | Anchor: `discretization.rst` § on closure spectrum + `slab_one_group.rst` § "Generic affine outflow reconstruction" | H |
| 309–321 | spatial_basis_per_axis attribute | LD carries 2 moments per axis; per-cell 2^d unknown; per-face 2^{d-1} transverse moments; per-axis gate on multi-moment wiring (#240 D5b) | MOVED | Tensor-product Kronecker structure teaching — belongs in discretization.rst or `_ubld` theory | M |
| 323–337 | diffusion_limit_consistent attribute | Full LD thick-diffusion limit IS consistent (LM-1989 Part II §IV Eqs. 4.16-4.19); REQUIRES slope SOURCE threaded (Increment C, now landed); test pinned; PAIR validity with angular closure | TWIN | Anchor: `discretization.rst` closure table + LMM-1987 citation + the note on diffusion-limit consistency | H |
| 339–348 | supports_curvilinear attribute | LD slab/Cartesian ONLY; curvilinear unpublished (#158 arm / #6); raise NotImplementedError on curvilinear visit | CONTRACT | Selection logic guard | H |
| 350–362 | has_transpose_kernel attribute | LD adjoint faces TYPED DEFERRAL (#280); reverse-mode VJP (UBLD Schur + moment-frame involution) unimplemented; False guards honest is_adjointable + reverse-walk entry | MOVED | Belongs in adjoint/reverse-mode operations theory (when that page lands) | M |
| 364–369 | theta attribute | Slope-moment weight θ = 1/3 (LM-1989 Eq. 4.3b); SN-exact LD; NOT the LLD mass-lumped variant | CONTRACT | Closure constant; documented value + warning against LLD | H |

### _schur_terms method (lines 373–417)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 373–417 | _schur_terms docstring | Build per-cell Schur intermediates; single algebra site (Pattern 2) consumed by BOTH update + residual; delegates to d1_closed_form | CONTRACT | Method contract + the Pattern 2 claim (single algebra source) | H |

### update method (lines 419–443)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 419–443 | update() docstring | Solve LD cell system; average + slope reconstruction; return downstream face ψ_out = ψ̄ + ψ̂; slab only | CONTRACT | Method contract | H |

### residual method (lines 445–463)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 445–463 | residual() docstring | Apply-direction Schur-reduced scalar system; at update's ψ̄ vanishes to FP; linear in cell_avg; affine in source; slab only | CONTRACT | Method contract + round-trip invariant | H |

### Comment block (DAG-family batched kernel) (lines 465–497)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 465–497 | Comment block (DAG kernel context) | FullFieldWavefront oracle (SAME contract as DD); ÷V Cartesian streaming g = |μ|/Δ; unified moment kernel (#240 D5b-S3); d=1 Schur is fast path (CumprodScan); ÷V dense system scale-free; d=1 PRODUCTION closed form; flat-source artifact retired | HISTORY | #240 D5b-S3 milestone (unified moment kernel); the architectural layering (fast-path vs dense) | L |

### _ubld_system method (lines 498–547)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 498–547 | _ubld_system docstring | Assemble batched ÷V UBLD dense system + source-moment RHS; shape-agnostic source (moment vs flat); average-moment slot single-sourced | CONTRACT | Method signature + input/output contract + discriminator (moment rank vs batch rank) | M |

### _ubld_inflow method (lines 549–591)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 549–591 | _ubld_inflow docstring | Sum per-axis upwind-inflow RHS contributions; d-generic transverse-moment weighting; d=1 axis-append for scalar face (walk convention difference) | MOVED | Tensor-product Kronecker d-generic structure — teaching belongs in `_ubld` theory page | M |

### _ubld_outgoing_faces (lines 593–624)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 593–624 | _ubld_outgoing_faces docstring | Reconstruct d downstream faces from 2^d moment vector; trace at downstream node (s_a = +1); Kronecker layout reshaping + trace order + flattening | MOVED | d-generic tensor-Legendre reconstruction — teaching belongs in `_ubld` theory page | M |

### cell_kernel_batch method (lines 626–658)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 626–658 | cell_kernel_batch docstring | Pure batched LD cell SOLVE; d-generic DAG wavefront kernel; UNIFIED moment path for every d; d=1 system is 2×2 whose Schur is d1_closed_form (proven equal); storage-free | TWIN | Anchor: `linear_discontinuous.rst` theory page (if exists) + the Schur form proof reference | M |

### residual_kernel_batch method (lines 660–705)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 660–705 | residual_kernel_batch docstring | Pure batched LD operator residual; apply-direction; d-generic moment UBLD system; mass-diagonal normalization (#240 D5b-S3); round-trip at solved state; moment-source consistency | TWIN | Anchor: reciprocity/matvec section; the mass-normalization clause is load-bearing for slope rows | M |

### Comment block (Scan-family coefficients) (lines 707–776)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 707–714 | Comment block (scan context) | ×V "denom" convention (CollisionCache); SAME LD as ÷V kernel, scaled by V; CumprodScan ≡ FullFieldWavefront gate pins equivalence; slab/Cartesian only | HISTORY | #240 D5b-S1 Branch 2 milestone; the equivalence gate (two-paths test) | L |

### affine_scan_coefficients method (lines 715–775)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 715–775 | affine_scan_coefficients docstring | LD Σ_t-epoch scan coefficients (×V convention); Schur-reduced 2×2 in flat-source form; slab-only guard (curvilinear unpublished); single-source d1_closed_form helper | TWIN | Anchor: `linear_discontinuous.rst` theory page (if exists) + d1_closed_form reference | M |

### moment_scan_closure method (lines 777–810)

| lines | symbol | content-id | verdict | destination / anchor | conf |
|-------|--------|------------|---------|----------------------|------|
| 777–810 | moment_scan_closure docstring | Per-cell d1_closed_form for 1-D moment SCAN (#240 D5b-S3 OWED-2); slope-source-aware companion of affine_scan_coefficients; scan_slope_face_source + scan_reconstruct; Pattern 2 — ONE LD algebra site | MOVED | The moment-scan extension (slope-source threading in 1-D) belongs in `_ubld` or a dedicated LinearDiscontinuous theory page | M |

---

## Summary

### Verdict counts by file

**scheme.py (1,356 lines):**
- CONTRACT: 11 blocks
- TWIN: 14 blocks
- MOVED: 4 blocks
- HISTORY: 3 blocks
- COMMENT-cut: 1 block (separator comments too short)

**diamond.py (785 lines):**
- CONTRACT: 3 blocks
- TWIN: 14 blocks
- MOVED: 3 blocks
- HISTORY: 5 blocks

**linear_discontinuous.py (811 lines):**
- CONTRACT: 8 blocks
- TWIN: 8 blocks
- MOVED: 6 blocks
- HISTORY: 3 blocks

**Aggregate:**
- CONTRACT: 22 blocks (call-site contracts, selection logic, method signatures)
- TWIN: 36 blocks (theory pages carry these concepts; code mirrors them)
- MOVED: 13 blocks (teaching content that should move to/extend theory pages)
- HISTORY: 11 blocks (campaign/milestone narration — verify issue links)
- COMMENT-cut: 1 block (narration-only comments; consider replacing with named constant or reference comment)

### Docstring lines proposed keep vs cut

- **Keep (CONTRACT + TWIN):** ~550 lines total across all three files (call-site contracts, protocol definitions, field documentation, trait specifications)
- **Move to theory (MOVED):** ~180 lines (tensor-product structure, advection-reaction interface, scan-march extensions, moment-threading, adjoint operations)
- **Verify/archive (HISTORY):** ~90 lines (link to issue/milestone records; these are progress narration, valid context but not evergreen documentation)

### Posing flags

**Dated operator-algebra posings found:**

1. **scheme.py line 75:** "the assembled cell balance **is** the loss operator :math:`L+C`" — This is the honest algebra posing (A = L + C − S − B); verified against operator_algebra.rst. ✓ CORRECT.

2. **diamond.py line 27–32:** Comment-stated close to "ψ^s_out = 2·ψ_avg − ψ^s_in" (the diamond relation). References Lewis & Miller §5.3; the symbol uses the modern ψ notation. ✓ CORRECT.

3. **linear_discontinuous.py line 58–67:** "2×2 system" with explicit matrix stated (LM-1989 Eqs. 4.3a-b verified against cited reference). ✓ CORRECT.

4. **linear_discontinuous.py line 83–94:** Schur-reduced "S·ψ̄ = eff_source + eff_numer_upstream" (closed form stated; SymPy-verified claim in docstring). ✓ CORRECT.

No dated/superseded posings (e.g., "A = L − (S+B)" or legacy "Σ_r" fold claims) found in these files.

### Confidence levels of uncertain blocks (marked L)

1. **scheme.py line 1–27 (Why abstraction, HISTORY):** Round 1 Wave C is in git history; Issue #157/#158 are the campaign trackers. Confidence is mechanical archive link. | L

2. **diamond.py line 187–193 (has_transpose_kernel, HISTORY):** Wave O / O.2b milestones are in git; the reverse-walk transpose implementation is testable. Confidence is moderate (campaign is real). | L

3. **linear_discontinuous.py line 12–30 (Module note, HISTORY):** Increment C (Issue #240 D5b-S3) is recent work; test_ld_thick_diffusive_limit pinning is green. Confidence is high for the landed claim. | L (conservative rating because "current implementation" status dates fast)

---

## Next Actions

1. **Archive HISTORY blocks:** Link to GitHub issue/commit records; don't move to theory (they are progress narration).
2. **Move MOVED blocks:** Extend or create theory pages:
   - `discretization.rst` § new "Advection-reaction interface" subsection (generic model-agnostic design)
   - New `_ubld.rst` theory page (tensor-product Kronecker structure, moment layouts, d-generic UBLD system)
   - New/extend `linear_discontinuous.rst` theory page (Schur form, diffusion-limit consistency, slope-source threading)
   - Adjoint/reverse-mode operations page (when #280 lands)
3. **Verify CONTRACT & TWIN blocks:** These stay in code; cross-check method signatures against Sphinx domain refs (no broken links on push).

