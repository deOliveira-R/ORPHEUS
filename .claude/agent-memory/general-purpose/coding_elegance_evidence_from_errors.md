---
name: coding-elegance-evidence-from-errors
description: Evidence base for a coding-elegance skill, distilled from every ERR-NNN in `.claude/skills/vv-principles/error_catalog.md` (ERR-001..ERR-047 + ERR-026 Phase F closeout). For each bug: the anti-elegance that enabled it, the elegance pattern that would have made it structurally impossible. Then aggregate patterns, anti-patterns, and implicit lessons.
metadata:
  type: project
  source: .claude/skills/vv-principles/error_catalog.md (3263 lines, 47 entries)
  cross_ref: .claude/agent-memory/method-implementer/issue_168_phase_f_closeout.md (twin-path fix incompleteness)
---

# Coding-elegance evidence base, distilled from L0 error catalog

This memo answers the question: **for each ERR-NNN in the catalog, what
architectural / API / type-system / abstraction choice would have made the
bug structurally unrepresentable, not merely easier to catch?**

The deliverable is organized in four sections:

- **§1** — Per-bug elegance reading (47 ERR entries, grouped by failure
  mode).
- **§2** — Aggregate pattern catalog (the elegance patterns the bugs
  collectively evidence).
- **§3** — Anti-pattern catalog (the code shapes that recurred across
  multiple bugs).
- **§4** — Implicit lessons the catalog teaches but does not name.

Coverage tags below: **AE** = the anti-elegance enabling the bug; **EL**
= the elegance pattern that would have made it structurally impossible.

---

## §1. Per-bug elegance reading

### Group A — Twin-path duplication (the same physics implemented twice on different code paths, fix lands on one, sister rots)

#### ERR-006 / ERR-007 — α recursion + ΔA/w missing in BOTH sweep AND BiCGSTAB operator

- **What:** Wrong α recursion and missing `ΔA/w` geometry factor in
  curvilinear sweep; the identical physics bug also existed in the FD
  transport operator (ERR-007 explicitly notes "When a bug is found in
  one code path, check ALL code paths that implement the same physics").
- **AE:** `_sweep_1d_spherical` and `build_transport_linear_operator_spherical`
  are two implementations of the SAME continuous integro-differential
  operator. They were authored independently and the per-ordinate
  redistribution term (`α_{n+1/2}·ψ_{n+1/2} − α_{n-1/2}·ψ_{n-1/2}) · ΔS_i
  / (2 𝒲_n V_i)`) was algebraically duplicated. Fixing one did not fix
  the other.
- **EL:** **Single source of truth for the discrete operator** — extract
  the per-cell recurrence kernel into a single helper that both the
  forward sweep and the matvec call. The Hébert §3.9.4 cell update
  IS a function; making it a function (rather than two inlined copies)
  forces both paths through the same code.

#### ERR-026 manifestations 1–7 (the ERR-026 saga in toto) — Curvilinear sweep WDD vs apply WDD

- **What:** Sweep's one-directional WDD closure
  `ψ_{n+1/2} = (ψ_n − (1−τ)·ψ_{n−1/2}) / τ`
  ≠ apply's symmetric closure `ψ_{n+1/2} = τ·ψ_{n+1} + (1−τ)·ψ_n`.
  Both consistent for flat flux; both diverge for non-flat. Six phases
  (Wave H A→F) were needed to reconcile.
- **AE:** **Twin-path duplication on a non-trivial discretisation**.
  Three separate inlined math implementations of the same continuous
  operator: (1) `_sweep_1d_spherical`, (2)
  `transport_operator_matvec_spherical` (Krylov apply), (3)
  `build_transport_linear_operator_spherical` (legacy BiCGSTAB). Each
  evolved independently. The Phase F closeout names this
  explicitly: *"Twin-path fix incompleteness: Phase D fixed the apply-
  matvec Carlson seed but NOT the SI/sweep twin"*.
- **EL:** **Operator-as-object with strategy-pattern closures.** Phases
  A–F arrived at this incrementally: `PoleAngularClosure` Protocol with
  `LegacyTauSymmetricInterpolation` / `BaileyFlatFluxRedist` /
  `MorelMontryAngularSweep` strategies; `BoundaryFaceFlux` Protocol;
  `PsiHalfAngleSeed` Protocol; `carlson_inward_sweep_from_source` as a
  reusable helper. The lesson: **the cell update IS a Protocol**. Once
  it is one, the sweep and the apply consume the SAME implementation
  via the same Protocol, and fixing one fixes both. The Phase F factored
  helper `psi_half_angle_seed.carlson_inward_sweep_from_source` is the
  archetypal example: one function, two callers.

#### ERR-039 — `apply_transpose` claimed Π* = R but actually multiplied by (2ℓ+1)

- **What:** `HarmonicMomentProjection.apply_transpose` advertised
  `CAP_APPLY_TRANSPOSE` and delegated to
  `HarmonicMomentReconstruction.apply` (the addition-theorem
  reconstruction WITH the `(2ℓ+1)` factor) — which is NOT the
  W-weighted Hilbert adjoint.
- **AE:** **Two adjacent transforms that share a name root** ("moment
  reconstruction" vs "moment adjoint") were written as if they were
  duals when in fact they live in different inner-product spaces. The
  `apply_transpose` was implemented by **delegating to the wrong sister
  operation**.
- **EL:** **Type-encoded inner-product space**. If `Π: V → C` carried a
  type tag for the inner products on both ends (`V_inner` =
  weighted-by-w, `C_inner` = standard), then the adjoint operator's type
  would be `Π*: C → V` with the W-weighted inner product on the right,
  and a delegation to the `R: C → V` reconstruction would be a type
  error. Failing that, the **adjoint identity must be a constructor-
  level contract**: `apply_transpose` cannot be declared without a test
  `⟨Πψ, c⟩ == ⟨ψ, Π*c⟩` (the catalog's L1 test
  `TestApplyTransposeIsWWeightedAdjoint` is precisely that, retrofitted).

---

### Group B — Convention drift across module boundaries (the "untyped numpy array crossing module boundaries" family)

#### ERR-003 — Octant batching breaks reflective BC ordering

- **What:** Batching ordinates by sweep direction silently changed
  the order in which reflective BC boundary fluxes were updated.
- **AE:** **Sequential processing order was an unstated contract**.
  Reflective BCs depend on a partner-ordinate's flux being computed
  first, but this dependency was encoded only by ordering of operations,
  not by data.
- **EL:** **Make data dependency explicit in the structure**, not in the
  iteration order. A dependency graph (DAG) over ordinates / cells, or
  a `creates_sweep_cycle` flag on the BC (which Wave 0 of the SN
  refactor introduced), surfaces the constraint as a checkable invariant
  instead of a hidden ordering assumption.

#### ERR-008 — Half-cell vs full-cell volume convention drift

- **What:** `CartesianMesh.volume` halved first/last cell volumes for
  the old MATLAB port convention; SN solver used these volumes in the
  keff numerator/denominator, where the reflective-BC physics wants
  FULL volume.
- **AE:** **The same geometry object served two consumers with different
  semantic needs.** "Mesh for sweeping" and "Volumes for integral
  quantities" wanted opposite conventions; the object exposed one
  `volume` attribute.
- **EL:** **Domain-language API** — separate `sweep_volumes` and
  `integral_volumes`, each named for its consumer. The aliasing forces
  the developer to think about which they want.

#### ERR-011 — `fuel.r` vs `fuel.dz` mixed in gap geometry update

- **What:** MATLAB: `gap.r_ = (clad.r(1) + fuel.dz)/2` — a cladding
  inner radius (4.22 mm) added to a fuel axial node height (1.5 m).
- **AE:** **Variable names that differ only by suffix** with mixed
  scalar/vector indexing in the same scope. Pure stringly-named
  attributes carry no unit / dimension information.
- **EL:** **Dimensional types** (e.g. `Length` vs `RadialPosition` vs
  `AxialHeight`). With distinct types, `clad.r + fuel.dz` would be a
  compile-time type error.

#### ERR-022 — Negative lethargy bin width flips spectrum plot sign at three sites

- **What:** Three independent call sites (`mc.solver`, `homogeneous.solver`,
  `plotting.plot_moc_spectra`) divided non-negative tallies by signed
  `du[g] = log(eg[g+1] / eg[g])`, producing uniformly negative
  `flux_per_lethargy`.
- **AE:** **Definition-site convention violated at three independent
  consumer sites**. The "lethargy bin width" should be non-negative by
  convention; the energy-grid ordering convention is "descending" and
  three places silently relied on different sub-conventions.
- **EL:** **Normalize at the definition site, not at every consumer.**
  Take `np.abs(np.log(eg[1:] / eg[:-1]))` once, in the helper that
  produces `du`. "Width" is a measure; sign is a direction. Conflating
  them at consumers is anti-elegance.

#### ERR-034 — Slab Variant α first-leg trajectory: missing μ in `x_traj = x - μ·s`

- **What:** `x_traj = x - s_pts_first` (s as arclength) used as if s
  were x-distance. Invariant under any uniform source.
- **AE:** **Arclength `s` and Cartesian position `x` are different
  geometric quantities** (related by `Δx = μ·s`) but both stored as
  `float`. The naming `s_pts_first` does not encode the geometric
  contract.
- **EL:** **Type-encoded geometric quantities**. `Arclength` and
  `Position` as distinct types, parametrised by `μ`. The conversion
  `Position(x) - μ * Arclength(s)` becomes the only way to write the
  trajectory, and `Position(x) - Arclength(s)` is a type error.

#### ERR-037 — Atalay z_0 quadrature: `μ = tanh(t)` substitution avoids endpoint pole

- **What:** Bracket `1/(1-ν²)` and `g_1(ν) ~ ln²(1-ν)` algebraically
  cancel but `mp.quad` was 1.5% under-evaluated at dps=35.
- **AE:** **The cancellation is an algebraic fact, but `mp.quad` does
  not know it.** The relationship between the integrand's nominal pole
  and its true behaviour was not captured by the integrand's *type*.
- **EL:** **Choose the integration variable that makes the regularity
  manifest.** This is the same elegance principle as ERR-008's half-
  volume one: pick the representation that makes the contract trivial.

---

### Group C — Inelegant arithmetic that destroys an invariant numerically

#### ERR-005 — DD recurrence rewrite introduces catastrophic cancellation

- **What:** Algebraically equivalent rewrite `2*(…) − psi_in` instead
  of `0.5*(psi_in + psi_out)` averaging.
- **AE:** **"Algebraic equivalence" was assumed to imply numerical
  equivalence.** The unstable form subtracts nearly-equal large numbers.
- **EL:** **Stable forms encoded once**. Codify the canonical numerically-
  stable form in a single helper. Future "optimizations" must go through
  that helper.

#### ERR-020 — `cbrt(x)**3 != x` destroys equal-volume bit-invariant

- **What:** `_subdivide_zone` made spherical edges via
  `cbrt(inner³ + k/n · (outer³-inner³))`, `compute_volumes_1d` then
  re-derived volumes as `(4/3)π·diff(edges**3)`. ULP-level drift.
- **AE:** **A round-trip through `cbrt ↔ **3` was assumed bit-exact**.
  The construction route ("place equal-volume edges") and the
  reconstruction route ("derive volumes from edges") encoded the same
  invariant twice — but only one was numerically clean.
- **EL:** **Compute the invariant once and propagate**. Have
  `_subdivide_zone` return `(edges, volumes)` together; volume is a
  one-scalar-per-zone broadcast. **"When `op(op_inverse(x))` shows up,
  ask if it's bit-exact — if not, refactor to avoid the round trip."**

---

### Group D — Tautological / structurally-circular tests

#### ERR-016 — Tautological inner-iteration convergence residual

- **What:** GS inner loop computed
  ```
  phi_new = transported / denom
  res = ||denom * phi_new - transported||
  ```
  which is identically zero by construction.
- **AE:** **The "residual" compared a quantity to its own algebraic
  re-derivation.** The two sides of the residual were not independent.
- **EL:** **Compare independent quantities** in any convergence check.
  The fix `||φ_new - φ_old|| / ||φ_new||` compares two iterates, not
  a quantity to its inverse-applied self.

#### ERR-032 — Two "independent" derivations both used `∫E_2 = 1 − E_3` (wrong; should be `½ − E_3`)

- **What:** A closed-form analytical reference and a fixed-point
  cross-check **both inherited the same wrong antiderivative** from the
  developer's working memory.
- **AE:** **Procedural independence (two code paths) was confused with
  structural independence (two derivations from different upstream
  identities).** They shared the wrong upstream identity table.
- **EL:** **Cross-checks must be structurally independent**, not
  procedurally so. The row-sum identity (from the *kernel*) caught the
  closed-form error; the *closed form* and its derivative could not.
  Pick cross-checks from different integrands or different identities
  whenever possible.

#### ERR-035 — Phase-3A heuristic closure assumed sphere/cylinder formula generalizes to slab by analogy

- **What:** `ψ_surf = α · B_period / (1 - α²·e^{-2τ})` was an analogical
  generalization of the sphere/cylinder rank-1 form — coincides with
  the first-principles result only at α∈{0,1}.
- **AE:** **Analogical reuse of a closure formula across geometries
  without re-derivation.** The "two-bounces-per-period" reasoning
  motivated `α²` in the denominator, but the first-principles
  derivation has only one `α` and a single-transit `B`.
- **EL:** **Delegate the symmetric case to the asymmetric implementation
  at the symmetric corner.** The Phase-3B fix made
  `solve_greens_function_slab` a thin wrapper that calls
  `solve_greens_function_slab_asymmetric` with `α_L = α_R = α`. One
  derivation, two API entrypoints, no possibility of drift.

---

### Group E — Stringly-typed APIs / positional arguments

#### ERR-031 — `compute_P_ss_cylinder(radii, sig_t)` called with swapped positional arguments

- **What:** `compute_P_ss_cylinder(sig_t, radii)` silently produced wrong
  P_ss; "radii" happened to be increasing so no guard fired. Test was
  meaningless for months.
- **AE:** **Positional API with two arguments of the same type
  (`np.ndarray`).** No way for the type system to distinguish them.
- **EL:** **Keyword-only arguments** + **input contract validation**.
  After Issue #134 the `chord_quadrature` recipe validates "radii must
  be strictly increasing" at construction, which caught the test bug
  as a side effect. Stronger still: distinct types
  (`Radii = NewType('Radii', np.ndarray)` and `MacroXS = NewType(...)`)
  would make the swap a type error.

#### ERR-023 — MC `_random_walk` silently ignored `Sig2` because the mixture field was never read

- **What:** Two-branch `sig_t = sig_a + sig_s_sum` did not include
  `Sig2.sum`; the (n,2n) majorant fraction was always rejected.
- **AE:** **A mixture field (`Sig2`) was silently dropped on the floor**
  because no transport-kernel code path consumed it. Zero default values
  meant existing tests passed.
- **EL:** **Exhaustive destructuring on construction.** When a
  `Mixture` is constructed, ANY field that the consumers do not consume
  must trigger an explicit "what should I do with this?" decision. In
  Python: dataclass with `__post_init__` doing a fingerprint
  cross-check, or sum-type representation that forces exhaustive
  pattern-matching.

#### ERR-024 — MC tally claimed `flux_per_lethargy` but used scattering estimator

- **What:** `tally[ig] += w / sig_s_sum` is a well-defined scattering
  estimator (`response weighted by Σ_s`), but field was named
  `flux_per_lethargy`. Label and content disagreed.
- **AE:** **Naming a result field by what we WANT, not by what the code
  COMPUTES.** "Unbiasedness alone is not a specification" — a
  scattering estimator is unbiased for the scattering response, not
  for flux.
- **EL:** **Specification-driven naming with contract test.** If the
  field is `flux`, a spectral cross-check against the analytical
  eigenvector must gate the field. The catalog's `test_2g_flux_ratio_homogeneous`
  is exactly that gate, retrofitted to fail the scattering estimator.

---

### Group F — Hidden hardcoded constants / convention-dependent magic numbers

#### ERR-004 — Hardcoded `4π` in BiCGSTAB RHS

- **What:** `build_rhs` hardcoded `4π` for angular normalization;
  correct for Lebedev (Σw = 4π), wrong for Gauss-Legendre (Σw = 2).
- **AE:** **A quadrature-dependent constant was hardcoded** in a function
  that takes a quadrature object as an argument. The constant should
  have been derived from the quadrature.
- **EL:** **Derive constants from the inputs**. Replace `4π` with
  `quadrature.weight_sum` at the use site. Hardcoded numerical constants
  in equations are anti-elegant; the equation belongs to a quadrature
  and the quadrature carries its own normalization.

#### ERR-017 — Wigner-Seitz pitch formula doubled: `pitch = r_cell · √π · 2`

- **What:** Extra factor of 2 quadrupled cell area; subcritical when
  it should be supercritical.
- **AE:** **Pitch formula written from the user's understanding** instead
  of inverted from the existing factory convention `pwr_pin_equivalent`,
  which uses `r_cell = pitch/√π`.
- **EL:** **Inversion via the factory convention**, not parallel
  derivation. If `pwr_pin_equivalent` knows `r_cell = pitch/√π`, the
  inverse `pitch = r_cell · √π` should be a helper alongside it (or
  the factory should accept either parametrisation).

#### ERR-018 — MATLAB direction sampling uses uniform θ instead of isotropic

- **What:** `theta = π · rng.random()` (uniform on [0, π]) instead of
  `arccos(1 - 2ξ)`. Classified as intentional MATLAB carryover.
- **AE:** **A "formula that looks wrong" was preserved without
  documenting whether it's intentional.** No way for a reader to know
  if it's a bug or a deliberate approximation.
- **EL:** **Pin intentional simplifications with a test that verifies
  the INTENDED formula**, with a docstring explaining the divergence
  from the physically-correct form. The catalog's
  `test_direction_sampling` verifying `E[dir_x²] = 1/4` does exactly
  this.

---

### Group G — Closed-form / specialization opportunities missed

#### ERR-025 — DD cumprod recurrence: missing `−Σ_t` in numerator AND missing `1/W` source normalization

- **What:** Two compensating errors that cancel for homogeneous problems
  but fail at material interfaces. A canonical SymPy derivation existed
  in `orpheus.derivations.sn_balance.derive_cumprod_recurrence` — and
  the implementation did not visibly reference it.
- **AE:** **Implementation derived independently from the
  symbolically-canonical derivation.** Two parallel sources of truth.
- **EL:** **The derivation IS the source of truth.** Cardinal pattern:
  the SymPy module is the algebra of record (Branch 1) and the numpy
  implementation (Branch 2) is L1-verified against it. The comment
  `# canonical DD recurrence — see derive_cumprod_recurrence` is the
  minimum elegance; better is to AST-extract or doctest the symbolic
  derivation into the production module's tests.

#### ERR-027 / ERR-028 / ERR-029 / ERR-033 — Peierls kernel assemblies missing adaptive / closed-form treatment

- **What:** ERR-027: one-point GL collocation where adaptive is needed
  for cross-panel `E_1` integral. ERR-028: GL fails on a derivative
  discontinuity at `x'=x_i`. ERR-029: ρ/ω integration with no
  subdivision at panel crossings / tangent angles. ERR-033: finite-N
  GL fall-through where `½·E_2(τ)` closed form applies in slab-polar.
- **AE:** **Same `(½) ∫ E_n` integral evaluated by different methods
  at different call sites** — adaptive in one, plain GL in another,
  closed-form in a third. The integration strategy was site-local, not
  kernel-class-local.
- **EL:** **Strategy on the kernel class, not on the call site.**
  Either the kernel knows its own integration strategy (`SlabPolarKernel`
  has `evaluate_against(basis)` that returns the closed form when τ is
  μ-independent), or a Nyström assembly recipe dispatches on the
  kernel type. **"Closed-form detection" is a kernel-class invariant,
  not a site-local optimisation.**

#### ERR-036 — Log-singular kernel diagonal truncation in plain GL

- **What:** Discrete operator
  `K[i,j] = (c/2) Σ_k (w_k / μ_k) exp(-|z_i - z_j| / μ_k)` is qualitatively
  wrong at the diagonal (truncates `+∞` to `~2 log(n_μ)`).
- **AE:** **Sampling a singular integrand at zero** with no detection
  of the singular point.
- **EL:** **Product-Nyström as the canonical Nyström pattern for
  weak singularities**: split the kernel `K = K_smooth + K_singular`,
  integrate `K_singular` analytically against the trial basis. The
  fix module `peierls_atkinson_nystrom` makes this a reusable recipe.

#### ERR-030 — Rank-N white-BC: mode-0 / mode-n≥1 normalization mismatch

- **What:** rank-1 Mark and rank-N Marshak modes use different mode-0
  normalizations; result is bit-exact at N=1 by construction and
  catastrophically wrong at N=2 on MR.
- **AE:** **Backward-compatibility hack** that pinned the historical
  rank-1 result by routing mode-0 through a different normalization
  than modes n≥1. Two different bases summed into the same series.
- **EL:** **Single basis at all ranks.** Either (a) the legacy rank-1
  is preserved by absorbing its normalization into a `mode_0_correction`
  factor that is applied UNIFORMLY at every rank, or (b) the rank-1
  result is re-validated post-flip. Mixing two bases is anti-elegant.

---

### Group H — Boundary-condition contracts (ERR-040..ERR-047, the Wave 3 typed-error suite)

These bugs were caught by the Wave 3 typed-error refactor and document
a single elegance pattern: **make illegal BC states unrepresentable via
typed contracts**.

#### ERR-040 — Tangential ordinate silently classified as inflow or outflow

- **AE:** Untyped "inflow / outflow / tangential" classification — a
  consumer assumes the partition is strict but the data allows a third
  state.
- **EL:** **Sum-type for ordinate classification**:
  `Inflow | Outflow | Tangential`, exhaustive pattern-matching.
  Tangential is a first-class case, not a "small ε".

#### ERR-041 — Vacuum BC constructed against an outgoing trace (`Γ_+` instead of `Γ_-`)

- **AE:** Untyped trace orientation — `OutflowTraceSpace` and
  `InflowTraceSpace` were both `np.ndarray` to consumers.
- **EL:** **Type-encoded trace orientation**. The BC's domain is
  typed `Γ_-`; constructing it with a `Γ_+` is a type error.

#### ERR-042 — Reflection-index table inconsistent with quadrature weights

- **AE:** Reflection table built against `mu_x` ordering but applied to
  a quadrature whose ordering treated `mu_x = -mu_y`. The "measure-
  preserving" property of the reflection was an assumption, not a check.
- **EL:** **Invariants asserted at construction**: the
  `assert_geometry_map_measure_preserving` invariant raises
  `BoundaryGeometryMapNotMeasurePreservingError` at construction time,
  not at the first solve.

#### ERR-043 — Boundary response kernel produces negative output (sign-flipped construction)

- **AE:** Sign of `R` (the response kernel) was an arithmetic property,
  not a constructional invariant. Negative R is geometrically
  meaningless for white / albedo.
- **EL:** **Sign invariants checked at construction.** White / Albedo
  inherit a `assert_response_positive_if_declared` invariant that
  fires `BoundaryResponseNotPositiveError` if any output is negative.

#### ERR-044 — Reflection permutation is not an involution

- **AE:** Reflection table is a permutation, and "an axis reflection is
  its own inverse" was an unchecked property.
- **EL:** **Constructional invariant on permutations.**
  `PermutationOperator.is_involution` flag computed once; BC layer
  consults it and raises `ReflectionNotInvolutiveError` on construction.

#### ERR-045 — Reflection maps an inflow ordinate to itself rather than to its outflow partner

- **AE:** **Multi-invariant contract** that was checked piecewise (only
  involution, only inflow partition, only image surjectivity) — the
  intersection of all three was where the bug lived.
- **EL:** **Composite invariant** that asserts all three contracts
  simultaneously, at construction.

#### ERR-046 — Albedo / white kernel with α > 1 (sub-Markov violation)

- **AE:** Albedo accepted as `float` with no range check. "α is a
  float" is not a contract.
- **EL:** **Range constraints encoded in the type** — `Albedo = Annotated[float, ValidRange(0, 1)]`, or a smart constructor that
  raises `SubmarkovViolationError`. **"α ∈ [0, 1]" is the contract.**

#### ERR-047 — Boundary source `q` has nonzero entries on outflow trace

- **AE:** `q` (boundary source) shaped like the full ordinate set,
  trusting the consumer to mask the outflow side. Bug: consumer
  doesn't mask.
- **EL:** **`q` lives on `Γ_-`** as a typed object — the constructor
  rejects any array with nonzero entries on the outflow ordinates.

The cluster ERR-040..ERR-047 is the densest evidence in the catalog
for the **"make illegal states unrepresentable"** pattern. A typed BC
contract layer (Wave 3) raised eight distinct typed errors AT
CONSTRUCTION; previously the same misconfigurations would have produced
silently-wrong results at solve time (negative flux, wrong eigenvalue,
broken adjoint).

---

### Group I — Reference / paper contamination, dead code, third-party API

#### ERR-010 — pyXSteam viscosity NaN at T > 900 °C

- **AE:** **Third-party library guard mis-attributed as physical
  validity**. The IAPWS 2008 correlation is fine; the library imposed a
  cutoff.
- **EL:** **Wrap the third-party call** in a project-internal helper
  (`_iapws_viscosity`) that knows the actual physical validity range and
  exposes it as the contract.

#### ERR-012 — Static heat-transfer areas in deformable TH modules; `clad_a_bnd_def` computed but never used

- **AE:** **Dead code + stale data.** A correctly-computed deformed area
  variable existed but was never plumbed to its consumer.
- **EL:** **Use Nexus-style call-graph checks** at refactor time to
  surface "computed but never read". Better: the deformable-geometry
  pipeline should be a single pure function whose return value is the
  consumer's input — no chance to forget the wiring.

#### ERR-013 — Closed-gap stress BC uses fabrication gap width

- **AE:** **A geometric parameter used past its validity regime.**
  `gap_dr0` (fabrication) lived in `params` and was used in BC3/BC4 even
  after the gap closed onto roughness scale.
- **EL:** **Phase-dependent state via a sum-type** —
  `OpenGap(dr0) | ClosedGap(roughness)`. The BC's matchstatement
  selects the right gap thickness per phase. Or, even simpler: rewrite
  the BC as displacement-based, eliminating the gap-width division
  entirely (the catalog's `fix` paragraph does this — the elegant move
  was to remove the dependency rather than parametrise it).

#### ERR-014 — `sigT` stored vs `sigT` recomputed disagree because truncation is non-bijective

- **AE:** **Storing a derived quantity** (`sigT = sum of components`) when
  the components are already stored. Floating-point format truncation
  destroys consistency.
- **EL:** **Single source of truth for derived quantities**: either
  store only the components and always recompute, or assert at load time
  that `stored == recomputed`.

#### ERR-038 — Atalay 1997 first-order Fredholm precision floor at small slab thicknesses

- **AE:** **Reference treated as exact** when the paper's own text
  states "first-order approximation… expect improvement at small slab
  thicknesses." Wave 2-B chased it as a code bug for 10+ days.
- **EL:** **Reference precision is part of the reference's metadata.**
  When importing a published value table, capture the paper's stated
  precision claims in the docstring AND in the test tolerance — and
  cross-check against an independent reference (or against the same
  paper in a regime where the approximation is exact).

---

### Group J — Single-purpose minor bugs (still teaching something)

#### ERR-001 — Z-ordinate weight loss in 2D Lebedev sweep

- **AE:** **Special-case `continue` in a sweep loop** silently dropped
  ordinates whose weight contributes to the integral.
- **EL:** **Avoid in-loop dispatch on the data**; the per-ordinate
  contribution should be uniform, and "ordinates whose `mu_x = mu_y = 0`"
  is a special case that the *quadrature class* handles (by ensuring
  none, or by giving them a degenerate-but-correct treatment).

#### ERR-002 — Scattering matrix double-transposed: `phi @ SigS^T` instead of `phi @ SigS`

- **AE:** **Transpose-by-thinking-twice in vectorized code.** The
  identity `(A^T v) = (v^T A)^T` was applied wrong.
- **EL:** **Test with asymmetric inputs always.** Better: encode
  scattering-matrix orientation in the type — `SigS_g_to_gp` vs
  `SigS_gp_to_g` as distinct types, with the dispatch on
  `@` overloaded.

#### ERR-009 — CP `phi = P_inf @ source` instead of `P_inf.T @ source`

- **AE:** **Convention of `P[i,j]` (birth-first vs collision-first)
  silently encoded in the matvec.** The `.T` is the convention's bridge.
- **EL:** **Operator-overloaded transport step** — `solve_step(P, source)`
  hides the `.T` from the caller; the operator object knows its own
  convention.

#### ERR-015 — `compute_keff` ignored (n,2n) net neutron production

- **AE:** **Two consistency-required quantities** (transport-source side
  and eigenvalue-numerator side) were computed in different functions,
  and one was updated when (n,2n) was added; the other was not.
- **EL:** **Single function that returns both** — or: a single
  `production / removal` calculation that uses the same `Sig2`
  consumption pattern at both sites. Reduces the surface area for
  drift.

#### ERR-019 — Missing `4π · sin(θ)` weight in MOC scalar flux update

- **AE:** **Multi-factor weight formula** with three separable
  contributions (`omega_a`, `omega_p`, `t_s`, plus implicit `4π`,
  `sin(θ_p)`). Each was knowable, none was named explicitly.
- **EL:** **Named intermediates**: `angular_integration_weight =
  4π · omega_a · omega_p · sin(theta_p)`, then `update = weight · t_s ·
  delta_psi`. Naming the factors makes "did we include 4π?" a textual
  question, not an algebraic one.

#### ERR-021 — Degenerate ray tangent to corner → `IndexError`

- **AE:** **Bare `list[1]` indexing** with no checked precondition that
  the list has ≥2 entries.
- **EL:** **Return `tuple | None`** (or an explicit `Option` type) from
  geometric primitives that can return zero crossings; callers
  pattern-match on `None`. **"This ray contributes nothing" is a
  first-class outcome of ray tracing**, not an exception.

---

## §2. Aggregate pattern catalog

The 47 bugs cluster around 7 elegance patterns. Listed in
decreasing weight of evidence (number of bugs each would have
prevented):

### Pattern P1 — **Single source of truth for the same computation**

**What it does:** When the same mathematical operation appears in
multiple places (sweep + matvec; closed-form + numerical;
construction + reconstruction), extract it into ONE function and call
it from every site. Twin / triplet implementations of the same
math are anti-elegant; a fix to one will not propagate to the others.

**Prevents:** ERR-006/007 (sweep + BiCGSTAB), ERR-014 (stored vs
recomputed `sigT`), ERR-020 (edges vs derived volumes), ERR-022
(three call sites of `du`), ERR-026 (the entire Phase A..F campaign
to reconcile three discretizations), ERR-027/028/029/033 (Peierls
integration strategy split across call sites), ERR-035 (Phase-3A
heuristic vs Phase-3B first-principles), ERR-039 (apply vs
apply_transpose using wrong sister), Phase-F twin-path fix
incompleteness.

**Invocation (ORPHEUS-style):**

```python
# Anti-pattern: two implementations of the same cell update
def _sweep_1d_spherical(...):
    psi_face_out = 2*psi_avg - psi_face_in  # WDD spatial closure (inlined)
    ...

def transport_operator_matvec_spherical(...):
    psi_face_out = 0.5*(psi_avg_left + psi_avg_right)  # arithmetic avg (drifts)
    ...

# Pattern: single helper, two callers
def wdd_face_out(psi_avg, psi_face_in):
    return 2 * psi_avg - psi_face_in

def _sweep_1d_spherical(...):
    psi_face_out = wdd_face_out(psi_avg, psi_face_in)
    ...

def transport_operator_matvec_spherical(...):
    psi_face_out = wdd_face_out(psi_avg, psi_face_in)
    ...
```

### Pattern P2 — **Make illegal states unrepresentable (via types + constructional invariants)**

**What it does:** Bake contracts into the type system or the
constructor of a domain object so that violating the contract is
impossible (type error) or fails immediately (constructor raises). The
contract becomes structural, not behavioral.

**Prevents:** ERR-008 (volume convention), ERR-011
(`fuel.r` vs `fuel.dz`), ERR-022 (signed widths), ERR-031 (positional
arg swap), ERR-034 (arclength vs position), ERR-040..ERR-047 (the
entire BC contract suite), ERR-046 (`α ∈ [0,1]`).

**Invocation:**

```python
# Anti-pattern: positional args of same dtype
def compute_P_ss_cylinder(radii: np.ndarray, sig_t: np.ndarray): ...
# silently broken when called as compute_P_ss_cylinder(sig_t, radii)

# Pattern: smart constructor + range check on every input
@dataclass(frozen=True)
class Radii:
    values: np.ndarray
    def __post_init__(self):
        if not np.all(np.diff(self.values) > 0):
            raise ValueError("radii must be strictly increasing")
        if not np.all(self.values > 0):
            raise ValueError("radii must be positive")

# Now Radii cannot be confused with MacroXS; mistakes are type errors.
```

### Pattern P3 — **Cross-check against a structurally-independent reference**

**What it does:** A correctness claim should be verified against a
reference that does NOT share an upstream identity, antiderivative, or
operator with the implementation. Procedural independence (different
code paths) is not enough.

**Prevents:** ERR-005 (algebraic rewrite vs original form — they shared
the algebraic identity; only numerical-stability test saw the drift),
ERR-016 (residual compared to its own inverse), ERR-018 (intentional
simplification vs physical-correct check), ERR-024 (scattering vs
collision estimator), ERR-032 (two derivations shared `∫E_2`),
ERR-035 (Phase-3A vs Phase-3B both used analogical reasoning), ERR-038
(reference precision floor as code bug), ERR-039 (apply_transpose
delegated to the wrong sister).

**Invocation:**

```python
# Anti-pattern: cross-check from the same identity
phi_analytical = (1/(2 * sig_t)) * (1 + ...)  # uses ∫E_2 = 1 - E_3 (WRONG)
phi_fixed_point = iterate_to_convergence(...)   # also uses ∫E_2 = 1 - E_3
assert np.allclose(phi_analytical, phi_fixed_point)  # passes at 1e-39, both wrong

# Pattern: cross-check from a DIFFERENT integrand
# Row-sum identity is derived from the kernel itself, not its volume integral
row_sum_expected = (1/(2*sig_t)) * (2 - E_2(tau_i) - E_2(tau_i_prime))
row_sum_actual = sum(K[i, j] for j in range(N))
assert np.allclose(row_sum_actual, row_sum_expected)  # catches the algebra bug
```

### Pattern P4 — **Named intermediates over inline algebra**

**What it does:** Give every conceptual term in a multi-factor formula
a name with a comment about its role. The bug-surface area of named
intermediates is the misnaming; the bug-surface area of inline algebra
is *every factor*.

**Prevents:** ERR-001 (z-ordinate skip continue), ERR-005 (algebraic
rewrite), ERR-019 (missing `4π · sin(θ)` weight), ERR-025 (DD recurrence
factor-of-two errors), ERR-030 (mode-0 vs mode-n normalisation).

**Invocation:**

```python
# Anti-pattern: dense inline algebra; "did we include 4π?" is unanswerable
delta_phi[i] += omega_a * omega_p * t_s * delta_psi

# Pattern: named factors with comments
angular_integration_weight = 4 * np.pi * omega_a * omega_p * sin_theta_p
# 4π    : full-sphere normalization (Boyd 2014 Eq 45 outer factor)
# sin_θ : 2D→3D segment-averaging Jacobian
segment_contribution = angular_integration_weight * t_s * delta_psi
delta_phi[i] += segment_contribution
```

### Pattern P5 — **Operator-as-object (Protocol / strategy pattern), not bare numpy**

**What it does:** Wrap the discrete operator (BC, cell update, kernel,
closure, projection) as a class implementing a Protocol with an explicit
`apply` / `apply_transpose` / `assert_invariants` interface. Untyped
numpy arrays crossing module boundaries are anti-elegant because they
carry no contract.

**Prevents:** ERR-002 (`SigS` transpose), ERR-009 (`P` transpose),
ERR-026 (entire Wave H reaches this conclusion via `PoleAngularClosure`,
`BoundaryFaceFlux`, `PsiHalfAngleSeed` Protocols), ERR-039 (adjoint
identity as Protocol contract), ERR-040..ERR-047 (typed BC contracts),
ERR-043 (sign invariant on `R`).

**Invocation:**

```python
# Anti-pattern: bare matvec
phi_new = P_inf @ source  # which convention? (ERR-009: should be .T)

# Pattern: operator-as-object
class CPTransportStep:
    """P_inf is birth→collision; apply moves source from birth to collision."""
    def __init__(self, P_inf: np.ndarray):
        self.P_inf = P_inf
    def apply(self, source: np.ndarray) -> np.ndarray:
        return self.P_inf.T @ source  # convention captured ONCE here
    def apply_transpose(self, response: np.ndarray) -> np.ndarray:
        return self.P_inf @ response
    def assert_convention(self): ...  # invariant probe
```

### Pattern P6 — **Closed-form / specialization detection at the kernel class, not the call site**

**What it does:** When a kernel admits a closed-form evaluation in some
regime (e.g. slab-polar with μ-independent τ → `E_n`), encode that
detection AT the kernel-class level. Plain finite-N quadrature on a
specialisable kernel is a bug, not a tolerance choice.

**Prevents:** ERR-027/028/029/033 (Peierls slab and curvilinear kernel
assembly), ERR-036 (log-singular product-Nyström), ERR-037 (tanh
substitution for endpoint pole).

**Invocation:**

```python
# Anti-pattern: site-local finite-N GL
def compute_P_esc(radii, sig_t, n_quad=64):
    # GL on every angular integral, regardless of structure
    return finite_n_gl(integrand, n_quad)

# Pattern: kernel-class strategy
class SlabPolarKernel:
    def evaluate_against(self, basis: LagrangeBasis) -> np.ndarray:
        # τ is μ-independent in slab-polar → closed form
        return 0.5 * E_2(self.tau)  # exact, no quadrature

class CurvilinearAngularKernel:
    def evaluate_against(self, basis):
        # τ depends on ω → adaptive quadrature with subdivision at tangent angles
        return adaptive_quad_with_hints(integrand, hints=self.tangent_angles)
```

### Pattern P7 — **Normalize at the definition site**

**What it does:** Convention-dependent values (signs, normalizations,
quadrature weight sums, energy-grid ordering) should be fixed at the
ONE place where the value is defined. Every consumer that re-applies
the convention is an opportunity for the bug to resurface.

**Prevents:** ERR-004 (hardcoded 4π vs quadrature weight sum), ERR-008
(half-cell volume), ERR-014 (truncated sigT), ERR-018 (intentional
direction sampling), ERR-022 (signed du / de), ERR-025 (1/W
normalization), ERR-031 (positional swap caught by upstream
validation).

**Invocation:**

```python
# Anti-pattern: convention re-applied at each consumer
du_mc      = np.log(eg[1:] / eg[:-1])  # signed in MC solver
du_homog   = np.log(eg[1:] / eg[:-1])  # signed in homogeneous solver
du_plot    = np.log(eg[1:] / eg[:-1])  # signed in plotting helper
# Three places, three bugs (ERR-022).

# Pattern: convention fixed at definition site
def lethargy_bin_widths(eg: np.ndarray) -> np.ndarray:
    """Returns non-negative lethargy bin widths. Convention-independent."""
    return np.abs(np.log(eg[1:] / eg[:-1]))
# Every consumer calls this and gets the right sign by construction.
```

---

## §3. Anti-pattern catalog

Anti-patterns that recurred across multiple bugs. For each: the code
shape, the architectural reasoning, the substitution.

### A1 — **Twin-path duplication of physics math**

**Shape:** Two (or more) modules implement the same continuous math
on different code paths (sweep vs apply matvec; direct vs Krylov;
analytical vs fixed-point). The implementations drift independently.

**Evidence:** ERR-006/007, ERR-026 (six waves of reconciliation),
ERR-035, ERR-039, Phase F twin-path fix incompleteness.

**Architectural reasoning:** Each duplicate is a maintenance liability;
when math changes, *all* must be updated; when one is fixed, the other
silently rots. Cardinal Rule 2 violation: shared math should imply
shared code.

**Substitution:** Pattern **P1 (Single source of truth)**. Factor the
shared math into a helper or Protocol; both paths consume the same
implementation. If the two paths legitimately need different
discretizations (e.g. WDD sweep vs symmetric apply), document the
difference at the Protocol level and verify they agree on the
problems where they SHOULD agree.

### A2 — **Bare numpy across module boundaries**

**Shape:** Functions accept and return `np.ndarray` with shape /
convention encoded only in docstring or test fixtures. No type-system
help to catch swapped arguments, wrong conventions, or unit drift.

**Evidence:** ERR-002 (SigS transpose), ERR-009 (P transpose),
ERR-011 (r vs dz), ERR-022 (signed du), ERR-031 (positional swap),
ERR-034 (arclength vs position), ERR-040..ERR-047 (BC contracts).

**Architectural reasoning:** Numpy's `ndarray` is a single type for all
quantities. Distinct physical quantities (radii, cross-sections,
positions, arclengths) have distinct semantic meaning but identical
representation. The bugs hide because there is no static check.

**Substitution:** Pattern **P2 (Make illegal states unrepresentable)**.
Use `NewType`, frozen dataclasses, or Protocols to encode the
distinctions. At minimum: keyword-only arguments + constructor input
validation.

### A3 — **Convention re-applied at every consumer**

**Shape:** A convention-dependent value (sign of a width, normalization
of a quadrature, transpose of a matrix, half-cell vs full-cell volume)
is computed at multiple consumer sites independently.

**Evidence:** ERR-004 (4π hardcoded across solvers), ERR-008 (volume
convention), ERR-022 (signed du at three sites), ERR-025 (1/W
normalization).

**Architectural reasoning:** N consumers ⇒ N opportunities to drift
from convention. The convention should be a property of the producer,
not a requirement on the consumer.

**Substitution:** Pattern **P7 (Normalize at the definition site)**.
Fix the convention at the single helper that produces the value;
consumers receive a quantity that already obeys the convention.

### A4 — **Tautological / circular cross-check**

**Shape:** A "convergence check" or "cross-check" compares quantities
that are derived from each other by construction. Two "independent"
references both inherit the same upstream identity.

**Evidence:** ERR-016 (residual is identically zero), ERR-032 (two
derivations shared `∫E_2`), ERR-035 (Phase-3A heuristic shared
analogical reasoning with itself).

**Architectural reasoning:** A check that compares `f(x)` to `g(f(x))`
where `g` is the inverse of the step that produced `f(x)` is testing
nothing. Procedural independence (different code paths) does not
imply structural independence (different mathematical derivations).

**Substitution:** Pattern **P3 (Cross-check against structurally-
independent reference)**. Pick the reference from a different integrand
(kernel row-sum vs closed form), a different identity, or a different
limit of the problem.

### A5 — **In-loop algebra without named intermediates**

**Shape:** A multi-factor formula written as one inline expression in a
hot loop. Each factor is implicit; "did we include factor X?" is an
algebraic question, not a textual one.

**Evidence:** ERR-001 (z-ordinate skip + weight integration tangled),
ERR-005 (DD recurrence rewrite), ERR-019 (missing `4π · sin(θ)`),
ERR-025 (DD coefficients), ERR-030 (rank-N mode-0 normalisation).

**Architectural reasoning:** Inline algebra hides assumptions in the
arithmetic. Named intermediates make the assumptions textual and
greppable.

**Substitution:** Pattern **P4 (Named intermediates)**. Break the
formula into commented stages. The cost is one named variable; the
benefit is reviewer-visible assumptions.

### A6 — **Hardcoded numerical constants in physics code**

**Shape:** `4π`, `1/2`, `1.0`, etc. appear as numeric literals in
formulas that depend on quadrature / convention / problem setup. The
value is correct for the original test case and silently wrong
elsewhere.

**Evidence:** ERR-004 (4π is quadrature-dependent), ERR-017 (Wigner-
Seitz pitch with extra factor 2), ERR-020 (cbrt round-trip ULPs).

**Architectural reasoning:** Hardcoded constants couple the formula to
a specific test setup. The constant's *meaning* is lost; only its
*value* survives.

**Substitution:** Derive constants from inputs (`quadrature.weight_sum`
instead of `4π`). When the constant is genuinely universal (e.g.
`np.pi`), name it in a comment. Reserve numeric literals for truly
dimensionless quantities.

### A7 — **Site-local optimisation instead of kernel-class strategy**

**Shape:** Integration strategy (finite-N GL vs adaptive vs closed-form)
chosen at the call site, not at the integrand-class level. Same
integrand assembled by different methods at different sites.

**Evidence:** ERR-027/028/029/033 (Peierls slab + curvilinear kernel),
ERR-036 (log-singular Atkinson), ERR-037 (tanh substitution).

**Architectural reasoning:** Integration strategy is a property of the
integrand's smoothness / singularity structure, not of the site that
invokes it. Decisions made at the site rot when new sites get added.

**Substitution:** Pattern **P6 (Closed-form / specialization at the
kernel class)**. The kernel knows whether it admits a closed form;
the consumer calls `kernel.evaluate_against(basis)` and trusts the
kernel.

### A8 — **Stored derived quantities**

**Shape:** A value `Y = f(X_1, X_2, ...)` is computed from components,
then stored alongside the components. Consumers receive both and trust
the stored `Y`; truncation / serialization / refactor can break the
invariant.

**Evidence:** ERR-014 (sigT stored vs recomputed), ERR-020 (edges and
volumes both stored).

**Architectural reasoning:** Duplicate state requires synchronisation
discipline. The stored value can desync from the components silently
under any non-bijective transformation.

**Substitution:** Either store only the components (recompute `Y` on
demand) or assert `stored == recomputed` at load. Pattern **P1
(Single source of truth)** applies to data as well as code.

### A9 — **Untyped capability advertisement**

**Shape:** An operator advertises a capability flag (e.g.
`CAP_APPLY_TRANSPOSE`) but the only test of that flag is "set is
non-empty"; no direct contract test exercises the production code path.

**Evidence:** ERR-039 (`apply_transpose` claimed `R` but delegated to
the wrong sister).

**Architectural reasoning:** Capability advertisement is a contract; an
untested contract is a lie.

**Substitution:** Every advertised capability must ship with a
direct, production-code-path test against an independent reference
(adjoint identity for `apply_transpose`, solve round-trip for `solve`).

### A10 — **Dead-code-adjacent stale references**

**Shape:** A correctly-computed value sits in scope unused while a
stale value is consumed downstream. Classic "computed but never read"
fingerprint.

**Evidence:** ERR-012 (`clad_a_bnd_def` computed and ignored), ERR-013
(`gap_dr0` used after gap closure).

**Architectural reasoning:** Dead code is a refactor hazard and a
red flag for incomplete plumbing. The computed-but-unused variable IS
the bug fingerprint.

**Substitution:** Use Nexus call-graph analysis at refactor time to
detect "produced but never consumed". Better: organise the dataflow as
a pure function pipeline so a missing wire becomes a missing argument
(type error).

---

## §4. Lessons the catalog implicitly teaches

These principles emerge across the bug history but are not stated as
named lessons in the catalog text.

### L1 — **Symmetry in math should imply symmetry in code (and absence of symmetry should imply absence of duplication)**

The catalog repeatedly demonstrates that a *fix* applied to one of two
structurally-mirrored production paths almost never propagates to the
sister (ERR-006/007 sweep vs operator, ERR-026 apply vs sweep across
6 phases, ERR-039 apply vs apply_transpose, the Phase F closeout's
explicit "Anti-pattern: Twin-path fix incompleteness"). The lesson:
**when math is symmetric, code should be symmetric — via a single
shared implementation, not via independent code paths that happen to
implement the same equation.** Conversely, when math is NOT symmetric,
trying to fold two paths into one (e.g. SI default vs Krylov default
for curvilinear) is wrong — the asymmetry is genuine and should be
expressed in the type system (e.g. as different Protocol strategies).

### L2 — **"It works on the homogeneous case" is necessary but not sufficient**

Eleven distinct bugs (ERR-002, ERR-006, ERR-009, ERR-015, ERR-019,
ERR-022, ERR-023, ERR-024, ERR-025, ERR-026, ERR-030) survived
homogeneous tests because the bug term cancelled, the flux was flat,
the eigenvalue was shape-invariant, or the symmetric case degenerated
into another problem. **Every solver should ship a heterogeneous /
multi-region / asymmetric-scattering test as a primary V&V gate, NOT
as a follow-up.** The catalog's Meta-Lesson 2 ("Homogeneous is
degenerate") and Meta-Lesson 4 ("20 passing tests don't mean correct")
state this in aggregate; the per-bug evidence is overwhelming.

### L3 — **Cross-cutting concerns belong in the type system, not in coding discipline**

The Wave 3 typed-error refactor (ERR-040..ERR-047) takes 8 distinct
contracts that previously lived as "things you have to remember" and
encodes each as a typed error raised at construction. The eight cases
all share a single pattern: **the property was always meant to hold;
discipline cannot enforce it; the type system can.** This generalises:
units, conventions, ranges, invariants, capability declarations should
all live in types or constructors, not in reviewer attention.

### L4 — **Probe the boundary between phases / regimes / corners**

Multiple bugs hid because tests only exercised the *corners*: α=0 or
α=1 (ERR-035); homogeneous or single-material (most curvilinear bugs);
zero or non-zero Sig2 (ERR-015, ERR-023). **The bug-rich region is
between the corners — partial reflection, multi-region, intermediate
values of any parameter.** A "rank-1 → rank-2" pin (ERR-030) and a
"α∈{0.5}" pin (ERR-035) catch what α∈{0,1} cannot. The catalog's
implicit lesson: design tests against the *interpolation* between
known-good corners, not just the corners themselves.

### L5 — **Tools that detect "computed but not consumed" reveal bug seeds**

ERR-012 (`clad_a_bnd_def` computed but never used) is the cleanest
example: a correctly-computed value was orphaned. ERR-023 (`Sig2`
field of `Mixture` never read in the transport kernel) is the dual: a
correctly-stored input was orphaned. Both bugs are detectable by
call-graph / use-def analysis (Nexus). **Static analysis of "is this
produced value consumed?" and "is this required input read?" is a
project-wide elegance practice**, not a one-off cleanup.

### L6 — **Catastrophic compensation is the most insidious bug class**

ERR-025 (DD recurrence) has TWO factor-of-two errors that compensate
on homogeneous problems and only break at material interfaces. ERR-006
combined a wrong α with a missing ΔA/w that gave correct homogeneous
results. ERR-035 had a closure formula that coincides with truth at
α=0 and α=1. **When two bugs interact to give correct results in
common cases, every test except the asymmetry-stressing one will
agree with the implementation.** The user-memory note "Fix bugs
immediately" applies here: deferring even ONE known bug risks letting
a second bug interact with the first, hiding both.

### L7 — **Documentation in the codebase IS code architecture**

Several bugs cite a missing or out-of-sync source-of-truth pointer
(ERR-025: derivation existed but production code did not reference it;
ERR-026: Bailey-Morel-Chang 2010 citation was wrong; ERR-018: a known
simplification was undocumented; ERR-038: paper precision floor was
under-read). **Sphinx docs / module docstrings / source comments are
part of the elegance contract** — they make the implicit explicit and
prevent re-derivation drift. Cardinal Rule 3 ("Sphinx IS the LLM's
brain") covers this from the documentation side; the catalog evidence
backs it from the bug-prevention side.

### L8 — **Anti-pattern not yet named in the catalog: "Specialised method that hides the general method"**

This is the ERR-030 / ERR-035 family pattern that the catalog does
NOT yet name. When a specialised closure (rank-1 Mark) is shipped
before the general one (rank-N Marshak), and the specialised one is
calibrated to "look right" by historical accident, the general method
is then constrained to reproduce the bit-exact specialised result —
which forces a mode-0 normalisation mismatch (ERR-030) or an
analogical heuristic (ERR-035). The **right** development order is
general-first: derive the rank-N case from first principles, then
specialise to rank-1 by collapsing modes. The catalog should add this
as a named anti-pattern: **"Premature specialisation that constrains
the general method"**. (Issue #131 / Issue #132 are tracking the fix
for this; the pattern is not yet documented.)

### L9 — **Anti-pattern not yet named: "Sequential ordering as silent contract"**

ERR-003 (octant batching breaks reflective BC ordering) and ERR-044/045
(reflection partnership in BC sweep) share the pattern that *order of
operations* encodes a data dependency. Order is fragile under
parallelisation, vectorisation, or refactor. **Encode the dependency
as a DAG, not as an iteration order.** The Phase C sweep-frame
rewrite (`SNMesh.iter_cells_by_direction`) is the elegant resolution —
the DAG ordering replaces "the loop happens to visit cells in a good
order".

### L10 — **Anti-pattern not yet named: "Approximation tolerance accepted as feature"**

ERR-036 (log-singular Atkinson, "loose tolerance 5e-2 because singularity-
aware is out of scope") and ERR-038 (Atalay precision floor accepted
as 5e-2 floor) both have the pattern: **a known approximation gap is
papered over by loosening the test tolerance, rather than flagged as
follow-up work.** The right elegance move: **the tolerance is a
contract; if the contract is loose, document WHY in the docstring AND
the test, and pin the gap with a structurally-independent reference
to prove it is a genuine approximation and not a bug**. The catalog
already pins both with explicit "this is the X approximation floor,
not the solver floor" tests — but the pattern itself is not named.

---

## Cross-references back to evidence

Every claim in §2 and §3 is backed by one or more ERR-NNN entries from
the catalog. The complete map:

- **P1 (single source of truth)**: ERR-006, ERR-007, ERR-014, ERR-020,
  ERR-022, ERR-026, ERR-027, ERR-028, ERR-029, ERR-033, ERR-035, ERR-039
- **P2 (illegal states unrepresentable)**: ERR-008, ERR-011, ERR-022,
  ERR-031, ERR-034, ERR-040, ERR-041, ERR-042, ERR-043, ERR-044, ERR-045,
  ERR-046, ERR-047
- **P3 (structurally-independent cross-check)**: ERR-005, ERR-016,
  ERR-018, ERR-024, ERR-032, ERR-035, ERR-038, ERR-039
- **P4 (named intermediates)**: ERR-001, ERR-005, ERR-019, ERR-025,
  ERR-030
- **P5 (operator-as-object)**: ERR-002, ERR-009, ERR-026, ERR-039,
  ERR-040..ERR-047
- **P6 (closed-form at kernel class)**: ERR-027, ERR-028, ERR-029,
  ERR-033, ERR-036, ERR-037
- **P7 (normalize at definition site)**: ERR-004, ERR-008, ERR-014,
  ERR-018, ERR-022, ERR-025, ERR-031

- **A1 (twin-path duplication)**: ERR-006/007, ERR-026, ERR-035, ERR-039
- **A2 (bare numpy across boundaries)**: ERR-002, ERR-009, ERR-011,
  ERR-022, ERR-031, ERR-034, ERR-040..ERR-047
- **A3 (convention re-applied at every consumer)**: ERR-004, ERR-008,
  ERR-022, ERR-025
- **A4 (tautological cross-check)**: ERR-016, ERR-032, ERR-035
- **A5 (in-loop algebra)**: ERR-001, ERR-005, ERR-019, ERR-025, ERR-030
- **A6 (hardcoded constants)**: ERR-004, ERR-017, ERR-020
- **A7 (site-local quadrature strategy)**: ERR-027, ERR-028, ERR-029,
  ERR-033, ERR-036, ERR-037
- **A8 (stored derived quantities)**: ERR-014, ERR-020
- **A9 (untyped capability advertisement)**: ERR-039
- **A10 (dead-code-adjacent stale refs)**: ERR-012, ERR-013

47 / 47 ERR entries cited at least once.
