# D — Flux ontology: torsor architecture vs cone proposal — blast radius + Phase-C substrate

Explorer memo, 2026-08-19, measured at HEAD `7aae9bf1` (main). All `file:line` are
at-HEAD and re-derivable; the durable claims are the structure. **Inventory only —
no ruling on which ontology is right.** Siblings: space-layer identity/metric facts
→ memo A; operator-arrow machinery → memo B; kernel binding → memo C. This memo
touches those seams only where the flux ROLE algebra reaches them.

**Headline measurements.** The plan's item 0.7 scope ("sweep the prose from
`coupled_system.py` and the member-family docs") is wrong by an order of magnitude,
confirming the brief's probe: `[M]` `grep -rl "torsor" … --include="*.py" --include="*.rst"`
→ **16 production files, 5 theory/verification pages, 16 test files** carry it, and
the torsor is a dunder-enforced runtime algebra (7 flux leaves × overridden
`__add__`/`__sub__`, a 7-leaf displacement package, a Rep-keyed registry, container
delegation through 3 composite layers, and 2 production consumers of the
displacement TYPE). Meanwhile the cone side is not empty either: a scheme-level
`is_positivity_preserving` flag already ships (DD honestly `False`, with a numerical
witness test), a cone-membership refusal ships in the BC realizer, coefficients
already have an explicit cone-as-predicate doctrine WITH its own test battery, and
`power_iteration` already implements ray normalization (unit production rate).

---

## 1. Torsor-architecture inventory (Q1)

### 1.1 The enforcing types

| type / mixin | file (at HEAD) | enforces | notes |
|---|---|---|---|
| `FluxRole` mixin | `orpheus/transport/fields/_flux_role.py:59` | `__sub__`: flux ⊖ flux → mints sibling displacement (the ONLY mint). `__add__`: flux ⊕ displacement → flux (torsor action; mesh-identity + shape checked by `_check_torsor_partner` :186); flux + same-class flux → `TypeError` naming `affine_combination` (:107-116); cross-class → Layer-1 gate. `affine_combination(pairs)` (:125): Σλᵢ=1 enforced to 1e-12, the only legal multi-flux blend. | **Scalar scaling deliberately UNTOUCHED** (:26-30): `ψ·c`, `ψ/k`, unary `−` inherited from `Field`. `_DISPLACEMENT_CLS` derived via `Displacement.sibling_of(Rep)` (:69-82). |
| `Displacement` marker + diagnostics | `orpheus/transport/displacements/_displacement.py:48` | Vector-space algebra INHERITED (V is a genuine vector space). Adds `contraction_ratio` (:92), `true_error_estimate` = ‖Δψ‖/(1−ρ) (:120), `where_largest` (:150). `_BY_REP` registry via `__init_subclass__` (:74-78). | Its own docstring (:64-66): the SI STOPPING test stays flat `np.linalg.norm(to_flat())`; the typed methods are diagnostics. |
| `CoefficientRole` mixin | `orpheus/transport/fields/_coefficient_role.py:73` | Overrides NOTHING — the *absence* of the affine gate is its content: coefficients (`CrossSectionField`) are a plain vector space, `σ+σ` legal, `Σ=0` IS the origin. | The complement table (:11-18) is the tree's own statement that the torsor is flux-specific. |
| `Field._check_partner` (Layer 1) | `orpheus/numerics/field.py:244` | `type(self) is type(other)` + equal `space`; cross-class arithmetic forbidden even at equal units ("same units grant permission to add in linear algebra, not meaning"). | The affine gate sits ON TOP of this; a cone overturn keeps or removes this layer independently (memo A's seam). |
| `Composite` (container) | `orpheus/transport/full_field.py:348-364` | `_check_partner` container-level only; member dunders are the single source of truth — deliberately does NOT pre-check member types "(that would BLOCK the legitimate composite torsor flux + displacement)" (:358). | Same discipline in `TimedFullField` (`_recombine` override only) and `RadialCharacteristicField` (`radial_characteristic_field.py:105-108` "the affine torsor propagated to both blocks"). |
| `CoupledField` | `orpheus/numerics/coupled_system.py:231, 292` | Member-wise delegation; container `_check_partner` arity-only, "pre-checking member types here would block the legitimate member torsor" (:292). | Module-docstring closing paragraph carries the verbatim prose the plan cites (§6 below). |

### 1.2 The leaves (the role grid's Flux and Displacement columns)

`[M]` `grep -rn "class .*FluxRole" orpheus/` → **7 flux leaves** mix in `FluxRole`:
`AngularFlux`, `ScalarFlux`, `HarmonicMomentFlux`, `AngularBoundaryFlux`,
`ScalarBoundaryFlux`, `RadialCharacteristicInteriorFlux`,
`RadialCharacteristicBoundaryFlux`. **7 displacement leaves** mirror them
(`orpheus/transport/displacements/__init__.py`): Angular, Scalar, Moment,
AngularBoundary, ScalarBoundary, RadialCharacteristicInterior/Boundary — Rep-keyed
1:1 (the `MomentDisplacement` ↔ `HarmonicMomentFlux` name asymmetry is why the
pairing is structural, not nominal). `AngularDisplacement.integrate_angular`
(`angular_displacement.py:49`) is the tangent map of the flux reduction — one shared
body, role-preserving on displacements; the consistent-DSA restriction consumes
exactly this (#2, landed `ec116e75` 2026-07-26).

The two `RadialCharacteristic*Displacement` leaves state their own forcing:
"Its existence is FORCED by the composite torsor algebra: System B's flux composite
ψ₂ − ψ₁ mints a displacement PER BLOCK" (interior :16-19, boundary :17). I.e. every
new block locus under the torsor mints a displacement sibling by construction —
the mechanism that grew 7 leaves from #208's original 4.

Source/sinks are deliberately BARE (vector space): `harmonic_moment_source_sink.py:27-42`
— "rate densities add vectorially … CLOSED, with no affine/torsor gate and no
displacement mint". The answer to the brief's "why?": the role split IS the design —
state (affine) vs rate-density (vector), same storage family.

### 1.3 Which tests gate it

- Foundation battery: `tests/numerics/test_affine_flux_algebra.py` — `[M]`
  `grep -c "def test_"` → **10 tests**, parameterized over leaves (mint bit-exact,
  torsor round-trip to 8 ULP, telescoping, gate, Σλ=1, diagnostics).
- Gate tests naming the forbidden add: `[M]` `grep -rln "add_flux_flux_forbidden\|flux_add_flux_forbidden\|flux + flux" tests/` →
  **7 files** (test_timed_full_field, test_angular_flux, test_angular_boundary_flux,
  test_affine_flux_algebra, test_harmonic_moment_flux, test_typed_source_sinks,
  test_2d_l2_matvec_correctness).
- `[M]` `grep -rn "raises(TypeError" tests/numerics/test_affine_flux_algebra.py tests/transport/ tests/sn/primitives/test_harmonic_moment_flux.py | wc -l` →
  **45 lines** (predicate: ALL TypeError raises in those trees — includes cross-class
  and container gates, not only the affine one). Narrower: `[M]`
  `grep -rn "raises(TypeError" tests/ -A2 | grep -c "affine\|origin"` → **16**
  raise-sites tree-wide whose following 2 lines name the affine/no-origin gate.
- `[M]` `grep -rln "Displacement" tests/ --include="*.py"` → **9 files / 36 lines**
  (incl. `test_flux_displacement_diagnostics.py` — SI wiring; `test_dsa_low_order.py`;
  `test_psi_half_coupling.py:3632` fails if `contraction_ratios` is empty on the
  coupled iterate — a wiring pin).
- Torsor-spelled operator-law gates: `tests/sn/operators/test_declared_law_is_linear.py`
  (§4.3 below), `test_operators_apply_typed.py`, `test_streaming/scattering/collision_operator.py`,
  `tests/sn/sweep/core/test_phase_c_gates.py:260` ("affine additivity in torsor form").
  ⚠ **Name collision**: `test_phase_c_gates.py` is issue **#168**'s Phase C
  (sweep-frame matvec gates), not the plan's cone-battery Phase C — the campaign
  should name its battery to avoid colliding with this shipped file.

---

## 2. Displacement-consumer map (Q2) — what structurally depends on the TYPE

**Production consumers outside `transport/{fields,displacements}`:** `[M]`
`grep -rn "Displacement" orpheus/ | grep -v transport/displacements | grep -v transport/fields` →
exactly **2 structural consumers** + 1 doc mention:

1. **`orpheus/numerics/iteration.py` (SourceIteration)** — the convergence machinery.
   - `_DisplacementLeaf` Protocol (:399): duck-typed on `l2` + `contraction_ratio`
     (numerics MUST NOT import transport — L1↛L2). "Deliberately consumer-minimal."
   - `_flux_displacement_leaf` (:421): unwraps `.interior` / recurses into
     `CoupledField.systems[0]` to find the leaf.
   - The SI loop (:781-793): `displacement = psi - psi_prev` — the FluxRole mint
     fires EVERY typed iteration; `contraction_ratios` list + `last_displacement`
     populated. **Comment: "DIAGNOSTICS only … the STOP rides the residual above."**
     The ρ-honest stop (:766-775) is `‖rhs_prev − rhs‖/‖q‖` — a SOURCE-typed
     subtraction (vector-space legal) + flat `_l2_norm`. So: the stopping logic
     consumes plain norms of source differences; the displacement TYPE feeds the
     ρ/‖Δψ‖/(1−ρ) diagnostic surface, not the stop.
   - The corrector step (:773): `psi = psi + self.corrector.apply(psi - psi_prev)` —
     corrector input AND output are displacement-typed; the update is the torsor add.
2. **`orpheus/sn/acceleration/dsa.py` (DSACorrection.apply :602-696)** — the hard
   structural consumer. `isinstance(interior, (AngularDisplacement, AngularFlux))`
   or TypeError (:635); the returned composite is MINTED as
   `AngularDisplacement.from_mesh` + `AngularBoundaryDisplacement.from_face_arrays`
   "so the torsor add ψ ⊕ Δψ is well-formed on both postures (a flux-typed boundary
   would trip the affine flux+flux gate)" (:680-684). ⭐ Its docstring records the
   Krylov leak: "GMRES vectors are Krylov directions whose role is erased at the
   scipy `from_flat` boundary (the template rebuilds them flux-typed), and **the
   swept vector IS the displacement from zero**" — Krylov's V-space algebra runs
   on flat ndarrays outside the type system, and V-vectors circulate flux-typed.

**The #340 convergence records do NOT consume the type**: `IterationRecord` /
`StoppingCriterion` / `IterationBudget` (`orpheus/numerics/convergence.py:910/418/691`)
carry `float` trajectories only. `power_iteration` (`numerics/eigenvalue.py:374`)
keeps `keff_history: list[float]` + named criteria — no field objects retained.

**`affine_combination` has ZERO production callers** — `[M]`
`grep -rn "affine_combination" orpheus/` → only the definition + the TypeError
message text. Its consumers are 5 test files. The relaxation it canonicalizes is
not spelled anywhere in production (SI has no relaxation step; the corrector path
uses the torsor form).

**Sites citing the doctrine as a constraint:** `orpheus/sn/solver.py:2960` — the
adjoint-eigenvalue warning path omits its balance-defect clause partly because
"the `1/k` scaling of a field crosses the affine-torsor arithmetic rules (a flux
state is not a vector …)", tracked as **#353**. (As written; note the shipped
algebra actually permits scalar scaling — the comment cites the doctrine, not the
dunder that would fire.)

---

## 3. Cone-machinery inventory at HEAD (Q4)

`[M]` `grep -rni "cone\|positivity\|is_positive\|nonnegative" orpheus/numerics/ orpheus/transport/ orpheus/sn/ | grep -vi test` → the hits sort into six clusters:

1. **Scheme-level positivity flag — EXISTS, DECLARED, UNREAD.**
   `DiscretizationScheme.is_positivity_preserving` (Protocol `transport/spatial/scheme.py:536`,
   ABC ClassVar :783; prose :446-448 "DD is **not** … Step is", :528 "gates
   negative-flux diagnostics"). Values: DD `False` (`diamond.py:129`, with :131-133
   "can produce negative outgoing face flux even when source and inflow are
   nonnegative"), LD `False` (`linear_discontinuous.py:233`). `[M]`
   `grep -rn "is_positivity_preserving" orpheus/` minus declarations → **0 production
   readers**; the ":528 gates negative-flux diagnostics" claim is aspirational — no
   such diagnostic exists. Tests: tag-pins (`test_diamond.py:139,145`,
   `test_linear_discontinuous.py:216`, `test_discretization_scheme_protocol.py:261-297`)
   **plus one genuine behavioral witness**: `test_diamond.py:793-830`
   `TestPositivityFailure.test_thin_cell_large_source_can_produce_negative_outgoing`
   constructs ΔxΣ_t ≫ 2|μ| and demonstrates negative ψ_out from positive inputs.
2. **Negative fluxes today: no fixup, no warning, nothing.** `[M]` grep "fixup|clip"
   over `orpheus/sn/ orpheus/transport/spatial/` → only prose mentions ("Step's
   positivity-fixup" in a trait docstring, "the positivity-fixup menu" in LD's
   references). DD output flows through unexamined.
3. **A production cone-membership REFUSAL at the BC realizer.**
   `orpheus/sn/boundary/realizer.py:792-812` — "REFUSAL AXIS: state-cone / sign":
   `ZeroFluxBoundary` (𝒜 = −1) raises `BoundaryError` because "the SN trace is an
   angular flux with ψ ≥ 0, and a negative inflow is outside that cone … (ψ ≥ 0 ⟹
   J± ≥ 0)". Also `sn/mesh/augmented_mesh.py:194` "a negative angular inflow is
   unrepresentable — use vacuum".
4. **The coefficient cone doctrine + its battery — the cone-as-predicate precedent.**
   `cross_section_field.py:35-37`: "**Nonnegativity is the cone, a property — not a
   type invariant.** … verify the cone is *closed* under the cone operations
   (+, λ≥0·)"; construction does NOT enforce Σ ≥ 0 (:95); `material_xs_field.py:694`
   "cone-valued — but the cone is a [property]". Tests:
   `tests/transport/fields/test_coefficient_fields.py` `TestCrossSectionConeAlgebra`
   (`test_cone_closed_under_addition:77`, `test_cone_has_origin:96`,
   `test_is_vector_space_not_torsor:107`). Same posture one method over:
   `scalar_boundary_flux.py:32` "positivity of J± is a property of the PHYSICAL laws
   𝒜∈[0,1], not a type invariant".
5. **Quadrature-weight positivity — a hard admission contract (different cone).**
   `numerics/quadrature/rules_sphere.py:578` "⛔ Positivity is not a preference",
   LP certificates of no-nonnegative-decomposition (:594,1024), the μ₁ positivity
   frontier (`registry.py:640,673`).
6. **Flux-positivity assertions in tests** — `[M]`
   `grep -rnE "(psi|phi|flux)[a-z_.]* *(>=|>) *0" tests/` (filtered) → **12 sites**:
   1 SN production-path gate (`tests/sn/operators/test_solver_components.py:394`
   `assert np.all(phi >= 0), "Negative flux from positive source"` — single fixture),
   2 diffusion (`test_properties.py:77`, `test_solver.py:210`), 7 derivations-tier,
   plus `tests/sn/sweep/curvilinear/test_psi_half_positivity.py` — a whole
   CHARACTERIZATION module on the M-M angular recurrence's ψ̂ sign (positive on the
   production value path, −12…−133 under inconsistent seeds; no `verifies` markers
   by design — "there is no equation whose truth they establish").

Also pre-existing vocabulary near the plan's Tier-1: the SI/Krylov linearity
batteries (§4.3), and `power_iteration`'s production-rate ray normalization (§7).

---

## 4. The superposition measurement (Q5 — the crux)

### 4.1 The legal additive/multiplicative set on a flux leaf (read off `_flux_role.py` + `field.py`)

| op | result | legality |
|---|---|---|
| `ψ₂ − ψ₁` (same class/space/mesh) | sibling `…Displacement` | legal — the only mint |
| `ψ + Δψ` (sibling displacement, same mesh, same shape) | flux | legal — torsor action |
| `ψ + ψ'` (same class) | — | `TypeError` (names `affine_combination`, `ψ+Δψ`, or rate-density adds) |
| `ψ + other-class` | — | `TypeError` (Layer-1, even at equal units) |
| `affine_combination([(λᵢ,ψᵢ)])`, Σλ=1±1e-12 | flux | legal — only multi-flux blend |
| `c·ψ`, `ψ/c`, `−ψ` | flux | **legal — scalar algebra deliberately untouched** ("eigenvalue normalisation ψ/k … all survive") |
| `zeros_for_mesh` / `from_mesh(0)` | flux | legal — zero fluxes constructed freely (initial iterates, `B(0)=0` gates) |
| `Δψ ± Δψ'`, `c·Δψ` | displacement | legal — V is a vector space |

### 4.2 Does ANY production path add two flux-TYPED fields?

**No — unspellable on typed paths by the gate, and `[M]` no site attempts it**
(the 45 TypeError test lines above are the negative controls; grep of `orpheus/`
for the gate's message fragment appears only at its definition). The three places
V-space linear algebra genuinely happens route AROUND the flux type:

- **Krylov/GMRES**: scipy iterates FLAT ndarrays (`to_flat`/`from_flat`,
  `full_field.py:410+`); roles are erased at the boundary and templates rebuild
  flux-typed (dsa.py docstring; issue #289's parse inventory is the same fact seen
  from pyright). The Arnoldi linear combinations — sums of V-vectors — happen
  outside the type system entirely.
- **The SI update**: spelled as torsor (`psi = psi + corrector.apply(psi − psi_prev)`),
  never as a flux sum.
- **Multi-source combination is at the SOURCE level, always**: the within-group rhs
  is `rhs = q_ext + Σ gᵢ.apply(ψ)` (`iteration.py:757-762`) — source/sink adds
  (vector space, closed); fission enters as `compute_fission_source(...)` handed to
  `solve_fixed_source` as q; prescribed inflow rides `q.boundary` (the affine-split
  convention — memory: no affine operator, affine law = linear op + typed source).

**"flux₁ + flux₂ of independent sources in the same medium" — zero occurrences in
production, tests, or examples as a typed add.** The nearest thing the tree does:
the double-delivery gate reads `γ₋ = 2.5 + 2.5 = 5.0` at the VALUE level
(`tests/sn/solve/test_declared_inflow_reaches_the_rhs.py:330-400`), and superposition
of PROBLEMS is exercised as solve-level linearity in q (below) — never by summing
two solved fluxes. So at HEAD, superposition is used exclusively at the source
level; no consumer is currently blocked wanting a flux sum. (What IS blocked is the
operator-on-V question — that is #331, §8.)

### 4.3 How the tree spells superposition when it needs it (verification surface)

- `tests/sn/operators/test_declared_law_is_linear.py` — additivity WITHOUT naming
  `B_V`: `B.apply(ψ₁+σ) − B.apply(ψ₂+σ) == B.apply(ψ₁) − B.apply(ψ₂)` with σ a
  displacement (:505-512; same for the full matvec :572-581). Operator OUTPUTS are
  source-typed so their `−` is legal; base-points shift by ⊕σ. The module header
  (:52-56) says this is "the affine-map law stated without ever naming B_V" and
  cites #331. The battery also includes `test_B_vanishes_at_zero` (:439) and
  `test_B_is_homogeneous` (:463) — vanishing-at-zero and homogeneity asserted on
  flux-typed inputs (both spellable because zero construction and scalar scaling
  are legal).
- `docs/theory/foundations/boundary_conditions.rst:2529-2540` — a superposition
  gate ("flux affine in the inflow amplitude") was **considered and rejected** for
  the delivery-count claim as a Mode-12 non-catcher (doubling is q→2q, still
  affine); only the coefficient against an independent reference has teeth.
- `tests/sn/sweep/core/test_ordinate_scan.py:379` — joint linearity of the scan
  kernel `scan(a,b₁+b₂,p₁+p₂) == scan(a,b₁,p₁)+scan(a,b₂,p₂)` on bare ndarrays.

### 4.4 Internal tensions of the shipped algebra (measured, not ruled)

1. **Scalar scaling is legal on flux** — a literal affine space admits no scalar
   action, and the doctrine keeps it "untouched" for ψ/k. So the shipped object is
   precisely: a vector space with binary same-class `+` disabled and `−` retyped.
2. **Zero fluxes are constructed everywhere** (initial iterates, `from_mesh` zeros,
   `B(0)=0` gates), while the doc says "the zero flux is a *chosen base point*, not
   an identity" (`field_algebra.rst` §Why-affine).
3. **DSA's own docstring uses origin language**: "the swept vector IS the
   displacement from zero" (dsa.py:622-624).
4. **`power_iteration` normalizes iterates to unit production rate**
   (`eigenvalue.py:438-451`, `flux_distribution / p`, convention ∫νΣ_f φ dV = 1) —
   operationally the plan's "ray in P(K), magnitude fixed by power normalization",
   already shipped as scalar scaling.
5. The cone ontology's needs vs the legal set: differences (have, as a distinct
   type), positive scaling (have), superposition of fluxes (unspellable — but §4.2:
   no production consumer wants it today; the demand is in the OPERATOR domain, #331).

---

## 5. The DD/step witness check (Q6)

`[M]` `grep -rni "step_characteristic\|class Step\|\"step\"\|'step'" orpheus/` →
**no Step and no step-characteristics realization exists.** The registry
(`grep 'key="' orpheus/transport/spatial/`) holds exactly TWO schemes:
`diamond_difference` and `linear_discontinuous` — both `is_positivity_preserving = False`.
`class Step(DiscretizationSchemeBase, key="step")` appears ONLY as a docstring
example of how registration works (`scheme.py:710-717`), plus prose claims that
Step *would be* positivity-preserving (:448, :701-702) and
`test_diamond.py:806-809` "Wave C-extension's :class:`Step` and
:class:`ExponentialCharacteristic` strategies will be positivity-preserving by
construction" (future-tense, unbuilt).

⟹ **The plan's C1 gate clause "DD documented non-cone-preserving … with step/SC
`True`" has no step/SC witness at HEAD** — a plan-authoring §6c gate-without-witness
defect unless C1 is fused with (or gated on) a step that lands a Step/SC
realization. What C1 does NOT need to build: the flag axis itself (exists on the
Protocol + ABC), DD's honest `False` (exists), and DD's negative-flux numerical
witness (`test_diamond.py:793-830` exists). What has no witness: any scheme for
which the flag is `True` — today the only `True` values live on test doubles
(`IdentityDiscretizationScheme` in `test_discretization_scheme_protocol.py:266`).

---

## 6. Doc surface (Q3)

**`docs/theory/foundations/field_algebra.rst`** — 602 lines, 9 sections (map:
:57 the affine-typed algebra + Key Facts; :112 two dimensional universes;
:164 the affine structure — eq `affine-torsor-algebra` :176; :242 **Why affine —
the design rationale** (the section correction #15 overturns: "flux states have no
natural zero (the zero flux is a *chosen base point*, not an identity)"; the
three-frames table affine/Banach/Krylov; the deliberate ≥2-consumer exception);
:313 displacement leaves; :382 convergence diagnostics (eqs
`affine-contraction-ratio` :398, `affine-true-error` :417); :467 typed residual
(eq `affine-typed-residual-eq` :475); :542 numerical evidence). Extracted verbatim
from `operator_algebra.rst` at #231 Phase 3; machine header marked PROVISIONAL.

**Graph view** (`nexus context std:file:theory/foundations/field_algebra`, degree
50): contains 4 equation nodes + 3 section nodes; `documents` **34 code objects**
(25 shown + 9 omitted — FunctionSpace, the flux leaves, Field, SourceIteration,
FluxRole + its dunders/affine_combination, Displacement + diagnostics, dsa module,
residual types, `_l2_norm`, `evaluate_residual`…). All four equation labels are
`vv-status: documented` entries in `docs/theory/verification/matrix.rst:873-877`
(within the 539-sentinel population), i.e. deliberately test-deferred.

**Citers**: `[M]` grep across docs → `operator_algebra.rst` (**11 ref sites**,
incl. :117, :136, :3242, :3400-3402, :3534, :3832-3861, :3955-3965, :4017, :4035,
:5824 — the Role-fiber reading, the "Flux earns a class because of the torsor"
argument, the phantom-type rejection); `wavefront_cochain.rst:35`;
`foundations/index.rst`; `methods/sn/acceleration.rst:887` (DSA update is the
torsor action); `methods/sn/curvilinear_one_group.rst:4170` (composite torsor);
`verification/matrix.rst:875`.

**`orpheus/numerics/coupled_system.py` closing prose** (the plan's cited site —
it is the module docstring's closing paragraph, ~:108, not the file end): "Role
semantics stay on the members: the machinery never adds arrays — every `+` it
evaluates is a member `+`, so **the affine flux torsor, displacement minting, and
units law of the member family apply unchanged** (the same delegation discipline
as `Composite._map_binary`)." Restated at `CoupledField` (:239) and
`_check_partner` (:292). Also `numerics/vector.py:55-58` ("the #208 affine torsor
gate") and `transport/__init__.py:31-32` ("the difference vector space V the
fields form a torsor over").

---

## 7. Phase-C hooks (Q7)

**(a) Power-iteration history.** `power_iteration` (`numerics/eigenvalue.py:374-500`)
retains per outer iterate: `keff_history: list[float]` and a dict of named
`StoppingCriterion` trajectories, accumulated per reading name via
`solver.measure_stopping_criteria(keff, keff_old, flux_distribution, flux_old)`
(SN impl at `sn/solver.py:1768`) — **the named-criteria mechanism is extensible
by construction** ("a solver may report any number of them and the loop never has
to know which"). Per-iterate FLUX is not retained, but at the loop top BOTH
`fission_source` (= Fφ, computed :420) and `flux_distribution` (= φ) are live —
the Collatz–Wielandt bracket inf/sup (Fφ)/φ is computable per iterate and could
ride the existing named-trajectory channel (two floats per outer) or the
`IterationRecord.criteria` tuple (`numerics/convergence.py:418/910`). No new
history object needed. Caveat: `compute_fission_source` uses the PREVIOUS iterate's
φ against the k of the same step — the bracket wants F applied to the CURRENT
iterate, so placement inside the loop (top of iteration n+1) matters.
Iterate normalization to unit production rate (:438-451) is already the
scale-invariant footing brackets want.

**(b) Unaccelerated SI is the DEFAULT, distinctly selectable.**
`solve_sn_fixed_source(acceleration: str | None = None)` (`sn/solver.py:3461,3645-3656`)
— `None` leaves both inner paths byte-untouched; `"dsa"` opt-in builds the ONE
`DSACorrection` consumed as `SourceIteration(corrector=)` (SI posture) or wrapped
around the swept vector (Krylov, :790). `inner_solver` defaults to
`"source_iteration"` (:2341,1205). So SI-MONOTONE's required substrate (a
selectable unaccelerated SI) exists as the default path; DSA exemption = the
already-existing `corrector is None` branch.

**(c) #365 / #366 — both convergence-verdict issues, adjacent not blocking.**
#365: `fully_converged`/`Solution.converged()` has 100 % recall / 33 % precision
(`[M]` 17/17, 17/51 at `0f5ca91c`) — the verdict a cone battery would sit beside
is currently over-conservative. #366: two public per-level `converged` surfaces
one attribute hop from the tree verdict (`IterationHistory.converged`,
`PowerIterationOutcome.converged` — the latter with zero production readers).
Relevant to WHERE bracket/monotonicity verdicts should surface (the record tree),
not to whether they can attach.

---

## 8. Issue map (Q8)

| issue | content (measured) | relation to the ontology fork |
|---|---|---|
| **#331** (OPEN) — "SN operator leaves disagree on whether their domain includes the DISPLACEMENT space" | The three leaves of ONE sum `A = (L+C) − S − B`: `L.apply(displacement)` OK → AngularSourceSink; `S.apply` TypeError (input-type enumeration excludes AngularDisplacement); `B.apply` TypeError (boundary must be AngularBoundaryFlux). So `build_within_group_system(...).loss` refuses too. A Krylov method iterates on V; today V-consumers must drop to `to_flat`. Not a correctness bug — every production path applies operators to fluxes. Asks for a user ruling: (1) every operator carries its tangent map `B_V` (the arm L already has), or (2) operators act on 𝔸 only and L's arm is the outlier. | **The sharpest statement of the tension.** Under the torsor, option 1 doubles every operator's domain contract (A-arm + V-arm); under the cone ontology (flux and difference = one tensor type in V, predicate apart) the disagreement is unspellable — `S.apply(Δψ)` type-checks by construction and the A/V ceremony dissolves. Its P5 gates already had to spell additivity without `B_V` (§4.3). Either ruling on #331 IS a partial ruling on the ontology. |
| **#289** (OPEN) — FullField generic over leaf types | The composite's `bulk: BulkField` erases the role narrowing; ~8 production sites carry local `isinstance` parses (solver drivers, boundary op, solution, fission). Candidate: `FullField(Generic[B, D])` (~14 files). Interlocks with #288; notes the driver conflation `SourceIteration[V]` (operator carrier) vs iterate carrier. | The parse inventory is role-multiplicity cost made visible: parses discriminate {Flux, SourceSink, Residual, Displacement} × window states. A cone overturn that merges Flux/Displacement into one type shrinks the discriminated set (flux-vs-displacement parses go away; flux-vs-source parses stay); a generics carve done FIRST would bake the 4-role grid into type parameters that an ontology change then reshapes. Ordering dependency, both directions. |
| **#288** (OPEN) — cross-class source-sink dunder statically unspellable | `ScalarSourceSink + AngularSourceSink → AngularSourceSink` (the #207 containment exception) cannot be typed under LSP; three candidate resolutions (named-composition rewire / NotImplemented protocol / accept-untyped). | Orthogonal to flux-vs-cone: it lives entirely in the RATE-DENSITY (vector-space) family, which both ontologies keep additive. Relevant only as precedent: it is the existing worked example of "a principled algebra the dunder type system cannot carry", the same class of cost any predicate-carrying flux type would meet at `__add__`. |

---

## 9. What a cone overturn would actually touch (verdict-free enumeration)

Counts carry their predicate; line numbers re-derivable.

**Production code — the torsor grep floor** (`[M]` `grep -rl "torsor" orpheus/ --include="*.py"` → 16 files):
`transport/fields/_flux_role.py` (the mixin — the object at stake),
`transport/fields/_coefficient_role.py` (defines itself as the complement),
`transport/fields/cross_section_field.py`, `transport/fields/radial_characteristic_interior_flux.py`,
`transport/displacements/` (whole package: `__init__`, `_displacement`, 2 RC leaves
carry the word; all 8 modules are the V-space realization),
`transport/full_field.py` (delegation prose + gate comments :107,:352,:358),
`transport/radial_characteristic_field.py:106`, `transport/__init__.py:32`,
`transport/source_sinks/harmonic_moment_source_sink.py:33`,
`numerics/coupled_system.py` (:108,:239,:292), `numerics/vector.py:58`,
`sn/solver.py:2960` (+ #353), `sn/acceleration/dsa.py` (:608,:680-696 — behavior,
not just prose: the minted return types).
Beyond the grep floor (no "torsor" spelling but structurally coupled): the 7
FluxRole leaf classes, `numerics/field.py` Layer-1 seam, `numerics/iteration.py`
(`_DisplacementLeaf`, `_flux_displacement_leaf`, SI diagnostics wiring),
`timed_full_field.py`, the #289 parse sites (~8), and every `…Displacement`
annotation (`grep -rn "Displacement" orpheus/` → 2 consumer modules + the package).

**Types**: 1 mixin (FluxRole) + 1 marker-with-registry (Displacement) + 7
displacement leaf classes + the sibling_of registry mechanics; plus the role-grid
docstrings on all 7 flux leaves and the residual/source packages' role-grid tables
(`transport/{displacements,residuals}/__init__.py`).

**Tests**: 16 files spell "torsor" (`[M]` grep -rl, §1.3 list); 7 files pin the
flux+flux TypeError by name; 10-test foundation battery
(`test_affine_flux_algebra.py`) asserts the algebra itself; 9 files / 36 lines
reference Displacement types; ~16 raise-sites name the affine gate;
`test_flux_displacement_diagnostics.py` + `test_psi_half_coupling.py:3632` pin the
SI diagnostics wiring; the torsor-form linearity gates
(`test_declared_law_is_linear.py`, `test_phase_c_gates.py:260`) would re-spell
trivially under either ontology (their assertions are value-level).

**Docs**: 5 pages spell "torsor" — `field_algebra.rst` (the 602-line page whose
§"Why affine" is the overturned argument; 4 equation labels + 3 section anchors,
all cross-referenced), `operator_algebra.rst` (11 ref sites incl. the
"Flux earns a class" §3854 and the phantom-type rejection §3977-4017),
`acceleration.rst:887`, `curvilinear_one_group.rst:4170`,
`verification/matrix.rst:875` (4-5 `affine-*` labels inside the 539-sentinel
documented population). Plus the `nexus-meta` role line and 34 `documents` edges
from the page to code objects. A retirement/rewrite owes the standard
four-search audit; note the equation label `affine-torsor-algebra` is `:eq:`-cited
from `operator_algebra.rst:3402` — a labelled equation is an API.

**Out-of-tree surfaces the plan should know exist** (referenced from docstrings):
`.claude/agent-memory/cross-domain-attacker/issue_208_delta_psi_affine_frames.md`,
`.claude/plans/issue_208_residual_typing_closeout.md`, the coding-elegance skill's
anti-pattern #18 (cites FluxDisplacement as a POSITIVE precedent — the skill text
itself would teach a retired ontology).

**What a cone overturn would NOT touch** (measured invariants): the stopping logic
(flat norms of source-typed residuals), the source/sink/residual vector-space
algebra and #288, the Layer-1 cross-class gate (memo A's seam), the Krylov flat
path (already role-erased), `is_positivity_preserving` and its witness test
(already cone-vocabulary), the realizer's cone refusal, quadrature-weight
positivity, and the production-rate normalization.

**Recent movement (Q9)**: `[M]` `git log --oneline --since=2026-08-08 -- orpheus/transport/fields orpheus/transport/displacements orpheus/transport/full_field.py` →
**empty**. Last commits on those paths: `b4984773` 2026-08-03 (docstring xref
sweep), `8d996f53` 2026-07-28 (xrefs), `ec116e75` 2026-07-26 (displacement tangent
map for the moment-0 reduction, #2 — the most recent structural change),
`167ac25b` 2026-07-19. The torsor layer is quiescent; no in-flight edits collide
with a campaign start.
