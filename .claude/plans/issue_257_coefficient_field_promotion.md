# CoefficientField + operator-as-promotion (#257) — implementation plan

> **POST-COMPACTION REGROUNDING (read in this order before any code):**
> 1. THIS plan (the roadmap + the embedded census facts below).
> 2. GH **issue #257** (`gh issue view 257`) — the durable thesis.
> 3. Frame memo: `.claude/agent-memory/cross-domain-attacker/coefficient_field_promotion_frames.md`
>    (the 4 structural frames + the multiplier-algebra law-suite).
> 4. The campaign anchor `.claude/plans/field_typed_operator_algebra_campaign.md` (the #256→#257 re-scope + the baseline reds).
> 5. Then the code, in the order the stages touch it.
>
> **Branch** `feature/field-typed-operator-algebra` (off `main`, NOT pushed). Foundation already
> committed: #256 step 1 (`Vector` Protocol, `cfb651b`) + step 3 (`apply(x:V)->V`, `41a92cb`).
> **On implementation start, FIRST copy this plan to `.claude/plans/issue_257_coefficient_field_promotion.md`
> (repo-durable; `~/.claude` is host-ephemeral) and commit it.**

## Context

`apply(self, x: V) -> V` over a bare `Generic[V]` was diagnosed (by the user) as unexplanatory —
"if type hinting doesn't help, it's compliance theatrics." That thread surfaced that grand report
§5.5–5.7 already specifies the target the code never reached: **cross-sections are `CoefficientField`s**
(a `Field` sibling), and **operators are those fields promoted** (`C = MultiplicationOperator(σ_t)`,
the multiplier-algebra embedding `M: L^∞ → B(L²)`). The opaque generic becomes the honest
**`TransportState`** bound. This re-scopes #256 (steps 1+3 are the foundation that stands); the rest
(scipy single-source, `Functional`, `BoundaryMomentField`) re-homes into the right taxonomy.

## The operator model (user-confirmed) — what is MODEL-SPECIFIC vs SHARED

- **Model-specific** (live in `sn/`, and each future model — CP/MoC — has its own):
  - `StreamingOperator` = **L** (pure streaming, `Ω·∇`; σ-free).
  - The **loss representation** = **L + C** = the `InvertibleOperator` / the sweep `(L+C)⁻¹`
    (`orpheus/sn/loss_representation.py`). This is where σ legitimately lives (the cell optical depth).
- **Shared across models** (live in `transport/`):
  - `CoefficientField` (+ `CrossSectionField`, `SpectrumField`), `MultiplicationOperator`,
    `C = CollisionOperator = M[σ_t]`, `ScatteringOperator`, `FissionOperator`, the `Functional` category,
    `TransportState`.

```
Field (numerics base = plain vector space; FluxRole adds the affine torsor)
├── CoefficientField (transport/; CoefficientRole = base algebra, NO affine gate)
│   ├── CrossSectionField  (cone σ≥0, units 1/cm)
│   └── SpectrumField      (probability simplex Σχ=1, dimensionless)
LinearOperator
├── MultiplicationOperator (transport/; = M[f], promotion of a CoefficientField)
│   └── CollisionOperator   (C = M[σ_t]; today ALREADY `σ[None]*ψ`)
├── StreamingOperator (sn/; = L, model-specific) ; loss rep = L+C (sn/, model-specific)
└── IntegralKernelOperator (transport/) → ScatteringOperator / FissionOperator
```

## Embedded census facts (so this plan regrounds without re-exploring)

- `CollisionOperator.apply` (`orpheus/sn/operator.py:592`) **is already** `self.sigma[None]*psi.bulk.values`;
  `.solve` (`:620`) `q/σ`; `apply_transpose`=`apply` (self-adjoint, currently untyped). The cleanest M[σ] target.
- The **C-fold** is in `StreamingOperator.apply` (`operator.py:417-418`): `loss_action(σ,ψ).bulk − σ[None]*ψ`
  — i.e. it computes (L+C) then subtracts C to *present* L. Target: StreamingOperator computes pure L
  (σ-free); the (L+C) action stays in the loss representation. `InvertibleOperator.apply` (`:884`)
  single-sources σ from `self.diagonal.sigma`.
- The **σ `np.ndarray` thread**: `LossRepresentation.loss_action(sigma, psi)` — 8 sites in
  `loss_representation.py` (`:252,272,429,442,773,933/947,1126/1147,1397/1418,1667`) + 4 callers in
  `operator.py` (`:417,461,884,903`) + `StreamingOperator.sigma_t` field (`:306`). σ in the loss-rep
  SWEEP is irreducible (cell discretization); σ in StreamingOperator (pure L) is removable.
- Primitives (`numerics/operator.py`): **only** `RankOneOperator` (`:1554`, fission, 1 prod site
  `fission.py:218`) and the σ-multiply are promotion targets. `DiagonalOperator` (`:1436`) has **NO prod
  caller** (`from_measure` `:1501` is the angular-W factory, no live caller). The 4 BC primitives
  (`Permutation`/`IncomingOrdinateMask`/`PeriodicWrap` + their `TensorProduct` folds, all in
  `boundary_realizer.py`) and `SumOfTensorProductsOperator` are basis-maps/projections → **stay flat**.
- `ScatteringOperator` already delegates to typed `mat_xs.apply_legendre_scattering_moments`/
  `apply_p0_in_scatter`/`apply_n2n` (`scattering.py`) — mostly a re-presentation, not a rewrite.
- `FissionOperator.kernel` (`fission.py:168`) = `RankOneOperator(χ,νΣf,axis=0) & IdentityOperator()`;
  `apply` is `@singledispatchmethod` over {TFF, ScalarFlux, ndarray}.
- `MaterialXSField` (`sn/material_xs_field.py`) = bare `@dataclass` (NOT a Field); one prod builder
  (`solver.py:858`). Views: `total_cross_section`/`absorption`/`fission_production`/`emission_spectrum`
  (ng,*spatial) + `sig_s_legendre`/`n2n_matrix` (group→group matrices). σ_s = a Kernel, not a Coefficient.
- `Functional` candidates: `iteration.py:251/272` estimators; the per-cell production rate is the
  anonymous `inner` einsum in `RankOneOperator.apply` (`operator.py:1671`).
- `units.py`: 4 signatures (`:130,136,140,144`); add `CROSS_SECTION_UNITS = UREG.cm**-1` (5th).
- Layering: `numerics ↛ transport` holds (zero `from orpheus.transport` in numerics) — `Vector`/`V` STAY
  in numerics; `TransportState` + coefficient leaves + `MultiplicationOperator` go in `transport/`.
- §4158 `.multiplication_operator()` is grand-plan VAPORWARE (`CriticalityEigenproblem` doesn't exist)
  — **no code rename needed**; reserve the name (the future eigen-iteration verb is `iteration_operator()`).
- `as_scipy_linop` (`operator.py:1680`): ZERO prod callers; live scipy boundary is the inline closure
  `iteration.py:744-766`.

## Staged implementation (each stage = one reviewed [elegance + qa] + gated commit)

**S1 — `CrossSectionField` + `SpectrumField` + `1/cm` unit + `CoefficientRole`.** Pure addition.
`units.py` 5th signature; `transport/fields/` `CoefficientRole` (base plain-vector-space algebra, NO
affine gate), `CrossSectionField(CoefficientRole, ScalarField)` (cone), `SpectrumField` (simplex Σχ=1,
`.mix` convex blend). **Tests of intrinsic properties (mandatory, per user):** the cone closure
(σ+σ ✓, λσ≥0 ✓, σ=0 origin), the simplex invariant (Σχ=1 enforced at construction + a failing-input
negative test), the role algebra differs from `FluxRole` (no `flux+flux→TypeError`-style gate).

**S2 — `MaterialXSField` typed accessors.** Bit-identical. Typed field accessors alongside the raw
ndarray views; `.values` bit-equal; existing consumers untouched this stage.

**S3 — `MultiplicationOperator` (transport/) + promote `CollisionOperator`.** Expect bit-identical.
`apply` = coeff broadcast-multiply on the matching subspace (einsum on `.values`, wrap once); `.solve`
= divide; `.H = self` (real coeff). Rewire `solver.py:217,925,1000`. **Law-suite tests (the math
concept's intrinsic properties):** `M_f@M_g=M_{fg}`, `M_1=I`, `M_0=ZeroOperator`, `M_{af+bg}=aM_f+bM_g`,
`M.H=M`, `spec=ess-range` (invertible iff min|f|>0). Gate vs old `CollisionOperator` (see §verification).

**S4 — `TransportState(Protocol)` (transport/, refines `Vector`).** Bit-identical annotations. Names
`.bulk`/`.boundary`; re-points the transport/SN operator-leaf annotations to read like the domain;
`numerics` keeps `Generic[V: Vector]` (L0-ndarray + scipy path). Intrinsic-property test: `np.ndarray`
is NOT a `TransportState`; `TimedFullField` IS.

**S5 — `Functional` category.** Bit-identical reductions. keff/production estimators + the per-cell
`ProductionRateFunctional` as a `Functional` `evaluate(x)->scalar/field` (§5.6 suffix; #256-step-5).
Test: a `Functional` is NOT a `LinearOperator`.

**S6 — `FissionOperator` as `IntegralKernelOperator` = `M_χ ∘ ProductionRateFunctional ∘ M_νΣf`.**
⚠ **May be principled-not-bit-identical** (the composition changes the reduction order vs the fused
einsum) — judge by the 3-criteria, keep `RankOneOperator` as the structurally-independent cross-check.

**S7 — scipy single-source.** Bit-identical. Route the inline closure (`iteration.py:744-766`) through
one adapter; retire/fold the zero-caller `as_scipy_linop`. #256-step-4. Gate: keff unchanged to solver
tol; `from_flat(to_flat(x))==x`.

**S8 — restore `StreamingOperator` to pure L; loss rep = L+C. BEHAVIORAL.** Drop the `(L+C)−C` fold:
StreamingOperator computes pure streaming (σ-free); `C` is the shared `MultiplicationOperator`; the
loss representation (model-specific) is `L+C` (the sweep, which keeps σ in its cell discretization).
Retire σ from `StreamingOperator`'s surface; σ stays in the loss-rep sweep. **Dispatch test-architect
FIRST** (operator-composition seam). Gate: SI rhs / Krylov matvec value-equal (the composition graph
changes — verify the field-level `L + C − S − B`); `InvertibleOperator` σ single-source consistent.

**S9 — `BoundaryMomentField` + close the moment-state boundary drop. BEHAVIORAL.** #256-step-6, fork
#1 = general moment-tail (the #251 `boundary_face_layout` lever). The `(moment_buf, None)` drop is in
`apply_windowed`, not public `solve_moments`. **Dispatch test-architect FIRST** (scalar↔angular↔moment).
Gate: moment tensor byte-identical; new boundary block provably == old `None`.

## Verification (Cardinal Rule 1)

- **Per stage, judge the gate principledly — NOT a blanket 0 ULP.** For each touched output ask the
  vv-principles 3 criteria: (a) is every intermediate a NAMED principled quantity? (b) is it verified
  against a STRUCTURALLY-INDEPENDENT reference (the multiplier-algebra law, k_∞=νΣf/Σa, the
  RankOne cross-check, mpmath), not just old-vs-new? (c) is any drift FP-non-associativity bounded
  (`reduction-depth × ULP`)? 0 ULP is the *expectation* for pure re-typing/wrapping (S1–S5,S7); a
  principled reduction-order change (S6) may legitimately differ — accept ONLY if all 3 hold, then
  narrow the regression contract (`assert_array_almost_equal_nulp`) for that output with the documented
  justification. Be prepared for surprises and classify them, don't force-fit 0 ULP.
- **Every mathematical concept gets a test of its intrinsic properties** (user directive): the cone /
  simplex (S1), the multiplier-algebra laws (S3), the Functional-≠-Operator (S5), the
  `TransportState`-≠-ndarray (S4). These are the L0/L1 gates, not afterthoughts.
- Regression subset (route around the 7 baseline reds #250/#232 + #212): `.venv/bin/python -O -m pytest
  -q tests/sn/operators tests/sn/spatial tests/sn/sweep/core tests/sn/solve tests/numerics
  --deselect tests/sn/solve/test_keff_slab.py::test_heterogeneous_absolute_keff`.
- CLI `npx pyright` per stage: the opaque-generic surface shrinks; **zero new `# type: ignore`**.
- Sphinx clean; archivist documents the CoefficientField/promotion algebra after S6/S9.

## Execution discipline

Commit per verified stage (no push without an ask); `Co-Authored-By: Claude Opus 4.8 (1M context)
<noreply@anthropic.com>`; L28 (never `git checkout/restore/stash` on tracked paths); explicit paths;
never `.claude/skills/*`/`CLAUDE.md`/`.claude/rules`/`.claude/hooks`. First implementation action:
copy this plan into the repo + update issue #257 with the staged roadmap (durability before compaction).
