---
name: nexus82-operator-algebra-implementers
description: Ground truth for `.. implements::` on the 40 labelled equations of docs/theory/foundations/operator_algebra.rst (nexus #82) — 32 declarable, 8 NONE, grouped by statement KIND, plus three doc-drift findings the declarations expose.
metadata:
  type: project
---

# nexus #82 — who implements `theory/foundations/operator_algebra`'s 40 equations

Determined 2026-08-17 against HEAD (graph build `8bb695d3`, `changed_since_build: false`).
Every dotted path below was existence-checked against `.nexus/graph.db` — **55 of 55
node IDs resolve**. Line numbers are re-derive-via-Nexus; the verdicts are durable.

**Why:** nexus currently mints 782 *inferred* `implements` edges into this page's 40
labels by shared name token. Declaring any implementer of an equation stands the
guessing down for that whole equation — so a wrong declaration is worse than a guess
(it ships at `confidence=1.0`).

**Headline: 32 of 40 declarable, 8 NONE.**

⛔ **CORRECTED 2026-08-17.** This line read *"21 of 40 declarable, 19 NONE"* — wrong
against **this file's own table**, whose kind breakdown sums to 8. The error was
inherited by the archivist's brief and by two GitHub comments before an archivist
re-counted while building to the table. ⟹ a headline is a SUMMARY of a table and
rots independently of it; re-derive it from the rows rather than carrying it.

**Headline: 32 of 40 declarable, 8 NONE.** The NONE set is not a gap — it is the
finding (§3): a fifth of this page's equations are *statements about* the algebra
(law / definition / identity / canonical-form) rather than *computations in* it.
`implements` is meaningless for that class, and the inference engine cannot tell
them apart from a label.

⛔ The ratio here read *"3 : 2"* under the wrong 21/19 headline. `[M]` it is
**8 of 40** — and the *typing rule* was wrongly listed among the non-implementable
kinds: `carrier-grid-operator-typing`, `harmonic-frame-is-galerkin` and
`product-solve-reroute` all read as identities and all have real carriers. The
usable split is **identity between QUANTITIES → no implementer; identity between
TYPES → a class declaration IS the implementer.**

---

## 1. The declaration table

| label | implementer (full dotted path) | file:line | one-line justification |
|---|---|---|---|
| `operator-apply` | `orpheus.numerics.operator.LinearOperator.apply` | `orpheus/numerics/operator.py:740` | The Protocol contract stub whose docstring IS the equation ("Return :math:`L\,x`"); declaration-site reading (§2.1). |
| `operator-solve` | **NONE** — definition (a verb name) | — | The page's own base-hosting rule forbids a declaration site; `[M]` `LinearOperator` declares only `apply`, and the narrowing Protocols declare `inverse()` / `apply_transpose` / `assemble` — never `solve`. |
| `operator-apply-transpose` | `orpheus.numerics.operator.SupportsAdjoint.apply_transpose` | `orpheus/numerics/operator.py:1168` | The one algebra-level declaration of the verb ("the raw Euclidean transpose verb"). Same reading as `operator-apply`. |
| `operator-fixed-source` | `orpheus.sn.solver.SNSolver._solve_source_iteration`<br>`orpheus.sn.solver.SNSolver._solve_krylov` | `orpheus/sn/solver.py:1798`<br>`orpheus/sn/solver.py:2024` | The two within-group realizations that POSE `Aψ = q`: each binds `A = L+C`, gains `(S, B)`, `F = 0_wg`, and puts the source on the right. If the posing were wrong, both assemble the wrong problem directly. |
| `operator-eigenvalue` | `orpheus.sn.solver.SNSolver.compute_fission_source`<br>`orpheus.numerics.iteration.KEigenvalue.compute_fission_source` | `orpheus/sn/solver.py:1424`<br>`orpheus/numerics/iteration.py:1403` | `return self.fission_op.apply(φ) / keff` — literally the `(1/k)Fψ` right-hand side, which is the *only* thing distinguishing this equation from `operator-fixed-source` ("they differ only in what sits on the right"). |
| `operator-within-group-composition` | `orpheus.sn.coupled_system.build_within_group_system` | `orpheus/sn/coupled_system.py:425` | `A_AA = LC - S - B_a` — literally `A = L+C−S−B`, and its own docstring calls itself the single production spelling. |
| `diagonal-operator-action` | `orpheus.numerics.operator.DiagonalOperator.apply` | `orpheus/numerics/operator.py:3352` | `self._broadcast(x.ndim) * x` — the elementwise multiply along the tagged axis. |
| `multiplication-operator-action` | `orpheus.transport.operators.multiplication_operator.MultiplicationOperator._apply_impl` | `orpheus/transport/operators/multiplication_operator.py:377` | Both registered arms compute `f_{g,r}·ψ_{n,g,r}` (composite → the `DiagonalOperator` broadcast engine; bare-array → `coefficient.values * phi`). ⚠ target `_apply_impl`, not `apply` — the public name is a runtime alias that exists only under `TYPE_CHECKING` and has no graph node. |
| `multiplication-operator-embedding` | `orpheus.transport.operators.multiplication_operator.MultiplicationOperator` | `orpheus/transport/operators/multiplication_operator.py:130` | The class IS the map `f ↦ M[f]`: it stores **only** a `CrossSectionField`, so "the cross section is the operator" is the type, not a comment. (The class = the embedding; its `_apply_impl` = the discrete action, previous row. No overlap.) |
| `inverse-as-operator` | `orpheus.numerics.operator.InverseOperator.apply`<br>`orpheus.sn.operators.sweep_operator.SweepOperator.apply` | `orpheus/numerics/operator.py:2269`<br>`orpheus/sn/operators/sweep_operator.py:~120` | Both bodies ARE the keystone `A.inverse().apply(b) == A.solve(b)`: `InverseOperator.apply` is `return self.inner.solve(x)`; `SweepOperator.apply` delegates to `inner.solve` bit-identically (the `(L+C)` instance the equation is written for). ⛔ do NOT add `GreenOperator.apply` — its forward `OperatorSum` has no `solve`, so it is the case the identity does not cover. |
| `carrier-grid-operator-typing` | `orpheus.numerics.operator.LinearOperator`<br>`orpheus.numerics.operator.Domain`<br>`orpheus.numerics.operator.Codomain` | `orpheus/numerics/operator.py:588`<br>`:117`<br>`:118` | The equation is a typing rule and `class LinearOperator(Protocol[Domain, Codomain])` with `apply(x: Domain) -> Codomain` is that rule written in code; the two TypeVars (PEP-696 `default=Domain`) are the grid coordinates it names. |
| `apply-solve-source-iteration-series` | `orpheus.numerics.iteration.SourceIteration.solve` | `orpheus/numerics/iteration.py:660` | The loop `rhs = q_ext + Σ g.apply(ψ)` → `ψ = A_inv.apply(rhs)` IS the partial sums of `Σ_k[(L+C)⁻¹S]^k(L+C)⁻¹q`. |
| `apply-solve-within-group-balance` | `orpheus.sn.operators.streaming.StreamingCollisionOperator` | `orpheus/sn/operators/streaming.py:445` | The class IS the discretised `Ω·∇ + Σ_t`: its `apply` computes the balance's LHS, its `solve` inverts it against `q`. |
| `apply-solve-cell-resolvent` | `orpheus.transport.spatial.diamond.DiamondDifference.update`<br>`…DiamondDifference.cell_kernel_batch`<br>`…DiamondDifference.affine_scan_coefficients`<br>`…DiamondDifference.cartesian_scan_coefficients` | `orpheus/transport/spatial/diamond.py:177`<br>`:388`<br>`:572`<br>`:694` | The four DD realizations of the **divide by the summed denominator**: `(source + numer_upstream)/denom`, `numer/denom`, and (scan form) `inverse_denom = 1.0/denom` twice. Discriminator vs the row below: this equation is the *division*; that one is the *diagonal assembly*. |
| `streaming-action-cell-balance` | `orpheus.transport.spatial.cell_balance.cell_balance_for_streaming`<br>`orpheus.transport.spatial.cell_balance.cell_balance_terms`<br>`…DiamondDifference._cartesian_streaming_diagonal`<br>`…DiamondDifference.affine_scan_coefficients`<br>`…DiamondDifference.residual_kernel_batch` | `orpheus/transport/spatial/cell_balance.py:123`<br>`:268`<br>`orpheus/transport/spatial/diamond.py:326`<br>`:572`<br>`:444` | The five sites that assemble `S = S_stream + σ_t V` and its upstream numerator (`denom, numer_upstream` / `denom = Σ_t + Σ_a 2g_a`), plus the residual form `r = Sψ̄ − (Q + Σ 2g_a ψ_in)`. |
| `harmonic-frame-is-galerkin` | `orpheus.transport.frames.harmonic_frame.HarmonicFrame`<br>`…HarmonicFrame.from_galerkin` | `orpheus/transport/frames/harmonic_frame.py:83`<br>`:110` | `class HarmonicFrame(GalerkinFrame)` IS the Liskov claim; `from_galerkin` (`cls(frame.basis, frame.measure)`) is what makes it hold with zero rebuild, which is the equation's operative content. |
| `streaming-action-pure-l` | `orpheus.sn.loss_representation._LossRepresentation.streaming_action` | `orpheus/sn/loss_representation/__init__.py:485` | One line: `return self.loss_action(self._zero_sigma_for(psi), psi)` — literally the equation's right half `streaming_action(ψ) = loss_action(0, ψ)`. ⛔ `StreamingOperator.apply` is a pure delegation to it — a consumer, not an implementer. |
| `product-solve-reroute` | `orpheus.numerics.operator.OperatorProduct.solve` | `orpheus/numerics/operator.py:1665` | `return self.b.inverse().apply(self.a.inverse().apply(b_vec))` — the equation transcribed. |
| `tensor-product-space-agreement` | `orpheus.numerics.operator._agreed_space`<br>`…TensorProductOperator.domain`<br>`…TensorProductOperator.codomain` | `orpheus/numerics/operator.py:338`<br>`:2872`<br>`:2887` | `_agreed_space` IS the three-way law (agree / silent→`None` / disagree→`IncompatibleOperatorComposition`), written once and shared with `OperatorSum`; the two properties are the equation's named subject. |
| `scattering-as-tensor-product-sum` | `orpheus.transport.operators.scattering.LegendreMomentScattering`<br>`…LegendreMomentScattering.apply` | `orpheus/transport/operators/scattering.py:115`<br>`:183` | The class docstring self-identifies with this label — "the diagonal spectrum of the sum-of-tensor-products form `Λ = Σ_ℓ P_ℓ ⊗ Σ_{s,ℓ}`" — and `apply` computes the per-ℓ block-diagonal group transfer. ⚠ see §4 finding C: the equation writes `S`, and the *full* per-ordinate `S` is explicitly NOT this shape. |
| `production-rate-functional` | `orpheus.transport.reaction_rate_functional.ReactionRateFunctional` | `orpheus/transport/reaction_rate_functional.py:76` | The typed co-vector `⟨Σ_x, ·⟩` (weight = the cross-section values, `axis=0` = groups). ⛔ NOT `InnerProductFunctional.evaluate`: that is the *generic* `⟨w,·⟩` and would stay correct if this equation were wrong. |
| `apply-distributes` | `orpheus.numerics.operator.OperatorSum.apply` | `orpheus/numerics/operator.py:1473` | `return self.a.apply(x) + self.b.apply(x)` — "applying a sum of linear maps is the sum of the applications" as a body, not a claim. ⚠ the *most contestable* declaration on the page; see §2.2 and §4 finding D. |
| `scattering-carrier-grid` | `orpheus.transport.frames.harmonic_frame.HarmonicFrame.analyse` (M)<br>`orpheus.transport.operators.scattering.LegendreMomentScattering.apply` (Λ)<br>`orpheus.transport.frames.harmonic_frame.HarmonicFrame.reconstruct` (R) | `orpheus/transport/frames/harmonic_frame.py:147`<br>`orpheus/transport/operators/scattering.py:183`<br>`orpheus/transport/frames/harmonic_frame.py:184` | The diagram's three edges are three materialized methods carrying exactly the typed endpoints the diagram draws (role-preserving M/R via `@overload`, role-CHANGING Λ). The four nodes are types, not code. |
| `scattering-aniso-composite` | `orpheus.transport.operators.scattering.ScatteringOperator.build_aniso_source`<br>`…ScatteringOperator.kernel`<br>`…ScatteringOperator._apply_impl` | `orpheus/transport/operators/scattering.py:910`<br>`:598`<br>`:1093` | `build_aniso_source` is literally `self.kernel.apply(ψ.values) / sum_w` = `(1/W)·(R∘Λ∘M)`; `kernel` is the `frame.conjugate(Λ)` composite (pre-`1/W`); `_apply_impl`'s `HarmonicMomentFlux` arm is the **second, deliberately-kept** explicit-typed realization `frame.reconstruct(Λ.apply(φ))/W`. The page states both arms are production by design. |
| `reaction-rate-kinf-oracle` | `orpheus.derivations.common.eigenvalue._infinite_medium_matrices`<br>`…eigenvalue.kinf_homogeneous`<br>`…eigenvalue.kinf_and_spectrum_homogeneous` | `orpheus/derivations/common/eigenvalue.py:39`<br>`:67`<br>`:98` | The `(A, F)` assembly and the `λ_max(A⁻¹F)` eig — in `orpheus/`, not `tests/` (the test *imports* it). ⚠ see §4 finding B: the page writes `A = diag(Σ_a) − Σ_s`; the code builds `A = diag(Σ_t) − (Σ_s + 2Σ_2)ᵀ`. |
| `fission-as-dyad` | `orpheus.transport.operators.fission.FissionOperator.kernel`<br>`…FissionOperator.production_rate` | `orpheus/transport/operators/fission.py:257`<br>`:313` | `return outer(self.chi, self.production_rate) & IdentityOperator()` — the equation transcribed; `production_rate` is the row co-vector `⟨νΣf|`. `[M]` verified the matvec really routes through it: both the `ScalarFlux` and `ndarray` arms call `self.kernel.apply(...)`, whose `RankOneOperator.apply` is `recon * functional.evaluate(x)` — no procedural twin. |
| `trace-half-decomposition` | `orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.inflow_indices_for_face` (Γ₋)<br>`…AngularTraceSpace.outflow_indices_for_face` (Γ₊) | `orpheus/numerics/spaces/angular_trace_space.py:474`<br>`:483` | The two directional selectors that realize the sign split — and, exactly as the equation's surrounding prose says, they are disjoint but NOT exhaustive (the tangential band is in neither). |
| `per-face-inflow-mask` | `orpheus.numerics.spaces.angular_trace_space.build_omega_dot_n`<br>`…AngularTraceSpace.inflow_indices_for_face` | `orpheus/numerics/spaces/angular_trace_space.py:224`<br>`:474` | `build_omega_dot_n` builds the `(n_faces, n_ordinates)` `Ω·n̂_f` table; `inflow_indices_for_face` applies `row < -TANGENTIAL_EPS` — the equation's `< −ε` predicate, with `TANGENTIAL_EPS = 4·eps` at `:186`. |
| `keff-as-integrated-rates` | `orpheus.sn.solver.SNSolver.compute_keff` | `orpheus/sn/solver.py:1552` | The only site computing `k` as a ratio of `IntegratedReactionRate`s. ⛔ **BLOCKED — fix the equation first**; see §4 finding A: the labelled equation is present-tense-FALSE against this body. |
| `integral-kernel-category` | `orpheus.transport.operators.integral_kernel_operator.IntegralKernelOperator.kernel` | `orpheus/transport/operators/integral_kernel_operator.py:195` | The Protocol member declaration whose docstring is verbatim this equation; declaration-site reading, same as `operator-apply`. The two concrete `kernel`s carry their own labels. |
| `tensor-product-adjoint-distributivity` | `orpheus.numerics.operator.TensorProductOperator.apply_transpose` | `orpheus/numerics/operator.py:2912` | The body loops `op.apply_transpose` over the factors — `(A⊗B)ᵀ = Aᵀ⊗Bᵀ` executed, not asserted. |
| `tensor-product-inverse` | `orpheus.numerics.operator.TensorProductOperator.inverse` | `orpheus/numerics/operator.py:2937` | `return TensorProductOperator(tuple(factor_inverses))` — the law constructed. |
| `tensor-product-action` | `orpheus.numerics.operator.TensorProductOperator.apply` | `orpheus/numerics/operator.py:2906` | The sequential per-factor loop `for op in self.ops: out = op.apply(out)`. |
| `apply-solve-parallel-identity` | **NONE** — identity (counter-example) | — | Unspellable in production: `StreamingOperator` declares no `inverse()`/`solve` at all, so `L⁻¹` does not exist as an object. Stated only to say what is *not* the coupled inverse. |
| `apply-solve-neumann-series` | **NONE** — identity (Lewis & Miller §3.2) | — | The splitting around **C** with `C⁻¹L` is never run. `GreenOperator` is a Neumann/splitting iteration, but around the sum's **leading term** (`A − B`), and production always spells `L + C` and promotes to the direct-sweep `StreamingCollisionOperator`. Declaring `GreenOperator` here would be a wrong declaration. |
| `apply-solve-neumann-expansion` | **NONE** — identity | — | The term-by-term expansion of the row above; same reason. |
| `apply-solve-denominator-inequality` | **NONE** — identity (an inequality) | — | Nothing computes a non-equality. |
| `solve-does-not-distribute` | **NONE** — law (a negative statement) | — | ENFORCED structurally (`OperatorSum` carries no `solve` verb; `StreamingOperator` declares no `inverse()`), but an absence is not an implementation. |
| `streaming-as-tensor-product-sum` | **NONE** — canonical form (not realized) | — | `[M]` `SumOfTensorProductsOperator` has **zero** production consumers (grep `orpheus/`: two `numerics/__init__.py` export lines + one `functional.py` docstring). Streaming is realized as a sequential walk; the page's own `operator_tensor_network` cross-reference records that it resists a clean factorisation. |
| `carrier-grid-cell` | **NONE** — definition (a taxonomy) | — | `Carrier = (Representation, Role)` neither computes a quantity nor performs an operation. Its code counterpart is the 16-leaf flat-MI hierarchy — a *structure*. Related-not-implementing: `orpheus.transport.fields._bases` (Representation ABCs) and `orpheus.transport.fields._flux_role.FluxRole` / `orpheus.transport.displacements._displacement.Displacement` (Role mixins). |

---

## 2. The two consequential judgement calls, and the readings taken

### 2.1 The three verb-definition labels — the **declaration-site** reading

`operator-apply` / `operator-solve` / `operator-apply-transpose` name the algebra's three
primitive verbs. Three readings were available:

| reading | consequence |
|---|---|
| (a) every concrete override | ~50 symbols each; reproduces the guessing it replaces, and is *wrong* — `DiagonalOperator.apply` implements `diagonal-operator-action`, not the verb |
| (b) **the declaration site, if one exists** | ← **taken** |
| (c) NONE for all three | loses the one symbol whose whole purpose is to BE the equation |

The reading is uniform and **falsifiable from the page's own text**: §"The base-hosting
rule" says a method lives on the base iff a universal realization exists, and that
`inverse()` / `solve` / `apply_transpose` are NOT base-hosted. So the asymmetry in the
answers is *predicted*, not arbitrary:

- `apply` — mandatory, base-hosted ⟹ `LinearOperator.apply` ✅
- `apply_transpose` — has a narrowing-target Protocol ⟹ `SupportsAdjoint.apply_transpose` ✅
- `solve` — has neither (`[M]` no `SupportsSolve` exists; the inverse axis's Protocol
  declares `inverse()`, not `solve`) ⟹ **NONE** ✅

Consequence for each concrete verb: it implements *its own* equation where the page has
one (`diagonal-operator-action`, `multiplication-operator-action`, `streaming-action-pure-l`,
`product-solve-reroute`, `apply-solve-cell-resolvent`, `inverse-as-operator`,
`tensor-product-action`), and nothing here otherwise.

### 2.2 `apply-distributes` — declared, against expectation

The brief predicted NONE, grouping it with `solve-does-not-distribute`. They are **not**
the same kind of statement: the second is a non-identity (nothing computes a `≠`), the
first is the **definition of operator addition's action**, and `OperatorSum.apply`'s body
is literally its right-hand side. Applying the brief's own test — *if this equation were
wrong, would this function compute the wrong thing directly?* — `OperatorSum.apply` fails
directly. Declared.

⚠ The nuance that makes it interesting is in §4 finding D: on the *specific* `(L+C)`
instance the equation is written for, the tree does **not** execute `Lψ + Cψ`.

### 2.3 `operator-fixed-source` / `operator-eigenvalue` — the **posing** reading

A posing equation states which operator is on the left and what is on the right. The code
that IS a posing is the code that *assembles* the two sides. For the fixed-source row that
is the two SN within-group arms (each binds `F = 0_wg` and routes the source to `q_ext`);
for the eigenvalue row it is the `(1/k)Fψ` right-hand side, which is the single thing
distinguishing the two equations.

Alternatives, deliberately NOT taken, with the structural reason each was rejected:
- `SourceIteration` / `KrylovAcceleration` — method-agnostic drivers; they solve
  `(A − Σ gains)ψ = q` without knowing what the leaves are. Consumers.
- `SNSolver.solve_fixed_source` — a three-line `if/elif` dispatcher; assembles nothing.
- `power_iteration` for `operator-eigenvalue` — it realizes the *generic standard form*,
  whose own labels (`eigen-standard-form`, `eigen-resolvent`, `eigen-k-posing`) are
  unclaimed on this page and are its proper home. Declaring it here would poach them.
  (It is the natural second declaration if a whole-loop target is wanted.)

---

## 3. The NONE answers, grouped by KIND — the narrowing signal for #82

19 of 40. **This grouping is the transferable output**: it names the statement classes an
`implements` inference should never guess at, because no code can implement them.

**Identity — a mathematical fact stated as a landmark or a counter-example (4)**
`apply-solve-parallel-identity`, `apply-solve-neumann-series`,
`apply-solve-neumann-expansion`, `apply-solve-denominator-inequality`.
All four live in one section (`§apply is linear in the operator; solve is not`) whose
*purpose* is to exhibit what the code deliberately does **not** do. Two of them contain a
`≠`. Their tell: the section they sit in argues a negative.

**Law — a property the algebra obeys, enforced by ABSENCE (1)**
`solve-does-not-distribute`. Enforced by `OperatorSum` carrying no `solve` verb and
`StreamingOperator` declaring no `inverse()`. An absence is not an implementation, and
there is no symbol to point at — which is exactly why the guessing engine attaches 9
random ones.

**Canonical form — a motivating decomposition the tree deliberately does not realize (1)**
`streaming-as-tensor-product-sum`. `[M]` its named type `SumOfTensorProductsOperator` has
zero production consumers. Tell: the page cross-references a *rejection* rationale
(`operator_tensor_network`, "why streaming resists a clean factorisation").

**Definition — a name, a taxonomy, or a verb with no declaration site (2)**
`operator-solve`, `carrier-grid-cell`.

⭐ **The cross-cutting tell, and it is greppable:** every NONE row's `.. (vv-status
rationale)` block (where it has one) says some variant of *"Mathematical identity"*,
*"Definitional"*, *"Not a solver claim"*, or *"a lens, not a type"* — while every
declarable row's rationale names a **verb, a value, or a file**. The page already labels
its own two classes; nothing reads the label. That is a cheap narrowing feature for #82,
independent of any name-token heuristic.

⚠ Counter-example that stops this from being a rule: `carrier-grid-operator-typing`,
`harmonic-frame-is-galerkin` and `scattering-carrier-grid` all carry
"Structural / representational identity — not a solver claim" rationales and are
**declarable anyway**, because a *typing rule* can have a materialized carrier (a class
declaration, a Protocol parameter list, a set of typed methods) where an *identity* or a
*law* cannot. So the split is `{identity, law, canonical-form} → NONE` and
`{typing-rule, definition} → check for a declaration site`.

---

## 4. Four doc-drift findings the declarations exposed

**A. ⛔ `keff-as-integrated-rates` is present-tense-FALSE — do not declare until fixed.**
The equation reads `k = (∫⟨νΣf,φ⟩dV + (n,2n)) / ∫⟨Σa,φ⟩dV`. The shipped
`SNSolver.compute_keff` computes `production / (absorption + leakage − emission_n2n)`:
`(n,2n)` moved to the **denominator with a minus** and the vacuum-boundary **leakage term
was added** (#291 + the R7 convention, both dated 2026-07-03 in the method's own
docstring, whose `.. math::` block is correct). The stale numerator
`∫⟨νΣf,φ⟩dV + (n,2n)` is now exactly `SNSolver.compute_production_rate` — the
renormalisation scale anchor, a different quantity. Declaring `compute_keff` against the
current text would ship a false claim at `confidence=1.0`; fix the equation in the same
change.

**B. `reaction-rate-kinf-oracle` spells its `A` two ways.** The page: `A = diag(Σ_a) − Σ_s`.
The algebra-of-record it is verified against (`_infinite_medium_matrices`):
`A = diag(Σ_t) − (Σ_s + 2Σ_2)ᵀ`. Equivalent in the infinite medium, but the **transpose**
and the `2Σ_2` term are in the code and absent from the equation. The declaration is
sound; the equation should adopt the code's spelling.

**C. `scattering-as-tensor-product-sum` writes `S`, implements `Λ`.** The label's body says
`S = Σ_ℓ P_ℓ ⊗ Σ_{s,ℓ}`, but ~900 lines earlier the same page records that this shape was
"**considered and rejected**" for `S` (R and M mix the ordinate and harmonic axes, so they
are not valid tensor-product factors; tracked on #260). What ships is the **moment-space
Λ** — and `LegendreMomentScattering`'s class docstring cites *this exact label* for
itself. Recommend re-spelling the equation's left side `Λ` (or `S_moment`) so the
declaration is unambiguous.

**D. `apply-distributes` holds by theorem, not by code path, on the operator it names.**
`StreamingCollisionOperator.apply` **overrides** `OperatorSum.apply` and computes
`loss_representation.loss_action(σ, ψ)` directly — so the production `(L+C)` path never
executes `Lψ + Cψ`. Its own docstring says the inherited leaf-sum is "value-equal only by
the forward-direction affine-in-σ coincidence", i.e. by `streaming-action-pure-l`. The
declaration therefore correctly targets the **generic** `OperatorSum.apply` and must NOT
be extended to `StreamingCollisionOperator.apply`.

**E. (minor) One dead Python-domain cross-reference.** The page cites
`:func:`~orpheus.transport.spatial.diamond.cell_balance_for_streaming``; the function
lives in `orpheus.transport.spatial.cell_balance` and is merely *imported* into
`diamond`. Renders as plain text with no `-W` warning (the standard silent class).

---

## 5. Two mechanical cautions for whoever writes the directives

- **Target `_apply_impl`, never `apply`, on the three `singledispatchmethod` operators**
  (`MultiplicationOperator`, `ScatteringOperator`, `FissionOperator`). `apply` exists only
  under `if TYPE_CHECKING`; at runtime it is an alias, and `[M]` the graph has a node for
  `_apply_impl` and none for `apply`. The per-carrier `@…register` arms are anonymous
  `_` functions with no nodes of their own, so `_apply_impl` is the finest addressable
  grain.
- **`orpheus/derivations/` is declarable.** It is inside the package, not under `tests/`;
  the `reaction-rate-kinf-oracle` implementers live there and the test *imports* them.
