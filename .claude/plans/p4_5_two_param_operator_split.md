# Phase 4.5 — the two-param `LinearOperator[Domain, Codomain]` split

> **Part of the Frame projection/reconstruction machinery campaign** (`.claude/plans/frame_projection_machinery.md`).
> Inserted as **P4.5**, to run BEFORE P5 (energy condensation) / P6 (adjoint homogenisation): the user's
> directive is "polish the architecture so we don't blindly build P5/P6 on an operator algebra that must
> further change." Branch `refactor/operator-inverse-algebra`. Surgical, main-agent-direct, user-steered;
> NO `method-implementer`. Host env → `.venv/bin/python` (runtime 3.14). **Cold-pickup-ready** (written
> because the session compacts before implementation — read this end-to-end, then `operator.py`/`vector.py`).

---

## ⚑ SCOPE EXPANSION + STATUS (2026-06-24) — read FIRST

> **✦✦ SESSION 2026-06-25 — COLD-PICKUP MARKER (read this first; supersedes the P4.5a-uncommitted note below).**
> **W-A, W-B, W-C, W-D LANDED + W-E RESOLVED (docs)** on `refactor/operator-inverse-algebra` (NOT uncommitted — reconcile
> against git, not the P4.5a paragraph below which is now archaeology):
> - `70d5d78` **W-A** base collapse (one invariant `LinearOperator(Protocol[Domain,Codomain])`, Mixin retired) + `53e54c6` plan.
> - `f55bdad` **W-B** projection ABCs generic `[Domain,Codomain]` + `b88ceec` plan.
> - `34f8eaa` **W-C** TimedFullField→FullField solve/residual boundary + `bf3adf0` plan.
> - `61e8ddd` **W-C follow-on** scattering/fission runtime register keys → `FullField` (the `@overload` STUBS stay
>   `TimedFullField` for W-F — see the W-F HEAD START note).
> - `0610b39` **W-D** cross-method recognition: C/S/F gain real composite `domain`/`codomain`; the `(L+C)−S` guard
>   is LIVE; de-SN-ified `sn_full_field` → `full_field`. (See the W-D section below for the 1-site/C-activator learnings.)
> - `d30d4a6` **W-E** resolution — docs/theory ONLY (no scattering/solver code): the phase-space-promotion premise was
>   REFUTED (the angular Frame is scattering's EIGENBASIS, Funk–Hecke; `S=R∘Λ∘M`=spectral theorem); the ruling + the two
>   relocation tripwires (independent-L consumer; cross-method-measure → #261) are archived. See the W-E CLOSED note below.
> **Net: pyright `orpheus/` 419 → 412 (Δ−7 cumulative; W-D Δ0), ZERO regressions; all gates green** (SN 1414 + the
> 7-and-only-7 pre-existing reds, 0-ULP scattering canary, tests/numerics 684, solve+eigenvalue 110, verification
> L0/L1 88, Sphinx -W clean).
> **W-E CLOSED — refuted + resolved (2026-06-25).** The phase-space-promotion premise was FALSIFIED: the angular
> Frame is the scattering operator's EIGENBASIS (Funk–Hecke; `S=R∘Λ∘M` = the spectral theorem `UΣU*`), so it is
> ALREADY correctly scattering-owned. Move 1 (promote → `SNMesh`) is DEAD; the resolvent borrow is NOT a smell; the
> typed M/R faces DEFER to P6 (their `S.H = M.H∘Λ.H∘R.H` consumer). The unifying ruling (operator-symmetry → frame
> ownership; Galerkin-angular vs PG-energy/spatial) is archived in the theory pages + the `frame.conjugate` 2-cell
> docstring. See the rewritten W-E section below. **NEXT = W-F** (the `@overload` retirement — register sigs already
> done in W-C; restructure the overload stack) OR **W-H** (independent `_DISPLACEMENT_CLS` derive, quick) → W-G
> (docs — the W-E eigenbasis ruling already landed there).

**P4.5a (the numerics typevar foundation) is BUILT + VERIFIED (uncommitted, lands WITH 4.5b):** `operator.py`
two-param `LinearOperator[Domain, Codomain]` / Mixin `Generic[V, W]` / composers threaded / `_AdjointOperator`
`Generic[V, W]` / `ZeroOperator` `cast` bridge; `vector.py` prose corrected; `pyproject` floor → 3.13 +
pyright `pythonVersion=3.13`. **CLI pyright: operator.py+vector.py 0 errors; full `orpheus/` 419→419 (Δ+0)**;
`tests/numerics` 684 passed; the 7 `sn/operators` reds proven pre-existing (byte-identical on HEAD). The
spike (`$CLAUDE_JOB_DIR/tmp/spike_two_param.py`) banked corrections C1/C2 into §2.4.

**THE USER RE-SCOPED P4.5 (2026-06-24) from a narrow S/F retype into a FUNDAMENTAL cross-subsystem carve.**
This session became EXPLORATION + PLAN REFINEMENT; implementation follows AFTER a refinement + re-compaction.
Standing decisions that OVERRIDE the §3.1–§3.4 / §6 drafts below:

- **D1 — base carrier pair is `[FullField, FullField]`, NOT `[TimedFullField, FullField]`.** The loss leaves
  are ALREADY `FullField`-typed; the flux→source role change is at the **bulk leaf** (`AngularFlux →
  AngularSourceSink`), NOT the composite. The composite stays `FullField` (same type both sides). **`TimedFullField`
  is reserved for a FUTURE `TimeOperator` (not yet built)** — so the operator boundary should speak `FullField`,
  with `TimedFullField` confined to where time/history is genuinely load-bearing. (`TimedFullField` has 225
  usages vs `FullField`'s 179 — the open question the transport explorer is mapping: how much of that 225 is
  incidental-at-the-boundary vs real time-carrying. **This supersedes the `[TimedFullField, FullField]` in
  §3.1/§3.2.**)
- **D2 — leverage EVERYTHING, retire EVERY awkwardness.** The `[Domain, Codomain]` axis is a fundamental
  change; inventory every `@overload` (exactly TWO modules: `sn/scattering.py` + `sn/fission.py` — fission is
  the "one more, made reluctantly"), every `cast`, every `.values` unwrap, every operator that changes carrier
  type, across numerics + transport + sn. Explicit > implicit.
- **D3 — challenge fundamentals, including `LinearOperatorMixin`'s existence.** "Why do we have an Operator
  Mixin? Does the two-param change bring a better way?" The Protocol(variant `[Domain,Codomain]`)/Mixin(invariant
  `[V,W]`) split is itself under review.
- **D4 — exploration-first.** 4 parallel agents dispatched 2026-06-24 (explorer×3 numerics/transport/sn +
  cross-domain-attacker on the Mixin challenge). Synthesise → **refine THIS plan's §3–§6** → re-compact →
  implement. Do NOT implement the §3/§4 sub-phases as written until the refinement lands.
- **D5 — C/F/S/B are CROSS-METHOD operators (2026-06-24 user).** Scattering / Collision / Fission / Boundary
  are method-agnostic (every method scatters/collides/fissions) — they sit in `orpheus/sn/` ONLY because the
  relocation effort (#261) hasn't happened. **Only `StreamingOperator` + `LossRepresentation` are method-specific**
  (the WDD sweep discretisation) and genuinely belong in `sn/`. Consequence: their full-field space is NOT
  `sn_full_field` — it is a TRANSPORT-level space; the SN-flavoured name (`full_field_space.py:170`) is a
  misplacement symptom. The two-param `Operator[Domain,Codomain]` algebra IS the cross-method abstraction, so P4.5
  makes C/F/S/B honest + relocation-READY; the file-move itself is #261 (core/adapter split: XS-core → transport
  L2, quadrature angular adapter → sn L3). **Don't tag C/F/S with an SN-specific space — recognise them as
  transport-level.**

---

## ✦✦ THE CATEGORICAL FOUNDATION (2026-06-25) — the WHY behind P4.5

A deep architecture exploration (4 explorers + 2 cross-domain-attacker digs + 2 pyright spikes, 2026-06-24/25)
established what P4.5 *is*, and it is load-bearing for every workstream. **Bank this; it is the rationale.**

**The carriers form a (Representation × Role) grid that is a DOUBLE CATEGORY** — equivalently, a category
**fibered over Representation, with Role the fiber coordinate, carrying a TORSOR on the Flux fiber**:
- **Objects** = the grid cells `(Representation, Role)`. Representation ∈ {Angular, Moment, Scalar, Trace}
  (sets SHAPE + a change-of-basis); Role ∈ {Flux, Source, Residual, Displacement} (sets the ARITHMETIC
  interface — the #208 affine torsor: Flux is an affine point, Source/Residual/Displacement are vectors).
- **Horizontal morphisms** = Representation-changes (M=analysis, R=reconstruction): a base change that FIXES
  the fiber ⟹ **role-generic**. Built from `(basis, measure)` = the Frame.
- **Vertical morphisms** = Role-changes (the cross-sections C=σ_t, Λ=Σ_{s,ℓ}, F=χ⊗νΣf): a fiber morphism
  identical over every base ⟹ **rep-generic**. The role-change IS the physics (the cross-section).
- **Scattering `S = (1/W)·R∘Λ∘M` is the 2-CELL** `frame.conjugate(Λ)` (a vertical morphism Λ conjugated by
  the horizontal adjoint pair M/R). The role-change is LOCALIZED at Λ; M/R preserve role. The bit-identical
  windowed-vs-full-angular crosscheck ALREADY in the code is its **interchange-law coherence witness**.

**THE KEY IDENTITY — a grid cell IS an operator's `(Domain, Codomain)`.** The parametrization belongs on the
OPERATOR, NOT the carrier: `LinearOperator[Domain, Codomain]` IS the typed grid traversal. `M :
LinearOperator[AngularFlux, HarmonicMomentFlux]` (a horizontal edge); `Λ : LinearOperator[HarmonicMomentFlux,
HarmonicMomentSourceSink]` (a vertical edge); `S : LinearOperator[AngularFlux, AngularSourceSink]` (the 2-cell).
**P4.5 (`LinearOperator[Domain, Codomain]`) is therefore the right and COMPLETE machinery for traversing the
grid** — every branch of the exploration (4-operators-vs-methods, frame ownership, Rep×Role, "type the grid")
converges back onto P4.5, now with a categorical reason it is correct, not merely expedient.

**THE CARRIERS ARE ALREADY THE ELEGANT NORMAL FORM — DO NOT REDESIGN THEM.** Both digs converge:
- Role changes the arithmetic interface (the torsor) ⟹ **Role MUST be a class** (a phantom `Generic[Role]`
  param is erased and cannot specialize `__add__` — refutes `Field[Rep,Role]` and rep-outer `Angular[Role]`).
- Representation changes shape ⟹ **Rep MUST be a class** (refutes role-outer `Flux[Rep]`).
- The only both-classes form with role-arithmetic-once + rep-shape-once is **the flat MI leaf
  `AngularFlux(FluxRole, AngularField)` that ALREADY EXISTS**. A parameterized `Carrier[Rep,Role]` is
  structurally impossible AND would break the **runtime units gate** (`type(self) is type(other)` — the units
  enforcement) via generic erasure. No novel encoding exists.
- The role-axis "asymmetry" (Flux/Displacement are mixins; Source/Residual are bare leaves) is **PRINCIPLED,
  not a defect** — the type-vs-property rule applied correctly: mint a role-object ONLY where a non-identity
  morphism lives (Flux has the torsor, Displacement has diagnostics; Source/Residual are plain vector roles
  with no special arithmetic ⟹ correctly carry nothing). **Do NOT uniformize the role axis — that is ceremony.**

**PRINCIPLED EXCEPTIONS — leave both (user-confirmed 2026-06-25):**
- **`HarmonicMomentResidual` is absent** — moment-space is never the subject of a balance equation (balances
  are bulk-angular `(L+C−S−F)ψ−q` and boundary), so there is NO `from_balance` consumer. Principled hole, not
  an oversight. Do not mint it.
- **The `iso+aniso → AngularSourceSink` source-injection** (`angular_source_sink.py`/`scalar_source_sink.py`):
  a hand-rolled subspace-containment Representation-traversal inside a dunder. **User personally endorsed it at
  creation.** Leave it.

**Memos (durable):** `.claude/agent-memory/cross-domain-attacker/rep_role_grid_double_category_frames.md`
(the double-category/fibration/torsor verdict; `lessons.md` L-004 sharpened with the class-vs-phantom decider);
`.claude/agent-memory/cross-domain-attacker/flux_to_sourcesink_operator_contract_frames.md` (the affine-map
operator structure); `operator_protocol_mixin_collapse_frames.md` (the Protocol/Mixin collapse, W-A).

**Naming (user, 2026-06-25): the operator's two parameters are `Domain` / `Codomain`, spelled in FULL** — not
`Domain`/`Codomain` (the abbreviation said nothing; "Domain" already reads as the input). Applied to P4.5a's
`operator.py` (the variant Protocol pair; the invariant Mixin still spells `V`/`W` + auxiliaries `Cmid`/`D2`
until W-A collapses the split and unifies them to the single invariant `[Domain, Codomain]`).

---

## ✦ REFINED WORKSTREAMS (2026-06-24) — AUTHORITATIVE; supersedes §3–§4 below

Synthesis of the 4 exploration maps + 2 spikes (`$CLAUDE_JOB_DIR/tmp/spike_two_param.py`,
`spike_collapse.py`). The old §3 (consumers) / §4 (sub-phases a–e) stay as superseded detail; THIS section is
the plan of record. Eight workstreams, dependency-ordered. Each: surgical, main-agent-direct; green +
CLI-pyright ≤ baseline (419 @3.13) + runtime-zero where typing-only; ff-merge.

**Convergent truths all four agents found:**
- The bulk operators (L/C/B/L+C) are `[FullField, FullField]` at the **composite** — the flux→source role
  change is at the **bulk leaf** (`AngularFlux→AngularSourceSink`), never the composite. The genuine `Codomain≠Domain`
  demonstrators are `projection.py` (mono-type ndarray, real space-change) and the **HarmonicFrame faces**
  (`[AngularFlux, HarmonicMomentFlux]`). Sequence so a low-risk demonstrator proves the pattern before the
  variance-heavy / 0-ULP-adjacent carves.
- **Two complementary levels stay**: the static carrier type (`Domain`/`Codomain`, `Vector`-collapsed at the numerics
  layer — types `apply`) and the runtime `(name,shape)` FunctionSpace guard (`domain`/`codomain` — value-based,
  no-op when `None`). Neither subsumes the other. The split's win: more operators populate `domain`/`codomain`
  honestly → strengthens the runtime guard's coverage (fewer `None`-skips).
- The `@overload` stacks are **convenience multi-dispatch, NOT the endomorphism lie**. `[Domain,Codomain]` retires the
  `apply(x)->Any` fallback + makes the primary arm honest; the secondary-carrier arms (ScalarFlux, bare
  ndarray, the TimedFullField composite) collapse only by routing `apply` through the shared typed `kernel`.
- `TimedFullField` history is **dead-on-arrival in production** (`advance`/`at_lag` called by zero solve code;
  every production build passes `_history=()`); it is the fuller-view ANCHOR for a future `TimeOperator` —
  **confine, don't delete** (coding-standards fuller-view-oracle exception).

### W-A — THE BASE COLLAPSE (keystone; REVISES P4.5a). ✅ **LANDED 2026-06-25 (uncommitted).**
**Result:** `operator.py` collapsed to ONE `@runtime_checkable LinearOperator(Protocol[Domain, Codomain])`
(invariant; `Codomain=Domain` default) carrying the algebra dunders as default-method bodies; the
`LinearOperatorMixin` class (270 lines) + the `TYPE_CHECKING` phantom `apply` stub + the `V`/`W` typevars
RETIRED. `_AdjointOperator` keeps the explicit `Generic[Domain, Codomain]` (C2 PEP-696 ordering — the
`[Codomain, Domain]` base swap still needs it). All 13 production operator classes + 7 test conformers +
10 doc `:class:` cross-refs migrated `LinearOperatorMixin→LinearOperator` (3 import dedups: numerics/__init__,
sn/scattering, geometry _bound_compat — the last dropped a now-redundant TYPE_CHECKING import). Variance
DROPPED (confirmed nothing relies on it; static type collapses to `Vector` at numerics). **GATE MET:**
CLI-pyright `orpheus/` **419→414 (Δ−5, ZERO regressions** — the 3 tightening-signature errors all pre-exist:
`solver.py:317` bare-`object`, `:748` resolvents that lack `domain`/`codomain` so they failed the OLD Protocol
too; capabilities-override + `.solve`/`.apply_transpose`-access errors are structurally identical to P4.5a's
base). `tests/numerics` 684✓, SN operators/solve/sweep/eigenvalue/verification **1414✓ + the 7-and-only-7
pre-existing reds (byte-identical ULP)**, **0-ULP scattering-kernel crosscheck canary GREEN**, `tests/transport`
286✓ (all 13 `isinstance` sites pass; the negative `not isinstance(func, LinearOperator)` is now MORE robust —
the tightening). Spikes: `spike_collapse.py` (the collapse) + `spike_mi.py` (the `(LinearOperator, ABC)` /
`(LinearOperator, RegistryMixin, ABC)` / bare / string-subscript MI combos — all 0 pyright + runtime clean).
NEXT = W-B. (Original spec retained below for the record.)

Collapse `LinearOperator(Protocol[V])` + `LinearOperatorMixin(Generic[V])` into ONE
`@runtime_checkable Operator(Protocol[Domain, Codomain])` — **invariant** Domain/Codomain (`Codomain=Domain` default for `[V]≡[V,V]`)
— carrying the algebra dunders as **default-method bodies** that explicit subclasses inherit. `spike_collapse.py`:
pyright 0 errors + runtime isinstance holds + composition runs. **Retires:** `LinearOperatorMixin` (whole class),
the `TYPE_CHECKING` phantom `apply` stub, the variant `Domain`/`Codomain` ceremony, the `_AdjointOperator` PEP-696
ordering workaround (re-evaluate — the `[Codomain,Domain]` swap may still need explicit ordering). The numerics
explorer's "Protocols can't carry impls" was FALSE (spike-disproved); the cross-domain-attacker's Frame-1 wins.
- **Name: KEEP `LinearOperator`** (fork 2 resolved — it encodes linearity via the adjoint/`solve`/Krylov
  contract; collapse the Mixin INTO `LinearOperator(Protocol[Domain,Codomain])`, do not rename).
- **Variance check:** confirm NO consumer relies on operator covariance/contravariance (15/18 leaves inherit
  bare; the static type collapses to `Vector` at the numerics layer ⟹ variance buys nothing).
- **Migration detail:** dunders-on-Protocol TIGHTEN `isinstance` to also require the dunders present — verify
  the 14 `isinstance(x, LinearOperator)` test sites (all current ones inherit the Mixin ⟹ pass; the negative
  functional test gets MORE robust). Gate: pyright ≤419; `tests/numerics`+operator/iteration green; the 14
  isinstance sites pass.

### W-B — `projection.py` honest faces. ✅ **LANDED 2026-06-25 (uncommitted).**
`AnalysisOperator`/`ReconstructionOperator` are now generic `LinearOperator[Domain, Codomain]` with
`apply(self, x: Domain) -> Codomain` (the role IS `M : V → W` — two-typed by definition, not speculation);
dropped the now-unused `import numpy as np`. **LEARNING (load-bearing for W-E):** the numerics frame faces
(`_FrameAnalysis`/`_FrameReconstruction`) STAY BARE (un-subscripted) — `np.ndarray`/`NDArray` does NOT satisfy
the `Vector` bound (numpy's `__add__` stubs don't return `Self`), so `AnalysisOperator[NDArray, NDArray]` is
REJECTED by pyright (`reportInvalidTypeArguments`) and cascades to the `frame.conjugate` `OperatorProduct`
composition. This confirms §2.5: ndarray-level numerics primitives stay bare; the V/W distinction lives in the
runtime `domain`/`codomain` FunctionSpaces. **⟹ W-E's transport wrappers MUST subscript with real `Field`
carriers (`AnalysisOperator[AngularFlux, HarmonicMomentFlux]`), NEVER ndarray.** Gate MET: pyright
projection+frame 0 errors, full `orpheus/` **414 (Δ0, pyright-neutral)**, tests/numerics 684✓. Runtime-zero
(subscript erased; faces unchanged). NEXT = W-C / W-D.

### W-C — `TimedFullField → FullField` boundary confinement. ✅ **LANDED 2026-06-25 (`34f8eaa`).**
The SN operator **solve/residual INPUT boundaries** speak `FullField`; `TimedFullField` confined to the driver
iterate. `InvertibleOperator.solve`/`solve_moments`/`_solve_timed_full_field` `rhs`/`initial_guess` → `FullField | None`
(dead `AngularFlux` arm dropped); iterate OUTPUT stays `TimedFullField`. The 2 `isinstance` guards (operator.py)
relax `TimedFullField → FullField` (timed iterate passes via inheritance; bare ndarray/`AngularFlux` still
rejected). `evaluate_residual` now returns bare `FullField` (residual = a depth-0 timed field = the FullField
degenerate); `boundary_vs_interior_split` takes `FullField`. The 2 dead doc-refs fixed + the dead-AngularFlux
narrative pruned (archivist). Gate: pyright `orpheus/` **414→412 (Δ−2, zero regressions)**; SN 1414 + the 7
pre-existing reds; 0-ULP canary green; tests/numerics 684; Sphinx -W clean.
**LEARNINGS (load-bearing):**
- **`history_depth=0` IS the `FullField` degenerate** (user insight): a `TimedFullField` with depth 0 has
  permanently-empty `_history` (`advance` trims to `[:0]`), so `at_lag(k≥1)` raises — only the current frame,
  identical content to `FullField`. The history_depth echo on a bare-`FullField` input is principled `else 0`.
  `evaluate_residual` already constructed `history_depth=0` ⟹ it was always "effectively FullField".
- **Blast radius was NOT funnelled through `_within_group_triple`** (that is the OPERATOR triple, no
  `TimedFullField` in it) — the field-type retype spread across **operator.py solve-family + solver.py
  residual/source-builders + loss_representation.py's 7 `sweep` strategies** (9 `initial_guess` annotations).
  ~13 signature edits, not a one-line funnel. (The explorer's "single source of truth" claim was about
  operators, not the field carrier.)
- **scattering/fission `@singledispatch` RUNTIME register-key → DONE (W-C follow-on, `61e8ddd`); `@overload`
  stubs → DEFERRED to W-F** (user directive 2026-06-25: change the runtime signatures now so W-F focuses on the
  overloads with the signatures already done). The register `def _(self, psi: FullField)` is runtime-safe (the
  `FullField` arm catches the `TimedFullField` iterate via MRO — verified: 86 scattering/fission tests + the 0-ULP
  crosscheck canary green). The `@overload` STUBS stay `TimedFullField` because a `FullField`-first `@overload`
  SHADOWS the bulk-field arms (`ScalarFlux`/`AngularFlux`/`ndarray`) — pyright's overlap detection treats the
  composite BASE `FullField` as ⊇ those (the narrow subclass `TimedFullField` did not), so they become "never used"
  and `S.apply(scalar)` → `Unknown` (`reportOverlappingOverload`). ⟹ W-F RESTRUCTURES the overload stack, not a key
  swap. (register/overload now intentionally disagree — the transient W-F closes.)
- **Pre-existing doc debt the archivist FLAGGED (future sweep, NOT W-C):** 4 dead `:class:`~orpheus.sn.angular_flux.AngularFlux``
  refs in `solver.py:1188/1366/1711/2366` (canonical = `orpheus.transport.fields.angular_flux`); `operator.py:304-305`
  "apply accepts ONLY TimedFullField" is stale vs the live `FullField` matvec guard (#257 S8a).
NEXT = W-D (cross-method recognition — needs proactive `test-architect`).

### W-D — cross-method recognition + §3.1 composition-guard close (D5). ✅ **LANDED 2026-06-25 (`0610b39`).**
`ScatteringOperator`/`FissionOperator`/`CollisionOperator` gained real `domain`/`codomain` returning the composite
full-field space; the `(L+C)−S` guard is LIVE and `IncompatibleOperatorComposition` is reachable on a mis-spaced
operand. De-SN-ified `sn_full_field` → `full_field` (user-approved 2026-06-25 — the operators that advertise it
are cross-method). **S/F thread the `FullFieldSpace` via `from_solver_data`** (a numerics `FunctionSpace`, NOT an
SN mesh — D5-honest, relocation-ready for #261); **C** uses its existing `self.sn_mesh.full_field_space`.
**LEARNINGS (load-bearing; test-architect re-validated its 2026-06-24 memo against HEAD):**
- **ONE production `OperatorSum` site, not two.** `evaluate_residual` (`solver.py`) has ZERO production callers
  (test-only); the SOLE production OperatorSum is `InvertibleOperator(L,C)` (`_within_group_triple`, built on every
  within-group solve). SI/Krylov/eigenvalue realise the matvec as `L.apply − Σgᵢ.apply` (no OperatorSum) ⟹ the
  guard never fires on the converging path — only at the `L+C` build.
- **Role split:** **C is the production gate-activator** (mis-naming C reds the `L+C` build ⟹ every solve); **S is
  the residual-test vehicle** (mis-naming S reds the test-built `(L+C)−S−B`); **F NEVER composes** into a production
  OperatorSum (always `F.apply` separately) ⟹ its space is the honest cross-method tag, gated at construction-level
  only. Minimal gate-activating set = {C, S}; F rides along (bias-to-completion).
- **Adjoint-unreachability of `(L+C−S−F−B).H` is preserved by the CAPABILITY lattice** (S/F advertise no
  `apply_transpose`), NOT by `None` spaces — W-D's spaces don't touch it. The stale `StreamingOperator.domain`
  "C/S/F report None" docstring + the `full_field_space`/`operator_algebra.rst` guard-claim prose were updated.
- **`LegendreMomentScattering` precedent** already landed (P2/#47): it carries real `frame.basis_space` (MOMENT
  space — NOT the composite; don't conflate) and retired its `cast`.
- **The de-SN-ify was scoped to the COMPOSITE name only.** The leaf names `sn_bulk` (`geometry.py`) / `sn_trace`
  (`trace_space.py`) stay (they are `compare=False` block metadata, out of the guard's `(name,shape)` identity, and
  belong to #261's relocation). The 3 remaining `sn_full_field` mentions in source are historical prose ("de-SN-ified
  from …"), not functional literals.
Gate MET: pyright `orpheus/` **412 (Δ0)**; 7-and-only-7 pre-existing reds; 0-ULP scattering canary; tests/numerics
684; solve+eigenvalue 110; verification L0/L1 88; Sphinx -W clean. Tests: T4.5 positive / negative-S (+Mode-11
`co_qualname` anchor to `OperatorSum.__init__`) / production-teeth-C in `test_typed_residual_evaluation.py`.
NEXT = W-E.

### W-E — angular Frame ownership. ✅ **CLOSED — REFUTED + RESOLVED (2026-06-25). The frame is already correctly placed.**
**The phase-space-promotion premise (the whole original W-E) is FALSIFIED.** A literature + structural investigation
(literature-researcher + cross-domain-attacker, 2026-06-25) tested the user's claim *"every method that needs the
HarmonicFrame needs it BECAUSE of scattering"* and CONFIRMED it on two independent legs:
- **Structural (Funk–Hecke — a theorem, not an analogy):** the scattering kernel `Σ_s(Ω·Ω')` is a zonal kernel on
  S² ⟹ the spherical harmonics are its EIGENFUNCTIONS, eigenvalues = the Legendre moments `Σ_{s,ℓ}` = the diagonal
  of `LegendreMomentScattering` Λ. So **`S = R∘Λ∘M` IS the spectral theorem `A = UΣU*`** (M=U analysis into the
  eigenbasis, Λ=Σ spectrum, R=U* synthesis). Streaming `Ω·∇` is the ℓ=1 tensor operator — does NOT commute with
  SO(3), couples ℓ↔ℓ±1 (the PN recurrence), is NOT diagonalized ⟹ the basis exists to diagonalize SCATTERING;
  streaming merely tolerates it.
- **Literature (no falsifier across SN/MoC/CP/PN/FCS/random-ray):** every documented flux→SH-moment projection is
  anisotropic-scattering-rooted (Hébert §3.3 Eq. 3.55 M used ONLY in Eq. 3.54 scattering source; the integral form
  wants isotropic sources Eq. 3.42; Brockmann Eq. 47; Ahrens LDO REMOVES M by reformulating the scattering source).
  The only structurally-independent falsifier — an output detector functional of order `L_d > L_scatter` — is
  ABSENT from ORPHEUS (the sole output functional is the ℓ=0 scalar flux, computed off-frame).
**⟹ THE FRAME BELONGS ON SCATTERING (its eigenbasis), NOT the phase-space. Move 1 (promote → `SNMesh`) is DEAD.**
The placement is ALREADY correct: generic machinery in `transport/frames`+`numerics` (shared with the homogenization
PG frame), `ScatteringOperator` holds the constructor at `L=scattering_order`, Λ on scattering, M/R on the frame.
The resolvent "borrow" (`_maybe_window` reads `scattering_op.frame`, `solver.py:617`) is **NOT a smell** — the
windowed resolvent's in-sweep M IS the scattering projection (reduces the iterate to moments BECAUSE scattering
consumes moments), same frame, same L. Nothing to de-borrow; no code change to scattering/solver.
**THE UNIFYING PRINCIPLE (the campaign prize):** *an operator owns its frame IFF the frame is its eigenbasis.*
Angular scattering = SO(3) symmetry → SH eigenbasis → **Galerkin, scattering-owned**; energy condensation / spatial
homogenization = NO symmetry → no eigenbasis → **Petrov-Galerkin, solution-weighted, owned by no operator** (the
projection verb / the test side). ONE structural cause for the campaign's Galerkin-vs-PG split (re-confirms fission
owns no angular frame — energy has no eigenbasis).
**RELOCATION TRIPWIRE (record, don't act):** the constructor relocates from `ScatteringOperator.frame` to the
neutral `Quadrature.angular_frame(L)` (already exists, with an anticipating naming tripwire at
`directional.py:439-447`) ONLY when a 2nd consumer with an L INDEPENDENT of scattering_order arrives (a detector
functional, or PN/SPN with `L_flux > L_scatter`). None today. **SECOND trigger (user 2026-06-25 → #261):**
cross-method use of `ScatteringOperator` — a `HarmonicFrame` needs an angular MEASURE (`Quadrature`) that CP/MoC
LACK, so the frame can't live as a field on the shared operator; resolve by (a) relocate to `Quadrature.angular_frame`
OR (b) specialize `ScatteringOperator` per method (SN subclass holds the frame; measure-free cross-method base carries
only Λ). Recorded on #261 + the theory anchor `frame-eigenbasis-relocation-tripwire`.
**LANDED close-out:** the ruling is archived in `docs/theory/foundations/frame.rst` (capstone, anchor
`frame-eigenbasis-ownership`) + `spherical_harmonics.rst` + `operator_algebra.rst` (archivist, Sphinx -W exit 0);
`frame.conjugate(Λ)` docstring NAMES the double-category 2-cell + the spectral-theorem reading (`frame.py:206`);
memos `harmonic_frame_ownership_funk_hecke.md` (+ lesson L-009) + `sh_flux_moment_projection_root_cause.md` +
`feedback_defer_only_when_architecture_vague.md` + `feedback_mesh_owns_machinery_not_storage_init.md`.
**The typed composable M/R faces + adjoints (the original move 2) DEFER to P6** — they are P6's `S.H = M.H∘Λ.H∘R.H`
tooling, built WHERE the adjoint-homogenization consumer exercises them (user decision 2026-06-25).

### W-F — `@overload` retirement + the shared `Operator[Flux, SourceSink]` emission abstraction (the DEEP one).
**HEAD START (W-C follow-on, `61e8ddd`):** the RUNTIME `@singledispatch` register keys for the composite
scattering/fission `apply` arm are ALREADY `FullField` (`def _(self, psi: FullField)`) — W-F starts with the
runtime signatures done and FOCUSES on the `@overload` static stubs. **CAVEAT (load-bearing, why this is the deep
one):** a naive `@overload def apply(self, psi: FullField)` as overload-1 SHADOWS the bulk-field arms
(`ScalarFlux`/`AngularFlux`/`ndarray`) — pyright treats the composite BASE `FullField` as ⊇ those (the narrow
subclass `TimedFullField` did NOT shadow them), so the bulk overloads become "never used" and `S.apply(scalar)` →
`Unknown` (`reportOverlappingOverload`). ⟹ W-F must RESTRUCTURE the overload stack (reorder broad-last, or
collapse to the typed-kernel `apply`), NOT swap the overload key. The register/overload now disagree
(register=`FullField`, overload stubs=`TimedFullField`) — an intentional transient W-F closes.

S and F ALREADY share the abstraction: both are a typed `kernel` (`OperatorProduct` `R∘Λ∘M` / `TensorProduct`
rank-1) satisfying `IntegralKernelOperator`; **F is the rank-1 degenerate** `M_χ ∘ ProductionRate ∘ M_νΣf` of
S's frame. Make C/F/S declare the honest leaf base `Operator[<Flux>, <SourceSink>]`, route `apply` through the
typed kernel for the primary carrier, retire the `apply(x)->Any` fallback overloads. **The secondary-carrier
arms** (ScalarFlux→ScalarSourceSink for diffusion/CP/kinetics; bare-ndarray for `KEigenvalue`/depletion) are
genuine distinct `(Domain,Codomain)` pairs — they collapse only by migrating those consumers to typed carriers
(decide in-session: retire-and-migrate now vs keep a NARROWED honest overload set). The `TimedFullField`
composite arm becomes the `BlockRole.BULK` lift of the leaf (structural, not a dispatch arm) once W-C lands.
Gate: 0-ULP canary; full scattering/fission/solve/eigenvalue green; pyright DOWN (overload confessions retired).

### W-G — docs (archivist) + the categorical naming.
`operator_algebra.rst`: the collapsed `LinearOperator[Domain,Codomain]` contract (one base, dunders-on-Protocol,
the affine-map / linear-part-of-affine structure), **the (Representation × Role) grid as a DOUBLE CATEGORY**
(horizontal = M/R Representation-morphisms, vertical = C/Λ/F Role-morphisms, scattering = the 2-cell, the
crosscheck = interchange coherence witness — see the Categorical Foundation above), the carrier-vs-FunctionSpace
two levels, FullField-as-composite vs leaf-role-change, fission-as-rank-1-degenerate-on-the-energy-axis, the
cross-method placement (D5), the principled exceptions (`HarmonicMomentResidual` absence, the source-injection).
Rewrite the endomorphism prose. Fix the W-C stale refs. Sphinx `-W` clean.

### W-H — minor: derive `_DISPLACEMENT_CLS` (the one carrier transcription smell).
The only carrier-grid smell the categorical dig found (shape-2): every flux leaf HAND-SETS `_DISPLACEMENT_CLS`
(the `flux⊖flux → Displacement[Rep]` codomain edge transcribed per leaf). Derive it from a `Rep → Displacement[Rep]`
lookup instead. Tiny, independent, bit-identical. Gate: the affine-torsor tests (`flux−flux`, `flux+displacement`)
green. (Do NOT touch anything else in the carrier grid — it is already the elegant normal form; see Categorical
Foundation. Do NOT mint `HarmonicMomentResidual`; do NOT uniformize the role axis; do NOT parameterize carriers.)

### Dependency order
W-A (foundation, revises P4.5a) → W-B (proves pattern) → { W-C, W-D } (boundary + cross-method, ~independent) →
W-E (Frame promotion + faces + kernel + name the 2-cell) → W-F (overload retirement, deepest) → W-G (docs) ;
**W-H** (the `_DISPLACEMENT_CLS` derive) is independent — land it any time. **Proactive `test-architect` before
W-D/W-E/W-F** (0-ULP-adjacent). The test-architect's prior plan (the §3.1 guard pair, the assert_type pins, the
Mode-11 sentinel) feeds W-D/W-E.

### Forks — RESOLVED (2026-06-24 user) + remaining
1. **W-A direction** — full collapse (recommended) vs keep the P4.5a split. *[the keystone — awaiting final nod]*
2. **Rename `LinearOperator`→`Operator`? — RESOLVED: NO, keep `LinearOperator`.** The name EARNS its adjective:
   the adjoint (`.H`/`apply_transpose` — the Hilbert adjoint `⟨Ax,y⟩=⟨x,A*y⟩` exists ONLY for linear `A`),
   `solve` (the linear inverse), and the Krylov consumer contract (GMRES/BiCGSTAB require linearity) all encode
   linearity. Every current leaf is linear in ψ. A future NONLINEAR operator would lack these ⟹ it is a DIFFERENT
   abstraction; mint an `Operator` supertype THEN (defer-until-≥2-instances; zero nonlinear ops today). W-A keeps
   the name; it collapses the Mixin INTO `LinearOperator(Protocol[Domain,Codomain])`.
3. **W-E — RESOLVED + REFRAMED (2026-06-25): the "four edges vs singledispatch" was the WRONG question.**
   The Categorical Foundation answers it: a face is ONE `LinearOperator[Domain,Codomain]` naming ONE grid edge
   (no singledispatch, no `@overload`, no `Any` — build-primitives), and the genericity lives on the operator,
   not the carrier. But the FIRST move is OWNERSHIP, not the face count: promote the Frame from scattering-private
   to phase-space (the resolvent-borrow smell), make it the M/R factory (role-agnostic core, adjoint-paired), and
   build only the grid edges consumers use. `S = (1/W)·R@Λ@M` is the 2-cell `frame.conjugate(Λ)`. See the rewritten
   W-E. (`Domain/Codomain` retired → `Domain/Codomain`.)
4. **W-D / #261 — RESOLVED: P4.5 does cross-method recognition + de-SN-ify the composite space; defers the
   file-move to #261** (core/adapter split is its own effort). C/F/S/B become honest + relocation-ready.
5. **W-F — RESOLVED: lean AGGRESSIVE retirement.** No cross-method stranding exists — CP/MoC/Diffusion
   reimplement C/F/S inline on bare `(N,ng)` arrays and reuse NONE of the SN operators (#261 map), so the
   secondary-carrier arms (`ScalarFlux`, bare `ndarray`) have NO cross-method consumers. The decision is a pure
   SN-INTERNAL audit: migrate the few live SN consumers (the bare-`ndarray` fission arm the `KEigenvalue`/
   depletion drivers touch) to typed carriers; keep at most the ONE `ndarray` arm that is genuinely the scipy
   `to_flat`/`from_flat` serialization boundary; retire the rest as `Any`-laundering convenience. Default: retire.

---

## §0. Why P4.5 exists (the finding that justifies it)

A post-P4 architecture review (explorer + cross-domain-attacker maps, 2026-06-24) established that **every
bulk transport operator is a map `flux STATE → rate-density SOURCE/SINK`** — the balance
`(L + C − S − F/k)ψ = q` is a *typed operator algebra* where each term `Oψ` is a SourceSink and they sum in
the **closed** source/sink vector algebra. This contract is **already ~80% built and uniform**:

- `ScatteringOperator.apply(AngularFlux)→AngularSourceSink`, `FissionOperator.apply(ScalarFlux)→ScalarSourceSink`,
  `(L+C).apply(...)→FullField(AngularSourceSink, BoundarySourceSink)` — all realised.
- A **named shared abstraction** exists: `IntegralKernelOperator` Protocol (S and F both satisfy it — both are
  nonlocal integral kernels; fission is the **rank-1 degenerate** of scattering's `R∘Λ∘M` frame:
  `F = M_χ ∘ ProductionRate ∘ M_νΣf`), the `BlockRole` markers, and the closed source/sink `__add__` algebra.

**But the operator TYPE lies.** `LinearOperator(Protocol[V])` declares `apply(x: V) -> V` — *endomorphic*. The
non-endomorphic operators (S/F mapping Flux→SourceSink) patch over this with `@overload` "confession" stacks
(`scattering.py:1388`, `fission.py:466`) that literally say *"NOT an endomorphism V→V; maps each input carrier
to a DISTINCT output carrier."* These are Python's weak spot — `@overload` + a singledispatch runtime + an
`Any/Any` fallback — and they exist ONLY because the base type is wrong. The `vector.py` `V` comment
(`:144-151`) even asserts the single-`V` endomorphism is "honest" — a *documented falsehood* (the module
docstring's own words about the prior `ndarray→ndarray` lie; this is the next correction in the same series).

**P4.5 makes the operator type honest: `LinearOperator[Domain, Codomain]`, `apply(x: Domain) -> Codomain`.** This:
1. **Retires the `@overload` endomorphism confessions** — the base is honestly non-endomorphic.
2. **Carrier-types the kernel end-to-end** — `kernel = R∘Λ∘M` becomes `LinearOperator[AngularFlux, AngularSourceSink]`,
   so `kernel.apply(AngularFlux)→AngularSourceSink` *directly*. This **dissolves the P4 option-2 asymmetry**
   (the windowed arm is carrier-typed, the full-angular arm flows ndarray through the fused kernel and types only
   at the boundary): both paths become uniformly typed, the only difference being "windowed skips M."
3. **Unifies S/F/L under one `Operator[Flux, SourceSink]` contract** — the genuine uniform abstraction the
   `@overload` stacks approximate.

The cross-domain-attacker's structural read: the operator is the **linear part of an affine map**
`A_flux → W_source` (the flux is an affine point — no origin, `flux+flux` forbidden; the source is a vector —
`source+source` closed), with the cross-section `1/cm` unit-gain as the **fiber transition**. The two-param
`[Domain, Codomain]` IS that affine-map type. (Memo: `.claude/agent-memory/cross-domain-attacker/flux_to_sourcesink_operator_contract_frames.md`.)

---

## §1. THE GATING DECISION — PEP-696 availability (resolve FIRST)

`LinearOperator[V] ≡ LinearOperator[V, V]` (keeping all ~7 existing single-param subscript sites working) needs
a **PEP-696 TypeVar default** `Codomain = TypeVar("Codomain", default=Domain)`. Native `typing.TypeVar(default=...)` exists
only in **Python ≥3.13**. Current state (explorer Map 0):
- `pyproject.toml:9` `requires-python = ">=3.11"`; runtime `.venv` = 3.14.3; pyright infers target = **3.11**
  (no `pythonVersion` pin); `typing_extensions` is **absent**.

**On 3.11, `default=` is a syntax-unsupported kwarg on `typing.TypeVar`.** So P4.5 MUST choose ONE:
- **(A) Raise `requires-python` to `>=3.13`** + add `[tool.pyright] pythonVersion = "3.13"`. Cleanest (native
  `typing`), but a public-floor bump — **needs explicit user sign-off** (does ORPHEUS commit to 3.13+?).
- **(B) Add a `typing_extensions` dependency**, `from typing_extensions import TypeVar` (the back-port supports
  `default=` down to 3.11). Keeps the 3.11 floor; adds one dep. Lower-commitment, reversible.

**DECISION — (A) RESOLVED 2026-06-24 (user): raise the floor to Python ≥3.13.** Native `typing.TypeVar(default=…)`,
no `typing_extensions` dependency, no import shim. First implementation actions: set `pyproject.toml`
`requires-python = ">=3.13"` AND pin `[tool.pyright] pythonVersion = "3.13"` (so the CLI oracle matches the new
floor), in the SAME commit as the typevar additions (P4.5a). The runtime is already 3.14; this commits ORPHEUS to
dropping 3.11/3.12 support — note it in the changelog / a `CHANGES` line if one exists. (Option B / `typing_extensions`
is retired — do not use it.)

(Aside: PEP-696 defaults also enable `W = TypeVar("W", default=V)` for the Mixin's codomain. Same gate.)

---

## §2. THE DESIGN — the two-param split mechanics

The design is the deferred `typed_carrier_grid_carve.md` Phase A, now grounded against the live `operator.py`
(explorer Maps 1-2). **`from __future__ import annotations` is active (`operator.py:47`)** — all annotations are
strings, so subscript-arity changes never break at import time; only pyright evaluates them.

### §2.1 The typevars (`operator.py`, near the current `V` import at `:62`; `vector.py` owns `V`)
```python
# vector.py already defines: V = TypeVar("V", bound=Vector)        # invariant — KEEP (composer + bodies need it)
W    = TypeVar("W",    bound=Vector, default=V)                    # invariant codomain (Mixin)
Domain  = TypeVar("Domain",  bound=Vector, contravariant=True)          # Protocol input  (variant interface)
Codomain = TypeVar("Codomain", bound=Vector, covariant=True, default=Domain) # Protocol output (variant interface)
Cmid = TypeVar("Cmid", bound=Vector)                              # OperatorProduct intermediate
D2   = TypeVar("D2",   bound=Vector)                              # __matmul__ other-operand domain
```
- `V`/`W` invariant: the generic composer `+` bodies and the Mixin need invariance (a `LinearOperatorMixin[V, W]`
  is its own `apply` consumer/producer). `Domain`/`Codomain` variant: the *interface* Protocol is contra/co-variant for
  substitutability. (This Protocol-variant + Mixin-invariant split is the standard read/impl pattern.)
- Defaults (`W=V`, `Codomain=Domain`) are the PEP-696 mechanism making `[V]≡[V,V]` — the ~7 existing single-param sites
  (explorer Map 6: `iteration.py:415-416,680-681`, `integral_kernel_operator.py:162`,
  `multiplication_operator.py:115`, `sn/operator.py:223,631`) keep compiling unchanged.

### §2.2 `LinearOperator` Protocol (`operator.py:251-325`)
| Member | Before | After |
|---|---|---|
| class | `Protocol[V]` (`:251`) | `Protocol[Domain, Codomain]` |
| `apply` | `def apply(self, x: V, /) -> V` (`:312`) — **the falsehood** | `def apply(self, x: Domain, /) -> Codomain` |
| `domain`/`codomain` | `-> Optional[FunctionSpace]` (`:291,304`) | **UNCHANGED** (FunctionSpace level — see §2.6) |
| `capabilities` | `frozenset[str]` (`:289`) | unchanged |
- `@runtime_checkable` stays. `H`/`block_role` are NOT on the Protocol (only the Mixin) — unchanged.

### §2.3 `LinearOperatorMixin` (`operator.py:344-592`) — `Generic[V]` → `Generic[V, W]`
The TYPE_CHECKING `apply` stub (`:409-417`, `apply(self, x: V) -> V`) → `apply(self, x: V) -> W`. Thread EVERY
dunder (explorer Map 1 table) so none assumes endomorphism — **this must be COMPLETE; the prior +16-red attempt
was a HALF-threaded cascade** (`typed_carrier_grid` Phase A note):
| Dunder | Before (`[V]`) | After |
|---|---|---|
| `__add__`/`__radd__`/`__rsub__` (`:435,438,444`) | `(other: LinearOperator[V]) -> OperatorSum[V]` | `(other: LinearOperator[V, W]) -> OperatorSum[V, W]` |
| `__sub__` (`:441`) | `-> OperatorSum[V]` | `-> OperatorSum[V, W]` (body `OperatorSum(self, ScaledOperator(-1.0, other))` unchanged) |
| `__mul__`/`__rmul__`/`__neg__`/`__truediv__` (`:447,452,455,470`) | `-> ScaledOperator[V]` | `-> ScaledOperator[V, W]` |
| `__matmul__` (`:487`) | `(other: LinearOperator[V]) -> OperatorProduct[V]` | `(other: LinearOperator[D2, V]) -> OperatorProduct[D2, W]` (intermediate `V`=self.domain=other.codomain; `D2`=other.domain — **precise, no Any**; this IS the `R @ Λ @ M` chain) |
| `__and__`/`__rand__` (`:490,504`) | `-> TensorProductOperator` (bare) | **UNCHANGED** (tensor product is carrier-agnostic — §2.5) |
| `__pow__` (`:520`) | `-> LinearOperator[V]` | `-> LinearOperator[V, V]` (only endomorphic operators are powerable — assert/doc it) |
| `adjoint`/`H` (`:551,575`) | `-> LinearOperator[V]` | `-> LinearOperator[W, V]` (the adjoint SWAPS domain/codomain) |
- `block_role` (instance attr, NOT ClassVar — see the rationale comment `:380-407`): unchanged.

### §2.4 The composers (`operator.py`) → two-param

> **SPIKE-CORRECTED 2026-06-24** (`$CLAUDE_JOB_DIR/tmp/spike_two_param.py`, pyright 1.1.410 @ 3.13, `0 errors`).
> The plan's first draft (below, struck) put the *variant* `Domain`/`Codomain` on `OperatorProduct`'s Mixin base
> — a **variance clash**: pyright rejects passing a contravariant TypeVar as an argument to the *invariant*
> Mixin parameter (`reportInvalidTypeArguments`). The two load-bearing corrections:
> - **C1 — variance ONLY on the Protocol.** `Domain`/`Codomain` (contravariant/covariant) appear ONLY in the
>   `LinearOperator(Protocol[Domain, Codomain])` definition. EVERY impl class — the Mixin AND all composers — uses
>   the **invariant** `V`/`W` (+ the invariant `Cmid`/`D2`). This is the standard variant-read-Protocol /
>   invariant-impl split. `OperatorProduct` is `LinearOperatorMixin[V, W]` with `__init__(a: LinearOperator[Cmid, W],
>   b: LinearOperator[V, Cmid])` (input `V`, output `W`, intermediate `Cmid`) — NOT `[Domain, Codomain]`.
> - **C2 — `_AdjointOperator` needs explicit `Generic[V, W]`.** Its base is the SWAPPED `LinearOperatorMixin[W, V]`,
>   so its inferred param list is `[W, V]` — but `W` carries `default=V`, and PEP-696 forbids a defaulted
>   typevar before a non-defaulted one (and `W`'s default `V` would be out of scope). Declaring
>   `class _AdjointOperator(LinearOperatorMixin[W, V], Generic[V, W])` pins the param order to `[V, W]`
>   (non-defaulted `V` first) while the base binds them swapped. This is the ONLY composer with `[W, V]`
>   order — all others are `[V, W]` (correct: non-defaulted first), so no other composer needs the explicit
>   `Generic`.

| Class (line) | Before | After (spike-corrected) |
|---|---|---|
| `OperatorSum` (`:680`) | `LinearOperatorMixin[V]`, `__init__(a: LinearOperator[V], b: LinearOperator[V])` | `LinearOperatorMixin[V, W]`, `__init__(a: LinearOperator[V, W], b: LinearOperator[V, W])` (both operands SAME `[V,W]`; body `a.apply(x)+b.apply(x)` unchanged) |
| `OperatorProduct` (`:766`) | `LinearOperatorMixin[V]`, `__init__(a: LinearOperator[V], b: LinearOperator[V])` | **`LinearOperatorMixin[V, W]`** (C1 — invariant, NOT `[Domain, Codomain]`), `__init__(a: LinearOperator[Cmid, W], b: LinearOperator[V, Cmid])` (input `V`, output `W`, **intermediate `Cmid` captured at construction**; body `a.apply(b.apply(x))` unchanged) |
| `ScaledOperator` (`:839`) | `LinearOperatorMixin[V]`, `__init__(scalar, op: LinearOperator[V])` | `LinearOperatorMixin[V, W]`, `__init__(scalar, op: LinearOperator[V, W])` |
| `ZeroOperator` (`:919`) | `LinearOperatorMixin[V]` | `LinearOperatorMixin[V, W]` (a zero map can be non-endomorphic — the `(L+C-S-F)` zero summand) |
| `IdentityOperator` (`:897`) | `LinearOperatorMixin[V]` | `LinearOperatorMixin[V]` (≡`[V,V]` — identity IS endomorphic; keep single-param) |
| `_AdjointOperator` (`:600`) | `LinearOperatorMixin[V]`, `__init__(inner: LinearOperator[V])` | **`LinearOperatorMixin[W, V], Generic[V, W]`** (C2 — explicit `Generic` pins param order), `__init__(inner: LinearOperator[V, W])` (presents the SWAPPED `[W, V]`; already swaps domain/codomain at `:634-641`, metric-aware body `:643-670` unchanged) |

### §2.5 Bare carrier-agnostic primitives — STAY BARE (do NOT subscript)
Explorer Map 3: `PermutationOperator` (`:981`), `IncomingOrdinateMaskTensor` (`:1085`), `PeriodicWrapOperator`
(`:1175`), `TensorProductOperator` (`:1223`), `SumOfTensorProductsOperator` (`:1344`), `DiagonalOperator`
(`:1437`), `RankOneOperator` (`:1660`). These act on **one tagged axis of a bare `ndarray`** and span multiple
carrier worlds (e.g. `TensorProductOperator` on the boundary trace AND the bulk group array) — by the
type-vs-property mint-iff rule they MUST stay unsubscripted (`→ [Unknown, Unknown]`, correct for an
implementation-layer primitive). The carrier role is supplied at the realisation seam (the `_as_boundary`
instance-stamp precedent). **`MultiplicationOperator` is NOT bare** (explorer correction): it is
`LinearOperatorMixin["FullField"]` (`multiplication_operator.py:115`), a typed-`FullField` flux→source operator
(its docstring `:185-241`: "The codomain is a *source*") — it belongs in the typed group (§3) with `CollisionOperator`.

### §2.6 THE TWO LEVELS — carrier type vs FunctionSpace (do not conflate)
There are TWO orthogonal typings, and P4.5 only adds the first:
- **Carrier type** (`Domain`/`Codomain`, `Vector`-bound): `AngularFlux`, `AngularSourceSink`, `TimedFullField`,
  `FullField`, … — the static `apply` in/out. **This is what P4.5 adds.**
- **FunctionSpace** (`domain`/`codomain`, runtime): the shape/metric object the `OperatorProduct`/`OperatorSum`
  **runtime guards** (`operator.py:796-806`, `:714-732`) check. **Unchanged by P4.5.**

These are complementary. A loss operator is carrier `TimedFullField → FullField` (role-changing) yet
FunctionSpace `full_field_space → full_field_space` (the flux and source composites SHARE the underlying space
structure — same shape/metric). So **a role-changing operator can be space-endomorphic.** P4.5 types the
carrier; the space guard is independent and already correct. **Do not try to make the space guard role-aware** —
that conflates the two levels.

---

## §3. THE CONSUMERS — activate honest types + retire the confessions

### §3.1 `ScatteringOperator` / `FissionOperator` — retire the `@overload` confessions
Both inherit bare `LinearOperatorMixin`, carry `None` spaces, and confess via `@overload`
(`scattering.py:1388-1406`, `fission.py:466-480`). After the split:
- Subscript them to their **primary balance carrier pair**. The balance matvec carrier is the composite:
  `ScatteringOperator(LinearOperatorMixin["TimedFullField", "FullField"])`, same for `FissionOperator`. The base
  `apply(TimedFullField) -> FullField` is then HONEST — the endomorphism-confession comment retires.
- **OPEN DESIGN DECISION (resolve in impl, §6 Q1):** the *multi-arm* dispatch (`ScalarFlux→ScalarSourceSink`,
  `AngularFlux→AngularSourceSink`, `HarmonicMomentFlux→AngularSourceSink`) is a SEPARATE concern from the
  endomorphism lie. Options: (a) keep a *narrowed* `@overload` set for the convenience arms (now WITHOUT the
  "not endomorphic" confession — the base is honest, the overloads are genuine multi-input polymorphism);
  (b) a type-level `Flux→Source` mapping so `apply` is generic (heavier); (c) split the convenience arms into
  named verbs. Lean (a) — the overloads stop being a *confession* and become honest multi-dispatch declarations.
- **CRITICAL RISK (explorer flag #2):** S/F currently carry `None` spaces, so the `OperatorSum`/`OperatorProduct`
  runtime guards are **no-ops for them**. The balance `(L+C) − S − F` sums them today only because all are `None`.
  Two sub-risks: (i) if you give S/F real spaces, the guards ACTIVATE — verify `(L+C)`, `S`, `F` all share the
  composite `full_field_space` so `OperatorSum` (equal domains AND codomains) PASSES; (ii) the carrier subscript
  `["TimedFullField", "FullField"]` must match `(L+C)`'s — confirm `StreamingOperator`/`InvertibleOperator` are
  the SAME `["TimedFullField", "FullField"]` (today they are `["FullField"]` ≡ endomorphic — see §3.2). **The
  balance `OperatorSum` operands must be uniformly `[TimedFullField, FullField]` or the split breaks composition.**

### §3.2 The loss operators `L` / `C` / `(L+C)` (`sn/operator.py`, `transport/multiplication_operator.py`)
Explorer Map 3 (loss): `StreamingOperator(LinearOperatorMixin["FullField"])` (`:223`),
`InvertibleOperator(OperatorSum["FullField"])` (`:631`), `CollisionOperator(MultiplicationOperator)`
(`MultiplicationOperator(LinearOperatorMixin["FullField"])`). Today these are `["FullField"]` ≡ `[FullField, FullField]`
(endomorphic on the composite). But the loss ACTION is `TimedFullField → FullField` (flux composite in, source
composite out): `loss_action` returns `FullField(bulk=AngularSourceSink, boundary=BoundarySourceSink)`. So the
HONEST loss type is `[TimedFullField, FullField]`, NOT `[FullField, FullField]`. **Decision (§6 Q2):** retype the
loss leaves to `[TimedFullField, FullField]` to match S/F, so the balance `OperatorSum` is uniformly
`[TimedFullField, FullField]`. NB the `.solve` resolvent FLIPS the role back (`apply`→sink, `solve`→flux:
`solve → TimedFullField`) — `solve`'s return type is the INVERSE direction `[FullField, TimedFullField]`-ish;
type `solve` separately from `apply` (it already is a distinct method).

### §3.3 The carrier-typed kernel — DISSOLVE the P4 option-2 asymmetry
Today `kernel = frame.conjugate(Λ) = OperatorProduct(R, OperatorProduct(Λ, M))` is an **ndarray composition**
(explorer Map 5): M/R (`_FrameAnalysis`/`_FrameReconstruction`) declare distinct FunctionSpaces but
`apply(NDArray)→NDArray`; Λ has a typed arm but is invoked on `.values` inside the product. The flux→source
carrier typing lives ONLY at the `S.apply` boundary. To make `kernel.apply(AngularFlux)→AngularSourceSink`
*directly* (so full-angular and windowed are uniformly typed):
- **Transport-level typed face operators.** The numerics frame faces CANNOT return transport carriers (layering —
  numerics is below transport; the P4 ruling). So the typed M/R operators live in **transport** (on/with
  `HarmonicFrame`): expose `HarmonicFrame.analysis_operator : LinearOperator[AngularFlux, HarmonicMomentFlux]` and
  `.reconstruction_operator : LinearOperator[HarmonicMomentSourceSink, AngularSourceSink]` — thin typed wrappers
  that unwrap `.values`, delegate to the generic numerics ndarray face, re-wrap. (The P4 `HarmonicFrame.analyse`/
  `reconstruct` METHODS already do exactly this cast; lift the cast into composable LinearOperator objects.)
- **Λ already has the typed arm** (`HarmonicMomentFlux → HarmonicMomentSourceSink`, P4 `ed1e14d`). With the
  two-param split it is honestly `LinearOperator[HarmonicMomentFlux, HarmonicMomentSourceSink]` (drop the
  `np.ndarray | …` union arm OR keep it as a separate ndarray-composition path — §6 Q3).
- **The typed kernel** = `R_typed @ Λ_typed @ M_typed : LinearOperator[AngularFlux, AngularSourceSink]` (the
  two-param `__matmul__` threads `AngularFlux → HarmonicMomentFlux → HarmonicMomentSourceSink → AngularSourceSink`
  precisely). Then: full-angular `build_aniso_source(ψ) = (1/W)·kernel.apply(ψ)` (typed, no `.values`); windowed
  = the `R_typed @ Λ_typed` sub-product (typed, M skipped). **Both uniformly typed; the asymmetry dissolves.**
- The §5.6 `IntegralKernelOperator` abstraction is PRESERVED (the kernel is still ONE composed operator, now
  carrier-typed). The 0-ULP `test_scattering_kernel_crosscheck` rewires to the typed kernel (its
  `_aniso_source_from_moment_values` ndarray oracle stays as the structurally-independent reference — verify the
  typed kernel == the oracle bit-for-bit).
- **Generalises to fission for free:** `F = M_χ @ ProductionRate @ M_νΣf` (rank-1 middle) types the same way —
  the same frame machinery covers both (cross-domain-attacker Q4). Fission's carrier-typed kernel is a follow-on,
  not blocking, but the design must NOT preclude it.

### §3.4 The `projection.py` ABCs + frame faces — the second confession
`AnalysisOperator`/`ReconstructionOperator` (`projection.py:90,120`) are bare `LinearOperatorMixin` with
`apply(ndarray)→ndarray` yet declare distinct domain/codomain spaces — the SAME confession shape as S/F, without
the overload stack. Subscript them: `AnalysisOperator[MeasureSpaceCarrier, TestSpaceCarrier]`,
`ReconstructionOperator[BasisSpaceCarrier, MeasureSpaceCarrier]`. **CAVEAT:** these are NUMERICS and operate on
ndarray — their `Domain`/`Codomain` are the ndarray/`Vector` level, NOT transport carriers (layering). The transport
carrier typing is the §3.3 transport wrappers. Keep the numerics ABCs ndarray-bound; do not leak transport types
into numerics.

---

## §4. SUB-PHASES (each: green + CLI-pyright ≤ baseline (target: DOWN) + runtime-zero; ff-merge)

- **P4.5a — the gating decision + the numerics split foundation.** Resolve §1 (AskUserQuestion: 3.13-floor vs
  `typing_extensions`). Add the typevars (§2.1); retype the Protocol (§2.2), the Mixin + all dunders (§2.3), the
  composers (§2.4); rewrite the `vector.py` `V` "endomorphism-is-honest" prose (`:144-151`) — it becomes a twin
  source of truth against the new Protocol if left (Cardinal Rule 2 / explorer flag #5). Bare primitives untouched
  (§2.5). **NO consumer yet — gated additive, lands WITH P4.5b** (do not merge alone). Gate: full
  operator/realizer/iteration/GAP nets green under `-O`; **zero runtime change**; pyright **no new reds** (the
  precise complete threading clears the half-threaded-cascade reds).
- **P4.5b — activate S/F/L honest types + retire the `@overload` confessions.** Retype `ScatteringOperator`/
  `FissionOperator` (§3.1) and the loss leaves (§3.2) to `[TimedFullField, FullField]`; resolve the multi-arm
  decision (§6 Q1); VERIFY the balance `OperatorSum` composition guards activate cleanly (§3.1 risk). Gate: the
  balance `(L+C−S−F−B)` matvec + SI + Krylov + eigenvalue nets green; pyright DOWN (confessions retired);
  0-ULP canary green.
- **P4.5c — the carrier-typed kernel (dissolve the option-2 asymmetry).** §3.3: the transport typed face
  operators on `HarmonicFrame`, the typed `kernel`, retype `build_aniso_source` + the windowed arm to the typed
  kernel/sub-product, rewire the 0-ULP crosscheck to the typed kernel (oracle preserved). Gate: 0-ULP crosscheck
  bit-identical to the ndarray oracle; Phase-5a guard green; the Mode-11 windowed sentinel green; full-angular ==
  windowed symmetric (both typed).
- **P4.5d — projection ABCs + frame faces** (§3.4): subscript the numerics ABCs (ndarray-bound). Gate: frame/
  numerics nets green; pyright DOWN.
- **P4.5e — docs.** `operator_algebra.rst`: the `Operator[Flux, SourceSink]` contract as the uniform balance
  algebra (affine-map / linear-part-of-affine-map structure), the `[Domain, Codomain]` split, the carrier-vs-space two
  levels, fission-as-rank-1-degenerate. Update the P4 carrier-grid section (the kernel is now typed). Rewrite any
  `vector.py`/operator docstrings asserting endomorphism. Sphinx `-W` clean. (Archivist.)

---

## §5. VERIFICATION
- `.venv/bin/python -O -m pytest tests/numerics tests/sn/operators tests/sn/solve tests/sn/sweep
  tests/sn/eigenvalue tests/sn/verification -q -rfE --timeout=300 -p no:cacheprovider` — SERIAL (xdist unstable);
  the **0-ULP `test_scattering_kernel_crosscheck`** is the load-bearing canary; the **7-and-only-7 pre-existing
  reds** baseline. P4.5a/b/d are RUNTIME-ZERO (pure typing + carrier wraps) — any value change = STOP. P4.5c is
  bit-identical (the typed kernel == the ndarray oracle; same Λ kernel + R face).
- `npx --no-install pyright --outputjson orpheus/` — the ORACLE (NOT the streamed `<new-diagnostics>` LSP lag,
  #226). Target: net reds **DOWN** (the `@overload` confessions + the `None`-space honesty retire; the
  half-threaded-cascade reds do NOT appear because the threading is complete). Pin `[tool.pyright] pythonVersion`
  to the chosen floor (§1).
- **Proactive `test-architect` before P4.5b and P4.5c** (operator-algebra carves crossing the 0-ULP path). The
  intrinsic-property tests: a non-endomorphic operator's `apply` carrier type is statically `Codomain≠Domain` (an
  `assert_type` pin); the composition guard activates for S/F (a negative test: composing incompatible-space
  operators raises `IncompatibleOperatorComposition`).

---

## §6. OPEN DESIGN QUESTIONS (resolve in-session, do NOT pre-bake)
- **Q1 — multi-arm dispatch vs `@overload`.** After the base is honest `[TimedFullField, FullField]`, do the
  convenience arms (`ScalarFlux`/`AngularFlux`/`HarmonicMomentFlux`) stay as a *narrowed* honest `@overload` set,
  or get a cleaner mechanism? Lean: narrowed overloads (no longer a confession). §3.1.
- **Q2 — loss carrier type.** Retype the loss leaves `["FullField"]` → `["TimedFullField", "FullField"]` to match
  S/F (so the balance OperatorSum is uniform)? Confirm `(L+C).apply` is genuinely `TimedFullField→FullField` at
  every call site, and `solve` (the inverse direction) types separately. §3.2.
- **Q3 — Λ's ndarray arm.** With Λ honestly `[HarmonicMomentFlux, HarmonicMomentSourceSink]`, keep its
  `np.ndarray | …` union arm (the endomorphic ndarray view the OLD ndarray kernel composed on) or retire it once
  the typed kernel replaces the ndarray composition? §3.3.
- **Q4 — does `ScatteringOperator`/`FissionOperator` carry FunctionSpaces now?** Today `None` (the guard is a
  no-op). Activating real composite spaces turns on the runtime guard — desired (catches mis-composition) but
  must be verified non-breaking across EVERY existing composition. §3.1 risk.
- **Q5 — fission's carrier-typed kernel.** Fold `F = M_χ @ ProductionRate @ M_νΣf` typing into P4.5c, or a
  follow-on? (Design must not preclude; implementation may defer.)

---

## §7. SCOPE / DISCIPLINE / CRITICAL FILES
- Surgical, main-agent-direct, user-steered; NO `method-implementer`. Per-sub-phase ff-merge; `main` always green.
  P4.5a is additive-only (lands with 4.5b). The pyright ratchet target FLIPS from "guide not gate" (the P1-P4
  campaign deferral) to **"net DOWN is the deliverable"** — P4.5 is precisely the #226 generic-Protocol cleanup
  the campaign deferred.
- **Critical files** (explorer maps; re-derive line numbers via Nexus): `orpheus/numerics/operator.py` (the split
  core — Protocol, Mixin, composers), `orpheus/numerics/vector.py` (typevars + the prose rewrite),
  `orpheus/numerics/projection.py` (the ABCs), `orpheus/numerics/frame.py` (the faces),
  `orpheus/transport/frames/harmonic_frame.py` (the typed face operators — NEW), `orpheus/sn/scattering.py` +
  `orpheus/sn/fission.py` (retire the @overload confessions + the carrier-typed kernels),
  `orpheus/sn/operator.py` + `orpheus/transport/multiplication_operator.py` (the loss leaves),
  `orpheus/transport/integral_kernel_operator.py` (the IntegralKernelOperator Protocol — refine to `[Domain, Codomain]`),
  `orpheus/numerics/iteration.py` (the SI/Krylov drivers — the variadic `*gains: LinearOperator[V, W]`),
  `pyproject.toml` (the §1 gating decision). Memos:
  `.claude/agent-memory/cross-domain-attacker/flux_to_sourcesink_operator_contract_frames.md`.
- **Durable principle banked:** an operator is the *linear part of an affine map* `Flux → SourceSink` (affine
  point domain, vector codomain, the cross-section unit-gain as fiber transition); the balance
  `(L+C−S−F)ψ=q` is a typed source/sink algebra; `[Domain, Codomain]` IS that type; `@overload` endomorphism stacks are
  the symptom of the missing second parameter.
