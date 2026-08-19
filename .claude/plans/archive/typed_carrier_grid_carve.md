# Typed-carrier operator carve — complete the (angular ⊗ moment) × (flux ⊗ source) grid

> **⚠️ SUPERSEDED (2026-06-23) by `.claude/plans/archive/frame_basis_carve.md`.** A deep
> architecture discussion reframed this: the grid's friction was a *symptom* of the
> numerics/transport partition + a missing `Frame` abstraction. The active carve builds a
> generic `Frame` + `Basis` ABC (the discrete-frame `analysis`/`reconstruction` pair).
> Retained for history only — do NOT execute this plan.

> **APPROVED PLAN (user-approved 2026-06-22).** Cold-pickup-ready: a fresh
> post-compaction session runs this end-to-end. Read §0 first, then the body.

---

## §0. Pickup-cold orientation (read FIRST)

**What this is.** Make the SN operators "speak typed carriers" (typed `Field`
Vectors) instead of raw `np.ndarray`, by completing the 2×2 carrier grid
`(angular ⊗ moment) × (flux ⊗ source)`. This is the user's deep-architecture
answer to "why do we need ndarray?" — and it gives the two-param domain/codomain
operator split its first real consumer (`HarmonicMomentProjection` /
`…Reconstruction` are genuinely non-endomorphic).

**Surgical / main-agent-direct / user-steered** (`.claude/rules/delegation.md`):
the main agent writes the code with the user steering; **do NOT dispatch
`method-implementer`**. `test-architect` / `explorer` / `qa` /
`elegance-enforcer` remain available.

**Environment.** Host → `.venv/bin/python` (Python **3.14** — PEP 696 TypeVar
defaults available). Temp files in `$CLAUDE_JOB_DIR/tmp`. Run the SessionStart
protocol (vv-principles, subagent-handoff-protocol, coding-elegance skills;
`.claude/lessons.md`; nexus `session_briefing`; `docs/development.rst`).

**Branch & state (reconcile against git — trust git, not this note).**
- Branch `refactor/operator-inverse-algebra` (off `main`). `git log --oneline -1`.
- **Landed on the branch:** `ef0631a` (the *original* inverse-as-operator carve
  plan), `ab75d02` (Phase 0 safety net — GAP-1/2/3 tests + mutation recipes),
  `272caa9` (**the foundation: RC1 `block_role` ClassVar→instance + RC3
  `_as_boundary`/`realize` retyped to `LinearOperatorMixin` + `α·B` dunder form;
  −12 pyright reds 420→408, zero runtime**). Working tree clean.
- This carve **replaces** the original `operator_inverse_algebra_carve.md` Phase 1
  (two-param split, now Phase A here, justified by the grid). The
  inverse-as-operator work (`.solve`→`.inverse`, retire `CAP_SOLVE`) is a SEPARATE
  follow-on that builds on Phase A.

**Oracle & constraints.** `npx --no-install pyright --outputjson orpheus/` is the
error ORACLE — the streamed `<new-diagnostics>` is #226 langserver LAG (import-
resolution cascades on hub-module edits); **never act on the stream, verify with a
CLI run.** **NO `# type: ignore`.** Canonical tests `python -O -m pytest`. Per-phase
ff-merge+push when green (user-authorized); `main` always green.

**Recommended execution order:**
1. **FIRST: dispatch `test-architect`** (proactive trigger — operator-algebra carve
   crossing subsystems). Brief it for: the verification plan for the three typed
   edges (M/Λ/R), the new `HarmonicMomentSourceSink` intrinsic-property test, and
   which existing gates pin the contracts (Galerkin `Π R = 4π·I`, the scattering
   0-ULP crosscheck, keff/MMS). It shapes the carve.
2. **B1 — the `HarmonicMoment*` renames** (behavior-preserving, lowest-risk; clean
   first commit, can ff-merge alone). Use Nexus `rename`/`impact`.
3. **Phase A — two-param operator foundation** (the composer cascade; iterate vs the
   CLI oracle — the earlier +16-red attempt was a HALF-threaded cascade; §A gives
   the PRECISE threading that clears it).
4. **B2 / B3 / B4** — type the edges, the kernel `(1/W)*(R@Λ@M)`, retire the casts,
   migrate the tests.
5. Verify (pyright down + net green) → re-baseline ratchet → Sphinx → **EXIT to the
   #226 pyright burndown** (`pyright_signal_cleanup.md`).
   Note: Phase A has no consumer alone — **do not ff-merge A without B**.

**Reading list:** this plan (§1–end); the three explorer maps' findings are folded
in below (carrier grid, the carrier-agnostic-primitive verdict, the boundary
raw-axis verdict); `operator_inverse_algebra_carve.md` (the inverse follow-on);
`.claude/plans/archive/issue_226_b4_operator_generics_verification.md` (the Phase-0 net).
Memory: `[[feedback-naming-consistency-greppable]]` (the family-naming principle
that drives the renames).

---

## Context

The #226 B4 work surfaced a deeper architectural question the user named directly:
**"operators should speak typed carriers, not `np.ndarray`."** Investigation (3
explorer maps + a CLI linchpin spike) established:

- **Typed Field carriers satisfy the `Vector` bound; raw `np.ndarray` does NOT** (proven:
  `AngularFlux`/`ScalarFlux`/`HarmonicMomentField`/`TimedFullField` → 0 pyright errors;
  `np.ndarray` → "incompatible with protocol Vector" because numpy's overloaded `__add__`
  can't match a `Self`-typed protocol). So typing operators on **typed carriers** sidesteps
  the bound friction entirely — the friction only ever came from trying to type the generic
  primitives on `np.ndarray`.
- **The genuine non-endomorphic operators are exactly `HarmonicMomentProjection` (Π, renamed
  from `MomentProjection` — see B1) and `HarmonicMomentReconstruction` (R)** — they map
  angular ↔ moment (different `Field` types). This is *why* the `cast(LinearOperator, …)`
  workaround exists at
  `scattering.py:648-660`: the single-`V` endomorphism `apply(x: V) -> V` cannot express
  them. They are the **first real consumer** of a domain/codomain split — it is no longer
  speculative.
- **The current type system is asymmetric:** angular space has role leaves
  (`AngularFlux`/`AngularSourceSink`), moment space has only one (`HarmonicMomentField`).
  Scattering's anisotropic in-scatter source `S = (1/W)·R·Λ·M·ψ` therefore can't be typed as
  a clean path — the flux→source role change has nowhere to live and leaks to the consumer.

**The design (user-confirmed): complete the 2×2 carrier grid.** Make every transformation a
typed edge that changes *either* the axis (angular↔moment) *or* the role (flux→source),
never both. The role change then lives exactly on Λ (scattering), where the physics puts it.

```
                M : AngularFlux ───────────────▶ HarmonicMomentFlux
   (project, role-preserving, Galerkin Π)                │
                                                          │ Λ : scatter, role-CHANGING,
                                                          │     axis-preserving  (Σ_{s,ℓ})
                                                          ▼
        AngularSourceSink ◀──────────────── HarmonicMomentSourceSink
                R : (reconstruct, role-preserving, Galerkin R)

   kernel  S = (1/W) * (R @ Λ @ M)  :  AngularFlux → AngularSourceSink
```

**Locked decisions:**
1. **1/W** is an explicit `(1/W) *` `ScaledOperator` factor in the composed kernel — out of
   the consumer, M/R stay the pure Galerkin pair (`Π R = 4π·I` intact), reads like the math.
2. **R's codomain is `AngularSourceSink`** — honest, because its *input* is
   `HarmonicMomentSourceSink` (Λ produced it) and R preserves role. The flux→source change is
   Λ's job, on its own typed edge — no casts, no implicit steps.
3. **Naming (Domain-first, matching every existing carrier):** rename
   `HarmonicMomentField` → **`HarmonicMomentFlux`** and `MomentProjection` →
   **`HarmonicMomentProjection`**; add **`HarmonicMomentSourceSink`**. Makes `grep Flux` /
   `grep SourceSink` / `grep HarmonicMoment` find every member (today `HarmonicMomentField`
   and `MomentProjection` are off-pattern). See `[[feedback-naming-consistency-greppable]]`.

## What stays UNTOUCHED (validated — do NOT type these)

- **The carrier-agnostic array primitives** (`TensorProductOperator`,
  `SumOfTensorProductsOperator`, `DiagonalOperator`, `RankOneOperator`,
  `PermutationOperator`, `IncomingOrdinateMaskTensor`, `PeriodicWrapOperator`). Confirmed they
  act on **one tagged axis of a bare ndarray** and span **multiple non-isomorphic carrier
  worlds** (`TensorProductOperator` is used on the boundary trace AND the bulk group array).
  By the project's type-vs-property mint-iff rule they MUST stay carrier-agnostic; the carrier
  role is supplied at the realization seam (the `_as_boundary` instance-stamp precedent). They
  stay bare `LinearOperatorMixin` (→ `[Unknown, Unknown]`, which is correct for an
  implementation-layer primitive).
- **The boundary operators.** The realized BC laws act on a raw per-face slice (`face_view`),
  with the typed `BoundaryFlux` one layer up at `SNBoundaryOperator`'s surface and the single
  `_reflect_trace` seam between — a genuine raw-axis layer, not a leak.

## Phase A — two-param operator foundation (the split, now justified)

`orpheus/numerics/operator.py` + `orpheus/numerics/vector.py`.

- `V` stays `bound=Vector` invariant (load-bearing — the generic composer `+` bodies need it).
- Add typevars: `W = TypeVar("W", bound=Vector, default=V)` (invariant codomain);
  `Din = TypeVar("Din", bound=Vector, contravariant=True)`,
  `Cout = TypeVar("Cout", bound=Vector, covariant=True, default=Din)` (variant interface);
  `Cmid = TypeVar("Cmid", bound=Vector)`, `D2 = TypeVar("D2", bound=Vector)` (composition
  intermediates). PEP 696 defaults (Python 3.14) keep endomorphic subscripts single-param —
  `LinearOperator[V]` ≡ `[V, V]`, `LinearOperatorMixin["FullField"]` ≡ `[FullField,
  FullField]` — so **every existing single-param subscript keeps working unchanged.**
- `LinearOperator(Protocol[Din, Cout])`, `apply(self, x: Din) -> Cout`,
  `H -> LinearOperator[Cout, Din]` (the adjoint swaps; trivial when endomorphic).
- `LinearOperatorMixin(Generic[V, W])`, `apply(self, x: V) -> W`, `adjoint/H ->
  LinearOperator[W, V]`. Thread the dunders so they no longer assume endomorphism:
  - `__add__/__sub__/__radd__/__rsub__(self, other: LinearOperator[V, W]) -> OperatorSum[V, W]`
  - `__mul__/__rmul__/__neg__/__truediv__ -> ScaledOperator[V, W]`
  - `__matmul__(self, other: LinearOperator[D2, V]) -> OperatorProduct[D2, W]` (the
    intermediate `V` = self's domain = other's codomain; `D2` = other's domain — **precise,
    no `Any`**; this is exactly the `R @ Λ @ M` chain).
- Composers → two-param: `OperatorSum(LinearOperatorMixin[V, W])` with
  `__init__(a: LinearOperator[V, W], b: LinearOperator[V, W])`;
  `ScaledOperator(LinearOperatorMixin[V, W])`; `ZeroOperator(LinearOperatorMixin[V, W])`;
  `OperatorProduct(LinearOperatorMixin[Din, Cout])` with `__init__(a: LinearOperator[Cmid,
  Cout], b: LinearOperator[Din, Cmid])` (intermediate `Cmid` captured at construction —
  precise); `IdentityOperator(LinearOperatorMixin[V])` (endomorphic); `_AdjointOperator`
  presents `[W, V]`.
- Bare primitives stay bare (unsubscripted) — carrier-agnostic.
- **Gate:** full operator/realizer/iteration net + GAP gates green under `-O`; **zero runtime
  change**; CLI pyright shows **no new reds** (the earlier +16 were from a half-threaded
  cascade — the precise threading above clears them). This phase has NO consumer yet; it is
  gated as additive and lands with Phase B (do not merge A alone — it is justified by B).

## Phase B — complete the grid + type the projection/scattering edges

### B1 — Carrier types + `HarmonicMoment*` family renames (behavior-preserving, first commit)
Two consistency renames + one new leaf. The whole moment/harmonic family becomes uniformly
`HarmonicMoment<Thing>`, so `grep HarmonicMoment` / `grep Flux` / `grep SourceSink` find every
member. Use Nexus `rename`/`impact` for each reference set; all are pure renames (no behavior
change) — land as a clean first commit.
- **Rename carrier** `HarmonicMomentField` → `HarmonicMomentFlux` (it is `(FluxRole,
  MomentField)`; the rename only makes the role token explicit). Wide reference set (the
  `… | None` threading + the in-sweep emit read it).
- **Rename operator** `MomentProjection` → `HarmonicMomentProjection` — restores symmetry with
  its dual `HarmonicMomentReconstruction` (today `MomentProjection` is the lone `HarmonicMoment`
  family member missing the qualifier). `grep HarmonicMomentProjection`/`…Reconstruction` now
  finds the Π/R pair as a matched set.
- **Add** `HarmonicMomentSourceSink(MomentField)` in `source_sinks/` mirroring
  `AngularSourceSink(AngularField)` (a `MomentField` source/sink leaf, no `FluxRole`).
  Ships its intrinsic-property test (per `feedback_test_intrinsic_properties`).

### B2 — Type the three edges (`orpheus/numerics/projection.py`, `orpheus/sn/scattering.py`)
- `HarmonicMomentProjection.apply(x: AngularFlux) -> HarmonicMomentFlux` — einsum body
  unchanged, `.values`-unwrap internally, re-wrap output via the input carrier's mesh/`L`.
  (`domain`/`codomain` already present.)
- `HarmonicMomentReconstruction.apply(x: HarmonicMomentSourceSink) -> AngularSourceSink` —
  same; **add** `domain`/`codomain` properties (currently absent — greenfield).
- `LegendreMomentScattering.apply(x: HarmonicMomentFlux) -> HarmonicMomentSourceSink` — the
  sole role-changing edge; the existing typed arm already wraps/unwraps, retype it to the
  role-typed in/out (the `np.ndarray` raw arm retires or stays as the serialization escape).

### B3 — Kernel + consumer migration + cast retirement (`orpheus/sn/scattering.py`)
- `kernel = (1/W) * (R @ Λ @ M)` typed `LinearOperator[AngularFlux, AngularSourceSink]`.
- **Retire** the `cast(LinearOperator, …)` workarounds (`scattering.py:648-660`) — the
  `OperatorProduct` now types end-to-end (Phase A's precise `__matmul__`/`__init__`).
- Migrate the 2 production `.apply` sites: `build_aniso_source` (`scattering.py:951` — pass
  the typed `AngularFlux`, not `.values`) and `_aniso_source_from_moment_values`
  (`scattering.py:540-541` — typed moment in/out); the manual `/sum_w` (`scattering.py:953/
  957`) folds into the kernel's `(1/W) *`.

### B4 — Test migration (verification net)
- `tests/numerics/test_projection_operators.py` (shape/ABC tests pass **raw ndarrays** → must
  build typed `AngularFlux`/`HarmonicMomentFlux` first — behavioral contract migration, per
  the retirement/test-migration rule).
- Keep green: `tests/numerics/test_spherical_harmonic_space.py` (Galerkin `Π R = 4π·I`,
  `Π* = g_C·S₀`), `tests/sn/operators/test_scattering_kernel_crosscheck.py` (0-ULP), the
  keff/MMS L1 net.

## Verification (end-to-end)

1. `npx --no-install pyright --outputjson orpheus/` — net reds **down** (cast retirements +
   the projection family typed), **no new reds**; re-baseline the ratchet.
2. `.venv/bin/python -O -m pytest tests/numerics tests/sn/operators tests/sn/eigenvalue
   tests/sn/verification` — Galerkin invariants, scattering 0-ULP crosscheck, keff (≥2G
   heterogeneous), MMS slab/2d all green. Tolerances unchanged (re-baseline only if a
   principled FP-reassociation is forced, per `vv-principles`).
3. The grid is real: `grep -rn "class HarmonicMomentFlux\|class HarmonicMomentSourceSink"`
   confirms the symmetric leaves; `grep Flux` / `grep SourceSink` now find every carrier.
4. Sphinx rebuild clean (Cardinal Rule 3) — update the operator-algebra + scattering theory
   pages with the grid (the 4-node carrier diagram, M/Λ/R as role-preserving/role-changing
   edges, the `(1/W)·R·Λ·M` kernel).

## Relationship to the broader carve

- This **replaces** the original `operator_inverse_algebra_carve.md` Phase 1 (the two-param
  split is now Phase A here, with a real consumer).
- The **inverse-as-operator** work (`.inverse()` returns an operator, `.solve`→`.inverse().
  apply`, retire `CAP_SOLVE`/`MissingCapability` solve-axis; the original detour's Phases
  2-5) is a SEPARATE follow-on that builds on Phase A's two-param `.H`/`.inverse` swap types.
- **EXIT back to the #226 pyright burndown** (per `pyright_signal_cleanup.md`) after the
  ratchet re-baseline: resume B5 (scattering/solver union dispatch), units.py, B6.
