# Phase 4.5 — the two-param `LinearOperator[Domain, Codomain]` split

> **Part of the Frame projection/reconstruction machinery campaign** (`.claude/plans/frame_projection_machinery.md`).
> Inserted as **P4.5**, to run BEFORE P5 (energy condensation) / P6 (adjoint homogenisation): the user's
> directive is "polish the architecture so we don't blindly build P5/P6 on an operator algebra that must
> further change." Branch `refactor/operator-inverse-algebra`. Surgical, main-agent-direct, user-steered;
> NO `method-implementer`. Host env → `.venv/bin/python` (runtime 3.14). **Cold-pickup-ready** (written
> because the session compacts before implementation — read this end-to-end, then `operator.py`/`vector.py`).

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

**P4.5 makes the operator type honest: `LinearOperator[Din, Cout]`, `apply(x: Din) -> Cout`.** This:
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
`[Din, Cout]` IS that affine-map type. (Memo: `.claude/agent-memory/cross-domain-attacker/flux_to_sourcesink_operator_contract_frames.md`.)

---

## §1. THE GATING DECISION — PEP-696 availability (resolve FIRST)

`LinearOperator[V] ≡ LinearOperator[V, V]` (keeping all ~7 existing single-param subscript sites working) needs
a **PEP-696 TypeVar default** `Cout = TypeVar("Cout", default=Din)`. Native `typing.TypeVar(default=...)` exists
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
Din  = TypeVar("Din",  bound=Vector, contravariant=True)          # Protocol input  (variant interface)
Cout = TypeVar("Cout", bound=Vector, covariant=True, default=Din) # Protocol output (variant interface)
Cmid = TypeVar("Cmid", bound=Vector)                              # OperatorProduct intermediate
D2   = TypeVar("D2",   bound=Vector)                              # __matmul__ other-operand domain
```
- `V`/`W` invariant: the generic composer `+` bodies and the Mixin need invariance (a `LinearOperatorMixin[V, W]`
  is its own `apply` consumer/producer). `Din`/`Cout` variant: the *interface* Protocol is contra/co-variant for
  substitutability. (This Protocol-variant + Mixin-invariant split is the standard read/impl pattern.)
- Defaults (`W=V`, `Cout=Din`) are the PEP-696 mechanism making `[V]≡[V,V]` — the ~7 existing single-param sites
  (explorer Map 6: `iteration.py:415-416,680-681`, `integral_kernel_operator.py:162`,
  `multiplication_operator.py:115`, `sn/operator.py:223,631`) keep compiling unchanged.

### §2.2 `LinearOperator` Protocol (`operator.py:251-325`)
| Member | Before | After |
|---|---|---|
| class | `Protocol[V]` (`:251`) | `Protocol[Din, Cout]` |
| `apply` | `def apply(self, x: V, /) -> V` (`:312`) — **the falsehood** | `def apply(self, x: Din, /) -> Cout` |
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
| Class (line) | Before | After |
|---|---|---|
| `OperatorSum` (`:680`) | `LinearOperatorMixin[V]`, `__init__(a: LinearOperator[V], b: LinearOperator[V])` | `LinearOperatorMixin[V, W]`, `__init__(a: LinearOperator[V, W], b: LinearOperator[V, W])` (both operands SAME `[V,W]`; body `a.apply(x)+b.apply(x)` unchanged) |
| `OperatorProduct` (`:766`) | `LinearOperatorMixin[V]`, `__init__(a: LinearOperator[V], b: LinearOperator[V])` | `LinearOperatorMixin[Din, Cout]`, `__init__(a: LinearOperator[Cmid, Cout], b: LinearOperator[Din, Cmid])` (**intermediate `Cmid` captured at construction** — precise; body `a.apply(b.apply(x))` unchanged) |
| `ScaledOperator` (`:839`) | `LinearOperatorMixin[V]`, `__init__(scalar, op: LinearOperator[V])` | `LinearOperatorMixin[V, W]`, `__init__(scalar, op: LinearOperator[V, W])` |
| `ZeroOperator` (`:919`) | `LinearOperatorMixin[V]` | `LinearOperatorMixin[V, W]` (a zero map can be non-endomorphic — the `(L+C-S-F)` zero summand) |
| `IdentityOperator` (`:897`) | `LinearOperatorMixin[V]` | `LinearOperatorMixin[V]` (≡`[V,V]` — identity IS endomorphic; keep single-param) |
| `_AdjointOperator` (`:600`) | `LinearOperatorMixin[V]`, `__init__(inner: LinearOperator[V])` | `LinearOperatorMixin[W, V]`, `__init__(inner: LinearOperator[V, W])` (presents the SWAPPED `[W, V]`; already swaps domain/codomain at `:634-641`, metric-aware body `:643-670` unchanged) |

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
- **Carrier type** (`Din`/`Cout`, `Vector`-bound): `AngularFlux`, `AngularSourceSink`, `TimedFullField`,
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
ndarray — their `Din`/`Cout` are the ndarray/`Vector` level, NOT transport carriers (layering). The transport
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
  algebra (affine-map / linear-part-of-affine-map structure), the `[Din, Cout]` split, the carrier-vs-space two
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
  intrinsic-property tests: a non-endomorphic operator's `apply` carrier type is statically `Cout≠Din` (an
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
  `orpheus/transport/integral_kernel_operator.py` (the IntegralKernelOperator Protocol — refine to `[Din, Cout]`),
  `orpheus/numerics/iteration.py` (the SI/Krylov drivers — the variadic `*gains: LinearOperator[V, W]`),
  `pyproject.toml` (the §1 gating decision). Memos:
  `.claude/agent-memory/cross-domain-attacker/flux_to_sourcesink_operator_contract_frames.md`.
- **Durable principle banked:** an operator is the *linear part of an affine map* `Flux → SourceSink` (affine
  point domain, vector codomain, the cross-section unit-gain as fiber transition); the balance
  `(L+C−S−F)ψ=q` is a typed source/sink algebra; `[Din, Cout]` IS that type; `@overload` endomorphism stacks are
  the symptom of the missing second parameter.
