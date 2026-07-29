# SN operator assembly + problem posing — ground truth map

Grounding for a design review of `.claude/plans/sn_operator_realization_and_posing_plan.md`
(dated 2026-07-24) against branch `main` @ `8654d348` (post-DSA-merge).
Every claim is read from the code BODY, not a docstring (L33). Where a docstring
makes a claim the body contradicts, it is flagged **DOC-DRIFT**. Plan mismatches
are flagged **⚠ PLAN-DRIFT**.

Line numbers current at `8654d348`; re-derive via Nexus if the tree moves.

---

## 0. Verdict summary on the plan's factual claims

| # | Plan claim | Verdict |
|---|---|---|
| 1 | `self.L = build_streaming_collision(...)`, i.e. `self.L` holds `L+C` | **TRUE** (solver.py:1040) — and worse than stated: `L` has 3 meanings in one `__init__`, and `self.L`/`self.S`/`self.F` are production-DEAD |
| 2 | `InvertibleOperator` exists in `sn/operators/streaming.py` | **TRUE**, unchanged (streaming.py:445); note it is NOT in the module's `__all__` |
| 3 | `WithinGroupSystem.resolvent` holds `M` (not `M⁻¹`); field `gains` | **TRUE** (coupled_system.py:327/328, 483-488, 531-539) |
| 3b | "regular splitting" appears in code | **TRUE** — 5 code sites + 4 theory-page sites; the VOCABULARY exists, the TYPE does not |
| 4 | `A = L + C − S − B` composed in production | **TRUE**, one site: `A_AA = LC - S - B_a` (coupled_system.py:464) |
| 5 | eager fusion to `A⁻¹F` vs held `(A, F)` pair | **BOTH exist, split by path**: `KEigenvalue` holds `(A,S,F)` un-fused; `direct_eigenvalue` fuses eagerly; the SN **production** path holds NO pair at all — posing is implicit in method bodies |
| 6 | `SystemRole.A`/`.B` collide with `A`=loss / `B`=boundary | **TRUE**, and demonstrable inside a single docstring (operator.py:248). ⚠ minor PLAN-DRIFT: only `SystemRole.B` is the radial-characteristic carrier; `SystemRole.A` is the ordinary transport system |
| 7 | `CoupledOperator`/`CoupledField` = explicit block grid over a direct sum | **TRUE** (numerics/coupled_system.py:617, 232) |
| 8 | `SNMethodSpace` is a plausible "reaction realizer" input | **FALSE / ⚠ PLAN-DRIFT** — it is a per-FACE boundary-realizer payload; carries no σ, no scheme, no function space |
| 9 | DSA already carries "posing" vocabulary | **PARTIALLY** — it carries `corrector`/preconditioner/low-order-system vocabulary and a `str`-flag `acceleration="dsa"`; it introduces NO splitting/implicit-operator/schedule TYPE |

---

## 1. `SNSolver` operator attributes — the plan's `self.L = L + C` claim is TRUE

`orpheus/sn/solver.py:1017-1042` (`SNSolver.__init__`, the "Operator triple" block):

```python
1025  self.scattering_op = ScatteringOperator.from_solver_data(...)
1031  self.fission_op   = FissionOperator.from_solver_data(...)
1040  self.L = build_streaming_collision(sn_mesh, self.mat_xs)
1041  self.S = self.scattering_op
1042  self.F = self.fission_op
```

**`self.L` holds `L + C` — CONFIRMED.** `build_streaming_collision`
(`orpheus/sn/coupled_system.py:356-377`):

```python
374  return StreamingOperator(sn_mesh) + MultiplicationOperator(
375      coefficient=mat_xs.total_cross_section_field,
376      space=sn_mesh.full_field_space,
377  )
```

σ-free streaming `+` the collision diagonal `C = M[σ_t]`; declared return type
`"InvertibleOperator"` (coupled_system.py:358).

**DOC-DRIFT (the misnomer is written out).** `orpheus/sn/solver.py:1035` says
*"The full transport operator :math:`L = \Omega\cdot\nabla + \Sigma_t`"* — that
object is `L + C`, not `L`. Meanwhile `rebind_cross_sections`'s comment states
the honest convention one line above its own re-bind:
`orpheus/sn/solver.py:1100-1102` — *"the streaming leaf L is σ-free since #257
S8b — only C carries σ"*, then `self.L = build_streaming_collision(...)`.
The codebase knows `L` is σ-free while the attribute named `L` is not.

**Extra collision the plan does not mention: three meanings of `L` in one
`__init__`.** `orpheus/sn/solver.py:990` binds the local `L` to the Legendre
scattering order:

```python
990  L = min(scattering_order, min(len(m.SigS) - 1 for m in materials.values()))
994  self.scattering_order = L
```

So within `SNSolver.__init__`: `L` = Legendre truncation order (local, :990),
`self.L` = the streaming+collision composite (:1040), and the honest σ-free
streaming leaf `L` = `StreamingOperator` (the thing the composite's left summand
actually is).

**All three of `self.L` / `self.S` / `self.F` are PRODUCTION-DEAD.** Repo-wide
grep for reads: `self.S` / `self.F` are assigned at solver.py:1041-1042 and read
**nowhere** (the only `self.S`/`self.F` reads in the tree are `KEigenvalue`'s own
attributes, `orpheus/numerics/iteration.py:1172-1322`). `self.L` is read by
exactly two TESTS and no production code:

* `tests/sn/operators/test_typed_residual_evaluation.py:233` — `(solver.L, "L+C")`
  (the test LABELS it `"L+C"`, i.e. the test already carries the honest name)
* `tests/sn/verification/mms/test_curvilinear_operator_admits_mms.py:84` —
  `LC, S = solver.L, solver.scattering_op` (the test RENAMES it to `LC` at the
  call site)

The production solve paths do not touch these attributes: `_solve_source_iteration`
(solver.py:1568), `_solve_krylov` (solver.py:1756) and the `solve_sn` finalize
(solver.py:2140) each call `build_within_group_system(...)`, which builds its own
`LC` internally at `orpheus/sn/coupled_system.py:462`. The live cached operators
are `self.scattering_op` / `self.fission_op` (used at solver.py:1132, 1569, 1757,
1854, 2141, …), not the `S`/`F` aliases.

⇒ The `(L, S, F)` "operator triple" on `SNSolver` is a **vestigial 2026-era
framing surface**, not the assembly path. Retiring/renaming it costs 2 test edits.

---

## 2. `InvertibleOperator` — still there, unchanged, 4 direct constructors

**Exists under that name**: `orpheus/sn/operators/streaming.py:445`.

**What it is**: a CLASS, subclassing the composer with pinned leg types —

```python
445  class InvertibleOperator(
446      OperatorSum[
447          "FullField", "FullField", "StreamingOperator", "MultiplicationOperator",
448      ],
449  ):
```

i.e. it wraps exactly two legs: a `StreamingOperator` (`self.a`, aliased
`.streaming`, :620) and a `MultiplicationOperator` (`self.b`, aliased
`.diagonal`, :625). Constructor guards (streaming.py:564-603): both leg types
(TypeError), **mesh-identity** (`streaming.sn_mesh is diagonal.coefficient.mesh`,
:585), and **σ > 0 everywhere** (:593).

**Domain / codomain**: inherited from `OperatorSum` (operator.py:1389-1397), which
takes the left summand's — i.e. `StreamingOperator.domain` =
`sn_mesh.full_field_space` (streaming.py:305/310), the composite bulk ⊕ trace
space. Endomorphism on `FullField`.

**Capabilities advertised** (each verified in the body):

| Capability | Where | What it actually does |
|---|---|---|
| `is_invertible` → `True` | streaming.py:608-615 | the SOLE invertible `OperatorSum` (base is `self.a.is_invertible`, operator.py:1445) |
| `inverse()` → `SweepOperator` | streaming.py:745-768 | shadows the generic `GreenOperator` by MRO |
| `solve(rhs, *, initial_guess=None)` | streaming.py:770-832 | the WDD sweep; `initial_guess` **accepted-and-dropped** (`del initial_guess`, :831) |
| `solve_transpose(b)` | streaming.py:1010-1109 | reverse-scan `(L+C)⁻ᵀ`; 1-D-scan-family only |
| `apply` OVERRIDE | streaming.py:659-678 | `loss_representation.loss_action(self.sigma, psi)` — NOT the leaf sum |
| `apply_transpose` OVERRIDE | streaming.py:680-703 | `loss_action_transpose(self.sigma, phi)` |
| `is_adjointable` | inherited `a ∧ b` (operator.py:1416) | |
| `__sub__` dispatch | streaming.py:707-741 | `(L+C) - SNMaskedBoundaryOperator → ScheduledInvertibleOperator` |
| `block_role` | DERIVED `join(FULL, BULK) = FULL` | operator.py:1369, streaming.py:604-606 (hand-stamp retired) |

**Nexus/grep consumer census** (four-search discipline, AGENT.md P4):

*Direct constructors* `InvertibleOperator(...)` — 4 production-ish sites, all inside
streaming.py itself (`:436` the `L + C` dispatch; `:489`/`:491` docstring examples;
`:445` the class line), then **8 test sites** in
`tests/sn/operators/test_invertible_operator.py`, **9** in
`test_removal_form_matvec_sweep.py`, **2** in `test_typed_residual_evaluation.py`,
and **3 derivations diagnostics** (`derivations/diagnostics/diag_curvilinear_seed_sensitivity.py:102`,
`diag_sphere_fixedpoint_consistency.py:82,123`). Nexus `calls`-edge in-neighbours
confirm exactly those three diagnostics + the tests + `StreamingOperator.__add__`.

*Production TYPE/annotation consumers* (13 files): `sn/operators/sweep_operator.py`
(16 mentions, incl. two runtime `isinstance` gates at :154 and :185),
`sn/operators/scheduled_invertible.py` (14, incl. an `isinstance` ctor guard :110),
`sn/coupled_system.py` (9 — the `WithinGroupSystem.resolvent` union type :327 and
`build_streaming_collision`'s return :358), `sn/solver.py` (8 — all annotations
`_select_si_resolvent` :708, `_within_group_si`, plus docstrings),
`numerics/operator.py` (5, docstrings only), `numerics/green_operator.py` (2, docs),
`numerics/iteration.py` (1, docs), `transport/operators/multiplication_operator.py`
(1, docs), the two `__init__.py` docstrings.

*Test-file census*: 27 test files mention it; 32 mentions in
`tests/sn/operators/test_invertible_operator.py` alone.

*Doc-node census* (the silent-blast-radius search): **20 live theory pages**
mention `InvertibleOperator` — heaviest `docs/theory/methods/sn/loss_representation.rst`
(27), `docs/theory/foundations/operator_algebra.rst` (11), plus
`operator_inverse_family.rst` (4), `boundary_conditions.rst` (6),
`curvilinear_numerics.rst` (6), `history.rst` (7), `slab_one_group.rst` (7),
`acceleration.rst` (1), `adjoint.rst` (1), and 11 more. (Plus a stale
`docs/_build/html_tau_stepA/` tree — build artifacts, not sources.)

**Note**: `orpheus/sn/operators/streaming.py:114-116` — `__all__ = ["StreamingOperator"]`.
`InvertibleOperator` is importable but not exported. Any rename must still fix
every `from .streaming import InvertibleOperator` (coupled_system.py:152,
solver.py:90, sweep_operator.py:73/151/183, scheduled_invertible.py:71).

---

## 3. `WithinGroupSystem` — field names match the plan; `resolvent` does hold `M`

`orpheus/sn/coupled_system.py:282-328`, `@dataclass(frozen=True)`, **exactly four
fields, no methods**:

```python
325  loss: "CoupledOperator"
326  space: "CoupledSpace"
327  resolvent: "CoupledOperator | InvertibleOperator"
328  gains: "tuple[LinearOperator, ...]"
```

**The plan is CORRECT: `resolvent` is `M`, the un-inverted forward operator.**
Proof in the two constructor arms:

* **seedless** (coupled_system.py:477-488): `resolvent=LC`, `gains=(S, B_a)`,
  `loss=CoupledOperator([[A_AA]], ...)` — a 1×1 grid.
* **carrying** (coupled_system.py:502-539):
  `resolvent = CoupledOperator([[LC, A_AB], [None, march]], ...)` (:531),
  `gains=(N,)` with `N = CoupledOperator([[N_AA, None], [emission, B_b]], ...)` (:521),
  `loss = CoupledOperator([[A_AA, A_AB], [-emission, A_BB]], ...)` (:514).

And every consumer inverts it explicitly:
`SourceIteration(system.resolvent.inverse(), *system.gains, ...)`
(solver.py:850-853); `_maybe_window(base_resolvent.inverse(), S, sn_mesh)`
(solver.py:866); `_within_group_krylov(system.resolvent, *system.gains, ...)`
(solver.py:1788); `final_resolvent.solve(final_rhs_a)` (solver.py:2216).
So the field named `resolvent` is the operator whose **inverse** is the resolvent
in the usual sense — the naming inversion the plan flags is real.

`WithinGroupSystem` is a **pure record**: no methods, no `__post_init__`, no
consistency check that `loss == resolvent − Σ gains`. That invariant is upheld
only by the two constructor arms being written correctly side-by-side
(coupled_system.py:451-539).

**"regular splitting" IS in the code** (5 sites) and in the theory corpus (4+):

* `orpheus/sn/coupled_system.py:77-79` — *"they consume its **regular splitting**
  `A = M − N` (Hackbusch 2016 §11) — `M` the sweepable part inverted every step,
  `N` the lagged coupling gains"*; citation at :120-121.
* `orpheus/sn/coupled_system.py:290` — the record docstring.
* `orpheus/sn/solver.py:172` (module comment), `:724` (`_select_si_resolvent`:
  *"the regular splitting `(L+C−B) = M − B_upper`"*).
* `orpheus/sn/loss_representation/sweep_schedule.py:173`.
* Docs: `docs/theory/foundations/coupled_block_operator.rst:413,444`,
  `docs/theory/foundations/boundary_conditions.rst:731`,
  `docs/theory/methods/sn/history.rst:218`,
  `docs/theory/methods/sn/cartesian_multid.rst:67,3809`.

⇒ The **vocabulary** of a named splitting is fully established in prose. What does
not exist is a splitting **object with behaviour** — the record cannot apply `A`,
cannot apply `M`, cannot check `A = M − N`; every consumer destructures it.

---

## 4. The operator algebra as actually written — one composition site

`orpheus/sn/coupled_system.py:461-467`:

```python
462  LC = build_streaming_collision(sn_mesh, mat_xs)     # L + C
463  B_a = SNBoundaryOperator(sn_mesh)
464  A_AA = LC - S - B_a
466  # C-fwd explicit stamp: ...
467  A_AA.system_role = SystemRole.A
```

**The honest algebra in production is `A_AA = (L + C) − S − B_a`.** Fission is
structurally absent (it enters as `q_ext` — stated coupled_system.py:421-423),
so the equation as posed is `A = L + C − S − B`, matching the AGENT.md durable
shape and `docs/theory/foundations/operator_algebra.rst`'s Option-Y spelling.

**The resulting TYPE is a left-nested `OperatorSum`, not a flat n-ary sum.**
Trace of the dispatch:

1. `LC - S`: `LC` is an `InvertibleOperator`; its `__sub__` (streaming.py:715-741)
   checks `isinstance(other, SNMaskedBoundaryOperator)` — `S` is a
   `ScatteringOperator`, so it falls to `super().__sub__(other)`, i.e.
   `LinearOperator.__sub__` (operator.py:695-698) =
   `OperatorSum(self, ScaledOperator(-1.0, other))`.
2. `(...) - B_a`: the left is now a plain `OperatorSum`, which does not override
   `__sub__`, so again `OperatorSum(OperatorSum(LC, −S), ScaledOperator(-1.0, B_a))`.

So `A_AA` = `OperatorSum(OperatorSum(InvertibleOperator, ScaledOperator(−1, S)),
ScaledOperator(−1, B_a))`. `A_AA.apply` therefore recurses through two
`self.a.apply(x) + self.b.apply(x)` levels (operator.py:1399-1400) and
`A_AA.is_invertible` is `self.a.is_invertible` chained down the LEFT SPINE to
`InvertibleOperator.is_invertible == True` (operator.py:1421-1445) — the
"canonical-ordering contract" documented there. **The invertibility of the whole
loss is therefore an artefact of spelling order**: `(-S) + LC - B_a` would report
non-invertible.

Then `A_AA` becomes a BLOCK of a `CoupledOperator` grid, never the top-level `A`:

* seedless → `loss=CoupledOperator([[A_AA]], domain=space, codomain=space)` (:484)
* carrying → `loss=CoupledOperator([[A_AA, A_AB], [-emission, A_BB]], ...)` (:514)
  with `A_BB = march - B_b` (:501)

**Sign asymmetry inside the grid** (documented coupled_system.py:24-51, verified
in the body): `(A,B) = +A_AB` POSITIVE (the `RadialCharacteristicSeeding`
internalizes its loss sign), `(B,A) = -emission` (a `ScaledOperator`),
`(B,B) = march − B_b`.

**Post-construction MUTATION of role**: `A_AA.system_role = SystemRole.A` (:467)
and `N_AA.system_role = SystemRole.A` (:520). The role is assigned by the
composition CONTEXT after the object is built, because `_join_system_roles`
(operator.py:317-335) propagates `None` from the model-generic leaves. So the
composite's system membership is an out-of-band `setattr`, not a constructor
argument and not derivable.

Also note `evaluate_residual` (solver.py:231-332) re-poses the same algebra
INDEPENDENTLY in its docstring (`"On a SEEDLESS mesh ``A = L + C - S - B``"`,
:247) and reaches the operator through `_bare_loss_arm(system)` =
`system.loss.blocks[0][0]` (solver.py:446-455) — i.e. it un-wraps the 1×1 grid to
get back the `A_AA` composition, because `evaluate_residual`'s bare arm takes the
composition, not the arity-guarded grid.

---

## 5. The operator-algebra surface, and where the eigenproblem is POSED

### 5a. `orpheus/numerics/operator.py` — the algebra surface

`LinearOperator` is a `Protocol[Domain, Codomain]` (:504) that **carries real
default-method bodies for the algebra dunders** (:675-807), so an explicit
subclass inherits both the `apply` contract and the algebra:

| Dunder / member | Line | Returns |
|---|---|---|
| `__add__` / `__radd__` | 685 / 690 | `OperatorSum(self, other)` |
| `__sub__` / `__rsub__` | 695 / 700 | `OperatorSum(self, ScaledOperator(-1.0, other))` |
| `__mul__` / `__rmul__` (scalar) | 705 / 710 | `ScaledOperator(float(other), self)` |
| `__neg__` | 713 | `ScaledOperator(-1.0, self)` |
| `__truediv__` (scalar) | 724 | `ScaledOperator(1/α, self)` |
| `__matmul__` | 739 | `OperatorProduct(self, other)` |
| `__and__` / `__rand__` | 748 / 763 | `TensorProductOperator` (⊗) |
| `__call__` | 766 | alias for `apply` |
| `__pow__` | 775 | `IdentityOperator` (n=0) / repeated `@`; **negative n raises** |
| `adjoint()` / `.H` | 813 / 856 | `_AdjointOperator(self)`, eager `MissingAdjoint` gate |
| `as_matrix()` | 868 | probing dense `Op → Mat` functor |

There is **no `solve` and no `inverse` on the Protocol base.** They are per-class
"realization verbs", declared only where a native realization exists — with three
runtime predicates as the advertisement: `is_invertible` (:605, default `False`),
`is_adjointable` (:623, default `False`), `is_assemblable` (:637, default `False`).
Static bridges `invertible()` / `adjointable()` / `assemblable()` (:1089/1113/1128)
narrow to `SupportsInverse` / `SupportsAdjoint` / `SupportsAssembly` Protocols.

Composers:

* **`OperatorSum`** (:1287) — generic over `[Domain, Codomain, SummandA, SummandB]`.
  `apply` = `a.apply(x) + b.apply(x)` (:1399). `apply_transpose` = the `(A+B)ᵀ`
  law with a `MissingAdjoint` guard (:1402). `is_invertible` = **LEFT-SPINE-HEAD
  only** (:1445). `inverse()` → `GreenOperator(self)`, the preconditioned-splitting
  iterative inverse (:1447-1465). Explicit comment at :1467-1470: *"NO `solve` on
  a generic sum … solving is `.inverse().apply(b)`"*. `assemble()` = the additive
  homomorphism `[A+B] = [A]+[B]` (:1478).
* **`OperatorProduct`** (:1504) — `apply`, `solve` (:1591, `B⁻¹(A⁻¹q)`),
  `inverse()` → `InverseOperator` (:1633), `assemble` (:1676).
* **`ScaledOperator`** (:1700) — `apply` (:1758, forwards `*extra/**kwextra`),
  `apply_transpose` (:1766), `inverse()` → `ScaledOperator` (:1782),
  `assemble` (:1812).
* **`_AdjointOperator`** (:1155) — private, constructed ONLY via `.adjoint()`/`.H`.
  Swaps `domain`↔`codomain` (:1195-1202). `apply` = the **G-adjoint**
  `G_V⁺ ⊙ apply_transpose(G_W ⊙ y)` delegated to `FunctionSpace.apply_metric` /
  `apply_inverse_metric` (:1204-1227). Carries the **swap law**:
  `is_invertible = invertible(inner) and adjointable(inner.inverse())` (:1249),
  `inverse()` returns `inner.inverse().H` as an OBJECT IDENTITY (:1251-1284).
  Its own `apply_transpose` raises `NotImplementedError` (:1229-1234).
* **`InverseOperator`** (:2103) — the GENERIC inverse-family member;
  `apply` delegates to `inner.solve(x)` bit-identically (:2153-2165), accepts and
  drops `initial_guess`. Named siblings with stronger invariants: `SweepOperator`
  (triangular sweep), `GreenOperator` (splitting), `MatrixInverseOperator` (LU).
  Back-half shared via `InverseWrapMixin` (:2007).

### 5b. Is there a pencil-like held `(A, F)` pair?

**Three different answers on three different paths — this is the load-bearing
finding for the plan's posing question.**

**(i) The SN PRODUCTION eigenvalue path holds NO pair at all.**
`solve_sn` (solver.py:1990) → `power_iteration(solver, max_iter=max_outer)`
(solver.py:2095). `power_iteration` (`orpheus/numerics/eigenvalue.py:203-272`)
never sees an operator: it consumes the `EigenvalueSolver` Protocol
(eigenvalue.py:82-166) = `initial_flux_distribution` / `compute_fission_source` /
`solve_fixed_source` / `compute_keff` / `converged`, plus the optional
`ProductionRateSolver.compute_production_rate` (:170-200). The loop body
(:238-266) is four method calls and a renormalising division. The eigenproblem's
posing — `A_loss = L+C−S−B`, `M = F`, `k = μ` — exists **only as prose** in the
module docstring (eigenvalue.py:22-36) and as the arithmetic scattered across
`SNSolver.compute_fission_source` (`fission_op.apply(φ)/keff`, solver.py:1132),
`SNSolver.solve_fixed_source` (a `str`-dispatch on `self.inner_solver`,
solver.py:1141-1146), and `SNSolver.compute_keff` (an assembled
production/net-removal RATIO of `IntegratedReactionRate`s + a boundary-leakage
term, solver.py:1250-1305). **No object anywhere in the SN production path
represents "the posed k-eigenproblem".**

Note `SNSolver.compute_keff` does NOT use `A`/`F` at all — it re-derives k from
reaction-rate functionals plus `_boundary_leakage_rate` (solver.py:1307-1401),
which reads the stashed `self._psi_typed` trace and rescales by a fission-production
ratio. That estimator's consistency with the posed problem is maintained by
comment and by a dedicated gate (`tests/sn/eigenvalue/test_keff_estimator_gate.py`,
cited solver.py:1427), not by construction.

**(ii) `KEigenvalue` DOES hold the triple un-fused — but the SN production path
does not use it for the forward solve.**
`orpheus/numerics/iteration.py:1049`, `__init__(A, S, F, *, …)` (:1138-1204):

```python
1164  if not A.is_invertible:  raise NotInvertible(...)
1171  self.A = A
1172  self.S = S
1173  self.F = F
1192  self._inner = SourceIteration(seeded_inverse(self.A), self.S, ...)
```

The triple survives as three attributes for the object's lifetime; `A` is never
fused with `F`. This IS the pencil-like object the plan is asking about, and the
posing is spelled at iteration.py:1050 (*"the k-eigenvalue problem
`(A − S)ψ = Fψ/k`, posed from an operator triple"*) and in the boundary comment
block iteration.py:1206-1220. Its `compute_keff` (:1286-1324) is the **operator-form
Rayleigh estimator** `Σ(Fψ) / (Σ(Aψ) − Σ(Sψ))` — a genuinely different (and
leakage-inclusive-through-`A`) estimator from `SNSolver.compute_keff`.
`solve()` (:1349) delegates to `power_iteration(self, ...)` (:1400).

Production consumer: **exactly one** — `orpheus/sn/solver.py:2472-2474`, inside
the **adjoint** entry family (`_adjoint_posing_parts` at solver.py:2322 →
`KEigenvalue(...)`). Nothing else in `orpheus/` constructs it (the other
`KEigenvalue` hits are docstrings + the `numerics/__init__.py` re-export).
6 test files use it. ⇒ **the un-fused-triple posing object exists and is
EXERCISED, but only on the adjoint path**; the forward SN eigenvalue solve
bypasses it entirely.

**(iii) `direct_eigenvalue` eagerly fuses — on dense ndarrays, not operators.**
`orpheus/numerics/eigenvalue.py:363-468`. Signature `(A: np.ndarray, F: np.ndarray)`
— **bare matrices, not `LinearOperator`s**. The last line (:468) is the eager fuse:

```python
468  return dominant_eigenpair(np.linalg.solve(A, F), imag_tol=imag_tol)
```

`dominant_eigenpair(M)` (:288-360) then takes the already-materialized resolvent
`M = A⁻¹F` and does `np.linalg.eig` + the Perron–Frobenius validation
(complex-dominant rejection :351, sign normalisation via `_sign_normalised` :275).
`rayleigh_quotient_iteration(A, F)` (:471-610) is the one algorithm that keeps the
pair un-fused throughout (the bordered form, `M = A⁻¹F` never formed — :575-593),
again on ndarrays.

⇒ The `(A, F)` pair exists as a **calling convention on three functions** and as
attributes on ONE class used only for adjoint. There is no `Pencil` type, no
`PosedEigenproblem`, no shared object the three engines agree on. `power_iteration`
and `direct_eigenvalue` "solve the SAME posed problem" only by the assertion in a
docstring (eigenvalue.py:392-394).

---

## 6. `BlockRole` / `SystemRole` — the collision is real and demonstrable

**`BlockRole`** (`orpheus/numerics/operator.py:208-231`), three members:

```python
229  BULK = "bulk"        # A_bb only — C, S, F
230  FULL = "full"        # off-diagonal coupling — L (the only irreducible one)
231  BOUNDARY = "boundary"  # A_ss only — a realized boundary law B
```

Join law `_join_block_roles` (:283-314): union lattice, `None` propagates.
`isinstance` markers `BulkOperator`/`FullOperator`/`BoundaryOperator` (:362/368/374)
via `_BlockRoleMeta.__instancecheck__` reading `op.block_role` (:358).

**`SystemRole`** (`orpheus/numerics/operator.py:234-280`), three members:

```python
278  A = "system_a"
279  B = "system_b"
280  COUPLED = "coupled"
```

Docstrings (verbatim, :246-265): *"**System A** — the transport bulk ⊕ trace (the
angular-flux `FullField` …), governed by `A_AA = L + C − S − B`"*; *"**System B**
— the ψ½ radial-characteristic ray …, governed by `A_BB`"*; *"`COUPLED` — maps
BETWEEN the systems"*.

**⚠ PLAN-DRIFT (minor, but worth correcting in the review):** the plan says
`SystemRole.A` / `SystemRole.B` "name two radial-characteristic carrier systems."
Only **B** is the radial-characteristic (ψ½) carrier. **A** is the ordinary
transport composite — i.e. every seedless slab/2-D solve is `SystemRole.A`, with
no radial-characteristic content whatsoever.

**The collision the plan flags is REAL and visible inside a single docstring.**
`orpheus/numerics/operator.py:248` reads *"governed by `A_AA = L + C − S − B`"* —
in that one line, `A` (as `A_AA`) means System A **and** the loss operator, and
`B` means System B's label **and** the boundary operator. Both letters carry two
meanings in the same sentence. Downstream, `orpheus/sn/coupled_system.py` uses
`A_AA`/`A_AB`/`A_BA`/`A_BB` for grid POSITIONS while `B_a`/`B_b` are BOUNDARY
operators (:463, :497) and `A_BB = march - B_b` (:501) — so `A_BB` contains `B_b`.

**Where each member is used** (complete census):

| Member | Assignment sites |
|---|---|
| `SystemRole.A` | `orpheus/sn/coupled_system.py:467` (`A_AA.system_role = …`), `:520` (`N_AA.system_role = …`) — both post-construction mutations |
| `SystemRole.B` | `orpheus/sn/operators/boundary.py:537` (`RadialCharacteristicBoundaryOperator`, class attr), `orpheus/sn/operators/radial_characteristic.py:270` (`RadialCharacteristicOperator`, class attr) |
| `SystemRole.COUPLED` | `orpheus/numerics/coupled_system.py:663` (`CoupledOperator`, class attr), `orpheus/sn/operators/radial_characteristic.py:684`, `:953`, `:1194` (the three off-diagonal/coupling operators) |

**`system_role` is WRITE-ONLY in the current tree.** Grep across `orpheus/` finds
only the assignments above plus the propagation machinery: the union join
`_join_system_roles` (operator.py:317-335, called from `OperatorSum.__init__`
:1375) and the passthroughs on `_AdjointOperator` (:1193) and `ScaledOperator`.
**No consumer branches on it, guards on it, or dispatches on it** — there is no
`if op.system_role ==` / `SystemRole.` read anywhere outside the enum's own
module and the assignment sites. (`orpheus/numerics/coupled_system.py:30` mentions
it in prose: *"`system_role` tags are …"*.) The `_join_block_roles` twin at
operator.py:304-310 documents the collapse trigger: *"a THIRD parallel role axis
(a DSA / multiphysics role) makes the shared abstraction pay — unify then."*

---

## 7. `CoupledOperator` / `CoupledField` — yes, an explicit block grid over a direct sum

`orpheus/numerics/coupled_system.py`.

**`CoupledField`** (:232, a dataclass) — the direct-sum STATE: `systems: tuple[...]`
of `SystemField` members (Protocol at :166, requiring `to_flat` :190 and `copy` :192).
API: `n_systems` (:280), `_check_partner` (:286, the cross-arity/class guard),
`_map_binary`/`_map_unary` (:306/:320), full vector algebra `__add__` `__sub__`
`__neg__` `__mul__` `__rmul__` `__truediv__` (:323-338), `copy` (:341),
**`to_flat`** (:347, concatenates the members — this is what sizes GMRES `restart`)
and **`from_flat`** (:361, classmethod-style rebuild via `_member_from_flat` :195).

**`CoupledSpace`** (:398, a `FunctionSpace["CoupledField"]`) — `systems` tuple
(:445), `zeros_factory` (:448), `from_systems` (:469), `zeros()` (:496),
`n_systems` (:535), **`system_slices`** (:540, the flat offsets),
`apply_metric` / `apply_inverse_metric` (:573/:589, member-wise — this is what
makes `.H` metric-correct across the direct sum), `inner_product` (:599).

**`CoupledOperator`** (:617, `LinearOperator["CoupledField","CoupledField"]`) —
the N×N typed block grid `A_ij : System_j → System_i`. Class attribute
`system_role = SystemRole.COUPLED` (:663).

Constructor (`blocks`, `*, domain: CoupledSpace, codomain: CoupledSpace`, :665-736)
validates, in order: row count vs codomain arity (:676), per-row column count
(:682), **no all-`None` row** (:689) and **no all-`None` column** (:696), then
**per-block space typing** — every block that DECLARES spaces must match
`domain_members[j]` / `codomain_members[i]` or
`IncompatibleOperatorComposition` raises at construction (:707-733). A `None`
block IS the structural zero map (no zero arithmetic runs, :627-628, :790-798).

API surface:

| Member | Line | Behaviour |
|---|---|---|
| `blocks` / `n_rows` / `n_cols` | 741 / 746 / 750 | row-major grid accessors |
| `domain` / `codomain` | 754 / 758 | the `CoupledSpace`s |
| `apply` | 780 | `y_i = Σ_j A_ij x_j`, arity-guarded |
| `apply_transpose` | 800 | the transposed grid `(Aᵀ)_ji = (A_ij)ᵀ` (Euclidean only) |
| `is_adjointable` | 832 | all present blocks adjointable |
| `_diagonal_blocks` / `_triangular_orientation` | 840 / 850 | structure detection |
| `_substitution_ready` / `_transpose_substitution_ready` | 877 / 898 | triangular-solve readiness |
| `_materialization_route` | 922 | the LU-EXTRACT fallback |
| `is_invertible` | 941 | **structure-keyed, two DIRECT routes**: triangular substitution, else LU materialize |
| `solve` / `solve_transpose` | 959 / 982 | block back/forward substitution |
| `_solve_triangular` | 1007 | the substitution body |
| `inverse` | 1077 | `CoupledSubstitutionOperator` or `MatrixInverseOperator` |
| `as_matrix` | 1106 | dense materialization via the space's zero exemplar |
| `is_assemblable` / `assemble` | 1154 / 1158 | recursive block assembly |

`CoupledSubstitutionOperator` (:1212) is the inverse-family member for the
triangular case: `apply` (:1255) runs the substitution, `apply_transpose` (:1277).

**Load-bearing structural note for the review:** this grid is the ONLY place in
the SN stack where a multi-block operator equation is a first-class typed object
with construction-time position checking. The `A_AA` composition (§4) has none of
that — it is an untyped left-nested `OperatorSum` whose invertibility depends on
spelling order.

---

## 8. `SNMethodSpace` — exists, but is a per-FACE boundary payload (⚠ PLAN-DRIFT)

**Exists under that name**: `orpheus/sn/mesh/method_space.py:72`,
`@dataclass(frozen=True)`, `__all__ = ["SNMethodSpace"]` (:68).

**What it holds** — five fields, four optional (:92-96):

```python
92  quadrature: "Quadrature"
93  face: Optional[str] = None
94  inflow_indices: Optional[np.ndarray] = None
95  mesh: "Optional[Mesh1D | Mesh2D]" = None
96  trace: "Optional[AngularTraceSpace]" = None
```

Factories: `minimal(quadrature)` (:99, everything else `None`) and
`for_face(*, mesh=None, quadrature, face, trace=None)` (:112-158, which derives
`inflow_indices = trace.inflow_indices_for_face(face)`). One method:
`inflow_indices_for_face` (:160), delegating to the held `trace`.

**Its ONE job** (module docstring :3-5, verified by the consumer census): *"the
SN-specific argument to `SNBoundaryRealizer.realize`"* — the payload that turns a
`BoundaryTraceLaw` into a 1-arg `LinearOperator`. Consumers:
`orpheus/sn/boundary/realizer.py:150` (the `realize(law, method_space)` parameter),
`orpheus/sn/mesh/augmented_mesh.py:420` (`SNMethodSpace.for_face(...)` — the
canonical construction site, per face, inside `realize_boundary_law`), and six
`orpheus/geometry/boundary/*.py` law modules that build a `minimal(quad)` for
their own realization examples/self-tests (`vacuum.py:39`, `white.py:61`,
`reflective.py:67`, `prescribed_inflow.py:64`, `_base.py:119,385`,
`_realizer.py:246`). Mirror type: `orpheus/diffusion/method_space.py`
(`DiffusionMethodSpace`, explicitly *"the mirror of `SNMethodSpace`"* :5).

**⚠ PLAN-DRIFT — it is NOT a plausible "reaction realizer" input**, on four
counts, all checkable:

1. **No cross-section data.** No `mat_xs` / `MaterialXSField` / σ field. A
   reaction operator (`C = M[σ_t]`, `S`, `F`) needs exactly that — compare
   `build_streaming_collision(sn_mesh, mat_xs)` (coupled_system.py:356) and
   `ScatteringOperator.from_solver_data(mat_xs=…, quadrature=…, scattering_order=…,
   full_field_space=…)` (solver.py:1025-1030).
2. **No function space.** No `full_field_space` / `FunctionSpace`. Every reaction
   leaf must advertise `domain`/`codomain` for the `OperatorSum` guard to validate
   the build (streaming.py:291-296 spells this requirement out).
3. **No scheme.** No `DiscretizationSchemeBase` — so nothing that depends on
   `spatial_basis_per_axis` / `is_affine_scannable` can be built from it.
4. **It is FACE-SCOPED, not domain-scoped.** `face: Optional[str]` and
   `inflow_indices` are codimension-1 data; the canonical factory is literally
   named `for_face`. A reaction operator is a bulk (codimension-0) object.
   Additionally its `mesh` field is documented as *"OPTIONAL metadata (C5.3, #225):
   nothing in the realizer chain reads it"* (:134-136) — i.e. the one field a
   reaction realizer would most need is the one declared unread.

Any "reaction realizer" would take `(sn_mesh, mat_xs)` — the pair
`build_within_group_system` already takes (coupled_system.py:380-385) — not this
type. Reusing the NAME "method space" for both would collide two very different
scopes.

---

## 9. DSA / acceleration — what landed 2026-07-27 and how it plugs in

`orpheus/sn/acceleration/` — two files: `__init__.py` (14 lines, re-exporting
`DSACorrection`, `DSALowOrderSystem`) and `dsa.py` (674 lines). The package's
existence on the SN side is the R4 ruling (`__init__.py:6-9`): *"an accelerator's
low-order coefficients are properties of the SN discretization … not of standalone
diffusion physics."*

**The two types:**

* **`DSALowOrderSystem`** (`dsa.py:137`, `@dataclass(frozen=True)`) — the
  per-group consistent low-order EDGE systems. Fields (:154-164): `a_low`
  `(ng, K+1, K+1)`, `g_map` `(ng, K+1, 2K)` (the displacement→source map), `_lu`
  (per-group `scipy.lu_factor` tuple), `_dh` `(ng,K)` = `D_i/h_i`, `_a_coef`
  `(ng,K)` = `σ_s1/(σ_t−σ_s1)`. Builders: `from_sn_mesh` (:169) —
  **loud admission gates**: 1-D Cartesian + `Mesh1D` (:191, else
  `NotImplementedError` citing #282 for curvilinear and Alcouffe corner moments
  for 2-D), `scheme.key == "diamond_difference"` (:204, the R5 seam for
  weighted-diamond), and `bc.kind ∈ {"vacuum","reflective"}` (:213-223) — then
  `_build` (:257) which assembles Larsen (23a–f)/(27)/(38)/(39) rows with TWO
  numerically-guarded convention boundaries (`Σw == 2` :279, `W2 == 1/3` :286)
  and a positivity guard `σ_t − σ_s1 > 0` (:297). Methods: `solve_correction`
  (:372, LU-solves for edge `f_0`), `cell_update` (:417, Larsen (28a)),
  `moment1_update` (:427, Larsen (28b)).
* **`DSACorrection`** (`dsa.py:453`, a `LinearOperator["FullField","FullField"]`)
  — the correction operator. `apply(displacement)` (:578-674) is
  `Δψ ↦ P [(28)] A_low⁻¹ G R Δψ`: restriction `R` = `interior.integrate_angular()`
  (:623) plus, on the P1 arm, the frame's ℓ=1 row `w·μ` (:634); then the LU solve,
  the (28a)/(28b) updates, and injection `P` — isotropic `/Σw` on the P0 arm
  (:645) or Larsen's (33) `Ψ = f₀ + (μ/W₂)f₁` on the P1 arm (:638-641). It
  corrects **both blocks**: bulk `AngularDisplacement` and a **trace**
  `AngularBoundaryDisplacement` from the wall-edge solutions (:660-670,
  load-bearing under the lagged reflective `B` gain — docstring :504-520 records
  `ρ > 1` observed without it). Input type gate at :611: `AngularDisplacement`
  (SI posture) or `AngularFlux` (Krylov swept vector); a `HarmonicMomentFlux`
  refuses loudly. Output is always the DISPLACEMENT role so the torsor add
  `ψ ⊕ Δψ` stays well-formed (:578-588). `from_sn_mesh` (:555) is the factory.

**How it plugs into the iteration — via a `corrector` PARAMETER, not a type:**

* `SourceIteration.__init__(A_inv, *gains, max_iter, tol, corrector=None)`
  (`orpheus/numerics/iteration.py:589-596`). The loop applies it at
  `iteration.py:730-731`:
  ```python
  730  if self.corrector is not None:
  731      psi = psi + self.corrector.apply(psi - psi_prev)
  ```
  Byte-inert when `None` (:728-729). Contract documented :535-564: correctness-safe
  because the correction vanishes with the displacement.
* Krylov posture: `_within_group_krylov(..., corrector=None)`
  (`orpheus/sn/solver.py:458-518`). With a corrector it builds the
  **transport-corrected left preconditioner** at :507-511:
  ```python
  507  sweep = seeded_inverse(LC)
  509  def preconditioner(q):
  510      swept = sweep.apply(q)
  511      return swept + corrector.apply(swept)
  ```
  Without one, an explicit identity `lambda q: q` (:505) — issue #200 still tracks
  the block-inverse face preconditioner.
* Threading: `_within_group_si(..., corrector=None)` (solver.py:781) forwards it
  to `SourceIteration` (:874) with **two loud refusals**: the coupled
  (curvilinear/carrying) arm raises `NotImplementedError` (:844-849, "#282 — no
  stability theory"), and the 2-D moment-windowed arm raises (:867-872).
* **Entry point**: `solve_sn_fixed_source(..., acceleration: str | None = None)`
  (solver.py:2928). Validated against `(None, "dsa")` at :3050-3054, and
  `corrector = DSACorrection.from_sn_mesh(sn_mesh, scattering_order=…)` at
  :3056-3061, then passed to the SI path (:3080) or the Krylov path (:3090).
  Also threaded through `_solve_fixed_source_si` (:3102 → :3177) and
  `_solve_fixed_source_krylov` (:3312 → :3414).
* **DSA is FIXED-SOURCE ONLY.** `solve_sn` (the eigenvalue entry, solver.py:1990-2003)
  has **no `acceleration` parameter**, and neither `_solve_source_iteration`
  (solver.py:1571-1575) nor `_solve_krylov` (solver.py:1788-1792) passes a
  `corrector` — both call the builders with the corrector defaulted to `None`.

**Does DSA carry the plan's "posing" vocabulary?** Partially, and instructively:

* ✅ **"low-order system"** is a first-class named TYPE (`DSALowOrderSystem`) that
  holds an assembled operator + its factorization + a source MAP (`g_map`) —
  i.e. it is closer to a "posed system" object than `WithinGroupSystem` is
  (which holds no factorization and no map, and has no methods).
* ✅ **"preconditioner" / "corrector"** vocabulary is live, but as a **parameter
  name and a closure**, not a type: the Krylov preconditioner is a bare `def`
  (`solver.py:509`) and the SI corrector is a `LinearOperator | None` slot.
* ❌ **No splitting type.** The `A = M − N` splitting is still the un-typed
  `WithinGroupSystem` record's two fields; DSA did not mint one. `_select_si_resolvent`
  (solver.py:707-775) still returns a bare `tuple[resolvent, gains]`.
* ❌ **No "implicit operator" type.** The DSA-preconditioned matvec `M⁻¹ ≈ (A − Σg)⁻¹`
  is realized as a Python closure captured by `KrylovAcceleration`, never as an
  object.
* ❌ **No "schedule" type from DSA.** A `SweepSchedule` type does exist, but
  pre-DSA (`orpheus/sn/loss_representation/sweep_schedule.py`, consumed by
  `_select_si_resolvent` at solver.py:758/773 and by `ScheduledInvertibleOperator`).
  DSA's own selector is the stringly-typed `acceleration: str | None`, which joins
  the existing `inner_solver: str` (solver.py:928) and `inner_schedule: str`
  (:937) family of string flags — i.e. **the DSA merge ADDED a third string-flag
  posing knob** rather than converting any of them to a type.

⇒ The DSA work is a strong data point *for* the plan's thesis (it shipped a real
posed low-order system as a type) and simultaneously *against* the premise that
the posing vocabulary is already in place (it reached for a `str` flag and a
closure at the two places a type would have been the elegant move).

---

## 10. Cross-cutting observations the review may want

1. **`build_within_group_system` is called fresh on every inner solve** —
   solver.py:1568 (SI, once per outer iteration), :1756 (Krylov, once per outer
   iteration), :2140 (finalize). Each call re-runs `build_streaming_collision`
   (coupled_system.py:462), re-constructs `SNBoundaryOperator` (:463), and on a
   carrying mesh re-constructs `Seeding`/`Emission`/`B_b`/`march` (:495-500) and
   two `CoupledSpace`s with fresh `zeros` lambdas. The record is frozen, but its
   construction is not memoized anywhere. (`self.L` — the one cached composite —
   is the dead attribute of §1.)
2. **Two independent "loss" spellings coexist by design and are kept in sync by
   comment**: the grid `system.loss` (used by the certificate on the coupled arm,
   solver.py:1630/1814) and `_bare_loss_arm(system)` = `system.loss.blocks[0][0]`
   (the un-wrapped `A_AA` composition, used on the seedless arm). The gate that
   they agree is the `test_typed_residual_evaluation` family.
3. **`build_coupled_system`** (coupled_system.py:163-219) is now a pure VIEW —
   `system = build_within_group_system(...); return (system.loss, system.space)`.
   Its docstring names its planned consumers as *"the d2 `evaluate_residual`
   re-type, the assembly arm, and the DSA substrate (#2)"* (:174-177) — but the
   landed DSA does NOT consume it (§9: `DSALowOrderSystem.from_sn_mesh` reads
   `sn_mesh` + `mat_xs` directly, dsa.py:225-254). **DOC-DRIFT**: that
   "PLANNED such consumers" list is now stale for the DSA item.
4. **`_reflect_outflow_into_inflow`** (solver.py:1909-1983) survives with
   **exactly ONE production call site**: solver.py:2162, the eigenvalue finalize
   sweep. Its own docstring is **DOC-DRIFT on both of its other claimed
   consumers** (an L33 catch — two docstrings in the same file contradict each
   other):
   * It claims (solver.py:1920-1924) that *"the DIRECT fixed-source SI loop
     (`_solve_fixed_source_si`) … call this helper"*. But
     `_solve_fixed_source_si`'s own docstring says the opposite —
     solver.py:3129-3131: *"The whole-trace `_reflect_outflow_into_inflow` route
     is **NOT needed** on this path"* — and there is no call in its body.
   * It claims (solver.py:1947-1950) the helper survives for *"Phase 3's
     octant-group Gauss-Seidel resolvent (which calls it face-restricted between
     octant-group sweeps — see `faces`)"*. The G-S resolvent actually passes
     `reflect=self.lower.reflect_rows_inplace`
     (`orpheus/sn/operators/scheduled_invertible.py:260`), a
     `SNMaskedBoundaryOperator` method (`boundary.py:770`) — **not** this helper.
   * Consequence: the `faces: Iterable[str] | None = None` parameter
     (solver.py:1910) has **zero callers passing it** — production or test. It is
     a dead parameter guarded by a stale docstring.
   AGENT.md's durable-shape section says the *"`_reflect_outflow_into_inflow`
   driver shim"* is RETIRED — true of the DRIVER route; the helper survives for
   the ONE driver-less finalize sweep. Worth a correction in AGENT.md and a
   docstring fix here.
5. **The `str`-flag posing family**, all on one entry point:
   `inner_solver ∈ {"source_iteration","krylov"}` (solver.py:928),
   `inner_schedule ∈ {"jacobi","gauss_seidel"}` (:937 and again at :747 in
   `_select_si_resolvent` — the SAME validation spelled twice),
   `acceleration ∈ {None,"dsa"}` (:3050), `boundary_condition: str` (:2920),
   `eigenvalue_method: str` (iteration.py:1149). Nexus `discriminations` territory.
