---
name: issue-256-apply-endomorphism-step3
description: #256 Step 3 — retype LinearOperator.apply to the V→V endomorphism via class-level Generic[V] bound to Vector; PASS-WITH-NITS, keep the TYPE_CHECKING apply-stub (self-typing is strictly worse), one real NIT = incomplete Generic[V] retype inside iteration.py drivers
metadata:
  type: project
---

# #256 Step 3 — `apply(x: V) -> V` endomorphism retype review

PASS-WITH-NITS (`feature/field-typed-operator-algebra`, uncommitted, on top of step-1 commit
`cfb651b`). Retypes `LinearOperator.apply` from `(x: np.ndarray)->np.ndarray` to the honest
endomorphism `(x: V, /)->V` via CLASS-LEVEL `Generic[V]` bound to the `Vector` Protocol:
`LinearOperator(Protocol[V])`, `LinearOperatorMixin(Generic[V])`, all composers `[V]`; SN leaves
bind `LinearOperatorMixin["TimedFullField"]`. Removes 10 `# type: ignore[arg-type]` from the
Mixin dunders. Annotation-only (bit-identical by construction; verified import OK + 614 passed).

## ⭐ Q1 ADJUDICATED — KEEP THE STUB (self-typing is strictly WORSE), and the brief's premise is WRONG.

The brief claimed the `TYPE_CHECKING` apply-stub on the Mixin "introduced 2 NEW
reportIncompatibleMethodOverride" at fission.py:223 / scattering.py:962. **FALSE — they
PRE-EXIST at HEAD.** Measured (combined `npx pyright` on the 4 brief files, which is the ONLY
honest scope — per-file runs report 0 because cross-module type resolution is stubbed out):

| state | total | reportArgumentType | override | reportOperatorIssue |
|-------|-------|--------------------|----------|---------------------|
| HEAD (`ff497d7`) | 48 | 23 | **2** | 0 |
| WORK (stub) | **39** | 14 (−9) | **2** | 0 |
| stub REMOVED | 59 | 33 | 2 | 1 |
| self-typed dunders (stub removed) | **50** | 23 | **2** | 1 |

Net of the PR = **−9 reportArgumentType, +0 override, +0 anything.** The override (2) is INVARIANT
across all four states → it is a pyright `@singledispatchmethod`→`singledispatchmethod[NoReturn]`
modelling limitation orthogonal to this PR; S/F's `@singledispatchmethod apply` "incompatibly
overrides" ANY typed `apply` the ancestry carries (at HEAD the conflict resolves against the
Protocol member; the stub does NOT add it). Neither the stub NOR self-typing clears it — out of
scope, file separately if wanted (NO ignore).

**Why self-typing FAILS (prototyped on a temp copy, measured, then discarded — working tree
UNTOUCHED):** `self: "LinearOperator[V]"` on the composer-building dunders fixes only the dunder
bodies' internal `self`-passing; it does NOT reconcile the concrete leaves' `apply(q_ext:
TimedFullField)` with the Protocol (that reconciliation is exactly what `Generic[V]` + the
`["TimedFullField"]` binding + the stub provide). Self-typing REGRESSES reportArgumentType back
to 23 (the −9 q_ext cascade re-opens) AND introduces NEW errors the stub avoids:
`__rmul__`→`self.__mul__` (line 442 `Cannot access attribute "__mul__"` — narrowing `__mul__`'s
self breaks the `__rmul__` delegate), `__pow__`'s `result @ self` (reportOperatorIssue), and 5
NEW `SNMesh`/`from_mesh` arg errors in sn/operator.py. **50 > 39.** The stub is the cleaner
technique by 11 errors and zero collateral. ⭐ The stub is also the RIGHT elegance shape: it
declares the endomorphism contract the mixin's dunders RELY on (every concrete leaf honours it),
costs nothing at runtime (`if TYPE_CHECKING`), and its docstring states the no-runtime-apply
invariant honestly.

## NITS

- **NIT-1 (CONCERN, do-now, the one real gap) — the `Generic[V]` retype is INCOMPLETE inside
  iteration.py's own drivers.** `SourceIteration(Generic[V])`/`KrylovAcceleration(Generic[V])`
  got `solve(q_ext: V)->tuple[V,...]`, but the private next-iterate producer
  `_solve_with_seed(self, rhs, psi_prev)` (iteration.py:443) is UNTYPED → returns implicit
  `Unknown`, which unions into `psi` as `Unknown | V | NDArray` and then fails `V` at
  `rhs + g.apply(psi)` (line 506) and `return psi` (line 547); same shape at Krylov 751/753.
  **+3 NEW errors in iteration.py** (HEAD iteration.py 11 → WORK 14). Bug habitat: a half-typed
  generic driver is a documented-falsehood of the same kind vector.py:13 indicts — the `V` flow
  is advertised but broken at the helper seam, and a future reader "fixes" it with an ignore (the
  user rejects ignores). Required: type `_solve_with_seed(self, rhs: V, psi_prev: V) -> V` (and
  check `_l2_norm`/`_flux_displacement_leaf` don't widen) so `psi` stays `V` end-to-end. This is
  the SAME incomplete-audit-scope lesson as step-1 (`eigenvalue.py` was the missed driver then;
  the iteration.py PRIVATE helpers are the missed surface now). Verify: WORK 5-file is still a net
  −3 (29→26) so it's not blocking, but it leaves iteration.py NOT clean under its own new generic.

- **NIT-2 (PASS, Q2) — `apply(self, x: V, /)` positional-only is CORRECT and elegant**, not a
  smell. The leaves spell the param `psi`/`moments`/`coefficients`/`phi`; `/` lets them satisfy
  the Protocol structurally without name-matching (a Protocol with a named param would force every
  leaf to rename to `x` — that WOULD be the smell). Applied consistently (Protocol + all composer
  applies + ProjectionOperator/ReconstructionOperator ABCs). Mirrors std-lib `__add__(self, value,
  /)`.

- **NIT-3 (PASS, Q3) — the pre-existing Mixin ignores left in place are CORRECTLY out of scope.**
  `__pow__` `[assignment]`/`[return-value]` (operator.py:540/541), `__call__` `[attr-defined]`
  (:517), `adjoint` `[arg-type]` (:572) are all `self`-is-a-bare-Mixin-not-yet-a-full-operator
  artefacts unrelated to the `V` retype; clearing them needs a separate self-type pass that (per
  the Q1 prototype) regresses other counts. Leaving them is the right carve boundary (don't expand
  scope mid-PR). They are tracked in [[issue-256-vector-protocol-step1]]'s campaign.

## Q4 (reads-like-domain) — PASS, exemplary.
`(L+C−S−F−B)ψ=q` over ONE Vector carrier is exactly what the `Generic[V]` single-`V` (not
per-operator `Generic[F]`) expresses; vector.py's module docstring names the endomorphism honestly
and the `apply` docstring cites the carriers. The single `V` is the honest call: every operator
maps a carrier to a same-typed carrier; the inner leaf changes ROLE (flux→source) at runtime but
the vector TYPE does not — so one `V`, not N. No new anti-pattern.

## Step-1 follow-through CONFIRMED LANDED.
The step-1 BLOCKING fix (add `Vector.__truediv__`) shipped in `cfb651b` (vector.py:134, `Fψ/k` +
`ψ/p` cited) AND all 4 Vector dunders were sharpened to `-> Self` (self-typed carrier narrowing —
the axis-4 refinement my step-1 note called a non-issue is now belt-and-suspenders correct).

## Gate state (annotation-only ⇒ reds are pre-existing baseline, NOT this PR).
614 passed; 7 failed = the documented baseline reds (SPHERE snapshots #232/#250 +
`ymin`/`mu_y` ny=1 2-D); diff touches ZERO numerical kernel (sn/operator.py change = 3 class-decl
generic bindings only). The +1 ULP DD DriftWarnings are within nulp=5 tol.
