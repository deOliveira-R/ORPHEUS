# The affine boundary source channel — `γ₋ψ = L γ₊ψ + q`

> **STATUS: PLANNED, not started.** Branch `refactor/operator-strategy-layers`.
> Arose from G6.3 **step 6** (`.claude/plans/geometric_transformation_consolidation.md`
> §7h order table) and outgrew it. That plan's step 6 now means "bind `ZeroOperator`'s
> spaces"; everything else about prescribed inflow lives here.
>
> **User ruling (2026-08-05):** build the channel — *"not because of trigger, but
> because if we build it, maybe next time we're working on BC, you will be able to
> connect the prescribed inflow machinery to the source at boundary. And let's use it
> in a non-trivial MMS, which is an excellent smash test for the numerical schemes.
> Test the matvec path as well."*

---

## 1. The finding — the design of record is already written; only half is built

The affine boundary law is

```
γ₋ψ = L γ₊ψ + q
```

and every law is one linear chain `L` plus one source `q`:

| law | `L` | `q` |
|---|---|---|
| vacuum | `0` (the zero morphism) | `0` |
| **prescribed** | **`0`** | **≠ 0** |
| specular | `P` | `0` |
| white / Lambertian | `B @ C` | `0` |

⭐ **Vacuum and prescribed differ ONLY in `q`.** That is the whole content of the
carve, and it is what makes `IncomingSourceOperator` unnecessary rather than merely
untidy.

**Two independent sites already state this as the design**, and both were written
before this session:

* `orpheus/diffusion/boundary_realizer.py:103-110` — *"`PrescribedInflow` is the
  rank-0 AFFINE law `J⁻ = q`: its realization is the boundary **source**
  `q.boundary`, not a linear boundary operator `B` (**exactly the SN split** — SN
  realizes it as an `IncomingSourceOperator` and does NOT stamp it BOUNDARY)."*
  Diffusion **REFUSES** the law outright on these grounds.
* `orpheus/numerics/operator.py:385-386` — the same claim, given as the reason the
  law carries no `BlockRole.BOUNDARY`.

`AngularBoundarySourceSink` is the typed `q`
(`orpheus/transport/source_sinks/angular_boundary_source_sink.py:1-17`: *"the `q`
term of the affine boundary law"*), packed into the `AngularTraceSpace` flat layout,
and `FullField` is `Composite[BulkField, BoundaryField]` — the bulk and boundary
source are one object, exactly as the user said.

### What is missing — the fence exists, the channel does not

`[M]` this session:

| claim | measured |
|---|---|
| the realized prescribed operator is AFFINE | `A(0) = 2.5`, `A(2x) − 2A(x) = 2.5` at `q = 2.5` |
| it is fenced out of the `B` block | `role=None`; vacuum/reflective are `BlockRole.BOUNDARY` |
| `q_ext.boundary` carries the source | ⛔ **NO** — `AngularBoundarySourceSink.zeros_on(...)` **unconditionally**, `solver.py:1590` and `:1781` |
| the comment above those two lines | ⛔ **present-tense-FALSE**: *"zero for vacuum/reflective; prescribed inflow otherwise"* |

⟹ **a nonzero prescribed inflow is built, returned, unstamped, and dropped.** The
source never reaches the channel the docs say it belongs in.

### ⚠ NOT a live bug — record this so it is not re-alarmed

The affine operator does **not** reach a Krylov matvec, fenced twice over:

1. `augmented_mesh.py:186-190` — `prescribed_inflow` is **not a registered `BC` kind**,
   so the law is *"declarable only by constructing the law directly, never from a
   `BC(...)` tag"* (#189). Only tests construct it.
2. The missing `BlockRole.BOUNDARY` stamp keeps it out of the `B` block.

So this carve is **prophylactic architecture**, not a bug fix. Say so in the commit;
a future reader finding "affine operator in a linear slot" will otherwise re-open it
as a defect.

---

## 2. Phases

⏸ **COMPACTION POINTS after P2 and P4.** Commit, then checkpoint this file.

### P1 — the bridge: `AngularBoundarySourceSink.from_spec(spec, mesh)`

Materialize an `InflowSourceSpec` (the lazy per-face *recipe*) onto the trace (the
eager whole-boundary *snapshot*), looping `spec.evaluate(face_shape)` per face and
packing the flat layout. The target docstring already specifies it
(`angular_boundary_source_sink.py:85-95`) and marks it deferred under
`feedback_unify_after_two_instances`.

⚠ **Retire that deferral note in the same commit** — it becomes present-tense-false
the moment the method exists. This is the campaign's most-repeated defect class.

Gates: round-trip `from_spec(ConstantInflowSource(v)) ≡ prescribed_inflow({face: v·ones})`
on the inflow slots and **zero elsewhere** (the outflow/tangential slots are the
claim that the mask dissolved into the type); `NoSource → zeros_on`.

### ⛔ P2 AS WRITTEN IS REFUTED — the channel already exists (2026-08-05)

P2 originally read *"route it: `q_ext.boundary` carries the source"*, targeting
`solver.py:1584` / `:1775`. **That premise was wrong, and this is the SIXTH time in
this campaign that a plan claim has been refuted by the realization.** The evidence
was in the corpus the whole time — `docs/theory/verification/sn.rst:3450-3458`:

> the inflow is carried in the `q.boundary` slot
> (`AngularBoundarySourceSink`) and consumed directly by `(L+C).solve` as the sweep
> inflow seed … the inflow is supplied as the boundary leaf of the **composite
> source** `q = q_bulk ⊕ q_∂` that `solve_sn_fixed_source` now accepts directly.

`[M]` `solver.py:3110` — `_build_fixed_source_rhs` normalises to `q = q_bulk ⊕ q_∂`
and **both** inner paths (`_solve_fixed_source_si`, `_solve_fixed_source_krylov`)
consume it. A non-vacuum MMS (§4.6) already drives it.

⟹ the two `zeros_on(...)` sites I flagged are the **eigenvalue** paths, where a zero
boundary source is CORRECT (a k-eigenvalue problem has no external inflow). Their
comments are misleading, but the code is right. **Fix the comments; touch nothing
else.**

**What is actually missing is narrower:** the **mesh-BC bridge** — a
`PrescribedInflow` *declared as a boundary law* never reaches the composite source.
`sn.rst:3900` says so in as many words: *"No `from_specs` / `PrescribedInflow`-BC
bridge is touched"* — the §4.6 MMS deliberately supplies `q_∂` directly instead.
That bridge is the one thing standing between the law tier and the channel, and it
is exactly what the user asked for: *"next time we're working on BC, you will be
able to connect the prescribed inflow machinery to the source at boundary."*

### P2′ — the mesh-BC bridge (REPLACES P2)

Materialise a face's declared `PrescribedInflow.source` into `q.boundary` via
`from_specs`, so a law-tier declaration reaches the existing channel. Scope this
against **#189** (prescribed_inflow is not a registered `BC` kind), which is the
reason the law is only constructible directly today — the bridge and the
registration are separable, and this plan owns only the bridge.

### P3 — collapse the operator; retire `IncomingSourceOperator`

`realize(PrescribedInflow)` returns the **zero morphism** — literally
`_narrowed_zero_operator`, the same object vacuum returns — stamped
`BlockRole.BOUNDARY` like every other law, because it is now genuinely linear.
Retire `IncomingSourceOperator` (three searches: graph callers, text-grep across
code/tests/docs, direct constructors).

Folds in G6.3 step 6's original scope: `ZeroOperator` gains `domain`/`codomain`,
bound in `_narrowed_zero_operator`.

⚠ **Two law-tier consumers must keep discriminating prescribed from vacuum** after
the operators become identical — both read the LAW, so both should survive, but
**verify rather than assume**: `sn/operators/boundary.py:175` (corner block, excludes
by FAMILY) and `sn/acceleration/dsa.py:227` (*"excluded by FAMILY, not by its current
q"*). Their own comments say family-not-value, which is exactly what makes them
survive — that is a prediction to check, not a reassurance to accept.

⚠ **`is_adjointable` flips.** `IncomingSourceOperator` reports `False` (measured);
the zero morphism is adjointable. Grep for gates pinning the `False`.

### P4 — the non-trivial MMS ⭐ the user's smash test

A manufactured solution with a **nonzero inflow on at least one face**, driven
through the new channel.

⛔ **vv-principles Mode 7 is the whole risk here.** The ansatz must ACTIVATE the
boundary source, not null it. State explicitly which terms the chosen ψ activates
and which it nulls; if ψ happens to vanish on the inflow trace, the MMS is blind to
exactly the thing it was built to test. Design the ansatz so `γ₋ψ ≠ 0` on the source
face and check that a mutation of `q` moves the measured error.

### P5 — the matvec/linearity gates ⭐ explicitly requested

* `B.apply` is **linear for every law**, prescribed included: `B(0) = 0`,
  `B(2x) = 2B(x)`, `B(x+y) = B(x)+B(y)`. Today prescribed fails all three at
  `q ≠ 0` — this is the gate that would have caught the affine weld.
* The affine content lives in `q` and only in `q`: same solve via
  (zero-operator + `q_ext.boundary`) ≡ the pre-carve (affine operator + zero
  `q_ext.boundary`), on a fixture where both are reachable.
* The Krylov matvec stays linear with a prescribed law installed.

---

## 3. Standing constraints

* ⛔ **NEVER** `git checkout <path>` / `git restore` / `git stash` / `git clean` on a
  tracked path — the tree carries irrecoverable uncommitted state. Revert by
  re-editing; compare via `git worktree add`.
* Mutation-test by in-process monkeypatch, never a git-level discard.
* `.venv/bin/python`; `python -O -m pytest`; **SERIAL**; always `-m "not slow"`.
* ⚠ `grep "^FAILED"` does NOT work: `-q --tb=no` emits no FAILED lines, and ANSI
  colour breaks the `^` anchor even with `-rf`. Use `--color=no -rf`. Cost two blind
  measurements on 2026-08-05.
* ⚠ Verify a test path exists before running it — a run collecting nothing exits 0
  in 0.01 s and looks green.
* Commit trailer: `Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>`.
* Known reds (4, none from this campaign): `TestWhiteXminPartial03GLSnapshot::
  test_matches_the_frozen_scaled_lambertian`; `test_cart2d_1g_vacuum_apply_principled_equiv`;
  `test_cart2d_2g_specular_apply_principled_equiv`;
  `TestBitIdenticalCurvilinear::test_spherical_inward_bit_identical`.
* Gate costs: `tests/geometry + tests/sn/operators` ≈24 s (`3 failed / 1747 passed`);
  wide slice ≈10 m 15 s (`4 failed / 4626 passed`).
