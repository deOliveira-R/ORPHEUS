# The affine boundary source channel — `γ₋ψ = L γ₊ψ + q`

> **STATUS: P1 + P2′ LANDED (HEAD `48657072`); P3–P6 remain. Read ⏸ COMPACTION
> POINT #1 at the END of this file FIRST — it carries the corrected red-set
> baseline (7, not 4 — the inherited wide slice omitted `tests/sn/solve`), the
> gate costs, and the durable lessons.** Branch `refactor/operator-strategy-layers`.
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

### ⛔⛔ REFUTED 2026-08-05 — it WAS a live bug, and fence (2) never existed

**The text below is kept as the falsified claim, because the way it was wrong is the
lesson.** It read:

> The affine operator does **not** reach a Krylov matvec, fenced twice over:
> (1) `augmented_mesh.py:186-190` — `prescribed_inflow` is not a registered `BC`
> kind (#189); only tests construct it. (2) The missing `BlockRole.BOUNDARY` stamp
> keeps it out of the `B` block.
> So this carve is **prophylactic architecture**, not a bug fix.

**Fence (1) is real. Fence (2) does not exist.**
`SNBoundaryOperator._face_laws` (`orpheus/sn/operators/boundary.py:311-324`) collects
`sn_mesh.bc[face]` for **every** face in the trace layout with **no `block_role`
filter**. `SNBoundaryOperator` carries the stamp itself and applies all of them, so
the unstamped affine leaf entered `B` exactly like any other law. `[M]`:

| claim | measured |
|---|---|
| `B` is linear with a declared inflow at `q = 2.5` | ⛔ `|B(0)| = 2.5`, `|B(2x) − 2B(x)| = 2.5` |
| SI, declared inflow, HEAD (post-P2′) | `γ₋(xmin) = 5.000000000000` — **delivered twice** |
| SI, declared inflow, pre-P2′ | `γ₋(xmin) = 2.500000000000` — worked, once, correct sign |
| Krylov, declared inflow, **both** sides of P2′ | ⛔ **RAISES** `ConvergenceCertificateError`, `‖Aψ−q‖/‖q‖ = 1.718` |
| vacuum control | `0.000000000000` |

⟹ **three corrections to this plan's own record:**

1. A declared `PrescribedInflow` was **never "silently inert."** On SI it was
   FUNCTIONAL — delivered once, with the correct sign, through `+Bψ` on the RHS.
2. **P2′ (`48657072`) is a double-delivery regression on the SI path**, ratio
   `2.000000` exactly against the vacuum control. On Krylov there was nothing left
   to regress: an affine `A(x) = A_lin(x) − c` breaks GMRES's Arnoldi relation
   `A V_k = V_{k+1} H_k`, so declared-prescribed × Krylov had been **unusable all
   along** and `_certify_within_group_exit` (`solver.py:435`) was catching it.
3. ⟹ **P3 is a BUG FIX, not prophylaxis.** The commit must say so. P2′ and P3 are
   two halves of one change and the tree is incorrect with only one of them landed.

⭐ **Why no gate saw any of it** — fence (1) is the whole answer. `prescribed_inflow`
is not a registered `BC` kind, so no production driver installs the law, and **no
test ran a full solve with a DECLARED law.** The P2′ gates stop at
`_build_fixed_source_rhs`: they correctly verify the RHS receives `q` and are
structurally blind to `B` also delivering it. That is `:ref:`verification-user-path``
failing on its own author one commit after landing it — the gates travelled part of
the user path and stopped before the solve.

⭐ **The transferable lesson: an "it is fenced, so it is not a bug" argument is a
claim about a CONSUMER, and must be measured at the consumer.** Both fences were
read off the *producer* (the realizer returns an unstamped operator; the `BC` registry
has no key). Fence (1) happened to be a genuine consumer fact. Fence (2) was an
inference about what `SNBoundaryOperator` *would* do with an unstamped leaf, never
checked against its 14 lines of code — and the check is one `B.apply(0)`.

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

### P3 — collapse the operator; retire `IncomingSourceOperator` ✅ LANDED 2026-08-05

**Verification plan of record: `scratch/p3_verification_plan.md`** (1033 lines,
`test-architect`) — read §0 for the measurements and §8 for the residual gap list.

**What landed:** `realize(PrescribedInflow)` returns `_narrowed_zero_operator(...)`
stamped `BlockRole.BOUNDARY`; `IncomingSourceOperator` retired; `ZeroOperator` gained
`domain`/`codomain` (G6.3 step 6 folded in), bound in `_narrowed_zero_operator` via
`checked_space_extent`. Post-carve: `|B(0)| = 0`, `|B(2x) − 2B(x)| = 0`,
`block_role = BOUNDARY`, `domain is Γ₊(f)`, `codomain is Γ₋(f)`.

⭐ **The keystone gate is `γ₋ψ|_f = q_f` on the converged answer**
(`tests/sn/solve/test_declared_inflow_reaches_the_rhs.py`), parameterized over
`{source_iteration, krylov}`. The boundary condition is a DEFINITION, so this needs no
reference solver and no discretization assumption, and it separates the three
outcomes exactly: `5.0` double / `2.5` single / `0.0` lost.

⚠ **Exactness differs by path, MEASURED.** SI writes the source into the inflow slot
and sweeps from it, so `γ₋ψ` is a COPY — bit-exact. Krylov reaches the trace through
the GMRES iterate and carries `2.500000000000008` = **18 ULP** at 2.5. The gate
asserts `array_equal` on SI and `assert_array_almost_equal_nulp(nulp=64)` on Krylov,
stated in ULP rather than `rtol` so the budget cannot quietly grow into the 1×/2× gap.

⭐⭐ **Why the leaf-linearity gate cannot replace the keystone — measured, not argued.**
Mutation battery over the 5 gate files (baseline `111 passed`):

| mutation | class | keystone | `B(0)==0` leaf gate | total reds |
|---|---|---|---|---|
| `q2` — double `q` in the source channel | **in-class** | 🔴 | 🟢 | 6 |
| `Lident` — realize `L := IdentityOperator` | **in-class** (linear!) | 🔴 | 🟢 | 10 |
| `affine` — reinstate the retired operator | out-of-class | 🔴 | 🔴 | 17 |
| `pcplus` — realize vacuum as identity (**control**) | control | 🔴 | 🟢 | 7 |

`Lident` is the discriminator: it is *perfectly linear*, so `B(0) = 0` still holds and
the leaf gate is blind, while the keystone reddens. And `affine`'s inflated 17 reds are
the vv anti-pattern #18 trap — it breaks LINEARITY, so most of its reds see the broken
law rather than the delivery count, and must not be credited as coverage.

⚠ Two `vv` rulings this battery confirms: candidate "φ is affine in `q` with the right
slope" is a **provable non-catcher** (a doubled delivery is `q → 2q`, still exactly
affine in `q` — Mode 12); and the `krylov` parameter is load-bearing because pre-P3 it
*raised* rather than answering wrongly.

⚠ **`_build_fixed_source_rhs`'s two arms take DIFFERENT bulk-source types** (a
per-ordinate array vs an `AngularSourceSink`). Every two-channel gate routes its bulk
through one helper (`_composite`) so the trap is unspellable; feeding one arm each way
gave `φ[0] = 3.083` vs `2.480`, a bulk difference that reads as a channel difference.

**Five docstrings stated the retired split as the design of record** and were corrected
(the campaign's own inventory found two): `numerics/operator.py:383`,
`diffusion/boundary_realizer.py:103`, `geometry/boundary/_realizer.py:124`,
`geometry/boundary/_bound_compat.py:218`, `sn/operators/boundary.py:283`.

**Also cleaned:** four imports in `sn/boundary/angular.py` (`AXIS_NAMES`, `face_name`,
`TANGENTIAL_EPS`, `build_omega_dot_n`) had been orphaned since **G6.3 step 3b** retired
`AngularAverageOperator` — a pre-existing retirement-audit miss, fixed while in the file.

**Open from `scratch/p3_verification_plan.md` §8** (not blockers for P3's landing):
the `SNBoundaryOperator.is_adjointable` widening row (§4.2, with its warning that a
reciprocity gate on prescribed alone is a provable non-catcher — the zero morphism is
metric-blind, so the honest gate pairs prescribed with a Lambertian face); ERR-047's
catalog catcher paths; and whether the double delivery earns its own `ERR-NNN`
(mode #6 convention drift: the definition site said "unstamped ⟹ not in `B`" and the
usage site never filtered). ⛔ `error_catalog.md` is forbidden to commit, so both
catalog items are handoffs.

#### P3 as originally specified

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

### ⭐⭐ THE RULING that orders the rest (user, 2026-08-05)

> **Tests must route through the machinery that a user would exercise without
> bypassing code functionality. Or else it's not testing the path the users go
> through.**

This is V&V doctrine, not a preference about this MMS, and it belongs in the
corpus (`docs/theory/verification/principles.rst`) alongside the ladder. It
disqualifies the §4.6 MMS's shape as a model for P4: that MMS hand-builds `q_∂`
and supplies it directly, *deliberately* bypassing the law tier
(`sn.rst:3900`). It verifies the CHANNEL and is silent about everything a user
actually touches to get an inflow.

**Three consequences, in order:**

1. ⛔ **P2′ BLOCKS P4.** An MMS cannot route through a bridge that does not
   exist. The old plan had P4 independent; it is not.
2. ⭐ **The MMS cannot use a `BC(...)` tag, and that is structural, not a gap.**
   `BC.params` is `dict[str, float]` (`geometry/mesh.py:59-64`), so a tag can
   carry `{"albedo": 0.7}` but never a manufactured solution restricted to a
   face. A non-trivial prescribed inflow is *inherently* not tag-expressible ⟹
   the user path here is constructing `PrescribedInflow(source=<spec>)` and
   installing it, which is a public surface (`orpheus.geometry.boundary`
   exports it). #189 (registering the kind) is a SEPARATE, weaker question
   about the constant case and is NOT a prerequisite.
3. ⭐ **The MMS fires the deferral trigger retroactively.** `from_specs` waited
   for *"the first real consumer that both declares a non-trivial
   `InflowSourceSpec` AND drives a sweep that consumes a typed boundary-source
   field"*. That is precisely P4. Built ahead of the trigger by ruling; the
   trigger then arrives two phases later and vindicates the ordering.

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

### P6 — promote what P4 hand-rolled into production ⭐ the doctrine's second half

**Added by user ruling, 2026-08-05**, refining the P4 ruling. (Numbered P6
because P5 — the matvec/linearity gates — was already in this plan.)

> If the MMS is exercising custom machinery that is not part of production, then
> the proper shape of the machinery should be implemented so that the MMS can
> use it. **It is a sign of a gap.**

P4's manufactured inflow needs an `InflowSourceSpec` that evaluates a
manufactured solution on a face. Implementing the Protocol in the test module is
*using* the machinery, not bypassing it — but it is a **stopgap with an owner**,
because production offering only `NoSource` and `ConstantInflowSource` is
precisely the gap. Nothing shipped can express an inflow that varies in angle,
space, or energy, which is the whole content of a non-trivial boundary
condition.

P6 promotes the shape P4 needed. **Do not design it before P4** — the honest
production shape is whatever the MMS turns out to require, and guessing it first
is how a speculative abstraction gets minted. Let the test state the need, then
build to it and retire the private version.

Landed doctrine: :ref:`verification-user-path` in
`docs/theory/verification/principles.rst` carries both halves, with this boundary
source as the worked example.

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

---

## ⏸ COMPACTION POINT #1 — P1 + P2′ landed; P3–P6 remain

A session picking up cold re-anchors from **this file + `git log`**, never from a
conversation summary. Companion checkpoint: `geometric_transformation_consolidation.md`
⏸ COMPACTION POINT #2 (that campaign's G6.3 steps 6–8 are still open and step 6
now means only "bind `ZeroOperator`'s spaces").

- **HEAD `48657072`** on `refactor/operator-strategy-layers`. **Tree clean.**
  Verify any hash with `git merge-base --is-ancestor <hash> HEAD`.
- **What landed:**
  | phase | commit | what |
  |---|---|---|
  | plan | `4a23c18a` | the plan of record |
  | P1 | `9dfddeaf` | `from_specs` — the recipe→snapshot bridge |
  | P2′ | `48657072` | `from_mesh_laws` + the RHS wiring; the user-path doctrine in `principles.rst` |
- **REMAINING: P3, P4, P5, P6** (§2). Order is forced: P3 → P4 → P6, with P5
  attachable any time after P3. **P2′ blocks P4** (an MMS cannot route through a
  bridge that does not exist) and **P4 blocks P6** (do not design the promoted
  shape before the MMS states what it needs).

### ⛔ THE RED-SET BASELINE WAS WRONG BY A WHOLE DIRECTORY — read this first

The "wide slice" this campaign inherited —
`geometry + numerics + sn/{operators,sweep,architecture} + diffusion` —
**OMITS `tests/sn/solve`**, and therefore reported **4** known reds when the true
count is **7**. The three it never saw:

* `tests/sn/solve/test_affine_carve_bit_identity.py::
  test_converged_flux_bit_identical_after_affine_carve` — `[si_2d_p1_aniso_het]`,
  `[krylov_2d_p1_aniso_het]`, `[si_slab_2g_het]`

⚠ **They read as a P2′ regression and are not one.** `[M]` with
`from_mesh_laws` monkeypatched back to `zeros_on` (the pre-P2′ behaviour) all
three STILL fail. Count check: 2 (operators) + 3 (solve) + 1 (sweep) = **6 in
`tests/sn`**, matching the "6 deliberate reds, documented in `81689a58`" that
project memory records, plus the 7th in `geometry`.

⟹ **Run `tests/sn` whole, not a hand-listed subset.** Measured cost:
`tests/{sn,transport,geometry,numerics,diffusion} -m "not slow"` = **17 m 35 s**
at `7 failed / 5882 passed`. The narrow slice is 10 m 15 s and under-covers.
⚠ Also note `tests/transport` was absent from the inherited slice too.

### Gate costs, measured

`tests/transport + tests/geometry + tests/sn/operators` ≈45 s
(`3 failed / 2205 passed`) · the full slice above ≈17 m 35 s
(`7 failed / 5882 passed`). Static: `tools/check_docstring_xrefs.py orpheus tests
docs` → `DEAD TARGETS 0` · `tests/test_pyright_ratchet.py` (NOT bare `npx
pyright`, which scans `scratch/`) · `sphinx -E -W` → 0 warnings.

### ⭐ Durable lessons from P1/P2′

1. ⭐⭐ **The user-path doctrine is now CORPUS, not a plan note** —
   :ref:`verification-user-path` in `docs/theory/verification/principles.rst`,
   both halves plus two corollaries. Do not re-derive it here; cite it.
2. ⭐⭐ **TWICE in one session my own scope claim was refuted BY THE CORPUS.** I
   proposed the affine carve as unbuilt — `diffusion/boundary_realizer.py:103`
   and `numerics/operator.py:385` already stated it as the design of record. I
   then wrote P2 as "route `q_ext.boundary`" — `sn.rst:3450` already described
   the composite channel, and `solver.py:3110` already built it. ⟹ **before
   proposing a carve, grep the corpus for the thing you are about to claim is
   missing.** The docs are the design of record and they were right both times.
3. ⭐ **A convenience factory can BE a bypass.** `prescribed_inflow` is public,
   production, correct — and a driver using it still skips the declaration tier.
   "It is a production symbol" does not establish that the user path was
   travelled. (Now recorded as a corollary in the corpus.)
4. ⭐ **A silently-inert declaration is invisible to a green suite** and the
   inertness is *caused by* the missing wiring the tests bypassed. The two
   defects protect each other: the bypass is why no test could see the
   inertness, and the inertness is why the bypass looked harmless.
5. **Refuse a double specification; do not resolve it.** Adding double-counts,
   overriding makes an input a silent no-op. Both wrong answers are quiet.
6. ⚠ **`grep "^FAILED"` failed TWICE, both from the anti-pattern-#17 list:** no
   FAILED lines are emitted under `-q --tb=no`, and ANSI colour breaks the `^`
   anchor even with `-rf`. Use `--color=no -rf`. Each empty result read as
   "no reds".
