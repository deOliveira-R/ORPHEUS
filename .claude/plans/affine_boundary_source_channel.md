# The affine boundary source channel — `γ₋ψ = L γ₊ψ + q`

> **STATUS: ⏹ COMPLETE — every phase landed 2026-08-05/06.**
>
> | phase | commit |
> |---|---|
> | P1 `from_specs` | `9dfddeaf` |
> | P2′ `from_mesh_laws` + RHS wiring | `48657072` |
> | **P3** retire the affine operator (ERR-075) | `8d552395` |
> | **declaration channel** — a law is a legal `bcs=` entry | `985497b5` |
> | **P4** the §4.6 MMS through a declared law | `c86df8fb` |
> | **P5** linearity of `B` and of the full matvec | `d3066188` |
> | **P6** step 1 — trace tiers carry their directions | `16835c68` |
> | **P6** — `evaluate(space)`; the probe retired | `6ecb637a` |
>
> Verify any hash with `git merge-base --is-ancestor <hash> HEAD`. Follow-ups
> live as GitHub issues (**#331** the displacement-domain disagreement,
> **#332** the `principles.rst` redirect backlog), not here.
>
> **Read ⏸ COMPACTION POINT #2 at the END of this file** for the measured
> baselines, gate costs and durable lessons. **#2 supersedes #1 where they
> disagree, and they do** — §1's "fenced twice over / prophylactic" premise was
> refuted, see the ⛔⛔ REFUTED block below.
>
> ⭐ **The campaign's one-line result:** a boundary law is `γ₋ψ = Lγ₊ψ + q`,
> the realizer tier carries `L` **only**, and vacuum and prescribed differ
> **solely in `q`** — now true at the leaf, at the assembled `B`, at the full
> matvec, and in the declaration a user writes.
> Branch `refactor/operator-strategy-layers`.
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

### P5 — the matvec/linearity gates ⭐ explicitly requested ✅ LANDED 2026-08-06

> ⛔ **The bullets below are the ORIGINAL brief and one of them is
> unspellable.** Kept because the correction is the phase's main lesson; read
> the "as built" block underneath for what actually landed.

* `B.apply` is **linear for every law**, prescribed included: `B(0) = 0`,
  `B(2x) = 2B(x)`, `B(x+y) = B(x)+B(y)`. Today prescribed fails all three at
  `q ≠ 0` — this is the gate that would have caught the affine weld.
* The affine content lives in `q` and only in `q`: same solve via
  (zero-operator + `q_ext.boundary`) ≡ the pre-carve (affine operator + zero
  `q_ext.boundary`), on a fixture where both are reachable.
* The Krylov matvec stays linear with a prescribed law installed.

#### ⛔ `B(x+y) = B(x)+B(y)` IS UNSPELLABLE — the carrier is an AFFINE space

Written as specified, the additivity row raises before it can assert::

    TypeError: cannot add two AngularFlux states: flux states form an affine
    space with no origin, so '+' between two fluxes is undefined
    (Σλ = 2 lands off the affine subspace).

Flux states are a torsor over a distinct displacement space (the #208 carve,
`orpheus/transport/fields/_flux_role.py`); `ψ ⊖ ψ → Δψ` and `ψ ⊕ Δψ → ψ` are
the legal moves. Scalar scaling is untouched, so homogeneity survives verbatim.
`B`'s CODOMAIN, by contrast, is a vector space — it returns rate densities
(`AngularSourceSink ⊕ AngularBoundarySourceSink`), so output differences are
ordinary subtraction. ⟹ **`B` is an affine map from an affine space into a
vector space**, and the honest third law is base-point independence of the
increment:

.. math::  B(\psi_1) - B(\psi_2) = B(\psi_1 \oplus \sigma) - B(\psi_2 \oplus \sigma)

which says the increment depends only on `ψ₁ ⊖ ψ₂` — the same content as
additivity, stated without ever naming `B`'s induced tangent map (unspellable:
`B.apply` refuses a displacement argument — **#331**, filed).

#### As built

**`tests/sn/operators/test_declared_law_is_linear.py`** (new, 14 rows) on
prescribed(`xmin`) + **reflective**(`xmax`), het 2G GL-8, declared on the
GEOMETRY:

* §0 the activation guard, **decomposed per face** — `xmin`'s whole slot
  exactly `0.0`, `xmax`'s inflow rows `1.8239310798528774`. A bare
  `|B(x)| > 0` would not have been attributable.
* §1 `B(0)=0` · homogeneity at `c ∈ {3.7, −2, 1e3, 1e−3}` · base-point
  independence. `[M]` **`B` is EXACT** — `0.0` at every scalar; the
  base-point row reads 4–8 nulp, which is the re-basing's own rounding
  (`fl(a+σ) − fl(b+σ) ≠ a − b`), not the operator.
* §2 the same three for `A = (L+C) − S − B` off `build_within_group_system`.
  `[M]` worst relative `3.7e-16` ⟹ `_MATVEC_RTOL = 1e-14` (27× headroom,
  twelve orders against the affine regression's `1.6e-1`).
* §3 ⭐⭐ **the campaign theorem at the assembled tier**: prescribed and vacuum
  give **bit-identical** `B` *and* `A` (`|Δ|_inf = 0.0` at `|A(ψ)| = 42.685`)
  while `q_∂` reads `2.5` vs `0.0`. Both halves in one row — without the
  second, a channel that dropped the declaration passes the first perfectly.

**`tests/sn/operators/test_g_adjoint_reciprocity.py`** (extended, +6 rows) —
the `is_adjointable` WIDENING, added to the family that already owns
reciprocity rather than as a second spelling. Two new `_BUILDERS` cases
(`slab_declared_prescribed_2g`, `…_white_2g`) + one `_FULL_LOSS_BUILDERS` case.
A partner face is MANDATORY: the zero morphism is metric-blind (`0ᵀ = 0` under
every metric), so prescribed alone is a provable non-catcher.

#### `[M]` The mutation battery — 66 rows, baseline green

| mutation | reds | reading |
|---|---|---|
| `control` — `B.apply` scaled by `\|ψ_trace\|_inf` (NONLINEAR) | **27** | positive control (vv #17); every homogeneity + both base-point rows + every reciprocity case |
| `affine` — the P3 regression (`+q` inside the operator) | **18** | `B(0)`, `A(0)`, all 8 homogeneity, §3, §0, the 3 NEW reciprocity cases, the P3 trace gate. **The 9 pre-existing reflective reciprocity cases stay GREEN** ⟹ the new ones are attributable |
| `identity` — `L := I` for prescribed (perfectly LINEAR) | **5** | every linearity row stays green **as predicted**; caught by §0, §3, and the P3 trace gate |

⭐ The `affine` column's **misses** are the informative part: neither
base-point row reddens, exactly as their docstrings say *in advance*. A
coverage audit counting them as ERR-075 catchers would have been wrong by two
rows — `vv` #18 working as intended.

#### ⭐ What the QA pass found that I had missed (and why)

A fresh-context `qa` review of the above returned 8 findings; all were acted on
in the same commit. The two that generalise:

* ⭐⭐ **My own correction pass left a third falsified claim standing, ~470 lines
  from one it had just corrected to say the opposite.** The
  `SNBoundaryOperator` class docstring's re-emission paragraph reads "a diffuse
  one to the same Lambertian white uses **(not adjointable)**" — `[M]` white,
  diffuse albedo and specular albedo ALL report `True`. It survived because
  "white" and "(not adjointable)" sit on ADJACENT lines and my grep was
  line-based. Sharpest detail: **that paragraph was itself a correction pass** —
  it ends "it is the enumeration in prose that had to stop naming classes" — and
  the fix had been applied to the subjects but not to the parenthetical verdicts
  beside them. Now `vv-principles` **#21**.
* ⭐⭐ **"+6 rows" over-counted the reciprocity work by 2×.** 3 cases × the
  module's parametrized consumers = 6 rows, but 3 have bodies that cannot see a
  boundary law at all: the metric-reference row builds `G` from
  `volumes`/`weights`/`omega_dot_n` and never reads `sn.bc` (`[M]` bit-identical
  to `slab_2g`), and the one-hot row zeroes the whole trace block where `B`'s
  range lives. Green under four mutations. Now `vv-principles` **#20**.
* ⭐ A third, **#19**: I argued in prose that a metric-loaded partner face is
  mandatory, then cited the TRUE-metric residual (`1.7e-15`) as evidence it was
  loaded — the one number that carries zero information about loading. The
  wrong-metric readings are `2.4e-01` / `1.4e-01`, and the two new cases were
  never added to the committed wrong-metric control. Fixed both.

Also fixed: a `[M]`-badged false rationale in `test_reemission_closure.py`
(`LambertianReemission.is_adjointable` is `True`, not `False` — and its early
return therefore stopped firing, silently taking that test's value leg live); an
invented `≈ 1.6e-1` regression magnitude that is actually the fixture's
boundary/interior ratio (true reading `4.3e-2`); a silently-dropped `bc_left` on
the non-Cartesian arm of `_full_loss_case` (now raises — accepting and ignoring
a declaration is this campaign's own opening defect); and three nit-level `[M]`
numbers. Plus one row ADDED, `test_an_all_prescribed_mesh_makes_B_the_zero_morphism`
— this module's premise for rejecting the MMS fixture was prose only.

⚠ The battery harness broke a **third** time during that review, in the
reviewer's own hands, the same zsh way (`pytest $VAR` → one bogus path → 0
collected, exit 0, 0.01 s, reads as "0 caught").

#### Two falsified doc claims fixed on sight (measured, not inferred)

`test_sn_boundary_operator.py`'s header ("the `is_adjointable` conjunction —
**white would drop it**", from `d7e13164`, 2026-06-03) and
`SNBoundaryOperator.apply_transpose`'s docstring ("**the white BC has no
Euclidean transpose**") were both present-tense FALSE: `[M]` a
`WhiteBoundary()` face reports `is_adjointable = True` and `apply_transpose`
returns, since B3.4b factored the Lambertian's re-emission. The second
contradicted `is_adjointable`'s own correct past-tense note **twelve properties
above it in the same class**. Discriminated by scope, NOT blanket-replaced:
`RadialCharacteristicBoundaryOperator` (`B_b`) genuinely does defer white — a
different predicate about the sphere's off-quadrature μ = ±1 ray corner — and
its wording stands.

### P6 — a directional inflow is declarable without smuggling ⭐ the doctrine's second half

> ⛔ **RETITLED 2026-08-06. This phase read "promote what P4 hand-rolled into
> production" and that title cost a design cycle.** It named a MECHANISM (move
> the class) when the goal is an OUTCOME, and the agent that wrote it was the
> agent it later misled — reading its own words back, after compaction, as a
> specification. Moving the class would have been *wrong*: a manufactured
> solution is legitimately test-owned, and promoting it would put a
> verification artefact in the production tree. The gap was never the class's
> location; it was that the INTERFACE forced the class to smuggle trace
> knowledge through its constructor. This is now the worked example in
> `.claude/rules/plan-authoring.md` §1.
>
> **Goal.** A user can declare `q(Ω)` on a face and have it delivered, without
> the source carrying trace knowledge it should not need.
>
> **Done when** `_ManufacturedFaceInflow` is an ordinary user-written source —
> concretely, when it no longer takes `mu_inflow=` or `n_ordinates=`, because
> the interface supplies both.

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

---

## ⏸ COMPACTION POINT #2 — P3 + the declaration channel + P4 landed; P5/P6 remain

A session picking up cold re-anchors from **this file + `git log`**, never from a
conversation summary. Supersedes ⏸ COMPACTION POINT #1 where they disagree —
and they DO disagree, because #1's §1 premise was refuted (see the "⛔⛔ REFUTED"
block in §1).

### State

**HEAD `c86df8fb`** on `refactor/operator-strategy-layers`. Verify any hash with
`git merge-base --is-ancestor <hash> HEAD`.

| phase | commit | what |
|---|---|---|
| P1 | `9dfddeaf` | `from_specs` — the recipe→snapshot bridge |
| P2′ | `48657072` | `from_mesh_laws` + the RHS wiring; the user-path doctrine |
| — | `ef4c3537` | checkpoint #1 |
| **P3** | **`8d552395`** | retire `IncomingSourceOperator` onto the zero morphism; `ZeroOperator` gains `domain`/`codomain` (G6.3 step 6 folded in); ERR-075 |
| **channel** | **`985497b5`** | a `BoundaryTraceLaw` is a legal declaration wherever a `BC` tag is |
| **P4** | **`c86df8fb`** | the §4.6 non-vacuum MMS re-routed through a DECLARED law |

**REMAINING: P5, P6.** Both specified below. Nothing blocks either.

### ⛔ Read these three corrections before trusting anything older

1. **P3 fixed a LIVE bug, not a hypothetical.** §1's "fenced twice over /
   prophylactic architecture" was wrong: `SNBoundaryOperator._face_laws` has no
   `block_role` filter, so the affine operator reached `B`. SI double-delivered
   (ratio `2.000000`); **Krylov RAISED** `ConvergenceCertificateError` on both
   sides of P2′ and had been unusable all along.
2. **The existing MMS suite is NOT uniformly vacuum-nulling.** 4 of 13 builders
   are non-vacuum by design (anisotropic `(A_g + μ_n B_g)/W`, SymPy provenance,
   live L1 consumer). P4 is therefore a **re-routing**, not a new ansatz. A
   claim generalised from the 1-D-1G case's docstring cost a whole design cycle.
3. **The sphere is NOT refused** for a declared prescribed inflow — declare,
   materialise and solve all succeed; `_reflect_corner` fires only when the
   ray-corner action is *invoked*. The real restriction is `is_adjointable=False`.

### ⭐ P5 — ✅ LANDED 2026-08-06. The spec below was RIGHT about the fixture and
### WRONG about one law; §2's "P5 — as built" block is the record.

⛔ The `B(x+y) = B(x)+B(y)` row it specifies **cannot be written**: the flux
carrier is an affine space and refuses `flux + flux`. The landed third law is
base-point independence of the increment. Full account in §2's P5 block; the
architecture gap it exposed is **#331**.

**The whole point: on the P4 MMS fixture every linearity row is a TAUTOLOGY.**
Both slab faces are prescribed ⟹ after P3 both realize to the zero morphism ⟹
`|B(x)|_inf = 0.0` for random `x`. `B(0)=0` and `B(2x)=2B(x)` then hold because
both sides are structurally zero; no input can red them (`vv` Mode 8).

⟹ **P5's fixture is prescribed(`xmin`) + REFLECTIVE(`xmax`)**, het 2G slab, GL-8.
`[M]` there: `|B(x1)|_inf = 1.3201645939238549` (the activation guard), `|B(0)| = 0`,
`|B(c·x) − c·B(x)| = 0` exactly; both drivers converge, `γ₋(xmin) = 2.500000000000`
exactly; SI 0.05 s, Krylov 0.24 s.

P5 owns, all on that fixture, each with the activation guard:
* `B.apply` linearity — `B(0)=0`, homogeneity, additivity;
* full matvec `A = (L+C) − S − B` linearity with a declared law. On the MMS
  fixture the prescribed contribution to `A` is structurally zero, so reds would
  come from `L+C−S` — unattributable (`vv` #18). The mixed fixture is what makes
  them attributable.

**ALREADY LANDED, do not duplicate** (design `scratch/p4_mms_design.md` §9.3):
the affine-content-lives-in-`q` claim (P3's `test_the_two_user_paths_reach_the_same_fixed_point`,
bit-identical); the leaf-level zero-morphism rows; the zero-morphism structural
row and the Krylov-reproduces-the-MMS row (both in P4's module).

Also open from `scratch/p3_verification_plan.md` §4.2: a row that
`SNBoundaryOperator.is_adjointable` is now `True` with a declared prescribed face
(a capability WIDENING — pre-P3 one prescribed face poisoned the whole block's
transpose). ⚠ Its warning: the zero morphism is the most metric-blind law there
is (`0ᵀ = 0` under every metric), so a reciprocity gate on prescribed **alone**
is a provable non-catcher. Pair it with a Lambertian/white face.

### P6 — two items, both with measured specifications

1. **The `InflowSourceSpec` shape.** `_ManufacturedFaceInflow`'s constructor list
   IS the spec (`test_mms_declared_inflow.py`, and `p4_mms_design.md` §10): the
   per-row `μ` **in trace order**, the face's own coordinate, the group axis, the
   `1/W`. ⟹ "a source that receives the TRACE and the FACE", not a bare shape.
2. **The ERR-047 presence probe is opt-out-able.**
   `assert_source_lives_on_incoming_trace` opens
   `probe = source.evaluate((N,)); if not np.any(probe): return` — `[M]` a spec
   returning zeros at rank-1 and `7.0` at rank-2 realizes cleanly with the
   certification SKIPPED. A presence predicate a source can decline is not a guard.

### Gate costs and the red baseline, measured

* **Wide** `tests/{sn,transport,geometry,numerics,diffusion} -m "not slow"` =
  **`7 failed / 5904 passed`, ≈17–18 min.** The 7 are pre-existing:
  2× `test_streaming_operator` `cart2d_*_principled_equiv`, 3×
  `test_affine_carve_bit_identity` (that file is **#208's** carve, not this
  campaign's — its `sha256` is measured unchanged by P3), 1×
  `test_diamond` `spherical_inward_bit_identical`, 1×
  `TestWhiteXminPartial03GLSnapshot` (task #33).
* `tests/sn/verification + tests/sn/solve + tests/transport` ≈ 7 min.
* `tests/geometry + tests/transport + tests/diffusion + tests/sn/operators` ≈ 46 s.
* Static: `tests/test_pyright_ratchet.py` (NOT bare `npx pyright` — it scans
  `scratch/`) · `tools/check_docstring_xrefs.py orpheus tests docs` → `DEAD
  TARGETS 0` · `sphinx -E -W` → 0 warnings · `python -m tests._harness.audit` →
  0 orphan equations, ERR coverage 74/75 (the gap is pre-existing ERR-074).

### ⭐⭐ Durable lessons from this stretch

1. ⭐⭐ **"It is FENCED, so it cannot happen" is a claim about a CONSUMER and must
   be measured AT the consumer.** Both fences were read off the producer; one was
   genuine, the other an unchecked inference about 14 lines nobody opened. The
   check was one line (`B.apply(0)`) needing no oracle. **When a design note says
   "X is excluded because it lacks Y", grep for who READS Y.**
2. ⭐⭐ **A mutation battery reporting 0 caught is a broken INSTRUMENT until
   proven otherwise — and it broke twice this stretch.** (a) `grep "^FAILED"`
   matches nothing under `-q --tb=no`, and ANSI colour breaks `^` even with `-rf`
   → use `--color=no -rf`. (b) A pytest plugin patching a module it imported
   itself in `pytest_configure` patches a DIFFERENT object from the one pytest
   collected: 0/6 caught while the mutations were provably live (ratio 2.0
   in-process). **Patch `item.module` in `pytest_runtest_setup`.** Fixing the
   harness alone took the P4 battery 0 → 5 reds with no test change.
3. ⭐⭐ **Read a sub-agent measurement's CONFIGURATION, not just its value.** I
   "corrected" a right statement (channels are bit-identical) using a **pre-carve**
   number (`1.998e-13`, operator-route vs channel). Post-carve it is exactly
   `0.0`. Three real numbers, three different comparisons; the wrong one would
   have justified an `rtol` gate blind to the defect.
4. ⭐ **zsh does NOT word-split unquoted variables.** `pytest $FILES` passes ONE
   bogus path → 0 collected, exit 0, 0.01 s, reads green. Always read the count.
5. ⭐ **A ULP budget coupled to `inner_tol` is not a floating-point claim.** The
   Krylov trace deviation scales with the solver tolerance: 18 ULP at `1e-13`,
   **27 580** at `1e-10`. Still spell it in ULP, not `rtol`.
6. ⭐ **pyright as architecture-smell detector, again.** Widening the declaration
   source without the axis SINK left 8 new `transport` errors — it found the
   half-done change before any test could.
7. ⚠ **`docs/theory/verification/matrix.rst` is GENERATED by the Sphinx build,
   so a commit that adds tests without committing a rebuild leaves it stale —
   and the staleness is invisible in a green tree.** `[M]` at P5 it was stale by
   **three** commits (`985497b5`, `c86df8fb`, `30637db7`, all of which added
   rows): 8774 → 8812 tests, +15 L1 and +23 foundation, five modules missing or
   under-counted. It reads as an authoritative V&V inventory while under-
   reporting the coverage that exists. Rebuild Sphinx and `git add` the matrix
   in the SAME commit as any test addition.
8. **A "forbidden to commit" note is a point-in-time snapshot.** The
   `vv-principles` prohibition was stale (committed at `34af8474`); verified with
   `git merge-base --is-ancestor` + a content grep before appending ERR-075.
