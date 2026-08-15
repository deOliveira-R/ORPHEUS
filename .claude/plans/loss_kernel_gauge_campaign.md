# The reflective trace has a canonical value, and the solver returns it

> **Plan of record, 2026-08-14.** Supersedes the Q6/track-B content this file
> previously held (that work is COMPLETE — `e17116f9`…`83aa324e`, see
> `.claude/plans/q6_selection_criterion_and_track_b_finish.md`).
> Branch `refactor/track-b-remainder`, HEAD `5c39226b`, 11 commits ahead of `main`.

---

---

## ⏸ COMPACTION POINT — 2026-08-14. Steps 0, 2, 3 LANDED. Step 1 WITHDRAWN.

**Re-anchor from THIS FILE + `git log`, never from a conversation summary.**
Branch `refactor/track-b-remainder`, **12 commits** ahead of `main`, ⚠ not
merged — verify with `git merge-base --is-ancestor <hash> HEAD`.

| step | state | where |
|---|---|---|
| **0** — does the finalize sweep carry `ker A`? | ✅ **YES, it does** | `scratch/probe_344_step0_finalize_sweep.py` |
| **1** — single-source `Solution` construction | ⛔ **WITHDRAWN**, premise refuted | see §Step 1 below |
| **2** — the derived predicate (both halves) | ✅ LANDED `5def63b0` | `transport/spatial/scheme.py`, `sn/mesh/augmented_mesh.py` |
| **3** — `InverseMetricOperator` | ✅ LANDED `5def63b0` | `numerics/operator.py` |
| **4** — `LossKernelGauge` | ▢ not started | — |
| **5** — fire it + the warning + the record | ▢ not started | — |
| **6** — the two operator-algebra holes | ▢ not started | user ruled: fix both |
| **7** — promote the characterization suite | ▢ designed + measured, NOT executed | task #78 |

`[M]` `5def63b0`: **25 gates, 6-leg mutation battery, all arms bite, control
clean.** Neighbourhood 281 passed / 1 skipped, 22.78 s; `npx pyright` 0 errors
on the three production files. **Pure addition — nothing consumes either
surface yet, so no behaviour changed.**

### ▶ Where work resumes — Step 4, and what it may now assume

The predicate is callable and gated. `Π = R ∘ G⁻¹ ∘ M` is now *spellable*:
`frame.conjugate(InverseMetricOperator(frame.gram))`, gated by
`test_the_frame_can_now_spell_its_own_projector`. What remains is the
`LossKernelBasis` (the closed-form modes as a `Basis`), the caching, and the
firing.

⚠ **Existence-check before designing** (`plan-authoring` §1): `LossKernelGauge`,
`LossKernelBasis` — `[M]` **0 hits in `orpheus/`** as of `5def63b0`. Neither
exists; both are to be written. ✅ Re-verified 2026-08-14 at `839c36e9`.

⛔ **REFUTED 2026-08-14, before any edit — the MODEL half of this pointer was
false.** Original text, kept per §3: *"The basis construction exists twice in
throwaway form — `scratch/probe_344_kernel_basis.py` (563 lines) and, derived
independently, inside
`derivations/diagnostics/diag_344_reflective_box_loss_nullspace.py` (~50 lines,
the compact re-derivation). **Use the compact one as the model.**"*

`[M]` the diagnostic contains **no closed-form basis at all**. All 8 of its
null-space uses route through `_null_basis(A)` (`:122-126`) — a **dense
`np.linalg.svd` of an `A` assembled column-by-column** by `_dense()`
(`:100-109`, `n` matvecs). That is the exact route a closed form exists to
replace, and the memo prices it at **23.0 s vs 0.00 s** at `ndof = 3744`.

⟹ **The only closed form in the tree is the file this pointer forbids**:
`scratch/probe_344_kernel_basis.py::closed_form_kernel_basis` (`:306-358`,
helpers `:162-305`). It is **~145 lines**, not 50, and it *does* call
`np.linalg.svd` — but on the tiny per-orbit **coefficient** matrix, to reduce
an over-complete generating set (§2.3 of the memo: `2Σn` generators → `2Σn−1`
rank, exactly one relation). That is a legitimate and cheap use; it is not an
SVD of `A`.

▶ **The real model is the MEMO, not either file**: `scratch/issue_344_kernel_basis.md`
(714 lines, status COMPLETE) — §1.4 is the closed-form basis, §2 the
construction recipe in 5 steps, §6.2 the blocked representation that makes it
scale, §6.3 the two hazards. The memo is what the probe was written from.

⚠ Mechanism, for the surprise log: this is `plan-authoring` §1's **PRECEDENT**
clause — *"model this on X"*, where the symbol exists and every property
claimed of it is wrong — repeating **the same day the clause was written**, in
a different campaign. The §1 existence-check I did run (does `LossKernelBasis`
exist?) passed and is what made the pointer *feel* verified; the model claim
was never checked, and nothing in the pointer's shape distinguished the two.

### Three findings from the landed work that change nothing in the plan but would cost a re-derivation

1. `ρ(Σ) ≥ 1` needs **no subspace construction** — the absorption-damped mode
   is `< 1` and the blind sawtooth is exactly `1`. This is why it survives LD's
   moment tail and why `ndim = 1` needs no special case.
2. **`cell_kernel_batch` is implemented by BOTH shipped schemes** — its own
   docstring says only Diamond overrides it, which is stale. That staleness is
   what made the derivation possible at all.
3. The **`Q_cells` shape is `(1, 1)` for both** closures regardless of moment
   tail; the face shape is `(1,1,1,m)` with `m = basis**(ndim-1)`. Two of my
   probe iterations died on guessing otherwise, and both failures LOOKED like
   "the scheme is unclassifiable".

### ⚠ Three instrument bugs hit while gating this, all flattering

Recorded because each cost a cycle and each would recur: **(a)** ANSI codes
swallowed the pytest summary — verbatim `vv #17`; **(b)** a quoted `"$T"`
collapsed three test paths into one, so **zero tests ran** and the run reported
clean; **(c)** the `pairs_counts_faces` mutation left all six geometry rows
GREEN because `faces // 2` equals the pair count on every fixture written —
the gate was right for the wrong reason, and only a *second* mixed axis
separates the two rules. Always `--color=no`, always check the collected count,
always mutate before believing a gate.

---

## Context — why

`A = L + C − S − B` is **exactly singular** on any `d ≥ 2` Cartesian
**diamond-difference** mesh with ≥ 2 reflective axis pairs — and all-reflective is
the default for the eigenvalue entry, so this is the standard `k_inf` lattice.
`[M]` `dim ker A = 12` at d=2 LS4 ng=2; `138` at d=3 (3,4,5).

The consequence a user sees: **the returned boundary trace is a function of the
cold start, not of the problem.** `[M]` three cold starts differing only inside
`ker A` give identical `n_iter = 252`, identical residual `8.54e-14`, all
`converged=True` — while the trace differs by up to **27.3 %** and the bulk is
bit-stable at `7e-16`. Both convergence functionals are blind (residual
`7.3e-16`, balance projection `1.8e-16`, against a positive control that moves
them 15 and 14 orders).

**Disposition (settled, #344): DETERMINISM, not a correctness bug.** The
splitting is coherent (`M − N ≡ A` bit-exactly `0.0`, 20/20), `ψ_exact` solves the
discrete system at every mesh, and the trace error converges at exactly `O(h)`.

**Outcome wanted.** A user who solves an all-reflective box twice gets the same
boundary trace, and it is the physically correct one — with the solver saying so,
and saying what would remove the freedom at the root.

**Why it is fixable and not merely pinnable:** `[M]` `ψ_exact` is G-orthogonal to
`ker A` (`1.27e-15`), and this is now a **theorem** — every kernel mode carries a
non-trivial sign character on every axis, so any mirror-even functional
annihilates it. So the minimum-`‖·‖_G` member **is** the physical answer. The
gauge is canonical, not conventional.

---

## The two rulings that shape the design (user, 2026-08-14)

**R1 — "gauge" is the right word; it needs SPECIFICITY.** The package already
gauge-fixes the ψ½ System-B corner metric
(`numerics/spaces/radial_characteristic_space.py:119-144`). A gauge should be
called a gauge — qualified by what it fixes. ⟹ the family is **`LossKernelGauge`**.

**R2 — ⭐ DERIVE the applicability, never tabulate it.** Do not hard-code
"DD and all-reflective". Ask two questions and let the answer fall out:

```
    gauge_freedom  ⟺  scheme admits undamped zero-mean trace modes
                      ∧  the problem is reflective on ≥ 2 axis pairs
```

This is the Q6-E *ask-don't-tabulate* ruling applied again, and it is why a
scheme change is a real alternative rather than a coincidence. The scheme
property is computable from the closure's own face-to-face transmission
`Σ = (2/D)·1wᵀ − I`:

| scheme | spectrum of `Σ` on `{v : wᵀv = 0}` | gauge freedom |
|---|---|---|
| **diamond** | `−1`, multiplicity `d − 1` — **undamped** | **yes** |
| **step** | `0`, multiplicity `d − 1` — maximally damped | no |
| **linear-discontinuous** | `[M]` `dim ker A = 0` on the identical box | no |
| any scheme, d = 1 | no zero-mean face mode exists | no |

⟹ the predicate reads the scheme's own coefficients. A future scheme answers
for itself with no edit here.

**R2b — the freedom must be AUDIBLE.** When gauge freedom is detected the solver
warns: the fixing was applied, AND changing to a scheme without the undamped
zero-mean trace mode removes the freedom **at the root**. Follows the #340
`ConvergenceWarning` family (`numerics/convergence.py`), which is the precedent
for "say what it cost and what to do".

---

## Measured findings that de-risk the build

- ⛔ **REFUTED 2026-08-14 — the DIAGONAL claim was a d=2 reading promoted to a
  universal, and acting on it would have shipped a silently WRONG projector.**
  Original text, kept per §3: *"⭐ **The Gram is EXACTLY diagonal.** `[M]`
  `max|offdiag|/max|diag| = 0.000e+00`, `cond = 1.60` (LS4) / `2.10`
  (lebedev(11)). The closed-form modes are already G-orthogonal (disjoint
  supports). ⟹ `GramStructure.DIAGONAL`, which `FrameBase.gram` accepts —
  **#275 (the dense-Gram seam) does NOT block this** — and **no QR is needed**,
  killing the second-largest cost (44 ms / 748 ms)."*

  `[M]` `$CLAUDE_JOB_DIR/tmp/probe_gram_diagonality.py`, on the raw pair-generator
  basis, `max|offdiag| / max|diag|` of `BᵀGB`:

  | configuration | per-orbit rank | WHOLE | worst BLOCK |
  |---|---|---|---|
  | d=2 (3,4) LS4 ng=2 | **1** | `0.000000e+00` | `0.000000e+00` |
  | d=3 (2,2,2) LS4 ng=1 | 11 | `1.723861e-01` | `1.914395e-01` |
  | d=3 (3,4,5) LS4 ng=2 | 23 | `4.045273e-01` | `4.312919e-01` |

  ⭐ **Why the original reading was vacuous:** at d=2 `κ({x,y}) = 1`, so each
  orbit carries exactly **ONE** mode — a 1×1 Gram is diagonal for free, and
  `cond = 1.60` is the spread *across* blocks, not evidence of orthogonality
  *within* one. At d=3 an orbit carries `2Σn − 1` modes and they are **not**
  G-orthogonal: the `{a,b}` and `{a,c}` pair generators both live on the `a`
  faces. `plan-authoring` §2's QUANTIFIER clause — the denominator was
  "the only shipped rank-1 case".

  ⚠ **The consequence is not a bad constant, it is a wrong answer.** `frame.gram`
  computes its diagonal by the row-sum probe `analysis(reconstruction(ones))`,
  which equals `Σ_j (MR)_{kj}` — the row sum, correct *only* if `MR` is
  diagonal. Declaring `DIAGONAL` on a 43 %-off-diagonal Gram normalises every
  coefficient by the wrong number, and nothing raises. This is precisely the
  refusal `GramStructure`'s own docstring exists to enforce.

  ✅ **FIX, measured — one `sqrt(G)`-weighted SVD per block, which REPLACES both
  of the memo's factorisations rather than adding a third.** The memo runs an
  SVD in the coefficient space (to cut the over-complete generators) *and* a
  QR in the state space (to G-orthonormalise). Fusing them: `U, s, _ =
  svd(√G · B_gen)`; keep `rank = #{s > 1e-10·s₀}`; `Φ = U[:, :rank] / √G`. Then
  `ΦᵀGΦ = I` — so `DIAGONAL` becomes **true by construction**, `G⁻¹` is the
  identity, and the rank the SVD finds IS the counting law (a free gate).
  `[M]` `$CLAUDE_JOB_DIR/tmp/probe_fused_orthonormal.py`:

  | configuration | rank == dim_R | `‖ΦᵀGΦ − I‖∞` | span vs `B_R` |
  |---|---|---|---|
  | d=2 (3,4) LS4 ng=2 | 12 ✅ | `4.441e-16` | `8.063e-16` |
  | d=3 (2,2,2) LS4 ng=1 | 33 ✅ | `2.220e-15` | `3.188e-15` |
  | d=3 (3,4,5) LS4 ng=2 | 138 ✅ | `1.132e-14` | `1.580e-14` |
  | d=2 (3,4) lebedev(11) ng=2 | 18 ✅ (T = 224) | `8.882e-16` | `1.182e-15` |

  ⟹ **`LossKernelBasis` IS the G-orthonormal basis, not the pair generators.**
  The generators are the *construction*; the orthonormal frame is the *object*.
  The `√G` inverse is well-posed only because R has no support on `G == 0` rows
  (the next bullet) — `[M]` re-verified as an explicit precondition on all four
  rows above, including the lebedev row where `T = 224`.
- ⭐ **`R` and `T` are separated by the metric's zero-set.** `[M]`
  `max|B_R|` on `G == 0` rows is **`0.000000e+00`** on LS4, product(4,4) AND
  lebedev(11). So `T = ker A ∩ ker G` and `R ⊆ (ker G)^⊥`. The gauge is
  well-posed on `R` at **every** quadrature, not just the tangential-free ones.
  `T` is genuinely ungaugeable (no metric to minimise) and stays out of scope.
- **The basis never reads a cross-section** (`[M]` fissile vs absorber
  bit-identical `2.799e-16`) ⟹ Stratum-1 cacheable, per
  `orpheus/sn/sweep/cache.py:121` (`GeometryCoefficients`).
- ⛔ **REFUTED 2026-08-14 — this is the source of the false model pointer
  above, and the claim itself does not survive.** Original text, kept per §3:
  *"The `test-architect` **independently re-derived the basis in ~50 lines**
  and reproduced the numbers bit-for-bit — a second, structurally independent
  construction exists."*

  `[M]` **no such construction is in the tree.** Searched every `.py` under
  `scratch/`, `derivations/`, `tests/`, `orpheus/` and every `.md` under
  `.claude/agent-memory/`: exactly ONE function builds a kernel basis without
  an SVD of `A` — `scratch/probe_344_kernel_basis.py::closed_form_kernel_basis`.
  What the diagnostic agent *did* produce independently is the dense-SVD
  **ground truth** the closed form was verified AGAINST (memo §3.1/§3.2: 13
  configurations, subspace gap `≤ 2.2e-14`), and the counting-law evaluator
  `predicted_dim` (`probe:370`) which never builds a vector. Those ARE two
  genuine independent checks — one numerical, one combinatorial — so the
  *verification* claim stands in a stronger form than written. The
  **construction** claim does not: there is one construction, not two.

  ⚠ The consequence for Step 4 is a real one, not bookkeeping: I believed a
  second implementation existed to cross-check the port against. It does not,
  so the port's oracle has to be **the dense SVD on the small fixtures** (which
  is exactly what `diag_344`'s `_null_basis` provides, at d=2/small-d=3 sizes)
  plus **the counting law** — never "the other implementation".

---

## ⛔ Premise corrections — three of mine, measured after I wrote them

1. **`CAP_APPLY`/`CAP_SOLVE`/`MissingCapability` DO NOT EXIST** — retired at #226.
   The tree uses three predicate/Protocol/TypeGuard axes
   (`is_invertible`/`SupportsInverse`/`invertible()`/`NotInvertible`, and the
   adjoint/assembly siblings), `operator.py:688-1211`. **"Apply but not solve" is
   spelled by ABSENCE of the method** (`TraceRestrictionOperator:2773` is the
   precedent); a raising stub is the named anti-pattern.
2. **Spaces CANNOT decide a sum's invertibility.** `[M]` proof by example:
   `I + (−1)·I` **is the zero map** and reports `is_invertible = True`, with
   identical spaces to `I + I`. The existing contract
   (`OperatorSum.is_invertible` returns the left-spine head, designating a
   `GreenOperator` preconditioner) is deliberate and correct for its purpose;
   invertibility is spectral, not structural. ⟹ Hole 2's fix is **a typed
   complement**, not a change to the global predicate.
3. **The cost figures I first relayed do not reproduce** (apply `0.094 ms` vs
   `398.6 µs`; `(8,8,8)` `0.703 ms` vs `10.95 ms`). ⚠ INCONCLUSIVE — three agents
   were competing for CPU. **Re-measure single-process before quoting any cost.**

---

## ⚠ The one unverified premise — measure it FIRST

**`solve_sn`'s returned trace is not the converged iterate's trace.** The
eigenvalue outer iterates on the SCALAR flux; the angular composite is a
side-stash (`solver._psi_typed`). The returned `angular_flux`/`boundary_flux`
come from a **final reconstruction solve** at `orpheus/sn/solver.py:2489`
(`final_implicit.solve(...)`) with the reflected converged trace as prescribed
inflow. **Whether that sweep preserves, damps, or re-injects the kernel component
is unmeasured** — my probe timed out before reporting.

### ✅ Step 0 — ANSWERED. The sweep PRESERVES it. Gauge fires at the eigenvalue entry too.

`[M]` `scratch/probe_344_step0_finalize_sweep.py`, d=2 all-reflective LS4 2g
**fissile**, both meshes, with both controls firing:

```
(3,3)  jacobi        keff=1.8750000000  n_outer=3  trace spread=4.190021e-11
       gauss_seidel  keff=1.8750000000  n_outer=3  trace spread=3.464537e-02
       ||P_ker (GS-J)|| / ||GS-J||  = 1.000000     <-- ENTIRELY a null vector
(3,4)  ||GS-J||/||J|| = 6.091755e-02, ||P_ker||/||·|| = 1.000000
CONTROL true null vector = 1.000000 (must be 1)   off-kernel = 4.16e-17 (must be 0)
```

Jacobi lands on the correct near-uniform member (spread `~2e-11`); G-S carries a
frozen `~3.3e-02` kernel component. `keff` is identical and equals the analytic
`k_inf = 1.875` — the mirror-even blindness, confirmed end-to-end.

⚠ **Instrument note, cost one relaunch:** the first run used
`get_mixture("B","2g")`, which is **not fissile** (`[M]` `SigF.sum = 0.0`) — an
eigenvalue entry on it cannot converge. `[M]` of the 2g library **only "A"
fissions**. With the right material the whole probe runs in **~7 s**
(`n_outer = 3`), so the "the eigenvalue path is slow, budget it" worry in the
original text was an artefact of the wrong mixture, not of the entry.

---

## Design

### ⛔ Step 1 — WITHDRAWN 2026-08-14, before any edit. Its premise is refuted.

*Original text:* "one `Solution` construction site — `_package_solution`
(`solver.py:2615`) documents itself as *'The ONE `SolutionBase` construction
convention'*; `[M]` present-tense FALSE, 2 of 3 forward routes bypass it at
`:3804`/`:3985`. Route both through the helper so the gauge is **one** call
site, not three."

The docstring IS false and that stays worth fixing — but routing the two arms
through the helper would be a **silent behaviour change outside this scope**,
and the "one call site" requirement it served is not a real requirement.

1. ⛔ **The bypass encodes a DELIBERATE, documented difference, not drift.**
   `solve_sn` packages `_cell_average_angular(...)` — the slot-0 cell-average
   view. The two fixed-source arms return `angular_out` **whole**, the full
   moment-tailed interior, and `solver.py:3768-3773` says why: *"a multi-moment
   closure's φ̂ slopes are internal DG structure, not the scalar flux the
   Solution reports (#240 D5b-S3)."* Unifying would discard DG slope structure
   from every fixed-source return.
   (The other split — `_history=()` vs the driver's history — IS latent: `[M]`
   nothing outside `timed_full_field.py` reads `_history` off a returned
   `Solution`. Only the moment-view split is load-bearing.)
2. ⛔ **"N call sites = twin path" is false here, and the tree says so.** Pattern 2
   forbids two IMPLEMENTATIONS of one computation, not N invocations of one
   function. `[M]` the house pattern is exactly N: `_exit_balance_defect` **4**
   call sites, `_certify_within_group_exit` **4**, `_package_solution` **2**.

⟹ **The gauge follows the `_exit_balance_defect` pattern**: one named function
plus one shared named predicate (Step 2), invoked from each entry that needs it.
⚠ With one sharpening, because a gauge MUTATES where those hooks only REPORT —
forgetting a site returns a silently ungauged answer. So the coverage is gated:
a test that enumerates the public entries and asserts each one gauges. That
replaces the structural guarantee unification would have given.
▸ Fix the false `_package_solution` docstring in passing; do not act on it.

### Step 2 — the scheme property, asked not tabulated
A derived predicate on the spatial scheme: *does the closure admit an undamped
zero-weighted-mean face mode?* Computed from the scheme's own transmission
coefficients (the `Σ` table above), beside the existing
`_residual_is_expressible` (`solver.py:396`), which is the tree's precedent for
a **shared, named, single-sourced applicability predicate** (its docstring at
`:407-412` says exactly why it exists). Plus `SNMesh.reflective_axis_pairs` —
`[M]` no such predicate exists today; the law is reachable per-face as
`sn_mesh.bc[face].law`, and it belongs on `SNMesh` beside
`radial_characteristic_field_space`, not inlined at call sites.

### Step 3 — the missing numerics primitive (ONE new thing)
`[M]` The projection algebra already exists, factored as frame theory —
`Π = R ∘ G⁻¹ ∘ M`, written verbatim at
`orpheus/numerics/basis/indicator_basis.py:67-69` and `numerics/frame.py:75`.
The **only** gap: `G⁻¹` lives as a *space metric* (`FrameBase.gram` returns a
`FunctionSpace`), not as an operator, so `frame.conjugate(G⁻¹)` cannot be
written, and `FrameBase.project` (`frame.py:310`) stops at coefficients.
⟹ add a `LinearOperator` view of a `FunctionSpace`'s (inverse) metric. Zero twin
exists (`grep inner_product_weights orpheus/numerics/operator.py` → nothing).
Reuse `apply_inverse_metric`'s existing Moore–Penrose masking
(`space.py:315-322`) rather than re-deriving it.

### Step 4 — `LossKernelGauge`
- `LossKernelBasis` — the closed-form modes as a `Basis`
  (`numerics/basis/base.py:117`), `gram_structure = DIAGONAL` (measured).
- The projector via the frame route, which `[M]` **avoids Hole 1** — both frame
  faces declare a real `apply_transpose`, so `(R @ M).H` builds, whereas
  `A.H @ A` raises `MissingAdjoint`.
- **Cached on the phase-space object that owns the geometry** (the trace space /
  `SNMesh`), following `AngularTraceSpace._face_restrictions`
  (`angular_trace_space.py:494-538`) — *not* on the solver, *not* on the operator.
- Scope: **`R` only.** `T` (tangential, `G ≡ 0`) has no minimum-norm
  representative — that is the separate "type it away" issue.

### Step 5 — fire it, and say so
One call from the single construction site of Step 1, guarded by Step 2's
predicate. Record it on `IterationHistory` following the `balance_defect`
precedent (`solution.py:197`) with the same `None`-means-*not-measured*
discipline. Emit R2b's warning.

### Step 6 — the two operator-algebra holes (user ruling: fix both)
- **Hole 1** — `_AdjointOperator.apply_transpose` raises `NotImplementedError`
  (`operator.py:1312-1317`) *"until a consumer demands it"*, and
  `_AdjointOperator` does not override `is_adjointable`. A self-adjoint projector
  is that consumer. ⚠ Shared numerics: needs its own gates and a wide re-run.
- **Hole 2** — per correction 2 above, the fix is a **typed complement** that
  declares `is_invertible = False` (Pattern 4), NOT a change to
  `OperatorSum.is_invertible`. File the general cancelling-sum defect separately
  — `I + (−1)·I` claiming invertible is real and predates this work.

### Step 7 — land the promotion (task #78, designed and measured)
`derivations/diagnostics/diag_344_reflective_box_loss_nullspace.py` →
`tests/sn/operators/test_loss_nullspace_reflective_box.py`. ⭐ `[M]`
`pyproject.toml` sets `testpaths = ["tests"]`, so the diagnostic **currently
contributes zero executed coverage** — this is a net gain, not a relocation.
Trimmed: **40 tests ≈ 53 s** projected, slowest single test 9.1 s, **not**
`slow`. ERR-056 control re-verified (`1.0000e+00` trace / `4.7891e-01` bulk vs
`1.1414e-12` baseline). ⚠ `docs/theory/verification/matrix.rst` is a committed
generated artifact counting foundation tests — regenerate and stage it.

---

## Verification

- **Step 0 gates everything.** Do not choose a placement before it reports.
- **The flagship**: on the #344 fixture, two cold starts differing only inside
  `ker A` return traces agreeing to solver tolerance, and the returned trace
  matches the analytic uniform answer. `[M]` today: 27.3 % apart.
- **The negative control that stops it being a universal absorber** (`vv` #19):
  Jacobi already lands on the correct member, so its deviation is pure round-off
  and must be measured **out of span** — `[M]` `1.0000e+00`. The gauge must
  remove *nothing* there.
- **Residual-neutrality**: `A(ψ − Πψ) = Aψ` by construction ⟹ no convergence
  certificate may move. Assert it, do not assume it.
- **Bulk bit-identity**: `[M]` the kernel is pure-trace (bulk share `1.1e-28`),
  so `scalar_flux`, `k_eff` and every bulk snapshot must be **unchanged**. This
  is the cheap regression check and it covers the DD snapshot set, which stores
  `scalar_flux`.
- ### ✅ Blast radius — ENUMERATED. **Zero artefacts, zero tests to re-baseline.**
  `[M]` a swap-and-run census over the wide gate (**8501 passed, 1:01:37, zero
  red**) with counters on `SNMesh.__init__` *and* all three `Solution`
  construction paths: **298 tests** build the singular configuration (46 files),
  **24** reach a `Solution` exit on it, and **exactly 1** asserts on a trace —
  `test_bc_extraction_2d.py:404` — which is **structurally invariant**, because
  the two operands it differences are the inflow-slot component of `Aψ`, and
  `A(ψ−Πψ) = Aψ`. Every frozen artefact in all 8 directories is classified in
  the agent's table: the look-alikes (`2d_octant_equivalence_02`'s face
  payloads, `walk_matvec_cart2d_2g`'s `fwd_trace`/`adj_trace`,
  `wave_t_t4`'s `*_specular_apply_boundary`) all carry a real trace on a real
  singular fixture but are **operator/sweep-level, never a solve exit**.
  `sha256` pins: none survive as frozen references (#333 ended that era).
- ⛔ **…and that is the finding: the change ships ENTIRELY UNGATED.** `[M]` zero
  production consumers of `boundary_flux`; `Solution.compare()` diffs
  `interior.values` only. Per `plan-authoring` §6c the change owes a **witness**,
  and the enumeration hands us its design: the `|Ω·n|^0` full-face moment is the
  *only* functional that moves, and only where the quadrature has **tangential**
  ordinates — and `[M]` **no shipped all-reflective 2-D fixture uses one** (all
  are `level_symmetric`). ⟹ the gate is an all-reflective 2-D DD solve under
  **`lebedev`**, asserting the full-face `|Ω·n|^0` moment against a pre-gauge
  literal, **paired with a `level_symmetric` control that must NOT move**.
- ⚠ **Decide before gauging the Krylov arm:** its returned `boundary` is *"the
  matvec's B1'' face residual"* by the code's own comment, not a flux trace.
  Projecting a **defect** off `ker A` is a different operation, and by
  residual-neutrality it is a **no-op**. Either it is the wrong thing to return,
  or it must be exempted — say which, in the code.
- ⭐ **All-reflective is the default more broadly than assumed:** `solve_sn` and
  `solve_sn_adjoint` have **no `boundary_condition` parameter at all**, and any
  bare `SNMesh(mesh, quad, mats)` resolves to all-reflective
  (`transport/method.py:257`). That is why the exposure is 298 tests, not ~20.
  ⚠ And `_apply_default_bcs` fills only when **all** faces are `None`, so a
  *partial* declaration silently keeps reflective on the rest.
- Two pre-existing defects found in passing, both worth a separate fix: the
  false `_package_solution` "ONE construction convention" docstring, and three
  non-existent `cyl_*` snapshot filenames at
  `docs/theory/methods/sn/curvilinear_numerics.rst:1810-1812` (renamed to
  `folded_4x8`/`folded_2x4` at Q5.6.3).
- **Cost, single-process, no competing agents** — re-measure build/apply and
  quote them with the configuration.
- Wide gate `python -O -m pytest -m "not slow"`, SERIAL, **budget ≥ 90 min**,
  detached via `Popen(start_new_session=True)` + a persistent `Monitor` writing
  to a LOG. Sphinx `-W` clean; `npx pyright` clean.

## Out of scope (file, do not build)
Typing away the tangential component `T`; the general cancelling-sum
`is_invertible` defect; #343's octant reorder; the ψ½ gauge's own naming.
