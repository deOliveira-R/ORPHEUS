# B3.0–B3.2 doc repair — `docs/theory/**` audit

Ground truth read this session (NOT the brief, NOT docstrings alone):

* `9e2139b4` (B3.0), `b39502f8` (B3.1), `7f02de15` (B3.2) — commit messages + diffs
* `.claude/plans/b3_domain_narrowing_crosswalk.md` §7 / §9 / §9.1 / §9.2
* `orpheus/geometry/boundary/_factors.py`, `.../reflective.py`, `.../periodic.py` (live)
* `orpheus/sn/boundary/realizer.py`, `orpheus/sn/boundary/angular.py` (live)
* `orpheus/sn/operators/boundary.py::_reflect_trace` (live) + pre-B3.2 via `git show 7f02de15^`
* `orpheus/numerics/operator.py::TraceRestrictionOperator` / `IncomingOrdinateMaskTensor`
* `orpheus/numerics/spaces/angular_trace_space.py` (`outflow_restriction` / `inflow_restriction`)
* `orpheus/diffusion/boundary_realizer.py` (whole file — the tier-2 / refusal-axis quotes)
* `tests/geometry/test_bc_equivalence_snapshot.py`, `tests/sn/sweep/core/test_phase_c_gates.py`
* `.claude/skills/vv-principles/error_catalog.md` ERR-042 / ERR-044 / ERR-045
* MEASURED: `law.geometry_map` / `law.response_kernel` on all seven laws;
  `dataclasses.fields(PeriodicBoundary) == ['axis']`; `TANGENTIAL_EPS = 8.88e-16`;
  `ZeroOperator.is_adjointable/is_invertible`; exports of the five factor types

Sphinx `-E -W`: **pre-edit baseline EXIT 0, 0 WARNING/ERROR/CRITICAL** → **post-edit EXIT 0, 0**.

Scope honoured: **no file outside `docs/theory/` touched.** No `git checkout/restore/stash/clean`
run on any path. No commit.

---

## 1. Verdict table — claims checked

Verdicts: **FALSE** = present-tense claim contradicted by live code · **STALE** = true history
presented as present tense · **TRUE** = verified, left alone.

### `docs/theory/foundations/boundary_conditions.rst`

| # | Claim (pre-edit) | Verdict | Action |
|---|---|---|---|
| 1 | Key Facts: "`G` is a geometric map (permutation, pushforward, **angular average**)" | **FALSE** (B3.0 — an average fails multiplicativity, so it is `R`) | Rewrote the bullet: `G : Γ₊→Γ₋` deck transformation carrying the crossing, `R : Γ₋→Γ₋` constitutive kernel, multiplicativity criterion stated, B3.0 named |
| 2 | Key Facts: "`VacuumInflow` realises to `IncomingOrdinateMaskTensor`, which zeroes only the inflow ordinates and preserves the outflow trace … the **sole** path to vacuum action" | **FALSE** (B3.2 — zero map `Γ₊→Γ₋`) | Replaced with two bullets: the Γ₊-domain narrowing + `B_face = ι₋∘law∘γ₊`, and vacuum-as-zero-map with the mask's fate; added a THIRD bullet on the three-way partition + `|Γ₊|=|Γ₋|` shape blindness |
| 3 | Layer-1: "the signed-projection table is what `SNBoundaryRealizer` reads to build the sparse **vacuum-mask** operator" | **FALSE** | Rewrote: the realizer reads BOTH half-traces — codomain `Γ₋` (ERR-041-checked) and domain `Γ₊` |
| 4 | Layer-3: "SN realizes vacuum as a sparse per-ordinate `IncomingOrdinateMaskTensor` on the inflow ordinates" | **FALSE** | Replaced with the zero map + both space hooks; added a diffusion row (`𝒜 = 0`, Marshak) so the "structurally different operators" list shows the two *shipped* methods |
| 5 | Design correction 2: "The realized per-face law **is** a full-face operator … the fix: `_apply_faces` **projects** the emission onto the codomain row — `apply` writes `inflow_indices_for_face`; `apply_transpose` writes `outflow_indices_for_face`" | **FALSE** as a contract. Pre-B3.2 the transpose wrote the **whole face** (`face_view[...] = law.apply_transpose(masked)`), and `boundary.py`'s own ⚠ names output-projecting `lawᵀ` onto the *law's codomain* as the WRONG spelling. So the doc documented neither the old code nor the right rule. | Past-tensed the diagnosis (kept — it is why the residual must carry no `B` term), replaced the remedy with B3.2's typing, and added a `.. warning::` stating the **surviving** trap precisely: the transpose scatters over `Γ₊`, never `Γ₋`; the wrong spelling extracts the diagonal block (spurious `+1` for vacuum), caught only by A2a grid-reciprocity on the het-VACUUM sphere because permutation laws are bit-identical under either spelling |
| 6 | `reflect_into_inflow` §: "The projection onto the inflow row is load-bearing: the realized law is a *full-face* operator … `_reflect_trace` projects the forward action onto `inflow_indices_for_face` and the Euclidean transpose onto `outflow_indices_for_face`" | **FALSE** | Rewrote to typing-not-projection; stated the transpose's input-mask / output-scatter asymmetry and cross-linked the warning |
| 7 | Same §: per-face action is "apply the law … and project onto the inflow row", "zero for vacuum" | **FALSE** | Rewrote as restrict → apply → scatter; specular is a permutation on the **reduced** axis |
| 8 | `bc-law-layer`: "`geometry_map` and `response_kernel` return the ABC's `None` on **every** law and are read by nothing" | **FALSE** (falsified by **B1** `a0fd17b4`, never doc-fixed; entangled with B3.0) | Replaced with a **measured** 7-row table of the live `(G, R)` per law, plus the B2.2 production reads and why the ABC keeps `None` defaults |
| 9 | Law table: periodic `G` = "spatial pushforward (**caller-supplied**)" | **FALSE** (B1 — the law carries `axis`; the realizer derives the partner) | Corrected |
| 10 | Descriptor-model bullet: "only `source` is ever overridden; `geometry_map` and `response_kernel` inherit the ABC's `None` on every law" | **FALSE** | Corrected + pointer to the measured table |
| 11 | B2.0 §: "collapsing diffusion's five-arm ladder — `bc[face].law.response_kernel.**scalar**`" | **FALSE** (B3.0 renamed `.scalar` → `.amplitude`) | Corrected, with the reason ("scalar" is actively wrong once a rank-one kernel exists) |
| 12 | `bc-trace-structure`: "The SN realizer's vacuum branch … indices are the constructor argument to `IncomingOrdinateMaskTensor`" | **FALSE** | Rewrote to the both-masks / domain+codomain reading incl. the `outflow_indices` field |
| 13 | `bc-tensor-primitives` table: vacuum = `IncomingOrdinateMaskTensor`; reflective = "bare `PermutationOperator` with `perm = reflection_index(axis)`" | **FALSE** (both) | Vacuum → two-hook zero map with the "coincidence not contract" note; reflective → TP on the reduced axis with `local_perm = γ₊.to_local(...)` and why `to_local` is mandatory; white/albedo/periodic marked **not yet narrowed**, prescribed marked *raises*; added a note that a shape assertion cannot tell the typings apart |
| 14 | Fast-path §: bit-identity rationale stops at the Wave-6 snapshot | **STALE** | Added the B3.2 bit-identity gate: reference materialised off the law *descriptor*, `np.array_equal`, 5 fixtures, falsified by the naive-`arange` mutation (6 rows red, 26 interceptions positive control) |
| 15 | `SNMethodSpace` attribute list: no `outflow_indices`; `minimal` described as a general fixture | **STALE** (B3.2) | Added the `outflow_indices` bullet (domain, sibling-of-codomain rationale, the realizer's error) and the "`minimal` is now a partial constructor" note |
| 16 | Worked example header: `… → IncomingOrdinateMaskTensor` chain | **FALSE** | → `ZeroOperator(Γ₊ → Γ₋)` |
| 17 | Worked example Step 3: "`for_face` derives inflow_indices" | **STALE** | Now derives BOTH; comment states neither is the other's complement |
| 18 | Worked example Step 4: realizer code block returning `IncomingOrdinateMaskTensor`; "self-adjoint, idempotent projector" | **FALSE** | Replaced with the live B3.2 arm (γ₊ + two `_zero_rows` hooks), plus a note on what it replaced and why the hooks are load-bearing (`0.0*x` echo would be right only by the `|Γ₊|=|Γ₋|` accident) |
| 19 | Worked example Step 6: "Shim forwards to `IncomingOrdinateMaskTensor.apply` … outflow rows pass through unchanged" | **FALSE** twice over (sweep no longer calls it; law no longer full-face) | Replaced with a then/now code contrast + a `.. warning::` recording the **measured** B3.2 gap: the narrowed law does NOT validate its own domain — unreachable via `_reflect_trace`, reachable via `sn.bc[face].apply`, 4 strict xfails, closes at B3.4 (one place, not seven copies) |
| 20 | `bc-vacuum-semantic-correction`: "The post-Wave-8 realizer's vacuum branch **returns** `IncomingOrdinateMaskTensor` … the right algebraic object" | **FALSE** | Section retitled "… and its dissolution at B3.2" + `.. important::` status banner; body past-tensed. Section KEPT (genuine history) |
| 21 | "Why this matters" — 3 consequences, incl. *"a future adjoint sensitivity path needs the outflow trace preserved"* | **FALSE** (#1 — consumer **measured not to exist**); **moot** (#2, #3) | `.. note:: **Retraction (2026-07-31, B3.2)**` + each argument preserved and past-tensed with an explicit **Disposition**: #1 = the declared-capability-no-consumer pattern; #2 = *inverted* (a non-endomorphism cannot be idempotent; `ZeroOperator` is now the right tag); #3 = right observation, mechanism one layer too shallow |
| 22 | "…quadrature adapters carry every tangential ordinate at μ = 0 so ψ = 0 on `T_f` for a properly-initialised flux" | **FALSE as stated** (μ≈0 is the *definition* of tangential, it does not make ψ vanish) | `.. warning::` with two measured corrections: (a) the honest weaker statement (no operator writes tangential slots, so a zero-initialised carrier keeps them zero — an initialisation property; the zero *metric weight* is the separate structural fact); (b) `outflow_indices_for_face` is strict `O_f` while the mask preserved `O_f ∪ T_f`, so `I − P_in ≠ P_out`; cyl `product(2,4)` = 4/8 tangential at `xmax` |
| 23 | Wave-8 call-site audit with `sweep.py:334` … line refs | **STALE** (files dissolved; Wave O removed `bc.apply` from the sweep) | Past-tensed + explicit "these file:line refs are frozen at that audit and do not resolve against the live tree"; added the second reading (no consumer ⇒ B3.2 acted on it) |
| 24 | "The curvilinear sweeps **now consume** the realizer-routed `IncomingOrdinateMaskTensor`" | **FALSE** | Past-tensed the C188.3 evidence; kept the uniformity claim and noted it survives B3.2 |
| 25 | "the `vacuum_lebedev17` case compares **inflow-row outputs only**" | **FALSE** (re-posed at B3.2) | Rewrote to the live assertion (index-set cross-check + `Γ₋` shape + all-zero) |
| 26 | V&V-status block: "intentional semantic capture of the §16A.5 inflow-only mask for vacuum" | **FALSE** | Corrected to the narrowed zero map |
| 27 | Option-(a) status: "the only path to vacuum action is … producing `IncomingOrdinateMaskTensor` output … the §16A.5 inflow-only-mask body is the **unique** vacuum semantics" | **FALSE** | Uniqueness claim kept (it survived); the operator behind it re-pointed to B3.2 |
| 28 | Option-A cost 2: "a documentation-burden landmine for future adjoint-sensitivity consumers that read outflow rows" | **FALSE** (same phantom consumer) | Kept as history, marked measured-not-to-exist, and noted Option A's retirement stands on costs 1 and 3 |
| 29 | β2 §: "realize-then-apply goes through the inflow-only-mask path … `vacuum_lebedev17` still pins inflow rows only" | **FALSE** | Corrected; the *uniqueness-by-design* argument preserved |
| 30 | Wave-6 snapshot table: vacuum tolerance "inflow-rows-only"; mixed row silent about its xfail | **FALSE / STALE** | Vacuum row → zero map; mixed row → **xfail(strict) pending B3.4**; added a note distinguishing the 3 re-posed-onto-`Γ₊` rows from the 4 still-full-face rows, and framing the Marshak xfail as an honest red that flips on B3.4 (never a suppression) |
| 31 | `bc-two-bc-applies-per-matvec` (whole §), incl. two present-tense `IncomingOrdinateMaskTensor` sentences | **STALE** (pre-dates B3 — Wave O made the matvec call `bc.apply` **zero** times; verified against the live `test_phase_c_gates.py` docstring) | Added an `.. important::` historical banner naming what changed and why the § is kept; past-tensed the two mask sentences. Deep rewrite of the § deliberately **not** done (out of B3 scope) — flagged below |
| 32 | `ordinate-partition-inflow-outflow` (3-way partition eq) | **TRUE** — and it is a `verifies` target (2 tests in `tests/numerics/test_angular_trace_space.py`) | Untouched (L-003) |
| 33 | `bc-universal-invariants` table (the B0 "*Intended* / *Implemented*" rows, ERR-042 row) | **TRUE** | Untouched; my new equivariance § is consistent with it |

### Other files

| # | File / claim | Verdict | Action |
|---|---|---|---|
| 34 | `methods/sn/boundary_conditions.rst` — resolution table: vacuum = `IncomingOrdinateMaskTensor`, reflective = bare `PermutationOperator(reflection_index)` | **FALSE** | Both rewritten (zero map / reduced-axis TP); white / periodic / albedo marked *not yet narrowed*; prescribed marked *raises*; added the B3.2 note incl. the shape-blindness caveat |
| 35 | same — `R = c₁ G_refl + c₂ G_diff` (Marshak) and "a linear operator **R** mapping outgoing → incoming" | **FALSE** (B3.0: two meanings of `R`; and `G_diff` = the Lambertian is a **response**, not a geometry) | Prose `R` → `B`; added a `.. warning::` typing the two decompositions (rank-N **sum** `Σ G_α ⊗ A_α` vs the affine **factorisation** `R G`, composite named `R ∘ G`, never `R`); Marshak line → `B = c₁B_refl + c₂B_diff` with what actually distinguishes the terms |
| 36 | same — primitives table: white `G_α` = "cosine-weighted hemispheric average"; vacuum `G_α` = `0`; periodic "caller-supplied" | **FALSE** (all three, B3.0/B1) | Added an *affine factors (G, R)* column with the live values; white row flagged response-not-geometry; trailing para on why vacuum's `G` is the identity deck element, not zero |
| 37 | same — 4 × `:meth:`SNMesh._resolve_one`` | **FALSE** (dead symbol; renamed `realize_boundary_law` at #290 P7b) | Repointed all 4 (grep-gated — `-W` is blind to these) |
| 38 | `foundations/operator_tensor_network.rst` — "`IncomingOrdinateMaskTensor(...) & IdentityOperator()` for vacuum" | **FALSE** | Rewrote the row: specular TP on the **reduced** axis, vacuum drops out of the TP form entirely (bare two-hook `ZeroOperator`) |
| 39 | `foundations/operator_algebra.rst` — "the SN realizer's vacuum branch consumes `inflow_indices_for_face` to construct an `IncomingOrdinateMaskTensor`" | **FALSE** | Rewrote to both-selectors / domain+codomain |
| 40 | same — "tangential band … (default `ε = 1e-12`)" | **FALSE** (measured `TANGENTIAL_EPS = 4·eps ≈ 8.9e-16`) | Corrected + added "so 'not inflow' is never 'outflow'" |
| 41 | `conventions/indexing_and_layout.rst` — `BoundaryTraceLaw` role given as `ψ⁺_Γ = T(ψ⁻_Γ)` | **FALSE** — direction reversed against the page's own `Γ₋`=inflow / `Γ₊`=outflow convention | Replaced with the affine form + "outflow trace in, inflow trace out" |
| 42 | same — `InflowTraceSpace / OutflowTraceSpace` row; `IncomingOrdinateMaskTensor` row | **FALSE / STALE** (#205/#201 collapse; B3.2 deprecation) | → `AngularTraceSpace` (one space, two selectors, three-way partition) and → `TraceRestrictionOperator` (sibling-not-subclass); realizer row gained the B3.2 typing |
| 43 | `methods/sn/curvilinear_one_group.rst` — "the operator … internally consumes the outflow selector … and writes only the inflow slots; **the outflow slots in the output are unspecified**" | **FALSE** (twice: it does not consume the selector — the *consumer* restricts; and there are no outflow slots in the output) | Added a note marking the code sketch as the Wave-8-era shape with the two changes named; rewrote the prose; kept the "ghost cell" idiom, which survives |
| 44 | `foundations/operator_adjoint.rst` — tangential slots never sourced; `B` carries the composite `full_field_space` | **TRUE** (verified: `SNBoundaryOperator.domain/codomain` still return `sn_mesh.full_field_space`; the narrowing is per-face, inside `_reflect_trace`) | Untouched |
| 45 | `verification/sn.rst` — `AngularBoundarySourceSink.prescribed_inflow` writes only inflow slots | **TRUE** (unrelated to the law's typing) | Untouched |
| 46 | Retired-name grep: `HemisphericalAverage`, `NullMap` | **absent from `docs/theory/**`** | No action needed |

---

## 2. What was ADDED (the brief's "one thing to add")

Three new `~`-level sections in `docs/theory/foundations/boundary_conditions.rst`, placed
immediately after `bc-factor-quotients` as requested:

* **`bc-method-realizability`** — a method is a projection `Π : Γ → Γ_h`; realizability is the
  commuting square `Π ∘ (R G) = (R_h G_h) ∘ Π` (new eq-label
  `bc-realizability-square`). Tier 1 exact-and-faithful (`R = α·I` for **any** α — why the
  diffusion realizer is one line); tier 2 exact-but-NOT-faithful (P1 identifies specular and
  Lambertian — the diffusion realizer's own module docstring quoted verbatim, and tied back to
  the rank-one theorem read "from the other side": the projection, not the response, is what
  destroys the distinction); tier 3 not exact. Closes with **scalar vs angular, NOT trivial vs
  non-trivial**.
* **`bc-equivariance`** — `G` is method-independent as a geometric object; `G_h` exists iff the
  discretization is equivariant, and that splits by which coordinate `g` touches. Specular acts
  on the ANGULAR coordinate ⇒ ERR-042/044/045 are **discretization-admits-the-symmetry** checks,
  not physics checks (all failure mode #5; independence stated with the two catalog-documented
  mutants). A spatial wrap acts only on the SPATIAL coordinate ⇒ every angular discretization is
  trivially equivariant — periodic is the more method-agnostic one, for that sharper reason.
* **`bc-refusal-axes`** — the three INDEPENDENT axes (angular resolution / spatial-topological /
  state-cone), as a `list-table`, each row citing the shipped guard that names it: diffusion
  refuses periodic ("the one geometry P1 cannot integrate away into a per-face albedo" — nothing
  to do with angle), SN refuses zero-flux (`𝒜 = −1` needs a signed current; `ψ ≥ 0`). Explicit
  closing point that the two methods are **incomparable**, not coarser/finer. `q` treated as the
  fourth thing: a **vector in `Γ₋`**, so diffusion's refusal is **plumbing** (#290 P5) and
  vanishes with no theory changing.

Plus one new `-`-level section carrying the narrowing itself:

* **`bc-domain-narrowing`** (+ sub-anchor `bc-narrowing-what-it-removed`) — `γ_S`/`ι_S` definition
  (new label `bc-trace-restriction-pair`), sibling-not-subclass rationale, the face action
  `B_face = ι₋∘law∘γ₊` (new label `bc-face-action-narrowed`), the three-spellings-one-pair table,
  a `.. warning::` on the two measured traps (`searchsorted`-not-`arange` with the complementary
  slab/2-D fixtures; the Mode-12 shape blindness), and the un-narrowed remainder (6 rows, 4 kinds,
  16 strict xfails, B3.4, plus the two carried consequences).

All three new eq-labels are tagged `.. vv-status: … documented` with rationale comments naming
the real gates. They land in the matrix's **Documented-only** bucket (auto-regen 521 → 524; no
orphan added). Self-check: every `vv-status: X` in `docs/theory/**` has a same-file `:label: X`.

---

## 3. Flagged, NOT fixed (outside `docs/theory/**` or outside B3 scope)

**Code / test docstrings (yours — I did not touch `orpheus/` or `tests/`):**

1. `orpheus/geometry/boundary/__init__.py:~135` — layer-2 bullet still says
   "`geometry_map` and `response_kernel` return the ABC's `None` on **every** law and are read by
   nothing — campaign phase **B1** mints … and populates them". **B1 landed** (`a0fd17b4`); this
   is the same falsehood I fixed at `bc-law-layer`.
2. `orpheus/geometry/boundary/__init__.py:43` defines `.. _bc-method-realizability:` **in a
   docstring**. The package is not `automodule`'d anywhere, so it is inert today — but it is a
   latent duplicate-label if the package is ever surfaced. The established pattern (see
   `_factors.py`, which `:ref:`s `bc-factor-roles`) is: **the theory page owns the label, the
   docstring references it.** I used the same label name on the page so the two agree; converting
   the docstring anchor to a `:ref:` is the right end-state.
3. `orpheus/diffusion/boundary_realizer.py` module docstring + `_partial_current_albedo` docstring
   — both say "a specular mirror, a Lambertian average, an identity and **a null map**". `NullMap`
   retired at B3.0.
4. `orpheus/sn/boundary/angular.py::IncomingSourceOperator` docstring — "for the rank-0 case where
   **`R = G = 0`**". B3.0 retired that spelling (`G` is the identity deck element; the zero map is
   not a bijection).
5. `tests/geometry/test_bc_equivalence_snapshot.py` — `TestVacuumLebedev17Snapshot`'s class
   docstring and its method docstring still describe the **pre-B3.2** semantics ("zeroes only the
   inflow ordinates; outflow rows pass through unchanged"; "The realized `IncomingOrdinateMaskTensor`
   masks only inflow rows") while the body asserts the zero map. `TestMixed30Spec70WhiteLS4Snapshot`
   still writes `R = 0.3 G_refl + 0.7 G_diff`.

**Pre-existing doc staleness NOT caused by B3 (reported, deliberately not deep-rewritten):**

6. `docs/theory/foundations/boundary_conditions.rst` `bc-two-bc-applies-per-matvec` (§ ~4340–4470)
   describes the Phase-D-era curvilinear matvec that Wave O deleted; the live gate asserts
   **zero** `bc.apply` calls. I added a historical status banner and past-tensed the two
   B3-touched sentences; a full rewrite is its own task.
7. `:meth:`SNMesh._resolve_one`` survives at 5 sites in that same file (lines ~3202, 4603, 4639,
   4841, 4850, 4930) — all inside dated archaeology where naming the then-current method is
   legitimate, so I left them. The **present-tense** ones on `methods/sn/boundary_conditions.rst`
   were repointed.
8. `docs/theory/foundations/boundary_conditions.rst` still describes `SNMethodSpace.minimal` as a
   general test fixture in places beyond the one I annotated; B3.4 makes it un-realizing and it
   becomes a retirement candidate.

---

## 4. Gate

```
.venv/bin/python -m sphinx -E -W docs docs/_build/html
EXIT=0
grep -c 'WARNING:|ERROR:|CRITICAL:'  ->  0
build succeeded.
```

Baseline (pre-edit, same command) was also EXIT 0 / 0 — **set unchanged, not merely count.**

Cross-ref grep gate (`-W` is blind to these): every new cross-doc `:ref:` verified as a real
`<a href>` in the built HTML from all five referring pages; all five new section anchors and all
three new `equation-*` ids present in the built page.

Files changed (all under `docs/theory/`):

```
 docs/theory/conventions/indexing_and_layout.rst       |   28 +-
 docs/theory/foundations/boundary_conditions.rst       | 1214 +++++++++++++---
 docs/theory/foundations/operator_algebra.rst          |   16 +-
 docs/theory/foundations/operator_tensor_network.rst   |   22 +-
 docs/theory/methods/sn/boundary_conditions.rst        |  144 +-
 docs/theory/methods/sn/curvilinear_one_group.rst      |   36 +-
 docs/theory/verification/matrix.rst                   |    5 +-   (AUTO-REGENERATED by the build)
```
