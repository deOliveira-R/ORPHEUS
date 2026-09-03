# The symmetry machine computes what it claims

**Groups are realized, invariance lives with the measure, the registry names what it
licenses.** The review-driven successor to `angular_spaces_derived_from_symmetry.md`
(#429). Tracking issue **#434**. Branch `fix/angular-phantom-support` (continues #429's
branch; base `c1fca8bd`, the #429 compaction record). Prior art, not restated here:
`.claude/plans/symmetry_quotient_plan.md` (2026-08-31, the theory/placement plan that
preceded #429's execution) and #429's Part XIV.

Governed by `.claude/rules/plan-authoring.md`. Every number below carries its
configuration; a bare number is a plan defect.

---

## ▶ RESUMES AT — read this block first, whole

**STATE 2026-09-03.** R1 IMPLEMENTED on the working tree, gated, documented — the
13-tree gate running at the time of writing (`scratch/_r1_full_gate.log`); the commit
lands when it reads 13 of 13 rc=0. ✅ R1 LANDED `f9d3b15b` (2026-09-03; gate 10737 collected, 13 of 13 rc=0, 59 min 35 s).
✅ R4 LANDED `13423a59` (2026-09-03; gate 10840 collected — predicted
10737 + 103 to the row: 102 architect rows + 4 elegance + 2 witness − 5 retired — 13 of
13 rc=0, 59 min 36 s).
Next act: **R2** (§II.R2) — on the clean tree run `scratch/_r2_edit_production.py`
then `scratch/_r2_edit_tests.py` (anchor fixes expected; the architect's §7 order
re-points `tests/_harness/references.py:31` first and re-collects) → scoped suites →
behaviour contract (`scratch/_r1_behaviour.py`; the two INTENDED moves: gauss_legendre
walks → `{O2_x, D_2h}`, the folded candidate set 20 → 18) → exit re-gate → pyright →
battery `scratch/_r2_mut.py` (19 arms + control; the architect is building it) →
elegance → archivist (17 `is_invariant` role xrefs across 6 files; manifolds.rst's
cycle paragraph; a page for the new module) → sphinx → `dead_references` → gate →
commit `Refs #434`.

**R1's record, for the reader who resumes here.** `[M]` behaviour contract
(`scratch/_r1_behaviour.py`): 9925 answers, 53 moved, all intended. Walk cost, min of 15
interleaved repeats, pristine → R1 (`scratch/_r1_walk_timing_3.log`): cold slab 342 →
344 ms (+0.4 %, the first walk closes the polyhedral groups either way), product(4,8)
155 → 39, lebedev(9) 371 → 72, folded 108 → 4, level_symmetric(8) 1774 → 459; warm
5.9 → 1.6, 137 → 12, 324 → 15, 107 → 4, 1421 → 33 ms — the containment answer is
memoised on the tag pair and membership is linear algebra on the matrices (without the
memo the slab walk was +450 %: `[M]` 3888 `RigidMotion` mints per walk, each paying the
orthogonality check). Elegance review: 4 violations + 6 concerns, all taken (the
realization's guards; the kernel takes the GROUP; `_orbit_space_of`; `orbit_certificate`
reads the realization; D∞h now answers `generic_images`; `images` → `generic_images`;
`conjugated_by` retired for lack of a consumer). Test-architect: 28 gates in 9 classes;
the 45-row α-renamed control retired; M1/M2/M3 denominators re-keyed (1090 → 919; 23 →
22 members, 131 → 108 strict, 23 → 20 touching). Battery 20 arms: 19 bite; the
two-generator arm is UNINSTALLABLE (refused by the guard) — its direct witness is
`TestR1TheRealizationRefusesIllegalStates`. Archivist: +841/−90 over 4 pages, sphinx 0/0,
`dead_references` 0/52, sentinels 587 → 591. Residue filed: #442 (the four-way naming
dispatch → a tag Protocol at the 6th family), #443 (SymPy derivations for the two theorems).
Surprise logged in place (not promoted — a one-off): the done-when `Cn(1) is Trivial` was
false as written; `==` is the contract (§II.R1).

**R4's record, for the reader who resumes here.** `[M]` the invariance kernel is
Mode-12 BLIND to the carve — the chart is the column selection the projector
re-writes — so 0 of 9925 behaviour answers moved and the gates assert at the ambient
tier; exit re-gate 135 of 135 arrays `array_equal`, slab `+4.000000000000` at L = 0..3,
census 57 of 60; pyright 0. Elegance: 1 violation + 5 should-fix + 3 nits, all taken —
`lift_codomain` COMPARED (an entry's identity includes where its lift lands);
`generic_orbit_dimension(points)` is the GROUP's method, asked of the base's probe set
(`manifold._generic_points`); the generic orbit's dimension is the MAXIMUM of
`rank{X p}` over that set (upper semicontinuity), not a single probe's — see the ⛔ in
§II.R4; `__post_init__` states four theorems in order (stabiliser, dimension, the
lift's ambient width, fundamental-domain agreement), each with its own witness.
Test-architect: 102 rows in 9 classes; battery 22 of 22 arms bite, and the
`dim_law_reads_dim_H` arm reddened **0 of 4597** catalogue rows (every shipped entry
has `rank[X p] = dim H`) — hence the two §6c witness rows, `S²/O(3)` is a POINT and
`ℝ³/O(3)` a RAY, the only constructible inputs on which the two laws disagree.
Archivist: +1127/−196 over manifolds.rst, error_catalog.rst, angular_quadrature.rst
(matrix.rst regenerated), every `[M]` re-measured (`embed ∘ select == P_H` from the
realized matrices, `array_equal` on 8 of 8 entries × 41 vectors; harness 31 of 33 rows
unchanged, the 2 movers the y-folds; equivariance breaks by exactly 1.0 under a quarter
turn between axes), sentinels 591 → 593, sphinx 0/0, `dead_references` 0/52; three false
docstring claims in manifold.py found and corrected (repr visibility ×2 — both chart
fields are `repr=False`; the forgery-ordering clause). Surprises, logged in place, not
promoted: (a) this plan's *"no test pinned the name `S^2/sigma_y/Trivial`"* — `[M]` 2
hits, both docstrings of the new gate, so *no assertion* pins it and the sentence as
written did not survive; (b) §II.R4's means spelled the generic point as ONE probe —
refuted by the elegance pass (⛔ at the bullet).

**ORDER (dependency-forced): R1 → R4 → R2 → R3.** R4 reads R1's `dim`; R2 moves the
kernel onto R1's realization and R4's lift; R3's coverage predicate reads R1's
membership test and R2's measure verbs.

**Per-carve protocol** (the #429 discipline, unchanged): phase opener re-measures its
own section's premises against the tree → `test-architect` designs the gates BEFORE the
edit (proactive trigger: the carves cross numerics ↔ quadrature ↔ measure) → main agent
implements → `elegance-enforcer` (every finding taken in the same commit) → `archivist`
(docs re-derived, not transcribed) → scoped suites → 13-tree gate detached to a LOG +
`sphinx -E -W` + `mcp__nexus__dead_references` → commit with `Refs #434` (the last one
`Closes #434`) → push. Predict the gate count before running it; reconcile the delta.

**Exit instruments (re-run at the end, must not move):** `scratch/_exit_regate.py
compare scratch/_exit_base.npz scratch/_exit_head.npz` (135 of 135 arrays
`array_equal` — re-capture HEAD first with `capture - scratch/_exit_head.npz`), `slab`
(LEG A `+4.000000000000` at L = 0..3), `census` (57 of 60 rows, the documented three
exceptions).

---

## 0. Provenance — the review of 2026-09-02

**Who.** The main agent read both modules whole (`orpheus/numerics/symmetry.py` 2455
lines, `orpheus/numerics/manifold.py` 1958) and the four touch points; four independent
reviewers: `elegance-enforcer` (`scratch/_rev_elegance_probe{1,2,3}.py`),
`cross-domain-attacker` (⚠ no Bash in that agent — derivation and reading only, every
claim `[D]`/`[M-read]`, none by execution), `qa` (`scratch/_rev_qa_{1..9}*.py`, a 17-arm
in-process mutation battery over 670 tests in 6 files — `_rev_qa_mut.py`, per-arm logs
`_rev_qa_arm_*.log`, pristine copies diffed clean), `explorer` (AST census + Nexus smell
sweep on a graph rebuilt at `c1fca8bd`; scripts `scratch/_rev_explorer_*.py`).

**The four rulings (user, 2026-09-02, each the recommended option):**
1. **Realize the group** — (identity component, coset representatives, membership
   predicate); every predicate computed; the tag survives for identity, naming and
   constructors; axes become vectors in the realization.
2. **Split** groups from invariance; **the invariance verbs live on the measure**.
3. **Three registry fields + coverage** — spent / unspent / owed; stage 0 =
   descent exists ∧ `X.by ⊆ unspent·spent` (total; retires `spent_group`'s refusal).
4. **The lift is the projector everywhere** — a builder-populated derivation output,
   spelled once; the hemisphere section retires; the dimension law lands.

**What the review did NOT find — a false certification.** `[M]` qa, brute force built
in plain numpy outside the module: `is_normalised_by` 5103/5103 (27 tags × 189
motions), `normalises` 729/729 (27²), `contains` 575/575 (24², after the probe's own
improper "SO(3)" sampler was fixed — 77 of 200 QR samples had det < 0), `orbit_stabiliser`
24/24 against a genuine maximum search (17 points × 360 stabiliser samples), the
vv-principles #15 compatibility law 1260/1260 (9 fixtures), barycentre equivariance
183/183 at 6.66e-16. The defects are all in the value types, the boundary, and claims
the code no longer honours.

**Withdrawn attacks (first-class, with the reason):** the six-key catalogue as a
"table" (2 procedures, fields pinned symbolically); `descending_slots` and `reference`
on `Quotient` (both derivation outputs); the hemisphere sign under-gated (provably
invisible to the kernel — every normalising `g` is block-diagonal); the orbit-space
reading over-certifying (its only divergence from the ambient reading is conservative:
`[M]` `C_4`, `D_4h`, `O_h` ambient-True / folded-False on an octahedral set); the
barycentre off `S²` (no production consumer needs a point of `S²`; the kernel is linear
and reads only the chart); `manifold` and `symmetry` both knowing about axes
(`AXIS_LETTER` single-sourced, the inverse derived).

---

## I. Findings ledger → owner carve

| id | finding | evidence `[M]` (configuration) | owner |
|---|---|---|---|
| A1 | the group is a tag dispatched on by family | 66 `isinstance(tag, T)` sites in 14 functions + 23 `is _NamedSubgroup.X` tests; 31 of them in `_contains` alone — explorer, AST over `symmetry.py` @ `c1fca8bd`, predicate = `isinstance` against the tag set ∪ `Compare` against `_NamedSubgroup.X` | R1 |
| A2 | the continuous half of the lattice is still hand-tabulated | `_NAMED_LATTICE` 8 edges; `_axial_contains`; `_continuous_decomposition`; `_identity_component_normalises`; the finite half is already computed (`_finite_contains`, `_fixes_axis`) | R1 |
| A3 | `SubgroupOfO3` is mutable | `__slots__ = ("_tag",)`, plain `__init__`; `g._tag = …` succeeds and moves `hash(Quotient)` (elegance probe1 P1/P2) | R1 |
| A4 | `Cn(1)` ≠ `Trivial` — one group, two values | each `contains` the other, unequal; `SPHERE.quotient(Cn(1))` → "no catalogue entry" (the message the door promises never to give) | R1 |
| A5 | `_maximal` tests `!=` where its docstring says "strictly contained"; `maximal_invariance_groups` returns `()` for a measure whose symmetry is `{e}` | elegance probe2 on an asymmetric 9-node sphere rule, walk and bruteforce both `()`; **0 production callers** of the walk (explorer) | R2 |
| A6 | `identity_component` returns the group itself for finite members | wrong on 12 of 22 members (qa probe 3); violates its own docstring's property on 11; **0 readers** tree-wide (3 filters + control) | R1 |
| A7 | stale `_contains` comments ("`Cn ⊂ Oh` out of scope") | `[M]` `O_h ⊇ C_4`, `D_4h`, `C_2`, `D_2h` True today (both finite → computed); the named-vs-finite arms' `inner is Trivial` clauses are unreachable | R1 |
| A8 | the "brute-conjugation CONTROL" test is production α-renamed | AST-unparsed + α-normalised: identical; its element list is production's `_group_elements` (qa) | R1 |
| A9 | tolerance twins | `_MEMBERSHIP_ATOL` = 1e-9 (element equality, symmetry) vs 1e-12 (point on manifold, manifold); `_X/_Y/_Z_AXIS` beside `_axis_vector`; `SubgroupTag` 0 readers; `is_subgroup_of` 0 production callers, 14 test sites | R1 |
| B1 | `Quotient.lift` is a three-arm branch on the group's tag; its `NotImplementedError` says "add … to `Quotient.lift`" | every OTHER derivation output is a field; 9 functions in 4 modules re-derive the family from the tag (elegance, AST) | R4 |
| B2 | the dimension law is unenforced | forged `S^2/O2_z` realized on `Ball(2)` and `S^2/sigma_x` on `[-1,1]` both construct (probe1 P5/P7); the law holds on 10 of 10 legitimate entries | R4 |
| B3 | `(M/H)/{e}` builds a second `Quotient` | `S^2/sigma_y/Trivial`, unequal to `S^2/sigma_y` — the class the door refuses for `G ⊆ H` | R4 |
| B4 | the lift is the Reynolds projector on all three arms | axial: `P p = (p·ê_a)ê_a`; mirror: zero column `a`; trivial: `I` (cross-domain, derived; equivariant under the normaliser because `g P_H g⁻¹ = P_{gHg⁻¹}`) — so `_hemisphere_section`, its `sqrt` and the `ρ² > 1 + 1e-12` refusal are unnecessary | R4 |
| B5 | `name == "Trivial"` string compares | 3 in manifold.py, 2 in `basis/descent.py` | R4 |
| B6 | `ambient_representatives` promises a representative and returns a barycentre; 8 locals named `section` in symmetry.py for the thing the entry says it lacks | elegance F7 | R4 |
| B7 | `_mod_trivial`'s docstring names a consumer that does not exist | `registry.support` returns `SPHERE`, never the entry; 0 production `SPHERE.quotient(Trivial)` | R4 |
| C1 | "ONE closure" is false in the code | `_invariance_on_orbit_space:1709-1716` inlines `_orbit_space_closure:1584-1594` (identical lambda); three docstrings claim otherwise (`symmetry.py:1567`, `:2390`, `directional.py:518`) | R2 |
| C2 | kernel step 2 is a theorem (0 reds deleted); step 1 duplicates `induced_action`'s guard (0 reds with the raise swallowed, 16 with it escaping) | qa battery M5/M5b/M6, 2290 kernel calls | R2 |
| C3 | `candidate_groups` reads node STORAGE width | one fold, two spellings: `{D_2h}` vs `{σ_x, σ_y, σ_z}` (qa); selection unaffected | R2 |
| C4 | `_distinct_azimuths` hard-codes `1e-9` shadowing its own `atol`; `_fixes_every_point` is fed the weight window (1e-13) while node matching uses 100× it | qa census | R2 |
| C5 | present-tense-false docstrings: module docstring ("fingerprint strategy", "static lookup table", "Issue 5 … will install"); `is_invariant`'s `O_h`/`I_h` bullets describe retired code, the `I_h` one an ERR-072-shaped sample | AST-unparsed call chain contains none of `sorted/fingerprint/multiset/icosahedron/vertex` (qa) | R1 (module), R2 (kernel) |
| C6 | ERR-045's message misdiagnoses the fold's identity permutation | `_specular.py:252` on `folded_product(4,8)`, axis y | R2 |
| C7 | the deferred / `TYPE_CHECKING` cycle notes are the boundary's symptom | manifold.py duck-types 7 attributes of `SubgroupOfO3`; measure.py imports the certificate at function scope | R2 |
| C8 | `_identity_component_normalises`'s five arms rest on 1–2 tests each; arm (c) invoked ONCE in 670 tests | qa per-arm table | R1 (dissolved) + gates |
| D1 | Γ recorded as a closure requirement, read as a fold licence | `folded_product(4,8)` admitted at both stages for `cartesian2d`; 2 of 4 `(sign μ_x, sign μ_y)` sweep quadrants empty (qa probe 9); bounded — not in the default registry, `select_quadrature("cartesian2d", 8)` → `lebedev(9)` | R3 |
| D2 | `continuous_isotropy` documented as `G⁰`; `O(2)_x` is disconnected | registry.py:834-847, 1052-1054 (its own comment concedes the Klein residual) | R3 |
| D3 | `spent_group` is not a functor; its coset arm raises where coverage is total | manifold.py:1422-1428 | R3 |
| D4 | two `vv-principles` anti-patterns owed (qa): a symmetry recorded for one job spent on another; the α-normalised AST check for a "brute-force control" | — | R3 docs pass |

---

## II. The carves

### R1 — Every question about a group is computed from its realization

**Goal.** A closed subgroup of `O(3)` is represented by what it IS — its identity
component and its coset representatives — and every predicate (`contains`,
`is_normalised_by`, `normalises`, `identity_component`, `dim`, `generic_images`, "G⁰
fixes every node") is ONE body reading that representation. Nowhere is a relation
between two groups written down.

**The theorem the design rests on** `[D]`: `so(3)` is simple and 3-dimensional; its
subalgebras are `{0}`, `ℝ·[â]ₓ` (one per axis `â`), and `so(3)` — none of dimension 2.
So a compact `G ≤ O(3)` is exactly the pair `(𝔤, {coset representatives of G⁰ in G})`.

**Means (design, 2026-09-02 — the opener re-measures nothing here beyond the §6b set,
because the review already measured the tree).**

```python
@dataclass(frozen=True, eq=False)
class IdentityComponent:
    """G⁰ = exp(𝔤): 𝔤 has dimension 0, 1 or 3 (so(3) has no 2-dim subalgebra)."""
    generators: tuple[NDArray, ...]        # skew-symmetric basis of 𝔤; () / ([â]ₓ,) / (Lx, Ly, Lz)
    dim; axis (the torus axis when dim == 1, else None)
    contains_element(h, atol)   # dim 0: h ≈ e; dim 1: det h > 0 and h â = â; dim 3: det h > 0
    contains(other)             # 𝔥 ⊆ 𝔤 — a 3×3 case table on dimensions + axis parallelism
    is_normalised_by(g, atol)   # Ad_g 𝔤 = 𝔤  (dim 1: g â = ±â)
    normalises(other: Realization, atol)   # EXACT: [𝔤, 𝔥] ⊆ 𝔥  ∧  X − Ad_s X ∈ 𝔥 for every rep s of other
    fixes(points, atol)         # X p = 0 for every generator, every point  (dim 1: p on the axis; dim 3: p = 0)
    images(points)              # exp(θ_k X) for the incommensurate θ (dim 1); () for dim 0; refuse for dim 3
    conjugated_by(g)

@dataclass(frozen=True, eq=False)
class Realization:
    """G = ⊔ r·G⁰ — the identity component and one representative per component, identity first;
    for a finite group every element."""
    component: IdentityComponent
    representatives: tuple[RigidMotion, ...]
    is_finite; elements (finite only)
    contains_element(h)   = any(component.contains_element(r⁻¹ @ h) for r in representatives)
    contains(other)       = component.contains(other.component) and all(self.contains_element(s) for s in other.representatives)
    is_normalised_by(g)   = component.is_normalised_by(g) and all(self.contains_element(g r g⁻¹) for r in representatives)
    normalises(other)     = component.normalises(other) and all(other.is_normalised_by(r) for r in representatives)
    conjugated_by(g)
```

*Why `component.normalises` is exact and not sampled* `[D]`: `G⁰` normalises `H` iff
`𝔤 ⊆ Lie(N(H))`, and `Lie(N(H)) = {X : [X, 𝔥] ⊆ 𝔥 and X − Ad_s X ∈ 𝔥 ∀ s ∈ H}`. Necessity is
the derivative at `t = 0` of `exp(tX) s exp(−tX) s⁻¹ ∈ H⁰`; sufficiency: with
`Y_s = X − Ad_s X ∈ 𝔥`, the curve `f(t) = exp(tX) s exp(−tX) s⁻¹` has
`f'(t) f(t)⁻¹ = Ad_{exp(tX)} Y_s ∈ 𝔥` (since `[X, 𝔥] ⊆ 𝔥`), so `f` stays in `H⁰`. On a
finite `H` (𝔥 = 0) this is "every element commutes with `X`" — today's
`_identity_component_normalises`; on `O(2)_b` it is `a ∥ b` (`[[â]ₓ,[b̂]ₓ] = [â×b̂]ₓ ∈
span[b̂]ₓ ⟺ â ∥ b̂`) plus `X − Ad_σ X = 2X ∈ 𝔥` — today's axis-equality arm. One body,
both answers.

- `SubgroupOfO3` becomes `@dataclass(frozen=True, slots=True)` with the tag as its one
  field; `__post_init__` normalises `Cn(1)` to the `Trivial` tag (Pattern 4: the second
  spelling is unrepresentable, not merely unconstructed). The generated `__eq__`/`__hash__`
  replace the hand-written pair. `realization` is a cached property derived from the tag by
  ONE `match` (`_realize(tag)`) — the sole remaining discrimination beside `name`/`__repr__`
  and the constructors.
- `dim` = `realization.component.dim`; `identity_component` → the NAMED member of the
  component (`Trivial` / `SO2(letter)` / `SO3`; refuse an unnamed axis); `orbit_stabiliser`
  structural: finite → self; component `SO(3)` → `O3`; torus about `a` → `O2(a)` if
  `O2(a) ⊇ self` else self (`[D]` the orbits are unions of circles about `a`; the
  partition's stabiliser is `O(2)_a` when `G` acts trivially on `μ`, `D_∞h(a)` otherwise —
  the only shipped instance of the latter is `Dinfh` itself).
- `rotation_axis` / `mirror_axis` stay as the tag-level readers the basis modules use
  (their retirement is #438's / the descent registry's business, not this carve's).
- **Name ruling (main agent, 2026-09-03; the test-architect's R-1).** `Realization` /
  `IdentityComponent` KEPT. `[M]` "realization" is taken 32× in the prose corpus by the
  operator campaign's third axis (an abstract operator realized as a kernel) and by
  `Quotient.realization` (the manifold an orbit space IS); but it is also this module's
  own prior vocabulary (`_realized_ops`, "every realization is in the standard setting",
  the tests' `realization` rows) and the textbook term for a matrix realization of a
  point group — the same concept, one tier down: an abstract object made concrete. The
  `SubgroupOfO3.realization` docstring disambiguates from the manifold's; the trigger
  to rename is R2 writing both in one frame (then `component_decomposition`).
- Tolerances: `_MEMBERSHIP_ATOL` (symmetry) → `_ELEMENT_ATOL = 1e-9`, the one band for
  every element-level comparison in the realization; `_X/_Y/_Z_AXIS` → `_axis_vector`.
- Kept, because 11 test files import them by name (`[M]` explorer §2c): `_realized_ops`
  (the generating sets — now what `_realize` closes), `_close_group`, `_group_elements`
  (`realization.elements`), `_cyclic_ops`, `_vertical_mirrors`, `_axis_vector`,
  `_INCOMMENSURATE_ANGLES`.

**Retired:** `_NAMED_LATTICE`, `_named_contains`, `_contains` (109 lines),
`_finite_contains`, `_fixes_axis`, `_axial_contains`, `_rotation_generator`,
`_maps_axis_to_itself`, `_continuous_decomposition`, `_fixes_every_point`,
`_identity_component_normalises`, `_is_axis_supported`, `_is_origin_supported` (the last
three's exact criteria become `IdentityComponent.fixes`), `SubgroupTag`, `is_subgroup_of`
(14 test sites → `contains`), the module docstring's "fingerprint" / "static lookup table"
paragraphs.

**Behaviour contract — captured BEFORE the edit, compared AFTER, bit-identical except
the three intended changes:** (a) the full containment table over every expressible
member (the 24 tags of qa probe 4 + `Cn/Dnh` at n = 1..8); (b) `is_normalised_by` over
qa probe 1's 27 × 189; (c) `normalises` over 27²; (d) `orbit_stabiliser` × 24;
(e) `is_invariant` over (every candidate group × the 12 shipped rules of
`scratch/_exit_regate.py`'s `ALL_RULES`); (f) `maximal_invariance_groups` on the 5 rules
of the review's probe (`[M]` today: lebedev(9) → `(Oh,)`, product(4,8) → `(Dnh(8),)`,
gauss_legendre(8) → `(Mirror('x'), O2('x'))`, folded_product(4,8) → `(Dnh(2),)`).
**Intended changes:** `Cn(1) == Trivial`; `identity_component` of a finite member is
`Trivial`; `SO2('x') ⊇ Cn(1)` stays True (was already computed). Walk cost: re-measure
with the min-of-15 interleaved protocol (`scratch/_p19_probe6.py`'s shape), not one draw.

**§6b set** `[M]` explorer: production readers of `SubgroupOfO3` outside the two modules
— 30 sites / 10 files (constructors and `contains`/`is_invariant` only; none touches the
tag); tests importing privates (§2c): 20 sites / 11 files, of which `_MEMBERSHIP_ATOL`
(`test_symmetry.py:2404`) and `is_subgroup_of` (14 sites, 2 files) change spelling.

**Done when:** `grep -c "isinstance(tag" symmetry.py` ≤ 3 (naming + `_realize`)
(⛔ REFUTED 2026-09-03 as a number, kept as intent: `[M]` 26 sites in 10 owners remain,
19 of them the four-way naming/constructor dispatch the ruling keeps — #442 tracks its
collapse; the STRUCTURAL dispatch is one `_realize`); `_NAMED_LATTICE` absent; the
before/after tables agree except the three intended rows;
`SubgroupOfO3.Cn(1) == SubgroupOfO3.Trivial` (⛔ the first draft said `is`; `[M]`
test-architect: `__post_init__` normalises the TAG, so the constructor returns a fresh
value, equal and same-hash, not the singleton — `is` would red a correct implementation);
`FrozenInstanceError` pinned; the qa
brute-force probes 1, 2, 4 promoted to gates with the SAME denominators (27 × 189, 27²,
24², 24) and an INDEPENDENT construction (not `_group_elements` — build the element sets
from the tag's definition in the test); the C8 arms each reddened by a gate that names
the arm; 13-tree gate green; sphinx `-E -W` 0; `dead_references` 0 dead.

**Docs (archivist):** `docs/theory/foundations/manifolds.rst`'s normaliser chapter
(`manifold-normaliser`) re-derived on the realization (the theorem above); the lattice
section; `discrete_measures.rst`'s "Symmetry groups" section; the symmetry module
docstring (module docstrings are read by fresh sessions — the "fingerprint"/"lookup table"
prose is present-tense false and is a bug).

### R4 — The lift is a derivation output, and an orbit space's dimension is a theorem

**Goal.** A `Quotient` carries its lift the way it carries its chart — as a field the
derivation populated — and cannot be constructed with a chart of the wrong dimension.

**Means.**
- `Quotient.lift_coordinates: Callable[[NDArray], NDArray]` and `lift_codomain: Manifold`
  (fields beside `orbit_coordinates`); `lift` (property) assembles the typed arrow.
- ONE helper spells the pair for a column-selection chart:
  `_coordinate_chart(columns, ambient) -> (select, embed)` — `select = points[:, columns]`
  (today's `_ambient_columns`, bit-identical), `embed = scatter the chart columns into a
  zero vector`. `embed ∘ select = P_H`, the Reynolds projector onto `(ℝ³)^H` (B4). The two
  sphere builders call it with their invariant columns; `_mod_trivial` passes the identity
  with codomain the base. A future entry whose invariants are not linear (`S²/C_n`, #439's
  `ℝP²`) supplies its own pair — the field is REQUIRED, so it cannot be forgotten.
- `ambient_representatives` → **`orbit_barycentres(points)`**: chart-width → `lift_coordinates`;
  ambient-width → `lift_coordinates(orbit_coordinates(points))` = `P_H p`. One concept,
  named for what it returns on every entry (B6); the `section` locals in symmetry.py →
  `barycentres`. `barycentre(entry)` (public, `functools.cache`) = `entry.lift` for EVERY
  entry — the named arrow, no longer axial-only; its docstring becomes the projector's.
- `_hemisphere_section` retired with its literal and its `sqrt` (B4). `[M]` its only
  reader is `lift`; the mirror entries' `fundamental_domain` still validates a fold's
  representatives (`Quotient.contains`), unchanged.
- **Dimension law in `__post_init__`:** `realization.dim == base.dim − dim(generic orbit)`,
  with `dim(generic orbit) = rank[X p : X ∈ 𝔤]` at a generic point `p` of the base
  (`_GENERIC_SPHERE_PROBE[0]` for the sphere; a generic vector for `RealSpace`/`Ball`;
  0 for a finite group without a point). ⛔ REFUTED 2026-09-03 (elegance pass): ONE
  probe is not generic — `[M]` with the probe on the axis the law both refused the
  honest `S²/O2_z` and admitted the disk forgery; shipped as the MAXIMUM of
  `rank{X p}` over the base's probe set (`manifold._generic_points`,
  `SubgroupOfO3.generic_orbit_dimension`), the maximum being the generic value by
  upper semicontinuity of rank. `[D]` `O(2)_a` on `S²`: rank 1 → 1 ✓;
  finite: 0 → 2 ✓; trivial on any base: 0 ✓; `O(3)` on `S²` (#440): rank 2 → 0 ✓;
  `SO(3)` on `ℝ³`: rank 2 → 1 ✓ — `dim H` alone is WRONG for the last two, which is why
  the law is stated on the orbit. Negative legs: the two forgeries of B2.
- `_catalogued_quotient`: a `Quotient` base under `Trivial` returns the base (B3); `group.is_trivial`
  replaces the three `name == "Trivial"` compares here and the two in `descent.py` (B5).
- `_mod_trivial`'s false consumer claim corrected (B7).

**§6b set** `[M]` explorer: `Quotient.lift` 0 production readers outside the two modules
(17 test sites, 1 file; 20 doc mentions in manifolds.rst); `ambient_representatives` 0
production outside (6 test sites; 5 doc roles); `barycentre` 0 production outside (6
test sites; 10 roles + 47 prose in manifolds.rst — archivist); `_ORBIT_CATALOGUE` and
`_sphere_mod_o2` imported by tests (`test_manifold.py:671,1834`, `test_slab_orbit_space.py:299`).

**Done when:** the six catalogue entries' `orbit_coordinates` are `array_equal` to
HEAD's on the exit instrument's 12 rules (chart bit-identity — the frame tables must not
move); the two forged entries are REFUSED at construction naming the law; `(S²/σ_y)/{e}
is S²/σ_y`; `_hemisphere_section` absent; `lift(chart)` for the mirror entry equals
`P_H(section)` for every fold node (the old section's projection — the re-keyed
`test_a_mirror_entrys_lift…` gate); the kernel's answers on (candidate × 12 rules)
unchanged from R1's table; 13-tree gate; sphinx; `dead_references`.

**Docs:** manifolds.rst's lift / barycentre / hemisphere-section passages (the "not a
section" caveats become the projector's one sentence); error catalogue II.11 wording;
`frame.rst`, `matrix.rst` mentions of `ambient_representatives`.

### R2 — Invariance is the measure's question; groups import geometry only

**Goal.** The import graph reads like the mathematics: `geometry.transformation` ←
`symmetry` (groups) ← `manifold` (a quotient is a manifold and a group) ← `measure` (a
measure lives on a manifold) ← `invariance` (a measure × a group). No deferred import,
no `TYPE_CHECKING` cycle note, no duck-typed group.

**Means.**
- New module `orpheus/numerics/invariance.py`: the ONE closure
  `orbit_space_closure(measure, motions, atol) -> OrbitCertificate | None` (normaliser
  check per motion, chart, induced action, `_orbit_closure`), `OrbitCertificate`,
  `_orbit_closure` (now REQUIRING the action — the `images_of=None` default is dead after
  the twin collapses), `_NODE_WINDOW_FACTOR`, `_embedded_nodes` → `orbit_barycentres_of(measure)`
  (tests import it from 6 files — keep the old name as the function's name or migrate
  the 7 sites; decide at the opener), `_ambient_orbit_space`, the kernel
  `is_invariant_under(measure, group, atol)` = `group.realization.normalises(H) ∧
  G⁰.fixes(barycentres) ∧ closure(reps) is not None` — three steps, no tag, no step 2
  (C2: a theorem the closure re-proves), the normaliser asked ONCE (C2's duplicated guard:
  the closure's per-motion check is the predicate; `induced_action`'s refusal stays as the
  type's contract), `candidate_groups` on the EMBEDDED nodes (C3), `_distinct_azimuths`
  taking its window from its argument (C4), the position test at the NODE window (C4),
  `_maximal` on strict containment (A5), `maximal_subgroups`, the walk.
- `DiscreteMeasure` gains the verbs (module-scope import of `invariance` — no cycle:
  `invariance` imports `measure` under `TYPE_CHECKING` only): `is_invariant_under(group,
  *, atol=1e-13)`, `certificate_under(group)`, `permutation_under(motion, *, atol)`,
  `singular_set_under(group)`, `symmetry_groups(*, atol, candidates, method)`;
  `quotient` reads `self.certificate_under`. `SubgroupOfO3.is_invariant` is DELETED (not
  a façade — a façade would keep one deferred import alive).
- `manifold.py` imports `symmetry` at module scope; `Quotient.by: SubgroupOfO3` typed;
  `_trivial_group` retired; the module docstring's cycle paragraph deleted (with the
  `[M]` import-order measurements kept as history in the docs page, not the module).
  ⛔ **R2 opener (test-architect, 2026-09-03): as written R2 does not import.**
  `symmetry.py:105` reads `AXIS_INDEX, AXIS_LETTER` from `manifold` at module scope at
  6 sites that do not move, so `manifold → symmetry` at module scope closes a 2-cycle —
  `[M]` shadow tree, one subprocess per (variant, entry point): 3 of 9 entry points
  clean, the rest dying `ImportError: cannot import name 'AXIS_INDEX' from partially
  initialized module`. **Ruling (main agent): the axis table moves BACK to `symmetry.py`**
  (its home until 2026-09-02; `manifold.py:1060`'s own comment records the move and its
  reason — the import direction — which R2 reverses); `manifold`, `basis/descent.py`
  and `basis/spherical_harmonic_basis.py` import it from there (2 production + 4 test
  readers). `[M]` with that move plus a real `invariance.py`: 10 of 10 entry points.
- ⛔ **C3 moves shipped answers the plan did not list** (`[M]` test-architect): with
  `candidate_groups` on the EMBEDDED nodes, `gauss_legendre(2/8/16)`'s walk goes
  `{O2_x, σ_x}` → `{O2_x, D_2h}` (correct: the embedded `(μ, 0, 0)` set has two
  z-azimuths, `D_2h ⊇ σ_x`, and `D_2h` acts on `S²/O(2)_x` by `μ ↦ ±μ`), and
  `folded_product(4,8)`'s candidate set shrinks 20 → 18. Ruling: ACCEPT (the walk has 0
  production callers; the answer is more informative) and re-key the two gates.
- Three of R2's four behaviour changes are INERT on every shipped rule (`[M]` position
  window 0 of 15, azimuth window 0 of 15, `_maximal` strictness 0 of 31) — each lands
  with a MANUFACTURED witness (the gates draft) or it lands unfalsifiable; the deleted
  step 2 is the opposite shape (32 shipped rows, a MUST-STAY-GREEN table).
- `_distinct_azimuths`' window: `atol · _NODE_WINDOW_FACTOR` (a position question, the
  node window — identical on all 15 shipped rules either way). `_embedded_nodes` keeps
  its name (5 importers incl. the harness). Continuous group ⇒ `certificate is None`
  confirmed as policy (`[M]` 6 of 135 rows are `invariant=True, certificate=None`).
- `Quadrature.ordinate_permutation` → `self.measure.permutation_under(motion3.linear_part,
  atol=_REFLECTION_ATOL)`; `AngularSymmetry.admits_symmetry` → `measure.is_invariant_under(…)`.
- `orbit_certificate` for a continuous group stays `None` by policy (its docstring says
  why: the discrete part's certificate is not the group's).
- C5 (the `is_invariant` strategy bullets), C6 (the ERR-045 message names the fold).

**§6b set** `[M]` explorer §1a/1b/2c: `SubgroupOfO3.is_invariant` — 1 production
caller (`registry.py:1012`), **91 `.is_invariant(` CALL sites / 5 files** (⛔ the
explorer's "170 attribute sites / 3 files" counted attribute LOADS incl. docstrings;
the test-architect's AST call census is the migration's denominator) (rewrite by AST with
source-span replacement, receiver ↔ first argument: `G.is_invariant(m, …)` →
`m.is_invariant_under(G, …)`); `orbit_certificate(` 1 production (`measure.py:1169`) + 8
test sites / 2 files; `singular_set(` 11 test sites / 4 files; `maximal_invariance_groups(`
15 / 2; `candidate_groups(` 11 / 1; `induced_permutation(` 1 production
(`directional.py:548`); `_embedded_nodes` 7 test sites / 6 files incl.
`tests/_harness/references.py:31`; `from orpheus.numerics import symmetry` at
`test_manifold.py:1359`; `numerics/__init__.py:43` re-exports `SubgroupOfO3` (unchanged).
`tests/test_layer_imports.py` permits every new edge (`numerics → geometry` allowed;
none of these modules is whitelisted). ⚠ §6d: inject and RUN the new module-scope edges in
all three import orders before trusting the AST — the last census that said "no cycle"
was blind to relative imports (plan-authoring 2026-08-31).

**Done when:** `grep -rn "TYPE_CHECKING" orpheus/numerics/manifold.py` finds only the
`exactness`/`transformation` type imports (or none); the CYCLE-MOTIVATED deferred imports
are gone — `measure.py`'s function-scope `orbit_certificate`, `manifold.py`'s
`_trivial_group`, `symmetry.py`'s `measure` import — (⛔ REFUTED as "no function-scope
imports at all": `[M]` `measure.py` legitimately defers 6 imports R2 does not touch); `is_invariant_under ≡ certificate_under is not None`
for every finite group on every shipped rule (`[M]` the R2 gate); the asymmetric rule's
`symmetry_groups()` is `(Trivial,)`; the fold reports `(Dnh(2),)` in both spellings; the
kernel's (candidate × 12 rules) table unchanged from R1; step 2's deletion reddens nothing
and the docstring SAYS it is an optimisation; 13-tree gate; sphinx; `dead_references`
(the `:func:`/`:meth:` roles of every moved name — 22 doc mentions of `lift` alone).

**Docs:** manifolds.rst / discrete_measures.rst "one closure" chapter re-pointed to
`invariance`; the module docstrings of all four modules; the layer page in
`development.rst` if it lists `numerics` modules.

### R3 — The registry names what it spends, what it leaves, and what it owes

**Goal.** `AngularSymmetry` records three different facts as three fields, and stage 0
asks the fold question of the field that licenses folds.

**Means.**
- `AngularSymmetry(spent, unspent, owed)` — `spent`: the stabiliser the reduction
  integrates away (the orbit space's group; was `continuous_isotropy`); `unspent`: the
  finite symmetry the SOLUTION still has in the local frame (licenses a fold; new);
  `owed`: the closure a reflecting face needs (was `discrete_residual`). ⚠ the name
  `unspent` is the proposal — confirm at the opener (candidates: `retained`, `even_under`);
  the docstring drops `G⁰` / `G/G⁰` (D2: the pair is a factorisation, and `O(2)_x` is
  disconnected).
- Table `[R]`, to be confirmed against each geometry's transport equation at the opener:
  slab `(O2('x'), Trivial, Mirror('x'))`; sphere the same; cylinder `(Trivial,
  Mirror('y'), Dnh(2))` — the local-frame azimuthal reflection the shipped
  `folded_product` folds by; cartesian2d `(Trivial, Mirror('z'), Dnh(2))` — z-uniform ⟹
  `ψ` even in `μ_z`; the σ_y fold is NOT licensed (D1).
- Stage 0: `descent = quotient_onto(self.support, measure.support)`; admitted iff
  `descent is not None and covered(measure.support.by, by=self.unspent, over=spent)` where
  `H ⊆ Γ·K ⟺ H⁰ ⊆ K ∧ ∀ rep r of H: ∃ γ ∈ Γ: K.contains_element(γ⁻¹ r)` (`[D]`: `ΓK` is a
  finite union of closed cosets; `H⁰` connected through `e` lands in `K⁰`). Total — no
  join, no coset space, no new registry field; `spent_group` retired (D3). Reproduces:
  slab rule on the slab `O2_x ⊆ {e}·O2_x` ✓; σ_y fold on the cylinder `σ_y ⊆ σ_y·{e}` ✓;
  σ_y fold on cartesian2d `σ_y ⊄ σ_z·{e}` ✗ (D1 closed); σ_z fold on cartesian2d ✓.
- The `vv-principles` skill gains qa's two anti-patterns (D4): "a symmetry recorded for
  one job spent on another" and the α-normalised AST check for a "brute-force control".

**§6b set** `[M]`: `continuous_isotropy` / `discrete_residual` — the 4 table rows
(`registry.py:1056-1071`), `support`/`reference`/`admits_domain`/`admits_symmetry`, the
registry tests (`test_registry.py` G12–G14 and the pre-existing rows), `spent_group` 1
production + 0 tests + 3 doc roles; the stage-0 message; the registry's Sphinx page.

**Done when:** the σ_y fold is REFUSED for cartesian2d at stage 0 naming `unspent` and
ADMITTED for the cylinder; a σ_z fold is admitted for cartesian2d; `admits_domain` is
total over (every shipped rule × every geometry) with no `NotImplementedError`; the slab
and sphere admissions unchanged; `spent_group` absent; 13-tree gate; sphinx;
`dead_references`; `Closes #434`.

---

## III. Verification instruments shared by every carve

- **The behaviour tables** (R1's (a)–(f)), captured at the R1 opener from HEAD into
  `scratch/_r1_before.npz`/`.json` and compared after EVERY carve — the machine's
  observable answers must not move except where the plan says so.
- **The #429 exit instruments** (`scratch/_exit_regate.py`): re-run at R4 (charts) and at
  the end.
- **The mutation batteries** (qa's `_rev_qa_mut.py` shape — in-process `pytest_configure`
  plugin, pristine copies diffed after): re-run per carve on the moved code; a green under
  a real mutant is a finding (C2's two arms are the standing example).
- **The 13-tree gate**: last measured 10752 at `a7c8de6d` (59 min 27 s, serial,
  `-m "not slow"`); predict each carve's delta from the gates it adds and retires.

## IV. Residue filed 2026-09-02

#435 engine scope by Chevalley–Shephard–Todd · #436 face pairings on `FundamentalDomain`
· #437 `Sym(μ)` from the node set · #438 character-sum gate for `descending_slots` ·
#439 `{e, −I}` and `ℝP²` · #440 `S²/O(3)` a point · #441 a non-polynomial exactness
reference. Also noted, not filed: `SubgroupOfO3` is 3-D while `RigidMotion` is
dimension-generic (a `Circle` quotient by a non-trivial group is a latent width
mismatch; R1's realization is dimension-generic by construction, the tag is not).

---

## ⏸ COMPACTION POINTS

**After R4 (2 of 4)** and **after R3**: the phase → commit table, corrections superseding
older text, the measured gate counts and costs, the surprise log rows owed to
`plan-authoring.md`.

## Landmines

1. `is_normalised_by`/`normalises` are EXACT today (qa 5103/5103, 729/729); R1 must not
   trade exactness for generality — the normaliser theorem in §II.R1 is the exact
   criterion, and the gates carry the same denominators.
2. The frame tables read `quotient_onto → target.orbit_coordinates`; R4 must keep the
   chart a column selection, bit-identical (the 135 arrays).
3. `tests/_harness/references.py:31` imports `symmetry._embedded_nodes` — the harness
   is loaded by every test module; break it and the whole suite dies at collection.
4. `numerics/__init__.py:43` re-exports `SubgroupOfO3`; `tests/test_layer_imports.py`
   gates `geometry → numerics` by SUBMODULE — a new `numerics/invariance.py` must not be
   imported from `geometry`.
5. Any timing (the walk's cost under the realization) is a stochastic measurement: min of
   ≥ 15 interleaved repeats, or it is one draw (plan-authoring 2026-09-02).
