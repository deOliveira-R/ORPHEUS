# Angular spaces are DERIVED from symmetry, and no object claims a property it does not have

**Status:** OPEN. Diagnosis COMPLETE and converged from four independent
sources. No production line changed yet.
**Tracking issue: #429** (OPEN — the work item, not a record). It is a
POINTER to this file; this file is the artifact. **#426 is BLOCKED by it**
(comment recorded on #426, 2026-08-31).
**Opened:** 2026-08-31. **Supersedes** the first draft of this file
(`angular_basis_and_measure_compatibility.md`, landed `2fbae24a`, renamed here);
that draft's content is absorbed in full — nothing was dropped, and its refuted
premises are kept in Part VII per `plan-authoring` §3.

**Sources merged, and who owns what:**

| source | contributes | where |
|---|---|---|
| the P_L defect investigation (main agent, this session) | the measured defect, the root cause, the leak inventory, the phantom census | Parts I–II |
| the user's **Symmetry & Quotient Layer** plan (2026-08-31) | the theoretical spine, the three gates, the placement ruling, phases 5–7 | Parts III–V, Phases 5–7 |
| `scratch/n2n_pl_consistency_literature.md` (555 lines, page-cited) | Morel/Galerkin quadrature; Bell & Glasstone's 1-D derivation | §III.6 |
| `scratch/n2n_pl_frames_attack.md` (858 lines) + `scratch/n2n_pl_blast_radius.md` (869 lines) | the well-posedness condition; the tree census | §III.5, Part VI |

> ⚠ **Read `plan-authoring.md` §1–§3 before editing.** Every claim carries an
> epistemic marker: `[M]` measured (with its configuration), `[R]` reasoned and
> unverified, `[D]` derived/proved, ⛔ refuted (kept in place), ⚡ land
> regardless of any ruling.
> ⚠ `scratch/` is UNTRACKED. Every load-bearing number from those memos is
> carried into this file so a `git clean` cannot destroy the finding.

---

# Part 0 — The decisions, stated first

**D0.1 — Build the catalog and the descent machinery. Refuse the general
orbit-space engine.** Computing an orbit space from scratch is mechanical
(invariant ring → syzygy ideal → Procesi–Schwarz PSD condition), but the groups
that occur in transport number about a dozen. A Gröbner engine is abstraction
without a consumer, and its failure mode is debugging elimination orderings
instead of transport. Each catalog entry is **derived once by the procedure**
(§III.1), then recorded with its derivation in the docstring and a symbolic
regression test.

**D0.2 — Symmetry reduction is a SECOND BINDING of the same kernel.** `Λ` is
representation-free; `.on(V)` supplies `(M, R)` from `V`'s frame; a quotient
space `V/H` is just another `V`. Nothing in the kernel layer changes and **no
operator is ever assembled-then-projected**.

**D0.3 — There is no `reduce(A, H)`.** Reduction is not an operation applied to
an assembled operator; it is a choice of which space to bind to, made at
Realization (stage 2). A `reduce(A, H)` function is the stage-inversion error
the operator-machinery report catalogues.

**D0.4 — Nothing declares what it can derive.** This is the discipline the whole
plan is built to satisfy, stated by the user 2026-08-31:
> *items carry their intrinsic properties, and concepts do not leak between
> objects or their combination; composability is checked exactly where the two
> objects meet; properties at the meeting place are either intrinsic or DERIVED
> from them — **if something has to be hard-coded, that is a signal that the
> individual objects are missing properties**; and this is applied
> retroactively, clean-before-extend.*

**D0.5 — Effort is a schedule, not a veto.** *"If something was more correct,
and we had to change 10 thousand sites, most likely we would do it anyway — so
effort is a measure of how many sessions something will take, not if it will be
done or not"* (user, 2026-08-31). Nothing in this plan is scoped out on size.

**D0.6 — What justifies the layer.** The user's plan listed three consumers and
said *"absent all three, do not build this."* There are now **four**, and the
fourth is remedial rather than prophylactic:

0. ⛔ **A live wrong-answer bug** — Part I. `solve_sn(scattering_order ≥ 2)` on
   any 1-D chart is wrong, publicly reachable, unguarded, with **zero** test
   witnesses. Gate G2 does not catch a hypothetical class; it catches a firing
   one. **This makes the justification unconditional.**
1. Curvilinear S_N — the quotient's connection term *is* the angular
   redistribution operator (`[R]`, §III.4).
2. The compatibility gates (Part IV) — a runtime check on quadrature/symmetry
   mismatch that catches a real error class no production code checks.
3. Isotypic decomposition — required by higher harmonics (Riesz contour) and
   reactor noise regardless; symmetry is where it comes from.

---

# Part I — The defect that forces this

**Reproducer:** `scratch/_pl_slab_defect_repro.py` (two self-checking legs).

## I.1 End-to-end, 1-group infinite medium `[M]`

Configuration: `Σ_t = 1.0`, `Σ_s0 = 0.5`, `Σ_s1 = 0.40`, `Σ_s2 = 0.45`;
`Quadrature.gauss_legendre(n_ordinates=8)`; 4 cells, all-reflective; uniform
isotropic unit source; `inner_solver="krylov"`, `inner_tol=1e-13`.
In an infinite medium the flux is isotropic, so every `φ_ℓ` with `ℓ ≥ 1` is
ZERO and `scattering_order` **cannot** matter. Analytic
`φ = W·Q/(Σ_t − Σ_s0) = 4.000000000000`.

| L | φ | rel err |
|---|---|---|
| 0 | +4.000000000000 | 1.2e-15 |
| 1 | +4.000000000000 | 4.4e-16 |
| **2** | **−3.764705882353** | **1.94** |
| **3** | **−3.764705882353** | **1.94** |

Wrong sign. 194 %.

## I.2 On real data `[M]`

421-group pure-Be slab (`BE009`, N = 0.1236 at/b-cm), 20 cm, 12 cells, vacuum
both faces, uniform U-235-χ source, S8 GL, Krylov: `sum(φ)` =
**2.149e3 / 1.978e3 / 2.267e4** at L = 0 / 1 / 2 — the L=2 answer is **10.6×**
the L=0 one.

## I.3 Reachable, unguarded, unwitnessed `[M]`

* Public: `solve_sn(..., scattering_order=2)`, `solve_sn_fixed_source(...)`.
* **Nothing validates `L` against the quadrature.** The sole clamp is
  `sn/solver.py:1357-1363`, against `min(len(m.SigS) − 1)` — the DATA's Legendre
  depth, not the rule's. The Q6 `exactness` / `degree_of_exactness` machinery is
  structurally disconnected: 0 occurrences in `frame.py`, `harmonic_frame.py`,
  `scattering.py`, `n2n.py`, `fission.py`.
* **No witness anywhere.** The corpus uses `scattering_order` = 0 (74 sites),
  1 (46), 3 (2) — and BOTH `L=3` sites
  (`tests/sn/operators/test_scattering_operator.py:1593`,
  `tests/sn/_fixtures/wave_t_t3/_capture_pre_t3_snapshots.py:134`) use
  `Quadrature.lebedev(17)` on a 2-D mesh, where it does not bite.

## I.4 Scope, on PRODUCTION's own bound basis `[M]`

`q.spherical_harmonics(L)` — the `ℓ ≥ 1` moment of an ISOTROPIC flux, which must
be 0:

| rule × L | live | rank | max offdiag | raw `ℓ≥1` of isotropic | verdict |
|---|---|---|---|---|---|
| `gauss_legendre(8)` L=1 | 2 | 2 | 8.8e-17 | 6.9e-18 | clean |
| **`gauss_legendre(8)` L=2** | 5 | **4** | 1.155 | **5.774e-1** | RANK-DEF SPURIOUS DENSE |
| **`gauss_legendre(8)` L=3** | 9 | **7** | 1.155 | **5.774e-1** | RANK-DEF SPURIOUS DENSE |
| **`gauss_legendre(16)` L=2,3** | 5, 9 | **4, 7** | 1.155 | **5.774e-1** | same |
| `level_symmetric(4)` L=3 | 16 | 16 | 1.9e-1 | 8.8e-17 | DENSE only |
| `folded_product(2,4)` L=3 | 7 | **4** | 2.793 | 1.1e-16 | RANK-DEF only |
| the other 18 rows | — | — | — | ≤ 1.2e-16 | clean |

**7 of 24 rows are not clean; 4 are the wrong-answer defect.** The spurious
value is exactly `1/√3 = 0.5773503` at every slab order tested (`n_ordinates` =
2, 4, 6, 8, 16) — a structural constant, not a quadrature-accuracy artifact.

The two DENSE-but-full-rank rows are a **separate class** — *degree L asked of a
rule too coarse to resolve it* — whose natural home is the Q6 exactness
machinery (Phase 4.4).

---

# Part II — Root cause, and the leak inventory

## II.1 The mechanism, read from the source `[M]`

`_evaluate_real_sh`, `orpheus/numerics/basis/spherical_harmonic_basis.py`, hard-
codes `ℓ = 1` in CARTESIAN form:

```python
Y[:, 1, 0] = mu_z ;  Y[:, 1, 1] = mu_x ;  Y[:, 1, 2] = mu_y
```

then for `ℓ ≥ 2` switches to the SPHERICAL route:

```python
cos_theta = mu_x
sin_theta = sqrt(1 - mu_x**2)      # [M] 0.279 .. 0.983 on GL8 -- NOT small
on_axis   = sin_theta < 1e-15      # [M] fires 0 of 8 -- never triggers
cos_phi   = mu_y / sin_theta       # 0 / 0.279 = 0
sin_phi   = mu_z / sin_theta       # 0
phi       = arctan2(0.0, 0.0)      # = 0.0  <-- FABRICATED
```

so `cos(m·φ) = 1` and **every `m > 0` slot is non-zero**.

⭐ **`ℓ = 1` is clean ONLY because it never takes this path** — which is exactly
the measured L boundary in §I.1. That coincidence is why the defect reads as an
"L ≥ 2 problem" rather than as what it is.

⭐ **The `on_axis` guard proves the case was anticipated and mis-aimed.** It
fires when `sin_theta ≈ 0`, i.e. `|μ| ≈ 1` — a direction ALONG the polar axis.
The slab's situation is the opposite: the TRANSVERSE part is zero while
`sin_theta` is large. `[M]` over 10 shipped rules the guard fires on **2 nodes
of 1 rule** (the Lebedev poles) while the case it misses is **100 % of the
ordinates on every slab rule**. ⚠ It is written with `np.where`, not `if`, so
**no coverage tool can ever report the branch unexercised.**

⭐⭐ And §V.5 below shows those same 2 nodes ARE the singular stratum of the
`SO(2)` action — so the guard is a *fixed-point detector* that already exists,
pointed at the wrong question.

## II.2 The convention mismatch, stated exactly

> `Quadrature` supplies `mu_y = mu_z = 0` meaning **"there is no azimuthal
> information"**. `_evaluate_real_sh` reads it as **"the azimuth is 0"**.
> Two meanings, one set of zeros.

AI-failure-mode #6 (convention drift) at the quadrature↔basis boundary. And
`Quadrature.spherical_harmonics`'s own docstring (`directional.py:538`) records
the FIRST meaning — *"only the m = 0 harmonics carry non-zero values; the other
slots are filled with zeros"* — while the code implements the second. **The
docstring states the invariant that would make the code correct.**
`[M]` that exact claim appears at **exactly 1 site tree-wide**.

## II.3 ⭐⭐ The redundant harmonic IS `det P` — a closed form, not a coincidence

The catalog derivation (§III.2) gives `det P = 4(1 − μ²)` as the **squared orbit
radius** of the `SO(2)` action. `[D]` the harmonic that makes the slab Gram
rank-deficient is that same polynomial:

```
Y₂^{+2} = 0.288675 · P₂²(μ) = 0.866025 · (1 − μ²)        [= det P, up to norm]
        = 0.577350 · (Y₀ − Y₂⁰)     since 1 − μ² = ⅔(1 − P₂)
```

⟹ predicted null direction `[−0.447214, 0, +0.447214, 0, +0.774597]` over the
live slots `{(0,0),(1,0),(2,0),(2,1),(2,2)}` — **identical to the SVD-measured
one** (`[M]`, agreement 4e-07, limited only by my hand-transcribed digits).

**So the rank deficiency is a theorem about the quotient, not a numerical
curiosity: the redundant slot is the orbit-radius polynomial, and on the
quotient the orbit radius is determined by μ.** This is a **third** appearance
of `(1 − μ²)` — orbit radius, angular-redistribution coefficient, and now the
redundant harmonic — and unlike the first two it is *proved*. It is live
evidence for Phase 1.3.

## II.4 The leak inventory

Applying D0.4 retroactively. Each row: the leak, its evidence, and the missing
property it signals.

| # | leak | evidence `[M]` | missing property |
|---|---|---|---|
| **L1** | `angular_frame` **fabricates** `support=SPACE_SPHERE` onto interval nodes | `directional.py:587-593`; rule support `[-1,1]` → frame support `S^2`; `‖node‖` 0.1834..0.9603 | `Basis` cannot declare its **domain**, so the binder forges a matching support |
| **L2** | `Basis` has **no `domain`** | `basis/base.py` declares `evaluate / synthesize / analyze / analyze_transpose / reconstruct / reconstruct_transpose / gram_structure / mass_matrix / space` — no domain. Its own docstring: *"the nodal/domain space comes from the measure"* | the intrinsic property itself |
| **L3** | no composability check where the two meet | `GalerkinFrame(basis, measure)` validates neither support nor domain | the check has no home even if L2 existed |
| **L4** | a 1-D rule exposes a 3-D direction tuple with fabricated zeros | see §II.5 — **ONE** production consumer, and it is L1 | the **chart map is not an object** (L6) |
| **L5** | `_harmonic_basis` selects the basis from `folded_by` alone | `directional.py:505-530`; a slab GL rule has `folded_by is None` | a NATIVE 1-D rule was never "folded" from anything. The question is the measure's **symmetry**, not its construction history |
| **L6** | `pushforward(φ, *, new_space)` takes the codomain as an ARGUMENT | `measure.py:740-760` | **the chart is an anonymous `Callable`** with no declared domain→codomain |
| **L7** | `product()` passes `support=SPACE_SPHERE` as a literal | `rules_product.py:526-534` — and it sits **between two fields whose own comments insist they are derived**: `invariance_group` *"COMPUTED from the factors … never a declared literal"*, `exactness` *"DERIVED from the two factors' own claims"* | same as L6; the constructor imports `gauss_legendre_on_mu` and then **discards the factor structure** |
| **L8** | `Space = str`, *"recommendations, not constraints"* | `measure.py:110-119` | nothing can enforce a declared domain today |

⭐ **The codebase already learned L1/L7's lesson one level down.**
`measure.py:120-128` records that `SPACE_CIRCLE` was `"[0,2π)"` until
2026-08-02, which *"made the tag assert a coordinate"*, and was renamed to name
the **manifold**. Writing `S^2` onto interval nodes is the same defect one level
up: **the tag asserting a membership the nodes do not have.**

## II.5 ⭐ The phantom census — ONE production consumer `[M]`

Probe: `scratch/_phantom_axis_probe.py`. In-process, make `axis_cosines(i)`
RAISE for `i >= dim`, then drive production 1-D paths.

```
[patched] axis_cosines(i >= dim) RAISES:
   PHANTOM    PRODUCTION slab solve, scattering_order=0   <- directional.py:589
   PHANTOM    PRODUCTION slab solve, scattering_order=1   <- directional.py:589
   survives   3-D control: lebedev(11).mu_z
   survives   3-D control: product(4,8).angular_frame(2)
```

`directional.py:589` is `column_stack([axis_cosines(0), axis_cosines(1),
axis_cosines(2)])` inside `angular_frame` — **the same expression that forges
`support=SPACE_SPHERE`.** Both 3-D controls survive, so the patch is correctly
scoped.

⟹ **The phantom components exist SOLELY to feed the fabricated S² measure.
L1 and L4 are one fix, not two** — which collapses an unusable
"171 occurrences of `mu_y`/`mu_z`/`axis_cosines` over 20 files" into one
expression.

⚠ Note it fires at **L = 0** too: the fabricated measure is built on *every*
slab solve. It is harmless at L=0 (no `m ≠ 0` slots exist) and at L=1 (the
Cartesian special case makes them genuinely zero), and harmful from L=2. **The
fabrication is universal; the symptom is not.**

⚠ **Denominator, stated honestly:** the probe covers 2 production solve paths,
4 direct accessor cases, 2 controls. It is **not** the whole tree — curvilinear
/ sphere paths, `geometry/boundary/_specular.py:162`, `sn/angular/closure.py`
and the MMS derivations also read `mu_y`/`mu_z` and were not exercised. Running
the same 20-line probe under the full suite is the honest census and is
Phase 2.4's first act.

## II.6 The measure machinery is already clean `[M]`

`gauss_legendre(8).measure.nodes.shape == (8,)` — **scalars**; `dim = 1`;
`support = "[-1,1]"`. **The measure never had three dimensions.** The phantoms
live entirely in the `Quadrature` façade. This is the single most encouraging
fact in the plan: the space layer is right, and the fabrication is a thin,
locatable crust on top of it.


## II.7 ⭐⭐ The tree ALREADY SHIPS this ontology — derived, tested, and production-unreachable

`[M]` 2026-08-31, at execution. Before deriving the `S²/SO(2)` catalog entry,
I looked for it. It exists, in `orpheus/numerics/quadrature/registry.py`:

```python
class AngularSymmetry:
    continuous_isotropy: SubgroupOfO3          # G⁰, the spent continuous part

    @property
    def support(self):                          # ← "derived, not declared"
        if spent == SubgroupOfO3.SO2:
            # S²/SO(2) — the polar marginal. The orbits of the axial
            # rotation are the constant-μ circles, so the quotient is
            # parameterised by μ alone.
            return SPACE_INTERVAL_M11
        if spent == SubgroupOfO3.Trivial:
            return SPACE_SPHERE

    @property
    def reference(self):                        # ← also derived
        return LEGENDRE if spent == SO2 else UNIFORM_ON_SPHERE

    def admits_domain(self, measure):           # ← Stage 0
        return measure.support == self.support
```

`[M]` the shipped table, printed:

| geometry | `G⁰` | derived `support` | derived `reference` |
|---|---|---|---|
| slab | `SO2` | `[-1,1]` | legendre |
| sphere | `SO2` | `[-1,1]` | legendre |
| cylinder | `Trivial` | `S^2` | uniform-on-sphere |
| cartesian2d | `Trivial` | `S^2` | uniform-on-sphere |

⟹ **Three consequences, and each shortens a phase.**

1. **Phase 1.1's `SO(2)/S²` entry is largely AUTHORED ALREADY** — the quotient,
   its chart (`μ` alone), its orbit description (constant-μ circles) and its
   exactness reference (Legendre) are written, with a docstring that already
   insists they are *derived, not declared*. 1.1 becomes *reconcile and extend*,
   not *derive*. ⚠ It is NOT complete: no `det P`, no named singular stratum,
   no action on the axis index set — those are still owed.
2. **Phase 2.2's G0 predicate EXISTS** as `admits_domain`. 2.2 is a WIRING
   problem, not a design problem.
3. ⛔ **And it is unreachable from production.** `[M]` `select_quadrature`
   (the only caller of `admits_domain` / `admits_symmetry`) has **21 callers,
   every one a test** — `grep -rn "select_quadrature" orpheus/` returns 12 hits
   and **all 12 are docstrings, the definition, or `__all__`**. The registry's
   own prose calls it *"opt-in convenience"*. So this is the
   `[[project-test-dependence-nexus-config]]` lesson verbatim: **authored
   knowledge may be INERT — look for it before deriving it.**

## II.8 ⭐ The rule tells the TRUTH; only the frame forges — and it falsifies TWO tags

`[M]` printed across 12 shipped rules:

| rule | `measure.support` | `angular_frame(2).measure_space` |
|---|---|---|
| `gauss_legendre(2,8)` | **`[-1,1]`** | `L2[S^2]` ⛔ |
| `level_symmetric(4,8,12)` | `S^2` | `L2[S^2]` ✅ |
| `lebedev(11,17)` | `S^2` | `L2[S^2]` ✅ |
| `product(4,4),(4,6),(8,8)` | `S^2` | `L2[S^2]` ✅ |
| `folded_product(2,4),(4,8)` | **`S^2/sigma_y`** | `L2[S^2]` ⛔ |

⭐ **The measure layer was never wrong.** `gauss_legendre(8).measure.support`
already reads `'[-1,1]'` — agreeing with `AngularSymmetry['slab'].support`
independently. §II.6 said the space layer is right; this says the *tag* layer
is right too. The forgery is **one expression**, in `angular_frame`.

⭐ **NEW leak, absent from §II.4's inventory: the frame destroys the FOLD's
quotient tag too.** `folded_product` declares `support='S^2/sigma_y'` — the
measure layer already speaks quotients — and `angular_frame` overwrites it with
`S^2`. So L1 falsifies **two** truths: the slab's `[-1,1]` AND the fold's
`S^2/sigma_y`. ⚠ This one is *behaviour-relevant*, not merely cosmetic, and `[M]` the
mechanism is now named exactly. Fixing it moves `frame.measure_space` from
`L2[S^2]` to `L2[S^2/sigma_y]`; `[M]` `measure_space` has **10 occurrences in
`orpheus/`, of which 3 are real production reads** — the property itself
(`frame.py:328`) and two OPERATOR FACES: `_FrameAnalysis.domain`
(`frame.py:582`) and `_FrameReconstruction.codomain` (`frame.py:627`). So the
tag is a **typed arrow's domain**, and `OperatorSum` refuses unequal domains
while `_AdjointOperator` branches on them — `plan-authoring` §8's worked example
verbatim (*"a field a consumer BRANCHES on is an input, not metadata"*).
⟹ **0.1c is a live behaviour change on any composition mixing a folded-rule
angular frame with an `S^2`-tagged operand. Measure it with the suite before
landing; it is NOT part of the bit-identical half, and it owes a gate at the
tier the change is observable.**

`[M]` **10 of 12 rows: the forged nodes are BIT-IDENTICAL to the rule's own**
(`np.array_equal`, exact shapes). The two exceptions are the 1-D rules, where
the rule's nodes are scalars `(N,)` and the forgery is `(N,3)`.

## II.9 ⛔ Why 0.1 cannot fully land without 3.4 — the lift does not exist

`[R]`, and it is the reason the coupling is essential rather than incidental.
`angular_frame` is trying to map `[-1,1] → S²`. **There is no such map**: a
point of `[-1,1]` is an ORBIT of the `SO(2)` action, not a point of `S²`, and
choosing a representative is exactly the "fabricate an azimuth" move
(§II.2). The arrow that exists runs the other way — the quotient `S² → [-1,1]`
— which is why the repair must change the **basis** (Phase 3.4's trivial
isotypic component) and cannot be done by fixing a tag.

⟹ **Phase 0.1 therefore SPLITS**, and the plan did not say so (it flagged this
obligation on 0.6 only):
* **0.1a — landable now, bit-identical.** A rule whose measure already lives on
  `S²` hands the frame **its own measure**; the `column_stack` + literal
  `support=SPACE_SPHERE` is deleted for those rows. `[M]` 10 of 12 shipped rules,
  byte-identical nodes; and it removes the phantom read the census attributes to
  `directional.py:589`.
* **0.1b — rides Phase 3.4.** The 1-D rows, where there is nothing honest to
  build until the basis changes. Until then the construction stays, but named,
  documented, and tagged ERR-080 with its retirement trigger — the
  `coding-standards` transitional-violation idiom, not an anonymous literal.
* **0.1c — the fold's tag**, gated on the §8 blast-radius measurement above.



## II.10 ⛔⛔ The ontology has a slot for the `G⁰` quotient and NONE for a `Γ` quotient — and the forgery is what HIDES that

`[M]` 2026-08-31, verified by me after an `explorer` survey raised it.
`AngularSymmetry.support` derives from `continuous_isotropy` **alone**. A rule
may ALSO have been quotiented by a *discrete* subgroup, and there is nowhere to
say so:

```
folded_product(4,8).measure.support = 'S^2/sigma_y'
  admits_domain(slab)        = False   (requires '[-1,1]')
  admits_domain(sphere)      = False   (requires '[-1,1]')
  admits_domain(cylinder)    = False   (requires 'S^2')      ← the shipped pairing
  admits_domain(cartesian2d) = False   (requires 'S^2')
```

⚠ **`folded_product` IS the shipped cylinder configuration.** `[M]` two
production error messages *instruct the user to build it*
(`sn/angular/closure.py:1052`, `:1228`); two MMS defaults use it
(`derivations/continuous/mms/sn.py:2104`, `:3873`); **51 test files** pair it
with `CoordSystem.CYLINDRICAL`. `[M]` over 20 (constructor × geometry) pairs,
**5 are production builds that stage 0 rejects.**

⟹ **Phase 2.2 cannot be "just wiring" after all** (this REFINES §II.7's
conclusion, which was measured before I had this): `admits_domain` is a correct
predicate over an **incomplete domain vocabulary**. Wiring it as-is would refuse
the cylinder. The ontology owes a second slot — the discrete quotient the rule
took — and `SubgroupOfO3` already exists to name it.

⭐⭐ **And here is the part that makes this a Phase-0 argument, not just a
Phase-2 one: the forgery is what CONCEALS the gap.** By overwriting
`S^2/sigma_y` with `S^2`, `angular_frame` makes the folded rule *look* like it
lives on the sphere — so nothing downstream ever presented `admits_domain` with
a quotient tag, and the missing slot could not be discovered. **0.1c does not
merely stop a lie; it surfaces the ontological gap the lie was hiding.** That is
the strongest available argument for landing 0.1c early rather than bundling it,
and it is why the user's ruling (*measure it, then land it in Phase 0*) is the
right call.

## II.11 ⛔ Collateral defect lead — `orbit_certificate` refuses 1-D by SHAPE, and says something false when it does

`[M]` verified: `SubgroupOfO3.Mirror("x").is_invariant(gauss_legendre(8).measure)`
is **`True`**, while `orbit_certificate(measure, Mirror("x"))` returns **`None`**
— refused by a shape test at `numerics/symmetry.py:1573-1574`, *before*
invariance is ever asked. Both refusal messages (`measure.py:1010`,
`symmetry.py:1612`) read *"not σ_x-invariant (or σ_x is continuous)"* —
**both halves are false for this input**, and the docstring names 2 of the 3
arms, omitting the operative one. `vv-principles` #17's multi-arm granularity
trap, in a production guard. **Not on the exit path** — file it, do not fix it
here.

## II.12 The nearest geometry-carrying frame, and the join that does not exist

`[M]` the production chain (`explorer`, grep — the graph could not answer this:
`callers(angular_frame)` returned `total: 0` with 25 unresolved spellings):

```
solve_sn:2531 → _as_sn_mesh:137 → SNMesh:171 → augmented_mesh.py:1186
  → full_field_space:1249 → solver.py:1413 ScatteringOperator.from_solver_data
  → scattering.py:811 HarmonicFrame.for_space → harmonic_frame.py:413,416
  → directional.py:592   ← the forge
```

⭐ `solve_sn` does not CHOOSE a quadrature — it is handed one (required
positional, no default), which is why `select_quadrature` is unreachable: **it
was never on the path.** The nearest frame carrying a geometry is
`SNMesh.__init__`'s `match self.coord` (`augmented_mesh.py:325`), which already
runs a (quadrature, coord) admission gate 31 lines later
(`assert_carrying_quadrature`, `:356`). The tag is dropped at `:1186`.
⚠ Two obstacles, both `[M]`: the tag is a `CoordSystem` (3 members) while
`GEOMETRY_ANGULAR_SYMMETRY` is keyed by 4 strings — **the `(coord, ndim) → key`
join exists nowhere in the tree**; and that precedent gate has exactly **one**
call site, inside `case CYLINDRICAL` (the `plan-authoring` §2 gate-as-denominator
trap, already logged for this very symbol).

## II.13 Phase 1.1 is RECONCILE — the surfaces that already state the ontology

`[M]` beyond §II.7's `registry.py`: `numerics/quadrature/rules_1d.py:1,:10-14,:127-129`
(module docstring naming `AngularSymmetry.continuous_isotropy` as the home and
stating `S²/SO(2) = [-1,1]`); a **live production guard** at
`numerics/generating_measure.py:398-402`; `numerics/symmetry.py:212,:805`; and
four authored sections in `docs/theory/foundations/discrete_measures.rst`
(`:739, :1114, :1236, :1308`), including a worked stage-0/1 table at
`:1205-1213`; plus 5 test modules.
⭐ And `MirrorEvenSphericalHarmonicBasis` (`spherical_harmonic_basis.py:401-426`)
**already states this campaign's diagnosis and DERIVES its parity mask by probe
rather than hand-listing it** — a genuine precedent for 3.4's mechanism. ⚠ But
its KEY cannot be borrowed: `folded_by` is a group object recording a *discrete*
fold, structurally unavailable for a *continuous* quotient. **The tree holds two
quotient records — a group object (wired, discrete) and a string tag (inert,
continuous) — and the slab is the unwired half.**


---

# Part III — The theoretical spine

## III.1 The derivation procedure (catalog generator, run by hand, once per entry)

For `G` compact acting orthogonally on `ℝⁿ`, `X ⊆ ℝⁿ` a `G`-stable real
algebraic set:

1. Minimal generators `p₁,…,p_k` of the invariant ring `ℝ[x]^G` (finite by
   Hilbert–Weyl).
2. `p : ℝⁿ → ℝᵏ` is proper and separates orbits ⟹ `ℝⁿ/G ≅ p(ℝⁿ)`, closed
   semialgebraic.
3. Syzygy ideal `I = ker(ℝ[y] → ℝ[x], yᵢ ↦ pᵢ)` by elimination.
4. **Procesi–Schwarz** (Invent. Math. 81, 1985): `p(ℝⁿ) = { y ∈ V(I) : P(y) ⪰ 0 }`
   with `P_ij = ⟨∇pᵢ, ∇pⱼ⟩` re-expressed in `y`. `V(I)` gives equalities;
   `P ⪰ 0` gives inequalities.
5. Intersect with the ideal of `X`.

## III.2 Worked entry — `S²/SO(2)`, and its three consequences

Invariants `p₁ = z`, `p₂ = x² + y²`; empty syzygy ideal;

```
P = [[1, 0], [0, 4p₂]],    det P = 4p₂
```

so `ℝ³/SO(2) = {p₂ ≥ 0}`. Adjoining `p₁² + p₂ = 1`:

```
S²/SO(2) = { μ ∈ ℝ : 1 − μ² ≥ 0 } = [−1, 1],    μ = Ω̂ · ẑ
```

Three facts, all load-bearing downstream:

- **The action is not free.** The poles have full stabilizer. `[−1,1]` is a
  manifold-with-boundary / orbifold, not a quotient manifold; `μ = ±1` is the
  **singular stratum**. Any design assuming quotients are smooth submersions is
  wrong exactly there.
- **`det P = 4(1 − μ²)` is the squared orbit radius**, and the orbit-space
  boundary is where it vanishes. `[M]` this is the same polynomial as the
  curvilinear angular-redistribution coefficient in `(1/r)∂_μ[(1−μ²)ψ]`, the
  same locus as the vanishing characteristic trace measure, **and** (§II.3) the
  redundant harmonic `Y₂^{+2}`.
- **`(SO(3), SO(2))` is a Gelfand pair**, and `S² = SO(3)/SO(2)`, so the object
  is really the double coset space `SO(2)\SO(3)/SO(2)`. Bi-invariant functions
  on it are the zonal spherical functions = Legendre polynomials.
  **Funk–Hecke is Schur's lemma on a Gelfand pair.** `RΛM` with `Λ` diagonal is
  a theorem of harmonic analysis on the quotient, not a chosen factorization.

## III.3 Encode the contravariant side, not the point set

`C(X/G) ≅ C(X)^G` and `L²(X/G, π_*μ) ≅ L²(X,μ)^G`. ORPHEUS never touches orbit
spaces; it touches function spaces over them. Reasons, in order of force:

1. Bases and quadratures live on function spaces, not orbit spaces.
2. The invariant side is finitely presented and computable.
3. **It is already in the taxonomy.** For `G` compact acting unitarily, the Haar
   average `ℛf = ∫_G f(g⁻¹·) dg` satisfies `ℛ² = ℛ` and `ℛ* = ℛ`: the orthogonal
   projection onto the invariant subspace. It factors as `ℛ = π* ∘ 𝔄` with
   `𝔄 π* = id`. **`(π*, 𝔄)` is a retraction/section pair.**

**Consequence — operator descent is retract conjugation:** `A_H = 𝔄 ∘ A ∘ π*`,
equal to the restriction of `A` to the invariant subspace **iff
`[A, U_g] = 0 ∀ g ∈ H`**. That is a per-operator condition, checked at
Declaration (gate G3), and it is where physics enters: geometry, material field
and BC must each be `H`-equivariant.

**Quotienting is an AXIS-level operation, not a space-kind-level one.** The
space kinds are untouched; the one genuinely new piece of axis data is the
stratification descriptor (§III.7).

⚠ **Vocabulary reconciliation, measured against THIS tree** (`plan-authoring`
§7.2 — the source plan was written against a target architecture):

| the source plan says | `[M]` this tree |
|---|---|
| `numerics/symmetry/groups.py` (a package) | `orpheus/numerics/symmetry.py` — a **module**. Phase 3.1 implies a re-home; `plan-authoring` §6d applies (check the import edge in BOTH directions) |
| *"it is a `Retract`; the lattice already has the structure"* | **no `Retract` type.** But the machinery EXISTS: `FunctionSpace.retraction()` / `.section()` (`space.py:452/487`) + `class AxisRetractionOperator` (`operator.py:3362`), a "split epi/mono pair" from CS4b S6.0b. ⟹ **the structural claim holds; the gate's vocabulary must be re-pointed** |
| Phase gate *"no new `Relation` subclass added"* | **no `Relation` type** — the gate is unrunnable as written |
| the per-axis `WithTrace` descriptor map | **0 hits in `orpheus/`** — forward-referenced, not yet built |

## III.4 `ℛ* = ℛ` is another instance of "adjoint = dual iff tight"

`ℛ² = ℛ` is algebraic and always holds. **Self-adjointness does not.** `ℛ* = ℛ`
requires the Gram to be `H`-invariant — the group must act by isometries of the
*discrete* inner product, not merely of the continuum measure. On a non-tight
frame the adjoint and the dual diverge and `ℛ` becomes an oblique projection in
the adjoint pairing. ⟹ the tightness gate and the symmetry gate are **the same
gate on a different operator**; share the implementation.

**Tightness descends for free.** If `{fᵢ}` is a tight frame with bound `A` and
`P` is an orthogonal projection, `{Pfᵢ}` is tight for `PV` with the same bound.
With `P = ℛ`: the `B/A` ratio is **inherited, not recomputed**.

## III.5 ⭐ The well-posedness condition for a per-mode operator

From `scratch/n2n_pl_frames_attack.md`, and it is the sharpest statement in the
plan:

> **`reconstruct ∘ multiply ∘ analyze` is well-defined iff `ker(synthesis)` is
> invariant under the multiplier** — a property of the (basis, measure,
> operator) **TRIPLE**, not of any pair.

If synthesis has a non-trivial kernel and `Λ` does not preserve it, two moment
vectors representing the SAME angular function produce DIFFERENT sources: **the
operator is not a function of the angular flux at all, only of the coefficient
representative.** A well-posedness failure, not an accuracy failure.

⟹ Corollary that kills a whole family of proposed repairs: **no left-, right-,
or pseudo-inverse can restore representation-independence that the retained
function space does not have.**

`[M]` and the practical discriminator: over the 8 shipped test frames,
**rank-over-LIVE-slots is a perfect detector** — slab-GL8-L2 is the only
rank-deficient row, while `product(4,4)-L2` and `LS4-L3` are DENSE for
legitimate reasons and full rank. ⟹ **rank is an excellent DETECTOR and a bad
REFUSAL predicate** (it over-rejects padding and harmless intra-block
degeneracy); a guard keyed on DENSE-ness would fire on healthy frames.

## III.6 What the literature settles

Full memo: `scratch/n2n_pl_consistency_literature.md`.

**The name.** *Galerkin quadrature* — **Morel, J. E. (1989), NSE 101(1), 72–87**,
local at `scratch/literature/`. The moment-to-discrete / discrete-to-moment
matrices `M`, `D`.
⚠ Morel's primitive is `D` and he writes `M = D⁻¹` (Eq. 27); the modern
convention makes `M` primitive and writes `D = M⁻¹`. Same statement — but quoted
side by side they read as a contradiction.

**The 1-D answer is a DERIVATION, not a convention.** Bell & Glasstone (1970)
§2.6a, printed p. 103: in the addition theorem the `cos m(φ − φ′)` terms
*"vanish upon integration over φ′"*, leaving Eq. (2.79), a **pure Legendre**
slab equation. §5.3a, printed p. 227: spherical harmonics appear in general
geometry *"because there is no azimuthal symmetry to eliminate terms from the
addition theorem of Legendre polynomials"*. Morel's 1-D `D`/`M` (Eqs. 39a/39b)
carry **no `m` index at all**; the `(ℓ,k)` double index first appears for the
CYLINDER (Eq. 45) with an `ℓ + k` EVEN restriction.

**The crux — invertibility is the PREMISE, not the output.** Larsen & Morel
(2010) Eqs. (1.104)/(1.106) declare `D`, `M` to be `N × N` **by definition**,
justified by unisolvence: *"One can uniquely define a polynomial of degree N−1
either in terms of N Legendre moments or in terms of N discrete function values
at N distinct points. Therefore D and M should be inverses of one another."*

**What the literature does instead.** Morel §VI: 2-D `S₂` has 4 directions but a
`P₁` source has 3 harmonics, so `M` is `4 × 3`. The remedy is not a
pseudo-inverse, not least-squares, not regularisation: *"one must find an
additional higher order spherical harmonic function that will provide a
moment-to-discrete matrix that can be inverted."* ⟹ **change the moment set
until `M` is square and invertible.** ORPHEUS's slab has too many AND a
dependent set, so the mirror move is to remove.

⚠ **An honest negative.** **UNFOUND**: no source states that 3-D real spherical
harmonics become linearly dependent on a 1-D ordinate set. The literature does
not pose the question **because no established formulation ever puts them
there.** The absence is the finding.

## III.7 The singular stratum

`μ = ±1` is the fixed-point set of the `SO(2)` action; there the orbit collapses
to a point and the pushforward trace measure vanishes. This is the
characteristic angular boundary already flagged as the forcing consumer for a
per-axis descriptor map. **The two are the same object**, and the degeneracy is
*per-axis and per-point, computed from the group* — not a global property of the
space. `[M]` §II.1's `on_axis` guard already detects exactly this stratum.

## III.8 Isotypic decomposition and the eigenvalue payoff

If `[A, U_g] = 0`, Schur block-diagonalizes `A` over the irreps of `H`; the
quotient space is the **trivial isotypic component**.

⭐ **This gives the correct NAME for the repair.** The sub-basis that fixes the
slab is not a "zonal special case" — it is the trivial isotypic component of the
`SO(2)` action, and that says *why those slots and not others*.
`[M]` probing invariance under SO(2)-about-x, using the same probe idiom
`MirrorEvenSphericalHarmonicBasis` uses, over **6 INCOMMENSURATE angles × 9
generic random directions** at `L = 4`: the invariant set is **exactly**
`{(0,0),(1,0),(2,0),(3,0),(4,0)}` — 5 of 25 real slots, **zero false positives,
zero false negatives**.
⚠ Incommensurate deliberately: a sample of right angles generates `C₄`, not
`SO(2)` (`vv-principles` #13's generator-set trap, which `MirrorEven` avoids
only because a mirror has order 2).

**Theorem (high confidence).** The `k`-eigenvalue fundamental mode lies in the
trivial isotypic component. Krein–Rutman gives a unique non-negative dominant
eigenvector `ψ`; `ℛψ` is non-negative, non-zero, and an eigenvector with the
same eigenvalue; uniqueness forces `ℛψ = ψ`. **So `k` may be computed entirely
on the quotient** — a theorem, not a modelling assumption. The same argument
applies to the dominant `α` mode when the prompt semigroup is positive
(*moderate* — needs the Perron-root argument for discrete `α`).

**Higher harmonics live in the nontrivial isotypic blocks**, giving the Riesz-
projection work an exact labelling: azimuthal modes indexed by `SO(2)` irrep.
**Ground field interacts** — over `ℝ` the nontrivial `SO(2)` irreps are
2-dimensional rotation blocks, over `ℂ` they are 1-dimensional (`e^{imφ}`), so
the isotypic decomposition is clean only over `ℂ`. The noise application and the
isotypic layer land together.

## III.9 ⭐ The support algebra of composed measures

The third leg, contributed this session. A measure's support is not a label to
be asserted; it is **derived by the verb that built the measure**:

| verb | support rule | `[M]` status |
|---|---|---|
| tensor `μ ⊗ ν` | `f"{μ.support} × {ν.support}"` | ✅ **derived** — `measure.py:588` |
| quotient `μ/G` | `f"{μ.support}/{G.name}"` | ✅ **derived** |
| pushforward `φ_*μ` | `new_space` — an **argument** | ⛔ **asserted** (L6) |

Two of three already work. **Only the pushforward leaks**, because `φ` is an
anonymous `Callable` carrying no codomain. The honest chain for a product rule:

```
GL[-1,1]  ⊗  RoU[S^1]   →  "[-1,1] × S^1"      (derived today)
   --pushforward(Archimedes chart)-->  "S^2"    (asserted today)
   --quotient(SO(2))-->                "[-1,1]" (the slab rule's TRUE provenance)
```

and the Archimedes cylindrical projection is not an implementation detail — **it
is the theorem that licenses product quadrature at all**, since `dμ dφ` pushes
forward to the uniform measure on `S²`.

⟹ **The missing property is that the chart map is not an object.** Give it a
type with a declared **domain → codomain** and:
* `support(φ_*μ) = codomain(φ)` becomes derived;
* `SPACE_SPHERE` stops being spellable as a literal in `product()` (L7);
* `Quadrature.dim`, `measure.support` and the node arity all become consequences
  of the construction rather than parallel assertions.

⭐ **And the same type kills the phantoms by refusing to exist.** `μ ↦ (μ,0,0)`
is not a chart into `S²` — `[M]` `‖(μ,0,0)‖ = |μ| ≠ 1`. The real relation runs
the OTHER way (`S² → [-1,1]` is the quotient), and a section of it does not
exist as a map because the fibre is a circle. So there is nothing legitimate for
`mu_y`/`mu_z` to return on a 1-D rule, and what a slab consumer actually needs
is not a direction but the streaming coefficient `Ω·n̂ₓ = μ` — which **is** the
node.

---

# Part IV — The well-posedness gates

Four distinct conditions on distinct objects. **They must not be fused into one
flag.** G1–G3 are the source plan's; **G0 is added this session** — G1–G3 are
all *group* conditions and none of them catches a category error.

| Gate | Statement | Object | Stage | Failure when violated |
|---|---|---|---|---|
| **G0 category** | `measure.support` lies in `basis.domain` | basis ⊗ measure | 2 | an angular basis bound to an energy measure typechecks; the binder forges a support (L1) |
| **G1 measure** | `H ≤ Isom(axis metric)` and `μ_axis` is `H`-invariant | axis | 2 | `π_*` ill-defined; `ℛ` not self-adjoint (§III.4) |
| **G2 discretization** | `H ≤ Sym(Q)` **and** the retained basis set is `H`-stable | frame / quadrature | 2 | discretization and quotienting do not commute; residue = **ray effects**. ⛔ **and the slab P_L defect of Part I** |
| **G3 commutation** | `[A, U_g] = 0 ∀ g ∈ H`, per operator term | operator | 3 | `𝔄 A π*` is not the restriction; silently wrong answer |

**G0 and G2 are complements, not duplicates.** G0 compares *manifolds* and
catches "wrong space entirely"; G2 compares *groups* and catches "right
manifold, wrong invariance". Neither sees the other's failure.

⚠ **G0 cannot be a tag-equality test.** The slab measure honestly says
`[-1,1]`; a trivial-isotypic basis honestly lives on `S²/SO(2)`; **those name
the same manifold.** String equality would refuse a legal pairing, and an alias
table would be exactly the hard-code D0.4 forbids. ⟹ **G0 compares the AMBIENT
manifold; G2 carries the quotient.** The principled predicate is:

> a function may be evaluated on a quotient **iff** it is invariant under the
> group being quotiented out ⟹ `measure.quotient_group ⊆ basis.invariance_group`

`[M]` and the lattice that answers it already ships —
`SubgroupOfO3.is_subgroup_of` gives exactly the right verdicts:

| basis invariance | measure quotient | verdict | meaning |
|---|---|---|---|
| `Trivial` (full SH) | `SO2` (slab) | **False** | ⛔ refuses the Part I bug, categorically |
| `SO2` (trivial isotypic) | `SO2` (slab) | **True** | ✅ the repair |
| `SO2` | `Trivial` (sphere) | **True** | ✅ a smaller space on a full rule is legal |
| `Mirror(y)` | `Mirror(y)` | **True** | ✅ the shipped fold |

**G2 catches a second, independent error class the source plan names**, and it
is worth shipping for that alone: level-symmetric `LQ_N` sets carry octahedral
symmetry, and `C₆` about a hexagonal axis is **not** a subgroup of `O_h`. A
1/6-core hexagonal reduction with a standard `LQ_N` set is not well-posed and
the residue is ray effects aligned with the mismatch. **If this layer ships
nothing else, ship G2.**

**G3's group is COMPUTED, not asserted.** `H_max = Sym(geometry) ∩ Sym(material
field) ∩ Sym(BC)`, an intersection over the catalog evaluated at setup. The user
does not declare "this is a quarter-core problem"; the code computes the largest
admissible `H` and reports it beside the acyclicity certificate. A user-supplied
`H` is validated against `H_max`, **never trusted**.

⚠ **Two measured obstacles before any gate can be armed:**
1. `[M]` `SubgroupOfO3.SO2.is_invariant(gauss_legendre(8).measure)` returns
   **`False`**. `[R]` a predicted axis-convention collision — `_embedded_nodes`
   uses **column 0** as the polar axis while the `SO2`/`Dinfh` realizations use
   **z** (`symmetry.py:800-819` vs `:889-899`, `:1014`). The measured `False` is
   fact; the "should be `True`" is reasoned and **unverified**. Phase 2.4.
2. `[M]` the rules populate their symmetry fields **inconsistently and with two
   different concepts**:
   `gauss_legendre(8)`: `invariance_group=Mirror('x')`, `folded_by=None`;
   `folded_product(4,8)`: `invariance_group=None`, `folded_by=Mirror('y')`.
   `invariance_group` = the symmetry of the NODE SET **within** its support;
   `folded_by` = the QUOTIENT group. **The slab's quotient group is recorded
   nowhere.** Anything built on these must not conflate them.

---

# Part V — Placement and the layer stack

## V.1 Pipeline

```
0  DATA          cross sections, geometry        → Sym(materials), Sym(geometry) computed here
1  KERNEL        Λ, representation-free          → UNCHANGED. Kernels do not know about H.
2  REALIZATION   kernel + frame + measure → arrow → QUOTIENT LIVES HERE.
                                                    OrbitAxis minted; G0/G1/G2 evaluated;
                                                    retraction/section pair minted;
                                                    stratification descriptor set
3  DECLARATION   A_loss = L + C − S − B          → G3 evaluated per term
4  POSING        Pencil(A, Π)                    → isotypic block selection
5  STRATEGY      splitting recipe                → computed on the quotient digraph
6  INSTANTIATION recipe applied to A − σΠ        → unchanged
7  SOLUTION      iterate / invert / normalize    → unchanged
```

**Ruling (D0.3 restated where it will be read): there is no "symmetry reduction"
step.** It is a choice of which space to bind to, made at stage 2.

## V.2 Layer stack (import direction preserved)

```
L1  numerics/symmetry/     group catalog, Haar/Reynolds, orbit charts,
                           pushforward measures, retraction/section minting.
                           Knows nothing of transport.
    numerics/basis/        H-stable sub-basis selection (Y_ℓm → Y_ℓ0 under SO(2))
    numerics/chart.py      NEW (§III.9): the typed chart, domain → codomain
L2  transport/symmetry.py  Sym(geometry), Sym(materials), Sym(BC), H_max
L3  sn/, moc/, cp/         method-specific G2 checks (quadrature symmetry groups)
```

`numerics/symmetry` must not import from `transport/`. "This reactor problem has
symmetry `H`" is L2; "`SO(2)` acting on `S²` yields `[−1,1]` with `dμ/2`" is L1.

## V.3 Naming (each name earned by a distinguishing invariant)

| Name | Invariant that earns it | Rejected alternative |
|---|---|---|
| `IsometryGroup` | acts by isometries of the axis metric — the precise G1 condition | `SymmetryGroup` (states no invariant) |
| `OrbitAxis` | its index set **is** the orbit space | `QuotientAxis` (names the construction) |
| `ReynoldsProjection` | idempotent + self-adjoint under the invariant Gram | `Symmetrizer` (verb, no invariant) |
| `IsotypicDecomposition` | Schur block structure indexed by irrep | `SymmetryBlocks` (collides with partition/block vocabulary) |
| `Chart` | a map carrying its own **domain → codomain**; §III.9 | `Embedding` (asserts injectivity we do not always have) |
| *(no new relational type)* | `(π*, 𝔄)` re-uses `AxisRetractionOperator` + `FunctionSpace.section` | a new `QuotientMap` — would duplicate CS4b S6.0b |

⚠ ⛔ **Do NOT name the sub-basis `ZonalSphericalHarmonicBasis`.** §III.8 shows
the earned name is the **trivial isotypic component**; "zonal" describes the
1-D case and would read as a special case of nothing.
**Compact groups only** — periodic lattices are isometries but non-compact;
they enter as finite cyclic groups on the discretized index set. State this in
the `IsometryGroup` docstring; it is a real limitation, not an oversight.

---

# Part VI — The action plan

Phases are ordered by dependency, not by size (D0.5). The `Q#` column maps each
item to the source plan's numbering so the two stay reconcilable.

## ⚡ Phase 0 — IMMEDIATE. Land regardless of any ruling below.

**Goal.** No object asserts a property its own data contradicts.
Every item here is true under every design choice in Phases 1–7.

| # | item | site | Q# |
|---|---|---|---|
| **0.1** | ⚡ **Stop fabricating the support, and delete the phantom read in the same expression.** `angular_frame` must not write `support=SPACE_SPHERE` over nodes with `‖Ω‖ ≠ 1`. §II.5 proves L1 and L4 are ONE fix. ⛔ **SPLIT at execution 2026-08-31 — see §II.9**: `0.1a` (3-D rules hand the frame their OWN measure; `[M]` bit-identical on 10 of 12 shipped rows) lands now; `0.1b` (the 1-D rows) **rides Phase 3.4**, because the lift `[-1,1] → S²` *does not exist as a map*; `0.1c` (the fold's `S^2/sigma_y` tag, §II.8's new leak) is gated on a `plan-authoring` §8 blast-radius measurement of `FunctionSpace` equality. | `directional.py:587-593` | — |
| **0.2** | ⚡ **A phantom component becomes unspellable.** *Proposed means as of 2026-08-31, SHARPENED at execution — see ⚠ below:* `axis_cosines(i)` raises for `i ≥ dim`. ⛔ a `raise`, **not** an `assert` — the canonical runner is `python -O`. ⚠ **BLOCKED on the full-suite census** (moved here from 2.4). | `directional.py` | — |

> ⚠ **0.2's means is not yet safe as written — `axis_cosines` is TWO
> functions sharing one name, and only one of them is lying.** `[M]`
> 2026-08-31, by reading all five production consumer clusters: the frame
> (`directional.py:589`) and `spherical_harmonics` (`:543-545`) ask *"what is
> the direction cosine along axis i?"* — a question with **no answer** for
> `i ≥ dim`, so their zeros are fabricated data. But
> `transport/mesh/axis.py:441` (`face_outflow_ordinates`),
> `sn/mesh/augmented_mesh.py:1837` (per-axis streaming) and
> `geometry/boundary/_specular.py:162` (the cosine measure `w·|μ_a|`) all key
> on the **MESH's** axis index and ask *"how much does this ordinate flow
> along axis i?"* — for which **zero is the correct answer**, and the
> accessor's own docstring says so (*"no quadrature data on this axis"* means
> *"no ordinate is outflowing on it"*).
>
> ⟹ a blanket `raise` would refuse three call sites that are asking a
> legitimate question. The elegant repair is Pattern 3 + Pattern 4: **split
> the accessor by its two meanings** — a strict `axis_cosines` whose domain is
> `i < dim`, and an honestly-named projection verb whose zero IS its answer —
> so the fabricated read is unspellable while the flow question keeps its
> home. ⭐ Note this is the SAME defect as L1 one level down: one set of zeros
> carrying two meanings (§II.2), here inside the accessor rather than across
> the quadrature↔basis boundary.
>
> ⚠ **The denominator is still open.** The three consumers above are keyed on
> the mesh axis, and whether any of them ever actually requests `i ≥ dim` at
> runtime is a question static reading cannot answer. The full-suite recording
> census (`scratch/_phantom_census_plugin.py`, non-invasive — records and
> returns the normal zeros) is running to settle it; **0.2 does not land until
> its table does.** `[M]` positive control: the instrument fires at exactly
> `directional.py:589`, `axis=[1,2] dim=[1] via=angular_frame`, with the 50
> tests of `test_frame.py` still passing.
| **0.3** | ⚡ **Correct `metric.py`'s "noise mode".** §VIII.1 — it is a null direction, not roundoff; the sentence is present-tense-FALSE and `metric.py` is read by every frame consumer. Re-derive the `_DENSE_METRIC_RCOND` cliff against the repaired frame. | `metric.py:105-119` | — |
| **0.4** | ⚡ **Correct `Quadrature.spherical_harmonics`'s docstring** — *"the other slots are filled with zeros"* is exactly the broken property. `[M]` this claim appears at **1 site tree-wide**. | `directional.py:538` | — |
| **0.5** | ⚡ **Correct the theory page's slab-frame passage** — a DIFFERENT claim from 0.4. It describes the slab frame as merely *"not even diagonal"* and cites P7's *"best diagonal candidate reads 1.806"*, inheriting the framing that omits RANK DEFICIENCY. Fix alongside 0.3. | `docs/theory/foundations/spherical_harmonics.rst:640-652` | — |
| **0.6** | ⚡ **A real `raise` in `_evaluate_real_sh`** refusing non-unit direction vectors, so the fabricated azimuth is unspellable at the source. ⚠ **sequencing**: this reddens the CURRENT slab pairing, so it lands WITH Phase 3.4 or immediately behind it. That is a `plan-authoring` §6b obligation, not a reason to weaken it. | `spherical_harmonic_basis.py` | — |
| **0.7** | ⚡ **Promote the reproducer to a strict-xfail gate.** `scratch/_pl_slab_defect_repro.py` LEG A becomes a committed `@pytest.mark.xfail(strict=True)` test, so the defect is VISIBLE in the suite and the marker self-retires the moment 3.4/4.1 land (an XPASS fails). ⛔ the marker form, **not** `pytest.xfail()` — the imperative call can never XPASS (`vv-principles` Mode 8, ninth class). ⚠ structure the body so **exactly one statement can fail and it is the documented one**. | `tests/sn/` | — |

**Done when:** `grep` finds no production site writing a `support` tag its own
nodes do not satisfy; and the §II.5 probe, run under the FULL suite, reports
zero phantom reads outside the sites 0.1/0.2 retire.

## Phase 1 — Derivations and rulings. No code. Prevents rework.

**Goal.** Every later phase builds on a written derivation rather than on a
remembered result.

| # | item | target | done when | Q# |
|---|---|---|---|---|
| 1.1 | ⛔ **RESCOPED 2026-08-31 — §II.7: the `SO(2)/S²` entry is largely AUTHORED ALREADY** in `registry.py` (`AngularSymmetry.support` → `SPACE_INTERVAL_M11` with the orbit description; `.reference` → `LEGENDRE`), so this is *reconcile and extend* (add `det P`, the named singular stratum, the action on the axis index set), NOT *derive*. **Catalog derivations.** Run §III.1 for `SO(2)/S²`, `ℤ₂` antipodal, `C_n`/`D_n` about an axis, the `O_h` sublattice for octant symmetry, `SO(3)` (diagonal, 1-D spherical), `SO(2)×ℝ_z` (1-D cylindrical). Record generators, syzygies, `P`-matrix, chart, pushforward measure. **Include §II.3's finding**: `det P` is not only the orbit radius but *the redundant harmonic*, in closed form. | `docs/theory/foundations/symmetry.rst` (new) | each entry has an explicit `det P` and a **named** singular stratum | Q0.1 |
| 1.2 | **The Gelfand-pair reading of Funk–Hecke.** `RΛM` with diagonal `Λ` is the spherical transform on a double-coset space. Supersedes any framing of harmonic projection as an optimization. | `frame.rst` / `scattering.py` | Funk–Hecke cited as a consequence of Schur, not as a special-function identity | Q0.2 |
| 1.3 | **Verify or kill the connection-term identification** (§III.2): is curvilinear `(1/r)∂_μ[(1−μ²)ψ]` the connection of the `SO(3)` phase-space quotient, with `(1−μ²) = det P`? Derive the reduction; do not argue from coincidence. ⭐ §II.3 supplies a **third**, *proved* appearance of the same polynomial — new evidence the source plan did not have. | `derivations/` script + `symmetry.rst` | a derivation, or an explicit *"coincidental in 1-D spherical"* ruling — **no third outcome ships** | Q0.3 |
| 1.4 | **Record the four gates** (Part IV) as a posing-table precondition block, with `C₆`/`LQ_N` as the worked negative and **the Part I defect as G2's worked positive**. | `operator_algebra.rst` | G0/G1/G2/G3 stated as distinct conditions on distinct objects | Q0.4 |
| 1.5 | **Record the placement ruling** (D0.3): no `reduce(A, H)`. Pre-commits against the fork an agent is primed to implement. | `symmetry.rst` + agent memory | ruling written before any code lands | Q0.5 |
| 1.6 | **Record `ℛ* = ℛ ⟺ H`-invariant Gram** as an instance of adjoint-vs-dual, and tightness descent as free inheritance. | `operator_adjoint.rst` | the quotient appears as a further entry in the adjoint-coincidence table | Q0.6 |
| **1.7** | ⭐ **NEW — the support algebra** (§III.9). Write the three verbs and their support rules; state that **`pushforward` is the only one that asserts**, and that the chart's missing codomain is the reason. Name Archimedes as the theorem licensing product quadrature. | `symmetry.rst` + `measure.py` | the three rules written in one place; `product()`'s `SPACE_SPHERE` literal cited as the worked leak | — |
| **1.8** | ⭐ **NEW — reconcile the vocabulary against the tree** (§III.3's table). Re-point every gate that names `Retract`, `Relation` or `WithTrace` at what exists (`AxisRetractionOperator`, `FunctionSpace.retraction`/`.section`), or mark it as depending on unbuilt machinery. | this plan | no phase gate names a type `[M]` absent from the tree | — |

## Phase 2 — Objects carry their intrinsic properties

**Goal.** Composability is decidable at the meeting place because both sides can
state what they are. **This is the phase D0.4 generates**, and it is the
prerequisite for every gate in Part IV.

| # | item | goal, separately from means | done when | Q# |
|---|---|---|---|---|
| 2.1 | **`Basis.domain`** (closes L2). *Proposed means, unverified:* a `domain: Space` property beside the existing `space`, giving the basis the same two-level structure the measure already has (`support` + `space`). | `SphericalHarmonicBasis.domain == SPACE_SPHERE`; every shipped `Basis` subclass answers | — |
| 2.2 | **G0 at frame construction** (closes L3). ⭐ **§II.7: the PREDICATE already exists** as `AngularSymmetry.admits_domain` (`measure.support == self.support`) — `[M]` reachable only from `select_quadrature`, whose 21 callers are **all tests**; and `[M]` `solve_sn` is HANDED a quadrature (required positional, no default), so the selector was never on the path. ⛔ **REFINED by §II.10 — NOT "just wiring".** The predicate is correct over an **incomplete domain vocabulary**: `support` derives from `continuous_isotropy` alone, so `folded_product`'s `S^2/sigma_y` matches **no** geometry and stage 0 would **refuse the shipped cylinder**. 2.2 owes the ontology a second slot (the *discrete* quotient the rule took — `SubgroupOfO3` already names it) **and** the `(CoordSystem, ndim) → key` join, which `[M]` exists nowhere. Do not hand-roll a second predicate; do extend this one's vocabulary. | `GalerkinFrame(SphericalHarmonicBasis(L=2), slab_measure)` **RAISES**, with a test witnessing it on the exact pairing that ships today | — |
| 2.3 | ⭐ **The typed `Chart`** (closes L6, L7). *Proposed means:* a map carrying `domain → codomain`; `pushforward` derives `support` from the chart rather than taking `new_space`; `product()` composes its two 1-D factors through the Archimedes chart instead of `column_stack` + a literal. | `grep` finds no `SPACE_SPHERE` literal at a measure constructor; `product(4,8).measure.support` is **derived**; the factor structure is recoverable from the built rule | — |
| 2.4 | **The slab rule declares its quotient group** (closes L5). Includes settling Part IV's obstacle 1 — the `SO2` axis-convention collision — and obstacle 2, the `invariance_group` vs `folded_by` conflation. ⛔ **MOVED to Phase 0 on 2026-08-31** — the full-suite phantom census was scheduled here, but **0.2 depends on it** (it is 0.2's denominator, not 2.4's). Scheduling a step's own precondition three phases downstream is the `plan-authoring` §6b defect with a *census* in place of a call site. Run and published at Phase 0. | `SubgroupOfO3.SO2.is_invariant(gauss_legendre(8).measure)` returns the derived answer, and the plan records whether that answer is `True` **with the derivation**, not the prediction | — |

## Phase 3 — The symmetry core

**Goal.** A quotient space, its retraction/section pair, and its invariant
sub-basis are all *derived from a catalog entry*.
⚠ **Phase 3 does not ship alone** — it lands with Phase 4.1, its first consumer.

| # | item | depends on | done when | Q# |
|---|---|---|---|---|
| 3.1 | `numerics/symmetry/` catalog: `IsometryGroup` entries from 1.1, each carrying its action on the axis index set, Haar (or counting) measure, fixed-point strata, and a symbolic test reproducing its own derivation. ⚠ `[M]` `numerics/symmetry.py` is a **module** today — this is a re-home; `plan-authoring` §6d (check the import edge BOTH ways, and look for the import-linter's declared contract first). | 1.1 | every entry's chart and measure reproduced by its own regression test | Q1.1 |
| 3.2 | `ReynoldsProjection` as a first-class operator with fixed domain/codomain, built via `.on(V)` like every other bound arrow, factored `π* ∘ 𝔄`. | 3.1 | `ℛ² = ℛ` and `𝔄π* = id` to machine precision | Q1.2 |
| 3.3 | `OrbitAxis` minting: `Quotient(axis, H) → (OrbitAxis, retraction/section)`. **Re-uses `AxisRetractionOperator` + `FunctionSpace.section`** (§III.3) — no new relational type. Stratification descriptor populated from `H`'s fixed-point set. | 3.1, 1.8 | the existing retraction/section tests pass with the new instance; no new relational type added | Q1.3 |
| **3.4** | ⭐ **`H`-stable sub-basis: the TRIVIAL ISOTYPIC COMPONENT.** `SphericalHarmonicBasis` → `{Y_ℓ0} ≅ {P_ℓ}` under `SO(2)`, mask **DERIVED by probe** (§III.8), generalising `MirrorEvenSphericalHarmonicBasis` so that it becomes an instance rather than a sibling. **THIS IS THE FIX for Part I.** | 3.1 | §VII's done-when 1 and 2 both hold | Q1.4 |

⚠ **Phase 3.4 leaves ONE encoding question deliberately open, and it should be
ruled here rather than assumed:**
* **(a) layout-preserving** — keep the rectangular `(L+1, 2L+1)` table and
  structurally zero the non-invariant columns, exactly as `MirrorEven` does
  (which already zeroes `|m| > ℓ` padding). Consumer-transparent; smallest diff.
* **(b) true `LegendreBasis` on `[-1,1]`** bound to the rule's own scalar
  measure — nothing fabricated anywhere, the basis's domain IS the measure's
  support. Changes the moment tensor to `(N, L+1)` and touches every consumer.
⟹ (b) is ontologically cleaner; (a) is precedent-shaped. **If (a) is taken its
docstring MUST say the layout is an ENCODING and name (b) as the successor**, or
the encoding will later read as the ontology. The `≅` in Q1.4's one line is
doing this work silently.

### ⭐ 3.4-R — the evidence, measured 2026-08-31 (was previously unmeasured)

**The blast radius of (b) is 15 sites, and only ONE is production.** `[M]`
by grep for explicit `(ℓ,m)` slot indexing into an SH/frame table, triaged by
meaning (the raw 22+19 hits include `values[:, 0, 0]` on FLUX arrays — a
different object — and the whole of `thermal_hydraulics`/`kinetics`, which index
unrelated `*_table` arrays):

| site | what it does | 1-D? |
|---|---|---|
| ⚠ `sn/acceleration/dsa.py:586` | `angular_frame(1).table[:, 1, 1]` — reads μ off the **ℓ=1 Cartesian slot** | **YES — DSA runs on slabs** |
| `basis/spherical_harmonic_basis.py:533-539` | the producer itself | — |
| `tests/numerics/test_spherical_harmonic_basis.py:52,60-63` | pins the ℓ=0/1 layout | 5 |
| `tests/sn/operators/test_solver_components.py:541,550,559-562` | orthogonality + layout pins | 6 |
| `tests/sn/acceleration/test_dsa_rate.py:747-748` | pins `table[:,0,0]`, `table[:,1,1]` | 2 |
| `tests/numerics/test_quadrature_directional.py:390` | `Y[:, 0, 0]` constant | 1 |

⚠ **This is a DIFFERENT inventory from §4.9's "12 sites"** — that one is *"consumes
the slab-GL-L2 frame as a fixture"*, this one is *"indexes the table by explicit
slot"*. Both are real; neither subsumes the other.

⭐ **The decisive argument is the project's own type-minting criterion**
(`coding-standards`, "Type vs property"): mint a type **iff** (a) there are ≥2
**non-isomorphic** realizations AND (b) a **non-identity morphism** is applied.
Here (a) holds — `{Y_ℓ^m}` on `S²` has `(L+1)²` members, `{P_ℓ}` on `[-1,1]` has
`L+1`; they are not isomorphic — and (b) holds: the quotient/Reynolds projection
and its section. **So (b) is what the rule selects**, and (a) is the criterion's
own named failure mode: *the same object wearing a label*.

⚠ And the technical point that reframes (a): **(a) does not remove the rank
deficiency, it only corrects the VALUES.** Structurally-zeroed columns still
contribute zero rows to the Gram, so `G` stays singular and `G⁺` keeps
discarding modes — the difference is that the discarded directions become honest
zeros instead of fabricated constants. `[R]` (a) is nonetheless *well-posed* by
§III.5: `ker(synthesis)` is then spanned by the dead slots and the per-ℓ
multiplier maps dead slots to dead slots, so the kernel is invariant. **Both
options give the right answer; only (b) makes the wrong one unspellable.**

## ⏸ COMPACTION POINT 1 — after Phase 3

Carry forward: the phase→commit table; the measured baseline; §VII's refuted
premises; and the Phase 3.4 encoding ruling once made.

## Phase 4 — The retrofit and its gates

| # | item | depends on | done when | Q# |
|---|---|---|---|---|
| 4.1 | **Retrofit the 1-D slab angular axis as `S²/SO(2)`.** The hard-coded `μ`-only path becomes the zero-configuration instance of the general machinery, **deleted, not shadowed**. | 3.1–3.4, 0.1 | see the corrected gate below | Q2.1 |
| 4.2 | **Normalization audit, forced by 4.1.** The pushforward `dΩ/4π ↦ dμ/2` fixes the `2π` exactly once. | 4.1, 2.3 | every angular normalization traced to a single pushforward statement; any discrepancy filed as an ERR, never silently patched | Q2.2 |
| 4.3 | **Realization/quotient commuting square:** `ℛ ∘ bind(K).on(V) == bind(K).on(V/H)` for `S`, `F`, `C`. The central correctness assertion, and cheap. | 3.2, 4.1 | passes for all three kernel kinds; a deliberately `H`-asymmetric `Σ_t` REDs it (negative control) | Q2.3 |
| 4.4 | **Quadrature symmetry annotation** + the **coarse-rule class**: `level_symmetric(4)` L=3 and `folded_product(2,4)` L=3 are DENSE/rank-deficient because `L` exceeds what the rule resolves. Route them to the Q6 `exactness` machinery, which `[M]` is currently disconnected (§I.3). | 3.1 | every quadrature carries a group; an unannotated set refuses to participate in a quotient; the two coarse rows have a stated verdict | Q3.1 |
| 4.5 | **G2 check** at realization. | 4.4, 3.4 | **`C₆` hexagonal + `LQ_N` REDs naming the missing subgroup relation**, AND **the Part I slab pairing REDs** — two independent positive controls | Q3.2 |
| 4.6 | `transport/symmetry.py`: `H_max = Sym(geom) ∩ Sym(mat) ∩ Sym(BC)`, reported beside the acyclicity certificate; user-supplied `H` validated against it. | 3.1 | an asymmetric material patch reports the REDUCED `H_max`, not the geometric one | Q3.3 |
| 4.7 | **G3 check** at Declaration: per-term commutation, randomized on a few `g` plus exact structural checks where available. | 4.6, 4.3 | each of `L, C, S, B, F` reports its own descent verdict; the sum descends iff all terms do | Q3.4 |
| 4.8 | **Ray-effect residue metric** when G2 fails and the user forces the reduction. | 4.5 | forced-`C₆` produces a measurable, reportable residue rather than a silent wrong answer | Q3.5 |
| **4.9** | ⭐ **NEW — re-key the collateral consumers.** `[M]` **12 sites repo-wide consume the slab-GL-L2 frame as a fixture** (7 tests, 4 docs, 1 production comment). In particular `tests/numerics/test_frame.py:906` asserts the live `ℓ=2` Gram diagonal is `[0.4, 0.8, 0.8]` at `rtol=1e-12` — **the two `0.8`s ARE the fabricated slots.** A green test PINS the defect. | 4.1 | every one of the 12 re-pinned, **not deleted** (`coding-standards`: retirement means test migration) | — |

### ⛔ Phase 4.1's gate, CORRECTED — the source plan's version is designed-green

The source plan's Q2.1 gate reads *"all existing 1-D results bit-stable"*.
`[M]` **that is wrong as stated and would pass for the wrong reason.** The
retrofit is bit-identical at `L ≤ 1` (`Δ = +0.000e+00`, measured in-process) and
**must change `L ≥ 2`**. It would "pass" only because nothing in the corpus
exercises `L ≥ 2` on a 1-D chart (§I.3) — i.e. it certifies the coverage hole,
and **as written it forbids the correction it should be demanding**
(`plan-authoring` §10: a gate that cannot detect its own campaign's success).

**The honest gate:**
> `L ≤ 1` bit-identical to pre-retrofit (`Δ = 0.000e+00`); **`L ≥ 2` MUST
> CHANGE**, from `−3.764705882353` to `+4.000000000000` on §I.1's fixture; and
> the 3-D rules are bit-identical at every `L` (`[M]` the candidate patch is
> **exactly inert** on a 3-D rule, `max|ΔY| = 0.0`).

## Phase 5 — Curvilinear S_N (the forcing consumer)

| # | item | depends on | done when | Q# |
|---|---|---|---|---|
| 5.1 | 1-D spherical angular axis as the diagonal `SO(3)` quotient; `OrbitAxis` carries the `μ = ±1` singular stratum with degenerate trace measure. | 4.1, 3.3, 1.3 | the stratum is degenerate-measure **by derivation from `H`**, not by hand | Q4.1 |
| 5.2 | The angular redistribution term as the quotient connection, **iff 1.3 verifies**. If 1.3 comes back negative this degrades to "hand-derived term, singular stratum still supplied by the quotient" — **the stratification payoff survives either way.** | 1.3, 5.1 | the α-recursion's degenerate starting direction lands on the singular stratum with no special case | Q4.2 |
| 5.3 | 1-D cylindrical as the `SO(2) × ℝ_z` quotient. | 5.1 | cylindrical and spherical share one code path differing only in catalog entry | Q4.3 |
| 5.4 | Curvilinear sweep admissibility on the quotient: the causal order must be compatible with the stratification. | 5.1 | acyclicity certificate emitted for curvilinear geometries | Q4.4 |

## ⏸ COMPACTION POINT 2 — after Phase 5

## Phase 6 — Spatial quotients (independent of Phase 5; may interleave)

| # | item | done when | Q# |
|---|---|---|---|
| 6.1 | **Reflective BC as chart boundary, not boundary operator.** A reflective face states "the domain is a fundamental domain of `H`"; on the quotient the identification is structural and there is **no `B` term** for that face. | quarter-core constructs with `B` absent on reflective faces; answer bit-consistent with the full-core solve averaged | Q5.1 |
| 6.2 | **Monodromy rank under quotient.** Prediction: closing-edge rank falls by `\|H\|` where the action is free, less where it is not. **Measure, do not assume** — fixed strata break freeness. | rank table extended with a MEASURED quotient column | Q5.2 |
| 6.3 | 2-D quarter symmetry derived as a quotient rather than a hard-coded acyclicity case. | the 4-orthant DAG produced by the quotient, one pass, zero boundary iterations | Q5.3 |
| 6.4 | Hexagonal `C₆`/`D₆` 1/6- and 1/12-core, **gated on 4.5 passing** — surfacing the quadrature requirement is the point. | either a compatible quadrature is supplied, or the configuration is refused with the reason | Q5.4 |

## Phase 7 — Isotypic decomposition

| # | item | depends on | Q# |
|---|---|---|---|
| 7.1 | `IsotypicDecomposition(A, H)`: Schur blocks over irreps; trivial block = the quotient space already built in Phase 3. | 3.2 | Q6.1 |
| 7.2 | `k` restricted to the trivial block by §III.8's theorem, cited in the posing table as a **precondition-free** reduction. | 7.1 | Q6.2 |
| 7.3 | Higher harmonics labelled by irrep; Riesz contour per block. | 7.1 | Q6.3 |
| 7.4 | Reactor noise: azimuthal mode index = `SO(2)` irrep; needs `𝔽 = ℂ` for 1-dimensional irreps. Lands **with** the ground-field tag, not before. | 7.1, `𝔽` tag | Q6.4 |

## Explicitly refused

| item | reason |
|---|---|
| General Gröbner / Procesi–Schwarz orbit-space engine | ~12 catalog entries ever; engine failure modes worse than the cost saved |
| `reduce(A, H)` post-assembly projection | stage inversion; pays full assembly cost; defeats the purpose |
| Non-compact groups as first-class | no Haar average; lattices enter as finite cyclic on the discrete index set |
| A new relational type for quotients | `[M]` `AxisRetractionOperator` + `FunctionSpace.section` already exist |
| A sixth space kind | quotienting is axis-level; the five kinds are untouched |
| Symmetry in Monte Carlo (tally symmetrization) | different mechanism — research tier, not this layer |
| ⛔ **A pseudo-inverse / Tikhonov / truncated-SVD repair of the rank deficiency** | §VII.1 — measured refuted, and independently confirmed by §III.6 |
| ⛔ **A rank-revealing "drop any dependent column" reduction** | would pass the isotropic control while keeping `Y₂^{+1} ∝ μ√(1−μ²)`, a non-polynomial that makes `c₀,c₁,c₂` stop being the Legendre coefficients — **designed-green** |

---

# Part VII — Done-when, for the plan as a whole

Checkable predicates, in order of strength:

1. `scratch/_pl_slab_defect_repro.py` LEG A reads `+4.000000000000` at
   `L = 0, 1, 2, 3` — today it reads `−3.764705882353` at `L ≥ 2`.
2. `[M]` a fresh census of the isotropic-flux `ℓ ≥ 1` moment over every shipped
   `(rule, L)` reads `0` to machine precision — today 4 of 24 rows do not.
3. Constructing a frame from a basis whose domain does not contain the measure's
   support RAISES, with a test witnessing it on the pairing that ships today.
4. No production site writes a `support` tag its own nodes do not satisfy, and
   no production site reads a direction component beyond the rule's `dim`.
5. `[M]` `C₆` + `LQ_N` REDs at realization, naming the missing subgroup relation.
6. The full 13-tree gate is green, with the per-tree deltas predicted before it
   is run (`reference_test_execution_env`).

---

# Part VIII — Collateral findings

## VIII.1 ⛔⛔ Campaign 1's P7 measured this null direction and called it "noise"

`orpheus/numerics/metric.py:111-112`, verbatim:

> `[M]` 2026-08-30 on the flagship slab Gram (singular values
> `2.71 / 1.42 / 4.92e-1 / 4.74e-2` live, one `~1e-16` **noise** mode)

`[M]` my measurement of that exact frame: `2.707550, 1.419220, 4.924502e-1,
4.744684e-2, 9.483994e-18` — **bit-for-bit the same four**. And the fifth is
**not noise**: its right-singular vector `[−0.4472, 0, +0.4472, 0, +0.7746]`
= `(1/√5)(−Y₀ + Y₂⁰) + √(3/5)·Y₂²` evaluates to **4.98e-16 at every node** —
an exact linear dependence, and §II.3 derives it in closed form.

⟹ The defect had been SEEN, one day before it was found, and mislabelled. The
comment then reasons carefully about a *"cliff's ANALYTIC edge"* treating
`σ₃ = 4.74e-2` as the smallest GENUINE live mode, and pins
`_DENSE_METRIC_RCOND = 1e-12` to sit *"~4 orders above the noise floor"* — i.e.
**the pin was tuned to discard the signal that reports this bug.**

⭐ The transferable half: the measurement was right and the INTERPRETATION
inverted it, and the rigour of the surrounding analysis is exactly what made it
unquestionable — `plan-authoring` §2's `[M]`-scope defect wearing a
carefully-derived cliff.

⚠ **Consequence for P7's own conclusion:** the DENSE-metric repair makes
Parseval a theorem AND **removes the singularity from the only property that
could report it.** P7 must not be read as having addressed this defect; the two
repairs are mutually exclusive in intent (fix the metric vs fix the basis).

## VIII.2 The diagnostics that exist are permissive, not inert

`[M]` `discrete_gram_structure` has exactly **one** production read
(`frame.py:286`) and it **routes a metric choice, never refuses**.
`Basis.gram_structure`'s `NotInvertible` cannot fire for SH (declared `DIAGONAL`
unconditionally, and it sits on `project`/`gram`, which the scattering path never
calls — 0 of 38 hits in `transport/operators/`). **Nothing anywhere computes a
rank or condition number of a frame Gram.** `Basis.mass_matrix` is a genuinely
inert authored diagnostic (0 production callers).

## VIII.3 The precedent, every adjective verified by reading

`MirrorEvenSphericalHarmonicBasis` is a genuine `@dataclass(frozen=True)`
**subclass** (frozen verified by attempted mutation); field
`mirror_axis: int = 1` **has a default**; **no registry, no factory
constructor**; keeps the full `(L+1, 2L+1)` layout and zeroes columns; the mask
is a `@cached_property` **DERIVED** by evaluating `_evaluate_real_sh` at 5
hardcoded generic directions and their mirrors, raising `RuntimeError` if a slot
classifies as neither even nor odd.
`Quadrature._harmonic_basis` selects on **one field, `folded_by`**, and its own
docstring records the seam this plan consumes: *"folding by any other group
refuses here until a consumer exists."* **The slab is that consumer.**

## VIII.4 Blast radius — measured, and far smaller than feared

`[M]` **4** `FrameBase` subclasses (denominator 457 `ClassDef`s / 340 files);
**11** production constructor sites in 7 files; the whole angular story is **2**,
chained `Quadrature.angular_frame(L)` → `HarmonicFrame.from_galerkin`, both
interned per `(rule, L)`. **4** production `HarmonicFrame.for_space` sites
(`scattering.py:811`, `:960`; `fission.py:305` and `n2n.py:200` **hardcode
`L = 0`**). One frame bypass: `dsa.py:585` reads `angular_frame(1).table` raw.
**The only ℓ-stacked per-mode apply is `material_field.py:244/270/284`.**
⭐ **HARD NEGATIVE:** CP, MoC, MC, diffusion, homogeneous, kinetics, fuel, TH —
**24 `.py` files, 0 angular-moment lines.**
**DSA is immune because it binds `L = 0`**, where `Λ = σI` and the kernel is
trivially invariant. **That is why nobody noticed.**

⟹ Today the defect has **one** consumer. `[R]` It acquires more the moment
`fission.py` / `n2n.py` stop hardcoding `L = 0` — which is exactly what #426
proposes for `(n,2n)`. **That is the coupling that paused Campaign 2.**
`[R, unconfirmed]` energy condensation on NON-NESTED grids is flagged as a
second consumer of the same shape (`OverlapBasis` declares
`PARTITION_OF_UNITY`, trial Gram measures DENSE, `frame.py:206-208`); the tree
census did not confirm an ℓ-stacked per-mode apply there. **Unresolved.**

---

# Part IX — ⛔ Refuted premises, kept in place (`plan-authoring` §3)

### IX.1 ⛔ "The measured Gram (pseudo-inverse) fixes it for every configuration"
**Whose:** the user's, ruled 2026-08-31 as the fix shape, with the explicit
instruction *"let's try it and see if it fixes"*.
**Refuted the same day, by measurement.** `G⁺` fixes every case where the Gram
is merely DENSE — including two that would otherwise have been missed
(`folded_product` L≤2, `level_symmetric(4)` L=3). It **cannot** fix rank
deficiency: `[M]` on slab GL8 L=2 the min-norm coefficients are
`c = 0.800·Y₀ + 0.200·Y₂⁰ + 0.346·Y₂²`, which reconstruct ψ **exactly**
(2.2e-16) yet are not `e₀`; `Λ` then scales the `ℓ=0` part by `Σ_s0` and the
`ℓ=2` parts by `Σ_s2`. Independently corroborated by §III.6.
**What survives, and it shaped the plan:** the instinct that ONE general
mechanism should cover every configuration was right — §III.8 supplies it, and
the DENSE-Gram route remains correct for the two dense-but-full-rank rows (4.4).

### IX.2 ⛔ "rank(G) == #live is the well-posedness predicate"
**Whose:** mine. **Refuted:** SUFFICIENT, not necessary — as a *refusal* it
over-rejects padding and harmless intra-block degeneracy. The necessary and
sufficient condition is §III.5's kernel-invariance.
**But** `[M]` rank-over-live is a **perfect DETECTOR** on the 8 shipped frames.
⟹ **excellent detector, bad refusal predicate.**

### IX.3 ⛔ "the redundancy is countable a priori: every even m≥2 slot is redundant"
**Whose:** mine. **Refuted out of sample:** `[M]` over
`n_ordinates ∈ {2,4,8,16,32} × L ∈ 1..6` — correct on **0 of 30 rows**. Two
errors: `(1,1)` is DEAD not live (the `ℓ=1` Cartesian special case sets
`Y[:,1,2] = mu_y = 0`), and the true rank saturates at `2L+1` for `L ≥ 3`
(measured 7, 9, 11, 13 at L = 3,4,5,6).
⟹ **I derived the rule from a basis I had not READ.** Reading
`_evaluate_real_sh` gave the mechanism in one pass and explained both errors.
§II.3 is the corrected, derived statement.

### IX.4 ⛔ "folded_product is spurious at ℓ≥1"
**Whose:** mine, published in a table earlier the same day.
**Refuted:** measured with a PLAIN `SphericalHarmonicBasis`, not the basis
production binds (`MirrorEven…`). On production's own path `folded_product` is
**clean at L ≤ 2**. ⟹ `vv-principles` #28: build the operand the way production
builds it.

### IX.5 ⛔ "all existing 1-D results bit-stable" as Phase 4.1's gate
**Whose:** the source symmetry plan's Q2.1. **Refuted:** see the corrected gate
in Part VI — it is designed-green and forbids the correction it should demand.

---

# Part X — V&V battery, by tier

**T0 (algebraic, machine precision):** `ℛ² = ℛ`; `𝔄π* = id`; `ℛ* = ℛ` under the
invariant Gram **and FAILS under a deliberately non-invariant one** (negative
control tying to §III.4).

**T1 (structural):** the realization/quotient commuting square (4.3);
`π_*Q = Q_quot` on the invariant basis with order of accuracy preserved;
sub-basis `H`-stability; `B/A` tightness **inherited, not recomputed**;
**G0 refuses the shipped slab pairing** (2.2).

**T2 (solve-level):** full-space solve averaged == quotient solve to iteration
tolerance (6.1); the **corrected** 4.1 gate; monodromy rank measured (6.2).

**T3 (physics):** `k` from quotient == `k` from full space; forced-G2 ray-effect
residue measured (4.8); curvilinear verification against a known solution;
**§I.1's infinite medium at `L = 2, 3`** — the plan's flagship red-to-green.

---

# Part XI — Confidence ledger

**High**
- The Part I defect, its mechanism, and its scope (all `[M]`, reproducer committed)
- `S²/SO(2) = [−1,1]`; `det P = 4(1−μ²)` as the squared orbit radius
- ⭐ `det P` **is** the redundant harmonic (§II.3 — derived AND matching the SVD)
- The phantom census: ONE production consumer, at the fabrication site (§II.5)
- Reynolds as retraction; `(π*, 𝔄)` re-using the existing retraction/section pair
- `ℛ* = ℛ` requiring an `H`-invariant Gram; tightness descent
- Gelfand pair `(SO(3), SO(2))`; Funk–Hecke as Schur
- `k`'s fundamental mode in the trivial isotypic component
- G2 as a genuine uncaught error class; `C₆ ⊄ O_h`
- The singular stratum `μ = ±1` coinciding with the vanishing trace measure
- Quotienting is axis-level; placement at stage 2

**Moderate — verify before the dependent item ships**
- The connection-term identification (blocks 5.2's *documentation*, not 5.1's code)
- Dominant `α` in the trivial isotypic block (needs the Perron-root argument)
- Monodromy rank dividing by `|H|` — stated as a measurement, deliberately not a theorem
- That the 4.2 normalization audit finds something. ⭐ **Partially discharged
  already**: the reconstruction uses the continuum `(2ℓ+1)` while the measured
  Gram is neither diagonal nor full rank — an inconsistency between the moment
  definitions and the scattering `Λ`, exactly as predicted

**I do not know**
- Whether a hexagonally-compatible high-order angular quadrature exists in usable
  form for 6.4, or whether that item terminates in "refuse the configuration"
- How `H_max` behaves on unstructured meshes where `Sym(geometry)` is only
  approximately realized by the discretization — the tolerance question
- `[R]` whether the `SO2` axis-convention prediction (Part IV obstacle 1) is right
- `[R]` whether energy condensation is a genuine second consumer (§VIII.4)
- Which of Phase 3.4's two encodings is right — **deliberately unruled**

---

# Part XII — ⏏ THE EXIT GATE — what returns us to Campaign 2

This plan is a SIDE-QUEST. It was opened mid-flight from
`.claude/plans/cs4c_binding_design.md` §18.6 step 1 (#426's missing (n,2n)
measurement), and the point of an exit gate is that **the whole plan is not the
prerequisite.** Most of it is genuinely separate work.

## XII.1 The coupling, stated precisely

#426 proposes that `(n,2n)` stop discarding its `ℓ ≥ 1` Legendre moments.
`[M]` `n2n.py:200` hardcodes `L = 0` — **that is the only reason the channel is
not already affected** by the Part I defect. Restoring the moments routes them
straight into the broken path. So the block is real, and it is narrow: it is a
block on the **P_L path being correct on a 1-D chart**, nothing more.

## ⚠ XII.1b The exit path is not a SEQUENCE — three items are one landing

`[M]`/`[R]` 2026-08-31, found at execution. The tracker lists 0.1, 2.2 and 3.4
as separate rows and the critical path reads them in order. **They cannot land
in order**, for one reason stated three ways:

* **0.1b** — `angular_frame` cannot stop forging `support=S^2` on a 1-D rule,
  because the honest alternative (the rule's own `[-1,1]` measure) is not
  consumable by `SphericalHarmonicBasis`. §II.9: *the lift does not exist*.
* **2.2** — G0 (`admits_domain`) refuses exactly the pairing production ships
  today, so arming it before the pairing changes reddens every slab solve.
  That is `plan-authoring` §6c's mirror: not *a gate with no witness*, but
  **a gate whose only witness is production**.
* **0.6** — `_evaluate_real_sh` raising on non-unit directions is the same
  refusal one layer down (the plan already flagged this one).

⟹ **0.1b + 0.6 + 2.2 + 3.4 are ONE commit**, and the §6b rule applies with a
*basis change* in place of a signature change. Everything else in Phase 0 and
Phase 2 is genuinely independent and lands first. Sequencing the four
separately would leave the tree red between them — or, worse, tempt a session
to weaken the gate so the interval stays green, which is how a real refusal
becomes a warning.

## XII.2 ⏏ The MINIMUM chain to return — the critical path

| # | why it is on the path |
|---|---|
| **Phase 0** (all 7) | ⚡ unconditional; 0.1/0.2 remove the fabrication the defect rides on |
| **Phase 1.1** — *the `SO(2)/S²` entry only* | 3.4's probe needs that one catalog row. §III.2 already derives it; the other entries serve Phases 5–6 and are **NOT on the path** |
| **Phase 2.1–2.4** | `Basis.domain`, G0, the typed `Chart`, and the slab rule's declared quotient group — the objects 3.4 reads |
| **Phase 3.1** — *scoped to the `SO(2)` entry* | the catalog home + the probe |
| **Phase 3.4** | ⭐ **the fix itself** — the trivial isotypic component |
| **Phase 4.1** | the retrofit, under the CORRECTED gate (Part VI) |
| **Phase 4.9** | re-key the 12 collateral sites, incl. the green test that pins the defect |

**Explicitly NOT on the path** (may land in any later session): Phase 1.2–1.8
except 1.1's one entry, Phase 3.2/3.3, Phase 4.2–4.8, and the whole of Phases
5, 6, 7. `[M]` §VIII.4 is what licenses this: today the defect has **one**
consumer, and CP/MoC/MC/diffusion/kinetics/fuel/TH carry **0 angular-moment
lines**, so nothing else is waiting on the layer.

## XII.3 ⏏ The exit predicate — all six, measured, before returning

1. `scratch/_pl_slab_defect_repro.py` LEG A reads **`+4.000000000000`** at
   `L = 0, 1, 2, 3`.
2. Phase 0.7's strict-xfail gate **XPASSes and its marker is deleted** — the
   self-retiring proof that the defect is closed.
3. `[M]` the isotropic-flux `ℓ ≥ 1` moment census reads `0` to machine precision
   on every shipped `(rule, L)` **except** the two coarse-rule rows, which carry
   a written verdict routed to Phase 4.4.
4. `[M]` 3-D rules bit-identical at every `L` (`max|ΔY| = 0.0`) — no 2-D/3-D
   result moves.
5. The **full 13-tree gate green**, with per-tree deltas PREDICTED before the
   run (`reference_test_execution_env`). ⚠ The baseline is UNMEASURED at this
   HEAD: step 4's `10106 / 0 / 19 sk / 66 xf` is the last full gate, `10116`
   predicted-not-run at CS4c §18.7.
6. `sphinx -W` clean and `dead_references` 0 dead — the four corrected prose
   surfaces (0.3/0.4/0.5 and `metric.py`) land with the code.

## XII.4 ⏏ On return

Resume `.claude/plans/cs4c_binding_design.md` **§18.6 step 1** — #426's missing
measurement, the Be-reflected fixture solved with and without the `ℓ ≥ 1`
(n,2n) source. It is now runnable at `ℓ = 1` **and** `ℓ = 2`, which it was not
when this side-quest opened.
⭐ And it gains a control it could not have had: `[M]` the same fixture at
`L ≤ 1` is bit-identical across the repair, so any `Δk` the (n,2n) moments
produce is attributable to them and not to the P_L machinery.

---

# Part XIII — Execution tracker

Update this table **in the same commit** as the work. `plan-authoring` §3: a
row that becomes false by being FIXED is edited here, not annotated elsewhere.
Status: `☐` not started · `▶` in flight · `✅ <hash>` landed · `⛔` refuted/dropped.

| # | item | on exit path? | status |
|---|---|:---:|---|
| 0.1a | 3-D rules hand the frame their OWN measure (`[M]` bit-identical, 10/12) | ⏏ yes | ☐ |
| 0.1b | the 1-D rows — ⛔ **rides 3.4**, the lift does not exist (§II.9) | ⏏ yes | ☐ (with 3.4) |
| 0.1c | the fold's `S^2/sigma_y` tag — gated on the §8 blast radius | ⏏ yes | ☐ |
| 0.2 | phantom component unspellable — ⚠ means SHARPENED, blocked on the census | ⏏ yes | ☐ |
| 0.2-census | full-suite phantom census (moved here from 2.4) | ⏏ yes | ▶ running |
| 0.3 | `metric.py` "noise mode" + rcond re-derivation | ⏏ yes | ✅ written (97-line re-derivation; value unchanged) |
| 0.4 | `directional.py:538` docstring | ⏏ yes | ✅ written |
| 0.7b | ⭐ NEW — gate the SECOND symptom (the `DenseMetric` RAISE) | ⏏ yes | ✅ written |
| 0.5 | `spherical_harmonics.rst` rank-deficiency correction | ⏏ yes | ✅ written |
| 0.6 | `_evaluate_real_sh` raises on non-unit directions | ⏏ yes | ☐ (with 3.4) |
| 0.7 | strict-xfail reproducer gate + **ERR-080** minted | ⏏ yes | ✅ written |
| 1.1 | catalog derivations — **`SO(2)/S²` entry is on the path** | ⏏ partial | ☐ |
| 1.2 | Gelfand-pair reading of Funk–Hecke | no | ☐ |
| 1.3 | connection-term verify-or-kill | no | ☐ |
| 1.4 | record the four gates | no | ☐ |
| 1.5 | record the placement ruling | no | ☐ |
| 1.6 | `ℛ* = ℛ` as adjoint-vs-dual | no | ☐ |
| 1.7 | the support algebra, written | no | ☐ |
| 1.8 | vocabulary reconciliation | no | ☐ |
| 2.1 | `Basis.domain` | ⏏ yes | ☐ |
| 2.2 | G0 at frame construction | ⏏ yes | ☐ |
| 2.3 | the typed `Chart`; pushforward derives support | ⏏ yes | ☐ |
| 2.4 | slab declares its quotient group; SO2 axis collision; **full-suite phantom census** | ⏏ yes | ☐ |
| 3.1 | `numerics/symmetry/` catalog (⏏ scoped to `SO(2)`) | ⏏ partial | ☐ |
| 3.2 | `ReynoldsProjection` | no | ☐ |
| 3.3 | `OrbitAxis` + retraction/section minting | no | ☐ |
| 3.4 | ⭐ **trivial isotypic sub-basis — THE FIX** | ⏏ yes | ☐ |
| 3.4-R | ⚠ the encoding ruling | ⏏ yes | ✅ **RULED 2026-08-31 (user): (b) true `LegendreBasis`** — see below |
| 4.1 | retrofit the slab axis (CORRECTED gate) | ⏏ yes | ☐ |
| 4.2 | normalization audit | no | ☐ |
| 4.3 | commuting-square test | no | ☐ |
| 4.4 | quadrature symmetry annotation + the coarse-rule class | no | ☐ |
| 4.5 | G2 check | no | ☐ |
| 4.6 | `H_max` | no | ☐ |
| 4.7 | G3 check | no | ☐ |
| 4.8 | ray-effect residue | no | ☐ |
| 4.9 | re-key the 12 collateral sites | ⏏ yes | ☐ |
| 5.1–5.4 | curvilinear S_N | no | ☐ |
| 6.1–6.4 | spatial quotients | no | ☐ |
| 7.1–7.4 | isotypic decomposition | no | ☐ |
| ⏏ | **the six exit predicates (XII.3)** | — | ☐ |

**Open rulings blocking execution:** ⭐ **NONE.** 3.4-R was ruled by the user on
2026-08-31 — **option (b), the true `LegendreBasis` on `[-1,1]`**, on the
measured evidence in Part VI (15 slot-indexing sites, only ONE of them
production; and the project's own type-minting criterion selects it, since
`{Y_ℓ^m}` and `{P_ℓ}` are non-isomorphic and a non-identity morphism — the
quotient projection — applies). The 0.1c fold tag was ruled the same day:
**measure the blast radius, then land it in Phase 0** (§II.10 gives the reason
this is right: the forgery is what CONCEALS the missing discrete-quotient slot).

---

# Part XIV — ▶ Resume surface

**Read in this order on pick-up:** Part XIII (the tracker — what is done),
Part XII (the exit gate — what returns us to Campaign 2 and what is NOT on that
path), then the phase you are executing. Part IX first if you are about to
re-propose something; four premises are already refuted there, three of them
mine.

**⚡ Phase 0 is unblocked and lands first**, independent of every ruling above.
**The only open ruling on the exit path is 3.4-R** (the encoding).

**This plan blocks** `.claude/plans/cs4c_binding_design.md` §18.6 step 1
(#426's missing (n,2n) measurement), PAUSED. The coupling is §VIII.4 and the
narrow form is Part XII.1: `n2n.py:200` hardcodes `L = 0`, which is the only
reason that channel is not already affected.

**Artifacts** (⚠ all in UNTRACKED `scratch/`; load-bearing content carried in
this file):
`_pl_slab_defect_repro.py` (the reproducer, two self-checking legs),
`_pl_slab_fix_probe.py` (in-process candidate fix; `[M]` L≤1 bit-identical,
L≥2 → `+4.000000000000`, inert on 3-D),
`_phantom_axis_probe.py` (the §II.5 census),
`n2n_pl_diagnosis_main.md`, `n2n_pl_consistency_literature.md`,
`n2n_pl_frames_attack.md`, `n2n_pl_blast_radius.md`.

**Standing constraints:** Host → `.venv/bin/python`; canonical runner
`-O -m pytest -p no:randomly -m "not slow" -q` SERIAL; branch + ff-merge; never
`git add -A`; never `git stash`; commit messages via `git commit -F -` with a
**quoted** heredoc; no source edits under running gates; `sphinx -W` +
`dead_references` at every phase exit; mutation batteries as subprocess-scoped
pytest plugins.
