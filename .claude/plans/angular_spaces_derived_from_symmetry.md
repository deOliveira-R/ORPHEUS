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
orbit-space ENGINE — but BUILD ITS EMBRYO.** Computing an orbit space from
scratch is mechanical (invariant ring → syzygy ideal → Procesi–Schwarz PSD
condition), but the groups that occur in transport number about a dozen. A
Gröbner engine is abstraction without a consumer, and its failure mode is
debugging elimination orderings instead of transport. Each catalog entry is
**derived once by the procedure** (§III.1), then recorded with its derivation in
the docstring and a symbolic regression test.

⭐ **SHARPENED 2026-08-31 (user).** *"Even if we do not build the entire engine
for operations on symmetry and quotient, the right solution is probably still to
build an embryo of the machinery, such that concepts are not stringly typed and
maybe some operations can still happen."* This corrects how D0.1 had been read —
including by me, twice this session. **What is refused is the DERIVATION engine,
not the machinery.** Those are different objects, and conflating them is what
left `Space = str` standing.

⭐⭐ **SHARPENED AGAIN, same day (user) — and the first version of this clause
was too weak.** The engine is **not ruled out**; it is ruled **not yet**:

> *"We're not outright ruling out building the engine. We're ruling that we will
> not prematurely build the engine. The embryo should be such that if the day
> ever arises that we decide building the engine is the right step, it should be
> a development of the embryo, instead of having to do the entire engine from
> scratch for a code base that was not ready to receive it."*

⛔ My first version said *"the catalog and the engine have the same INTERFACE …
a second backend behind an unchanged signature."* **That is not enough, and it
is the twin-path risk wearing a compliment.** A shared signature guarantees only
that the *call site* survives; the engine would still be an independent
implementation with its own representation of polynomials, ideals and PSD
conditions, plus a translation layer to whatever the catalog happened to store.
That is a from-scratch build with a seam, which is exactly what the ruling
forbids.

⟹ **The binding requirement is on the DATA MODEL, not the signature: a catalog
entry must BE the derivation procedure's OUTPUT, not a human summary of its
answer.** §III.1's procedure emits, per entry: the invariant generators
`p₁…p_k`; the syzygy ideal `I`; the matrix `P_ij = ⟨∇pᵢ, ∇pⱼ⟩` and `det P`; the
chart; the pushforward measure; and the stratum where `det P` vanishes. **Those
are the entry's FIELDS.** Then the engine is *"compute these fields instead of
reading them from a literal"* — a development, with no new vocabulary and no
seam.

⭐ **And the acceptance suite for that engine already exists, written years
before it.** D0.1 already requires every entry to ship *"a symbolic regression
test reproducing its own derivation"*. If the fields are the derivation's
outputs, those ~12 tests **are the engine's specification and its acceptance
gate** — the engine ships when it reproduces them by computation. `vv-principles`
#17: a spec written before the implementation cannot be shaped to flatter it.

⚠ Two further obligations that follow, both cheap and both easy to lose:
* **Provenance per entry** (`derived_by: "hand" | "engine"`), so a *mixed* state
  is expressible and visible. An incremental engine rollout is exactly a mixed
  state, and without the field it is unrepresentable — the migration would have
  to be all-or-nothing.
* **The refusal must name the missing WORK, not just the gap.**
  `M.quotient(C₆)` raises with *"no entry for `S²/C₆`; derive it per §III.1
  (invariants → syzygy → `P` ⪰ 0) and register it, or implement the engine"* —
  which turns the engine's absence into a work item a fresh session can pick up,
  rather than a wall.

⟹ **the falsifiable form of this ruling, for a future session to check:** *given
a catalog entry, could an engine populate its fields without introducing a
single new type?* If the answer is no, the embryo has drifted from being a seed
and the ruling has been violated — regardless of how clean the interface looks.

**What the embryo contains** (all of it operations that need no Gröbner basis):

| capability | how, without the engine |
|---|---|
| `M.quotient(H)` | catalog lookup; `raise` naming `M` and `H` when absent |
| `M.contains(points)` | a membership predicate per member (`‖Ω‖ = 1` for `S²`) — `[M]` this alone refuses the ERR-080 forged measure AT CONSTRUCTION |
| `M × N` | the product manifold; already spelled by `measure.__mul__`'s support rule |
| `chart : M → N` | a typed arrow (§III.9); makes pushforward's codomain derived |
| the singular stratum | `H`'s fixed-point set, carried by the quotient it produces |
| the group lattice | ✅ **already ships** — `SubgroupOfO3.is_subgroup_of`, and Part IV's G0 verdict table runs on it today |

⟹ the honest statement of what is NOT built: **no orbit space is ever COMPUTED**
from a group and a manifold. Everything above is bookkeeping over derivations a
human did once, plus predicates. `plan-authoring` §2: the refusal has a
denominator — `[M]` ~12 entries, enumerated in §III.1/1.1.

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

**D0.7 — ⭐ The manifold is a TYPE. `Manifold` is minted.** (User ruling,
2026-08-31, after the naming discussion below.) Today `Space = str` with
`measure.py:114` calling the tags *"recommendations, not constraints"* (L8) —
so nothing can express containment, quotient, membership, or composition, and
the "check" degenerates to string equality.

*Against the project's own minting criterion* (`coding-standards`, Type vs
property — mint **iff** ≥2 non-isomorphic realizations AND a non-identity
morphism is applied):
* **(a)** `S²`, `S¹`, `[-1,1]`, `[0,1]`, `[0,∞)`, `ℝ`, `S²/σ_y`, `spatial_Rᵈ`,
  `energy`, `index(axis)` — non-isomorphic, and `[M]` all shipped.
* **(b)** quotient (`measure.quotient`, shipped), product (`product()`'s
  Archimedes composition, shipped), pushforward (`measure.pushforward`,
  shipped). §III.9's support algebra IS the morphism list.

Both hold ⟹ earned, not speculative.

⚠ **The prior ruling this must not contradict, read rather than paraphrased**
(`plan-authoring` §1's PRECEDENT clause): `sn_reshape.md` Issue 2 —
*"Don't try to enforce `Space` types via Python generics — not expressive
enough without runtime overhead"*, quoted verbatim at `measure.py:106-109`.
`[R]` **it does not cover this.** It rejects `Generic[Tag]` **phantom type
parameters** — which `coding-standards` also rejects, since they are erased at
runtime and do not specialize dunders. A first-class **value** with real methods
(`quotient`, `__mul__`, `contains`, `dim`) is a different proposal, and the one
the morphisms above require.

⭐ **The NAME, checked against the prose corpus** (§1's clause: a free name can
be free *because it was rejected*). `[M]` `Manifold` — **2** hits in
`.claude/`, neither a rejection; **22** in `docs/theory/` as ordinary
mathematics; and it is already the corpus's own word for exactly this concept
(`measure.py:120` — *"The circle names the MANIFOLD, not a chart of it"*).
⛔ `Domain` was considered and REFUSED: `[M]` **137** prose hits, and it is
already the operator-algebra's TYPE vocabulary
(`LinearOperator[Codomain, Domain]`) — minting it for manifolds is a genuine
collision, unlike the *attribute* `domain`, which is not (see §III.10).

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
| **L1** | `angular_frame` **fabricates** `support=SPACE_SPHERE` onto interval nodes | `directional.py:624` (⛔ **was `:587-593`; re-measured 2026-09-01** — `angular_frame` now begins at `:580` and the fabrication sits at `:624`); rule support `[-1,1]` → frame support `S^2`; `‖node‖` 0.1834..0.9603 | `Basis` cannot declare its **domain**, so the binder forges a matching support |
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

⛔ **REFUTED 2026-08-31 (§II.15 R2/R3) — measured, and the traffic is ZERO.**
`[M]` **0 of 145 089** two-sided composability checks put an `L2[S^2…]` space on
either compared side; a positive control STRONGER than 0.1c reddens exactly one
test, on an *unfolded* rule; 0.1c itself reds nothing. The static chain above is
real code and **carries no traffic**. The paragraph stays because its *shape* is
the transferable lesson: §8's grep-the-branch check is a WIRING check and cannot
tell load-bearing from inert, and the error direction — over-caution — reads as
rigor.

`[M]` **10 of 12 rows: the forged nodes are BIT-IDENTICAL to the rule's own**
(`np.array_equal`, exact shapes). The two exceptions are the 1-D rules, where
the rule's nodes are scalars `(N,)` and the forgery is `(N,3)`.

⛔ **INCOMPLETE — corrected at 0.1a's execution, 2026-09-01. The forgery
destroys THREE truths, not one.** This section and §II.4's leak inventory both
enumerate only the **support tag**, because that is the field the campaign is
*about*. `[M]` by `dataclasses.fields(DiscreteMeasure)` — the type carries five
fields, and the rebuild passed three:

| field | rules carrying it | frames carrying it (pre-0.1a) | what it gates |
|---|---|---|---|
| `support` | 12 of 12 | 12 of 12, **4 FALSIFIED** — 2 slab (`[-1,1]`) + 2 fold (`S^2/sigma_y`), i.e. two distinct tags over four rules | `measure_space`, i.e. an operator face's `domain` |
| `invariance_group` | **10** of 12 | **0** of 12 | `DiscreteMeasure.phase` (`measure.py:409`) |
| `exactness` | **10** of 12 | **0** of 12 | `degree_of_exactness`, `generating_measure` |

⭐⭐ **The consequence nobody had looked at: `frame.measure.phase` RAISED
`NotImplementedError` on ALL 12 rules — the ANGULAR frame's own measure could
not say it was angular.** `phase` keys on `invariance_group` first and falls
through to `support.startswith("spatial"|"energy")`, and `'S^2'` matches
neither. `[M]` `lebedev(17).measure.phase == 'angular'` while
`lebedev(17).angular_frame(2).measure.phase` raised
*"phase is undetermined for support 'S^2'"*.

⚠ And the gate that should have caught it names the claim in its **docstring**
and asserts something weaker in its **body**:
`test_phase_is_a_closed_category_consistent_with_angular_frame`
(`tests/numerics/test_measure_phase.py:96`) says *"the angular quadrature's
frame and its measure both say angular"* and then asserts
`frame.measure.support == SPACE_SPHERE` — the support TAG, not the phase.
`[R]` the phase assertion would have raised, so the weaker one is what
survived. 0.1a makes the docstring's own claim assertable, and the strengthened
row lands with it (`coding-elegance` anti-pattern #20 — a docstring naming a
claim the body does not make).

⟹ **the transferable half** (`plan-authoring` §2's quantifier clause, applied
to a STRUCTURAL enumeration rather than a numeric one): *"the rebuild loses X"*
is a completeness claim over the **source type's field list**, and its
denominator is `dataclasses.fields(T)` — not the concept you are chasing. I
enumerated against the concept and got 1 of 3. One line answers it.

⭐ **`.phase` has ZERO production consumers** (`[M]` by AST over `orpheus/` +
`tests/`: 6 reads, all in `test_measure_phase.py`; positive control
`.weights` → 106 in `orpheus/`). So this is a **correctness** repair to an
authored property, not a live behaviour change — §8's own 2026-08-31
sharpening, and this time I priced it before writing the warning.

⭐ **A claim 0.1a DEMOTES, found by self-review and closed in the same commit.**
`angular_frame`'s docstring asserts *"``frame.table`` equals
``spherical_harmonics(L)`` bit-identically (both route through
``SphericalHarmonicBasis.evaluate`` **on these cosines**)"*. The parenthetical
was the *mechanism*: both spellings literally shared one
`column_stack(axis_cosines(0..2))` expression, so agreement was true **by
construction**. After 0.1a the frame routes `self.measure.nodes` while
`spherical_harmonics` still column-stacks, so the claim now holds **because two
arrays happen to be equal** — strictly weaker, and `coding-standards`' silent
claim-class demotion with *no gate at all*, only prose.
`[M]` re-measured across 12 rules × L ∈ {0,1,2}: **36 of 36 bit-identical**,
with a positive control that reads `False`. `[M]` and **nothing pinned it** —
0 tests compare the two. ⟹ the docstring's reason is corrected and the claim
gains the gate it always deserved (`test_q8_6`). It is not a tautology:
both sides call the same basis object, and independence lives in the **input
assembly** (measure nodes vs column-stacked cosines) — the
[[feedback-verify-shared-primitive-pure-math]] criterion.

⚠ **Collateral, and it belongs to 2.0c — not to 0.1a.** After 0.1a the phase
answers for 8 of 12 rules and still raises for 4: the two 1-D rules (whose
frame measure is ERR-080's fiction, correct to keep raising until 3.4) **and
the two FOLDED rules**, whose own measure carries `invariance_group=None`.
That last one is *deliberate and gated* — `test_measure.py:353` and `:949` pin
it, and rightly: a σ_y-quotient of an `O_h`-invariant measure is not
`O_h`-invariant, so inheriting the parent's group would be a false claim.
The real gap is that `phase`'s fallback is a **string-prefix** test
(`support.startswith("spatial"|"energy")`), which `'S^2/sigma_y'` cannot match
— stringly-typed dispatch, exactly what **2.0c** retires by making `support` a
`Manifold`. ⟹ **2.0c acceptance item: `folded_product(4,8).measure.phase ==
"angular"`, derived from the manifold rather than from a prefix.**

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
fold.

⛔ **The paragraph that followed here is REFUTED by 2.0c, and the refutation is
the useful part.** It read: *"`folded_by` … is structurally unavailable for a
continuous quotient. The tree holds two quotient records — a group object
(wired, discrete) and a string tag (inert, continuous) — and the slab is the
unwired half."* `[M]` 2026-09-01 at `8cc69e7f`: there is **no string tag any
more**, and `Quotient.by` is a group object that handles the *continuous* case
identically — `SPHERE.quotient(SubgroupOfO3.SO2).by` **is** `SubgroupOfO3.SO2`
(spelled `SO2("x")` since 2.4 — the group names its axis),
and `DiscreteMeasure.quotient_group` derives from it for discrete and continuous
folds alike. So the "two records" are now ONE mechanism, and the structural
obstacle this paragraph named **does not exist**.

✅ **What survives, and it is 2.4's whole content: the slab is still the unwired
half.** Not because the machinery is missing — it is present and exercised — but
because the slab rule declares `Interval(-1,1)` instead of the quotient. That is
now a one-line declaration against a working catalogue, which is why 2.4 stopped
being a research question (§V.5g). ✅ **WIRED at 2.4 (`17501245`)** — and the
"one line" was three: the group had to name its axis first, and the
declaration had to live on the slab's ADOPTER rather than on the chart-level
rule the product factor shares (2.4 EXECUTED).



## II.14 ✅ The full-suite phantom census — DISCHARGED 2026-08-31

§II.5's ⚠ denominator warning is now closed. Instrument:
`scratch/_phantom_census_plugin.py`, a **non-invasive recorder** (it records
the `axis_index >= dim` reads and returns the normal zeros, so it cannot
redden anything — a census that could fail would be measuring the plugin).
Positive-control validated before use: it fires at exactly
`directional.py:589`, `axis=[1,2] dim=[1] via=angular_frame`, with the 50
tests of `test_frame.py` still passing.

Run: all 13 trees, `-O -m "not slow" -p no:randomly`, serial, **rc=0
throughout** (numerics 2538 · transport 645 · sn 3375 · geometry 727 ·
diffusion · data · homogeneous · cross_method · moc · mc 39 · cp ·
derivations 1637 · rootfiles 414), 65 min wall-clock.

**`[M]` 11 phantom call sites, 1604 reads. 5 in `orpheus/`, 6 in `tests/`.**

| site | reads | the question it asks | verdict |
|---|---|---|---|
| `numerics/quadrature/directional.py:589` (`angular_frame`) | **1478** | *"what is the direction cosine along axis i?"* | ⛔ **FABRICATION** — no answer exists for `i ≥ dim` |
| `numerics/quadrature/directional.py:330` (`mu_z`) | 31 | reached from `spaces/angular_trace_space.py:220` (`_quadrature_axis` → `Ω·n̂_f`) — *"how much flows along i?"* | ✅ **zero is correct** |
| `numerics/quadrature/directional.py:320` (`mu_y`) | 29 | — | `[M]` **0 production outer-frames**; only tests reach it |
| `geometry/boundary/_specular.py:162` | 5 | the cosine measure `w·\|μ_a\|` on a face | ✅ **zero is correct** |
| `transport/mesh/axis.py:441` | 4 | `face_outflow_ordinates` | ✅ **zero is correct** |

⟹ **§II.5's "ONE production consumer" is CONFIRMED — but only for the
FABRICATION class, and only after triaging by MEANING.** The raw site count is
5, and a reader who stopped at the count would either believe the claim was
5× understated or, worse, act on it: `[M]` **0.2's blanket-`raise` means would
have broken three legitimate consumers**, all three asking the flow question the
accessor's own docstring documents. §II.5's number was right and its predicate
was doing all the work.

⭐ **And the sequencing this settles: 0.1 comes FIRST and shrinks 0.2 to a
clarification.** 0.1a deletes the `column_stack` — so the fabrication read
disappears with it, taking **92 %** of all phantom reads. What remains is three
consumers asking a legitimate question through an accessor named for a different
one. 0.2 is then a pure naming/typing split with no defect behind it, and its
done-when is re-runnable: **re-run this census after 0.1 and the
`directional.py:589` row must be GONE.**

⛔⛔ **REFUTED 2026-09-01 at 0.2's opener, and it is MY row — 0.1a removes
`0 %` of the phantom reads, not 92 %.** A phantom read is `axis_cosines(i)`
with `i >= dim`, which **can only happen when `dim == 1`**: a 3-D rule's
indices run 0..2 against `dim = 3` and never qualify. 0.1a routes exactly the
**3-D** rules — the ones that were never phantom-reading — so it cannot have
removed any. `[M]` measured with the census's own recorder design at
`2c1a06b1`: `lebedev(11)` / `product(4,8)` / `folded_product(4,8)` read **0**
phantom components, before and after; `gauss_legendre(8)` still reads **6**
(axes 1 and 2 × `L ∈ {0,1,2}`), merely relocated `:589 → :682`.
⟹ the 92 % belongs to **0.1b**, which rides **3.4** — the FUSED ORDER step (renumbered to #5 at the 2026-09-01 compaction), not step 1.
⭐ The mechanism of the mistake is worth more than the number: I wrote "0.1
deletes the column_stack" and then split 0.1 into 0.1a/0.1c (**now**) and 0.1b
(**at 3.4**) *in a later section*, without returning to re-attribute the 92 %.
A quantity attributed to a step **survives that step being split**, and it
attaches to whichever half the reader is looking at. ⟹ when a step splits,
re-read every number credited to its old name.

⭐⭐ **And the sharpest thing 0.2's opener found — the forged node is the orbit
BARYCENTRE, exactly.** `[M]` 2026-09-01, against an explicit
periodic-trapezoid quadrature of the `SO(2)` orbit of each
`gauss_legendre(8)` ordinate: the mean of the circle of unit vectors at polar
cosine :math:`\mu` is :math:`(\mu, 0, 0)` to **1.1e-16** in the constant
column and **1.7e-17** in the two averaged ones.

⛔ **This row first quoted "0.000e+00 on all three columns", and that figure
measured the wrong pair.** The probe had compared the FORGERY to
:math:`(\mu,0,0)`, which is exact by construction — those columns are literal
zeros — so it carried no information. The informative comparison is the ORBIT
MEAN against the forgery, whose precision is the reference quadrature's own.
Same shape as the 2026-08-31 transcription-precision row: the number that read
as the strongest evidence was the one wholly determined by how the probe was
written. ⭐ And the correction turned up a second thing worth carrying:
**refining the orbit rule makes it WORSE.** An equispaced rule IS the periodic
trapezoid, which is exact for trigonometric polynomials, so `n = 8` already
integrates :math:`\cos\varphi` / :math:`\sin\varphi` to machine precision;
`[M]` `|col0 − μ|` reads `1.1e-16` at n=8, `9.5e-15` at n=1024 and `1.1e-12`
at n=200 000 — pure summation round-off in the constant column. The gate's
first draft used 200 000 points and failed its own `1e-12` tolerance against
its own reference, not against the code.

⟹ **ERR-080's numbers were never wrong; its TYPE is.** The zeros are not
garbage and not a fallback — they are the correct value of
:math:`\langle \Omega_i \rangle` over the orbit. The defect is that a **mean**
was handed to a basis that needs a **point**, and that is precisely why
:math:`\lVert(\mu,0,0)\rVert = |\mu| \ne 1`: the barycentre of a circle of
unit vectors is its centre, which lies strictly inside the ball. It is the one
point of the convex hull guaranteed *not* to be on the orbit (except at
:math:`\mu = \pm 1`, where the orbit degenerates to a point).

⭐ This is the same verdict the σ_y derivation reached by a different route —
*"ERR-080 is a botched SECTION"* — now with its mechanism named: a section must
choose a **representative on the orbit**, and the barycentre is the canonical
way to fail that. It also settles 0.2's design question, because it says what
the flow accessor actually returns: **the orbit mean**, whose zero is DERIVED
rather than defaulted, and which coincides with the plain cosine whenever the
orbit is a single point.

### ▶ 2.1 — PRE-FLIGHT, measured 2026-09-01 at `a2befd9b` (before any edit)

Run so the next session inherits checked premises, not this section's prose.

**✅ The defect REPRODUCES, exactly as claimed.**
`IndicatorBasis((array([0.,1.,2.]),))` — read as 2 ENERGY groups — and
`IndicatorBasis((array([0.,0.5,1.]),))` — 2 SPATIAL cells — both mint
`FunctionSpace('L2[coarse_cells_R1]', shape=(2,))`, which is `==`-EQUAL **and
hash-equal**. `[M]` negative control (3 cells vs 2) correctly `False`, so the
equality is not a degenerate always-true. This is an illegal state that IS
representable, and `FunctionSpace.__eq__` is `(name, shape)` — so the false
name is doing all the damage.

**✅ The §6b sizing is RIGHT, for a reason the plan does not state.**
`[M]` `IndicatorBasis(` construction sites: **18 — 4 in `orpheus/`
(`data/energy_grid.py:220`, `geometry/mesh.py:444`, `:753`,
`numerics/frame.py:755`), 14 in `tests/`** — matching the plan exactly. ⚠ But
the indicator FAMILY is **39** sites (11 `orpheus/`, 28 `tests/`), and that
does not widen the step, because:

| class | relation | does a ctor field on `IndicatorBasis` reach it? |
|---|---|---|
| `OverlapBasis` | **subclass** of `IndicatorBasis` | ✅ by inheritance |
| `WeightedIndicatorBasis` | ⛔ **SIBLING** — direct subclass of `Basis` | ✅ anyway, by **COMPOSITION**: it holds an `indicator` field and derives its space from it |

`[M]` verified both ways: `WeightedIndicatorBasis` wrapping the energy and the
spatial indicator gives the SAME `L2[coarse_cells_R1]` pair, `==` and hash
equal — it inherits the defect *through the field*, so it inherits the repair
through the field too. ⟹ **18 sites is the right §6b set for the repair.**

**⚠ What 18 does NOT cover — the DONE-WHEN is wider than the defect.** The row
promises *"every shipped `Basis` answers; `SphericalHarmonicBasis.domain is
SPACE_SPHERE`"*. `[M]` by runtime introspection the shipped set is **6**
(`IndicatorBasis`, `OverlapBasis`, `WeightedIndicatorBasis`,
`SphericalHarmonicBasis`, `MirrorEvenSphericalHarmonicBasis`,
`LossKernelBasis`), and **0 of 6** defined `domain` when this was measured
(✅ 6 of 6 since `c461fe8d`). Three get it from the
indicator repair; the other three need their own answer, and two of those are
easy (`SphericalHarmonicBasis` → the sphere; `MirrorEven…` → its quotient,
which `Manifold` can now spell as `S^2/sigma_y`). **`LossKernelBasis` is the
one to look at first** — it lives in `sn/operators/loss_kernel_gauge.py`, not
in the basis package, and nothing in the plan says what manifold it is on.

⚠ **And the done-when's own spelling needs re-deriving before it is built.**
It says `SphericalHarmonicBasis.domain is SPACE_SPHERE` — but `SPACE_SPHERE`
is the **`str`** `"S^2"` (`measure.py:130`), while 2.1's stated type is
`Manifold`. After 2.0a the sphere object is `Sphere()`. Those are different
things and an `is` comparison against the string would either fail or, worse,
pass while typing the field as `str` and quietly re-introducing `Space = str`
— the exact thing 2.0c retires. ⟹ **decide the target type first**; the row's
`is SPACE_SPHERE` is pre-`Manifold` text that survived the mint.

### ▶ 2.1 — the DESIGN, measured 2026-09-01 at `ba05f773` (before any edit)

**Goal (`plan-authoring` §5, the domain's terms).** A basis function is a map,
and a map is not defined until you say what it EATS. Every basis states the
manifold its functions are defined ON, so that a space cannot claim an identity
it does not have.

⭐⭐ **The sharpest statement of the defect, found by reading the four
construction sites: the measure and the basis are built in the SAME FUNCTION,
three to five lines apart, and the measure names the manifold CORRECTLY while
the basis fabricates a spatial one.** Nothing can see the disagreement because
one of the two is a hard-coded f-string. `[M]` at `ba05f773`:

| site | the MEASURE's `support`, lines away | the BASIS's invented name | the `Manifold` |
|---|---|---|---|
| `frame.py` `:759` / `:755` | `f"index({axis_label})"` | `L2[coarse_cells_R1]` | `IndexSet(axis_label, n)` |
| `energy_grid.py` `:206` / `:220` | `"energy"` | `L2[coarse_cells_R1]` | `EnergyGroups(n_groups)` |
| `geometry/mesh.py` `:418` / `:444` | `"spatial_R1"` | `L2[coarse_cells_R1]` | `RealSpace(1)` |
| `geometry/mesh.py` `:731` / `:753` | `"spatial_R2"` | `L2[coarse_cells_R2]` | `RealSpace(2)` |
| `loss_kernel_gauge.py` `:1169` / `:1164` | `f"sn_trace_orbit{orbit}_g{group}"` | `loss_kernel_coeff[orbit..._g...]` | `IndexSet(f"sn_trace_orbit...")` |

⭐ The last row is the **positive control**: the one basis that does NOT
fabricate is the one whose author named its space by the block's own identity.
The defect is not "indicator bases are careless" — it is that four of five
sites had the answer in scope and re-invented it.

`[M]` the live collision, over the four production spaces:
`('L2[coarse_cells_R1]', (2,))` is minted by **two different manifolds**
(energy and spatial-1D); the `frame.py` one escapes only on `shape`.

**Proposed means, and the two things injection REFUTED before any edit.**

⛔ **(a) `domain` CANNOT be a dataclass field satisfying an abstract property.**
The obvious shape — `Basis.domain` abstract, `IndicatorBasis.domain` a ctor
field — does not run. `ABCMeta` re-checks `getattr(cls, name)` against the
base's `__abstractmethods__`, and an annotation-only field puts **nothing** in
`__dict__`, so the class stays abstract: `[M]` injected,
`TypeError: Can't instantiate abstract class ... without an implementation for
abstract method 'domain'`. Note it is order-dependent and would have "worked"
had the field carried a default -- i.e. the shape that also re-admits `None`.
The tree's own PRECEDENT (`plan-authoring` §1 -- read, not recalled) is
unanimous the other way: `[M]` across **4 ABCs and ~40 subclasses**
(`Basis.space`, `FrameBase.test`, `LinearOperator.domain` x17,
`Manifold.dim`/`name` x9) **every** abstract property is answered by a
`@property`, never by a field -- including `Manifold` itself, where `Interval`
stores `a`/`b` and exposes `name`. => the field carries the basis's own
vocabulary (`partition_of`) and `domain` is the uniform question, derived. That
is not two spellings of one thing: *what these edges partition* and *what these
functions eat* are the same manifold **by a theorem about indicator bases**,
and the property is where that theorem is written down.

⛔ **(b) The attractive STRONGER invariant would reject a shipped production
site.** `domain.contains(cell_centres)` reads as the better construction check
than a bare dimension match -- and `[M]` it REFUSES `frame.py:755`, whose
single-region partition `[-0.5, n-0.5]` has centre `(n-1)/2`, not an integer,
on an `IndexSet` whose `contains` requires integers. The honest invariant is
the ambient-dimension match (a `d`-axis tensor partition partitions a manifold
whose points carry `d` coordinates), which is `vv-principles` #16: do not
assert tighter than the type promises.

**The §6b set is 22 construction sites + 1 factory forward -- NOT the
pre-flight's 18.** `[M]` by AST over 895 files, unclipped:
`IndicatorBasis(` **18** (4 `orpheus/` + 14 `tests/`; **7 positional**),
`OverlapBasis(` **3** (all `tests/`), `OverlapBasis.from_indicator` **1**
(`energy_grid.py:270`, forwards `edges_per_axis=` and must forward the new
field too). ⚠ The pre-flight's table was right that `OverlapBasis` inherits the
FIELD; it never asked whether `OverlapBasis` has construction sites **of its
own**, and `OverlapBasis(` contains no `IndicatorBasis` substring -- the
2026-08-24 / 2026-08-29 §6b row, "members spelled without the symbol", fourth
instance in this campaign. `WeightedIndicatorBasis` (18 sites) and
`SphericalHarmonicBasis` (25) need **no** ctor change: both derive.

**Name-pin blast radius, `[M]` validated filter with a positive control:**
`coarse_cells` is **1 production line, 0 test pins, 4 doc lines**
(`manifolds.rst:468,480`, `spaces.rst:1900,3547` -- all four quote the false
name as a live example, so all four are Cardinal-Rule-3 MUST-FIX).

**Done-when** -- ⛔ NOT the tracker row's `SphericalHarmonicBasis.domain is
SPACE_SPHERE`, which is pre-`Manifold` text: `[M]` `SPACE_SPHERE` is the `str`
`"S^2"`, and `Sphere() == SPHERE` is `True` while `Sphere() is SPHERE` is
**`False`**, so the row's spelling would either fail or pass while typing the
field `str`. The predicate is:

* every one of the 6 shipped `Basis` subclasses answers `domain` with a
  `Manifold` (structural -- the ABC refuses to instantiate otherwise);
* the energy and spatial 2-cell indicator spaces are **not** `==` and **not**
  hash-equal, with the 3-cell negative control still unequal;
* `SphericalHarmonicBasis(L).domain == SPHERE` and
  `MirrorEvenSphericalHarmonicBasis(L, axis).domain == SPHERE.quotient(Mirror(axis))`
  (`[M]` names `S^2` and `S^2/sigma_y`, matching 0.1c's shipped fold tag);
* an `IndicatorBasis` whose partition rank disagrees with its manifold's
  ambient dimension cannot be constructed.

⚠ **What 2.1 does NOT do, deliberately** -- `DiscreteMeasure.support` stays a
`str`, so `f"L2[{support}]"` (`measure.py:331`) and the basis's new
`f"L2[coarse_cells({domain.name})]"` remain two spellings of "the L2 name over
a manifold". That twin is **pre-existing** (§III.10's two-site table) and is
inventory item #10, owned by **2.0c** (`FunctionSpace.manifold`); 2.1 makes the
second one HONEST without hoisting. Recorded here so the next reader does not
read the twin as new.

---

### ✅ 2.1 — EXECUTED 2026-09-01, landed `c461fe8d`. What building it actually found

**What shipped.** `Basis.domain` is an abstract `@property` returning a
`Manifold`; all 6 shipped subclasses answer, and the ABC refuses to
instantiate one that does not. `IndicatorBasis` takes a required
`partition_of: Manifold` and derives `domain` from it; its coefficient-space
name becomes `f"L2[coarse_cells({domain.name})]"`. `SphericalHarmonicBasis`
answers `SPHERE`; `MirrorEven...` answers `SPHERE.quotient(Mirror(axis))`
(**derived**, so it is the same object the folded measure carries — `[M]`
`S^2/sigma_y`, matching 0.1c); `WeightedIndicatorBasis` and `OverlapBasis`
delegate; `LossKernelBasis` answers its block's `IndexSet`.
`manifold._ambient` was promoted to the public **`ambient_dim`** — the
construction invariant is its first cross-module consumer.

⭐⭐ **The type did work on its FIRST DAY: it forced a distinction the string
tag could not make.** Assigning manifolds to the 17 test sites was supposed to
be transcription. `[M]` two `tests/data/test_energy_grid.py` sites partition
by **energy VALUE in eV** (`[1.0e7, 1.0e3, 1.0e-5]`) while
`tests/data/test_mixture_condense.py` and production `EnergyGrid.as_basis()`
partition by **group INDEX** (`arange(n+1) - 0.5`). Those are two different
manifolds — the continuous half-line and a counting set — and the old
vocabulary called both of them `"energy"`. They are now `HALF_LINE` and
`EnergyGroups(n)`, and `ambient_dim` would not have caught the confusion
(both are 1): only naming the point set does. ⟹ the mint's value is not
only that it makes the ERR-class unspellable; it is that **assigning the type
is itself a census** that reads every site's actual meaning.

⛔ **The §6b set was 23, not the pre-flight's 18 — and the 5 extras are three
distinct spellings, all invisible to a census of the changed class's name.**

| the 5 | why a `IndicatorBasis(` census cannot return it |
|---|---|
| `OverlapBasis(` ×3 (`test_frame.py:427/508/788`) | a SUBCLASS's own construction sites — the pre-flight's table correctly said `OverlapBasis` inherits the FIELD, and never asked whether it has constructors of its own |
| `OverlapBasis.from_indicator` (`overlap_basis.py:108`) | a FACTORY that forwards `edges_per_axis=` and must now forward `partition_of=` |
| `_DenseTrial((edges,))` (`test_frame.py:462`) | a TEST-LOCAL subclass, on the SAME LINE as a counted `IndicatorBasis(` call |

⟹ fifth instance of `plan-authoring` §6b's *"members spelled without the
symbol"* in this campaign, with a spelling not in the 2026-08-24/29
inventory: **a subclass's own constructor**. Cost: none — all three classes
fail loudly (`TypeError: missing required positional argument`), and the
sites were enumerated before any edit rather than by a red loop.

⛔ **Two designs REFUTED by injection before any edit** (both recorded in the
design block above): a dataclass FIELD cannot satisfy an abstract property,
and `domain.contains(cell_centres)` — the invariant that reads stronger —
**refuses a correct production caller**.

⚠ **One known divergence, PINNED rather than deferred silently.**
`LossKernelBasis`'s measure tags the bare label `sn_trace_orbit…_g…` while
`IndexSet` spells the same point set `index(sn_trace_orbit…_g…)`. `test_d6`
asserts both spellings explicitly, so **2.0c** (which retypes `support` to a
`Manifold`) must return here and cannot resolve it by accident. The other
four production pairs agree exactly today.

⭐ **`coding-elegance` #20 caught in my OWN gate, before it landed.** `test_d6`'s
docstring said the loss-kernel divergence *"is asserted below"* and the body
did not assert it — the exact anti-pattern (a docstring naming a claim the
body does not make) that §"the rebuild loses X" row logged three days ago,
committed by the author of that row. Found by re-reading the file before
installing it; the fix was the assertion, not the sentence.

**Gates** — `tests/numerics/test_basis_domain.py`, 13 rows, all `foundation`:
`d1` the defect keystone (energy vs spatial spaces separate, with a
negative control), `d2` the completeness claim by **runtime** `__subclasses__`
(landmine #4; filtered to `orpheus.`-defined classes so a test-local subclass
cannot make it import-order-dependent), `d3` the ABC's structural refusal
(the stub implements every OTHER abstract member, so the refusal cannot be
credited to the wrong arm — `vv-principles` #17), `d4` the invariant both
legs, `d5`/`d5b` the harmonic domains swept over 4 degrees × 3 mirror axes,
⭐ `d6` **the flagship — a frame's two halves name ONE manifold**, and `d7`
the delegation drift guard asserted by identity.

⭐ Why `d6` is the keystone rather than a name assertion: a gate on the basis's
own name is satisfied by any self-consistent lie, and the pre-2.1 name WAS
self-consistent. `d6` compares two independently-authored halves of one frame,
so it reddens the day either side re-invents the manifold — which is the
failure that actually happened.

**The 7-arm mutation battery — every arm a DISTINCT red set** (in-process
monkeypatch, so there is no file to leave mutated if the harness kills the run;
`scratch/_p21_mut/`). Run over 203 rows across 8 suites:

| arm | calls | reds | which |
|---|---|:--:|---|
| *(none — baseline)* | 0 | **0** | — |
| `m1` `domain` guessed from `ndim` (the original lie, one level up) | 102 | 3 | `d1`, `d6`, **`d7`** |
| `m2` `space` name hard-coded to the pre-2.1 literal | 97 | 2 | `d1`, `d6` |
| `m3` construction invariant removed | 163 | 1 | `d4` |
| `m4` the fold forgets its quotient (answers `SPHERE`) | 3 | 3 | `d5b` ×3 |
| `m5` the wrapper answers for itself instead of delegating | 1 | 1 | `d7` |
| `m6` harmonics eat a circle | 7 | 7 | `d5` ×4, `d5b` ×3 |
| ⭐ **POSITIVE CONTROL** — every basis eats a sphere | 105 | 6 | `d1`, `d4`, `d5b` ×3, `d6` |

⭐ `m1` and `m2` separate, which is the discrimination that matters: a wrong
NAME (`m2`) and a wrong DOMAIN (`m1`) are different failures, and only `m1`
reds `d7`, because a wrapper faithfully delegating a wrong answer is a third
thing again.

⚠ **Two gates no arm reddens, stated rather than left to read as unwitnessed.**
`d3` builds BOTH its legs itself (an abstract stub that must refuse and a
completed one that must construct), so it is its own witness and no production
mutation can reach it. `d2` is a **forward** guard — all 6 shipped bases answer
today, so it is unfalsifiable now by construction; its value is that a basis
added tomorrow is in scope without editing the file. Neither is a coverage gap;
both would be, if this note were not here.

⛔ **pyright caught 8 errors of mine, and the fix taught something about the
`type: ignore` rule.** `test_d3`'s abstract stub was first written with bare
`...` bodies and untyped signatures — so all 7 members "override `Basis` in an
incompatible manner", plus `space` returning `None`. Retyped with full
signatures and `raise NotImplementedError` bodies (a stub that silently returns
`None` is not a faithful stand-in for the ABC, and the checker says so).

⭐ That left **one** error, and it is the one the test EXISTS to assert:
`Cannot instantiate abstract class "_Incomplete" — "Basis.domain" is not
implemented`. I had removed the `# type: ignore[abstract]` on
`coding-elegance` #19 grounds ("no `type: ignore`") — wrong here, and the tree
settles it: `[M]` **228** `type: ignore` comments ship, and the dominant idiom
is exactly this one, a negative leg spelling the illegal thing inside
`pytest.raises` (`TimedFullField(interior="not a field", ...)`
`# type: ignore[arg-type]`). ⟹ #19 targets an ignore that **hides** a defect;
an ignore on a line whose whole purpose is to be refused **documents the
claim**, and removing it makes the file dirty rather than honest. Restored with
that reasoning written beside it. `[M]` the file is then `0 errors, 0 warnings`.

⚠ **And a scope note on this campaign's pyright figures, so a later reader does
not compare two populations** (`plan-authoring` §2). `[M]` a bare `npx pyright`
over the whole project reads **2533 errors** — it includes
`derivations/diagnostics/`, `tools/`, and other trees carrying long-standing
pre-existing errors. The campaign's "pyright 0" claims are DELTA claims on the
files a step touches, which is the right predicate and was never written down.
2.1's honest figure: **0 new errors across the 19 files it touches**; 15
pre-existing errors sit in 4 of those files at lines this work does not touch
(`test_mixture_condense.py`, `test_axis_marginal.py`, `test_frame.py`,
`test_harmonic_frame.py`).

⛔ **The battery's FIRST run was broken, and the call counter is the only reason
that was legible.** Every arm reported `0 failed` — including the positive
control, which cannot be right. `[M]` the cause: **zsh does not word-split an
unquoted scalar**, so `$SUITES` reached pytest as ONE nonexistent path and 0
tests were collected — `rc=4`, §II.16(c)'s recorded signature, third instance in
this campaign. Without `CALLS`, "0 reds" reads as *the gates are blind* and
would have sent me to strengthen gates that were already fine
(`vv-principles` #17: the harness lies before the code does, in the
safe-looking direction). Fixed with a zsh array + an `rc` column.

### ✅ 2.1b — EXECUTED 2026-09-01 `9b4a4d9c`. What building it actually found

**What shipped.** `Basis.invariance_group -> SubgroupOfO3 | None` — a
**concrete, `@final`, DERIVED** property on the ABC (`base.py:291`), right
after the abstract `domain` (`:247`). A `match` on the TYPE of `self.domain`:
`Quotient(base=Sphere(), by=group)` → `group`; `Sphere()` → `Trivial`;
anything else → `None`. `MirrorEvenSphericalHarmonicBasis` answers
`Mirror(axis)`, `SphericalHarmonicBasis` answers `Trivial` at every `L`
(including `L = 0`, which is `O(3)`-invariant — the reading is the LOWER
BOUND the domain guarantees, and under-declaring is legal but lossy), the
three indicator bases and `LossKernelBasis` answer `None`. Module-scope
imports of `Manifold, Quotient, Sphere` and `SubgroupOfO3` into `base.py`,
injected and run in **8** fresh import orders (§6d) — no cycle. Zero
subclass edits; zero test-local subclass edits.

⭐⭐ **The FOURTEENTH consecutive phase opener to correct its own section, and
the second to dissolve a requested FIELD into a derivation.** The tracker
read *"`Basis.invariance_group` — absent; derivable for every shipped
basis"* and §V.5h(b) measured *"0 of 6 subclasses answer it"*. Both true —
and the answer was already in the tree: `[M]`
`MirrorEvenSphericalHarmonicBasis(L=2, mirror_axis=1).domain.by` is
`Mirror('y')`, the SAME object `folded_product(4, 8).measure.quotient_group`
returns, because `domain` is `SPHERE.quotient(Mirror(axis))` and the quotient
is memoised. A function on `M/H` IS an `H`-invariant function, so for
FUNCTIONS "has `H`" and "descends to `M/H`" are one property, and the group a
basis HAS is the group its domain SPENT. Six overrides would have been six
second homes for a fact the domain already states — exactly the shape 2.0d
(*a `quotient_group` FIELD*) had when it dissolved into `Quotient.by` at
2.0c. ⟹ 2.1b landed as a **no-op extension through the single generic body**
(`coding-standards`, clean-before-extend), which is what made its §6b set
empty.

**The category ruling (§V.5h(b)'s open question): `None` off a non-angular
domain, `Trivial` off the bare sphere.** Three answers, not two. `Trivial` on
`RealSpace`/`EnergyGroups`/`IndexSet` would assert an `O(3)` action that
does not exist; the measure's `phase` gives the same category answer for the
same manifolds; and the sphere is NOT `None` because `O(3)` acts on it and
the full harmonics have spent none of it — *"the invariance group of the
degree-`L` real SH basis is trivial"* is the textbook statement, and `None`
would read as *not applicable*, which is false for an `S²`-basis. §8's
branchable-`None` hazard does not bite: G2 (2.2) must special-case the
MEASURE's `None` (spent nothing) regardless, so the basis's `None` costs the
predicate nothing new. `@final` because an override would let `domain` and
this answer disagree — the two-tags-that-drift state the derivation exists to
make unspellable; the remedy for a basis invariant under more than its domain
shows is to declare the finer domain (3.4-R: the TRUE `LegendreBasis` on
`S²/SO2_x`, not an `m = 0` sub-basis on `S²`).

**HAS vs SPENT, settled for the basis side.** A measure carries two group
slots because for a POINT SET they come apart — a rule can be closed under a
mirror (HAS) without having been folded by it (SPENT), and once folded it is
no longer closed. For functions the two coincide, so a basis carries ONE slot,
named for what it HAS (`invariance_group`, matching Part IV's predicate
`measure.quotient_group ⊆ basis.invariance_group`) and read off what its
domain SPENT. The two properties' `None`s therefore mean different things —
`quotient_group` `None` = spent nothing (Lebedev), `invariance_group` `None`
= no `O(3)` action (energy) — and both docstrings say so.

`[M]` **The pairings the tree ships, at `a4166ef2` + this step:**

| pairing | SPENT (`measure.quotient_group`) | HAS (`basis.invariance_group`) | `has.contains(spent)` |
|---|---|---|---|
| `folded_product(4,8)` vs its frame basis | `Mirror('y')` | `Mirror('y')` | **True** — and `has is spent`, `basis.domain is measure.support`: ONE object |
| `gauss_legendre(8)` (`S^2/SO2_x`) vs the full-sphere SH its frame binds today | `SO2('x')` | `Trivial` | **False** — ERR-080's pairing as a lattice verdict, for the first time |
| `lebedev(17)` vs SH | `None` | `Trivial` | spent nothing ⟹ admitted |
| `product(4,8)` vs SH | `None` | `Trivial` | same |

⚠ **What is NOT claimed.** Nothing REFUSES on the slab's verdict: the frame
gate is 2.2 (fused with 0.1b + 0.6 + 3.4), and `[M]`
`gauss_legendre(8).angular_frame(2).measure.support.name == 'S^2'` (the
forged one) while `gauss_legendre(8).measure.support.name == 'S^2/SO2_x'` —
§V.5h(c) re-confirmed, which is why `test_e1` reads the verdict off the
QUADRATURE's measure. ERR-080 stays OPEN (the strict-xfail gate still reads 3
xfailed). G0's form (equality vs ambient-base + lattice) is still the ruling
§V.5h(a) owes 2.2's opener.

**Gates** — `tests/numerics/test_basis_domain.py` section E, **11 rows**
(13 → 24 in the file), all `foundation`: ⭐ `e1` **the keystone — the fold's
two halves read ONE group object, and the slab's pairing is refusable**
(positive leg by lattice, equality AND identity; negative leg the ERR-080
verdict, stated as a measurement not a refusal, §6c); `e2` the full-sphere
harmonics HAVE `Trivial` over `L ∈ {0,1,3,7}`; `e2b` mirror-even HAS its
mirror **by identity** (`is basis.domain.by` — `[M]` `Mirror('y') is
Mirror('y')` is False for two constructions, so `is` separates *derived*
from *duplicated*), axis load-bearing; `e3` the category leg with a positive
control in the same test; `e4` the three arms on ONE class shape differing
only in `domain` (separates *reads the domain* from *knows the class*); `e5`
Part IV's four-row table on real objects (rows 2–3 via a stub declaring
`S²/SO2_x`, retirement trigger for 3.4). `test_d3`'s stub was hoisted to
module scope (`_AbstractStub`) so section E can declare domains on it.

**The 8-arm battery** (in-process monkeypatch, `scratch/_p21b_mut/`; **156
rows across 7 suites** — `test_basis_domain`, `test_spherical_harmonic_basis`,
`test_frame`, `test_measure_phase`, `test_slab_orbit_space`,
`test_quadrature_fold`, `test_harmonic_frame`):

| arm | calls | reds | which |
|---|---|:--:|---|
| *(none)* | 0 | **0** | — |
| `m1` always `None` | 19 | 11 | `e1`, `e2`×4, `e2b`×3, `e3`, `e4`, `e5` |
| `m2` always `Trivial` | 13 | 7 | `e1`, `e2b`×3, `e3`, `e4`, `e5` |
| `m3` the sphere arm answers `None` (the two-way design the opener rejected) | 29 | 8 | `e1`, `e2`×4, `e3`, `e4`, `e5` |
| `m4` the quotient arm ignores the base | 30 | **0** | ⚠ **BLIND BY CONSTRUCTION** — no non-sphere quotient ships, so no input reaches the difference; predicted before the run and stated in the arm |
| ⭐ `m5` the fold basis STORES a copy (equal, not identical) | 8 | 4 | `e1`, `e2b`×3 — the derived-vs-duplicated discriminator, and only the identity legs see it |
| `m6` the fold forgets its quotient | 8 | 8 | `d5b`×3, `e1`, `e2b`×3, `e5` |
| `m7` the quotient arm loses the axis | 21 | 6 | `e1`, `e2b`×3, `e4`, `e5` |
| **POSITIVE CONTROL** — every basis HAS `SO2('z')` | 13 | 11 | `e1`, `e2`×4, `e2b`×3, `e3`, `e4`, `e5` |

⚠ Read honestly: `m1` and the control red the SAME 11 (both wrong
everywhere), so "every arm distinct" is false for that pair; the
discriminating arms are `m2`/`m3`/`m5`/`m7`, which separate *always trivial*,
*two-way*, *copied* and *axis-blind* from one another. And **no arm reddens
anything outside `test_basis_domain.py`** — expected and worth stating: the
property has no consumer yet (§6c), so its only witnesses are its own, and the
first production consumer (2.2's G2) must bring its own.

**Seven-tree gate: not re-run, and why.** numerics **2690** rc=0
(2679 + 11, reconciled unit-for-unit against section E); the slice most likely
to notice — `test_quadrature_fold` + `transport/frames` + `sn/angular`, **79**
rows — rc=0; root docs gates (`test_docstring_xrefs`,
`test_error_catalogue_reconciles`, `test_vv_harness_audit`,
`test_layer_imports`, `test_pending_ports`) **415 + 5 xfail** rc=0 (⚠ a
different file set from 2.4's *"50 passed"*). data/geometry/transport/sn NOT
re-run: 2.1b changes no behaviour any tree outside `numerics` reads — a
derived property with zero consumers, and the battery reddens nothing outside
its own file. pyright **0** on the 4 touched files (⚠ it caught 2 of mine:
pyright does not narrow a PROPERTY re-read across statements — `assert x.g is
not None; x.g.contains(...)` still reads `Optional`; bind to a local).

⭐ **Three present-tense-false docstrings fixed on sight** (the articulation
standard's enforcement half), all residue of earlier steps and all OUTSIDE
this step's diff until it touched the file: (1) `test_basis_domain.py`'s own
module docstring — *"`DiscreteMeasure.support` is still a `str` (that retype
is tracker 2.0c)"*, false since `025834f5`; the file's `test_d6` docstring had
been re-tensed at 2.0c while the module header was not; (2) `measure.py:222`
— *"`SubgroupOfO3.SO2` for an azimuthally-symmetric product/slab rule"*,
false twice (a bare `SO2` is landmine #6 since 2.4; product rules carry
`Dnh(n_phi)`, the slab's marginal `Mirror('x')` — ERR-072); (3)
`symmetry.py:450` — a comment listing `SubgroupOfO3.SO2` among the ClassVar
singletons. ⚠ (2) and (3) sit in files 2.4's own sweep edited; its filter
looked for CODE spellings of the bare name and both survivors are PROSE
(`coding-standards`' tense-discrimination clause — grep the retired name
across prose too, and sort by tense).

**Docs** — the archivist's pass is recorded in the commit (manifolds §(b)
continuation + the inventory row; discrete_measures HAS/SPENT for the basis
side; ERR-080 dated status; matrix regenerated).

### ✅ 2.3 — EXECUTED 2026-09-02 `5ec3a00a`. What building it actually found

**What shipped.** `ManifoldMap` — ONE frozen value type in `manifold.py`
(`domain: Manifold`, `codomain: Manifold`, `apply`), with `__call__` on
`(n, d)` arrays and `@` as composition refused across mismatched endpoints;
two memoised named maps beside it, `archimedes(axis)` (`[-1,1] × S¹ → S²`,
:eq:`product-mu-phi-cosines` verbatim for `"z"`) and `barycentre(orbit_space)`
(`S²/SO(2)_a → D³`, `μ ↦ μ ê_a`); `DiscreteMeasure.pushforward(phi: ManifoldMap)`
with `new_space=` RETIRED — the support is READ off `phi.codomain` and a map
out of the wrong point set is refused; `quotient()` builds its orbit
retraction as a typed map; `spherical_product` IS
`(polar * azimuthal).pushforward(archimedes("z")).with_metadata(...)` — the
hand loop and the `support=SPHERE` literal (L7) are gone; `_embedded_nodes`'s
polar-marginal arm READS `barycentre`; `symmetry._AXIS_INDEX` moved to
`manifold.AXIS_INDEX` (the chart needed it and the runtime graph already runs
`symmetry → manifold`). The forgery arm (`directional.py:786`) is untouched
except a comment naming it as `barycentre` with a forged codomain.

⭐ **The SIXTEENTH consecutive phase opener — and the first whose section
survived it intact.** §V.5i(a)–(f) re-ran at `f2682869` and reproduced 6 of 6
(the production tree was identical to the pre-flight's; only the plan had
moved). What execution corrected was the row's **name** and one **home**:
the tracker's *"the typed `Chart`"* names a mechanism a differential geometer
would refuse — a chart is `M ⊃ U → ℝⁿ`, and of the three maps only the
*inverse* of the Archimedes one is a chart (the retraction lands on a
`Quotient`, the barycentre on a `Ball`) — ruled `ManifoldMap`
(`feedback_high_signal_names`: the precise term, Domain-first); and the
axis table the chart reads lived in `symmetry.py`, which `manifold` cannot
import at module scope, so it moved to the module whose curved members own
the convention. Two more found on sight, neither the opener's: the
`DiscreteMeasure` class docstring still spelled the product and quotient
supports as f-strings (`f"{μ.support} × {ν.support}"`, present-tense-false
since 2.0c — the far-from-the-diff tense survivor, third time this
campaign), and my own first draft of the domain-guard message offered
`on_orbit_space` as the remedy for a map out of the wrong product, which
runs the other way (chart → orbit space); re-worded before any test pinned
it.

⚠ And one of MINE the archivist refuted on landing: the `ManifoldMap`
docstring (and §V.5i(b)) said `Ball` had *"0 production consumers"* — the
census behind it ran `grep -v manifold.py` and so could not see `Ball(2)`,
the σ_y entry's `realization` at `manifold.py:922`. Corrected to what was
measured (0 consumers OUTSIDE its module; `Ball(3)` never constructed) in
the docstring, §V.5i(b) and the theory page. Three more present-tense-false
spots the archivist found and I fixed in the same landing:
`axis_cosines`'s docstring (`directional.py:298`) still promised the zero
fallback phase 0.2 retired for a `raise`; the ERR-080 gate module's
docstring counted *"two"* strict-xfail rows where `[M]` by AST there are
three; and the level-membership comment in `rules_product.py` kept the 8-space
indent of the loop it no longer lives in.

**The four rulings (§V.5i(g)), all taken at the opener:**

| question | ruled | why |
|---|---|---|
| shape | ONE frozen value type; named maps are FACTORIES (`archimedes`, `barycentre`), the way `SPHERE`/`LEGENDRE` are values | three one-line maps do not earn three classes; the retraction is data-dependent (built per certificate) so it is a constructor call either way |
| scope | all THREE maps in this step | the retraction is the ONE production `pushforward` caller (free once the signature changed); the barycentre is the Pattern-2 collapse §V.5i(f) named — and `Ball(3)`'s first production consumer |
| the pushforward REFERENCE measure (the `AngularSymmetry.reference` twin) | **DEFERRED to 3.1**, mechanism recorded | `[M]` 2026-09-02 a module-scope `manifold → exactness` import kills **5 of 5** fresh import orders (`ImportError: cannot import name 'Manifold' from partially initialized module`); the viable home is a `Quotient` field populated INSIDE the derivation function via function-scope import — 3.1's row already lists "pushforward measure" among the entry's fields |
| name | `ManifoldMap` | above |

⭐ **ERR-080, restated as a property of a MAP.** `barycentre`'s image satisfies
`1 − ‖μ ê_a‖² = 1 − μ² = det P / 4` — the squared orbit radius the catalogue
records as `Quotient.det_gram` — so it lies ON the sphere exactly on the
singular stratum (the poles) and INSIDE the ball everywhere else. The
forgery is this map with `codomain=SPHERE`; `test_barycentre_lands_in_the_ball…`
asserts `not SPHERE.contains(image[1:-1])` as a fact about the map, and the
forgery arm's comment now says which map it is. Through `pushforward` the
sentence has no spelling; through the raw constructor it still does, and
0.1b (fused step) closes that door — **2.3 did not move the ERR-080 gate
(still 3 xfailed) and was not expected to (§V.5i(b), §6c).**

`[M]` **Derivation agreements, at the carve** (`scratch/_p23_bitid.py`,
pristine copies loaded under alias module names):

| object | configurations | differ |
|---|---|---|
| `spherical_product` vs the pre-2.3 loop — every `dataclasses.fields` of the measure AND of the `LevelStructure` | **60** (n_μ ∈ {2,3,4,5,6} × n_φ ∈ {6,8,10,16,24,32} × {NODE_ALIGNED, STAGGERED}) | **0** |
| `Quadrature.folded_product` (the retraction path) vs pristine product + fold, incl. `support is` the memoised `S^2/sigma_y` | 5 ((2,8),(4,8),(4,16),(6,10),(3,24)) | **0** |
| `_embedded_nodes` vs pristine — 12 orbit-space rows (n ∈ {2,4,8,16} × x,y,z) + bare interval, circle, product | 15 | **0** |

Squaring, measured before spelling the chart: on every shipped Gauss node set
(n ∈ {2,3,4,5,6,8,12,16}) the scalar `mu_val**2` of the old loop and the
array `mu**2` of the chart are bit-equal (56 of 56 nodes; on *random*
doubles they are not, so this was a measurement, not an assumption).

**Gates** (`scratch/_p23_driver.log`; serial, `-O`, `-m "not slow"`):
numerics **2690 → 2711** (+21: `TestManifoldMap` 14 rows, `test_rules_product`
+7, `test_measure` net 0 — the REFUSES gate re-posed as three legs) ·
sn primitives+mesh+angular + transport/frames **579** rc=0 · sn/operators
**1274** + 5 xfailed rc=0 · sn/sweep **919 + 1 skipped + 31 xfailed rc=0** · geometry + _harness +
mms-ordering **736 + 4 skipped + 1 xfailed rc=0** · pyright **0** on the 8 touched files · 8 files
compile under `-W error` · 6 import orders clean.

**Mutation battery** (`scratch/_p23_mut/mutplugin23.py`, in-process, counted;
526 rows over 10 suites; names bound by `from … import` patched at
BOTH homes):

| arm | what it breaks | reds | reading |
|---|---|---|---|
| none | — | 0 / 526 | baseline |
| m1 `archimedes` left-handed (b ↔ c) | the chart's frame | **26** in 7 files | product bit-identity, group derivation, level structure, admission, the fold |
| m2 `archimedes` codomain = its own domain | the derived support | **58** in 8 files | everything that reads `support` downstream of `product()` — the widest arm |
| m3 `pushforward` domain guard OFF | the refusal | **2** | exactly the two refusal legs (`test_measure`, `test_rules_product`) — the guard's ONLY witnesses, both constructed (§6c) |
| m4 `pushforward` keeps its OWN support (ignores the codomain) | support READ off the map | **57** in 8 files | the fold lands on `S²` instead of `S²/σ_y`; the product on `[-1,1] × S¹` |
| m5 `barycentre` ignores the axis (always x) | the axis read | **7** | the y/z rows of three suites; x is blind BY CONSTRUCTION (the fixed axis is x) |
| m6 `barycentre` codomain = `SPHERE` (the forgery) | the honest codomain | **3** | ⚠ ONLY `TestManifoldMap`'s three barycentre rows — nothing else in the tree reads `Ball(3)` off a measure yet; the honest codomain's sole pin is the new gate |
| m7 `@` guard OFF | composition endpoints | **1** | the composition test |
| m8 `@` composes in the wrong order | functoriality | **2** | composition + functoriality |
| m9 `_embedded_nodes` pads by hand (the pre-2.3 body) | Pattern-2 collapse | **0**, CALLS 1152 | **BLIND by construction, predicted**: the hand padding and the map agree bit-exactly, so no value-level gate can tell a duplicated spelling from a derived one |
| control: `archimedes` negates μ | — | **9** | ⚠ a WEAK control, and the reason is `vv-principles` Mode 12: μ → −μ is σ_h, in the product rule's own invariance group, so the node SET is unchanged and only the 9 order-sensitive pins moved. m2/m4 (58/57) are the effective positive controls |

⚠ **The harness lied first, in the safe-looking direction** (`vv-principles` #17): the driver's inline `$(...)` capture returned EMPTY summaries for m1, m2, m4 and the control — the four widest arms — which reads as *"the arm did not run"*; re-run to log FILES they are the 26/58/57/9 above. Read from a file, never from a shell capture, when an arm can redden fifty gates.

**§6b, as executed.** `pushforward`'s set was the pre-flight's 1 + 8 exactly
(⚠ AFTER the step it has **2** production callers — `quotient()` and
`spherical_product` — so do not re-quote the 1)
— the red loop's first pass reddened the 6 tests holding those 8 sites and
nothing else across 489 tests of the directly affected suites; no
spelled-without-the-symbol member surfaced (the string form `"pushforward"`
has no other consumer). `spherical_product`'s 20 callers (3 production, 17
test) needed no edit: the signature is unchanged and the output bit-identical.

**Docs** (archivist, reviewed): `manifolds.rst` +926/−42 — a new `=`-level
section `manifold-arrows` with six subsections (the type; the three arrows the
tree was already spelling; Archimedes; the barycentre; composition and the
labelled `manifold-map-functoriality`; what 2.3 did NOT build — neither of
the entry's own chart/section ships), Key Facts, status YAML, seams,
Development-history row, and re-derived — not transcribed — `[M]` rows
including `π ∘ archimedes_a = pr₁` bit-exactly (500 random points per axis)
and the shipped fold as `retraction ∘ archimedes` on 5 of 5;
`discrete_measures.rst` (pushforward row reads `φ.codomain`; two stale
f-string cells from 2.0c fixed); `spherical_harmonics.rst` (its ERR-080
mechanism paragraph described the pre-0.2 `column_stack` — re-tensed with the
old text kept); `error_catalog.rst` ERR-080 dated ✅ Progress block, stays
OPEN, gate described as unmoved; `matrix.rst` regenerated **10595 → 10616**.
`sphinx -E -W` 0/0/0; `dead_references` **0 / 52**; root docs gates **415 +
5 xfailed** rc=0.

**What is NOT claimed.** ERR-080 stays OPEN. The forgery arm is a raw
constructor by design until 3.4. The reference measure still lives on the
registry (3.1). `archimedes` is consumed in production for `"z"` only — the
x/y arms are exercised by tests alone. `pushforward` still runs no membership
check on its image (0.1b's job, at construction).

**Residue for 3.1 (recorded here so the deferral is a work item, not a
loss):** a `Quotient.reference` field, typed under `TYPE_CHECKING`, populated
by `_sphere_mod_so2` with a function-scope `LEGENDRE` import and by
`_mod_trivial` with `UNIFORM_ON_SPHERE`; `AngularSymmetry.reference` then
reads `self.support.reference` — and the σ_y entry ships `None`, because its
pushforward is a *weighted* disk measure (`2 dx dz / √(1−x²−z²)`) no shipped
`ReferenceMeasure` realization can spell, which is honest: nothing consumes
it (the fold drops `exactness`).

### ✅ 3.1 — EXECUTED 2026-09-02 `67e38605`. What building it actually found

**Outcome.** An orbit-space catalogue entry carries ALL SEVEN of D0.1's
derivation outputs. `[M]` `dataclasses.fields(Quotient)` = **12** (was 10):
`orbit_coordinates` (required; the quotient map's action on the base's
ambient coordinates, returning the realization's coordinates — the
invariants that SURVIVE the base's ideal) and `reference` (`ReferenceMeasure
| None`, default `None`). The typed arrow is the DERIVED property
`Quotient.quotient_map` — `ManifoldMap(base, self, orbit_coordinates)`,
codomain the ENTRY, never the realization. `AngularSymmetry.reference`
READS the entry (the same collapse `support` underwent at 2.4); its
`LEGENDRE` import is gone. The opener (the **eighteenth**) is the second
whose pre-flight (§V.5j) survived intact: (a)–(e) reproduced 5 of 5 at
`3623adc2`, with two user rulings on NAME and MECHANISM and one on the
two honest `None`s.

**The three rulings (AskUserQuestion, 2026-09-02).**

| fork | ruled | why |
|---|---|---|
| the map's name + codomain | **`quotient_map`, codomain the ENTRY** (not D0.1's "chart" into the realization) | a chart is `M ⊃ U → ℝⁿ` (2.3's ruling) and `Ω ↦ Ω·ê_a` is not injective on `S²`; the exit-path consumer (3.4b's `Descent`) pulls back along π; landing on the realization is exactly the reading 2.4 made refusable |
| how it is carried | **stored apply + derived property** | a frozen dataclass cannot hold a map whose codomain is itself; the 2.1b pattern (`Basis.invariance_group` derived from `domain`); `field(compare=False, repr=False)` — a function has no value equality and the map is a function of `(base, by)` |
| the reference on the two non-axial entries | **`None` on both; the registry keeps its bare-sphere arm** | σ_a's pushforward is the WEIGHTED disk measure `2 dx dz/√(1−x²−z²)` no shipped `ReferenceMeasure` spells; `M/{e}`'s is Lebesgue on the BASE, whose orthogonal system `Manifold` does not carry; the registry's Trivial support is the bare `SPHERE`, never a quotient |

**Measured at execution** (`scratch/_p31_{census,witness,inject}.py`, logs beside them):

| claim | `[M]` |
|---|---|
| `π_a ∘ archimedes_a = pr₁` | bit-exact **12 of 12** (3 axes × n ∈ {2,4,8,16}, random φ), through the typed `@` |
| `barycentre_a ∘ π_a` = axial projection | bit-exact **3 of 3**, 1000 directions; `Q.contains(π(v))` |
| π is H-invariant (its DEFINING law) | bit-exact on **7 of 7** entries (3 rotations each about the axis; the reflection; the identity); negative leg: a rotation about ANOTHER axis moves it — ⚠ chosen, since σ_b (b ≠ a) does NOT move π_a |
| the numeric map IS the recorded symbolic invariants | `lambdify` of the surviving generators == `orbit_coordinates`, **7 of 7** entries (the D0.1 tie; the roster is the six keys + the trivial one) |
| pushforward identity on a REAL rule | `level_symmetric_sn(4)` along π_x: support **is** the entry, `∫μ² d(π_*μ) == ∫(Ω·x)² dμ` bit-exact, ≈ 4π/3; a rule on `[-1,1]` REFUSED |
| §6d, function-scope import, SHADOW copy | **alive 7 of 7** fresh orders; the same line at module scope **dead 7 of 7** (`cannot import name 'DiscreteMeasure' from partially initialized module 'orpheus.numerics.measure'` — the cycle closes through `measure`, which `generating_measure` imports) |
| import-time quotient construction | **0** module-scope `.quotient()`/`.on_orbit_space()` calls in `orpheus/` (positive control 7 anywhere) |
| §6b for the required field | 3 derivation returns + **1** test constructor (as pre-flighted); 3 `match Quotient(...)` patterns untouched (`basis/base.py:361` is a third the pre-flight did not list — a non-member for a keyword field) |
| the registry twin's readers | **1** production (`registry.py:1302`, the V stage); 2 tests (`test_registry.py:345`, `test_slab_orbit_space.py:206`) — both still green by IDENTITY |
| pickle | a `Quotient` round-trips equal (`partial` of a module-level function, not a lambda) |

**Gates** (serial `-O -m "not slow" -p no:randomly`, all rc=0, logs
`scratch/_p31_gate_*.log`): numerics **2711 → 2739** (+28: 21 new `test_manifold.py` rows + 1
`test_registry.py` row in the first pass, then +6 after the archivist's
findings — 3 invariance rows and 2 symbolic rows on the widened 7-entry
roster, and the module-scope mint gate; unit-for-unit; the 8 warnings are
2.3's 8, unchanged) · ERR-080 gate **2 passed, 3 xfailed**
(strict — OPEN) · sn primitives+mesh+angular+frames **579** · geometry +
harness + mms-ordering **736** + 4 sk + 1 xf · root docs gates
(`test_error_catalogue_reconciles` + `test_docstring_xrefs` +
`test_layer_imports`) **399**, re-run green after the docs pass (⚠ a different file set from 2.3's "415 + 5
xf" — not comparable) · sn/operators **1274 + 5 xf (identical to 2.3)** · sn/sweep
**919 + 1 sk + 31 xf (identical to 2.3; 4 min 39 s)**. pyright **0** on the 4 touched files; 4 files compile under
`-W error`. `sphinx -E -W` **0 warnings / 0 errors, rc=0**, run twice (before and after the post-archivist test edits; `scratch/_p31_sphinx_final.log`); `dead_references` **0 dead / 0 sites / 52 checked, 52 rescued** on the freshly rebuilt graph; matrix
**10616 → 10644** (+28 unit-for-unit: `foundation` 7592 → 7620, `numerics/test_manifold` 70 → 97, `numerics/test_registry` 74 → 75; the new documented label `manifold-quotient-pushforward` takes the sentinel count 577 → 578).

**The battery** (`scratch/_p31_mut/mutplugin31.py`, in-process, per arm, a
CALL counter; the memo verified EMPTY at configure; 10 files, 612 rows (final tree, after the post-archivist test edits);
per-arm logs `scratch/_p31_mut/_arm_*.log` — read from files, never a
shell capture):

| arm | reds | CALLS | reads as |
|---|---|---|---|
| none | 0 / 612 | 0 | the baseline |
| m1 SO2 `orbit_coordinates` reads the NEXT column | **14** | 31 | every new law: π∘φ ×3, β∘π ×3, invariance SO2 ×3, symbolic tie ×3, the pushforward row, the pickle row |
| m2 `quotient_map` codomain = the REALIZATION | **11** | 18 | the arrow-into-entry ×3, the two compositions ×6 (refused at the typed `@`), the pushforward row, the replace/pickle row |
| m3 axial entries carry NO `reference` | **18** (16 + 2 harness) | 6 | ⭐ **the registry READS the entry**: 14 registry/slab-orbit reds (the V stage raises for slab AND sphere) + the manifold reference row; ⚠ 2 are artefacts of wrapping the catalogue dict (`test_all_three_mirrors_share_ONE_derivation` pins builder identity; `test_d2` pins the memo) |
| m4 the registry TABULATES `LEGENDRE` again (the old twin) | **1** | 48 | ⭐ **the rewrite is behaviour-preserving on every shipped geometry BY DESIGN** — the tabulated answer coincides with the entry's; the one red is the `None` arm's missing-work message (§6c: the witness the gate landed with). Read with m3: together they say *the registry reads the entry, and the entry says what the table said* |
| m5 the quotient arm answers the SPHERE measure | **14** | 39 | the V stage refuses Gauss–Legendre for slab/sphere |
| m6 σ_a `orbit_coordinates` KEEP the dropped column | **7** (6 + 1 harness) | 52 | all three mirrors' invariance rows (the reflection now moves the image) + all three symbolic ties; ⚠ the pickle row reds because `pickle` resolves the patched module function by NAME |
| m7 trivial `orbit_coordinates` are ZERO | **2** | 4 | the Trivial symbolic tie + the invariance row's NEGATIVE leg (an outsider no longer moves zeros) |
| m8 `quotient_map` domain = the REALIZATION (expected weak) | **10** | 20 | stronger than expected: both compositions refuse at the typed `@`, the pushforward's domain guard refuses, the arrow row's `domain is SPHERE` |
| POSITIVE CONTROL every quotient named `S^2/Trivial` | **13** | 76 | across 6 files, none of them 3.1's gates — a control OUTSIDE the SUT's invariance group (vv-principles #17 as sharpened at 2.3) |

Nine of nine mutation arms red with distinct red sets; no arm blind; the
memo verified EMPTY at `pytest_configure` so every entry in every run was
built by the patched builder.

**§6b as executed.** The required field reached exactly the four sites the
pre-flight named. The registry's `reference` rewrite reached its one
production reader and two test readers with **no edit** — both compare by
identity and the entry's field IS the object the twin used to tabulate.

**NOT claimed.** The section (`fundamental_domain=None`) stays unchosen —
a section is a CHOICE and 3.1 shipped the quotient map, which is
derivation output, instead. ERR-080 is not repaired; the forgery arm is
untouched; the probe is 3.4's; `Descent` does not exist (`class Descent`
0); no weighted-disk `ReferenceMeasure` was added (the registry's raise
names it as the missing WORK); the `symmetry.py` package re-home is off
the path. The 13-tree gate baseline is still unmeasured (best spent before
the fused step).

**Residue for the fused step (⭐ 0.1b + 0.6 + 2.2 + 3.4 + 3.4b).** The
quotient map is what `Descent` pulls `P_ℓ` back along; `Quotient.reference`
is what a claim on the descended basis is against. 2.2 still owes its G0
ruling (§V.5h(a)).

**Docs** (archivist, reviewed): archivist, reviewed — `manifolds.rst` +~1200 (status YAML; a new Key Facts bullet; the engine-data-model table's two *not a slot* rows tombstoned into SLOTS, `9 of 9` over twelve fields; a new `~` subsection `manifold-second-twin-reference` with a re-measured 5-row table; arrows-not-built (1) half-remedied — the map ships, the SECTION does not — and (2) DISCHARGED with the 7/7-vs-0/7 measurement beside 2.3's 5/5 (a different edge); three new sections `manifold-quotient-map`, `manifold-pushforward-reference` (new documented label `manifold-quotient-pushforward`, deriving the hat-box AND the σ_a weighted disk from scratch, mass-checked to 4π symbolically and against `lebedev(11)`) and `manifold-value-at-function-scope`; the import-cycle section retitled MODULE scope and its edge table widened to five modules with four stale line numbers fixed; seams row; changelog row), `discrete_measures.rst` +58, `error_catalog.rst` +55 (ERR-080 Progress block, entry OPEN). **It refuted the brief four times**: *pickle round-trip equal* holds for the ENTRY (7/7) and not the callable (1/7 — `partial` has identity equality; the field's `compare=False` is what makes the entry equal); `.reference` reads are 10 not 9 (the new read is the tenth); the third `match Quotient` is `ambient_dim`, not `barycentre` (which uses `isinstance`); `:964` is the guard, `:965` the read. It widened two claims (symbolic tie 7 of 7 entries, not 5; `∫μ²` is 1 ULP from 4π/3) and found one mathematical fact: π_a is ALSO invariant under σ_b, b ≠ a — `O(2)_a` and `SO(2)_a` induce the same orbit partition, so a quotient map fixes the PARTITION and the partition does not fix the group; `Quotient.by` is a declaration. Four code findings, all acted on before landing: a future-tense word in `archimedes`'s docstring (fixed); `reference.support == realization` had no pin (now gated on the 7-entry roster — it is the chart codomain, never the entry); the function-scope import's safety condition — 0 module-scope quotient mints — had no gate (now an AST gate with a positive control); and the σ_y derivation's algebra of record is an UNTRACKED memo (filed as an issue — the corpus's derivation-source debt, its weakest dimension at 3/5)

---

### ✅ 2.5 — EXECUTED 2026-09-02 `58ad7c28`. What building it actually found

**The outcome.** The angular moment space has ONE home — the basis the
quadrature's frame bound — and every consumer READS it. `TruncatedBasis`
(`numerics/basis/base.py`, a `runtime_checkable` Protocol: `L` + `space`) is
the harmonic family's shared surface; `HarmonicFrame`'s two doors ask for
THAT surface (`_admit_truncated`, message *"truncation order"*), not for one
class, and the mints read `frame.truncation_order`; the Λ ends read
`self.flux_analysis.frame.basis.space` (`ScatteringOperator._moment_space`,
was `_sh_space`), the two tier-2 classmethods take a `basis: TruncatedBasis`,
the fission and (n,2n) ℓ=0 ends read `self.frame.basis.space`
(`_sh_space_l0` deleted), the moment-flux field's head is
`mesh.quad.angular_frame(L).basis.space` (`_angular_head_space`, a
`_CarriesQuadrature` Protocol; a quadrature-less `MaterialMesh` refused
with a typed message), and `truncate` asks the head for ITS OWN family one
order down (`SphericalHarmonicSpace.truncated`, name-preserving) and slices
by the new head's shape — head-layout-agnostic, so B13's `truncate` row is
already discharged. Minted at the fused step's opener (§V.5l) and executed
the same session — no separate opener; the test-architect's memo
(`scratch/_fused_verification_plan.md`) was read before the first line.

**Two rulings I made (the memo's A-R1 / A-R2), recorded here because a plan
must not hide a decision in a diff:**

1. **A-R1 — the ends bind `basis.space`, the CONTINUUM space, not the
   frame's Parseval-dressed `basis_space`.** `[M]` the fork, re-derived by
   the archivist draw-free (the test-architect's *"`apply_metric` moves
   96–161 %"* did NOT reproduce under any norm — a relayed number, §4):
   the two spaces are `(name, shape)`-equal and metric-DIFFERENT on **33 of
   33** (rule, L) rows (11 constructions from all 5 factories × L ∈
   {0,1,2}); the per-ℓ metric ratio is exactly `[(2ℓ+1)/4π]²` =
   `6.33e-3 / 5.70e-2 / 1.58e-1`; and Λ's Hilbert adjoint under the
   CONTINUUM end is `Λ* = Λᵀ` **exactly (0.0) on 33/33**, while under the
   dressed end it MOVES on 10/33 rows — every `gauss_legendre` /
   `folded_product` row at L ≥ 1 (`9.7e-2 … 1.37e-1` on the DIAGONAL-verdict
   ones, `1.08e-1 … 1.58` on the DENSE ones). ⛔ My first argument for the
   ruling — *"Λ is diagonal in ℓ so its adjoint agrees under either metric
   wherever the Gram is diagonal"* — is REFUTED by the six DIAGONAL movers
   and is struck. The ruling stands on the measurement: binding the basis's
   own space keeps every number and every `.H` bit-identical to the
   `from_L(L)` mint it replaces, which is A's whole acceptance. Gated with
   the dressed space as the NEGATIVE control (A2, and the m10 arm).
2. **A-R2 — `_space_for_mesh_and_L` RE-DERIVES from the frame, through the
   mesh's quadrature.** The memo's dichotomy (*each choice reds a
   population*) is a POST-B statement; with the SH basis bound everywhere
   both populations are green and bit-identical. `from_mesh_and_L(values,
   mesh, L)` keeps its signature — the mesh carries the quadrature, so it
   is frame-derived and not an `L`-mint; `[M]` 51 + 22 test constructions
   untouched.

**The §6b set, as executed** (predicate: every site that minted the space
or narrowed the basis, by the explorer's AST census + the red loop): 7
producers (`scattering.py` ×3, `fission.py`, `n2n.py`, `_bases.py`,
`harmonic_moment_flux.py`), 2 doors + 1 annotation (`harmonic_frame.py`),
the SH space's own `truncated` (new), the Protocol + export — 9 production
files; **26 test call sites in 8 files** re-keyed from `L=` to
`basis=SphericalHarmonicBasis(L=…)` (the explorer's grep read 28 LINES; two
were the same site counted twice), 1 refusal test re-pinned
(`test_harmonic_frame.py`: *"spherical-harmonic trial"* → *"truncation
order"*), 1 canary re-described + 1 external literal pin added. pyright
**0** on the 9 production files; **20** diagnostics on the 9 touched test
files, `[M]` **0 on a touched line** (intersected by script — all
pre-existing).

**Bit-identity, measured.** The slab GL8 fixed-source flux at
`L = 0, 1, 2` (`scratch/_fused_g0_census.py`'s fixture) is `array_equal`
before/after, `max|Δ| = 0.0` on all three (`scratch/_fused_baseline_slab_L01.npz`
vs `_fused_after_A_slab_L01.npz`); the sphere GL8 and folded(2,4)
cylinder solves reproduce their values; 21 frame constructions, 3 unique
pairings, unchanged. Every operator end and field head is
`(name, shape)`-equal AND metric-`array_equal` to `from_L(L)` on 11 shipped
rules × L ∈ {0,1,2} (gate A2, 33 rows).

**Gates** (`tests/transport/frames/test_moment_space_is_read_off_the_frame.py`,
36 rows, + 2 in `test_harmonic_frame.py`):

| gate | asserts | §6c witness it rejects |
|---|---|---|
| **A1** the ROUTE gate | a FOREIGN truncated basis (`_ForeignTruncatedBasis`: NOT an SH subclass, delegates the numerics, RENAMED coefficient space) bound into the quadrature's cache makes every operator end and the field head MOVE; `truncate` keeps the foreign name; the tier-2 mints take the basis | an end minted from `L` fails `OperatorProduct`'s guard `A.domain == B.codomain` — the negative leg, and every producer revert's red |
| **A2** metric identity (33 rows) + A2b on a posed composite | `frame.basis.space == from_L(L)` AND metric `array_equal`; the dressed `basis_space` is `==` and metric-different (`[M]` archivist: 33/33, ratio `[(2ℓ+1)/4π]²`) | the dressed space (the fork) |
| **A3** the door | an indicator trial refused at BOTH doors with *"truncation order"*; the foreign basis admitted, `truncation_order` read | the pre-2.5 `isinstance` narrowing |
| **A5** the demoted canary, re-described | both sides now derive from ONE source; keeps the discovery-path `is`-identity | — (a tautology named as one, per coding-standards) |
| **A6** the external pin | the moment codomain's shape as a hand-written literal `(3, 5, ng, nx)` | a wrong-`L` frame |

**Red loop.** Scope: `tests/transport`, `tests/numerics/{test_frame,test_basis_domain}.py`,
`tests/sn/operators`, `tests/sn/primitives`, the ERR-080 gate,
`test_dsa_rate.py`: **2403 passed, 2 failed** — both MY fixture defects in
the new gate (the foreign frame built on the RAW 1-D measure, which the SH
numerics refuse — it must ride the production frame's measure, today the
forged padding, tomorrow the rule's own; and a `scattering_order=2` request
against two-moment synthetic data), fixed; then **113 passed** on the five
gate files. The first launch died at COLLECTION with a `SyntaxError`: my
re-key script inserted the added import INSIDE a multi-line import block —
re-placed after the last top-level import by AST, 8 of 8 compile.

**Mutation battery** (`scratch/_fused_25_mut/mutplugin25.py` +
`_fused_25_battery.sh`, in-process, `pytest_configure`-time, a CALL COUNTER
per arm; scope 113 tests; `scratch/_fused_25_battery.log`):

| arm | mutant | reds | verdict |
|---|---|---|---|
| none | — | 0 / 113 | baseline |
| POSITIVE CONTROL | `SphericalHarmonicBasis.space` renamed | **34** (A2's 33 rows + 1), CALLS 109 | the harness sees |
| m1 `_moment_space` → `from_L` | | **1** — A1, CALLS 3 | caught |
| m2 fission ends → `from_L(0)` | | **1** — A1, CALLS 1 | caught |
| m3 (n,2n) ends → `from_L(0)` | | **1** — A1, CALLS 1 | caught |
| m4 field head → `from_L(L)` | | **1** — A1, CALLS 40 | caught |
| m5 `truncated` drops the name | | **1** — A1, CALLS 5 | caught (the name-preservation is what makes it visible) |
| m6 / m7 tier-2 mints → `from_L(basis.L)` | | **1** each — A1, CALLS 13 / 2 | caught |
| m8 / m9 a door narrows to SH again | | **3** each — A3, A1, the re-pinned door test, CALLS 174 / 173 | caught |
| m10 THE FORK: ends read `basis_space` | | **1** — A2b, CALLS 5 | caught — by the metric assertion only, as predicted: `==` is blind to it |

No arm blind; every mutant ran (CALLS > 0); integrity: only the 19 intended
files modified, nothing else in `orpheus/` or `tests/`.

**Gate** (13 trees, predicted first, `scratch/_fused_25_full_gate.sh` →
`_fused_25_full_gate.log`, serial canonical runner): predicted transport
646 → **683** (+36 in the new file + 1 in `test_harmonic_frame.py`), every
other tree unchanged, **10417 → 10454**. `[M]` **MEASURED, 13 of 13 rc=0:
10366 passed · 19 skipped · 69 xfailed · 0 failed = 10454, every tree equal
to its predicted denominator** — transport 682 + 1 sk = 683; the other
twelve trees bit-for-bit the pre-2.5 rows (numerics 2739 · geometry
727+4sk+1xf · data 237 · homogeneous 50 · diffusion 113 · cp 141 · moc 121
· mc 39+2xf · cross_method 81 · sn 3384+1sk+50xf · derivations
1637+13sk+11xf · root+harness 415+5xf); 61 min 55 s on a machine also
running the archivist. Delta against the pre-2.5 baseline: **+37 passed,
nothing else moved** — the 36 + 1 new gates, reconciled unit-for-unit.
⚠ The run started on the pre-docs tree and `root+harness` ran AFTER the
archivist's docs edits landed, so the docs gates tested the final corpus;
the five post-gate docstring corrections (the archivist's findings + my
refuted argument) are docstring-only and were re-run on the transport /
frame / primitives slices before the commit.

**Docs (archivist, same session).** `frame.rst` +452 — the SSOT section
`frame-moment-space-single-home` (eight-homes census, the fork, why the
continuum end is right for Λ, the `TruncatedBasis` door, the field head +
truncation, the three gates with their §6c witnesses; one new documented
label `moment-space-read-off-the-frame`, a Key Facts bullet, a
Development-history entry); `operator_algebra.rst` (the `from_galerkin`
guard's *"SH-only truncation order"* claim corrected + tombstoned);
`manifolds.rst` (status YAML + a paragraph under
`manifold-invariance-lower-bound`); `spaces.rst` (metric-blind identity's
gate consequence); `cartesian_multid.rst` (`two-moment-carrier-space`
scoped to the page's full-sphere setting); `error_catalog.rst` (ERR-080
Progress (2.5), entry OPEN); `matrix.rst` regenerated (sentinels 578 →
579). `sphinx -E -W` 0/0/0; `dead_references` 0 / 52; docstring-xref
checker 0. ⛔ **Three claims of my brief refuted by re-derivation**: the
mint census is 8 executable calls (7 re-mints + the basis's own), the fork
denominator is 33/33 not 12/12, and *"`apply_metric` moves 96–161 %"* did
not reproduce (see A-R1). Four docstring defects it could not edit (sources
frozen under the gate) — `harmonic_moment_source_sink.py:40-41` and `:113`
(the flux twin was re-tensed, the source/sink sibling was skipped), the
prose idiom `find_factor(SphericalHarmonicSpace).L` (5 sites, 0 production
calls), `_bases.py:628`'s two-axis shape contract, and
`directional.py:596`'s `-> SphericalHarmonicBasis` annotation (3.4's) —
applied after the gate, in the follow-up.

**NOT claimed.** No production number moved (by construction and by
measurement); no 1-D rule binds a Legendre basis yet; `HarmonicMomentFlux`'s
`[0, 0]` head reads (`isotropic_part`, `scalar_flux`, `anisotropic_part`,
`moment(ℓ)`) and `FissionMomentOperator`'s `m[0, 0]` are STILL two-axis
reads (the memo's B13 minus `truncate`) — they are the fused commit's, where
a rank-1 head first exists. `test_binding_tightness._equispaced_frame` is
untouched (G0 is the fused commit's).

**Residue for the fused commit** (so it is not re-derived): the foreign-basis
idiom in A1 is the shape of the Legendre gate's fixture (a non-SH basis
through the production chain); `_angular_head_space` and `_moment_space`
are the two sites a Legendre-bound frame flows through with no edit; the
`_TruncatesByOrder` Protocol is what the Legendre coefficient space must
satisfy (`L`, `truncated`, plus `TruncatedBasis`'s `space` on the basis).

### ✅ FUSED — EXECUTED 2026-09-02, landed `5436184e` (its gate and the typing narrowing rode the follow-up commit). 0.1b + 0.6 + 2.2-G0 + 3.4 + 3.4b: what building THE FIX actually found

**The outcome.** `solve_sn(scattering_order ≥ 2)` on a 1-D chart is RIGHT:
`[M]` LEG A (`scratch/_pl_slab_defect_repro.py`) reads `+4.000000000000` at
`L = 0, 1, 2` (rel `1.2e-15 / 4.4e-16 / 6.7e-16`; was `-3.764705882353` at
`L = 2`), and the ERR-080 gate's three strict-xfail rows (L = 2, L = 3,
GL16/L = 4 raising) XPASS. The slab's flux at `L = 0` and `L = 1` is
`array_equal` to the pre-fix baseline (`max|Δ| = 0.0`); every one of 135
sphere-rule arrays (9 rules × L ∈ 0..4: table, Gram, analysis) is
bit-identical to HEAD (a detached-worktree capture); the isotropic flux's
`ℓ ≥ 1` moments read ≤ `7.7e-16` on every shipped (rule, L) except six
labelled rows over three rules (`gauss_legendre(2)` L ≥ 4: `7.78e-01`, the
dead-slot theorem; `product(4,4)` / `folded_product(2,4)` L ≥ 4: `3.67` /
`4.89`, the two-pole aliasing under the harmonics, PRE-EXISTING) — the memo's
H-9 population exactly.

**What landed** (one commit, as XII.1b ruled; the memo's B-p1 … B-p10, plus
one §6b member the census had not listed):

| piece | where | the sentence |
|---|---|---|
| `LegendreBasis` | `numerics/basis/legendre_basis.py` | flat `(N, L+1)`; the MATCHED spelling (`1.0` / `μ` / `lpmv`) — `[M]` no single scipy routine reproduces the harmonics' m=0 column at the bit (H-1); dual factor `2ℓ+1`; continuum Gram `4π/(2ℓ+1)` = the pushforward Gram (NOT `2/(2ℓ+1)`); both coordinate systems of the orbit space accepted |
| `LegendreSpace` + `MomentHead` | `numerics/spaces/` | the flat head; the Protocol every carrier reads its layout through (`isotropic_slot`, `degree_block`, rank, `truncated`) |
| the entry's probe | `Quotient.descending_slots` (manifold.py) | fibre-constancy along `quotient_map`; generic sphere points × `SubgroupOfO3.generic_images` (every finite element; six INCOMMENSURATE angles for SO(2)) |
| `Descent` | `numerics/basis/descent.py` | the two realizations, the isomorphism at the BIT (`[M]` `0.0` on 7/7 sphere rules at L=4, 5 real slots of 25), the discriminator as `frame_basis`; the upstairs face axis-keyed (H-4) |
| G0 | `quotient_onto` (manifold.py) + `_descent_arrow` (frame.py) | ONE predicate — the arrow exists — at `FrameBase.__post_init__` and `GalerkinFrame.__init__`; the table pulled back along it (`FrameBase.descent` / `test_descent`) |
| 0.1b | `directional.py` | `_harmonic_frame_measure` DELETED; `angular_frame` binds `self.measure`; `_harmonic_basis` dispatches on the SUPPORT through `Descent.for_entry(...).frame_basis`; `folded_by` is provenance |
| 0.6 | `spherical_harmonic_basis.py` | `evaluate` refuses off-sphere directions (a real `raise`); `even_slot_mask` reads the entry; `_PARITY_PROBE_DIRECTIONS` retired (`[M]` off the sphere) |
| the carrier's head reads | `material_field.py`, `scattering.py`, `fission.py`, `n2n.py`, `_bases.py`, `harmonic_moment_flux.py` | `moment_source`/`moment_emission` take `head=` and pick the contraction by the head's RANK (⭐ the material field's einsum `mfc...,fg->mgc...` spelled the m-axis — a §6b member the census had not listed); `FissionMomentOperator`, `scalar_flux`, `isotropic_part`, `anisotropic_part`, `l_block`, `ng`, `zeros_for_mesh_and_L` read the head |
| the DSA μ-read | `dsa.py` | `axis_cosines(0)` |

**Rulings applied** (§V.5l(c) 1–7): flat `{P_ℓ}` (after 2.5); the probe on
the entry; G0 = the descent arrow with `P_ℓ` on a full sphere BUILT (`[M]`
Legendre on `lebedev(17)`: DIAGONAL Gram, isotropic `ℓ ≥ 1` ≤ `7.4e-16`; on
`level_symmetric(4)` DENSE at L=4, isotropic-clean; on `product(4,4)` the
aliasing row); the Γ-slot deferred to 2.2b; the O(2) refusal deferred to
1.9 (#432) — `[M]` Legendre-on-the-fold is refused by the lattice, with the
docstring naming the issue; the dead-slot frames keep the DenseMetric arm.

**The red loop** (`scratch/_fused_B_red1.log`, scope numerics + transport +
sn/{operators,primitives,acceleration,mesh,angular} + the ERR-080 gate):
**105 failed / 5187 passed**, classified with `[M]` frames — ~57 hand-built
RECTANGULAR moment arrays on 1-D fixtures; 12 material-field verb calls
without `head=`; 7 G0 refusals of a SIZED-vs-UNSIZED same-named manifold
pairing (`EnergyGroups(n)` vs `ENERGY`, `IndexSet(label, n)` vs
`IndexSet(label)`) — fixtures, `[M]` production names one object on both
halves (homogeneous + diffusion 163/163); 4 `_equispaced_frame`; 12 in 2.5's
own gate (its metric-identity rows compared 1-D rules to the SH space); the
D-slab/D3/D4/A2 DENSE-witness re-poses; `test_e1`; `test_q8_4`; the DSA slot
pin; the ERR-080 XPASS ×3; two 0.6 refusals of hand-typed non-unit
directions. **Zero production defects.** Handed to the test-architect
(`scratch/_fused_B_test_architect_brief.md`) as the gate-writing pass.

**Gates** (test-architect, `scratch/_fused_B_test_architect_brief.md` → 2 new
files + 23 edited, **579** test defs in the touched set; pyright ratchet delta
0): `tests/numerics/test_legendre_basis.py` (32 rows — B1 with the GL
dead-slot THEOREM row, B2 both coordinate systems, B3 the matched spelling
with the pure-`lpmv` re-spelling as an in-test mutation that reds, B3b the
independent three-term recurrence, B4 the dual factor, `P_ℓ(±1)`, the
`4π/(2ℓ+1)` Gram equal to the harmonics' and NOT `2/(2ℓ+1)`);
`test_descent.py` (20 — the isomorphism at the BIT on 7 sphere rules, the
axis-keyed refusal, the discriminator, `for_entry` off the sphere);
`test_manifold.py` +38 (B5 the entry's probe with the right-angle NEGATIVE
control at `L ≥ 4` — the memo's H-3: blind below 4 — B5b the padding
denominator, the mirror leg `array_equal` to the old mask, `quotient_onto`'s
arms); `test_frame.py` (the D-slab witness INVERTED — the slab is DIAGONAL
now — D2/D3/D4 re-posed DRAW-FREE on the generalized-eigenvalue range,
`_DENSE_FRAME_CASES` 4 → 6 mechanisms, `TestG0` 5 rows incl. the pullback
identity); `test_spherical_harmonic_basis.py` (B6 the mask through the entry,
B11 the 0.6 refusal with a SUBPROCESS row proving the `raise` is `-O`-live);
`test_quadrature_directional.py` (`test_q8_4` → the 1-D ROUTE gate; ⭐
`test_q8_6` RETIRED by me — it pinned `spherical_harmonics(L) is
angular_frame(L).table`, and the accessor is gone, so its claim is
structural); `test_dense_metric.py` (A2 re-keyed); `test_basis_domain.py`
(`test_e1` inverted to the ADMIT verdict, the Part I refusal kept as G0's
row-5 witness); `test_dsa_rate.py` (B12 + the REFUSAL leg, the only
discriminating one — H-8); `test_harmonic_moment_flux.py` ×2 (B13 the
rank-1-head field through every head read); `test_material_field.py`
(`TestFlatAngularHead`); `test_binding_tightness.py` (`_equispaced_frame`
re-declared on `S^2/SO2_x` + a `COSINE_INTERVAL` twin as G0's refusal
witness, the FIRST message present and the other two absent); the ERR-080
gate's three strict-xfail markers DELETED, `catches("ERR-080")` kept and
**mutation-verified**: the catalogued defect re-forged in-process (G0 off,
0.6 off, the `(μ,0,0)` padding with the full harmonics) reds **16 of 170**
rows with **0 of the 3-D rows** moving. ⛔ **Two of the memo's own claims
refuted at execution**: the A2 witness `folded_product(2,4).angular_frame(2)`
has naive-pinv asymmetry `9.0e-17`, BELOW A2's own `1e-15` floor (an inert
gate) — re-keyed to `equispaced(8)-L4` (`1.6e-30` honest vs `7.1e-15` naive)
+ `LS4-L3`; and `folded_product(2,4)` at L = 2 is DENSE only by rank
deficiency (live-block offdiag `8.6e-17`) — the genuine dense witness is
L = 3 (`1.000`), used for D1–D4. Also: `verifies("manifold-fibre-constancy")`
attached — `[M]` that label was an ORPHAN equation carrying no `documented`
sentinel.

⭐⭐ **G0 CAUGHT A PRODUCTION DEFECT, outside the red loop's scope.**
`OverlapBasis` (the fractional-overlap trial of energy condensation)
inherited `domain = partition_of` — the COARSE partition it SPANS — while
its `evaluate` validates the FINE row count it EATS: a basis whose domain
named the wrong end of its own map, and nothing compared the two until the
frame's G0 did. `[M]` every non-degenerate `Mixture.condense` (public;
reached from `Solution.condense`) raised on the first G0 — **45 tests / 3
files** (`tests/data/test_mixture_condense.py` 30, `tests/sn/test_condensation.py`
11, `tests/transport/test_kernels.py` 4), none in the red loop's scope,
which is why the loop reported zero production defects (a scope is a
denominator — §2). Fixed at the type: `OverlapBasis.fine: Manifold`
(kw-only) is what the functions EAT, `domain` returns it, `from_indicator`
takes `fine=`, `EnergyGrid.overlap_to` passes `EnergyGroups(self.n_groups)`.
Left RED by the test-architect rather than xfailed — the right call. Two
non-blocking findings of its, both fixed: `quotient_onto(S²/{e}, S²)` returned
`None` for an isomorphism (the trivial arm added, unreachable by any shipped
rule); `_moment_blocks`' head guard was one-sided (a rectangular tensor with
a flat head passed it and died later) — the group axis is guarded too.

**Retired with the commit** (coding-standards, 3-search audit): the
pass-through accessor `Quadrature.spherical_harmonics(L)` — its NAME became
false the day a 1-D rule bound its Legendre basis (the archivist's finding
#4), `[M]` 0 production consumers, 10 test reads rewired to
`angular_frame(L).table`, 1 docs mention re-pointed, `test_q8_6` (the `is`
pin between the two spellings) retired with it; `_PARITY_PROBE_DIRECTIONS`
(`[M]` its five directions were OFF the sphere, norms 0.83–0.998);
`Quadrature._harmonic_frame_measure` (the forgery arm).

**Docs (archivist, same session).** 12 pages + `matrix.rst` (+1756/−240):
ERR-080 **✅ FIXED 2026-09-02** with its before/after table against a pinned
pre-repair tree, the three aliasing rows and the dead-slot theorem;
`manifolds.rst` a new chapter *What descends* (`manifold-descending-slots`,
`manifold-descent`, `manifold-g0-descent-arrow`) + status YAML + Key Facts ×4
+ seams row + changelog; `frame.rst` `frame-g0-descent-arrow` (the seven-row
table, the pull-back, the MomentHead) and the DENSE roster re-measured;
`spaces.rst` `spaces-moment-head`; `spherical_harmonics.rst` the chapter
`sh-legendre-is-the-1d-family`; `angular_quadrature`, `adjoint`,
`operator_algebra`, `acceleration`, `indexing_and_layout`, `slab_multigroup`,
`cartesian_multid` re-tensed — ⭐ the moment layout `(L+1, 2L+1, ng, …)`
had been asserted as THE contract at **9 sites / 7 pages**. Three new
`documented` labels (sentinels 579 → 582). `sphinx -E -W` 0/0/0;
`dead_references` 0 / 52; the docstring-xref checker 0. Four production
docstring defects reported and applied after it returned (the two relayed
precision numbers in `legendre_basis.py` re-scoped — `[M]` pure `lpmv` moves
the TABLE by `4.4e-16` and the converged FLUX at L = 1 by `2.8e-14`;
`descent.py`'s dead 1-D-mask arm removed; the accessor retirement above).
⛔ Two of MY brief's claims refuted by its re-derivation: the isotropic
census reads ≤ `4.3e-16` on 49 of 52 rows (not "≤ 8e-16 on every"); the
135-array 3-D bit-identity was published as a THEOREM (same basis class,
same measure object, the identity arrow) rather than as a sample.

**Mutation battery** (`scratch/_fused_B_mut/mutplugin_B.py` + `_fused_B_battery.sh`,
in-process at `pytest_configure`, a CALL COUNTER per arm; scope the 16 gate
files, 585 rows; `scratch/_fused_B_battery2.log`, baseline 0 reds):

| arm | mutant | reds | verdict |
|---|---|---|---|
| POSITIVE CONTROL | `legendre_table` → ones | **43**, CALLS 163 | the harness sees |
| m0b `μ → −μ` in the Legendre table | (the memo expected NULL — GL nodes are symmetric) | **27** | LOADED: `P_1 = −μ` flips every odd column, and the descent's bit-identity sees it |
| m1 pure-`lpmv` spelling at ℓ = 1 | | **25** | caught (B3 and every bit-identity row) |
| m2 `eval_legendre` at ℓ ≥ 2 | | **22** | caught |
| m3 dual factor `2ℓ+1 → 1` | | **13** | caught |
| m4 `gram_structure` → DENSE | | **2** | caught |
| m5 probe: incommensurate → right angles | | **12** | caught — the `L ≥ 4` control (H-3) |
| m6 probe directions off the sphere | | **37** | caught (0.6 refuses the probe) |
| m7 `even_slot_mask` unrestricted | | **5** | caught |
| m8 upstairs refusal keyed on measured alignment | | **1** | caught — exactly the axis-keyed refusal row (H-4) |
| m9 G0 always admits | | **6** | caught |
| m10 G0 compares the REALIZATION | | **1** | caught — ⛔ BLIND in the first battery: the dispatch refuses a chart-level rule before G0 sees it, so no test reached the predicate; a DIRECT `GalerkinFrame(LegendreBasis, rule-on-[-1,1])` witness was added |
| m11 the forgery restored | | **81 + 4 errors** | caught (B-K, the route gate, 0.6, G0 rows) |
| m12 dispatch on `folded_by` | | **121 + 4 errors** | caught |
| m13 0.6 as a bare `assert` | | **3** | caught — the SUBPROCESS `-O` row is one of them |
| m14 DSA reads `mean_axis_cosine(0)` | DECLARED weak | **0** | blind, as the memo predicted (H-8): the value legs are `array_equal` on 5/5 GL rules and the refusal leg pins the ACCESSORS, not DSA's read; declared in B12's docstring rather than counted |
| m15 `moment_space_on` mints from `L` | | **11** | caught |
| m16 `scalar_flux` reads `[0, 0]` | | **6** | caught |

One declared-blind arm (m14), every other arm loaded, every mutant ran.
⚠ The first battery run (`_fused_B_battery.log`) ran on the tree BEFORE the
last seven fixture fixes (its `none` arm read 7) and exposed m10's blindness;
the table above is the re-run on the corrected tree.

**The 13-tree gate — exit predicate 5, MEASURED at `5436184e`** (`[M]` serial
`-O -m "not slow" -p no:randomly`, per-tree denominators PREDICTED first by
collect-only; `scratch/_fused_B_full_gate.sh` → `_fused_B_full_gate.log`,
per-tree `_fused_B_gate_<tree>.log`; launched on the working tree that became
`5436184e`, finished 13:00 — the commit was made mid-run, per the user's
pre-emptive-compaction ruling): **13 of 13 rc=0; 10471 passed · 19 skipped ·
66 xfailed · 0 failed = 10556 collected = the prediction, unit-for-unit on
every tree.**

| tree | predicted (collected) | measured | rc | wall |
|---|---|---|---|---|
| numerics | 2809 | **2809** passed | 0 | 30 s |
| transport | 707 | 706 + 1 sk | 0 | 23 s |
| geometry | 732 | 727 + 4 sk + 1 xf | 0 | 9 s |
| data | 237 (+2 desel) | 237 | 0 | 42 s |
| homogeneous | 50 | 50 | 0 | 2 s |
| diffusion | 113 | 113 | 0 | 2 s |
| cp | 141 (+13 desel) | 141 | 0 | 66 s |
| moc | 121 (+3) | 121 | 0 | 105 s |
| mc | 41 (+16) | 39 + 2 xf | 0 | 109 s |
| cross_method | 81 (+8) | 81 | 0 | 13 s |
| sn | 3439 (+116) | **3391** + 1 sk + **47** xf | 0 | 15 min 43 s |
| derivations | 1661 (+67) | 1637 + 13 sk + 11 xf | 0 | 37 min 16 s |
| root + harness | 424 (+2) | 419 + 5 xf | 0 | 20 s |
| **total** | **10556** | **10471 passed · 19 skipped · 66 xfailed · 0 failed = 10556** | 0 | 61 min |

Against the 2.5 row (`10454`): **+102 collected = +105 passed − 3 xfailed**,
reconciled FILE BY FILE by collect-only against the `2ebdbc05` worktree
(`scratch/_fused_B_collect_reconcile.log`, same runner and filter):
numerics **+70** = `test_legendre_basis` +32 (new) · `test_descent` +20 (new)
· `test_manifold` +11 · `test_frame` +10 · `test_spherical_harmonic_basis` +7
· `test_quadrature_directional` −10 (`test_q8_6` retired); transport **+24** =
`fields/test_harmonic_moment_flux` +13 · `test_material_field` +8 ·
`test_binding_tightness` +3; sn **+4** = `primitives/test_harmonic_moment_flux`
+3 · `acceleration/test_dsa_rate` +1 — and inside sn the ERR-080 gate's three
strict-xfail rows are now three passes (xf 50 → 47, passed +7 for +4
collected); root+harness **+4** = `test_layer_imports` 349 → 353, the layer
gate parametrizing over every production module and the commit adding four
(`legendre_basis`, `descent`, `legendre_space`, `moment_head`). ⚠ The
pre-compaction checklist had guessed two of these terms — `test_manifold`
"+38" and the harness +4 as "the matrix's rows" — from diff impressions, not
collect counts; the list above is the measurement. Wall 61 min on the machine
that was also running the session.

### 1.9 EXECUTED — `SubgroupOfO3.O2(axis)`: an orbit space is named by its STABILISER (2026-09-02, landed `c1d53206`)

**Goal (the outcome).** A basis on the slab's orbit space declares the FULL group its
functions have, so the frame's lattice gate G0 admits every mathematically admissible
pairing — in particular the Legendre basis on a σ_b-folded rule — and refuses nothing
that the mathematics admits.

#### §V.5m — the opener (the 20th consecutive opener to correct its own section)

Premises re-measured at `e590623f` before designing (`scratch/_19_census.py`, re-runnable):

* §1 existence, with positive controls (`class SO2` 1, `"SO2_x"` 3): `class O2` **0**,
  `def O2(` **0**, `"O2_"` keys **0**, word-bounded `O2(` **0** — landmine #14 held; the
  only `O2` in the tree was `Dinfh`'s rejected old name (`symmetry.py:136-148`).
* The prose corpus (§1's rejected-name clause): `.claude/agent-memory/cross-domain-attacker/
  quadrature_symmetry_quotient_frames.md:225` had ALREADY concluded, during the quadrature
  campaign, *"`GaussLegendre1D`'s tag must move `SO2 → O2` in the same commit"*; the
  test-architect's `L69f` records the over-refusal. The name was free because it had been
  rejected for the WRONG group and never minted for the right one.
* ⛔ **The issue's row was under-specified, and the opener corrected it.** #432 offered two
  designs and asked for a ruling: (1) ONE catalogue entry named by the stabiliser
  (`by = O2(a)`), or (2) TWO entries (the slab keeps `S²/SO2_a`, the Legendre basis moves to
  `S²/O2_a`, and G0's third arm learns an "identical invariants ⟹ identity map" case for the
  slab's `(N,)` μ-nodes). `[M]` the invariant rings coincide
  (`R[x]^{SO(2)_a} = R[x]^{O(2)_a} = R[x_a, x_b²+x_c²]`), so (2) is two unequal `Quotient`
  objects for one manifold — the defect class 2.1 and 2.4 removed — plus a special arm.
  **User-ruled (1)**, with two riders: `SPHERE.quotient(SO2(a))` REFUSES with the theorem
  rather than normalising (the user's rule: *"if option 2 is more principled but requires
  more effort, option 2 is the right option"* — a call site naming a non-stabiliser group is
  the state to make unspellable, and one orbit space then has one spelling in the tree), and
  the registry's slab/sphere `continuous_isotropy` flips to `O2("x")` (with O(2)_x spent the
  recorded residual `Mirror("x")` is exactly G/G⁰; under SO(2)_x the true residual is the
  Klein four-group and the recorded mirror is half of it).
* The naming blast radius, `[M]` by regex per tree BEFORE the edit: the token `SO2_[xyz]`
  in 5 production files (10 hits) / 6 test files (27) / 5 theory pages (50) / this plan
  (30); `SubgroupOfO3.SO2(` at 9 production sites (5 orbit-space callers + 3 in
  `symmetry.py` + 1 docstring) and 112 test sites in 11 files (most are LATTICE tests that
  keep SO2; the orbit-space ones flip).
* `[M]` the walk's cost before: slab 6.0 ms · product(4,8) 114 ms · level_symmetric(4)
  348 ms · folded_product(4,8) 110 ms · lebedev(11) 310 ms (warm caches, loaded machine).

#### What landed (one commit)

| piece | where | the sentence |
|---|---|---|
| `O2(axis)` | `symmetry.py` (`class O2`, `SubgroupOfO3.O2`, repr/name `O2_x`, `rotation_axis` answers) | the stabiliser `{g ∈ O(3): g ê_a = ê_a}` = C_∞v; identity component `SO2(a)`; `D_∞h = O(2)_z × {e, σ_h}`; the rejected-name note reconciled in place |
| computed axial relations | `_fixes_axis` (memoised) + `_axial_contains` (replaces `_so2_contains`) | a finite group ⊆ `O2(a)` iff every element fixes `ê_a`; ⊆ `SO2(a)` iff also proper — `[M]` the tabulated arm said `SO2('x') ⊉ C_1` while `SO2('x') ⊇ {e}`, and `test_cn_in_so2[1]` pinned the wrong one |
| exact invariance | `_invariance_on_points` | `O2(a)`-closed ⟺ axis-supported (the SO2 criterion; axis points are fixed by the vertical mirrors) |
| the probe's images | `generic_images` | both components: 6 incommensurate rotations + each ∘ one vertical mirror (12) |
| the walk | `candidate_groups` | offers `O2` on three axes; `[M]` slab → `{O2_x, sigma_x}` (walk = brute), the generated group D_∞h about x unspellable |
| ONE entry, refusal | `manifold._sphere_mod_o2` (was `_sphere_mod_so2`), catalogue keys `O2_*` (+ `SO2_*` routed to the refusal), `Quotient.by` documented | `by = O2(a)`; `SPHERE.quotient(SO2(a))` → `ValueError` *"S^2/SO2_x is the orbit space S^2/O2_x … spell SubgroupOfO3.O2('x') (#432)"* |
| the spenders | `registry.py` slab/sphere `O2("x")`; `rules_1d.gauss_legendre_on_polar_orbit` | `[M]` `gauss_legendre(8).measure.support.name == 'S^2/O2_x'`, `quotient_group == O2('x')`, space `L2[S^2/O2_x]` |
| the basis | `LegendreBasis.domain` → the O2 entry; `LegendreSpace.from_L` names itself off `basis.domain.name` | `invariance_group == O2(axis)` — the FULL group; one spelling of the entry's name |
| G0 | unchanged in code | `[M]` `GalerkinFrame(LegendreBasis(L), folded_product(4,8).measure)` constructs, arrow `S^2/sigma_y -> S^2/O2_x`, table `(16, L+1)`; isotropic ℓ≥1 moments ≤ 1.4e-15 at L = 2, 4, 6, ℓ = 0 = 4π (no aliasing: the φ-rule is exact to trig degree 7 and the fold is σ_y-even); the x-fold negative leg: `LegendreBasis(axis="x")` on `product(4,8).quotient(Mirror("x"))` REFUSED, `axis="z"` admitted |
| docstrings | 14 production modules | every entry-level spelling of the slab's orbit space reads `S^2/O(2)_a`; the math about the rotation half's orbits stays SO(2) where that is what the sentence is about |

`[M]` after: walk cost slab 5.7 ms · product 133 ms · LS4 366 ms · fold 115 ms · lebedev
328 ms. ⛔ My reading of that as *"+3 candidates ≤ 5 %"* was ONE warm-cache draw on a loaded
machine, and the archivist REFUTED it the same afternoon with the right instrument (min of 15
interleaved repeats, `scratch/_p19_probe6.py`): **+11.3 – 26.2 %** — slab 5.0→5.6 ms, fold
103.5→115.2, product 107.5→135.6, lebedev 283.8→329.9, LS4 320.8→367.1. Same cause as the
2026-08-28 draw row of `plan-authoring` §4: a timing is a stochastic measurement and its
configuration includes the repeat protocol. No new rebuilds (memo lesson 6). `[M]` the
probe on the O2 entry at L = 4: 5 real slots of 25 (the m = 0 column) — and under a
right-angle sample it generates C_4v, not C_4, admitting only the σ_v-EVEN member of each
±4 pair (6 at L = 4, 8 at L = 5; the SO(2)-only probe read 7 / 10): the mirrored half of
the probe does real work. pyright over the 14 touched modules: **0**.

#### The §6b set (the red loop is the enumeration)

Red loop 1 (numerics / transport / harness, `_p19_red1.log`): numerics collection KILLED
by `test_registry.py` building the old spelling at import (the harness lesson: production
in a parametrize list); transport 6 reds, all the refusal. Mechanical re-spell of the
tests' CODE patterns (`_p19_respell_tests.py`, 10 files: `SPHERE.quotient(SO2(` ×41,
names ×7, `continuous_isotropy` ×3, `.quotient_group ==` ×3, `.by ==` ×1,
`invariance_group ==` ×2, `has == spent ==` ×1). Red loop 2 (`_p19_red2.log`): transport
706 + 1 sk green, geometry 727 green, harness 403 green, numerics **17 reds / 2792** —
every one a SEMANTIC member: the G0 row that asserted the refusal (now the admitted
witness + the x-fold negative leg), the frame's and `quotient_onto`'s order-relation
walks (SO2 → O2 in the `has`/`groups` sets), `_ROSTER` and the registry rows dict, the
barycentre message pin (*"axial rotation group"* → *"axial group"*), the un-memoised
derivation's name + the refusal added, the right-angle control's numbers (7/10 → 6/8),
the catalogue listing, `test_cn_in_so2[1]` (the corrected edge), the walk's slab answer
(`_p19_fix_tests.py`, 4 files). Red loop 3 (`_p19_red3.log`): numerics 2809 green; the sn
subtrees 2360 passed + 16 xfailed / 0 (8 min 14 s).

#### The elegance round — the review's ONE absence, and what it cost to leave it

The elegance-enforcer's review of the production diff returned four VIOLATIONS and one
CONCERN that explained all four as one absence: **the tree had no *stabiliser* concept on the
group**, so the fact the whole step turns on — *the largest subgroup of O(3) with `SO(2)_a`'s
orbits is `O(2)_a`* — was spelled three times (three `SO2_a` decoy catalogue keys routing to a
refusal; `_sphere_mod_o2` rebuilding the group it was HANDED from an integer through
`"xyz"[axis]`; `Quotient.by`'s docstring declaring an invariant nothing enforced) and lived in
the wrong tier. `[M]` the review's forgery: `replace(SPHERE.quotient(O2('x')), by=SO2('x'))`
constructed, compared unequal to the entry under `==`, and was accepted by `barycentre` and the
polar-axis reader — ERR-080's own defect class reopened by the idiom Pattern 4 blesses, because
the refusal lived in ONE builder. Accepted in full, in the same commit (bias to completion):

| finding | fix |
|---|---|
| no home for the stabiliser | `SubgroupOfO3.orbit_stabiliser` — `SO2(a) → O2(a)`, `SO3 → O3`, every other member itself (`[R]` a finite group's orbit-preserving isometry agrees with one element on an open set, hence everywhere; `D_∞h`/`O(2)_a`/`O(3)` are the stabilisers of their own orbit families; `[M]` 24 of 24 lattice members gated). `Quotient.__post_init__` checks `by == by.orbit_stabiliser`; the door `_catalogued_quotient` checks it BEFORE the lookup — one helper, `_assert_named_by_stabiliser`, two callers, the message fragment *"S^2/SO2_x is the orbit space S^2/O2_x … spell SubgroupOfO3.O2('x')"* kept verbatim. `_sphere_mod_o2` is a pure derivation again with ONE function-scope import; the three decoy keys are gone (`[M]` the `NotImplementedError` listing named 9 entries of which 3 raised; now 6 of 6 obtainable) |
| the `manifold → symmetry` runtime edge the refusal had added, documented as *"the one runtime edge"* while there were two | gone with the import; the module docstring is true again |
| `D_1h = {e, σ_h}` in fresh prose, and its pre-existing twin in `_contains` | `[M]` order **4**, `{e, C_2^x, σ_y, σ_z}` — every element fixes ê_x, which is WHY it sits in `O(2)_x` and in no other axial stabiliser (the archivist measured the same, independently) |
| `_so2_contains` named in a comment after the rename | `_axial_contains` |
| `np.allclose(atol=1e-12)` — the default `rtol=1e-5` made the effective tolerance **1e-5** on a unit component | absolute, `_MEMBERSHIP_ATOL = 1e-9`, shared with `_finite_contains` (ONE band); `_fixes_axis` lost its `proper_only` flag and `_axial_contains` spells `SO(2)_a = O(2)_a ∩ SO(3)` by composing with the lattice's own `SO(3)` relation |
| the registry's geometry table still teaching `SO(2)`, and the chart `[-1,1]` as the domain | `O(2)_x`, `S^2/O(2)_x` (chart `[-1,1]`), *"1-D in x"* |
| `descent.upstairs_columns`' message naming SO(2); `measure.py`'s half-flipped sentence; `barycentre`'s *"axial rotation group"*; the `D_∞h` note's *"the group a cylinder actually carries"* (`[M]` archivist: the cylinder row spends `Trivial`; `Dinfh` is consumed nowhere outside `symmetry.py` and its tests) | the message reads the entry's group; the prose reconciled; the cylinder claim demoted to what is measured |
| `"xyz"[i]` open-coded at 4 sites, plus a literal forward map inside `mirror_axis` | `AXIS_LETTER` beside `AXIS_INDEX`, the index map DERIVED from the letters |
| rename `rotation_axis → fixed_axis` | ⛔ **NOT taken**: `C_n` fixes ẑ pointwise and answers `None`, so the proposed name lies for a family the old one merely omits; the precise term is the **principal axis**, a WIDER accessor with different answers on three families — a design step, **#433** |

pyright 0 on the six modules; red loop `_p19b_red1.log`: **3563 passed / 1 skipped** (56 s).
⚠ Two of the review's own attacks it withdrew are worth keeping: the descent under `O2` is
bit-identical to the pre-diff `SO2` answer at every axis and L ≤ 4 (the invariant-ring theorem
REALIZED, not asserted), and the function-scope `manifold → symmetry` import was safe in 7 of 7
import orders — the defect was honesty, not breakage.

#### Gates (test-architect) — 43 rows + 5 of mine; a 17-arm battery with two net-new witnesses

`[M]` collect delta numerics **+43** (`test_symmetry` +23 — the axial stabiliser section: axis
validated at construction, `rotation_axis` answers, strictly above its rotation half per axis,
three incomparable axes, exactly the mirrors THROUGH the axis, the finite relations computed,
the whole spellable lattice obeys the order-relation laws (23 members / 131 strict edges), the
downward-closure law exercises the new member (1090/166/138/28), invariance is axis support and
needs NO μ-reflection, no shipped cubature is O2-invariant on any axis, `generic_images` carries
both components and every image fixes the axis; `test_manifold` +8 — the entry records the
stabiliser and its two derivation outputs, the rotation half refused with the theorem on three
axes **as a `ValueError`, not the catalogue's `NotImplementedError`** (the placement leg that
gained teeth when the refusal moved to the door), named by the LARGEST group (140/33/107),
strictly bigger than the rotation half (7 vs 3 of 20); `test_frame` +5 — the axis × fold
pairings split exactly on the lattice edge, the admitted fold reproduces the isotropic moments
at L = 2/4/6, the on-axis refusal names both orbit spaces; `test_registry` +5; `test_legendre_basis`
+2 — the space name is READ off the basis's domain). Plus mine after the elegance round (+5):
`TestOrbitStabiliser` (closure + idempotence and the census of which members grow, 24 of 24)
and three rows on `TestTheOrbitSpaceIsNamedByItsStabiliser` (the `replace` forgery refused; the
door's spelling and `S^2/SO3 → S^2/O3` by the same law; every catalogue key OBTAINABLE).
Two existing docstrings were present-tense-false against their own bodies and fixed in place.

**Battery** (`scratch/_p19_mut/`, in-process, numerics + transport, 17 arms, every arm with a
BITE assertion): the P2 control (`O2.name` drops the axis) reddens **624** across 57 files — the
name is a `FunctionSpace` identity component; **two arms had ZERO pre-existing catchers** and now
have witnesses — `O(2)_a ⊄ SO(3)` (5 reds, 0 before) and `LegendreSpace.from_L` hard-coding the
old name (2, 0 before); one arm is declared INERT and gated as a named blindness row: dropping
the mirrored half of `generic_images` cannot move any `descending_slots` mask at incommensurate
angles, because `R[x]^{SO2_a} = R[x]^{O2_a}` — the theorem itself (18 of 18 rows `array_equal`;
the discriminating regime is the degenerate right-angle sample, C_4 vs C_4v, which is the
right-angle control). ⛔ Three refutations of MY brief, all measured: the "σ_v-odd control" is
mathematically impossible (same theorem); `O2(x) ⊇ D_1h` holds on **x alone**; the "9 fold
pairings" has only **two** shipped folds (`product(4,8).quotient(Mirror("z"))` RAISES — σ_z
permutes the polar levels — so the z row is a constructed half-sphere measure, labelled so).
Production observation, NOT a defect: `hash(SO2(a)) == hash(O2(a)) == hash(Mirror(a))` — a
frozen dataclass hashes its field tuple; `__eq__` discriminates the class; sets and both
`functools.cache`s stay correct. It matters for test design only (a `hash(a) != hash(b)` leg
reds on correct code). Post-elegance battery, 19 arms (three re-keyed + `d1` the stabiliser
check inert + `d2` `orbit_stabiliser` the identity): **19 of 19 arms bite** on a 0/0 baseline
(`scratch/_p19b_battery.log`, ~52 s per arm) — `d1` reddens **6** (the `replace` forgery gate
and the placement leg: with the check inert, `SO2` falls through to the catalogue's
`NotImplementedError`), `d2` **7**; the re-keyed `a` (drop the `SO(3)` intersection) **11**,
`k` (axis-blind) **20**, `d` (re-add a decoy key AND bypass the invariant, the two moves the
round made unnecessary) **7**; the controls P1 24, P2 586 + 39 collection errors; `h` 139 + 5
by raising; the thinnest arms `g` and `i` at **2**. In-process, nothing written: the same 14
production files after as before.

#### Docs (archivist) — written to the ACCEPTED design; 7 present-tense-false claims fixed

`+1144 / −355` across 9 theory pages + regenerated `matrix.rst`: a new `manifolds.rst` chapter
(`manifold-orbit-space-stabiliser`; two new eq-labels, `manifold-axial-stabiliser` and
`manifold-axial-invariant-rings`, the invariant-ring identity re-derived in SymPy with two
controls) stating the naming law, `orbit_stabiliser`'s three-row rule, the construction
invariant and the door, six keys / six entries, and the rejected nine-key design kept as a note
with its four costs; the lattice section rewritten for the computed axial relations; `frame.rst`
G0 row 7 ✅ + a negative row, the `TruncatedBasis` census 5/2 → **6/3**; ERR-080 and ERR-072
re-tensed. False claims found that were NOT re-spellings: `D_1h` order 4 (my brief said 2); the
walk cost (above); *"ERR-080 remains open"* ×4 in body prose; `frame.rst`'s "the third member
does not ship"; two future-tense "once tracker 3.4 lands"; `rotation_axis` described as
"exactly the rotations about one axis" (false once `O2` answers it; re-described as *the axis
whose polar interval the orbit space IS*, 6 of 21 groups); the lower-bound remedy sentence.
`[M]` stage 0 identical on **24 of 24** (geometry × rule) against a pinned pre-change tree.
Sphinx `-E -W`: 0 before, 0 after (the set, not the count). `dead_references` on the archivist's
tree: **1 dead / 53 checked** — `orbit_stabiliser`, unlanded at the time; on the final tree
(`scratch/_p19b_sphinx.log`, `-E -W` rc 0, 0 warnings) **0 dead / 52 checked**.

#### Gate — **13 of 13 rc=0, 10604 = predicted** — the fused row's 10556 + 48, all of it numerics (2809 → 2857: the 43 test-architect rows + my 5), every other tree unit-identical to the fused row, the harness 424 unchanged (the archivist's two new equation labels added no parametrized row this time; predicted from collect counts BEFORE the run, phase 1 of `scratch/_p19_full_gate.sh`). Wall **60 min 27 s** on the session's machine (`scratch/_p19_full_gate.log`, per-tree `_p19_gate_<tree>.log`), started at `2ea2da73` with the 52 files dirty.

### 2.2b EXECUTED — the Γ-slot: a geometry admits a FOLD of its domain, and asks the owed symmetry ON the orbit space (2026-09-02, landed `<HASH>`)

**Goal (the outcome).** The shipped cylinder configuration — `folded_product`, a rule on
`S²/σ_y` — is admissible for the cylinder by the selector's own two stages, because the
ontology can say that a rule was quotiented by part of the owed residual, and because the
owed residual is asked of the rule where the rule lives.

#### §V.5n — the opener (the 21st consecutive opener; this one CONFIRMED its section)

Premises re-measured at `4b7d24c3` before designing: `cylinder.admits_domain(folded)` **False**
(stage 0 was equality of supports), `admits_symmetry(folded)` **False** (Γ = D_2h asked on the
ambient representatives: σ_x and σ_z permute them, σ_y maps a `y ≥ 0` node to its absent
mate), `D_2h ⊇ σ_y` True, `select_quadrature` 0 production callers, `orbit_certificate` on the
fold under σ_y `None`. The §6c witness exists alone (the folded cylinder rule, refused at both
stages). The tracker row's four `[M]` claims all reproduce.

**Three rulings (user, all the recommended option):** (R1) the induced action of an isometry
on an orbit space is MACHINERY on `Quotient`, refusing a motion that does not normalise `by`,
and `SubgroupOfO3.is_invariant` asks a quotient-supported measure THROUGH it — one invariance
notion; (R2) stage 0 reads the existing descent arrow `quotient_onto` plus `Γ ⊇ X.by`; (R3)
the axial entries take the same route with the barycentre as the equivariant lift.

⭐ **What the algebra needed that the ruling did not spell.** An isometry `g` descends to
`M/H` iff it normalises `H` — for a FINITE `H` a conjugation check; for a CONTINUOUS group
asking to act, an exact criterion rather than a sampled one (ERR-072's shape): conjugation by
a connected group is continuous, so into a finite group it is constant, i.e. `SO(2)_a`
normalises a finite `H` iff every element commutes with the rotation generator `[ê_a]_×`, and
`SO(3)` iff `H ⊆ {e, −I}`; the discrete coset representatives (`σ_v`; `σ_h, σ_v, C_2'`;
`−I`) are then finite elements. And a continuous group's identity component must FIX every
node of a finite set on the orbit space (its orbits are connected), which for a finite `H` is
the axis-support / origin rule on the SECTION nodes — the criterion the ambient arm applied
since ERR-072, now stated once.

#### What landed (one commit)

| piece | where | the sentence |
|---|---|---|
| the point-set match, free of the motion | `geometry/transformation.py` `permutation_between`, `permutation_preserving`; `RigidMotion.permutes`/`preserves` call them | an induced action is not an ambient isometry, and its images must be matched by the same rule (ERR-073's bijection guard included) |
| the lift | `Quotient.lift` | barycentre (axial; not a section — equivariant under the normaliser, which is all an induced action needs, and why the axial entry needs no section) / hemisphere section (mirror) / identity (trivial); `π ∘ lift = id` |
| the induced action | `Quotient.induced_action(motion)` + `section_coordinates` | `[p] ↦ [g p]`, section coordinates in (chart-width lifted), chart coordinates out; REFUSED unless `by.is_normalised_by(motion)` — `[M]` `C_4` about z on `S²/σ_y` refused; `σ_x` acts as `(x,z) ↦ (−x,z)` on the disk and `μ ↦ −μ` on the interval |
| the normaliser | `SubgroupOfO3.is_normalised_by`, `normalises`; `_rotation_generator`, `_continuous_decomposition`, `_identity_component_normalises` | `[M]` `O2(x) ⊳ σ_x` True, `⊳ σ_y` False; `D_∞h ⊳ O2(z)` True; `SO(3) ⊳ σ_x` False; `O(3) ⊳ O2(x)` False; `D_2h ⊳ σ_y` True; `C_4 ⊳ σ_y` False |
| ONE kernel | `_invariance_on_orbit_space`; `_check_invariance` routes a bare support to `S²/{e}`; `_invariance_on_points` RETIRED; `_polar_axis_of` RETIRED; `_embedded_nodes` reads `section_coordinates` | normalise → `G ⊆ H` trivial → identity component fixes every node → coset representatives / elements permute the CHART set through the induced action (`_orbit_closure(images_of=)`) |
| the certificate | `orbit_certificate` on the orbit space | agrees with `is_invariant` by construction; II.11's refusal of 1-D nodes by SHAPE is gone as a by-product (`[M]` `orbit_certificate(gauss_legendre(8), σ_x)` exists) |
| the spent-group door | `_catalogued_quotient` | `(M/H)/G` for `G ⊆ H` is `M/H`: refused with the theorem, never "derive `S^2/sigma_y/sigma_y`" |
| stage 0 | `AngularSymmetry.admits_domain` | `quotient_onto(S²/G⁰, X)` exists AND `Γ ⊇ X.by` — the frames' G0 and the selector read ONE arrow; the rejection message names the missing arrow |

`[M]` the smoke on the final tree: cylinder ADMITS the fold at both stages; the fold under
`σ_x/σ_y/σ_z/D_2h` True, `C_4`/`O2_x` False; the slab rule's nine answers unchanged; product
(bare sphere) unchanged; `walk(fold) = {D_2h}`, `walk(slab) = {O2_x, σ_x}`, `walk(product) =
{D_8h}`; `orbit_certificate(fold, σ_y)` exists. pyright 0 on the four modules.

#### The §6b set (the red loop is the enumeration)

Red loop 1 (numerics + transport + geometry + the two sn `_embedded_nodes` consumers +
harness, `_22b_red1.log`): **1 red / 4367** — `test_the_fold_consumes_the_symmetry_idempotent_only_on_a_trivial_action`,
whose first leg pinned *"the fold is no longer σ_y-invariant"* (the ambient reading) and whose
second leg pinned the catalogue's `NotImplementedError` naming `S^2/sigma_y/sigma_y` as work to
derive. Both are the same fact now — a quotient by the SPENT group is refused at the door with
the theorem — and the test says so. Red loop 2 (numerics, `_22b_red2.log`): **2857 / 0**.

#### Gates (test-architect) — 14 classes / 148 rows, pasted; the §6b set READ by an in-process call census

The test-architect worked against a pristine `4b7d24c3` snapshot (the tree went unimportable
under it for twenty minutes — the decorator slip) plus a SHADOW of the design outside `orpheus/`
whose validity control is that it reproduces the shipped answers exactly on 395 of 415
(rule × candidate) rows. The §6b set was measured, not read: a `-p` plugin recorded every
`is_invariant` / `admits_domain` / `admits_symmetry` call over six trees (5725 tests) —
**1636 `is_invariant` calls / 74 tests, 1 verdict moves** (the fold-idempotence test); 107 + 59
stage calls, **0 move, not one ever sees a fold**; 1 internal pin (`_embedded_nodes ==
barycentre`, which survives because the lift was re-pointed rather than deleted). Pre-change
literals pinned: the slab and polar-orbit rules' full 15-candidate verdict dicts, the fold's
21-candidate dict before and after, the 16-row selector table, the lattice rows. Gates, one
per named input they REJECT at landing: the lift is a right inverse of the quotient map on all
three families; the induced action is well-defined (`C_4` about z and a rotation about y REFUSED
on `S²/σ_y` and `S²/O2_x`) and acts as the theorem says; `ambient_representatives` dispatches
on width and refuses anything else; the embedding reads the lift; the normaliser criterion
agrees with a brute conjugation CONTROL at incommensurate angles (the control may sample; the
production code never does); the trivial quotient answers exactly like the bare sphere; the
slab and sphere families answer exactly as before the carve (the pinned dicts); the fold is
invariant under exactly the groups that descend and close; the compatibility law holds on the
orbit-space route; the walk on a fold reports `D_2h`; `orbit_certificate` and `is_invariant`
AGREE (re-posed from "now answer different questions" once the certificate moved too);
stage 0 is the descent arrow plus the owed residual (cylinder × fold ADMITTED, slab × fold
refused, slab × slab identity); stage 1 on a fold asks the orbit space (a fold with a dropped
node REFUSED); the selector is UNMOVED (16 of 16 rows pick the same spec — the Γ leg is inert
at that tier because stage 2's V conjunct refuses a fold first). Refutations of MY memo, all
measured and kept: 4 flips per fold not 2 (`σ_y, C_2, D_1h, D_2h`; 20 of 415 rows); `σ_y`'s
flip is the least falsifiable (the `H ⊇ G` theorem never reads a node); "a fold by σ_x of the
σ_y-fold works today" was FALSE; `_embedded_nodes` also feeds `ordinate_permutation`, so R3
read as DELETE would have embedded an `S²/O2_y` rule along x (max |Δ| 0.96) — the re-point
route was taken. `[M]` the 148 rows: 147 green against the first carve (the one red the
`orbit_certificate` docstring, fixed in the elegance round), 592 of 592 across the five
affected files after it, 4882 / 0 on the wide loop.

#### The elegance round — 3 violations + 7 concerns, all taken

The reviewer's own control: the retired ambient kernel re-instated in-process and replayed
over 23 groups × 10 measures — **3 of 230 answers changed**, `σ_y`, `C_2`, `D_2h` on the fold,
all intended; every non-fold family bit-identical. The exact normaliser held against a
brute-force conjugation reference on **529 pairs and 138 motions, 0 disagreements** (its own
first run reported 46 — its reference's fault, `plan-authoring` §4's verify clause working).

| finding | what was wrong | fix |
|---|---|---|
| V1 two camps | `Quadrature.ordinate_permutation` matched the AMBIENT nodes: on the fold `σ_y` → `None` while `is_invariant` → True; the boundary realizer would REFUSE a reflection that acts trivially on the orbit space while the registry admits the rule — three subsystems away | `symmetry.induced_permutation(measure, motion, atol)`, the single-motion face of ONE closure `_orbit_space_closure` read by `is_invariant`, `orbit_certificate` and `ordinate_permutation`; `[M]` the fold's `σ_y` permutation is the identity of the stored ordinates |
| V2 a lift that lies by its name | `section_coordinates` returned the barycentre (inside the ball) under a name promising a point of the base — ERR-080's shape | RENAMED `ambient_representatives`, docstring says what it returns and why; `lift`'s arms typed consistently (domain = the entry; codomain `D^3` / `S^2` / `S^2`) |
| V3 a guard unable to fail | the docstring claimed a continuous `H` was "answered at step 2" (`H ⊇ G`) — false for `D_∞h` on `S²/O2_z` (`SO2_z ⊆ O2_z`, `D_∞h ⊄ O2_z`); the fall-through position test on the barycentre was a tautology (on the axis by construction) | step 2 tests `H ⊇ G⁰` too; `SubgroupOfO3.identity_component` (a GROUP; `_continuous_decomposition` returns it, the `"SO3"` sentinel gone); the position test runs ONLY for a finite `H`, where it is exact |
| C1 | the width validation re-duplicated (1 → 3 copies) by the very extraction meant to single-source | one home: `on_points` + `permutation_between`; one test pin re-pointed |
| C5 | `is_normalised_by` read the whole motion on the finite arm and the linear part on the continuous ones (a translation normalised `D_∞h` and not `O_h`) | the LINEAR part on both — the convention `ordinate_permutation` states (a point group acts on directions) |
| C6 | the spent-group door refused `(S²/σ_y)/{e}` while `S²/{e}` is an entry; and `S²/{e}` was the container for padded points OFF the sphere | `{e}` admitted on every base first, named as the one exception; the ambient container is `RealSpace(3)/{e}` |
| C7 | `admits_domain` read `X.by` relative to the WRONG base and needed an equality branch (`[M]` the slab's `σ_x ⊉ O2_x`) | `spent_group(D, X)` beside `quotient_onto` — what the ARROW spends; one expression |
| C4 | the hemisphere lift's clamp turned an off-chart point into a plausible wrong one | refuses `ρ > 1` past a round-off band; the draft gates witness all three lift arms |
| C2 / C3 / N1 | `_as_columns` at four sites; the `"SO3"` string; finiteness spelled twice in `orbit_certificate` | two sites inside the one helper (the producer's `(n,)` for a 1-D chart is the tree's node convention, left); the group; one spelling |

Withdrawn by the reviewer itself, worth keeping: `by.name == "Trivial"` is the house idiom (5
prior sites); the barycentre is an isometry of the chart onto the axis, so the match window is
unchanged; the four `assert elements is not None` are type-narrowing over a total dispatch.

**Battery** (`scratch/_22b_mut/`, in-process, numerics 3005 rows, 14 arms; baseline 0/0):
the positive control (`is_invariant` always True) reddens **60**; the normaliser made inert
**37**; the axial lift as the zero-pad forgery **28**; the induced action applied to the
AMBIENT points **90**; stage 0 back to equality **3** (the fold rows only — `[M]` equality and
the arrow+Γ predicate agree on 24 of 28 geometry × rule rows); stage 0 without the Γ leg **4**;
stage 0 without the identity case (every rule asked `Γ ⊇ X.by`) **44**; stage 1 asked on the
ambient nodes **18** (re-keyed after the round: the trivial orbit space for every support);
the identity component's fixedness dropped **23**; the chart closure without its weight leg
**4**; `_embedded_nodes` reading column 0 instead of the lift **24**. Two arms are the
findings rather than the reds: the bijectivity leg is **UNINSTALLABLE** — the `Permutation`
TYPE refuses a non-bijection at construction ("indices are not a permutation of range(17): 1
value(s) repeat"), so the ERR-073 guard is held twice and no battery can redden it (Pattern 4,
the strongest witness); and the kernel with BOTH containment short circuits removed reddens
**0 of 3005** with a real mutant installed (the original arm set a flag nothing read; re-keyed
to a copy of the kernel without steps 2) — the short circuits are theorems the closure
re-proves, not branches that decide an answer. Tree integrity after: the same 5 production
files (in-process, nothing written).

#### Docs (archivist) — a new chapter on who ACTS on an orbit space; +1404 over four pages

`manifolds.rst` **+1232**: a new `=` chapter *Who ACTS on an orbit space — the normaliser, the
lift, the induced action* (8 subsections: the descent criterion, eq-label
`manifold-normaliser-descent`; the sampling control; the lift; `induced_action`; the ONE
four-step kernel; the Γ-slot, eq-label `manifold-gamma-slot-stage-zero`;
`ordinate_permutation`; `orbit_certificate` and II.11; the spent-group door; what moved) plus
12 in-place repairs (a ⛔-tombstoned warning, the `_polar_axis_of` import-census row, a
DISCHARGED seam, the barycentre tables, the machine header, a dev-history row);
`discrete_measures.rst` **+112** (stage 0 through the arrow, eq-label
`quadrature-stage-zero-descent`; stage 1; a witness row; the 24-of-24 note);
`error_catalog.rst` **+42** (II.11 now WHOLE: `orbit_certificate(gauss_legendre_on_mu(8), σ_x)`
went `None` → 2 permutations); `frame.rst` +18; `matrix.rst` regenerated (584 → 587
documented labels, exactly the three new ones). Every "before" reading from a pinned `git
archive 4b7d24c3` copy: stage 0/1 — cylinder AND cartesian2d admit the fold at both stages
(both False before); `walk(fold)` `{σ_x, σ_z}` → `{D_2h}`; the compatibility law 0
violations at 342 and 450 pairs both sides; `π ∘ lift = id` exactly on all four entries;
bare-`S²` vs trivial-quotient agreement 150 of 150; the right-angle sample over-certifies the
normaliser on 2 of 8 pairs while ten incommensurate angles agree 8 of 8 (the control);
stage-0 refusals over the 5 factories × 4 geometries **12 → 10 of 20**, no pair `True →
False`. Sphinx `-E -W` 0 → 0 as a set; `dead_references` 0 / 52; the L-062-patched xref gate
0 dead with an end-to-end positive control (2 / 2). Claims of mine it refuted: the memo's
"fold by σ_x of the σ_y-fold works today" (false both sides); the delta's "refuses a translated
motion" (it takes the linear part) and "`Ball(3)`" (it is `RealSpace(3)`); and it found the
certificate's refusal message a THREE-arm disjunction wearing two-arm text (fixed in the same
commit: continuous / does not act on the orbit space / acts without permuting). Its one owed
follow-up: #370's gap 2 is CLOSED the way that issue demanded (the lattice arrow, not a
widened tag); gap 1 stands (`folded_product(4,8).measure.exactness is None`).

#### Gate — **13 of 13 rc=0, 10752 = predicted** — the 1.9 row's 10604 + 148, all of it numerics (2857 → 3005: the test-architect's 14 classes; my own test edits add no row), every other tree unit-identical to the 1.9 row, the harness 424 unchanged (the archivist's three new equation labels added no parametrized row); predicted from collect counts BEFORE the run (phase 1 of `scratch/_22b_full_gate.sh`). Wall **59 min 27 s** on the session's machine with Sphinx building alongside for the first minutes (`scratch/_22b_full_gate.log`, per-tree `_22b_gate_<tree>.log`), started at `4b7d24c3` with the 26 files dirty. The sn tree's 3439 is the witness that the moved `ordinate_permutation` changed nothing for the boundary realizer and the deck consumers on every shipped configuration.

Predicted BEFORE the run from collect counts: numerics **3005** (2857 + the 148 drafted rows;
my own test edits add no row), geometry 732, transport 707, every other tree as the 1.9 row
— **10752** against the 1.9 row's 10604, +148 all numerics, the harness expected at 424 unless
the docs pass adds a parametrized row.

### ✅ 0.2 — EXECUTED 2026-09-01, landed `ce46181c`. What the split actually found

**The shape that shipped.** `axis_cosines(i)` is now the COORDINATE question
and refuses a suppressed axis with a real `raise` (⛔ not an `assert` — the
canonical runner is `-O`). `mean_axis_cosine(i)` is the ORBIT MEAN, delegating
to `axis_cosines` in-domain so there is one source for the coordinate, and
returning the derived zero on a suppressed axis. The three census-confirmed
flow consumers moved to it; `mu_y`/`mu_z` became views over it (the census says
both are read for the flux question); `mu_x` and `eta`/`xi` stay strict —
`[M]` `xi` is never read on a 1-D rule (all 24 `.xi`/`.eta` reads are
cylindrical), and `augmented_mesh.py:1837` provably cannot phantom-read
(a mesh's axis count never exceeds the quadrature's `dim`) **and wants the
POINT anyway**, because the sweep marches individual ordinates, not orbits.

⭐ **THREE meanings shared those zeros, not two.** The plan's triage found the
coordinate/flux split. Execution found a third the old accessor also swallowed:
`i >= 3` is not a suppressed axis, it is **no axis at all**, and the old code
returned zeros for that too — so a genuine indexing bug was indistinguishable
from a legitimate 1-D read. `mean_axis_cosine` therefore bounds at the AMBIENT
dimension (`_DIRECTION_AMBIENT_DIM`) and its two refusals say which is which.
`[M]` gated by `test_q9_4`; the mutation that drops the bound reds it.

⭐ **`spherical_harmonics` single-sourced onto the frame.** It carried the
SECOND copy of the 1-D fabrication, so making `axis_cosines` strict forced the
question. `[M]` it has **0 production consumers** (9 test call sites) and `[M]`
was already bit-identical to `angular_frame(L).table` (36 of 36) — the SAME
math procedurally rearranged, which is `coding-standards`' *retire*
discriminator, not its *keep* one. It is now `return self.angular_frame(L).table`.
⚠ That DEMOTES `test_q8_6` from a two-path bit-identity gate to a tautology,
which the retirement rule says to handle rather than revert: the gate is kept,
re-scoped in its own docstring to pin the **single-sourcing** (an `is`
identity, which reds if a second evaluation path returns), and `[M]` the value
coverage survives externally in `test_q5_1`/`test_q5_2`, whose literals were
authored independently of both spellings. `[M]` mutation `m5` (re-introduce the
twin) reds 10 rows.

⭐ **A §6b member no census could return, found by the red loop.**
`test_q4_5` built `np.column_stack([q.axis_cosines(a) for a in range(3)])` to
apply a `RigidMotion` — i.e. it hand-rolled the canonical R³ **embedding**, a
Pattern-2 duplicate of `_embedded_nodes`, and it only worked because the
accessor silently padded. It is neither a coordinate consumer nor a flux one;
it is an *embedding* consumer, a category the triage did not have. Now calls
`_embedded_nodes(q.measure)`. ⟹ the transferable half: **a strictness change
surfaces duplicates of whatever the laxness was standing in for**, and those
are invisible to a census of the lax call itself.

⛔ **And `_embedded_nodes`' own docstring was made false by this step** — it
said the embedding is *"written down in `Quadrature.axis_cosines` and used by
`spherical_harmonics` internally"*, and 0.2 falsified **both halves** in one
commit. Corrected in place with the date and the surviving site.

⛔⛔ **The §6b set was short by four, and the red loop cost 29 tests across two
trees — the miss is instructive because the census had ALREADY FOUND them.**
I enumerated `orpheus/` consumers by AST, triaged each by meaning, and treated
`tests/` as "the red loop will catch it". It did. But §II.14's own census table
carries the line *"⚠ Note `tests/_harness/references.py:57` (48 reads) — a
harness-level phantom read"*, which I read in this same session and did not
connect to 0.2's blast radius. `[M]` geometry 21 red + sn 8 red, every one
reaching `axis_cosines` through that ONE harness helper.

The four, and three of them are the same shape:

| site | what it was | what it wanted |
|---|---|---|
| `tests/_harness/references.py:57` | `column_stack([axis_cosines(a) for a in range(3)])` | the R³ **embedding** (drives all 29) |
| `tests/sn/sweep/curvilinear/test_coupled_pole_mu_level_invariant.py:160` | the same three-call `column_stack` | the R³ **embedding** |
| `tests/numerics/test_quadrature_directional.py` (`q4_5`) | the same again | the R³ **embedding** |
| `tests/geometry/test_bc_universal_invariants.py:894` | `axis_cosines(AXIS_NAMES.index(axis))` | the **flux** verb (mirrors `_specular.py:162`) |

⟹ **the triage had two categories and the tree has three.** Coordinate, flux —
and **embedding**, which is neither: it wants the orbit's image as a POINT of
:math:`\mathbb{R}^3` in order to apply a `RigidMotion`. That category was
invisible to a census of the accessor, because all three sites spell it as a
`column_stack` of three accessor calls rather than as any symbol. It is
`plan-authoring` §6b's *"members spelled without the symbol"* with a new
member: **a hand-rolled composition of the changed call**.
⭐ And the repair is principled rather than mechanical: a mirror is LINEAR, so
it commutes with the orbit mean — matching embedded barycentres is equivalent
to matching orbits, which is why the 1-D arm stays correct through
`_embedded_nodes`.

⚠ **A fifth member, and the one a grep genuinely cannot return: a duck-typed
SURROGATE.** `_MutantQuadrature` (`test_bc_universal_invariants.py:158`)
stands in for a `Quadrature` while returning a deliberately-broken
`ordinate_permutation`. It delegated `axis_cosines` and nothing else, so
re-pointing the reference at `_embedded_nodes(quad.measure)` broke it on a
member it had never been asked for. It now passes through `measure` and
`mean_axis_cosine` too. Exactly the 2026-08-24 surprise-log row, third
instance.

⟹ **the transferable half, and it is cheap:** when a step makes a lax accessor
STRICT, the §6b set is not the accessor's call sites — it is the call sites
**plus every hand-rolled composition of them**, and those are found by asking
*"what was the laxness standing in for?"* (here: an embedding) and grepping for
THAT. ⛔ **And I first wrote "one grep would have found all four" — `[M]` it
finds THREE.** Re-run against `HEAD`, `column_stack` within four lines of
`axis_cosines` returns 6 hits covering the production site and the three
EMBEDDING duplicates, and misses `bc_universal_invariants:894` entirely,
because that one spells it `axis_cosines(AXIS_NAMES.index(axis))`. The fourth
needs a different question — *"which test mirrors the production line I just
re-pointed?"* — since it is the independent reference for `_specular.py:162`.
⟹ **two greps, not one**: one for the composition the laxness stood in for,
one for the twins of every production site the step touches. Cost as it
actually played out: one red loop (29 tests, both trees, all loud).

⚠ **A scanner failure, logged because it is `vv-principles` #17 verbatim.**
The mutation table's first run reported `m3_mean_diverges` as
**0 passed / 0 failed / 0 calls** — the "zero on both scanners" signature. The
arm was fine: re-run, it reds **11** rows (`test_q9_3` ×8 plus three siblings).
The broken instrument was my summary parser (`grep -oE '[0-9]+ passed' | head -1`
on a captured string); anchoring on the summary LINE fixed it. ⟹ the discipline
that caught it is the one #17 prescribes — a per-arm CALL COUNTER, so
"no reds" and "never ran" cannot be confused.

### ✅ 0.2-R — RULED by the user, 2026-09-01

**(a) The flow verb is `mean_axis_cosine`.** Chosen over `axis_flow`
(operational — names the USE, and is only honest on the 1-D side, since for a
3-D rule it returns a plain cosine and not a flow) and over `axis_projection`
(the plan's own placeholder — accurate for both, and silent on the mechanism,
which is the ambiguity that let ERR-080 hide). The ruled name is the only one
that states WHY the 1-D answer is zero: it is
:math:`\langle \Omega_i \rangle` over the orbit, so **the zero is DERIVED,
not defaulted**, and the SAME definition covers a 3-D rule, where the orbit is
a single point and the mean is the cosine itself. One definition, no branch at
the concept level.

**(b) 0.2 lands NOW, decoupled from 3.4** — the 1-D arm of
`_harmonic_frame_measure` spells its own zeros inline rather than obtaining
them from an accessor named *"direction cosine along axis i"*. Rejected:
folding 0.2 into the FUSED landing (ORDER step #5), which would make an already
highest-risk commit (six items, its only witness PRODUCTION) larger.

⟹ **What this does NOT change:** 0.2 still runs now, and the census's triage
still stands. What changes is its FRAMING — it is not *"a pure naming/typing
split with no defect behind it"*, because the fabrication consumer is still
live and stays live until 3.4.

⭐⭐ **And that turns out to be the design, not an obstacle.** If the direction
accessor is to RAISE for `i >= dim`, the 1-D arm of `_harmonic_frame_measure`
cannot keep calling it — so the arm must **spell its own zeros inline**:

```python
nodes=np.column_stack([self.measure.nodes, np.zeros(N), np.zeros(N)])
```

which is strictly better *on its own terms*, independently of 0.2: the zeros
ARE the fiction, and laundering them through an accessor named
*"direction cosine along axis i"* is precisely what made the fiction invisible
for as long as it was. Spelling them locally makes the invention visible at the
one site that commits it, and **decouples 0.2 from 3.4** — 0.2 no longer has to
wait for the fiction to be repaired, only for it to stop hiding.

⭐ **Collateral, found while enumerating the accessor's callers:
`Quadrature.spherical_harmonics` has ZERO production consumers.** `[M]` by AST
**0** call sites in `orpheus/`, **9** in `tests/` (Nexus `callers` returns an
empty list *with* an `unresolved: 9` warning — the AST census is what settles
it, and the warning is why the empty list was not banked). It is a Pattern-2
twin of `frame.table` — `[M]` 36 of 36 bit-identical — differing only in the
input assembly, i.e. the SAME math procedurally rearranged, which is the
`coding-standards` fuller-view-oracle rule's own **retire** discriminator
rather than its keep one. It also carries the second copy of the 1-D
fabrication (`:575-577`). ⟹ **retirement candidate, scheduled against 3.4**
(its docstring already names 3.4 as its repair); not taken here, because
retiring a method with 9 test consumers is its own carve and 0.2 is not it.

⚠ Note `tests/_harness/references.py:57` (48 reads) — a harness-level phantom
read, part of §4.9's migration set, not covered by §4.9's current inventory.

⚠ And a hazard for 0.2's eventual split: `angular_trace_space.py:220` reaches
the accessor as `getattr(quadrature, "mu_y", np.zeros_like(mu_x))` — the
`coding-standards` string-form idiom that **fails in the DEFAULT's direction**.
If `mu_y`/`mu_z` are ever retired, that call does not raise; it silently returns
the default and every branch keyed on it flips.



## II.15 ⛔⛔ FIVE of my execution-time claims, REFUTED by the verification plan (2026-08-31)

`test-architect`, 106 tool calls, monkeypatch-only batteries. Kept in place per
`plan-authoring` §3 — the refutations are worth more than the claims were.

### ⛔ R1 — *"all 6 `Basis` subclasses can DERIVE their domain from their own fields"* is FALSE for 3 of 6

`[M]` over every shipped `(basis, measure)` pairing production constructs, **4 of
10 agree** under my proposed derivations; one mismatch is the intended ERR-080
one; **five are not**, all in the indicator family:

| mint site | my proposed value | the REAL `measure.support` |
|---|---|---|
| `energy_grid.py:220` `as_basis()` | `spatial_R1` | **`energy`** |
| `energy_grid.py:270` `overlap_to()` | `spatial_R1` | **`energy`** |
| `mixture.py:424` | `spatial_R1` | **`energy`** |
| `frame.py:754` `_collapse_pair('angular')` | `spatial_R1` | **`index(angular)`** |
| `frame.py:754` `_collapse_pair('energy')` | `spatial_R1` | **`index(energy)`** |

⭐ The reason is real and not a slip: `edges_per_axis` is an **index partition
with no manifold identity**, and the tree legitimately mints `IndicatorBasis`
against **three** manifolds. `f"spatial_R{ndim}"` hard-codes one of three — the
exact defect this campaign exists to remove, committed by me while removing it.
⟹ **`IndicatorBasis` needs a CONSTRUCTOR-SUPPLIED tag** (5 mint sites), which is
the honest form: the indicator genuinely does not know its manifold until told,
and forcing the mint site to say so IS "no hard-coding at the meeting place".
⚠ This also kills the candidate LAW `basis.<tag> == measure.support` as a
universal — it would false-red on five correct pairings. **It survives only as
the ANGULAR-path predicate, which is all `angular_frame` needs.**

### ⛔ R5 — the property must be called `support`, NOT `domain`

`[M]` `domain` is already a loaded protocol name: **58** definers in `orpheus/`
and **13** `getattr(x, "domain", None)` readers (`operator.py:1168,:1485`,
`flat_operator.py:114`, `coupled_system.py:767`, …) — **all expecting a
`FunctionSpace`**. Part A would hand them a `str`, and the `getattr`-with-default
idiom **fails silently in the default's direction** (`coding-standards`).
`[M]` `support` has 15 definers and `DiscreteMeasure.support` is *the same
concept with the same type*. ⟹ **RULED: `Basis.support`.** Matching spellings
across the two objects is also what makes the composability check read as a
comparison of like with like.

### ⛔ R2/R3 — 0.1c is NOT "the behaviour-bearing half", and my §8 hazard is empirically vacuous

`[M]` a positive control **strictly stronger** than 0.1c (corrupt `support` on
*every* rule) reddens **exactly one** test over 2634 measured rows — and that
one is on `lebedev(17)`, an **unfolded** rule. Simulating 0.1c itself: **0 reds**
across all three scopes (2052 / 351 / 231).
`[M]` and the mechanism I named does not exist: over `OperatorProduct.__init__`
+ `_agreed_space`, **0 of 145 089** two-sided composability checks put an
`L2[S^2…]` space on either compared side. `frame.conjugate` checks
`basis_space`, never `measure_space`; the SN faces are
`HarmonicFrame.flux_analysis_on(angular_space)`, carrying the **caller's** space.

⚠⚠ **The lesson for me, and it is a §8 sharpening**: I traced a real static
chain — `measure_space` → `_FrameAnalysis.domain` → `OperatorSum` refuses
unequal domains — and published it as a measured hazard. Every link is real
code. **The traffic on it is zero.** §8's own prescribed check (*grep the field
against `is None` / equality and read what each hit decides*) is a **wiring**
check by construction, and cannot distinguish load-bearing from inert. And the
error direction is the flattering one: an over-warning reads as rigor, so
nothing prompts a re-check.

### ⛔ R4 — my proposed bit-identity gate CANNOT be the keystone

`[M]` inverting the new predicate outright (*"always take the bridge"*) is
**bit-identical on end-to-end `keff`** across slab/sphere/cylinder and reds
**0 of 120** (Tier 0) and **0 of 1913** (Tier 1). So a value gate is blind to
the routing decision the change IS. `plan-authoring` §10, third shape: the gate
I would have written cannot detect its own campaign's success OR failure.
⟹ **the keystone is a ROUTE gate, not a value gate**: `frame.measure is
q.measure` — `[M]` **red today on all 8 rows** (`False`), green after on the 6
routed rows. Identity, not equality: equality is what the bridge already
satisfies.

✅ **REMEDIED 2026-09-01 by 0.1a** (`plan-authoring` §3's third case — this row
did not become *wrong*, it became *past tense*, which is the class nothing
prompts you to edit). "Red today" is no longer today. Measured over the **12
shipped rules** rather than R4's own 8-row fixture set: the gate read `False`
on **12 of 12** before and reads `True` on the **10** that route.

⭐⭐ And R4's ruling was **re-derived independently on the gates it produced**,
which is stronger than inheriting it: under the exact pre-carve behaviour
(`m1_always_rebuild`), the value gate `test_q8_5` passes **8 of 8** while the
route gate `test_q8_1` reds **8 of 8**. The blindness R4 predicted is therefore
a measurement on this campaign's own committed gates, not a claim carried
forward from a probe nobody can re-run.

### ⛔⛔ R4b — the ENTIRE fold-basis machinery has zero solve-path witnesses

`[M]` **179 folded `angular_frame` calls across 72 tests**, reported as *"100 %
at `L = 0`"* — ⛔ **that half is REFUTED**: my own spy reads **31 calls at
`L = 1`** on a subset (Part XIV). What survives, and is the claim that matters,
is that `[M]` at `L = 0` the `MirrorEven` table is **bit-identical** to the plain
SH table (`max|Δ| = 0.000e+00`; it diverges at `8.688e-01` from `L = 1`). So every
one of those 179 calls witnesses nothing about the fold. `[M]` **deleting
`MirrorEvenSphericalHarmonicBasis` outright reds 0 of 1913.**
⟹ **PRE-CARVE OBLIGATION**: one parametrize row — 2-group heterogeneous
cylinder, `folded_product(4,8)`, `scattering_order=1` — `[M]` takes the fold
basis from **0** witnesses to an **18.5 %** catcher
(`0.7475027…` → `0.6095089…`). It must land **before** any carve touches the
fold path, or the carve is unfalsifiable there. `vv-principles` #17 / §6c.

### ⛔ And my §3 "two independent producers" check is a TAUTOLOGY

`Quadrature.quotient` sets `folded_by=group`; `_harmonic_basis` reads
`folded_by.mirror_axis`; `DiscreteMeasure.quotient` reads `group.name`. **Both
sides derive from ONE `SubgroupOfO3` instance**, so "the two agree" is `X == X`
— worse than a shared-token check, because it looks like corroboration. The
falsifiable content is the ROUND TRIP: `AXIS_NAMES[Mirror(a).mirror_axis] == a`
(`AXIS_NAMES` is at `numerics/face_layout.py:100`). Swap two entries and it reds.

### ✅ What SURVIVED

`[M]` Part A breaks nothing in-tree: 7 transitive subclasses (6 production +
`_DenseTrial(IndicatorBasis)`, which inherits), **0** duck-typed classes, **0**
`SimpleNamespace`/`Mock` surrogates, 0 aliased-base subclasses, 9 Basis/Frame
monkeypatches (all patching methods on real subclasses) — confirmed by an exact
`ABCMeta` simulation (720 + 231 passed, `still-abstract: []`). So **abstract is
safe**, and the §6b surrogate classes that bit twice this month are absent here.
⚠ The census is structurally blind to out-of-tree subclasses (§6.4 of the memo).


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

## III.4b Descent of the four kernel kinds — only `L` is nontrivial

⛔ **RESTORED 2026-08-31.** `[M]` this table is in the source plan
(`symmetry_quotient_plan.md` §I.4) and was **dropped in the merge** — found by a
whitespace-insensitive probe over both files (3 of 3 distinctive phrases absent
here, present there). It is the only such drop; 11 other probes survived. ⚠ The
lesson is the instrument: a **line-based grep cannot see a phrase that wraps**,
and my first pass reported a *second* false drop for exactly that reason
(*"the same\ngate on a different operator"*, which is present at §III.4).
Flatten whitespace before any merge-fidelity claim.

| term | kernel type | descends **iff** | note |
|---|---|---|---|
| `S` | integral, `σ_s(Ω̂·Ω̂′)` | **always** — bi-invariant under all of `SO(3)` | descent *is* Funk–Hecke; `Λ` diagonal in `P_ℓ` |
| `F` | rank-1, isotropic emission | the material field is `H`-invariant | `\|χ/4π⟩⟨νΣ_f\|`; trivially bi-invariant in angle |
| `C` | multiplier `Σ_t(r)` | `Σ_t ∘ g = Σ_t ∀ g ∈ H` | pure material-symmetry condition |
| `L` | differential | mesh and geometry `H`-equivariant | ⭐ **acquires a CONNECTION TERM** |

Reducing a differential operator along a **non-free** action produces a term
supported on the quotient's stratification — which is why `L` is the only hard
row, and why Phase 1.3 (verify-or-kill the connection identification) exists.

⭐⭐ **And the row this campaign was born in is `S`.** The table says scattering
descends *always*, and that its descent **is** Funk–Hecke with `Λ` diagonal in
`P_ℓ`. That is exactly the machinery Part I finds broken — so ERR-080 is not a
counterexample to the row, it is a **discretization** failure of a statement that
is true in the continuum: the theorem descends, and the *discrete* realization
does not, because the retained harmonics are not `H`-stable on the slab node set.
`[R]` this is the sharpest available statement of what G2 checks and why G0/G1
cannot catch it — the continuum row is unconditional, so any gate reading the
continuum will pass.

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


## III.10 ⭐⭐ Three levels, and the code names two — why `FunctionSpace` is not the domain

The user's question, 2026-08-31: *"What was wrong with `FunctionSpace` as a
domain? Is that not the right domain?"* Answering it is what earned D0.7.

A basis function is a map. What it EATS is a point:

```
Y_ℓ^m : S² → ℝ          the argument is a POINT of S²,
                        not a vector of degrees of freedom
```

`[M]` `FunctionSpace` is documented as *"a finite-dimensional vector space of
**discrete fields**"*, carrying a `shape` — the tensor shape of the DOFs. So
`L²(S²)` is the space the basis functions are **elements of**, never the space
they are **maps from**. `FunctionSpace` answers *"what do these live in"*; the
domain question is *"what do these eat"*. Different arrows.

There are **three** levels, and the tree currently names two:

| level | the object | `Basis` has | `DiscreteMeasure` has |
|---|---|---|---|
| the manifold `M` | `S²`, `[-1,1]`, `S²/σ_y`, `ℝᵈ`, energy, index | ⛔ **nothing** | `support` (a bare `str`) |
| fields on `M`, discretized | `L²(M)` at `N` nodes | — | `.space` (`FunctionSpace`) |
| coefficients | `ℝ^K` | `.space` (`FunctionSpace`) | — |

⭐ **And the `FunctionSpace`-level check already passes** — which is why the
defect survived. `[M]` the frame's arrow `analysis : measure.space →
basis.space` type-checks on the slab: shapes `(8,3) → (8,)` are fine. **The
check that fails is one level down, on the manifold, where no object exists.**

⚠ **The measure's own tag is mis-NAMED the same way.** `[M]`
`gauss_legendre(8).measure.support == '[-1,1]'` while `supp(μ)` is 8 points — so
`support` already names the **ambient manifold**, not the support. That is what
the `SPACE_CIRCLE` rename lesson (`measure.py:120-128`) says it learned; the
name simply did not follow.

### ⛔ The `support` rename of 2026-08-31, REFUTED the same day

§II.15 R5 ruled `Basis.support` on a measured name collision (58 `domain`
definers, 13 `getattr(x,"domain",None)` readers). **Both halves were wrong.**

* `[M]` **the collision is unreachable.** All 13 readers read off OPERATORS
  (`self.inner`, `self.a`, `self.b`, `self.op`, `block`), and `[M]` **0 `Basis`
  subclasses are operators**. ⚠ This is `plan-authoring` §8's brand-new
  sharpening — *a name's EXISTENCE is not its RELEVANCE* — applied by me to the
  agent's hazard and **not** to its naming evidence, one commit later.
* And `support` is **mathematically false for a basis**: `supp(f)` means *where
  `f` is non-zero*. For `IndicatorBasis` that is exactly ONE cell per function,
  so `IndicatorBasis.support = "spatial_R1"` states something untrue. For the SH
  basis it is accidentally near-right (almost all of `S²`), which is what let it
  past.

⟹ **RESTORED: `Basis.domain`.** Category-theoretically it is not even a
collision — `dom` is the source of a morphism in both cases, in **Man** for a
basis function and in **Vect** for an operator. Same functor, different
categories. The two objects legitimately keep different words for the same
manifold: *support of a measure* and *domain of a function* are both standard.

### ⭐ The manifold is currently smuggled through a NAME STRING

`[M]` two sites encode `M` inside a `FunctionSpace.name`:

| site | name built as | |
|---|---|---|
| `measure.py:331` | `f"L2[{self.support}]"` | derived — but from a `str` |
| `basis/indicator_basis.py:284` | `f"L2[coarse_cells_R{self.ndim}]"` | ⛔ **hard-coded, and already FALSE** |

⭐⭐ The second is **§II.15 R1's defect, live in the shipped tree**:
`EnergyGrid.as_basis()` (`energy_grid.py:220`) builds
`IndicatorBasis(edges_per_axis=(arange(n_groups+1) - 0.5,))` — an **index**
partition — and that basis then called its own coefficient space
`L2[coarse_cells_R1]`, naming a **spatial** manifold it has nothing to do with.
So R1 was not a hypothetical about a property I was about to add; the same lie
already shipped one level over, in a name string.
✅ **REMEDIED 2026-09-01 by tracker 2.1 @ `c461fe8d`**: the constructor takes
`partition_of=EnergyGroups(n_groups)` and the name derives —
`L2[coarse_cells(energy)]`. ⚠ The MEASURE's `f"L2[{support}]"` still
interpolates a `str`; `FunctionSpace` carrying its own `Manifold` (which
collapses both spellings) is **2.0c**.


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

> ⛔ **VOID as written since 2.4 (2026-09-01) — both premises moved.** The
> vocabulary is `Manifold` OBJECT equality (2.0c), and the slab measure says
> `S^2/SO2_x` (2.4), so "those name the same manifold" is now literally one
> memoised object. `[M]` at `8fc781d2` an equality G0 refuses exactly the Part
> I pairing (slab rule `S^2/SO2_x` vs SH `S^2`) and admits both shipped angular
> pairings (the fold's two halves are both `S^2/sigma_y`; Lebedev vs SH both
> `S^2`); the one legal pairing it would refuse — `P_ℓ` on a full sphere rule
> — ships nowhere. Whether G0 is equality or ambient-base + the lattice below
> is a RULING owed at 2.2's opener; see **§V.5h(a)**. The lattice predicate
> below stands as G2, and its MEASURE side is real since 2.4.

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

> ✅ `[M]` 2026-09-01 (2.1b): this table RUNS as
> `tests/numerics/test_basis_domain.py::test_e5` on real objects, and
> reproduces all four verdicts. ⚠ Two spellings differ from the table's
> vocabulary: (i) the MEASURE side spells row 3's `Trivial` as **`None`** — a
> full-sphere rule has SPENT nothing, and `quotient_group` answers `None`
> rather than `Trivial` for a point set (HAS ≠ SPENT there; for FUNCTIONS
> they coincide, which is why a basis carries one slot); (ii) rows 2–3's
> `SO2` basis does not ship (3.4-R's `LegendreBasis`, `class LegendreBasis`
> **0**), so the gate uses a test-local stub declaring the domain
> `S²/SO2_x`, with a retirement trigger in its docstring.

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
   ✅ **REMEDIED 2026-09-01 by 2.4** — the `[R]` was right about the mechanism
   and the fix was not a convention: `SO2` is axis-parameterised (user
   ruling), `[M]` `SO2("x").is_invariant(slab)` is **True** and `SO2("z")`
   **False**, both DERIVED (a finite set is `SO(2)_a`-closed iff on axis `a`;
   the marginal embeds along the axis its orbit space names). The 2.4
   EXECUTED section (c) carries the table.
2. `[M]` the rules populate their symmetry fields **inconsistently and with two
   different concepts**:
   `gauss_legendre(8)`: `invariance_group=Mirror('x')`, `folded_by=None`;
   `folded_product(4,8)`: `invariance_group=None`, `folded_by=Mirror('y')`.
   `invariance_group` = the symmetry of the NODE SET **within** its support;
   `folded_by` = the QUOTIENT group. **The slab's quotient group is recorded
   nowhere.** Anything built on these must not conflate them.
   ✅ **REMEDIED — `folded_by` retired at 2.0c, the slab's quotient group
   recorded at 2.4 (2026-09-01).** Both concepts now live in ONE mechanism:
   the quotient group is `Quotient.by`, read through `measure.quotient_group`
   for the discrete fold (`Mirror('y')`) and the continuous marginal
   (`SO2('x')`) alike, and `invariance_group` keeps its one meaning.

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


# Part V.5 — ⭐ The minting inventory (reviewed 2026-08-31)

What the math needs as an OBJECT, what it needs only as a FIELD, and what is
still an open architectural question. Written at the user's request after D0.7.

## A. Types to mint

| # | type | why it is EARNED (≥2 realizations + a non-identity morphism) | phase |
|---|---|---|---|
| 1 | ⭐ **`Manifold`** | D0.7 — **RULED**. 10 shipped realizations; quotient / product / pushforward all shipped as morphisms | 2.0a |
| 2 | **`Chart`** — a typed arrow `M → N` | §III.9. Makes `support(φ_*μ) = codomain(φ)` derived, kills the `SPACE_SPHERE` literal in `product()` (L7), and ⭐ **makes the phantom lift unspellable**: there IS no chart `[-1,1] → S²`, because a point of `[-1,1]` is an ORBIT | 2.3 |
| 3 | **`ReynoldsProjection`** `ℛ = π* ∘ 𝔄` | already planned; `ℛ² = ℛ` is its own law | 3.2 |
| 4 | **`OrbitAxis`** + retraction/section | already planned; re-uses `AxisRetractionOperator` | 3.3 |
| 5 | **`LegendreBasis`** on `[-1,1]` | 3.4-R **RULED (b)** — non-isomorphic to `{Y_ℓ^m}` (`L+1` vs `(L+1)²`), and the quotient projection is the non-identity morphism | 3.4 |
| 6 | **`IsotypicDecomposition`** | already planned; the trivial block IS 3.4 generalised | 7.1 |
| 7 | ⭐ **`Descent`** — the witness of `Funcs(M)^H ≅ Funcs(M/H)` | **3.4-R2 RULED (i)**, §C below. Two realizations of one space ship after 3.4-R; without this they are a Cardinal-Rule-2 twin. `[R]` it is Phase 7's object arriving early — the trivial isotypic block IS a descent, and higher irreps generalise it | 3.4b |

⚠ **`Manifold` has one requirement the plan must not lose: a QUOTIENT knows its
SINGULAR STRATUM.** §III.7 and Phase 5.1 want `OrbitAxis` to carry the
`μ = ±1` stratum *"by derivation from `H`"*, not by hand. That is a property of
`M/H` — the image of `H`'s fixed-point set — so `Manifold.quotient(H)` must
return a manifold that knows it. ⭐ And `[M]` the detector already exists,
pointed at the wrong question: `_evaluate_real_sh`'s `on_axis` guard fires on
exactly those 2 nodes of 1 rule (§II.1).

## B. Fields, not types — the type exists and the SLOT does not

| # | field | today | phase |
|---|---|---|---|
| 7 | `Basis.domain: Manifold` | absent (L2) | 2.1 |
| 8 | `Basis.invariance_group: SubgroupOfO3 \| None` | ~~absent~~ ✅ **REMEDIED 2026-09-01 by 2.1b** — DERIVED from `Basis.domain` (a `Quotient`'s `by`; `Trivial` off the bare sphere; `None` off a non-angular manifold) by one concrete `@final` property on the ABC. The requested FIELD dissolved into `domain.by` exactly as row 9's dissolved into `Quotient.by` | 2.1b |
| 9 | `measure.quotient_group: SubgroupOfO3 \| None` | ~~`[M]` **absent** — `folded_by` is the discrete case on the wrong object, and Part IV records *"the slab's quotient group is recorded nowhere"*~~ ✅ **REMEDIED** — derived from `Quotient.by` at 2.0c; the slab's answers `SO2('x')` since 2.4 (`17501245`) | 2.0d → 2.0c + 2.4 |
| 10 | `FunctionSpace.manifold` | smuggled through a NAME STRING at 2 sites, one of them `[M]` already FALSE (§III.10) | 2.0c |
| 11 | ⭐ **`AngularSymmetry` needs its Γ slot** | `[M]` `support` derives from `continuous_isotropy` ALONE, so `folded_product`'s `S^2/sigma_y` matches **no** geometry and stage 0 would refuse the shipped cylinder (§II.10) | **2.2 — NEW, was unowned** |

⭐ Note 7+8+9 together are exactly Part IV's G0 predicate
(`measure.quotient_group ⊆ basis.invariance_group`), which the plan has stated
since it was written and for which **neither side exists in the tree**. The
predicate was never wrong; it had no operands.

## C. ⚠ THE OPEN ARCHITECTURAL QUESTION — the descent isomorphism

`Funcs(M)^H ≅ Funcs(M/H)`: the `H`-invariant functions upstairs ARE the
functions downstairs. After 3.4-R the tree will carry **both realizations of
that one space**, and Cardinal Rule 2 says two spellings of one concept is a
flag until the relation is explicit:

| realization | example | domain | span |
|---|---|---|---|
| **upstairs** — invariant SUBSPACE | `MirrorEvenSphericalHarmonicBasis ⊂ SH(S²)` | `S²` (unchanged) | restricted, σ-odd columns zeroed |
| **downstairs** — the quotient's OWN basis | `LegendreBasis` on `[-1,1]` (3.4-R) | `S²/SO(2)` | full |

`[R]` **why they differ is practical, not mathematical**: `S²/SO(2)`'s function
space has a classical named basis (Legendre); `S²/σ_y`'s does not, so the
invariant-subspace realization is the only spellable one there. The `≅` is what
Q1.4's one-line *"`{Y_ℓ0} ≅ {P_ℓ}`"* has been carrying silently all along.

**Three options, unruled:**
* **(i)** mint a `Descent` / `InvariantRestriction` witnessing the iso, so both
  realizations are one object's two faces. Phase 7 needs the general form
  anyway — the trivial isotypic block IS this, and higher irreps generalise it.
* **(ii)** keep both, unrelated, with the discriminator RECORDED (*"downstairs
  when the quotient has a classical basis; upstairs otherwise"*) — cheapest,
  and honest only if written where both classes can see it.
* **(iii)** force one realization everywhere — ⛔ `[R]` refused on sight: (a) is
  option 3.4-R(a), already ruled out; (b) requires inventing a basis for
  `S²/σ_y` that has no classical form.

⟹ ✅ **RULED 2026-08-31 (i) — mint `Descent`.** It is Phase 7's object arriving
early: the trivial isotypic block IS a descent, so building it here is not a
detour but the first instance of machinery Phase 7 needs regardless. Both
realizations then become one object's two faces rather than siblings with no
stated relation.

**What `Descent` must carry** (`[R]`, to be checked against §III.8 when 3.4b is
designed):
* the pair `(M, H)` and the resulting quotient `M/H`;
* the two realizations — the **invariant subspace** upstairs and the
  **quotient's own basis** downstairs — and which one a given consumer gets;
* the isomorphism itself, checkable: analysing upstairs then descending equals
  descending then analysing, to machine precision;
* ⚠ the **discriminator**, written where both classes can see it: *downstairs
  when the quotient has a classical named basis (`S²/SO(2) → {P_ℓ}`), upstairs
  otherwise (`S²/σ_y` has no classical family)*. That sentence is the whole
  reason two realizations exist and it must not live only in this plan.

⚠ Scheduling: **not on the exit path** (Part XII), so it does not block the
return to Campaign 2 — but it must land **with or before 3.4**, since 3.4 is
what creates the second realization.

## D. Checked and NOT needed

`SubgroupOfO3` (exists; carries BOTH continuous `SO2` and discrete `Mirror`
members — 2.0d needs no new group type) · `GeneratingMeasure`/`ReferenceMeasure`
(exists) · the well-posedness triple (a PREDICATE — gates G0–G3 — not a type) ·
⛔ `Space` as a phantom generic (refused by `sn_reshape.md` Issue 2 **and** by
`coding-standards`' phantom-parameter corollary; D0.7 proposes a value type
instead, which is a different thing).

## E. Placement

`Manifold` is the most primitive object in the stack — `measure`, `basis` and
`space` all depend on it and it depends on none of them. ⟹ its own module with
no intra-`numerics` imports. ⚠ `plan-authoring` §6d: check the import edge BOTH
ways and read `tests/test_layer_imports.py`'s `FORBIDDEN_EDGES` before placing.



## §II.16 ⛔ 2.1-W's own premises, corrected at execution (2026-08-31)

The witness LANDED, and three of the claims that motivated it did not survive
being run. All three were mine or inherited; the *obligation* was right and its
*evidence* was wrong in three separate ways.

### ⛔ (a) "the fold basis has ZERO committed witnesses" — FALSE, there were four

`[M]` under the defining rebind (`_harmonic_basis` → `SphericalHarmonicBasis`),
**four pre-existing tests red**:
`test_quadrature_fold.py::TestFoldedHarmonics::{test_flat_moments_are_the_isotropic_moment_alone,
test_the_folded_frame_analysis_is_isotropic_on_a_flat_flux}` and
`tests/numerics/test_frame.py::{test_parseval_frame_square_closes,
test_parseval_dressing_installed_on_diagonal_frames}[folded8x8-L2]`.

⟹ **2.1-W stands on a narrower and better-stated premise.** All four are
**OBJECT-tier** gates — on the basis and on its frame — and `[M]` **none is
reached through `solve_sn`**. What had no witness was the **eigenvalue tier**:
nothing connected the sub-basis to an answer a user reads. That is the gap the
new gate closes, and it is a smaller claim than "zero witnesses".

### ⛔ (b) The `0.7475027… → 0.6095089…` (18.5 %) catcher figure does not reproduce

`[M]` on the fixture the row itself names, the measurement is
`0.9726641733732218 → 0.4159228684117852` — **57 %**, not 18.5 %. The original
configuration is unknown and was never written down (`plan-authoring` §4: a
relayed number travels with its configuration or it does not travel). **Neither
pair of numbers was published into the corpus**; the gate carries its own.

### ⛔⛔ (c) "deleting the class reds 0 of 1913" is not a collection-ERROR artifact — **pytest never ran at all**

I suspected `vv` Mode-8's third pipeline class (a mutation that raises kills
COLLECTION, and a scanner counting `^FAILED` reads 0). `[M]` it is one notch
worse. `directional.py:83` imports the class at module scope and
`tests/sn/primitives/conftest.py:7` reaches it transitively
(`conftest → tests.sn._test_helpers:30 → orpheus.transport:77 → …fields:49 →
cross_section_field.py:68 → orpheus.numerics:37 → …quadrature:71 →
.directional:83`). Reproduced with a `sitecustomize` meta-path hook: **rc=4,
0 collected, 0 `^FAILED` AND 0 `^ERROR`** — *both* scanners read zero, and
`--continue-on-collection-errors` does not help, because there is nothing to
continue from.

⟹ the honest re-measure, with its scope and filter written down. Population:
**80 test files** that can reach a folded rule — 74 by direct
`folded_product`/`.quotient` construction, **plus 6 indirect-only**, because
production's own MMS case builders (`derivations/continuous/mms/sn.py:2104`,
`:3873`) *default* to `folded_product`; one of the 6 is outside `tests/sn` and
`tests/numerics` entirely. Filter
`-O -m "not slow" -q --continue-on-collection-errors`:

| run | result |
|---|---|
| control | **1827 passed**, 35 deselected, 46 xfailed, 0 failed, 0 error, rc=0, 354 s |
| rebind | **7 failed**, 1820 passed, 0 error, rc=1, 358 s |

### ⭐ And three findings worth carrying past this step

* **A folded-vs-unfolded equivalence gate is UNWRITABLE on a cylinder.**
  `[M]` `Quadrature.product(4,8)` raises *"A cylindrical SNMesh admits only a
  quadrature whose every mu-level is CARRYING"* (the Q5.6.3 flip) — the
  unfolded parent cannot be posed on the chart its child serves.
* ⛔ **A CACHE can mask this mutation into a plausible bit-identical green.**
  `Quadrature._angular_frames[L]`: a rule built *after* the mutation reads
  `0.4159…`; the same rule warmed by an unmutated solve first reads `0.9727…`,
  `array_equal` to honest. ⟹ build after, or install at `pytest_configure`.
* ⭐ **A green mutation arm can be a GEOMETRY THEOREM, not a blind gate.**
  Over-masking slot `[1,0]` bit (`max|ΔY| = 8.611e-01`) and moved nothing:
  that slot carries `μ_z` (`corr = +1.000`), and a 1-D cylinder is symmetric
  under `μ_z → −μ_z`, so the basis function is annihilated by the chart at any
  refinement. A third hypothesis beside *insufficient mutation* and *blind
  gate* — and the only one of the three that is a fact about the geometry
  rather than about the test.

## §V.5b ⛔ 2.0a's migration sizing, re-censused at the phase opener (2026-08-31)

The 2.0a row priced the `Space = str` retirement as *"29 `Space` refs, 39
`support=` construction sites, 45 `.support` reads"*. **No predicate and no
tree was written beside any of them**, and re-derived they do not reproduce.

`[M]` by **AST** over `orpheus/` + `tests/` (a membership question is parsed,
never grepped — `plan-authoring` 2026-08-27; positive control: the known
`directional.py` forge site is found; 0 files unparseable):

| predicate (the AST node, exactly) | `orpheus/` | `tests/` | total |
|---|---:|---:|---:|
| `ast.keyword` with `arg == "support"` | 29 | **58** | **87** |
| `ast.Attribute` with `attr == "support"` | 31 | **31** | **62** |
| `ast.Name` with `id == "Space"` | 7 | 1 | **8** |
| `ast.Name` with `id.startswith("SPACE_")` | 25 | 14 | **39** |

⟹ **the real surface is 196 name-level touches, of which 104 — 53 % — are in
`tests/`.** The original row is not merely low; its three numbers answer three
*different* questions, and only one of them (`support=` at 29) reproduces at
all, and only on `orpheus/`.

⚠ **What actually changes.** Size is not a veto (D0.5), so this is not a
blocker — but the SHAPE moved, and the shape is what a step decomposition is
built from (`plan-authoring` §6b: the unit of work is the call-site set, not
the tidiness of the description). The dominant cost of 2.0a is **58 test
construction sites**, not the production tree. A decomposition that plans
`orpheus/` first and leaves `tests/` as "follow-up" would leave the tree
un-compiling between steps.

⭐ **The transferable tell, and it is new: a plan row that lists SEVERAL counts
in one breath implies a shared denominator that was never stated.** §2's
quantifier clause guards a *single* number's missing denominator; a *list*
(`29 … 39 … 45`) reads as one census reported three ways, so the reader
supplies the shared scope by assumption. Here there was none — and one of the
three (`39`) is *exactly* the answer to a neighbouring predicate the row does
not name (`SPACE_*` names, 25 + 14). `[R]` that is consistent with a label and
a measurement paired by proximity rather than by derivation; I cannot prove it,
and the coincidence is what made the row look internally consistent.

⟹ when a row prices work with more than one count, each count carries its own
predicate and its own tree, or the row is one census wearing three hats.


## §V.5c ⭐ 2.0a's REAL ground measure — the algebra is already there, spelled as string concatenation

§V.5b corrected the *size*. This is the *shape*, and it is what the mint should
be argued from — because a count of `support=` sites measures the migration's
bulk, while the sites below measure whether the type is EARNED.

`[M]` all by AST, 2026-08-31, over `orpheus/` + `tests/`.

### (a) A `SPACE_*` constant retypes for free; a raw string does not

`support=` kwargs, split by the ast node of the VALUE:

| the value passed | `orpheus/` | `tests/` | total |
|---|---:|---:|---:|
| raw string literal | 3 | **51** | 54 |
| `SPACE_*` constant | 13 | 3 | 16 |
| attribute (`m.support`, forwarded) | 5 | 3 | 8 |
| other `Name` (forwarded) | 4 | 1 | 5 |
| **f-string** — a tag BUILT by interpolation | **4** | 0 | 4 |

⟹ **`Space = str → Manifold` leaves the 16 constant sites and the 13 forwarding
sites untouched**; the work is 54 literals (94 % of them in `tests/`) and the 4
f-strings. The production tree is nearly free. This is the number the step
decomposition needs and §V.5b's `87` is not.

### (b) The member list in the 2.0a row is INCOMPLETE — the families are parametric

The row enumerates `S²`, `S¹`, `[-1,1]`, `[0,1]`, `[0,∞)`, `ℝ`, `spatial_Rᵈ`,
`energy`, `index(axis)`. `[M]` the four f-string producers show the shipped
members are **parametric families**, and one family is missing from the row:

| site | builds | family |
|---|---|---|
| `generating_measure.py:420` | `f"[{a},{b}]"` | the interval `[a,b]` — `SPACE_INTERVAL_M11`/`_01` are two MEMBERS of it, not two manifolds |
| `material_mesh.py:542` | `f"spatial_R{self.ndim}"` | `ℝᵈ`, indexed by `d` |
| `frame.py:759` | `f"index({axis_label})"` | a finite index set, by label |
| `loss_kernel_gauge.py:1169` | `f"sn_trace_orbit{orbit}_g{group}"` | ⭐ **the same index-set family, second spelling** — `nodes = indices.astype(float)`, a subset of trace DOFs per (reflection orbit, group). Not in the row's list |

⟹ two shipped sites mint the SAME kind of object under incompatible naming
schemes — Pattern 2, in the string vocabulary. And `Manifold` cannot be an enum
or a fixed set of singletons: it must be a value type with **parametric
constructors**, which is what D0.7 ruled and this measures.

### (c) ⭐⭐ Three of the plan's proposed VERBS already exist, as string ops

`[M]` **18** `.support` reads do string manipulation — **all 18 in `orpheus/`,
0 in `tests/`.** Three of them are the morphisms 2.0a proposes to mint:

| site | today | the verb it IS |
|---|---|---|
| `measure.py:588` | `new_space = f"{self.support} × {other.support}"` | `Manifold.__mul__` |
| `measure.py:1022` | `new_space = f"{self.support}/{group.name}"` | `Manifold.quotient(H)` — ⛔ **my "no lookup and no check" is HALF FALSE, corrected 2026-08-31.** The MEASURE *is* gated: `orbit_certificate` raises unless the nodes are `group`-invariant. What is ungated is the **TAG** — no catalogue lookup, so `[M]` `DiscreteMeasure(support="S^2")` accepts a node of norm `√2`. ⭐ The sharper statement is the campaign's own thesis: the measure's NODES are checked, its MANIFOLD CLAIM is not — level 2 gated, level 1 forged |
| `measure.py:802` | `f"φ_*({self.support})"` | the pushforward's codomain |
| `measure.py:331` | `name = f"L2[{self.support}]"` | §III.10 / 2.0c's derived name |

⟹ **the mint is not adding an algebra; it is giving a name to one that already
runs.** `coding-elegance` Pattern 1: the string interpolation IS the operation.
That is the strongest available form of the D0.7 criterion's clause (b) — the
morphisms are not merely *applicable*, they are *shipped*.

### (d) ⛔ A production branch arm with no producer — found by the same census

`measure.py:411`: `if self.support.startswith("spatial") or self.support == "cells":`

`[M]` **`"cells"` is produced by NO static producer.** Denominator, enumerated:
58 literal `support=`/`new_space=` sites (13 distinct values), the 6 `SPACE_*`
constants, and 36 non-literal producers — of which 4 are the f-strings in (b)
and 32 forward an existing support. `"cells"` appears in none. The only
`"cells"` elsewhere in the tree is `radial_characteristic_space.py`'s
`part="cells"`, a **walk-part label, not a measure support**.

So the arm is unreachable, and **three surfaces advertise it**: the branch, the
docstring at `:387` (*"(`"spatial_…"` / `"cells"`)"*), and the `raise` message
at `:416`. ⚠ `vv-principles` #17's granularity trap is why nothing caught it —
mutating the `phase` guard as a unit reddens on the `spatial` arm, which has
producers, and certifies both.

⟹ **acceptance item for 2.0a**: the mint must make this arm unspellable (a
manifold's phase is a property of the object, not an OR of two string prefixes,
one fictional), and the docstring and message must stop naming `"cells"`.

### (e) Two spellings of one quotient, both shipped

`[M]` the literal inventory contains **`'S^2/<sigma_y>'` and `'S^2/sigma_y'`** —
the same quotient, two tags, unequal under `==`. Also shipped: `'img'`,
`'probe'`, `'[-1,1]^slab'`. This is `measure.py:114`'s *"recommendations, not
constraints"* (L8) with its bill arriving: the type permits any string, so the
tree contains manifolds that are typos of each other and tags that name nothing.


## 2.0c EXECUTED (2026-09-01) — a measure names the point set it lives on

**Goal (outcome).** A measure states *where its atoms live* as an object with
structure, so the phase-space factor it belongs to, the space it induces, and
the orbit space it folds onto are all DERIVED rather than spelled.

### What landed

| # | the change | why it is not merely a retype |
|---|---|---|
| 1 | `Space = str` and the six `SPACE_*` tag constants **RETIRED**; `support: Manifold` on all six implementors (`DiscreteMeasure`, `GeneratingMeasure`, `UniformMeasure`, `ProductMeasure`, the `ReferenceMeasure` **Protocol**, and `AngularSymmetry` — the one spelled `-> str`) | the tree carried **two parallel vocabularies** for one concept (§V.5f(b)); this collapses them |
| 2 | ⭐ `DiscreteMeasure.phase` — stringly-typed dispatch → a `match` on the manifold's TYPE (user-ruled, §V.5f(e)) | `startswith("spatial")`/`== "cells"`/`startswith("energy")` gone; `manifold.py` stays pure geometry |
| 3 | ⛔⛔ **the LIVE defect fixed** — `folded_product(4,8).measure.phase` RAISED, now answers `angular` | a `Quotient(base=Sphere())` arm; the fold's `invariance_group` is legitimately `None` |
| 4 | ⭐ `μ.quotient(G)` routes through `Manifold.quotient` | the fold's target is the **catalogue's** orbit space, not a fabricated `f"{support}/{G.name}"` |
| 5 | ⭐ `μ ⊗ ν` and `ProductMeasure.support` route through `Manifold.__mul__` | the product algebra shipped as `" × ".join(...)` and `f"{a} × {b}"` — two hand-rolled spellings of one operator |
| 6 | ⭐ `pushforward(new_space=...)` is **REQUIRED**; the fabricated `f"φ_*({support})"` default is gone | the image of a manifold under an arbitrary map is a manifold nobody derived — naming it does not make it known. `[M]` 7 of 8 call sites already passed it |
| 7 | ✅ **2.1's handoff discharged** — `LossKernelBasis`'s measure now reads `support=basis.domain` | the one production frame whose halves disagreed in spelling; closed by construction, not by keeping two strings equal |
| 8 | `Mesh1D`/`Mesh2D`/`EnergyGrid`/`material_mesh`/`frame.py` — every f-string producer takes the manifold its own basis partitions | four of five build basis and measure in one function; `frame.py` now names the manifold ONCE and hands it to both halves |
| 9 | 2.0d absorbed: `DiscreteMeasure.quotient_group` derived from `Quotient.by` | `[M]` it agrees with `Quadrature.folded_by` — two independent records of one fact |

### ⭐ Four things the TYPE found that the string could not

Each was green before the retype and is a real defect, not a migration cost.

1. ⛔ **`EnergyGrid.as_measure` and `.as_basis` named DIFFERENT manifolds.**
   `[M]` the measure said `EnergyGroups(ng=None)`, the basis `EnergyGroups(ng=2)`.
   Under the tag both said `"energy"` and **agreed by being equally
   uninformative** — the string was lossy, and the loss is what hid the
   disagreement. This is precisely `test_d6`'s claim, caught on a pair 2.1 had
   already touched.
2. ⛔ **A test fixture called `(μ, 0, 0)` on `"[-1,1]^slab"` "the slab
   embedding" and folded it by σ_y.** Those points lie on **no** manifold with
   a σ_y action (`(-0.7,0,0)` is not on `S²`; an interval has no `O(3)` action).
   The string fabricated `"[-1,1]^slab/sigma_y"`; `Manifold.quotient` refuses.
   ⟹ ERR-080's defect class *inside a test that was green*. Replaced by the
   genuine σ_y-fixed set — the `xz` great circle — which is a stronger claim.
3. ⛔ **That same test contradicted its own docstring.** The header says *"the
   quotient CHANGES the space, so re-folding is ill-posed"*; the body then
   re-folded and asserted idempotence. `[M]` `(S²/σ_y).quotient(σ_y)` is
   correctly refused. `coding-elegance` #20 — a docstring naming a claim the
   body does not make — resolved by making the body assert the docstring.
4. ⛔ **`test_angular_phase_is_symmetry_derived_not_support_tagged` INVERTED.**
   It asserted that stripping `invariance_group` from a Lebedev rule makes
   `phase` RAISE, because `"S^2"` was an inert string no branch read. `Sphere`
   is a *reason*, so a measure on one is angular regardless — which is exactly
   what lets the shipped fold answer. Renamed and re-argued in place, per §3.

### ⚖ Gates: one promoted, one retired, one re-based

* ⭐ **PROMOTED** — `test_manifold.py`'s two-producer agreement gate compared
  `Manifold.quotient(...).realization` against `AngularSymmetry.support`. Both
  sides are now `Manifold`, so it asserts **object equality** where it could
  only assert **names** before. `coding-standards`' mirror clause: a retirement
  can silently strengthen a gate's claim class, and the docstring must move
  with it. ⚠ It is **not** tautological — `AngularSymmetry.support` still
  hand-writes its table while `Manifold.quotient` derives the orbit space from
  invariant theory; neither reads the other.
* **RETIRED with its arm** — `test_spatial_supports_have_spatial_phase`'s
  `"cells"` parameter. `[M]` §V.5e(e): zero production producers, and its only
  producer was *that parametrize list* — a synthetic witness for the arm it
  witnessed. The replacement runs `d = 1, 2, 3`; `d = 3` passes with no
  production change, where the prefix test needed a new tag per dimension.
* **RE-BASED, principled** — `LEGENDRE.on(0, 1).support` was the string
  `"[0.0,1.0]"` (an f-string over floats) and is now `Interval(0.0, 1.0)`,
  whose `name` normalises to `"[0,1]"`. Comparing the OBJECT says the interval
  is the same one whatever the float repr — which the string could not.

### `[M]` The mutation battery — 9 arms, and the one that came back BLIND

Crash-safe by construction (every touched file copied aside *before* the first
mutation and diffed against the copy afterwards — a `finally` does not survive a
SIGTERM, a copy does; all 5 files verified clean).

| arm | reds | what it pins |
|---|---|---|
| A | 1 | `phase`'s `Quotient(base=Sphere())` arm — the LIVE defect's fix |
| B | 1 | `phase`'s `Sphere()` arm |
| C | 3 | `phase`'s `RealSpace()` arm |
| D | 3 | the fold routes through `Manifold.quotient` |
| E | 1 | the tensor product routes through `Manifold.__mul__` |
| ⛔ **F** | **0 → 6** | the loss-kernel measure reads `basis.domain` |
| G | 1 | the energy measure carries its group count |
| H | 1 | `AngularSymmetry`'s table vs the derived orbit space |
| ★ | 7 | **positive control** — `Manifold.name` lies |

⛔⛔ **Arm F came back BLIND, and it is the one thing 2.1 explicitly handed this
step.** 2.1's `test_d6` pinned the `LossKernelBasis` spelling divergence and its
docstring said 2.0c *"must come back here and **cannot resolve it by
accident**"*. 2.0c resolved it — `support=basis.domain` — and **that is exactly
what resolving it by accident looks like**: the divergence became unspellable,
`test_d6` went green for a different reason, and replacing the expression with a
wrong `IndexSet` reddened **nothing tree-wide** (re-measured against the
loss-kernel suite too, to separate *no witness* from *not measured* — genuinely
zero). `test_d6` asserts the BASIS's half; no gate reached the production
MEASURE's.

⟹ the witness now exists —
`test_each_blocks_frame_names_ONE_manifold_on_both_halves`, 6 rows over the
singular fixtures — and it compares the two halves **through production
plumbing** (`block.gather.codomain.name == f"L2[{basis.domain.name}]"`, the
gather's codomain being the measure's induced space) rather than by re-reading
the constructor. `[M]` the arm now reddens 6 of 6.

⭐ The transferable shape, and it is `vv-principles` #17 read in an unfamiliar
direction: the battery is normally run to find gates with no teeth. Here it
found a **repair with no gate** — a change that is correct, that makes a defect
unspellable, and that nothing would notice being reverted. A prediction ("2.0c
cannot resolve this by accident") does not enforce itself; only a mutation does.

### `[M]` Bit-identity of every shipped space name

Before touching a call site: **10 of 10** tag constants reproduce their
spelling exactly through `Manifold.name` (`S^2`, `S^1`, `[-1,1]`, `[0,1]`,
`[0,inf)`, `R`, `energy`, `spatial_R1`, `spatial_R2`, `index(angular)`), and so
do the two derived spellings (`S^2/sigma_y` via `Quotient.name`,
`[-1,1] × [-1,1]` via `Product.name`). The single divergence is the affine
remap above.

### §6d — the new module-scope edge, injected and RUN

`measure.py` now imports `manifold.py` at module scope. `[M]` safe **by design,
not by luck**: `manifold → symmetry` is `TYPE_CHECKING`-only, so the cycle
`measure → manifold → symmetry → measure` never closes at runtime. Injected and
run in **five** import orders (measure / symmetry / geometry / registry /
sn.solver first) — all clean. Both modules are `orpheus.numerics.*`, so no
package-to-package edge is created and the layer contract is untouched.

### ⚠ Two process notes worth carrying

1. ⛔ **L25, at tree scale.** The vocabulary sweep was a tree-wide regex, and it
   corrupted **history prose in `manifold.py` itself** — `"was
   SPACE_INTERVAL_M11"` became `"was COSINE_INTERVAL"`, i.e. the file defining
   the replacements had its own record of what it replaced overwritten. Caught
   by reading the diff, not by any test (prose has no gate). L25 is written for
   a *file-internal* `replace_all`; the same premise — "every occurrence here is
   the target concept" — fails identically across a tree, and the file most
   likely to hold the counterexample is the one that documents the rename.
2. ⭐ **The census's own positive control failed, correctly, at the end.** It
   asserted *"an f-string producer exists"* — true when written, and false once
   the retype retired all four. `lessons` **L62** exactly: a control drawn from
   the TREE validates a filter only for shapes the tree still exhibits. A
   control that fails because the work succeeded is the control working.

---

## §V.5g ⭐⭐ 2.4's PRE-FLIGHT (2026-09-01, at `8cc69e7f`) — it repairs a LIVE defect, and it is the one step that is NOT name-preserving

The **12th** consecutive opener to correct its own section, and it changes 2.4's
CLASS: the tracker files it under wiring (*"the slab rule declares its quotient
group"*), and it is a defect fix of exactly 2.1's kind.

### (a) ✅ The premise holds, measured rather than inherited

`[M]` at `8cc69e7f`: `SPHERE.quotient(SubgroupOfO3.SO2)` →  name `S^2/SO2`,
`realization` **`Interval(-1.0, 1.0)`**, and `== COSINE_INTERVAL` is `True`;
`dim` 1 on both. `Quotient.contains(μ-nodes)` is `True`, so the declaration
changes **no coordinate** — only what the object knows about itself.

### (b) ⛔⛔ THE FINDING: the `[-1,1]` space collision is LIVE, and it is 2.1's defect one level up

`[M]` an 8-node slab **angular** rule and an 8-node **spatial** measure on
`Interval(-1,1)`:

| | support | space | |
|---|---|---|---|
| angular (`gauss_legendre(8)`) | `[-1,1]` | `L2[[-1,1]]`, shape `(8,)` | |
| spatial (a measure on `[-1,1]`) | `[-1,1]` | `L2[[-1,1]]`, shape `(8,)` | |
| | **supports `==`** | ⛔ **spaces `==` AND hash-equal** | |

`FunctionSpace.__eq__` is `(name, shape)`, so **an angular space and a spatial
space are the same object to every composability check.** This is precisely the
illegal state 2.1 made unrepresentable for the energy/spatial pair
(`L2[coarse_cells_R1]`) — and it survived 2.0c because both sides genuinely live
on `Interval(-1,1)`: the retype made the supports *honest*, and they are honestly
identical. Typing the support was necessary and is not sufficient; the slab's
support is still named by its **chart codomain** instead of by the orbit space it
is, which is ERR-080's defect class stated at the level of spaces.

`[M]` after the declaration: `L2[S^2/SO2]` vs `L2[[-1,1]]` — `==` is **`False`**,
and the collision is unspellable.

⟹ **2.4 is not wiring.** Like 2.1 it repairs something ALREADY WRONG, which is
the property that earned 2.1 its place at the head of the queue.

### (c) ⭐ It discharges 2.0c's surviving fallback arm — check it per §6c

`[M]` a measure on `SPHERE.quotient(SO2)` with `invariance_group=None` answers
`phase == "angular"` from the **manifold alone**. So 2.0c's ⚠-marked
`case _ if self.invariance_group is not None` arm — kept because `Interval` is
genuinely ambiguous — becomes unreachable for the slab. ⛔ **Do not delete it on
that basis alone**: `phase` is public and an `Interval`-supported measure is
still constructible by a caller. The §6c question is whether any *shipped* input
still reaches it; if none does, the honest move is to say so in its docstring,
not to remove the arm.

### (d) ⛔ It is the ONE step that RENAMES a shipped space — and 2.0c's bit-identity does not extend to it

Every previous step preserved names exactly (`[M]` 2.0c: 10 of 10). 2.4
**deliberately** does not: `L2[[-1,1]]` → `L2[S^2/SO2]` on every slab angular
space. `[M]` the literal-pin surface is small — **3** occurrences of
`L2[[-1,1]]` tree-wide, of which 1 is production's own f-string and 2 are
docstrings — but the real surface is **23** space-comparison sites in
`orpheus/`, because a `FunctionSpace` compares by `(name, shape)`. ⟹ size 2.4
by the **comparison** surface, never by the literal one; a literal census here
is the same wrong-side error §V.5f(d) caught in 2.0c's own row.

### (e) ⛔ Part IV's obstacle 1 is REAL and unresolved: `SO2` names no axis

`[M]` `symmetry.py:121` — `SO2 = "SO2"`, a bare named member, comment *"axial
rotations C_∞"*. It has **no `.axis`** (checked: absent), while `Mirror` IS
parameterised (`Mirror('x')`, `Mirror('y')`). So `S^2/SO2` does not say **which**
axis the rotation is about, and the slab's μ is a cosine against a specific pole.
`[M]` the shipped slab rule meanwhile declares `invariance_group=Mirror('x')`.
⟹ 2.4 must settle whether `SO2` becomes axis-parameterised (making
`S^2/SO2_z` ≠ `S^2/SO2_x`, and changing the catalogue key) or whether the
convention is fixed by fiat and documented. **This is a design ruling, not a
measurement** — surface it before building.

✅ **RULED 2026-09-01 (user): parameterise `SO2(axis)` like `Mirror`.** Put
to the user with the numbers below, and the recommended option chosen. What
made it a ruling and not a fiat is that the tree carries **two** poles, so no
single convention could be right: `[M]` the group catalogue's `SO2` tested
`hypot(x, y) = 0` (about **z**); the slab embeds its marginal as `(μ, 0, 0)`
(**x**); `_evaluate_real_sh` takes `cos_theta = mu_x` (the spherical-harmonic
pole is **x**); the product rules' polar factor is `μ = μ_z` (**z**); and ONE
function, `gauss_legendre_on_mu`, serves both the slab's rule (an orbit space
about x) and the product factor (about z) — `[M]` `rules_product.py:361/:696`
and `directional.py:1002` consume it, and `spherical_product` reads neither
its `support` nor its `invariance_group`. The two refused options and why:
*fix by fiat* re-ships the retired-`Z2` defect (an unnamed group "realized" on
one axis while a consumer needs another — the tree already ruled this once, on
2026-08-02); *move the slab to z* has no single pole to standardise on and
would touch the sweep, the SH basis, every `Mirror('x')` row and every slab
snapshot.

---

## 2.4 EXECUTED (2026-09-01) — the slab says what space its ordinates live on

**Landed as ONE commit with its gates, per §6b/§6c — `17501245`, pushed.** The step's own title
(*"the slab rule declares its quotient group"*) named a one-line declaration;
what shipped is the declaration **plus the two things it was silently
assuming**: that the group it declares names an axis, and that the rule it
declares on is the slab's and not the product rules' factor.

### (a) What changed, by file — production

| file | what |
|---|---|
| `numerics/symmetry.py` | `SO2(axis)` frozen tag beside `Mirror(axis)`; the bare `_NamedSubgroup.SO2` member and `SubgroupOfO3.SO2` singleton **retired**; `SubgroupOfO3.SO2(axis)` constructor, `rotation_axis` (the continuous dual of `mirror_axis`); `name` → `SO2_{axis}`, `repr` carries the axis (the walk and `_GROUP_CACHE` key on it); `_is_axis_supported(nodes, axis, atol)` — the exactness test about the NAMED axis; `_so2_contains` — the axis-dependent lattice rows, out of the enum table (`SO2(a) ⊆ SO3, O3` ∀a; `SO2(z) ⊆ D_∞h` for z alone; `C_n ⊆ SO2(z)` alone; distinct axes incomparable; no finite group but `{e}` inside); `candidate_groups` offers all three axes; ⭐ `_embedded_nodes` **reads the axis off the orbit space** — a 1-D measure on `S^2/SO2_a` embeds along `a`, column 0 only for a bare interval; `_octahedral_ops`/`_icosahedral_ops` memoised (see (e)) |
| `numerics/manifold.py` | `_sphere_mod_so2` reads the axis off the group (invariants `p1 = x_a`, `p2 = x_b² + x_c²`) — ONE derivation, three catalogue keys `SO2_x/y/z`, exactly as the mirrors; `Manifold.quotient` **memoised** (`_catalogued_quotient`, `functools.cache` on `(base, group)`); `gram` frozen to `ImmutableMatrix` in all three builders so a `Quotient` is hashable |
| `numerics/measure.py` | NEW `DiscreteMeasure.on_orbit_space(quotient)` — the same atoms read as chart coordinates of `M/H`; refuses unless `quotient.realization == self.support`; drops `invariance_group` (an embedding statement), keeps `exactness` (a chart statement). `phase`'s docstring re-tensed (it still described the **string** dispatch 2.0c retired) and the fallback arm's comment corrected — see (c) |
| `numerics/quadrature/rules_1d.py` | NEW `gauss_legendre_on_polar_orbit(n, axis)` = `gauss_legendre_on_mu(n).on_orbit_space(SPHERE.quotient(SO2(axis)))` re-tagged `Mirror(axis)`; `gauss_legendre_on_mu` stays on the chart (it is the product factor); module docstring rewritten around the two objects |
| `numerics/quadrature/directional.py` | `Quadrature.gauss_legendre` builds the declared x-marginal |
| `numerics/quadrature/registry.py` | `AngularSymmetry.support` = `SPHERE if Trivial else SPHERE.quotient(spent)` (the hand-written `COSINE_INTERVAL` row retired); `.reference` keys on `rotation_axis is not None`; slab/sphere rows spend `SO2("x")`; the `GaussLegendre1D` spec's factory is `partial(gauss_legendre_on_polar_orbit, axis="x")`; the stage-0 message drops its now-redundant `(= S^2/…)` |
| `numerics/quadrature/__init__.py`, `numerics/__init__.py` | export the adopter |

### (b) ✅ The keystone — the live defect is unspellable, and no coordinate moved

`[M]` `Quadrature.gauss_legendre(8).measure`: `support` **`S^2/SO2_x`**, space
**`L2[S^2/SO2_x]`**, `phase` angular, `quotient_group` `SO2('x')`,
`invariance_group` `Mirror('x')`; nodes and weights **bit-identical** to
`gauss_legendre_on_mu(8)`; `Quotient.contains(nodes)` True. Against an 8-node
spatial rule on `Interval(-1,1)`: spaces **unequal, hashes unequal**. Positive
control in the gate: two spatial 8-node rules on the interval ARE equal — the
ROLE is the discriminator, not the shape. Gate `test_a1` (3 rows) +
`test_a2` (3 rows), `tests/numerics/test_slab_orbit_space.py`.

### (c) ⭐ Part IV obstacle 1 — answered by DERIVATION, and the axis is load-bearing at stage 0

`[M]` `SO2('x').is_invariant(slab GL8)` **True**; `SO2('z')` **False**;
`Mirror('x'/'y'/'z')` True; `Dinfh` False. The derivation, one line: a finite
set is `SO(2)_a`-closed iff every node lies ON axis `a`; the marginal embeds
along the axis its orbit space names; hence True for its own axis only. And
the z-declared rule `gauss_legendre_on_polar_orbit(8, "z")` answers the
**opposite** way (`SO2('z')` True, `SO2('x')` False), which is what makes the
axis a parameter rather than a label. `maximal_invariance_groups(slab)` →
`{SO2_x, σ_x, σ_y, σ_z}`, walk == bruteforce.

`[M]` stage 0 for `slab` and `sphere`: the declared x-marginal **admitted**;
the chart-level `gauss_legendre_on_mu(8)` **REFUSED**; the z-declared and
y-declared marginals **REFUSED**. ⟹ §6c is satisfied by inputs that exist in
the tree at landing — before 2.4 the chart-level rule *was* the registered slab
rule. `select_quadrature("slab", target_degree=7)` still chooses
`GaussLegendre1D`, now on `L2[S^2/SO2_x]`.

⛔ **The pre-flight's (c) was half wrong, and my own 2.0c comment was the
wrong half.** §V.5g(c) said the `phase` fallback arm "becomes unreachable for
the slab" — true — and the comment I wrote INTO `measure.py` at 2.0c said
*"this arm becomes unreachable"* — **false**. `[M]` `gauss_legendre_on_mu(8)`
is a shipped, public rule on a bare `Interval` carrying `Mirror('x')`, and it
reaches the arm and answers `"angular"` there; it cannot declare an orbit
space because it is the product rules' polar factor about z. The §6c check the
pre-flight prescribed found it; the arm stays, its comment now says what
reaches it, and `test_a3` witnesses both the slab answering from its manifold
alone (`invariance_group=None` → angular, `quotient_group == SO2('x')`) and
the chart rule answering via the arm. *A prediction written into a code
comment does not enforce itself* — the 2.1 `test_d6` row, one file over.

### (d) The §6b set, as the red loop enumerated it

**29** comparison sites in `orpheus/` by AST (predicate: any `Compare` with an
`Attribute` named `support`/`space`/`domain`/`codomain`/`realization`/
`measure_space` on either side; positive control `admits_domain` found;
`scratch/_p24_compare_census.py`). ⛔ The pre-flight's *23* was a different,
unstated predicate — and **neither number was the set**: 26 of the 29 compare
a FIELD's space against a space derived from the SAME measure, so they rename
together; the 3 that can see a declaration are `registry.py:938`
(`admits_domain`, stage 0), `generating_measure.py:399` (the affine remap,
untouched because `LEGENDRE` stays on the chart) and `rules_product.py:435`
(the azimuthal factor, unaffected). ⛔ The pre-flight's *"3 literal
`L2[[-1,1]]` pins"* is **0** in `orpheus/`, `tests/` and `docs/` by
fixed-string grep; the 3 were in this plan.

The set that actually reddened (`[M]` loop 1: 12 failed + 7 errors + 1
collection error over 11 files): `test_symmetry.py` (**21** `SO2` sites, one
collection error), `test_manifold.py` (4 code sites + the name pin + the
two-producer gate), `test_registry.py` (the table, the derived-support gate,
and **5 selection tests whose test-local specs built chart-level rules** —
`CHEBYSHEV_T.gauss`, `LEGENDRE.gauss` — now refused at stage 0 BEFORE the
stage they exercise; each now reads its rule on the slab's orbit space through
one helper), `test_quadrature_directional.py` (the `q8_4` support pin),
`tests/sn/_test_helpers.py` (a **false** constant `invariance_group=SO2` on a
synthetic product rule — now the honest `Dnh(n_phi)`), and one stale prose
site in `test_sweep_regression.py`. ⭐ The registry five are the §6b spelling
that no census could return: **a test-local FACTORY producing an input the
production gate now refuses** — the gate got stricter, and every fixture that
relied on the old laxness had to say what it is.

### (e) ⚠ A pre-existing hot spot, amplified 1.2× and then removed

Loop 2 did not finish in 300 s. `[M]` `test_computed_symmetry_matches_the_construction`:
one `maximal_invariance_groups` walk over `product(4,3)` took **9.4 s**, of
which **9.3 s** was `_icosahedral_ops()` rebuilding its 120-element closure
**41 times** — `_contains` asks `_realized_ops` for the generating set in its
finite-vs-finite GUARD, before the element cache in `_group_elements` is ever
consulted, so the memo on the closure never covered the guard. `[M]` with the
pre-2.4 13-candidate set the same walk paid **35** rebuilds (≈8 s): the cost
was already there; three more axial candidates raised it to 41 and the two
lattice-wide tests (the construction walk over 9 rules, the compatibility law
over every rule × 256 pairs) over the harness timeout. Fixed at the source —
`_octahedral_ops`/`_icosahedral_ops` are constants and are now `functools.cache`d;
`[M]` the walk is **36 ms**, `test_symmetry.py` runs in **15 s**. ⟹ the
transferable half: *a memo on a derived object does not cover a guard that
asks for its inputs*, and a candidate-set growth multiplies the guard, not
the memo. (Companion to `lessons` L16 — the regression hid behind green.)

### (f) `[M]` Gates at this tree state (serial, `-O -m "not slow" -p no:randomly`)

| tree | at 2.0c | at 2.4 | delta | rc |
|---|---|---|---|---|
| numerics | 2657 | **2679** | +22 | 0 |
| transport | 645 +1sk | 645 +1sk | — | 0 |
| geometry | 727 +4sk +1xf | 727 +4sk +1xf | — | 0 |
| data | 237 | 237 | — | 0 |
| sn | 3384 +1sk +116desel +50xf | 3384 +1sk +116desel +50xf | — (15 min 24 s) | 0 |
| cross_method | (not gated at 2.0c) | 81 +8desel | — | 0 |
| diffusion | (not gated at 2.0c) | 113 | — | 0 |

**Seven of the thirteen trees gated at this tree state**, all rc=0. Not
measured here: cp, derivations, homogeneous, mc, moc, and the root-level
modules other than the two docs gates — exit predicate 5 still stands and is
still best spent immediately before the fused step.

⭐ **The +22 reconciles against a second, independently counted instrument**:
the V&V matrix (`docs/theory/verification/matrix.rst`, regenerated by the
Sphinx build) went **10562 → 10584**, `foundation` 7538 → 7560, with
`numerics/test_slab_orbit_space` entering at 17 and `test_symmetry` moving
105 → 110 — the same 22 units the tree gate counted. `dead_references`:
**0 dead / 52 checked / 0 undecidable** after the corpus pass.

⭐ The sn column is the one that matters for a RENAME step, and it moved by
nothing: `[M]` 3384 of 3384, no snapshot re-baselined, no space-name pin
anywhere in `tests/sn` — the 26 field-vs-mesh comparison sites really do
rename together, and the `.npz` snapshots carry arrays, not spaces.

The +22 reconciles unit for unit: **17** rows in `test_slab_orbit_space.py`
(a1 ×3, a2 ×3, a3, b1 ×3, b2 ×2, b3, c1, c2, d1, d2), **+4** new
`test_symmetry.py` gates (rotation axis validated at construction; three axes
incomparable; `rotation_axis` the dual of `mirror_axis`; the walk offered all
three axes), **+1** parametrize row (`test_o3_contains_every_proper_subgroup`
now lists `SO2("x")` and `SO2("z")` where it listed one bare `SO2`).
pyright **0** on the six touched production files and the new gate file
(⚠ the four errors it DID report were mine — sympy's stubs type
`ImmutableMatrix.det()` loosely — fixed by taking the determinant before
freezing; the 12 remaining hits are pre-existing lines in
`tests/sn/_test_helpers.py:671–865,1175` and `test_symmetry.py:1164`, outside
every hunk of this step).

### (g0) `[M]` The per-arm mutation battery — the verdict is a TABLE (`vv-principles` #17)

`scratch/_p24_mutation.py`, crash-safe by copy-aside (5 files, all byte-identical
after), scoped to the 7 suites that can redden; baseline green first. Every
arm mutates ONE claim; every arm reddened, each with a **distinct** red set:

| arm | mutation | red |
|---|---|---|
| A | the exactness test ignores the axis (always z — the retired behaviour) | **9** — the marginal's own-axis gate, `test_b1`, the candidate walk, two SO3-on-marginal gates |
| B | the embedding ignores the orbit space's axis (column 0 always) | **5** — the two gates that declare about y/z |
| C | the derivation ignores the axis (`p1 = z` whatever the group) | **1** — `test_d1` (the generators name the coordinate) |
| D | the slab adopter declares the WRONG axis | **26** — keystone, phase, registry selection, stage 0 |
| E | `on_orbit_space`'s chart guard removed | **1** — `test_c2` |
| F | the registry answers the CHART, not the orbit space (the pre-2.4 row) | **18** — every selection gate, stage 0 |
| G | lattice: `C_n ⊆ SO2` about EVERY axis | **7** — `test_cn_in_so2` + the downward-closure law on marginals |
| H | lattice: `SO2(any) ⊆ D_∞h` | **1** — `test_so2_chain` |
| I | the walk offered only z (the retired bare `SO2`) | **1** — the candidate-groups gate |
| J | the residual re-tagged on the WRONG axis | **5** — `test_a3`, `test_b1`, two selector gates that ask the nodes |
| K | `Quotient.name` drops the axis (three quotients, one name) | **58** — 51 failed + 7 collection errors (the catalogue keys collide) |
| ★ | positive control — `Mirror` containment inverted | **6** |

⭐ Two readings worth keeping. Arms **C** and **E** redden exactly one gate
each: the derivation's axis and the chart guard have precisely one witness
apiece, which is enough and is stated. Arm **B** is the one that separates
*the axis is stored* from *the axis is used*: with the embedding blind to it,
every x-declared assertion still passes and only the y/z-declared ones fail —
a battery over the shipped slab alone would have read B as BLIND, which is why
the gates declare about all three axes.

---

## §V.5l ⭐⭐ The FUSED STEP's OPENER (2026-09-02, at `17e0eb0e`) — §V.5k reproduced 6 of 6, the RULING's evidence was one LAYER too narrow, and a PRE-STEP (2.5) is minted

The **19th** opener, run at execution. The code tree at `17e0eb0e` is
`2f294ef1`'s (`[M]` `git diff --stat 2f294ef1 HEAD` touches only this plan and
one skill), so §V.5k's numbers were expected to reproduce, and they did — the
opener still corrected its own section, for the nineteenth time. Every `[M]`
is at `17e0eb0e`: `scratch/_p31c_reground.py` re-run (`scratch/_fused_reground.log`),
`scratch/_fused_isotypic_probe.py` (+ `.log`), `scratch/_fused_g0_census.py`
(+ `.log`, and the baseline `scratch/_fused_baseline_slab_L01.npz`), the
explorer's census `scratch/_fused_6b_census.md` (+ `.py` re-runnable, `.out`;
AST with 8 asserted positive controls), and the moment-space AST census
(`SphericalHarmonicSpace` mints, in `_fused_6b_census.py`'s P3-ext).

### (a) ✅ §V.5k(a)–(f) reproduce — each with a refinement

* **(a) the fusion premise** — `support=SPHERE` by regex 7 / AST 4, the
  forgery's ONE caller, `test_q8_4` at `:568`, 3 strict-xfail rows by AST,
  `_evaluate_real_sh` 4 in-file call sites, `select_quadrature` 0 production
  callers, `solve_sn`'s rule a required positional: all as pre-flighted. ⭐
  Refinement on REAL solves (`_fused_g0_census.py`: a spy on every frame
  construction over a slab GL8 fixed-source solve at `L = 0, 1, 2`, a sphere
  GL8 at `L = 0, 1`, a folded(2,4) cylinder at `L = 0, 1`): **21 frame
  constructions, 3 unique pairings** (`SphericalHarmonicBasis` on `S^2`,
  `IndicatorBasis` on `index(angular)`, `MirrorEven…` on `S^2/sigma_y`), and
  object-equality G0 refuses **0 of 3** today — the forgery keeps the slab
  frame's measure on `S^2`. After 0.1b the slab-SH pairing is the ONE refusal.
  The slab fluxes at `[0, 0]`: `+3.999999999999995e+00` (L=0),
  `+4.000000000000002e+00` (L=1), `-3.764705882352964e+00` (L=2) — the L≤1
  rows are the bit-identity baseline the fix is measured against.
* **(b) THE FIX's blast radius** — ⛔ the ruling's evidence was one layer too
  narrow; see (b) below.
* **(c) the probe** — `[M]` REPRODUCED and re-runnable: about x, 6
  incommensurate angles × 9 generic directions at `L = 4` → exactly
  `{(ℓ, 0)}`, **5 of 25**; four RIGHT angles falsely admit `(4, ±4)` (7 of 25)
  — vv-principles #13 as measured; the second spelling through
  `S^2/SO2_x.quotient_map` agrees (`π(R d) == π(d)` bit-exact on every angle,
  the same 5 slots); at `L = 2` 3 of 9, at `L = 6` 7 of 49. ⭐ NEW —
  **the slot idiom is x-ONLY**: about y and z the mask returns only `(0,0)`
  and one `(1, ·)` slot, because the basis's polar axis is x
  (`_evaluate_real_sh`: `cos θ = μ_x`), while a rank test shows the
  `SO(2)_a`-invariant subspace is **1-dimensional per degree about EVERY
  axis** (dims `[1, 1, 1, 1, 1]` at `L = 4` for x, y and z). The trivial
  isotypic component exists about every axis and is slot-aligned only about
  x ⟹ an `SO2(a ≠ x)` UPSTAIRS face has no slot spelling and is refused
  naming the work (no shipped rule lives on `S^2/SO2_y` or `S^2/SO2_z`).
  ⭐ NEW — `_PARITY_PROBE_DIRECTIONS` (`spherical_harmonic_basis.py:410`)
  are **OFF the sphere**: norms `0.998, 0.831, 0.954, 0.958, 0.989`,
  `SPHERE.contains` **False** — so 0.6's refusal reaches the mirror probe;
  `[M]` normalising them leaves `even_slot_mask` bit-identical for axes
  0, 1, 2 × `L ∈ {1, 2, 3, 4, 5}`.
* **(d) existence** — `class Descent` / `LegendreBasis` / `ReynoldsProjection`
  / `OrbitAxis` / `IsometryGroup` / `Chart` all **0**, controls 1. Hashes
  37 / 36 / 1 dangling (`de29bcc6`); supersession 60 / 56; surprise log 72.
* **(e) §V.5h(a)'s G0 table** reproduces (`==` False / True / True;
  ambient-base all True; `quotient_group` `SO2('x')` / `Mirror('y')` /
  `None`). The Γ-slot: `[M]` `cylinder.admits_domain(folded_product(4,8))`
  **False** AND `admits_symmetry` **False** (the folded nodes are not
  `D_2h`-closed — σ_y maps a representative to its non-representative
  mate), `discrete_residual.contains(Mirror('y'))` **True**. The slot needs
  stage 0 AND stage 1; since 3.1 the stage-1 predicate is spellable — Γ acts
  on the orbit space through the entry's `quotient_map` (`π(g·nodes)` must
  permute `π(nodes)`).
* **(f) exit predicate 5** — ✅ **MEASURED, 13 of 13 rc=0**: `[M]`
  **10329 passed · 19 skipped · 69 xfailed · 0 failed = 10417 collected, the
  predicted denominator exactly, tree by tree** (numerics 2739 · transport
  645+1sk · geometry 727+4sk+1xf · data 237 · homogeneous 50 · diffusion 113
  · cp 141 · moc 121 · mc 39+2xf · cross_method 81 · sn 3384+1sk+50xf ·
  derivations 1637+13sk+11xf · root+harness 415+5xf; 61 min 34 s on a loaded
  machine). The table is in Part XIV's baseline section. Against CS4c step
  4's `10106 / 0 / 19 sk / 66 xf`: +223 passed, +3 xfailed (ERR-080's rows).

### (b) ⛔ THE RULING'S EVIDENCE WAS ONE LAYER TOO NARROW — the rectangular layout is the transport layer's angular MOMENT-SPACE CONTRACT

3.4-R was ruled on *"15 slot-indexing sites, ONE production"* (2026-08-31). The
explorer re-derived that predicate by AST: production **8** true slot-index
lines (7 are the PRODUCER, `_evaluate_real_sh` — the inherited 4 missed the
ℓ ≥ 2 loop; **1** consumer, `dsa.py:584`), tests **21 / 28 lines in 6 files**
(the inherited 15-in-5 counted a COMMENT, `test_keff_curvilinear.py:982`, and
missed two files), of which **2 lines in ONE test** (`test_dsa_rate.py:747-748`)
are on a 1-D rule. On the `.table` surface the 1-D blast radius is **one site**
(`dsa.py:585`; the 13 reads / 5 files reproduce — 2 layout-assuming reads are
2-D-windowed only, 5 are pass-throughs, 5 read `LossKernelBasis`'s own table).
The ruling's number was RIGHT for its predicate — and the predicate could not
see the layer the fix breaks:

`[M]` **8 production `SphericalHarmonicSpace.from_L(L)` mints** (AST, positive
control the scattering site): `transport/operators/scattering.py:207 / :378 /
:521`, `fission.py:510` (**L = 0**), `n2n.py:327` (**L = 0**),
`transport/fields/_bases.py:682` (`_space_for_mesh_and_L` — the moment-flux
field's space `SH(L) ⊗ cells`), `harmonic_moment_flux.py:334` (`truncate`), and
the basis's own `space` (`spherical_harmonic_basis.py:403`). `HarmonicMomentFlux`
reads `L` back off that factor (`find_factor(SphericalHarmonicSpace).L`; 4
test sites), `OperatorProduct`'s guard compares operator ends to
`frame.basis_space` by `(name, shape)` (`operator.py:388-401`, pinned vocabulary
`"requires equal domains"`), and `SphericalHarmonicSpace.__post_init__`
(`spaces/spherical_harmonic_space.py:153`) FORBIDS an `(L+1,)` shape. Two more
doors: `HarmonicFrame.__init__` (`transport/frames/harmonic_frame.py:331`) and
`from_galerkin` (`:368`) narrow to `SphericalHarmonicBasis` by `isinstance` —
`[M]` **0 mentions in this plan until now** — so a non-subclass `LegendreBasis`
is refused at the door of every 1-D solve, at **L = 0** (`fission.py:305` and
`n2n.py:200` call `for_space` at 0): the ISOTROPIC 1-D solve is inside the §6b
set. The moment CARRIER's own head reads are a fourth reader class (`[M]` 10
two-index constant subscripts: `values[0, 0]` in `isotropic_part` /
`scalar_flux`, `m[0, 0]` in `FissionMomentOperator.apply`, `ng = space.shape[2]`,
`truncate`'s `(L+1, 2L+1, …)`).

⟹ a flat `(L+1,)` basis mismatches every moment field and every operator end
minted from `L` alone. That is not a defect the fix introduces — it is the
Stage-2 generator ruling's own defect (*a Factory Object induces space and
operator TOGETHER at one site; forgetting = retaining the induced parts*): the
frame induces the moment space, and seven sites retain a copy minted from `L`.
Pattern 2, one tier above the table. The surprise is logged in
`plan-authoring`'s table (a shape minted independently by consumers, and a
type guard at a door — two more §6b spellings).

### (c) ✅ The four RULINGS (user, 2026-09-02, at this opener)

1. **3.4's realization — flat `{P_ℓ}`, AFTER a moment-space pre-step.** Ruled
   among (i) flat after a pre-step, (ii) a layout-keeping `LegendreBasis` (the
   fold's idiom: `P_ℓ` in the m=0 column, other columns structurally zero —
   every consumer flows unchanged; `[M]` the m=0 numbers are bit-identical to
   today's slot), (iii) flat with the transport carve inside the one fused
   commit. ⟹ **NEW tracker row 2.5** — *the angular moment space is READ off
   the frame, never minted from `L`* — its OWN commit, bit-identical (with the
   SH basis still everywhere, every space object is `(name, shape)`-identical
   and every number unchanged), landing BEFORE the fused commit. 3.4-R stands
   as ruled; its evidence row is corrected in place.
2. **The probe lives on the catalogue ENTRY** (`Quotient`): *which slots of a
   basis on my base descend to me* = fibre-constancy along `quotient_map`,
   with `SubgroupOfO3` supplying generic orbit samples (every element of a
   finite group; INCOMMENSURATE angles about the axis for `SO2`). `Descent`
   and `MirrorEvenSphericalHarmonicBasis` both READ it; the latter's private
   probe retires (Pattern 2). Ruled over (ii) on `Descent` (a basis → descent →
   basis dependence inside `numerics/basis`) and (iii) on the basis (the probe
   spelled once per group; a `LegendreBasis` carrying an `SO(2)` mask
   production never reads).
3. **G0 — ONE predicate, and `P_ℓ` on a full-sphere rule BUILT NOW.** The
   user's answer, verbatim: *"Can we make one and prepare now the Pl on a full
   sphere? This seems legitimate, so maybe we should handle it now."* Read as:
   (a) ONE predicate — the lattice one, of which object equality is the special
   case `K = H`: a frame binding `basis` to `measure` is admissible iff **the
   descent arrow `measure.support → basis.domain` EXISTS** (one base; the group
   the measure SPENT is contained in the group the basis HAS —
   `measure.quotient_group ⊆ basis.invariance_group`), and `FrameBase.table`
   tabulates the basis THROUGH that arrow (the pullback `f ↦ f ∘ π` that
   `ManifoldMap`'s docstring already names as `Descent`'s job); (b) the
   full-sphere pairing — `LegendreBasis(L, "x")` bound to a Lebedev /
   level-symmetric rule (`Trivial ⊆ SO2('x')`), tabulated as `P_ℓ(Ω·ê_x)` —
   is BUILT and WITNESSED in the fused commit (its Gram against a 3-D rule is
   `4π/(2ℓ+1)` diagonal to the rule's exactness; an isotropic field analyses to
   `φ_{ℓ ≥ 1} = 0`). Home: `FrameBase`, every frame. ⚠ **If this reading is
   wrong, correct it HERE before the fused commit lands.**
4. **The Γ-slot (2.2's second half) is its OWN commit, right AFTER the fused
   one** — stage 0 admits a rule on `S²/H` for `H ⊆ Γ` (lattice, DERIVED, no
   declared column); stage 1 asks Γ-invariance through the entry's
   `quotient_map`. Its §6c witness exists alone (the folded cylinder, refused
   today → admitted). Tracker row **2.2b**.

⭐ **Three more, taken after 2.5 landed (user, 2026-09-02, before a line of
the fused commit):**

5. **The G0 reading above is CONFIRMED** — *"Yes — build it as recorded"*:
   one lattice predicate at `FrameBase`, the table pulled back along the
   descent arrow, the Legendre-on-a-full-sphere pairing BUILT and witnessed
   on level-symmetric / Lebedev rules, with `product(4,4)` and
   `folded_product(2,4)` at L = 4 as LABELLED aliasing controls (`[M]` the
   memo's J3: 3.7 / 4.9 — the two-pole convention, not a defect).
6. **The O(2)_a refusal (memo H-6) is NOT the fused commit's.** The derived
   lower bound stands (`SO2('x')`), the lattice G0 refuses a Legendre basis
   on a σ_y-folded rule (no shipped consumer), and *"schedule this
   development for the next step after a context compaction. It's another
   commit."* ⟹ NEW tracker row **1.9** — *a basis states the FULL group its
   functions have: `SubgroupOfO3.O2(axis)`* — first step AFTER the
   compaction that follows the fused commit; GitHub issue filed.
7. **A Legendre frame whose order out-resolves its 1-D rule keeps the
   `DenseMetric` arm** (status quo — every over-resolved SH frame gets it
   today); the GL dead-slot theorem (memo H-5: `GL_n` is diagonal-and-exact
   for L ≤ n−1 and dead at ℓ = n) is PINNED as a gate and is the closed-form
   verdict exit predicate 3 owes the coarse-rule rows.

### (d) ⟹ the honest shape now: THREE landings, in this order

**2.5 — the pre-step (bit-identical, its own gate):** `HarmonicFrame`'s two
doors accept any `Basis`; the seven production `from_L` mints read the frame's
basis space (`quadrature.angular_frame(L).basis.space` / `frame.basis_space`);
the moment-flux field gains a frame-derived constructor (`from_mesh_and_L`
stays as the full-SH family's honest `L`-only mint for the `[M]` 51 + 22 test
constructions that model no quadrature — priced by the test-architect's memo);
the carrier's head reads become basis-mediated (the basis names its isotropic
slot). Gate: a mutation of the frame basis's `space` must MOVE the field's
space and the operator ends (today it does not — that IS the twin); the two
canaries the explorer named (`test_harmonic_frame.py:162`,
`test_harmonic_moment_flux.py:308`) re-described for what they still separate
(coding-standards: single-sourcing a duplicate demotes every gate that compared
its copies — find the EXTERNAL pin or write one).

**The fused commit (0.1b + 0.6 + 2.2-G0 + 3.4 + 3.4b):** §V.5k(g), plus — the
descent arrow at `FrameBase` with the table pulled back along it;
`LegendreBasis.evaluate` accepting BOTH honest coordinate systems of the orbit
space (the realization's `(N,)` μ — the slab rule; the base's `(N, 3)` section
coordinates — a full-sphere rule, through `quotient_map`);
`Quadrature._harmonic_basis` dispatching on the MEASURE'S SUPPORT (`match`:
`Quotient(base=Sphere(), by=SO2(a))` → Legendre, `Quotient(by=Mirror)` →
MirrorEven, `Sphere()` → SH; never on `folded_by`); the Legendre coefficient
space a NEW class (`SphericalHarmonicSpace` forbids its shape);
`_PARITY_PROBE_DIRECTIONS` normalised; the DSA μ-read → `axis_cosines(0)`
(`[M]` explorer P7: bit-equal on 1-D rules — and the read is UNCONDITIONAL in
`L`, so every slab DSA solve is in the set); the DENSE-Gram witness of
campaign-1 P7 RE-POSED (`[M]` slab GL8 L=2 is its flagship — `test_frame.py`
D-slab / D3 / D4, `test_dense_metric` A2 — and a Legendre Gram on GL8 is
diagonal, so D3 needs a re-pose and A2 goes silently inert);
`test_binding_tightness._equispaced_frame` re-pinned (`[M]` a raw
`Quadrature(measure=… support=COSINE_INTERVAL)` outside the registry — G0
refuses it); the L = 1 solve-level pins (the two 1-D P1 regression snapshots,
`test_affine_carve_bit_identity`'s slab case) measured against the saved
baseline — bit-identical if the reconstruction's contraction order permits,
otherwise the three-criteria re-baseline, SAID rather than tolerated. Then exit
predicates 1–4.

**2.2b — the Γ-slot:** separate, after.

⭐ The **test-architect** was dispatched at this opener
(`scratch/_fused_test_architect_brief.md` → `scratch/_fused_verification_plan.md`)
per the proactive trigger (a carve crossing numerics → transport → sn); its
memo is read before a line of 2.5 is written.

---

## §V.5k ⭐⭐ The FUSED STEP's PRE-FLIGHT (2026-09-02, at `2f294ef1`) — the fusion holds, THE FIX's blast radius is the frame's TABLE LAYOUT, and three rulings are owed at its opener

The **19th** opener's raw material, measured at the compaction point after 3.1
landed. Every `[M]` below is at `2f294ef1` (`scratch/_p31c_reground.py` and
two AST censuses beside it; positive controls in the same pass — `class
Quotient` 1, `class ManifoldMap` 1). `plan-authoring` §1: everything under (g)
is a HYPOTHESIS until the fused step's own opener re-measures (a)–(f).

### (a) ✅ XII.1b's fusion premise re-confirmed by direct reading

`_harmonic_frame_measure` (`directional.py:730`) still forges the 1-D arm as
a RAW `DiscreteMeasure(…, support=SPHERE)` at `:794` — **1 of 4**
`support=SPHERE` keyword literals in `orpheus/` by AST (the other three are
honest: `exactness.py:237` is `UNIFORM_ON_SPHERE`, `rules_sphere.py:127`
and `:1113` are Lebedev and level-symmetric, 3-D rules on the sphere; ⚠ a
REGEX census reads **7** because three docstrings spell the kwarg — state
the predicate). The arm has exactly **1** caller (`angular_frame`, `:678`,
through `_harmonic_basis` + `GalerkinFrame` at `:726`). `test_q8_4`
(`tests/numerics/test_quadrature_directional.py:568`, parametrized over
`_LIFTED_RULES`) is the self-retiring gate that goes RED when the arm dies.
The ERR-080 gate declares **3** `xfail(strict=True)` rows (`:228`, `:255`,
`:277`). `_evaluate_real_sh` has **4** call sites, all inside
`spherical_harmonic_basis.py` (`:208`, `:236`, `:509`, `:512`) — 0.6's raise
is local, and it would FIRE on the forged frame's `(μ, 0, 0)` nodes
(`‖·‖ = |μ| < 1` off the poles), which is exactly why 0.6 rides 3.4. 2.2:
`select_quadrature` has **0** production callers (26 test sites; §II.7
counted 21) and `solve_sn(materials, mesh, quadrature, …)` takes its rule as
a required positional — so G0's home is the FRAME, and §V.5h(c) stands: the
frame's measure still says the forged `S^2` for the slab, so a G0 written on
`(frame.measure, basis)` is INERT until 0.1b and one written on
`(quadrature.measure, basis)` reddens every slab solve today. **Land the
four together, as ruled.**

### (b) ⛔ THE FIX's blast radius is the frame's TABLE LAYOUT — and the ruling's numbers are INHERITED

3.4-R ruled (b), a true `LegendreBasis`: `{P_ℓ}` on `S²/SO(2)_a`, `L+1`
members. So on a 1-D rule the rectangular `(L+1, 2L+1)` SH table that every
consumer indexes is **not kept** — unlike `MirrorEvenSphericalHarmonicBasis`,
which keeps the layout and structurally zeroes the σ-odd columns so its
consumers flow through unchanged. XII.1b's §6b rule therefore applies with a
BASIS CHANGE in place of a signature change, and its set is *every reader of
the table's layout on a 1-D rule*. `[M]` the raw material at this HEAD:

| what | `[M]` |
|---|---|
| `angular_frame(` production call sites | **3** — `directional.py:676` (the `spherical_harmonics` wrapper; 0 production consumers of its own), `sn/acceleration/dsa.py:585` (reads `angular_frame(1).table[:, 1, 1]`, μ off the ℓ=1 Cartesian slot — **DSA runs on slabs**), `transport/frames/harmonic_frame.py:416` (`HarmonicFrame.from_galerkin`) |
| `HarmonicFrame.moment_space_on(` | **2** production readers in `harmonic_frame.py` — derives the moment codomain from the basis's shape; consumed by `for_space` at **4** production sites (`transport/operators/{fission,n2n,scattering}.py`) |
| `.table` attribute reads in `orpheus/` | **13**, in 5 files: `numerics/frame.py`, `quadrature/directional.py`, `sn/acceleration/dsa.py`, `sn/loss_representation/__init__.py`, `sn/operators/loss_kernel_gauge.py` |
| explicit `[:, ℓ, m]` slot indexing on a table/`Y`/frame line | **4** production lines, all the PRODUCER (`spherical_harmonic_basis.py:578-584`) + **15** test lines in 5 files (`test_quadrature_directional`, `test_spherical_harmonic_basis`, `test_dsa_rate`, `test_keff_curvilinear`, `test_solver_components`) |
| `SphericalHarmonicBasis(` constructors | 3 production (`directional.py:611`, `spaces/spherical_harmonic_space.py:203`, `:220`) + 26 tests; `MirrorEven…(` 1 production (`directional.py:621`) + 6 tests; `GalerkinFrame(` 3 production (`frame.py:760`, `directional.py:726`, `loss_kernel_gauge.py:1199`) + 17 tests |
| `Basis` subclasses at runtime | 4 direct / 6 recursive, unchanged; `LegendreBasis` would be the 7th |

⚠ **The 3.4-R evidence table's *"15 sites, ONE production"* and §4.9's *"12
collateral sites"* are 2026-08-31 numbers under two different predicates**
(explicit slot indexing vs. *consumes the slab-GL-L2 frame as a fixture*);
the rows above are raw material, not the set. Re-census both at the opener,
by AST, and remember `plan-authoring` §6b's inventory of members spelled
WITHOUT the symbol: the §4.9 pins that assert the fabricated ℓ=2 Gram diagonal
`[0.4, 0.8, 0.8]` (`test_frame.py:906`) go RED by design — **re-pin, never
delete** (`coding-standards`: retirement means test migration).

⛔ **RE-CENSUSED at the opener, 2026-09-02 (§V.5l(b)):** the slot-index predicate reproduces at 8 production lines (7 producer + 1 consumer) and 21/28 test lines in 6 files — and it is the WRONG LAYER. The fix's blast radius is the angular MOMENT-SPACE contract: 8 production `SphericalHarmonicSpace.from_L` mints, two SH `isinstance` doors on `HarmonicFrame`, and the `(name, shape)` guards between them and the frame. A pre-step (tracker **2.5**) was minted for it, user-ruled.

### (c) The isotypic probe — §III.8's specification stands, its idiom ships, and 3.1 gave it a SECOND spelling

`[M]` `MirrorEvenSphericalHarmonicBasis.even_slot_mask`
(`spherical_harmonic_basis.py:498`, a `cached_property`; 2 production
readers, both in-file) evaluates the parent at `_PARITY_PROBE_DIRECTIONS`
(`:410`, 5 generic directions) and at their σ-images and classifies each
slot even/odd by `isclose(atol=1e-12)`, raising on an unclassifiable slot.
The `SO(2)` version rotates the probe directions about axis `a` by
INCOMMENSURATE angles and keeps the slots invariant under ALL of them —
§III.8 measured 6 angles × 9 directions at `L = 4` → exactly `{(ℓ, 0)}`,
5 of 25, zero false positives or negatives (`vv-principles` #13: right
angles generate `C₄`). ⚠ **No re-runnable reproducer of that measurement
survives** (`grep -l isotypic scratch/*.py` → only the σ_y probe); the opener
re-measures it. ⭐ And since 3.1 the probe has a second, more honest
spelling: a slot descends to `S²/SO(2)_a` iff its function is constant on
the fibres of `Quotient.quotient_map` — the pullback along π — which is
ALSO the isomorphism check `Descent` must carry. **RULING owed (§V.5j(d)):
the probe's HOME** — a method of the basis (the `MirrorEven` idiom, one per
group) or of the catalogue ENTRY (*which slots of a basis on my base descend
to me*, where `Descent` reads it). Two homes, one ruling; do not build it
twice.

### (d) 3.4b `Descent` — ruled (i) mint, *with or before 3.4*; since 3.1 it carries the ENTRY, not a pair

Part V.5 §C: `Descent` carries `(M, H, M/H)`, the two realizations (the
invariant subspace upstairs, the quotient's own basis downstairs), the
checkable isomorphism (analyse-then-descend == descend-then-analyse to
machine precision) and the DISCRIMINATOR written where both classes see it
(*downstairs when the quotient has a classical named basis*). ⭐ After 3.1
the pair and the arrow are ONE object — the `Quotient` entry with its
`quotient_map` and its `reference` — so `Descent` reads the entry rather
than carrying `(M, H)` again (Pattern 2). `[M]` existence at `2f294ef1`:
`class Descent` **0**, `class LegendreBasis` **0**, `class
ReynoldsProjection` **0**, `class OrbitAxis` **0**. ⛔ Do NOT name the
sub-basis `ZonalSphericalHarmonicBasis` (Part V.5's ban: "zonal" is the 1-D
case wearing the name of nothing).

### (e) 2.2 owes TWO rulings, not one

§V.5h(a): G0 as object equality (refuses the Part I bug, admits both shipped
angular pairings, refuses the unshipped `P_ℓ`-on-a-full-sphere pairing until
someone needs it) **or** ambient-base equality + lattice G2. And §II.10's
Γ-slot: `AngularSymmetry.support` derives from `continuous_isotropy` ALONE,
so stage 0 refuses the shipped CYLINDER (`folded_product` on `S²/σ_y`;
`[M]` §II.10: 5 production builds rejected) — 2.2 owns the discrete slot the
ontology lacks. Both rulings at the opener, before a line is written.

### (f) Exit predicate 5 — the 13-tree baseline is being MEASURED at this compaction

Launched detached at `2f294ef1` (`scratch/_p31c_full_gate.sh` →
`scratch/_p31c_full_gate.log` + `_p31c_gate_<tree>.log`), serial canonical
runner, with the per-tree denominators PREDICTED first by collect-only:
numerics **2739** · transport 646 · geometry 732 · data 237 (+2 desel) ·
homogeneous 50 · diffusion 113 · cp 141 (+13) · moc 121 (+3) · mc 41 (+16)
· cross_method 81 (+8) · sn 3435 (+116) · derivations 1661 (+67) ·
root+harness 420 (+2) — **10417 collected** under `-m "not slow"` (⚠ the
matrix's 10644 is a DIFFERENT population: every collected test, slow tier
included — §2's filter clause). The measured rows land in Part XIV's
baseline section, in a follow-up commit if the run outlives this
compaction. Predicates 1–4 and 6 have their reproducers:
`scratch/_pl_slab_defect_repro.py` LEG A, the ERR-080 gate's XPASS,
`_p01a_ground.py`'s 12-rule route table, and `sphinx -W` +
`dead_references`.

### (g) ⟹ the fused step's honest shape, as a hypothesis for its opener

ONE commit: `LegendreBasis` (domain the `S²/SO(2)_a` entry — so
`invariance_group` is `SO2(a)` for free through `domain.by`, 2.1b's
derivation; its slot set DERIVED by the probe, never hand-listed);
`Descent` reading the entry; the 1-D arm of `_harmonic_frame_measure`
DELETED (0.1b) with `test_q8_4` rewritten as the ROUTE gate for 1-D rules;
`_evaluate_real_sh` refusing non-unit directions (0.6); G0 and the Γ-slot
armed at the frame (2.2); the DSA μ-read re-pointed off the ℓ=1 Cartesian
slot; `HarmonicFrame.moment_space_on` deriving the moment codomain from the
basis's OWN shape; every §4.9 pin re-keyed. Then exit predicates 1–4 flip:
LEG A reads `+4.000000000000` at `L = 0..3`, the strict-xfail rows XPASS and
their markers are deleted, the isotropic-flux `ℓ ≥ 1` census reads `0`, and
every 3-D rule is bit-identical at every `L`.

---

## §V.5j ⭐⭐ 3.1's PRE-FLIGHT (2026-09-02, at `14c37aa9`) — the catalogue HOME is already landed, the row's re-home is not on the path, and what 3.1 still owes is TWO fields and one map

The **17th** opener, run at the compaction point after 2.3 landed. Every `[M]`
below is at `14c37aa9` (`scratch/_p23_reground.py`, AST/regex with positive
controls; the existence pass is the same instrument Part XIV's §1 uses).

### (a) ⛔ The row's premise is pre-2.0a: the "`numerics/symmetry/` catalog" LIVES IN `manifold.py`, ruled at E. Placement

3.1's row reads *"`numerics/symmetry/` catalog: `IsometryGroup` entries from
1.1 … `[M]` `numerics/symmetry.py` is a module today — this is a re-home"*.
`[M]` the catalogue is `manifold._ORBIT_CATALOGUE` (`manifold.py:1270`): six
keys, two procedures (`_sphere_mod_so2` `:1036`, `_sphere_mod_mirror` `:1132`)
plus the derived trivial quotient (`_mod_trivial` `:1229`), memoised through
`_catalogued_quotient` `:1005` — placed there by E. Placement (*"`Manifold` …
its own module with no intra-`numerics` imports"*) and consumed in production
since 2.4. `class IsometryGroup` **0** — the tree's group type is
`SubgroupOfO3` (`symmetry.py:414`, 7 classes in the module) and minting a
second would be the twin the plan forbids. ⟹ **the "catalog home" half of 3.1
is LANDED** (2.0a + 1.1 + 2.4), and XII.2's *"the catalog home + the probe"*
is half discharged. The re-home of `symmetry.py` (1847 lines, 7 classes, 29
functions, **30** importer files — 10 `orpheus/`, 20 `tests/`) into a package
serves Phases 5–6 (the other entries), which XII.2 says are **NOT on the
path**; and ⚠ §6d: `numerics/__init__.py:43` re-exports `SubgroupOfO3` and
`tests/test_layer_imports.py:172-190` documents the `numerics → geometry`
back-edge that runs THROUGH `symmetry` (it imports `geometry.transformation`),
tolerated only because geometry imports numerics by SUBMODULE. A package
re-home must keep that discipline and is its own §6d step — not 3.1's.

### (b) The data model, measured against D0.1's field list — 5 of 7 derivation outputs ship; the two missing are the ones 2.3 named

D0.1: an entry's fields ARE the procedure's outputs — *generators, syzygy
ideal, `P`, `det P`, the chart, the pushforward measure, the stratum*. `[M]`
`dataclasses.fields(Quotient)` = `base, by, realization, fundamental_domain,
generators, syzygy, gram, det_gram, derived_by, singular_stratum` (10).
Present: generators ✓ syzygy ✓ `P` ✓ (`gram`) `det P` ✓ stratum ✓ (+
`derived_by`, the provenance D0.1 also asks for, and the symbolic regression
tests in `test_manifold.py::TestQuotient`). **Absent: the chart, as a MAP —
and the pushforward measure.** `realization` is the chart's CODOMAIN; the
map `S² → S²/SO(2)_a`, `Ω ↦ Ω·ê_a`, *"is still not a value anywhere"*
(`manifolds.rst`, the arrows-not-built section). 2.3 minted the TYPE it
needs (`ManifoldMap`) and its witness already exists: the archivist measured
`π ∘ archimedes_a = pr₁` bit-exactly. The row's other items — *the action on
the axis index set, Haar (or counting) measure* — have **0** consumers on the
exit path (3.4 reads the mask, `Descent` reads the quotient map).

### (c) The pushforward REFERENCE measure — the residue 2.3 deferred, its §6b/§6c/§6d measured

*Twin:* `AngularSymmetry.reference` (`registry.py:911`) answers `LEGENDRE`
for any rotation axis, `UNIFORM_ON_SPHERE` for `Trivial`, raises otherwise.
`[M]` its ONE production reader is the selector's V stage,
`registry.py:1272` (`wanted = angular_symmetry.reference`). ⚠ Landmine: a
grep for `.reference` returns **46** hits because `ExactnessClaim.reference`
shares the name (`exactness.py:343/391/399`) — discriminate by the receiver.
*§6b for a new `Quotient` field with a default:* `[M]` **3** constructor
sites, all in `manifold.py` (the three derivation returns `:1106`, `:1201`,
`:1248`) + 1 in `tests/`; the other `Quotient(` hits (`measure.py:489`,
`manifold.py:1319`) are `match` patterns, untouched by a keyword field.
*§6d:* module scope is a **5 of 5** cycle (measured at 2.3); the viable
mechanism is a function-scope import inside the derivation function — ⚠ and
that is safe only if no quotient is built at IMPORT time. `[M]` today
`AngularSymmetry.support`/`.reference` are properties, `gauss_legendre_on_polar_orbit`
is a function, and the registry's static table stores specs, not measures —
so the first `quotient()` runs at rule construction, after every module is
loaded. **Inject and run at 3.1's opener anyway** (the 2026-08-26 §6d row:
order-dependent cycles pass smoke tests in the reassuring direction).
*§6c:* the σ_y entry's pushforward is the WEIGHTED disk measure
`2 dx dz / √(1 − x² − z²)`, which no shipped `ReferenceMeasure` realization
can spell (`UniformMeasure`, `GeneratingMeasure` — 1-D recurrences —,
`ProductMeasure`); it ships `None`, and honestly: `[M]` the fold drops
`exactness`, so nothing reads it. The `S²/Trivial` entry's is
`UNIFORM_ON_SPHERE` — but ⚠ `AngularSymmetry.support` returns the bare
`SPHERE` for `Trivial`, not the trivial quotient, so the registry's `Trivial`
arm cannot read it off the entry without a design choice (route `support`
through `_mod_trivial`, or keep the arm). Rule it at the opener.

### (d) "The probe" — it is 3.4's, and its idiom ships for a MIRROR

XII.2's *"the catalog home + the probe"* names §III.8's isotypic probe:
invariance of each real-SH slot under `SO(2)`-about-`x` over incommensurate
angles × random directions (`[M]` §III.8: exactly `{(ℓ,0)}`, 5 of 25 slots at
`L = 4`). The idiom ships as `MirrorEvenSphericalHarmonicBasis.even_slot_mask`
(`spherical_harmonic_basis.py:498`, `_PARITY_PROBE_DIRECTIONS`) for a group
of order 2; the `SO(2)` version must sample INCOMMENSURATE angles
(`vv-principles` #13 — right angles generate `C₄`). It is the input 3.4's
`LegendreBasis` mask is derived from, so it lands with 3.4 (fused step),
not as 3.1's deliverable — unless 3.1's opener rules the probe a method of
the catalogue ENTRY (*"which slots of a basis on my base descend to me"*),
which is where a `Descent` (3.4b) would read it from. Two homes, one
ruling; do not build it twice.

### (e) ⟹ the honest scope of 3.1, as a hypothesis for its opener to rule

**3.1 = the entry carries its two remaining derivation outputs**: (i) `chart:
ManifoldMap` — `S² → S²/SO(2)_a`, `Ω ↦ Ω·ê_a`, codomain the entry itself,
witnessed by `π ∘ archimedes_a = pr₁` and by `barycentre ∘ π` being the
axial projection; (ii) `reference: ReferenceMeasure | None` — `LEGENDRE`
for the three axial entries, populated inside `_sphere_mod_so2` by
function-scope import; `AngularSymmetry.reference` reads it. Small — three
constructor sites, one twin absorbed, one new field pair. The package
re-home and the other entries stay off the path (XII.2). §6c for (i):
nothing shipped is refused by a chart field; its witnesses are derivation
agreements, like 2.3's. `plan-authoring` §1: *this scope is a hypothesis
until the opener re-measures (a)–(d) at its own HEAD.*

✅ **CONFIRMED at the opener and LANDED 2026-09-02** — (a)–(e) reproduced
**5 of 5** at `3623adc2` (the one new census fact: a THIRD `match
Quotient(...)` pattern at `basis/base.py:361`, a non-member for a keyword
field). Three user rulings refined (e) without refuting it: the map is
`quotient_map` with codomain the ENTRY (not *chart*), carried as a stored
`orbit_coordinates` + a derived property (a frozen dataclass cannot hold a
map whose codomain is itself), and the trivial entry's reference is `None`
(not `UNIFORM_ON_SPHERE` — Lebesgue on the BASE is the base's to know).
See the *3.1 EXECUTED* section.

---

## §V.5i ⭐⭐ 2.3's PRE-FLIGHT (2026-09-02, at `4d946988`) — the SECTION has no consumer, the tree already spells THREE maps, and one of them TWICE

The **15th** opener, run at the compaction point after 2.1b landed so the next
session designs 2.3 against the tree, not against Part XIV's pre-2.4 pointer.
Every `[M]` below is at `4d946988` (the 2.1b tree plus two docstring lines).

### (a) ✅ §V.5h(e)'s hypothesis CONFIRMED by census — nobody reads a section

`[M]` `Quotient.fundamental_domain` (the section's IMAGE, in the base's
coordinates) has **0** production readers outside `manifold.py` — 7 sites,
all in `tests/numerics/test_manifold.py`; `Quotient.realization` has **1**
(`DiscreteMeasure.on_orbit_space`'s equality guard, `measure.py:1027`).
The only section-LIKE production consumer is the forgery arm
(`directional.py:786`), which 3.4 retires. ⟹ **2.3 is the QUOTIENT MAP
and the ARCHIMEDES CHART.** `_sphere_mod_so2`'s `fundamental_domain=None`
stays `None` on purpose, and Part XIV's pre-2.4 pointer (*"the SECTION MAP's
natural home"*) is ⛔ REFUTED as a scope for 2.3 (kept in place per §3).

### (b) ⭐ THREE maps around the quotient — two are spelled already, one of them TWICE

| map | where the tree spells it today | consumer after 3.4 |
|---|---|---|
| the quotient map `π: S² → S²/H` | `DiscreteMeasure.quotient` (`measure.py:1161`): the orbit-representative retraction `nodes ↦ nodes[representative]` pushed with `new_space=support.quotient(group)` — π as a node permutation plus a relabel | ⭐ `Descent` (3.4b) pulls `P_ℓ` back along it |
| the Archimedes chart `[-1,1] × S¹ → S²` — ⭐ already a LABELLED equation, `:eq:product-mu-phi-cosines` | `spherical_product` (`rules_product.py:527-531`): `column_stack([mu_x, mu_y, mu_z])` + `support=SPHERE` — L7's literal | `product()` and every folded product |
| ⭐ the BARYCENTRE map `S²/SO(2)_a → B³`, `μ ↦ μ·ê_a` | **TWICE.** Honestly by `_embedded_nodes` (`symmetry.py:951`, axis-aware since 2.4; **2** production callers — `is_invariant`'s dispatch `:1052` and the motion-preservation check `directional.py:536`), where an image OFF `S²` is CORRECT: the barycentre of an `SO(2)`-orbit lies on the axis inside the ball, and a closure test wants exactly that point. Dishonestly by the forgery arm (`directional.py:786`): the SAME `(μ, 0, 0)` padding declared `support=SPHERE` — its own comment says *"the barycentre is off S²"* | `_embedded_nodes` keeps it; the forgery retires at 3.4 |

⟹ **ERR-080 restated in 2.3's vocabulary: the barycentre chart with a
FORGED codomain.** `[M]` `Ball` (`manifold.py:454`, minted at 1.1) has **0**
production consumers — the barycentre map's honest codomain exists and
nothing names it. ⛔ **FILTER defect, caught by the archivist at 2.3's
landing (2026-09-02):** that census ran `grep -v manifold.py`, so it could
not see `Ball(2)` at `manifold.py:922` — the σ_y entry's own `realization`.
The honest claim is *0 consumers OUTSIDE its module, and `Ball(3)` never
constructed*; `barycentre` is the first arrow whose CODOMAIN is a ball.
(§2's FILTER clause: excluding the defining module from a consumer census
hides every in-module consumer, and reads as "nothing names it".) If `pushforward` derives `support = chart.codomain`, a
measure pushed along the barycentre chart lands on `Ball(3)` and the forgery
is unspellable THROUGH THAT VERB. ⚠ Not through the raw constructor: the
forgery arm calls `DiscreteMeasure(...)` directly, so 0.1b's refusal at
construction (`Sphere.contains`) is what closes that door — complementary,
not redundant — and **2.3 alone must not be expected to red the ERR-080
gate** (§6c: said here so the fused step is not credited to 2.3, or 2.3 to it).

### (c) ⛔ The tracker's done-when is §10's third shape — scope it to MAP-BUILT measures

The row's done-when reads *"`grep` finds no `SPACE_SPHERE` literal at a
measure constructor"*. `SPACE_SPHERE` retired at 2.0c; the literal is now
`support=SPHERE`, and `[M]` there are **5** in production: `exactness.py:238`
(the reference `UniformMeasure` — the DEFINITION of the sphere's measure, not
a rule), `rules_sphere.py:130` (Lebedev, tabulated ON `S²`),
`rules_sphere.py:1116` (level-symmetric — ⚠ the SAME two lines as the
product, `column_stack([mu_x, mu_y, mu_z])` at `:1112` + the literal, but
its cosines come from `_build_level_symmetric_arrays`, tabulated on `S²`
natively — NOT a map from a product chart), `rules_product.py:531`
(the product — map-built) and `directional.py:788` (the forgery —
map-built, retires at 3.4). **Only 2 of 5 are map-built.** The honest tell:
*no measure BUILT BY A MAP names its support by literal* —
`product(4, 8).measure.support` derived from the chart's codomain, nodes
bit-identical. A tell over all five is pinned red forever by three honest
tabulations the step must not touch.

### (d) `pushforward`'s §6b set — and a present-tense-false docstring beside it, fixed

`pushforward(phi, *, new_space: Manifold)` (`measure.py:836`) — `[M]` **1**
production caller (`:1161`, inside `quotient()`) and **8** sites in
`tests/numerics/test_measure.py` (`:261, :272, :303, :322, :344` — the
no-`new_space` negative leg — `:347, :716, :792`). If the signature takes a
chart, that is the set, 9, by grep on `.pushforward(`; the string form
`"pushforward"` returns nothing else (§6b's spelled-without-the-symbol
inventory, checked). ⛔ Found on the way and FIXED at this compaction
(`4d946988`): `measure.py`'s module and class docstrings still spelled
`μ.pushforward(φ)` with no codomain — the pre-2.0c call, present-tense-false
since `025834f5` retired the fabricated default. The method's own docstring
was re-tensed then; the two overviews 800 lines above it were not
(`coding-standards`' tense clause, the far-from-the-diff variant).

### (e) The pushforward MEASURE is a twin on the registry — and attaching it to the seed has a §6d hazard, measured

`[M]` `AngularSymmetry.reference` (`registry.py:911`) answers `LEGENDRE`
for any `rotation_axis` — with the Archimedes argument written in a COMMENT
(*"the pushforward of dΩ under μ = Ω·ê_a is 2π dμ"*) — and
`UNIFORM_ON_SPHERE` for `Trivial`, raising otherwise. That is D0.1's 8th
field (the orbit space's pushforward measure) living on the registry;
`Quotient` has **10** fields and none is a measure. ⛔ §6d, before anyone
designs it: `exactness.py:115-116` imports `manifold` at module scope and
`manifold.py` imports nothing from `exactness`, so a `ReferenceMeasure` field
on the catalogue entry populated at module scope creates
`manifold → exactness → manifold`. Inject-and-run at 2.3's opener before
choosing between a lazy field, a key the registry resolves, or the reverse
dependency.

### (f) §6c for 2.3, answered before it is asked

With `support = chart.codomain` derived, **no shipped input is refused**: the
fold's chart codomain is the quotient it already names; the product's is
`SPHERE`. A domain guard (`chart.domain == measure.support`) refuses only
test-built mismatches. ⟹ 2.3 is an enabler like 2.1b, and its honest
witnesses are derivation AGREEMENTS (`product(4, 8).measure` bit-identical
with its support derived; the fold's derived support `is` the memoised
quotient), refusal legs constructed in tests, and the Pattern-2 collapse of
(b)'s twice-spelled barycentre map onto ONE chart object. `[M]`
`DiscreteMeasure.__mul__` (`measure.py:608`) is the tensor product and
`COSINE_INTERVAL * CIRCLE` is a `Product` manifold, so
`spherical_product(polar, azimuthal)` is `(polar * azimuthal).pushforward(ARCHIMEDES)`
in everything but the chart object.

### (g) Order

✅ **RULED 2026-09-02 at 2.3's opener** — (a)–(f) reproduced 6 of 6 at `f2682869`; the import injection: **5 of 5** orders die. Rulings: ONE frozen value type named **`ManifoldMap`** (not `Chart`); ALL THREE maps typed in 2.3; the pushforward measure DEFERRED to 3.1 (populated inside the derivation function, never at module scope). See the *2.3 EXECUTED* section.

*The question as asked:* 2.3's opener re-runs (a)–(f) at its own HEAD (all greps and one import
injection — minutes) and then RULES the chart's shape with the user: a frozen
map type with `domain`, `codomain` and `__call__`; whether `quotient()`'s
retraction and `_embedded_nodes`'s barycentre become charts in the SAME step
or only the Archimedes one; and where the pushforward measure lives (e).
Then 3.1 scoped to `SO(2)`, then the fused step, where 2.2 owes its G0 ruling
(§V.5h(a)).

---

## §V.5h ⭐⭐ 2.1b + 2.3's PRE-FLIGHT (2026-09-01, at `8fc781d2`) — G0 is now a MANIFOLD question, and half of 2.3 may have no consumer

The **13th** opener, run at the compaction point so the next session designs
against the tree rather than against Part IV's 2026-08-31 vocabulary. Every
`[M]` below is at `8fc781d2`.

### (a) ⛔ Part IV's *"G0 cannot be a tag-equality test"* is VOID — both of its premises moved

Part IV argued: *the slab measure honestly says `[-1,1]`; a trivial-isotypic
basis honestly lives on `S²/SO(2)`; those name the same manifold; string
equality would refuse a legal pairing.* Since 2.0c the vocabulary is
`Manifold` **object** equality, and since 2.4 the slab measure says
`S^2/SO2_x`. Measured on the three pairings the tree ships, with two candidate
G0 predicates side by side:

| pairing | `measure.support` | `basis.domain` | G0 as `==` | G0 as ambient-base `==` | `measure.quotient_group` |
|---|---|---|---|---|---|
| slab rule vs `SphericalHarmonicBasis(L=2)` | `S^2/SO2_x` | `S^2` | **False** | True | `SO2('x')` |
| `folded_product(4,8)` vs its own frame basis (`MirrorEven…`) | `S^2/sigma_y` | `S^2/sigma_y` | **True** | True | `Mirror('y')` |
| `lebedev(17)` vs `SphericalHarmonicBasis(L=2)` | `S^2` | `S^2` | **True** | True | `None` |
| energy measure vs `SphericalHarmonicBasis(L=2)` | `energy` | `S^2` | False | False | — |

⟹ **An equality G0 refuses exactly the Part I bug and admits both shipped
angular pairings** — it is the `test_d6` keystone's own claim (*a frame's two
halves name ONE manifold*), now spellable for the slab. The ONE legal pairing
equality would refuse and the ambient-base + lattice form would admit is
row 3 of Part IV's table — `P_ℓ(μ)` evaluated on a full sphere rule — and
`[M]` it ships **nowhere** (`class LegendreBasis` is 0; no basis today has a
quotient domain except the fold's, which pairs with its own quotient rule).
⟹ **2.2's opener owes a RULING, not a measurement**: G0 = object equality
(stronger, one predicate, refuses the unshipped legal pairing until someone
needs it) **or** G0 = ambient-base equality + G2 = lattice. Part IV's
paragraph is annotated in place; it must not be built to as written.

### (b) ✅ G2's MEASURE side became real at 2.4; the BASIS side is all of 2.1b

`[M]` `gauss_legendre(8).measure.quotient_group == SO2('x')` (was `None`
before 2.4) and `folded_product(4,8).measure.quotient_group == Mirror('y')`.
The lattice answers Part IV's four rows on the measure side today:
`Trivial.contains(SO2('x'))` **False** (refuses the bug once a full-SH basis
says `Trivial`), `SO2('x').contains(SO2('x'))` True (the repair),
`SO2('x').contains(Trivial)` True (a smaller space on a full rule),
`Mirror('y').contains(Mirror('y'))` True (the fold). `[M]` by runtime
`__subclasses__`: **0 of 6** bases answer `invariance_group`. Derivable:
`SphericalHarmonicBasis` → `Trivial`; `MirrorEvenSphericalHarmonicBasis` →
`Mirror(axis)` (it carries `mirror_axis`); the future `LegendreBasis` →
`SO2(axis)`. ⚠ **A category question 2.1b's opener must rule:** the three
indicator-family bases and `LossKernelBasis` live on `RealSpace` /
`EnergyGroups` / `IndexSet`, where no subgroup of `O(3)` acts — `None`
(the `mirror_axis` tax again, and §8's branchable value) or `Trivial` (true,
and vacuous)? Read `IndicatorBasis.domain` before choosing.

> ✅ **RULED at 2.1b's opener (2026-09-01): `None` — and `Trivial` for the
> BARE sphere, so three answers by the domain's TYPE.** `Trivial` on
> `RealSpace`/`EnergyGroups`/`IndexSet` would assert an `O(3)` action that
> does not exist; the measure's `phase` already gives that category answer
> for the same manifolds; and §8's hazard does not bite, because G2 (2.2)
> must special-case the MEASURE's `None` (spent nothing) anyway, so the
> basis's `None` costs the predicate nothing new. ✅ And *"0 of 6 answer"* is
> **6 of 6 by inheritance** — the "derivable" list above was right, and the
> derivation is `domain.by`: ONE concrete `@final` property on the ABC,
> zero subclass edits (the requested FIELD dissolved the way 2.0d's did).
> Read `Basis.invariance_group`'s docstring for the HAS/SPENT argument.

### (c) ⛔ 2.2's gate at the FRAME is INERT until 0.1b — XII.1b's fusion re-confirmed at this HEAD

`[M]` `gauss_legendre(8).angular_frame(2).measure.support.name == 'S^2'` — the
FORGED measure — while `gauss_legendre(8).measure.support.name ==
'S^2/SO2_x'`. The frame never sees the rule's declaration, so a G0/G2 written
on `(frame.measure, basis)` admits the slab and catches nothing; one written
on `(quadrature.measure, basis)` refuses the slab **today**, i.e. reddens
every slab solve until 3.4 — the *gate whose only witness is production*
hazard XII.1b names. Land 2.2 with 0.1b + 0.6 + 3.4, as ruled.

### (d) 2.1b's §6c question, answered before it is asked

> ✅ **EXECUTED as `test_e1`, and the witness is STRONGER than an agreement:**
> the fold's two halves are the SAME object — `has is spent` and
> `basis.domain is measure.support` — because both read one memoised
> `Quotient`. The negative leg reads the slab's verdict off the
> QUADRATURE's measure (the frame's still says the forged `S^2`, per (c)).

Alone, 2.1b refuses nothing: it is a field with one consumer (2.2's G2). Its
honest witness at landing is not a refusal but an AGREEMENT on a shipped
pairing: `[M]` the fold's two halves — `folded_product(4,8).measure.quotient_group
== Mirror('y')` and its basis's `invariance_group` (to be `Mirror('y')`) —
satisfy `⊆`, plus the negative leg spelled as the lattice's verdict,
`SO2('x') ⊆ Trivial` is False for the slab rule against a full-SH basis. That
is a real witness (the fold ships) and not a frame refusal (which would be
§6c-vacuous until the fused step).

### (e) ⛔ 2.3 re-measured: the CHART half has consumers, the SECTION half may have NONE

`[M]` `dataclasses.fields(Quotient)` is **10**: `base, by, realization,
fundamental_domain, generators, syzygy, gram, det_gram, derived_by,
singular_stratum`. D0.1's chart (the MAP) and pushforward measure are still
absent (`class Chart` **0**). Consumers a typed arrow would have today:
`pushforward` has **1** production call site (inside `quotient()`,
`measure.py:1156`) and 1 test file; `spherical_product` composes `S²` nodes as
`column_stack([mu_x, mu_y, mu_z])` + a `support=SPHERE` literal
(`rules_product.py:531`) — the Archimedes parametrisation, unnamed; and the
fold's `quotient()`.

⛔ `[R]` **The SECTION half of 2.3 has no consumer once 3.4 lands.** A section
`S²/SO(2)_a → S²` exists to evaluate an `S²`-function at an orbit-space point
— which is ERR-080's forgery, the thing 0.1b retires. What 3.4b's `Descent`
needs is the **quotient map** `π: S² → S²/SO(2)_a`, `Ω ↦ Ω·ê_a`, to pull
`P_ℓ` back (`P_ℓ ∘ π` is the `m = 0` harmonic) — an arrow in the OTHER
direction. Part XIV's own pointer (written at the 2.0c compaction) calls 2.3
*"the SECTION MAP's natural home"*; treat that as a **hypothesis** and settle
it at 2.3's opener by enumerating who would call a section after 3.4. If the
answer is nobody, 2.3 is the quotient map + the product parametrisation, and
`_sphere_mod_so2`'s `fundamental_domain=None` stays `None` on purpose.

`[R]` **A Pattern-2 twin for 2.3/3.1 to absorb:** D0.1's *pushforward
measure* of a quotient is today answered by `AngularSymmetry.reference`'s
dispatch (`LEGENDRE` for any `SO2`, `UNIFORM_ON_SPHERE` for `Trivial`,
raise otherwise). That is the catalogue entry's 8th field living on the
registry; the seed should carry it and the registry should read it.

### (f) Order, and a §6b note

2.1b first (small, one field, the fold as witness); then 2.3 with (e) settled
at its opener. Both are enablers whose consumer is the fused step. ⚠ 2.1b
re-signatures nothing, but it ADDS an abstract member to `Basis` if made
abstract like `domain` — then the §6b set is every subclass constructor
**including test-local ones** (`_DenseTrial` at 2.1: 23, not 18), enumerated by
runtime `__subclasses__` after importing every module (landmine #4).

> ⛔ **VOID at execution (2026-09-01).** It landed CONCRETE and derived, so the
> §6b set was EMPTY: zero subclass edits, zero test-local subclass edits
> (`test_d3`'s stub was hoisted to module scope for section E to reuse, not
> re-signatured). The sizing was right about the abstract design and the
> opener chose the other one.

---

## §V.5f ⭐⭐ 2.0c's THIRD opener (2026-09-01, at `4921d737`) — `support` is an INTERFACE, not a field

> ✅ **BANNER: 2.0c LANDED at `025834f5`.** This opener's central finding — that
> the step named ONE implementor of SIX — shaped the landing: all six were
> retyped in one commit. Kept for its reasoning and its §6b spelling.

The **11th consecutive** phase opener to correct its own section, and the first
to change the step's *unit of work* rather than its sizing. `[M]` all by AST,
each pass carrying positive controls; probes `scratch/_p20c_{census,space_alias,
space_exact,compare,stringops}.py` (untracked, re-runnable).

### (a) ⛔ The tracker names ONE implementor of a SIX-implementor contract

The row reads *"`DiscreteMeasure.support: str → Manifold`"*. `[M]` by AST, the
bare `Space` alias is annotated at **7** sites across **5 classes** — and a
**6th implementor is spelled `str`**, so no alias census can return it:

| class | site | how it is spelled |
|---|---|---|
| `ReferenceMeasure` | `exactness.py:167` | ⭐ a **`Protocol`** — `support` is a structural CONTRACT |
| `UniformMeasure` | `exactness.py:210` | `support: Space` |
| `ProductMeasure` | `exactness.py:284` | `-> Space`, **derived** |
| `GeneratingMeasure` | `generating_measure.py:217` | `support: Space` |
| `DiscreteMeasure` | `measure.py:267` (+ `:744` `new_space: Space \| None`) | `support: Space` |
| ⛔ `AngularSymmetry` | `registry.py:869` | **`-> str`** — the alias's EXPANSION |
| *(test-local)* | `tests/numerics/test_exactness.py:58` `_FourierRef` | `support: Space` |

⭐⭐ **And the sixth is not a bystander — it is `DiscreteMeasure.support`'s
comparison partner.** `[M]` `registry.py:930`, stage 0 of the admission ladder:

```python
def admits_domain(self, measure: DiscreteMeasure) -> bool:
    return measure.support == self.support        # DiscreteMeasure vs AngularSymmetry
```

Retyping one side and not the other makes this `Manifold == str`, which is
**`False`, not an error** — Python compares unequal types silently. ⟹ the
`Space` retype **cannot be split by class**; the six move in one commit or
stage 0 silently refuses every rule (`plan-authoring` §6b).

⟹ **§6b's seventh spelling, for the 2026-08-24 / 08-29 / 09-01 inventory:
an implementor annotated with the alias's EXPANSION rather than the alias.**
`Space = str` means `-> str` and `-> Space` are the SAME annotation to the type
checker and DIFFERENT strings to every census. A `TypeAlias` is invisible at the
only tier a grep or an AST-annotation pass can see.

### (b) `[M]` The two vocabularies both ship, and `manifold.py` says otherwise

| the tag vocabulary | the type vocabulary |
|---|---|
| `measure.py:115-130` — `SPACE_R`, `SPACE_HALF_LINE`, `SPACE_INTERVAL_M11`, `SPACE_INTERVAL_01`, `SPACE_CIRCLE`, `SPACE_SPHERE` | `manifold.py:708-714` — `REAL_LINE`, `HALF_LINE`, `COSINE_INTERVAL`, `UNIT_INTERVAL`, `CIRCLE`, `SPHERE`, `ENERGY` |

⚠ `manifold.py:710-713` comments each new constant *"was `SPACE_INTERVAL_M11`"*
etc., and its section header reads *"under their **retired** tag names"*. `[M]`
**none of them is retired** — `generating_measure.py:398` compares against
`SPACE_INTERVAL_M11` today, and 9 production sites still import the family.
The comments are ASPIRATIONAL wearing the past tense (`coding-standards`'
present-vs-past-tense rule, inverted: here the *past* tense is the false one).
**2.0c is what makes them true**, and it owes them their re-tensing.

### (c) ⭐ `ProductMeasure.support` IS the product algebra, transcribed as string concatenation

`[M]` `exactness.py:285`: `" × ".join(f.support for f in self.factors)`, and
`measure.py:588`: `f"{self.support} × {other.support}"`. `Manifold.__mul__`
(`manifold.py:156`) already builds a `Product`. ⟹ two hand-rolled string
spellings of one operator that ships — Pattern 1 + Pattern 2 in one line each.

### (d) Predicates for the three counts, so none of them is re-read as another

`plan-authoring` §2's FILTER clause: each number below names what it counts.

| predicate | `[M]` at `4921d737` | reconciles with |
|---|---|---|
| `support=` **keyword arguments** (producers) | **87** (32 `orpheus/`, 55 `tests/`) | ✅ §V.5e's 87, exactly |
| ...of which **string literals** | **54**, 51 in `tests/` = 94 % | ✅ §V.5e, exactly |
| **all** `.support` Load reads | **71** (31 `orpheus/`, 40 `tests/`) | *(never counted before)* |
| ...feeding a **string operation** (14 f-string, 2 `.startswith`, 1 `str.join`) | **17** (16 `orpheus/`, 1 `tests/`) | ⚠ §V.5e said **16**; different filter, both stated |
| `support` **comparisons** (`==`/`!=`) | **33** (5 `orpheus/`, 28 `tests/`) | *(never counted before)* |

⭐ The 17 are the quiet ones and they are almost all **repairable mechanically**:
`f"L2[{self.support}]"` → `…{self.support.name}]`. The two `.startswith()` calls
are inside `phase` and the retype **dissolves them**. The 33 comparisons fail in
mixed directions — `admits_domain` goes silently `False`, `generating_measure.py:398`
raises always — which is why the set must land whole.

### (e) ⚖ The design fork this opener opens: WHERE does `phase` live?

§V.5e(c) proposed *"`phase` becomes a property of the MANIFOLD, and the measure's
version is a one-line forward"* — and §V.5e(f) then refuted its totality
(`Interval(-1,1)` is genuinely both the slab's μ-axis and a 1-D spatial
interval). **Both remain true**; what neither settled is the layering:

* `Phase` is a **transport** taxonomy (position × direction × energy) and is
  defined in `numerics/measure.py:132`.
* `Manifold` is **pure geometry** in `numerics/manifold.py`. A `Sphere` is not
  intrinsically "angular" — it is angular *because transport's direction
  variable lives on it*.

✅ **RULED by the user, 2026-09-01: dispatch on the manifold TYPE, `phase`
stays on the measure.** `manifold.py` stays pure geometry — no transport
vocabulary crosses into it. The reasoning below is the ruling's argument;
§V.5e(c)'s *"`phase` becomes a property of the MANIFOLD"* is ⛔ **SUPERSEDED**.

⟹ **the ruled means:** keep `phase` on
`DiscreteMeasure` and dispatch on the manifold's **TYPE** (`Sphere` /
`Quotient(base=Sphere)` → angular; `EnergyGroups` → energy; `RealSpace` →
spatial; `Interval` → ambiguous, defer to `invariance_group`). That kills the
stringly-typed dispatch without moving transport vocabulary into a pure-math
module. The alternative — `Manifold.phase` — is what §V.5e(c) proposed and is
recorded here so the fork is visible rather than silently resolved.

---

## §V.5e ⭐ 2.0c's SECOND phase opener (2026-09-01, at `10e48312`) — unblocked, and the prize is bigger than the retype

> ✅ **BANNER: 2.0c LANDED at `025834f5` (2026-09-01).** This is the opener that
> re-scoped it, kept for its reasoning. Its findings are all DISCHARGED: the live
> `folded_product(…).measure.phase` defect is fixed, the stringly-typed dispatch
> is gone, the `"cells"` arm retired with its synthetic witness, and the 16
> consumers were re-pointed. ⚠ Its `[M]` readings of `.support` quote STRINGS —
> that field is a `Manifold` now.

The 10th consecutive phase opener to correct its own section. `[M]` all by AST
or runtime, each with a positive control (§V.5d(f)4 records a `grep` returning
**zero** on this exact question while AST found all five — the ugrep anchor
hazard). Probe: `scratch/_p20c_census*.py` (untracked).

### (a) ✅ THE BLOCKER IS DISCHARGED

§V.5d(a) blocked 2.0c on the missing `S²/σ_y` catalogue entry. `[M]` at
`10e48312`, `SPHERE.quotient(SubgroupOfO3.Mirror("y"))` returns
`'S^2/sigma_y'` — 1.1's σ_y entry landed at `b55bba56`. **2.0c is open.**

### (b) ⛔ The tracker's sizing counted the PRODUCERS; a retype's risk is in the CONSUMERS

The row reads *"`[M]` 54 literals + 5 f-strings, 94 % of the literals in
`tests/`"*. Re-measured:

| the row's claim | measured at `10e48312` | |
|---|---|---|
| 54 literals | **54** | ✅ reproduces exactly |
| 94 % of literals in `tests/` | **51 of 54 = 94 %** | ✅ reproduces exactly |
| "5 f-strings" | **4** `support=f"..."` construction sites | ⚠ off by one — and it is the wrong question |
| *(never counted)* | ⛔ **16** sites CONSUMING `.support` as a string — 15 `orpheus/`, 1 `tests/` | the set that actually breaks |

⭐ Both halves of the row are about **producing** a support tag; the retype
breaks **consumers** — 12 f-strings interpolating `.support`, 2
`.startswith()` tests, and the `L2[...]` name. A producer that hands a
`Manifold` where a `str` was expected fails loudly at the call; a consumer that
interpolates one silently produces `L2[Sphere()]`. ⟹ the §6b set is
**87 producers + 16 consumers**, and the 16 are where the damage would be
quiet. Same shape as the 2026-08-16 surprise row (a fraction over *inputs*
retargeted at *outputs*) and as 2.1's own §6b miss.

### (c) ⭐⭐ THE PRIZE: `DiscreteMeasure.phase` is stringly-typed dispatch, and the retype dissolves it

`[M]` `measure.py:410-413`:

```python
if self.invariance_group is not None:      return "angular"
if self.support.startswith("spatial") or self.support == "cells":
                                           return "spatial"
if self.support.startswith("energy"):      return "energy"
raise NotImplementedError(...)
```

That is `coding-elegance`'s named anti-pattern — **stringly-typed dispatch** —
and `discriminations`' *"a repeated conditional is a missing type"*, on a tag
whose type is exactly what 2.0c mints. A `RealSpace` knows it is spatial; an
`EnergyGroups` knows it is energy; a `Sphere` or a `Quotient` of one knows it is
angular. ⟹ **`phase` becomes a property of the MANIFOLD**, and the measure's
version is a one-line forward. The three prefix tests, the `"cells"` special
case and the raising fallthrough all go away — replaced by a total function on
a closed sum, which is the same move `Manifold` itself was minted by.

⚠ **This is a GOAL, not yet a design.** Two things must be checked before it is
built: whether `invariance_group`'s angular arm still needs to precede the
manifold's answer (it is a property of the MEASURE, not of the point set), and
what a `Product` manifold's phase is (`[R]` a product of factors with different
phases has none — which may be the honest answer, and is a case the string
version could not even pose).

### (d) `[M]` `.phase` has **ZERO** production consumers — 7 reads, all in `tests/`

So the dispatch above is exercised only by its own gates today. ⚠ That is not a
licence to skip it: 0.1a's landing found `frame.measure.phase` **RAISING on all
12 rules** because the forgery had destroyed `invariance_group`, and it was a
test that caught it. It is a live surface with a latent consumer
(`process-discipline`: "no current consumer" ≠ speculative), and its dispatch is
the thing 2.0c makes honest.

### (e) Landmine #3 holds — but ⛔ the reasoning below was WRONG, and the way it was wrong is the transferable part

*(Original, kept per §3:)* "✅ Landmine #3 CONFIRMED by an independent census.
§V.5c(d) recorded that `measure.py:411`'s `self.support == "cells"` has **no
producer**. `[M]` the 10 distinct shipped literals are `'[-1,1]'` ×14, `'S^2'`
×13, `'spatial_R1'` ×10, `'R'` ×9, `'energy'` ×2, `'[0,1]'` ×2,
`'spatial_R2'`, `'index(angular)'`, `'[-1,1]^slab'`, `'probe'` — **`'cells'` is
not among them**, and no f-string site can produce it. It is a dead arm."

⛔ **REFUTED as reasoning, same day, by its own author.** `'cells'` appears at
**7** sites, 2 of them in `orpheus/` — and a search for those turns up
`radial_characteristic_space.py:293` (`_RayLeg(part="cells", …)`) and `:488`
(`_PART: ClassVar[str] = "cells"`), neither of which is a `support`. So the
CONCLUSION survives for production. But the one real producer is
`tests/numerics/test_measure_phase.py:67` — a `parametrize` list
`["spatial_R1", "spatial_R2", "cells"]` feeding `support=support`.

⭐⭐ **And my own census had already found it.** It is the row
`tests/numerics/test_measure_phase.py:72 (Name)` in the NON-LITERAL half of the
same table — I answered the question from the LITERAL half and stated the
result over the whole population. ⟹ the mechanism, which generalises past this
campaign: **a census that helpfully separates "literal" from "non-literal"
invites you to answer from the literal half, and a parametrized producer lives
in the other one.** The split is a convenience for reading and a trap for
concluding; a universal must be stated over the union or not at all.

⟹ **the corrected claim**: `support="cells"` has **zero production producers**
and **exactly one test producer, which exists to assert the arm**. That is not
a dead arm — it is a *production-unreachable* arm with a synthetic witness,
which `coding-standards` treats differently: retiring it means retiring its
gate row too, not silently transcribing both into the new type.

### (f) ⛔⛔ A LIVE DEFECT the opener found: the shipped FOLD cannot say it is angular

`[M]` at `10e48312`, over the five shipped rule families:

| rule | `support` | `invariance_group` | `phase` |
|---|---|---|---|
| `gauss_legendre(8)` | `[-1,1]` | `sigma_x` | `angular` |
| `lebedev(17)` | `S^2` | `Oh` | `angular` |
| `level_symmetric(4)` | `S^2` | `Oh` | `angular` |
| `product(4,8)` | `S^2` | `D_8h` | `angular` |
| ⛔ **`folded_product(4,8)`** | `S^2/sigma_y` | **`None`** | **RAISES** |

The fold's `invariance_group` is legitimately `None` — the mirror was
*quotiented away*, so the measure is not invariant under it — and
`'S^2/sigma_y'` matches neither `startswith("spatial")` nor
`startswith("energy")`. So the string dispatch has **no arm for a quotient**,
and the production folded-cylinder measure raises when asked what phase-space
factor it belongs to. ⚠ Latent: `[M]` 0 production consumers of `.phase`. It is
still a defect, and it is the SECOND time this campaign has found `phase`
raising on an angular object (0.1a found `frame.measure.phase` raising on all
12 rules, from a different cause).

⭐ It also answers §V.5e(c)'s first open question, and answers it **against** the
simple design: `phase` **cannot** be a total function of the manifold alone
while the slab's support is `[-1,1]`, because `Interval(-1,1)` is genuinely
ambiguous — it is the slab's angular μ-axis AND exactly how a 1-D spatial
interval would spell itself. The current code resolves that with
`invariance_group`, a property of the MEASURE, not of the point set.

⭐⭐ **And the deeper reading: `support='[-1,1]'` on a slab angular rule is
ERR-080's defect class again — a manifold named by its CHART CODOMAIN instead
of by the orbit space it is.** `[M]` §V.5d(b) already established that
`S²/SO(2)`'s `realization` **is** `[-1,1]`. So the honest support is
`SPHERE.quotient(SO2)`, whose realization is that interval — and with it, every
one of the five rows above becomes angular *by derivation from the manifold*,
the `Interval` ambiguity disappears, and the `invariance_group` arm is no
longer load-bearing. ⚠ That declaration is tracker **2.4**, NOT 2.0c. 2.0c
must therefore keep the `invariance_group` arm and add a quotient arm; the
collapse to a single manifold-derived answer is 2.4's to finish, and 2.4 now
has a measured reason to exist beyond tidiness.

## §V.5d ⛔ 2.0c/2.0d's phase opener (2026-08-31) — ✅ ITS BLOCKER IS DISCHARGED; 2.0d dissolves into 2.0c

> ✅ **BANNER, read before the section. Everything below about 2.0c being
> BLOCKED is HISTORY** (kept per §3, which requires the refuted text to stay).
> 1.1's σ_y catalogue entry landed at `b55bba56`, §V.5e(a) discharged the
> blocker, and **2.0c LANDED at `025834f5`** on 2026-09-01. `[M]` at
> `8cc69e7f`, `SPHERE.quotient(SubgroupOfO3.Mirror("y"))` answers
> `'S^2/sigma_y'` rather than raising, and `_ORBIT_CATALOGUE` is no longer a
> one-entry table. ⚠ Its `folded_product(4,8).measure.support` readings are
> also pre-2.0c: that field is a `Quotient` OBJECT now, not the string this
> section quotes. The section is retained for its *reasoning*, which is what
> re-scoped the step — not for its status claims.

The 9th consecutive phase opener to correct its own section. `[M]` all by AST
or by runtime probe, 2026-08-31 evening, at `a9bfebab`; probes in
`scratch/p20c_*.py` (untracked).

### (a) ⛔⛔ 2.0c is BLOCKED — the shipped fold has no catalogue entry

`[M]` `_ORBIT_CATALOGUE` (`manifold.py:633`) has **exactly ONE** entry,
`(Sphere, "SO2")`, plus the `Trivial` group handled by derivation. The shipped
cylindrical fold's support is `'S^2/sigma_y'`
(`[M]` `folded_product(4,8).measure.support`), and

```
SPHERE.quotient(SubgroupOfO3.Mirror("y"))
  -> NotImplementedError: no catalogue entry for S^2/sigma_y
```

⟹ **`DiscreteMeasure.support: str → Manifold` cannot land until `S²/σ_y` is
derived and registered.** That is tracker **1.1**, which the tracker orders
*after* Phase 2 and marks "⏏ partial". The dependency is real, forced, and was
unrecorded. **1.1's σ_y entry precedes 2.0c.**

### (b) ⛔ And 1.1's σ_y entry is itself blocked on a REPRESENTATION ruling

`[M]` the two shipped quotients use **incompatible node representations**:

| quotient | shipped `measure.nodes` | what it is |
|---|---|---|
| `S²/SO(2)` — `gauss_legendre(8)` | shape `(8,)`, 1-D | the **invariant coordinate** `p₁ = z = μ`; realization `[-1,1]` |
| `S²/σ_y` — `folded_product(4,8)` | shape `(16,3)`, unit norm, `μ_y ∈ [+0.1945, +0.8688]` | a **SECTION** — one representative per mirror pair, in the ambient sphere |

Procesi–Schwarz returns the orbit space in *invariant* coordinates; `[R]` for
σ_y those are `(x, z, y²)` and the sphere ideal cuts them to a closed 2-disk.
The tree ships the *section* (the closed upper hemisphere in ℝ³). `Quotient`
has one `realization` field and `Quotient.contains` delegates to it
(`manifold.py:509`), while `ambient_dim(Quotient) = ambient_dim(realization)`
(renamed public at 2.1) — so
the choice decides whether the shipped 16×3 folded nodes are ADMITTED or
REFUSED. A disk realization demands 2 columns; the nodes have 3.

⟹ this is tracker **3.3**'s retraction/section question (`OrbitAxis`), arriving
early because 2.0c needs it.

⛔ **REFUTED 2026-09-01 by the docs pass, and the refutation is sharper than the
claim.** *(Original: "⚠ It is not visible from the SO(2) entry: there the
invariant chart and the section coincide in dimension, which is exactly why one
worked entry did not expose the fork.")* `[M]` the dimensions coincide in
**BOTH** entries — `[-1,1]` 1 against an `SO(2)` half-meridian 1, `D²` 2 against
the σ_y hemisphere 2 — so dimension cannot discriminate, and **it cannot, by my
own construction**: `Quotient.__post_init__` now GATES that agreement, so it is
true of every entry that can be built.
⭐ The mechanism worth carrying: the sentence was written *before* the gate
existed, was diagnostic then, and my own step made it vacuous — a claim
invalidated by its campaign's success, which is `plan-authoring` §10's shape
moved from a metric onto a sentence. Nothing prompts the re-check, because the
sentence still reads true.
⟹ the two reasons that DO reproduce: (i) no canonical section exists for a
positive-dimensional group, so `SO(2)` has nothing to put in the second slot;
(ii) the tree's `SO(2)` data is **already chart coordinates** (`(8,)`) while the
fold's is a section (`(16,3)`), so only the fold's data can contradict a
chart-only reading.

### (c) ⭐⭐ 2.0d DISSOLVES — the field it proposes is a second home for a fact `Manifold` already carries

2.0d asks for `DiscreteMeasure.quotient_group`. `[M]` `Quotient.by` **IS** that
group (`SPHERE.quotient(SubgroupOfO3.SO2).by is SubgroupOfO3.SO2`). Once
`support: Manifold`, the field is a one-line derived property with **zero new
state**:

```python
@property
def quotient_group(self):
    return self.support.by if isinstance(self.support, Quotient) else None
```

⭐ And the derived form is *strictly better than a stored field*, because
`quotient_group` is a property of the **support**, not of the nodes — so it must
travel exactly where `support` travels. `[M]` the four verbs disagree on what
they forward: `consolidate` preserves `support` (`:882`), `restrict` preserves
it (`:1070`), `pushforward` REPLACES it (`:805`), `quotient` derives a new one
(`:1022`). A stored field needs that table maintained by hand in four places and
kept in agreement with `support`'s; the derived property cannot disagree with
`support` because it *is* `support`. Pattern 2, caught before it was minted.

⟹ **2.0d is not a separate landing.** It is one property + its witness, riding
2.0c. The tracker row is re-scoped, not dropped.

### (d) ⛔ 2.0d's done-when was unreachable as written

Stated: *"`gauss_legendre(8).measure.quotient_group is SubgroupOfO3.SO2`,
derived by construction, not declared"*. `[M]` `gauss_legendre(8).measure`
carries `support='[-1,1]'` and is built **directly**, never through
`.quotient()` — so "derived by construction" gives `None`, at 2.0d and at 2.0c
alike. The `SO2` answer is a **declaration**, and the plan already owns it as
tracker **2.4** (*"the slab rule declares its quotient group"*). The done-when
silently conflated the two steps.

⟹ corrected, and now split by what would falsify each:
* **2.0c (derived):** `folded_product(4,8).measure.quotient_group is
  SubgroupOfO3.Mirror('y')` — `[M]` reachable, since `folded_product` ends in
  `full.quotient(SubgroupOfO3.Mirror("y"))` (`directional.py:845`).
* **2.4 (declared):** `gauss_legendre(8).measure.quotient_group is
  SubgroupOfO3.SO2`, and `[M]` at 2.0c it is `None`.

### (e) ⭐ The `L2[coarse_cells_Rd]` false name WAS an ILLEGAL STATE THAT IS REPRESENTABLE

✅ **REMEDIED 2026-09-01 by tracker 2.1 @ `c461fe8d`** — everything below
describes the tree BEFORE that commit and is kept per §3, because the
measurement is what justified the repair. `[M]` at HEAD the two spaces are
neither `==` nor hash-equal
(`tests/numerics/test_basis_domain.py::test_d1`).

§III.10 records the name as *false*. `[M]` it is worse than false — measured
with both legs:

| | |
|---|---|
| a 2-group **ENERGY** indicator space (`EnergyGrid.as_basis().space`) | `name='L2[coarse_cells_R1]'`, `shape=(2,)` |
| a 2-cell **SPATIAL 1-D** indicator space | `name='L2[coarse_cells_R1]'`, `shape=(2,)` |
| ⛔ `energy == spatial` | **`True`**, and `hash` equal |
| ✅ NEGATIVE CONTROL `energy == FunctionSpace("L2[S^2]", (2,))` | `False` — the predicate does discriminate when names differ |

`FunctionSpace.__eq__`/`__hash__` are `(name, shape)` (`space.py:124`'s class
docstring: *"``name`` is the **identity** of the space"*), and `[M]` **26** sites
in `orpheus/` compare a space attribute with `==`/`!=`. So every composability
guard in the tree admits an energy↔spatial confusion on these two spaces.

⚠ **Honest scope, per `plan-authoring` §8's own sharpening (a branch's EXISTENCE
is not its RELEVANCE):** what is measured is that the TYPE cannot distinguish
them. Whether any production path puts these two on either side of one of the 26
guards is **NOT measured**, and no wrong answer is claimed. The deliverable is a
WITNESS, not a warning.

⟹ this clause needs **2.1**'s constructor field to be fixed *by construction*
(§II.15 R1: 18 `IndicatorBasis` sites over three families; any derived value
hard-codes one of three). ⛔ **2.0c's done-when cannot own it** — `plan-authoring`
§6b, with a *false name* in place of a call site. The clause moves to **2.1**.

### (f) Four smaller corrections from the same census

1. ⛔ §V.5c(e)'s *"two spellings of one quotient, **both shipped**"* — `[M]` all
   five odd literals (`'img'`, `'probe'`, `'[-1,1]^slab'`, `'S^2/<sigma_y>'`,
   `'S^2/sigma_y'`) are in **`tests/`**; `orpheus/` has **zero**. Production has
   ONE fold spelling, and `Quotient.name` reproduces it **bit-identically**
   (`[M]` `'S^2/sigma_y' == f"{SPHERE.name}/{Mirror('y').name}"`). The
   identity-collapse hazard is a test-fixture matter, not a shipped one.
2. ⭐ **Name identity is preserved for every production support tag.** `[M]`
   `Manifold.name` reproduces `[-1,1]`, `S^2`, `S^2/sigma_y`, `spatial_R{d}`,
   `index({label})`, `R`, `[0,1]`, `energy` verbatim — so the retype does not
   move `FunctionSpace` identity anywhere it matters. **Two exceptions**, both
   real name changes: `indicator_basis.py`'s `L2[coarse_cells_R{d}]` (which
   is the defect, per (e); ✅ **landed at 2.1 @ `c461fe8d`**, ahead of 2.0c)
   and `loss_kernel_gauge.py:1169`'s
   `sn_trace_orbit{o}_g{g}`, which under `IndexSet` becomes
   `index(orbit{o}_g{g})` — a space-identity change on the SN trace-gauge path,
   owed a gate.
3. ⭐ **The `φ_*(M)` pushforward fallback has exactly ONE consumer, and it is a
   test.** `[M]` by AST, 8 `pushforward` call sites: 7 pass `new_space`
   explicitly (including the sole production one, `measure.py:1020` inside
   `quotient`); the 1 that omits it is `tests/numerics/test_measure.py:331`,
   whose only job is to assert the fallback. ⟹ making `new_space` REQUIRED is a
   one-site change, and it is the principled one: a pushforward cannot know its
   codomain without the chart, so the fabricated tag `φ_*(M)` names nothing —
   §V.5c(e)'s own defect, in the production verb. `Chart` (2.3) then supplies it
   structurally.
4. ⛔ `[M]` a `grep` for these five literals returned **zero** while the AST
   census found all five — the L61(c)/`nexus-tools` ugrep hazard, live again.
   Every count in this section is AST or runtime, each with a positive control.

### ⭐ RULED 2026-08-31 (user) — `Quotient` carries TWO coordinate systems

§V.5d(b) put the fork to the user; the ruling is **two slots**, and the
reasoning is worth more than the choice.

`realization` **keeps** its documented meaning — the chart codomain, in the
INVARIANTS' language, the same language as `generators`, `gram` and `det_gram`,
so it stays consistent with its six neighbours. A new
`fundamental_domain: Manifold | None` carries the **section's image**, in the
BASE's coordinates — which is what `DiscreteMeasure.quotient` actually emits,
since it keeps `nodes[representative]`, a selection applying no chart.
`Quotient.contains` accepts either and dispatches on ambient width;
`ambient_dim` (public since 2.1) still reports the realization's, because a `Product` factor needs one canonical
width.

**⛔ The three alternatives, all measured, all refused** (`process-discipline`:
a refuted candidate is first-class output, with its structural reason):

| candidate | why it fails |
|---|---|
| `realization = SPHERE` | refuses the ERR-080 forgery but **ADMITS the mirrored orbit twins**; `[M]` `Quotient.contains` becomes bit-for-bit `SPHERE.contains` — a predicate that cannot distinguish `M/H` from `M`. Also topologically false: `D² ≇ S²` (`χ = 1` vs `2`) |
| the disk ALONE (chart only) | ⛔ **Mode-12 BLIND to ERR-080** — the chart `(x,y,z) ↦ (x,z)` drops exactly the coordinate the forgery corrupts, so `(μ,0)` is a legal disk point. `[M]` `max \|(μ,0)\|² = 0.9221566084920586` over the 8 forged rows, all inside |
| `Product(COSINE_INTERVAL, COSINE_INTERVAL)` | the bounding **square**, not the disk. `[M]` witness `(0.9, 0.9)` is in the square, not the disk, and is **no direction at all** (`η²+μ² = 1.62 > 1` forces `ξ² < 0`) |

**⭐⭐ And the finding the fork produced, which is larger than the fork:
ERR-080 is a botched SECTION of `S²/SO(2)`.** `[M]` the tree's 3-D embedding of
a 1-D rule is `(μ, 0, 0)` — norms `0.960, 0.797, 0.526, 0.183`,
`SPHERE.contains → False` — while an honest `φ = 0` half-meridian
`μ ↦ (√(1−μ²), 0, μ)` is on the sphere to `0.0`. The realization is a **chart**;
a consumer needed a **section**; the tree fabricated one by zero-padding. With
the second slot that padding has nowhere to live.
⚠ **Level-1 only.** It makes the nodes points of the manifold. The level-2 half
is untouched and is still the trivial isotypic sub-basis (tracker 3.4).

**The structural asymmetry** — ⛔ and note it is NOT a dimension mismatch
(refuted above, §V.5d(b)): what the first entry lacked is a *section to put in
the second slot* and *data that could contradict a chart-only reading*. For a
**positive-dimensional** group the chart is a genuine reduction
(`3 → 1` floats) and no section is canonical; for a **finite reflection** the
chart is no reduction (`3 → 2`, the third recoverable) and the section IS
canonical — a *strict* fundamental domain. Both choices are locally correct;
the TYPE was what had to serve both.

### ✅ What landed with the ruling

* **`Ball(d)`** — the closed unit ball, `name = "D^d"`. Minted because `S²/σ_y`
  IS the closed 2-disk in invariant coordinates and nothing could say so.
* **`FundamentalDomain(base, normals, label)`** — `base ∩ ⋂{⟨p,n⟩ ≥ 0}`.
  ⭐ An antipodal PAIR spells an equality, so ONE tuple expresses both the σ_y
  hemisphere (`dim 2`) and an `SO(2)` half-meridian (`dim 1`), and
  `dim = base.dim − (antipodal pairs)` is what the new `Quotient.__post_init__`
  gates against the realization's `dim` — a real check, since the two derive
  their dimension by different routes.
* ⛔ **The half-spaces are CLOSED, not strict**, and the reason is measured, not
  stylistic: `[M]` the cylindrical march seeds every level at `ξ = 0` exactly —
  ON the stratum — so `> 0` would refuse a direction production marches from.
  `coding-elegance` #18's half (ii), which is a claim about the PRODUCERS and is
  the half that gets skipped.
* **`_sphere_mod_mirror`** — one derivation, three catalogue keys
  (`sigma_x/y/z`); it reads the axis off the group. Syzygy empty **by theorem**
  (Chevalley–Shephard–Todd: a reflection group's invariant ring is polynomial),
  `P = diag(1,1,4p₃)`, `det P = 4p₃`, realization `Ball(2)`.
* ⛔ **`singular_stratum` RETYPED** `tuple[float, ...] → Any | None` — a SymPy
  locus in the realization's coordinates, `None` ⟺ free. `S²/σ_y`'s stratum is
  the disk's boundary **circle**; the first entry's SHAPE had become the field's
  type. It stays STORED (not derived) for a reason that survives review:
  recovering it needs the BASE's ideal (`4p₂ → 4(1−μ²)` only after
  `p₁²+p₂ = 1`), which a `Quotient` does not carry — the type cannot recompute
  it, which is exactly when storing is right. Contrast 2.0d's `quotient_group`,
  which the type CAN recompute and therefore must not store.

**`[M]` the acceptance measurement, on production data**, all four legs:

| input | verdict |
|---|---|
| `folded_product(4,8).measure.nodes` `(16,3)` | **ADMITTED** (section) |
| the mirrored orbit twins (`μ_y → −μ_y`) | **REFUSED** |
| the ERR-080 forgery `(8,3)`, not unit-norm | **REFUSED** |
| the charted form `(16,2)` | **ADMITTED** (chart) |

Mutation battery, per arm, restored byte-identical against a copy-aside:
whole-predicate control **3 reds** (both arms) · section dispatch disabled **2**
· half-spaces dropped **1** · strict-`>` **1** · `Ball` opened **1** ·
construction invariant off **1** — each the gate that names it. ⚠ The first
control attempt reddened only 1 and was **not** a control: it mutated one arm of
the new two-language dispatch (`vv-principles` #17's granularity trap, live).


### The re-ordered exit path

⟹ **1.1(σ_y) → 2.0c (retype + the derived `quotient_group`) → 2.1 → 2.1b →
2.2**, with 2.4 owning the slab's *declared* group. 2.0d is folded into 2.0c;
2.0c's `indicator_basis` clause moves to 2.1.


# Part VI — The action plan

Phases are ordered by dependency, not by size (D0.5). The `Q#` column maps each
item to the source plan's numbering so the two stay reconcilable. `[M]` the source plan is now IN-REPO at `.claude/plans/symmetry_quotient_plan.md` (399 lines) — read it for the Procesi–Schwarz derivation, the `Retract`-pair framing, and the confidence ledger this plan condenses.

## ⚡ Phase 0 — IMMEDIATE. Land regardless of any ruling below.

**Goal.** No object asserts a property its own data contradicts.
Every item here is true under every design choice in Phases 1–7.

| # | item | site | Q# |
|---|---|---|---|
| **0.1** | ⚡ **Stop fabricating the support, and delete the phantom read in the same expression.** `angular_frame` (`[M]` 2026-09-01: the site is `directional.py:624`, not the `:587-593` this row cited) must not write `support=SPACE_SPHERE` over nodes with `‖Ω‖ ≠ 1`. §II.5 proves L1 and L4 are ONE fix. ⛔ **SPLIT at execution 2026-08-31 — see §II.9**: `0.1a` (3-D rules hand the frame their OWN measure; `[M]` bit-identical on 10 of 12 shipped rows) lands now; `0.1b` (the 1-D rows) **rides Phase 3.4**, because the lift `[-1,1] → S²` *does not exist as a map*; `0.1c` (the fold's `S^2/sigma_y` tag, §II.8's new leak) is gated on a `plan-authoring` §8 blast-radius measurement of `FunctionSpace` equality. | `directional.py:587-593` | — |
| **0.2** | ⚡ **A phantom component becomes unspellable.** *Proposed means as of 2026-08-31, SHARPENED at execution — see ⚠ below:* `axis_cosines(i)` raises for `i ≥ dim`. ⛔ a `raise`, **not** an `assert` — the canonical runner is `python -O`. ✅ **census DISCHARGED (§II.14) — and it REFUTED this means**: `[M]` a blanket raise breaks **3 legitimate production consumers** asking the FLOW question. ⟹ runs AFTER 0.1 (which deletes the one fabrication site, 92 % of all reads); 0.2 is then the accessor SPLIT. Done-when re-runnable: re-run the census, `directional.py:589` must be GONE. ⛔⛔ **BOTH clauses of that last sentence are REFUTED (2026-09-01, §II.14):** `[M]` 0.1a removes **0 %** of the phantom reads — only a 1-D rule can phantom-read, and 0.1a routes the 3-D ones — so the 92 % rides **0.1b** at the fused step; and the `:589` row did not go away, it MOVED (to `_harmonic_frame_measure`) and is retired by 3.4, not here. ✅ 0.2 landed anyway, decoupled, at `ce46181c` — see §0.2-R and the 0.2 execution section. | `directional.py` | ✅ `ce46181c` |

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
| **2.0a** | ⭐ **MINT `Manifold`** (D0.7). The thing functions are defined on and measures live on: `S²`, `S¹`, `[-1,1]`, `[0,1]`, `[0,∞)`, `ℝ`, `spatial_Rᵈ`, `energy`, `index(axis)`, and quotients/products of those. Verbs: `quotient(H)` (⭐ **catalog lookup, and a `raise` naming the missing WORK** — *"no entry for `S²/C₆`; derive it per §III.1 and register it, or implement the engine"* — so the engine's absence is a pick-up-able work item, not a wall; D0.1), `__mul__`, `contains(points)`, `dim`, and the **singular stratum** a quotient inherits from `H`'s fixed-point set. ⛔ a VALUE type with methods — **not** a `Generic[Tag]` phantom parameter (`sn_reshape.md` Issue 2 stands; D0.7 explains why it does not bar this). | `Space = str` is retired. ⛔ **The three counts below are REFUTED as written — see §V.5b.** *(Original text, kept per `plan-authoring` §3: "`[M]` today 29 `Space` refs, 39 `support=` construction sites, 45 `.support` reads — that is the migration's size, not a veto (D0.5)".)* `[M]` re-censused by AST 2026-08-31: **87** `support=` kwargs / **62** `.support` reads / **8** `Space` names / **39** `SPACE_*` names, and **104 of the 196 are in `tests/`** — the tree the original counts never ranged over. Size is still not a veto (D0.5); the SHAPE is what moved | — |
| **2.0b** | ⭐ **`contains` is the membership PREDICATE the string never had.** `SPACE_SPHERE.contains(nodes)` is `‖Ω‖ = 1`. This is what makes the ERR-080 fabrication refusable **at the measure constructor**, three hops before the symptom. | `DiscreteMeasure.__post_init__` refuses nodes its declared manifold does not contain; the Part I forged measure is unconstructable | — |
| **2.0c** | ⭐ **`FunctionSpace` carries its `Manifold`; the two `L2[...]` NAME STRINGS become derived** (§III.10). ✅ **The BASIS half landed early at 2.1 @ `c461fe8d`**; what remains is the MEASURE's `f"L2[{support}]"`. | `grep` finds no `f"L2[...]"` literal; `[M]` `indicator_basis.py:284`'s live false name (`L2[coarse_cells_R1]` on an ENERGY basis) is gone by construction | — |
| **2.0d** | ⭐ **`measure.quotient_group`** — the field Part IV's own predicate needs and that `[M]` **exists nowhere**: `folded_by` records the DISCRETE fold on the `Quadrature` façade, and *"the slab's quotient group is recorded nowhere"* (Part IV obstacle 2). Derived by `quotient()`, `None` for an unquotiented measure. | `gauss_legendre(8).measure.quotient_group is SubgroupOfO3.SO2`, derived by construction, not declared | — |
| 2.1 | **`Basis.domain: Manifold`** (closes L2). ⛔ the 2026-08-31 `support` rename is REFUTED — §III.10: the collision is unreachable (`[M]` 0 `Basis` subclasses are operators) and `support` is mathematically false for a basis (`supp(f)` = where `f` is non-zero; for an indicator that is ONE cell). *Means:* a `domain` property beside the existing `space`, giving the basis the same two-level structure the measure has. ⚠ **`IndicatorBasis` takes it as a CONSTRUCTOR FIELD** (§II.15 R1: `[M]` ⛔ *"over 5 sites"* CORRECTED 2026-08-31 — by AST there are **18** `IndicatorBasis(...)` sites (**4** in `orpheus/`, 14 in `tests/`) over **three** manifold families; the *three* was right, the *five* was not. Any derived value hard-codes one of three — and `[M]` that lie already ships in its space NAME). | every shipped `Basis` answers; `SphericalHarmonicBasis.domain is SPACE_SPHERE`; the energy indicator no longer claims a spatial manifold | — |
| **2.1b** | ⭐ **`Basis.invariance_group`** — the OTHER side of Part IV's G0 predicate (`measure.quotient_group ⊆ basis.invariance_group`), ~~which today has **neither** side~~ ✅ both sides real — the measure's at 2.4, the basis's at 2.1b. Derivable for every shipped basis: full SH → `Trivial`; `MirrorEven` → `Mirror(axis)`; the trivial isotypic sub-basis → `SO2`. | the four-row lattice table in Part IV runs as a test and reproduces its verdicts | ✅ landed 2026-09-01 `9b4a4d9c` — read the *2.1b EXECUTED* section. ⛔ the row's MEANS (*a field answered by six subclasses*) DISSOLVED at the opener: the group is `domain.by`, read by ONE concrete `@final` property on the ABC — zero subclass edits, the way 2.0d dissolved into 2.0c. Category ruled: `None` off a non-angular domain, `Trivial` off the bare sphere. `[M]` the four-row table runs as `test_e5`; the fold's two halves read ONE object; ERR-080's pairing is the verdict `Trivial ⊉ SO2('x')`, refused by nothing yet (2.2). numerics 2679 → **2690**, reconciled |
| 2.2 | **G0 at frame construction** (closes L3). ⭐ **§II.7: the PREDICATE already exists** as `AngularSymmetry.admits_domain` (`measure.support == self.support`) — `[M]` reachable only from `select_quadrature`, whose 21 callers are **all tests**; and `[M]` `solve_sn` is HANDED a quadrature (required positional, no default), so the selector was never on the path. ⛔ **REFINED by §II.10 — NOT "just wiring".** The predicate is correct over an **incomplete domain vocabulary**: `support` derives from `continuous_isotropy` alone, so `folded_product`'s `S^2/sigma_y` matches **no** geometry and stage 0 would **refuse the shipped cylinder**. 2.2 owes the ontology a second slot (the *discrete* quotient the rule took — `SubgroupOfO3` already names it) **and** the `(CoordSystem, ndim) → key` join, which `[M]` exists nowhere. Do not hand-roll a second predicate; do extend this one's vocabulary. | `GalerkinFrame(SphericalHarmonicBasis(L=2), slab_measure)` **RAISES**, with a test witnessing it on the exact pairing that ships today | — |
| 2.3 | ⭐ **The typed `Chart`** (closes L6, L7). *Proposed means:* a map carrying `domain → codomain`; `pushforward` derives `support` from the chart rather than taking `new_space`; `product()` composes its two 1-D factors through the Archimedes chart instead of `column_stack` + a literal. ⛔ **The NAME was refuted at the ruling (2026-09-02)** — a chart is `M ⊃ U → ℝⁿ` and only the inverse of the Archimedes map is one; user-ruled **`ManifoldMap`**, one frozen value type, with `archimedes(axis)` and `barycentre(orbit_space)` as memoised factories. | `grep` finds no `SPACE_SPHERE` literal at a measure constructor (⛔ scoped at §V.5i(c) to MAP-BUILT measures — `[M]` an AST tell over `spherical_product`'s source runs as a test); `product(4,8).measure.support` is **derived** (`is archimedes("z").codomain`); the factor structure is recoverable from the built rule | ✅ landed 2026-09-02 — read the *2.3 EXECUTED* section. All three maps typed; `new_space=` retired (1 + 8 sites); `[M]` bit-identical on 60 product / 5 folded / 15 embedding configurations; the reference measure DEFERRED to 3.1 with the §6d cycle measured (5 of 5 orders); numerics 2690 → **2711** |
| 2.4 | **The slab rule declares its quotient group** (closes L5). Includes settling Part IV's obstacle 1 — the `SO2` axis-convention collision — and obstacle 2, the `invariance_group` vs `folded_by` conflation. ⛔ **MOVED to Phase 0 on 2026-08-31** — the full-suite phantom census was scheduled here, but **0.2 depends on it** (it is 0.2's denominator, not 2.4's). Scheduling a step's own precondition three phases downstream is the `plan-authoring` §6b defect with a *census* in place of a call site. Run and published at Phase 0. | `SubgroupOfO3.SO2.is_invariant(gauss_legendre(8).measure)` returns the derived answer, and the plan records whether that answer is `True` **with the derivation**, not the prediction | — |

## Phase 3 — The symmetry core

**Goal.** A quotient space, its retraction/section pair, and its invariant
sub-basis are all *derived from a catalog entry*.
⚠ **Phase 3 does not ship alone** — it lands with Phase 4.1, its first consumer.

| # | item | depends on | done when | Q# |
|---|---|---|---|---|
| 3.1 | ⛔ **RE-SCOPED at its opener (§V.5j) and LANDED 2026-09-02 as the entry's quotient map + pushforward reference — see the *3.1 EXECUTED* section. The catalogue HOME in this row's premise landed at 2.0a/1.1/2.4 in `manifold.py`, `class IsometryGroup` is `SubgroupOfO3`, and the package re-home is OFF the exit path.** *(Original row, kept per §3:)* `numerics/symmetry/` catalog: `IsometryGroup` entries from 1.1. ⭐ **D0.1's SEED requirement is binding here**: an entry's FIELDS are §III.1's OUTPUTS — invariant generators `p₁…p_k`, syzygy ideal `I`, `P_ij = ⟨∇pᵢ,∇pⱼ⟩` and `det P`, chart, pushforward measure, and the stratum where `det P` vanishes — **not** a human summary of the answer, so a future engine populates them without minting a single new type. Plus the action on the axis index set, Haar (or counting) measure, a `derived_by: "hand" | "engine"` provenance field, and a symbolic test reproducing its own derivation (⭐ those ~12 tests ARE the engine's acceptance suite, written before it). ⚠ `[M]` `numerics/symmetry.py` is a **module** today — this is a re-home; `plan-authoring` §6d (check the import edge BOTH ways, and look for the import-linter's declared contract first). | 1.1 | every entry's chart and measure reproduced by its own regression test | Q1.1 |
| 3.2 | `ReynoldsProjection` as a first-class operator with fixed domain/codomain, built via `.on(V)` like every other bound arrow, factored `π* ∘ 𝔄`. | 3.1 | `ℛ² = ℛ` and `𝔄π* = id` to machine precision | Q1.2 |
| 3.3 | `OrbitAxis` minting: `Quotient(axis, H) → (OrbitAxis, retraction/section)`. **Re-uses `AxisRetractionOperator` + `FunctionSpace.section`** (§III.3) — no new relational type. Stratification descriptor populated from `H`'s fixed-point set. | 3.1, 1.8 | the existing retraction/section tests pass with the new instance; no new relational type added | Q1.3 |
| **3.4** | ⭐ **`H`-stable sub-basis: the TRIVIAL ISOTYPIC COMPONENT.** `SphericalHarmonicBasis` → `{Y_ℓ0} ≅ {P_ℓ}` under `SO(2)`, mask **DERIVED by probe** (§III.8), generalising `MirrorEvenSphericalHarmonicBasis` so that it becomes an instance rather than a sibling. **THIS IS THE FIX for Part I.** | 3.1 | §VII's done-when 1 and 2 both hold | Q1.4 |

| **3.4b** | ⭐ **`Descent`** — witness `Funcs(M)^H ≅ Funcs(M/H)` (3.4-R2, RULED (i)). Lands **with or before 3.4**, which is what creates the second realization. Carries the discriminator (*downstairs when the quotient has a classical basis; upstairs otherwise*) where both classes can see it. | 3.1, 3.4 | analyse-then-descend ≡ descend-then-analyse to machine precision, on BOTH shipped pairs (`S²/SO(2)`, `S²/σ_y`) | — |

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
| ⛔ **CORRECTED 2026-08-31 — this row was never a refusal.** General Gröbner / Procesi–Schwarz orbit-space engine | **DEFERRED, not refused** (D0.1): ~12 catalog entries ever, so it is not yet worth its failure modes — but the embryo is built as its **seed**, and the day it is worth building it must be a *development* of that seed. The falsifiable check is in D0.1. Listing it here as "refused" is the tense error `plan-authoring` §3 exists to prevent. |
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
7. ⭐ **The SEED check (D0.1), and it is the one that outlives this campaign**:
   pick any catalog entry and ask — *could an engine populate its fields without
   minting a single new type?* If no, the embryo has drifted from being a seed
   and D0.1 is violated, **however clean the interface looks**. `[R]` the cheap
   form: the entry's fields are §III.1's outputs (`p₁…p_k`, `I`, `P`, `det P`,
   chart, pushforward, stratum) and its `derived_by` field admits `"engine"`
   as a value nothing has to be rewritten to produce.
8. ⭐ **`Space = str` is retired** — `grep` finds no `support=` taking a bare
   string literal, and `Manifold.contains` refuses the Part I forged measure at
   construction rather than three hops downstream.

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

⭐ **2026-09-02, at the fused step's opener (§V.5l): a FIFTH item PRECEDES the four and is NOT fused.** Tracker **2.5** makes the angular moment space READ off the frame while the SH basis is still bound everywhere — bit-identical alone, so it lands green as its own commit; without it the flat basis mismatches every operator end and moment field minted from `L` (`[M]` 8 production mints). The four above remain ONE commit after it, and the Γ-slot (**2.2b**) is a sixth, separate, after them.

## XII.2 ⏏ The MINIMUM chain to return — the critical path

| # | why it is on the path |
|---|---|
| **Phase 0** (all 7) | ⚡ unconditional; 0.1/0.2 remove the fabrication the defect rides on |
| **Phase 1.1** — *the `SO(2)/S²` entry only* | 3.4's probe needs that one catalog row. §III.2 already derives it; the other entries serve Phases 5–6 and are **NOT on the path** |
| **Phase 2.1–2.5** | `Basis.domain`, G0, the typed `Chart`, the slab rule's declared quotient group — the objects 3.4 reads — and, since 2026-09-02, **2.5**: the angular moment space READ off the frame (§V.5l), which must land BEFORE the fused commit |
| **Phase 3.1** — *scoped to the `SO(2)` entry* | ~~the catalog home + the probe~~ ⛔ the home landed at 2.0a/1.1/2.4 and the probe is 3.4's (§V.5j); ✅ **LANDED 2026-09-02** as the entry's quotient map + pushforward reference — what 3.4b's `Descent` pulls back along, and what its claim is against |
| **Phase 3.4** | ⭐ **the fix itself** — the trivial isotypic component |
| **Phase 4.1** | the retrofit, under the CORRECTED gate (Part VI) |
| **Phase 4.9** | re-key the 12 collateral sites, incl. the green test that pins the defect |

**Explicitly NOT on the path** (may land in any later session): Phase 1.2–1.8
except 1.1's one entry, Phase 3.2/3.3, Phase 4.2–4.8, and the whole of Phases
5, 6, 7. `[M]` §VIII.4 is what licenses this: today the defect has **one**
consumer, and CP/MoC/MC/diffusion/kinetics/fuel/TH carry **0 angular-moment
lines**, so nothing else is waiting on the layer.

## XII.3 ⏏ The exit predicate — all six, measured, before returning

1. ✅ **MET 2026-09-02** — `scratch/_pl_slab_defect_repro.py` LEG A reads **`+4.000000000000`** at
   `L = 0, 1, 2` (`[M]` rel 1.2e-15 / 4.4e-16 / 6.7e-16) and the gate's L = 3 row passes.
2. ✅ **MET** — the three strict-xfail rows XPASSed and their markers are deleted; `catches("ERR-080")` kept and mutation-verified (16 of 170 rows red under the re-forged defect, 0 of the 3-D rows).
3. ✅ **MET** — `[M]` ≤ 4.3e-16 on 49 of 52 (rule, L) rows (13 rules × L ∈ 1..4); the three exceptions carry their verdict: `gauss_legendre(2)` at L ≥ 4 (the GL dead-slot theorem, pinned), `product(4,4)` and `folded_product(2,4)` at L ≥ 4 (the two-pole aliasing under the harmonics, PRE-EXISTING) — *up to four rows, not two* (the memo's H-9), routed to Phase 4.4.
4. ✅ **MET** — `[M]` 135 arrays (9 sphere rules × L ∈ 0..4: table, Gram, analysis) `array_equal` to HEAD via a detached-worktree capture; published in the docs as a THEOREM (same basis, same measure object, the identity arrow).
5. The **full 13-tree gate green**, with per-tree deltas PREDICTED before the
   run (`reference_test_execution_env`). ✅ **BASELINE MEASURED 2026-09-02 at
   `2f294ef1` (code = `17e0eb0e`): 10329 / 19 sk / 69 xf / 0 failed = 10417
   collected, 13 of 13 rc=0, every tree equal to its predicted denominator**
   (Part XIV's baseline table). The predicate for the fused commit's re-gate is
   the DELTA against this row, predicted first. *(Kept per §3: before this
   date the baseline was UNMEASURED — step 4's `10106 / 0 / 19 sk / 66 xf`
   was the last full gate, `10116` predicted-not-run at CS4c §18.7.)*
6. ✅ **MET** — `sphinx -E -W` 0/0/0 (twice: after the archivist, after the accessor retirement); `dead_references` 0 / 53 after the API page's stale `:meth:` was re-pointed.

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
| 0.1a | 3-D rules hand the frame their OWN measure — ⚠ keystone is a **ROUTE** gate (`frame.measure is q.measure`), not bit-identity (§II.15 R4) | ⏏ yes | ✅ `2c1a06b1` — `Quadrature._harmonic_frame_measure`; `[M]` route gate **0 of 12 → 10 of 12**, nodes/weights bit-identical on all 10, and the forgery's **other two** losses reversed (§II.8's correction). Gates Q8.1–Q8.5, per-arm mutation battery 5/5 caught with a call counter |
| 0.1b | the 1-D rows — ⛔ **rides 3.4**, the lift does not exist (§II.9) | ⏏ yes | ✅ **the fused commit** (2026-09-02) — `_harmonic_frame_measure` DELETED, `angular_frame` binds the rule's own measure on 12 of 12; `test_q8_4` rewritten as the 1-D ROUTE gate |
| 0.1c | the fold's `S^2/sigma_y` tag — ⛔ §8 hazard REFUTED (0 reds; 0/145 089 checks); was blocked on **2.1-W** | ⏏ yes | ✅ `2c1a06b1` (with 0.1a) — `[M]` `folded_product` frames move `L2[S^2]` → `L2[S^2/sigma_y]`; gated by `test_q8_3` with an unfolded discriminating control |
| 0.2 | phantom component unspellable — ⛔ means REFUTED by the census; ⛔ *"runs AFTER 0.1"* also refuted (0.1a removes **0 %** of the phantom reads — only 1-D rules can phantom-read) | ⏏ yes | ✅ `ce46181c` — `axis_cosines` strict (the coordinate) + `mean_axis_cosine` (the orbit mean, zero DERIVED). ⭐ Found a **third** meaning in those zeros (`i >= 3` is no axis at all, not a suppressed one), single-sourced `spherical_harmonics` onto the frame (0 production consumers), and collapsed `_quadrature_axis`'s 3-arm ladder + its defaulted-`getattr` hazard. Gates Q9.1–Q9.5 + a re-scoped Q8.6; 5-arm mutation battery, every arm a distinct red set |
| 0.2-census | full-suite phantom census (moved here from 2.4) | ⏏ yes | ✅ **DISCHARGED** — §II.14, 13 trees rc=0, 11 sites / 1604 reads |
| 0.3 | `metric.py` "noise mode" + rcond re-derivation | ⏏ yes | ✅ `ae4dbc1f` |
| 0.4 | `directional.py` `spherical_harmonics` docstring | ⏏ yes | ✅ `ae4dbc1f` |
| 0.7b | ⭐ NEW — gate the SECOND symptom (the `DenseMetric` RAISE) | ⏏ yes | ✅ `ae4dbc1f` |
| 0.5 | `spherical_harmonics.rst` rank-deficiency correction | ⏏ yes | ✅ `ae4dbc1f` |
| 0.6 | `_evaluate_real_sh` raises on non-unit directions | ⏏ yes | ✅ **the fused commit** — `SphericalHarmonicBasis.evaluate` refuses off-sphere directions (a real `raise`, proven `-O`-live by a subprocess row); `_PARITY_PROBE_DIRECTIONS` retired (they were off the sphere) |
| 0.7 | strict-xfail reproducer gate + **ERR-080** minted | ⏏ yes | ✅ `ae4dbc1f` |
| 1.1 | catalog derivations — the `S²/σ_y` entry (the shipped fold), which §V.5d(a) showed PRECEDES 2.0c | ⏏ **yes** | ✅ **σ_y LANDED** — `_sphere_mod_mirror` (one derivation, 3 mirror keys), `Ball`, `FundamentalDomain`, `Quotient`'s two slots + construction invariant, `singular_stratum` retyped. `[M]` 4-leg acceptance on production nodes; 6-arm mutation battery, each arm its own gate; pyright 0. ⭐ The user's two-slot ruling is in §V.5d. **Remaining under 1.1:** the other catalogue entries (not on the exit path) |
| 1.2 | Gelfand-pair reading of Funk–Hecke | no | ☐ |
| 1.3 | connection-term verify-or-kill | no | ☐ |
| 1.4 | record the four gates | no | ☐ |
| 1.5 | record the placement ruling | no | ☐ |
| 1.6 | `ℛ* = ℛ` as adjoint-vs-dual | no | ☐ |
| 1.7 | the support algebra, written | no | ☐ |
| 1.8 | vocabulary reconciliation | no | ☐ |
| 1.9 | ⭐ **NEW 2026-09-02 (user-ruled after 2.5) — a basis states the FULL group its functions have: `SubgroupOfO3.O2(axis)`**, so `P_ℓ(Ω·ê_a)` declares `O(2)_a` (the stabiliser fact of 3.1: `O(2)_a` and `SO(2)_a` induce the same orbit partition) and the lattice G0 ADMITS a Legendre basis on a σ_b-folded rule (b ≠ a). Memo H-6. Its own commit, FIRST after the compaction that follows the fused commit | no | ✅ **LANDED 2026-09-02 `c1d53206`** — one entry named by its stabiliser, the refusal on the TYPE and at the door (`orbit_stabiliser`, the elegance round's one absence); read the *1.9 EXECUTED* section; the deferred accessor design is #433 |
| 2.0a | ⭐ **MINT `Manifold`** (D0.7, user-ruled) — retires `Space = str` | ⏏ yes | ✅ `b8c05d16` — 9 variants, 40 tests, pyright 0. ⚠ the TYPE only; `support` is still `str` (that is 2.0c/2.1) |
| 2.0b | `Manifold.contains` — the membership predicate; refuses the forged measure at construction | ⏏ yes | ◐ **HALF** — the predicate ships and is gated both legs (`[M]` refuses the `gauss_legendre(8)` forgery, norms `[0.1834, 0.9603]`; admits the normalised control) @ `b8c05d16`. The **refusal AT CONSTRUCTION** is unbuilt — that is the wiring, and it rides 0.1. ⛔ no `catches("ERR-080")` marker until then: refusing a forged ARRAY is not the production path refusing it |
| 2.0c | `DiscreteMeasure.support: str → Manifold`; the `L2[...]` name derives; `quotient_group` becomes a derived property | ⏏ yes | ✅ `025834f5` — all **six** implementors retyped in one commit (the row named one; §V.5f(a)). `phase` dispatches on the manifold's TYPE; the LIVE `folded_product` defect fixed; `pushforward`'s fabricated default retired; 2.1's handoff discharged **and witnessed** (⛔ the repair measured BLIND first). `[M]` 5 trees rc=0, +7 tests, matrix reconciling unit-for-unit. ~~▶ NEXT.~~ ~~⛔ BLOCKED on 1.1's σ_y entry (§V.5d(a))~~ ✅ **REMEDIED 2026-09-01** — 1.1 landed `b55bba56`; `[M]` `SPHERE.quotient(SubgroupOfO3.Mirror("y"))` → `'S^2/sigma_y'` (§V.5e(a)). ⭐ **Read §V.5e before designing** — the second opener re-shaped the step: the row's sizing counts PRODUCERS (`[M]` §6b ≈ 87 producers + **16 consumers**, and the 16 are where a retype fails QUIETLY), the prize is that `phase` is stringly-typed dispatch, and it carries a LIVE defect (`folded_product(4,8).measure.phase` **RAISES**). Re-scoped at the opener: **absorbs 2.0d** (§V.5d(c)); its `indicator_basis` clause **moves to 2.1** (§V.5d(e), §6b). `[M]` the retype preserves every production space name bit-identically bar two (§V.5d(f2)) |
| 2.0d | `measure.quotient_group` | ⏏ yes | ✅ **LANDED INSIDE 2.0c** `025834f5` as a derived property (`Quotient.by`), with its witness. ⛔ **DISSOLVED as a separate item** (§V.5d(c)) — `[M]` `Quotient.by` IS the group, so a stored field is a second home for it, and `support`'s four forwarding rules would have to be hand-kept in agreement. Its done-when was also unreachable as written (§V.5d(d)): the `SO2` answer is **2.4**'s declaration, not a derivation |
| 2.1 | `Basis.domain: Manifold` (⛔ `support` rename REFUTED — §III.10) — `IndicatorBasis` takes a ctor field. ⭐ **Now also owns 2.0c's `indicator_basis` clause** (§V.5d(e)): `[M]` the false `L2[coarse_cells_R1]` name makes a 2-group ENERGY space and a 2-cell SPATIAL space `==`-EQUAL and hash-equal, with a passing negative control — an illegal state that IS representable | ⏏ yes | ✅ `c461fe8d` — `Basis.domain` abstract; `IndicatorBasis` takes `partition_of`; `ambient_dim` promoted public as its first cross-module consumer. ⭐ Keystone is `test_d6` — **a frame's two halves name ONE manifold** — not a name assertion, because the pre-2.1 name was self-consistent. `[M]` §6b set was **23**, not 18 (3 `OverlapBasis(` + `from_indicator` + a test-local subclass); 13 gates; 7-arm mutation battery. ⛔ the *"26 space-comparing sites"* figure was never re-measured and is **struck** — the repair is at the 5 producers, not at the comparers |
| 2.1b | `Basis.invariance_group` — G0's other side; ~~today NEITHER side exists~~ both sides real (2.4 / 2.1b) | ⏏ yes | ✅ `9b4a4d9c` |
| 2.1-W | ⭐ **PRE-CARVE**: the fold-basis witness | ⏏ yes | ✅ `TestFoldedCylinderP1BindsTheQuotientBasis` — 3 rows, L1, `verifies("pn-scatter", "discrete-measure-quotient")`. Keystone is k_inf-invariance on a HOMOGENEOUS folded cylinder, not the 2g-het solve the row proposed (§II.16). `[M]` rebind reds all 3 at nine orders of separation. **0.1c is unblocked** |
| 2.2 | G0 at frame construction — ~~⚠ ALSO owns the `AngularSymmetry` Γ-slot (§II.10, inventory #11)~~ ⛔ the Γ-slot is **2.2b** since 2026-09-02 (user-ruled, its own commit after the fused one). ✅ **RULED 2026-09-02 (§V.5l(c)3): ONE predicate — the descent arrow `measure.support → basis.domain` EXISTS (the lattice; equality is `K = H`), the table pulled back along it, on `FrameBase`; and `P_ℓ` on a full-sphere rule BUILT and witnessed** | ⏏ yes | ✅ **the fused commit** — G0 = the descent arrow (`quotient_onto` + `_descent_arrow` at both construction paths, the table pulled back along it); the seven-row table gated; ⭐ it CAUGHT `OverlapBasis`'s misdeclared domain on landing. The Γ-slot is 2.2b |
| 2.3 | the typed `Chart` ⛔ ruled **`ManifoldMap`**; pushforward derives support | ⏏ yes | ✅ landed 2026-09-02 `5ec3a00a` — read the *2.3 EXECUTED* section |
| 2.4 | slab declares its quotient group (`SPHERE.quotient(SO2("x"))`); the SO2 axis-convention ruling. ~~full-suite phantom census~~ ⛔ **that MOVED to Phase 0 and is DISCHARGED** (`ee60010e`+`ccda0e61`) — it was 0.2's denominator, not 2.4's | ⏏ yes | ✅ `17501245` — read the *2.4 EXECUTED* section. **User-ruled:** `SO2(axis)` parameterised like `Mirror` (the tree carries TWO poles — slab/SH about x, product/`C_n`/`D_nh` about z — so no fiat could be right); `gauss_legendre_on_polar_orbit(n, axis)` declares, `gauss_legendre_on_mu` stays on the chart as the product factor. `[M]` the `[-1,1]` angular/spatial collision is unspellable (`L2[S^2/SO2_x]` ≠ `L2[[-1,1]]`, hashes differ); Part IV obstacle 1 answered by DERIVATION (`SO2("x")` True / `SO2("z")` False); stage 0 refuses the chart-level rule AND the wrong axis; 4 trees rc=0, +22 reconciling unit-for-unit; pyright 0. ⛔ the pre-flight's *23 sites / 3 literals* were the wrong instruments — the set was the **5 test-local factories** the stricter gate refused, found by the red loop. ⚠ Also removed a pre-existing 9 s hot spot the extra candidates amplified (§(e)) |
| 3.1 | ~~`numerics/symmetry/` catalog~~ ⛔ the catalogue HOME landed at 2.0a/1.1/2.4 (`manifold._ORBIT_CATALOGUE`, §V.5j(a)); 3.1 = the entry's quotient map + pushforward reference (⏏ scoped to `SO(2)`) | ⏏ partial | ✅ landed 2026-09-02 — read the *3.1 EXECUTED* section (`Quotient.orbit_coordinates` + derived `quotient_map`, `Quotient.reference`; `AngularSymmetry.reference` reads the entry) |
| 2.5 | ⭐ **NEW 2026-09-02 (the fused step's opener, §V.5l) — the angular moment space is READ off the frame, never minted from `L`.** `[M]` 8 production `SphericalHarmonicSpace.from_L(L)` mints build operator ends and the moment-flux field's space from `L` alone; `HarmonicFrame` narrows to `SphericalHarmonicBasis` at two doors; the `(name, shape)` guards compare the mints to `frame.basis_space`. The Stage-2 generator ruling's own defect, found by re-censusing 3.4-R's evidence. Bit-identical alone (the SH basis is still everywhere); PRECEDES the fused commit (§6b: the flat basis mismatches every retained mint) | ⏏ yes | ✅ `58ad7c28` — read the *2.5 EXECUTED* section: 7 producers + 2 doors read the frame's basis; `TruncatedBasis`; A-R1 (`basis.space`, not the dressed one) and A-R2 (re-derive through the mesh's quadrature) ruled and gated; slab L=0/1/2 `array_equal`; 11-arm battery, none blind |
| 2.2b | ⭐ **NEW 2026-09-02 — the Γ-slot, SPLIT off 2.2** (user-ruled): stage 0 admits a rule on `S²/H` for `H ⊆ Γ` (lattice, derived); stage 1 asks Γ-invariance THROUGH the entry's `quotient_map`. `[M]` today `cylinder.admits_domain(folded)` False AND `admits_symmetry` False; `select_quadrature` 0 production callers | ⏏ yes | ✅ **LANDED 2026-09-02 `<HASH>`** — stage 0 through the descent arrow + what it SPENDS; the owed symmetry asked ON the orbit space through `Quotient.induced_action` (one closure for `is_invariant`, `orbit_certificate` and `ordinate_permutation`); read the *2.2b EXECUTED* section |
| 3.2 | `ReynoldsProjection` | no | ☐ |
| 3.3 | `OrbitAxis` + retraction/section minting | no | ☐ |
| 3.4 | ⭐ **trivial isotypic sub-basis — THE FIX** | ⏏ yes | ✅ **the fused commit** (2026-09-02) — `LegendreBasis` + `LegendreSpace` + `MomentHead`; the entry's probe; `_harmonic_basis` dispatches on the SUPPORT; LEG A `+4.000000000000` at L = 0..2, the ERR-080 gate XPASSed, its markers deleted. Read the *FUSED EXECUTED* section |
| 3.4-R | ⚠ the encoding ruling | ⏏ yes | ✅ **RULED 2026-08-31 (user): (b) true `LegendreBasis`** — see below. ⛔ its EVIDENCE (*15 sites, ONE production*) was re-censused 2026-09-02 and found one layer too narrow (§V.5l(b)); **RE-AFFIRMED by the user at the opener, with 2.5 minted in front of it** |
| 3.4-R2 | ⚠ the DESCENT ruling | no | ✅ **RULED 2026-08-31 (user-endorsed): (i) mint `Descent`** — Part V.5 §C |
| 3.4b | ⭐ **`Descent`** — the iso witness; lands WITH or BEFORE 3.4 | no | ✅ **the fused commit** — `numerics/basis/descent.py`: the two realizations, the isomorphism at the BIT (7/7 sphere rules), the discriminator as `frame_basis` |
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

**Open rulings blocking the EXIT PATH:** ⭐ **NONE.** (One ruling is open OFF the exit path: **3.4-R2**, the descent isomorphism — Part V.5 §C. It must be ruled before 3.4 lands, but it does not block the return to Campaign 2.)

**Ruled — 2.0a-R, the `Manifold` SHAPE (user, 2026-08-31).** A **closed sum
type, split by TOTALITY**: `dim` / `name` / `contains` / `__mul__` are total, so
they live on an abstract base; the derivation fields are partial, so they live on
`Quotient` alone. Decided against a polymorphic hierarchy on measured evidence,
not taste — `[M]` D0.1's derivation fields belong to `Quotient` alone, so a
hierarchy must put them on the base returning `None` for every non-quotient, which
is the exact tax `SubgroupOfO3.mirror_axis` already pays and `directional.py:522`
already branches on. Both shipped precedents were read rather than recalled
(`plan-authoring` §1): `BoundaryTraceLaw` is 7 registered sibling classes;
`SubgroupOfO3` is ONE class over a `_tag` sum, dispatching with `[M]` **10
`isinstance` calls in 13 methods, 0 `match`** — so the ruled shape is that
precedent's own data model with the dispatch it would get today, keeping the two
halves of `M/H` structurally parallel.

⛔ **Corrected 2026-08-31 (docs pass): the shipped seed is `6 of 8`, not the 8
this ruling first claimed.** `[M]` `dataclasses.fields(Quotient)` carries
`generators`, `syzygy`, `gram`, `det_gram`, `singular_stratum`, `derived_by`.
D0.1's list also names the **chart** and the **pushforward measure**; the chart
ships only as its CODOMAIN (`realization`) and the pushforward not at all.
⟹ this does **not** violate D0.1's falsifiable test (*could an engine populate
the fields without introducing a new type?*) — but only because `Chart` is
already a SCHEDULED mint (V.5 A #2, tracker **2.3**). The seed is complete *up
to* `Chart`, and 2.3 completes it; writing the count as 8 read as *done*.

⛔ **REFUTED 2026-09-01 — it is NOT a twin, and the correction matters for
2.2.** `[M]` `AngularSymmetry.support` computes :math:`S^2/G^0`, the
**CONTINUOUS isotropy**, while a mirror lives in the *discrete residual* `Γ` —
so `σ_y` is a row that property **structurally cannot answer**, not a row it
answers in another vocabulary. The two do not compute the same function, so
`Manifold.quotient` does not absorb it; **2.2 must give the ontology its second
slot**, which §II.10 already says and this now explains. ⭐ The corroborating
break, re-measured independently: `[M]`
`GEOMETRY_ANGULAR_SYMMETRY["cylinder"].admits_domain(folded_product(4,8).measure)`
is **`False`** (`'S^2'` vs `'S^2/sigma_y'`) — stage 0 refusing the shipped
cylinder, inert today only because `folded_product` is not a registered spec.
⚠ What SURVIVES of the original claim below: the `Trivial` row really is the
same answer re-derived, and the mint really is a re-typing for the continuous
half. *(Original text kept per §3:)*

⭐ **An un-briefed Pattern-2 twin, found by the docs pass:
`AngularSymmetry.support` (`quadrature/registry.py:869`) is ALREADY an
orbit-space catalogue** — in the string vocabulary, with the same lookup
(`SO2 → "[-1,1]"`, `Trivial → "S^2"`) and the same *"extend it when a geometry
first spends it"* refusal. Strong evidence the mint is a **re-typing** rather
than an invention, and it is the twin tracker **2.2** must absorb. `[M]`
comparing the two exposed a real gap in mine: `S²/{e}` is legal and trivially
derivable and my catalogue RAISED. Fixed by deriving it — `M/{e} ≅ M` is a
theorem, not a table row, and the same procedure on the trivial group gives
`P = I`, `det P = 1` vanishing nowhere, hence no stratum, hence a free action
(right vacuously, the only element being the identity). It doubles as a
positive control on the machinery, and a committed row now pins the two
implementations in agreement with a ⚠ *do not re-baseline* note.

⛔ **And a §6d finding that inverted mid-check.** The placement note says
*"its own module with no intra-`numerics` imports"*. I questioned it, ran an AST
import census, and got *"`symmetry` imports only `geometry.transformation`, so
`symmetry → measure` = 0, no cycle"* — **wrong, and in the reassuring
direction.** My filter tested `node.module.startswith("orpheus")`, and an
`ImportFrom` with `level > 0` carries an UNQUALIFIED `.module`, so every
**relative** import was silently dropped — including
`symmetry.py:98`'s `from .measure import DiscreteMeasure`, a module-scope
RUNTIME edge. `[M]` re-run with relative resolution and a positive control:
`symmetry → measure` is real, and `measure`'s own `symmetry` import is
`TYPE_CHECKING` at `:91` plus a **function-scope** one at `:1005` — deferred for
exactly this reason. So `measure → manifold → symmetry → measure` is a genuine
cycle and the placement note was right, for a reason it never wrote down.
`Manifold` therefore references `SubgroupOfO3` under `TYPE_CHECKING` only;
verified by importing the trio in both orders and in a fresh interpreter.

**Ruled:** 3.4-R was ruled by the user on
2026-08-31 — **option (b), the true `LegendreBasis` on `[-1,1]`**, on the
measured evidence in Part VI (15 slot-indexing sites, only ONE of them
production; and the project's own type-minting criterion selects it, since
`{Y_ℓ^m}` and `{P_ℓ}` are non-isomorphic and a non-identity morphism — the
quotient projection — applies). The 0.1c fold tag was ruled the same day:
**measure the blast radius, then land it in Phase 0** (§II.10 gives the reason
this is right: the forgery is what CONCEALS the missing discrete-quotient slot).

---

# Part XIV — ⏸ COMPACTION **RECORD** (2026-09-02, 4th of the day — the FUSED commit LANDED, its follow-up owed) + ▶ Resume surface

⚠ **Naming, because it has been ambiguous since the plan was written.** Part VI
schedules two *planned* pauses — "COMPACTION POINT 1 — after Phase 3" and
"COMPACTION POINT 2 — after Phase 5". **This Part is neither.** It is the
running RECORD of where the campaign actually is, rewritten at each real
compaction. It has been rewritten eleven times (2026-08-31, six times on
2026-09-01, four times on 2026-09-02); do not try to match its number against Part VI's.

**Read in this order on pick-up:** this Part → Part XIII (the tracker) → Part XII
(the exit gate, esp. **XII.1b**) → the phase you are executing. **Part IX,
§II.15, §II.16, §V.5d, §V.5e and the 2.1 / 0.2 execution sections first if you
are about to re-propose anything** — the supersession table below is the count
that is kept current (`[M]` **60 rows, 56 of them ⛔**, counted 2026-09-02
at the fused step's pre-flight, `2f294ef1` (was 58 / 54 at 3.1's; 3.1 added exactly 2) — predicate: every `| `-led line of that table except the
header and separator, ⛔ anywhere in the status column; re-derived by script,
never incremented by hand. ⚠ The previous record read *48 / 47* under an
UNSTATED predicate and 2.1b added exactly 2 rows, so the two predicates
differ by ONE row — §2's FILTER clause, caught by writing the predicate
down), and
most of those are mine, several written the same day they were refuted, one
refuted *by my own fix*, and one a denominator I "corrected" in the WRONG
direction.

⭐ **The single most useful habit this campaign has, stated so it survives:**
every phase begins by RE-MEASURING its own section's premises before designing.
`[M]` **NINETEEN consecutive phase openers have corrected their own section**, and
two of them (P6, and 2.0a's sizing) dissolved or re-scoped the phase entirely.
The eleventh (2.0c's third opener, §V.5f) found that the step's own title named
**one implementor of a six-implementor contract**, one of them spelled so that
no census of the alias could return it; the twelfth (2.4's pre-flight, §V.5g)
found that the step filed as *wiring* repairs a **live** defect, and is the one
step in the campaign that RENAMES a shipped space; the thirteenth (2.1b + 2.3's
pre-flight, **§V.5h**, run at this compaction) found Part IV's G0 argument VOID
under the post-2.4 vocabulary and half of 2.3 possibly consumer-less; the
fourteenth (2.1b's own opener, at execution) found the tracker's MEANS — *a
field answered by six subclasses* — already sitting in the fold basis's
`domain.by`, so the step landed as ONE derived property and zero subclass
edits, the way 2.0d dissolved into 2.0c; the fifteenth (2.3's pre-flight,
**§V.5i**, run at this compaction) confirmed by census that a SECTION has no
production consumer, found the barycentre map spelled TWICE (once honestly,
once as ERR-080's forgery), and scoped a done-when that three honest
tabulations would have pinned red forever; the sixteenth (2.3's own opener,
at execution) is the first whose pre-flight survived intact — (a)–(f)
reproduced 6 of 6 — and it still corrected the row's NAME (*"the typed
`Chart`"* → `ManifoldMap`, user-ruled: a chart is `M ⊃ U → ℝⁿ`) and one
symbol's HOME (the axis table), and then had one of its OWN claims refuted by
the archivist on landing (the `Ball` census, §V.5i(b)); the seventeenth
(3.1's pre-flight, **§V.5j**, run at this compaction) found the row's
"catalog home" half already LANDED and its package re-home OFF the exit path,
leaving 3.1 two fields and one map; the eighteenth (3.1's own opener, at
execution, 2026-09-02) is the second whose pre-flight survived intact —
(a)–(e) reproduced 5 of 5 at `3623adc2` — and it still took three user
rulings (the map's NAME, its MECHANISM, the two honest `None`s) before a
line was written; the nineteenth (the fused step's own opener, at
execution, 2026-09-02, **§V.5l**) reproduced §V.5k 6 of 6 and still found
the RULING's evidence one layer too narrow — the rectangular layout is the
transport layer's moment-space CONTRACT, minted from `L` at 8 production
sites behind two `isinstance` doors the plan never mentioned — so the user
minted a bit-identical pre-step (2.5) in front of the fix and widened G0 to
the lattice with `P_ℓ` on a full sphere built now. Budget the opener; it has
never once been wasted.

⚠ This supersedes the two earlier versions of this RECORD. Their content is
folded in below; nothing from either was dropped.

---

## ▶ RESUMES AT — stated as OUTCOMES (`plan-authoring` §1)

✅ **1.9 (#432) LANDED 2026-09-02 as `c1d53206`** — the commit titled `feat(numerics): an orbit
space is named by its stabiliser — SubgroupOfO3.O2(axis)`, its hash stamped by the plan-only
commit that followed. **Nothing is owed on 1.9.** The record is the *1.9 EXECUTED* section in
Part XIII; the short form, in the order it happened:

1. **The rulings** (user, at the opener): `S²/SO(2)_a` and `S²/O(2)_a` are ONE catalogue entry,
   named by its STABILISER `O2(a)` (`[M]` the invariant rings coincide);
   `SPHERE.quotient(SubgroupOfO3.SO2(a))` REFUSES with the theorem — the more principled option
   regardless of the re-spelling it cost (now a memory rule: effort is never the tie-break); the
   registry's slab/sphere `continuous_isotropy` is `O2("x")`.
2. **The elegance round found the step's ONE absence** — the tree had no *stabiliser* concept
   on the group, so the theorem was spelled three times and the refusal lived in a builder
   (`[M]` `replace(entry, by=SO2('x'))` forged a second spelling and `barycentre` accepted
   it). Landed in the same commit: `SubgroupOfO3.orbit_stabiliser` (`SO2(a) → O2(a)`,
   `SO3 → O3`, every other member itself, 24 of 24 gated); `Quotient.__post_init__` and the
   catalogue door both refuse a non-stabiliser group; the three `SO2_a` decoy keys are gone;
   `_fixes_axis` spells `SO(2)_a = O(2)_a ∩ SO(3)` by composition with ONE absolute tolerance;
   `D_1h` is order 4; `AXIS_LETTER`. ⛔ The proposed `rotation_axis → fixed_axis` rename was
   REFUSED (`C_n` fixes ẑ pointwise and answers `None`; the precise term is the principal
   axis, a wider accessor) — **#433**.
3. **The gate**: **13 of 13 rc=0, 10604 = predicted** — the fused row's 10556 + 48, all of it numerics (2809 → 2857: the 43 test-architect rows + my 5), every other tree unit-identical to the fused row, the harness 424 unchanged (the archivist's two new equation labels added no parametrized row this time; predicted from collect counts BEFORE the run, phase 1 of `scratch/_p19_full_gate.sh`). Wall **60 min 27 s** on the session's machine (`scratch/_p19_full_gate.log`, per-tree `_p19_gate_<tree>.log`), started at `2ea2da73` with the 52 files dirty.
4. **Two of my own numbers were refuted by the agents and corrected in place before shipping**:
   the walk cost (*"+3 candidates ≤ 5 %"* was ONE warm-cache draw; `[M]` +11.3–26.2 % by the
   min of 15 interleaved repeats — the `plan-authoring` surprise log carries the row: a timing
   is a stochastic measurement whose configuration is the repeat protocol) and `D_1h = {e, σ_h}`
   (order 4, `{e, C_2^x, σ_y, σ_z}`).
5. **Gates and docs**: +48 rows (43 test-architect, 5 mine); the 19-arm in-process battery all
   biting on a 0/0 baseline, two of its arms the FIRST witnesses of `O(2)_a ⊄ SO(3)` and of the
   space name being READ off the domain, one arm theorem-inert and gated as a named blindness
   row; 9 theory pages + a new `manifolds.rst` chapter on the stabiliser law, 7 present-tense-
   false claims fixed; sphinx `-E -W` 0 before and after; `dead_references` 0 / 52.

✅ **2.2b (the Γ-slot) LANDED 2026-09-02 as `<HASH>`** — the commit titled `feat(numerics): a
geometry admits a fold of its domain, and asks the owed symmetry ON the orbit space`, its hash
stamped by the plan-only commit that followed. **Nothing is owed on 2.2b.** The record is the
*2.2b EXECUTED* section in Part XIII; the short form:

1. **The rulings** (user, three, all the recommended option): the induced action of an
   isometry on an orbit space is MACHINERY on `Quotient`, refusing a motion that does not
   normalise the quotienting group, and `is_invariant` asks a quotient-supported measure
   through it — one notion of invariance; stage 0 reads the descent arrow `quotient_onto`
   plus what it SPENDS in Γ; the axial entries take the same route with the barycentre as the
   equivariant lift.
2. **The algebra the ruling did not spell**: `g` descends to `M/H` iff it normalises `H`; for
   a continuous group the normaliser is EXACT through the rotation generator `[ê_a]_×` and one
   representative per coset, never sampled (ERR-072); a continuous group's identity component
   must FIX every node of a finite set on the orbit space (connected orbits) — the
   axis-support / origin rule on the section nodes, run only where the quotienting group is
   finite. `[M]` the exact normaliser against brute conjugation: 529 pairs, 138 motions, 0
   disagreements.
3. **What shipped**: `Quotient.lift` / `ambient_representatives` / `induced_action`;
   `SubgroupOfO3.is_normalised_by` / `normalises` / `identity_component`; ONE closure
   `_orbit_space_closure` read by `is_invariant`, `orbit_certificate` and the new
   `induced_permutation` (which `Quadrature.ordinate_permutation` now reads — on a fold σ_y is
   the identity permutation of the stored ordinates, where the ambient reading found none and
   a reflecting face would have been refused); `spent_group(D, X)` beside `quotient_onto`;
   `admits_domain` one expression; the spent-group door with `{e}` as its named exception; the
   ambient container `R³/{e}`; `_invariance_on_points` and `_polar_axis_of` RETIRED;
   `permutation_between` / `permutation_preserving` factored out of `RigidMotion`; II.11's
   1-D shape refusal gone as a by-product.
4. **Measured**: cylinder and cartesian2d ADMIT the folded rule at both stages (both False
   before); the kernel swap moved 3 of 230 (group × measure) answers, all on the fold, all
   intended; the slab and product families bit-identical; the selector UNMOVED on 16 of 16
   rows (stage 2's V conjunct refuses a fold first — the Γ leg is load-bearing at
   `admits_domain`, 4 of 28 rows). 148 gate rows (14 classes) from the test-architect, whose
   in-process call census read the §6b set (1636 `is_invariant` calls, 1 verdict moves); the
   14-arm battery: 12 bite, the bijectivity leg UNINSTALLABLE (the `Permutation` type
   refuses), the containment short circuits theorems (0 reds with a real mutant).
5. **The gate**: **13 of 13 rc=0, 10752 = predicted** — the 1.9 row's 10604 + 148, all of it numerics (2857 → 3005: the test-architect's 14 classes; my own test edits add no row), every other tree unit-identical to the 1.9 row, the harness 424 unchanged (the archivist's three new equation labels added no parametrized row); predicted from collect counts BEFORE the run (phase 1 of `scratch/_22b_full_gate.sh`). Wall **59 min 27 s** on the session's machine with Sphinx building alongside for the first minutes (`scratch/_22b_full_gate.log`, per-tree `_22b_gate_<tree>.log`), started at `4b7d24c3` with the 26 files dirty. The sn tree's 3439 is the witness that the moved `ordinate_permutation` changed nothing for the boundary realizer and the deck consumers on every shipped configuration.
6. **Refuted, corrected in place**: my memo's 2 flips per fold were 4; "a fold by σ_x of the
   σ_y-fold works today" was false; two claims in my delta to the archivist (translations
   "refused" — they act by their linear part; the container "the ball" — it is `R³`).
7. **Docs**: a new `manifolds.rst` chapter on who acts on an orbit space (+1232), stage 0/1
   re-stated, II.11 made whole; sphinx `-E -W` 0 → 0; `dead_references` 0 / 52. #370's gap 2
   is closed by this landing; its gap 1 stands.

The ⏏ ORDER table rules what follows; its first unlanded rows are **4.1** and **4.9**, then
the exit predicates' re-gate (XII.3). #426/#428 stay blocked behind #429.

✅ **THE FOLLOW-UP LANDED 2026-09-02 as `26151a8d`** — the commit titled
`docs(plans): the fused commit carries its landing hash and gate` (the hash
stamped by the compaction-prep commit that followed it). **The pre-emptive compaction's owed second
half is paid; nothing is owed on the FUSED landing.** What the five owed steps
found, in the order they ran:

1. **The gate** — **13 of 13 rc=0, 10556 = the prediction**, unit-for-unit on
   every tree; the +102 over 2.5 reconciled FILE BY FILE (the table in the
   *FUSED EXECUTED* section). ⚠ Two of the checklist's own guesses at the
   +102's terms below were wrong (`test_manifold` "+38" is **+11** by collect;
   the harness "+4" is `test_layer_imports` seeing four new production
   modules, not the matrix) — diff impressions, not counts.
2. **The typing narrowing** — applied. ⚠ The script's `ensure_import` was
   satisfied by a FUNCTION-LOCAL `MomentHead` import inside one test of
   `tests/sn/primitives/test_harmonic_moment_flux.py` and added no
   module-level one: **21 reds** (`NameError` — pytest REWRITES test-module
   asserts, so they run under `-O`; a bare `assert` is inert only in
   `orpheus/`), fixed by hoisting the import and dropping the local twin.
   The five files: **197 passed**. pyright over the fused commit's `.py`
   files: **33 → 31 diagnostics, 0 on any line the commit or the narrowing
   touched** (`[M]` every diagnostic intersected against `git diff -U0
   HEAD~1`; the 31 are pre-existing, on untouched lines).
3. **The record** — `5436184e` stamped at the four placeholder sites; the
   gate table in the *FUSED EXECUTED* section, the delta row in *Measured
   baselines*.
4. **#429** — the landing comment posted (`_fused_B_issue_comment.md`, filled).
5. The `_head_tree` worktree (the bit-identity oracle at `2ebdbc05`) removed.

▶ **NEXT — the user-ruled order:** ⏸ compaction, then **1.9** (#432,
`SubgroupOfO3.O2(axis)`, its own commit), then **2.2b** (the Γ-slot), then
4.1, 4.9, and the exit predicates' re-gate (XII.3: 1–4 and 6 met at
`5436184e`, 5 met by the gate above; the exit re-gate re-measures all six).
**Nothing on the exit path is blocked.**

*The checklist as it stood before the compaction — EXECUTED as recorded
above, kept per `plan-authoring` §3:*

⛔⛔ **READ THIS FIRST — a PRE-EMPTIVE compaction (the user asked for it
with the fused commit's gate mid-run), so the landing is SPLIT across two
commits and the second is OWED.** The FUSED commit (0.1b + 0.6 + 2.2-G0 +
3.4 + 3.4b — THE FIX, ERR-080 closed) was committed from a working tree
whose 13-tree gate was still running (`scratch/_fused_B_full_gate.sh` →
`scratch/_fused_B_full_gate.log`, predicted **10556** collected: numerics
2809 · transport 707 · geometry 732 · data 237 · homogeneous 50 · diffusion
113 · cp 141 · moc 121 · mc 41 · cross_method 81 · sn 3439 · derivations 1661
· root+harness 424; the run started ~11:35 and takes ~62 min). Everything
ELSE was green before the commit: the scoped red loop (14 trees' worth of
the affected scope, `_fused_B_red3.log`: 5809 passed / 7 fixture reds since
fixed, then 439 passed on the touched files + the condense consumers), the
19-arm battery (`_fused_B_battery2.log`), `sphinx -E -W` 0/0/0 (twice),
`dead_references` 0. **The follow-up commit, in this order, is the FIRST
thing a resumed session does:**

1. **Read the gate** — `grep -E '^tree=.*rc=|@@@' scratch/_fused_B_full_gate.log`.
   Expected: 13 of 13 rc=0, every tree equal to its predicted denominator
   above, total 10556 (= 10454 + 102: `test_legendre_basis` 32 +
   `test_descent` 20 + `test_manifold` +38 + `test_frame` +… − `test_q8_6`'s
   12 retired params + the matrix's harness rows). If a tree is red, its
   per-tree log is `scratch/_fused_B_gate_<tree>.log`; characterise before
   fixing (⚠ the gate ran on the tree as committed, so any red is real).
2. **Apply the typing-only narrowing on the test-architect's lines** —
   `.venv/bin/python scratch/_fused_B_pyright_tests.py` (24 pyright errors on
   NEW test lines: `space.factors[0].isotropic_slot` on a `FunctionSpace`-typed
   value, `metric.matrix` on an Optional, `S.frame.basis` as `TruncatedBasis`;
   behaviour-neutral, held back only because the gate was running — L37).
   Then re-run the five files it touches
   (`tests/numerics/test_frame.py tests/sn/operators/{test_n2n_operator,test_scattering_operator}.py
   tests/sn/primitives/test_harmonic_moment_flux.py tests/transport/fields/test_harmonic_moment_flux.py`)
   and `npx pyright` over `git diff --name-only HEAD~1 -- '*.py'` — expect 0
   on touched lines (the ~31 remaining diagnostics are pre-existing, on
   untouched lines — `[M]` intersected by script at commit time).
3. **Record the gate** in this Part's baseline section (the table pattern of
   the 2.5 row) and in the *FUSED EXECUTED* section's `<GATE — recorded in
   the follow-up>` placeholder; **stamp the fused commit's hash** at the three
   `(hash stamped in the follow-up)` sites [✅ stamped `5436184e`] (the section header, the tracker
   row, the ledger row) and the two ✅ pointers; commit as
   `docs(plans): the fused commit carries its landing hash and gate` with the
   two trailers; push.
4. **Post the #429 record** — `scratch/_fused_B_issue_comment.md` with
   `<HASH>` and `<GATE>` filled (`gh issue comment 429 -F`). ERR-080 is
   CLOSED in the catalogue; #429 stays open for 1.9, 2.2b, 4.1, 4.9, the
   exit predicates' re-gate.
5. Then ⏸ the compaction the user ruled — and after it **1.9** (#432,
   `SubgroupOfO3.O2(axis)`) then **2.2b** (the Γ-slot).

⭐ **The one-line state for a session that trusts nothing above:** `git log
--oneline -3` shows the fused commit; `git status --porcelain | grep -v
scratch/` is EMPTY after the follow-up and non-empty (the five typing-fixed
test files + this plan) before it. The scratch reproducers are listed under
*Artifacts* (`_fused_*`).


**Nothing on the exit path is blocked.** The σ_y orbit space landed, which was
the last item gating anything.

✅ **Phase 0 is DONE** — 0.1a + 0.1c @ `2c1a06b1`, 0.2 @ `ce46181c`. No
quadrature asserts a support its own nodes contradict, and the one surviving
fiction (the 1-D lift) is NAMED, local, and self-retiring.

✅ **A basis states the manifold its functions eat** — tracker **2.1** LANDED
`c461fe8d` (2026-09-01). `Basis.domain` is abstract on the ABC, all 6 shipped
bases answer, and the **energy/spatial** collision (`L2[coarse_cells_R1]`) is
unspellable.
⛔ **Its headline used to read *"no space claims an identity it does not
have"*, and that universal is FALSE** — `[M]` 2026-09-01 (§V.5g) an 8-node slab
ANGULAR space and an 8-node SPATIAL space on `[-1,1]` are still `==` AND
hash-equal. 2.1 closed one pair; the same defect class survives wherever two
roles honestly share a point set, and closing THIS one is **2.4**. The claim is
re-tensed to what was measured, per §2's quantifier clause.

✅ **A measure names the point set it lives on** — tracker **2.0c** LANDED
`025834f5` (2026-09-01), absorbing 2.0d. `Space = str` and its six `SPACE_*`
tags are retired across all six implementors; `phase`, the `L2[…]` name, the
fold and the tensor product are DERIVED from the manifold. Read the *2.0c
EXECUTED* section — it records four defects the TYPE found that the string
could not, and one repair the mutation battery caught with **no witness**.

✅ **A measure names the point set it lives on, and the phase-space factor it
belongs to is DERIVED from that** — tracker **2.0c** LANDED 2026-09-01,
absorbing 2.0d. Read the *2.0c EXECUTED* section: it records four defects the
TYPE found that the string could not, and one repair the mutation battery caught
with no witness.

✅ **The slab says what space its ordinates live on** — tracker **2.4**
LANDED `17501245` (2026-09-01). The slab's rule is
`gauss_legendre_on_polar_orbit(n, "x")`, on `S^2/SO2_x`; `SO2` is
axis-parameterised by user ruling; the `[-1,1]` angular/spatial collision is
unspellable; Part IV obstacle 1 is answered by derivation. **Read the *2.4
EXECUTED* section** — it records the ruling's evidence, the §6b set the red
loop actually found (five test-local factories, not the 23 sites), a
prediction of mine in a production comment that the §6c check refuted, and a
pre-existing 9 s hot spot the step amplified and then removed.

✅ **A basis states the symmetry group its functions are invariant under,
by naming what they eat** — tracker **2.1b** LANDED 2026-09-01 (hash stamped in
the follow-up). `Basis.invariance_group` is DERIVED (`domain.by`; `Trivial`
off the bare sphere; `None` off a non-angular domain), `@final`, zero subclass
edits; the fold's two halves read ONE object; ERR-080's pairing is the lattice
verdict `Trivial ⊉ SO2('x')`, refused by nothing yet. **Read the *2.1b
EXECUTED* section** — it records the fourteenth opener correction (the
requested FIELD dissolved into the domain, like 2.0d's), the category ruling
and its reason, an 8-arm battery with one arm blind by construction, and three
stale docstrings from earlier steps fixed on sight.

✅ **A measure's support is DERIVED from the map that built it — the product
rule's `S²` and a pushforward's target stop being literals** — tracker **2.3**
LANDED 2026-09-02 `5ec3a00a`. `ManifoldMap` (user-ruled name;
the row's *"typed `Chart`"* was refuted as imprecise), `archimedes(axis)`,
`barycentre(orbit_space)`; `pushforward(phi)` reads `phi.codomain` and refuses a
map out of the wrong point set; `spherical_product` IS
`(polar * azimuthal).pushforward(archimedes("z"))`; the fold's retraction and
`_embedded_nodes` read typed maps. **Read the *2.3 EXECUTED* section** — it
records the four rulings, the 60/5/15 bit-identity table, the eleven-arm
battery (one arm blind by construction, one control inside the rule's own
symmetry group, four arms the harness dropped), and the reference measure's
deferral to 3.1 with the §6d cycle MEASURED.

✅ **An orbit space is a catalogue entry that carries ALL of its derivation —
the quotient map and the measure it pushes forward included** — tracker
**3.1** LANDED 2026-09-02 `67e38605`. `Quotient` has 12 fields:
`orbit_coordinates` (the surviving invariants as a map on the base's
ambient coordinates) with the DERIVED arrow `quotient_map` (codomain the
ENTRY — user-ruled over D0.1's word *chart*, which 2.3 had already ruled
means `M ⊃ U → ℝⁿ`), and `reference` (`LEGENDRE` on the three axial
entries by the hat-box theorem; honest `None` on σ_a — a weighted disk
measure no shipped realization spells — and on `M/{e}` — Lebesgue on the
BASE, whose orthogonal system `Manifold` does not carry; both user-ruled).
`AngularSymmetry.reference` READS the entry; the registry's bare-sphere arm
stays because its Trivial support is `SPHERE`, never a quotient. **Read the
*3.1 EXECUTED* section** — the three rulings, the witness table (12/12,
3/3, 4 groups, 5/5), the §6d shadow-copy measurement (7/7 vs 0/7), and a
nine-arm battery whose m3/m4 pair triangulates *the registry reads the
entry, and the entry says what the table said*.

✅ **The angular moment space is READ off the frame, never minted from `L`** —
tracker **2.5**, MINTED 2026-09-02 at the fused step's opener (§V.5l) and
LANDED the same session as its own bit-identical commit BEFORE the fused one
(read the *2.5 EXECUTED* section). `[M]` 8 production
`SphericalHarmonicSpace.from_L` mints, two SH `isinstance` doors on
`HarmonicFrame`, `(name, shape)` guards between them and the frame — the
Stage-2 generator ruling's own defect, found by re-censusing 3.4-R's evidence.
Read §V.5l(b)–(d) and the test-architect's memo
(`scratch/_fused_verification_plan.md`) before designing.

✅ **The slab's P_L scattering is right because the basis its measure pairs
with is the trivial isotypic sub-basis of the spent group** — the fused
step ⭐ **0.1b + 0.6 + 2.2 + 3.4 (+ 3.4b)** LANDED 2026-09-02 `5436184e`
(read the *FUSED EXECUTED* section — ERR-080 CLOSED; G0 caught a
production defect on landing).

⏸ **NEXT: a context compaction** (user-ruled 2026-09-02), then **1.9**
(`SubgroupOfO3.O2(axis)`, #432) and **2.2b** (the Γ-slot), then 4.1, 4.9
and exit predicate 5's re-gate.
⭐ **Its opener RAN 2026-09-02 (§V.5l): four rulings taken, all recorded there —
read §V.5l first, then §V.5k.** It found
the fusion premise intact by direct reading, THE FIX's blast radius to be
the frame's TABLE LAYOUT on a 1-D rule (`[M]` 3 production `angular_frame`
readers, 13 `.table` reads in 5 files, the DSA μ-read off the ℓ=1 slot,
`moment_space_on` deriving the moment codomain from the basis's shape —
and the ruling's *15 sites / 12 sites* figures INHERITED from 2026-08-31
under two predicates, to be re-censused), the isotypic probe's idiom
shipping with no surviving reproducer, THREE rulings owed at the opener
(the probe's HOME; G0's predicate; the Γ-slot), and exit predicate 5 being
MEASURED (10417 collected, predicted first). ⚠ Open
with the opener (every one so far has corrected its own section): §III.8's
probe with INCOMMENSURATE angles (vv-principles #13; the mirror idiom ships
at `spherical_harmonic_basis.py:498`), `Descent` pulling `P_ℓ` back along
`Quotient.quotient_map`, its claim against `Quotient.reference`, 2.2's G0
ruling owed (§V.5h(a)), and the 13-tree gate baseline measured FIRST (exit
predicate 5 — best spent immediately before this step, where it is
load-bearing).

⛔ **The text below is the PRE-3.1 pointer, kept per §3.** Every claim in it
held at execution; the map's NAME (*chart* → `quotient_map`) and MECHANISM
(a stored apply + a derived arrow, not a stored `ManifoldMap`) were
user-ruled at the opener, and the trivial entry's reference was ruled
`None` rather than `UNIFORM_ON_SPHERE`.

▶ **An orbit space is a catalogue entry that carries ALL of its derivation —
including the measure the quotient map pushes forward** — tracker **3.1**,
scoped to `SO(2)`. ⭐ **Pre-flighted at this compaction — read §V.5j before
designing.** It found the row's premise pre-2.0a (`[M]` the catalogue lives in
`manifold.py`, six keys, ruled at E. Placement; `class IsometryGroup` 0 — do
not mint a second group type), the `symmetry.py` package re-home OFF the exit
path (XII.2) and a §6d step of its own (`numerics/__init__.py:43` re-exports
`SubgroupOfO3`; the `numerics → geometry` back-edge runs through it), and the
honest remaining scope to be TWO fields and ONE map: the entry's chart
`Ω ↦ Ω·ê_a` as a `ManifoldMap` (witness `π ∘ archimedes_a = pr₁`, measured
bit-exact) and its pushforward reference (`LEGENDRE`, populated INSIDE the
derivation function; `[M]` 3 constructor sites, 1 production reader of the
registry twin at `registry.py:1272`; the σ_y entry ships `None`, honestly).
"The probe" is 3.4's (§III.8; the mirror idiom ships at
`spherical_harmonic_basis.py:498`). ⚠ Open with the opener: re-measure 3.1's row (it names a
`numerics/symmetry/` PACKAGE re-home — §6d both ways, and the import-linter's
declared contract first) and the 2.3 residue (a `Quotient.reference` field
populated INSIDE `_sphere_mod_so2` / `_mod_trivial` via function-scope import;
`AngularSymmetry.reference` then reads it; the σ_y entry ships `None` because
its pushforward is a weighted disk measure no shipped `ReferenceMeasure` can
spell). Then the fused step, where 2.2 owes its G0 ruling (§V.5h(a)).

⛔ **The text below is the PRE-2.3 pointer, kept per §3.** Every claim in it
held at execution; the name it uses did not.

▶ **A measure's support is DERIVED from the map that built it — the product
rule's `S²` and a pushforward's target stop being literals** — tracker
**2.3** (the typed `Chart`). ⭐ **Pre-flighted at this compaction — read
§V.5i before designing.** It settled §V.5h(e) by census (`[M]` 0 production
readers of `fundamental_domain`), found the tree already spelling THREE maps
around the quotient — the quotient map in `quotient()`, the Archimedes chart
in `spherical_product`, and the BARYCENTRE map `S²/SO(2)_a → B³` spelled
TWICE (honestly in `_embedded_nodes`, dishonestly as ERR-080's forgery; `Ball`
has 0 production consumers) — scoped the done-when to MAP-BUILT measures (2 of
5 `support=SPHERE` literals), enumerated `pushforward`'s §6b set (1 + 8), and
measured a §6d cycle (`manifold → exactness → manifold`) waiting for anyone
who puts the pushforward measure on the seed at module scope. ⚠ The
§V.5h(e) question it answered: *"the SECTION MAP's home"* (this Part's pre-2.4
pointer) was a HYPOTHESIS — a section `S²/SO(2)_a → S²` exists only to do what ERR-080 does,
and after 3.4 it may have NO consumer, while the QUOTIENT MAP
`π: S² → S²/SO(2)_a` has one in 3.4b (`Descent` pulls `P_ℓ` back along it).
Enumerate who would call a section after 3.4; if nobody, 2.3 is the quotient
map + the Archimedes parametrisation of `product()`, and `_sphere_mod_so2`'s
`fundamental_domain=None` stays `None` on purpose. `[M]` at `8fc781d2`:
`Quotient` has **10** fields, `class Chart` **0**, `pushforward` has **1**
production caller (`measure.py`, inside `quotient()`), `spherical_product`
column-stacks `S²` nodes against a `support=SPHERE` literal
(`rules_product.py:531`). ⭐ A Pattern-2 twin to absorb: the pushforward
MEASURE of a quotient is answered today by `AngularSymmetry.reference`'s
dispatch (`LEGENDRE` for any `SO2`, `UNIFORM_ON_SPHERE` for `Trivial`) — the
catalogue entry's field living on the registry; the seed should carry it and
the registry should read it. Then 2.2 owes its G0 ruling at the fused step
(§V.5h(a)).

⛔ **The text below is the PRE-2.1b pointer, kept per §3.** Its 2.1b half
held at execution except for the sizing in §V.5h(f) (the §6b set was empty,
not 23) and the MEANS (a derived property, not a field); its 2.3 half is
carried forward above unchanged:

▶ **A basis states the symmetry group its functions are invariant under, and
the frame checks its two halves agree** — trackers **2.1b** (`Basis.invariance_group`,
G0's other side: `[M]` 0 of 6 subclasses answer it) and **2.3** (the typed
`Chart` — ⚠ *"the SECTION MAP's home"* is a HYPOTHESIS, see §V.5h(e): the
section may have no consumer after 3.4, while the QUOTIENT MAP has one in
3.4b). ⭐ **Pre-flighted at this compaction — read §V.5h before designing.**
It found Part IV's *"G0 cannot be a tag-equality test"* VOID (an equality G0
refuses exactly the Part I bug and admits both shipped angular pairings — a
RULING for 2.2's opener), G2's measure side already real (the slab's
`quotient_group` is `SO2('x')`), 2.2's frame gate inert until 0.1b (the frame
still sees the forged `S^2`), and 2.3's *"6 → 8"* premise pre-dating
`Quotient.by` and `on_orbit_space` (`Quotient` has 10 fields; chart and
pushforward measure absent; `pushforward` has ONE production caller).

⛔ **The text below is the PRE-2.4 pointer, kept per §3.** It described the
tree at `8cc69e7f`; its five findings all held at execution except the second
half of #4 (the fallback arm is reachable by the chart-level μ-rule) and the
sizing in #3 (neither number was the set):

~~▶ **The slab says what space its ordinates live on, and stops naming it by the
interval a chart happens to map onto.**~~ Tracker **2.4** — the slab's angular
support becomes `SPHERE.quotient(SubgroupOfO3.SO2)`. ⭐ **Pre-flighted
2026-09-01 at `8cc69e7f` — read §V.5g before designing; it re-classed the
step.** Five things a fresh session must not re-derive:

1. ✅ **The premise holds.** `[M]` `S²/SO(2)`'s `realization` **is**
   `Interval(-1,1) == COSINE_INTERVAL`, and `Quotient.contains(μ-nodes)` is
   `True` — the declaration changes **no coordinate**, only what the object
   knows about itself.
2. ⛔⛔ **It repairs a LIVE defect, so it is NOT wiring — it is 2.1's class.**
   `[M]` an 8-node slab **angular** space and an 8-node **spatial** space on
   `[-1,1]` are `==`-equal **AND hash-equal** (both `L2[[-1,1]]`, shape `(8,)`).
   2.0c could not fix this: it made the supports *honest*, and they are honestly
   identical. The slab's support is still named by its CHART CODOMAIN — ERR-080's
   defect class, stated at the level of spaces. `[M]` after the declaration,
   `L2[S^2/SO2]` vs `L2[[-1,1]]` are unequal and the collision is unspellable.
3. ⛔ **It is the ONE step that RENAMES a shipped space.** Every prior step
   preserved names exactly (`[M]` 2.0c: 10 of 10). Size it by the **23**
   space-comparison sites in `orpheus/`, **not** by the 3 literal `L2[[-1,1]]`
   pins — a literal census is the same wrong-side error §V.5f(d) caught in
   2.0c's own row.
4. ⭐ **It discharges 2.0c's surviving fallback arm.** `[M]` a measure on the
   quotient answers `phase == "angular"` with `invariance_group=None`, from the
   manifold alone. ⚠ Per §6c, do **not** delete the arm on that basis — `phase`
   is public and an `Interval` measure is still constructible; check whether any
   SHIPPED input reaches it, and if none does, say so in the docstring.
5. ⛔ **Part IV's obstacle 1 is real and needs a RULING, not a measurement.**
   `[M]` `symmetry.py:121` — `SO2 = "SO2"` is a bare member with **no `.axis`**,
   while `Mirror` is parameterised. So `S²/SO2` does not say which axis, and the
   shipped slab rule declares `invariance_group=Mirror('x')`. Decide:
   axis-parameterise `SO2` (changing the catalogue key) or fix the convention by
   fiat and document it. **Surface this before building.**

⭐ And ERR-080's level-1 half still needs the honest `φ = 0` half-meridian; the
tree fabricates it by zero-padding to `(μ, 0, 0)`, which is off `S²`.

⛔ **The text below is the PRE-2.0c pointer, kept per §3.** Its opener ran
2026-09-01 and is **§V.5e**; **§V.5f** is the third opener, which re-scoped the
step again. Four things it established:

1. ✅ **Unblocked.** `[M]` `SPHERE.quotient(SubgroupOfO3.Mirror("y"))` answers
   `'S^2/sigma_y'` at this HEAD; 1.1's entry landed at `b55bba56`.
2. ⛔ **The tracker row's sizing is the wrong side.** Its two numbers reproduce
   exactly (54 literals, 94 % in `tests/`) and both count PRODUCERS. A retype
   breaks CONSUMERS: `[M]` **16** sites read `.support` into a string operation
   (15 `orpheus/`, 1 `tests/`). A producer handing a `Manifold` where a `str`
   was expected fails LOUDLY; a consumer interpolating one silently emits
   `L2[Sphere()]`. §6b set ≈ **87 producers + 16 consumers**, one commit.
3. ⭐⭐ **The prize is not the retype — it is that `DiscreteMeasure.phase`
   (`measure.py:369`) is stringly-typed dispatch** on the tag 2.0c types
   (`startswith("spatial")` / `== "cells"` / `startswith("energy")`).
4. ⛔⛔ **A LIVE defect to fix inside it:** `[M]`
   `folded_product(4,8).measure.phase` **RAISES** — the fold's
   `invariance_group` is legitimately `None` and the dispatch has no arm for a
   quotient. Latent (`[M]` **0** production `.phase` consumers) and still wrong.
   ⚠ It also refutes the obvious design: `phase` **cannot** be a total function
   of the manifold while the slab's support is `[-1,1]`, which is ambiguous
   between the angular μ-axis and a spatial interval. The collapse to one
   manifold-derived answer needs the slab to declare `SPHERE.quotient(SO2)` —
   that is **2.4**, not 2.0c. Keep the `invariance_group` arm; add a quotient arm.

⚠ And one obligation 2.1 HANDED it: `test_d6` pins the one production pair
whose two halves do not agree in spelling — `LossKernelBasis`'s measure tags
the bare label `sn_trace_orbit…` while `IndexSet` spells it `index(…)`. 2.0c
must resolve that, and cannot do so by accident.

*Then* the rest of G0, and the fix.

⛔ **The text below is the PRE-2.1 pointer, kept per §3 (2.1 has since
LANDED — this described the tree before `c461fe8d`).** Tracker **2.1** —
`Basis.domain: Manifold`. It is the only remaining exit-path item that repairs
something ALREADY WRONG rather than enabling a later step: `[M]` re-verified at
`a2befd9b`, a 2-group ENERGY space and a 2-cell SPATIAL space are `==`-equal
AND hash-equal, both minting `L2[coarse_cells_R1]`, with a passing negative
control. **Read the 2.1 PRE-FLIGHT section before designing** — it checks the
sizing (18 sites is right, and why the 39-site family does not widen it), and
it flags that the row's own done-when (`domain is SPACE_SPHERE`) is
pre-`Manifold` text that would re-introduce `Space = str`.
⛔ Both of that paragraph's own claims were later refuted at execution: the §6b
set was **23**, not 18, and the done-when's spelling would have re-introduced
`Space = str`. See the supersession table.

### ⏏ The ORDER, and why it is this order

Part XII.2 says WHICH items are on the critical path; this says in what
sequence, which is not the tracker's numbering. Per-item status lives in
Part XIII — do not copy it here (`plan-authoring` §9).

| # | items | why here |
|---|---|---|
| ✅ | **0.1a**, **0.1c** | LANDED `2c1a06b1` |
| ✅ | **0.2** | LANDED `ce46181c` |
| ✅ | **2.1** | LANDED — see the 2.1 EXECUTED section. ⚠ It leaves ONE thing for 2.0c that did not exist before: `test_d6` pins the `LossKernelBasis` spelling divergence (bare label vs `index(...)`), so 2.0c must resolve it rather than inherit it |
| ~~1~~ | ~~**2.1**~~ | ⭐ the only remaining exit-path item that repairs something ALREADY WRONG. `[M]` **re-verified at `a2befd9b`**: a 2-group ENERGY space and a 2-cell SPATIAL space are `==`-equal AND hash-equal (both `L2[coarse_cells_R1]`, shape `(2,)`), negative control passing. Sizing also re-verified — **18 `IndicatorBasis(` sites, 4 `orpheus/` + 14 `tests/`**, and the 39-site family does NOT widen it (siblings inherit the repair by composition). ⚠ the *"26 space-comparing sites"* figure is INHERITED and was **not** re-measured — treat as unverified. **Read the 2.1 PRE-FLIGHT section** |
| ✅ | **2.0c** (+ **2.0d** absorbed) | LANDED — see the *2.0c EXECUTED* section. `Space = str` and its six tags retired across all six implementors; `phase` dispatches on the manifold's TYPE (user-ruled); the LIVE `folded_product(…).measure.phase` defect fixed; `pushforward`'s fabricated default retired; 2.1's `LossKernelBasis` handoff discharged **and witnessed** (its repair was measured BLIND first) |
| ✅ | **2.4** | LANDED `17501245` — see the *2.4 EXECUTED* section. The defect fix of 2.1's kind: `[M]` the `[-1,1]` angular/spatial collision is unspellable; `SO2(axis)` user-ruled; stage 0 now refuses a chart-level rule and a wrong-axis rule. ⛔ its pre-flight's sizing (*23 sites / 3 literals*) was the wrong instrument twice — the §6b set was the five test-local factories the stricter gate refused |
| ✅ | **2.1b** | LANDED `9b4a4d9c` — see the *2.1b EXECUTED* section. G0's other side by DERIVATION (`domain.by`), one `@final` property, zero subclass edits; the four-row table runs as `test_e5`; the fold's halves read ONE object. ⛔ its §V.5h(f) §6b sizing was void (the set was empty) |
| ✅ | **2.3** | LANDED 2026-09-02 `5ec3a00a` — see the *2.3 EXECUTED* section: `ManifoldMap` + `archimedes` + `barycentre`, `new_space=` retired, 60/5/15 bit-identical, reference measure deferred to 3.1. ⏬ *The pre-landing row, kept per §3:* **NEXT — read §V.5i (pre-flighted 2026-09-02).** ~~G0's other side (`[M]` `Basis.invariance_group` 0 of 6)~~ (✅ 2.1b, by derivation); the typed `Chart` (also the SECTION MAP's natural home — `_sphere_mod_so2`'s `fundamental_domain=None` comment says exactly what 2.3 must choose). ⚠ Open with the re-measure: 2.3's *"6 → 8"* premise predates `Quotient.by` (2.0c) and `on_orbit_space` (2.4) |
| ✅ | **3.1** | LANDED 2026-09-02 — see the *3.1 EXECUTED* section: `Quotient.orbit_coordinates` + the derived `quotient_map` (codomain the ENTRY, user-ruled), `Quotient.reference` (`LEGENDRE` on the axial entries, honest `None` on σ_a and `M/{e}`, user-ruled), `AngularSymmetry.reference` reads the entry; §6d measured on a shadow copy (function scope alive 7/7, module scope dead 7/7); nine-arm battery, no arm blind. ⏬ *The pre-landing row, kept per §3:* **NEXT — read §V.5j (pre-flighted 2026-09-02).** ~~the catalogue home + the probe~~ ⛔ the HOME is landed (`manifold._ORBIT_CATALOGUE`) and the PROBE is 3.4's; what remains is the entry's chart map + its reference measure. ⭐ It inherits 2.3's deferral: the pushforward REFERENCE measure as a `Quotient` field populated INSIDE the derivation function (function-scope import — `[M]` module scope is a 5-of-5 cycle), read by `AngularSymmetry.reference`; and §6d BOTH ways before the `symmetry.py` → package re-home |
| ✅ | **2.5** | LANDED 2026-09-02 `58ad7c28` — see the *2.5 EXECUTED* section (bit-identical; the moment space's ONE home is the frame's basis) |
| ✅ | ⭐ **0.1b + 0.6 + 2.2 + 3.4 (+ 3.4b)** | LANDED 2026-09-02 `5436184e` — THE FIX; read the *FUSED EXECUTED* section |
| ✅ | ⭐ **1.9** | LANDED 2026-09-02 `c1d53206` — an orbit space is named by its STABILISER (`SubgroupOfO3.O2(axis)`, `orbit_stabiliser`); read the *1.9 EXECUTED* section |
| ✅ | ⭐ **2.2b** | LANDED 2026-09-02 `<HASH>` — the Γ-slot: a geometry admits a FOLD of its domain (stage 0 through the descent arrow), and the owed symmetry is asked ON the orbit space (`Quotient.induced_action`, the normaliser, one closure); read the *2.2b EXECUTED* section |
| ⏸ | **COMPACTION** | user-ruled 2026-09-02: a context compaction follows the fused commit; the two commits below are the first steps after it |
| ⏳ | **1.9** — `SubgroupOfO3.O2(axis)` | IN FLIGHT 2026-09-02 (built, green on its trees, uncommitted — the ⏳⏳ block at the top of the resume section); its own commit, FIRST after the compaction (user-ruled, memo H-6) |
| **3** | **2.2b** — the Γ-slot | ✅ LANDED 2026-09-02 `<HASH>`, its own commit |
| 4 | **4.1**, **4.9**, then the six exit predicates (XII.3) | |

⚠ **2.1 and 2.0c were independent** — both needed only `Manifold`. 2.1 went
first because it was smaller and fixed a live defect. ✅ That is settled;
2.0c is now the head of the queue, and it inherits one obligation from 2.1
(the `LossKernelBasis` spelling divergence, pinned by `test_d6`).

### ⛔ Why the FUSED step (now #5) is the campaign's highest-risk landing

`plan-authoring` §6c's MIRROR: not a gate with no witness, but **a gate whose
only witness is PRODUCTION.** 2.2's G0 refuses exactly the pairing the tree
ships today, so arming it before 3.4 changes that pairing reddens every slab
solve. The temptation the interval creates is to weaken the gate so the tree
stays green — which is how a real refusal decays into a warning. Land the four
together or not at all.

### ⭐ What this session's ruling changed about a LATER step

⛔ **SUPERSEDED — kept per §3.** 2.3 LANDED (`5ec3a00a`) without the section
map (§V.5i(a) refuted it by census) and under the name `ManifoldMap`; the
paragraph below is the 2026-09-01 reading.

**2.3 (the typed `Chart`) gained weight.** It was *"completes the engine seed
6 → 8"*; after the two-slot ruling it is also the natural home for the SECTION
MAP — and `[M]` the honest `φ = 0` half-meridian is exactly what ERR-080's
level-1 half needs, since the tree currently fabricates it by zero-padding to
`(μ, 0, 0)`, which is off `S²`.

### ⛔ Landmines for the pick-up session's own greps and censuses (numbered, newest LAST)

1. **`invariance_group` returns ~55 grep hits and the plan's claim still holds.**
   `[M]` exactly ONE class defines it — `DiscreteMeasure` — meaning *the group
   this measure IS invariant under*, **not** the group it was quotiented BY.
   `[M]` **0 of 6** `Basis` subclasses define it. §3's ambiguous-name hazard.
2. **`Basis.domain` returns ~45 definers and they are ALL OPERATORS.** §III.10
   ruled the collision unreachable, and this is why a bare grep says
   otherwise. ⚠ **Updated 2026-09-01:** `[M]` since `c461fe8d` **6 of 6**
   `Basis` subclasses answer `domain` too, so the grep is now genuinely
   ambiguous between two LIVE meanings — `dom` in **Man** for a basis, `dom`
   in **Vect** for an operator. Discriminate by the class, not by the name.
3. ✅ **DISCHARGED at 2.0c** — `measure.py`'s `self.support == "cells"` arm is
   gone with the string dispatch, and its only producer (a `parametrize` list)
   retired with it. *(Kept per §3: it was a **production-unreachable arm with a
   synthetic witness**, not a dead one — §V.5e(e).)*
4. ⭐ **NEW — a `Basis` subclass census by AST UNDERCOUNTS.** `[M]` AST said
   3 direct / 5 recursive; runtime says **4 / 6**, and the runtime answer is
   the one that matches §II.15's independent count. Inheritance is a runtime
   relation. Use `__subclasses__` after importing the package.
5. ⭐ **`axis_cosines` RAISES on a suppressed axis** (since `ce46181c`). Any
   probe, diagnostic or scratch script that does
   `column_stack([q.axis_cosines(a) for a in range(3)])` will die on a 1-D
   rule. That composition is the R³ **embedding**: call
   `_embedded_nodes(q.measure)`. For the flux question, `mean_axis_cosine`.
6. ⭐⭐ **NEW at 2.0c — `support=` takes a `Manifold`, never a string.** Any
   scratch script, probe or diagnostic written before 2026-09-01 that does
   `DiscreteMeasure(..., support="S^2")` now stores a `str` in a field the rest
   of the tree calls `.name`/`.quotient()` on, and fails at the *consumer*, not
   at construction (nothing validates the field's type). The members are in
   `orpheus.numerics.manifold`: `SPHERE`, `CIRCLE`, `COSINE_INTERVAL`,
   `UNIT_INTERVAL`, `HALF_LINE`, `REAL_LINE`, `ENERGY`, plus `RealSpace(d)`,
   `IndexSet(label=…)`, `EnergyGroups(n)`, `Interval(a, b)`.
   ⚠ ~~And `pushforward` now **requires** `new_space`.~~ ⛔ superseded at
   2.3 — it takes a `ManifoldMap`; see landmine #10.

7. ⭐ **`IndicatorBasis` needs `partition_of`** (since `c461fe8d`).
   Every construction site takes a `Manifold`; a scratch script written before
   that dies with `TypeError: missing required positional argument`. The
   manifold is whatever the edges partition — `RealSpace(d)`, `EnergyGroups`,
   `IndexSet`, or `Interval` for a partition by physical VALUE.
8. ⭐⭐ **A CENSUS-READING trap, not a grep one.** A census that
   helpfully splits its hits into "literal" and "non-literal" invites you to
   answer a universal from the LITERAL half. `[M]` 2026-09-01 I concluded
   *"`support='cells'` has no producer"* from exactly that split; the real
   producer is a `parametrize` list feeding `support=support`, which my own
   table had already printed as a non-literal row. **State a universal over
   the union, or not at all.**
9. ⚠ **A hash in this plan can be AMENDED AWAY.** `[M]` `de29bcc6` was
   cited by two tracker rows and is orphaned (amended into `b8c05d16` 39 s
   later; no branch contains it, so it will be GC'd). Corrected 2026-09-01.
   ⟹ on pick-up, run `git merge-base --is-ancestor <h> HEAD` over EVERY hash
   the plan cites, not just the ones you intend to use — a dangling hash reads
   exactly like a live one.
10. ⭐⭐ **NEW at 2.3 — `pushforward` takes ONE argument, a `ManifoldMap`.**
    `μ.pushforward(fn, new_space=N)` (2.0c's spelling, alive for one day) is a
    `TypeError`; spell `μ.pushforward(ManifoldMap(μ.support, N, fn))`, and
    the map's `domain` must `==` the measure's support or it is refused by
    VALUE (the slab's `S^2/SO2_x` rule is refused by a map out of `[-1,1]`).
    `symmetry._AXIS_INDEX` is gone — `manifold.AXIS_INDEX`. A `grep` for
    `.reference` is ambiguous between `AngularSymmetry.reference` (the
    registry twin, 1 production reader) and `ExactnessClaim.reference`
    (46 hits tree-wide). And a mutation battery whose per-arm summary is
    read through a shell `$(...)` capture can return an EMPTY summary for
    the widest arms (`[M]` 4 of 11 at 2.3) — read from a log file.
11. ⭐ **NEW at 3.1 — `Quotient(...)` REQUIRES `orbit_coordinates`.** Any
    scratch script that builds a raw `Quotient(base=, by=, realization=)`
    dies with `TypeError: missing required keyword-only argument`; go
    through `M.quotient(H)`. `Quotient` compares/hashes WITHOUT that field
    (`compare=False`), so two honest constructions of one entry are still
    equal. `AngularSymmetry.reference` no longer imports `LEGENDRE` —
    a grep for `LEGENDRE` in `registry.py` finds only docstrings. And the
    manifold module docstring's *"imports nothing from numerics"* is now
    a MODULE-SCOPE claim: `_sphere_mod_so2` imports `generating_measure`
    at function scope, measured safe (7/7) on a shadow copy.
12. ⭐⭐ **NEW at the fused step's opener — the angular moment space has EIGHT
    homes, and the frame is not the one production reads.** `[M]`
    `SphericalHarmonicSpace.from_L(L)` is minted at 8 production sites
    (`scattering.py:207/378/521`, `fission.py:510`, `n2n.py:327`,
    `_bases.py:682`, `harmonic_moment_flux.py:334`, the basis's own `space`);
    `HarmonicFrame.__init__` (`harmonic_frame.py:331`) and `from_galerkin`
    (`:368`) REFUSE any basis that is not a `SphericalHarmonicBasis`;
    `SphericalHarmonicSpace.__post_init__` (`:153`) refuses any shape but
    `(L+1, 2L+1)`; `OperatorProduct` compares ends by `(name, shape)`. A
    census of `Y[:, l, m]` READS cannot return any of these. Until 2.5 lands,
    a basis with a different `space.shape` fails at the DOOR of the first
    `for_space` call (L = 0, fission), not at the table. And `angular_frame`
    is called at L = 0 on EVERY solve (fission + n2n through `for_space`),
    and by DSA unconditionally in `L`.
13. ⭐⭐ **NEW at the fused commit — a 1-D rule's frame is LEGENDRE.**
    `Quadrature.gauss_legendre(n).angular_frame(L).table` is `(N, L+1)`, its
    `basis` a `LegendreBasis`, its `basis.space` a `LegendreSpace` (name
    `legendre_space(S^2/SO2_x)`, field `spent_axis`); every scratch script
    that indexes `[:, 1, 1]`, slices `[l, :2l+1]`, or builds `(L+1, 2L+1, …)`
    moment arrays on a slab dies or reads garbage. Read the head:
    `field.head.isotropic_slot` / `.degree_block(l)`. `Quadrature.spherical_harmonics`
    is GONE (read `angular_frame(L).table`); `_harmonic_frame_measure` is
    GONE; `_PARITY_PROBE_DIRECTIONS` is GONE (the entry's probe is
    `Quotient.descending_slots`). `ScatteringMaterialField.moment_source` /
    `N2NMaterialField.moment_emission` REQUIRE `head=` (the operator's
    domain). `OverlapBasis.from_indicator` REQUIRES `fine=` (what the basis
    EATS — G0 caught the misdeclaration). `SphericalHarmonicBasis.evaluate`
    REFUSES non-unit directions (normalise probe directions). A
    `GalerkinFrame(basis, measure)` REFUSES when no quotient map
    `measure.support → basis.domain` exists (message fragment *"no quotient
    map"*); a sized `EnergyGroups(n)` and the unsized `ENERGY` are UNEQUAL
    manifolds under it — name one object on both halves.
14. ⭐ **NEW at the compaction prep after the follow-up — `O2` is a name this
    module REJECTED once, and a bare `O2(` grep lies.** 1.9 (#432) mints
    `SubgroupOfO3.O2(axis)`. `[M]` `orpheus/numerics/symmetry.py:136-148`:
    `Dinfh` was NAMED `O2` until 2026-08-02 and was renamed because that entry
    realizes axial rotations + σ_h (C_∞h), while true O(2) embedded in 3-D is
    C_∞v — rotations about the axis + the VERTICAL mirrors — which is exactly
    the group #432 wants (`P_ℓ(Ω·ê_a)` is invariant under rotations about `a`
    and mirrors through planes CONTAINING `a`). So 1.9's member is the C_∞v
    realization, NOT a re-parameterization of `Dinfh`, and that note's history
    sentence (*"named `O2` … wrong in a load-bearing way"*) must be reconciled
    IN THE SAME COMMIT so it cannot read as condemning the new member
    (`plan-authoring` §3's ambiguous-name hazard; §1's rejected-name clause —
    the name is free BECAUSE it was rejected). ⚠ The census: `grep 'O2('`
    returns **14** hits in `symmetry.py`, 12 of them `SO2(`; use
    `grep -P '(?<![A-Za-z_])O2\b'` (**2** hits, both the note). `[R]` the
    lattice edge is `SO2(a) ⊆ O2(a)` for every `a`, and `O2(z) ⊆ Dinfh` only —
    `Dinfh` is parameter-free about z; for `a ≠ z` the only shipped supergroup
    is `O3`. Reason it from the realizations, not from this row.
15. ⭐⭐ **NEW at 1.9 — `SPHERE.quotient(SubgroupOfO3.SO2(a))` RAISES, at the door AND at
    construction.** The orbit space is `SPHERE.quotient(SubgroupOfO3.O2(a))`, name `S^2/O2_x`,
    `by == O2(a)`; every scratch script spelling the SO2 quotient dies with `ValueError:
    S^2/SO2_x is the orbit space S^2/O2_x …`. The refusal is `_assert_named_by_stabiliser`
    reading `group.orbit_stabiliser` — called by `Quotient.__post_init__` and by
    `_catalogued_quotient` BEFORE the lookup, never by a builder; so `replace(entry,
    by=SO2(a))` raises too, and `SPHERE.quotient(SO3)` raises naming `O3`.
    `_ORBIT_CATALOGUE` has SIX keys (`O2_*`, `sigma_*`) and NO `SO2_*` key.
    `_sphere_mod_so2` is GONE (`_sphere_mod_o2`, a pure derivation); `_fixes_axis(tag, axis)`
    has no `proper_only` (the SO(2) half is `∩ SO(3)` in `_axial_contains`);
    `AXIS_LETTER` lives beside `AXIS_INDEX` in `manifold.py`. `quotient_group`
    of a slab rule is `O2('x')`; `LegendreBasis.invariance_group` is `O2(axis)`; the space is
    `L2[S^2/O2_x]` / `legendre_space(S^2/O2_x)`. The lattice test `_NAMED`/`groups` lists
    that build quotients must spell O2; the SO2 GROUP itself still exists for lattice
    questions (`Cn ⊆ SO2_z`, `SO2_z ⊆ Dinfh`). A `grep 'O2('` now returns the real member
    beside the 12 `SO2(` hits — use the word-bounded form.

### §1 existence-checks — re-run 2026-09-02 at `2f294ef1` (post-3.1)

`[M]` by Python regex over `orpheus/**` + `tests/**` (`scratch/_p31c_reground.py`),
positive controls `class Quotient` **1**, `class ManifoldMap` **1**,
`class SphericalHarmonicBasis` **1**, `class MirrorEvenSphericalHarmonicBasis`
**1**: `class Chart` **0** · `class LegendreBasis` **0** · `class Descent`
**0** · `class ReynoldsProjection` **0** · `class OrbitAxis` **0** · `class
IsometryGroup` **0** — nothing the fused step promises exists yet. Hash
check, same predicate as below: **36** distinct 8-hex hashes in this Part
before this rewrite, **35** ancestors of HEAD, 1 dangling (`de29bcc6`,
known). `[M]` by RUNTIME `__subclasses__`: `Basis` 4 direct / 6 recursive —
unchanged. `support=SPHERE` keyword literals in `orpheus/` by AST: **4**
(the forgery `directional.py:794`, `UNIFORM_ON_SPHERE`, Lebedev,
level-symmetric); ⚠ a regex reads 7 (three docstrings). `dataclasses.fields(Quotient)`
**12**. `select_quadrature` production callers **0** (26 tests). The
ERR-080 gate: **3** strict-xfail rows by AST.

⭐ **3.1 added these landmarks** (grep at landing, 2026-09-02): `Quotient`
`manifold.py:574` · `orbit_coordinates` `:612` · `reference` `:667` ·
`quotient_map` `:704` · `ManifoldMap` `:822` · `archimedes` `:940` ·
`barycentre` `:1010` · `_ambient_columns` `:1084` · `_all_coordinates`
`:1099` · `_sphere_mod_so2` `:1136` (function-scope `LEGENDRE` import
`:1194`) · `_sphere_mod_mirror` `:1255` · `_mod_trivial` `:1378` ·
`_ORBIT_CATALOGUE` `:1435` · `AngularSymmetry.support` `registry.py:875`,
`.reference` `:910` (guard `:964`, the read `:965`, the bare-sphere arm
`:979`), the V-stage reader `:1302`, `admits_domain` `:981`,
`GEOMETRY_ANGULAR_SYMMETRY` `:1019`, `select_quadrature` `:1130`. ⭐ **The
fused step's landmarks**: `_harmonic_frame_measure` `directional.py:730`
(the forgery `:794`), `angular_frame` `:678`, `spherical_harmonics` `:623`,
`axis_cosines` `:294`, `mean_axis_cosine` `:327`, `Quadrature.quotient`
`:548`, `gauss_legendre` `:870` · `SphericalHarmonicBasis`
`spherical_harmonic_basis.py:121`, `_PARITY_PROBE_DIRECTIONS` `:410`,
`MirrorEvenSphericalHarmonicBasis` `:422`, `even_slot_mask` `:498`,
`_evaluate_real_sh` `:550` · `Basis` `basis/base.py:120`, `.domain` `:247`,
`.invariance_group` `:291` · `test_q8_4` `test_quadrature_directional.py:568`
· the ERR-080 gate rows `:228`/`:255`/`:277` · the DSA μ-read `dsa.py:585` ·
`HarmonicFrame.from_galerkin`'s caller `harmonic_frame.py:416`.

⛔ **The block below is the 2.3-era re-run, kept per §3.**

### §1 existence-checks — re-run 2026-09-02 at `14c37aa9` (post-2.3)

`[M]` by Python regex over `orpheus/**` + `tests/**` (`scratch/_p23_reground.py`),
positive control `class Quotient` **1** (`manifold.py:562`): `class Chart`
**0** (its ruled successor `class ManifoldMap` **1**, `manifold.py:742`) ·
`class LegendreBasis` **0** · `class Descent` **0** · `class ReynoldsProjection`
**0** · `class OrbitAxis` **0** · `class IsometryGroup` **0** (3.1's row's
name for what ships as `SubgroupOfO3` — do not mint it). ⛔ The `Ball`
fact below was WRONG as written (its census excluded `manifold.py`; `Ball(2)`
is the σ_y realization at `:922`): the honest fact is *no consumer outside its
module, and `Ball(3)` never constructed* — `barycentre` (`:930`) is now the
first arrow with a ball codomain. Hash check, same predicate as below: **34**
distinct 8-hex hashes in this Part, **33** ancestors of HEAD, 1 dangling
(`de29bcc6`, known). `[M]` by RUNTIME `__subclasses__`: `Basis` 4 direct / 6
recursive, 6 of 6 answer `domain` AND `invariance_group` — unchanged.
`support=SPHERE` code literals in `orpheus/`: **4** (was 5; the map-built one
in `rules_product.py` is gone; the forgery's stays by design until 3.4).

⛔ **The block below is the 2.1b-era re-run, kept per §3.**

### §1 existence-checks — re-run 2026-09-02 at `4d946988`

`[M]` in `orpheus/`, by AST, **with positive controls in the same pass**
(`class Manifold` → `manifold.py:113`, `Quotient` → `:556`, `SO2` →
`symmetry.py:265`, `Mirror` → `:200` — all found):
`class Chart` **0** · `class LegendreBasis` **0** · `class Descent` **0** ·
`class ReynoldsProjection` **0** · `class OrbitAxis` **0**.
⟹ nothing the pointer promises exists yet. `[M]` re-run at `4d946988` by
Python regex over `orpheus/**` + `tests/**`: `class Chart` **0** ·
`class LegendreBasis` **0** · `class Descent` **0** · `class ReynoldsProjection`
**0** · `class OrbitAxis` **0**, positive control `class Quotient` **1**. ⭐ And
one existence fact the other way: `Ball` exists (`manifold.py:454`) with **0**
production consumers — §V.5i(b) names its first. Hash check, predicate stated:
**31** distinct 8-hex hashes cited in backticks within this Part, **30**
ancestors of HEAD, 1 dangling (`de29bcc6`, known).

`[M]` by **RUNTIME** `__subclasses__` after importing every `orpheus` module
(⛔ *not* AST — landmine #4): `Basis` has **4 direct / 6 recursive**, all six
`orpheus`-defined — `IndicatorBasis`, `OverlapBasis`, `WeightedIndicatorBasis`,
`SphericalHarmonicBasis`, `MirrorEvenSphericalHarmonicBasis`, `LossKernelBasis`.
✅ **6 of 6 answer `domain`**; ~~`invariance_group` is still **0 of 6**, so G0's
other side (tracker **2.1b**) remains absent~~ ✅ **6 of 6 answer
`invariance_group`** since 2.1b (2026-09-01), by inheritance from ONE concrete
property (`basis/base.py:291`).

⭐ 2.1b added one landmark file: `orpheus/numerics/basis/base.py` — `class Basis` :120 · `domain` :247 · `invariance_group` :291 (`[M]` grep, 2026-09-01; the module-scope imports it gained shift nothing below them but these).

⚠ **3.1 (2026-09-02) grew `manifold.py` by ~180 lines below `Quotient`'s
fields and `registry.py` by ~40 inside `reference`** — everything in the 2.3
block below from `ManifoldMap` down moved. `[M]` by grep at landing:
`Quotient` :574 · `orbit_coordinates` :612 · `reference` :667 ·
`Quotient` :574 · `orbit_coordinates` :612 · `reference` :667 ·
`quotient_map` :704 · `ManifoldMap` :822 · `archimedes` :940 · `barycentre`
:1010 · `_ambient_columns` :1084 · `_all_coordinates` :1099 ·
`_sphere_mod_so2` :1136 (the function-scope `LEGENDRE` import :1194) ·
`_sphere_mod_mirror` :1255 · `_mod_trivial` :1378 · `_ORBIT_CATALOGUE`
:1435 · `AngularSymmetry.support` `registry.py:875`, `.reference` :910 (the
`isinstance` guard :964, the read :965, the bare-sphere arm :979), the V-stage reader :1302 ·
`TestTheEntryCarriesItsQuotientMapAndItsReference` `test_manifold.py:949` (the module-scope
mint gate :1164) · `test_the_reference_is_READ_off_the_catalogue_entry` `test_registry.py:1315`.

**Landed and callable** (`file:line` re-derived by AST at `14c37aa9`,
2026-09-02, after 2.3. ⚠ `manifold.py` grew by ~290 lines at 2.3 — the arrows
section sits BETWEEN the constants and the catalogue, so everything from
`_catalogued_quotient` down moved by that much; `measure.py`, `symmetry.py`,
`directional.py`, `rules_product.py` moved by a few lines each):
`class Manifold` `manifold.py:119` · `Sphere` `:231` · `Interval` `:284` ·
`Ball` `:460` · `FundamentalDomain` `:492` · `Quotient` `:562`, its `by`
field `:580` · ⭐ `AXIS_INDEX` `:738` **NEW** · ⭐ `class ManifoldMap` `:742`
**NEW** (`__matmul__` `:835`) · ⭐ `archimedes` `:860` **NEW** · ⭐
`barycentre` `:930` **NEW** · `_catalogued_quotient` `:1005` ·
`_sphere_mod_so2` `:1036` · `_sphere_mod_mirror` `:1132` · `_mod_trivial`
`:1229` · `_ORBIT_CATALOGUE` `:1270` · `ambient_dim` `:1289` ·
`Basis.domain` `basis/base.py:247` · `Basis.invariance_group` `:291` ·
`Basis.space` `:370` · `IndicatorBasis` `basis/indicator_basis.py:108` ·
`_evaluate_real_sh` `basis/spherical_harmonic_basis.py:550` ·
`DiscreteMeasure.support` `measure.py:311`, `.space` `:353`,
`.quotient_group` `:414`, `.phase` `:437`, `.__mul__` `:612`, ⭐
`.pushforward` `:840` (takes a `ManifoldMap`; production callers
`:1192` inside `quotient()` and `rules_product.py:496`), `.on_orbit_space`
`:1013`, `.quotient` `:1056` · `class SO2` `symmetry.py:265`,
`SubgroupOfO3.SO2` (ctor) `:489`, `rotation_axis` `:496`, `_so2_contains`
`:907`, `_embedded_nodes` `:951` (⭐ reads `barycentre`), `_polar_axis_of`
`:1004`, `_is_axis_supported` `:1179`, `candidate_groups` `:1642` ·
`AngularSymmetry` `registry.py:826`, `.support` `:876`, `.reference` `:911`
(⚠ the twin 3.1 absorbs; reader `:1272`), `.admits_domain` `:951` ·
`gauss_legendre_on_polar_orbit` `rules_1d.py:155` · `Quadrature.quotient`
`directional.py:548`, `axis_cosines` `:294`, `mean_axis_cosine` `:327`,
`spherical_harmonics` `:623`, `angular_frame` `:678`,
`_harmonic_frame_measure` `:730` (the forgery line is now `:797`),
`Quadrature.gauss_legendre` `:870` · `spherical_product` `rules_product.py:332`
(⭐ the chart pushforward at `:496`), `_derived_product_group` `:233` ·
`ReferenceMeasure` `exactness.py:152`, `UniformMeasure` `:181`,
`UNIFORM_ON_SPHERE` `:237` · `LEGENDRE` `generating_measure.py:480`.

⛔ **The block below is the 2.1b-era list, kept per §3 — every `manifold.py`
line below `Ball` in it is STALE by ~290 lines.**

**Landed and callable** (`file:line` re-derived by AST at `4d946988`,
2026-09-02. ⚠ `measure.py` and `symmetry.py` moved at 2.1b and at this
compaction — docstring fixes only, every symbol below re-derived; `manifold.py`
did not move. ⚠ **`manifold.py` moved TWICE on 2026-09-01** — the `quotient`
body moved into a memoised module function at 2.4 (everything from `Ball`
down shifted up), then the module docstring grew by 8 lines at the archivist's
findings (everything shifted down). `measure.py`'s `phase` prose grew;
`registry.py` grew; `symmetry.py` gained a class and two functions. Do not
trust a copy from an earlier record):
`class Manifold` `manifold.py:113` · `Sphere` `:225` · `Interval` `:278` ·
`Ball` `:454` · `FundamentalDomain` `:486` · `Quotient` `:556`,
its `by` field `:574` · ⭐ `_catalogued_quotient` `:722` **NEW** (the memo) ·
`_sphere_mod_so2` `:753` (⭐ now reads the axis) · `ambient_dim` `:1006` ·
`Basis.domain` `basis/base.py:247` · ⭐ `Basis.invariance_group` `:291` **NEW** (2.1b) · `Basis.space` `:370` ·
`IndicatorBasis` `basis/indicator_basis.py:108` ·
`_evaluate_real_sh` `basis/spherical_harmonic_basis.py:550` ·
`DiscreteMeasure.support` `measure.py:307`, `.space` `:349`,
`.quotient_group` `:410`, `.phase` `:433`, `.__mul__` `:608` (the tensor product), `.pushforward` `:836`
(⚠ `new_space` **required**; its one production caller `:1161`), ⭐ `.on_orbit_space` `:994` (2.4) ·
⭐ `class SO2` `symmetry.py:265`, `SubgroupOfO3.SO2` (ctor) `:489`,
`rotation_axis` `:496`, `_so2_contains` `:907`, `_embedded_nodes`
`:951` (⭐ axis-aware — the BARYCENTRE map, §V.5i(b)), `_polar_axis_of` `:993`,
`_is_axis_supported` `:1173` (⭐ takes `axis`), `candidate_groups` `:1636` ·
`AngularSymmetry.support` `registry.py:876` (⭐ now `SPHERE.quotient(spent)`),
`reference` `:911` (⚠ the pushforward-measure TWIN, §V.5i(e)), `admits_domain` `:951` ·
⭐ `gauss_legendre_on_polar_orbit` `rules_1d.py:155` **NEW** (the slab's
adopter) · `Quadrature.gauss_legendre` `directional.py:861` ·
`Quadrature.axis_cosines` `directional.py:294` (STRICT),
`mean_axis_cosine` `:325`, `spherical_harmonics` `:621`,
`angular_frame` `:676`, `_harmonic_frame_measure` `:728` (the ERR-080 fiction,
1-D arm only — the forgery line is `:786`) · `spherical_product`
`rules_product.py:332` (the Archimedes composition at `:527`).

⭐ **Two `SO2` spellings are now landmines of their own**, for any scratch
script written before 2026-09-01: `SubgroupOfO3.SO2` is a **constructor** —
`SubgroupOfO3.SO2("x")` — and a bare `SubgroupOfO3.SO2` is a bound method
that fails at the consumer (`.name`, `.is_invariant`); and a catalogue key
`(Sphere, "SO2")` no longer exists — the keys are `SO2_x`/`SO2_y`/`SO2_z`.
`gauss_legendre_on_mu(n).support` is still `COSINE_INTERVAL` (it is the
product factor); the slab's rule is `Quadrature.gauss_legendre(n).measure`
or `gauss_legendre_on_polar_orbit(n, "x")`, on `S^2/SO2_x`.

## Landing ledger

| item | what | commit |
|---|---|---|
| 0.7 / 0.7b | the strict-xfail gate (3 rows) + **ERR-080** minted | `ae4dbc1f` |
| 0.3 | `metric.py`'s rcond comment, re-derived | `ae4dbc1f` |
| 0.4 / 0.5 | the two prose surfaces that argued against suspecting the defect | `ae4dbc1f` |
| 0.2-census | the full-suite phantom census, **discharged** | `ee60010e` + `ccda0e61` |
| — | five execution-time claims of mine, refuted and kept in place | `324137ec` |
| D0.7 / 2.0 / V.5 | `Manifold` ruled; `Basis.domain` restored; the minting inventory | `a0783c2a` |
| D0.1 / 3.4b / III.4b | the embryo; `Descent` ruled; the four-kernel table restored | `43dd7ee5` |
| D0.1 | the engine is **deferred, not refused** — the embryo is its SEED | `60c370ae` |
| XIV | compaction point 0 | `ac30598e` |
| — | the fold basis's "100 % at `L = 0`" refuted and reconciled | `52880c7a` |
| §V.5b | **2.0a's sizing REFUTED** at its own opener — three counts, three questions | `cfd98cc8` |
| §V.5c | 2.0a's real ground measure — the algebra ships already, as string concatenation | `212dfa0b` |
| **2.0a** | ⭐ **`Manifold` MINTED** — 9 variants, the S²/SO(2) catalogue entry, 40 tests | `b8c05d16` |
| — | `plan-authoring`: the AST import census drops relative imports | `fc71a84e` |
| — | numerics green as a gate; the delta reconciles exactly | `eeccc51d` |
| **2.1-W** | ⭐ **the fold basis is falsifiable at the EIGENVALUE tier** — 3 L1 rows | `daa783c5` |
| **2.0a docs** | the point-set layer's corpus page + three corrections to the mint | `fba4205a` |
| §V.5d | the 2.0c/2.0d opener — the step is BLOCKED, 2.0d dissolves, the order changes | `991097fb` |
| **1.1 σ_y** | ⭐ **`S²/σ_y` derived + `Quotient`'s TWO coordinate systems** (user-ruled) — `Ball`, `FundamentalDomain`, the construction invariant, `singular_stratum` retyped | `b55bba56` |
| **1.1 docs** | the corpus pass — and it REFUTED my dimension-coincidence claim, plus 3 code defects | `bf9296a1` |
| XIV | compaction point 2 — the exit path unblocked, and its ORDER | `4a5ac108` |
| **0.2** | ⭐ **the accessor SPLIT** — `axis_cosines` (coordinate, strict) vs `mean_axis_cosine` (the orbit mean, zero DERIVED). Found a THIRD meaning in the zeros, single-sourced `spherical_harmonics` onto the frame, collapsed `_quadrature_axis`'s ladder + its defaulted-`getattr`. `[M]` red loop 29 tests / 2 trees — the §6b set missed the EMBEDDING category | `ce46181c` |
| §V.5e | ⭐ **2.0c's opener** — unblocked; the tracker's sizing counts the wrong side; `phase` is stringly-typed dispatch; and `[M]` a LIVE defect (`folded_product(…).measure.phase` RAISES) | `765757a0` |
| **2.0c** (+2.0d) | ⭐⭐ **a measure names the point set it lives on** — `Space = str` and its six tags RETIRED across all six implementors; `phase`, the `L2[…]` name, the fold and the tensor product all DERIVED from the manifold. `[M]` 10 of 10 names bit-identical; a LIVE defect fixed (`folded_product(…).measure.phase` raised); 4 defects found that the string could not express, all green before; ⛔ the mutation battery read 2.1's handoff repair as BLIND and its witness now reds 6 of 6 | `025834f5` |
| **2.1** | ⭐ **no space claims an identity it does not have** — `Basis.domain` abstract, `IndicatorBasis` takes `partition_of`, `ambient_dim` promoted public. `[M]` the energy/spatial `L2[coarse_cells_R1]` collision (`==` AND hash-equal) is unspellable. 13 gates; 7-arm battery, every arm a distinct red set; §6b set was **23**, not 18 | `c461fe8d` |
| **0.1a + 0.1c** | ⭐ **the frame stops forging its own domain** — `_harmonic_frame_measure` routes the rule's measure; the 1-D fiction is named + self-retiring. `[M]` route gate 0→10 of 12; §II.8's *three* losses reversed; 5-arm mutation battery, every arm's red count reconciling to zero unexplained units | `2c1a06b1` |
| **2.4** | ⭐⭐ **the slab says what space its ordinates live on** — `SO2(axis)` parameterised (user-ruled: the tree has TWO poles), `S^2/SO2_x` declared by `gauss_legendre_on_polar_orbit`, `on_orbit_space` minted, the registry derives its domain. `[M]` the `[-1,1]` angular/spatial collision unspellable; Part IV obstacle 1 answered by derivation; stage 0 refuses the chart-level and wrong-axis rules; 4 trees rc=0, +22 reconciling; a pre-existing 9 s lattice-walk hot spot removed | `17501245` |
| **2.1b** | ⭐ **a basis states the symmetry its functions HAVE, by naming what they EAT** — `Basis.invariance_group` DERIVED from `domain` (`Quotient.by` / `Trivial` / `None`), `@final`, zero subclass edits; the fold's two halves read ONE object; ERR-080's pairing is a lattice verdict (`Trivial ⊉ SO2('x')`), refused by nothing yet. 11 gates, 8-arm battery (1 blind by construction), 3 stale docstrings fixed on sight; numerics 2679 → 2690 | `9b4a4d9c` |
| **2.3** | ⭐⭐ **a measure's support is READ off the map that built it** — `ManifoldMap(domain, codomain, apply)` with `@` composition; `archimedes(axis)` and `barycentre(orbit_space)` memoised; `pushforward(phi)` reads `phi.codomain`, refuses a map out of the wrong point set, `new_space=` retired (1 + 8); `spherical_product` IS the tensor product pushed along the chart; the fold's retraction and `_embedded_nodes` typed; ERR-080 restated as the barycentre map with a FORGED codomain (arm untouched, retires at 3.4). `[M]` bit-identical 60/5/15; eleven-arm battery; numerics 2690 → 2711; matrix 10595 → 10616; reference measure deferred to 3.1 with the §6d cycle measured (5 of 5) | `5ec3a00a` (+ stamp `14c37aa9`) |
| **3.1** | ⭐ **an orbit-space entry carries ALL of its derivation, the quotient map and the pushforward reference included** — `Quotient.orbit_coordinates` (required; the surviving invariants as a map on ambient coordinates) + the DERIVED `Quotient.quotient_map` (codomain the ENTRY), `Quotient.reference` (`LEGENDRE` on the three axial entries; honest `None` on σ_a and on `M/{e}`), populated INSIDE `_sphere_mod_so2` by a function-scope import; `AngularSymmetry.reference` READS the entry, its `LEGENDRE` import gone. `[M]` π∘φ_a = pr₁ bit-exact 12/12; β_a∘π_a the axial projection 3/3; π H-invariant on 4 groups; the numeric map == `lambdify` of the recorded invariants 5/5; §6d on a shadow copy 7/7 alive vs 0/7; nine-arm battery, none blind; numerics 2711 → 2739, matrix 10616 → 10644 | `67e38605` |
| **2.5** | ⭐ **the angular moment space is READ off the frame, never minted from `L`** — `TruncatedBasis` (Protocol) is the harmonic family's surface and `HarmonicFrame`'s door; the Λ ends, the fission/(n,2n) ℓ=0 ends, the moment-flux head and `truncate` all read the bound basis's space (7 producers retired, 2 doors widened); A-R1 binds the CONTINUUM space (`[M]` `Λ* = Λᵀ` exactly under it, 33/33; the dressed one moves Λ* on 10/33), A-R2 derives the field head through the mesh's quadrature. `[M]` slab L=0/1/2 flux `array_equal`; 33-row metric identity; route gate with a FOREIGN basis; 11-arm battery (control 34, each producer revert → A1, the fork → A2b only), none blind; 13 trees rc=0, 10454 = predicted | `58ad7c28` |
| **FUSED** | ⭐⭐ **THE FIX — a 1-D rule binds the Legendre basis on its orbit space; ERR-080 CLOSED** (0.1b + 0.6 + 2.2-G0 + 3.4 + 3.4b): `LegendreBasis` (the bit-matched spelling), `LegendreSpace` + `MomentHead`, the entry's probe, `Descent`, G0 as the descent arrow with the table pulled back along it, the forgery deleted, 0.6's refusal, the carriers reading their head, `OverlapBasis`'s misdeclared domain caught by G0 and fixed. `[M]` LEG A `+4.000000000000` at L = 0..2; slab L ≤ 1 `array_equal`; 135 sphere arrays bit-identical; 19-arm battery, one declared-weak arm; the 13-tree gate finished after the commit: **13 of 13 rc=0, 10556 = predicted**, reconciled file-by-file in the *FUSED EXECUTED* section | `5436184e` + follow-up `26151a8d` (gate record, typing narrowing) |
| **1.9** | ⭐ **an orbit space is named by its STABILISER — `SubgroupOfO3.O2(axis)`, and a basis on it declares the FULL group** (#432): `O2(a)` = the stabiliser `{g : g ê_a = ê_a}`, the axial lattice relations COMPUTED from realizations (`SO(2)_a = O(2)_a ∩ SO(3)` by composition), exact O(2) invariance, both-component `generic_images`; ONE catalogue entry `S^2/O2_a`, `SPHERE.quotient(SO2(a))` REFUSED with the theorem (user-ruled over normalising); the elegance round's `orbit_stabiliser` puts the refusal on `Quotient.__post_init__` and at the door, retires the three decoy keys and the second runtime import; the registry spends `O2("x")`; `LegendreSpace` names itself off the basis's domain. `[M]` the Legendre basis is ADMITTED on a σ_y-fold (ℓ≥1 isotropic moments ≤ 1.4e-15 at L = 2/4/6) and refused on the x-fold; 19-arm battery all biting; 13 of 13 rc=0, 10604 = predicted (fused 10556 + 48 numerics), 60 min 27 s | `c1d53206` |
| **2.2b** | ⭐ **the Γ-slot — a geometry admits a FOLD of its domain, and the owed symmetry is asked ON the orbit space** (#429): `Quotient.induced_action` (well-defined iff the motion normalises the quotienting group — `SubgroupOfO3.is_normalised_by` / `normalises`, EXACT for the continuous members through the rotation generator and one representative per coset), `Quotient.lift` / `ambient_representatives`, ONE closure `_orbit_space_closure` read by `is_invariant`, `orbit_certificate` and `induced_permutation` (which `Quadrature.ordinate_permutation` reads — on a fold σ_y is the identity permutation of the stored ordinates); stage 0 = `quotient_onto(D, X)` exists AND `Γ ⊇ spent_group(D, X)`; the spent-group door; `_invariance_on_points`/`_polar_axis_of` retired; II.11's 1-D shape refusal gone. `[M]` cylinder and cartesian2d ADMIT `folded_product` at both stages (both False before); 3 of 230 answers moved, all on the fold; 148 gate rows; 14-arm battery (12 bite, the bijectivity leg uninstallable, the short circuits theorems); 13 of 13 rc=0, 10752 = predicted (1.9's 10604 + 148 numerics), 59 min 27 s | `<HASH>` |

**Branch** `fix/angular-phantom-support`, pushed, ⚠ **nothing merged.**
⛔ No commit COUNT is recorded here — it is the field guaranteed to rot. Run
`git rev-list --count main..HEAD`, and `git merge-base --is-ancestor <hash> HEAD`
for merge status.

## ⛔ Corrections that SUPERSEDE older text in this file

Read these before quoting any earlier section:

| what was said | status |
|---|---|
| 0.2: *"`axis_cosines(i ≥ dim)` raises"* | ⛔ **REFUTED by the census** — breaks 3 legitimate flow-question consumers. §II.14 |
| 0.1: a single landable item | ⛔ **SPLIT** into 0.1a/b/c; the lift `[-1,1] → S²` does not exist as a map. §II.9 |
| 2.2: *"just wiring"* | ⛔ **REFINED** — the predicate is correct over an INCOMPLETE vocabulary; it would refuse the shipped cylinder. §II.10 |
| §II.8's 0.1c hazard | ⛔ **REFUTED** — `[M]` 0 of 145 089 composability checks carry such a space. §II.15 R2/R3 |
| `Basis.support` (ruled 2026-08-31 am) | ⛔ **REFUTED same day** — collision unreachable, and `support` is false for a basis. §III.10 |
| *"the engine is refused"* | ⛔ **DEFERRED, not refused** — and the embryo must be its seed. D0.1 |
| *"same INTERFACE, second backend"* | ⛔ **too weak** — the requirement is on the DATA MODEL. D0.1 |
| 2.0a: *"29 `Space` refs, 39 `support=`, 45 `.support`"* | ⛔ **REFUTED at the opener** — by AST 87 / 62 / 8, and 53 % of the surface is `tests/`. §V.5b |
| 2.1-W: *"zero committed witnesses"* | ⛔ **FALSE — there were four**, all OBJECT-tier, none through `solve_sn`. §II.16(a) |
| 2.1-W: the *18.5 %* catcher figure | ⛔ **does not reproduce** — 57 % on its own stated fixture. §II.16(b) |
| 2.1-W: *"reds 0 of 1913"* | ⛔ **pytest never RAN** — rc=4, 0 collected, zero on BOTH scanners. §II.16(c) |
| 2.0a-R: *"all 8 derivation fields"* | ⛔ **6 of 8** — chart + pushforward wait on `Chart` (2.3). The ruling block |
| §V.5c: *"`quotient` does no lookup and no check"* | ⛔ **HALF false** — the MEASURE is gated; the TAG is not. §V.5c(c) |
| 2.0d: *"`gauss_legendre(8).measure.quotient_group is SO2`, derived"* | ⛔ **UNREACHABLE as written** — that measure never goes through `.quotient()`; the `SO2` answer is **2.4**'s declaration. §V.5d(d) |
| 2.0d as a separate landing | ⛔ **DISSOLVED into 2.0c** — `Quotient.by` already IS the group. §V.5d(c) |
| 2.0c as the next step | ⛔ **BLOCKED** — no `S²/σ_y` catalogue entry; 1.1 precedes it. §V.5d(a) |
| §V.5c(e): *"two spellings of one quotient, both shipped"* | ⛔ **all five odd literals are in `tests/`**; `orpheus/` has zero, and `Quotient.name` reproduces the production tag bit-identically. §V.5d(f1) |
| 2.0c's `indicator_basis` done-when clause | ⛔ **MOVED to 2.1** (§6b) — it needs 2.1's ctor field to be fixed *by construction*. §V.5d(e) |
| 2.1: *"three manifolds over 5 sites"* | ⛔ **18 sites** (4 `orpheus/`, 14 `tests/`); the *three* was right. The 2.1 row |
| §V.5d(b): *"the SO(2) entry cannot expose the fork because chart and section coincide in dimension"* | ⛔ **REFUTED BY MY OWN FIX** — `[M]` they coincide in BOTH, necessarily, since `Quotient.__post_init__` gates it. §V.5d(b) |
| 2.0a-R: *"`AngularSymmetry.support` is an un-briefed Pattern-2 twin"* | ⛔ **NOT a twin** — it computes `S²/G⁰`, the CONTINUOUS isotropy; a mirror is in the discrete residual Γ, a row it cannot answer. **2.2 still owes the second slot.** The 2.0a-R block |
| 2.1 pre-flight: *"18 `IndicatorBasis(` sites is the right §6b set"* | ⛔ **23** — it cleared `OverlapBasis` as *"inherits the FIELD"* (true) and never asked whether it has CONSTRUCTORS of its own: `OverlapBasis(` ×3, `from_indicator`, and `_DenseTrial((edges,))`, a test-local subclass on the same LINE as a counted call. None contains the substring `IndicatorBasis`. §6b's sixth spelling |
| 2.1 row: *"across 26 space-comparing sites"* | ⛔ **STRUCK** — inherited, never re-measured, and the wrong predicate anyway: the repair is at the 5 PRODUCERS of the false name, not at the comparers |
| the 2.1 done-when *"`domain is SPACE_SPHERE`"* | ⛔ pre-`Manifold` text — `SPACE_SPHERE` is the `str` `"S^2"`, and `[M]` `Sphere() == SPHERE` is `True` while `is` is **`False`**. Built as written it would have re-introduced `Space = str` |
| 2.0a / 2.0b: `✅ de29bcc6` | ⛔ **DANGLING** — `[M]` amended 39 s later into `b8c05d16` and now orphaned (no branch contains it, so it will be GC'd). The ledger already said `b8c05d16`; the plan disagreed with itself. Re-pointed 2026-09-01 |
| 2.0c row: *"`DiscreteMeasure.support: str → Manifold`"* | ⛔ **names ONE implementor of SIX.** `[M]` the `Space` alias is annotated at 7 sites over 5 classes — incl. the `ReferenceMeasure` **Protocol** — and a 6th implementor (`AngularSymmetry.support`) is spelled **`-> str`**, the alias's EXPANSION, so no alias census can return it. It is `DiscreteMeasure.support`'s comparison partner in stage 0, and `Manifold == str` is silently `False`. §V.5f(a) |
| §V.5e(c): *"`phase` becomes a property of the MANIFOLD"* | ⛔ **SUPERSEDED by the user's ruling** (2026-09-01): dispatch on the manifold's TYPE, `phase` stays on the measure. `Phase` is a transport taxonomy and `Manifold` is pure geometry; a `Sphere` is angular *because transport's direction variable lives on it*, not intrinsically. §V.5f(e) |
| `manifold.py`'s *"the shipped members, under their **retired** tag names"* | ⛔ **ASPIRATIONAL when written** — `[M]` none of the six was retired; `generating_measure.py:398` compared against `SPACE_INTERVAL_M11` and 9 production sites imported the family. ✅ **True as of 2.0c**, and the comments now say *replaced*, in the past tense they had already claimed |
| 2.1's `test_d6`: *"2.0c must come back here and cannot resolve it by accident"* | ⛔ **it DID resolve it by accident.** The repair (`support=basis.domain`) makes the divergence unspellable, so `test_d6` went green for a different reason — and `[M]` the mutation battery read the repair as **BLIND tree-wide**. A prediction does not enforce itself. Witness written; arm now reds 6 of 6 |
| `pushforward(new_space=None)` defaulting to `f"φ_*({support})"` | ⛔ **RETIRED** — the image of a manifold under an arbitrary map is a manifold nobody derived, and naming it does not make it known. `[M]` 7 of 8 call sites already passed the target; the one that did not was a test asserting the fabrication. Its gate flips to the refusal |
| `test_measure.py`'s *"the slab embedding `(μ,0,0)`"* fixture on `"[-1,1]^slab"` | ⛔ **ERR-080's defect class INSIDE A GREEN TEST.** Those points are on no manifold with a σ_y action; the string fabricated `"[-1,1]^slab/sigma_y"` and `Manifold.quotient` refuses. Replaced by the genuine σ_y-fixed set (the `xz` great circle). The same test also **contradicted its own docstring** — header: *"re-folding is ill-posed"*, body: re-folds and asserts idempotence |
| `test_angular_phase_is_symmetry_derived_not_support_tagged` | ⛔ **its claim INVERTED.** It asserted that stripping `invariance_group` makes `phase` RAISE, because `"S^2"` was an inert string no branch read. `Sphere` is a *reason*, so the measure stays angular — which is what lets the shipped fold answer at all. Renamed + re-argued in place |
| §V.5e(e): *"`support='cells'` has no producer — a dead arm"* | ⛔ **REFUTED same day, by its own author.** Production yes; the real producer is a `parametrize` list at `tests/numerics/test_measure_phase.py:67`, which MY OWN census had found as a non-literal row and I failed to read. ⟹ *production-unreachable with a synthetic witness*, not dead |
| §II.13: *"`folded_by` is structurally unavailable for a continuous quotient; the tree holds TWO quotient records, a group object and a string tag"* | ⛔ **REFUTED by 2.0c** — `[M]` there is no string tag any more, and `Quotient.by` handles the continuous case identically (`SPHERE.quotient(SO2).by` **is** `SubgroupOfO3.SO2`), with `measure.quotient_group` deriving from it for discrete and continuous folds alike. The two records are ONE mechanism and the structural obstacle does not exist. ✅ What survives is 2.4's whole content: the slab is still the unwired half — now a one-line declaration against a working catalogue, not a research question |
| 2.1's headline: *"no space claims an identity it does not have"* | ⛔ **FALSE as a UNIVERSAL** — `[M]` 2026-09-01 (§V.5g) an 8-node slab ANGULAR space and an 8-node SPATIAL space on `[-1,1]` are still `==` AND hash-equal. 2.1 closed the energy/spatial pair; the class survives wherever two ROLES honestly share a point set. ⭐ And 2.0c could not close it — the retype made the supports *honest*, and they are honestly identical. Re-tensed to what was measured; closing it is **2.4** |
| 2.4's summary tracker row: *"…; **full-suite phantom census**"* | ⛔ **STALE — the census MOVED to Phase 0 on 2026-08-31 and is DISCHARGED** (`ee60010e`+`ccda0e61`); it was 0.2's denominator, not 2.4's. The DETAILED row (Part XI) recorded the move and the summary row did not, so the plan disagreed with itself — the same class as the `de29bcc6` dangling hash, found by the same sweep |
| 2.4 filed as wiring (*"the slab rule declares its quotient group"*) | ⛔ **RE-CLASSED at its pre-flight** — it repairs a **LIVE** defect of 2.1's kind (the `[-1,1]` space collision above), and it is the ONE step that RENAMES a shipped space, so its §6b set is the **23** comparison sites, not the 3 literal pins. §V.5g |
| 2.0c row: *"54 literals + 5 f-strings, 94 % in `tests/`"* | ⛔ **counts the PRODUCERS.** `[M]` the two headline numbers reproduce EXACTLY (54, 94 %) and "5 f-strings" is **4**; the set a retype breaks is the **16 CONSUMERS** (15 `orpheus/`), never counted. §V.5e(b) |
| `DiscreteMeasure.phase` is total on the shipped rules | ⛔ **FALSE — `folded_product(4,8).measure.phase` RAISES** (`[M]` support `S^2/sigma_y`, `invariance_group=None`, no quotient arm). A live defect, latent only because `[M]` `.phase` has **0** production consumers. §V.5e(f) |
| *"`Basis.domain` 0 of 6 subclasses"* | ⛔⛔ **the ORIGINAL 6 WAS RIGHT — my "correction" to 5 was itself wrong, and shipped at `4a5ac108`.** `[M]` 2026-09-01 by RUNTIME `__subclasses__` after importing every `orpheus` module: **4 direct, 6 recursive**, and `0 of 6` define `domain` or `invariance_group`. My static AST census undercounted (it matched only `ast.Name` bases). §II.15's own *"7 transitive subclasses (6 production + `_DenseTrial`)"* already said 6, so the plan was self-inconsistent and I resolved it toward the wrong side. ⟹ **for a subclass census use RUNTIME introspection, not AST** — inheritance is a runtime relation and an AST pass cannot see re-exports, aliased bases or qualified base names |
| 0.1/L1's site `directional.py:587-593` | ⛔ **DRIFTED** — `[M]` 2026-09-01 the fabrication is at **`:624`**; `angular_frame` begins at `:580` |
| §II.4 + §II.8: the forgery destroys the **support tag** | ⛔ **INCOMPLETE — it destroys THREE fields.** `[M]` `invariance_group` (10 of 12 rules → 0 of 12 frames) and `exactness` (same split) were dropped too, so `frame.measure.phase` RAISED on all 12. §II.8's correction block |
| §II.15 R4 (*a value gate cannot be the keystone*) | ✅ **RE-DERIVED INDEPENDENTLY at 0.1a**, on my own gates rather than inherited: under the exact pre-carve mutation the value gate `test_q8_5` passes **8 of 8** while the route gate `test_q8_1` reds **8 of 8** |
| §V.5g(d): *"size 2.4 by the **23** comparison sites"* | ⛔ **the 23 had no stated predicate, and no comparison count was the set.** `[M]` by AST under a stated predicate there are **29**; 26 compare a field's space against one derived from the SAME measure and rename together; the 3 that can see a declaration were all handled by design. The set the red loop found was **5 test-local FACTORIES** in `test_registry.py` that built chart-level rules the stricter stage 0 now refuses — a §6b spelling no comparison census returns. 2.4 EXECUTED (d) |
| §V.5g(d): *"3 literal `L2[[-1,1]]` pins, 1 production f-string + 2 docstrings"* | ⛔ **0** in `orpheus/`, `tests/` and `docs/` by fixed-string grep — the three were in THIS PLAN. The earlier grep used a regex with `\[` under ugrep, which is the nexus-tools landmine; re-run with `-F` |
| `measure.py`'s 2.0c comment on the `phase` fallback arm: *"tracker 2.4 makes the slab declare `SPHERE.quotient(SO2)` … and this arm becomes unreachable"* | ⛔ **FALSE, by my own hand, and the §6c check §V.5g(c) prescribed is what caught it.** Unreachable for the SLAB, yes; `[M]` the chart-level `gauss_legendre_on_mu(8)` — public, shipped, the product rules' polar factor, which cannot declare an axis — reaches it and answers `"angular"`. Arm kept, comment corrected, witnessed in `test_a3`. A prediction in a production comment enforces itself no better than one in a gate's docstring (the `test_d6` row) |
| `test_symmetry.py`: *"a polar marginal is NOT SO(2)-invariant"* (the 2026-08 inversion of *"trivially SO(2)-invariant"*) | ⛔ **INVERTED A SECOND TIME, and both earlier forms were right about what they measured.** With the axis named, `[M]` the marginal IS invariant under `SO2("x")` (its own spent group) and NOT under `SO2("z")`; the bare tag could only ever report one of the two. Re-argued in place, with the z-declared marginal answering the opposite way |
| the 2.0c `phase` docstring bullets (*"spatial — iff support is a spatial tag (`"spatial_…"` / `"cells"`)"*, *"angular — iff `invariance_group` is set"*) | ⛔ **present-tense-false since 2.0c** — they described the STRING dispatch 2.0c retired while the body underneath matched on manifold types. Re-tensed at 2.4 in passing; nothing tested prose |
| Part IV obstacle 2 (*"the slab's quotient group is recorded nowhere"*) | ✅ **REMEDIED** — `folded_by` retired at 2.0c, the slab's group recorded at 2.4 as `Quotient.by` read through `measure.quotient_group`, one mechanism for the discrete fold and the continuous marginal |
| §V.5h(f): *"2.1b ADDS an abstract member … then the §6b set is every subclass constructor including test-local ones"* | ⛔ **VOID at execution** — it landed as a CONCRETE derived property; the §6b set was EMPTY (zero subclass edits). The tracker's own MEANS (*a field, answered by six subclasses*) dissolved into `domain.by` the way 2.0d's did into `Quotient.by`. *2.1b EXECUTED* |
| tracker 2.1b's row: *"which today has **neither** side"* | ✅ **REMEDIED, and was already half-false when the 2.4 record was written** — the measure side became real at 2.4 (`quotient_group`), the basis side at 2.1b; corrected in place |
| 2.3's done-when: *"`grep` finds no `SPACE_SPHERE` literal at a measure constructor"* | ⛔ **§10's third shape** — `[M]` 5 `support=SPHERE` literals in production, **3 honest tabulations** (Lebedev, level-symmetric, the reference measure) the step must not touch; scope it to MAP-BUILT measures (2 of 5). §V.5i(c) |
| Part XIV's pre-2.4 pointer: *"2.3 — the typed `Chart` (also the SECTION MAP's natural home)"* | ⛔ **REFUTED by census** — `[M]` `fundamental_domain` has 0 production readers; the only section-like consumer is the forgery arm 3.4 retires. 2.3 is the quotient map + the Archimedes chart. §V.5i(a) |
| 2.3's row and every pointer: *"the typed `Chart`"* | ⛔ **the NAME was refuted at the ruling (2026-09-02)** — a chart is `M ⊃ U → ℝⁿ`; of the three maps only the inverse of the Archimedes one is a chart. User-ruled `ManifoldMap`, one frozen value type, factories for the named maps. *2.3 EXECUTED* |
| §V.5i(b), the `ManifoldMap` docstring and the archivist's brief: *"`Ball` has 0 production consumers"* | ⛔ **FALSE — the census excluded `manifold.py`**, where `Ball(2)` is the σ_y entry's `realization` (`:922`). Honest: 0 consumers OUTSIDE its module, `Ball(3)` never constructed. Caught by the archivist re-deriving it; corrected in code, plan and page. `plan-authoring` surprise log, 2026-09-02 |
| 2.0a-R's supersession row: *"chart + pushforward wait on `Chart` (2.3)"* | ⛔ **HALF discharged** — 2.3 shipped the chart's TYPE and three arrows, but `[M]` the ENTRY's own chart map `Ω ↦ Ω·ê_a` and its pushforward measure are still not fields (`Quotient` has 10 fields, 5 of D0.1's 7 derivation outputs); both wait on **3.1**, §V.5j(b)–(c) |
| §V.5i(d) and this Part's 2.3 pointer: *"`pushforward` has ONE production caller"* | ✅ true when written; **2** after 2.3 (`measure.py:1192`, `rules_product.py:496`) — do not re-quote the 1 |
| 3.1's row: *"`numerics/symmetry/` catalog … this is a re-home"* and XII.2's *"the catalog home + the probe"* | ⛔ **pre-2.0a premise** — `[M]` the catalogue is `manifold._ORBIT_CATALOGUE` (E. Placement), the probe is 3.4's (§III.8), and the package re-home is OFF the exit path. 3.1's honest scope is two fields + one map. §V.5j |
| §V.5j(e): *"(i) `chart: ManifoldMap` — codomain the entry itself"* as a stored field | ⛔ **NAME and MECHANISM refined at 3.1's opener (user-ruled)** — `quotient_map`, DERIVED from a stored `orbit_coordinates`; a frozen dataclass cannot store a map whose codomain is itself. The 3.1 EXECUTED section |
| §V.5j(c): *"The `S²/Trivial` entry's is `UNIFORM_ON_SPHERE`"* | ⛔ **user-ruled `None`** — Lebesgue on the BASE is a property of the base (`Manifold` carries no orthogonal system); the registry's Trivial arm reads the bare sphere, never the entry. The 3.1 EXECUTED section |

## Measured baselines and costs

⭐⭐ **RE-GATED at 2.4** (2026-09-01), serial `-O -m "not slow" -p no:randomly`,
**all rc=0** — the sn column is filled in by the 2.4 EXECUTED section's gate
row once its detached run reports (it is the long one; read that section
rather than this table for it):

| tree | at 2.1 | at 2.0c | at 2.4 | delta 2.0c → 2.4 |
|---|---|---|---|---|
| numerics | 2656 | 2657 | **2679** | +22 (17 new-file rows, +4 symmetry gates, +1 parametrize row — reconciled) |
| data | 237 (+2 desel) | 237 | 237 | — |
| geometry | 727 +4sk +1xf | 727 +4sk +1xf | 727 +4sk +1xf | — |
| transport | 645 +1sk | 645 +1sk | 645 +1sk | — |
| sn | 3378 +1sk +50xf | 3384 +1sk +116desel +50xf | 3384 +1sk +116desel +50xf | — (the rename moved nothing; 15 min 24 s) |

⭐⭐ **RE-GATED at 3.1** (2026-09-02, `67e38605`, same runner, all rc=0,
logs `scratch/_p31_gate_*.log` and `_p31_post_*.log`): numerics **2739**
(2711 + 28: 21 + 1 rows in the first pass, then 6 after the archivist's
findings — reconciled unit-for-unit) · ERR-080 gate **2 + 3 xf** (strict,
OPEN) · sn primitives + mesh + angular + transport/frames **579** ·
sn/operators **1274** + 5 xf · sn/sweep **919** + 1 sk + 31 xf · geometry +
`_harness` + mms-ordering **736** + 4 sk + 1 xf · root docs gates
(`test_error_catalogue_reconciles` + `test_docstring_xrefs` +
`test_layer_imports`) **399** — ⚠ a different file set from 2.3's *415 + 5
xf* row, not comparable. V&V matrix **10616 → 10644** (+28, matching the
numerics delta; sentinel count 577 → 578 for the new documented label).
`sphinx -E -W` 0/0 twice; `dead_references` 0 / 52; pyright 0 on the 4
touched files; the 8 numerics warnings are 2.3's 8. Nine-arm battery, no
arm blind (the *3.1 EXECUTED* section).

⭐⭐ **RE-GATED at 2.3** (2026-09-02, `5ec3a00a`, same runner, all rc=0):
numerics **2711** (2690 + 21: `TestManifoldMap` 14, `test_rules_product` +7,
`test_measure` net 0) · sn primitives + mesh + angular + transport/frames
**579** · sn/operators **1274** + 5 xf · sn/sweep **919** + 1 sk + 31 xf ·
geometry + `_harness` + mms-ordering **736** + 4 sk + 1 xf · root docs gates
**415** + 5 xf. V&V matrix **10595 → 10616** (+21, matching the numerics
delta unit-for-unit). `sphinx -E -W` 0/0/0; `dead_references` 0 / 52;
pyright 0 on the 8 touched files. ⚠ These are SLICES, not the tree table
above: `sn/operators` + `sn/sweep` + the primitives set are three of `sn`'s
subtrees, and the geometry number bundles two extra paths — do not compare
them against the per-tree column. data / transport / cp / derivations …
**not re-run** at 2.3 (the change is bit-identical on every measure any of
them builds, `[M]` 60/5/15).

⛔ *Prior, kept per §3:* ⚠ **No gate at the 1st 2026-09-02 compaction** — the tree was 2.1b's plus two docstring lines in `measure.py` (`4d946988`); 2.1b's row below was the baseline then.

⭐ **RE-GATED at 2.1b** (2026-09-01, same runner): numerics **2690** rc=0
(2679 + 11, reconciled unit-for-unit against section E's rows: `e1` 1, `e2` 4,
`e2b` 3, `e3` 1, `e4` 1, `e5` 1). data/geometry/transport/sn **not re-run** —
2.1b changes no behaviour any tree outside `numerics` reads (a derived
property with zero consumers; `[M]` the 8-arm battery reddens nothing outside
`test_basis_domain.py`); the slice most likely to notice
(`test_quadrature_fold` + `transport/frames` + `sn/angular`, **79**) is rc=0,
and the root docs gates read **415 + 5 xfail** rc=0. The matrix delta is in
the *2.1b EXECUTED* section.

⛔ **Prior table, kept per §3** — *"RE-GATED at 2.0c (`025834f5`): numerics
2657, data 237, geometry 727, transport 645, sn 3384"*; superseded by the
column above.

⭐ **The delta reconciles unit-for-unit against a second, independently counted
instrument.** The V&V matrix went **10555 → 10562** (+7): `+1` the new
`test_measure_phase` row (the live-defect witness) and `+6` the new
`test_loss_kernel_gauge` rows (the BLIND repair's witness) — matching the tree
gate's `+1` numerics and `+6` sn exactly. `foundation` 7531 → 7538.
`sphinx -E -W --keep-going` **build succeeded, 0 warnings**; `dead_references`
**0 dead / 52 checked / 0 undecidable** (⚠ it found **1** on the first pass — a
`:data:` xref my rename re-pointed at the retired module; nothing else in the
toolchain can see that surface); pyright **0** across all 18 touched files
(10 production + 8 test).

⛔ **Prior baseline, kept per §3** — *"RE-GATED at 2.1: numerics 2656
(`ce46181c`'s 2643 + 13 new `test_basis_domain.py` rows), … V&V matrix
10542 → 10555, all 13 in `foundation` (7518 → 7531)."* Correct when written and
superseded by the table above.
⚠ **pyright scope**, so a later reader does not compare two populations: a bare
`npx pyright` over the whole project reads `[M]` **2533** errors (it includes
`derivations/diagnostics/`, `tools/`, …). Every "pyright 0" in this plan is a
DELTA claim on the files a step touches — the right predicate, never written
down until now.

⭐ **The five trees below were gated at `ce46181c`'s tree state**, serial
`-O -m "not slow" -p no:randomly`. Every number is a GATE (rc checked), not a
census — do not compare them against the phantom-census figures, which came
from a recording plugin and are a different instrument.

| tree | result | rc |
|---|---|---|
| numerics | **2643** passed | 0 |
| transport | 645 passed, 1 skipped | 0 |
| geometry | **727** passed, 4 skipped, 1 xfailed | 0 |
| sn | **3378** passed, 1 skipped, 116 deselected, 50 xfailed | 0 |
| derivations | 1637 passed, 13 skipped, 67 deselected, 11 xfailed | 0 |
| rootfiles (3 modules) | 353 passed, 5 xfailed | 0 |

⚠ `derivations` (36 min) was gated on the FIRST 0.2 run, which already carried
every production change; the repairs after it were test-only and `[M]`
`tests/derivations` has **0** references to the touched harness module. Stated
so a re-reader does not mistake it for a stale number.

* **numerics reconciles to zero unexplained units across both landings**:
  2594 (pre-session) + 37 (0.1a/0.1c: Q8.1 ×8, Q8.2 ×8, Q8.3, Q8.4 ×2,
  Q8.5 ×8, Q8.6 ×10) = 2631, + 12 (0.2: Q9.1, Q9.2, Q9.3 ×8, Q9.4, Q9.5)
  = **2643**.
* **Verification matrix**: 10493 → 10530 → **10542**, i.e. +37 then +12, all
  `foundation` (7469 → 7506 → 7518) — the same units, counted by a second
  instrument.
* **Sphinx `-E -W --keep-going`**: `build succeeded`, **0 warnings**, run
  TWICE (the second after the ERR-080 edit, because the first predated it).
  **`dead_references`**: **0 dead / 52 checked / 0 undecidable**.
* **pyright (CLI, not the streamed LSP)**: **0 new errors**. 0.1a fixed 2
  pre-existing ones in `test_q7_4` by narrowing; 9 remain in
  `test_bc_universal_invariants.py`, `[M]` identical at HEAD and on lines this
  work does not touch.
* ⚠ **The 13-tree GATE baseline is still UNMEASURED.** Exit predicate 5 stands.
  Last full 13-tree gate: CS4c step 4's `10106 / 0 / 19 sk / 66 xf`.
  `[M]` measured green at `025834f5` (2.0c's exit gate, serial, rc=0 each):
  **numerics 2657 · data 237 · geometry 727 · transport 645 · sn 3384** —
  **five** of the thirteen. ⭐ **At 2.4: SEVEN of thirteen** (the five above
  re-gated — numerics **2679** — plus cross_method 81 and diffusion 113, all rc=0). ⛔ The previous version of this bullet said "six …
  the remaining seven (diffusion, **data**, …)" and listed `data` on BOTH
  sides; corrected by enumerating `tests/*/` rather than recalling it.
  Not measured at this HEAD: **cp, cross_method, derivations, diffusion,
  homogeneous, mc, moc**, plus the eight root-level `tests/test_*.py` modules
  (of which `test_error_catalogue_reconciles` + `test_docstring_xrefs` **were**
  run green at 2.0c — 50 passed — because the corpus pass touched them).
  ⟹ **best spent immediately before the fused step #5**, where it becomes
  load-bearing — not now. ⭐ **And that is NOW: launched detached at
  `2f294ef1` (`scratch/_p31c_full_gate.sh`), denominators PREDICTED first —
  10417 collected across the thirteen (§V.5k(f)); the measured rows are
  recorded below this bullet when the run reports, or in
  `scratch/_p31c_full_gate.log` if it outlives the session.** ✅ **It
  reported 2026-09-02 08:45 — 13 of 13 rc=0, 10417 = predicted; the table is
  the RE-GATED block below.** ⭐ At 2.3 (`5ec3a00a`) the slices in the RE-GATED
  block above were run (numerics whole; sn by three subtrees; geometry +
  harness) — still not the 13-tree gate.

⭐⭐ **RE-GATED — THE FULL 13-TREE GATE, exit predicate 5, MEASURED 2026-09-02
at `2f294ef1`** (`[M]` the code tree is `17e0eb0e`'s — the two commits between
touch only this plan and one skill), serial `-O -m "not slow" -p no:randomly`,
per-tree denominators PREDICTED first by collect-only (`scratch/_p31c_full_gate.sh`
→ `scratch/_p31c_full_gate.log`, `_p31c_gate_<tree>.log`), **all 13 rc=0**:

| tree | predicted (collected) | measured | rc | wall |
|---|---|---|---|---|
| numerics | 2739 | **2739** passed | 0 | 29 s |
| transport | 646 | 645 + 1 sk | 0 | 23 s |
| geometry | 732 | 727 + 4 sk + 1 xf | 0 | 9 s |
| data | 237 (+2 desel) | 237 | 0 | 41 s |
| homogeneous | 50 | 50 | 0 | 2 s |
| diffusion | 113 | 113 | 0 | 2 s |
| cp | 141 (+13 desel) | 141 | 0 | 66 s |
| moc | 121 (+3) | 121 | 0 | 104 s |
| mc | 41 (+16) | 39 + 2 xf | 0 | 109 s |
| cross_method | 81 (+8) | 81 | 0 | 13 s |
| sn | 3435 (+116) | **3384** + 1 sk + 50 xf | 0 | 15 min 12 s |
| derivations | 1661 (+67) | 1637 + 13 sk + 11 xf | 0 | 38 min 13 s |
| root + harness | 420 (+2) | 415 + 5 xf | 0 | 19 s |
| **total** | **10417** | **10329 passed · 19 skipped · 69 xfailed · 0 failed = 10417** | 0 | 61 min 34 s |

Every tree's measured count equals its predicted denominator unit-for-unit
(passed + skipped + xfailed = collected, 13 of 13). Against the last full
gate (CS4c step 4, `10106 / 0 / 19 sk / 66 xf`): **+223 passed, +3 xfailed**
— the campaign's tests, and the ERR-080 gate's three strict-xfail rows,
which are the rows the fused commit turns into XPASS-then-deleted. ⚠ The
wall time (61 min) is on a machine concurrently running two sub-agents and
the opener's probes; the `derivations` tree alone was 36 min at 0.2 on a
quiet machine. Predicate for the next re-gate: the same 13 trees, the same
filter; the fused commit's predicted deltas are +3 XPASS→pass in `sn/solve`
and whatever 2.5 and the fix add — write the prediction BEFORE the run.

⭐⭐ **RE-GATED at 2.5** (2026-09-02, the working tree that became the 2.5
commit; same runner, same 13 trees, `scratch/_fused_25_full_gate.log`):
predicted **10454** (transport 646 → 683, nothing else), measured **10366
passed · 19 sk · 69 xf · 0 failed = 10454, 13 of 13 rc=0**, every tree equal
to its predicted denominator; against the `2f294ef1` table above **+37
passed** (transport 682 + 1 sk) and every other row identical.

⭐⭐ **RE-GATED at the FUSED commit `5436184e`** (2026-09-02, launched on the
working tree that became the commit, finished 62 min after it;
`scratch/_fused_B_full_gate.log`): predicted **10556** (numerics 2809 ·
transport 707 · sn 3439 · root+harness 424, the other nine unchanged),
measured **10471 passed · 19 sk · 66 xf · 0 failed = 10556, 13 of 13 rc=0**,
every tree equal to its predicted denominator; against the 2.5 row **+102
collected = +105 passed − 3 xfailed** (numerics +70, transport +24, sn +4
with xf 50 → 47, root+harness +4), reconciled FILE BY FILE in the *FUSED
EXECUTED* section's gate table — read the numbers there, not here. The next
re-gate's delta is against THIS row; its own prediction is written before it
runs.

## Durable lessons — promoted OUT of this file

⛔ **Do not keep a COUNT here** (`plan-authoring` §9 — the file re-measures
itself, and a copied count is guaranteed to go stale while looking
authoritative). `.claude/rules/plan-authoring.md`'s surprise log is the record;
`[M]` it holds **72** dated rows at `14c37aa9`, and the campaign's own
contributions are the ones a `grep` for this plan's vocabulary finds. *(An
earlier version of this paragraph asserted "seven rows across the campaign" —
un-countable, since no mechanical predicate separates this campaign's rows from
the rest.)*

The most recent, 2026-09-02 (this compaction), and the two before it:

* **A consumer census that EXCLUDES the defining module cannot see
  in-module consumers** — and reads as "nothing names it". Founding failure:
  §V.5i(b)'s *"`Ball` has 0 production consumers"*, relayed into a docstring
  and a dispatch brief; `grep -v manifold.py` hid `Ball(2)`, the σ_y entry's
  own realization. (Promoted alongside it, to `vv-principles` #17: **a mutation
  control inside the SUT's own invariance group is a null control** — 2.3's
  μ → −μ control sat inside the product rule's σ_h and reddened 9 where two
  ordinary arms reddened 58/57.)

* **"The rebuild loses X" is a completeness claim over the source TYPE'S FIELD
  LIST** — and the denominator, `dataclasses.fields(T)`, is always available
  for free. Founding failure: two sections of this plan enumerated what a lossy
  `DiscreteMeasure` rebuild destroyed and both named 1 of 3, because they
  enumerated against the CONCEPT the campaign was about.
* **A §6b census keyed on the changed class's NAME cannot return its
  SUBCLASSES' constructors** — the most mechanical member of the
  "spelled without the symbol" family. Founding failure: 2.1's pre-flight sized
  the set at 18 having explicitly *considered and cleared* `OverlapBasis`; the
  set was 23.

⭐ **Five lessons from this campaign that are NOT yet promoted**, kept here
because it is not yet clear they generalise past it:

1. **A quantity attributed to a step survives that step being SPLIT.** The
   "92 % of phantom reads" was credited to 0.1; 0.1 then split into 0.1a/0.1c
   (now) and 0.1b (at 3.4) in a later section, and the 92 % silently attached
   to the half the reader was looking at. `[M]` 0.1a removes **0 %**. ⟹ when a
   step splits, re-read every number credited to its old name.
2. **A strictness change surfaces duplicates of whatever the laxness stood in
   for** — and those are invisible to a census of the lax call. `[M]` three
   sites hand-rolled the R³ embedding as `column_stack` of the accessor; a
   fourth was a duck-typed surrogate. 29 reds. ⟹ ask *"what was the laxness
   standing in for?"* and grep for THAT, plus the twins of every production
   site the step touches. `[M]` two greps, not one — one finds 3 of the 4.
3. **Refining a periodic-trapezoid reference makes it WORSE.** `[M]` the orbit
   rule is exact for trig polynomials at `n >= 2`; `|col0 - mu|` is 1.1e-16 at
   n=8 and 1.1e-12 at n=200 000. A gate's first draft failed its own tolerance
   against its own reference.
4. ⭐ **A `# type: ignore` on a line whose whole PURPOSE is to be refused
   DOCUMENTS the claim; it does not hide a defect.** `coding-elegance` #19 says
   no ignores, so at 2.1 I removed one from a `pytest.raises` leg that
   deliberately instantiates an abstract class — and pyright then correctly
   reported exactly the refusal the test asserts. `[M]` the tree settles it:
   **228** ignores ship and the dominant idiom is precisely this one
   (`TimedFullField(interior="not a field", …)  # type: ignore[arg-type]`).
   ⟹ #19 targets suppression, not documentation; the discriminator is whether
   the flagged line is the SUT or the assertion.
5. ⭐⭐ **Assigning a newly-minted type to its call sites is itself a CENSUS,
   and it finds things no count could.** 2.1's 17 test-site edits were meant to
   be transcription. `[M]` they separated two manifolds the string tag
   `"energy"` had conflated — a partition by energy VALUE in eV is an
   `Interval`, a partition by GROUP INDEX is `EnergyGroups`. Both have ambient
   dimension 1, so no structural or dimensional check could have found it;
   only being forced to NAME each site's point set does. ⟹ budget the
   site-by-site pass as investigation, not as chore.
7. ⭐ **A test that parses another test's OBJECT by NAME is a filter validated
   against one shape.** 3.1's symbolic-tie gate read the entry's SymPy
   generators assuming the sphere entries' `x y z` symbols; the generic trivial
   derivation spells `x0:n`, and the gate reddened on it at the first run.
   Found structurally (the slot is the one symbol with derivative exactly 1),
   the same gate then ran on all 7 entries. `plan-authoring` §2's FILTER
   clause, one tier down: a filter is validated against a known member of
   EVERY shape, and a naming convention is a shape.
6. ⭐ **A memo on a DERIVED object does not cover a guard that asks for its
   INPUTS — and growing a candidate set multiplies the guard, not the memo.**
   2.4 added three axial candidates to the lattice walk and two lattice-wide
   tests went from seconds to past the harness timeout. `[M]` the walk's
   time was 99 % one un-memoised generating-set builder (`_icosahedral_ops`,
   itself a 120-element closure) called 41 times from `_contains`'s
   finite-vs-finite guard — while `_group_elements`, one call downstream,
   had memoised the closure all along. The pre-2.4 set already paid 35
   rebuilds; nothing noticed because 8 s per walk hid inside a green suite
   (`lessons` L16). ⟹ when a change widens a set the lattice walks, profile
   ONE walk before and after; and memoise at the layer a guard touches, not
   only at the layer the result lives.

## Artifacts

⚠ **All in UNTRACKED `scratch/` — a `git clean` destroys them.** Load-bearing
content is carried in this file; the reproducers are not.
⭐ **`sigma_y_orbit_derivation.md` (758 lines) + `probe_sigma_y_orbit.py`** —
the Procesi–Schwarz derivation with verbatim SymPy, the Molien completeness and
`dim(𝔪/𝔪²)` minimality checks, and the chart-vs-section analysis the ruling was
made on. `[M]` the probe is re-runnable: **57 checks, rc=0**, 9 of 9 mutations
red. Its load-bearing content is carried into `manifolds.rst`, which is tracked;
the memo is not. · `p20c_*.py` (the 2.0c opener's censuses).
`_pl_slab_defect_repro.py` · `_pl_slab_fix_probe.py` · `_phantom_axis_probe.py`
· `_phantom_census_plugin.py` + `_driver.sh` + `.json` (**re-runnable — it is
0.2's done-when**) · `_fold_L_spy{,2}.py` · ⭐ **`_p01a_ground.py`** (the 12-rule route/tag/field table — re-runnable, and the fastest way to see what `angular_frame` carries) · `_p01a_mut/mutplugin.py` + `_p02_mut/mutplugin2.py` (the two per-arm mutation batteries, each with a CALL COUNTER — the thing that separated "no reds" from "never ran") · `_manifold_selfcheck.py`
(**re-runnable RST self-check: underlines, ladder, table widths, label
uniqueness, xref resolution with a positive control**) · `rcond_rederivation.md`
+ `probe_rcond_01…17_*.py` · `angular_symmetry_wiring_survey.md` ·
`basis_domain_verification_plan.md` · `n2n_pl_*.md` (4). ⭐ **NEW at 2.1/2.0c:** `_p21_mut/` (the 7-arm battery + its driver — ⚠ the driver's first version handed pytest ONE nonexistent path because **zsh does not word-split an unquoted scalar**; the fixed one uses an array and prints `rc` and `CALLS` per arm) · `_p21_redloop.log` (the 5-tree gate) · `_p21_numerics_final.log` · `_p21_pyright.log` (854 KB, whole-project) · `_p21_sphinx.log` · `_p20c_census.py` + `_p20c_census2.py` (**re-runnable — they ARE 2.0c's opener**, AST with positive controls). ⭐ **NEW at 2.4:** `_p24_compare_census.py` (the 29-site comparison census, AST, positive control) · `_p24_mutation.py` (**the 12-arm battery — re-runnable**, copy-aside integrity, prints its own `BATTERY COMPLETE`) · `_p24_gate_{numerics,transport,geometry,data,sn,cross_method,diffusion}.log` · `_p24_redloop{1,2}.log` (loop 1 is the §6b enumeration) · `_p24_commit.txt` · `_p24_pristine/` (the pre-2.4 production files, byte-exact). ⚠ `_p24_preflight.py` and `_p24_reground.py` predate the ruling and spell `SubgroupOfO3.SO2` bare — they no longer run as written; the reground script's landmark half still does. ⭐ **NEW at 2.3:** `_p23_bitid.py` (**re-runnable — the 60/5/15 bit-identity table**, loading the pristine copies from the session scratchpad; copy them into `scratch/` if they are wanted past the session) · `_p23_mut/mutplugin23.py` + `_p23_gate_driver.sh` (the eleven-arm battery; ⚠ its `$(...)` capture dropped four arms — `_p23_mut/_arm_*.log` are the re-runs to files) · `_p23_gate_{numerics,sn_prim_mesh_ang_frames,sn_operators,sn_sweep,geom_harness_mms,rootdocs,post}.log` · `_p23_red1.log` (the first red loop: exactly the 8 `new_space=` sites) · `_p23_executed_section.md` (the section's draft with placeholders) · `_p23_reground.py` (**re-runnable — this compaction's existence/hash/landmark/count pass**) · `_p23_issue_comment.md` · `_arch_23_sphinx.log`. ⭐ **NEW at 3.1:**
`_p31_census.py` (**re-runnable — 3.1's opener**: existence, fields,
constructor sites, `.reference` readers, importers, module-scope mints) ·
`_p31_witness.py` (π∘φ = pr₁, β∘π, on the real tree) · `_p31_inject.py` +
`_p31_inject.log` (**the §6d inject-and-run on a SHADOW copy** — 7/7 vs 0/7;
the shadow lives in the session scratchpad, rebuild it with `cp -R
orpheus <dir>/`) · `_p31_mut/mutplugin31.py` + `_p31_battery.sh` +
`_p31_battery{,2}.log` + `_p31_mut/_arm_*.log` (the nine-arm battery, per-arm
logs; `2` is the final tree) · `_p31_gate_driver.sh` + `_p31_gate_*.log` +
`_p31_post.sh` + `_p31_post_*.log` · `_p31_red1.log` · `_p31_sphinx{,_final}.log`
· `_p31_archivist_brief.md` (the brief the archivist refuted four times) ·
`_p31_executed_section.md` · `_p31_issue_comment.md` ·
`_p31_issue_derivations.md` (#431) · `_p31_pristine/` · `_arch31_*.py` +
`_arch31_shadow/` (the archivist's own probes and cycle reproduction) ·
`_arch_31_sphinx{,_baseline}.log`. ⭐ **NEW at this compaction:**
`_p31c_reground.py` (**re-runnable — the existence/landmark/caller/hash
pass**) · `_p31c_full_gate.sh` + `_p31c_full_gate.log` + `_p31c_gate_<tree>.log`
(**the 13-tree gate, prediction phase + run**). ⭐ **NEW at the FUSED commit:** `_fused_isotypic_probe.py` (**re-runnable — §III.8's probe, legs A–D**) · `_fused_g0_census.py` (the frame spy on real solves) + `_fused_baseline_slab_L01.npz` / `_fused_after_A_…` / `_fused_after_B_…` (the slab bit-identity baselines) · `_fused_iso_census.py` (**re-runnable — exit predicate 3**) · `_fused_6b_census.{md,py,out}` (the explorer's census) · `_fused_verification_plan.md` + `_fused_va_*.py` (the test-architect's memo and probes) · `_fused_edit_numerics{,_b}.py`, `_fused_edit_transport{,_b}.py`, `_fused_edit_bases_c.py` (the production edit scripts, anchored) · `_fused_B_mut/mutplugin_B.py` + `_fused_B_battery.sh` + `_fused_B_battery{,2}.log` (**the 19-arm battery — re-runnable**) · `_fused_B_red{1,3}.log` · `_fused_B_full_gate.{sh,log}` + `_fused_B_gate_<tree>.log` · `_fused_B_test_architect_brief.md`, `_fused_B_archivist_brief.md`, `_fused_B_executed_section.md`, `_fused_B_issue_comment.md`, `_fused_B_commit_msg.txt` · `_fused_B_pyright_tests.py` (⚠ its `ensure_import` is fooled by a function-local import — see the follow-up record) · `_fused_B_followup_tests.log` (the five-file re-run + pyright) · `_fused_B_collect_reconcile.log` (**re-runnable — the per-file collect diff against a `2ebdbc05` worktree**) · `_fused_B_followup_plan_edit.py` (**the follow-up's typing script**) · `_fused_B_edit_plan.py` · the HEAD-capture worktree `<scratchpad>/_head_tree` (`git worktree list` — remove with `git worktree remove` when done) + `_capture_tables.py` / `_tables_{HEAD,work}.npz`.

## Standing constraints

Host → `.venv/bin/python`; canonical `-O -m pytest -p no:randomly -m "not slow" -q`
SERIAL; branch + ff-merge, never squash; **never `git add -A`** (untracked
`scratch/` sits one flag away from the history); never `git stash`; never
`git checkout`/`restore` a tracked path; commit via `git commit -F -` with a
**quoted** heredoc; no source edits under a running gate, and no production edit
while a sub-agent holds the tree (L38); long gates = detached + a LOG, never a
2-minute foreground `pytest`; `sphinx -W` + `dead_references` at every phase
exit; a bare `assert` in `orpheus/` is inert under `-O` — use a real `raise`.

**This plan blocks** `.claude/plans/cs4c_binding_design.md` §18.6 step 1 (#426's
missing (n,2n) measurement), PAUSED. Tracking issue **#429**. The source plan is
in-repo at `.claude/plans/symmetry_quotient_plan.md`.
