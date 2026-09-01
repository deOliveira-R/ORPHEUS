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
partition — and that basis then calls its own coefficient space
`L2[coarse_cells_R1]`, naming a **spatial** manifold it has nothing to do with.
So R1 was not a hypothetical about a property I was about to add; the same lie
already ships one level over, in a name string. `FunctionSpace` carrying a
`Manifold` makes both names derived and the lie unspellable.


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
| 8 | `Basis.invariance_group: SubgroupOfO3` | absent | 2.1b |
| 9 | `measure.quotient_group: SubgroupOfO3 \| None` | `[M]` **absent** — `folded_by` is the discrete case on the wrong object, and Part IV records *"the slab's quotient group is recorded nowhere"* | 2.0d |
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


# Part VI — The action plan

Phases are ordered by dependency, not by size (D0.5). The `Q#` column maps each
item to the source plan's numbering so the two stay reconcilable. `[M]` the source plan is now IN-REPO at `.claude/plans/symmetry_quotient_plan.md` (399 lines) — read it for the Procesi–Schwarz derivation, the `Retract`-pair framing, and the confidence ledger this plan condenses.

## ⚡ Phase 0 — IMMEDIATE. Land regardless of any ruling below.

**Goal.** No object asserts a property its own data contradicts.
Every item here is true under every design choice in Phases 1–7.

| # | item | site | Q# |
|---|---|---|---|
| **0.1** | ⚡ **Stop fabricating the support, and delete the phantom read in the same expression.** `angular_frame` must not write `support=SPACE_SPHERE` over nodes with `‖Ω‖ ≠ 1`. §II.5 proves L1 and L4 are ONE fix. ⛔ **SPLIT at execution 2026-08-31 — see §II.9**: `0.1a` (3-D rules hand the frame their OWN measure; `[M]` bit-identical on 10 of 12 shipped rows) lands now; `0.1b` (the 1-D rows) **rides Phase 3.4**, because the lift `[-1,1] → S²` *does not exist as a map*; `0.1c` (the fold's `S^2/sigma_y` tag, §II.8's new leak) is gated on a `plan-authoring` §8 blast-radius measurement of `FunctionSpace` equality. | `directional.py:587-593` | — |
| **0.2** | ⚡ **A phantom component becomes unspellable.** *Proposed means as of 2026-08-31, SHARPENED at execution — see ⚠ below:* `axis_cosines(i)` raises for `i ≥ dim`. ⛔ a `raise`, **not** an `assert` — the canonical runner is `python -O`. ✅ **census DISCHARGED (§II.14) — and it REFUTED this means**: `[M]` a blanket raise breaks **3 legitimate production consumers** asking the FLOW question. ⟹ runs AFTER 0.1 (which deletes the one fabrication site, 92 % of all reads); 0.2 is then the accessor SPLIT. Done-when re-runnable: re-run the census, `directional.py:589` must be GONE. | `directional.py` | — |

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
| **2.0c** | ⭐ **`FunctionSpace` carries its `Manifold`; the two `L2[...]` NAME STRINGS become derived** (§III.10). | `grep` finds no `f"L2[...]"` literal; `[M]` `indicator_basis.py:284`'s live false name (`L2[coarse_cells_R1]` on an ENERGY basis) is gone by construction | — |
| **2.0d** | ⭐ **`measure.quotient_group`** — the field Part IV's own predicate needs and that `[M]` **exists nowhere**: `folded_by` records the DISCRETE fold on the `Quadrature` façade, and *"the slab's quotient group is recorded nowhere"* (Part IV obstacle 2). Derived by `quotient()`, `None` for an unquotiented measure. | `gauss_legendre(8).measure.quotient_group is SubgroupOfO3.SO2`, derived by construction, not declared | — |
| 2.1 | **`Basis.domain: Manifold`** (closes L2). ⛔ the 2026-08-31 `support` rename is REFUTED — §III.10: the collision is unreachable (`[M]` 0 `Basis` subclasses are operators) and `support` is mathematically false for a basis (`supp(f)` = where `f` is non-zero; for an indicator that is ONE cell). *Means:* a `domain` property beside the existing `space`, giving the basis the same two-level structure the measure has. ⚠ **`IndicatorBasis` takes it as a CONSTRUCTOR FIELD** (§II.15 R1: `[M]` ⛔ *"over 5 sites"* CORRECTED 2026-08-31 — by AST there are **18** `IndicatorBasis(...)` sites (**4** in `orpheus/`, 14 in `tests/`) over **three** manifold families; the *three* was right, the *five* was not. Any derived value hard-codes one of three — and `[M]` that lie already ships in its space NAME). | every shipped `Basis` answers; `SphericalHarmonicBasis.domain is SPACE_SPHERE`; the energy indicator no longer claims a spatial manifold | — |
| **2.1b** | ⭐ **`Basis.invariance_group`** — the OTHER side of Part IV's G0 predicate (`measure.quotient_group ⊆ basis.invariance_group`), which today has **neither** side. Derivable for every shipped basis: full SH → `Trivial`; `MirrorEven` → `Mirror(axis)`; the trivial isotypic sub-basis → `SO2`. | the four-row lattice table in Part IV runs as a test and reproduces its verdicts | — |
| 2.2 | **G0 at frame construction** (closes L3). ⭐ **§II.7: the PREDICATE already exists** as `AngularSymmetry.admits_domain` (`measure.support == self.support`) — `[M]` reachable only from `select_quadrature`, whose 21 callers are **all tests**; and `[M]` `solve_sn` is HANDED a quadrature (required positional, no default), so the selector was never on the path. ⛔ **REFINED by §II.10 — NOT "just wiring".** The predicate is correct over an **incomplete domain vocabulary**: `support` derives from `continuous_isotropy` alone, so `folded_product`'s `S^2/sigma_y` matches **no** geometry and stage 0 would **refuse the shipped cylinder**. 2.2 owes the ontology a second slot (the *discrete* quotient the rule took — `SubgroupOfO3` already names it) **and** the `(CoordSystem, ndim) → key` join, which `[M]` exists nowhere. Do not hand-roll a second predicate; do extend this one's vocabulary. | `GalerkinFrame(SphericalHarmonicBasis(L=2), slab_measure)` **RAISES**, with a test witnessing it on the exact pairing that ships today | — |
| 2.3 | ⭐ **The typed `Chart`** (closes L6, L7). *Proposed means:* a map carrying `domain → codomain`; `pushforward` derives `support` from the chart rather than taking `new_space`; `product()` composes its two 1-D factors through the Archimedes chart instead of `column_stack` + a literal. | `grep` finds no `SPACE_SPHERE` literal at a measure constructor; `product(4,8).measure.support` is **derived**; the factor structure is recoverable from the built rule | — |
| 2.4 | **The slab rule declares its quotient group** (closes L5). Includes settling Part IV's obstacle 1 — the `SO2` axis-convention collision — and obstacle 2, the `invariance_group` vs `folded_by` conflation. ⛔ **MOVED to Phase 0 on 2026-08-31** — the full-suite phantom census was scheduled here, but **0.2 depends on it** (it is 0.2's denominator, not 2.4's). Scheduling a step's own precondition three phases downstream is the `plan-authoring` §6b defect with a *census* in place of a call site. Run and published at Phase 0. | `SubgroupOfO3.SO2.is_invariant(gauss_legendre(8).measure)` returns the derived answer, and the plan records whether that answer is `True` **with the derivation**, not the prediction | — |

## Phase 3 — The symmetry core

**Goal.** A quotient space, its retraction/section pair, and its invariant
sub-basis are all *derived from a catalog entry*.
⚠ **Phase 3 does not ship alone** — it lands with Phase 4.1, its first consumer.

| # | item | depends on | done when | Q# |
|---|---|---|---|---|
| 3.1 | `numerics/symmetry/` catalog: `IsometryGroup` entries from 1.1. ⭐ **D0.1's SEED requirement is binding here**: an entry's FIELDS are §III.1's OUTPUTS — invariant generators `p₁…p_k`, syzygy ideal `I`, `P_ij = ⟨∇pᵢ,∇pⱼ⟩` and `det P`, chart, pushforward measure, and the stratum where `det P` vanishes — **not** a human summary of the answer, so a future engine populates them without minting a single new type. Plus the action on the axis index set, Haar (or counting) measure, a `derived_by: "hand" | "engine"` provenance field, and a symbolic test reproducing its own derivation (⭐ those ~12 tests ARE the engine's acceptance suite, written before it). ⚠ `[M]` `numerics/symmetry.py` is a **module** today — this is a re-home; `plan-authoring` §6d (check the import edge BOTH ways, and look for the import-linter's declared contract first). | 1.1 | every entry's chart and measure reproduced by its own regression test | Q1.1 |
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
| 0.1a | 3-D rules hand the frame their OWN measure — ⚠ keystone is a **ROUTE** gate (`frame.measure is q.measure`), not bit-identity (§II.15 R4) | ⏏ yes | ☐ |
| 0.1b | the 1-D rows — ⛔ **rides 3.4**, the lift does not exist (§II.9) | ⏏ yes | ☐ (with 3.4) |
| 0.1c | the fold's `S^2/sigma_y` tag — ⛔ §8 hazard REFUTED (0 reds; 0/145 089 checks); was blocked on **2.1-W** | ⏏ yes | ☐ ✅ **UNBLOCKED** — 2.1-W landed, so a carve on the fold path is now falsifiable at the eigenvalue tier |
| 0.2 | phantom component unspellable — ⛔ means REFUTED by the census; **runs AFTER 0.1**, then it is the accessor SPLIT (direction vs flow) | ⏏ yes | ☐ |
| 0.2-census | full-suite phantom census (moved here from 2.4) | ⏏ yes | ✅ **DISCHARGED** — §II.14, 13 trees rc=0, 11 sites / 1604 reads |
| 0.3 | `metric.py` "noise mode" + rcond re-derivation | ⏏ yes | ✅ `ae4dbc1f` |
| 0.4 | `directional.py` `spherical_harmonics` docstring | ⏏ yes | ✅ `ae4dbc1f` |
| 0.7b | ⭐ NEW — gate the SECOND symptom (the `DenseMetric` RAISE) | ⏏ yes | ✅ `ae4dbc1f` |
| 0.5 | `spherical_harmonics.rst` rank-deficiency correction | ⏏ yes | ✅ `ae4dbc1f` |
| 0.6 | `_evaluate_real_sh` raises on non-unit directions | ⏏ yes | ☐ (rides 3.4 — XII.1b) |
| 0.7 | strict-xfail reproducer gate + **ERR-080** minted | ⏏ yes | ✅ `ae4dbc1f` |
| 1.1 | catalog derivations — **`SO(2)/S²` entry is on the path** | ⏏ partial | ☐ |
| 1.2 | Gelfand-pair reading of Funk–Hecke | no | ☐ |
| 1.3 | connection-term verify-or-kill | no | ☐ |
| 1.4 | record the four gates | no | ☐ |
| 1.5 | record the placement ruling | no | ☐ |
| 1.6 | `ℛ* = ℛ` as adjoint-vs-dual | no | ☐ |
| 1.7 | the support algebra, written | no | ☐ |
| 1.8 | vocabulary reconciliation | no | ☐ |
| 2.0a | ⭐ **MINT `Manifold`** (D0.7, user-ruled) — retires `Space = str` | ⏏ yes | ✅ `de29bcc6` — 9 variants, 40 tests, pyright 0. ⚠ the TYPE only; `support` is still `str` (that is 2.0c/2.1) |
| 2.0b | `Manifold.contains` — the membership predicate; refuses the forged measure at construction | ⏏ yes | ◐ **HALF** — the predicate ships and is gated both legs (`[M]` refuses the `gauss_legendre(8)` forgery, norms `[0.1834, 0.9603]`; admits the normalised control) @ `de29bcc6`. The **refusal AT CONSTRUCTION** is unbuilt — that is the wiring, and it rides 0.1. ⛔ no `catches("ERR-080")` marker until then: refusing a forged ARRAY is not the production path refusing it |
| 2.0c | `FunctionSpace` carries its `Manifold`; the two `L2[...]` name strings become derived | ⏏ yes | ☐ |
| 2.0d | `measure.quotient_group` — `[M]` the slab's is recorded NOWHERE today | ⏏ yes | ☐ |
| 2.1 | `Basis.domain: Manifold` (⛔ `support` rename REFUTED — §III.10) — `IndicatorBasis` takes a ctor field | ⏏ yes | ☐ |
| 2.1b | `Basis.invariance_group` — G0's other side; today NEITHER side exists | ⏏ yes | ☐ |
| 2.1-W | ⭐ **PRE-CARVE**: the fold-basis witness | ⏏ yes | ✅ `TestFoldedCylinderP1BindsTheQuotientBasis` — 3 rows, L1, `verifies("pn-scatter", "discrete-measure-quotient")`. Keystone is k_inf-invariance on a HOMOGENEOUS folded cylinder, not the 2g-het solve the row proposed (§II.16). `[M]` rebind reds all 3 at nine orders of separation. **0.1c is unblocked** |
| 2.2 | G0 at frame construction — ⚠ ALSO owns the `AngularSymmetry` Γ-slot (§II.10, inventory #11): `support` derives from `continuous_isotropy` alone, so stage 0 refuses the shipped cylinder | ⏏ yes | ☐ |
| 2.3 | the typed `Chart`; pushforward derives support | ⏏ yes | ☐ |
| 2.4 | slab declares its quotient group; SO2 axis collision; **full-suite phantom census** | ⏏ yes | ☐ |
| 3.1 | `numerics/symmetry/` catalog (⏏ scoped to `SO(2)`) | ⏏ partial | ☐ |
| 3.2 | `ReynoldsProjection` | no | ☐ |
| 3.3 | `OrbitAxis` + retraction/section minting | no | ☐ |
| 3.4 | ⭐ **trivial isotypic sub-basis — THE FIX** | ⏏ yes | ☐ |
| 3.4-R | ⚠ the encoding ruling | ⏏ yes | ✅ **RULED 2026-08-31 (user): (b) true `LegendreBasis`** — see below |
| 3.4-R2 | ⚠ the DESCENT ruling | no | ✅ **RULED 2026-08-31 (user-endorsed): (i) mint `Descent`** — Part V.5 §C |
| 3.4b | ⭐ **`Descent`** — the iso witness; lands WITH or BEFORE 3.4 | no | ☐ |
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

# Part XIV — ⏸ COMPACTION POINT 1 (2026-08-31, evening) + ▶ Resume surface

**Read in this order on pick-up:** this Part → Part XIII (the tracker) → Part XII
(the exit gate) → the phase you are executing. **Part IX, §II.15 and §II.16
first if you are about to re-propose anything** — twelve premises are refuted
there, most of them mine, several the same day they were written.

⚠ This supersedes COMPACTION POINT 0. Its content is folded in below; nothing
from it was dropped.

---

## ▶ RESUMES AT — stated as OUTCOMES (`plan-authoring` §1)

**A function space knows the manifold its functions are defined on, and a
measure knows the group it was quotiented by.** Trackers **2.0c** and **2.0d**.
Today both facts are smuggled: the manifold through a NAME STRING at 2 sites
(one of them `[M]` already FALSE — §III.10), and the quotient group **nowhere at
all**. They are the two operands the G0 predicate has always named and never
had.

*Then* **2.1** (`Basis.domain: Manifold`) and **2.1b**
(`Basis.invariance_group`), which complete G0's other side.

### The migration this opens, measured — read §V.5c before decomposing it

`[M]` `support = str → Manifold` does **not** touch 87 sites. Split by the AST
node of the value passed:

| the value passed | `orpheus/` | `tests/` | what the retype costs |
|---|---:|---:|---|
| `SPACE_*` constant | 13 | 3 | **nothing** — the constant is retyped once |
| attribute / other Name (forwarded) | 9 | 4 | **nothing** — it forwards |
| raw string literal | 3 | **51** | one edit each; 94 % of it in `tests/` |
| **f-string** — a tag BUILT by interpolation | **4** | 0 | ⭐ **the real work** |

⟹ the production tree is nearly free; the bulk is 51 test literals; and the
**4 f-strings are the architecture**, because each is a parametric family the
string type hid (`[{a},{b}]`, `spatial_R{d}`, `index({label})`, and
`sn_trace_orbit{orbit}_g{group}` — the last a SECOND spelling of the
index-set family, which the 2.0a row's member list never named).

`[M]` and **18** `.support` reads do string manipulation, **all 18 in
`orpheus/`, 0 in `tests/`** — of which three ARE the verbs (`measure.py:588`
`__mul__`, `:1022` `quotient`, `:802` the pushforward codomain) and one is the
`L2[...]` name at `:331` that 2.0c retires.

### ⛔ Two landmines for the pick-up session's own greps

1. **`invariance_group` returns 55 hits and the plan's claim still holds.**
   `[M]` by AST exactly ONE class defines it — `DiscreteMeasure` — and it means
   *the group this measure IS invariant under*, **not** the group it was
   quotiented BY. `[M]` `hasattr(SphericalHarmonicBasis(L=1), "invariance_group")`
   is `False`; 0 of 6 `Basis` subclasses define it. So both sides of G0
   (`measure.quotient_group ⊆ basis.invariance_group`) are genuinely absent, and
   a bare grep says otherwise. §3's ambiguous-name hazard, live.
2. **`measure.py:411`'s `self.support == "cells"` has NO producer.** Denominator
   enumerated in §V.5c(d): 58 literal `support=`/`new_space=` sites (13 distinct
   values), 6 `SPACE_*` constants, 36 non-literal producers — `"cells"` in none.
   The only `"cells"` in the tree is `radial_characteristic_space.py`'s
   `part="cells"`, a walk-part label. Three surfaces advertise the dead arm (the
   branch, the docstring at `:387`, the `raise` at `:416`). **2.0a acceptance
   item**: the mint must make it unspellable and the prose must stop naming it.

### §1 existence-checks on every symbol this pointer names — run 2026-08-31 at `fba4205a`

`[M]` in `orpheus/`: `class Chart` **0** · `FunctionSpace.manifold` **0** ·
`measure.quotient_group` **0** · `Basis.domain` **0 of 6** subclasses ·
`class LegendreBasis` **0** · `class Descent` **0** · `class ReynoldsProjection`
**0** · `class OrbitAxis` **0**. ⟹ nothing the pointer promises already exists.

⭐ **And the DELIVERABLE-level check, which a symbol grep cannot answer** (§1's
clause): `AngularSymmetry.support` (`quadrature/registry.py:869`) **already
performs the orbit-space lookup 2.0a's `Manifold.quotient` mints** — same rows,
same refusal shape, in the string vocabulary. It is a Pattern-2 twin, not a
duplicate deliverable: tracker **2.2** absorbs it, and a committed row already
pins the two in agreement.

---

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
| 2.1: *"three manifolds over 5 sites"* | ⛔ **18 sites** (4 `orpheus/`, 14 `tests/`); the *three* was right. The 2.1 row |

## Measured baselines and costs

* **numerics tree, as a GATE** (`-O -m "not slow"` serial, `de29bcc6`):
  `[M]` **2578 passed / 0 failed**, 337 s. Reconciles exactly: census `2538` + 40.
* **Verification matrix**: 10433 → **10481**. Every unit attributed — +40 (2.0a
  foundation) +3 (2.1-W, L1) +1 (`test_layer_imports.py` builds its parametrize
  list by `rglob`, so a new module adds a row — and it PASSES, so the §6d edge
  I checked by hand is now enforced automatically) +4 (the trivial-quotient and
  twin-agreement rows).
* **Sphinx `-E -W --keep-going`**: **0 warnings / 0 errors**, ~2 min.
  **`dead_references`**: **0 dead / 52 checked**. Theory labels 905 → 907.
* **Phantom census** (a recording plugin, rc ignored — NOT a gate): 13 trees,
  65 min, numerics 2538 · transport 645 · sn 3375 · geometry 727 · mc 39 ·
  derivations 1637 · rootfiles 414.
* ⚠ **The 13-tree GATE baseline is still UNMEASURED at this HEAD.** Exit
  predicate 5 stands. Last full gate: CS4c step 4's `10106 / 0 / 19 sk / 66 xf`.

## Durable lessons — promoted OUT of this file

`plan-authoring` gained five rows across the campaign; they are the record, not
this plan: the transcription-precision row; the §8 sharpening (*a branch's
EXISTENCE is not its RELEVANCE*); the §2 VIEWPORT re-commit; **the LIST-of-counts
row** (a row pricing work with several counts implies a shared denominator that
was never stated); and **the §6d relative-import row** (the AST import census
§6d prescribes drops `from .x import y`, and fails toward *"no cycle, proceed"*).

## Artifacts

⚠ **All in UNTRACKED `scratch/` — a `git clean` destroys them.** Load-bearing
content is carried in this file; the reproducers are not.
`_pl_slab_defect_repro.py` · `_pl_slab_fix_probe.py` · `_phantom_axis_probe.py`
· `_phantom_census_plugin.py` + `_driver.sh` + `.json` (**re-runnable — it is
0.2's done-when**) · `_fold_L_spy{,2}.py` · `_manifold_selfcheck.py`
(**re-runnable RST self-check: underlines, ladder, table widths, label
uniqueness, xref resolution with a positive control**) · `rcond_rederivation.md`
+ `probe_rcond_01…17_*.py` · `angular_symmetry_wiring_survey.md` ·
`basis_domain_verification_plan.md` · `n2n_pl_*.md` (4).

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
