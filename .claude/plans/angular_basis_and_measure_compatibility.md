# P_L scattering is correct on every chart, and an incompatible basis–measure pairing is unspellable

**Status:** OPEN. Diagnosis COMPLETE and converged from three independent
perspectives; no production line changed yet.
**Opened:** 2026-08-31, during CS-ladder Campaign 2 (§18.6 step 1 — #426's
missing measurement), which is PAUSED behind this.
**Authority for the pause:** the user, on being shown the measurement below.

> ⚠ **Read `plan-authoring.md` §1–§3 before editing this file.** Every claim
> below carries an epistemic marker. `[M]` = measured, with its configuration.
> `[R]` = reasoned, unverified. ⛔ = refuted, kept in place per §3.
> ⚡ = LOW-HANGING FRUIT — land regardless of any later ruling.

---

## 0. GOAL, separately from any proposed means

**Goal.** An anisotropic (P_L) scattering source is correct on every chart
ORPHEUS ships, at every Legendre order the solver admits — and a basis bound to
a measure it is not compatible with cannot be constructed at all.

**Non-goal.** This plan does not chase performance, does not re-pose the
scattering operator, and does not touch the energy or spatial axes except where
§6 measures a shared consumer.

**Done when** (checkable predicates, in order of strength):
1. `scratch/_pl_slab_defect_repro.py` LEG A reads `+4.000000000000` at
   `L = 0, 1, 2, 3` — today it reads `-3.764705882353` at `L >= 2`.
2. `[M]` a fresh census of the isotropic-flux `l >= 1` moment over every shipped
   `(rule, L)` reads `0` to machine precision — today 4 of 24 rows do not.
3. Constructing a frame from a basis whose domain does not contain the measure's
   support RAISES, and a test witnesses the raise on the exact pairing that
   ships today.
4. No production site writes a `support` tag that its own nodes do not satisfy.

---

## 1. ⛔ THE DEFECT — measured

**Reproducer:** `scratch/_pl_slab_defect_repro.py` (two legs, both self-checking).

### 1.1 End-to-end, 1-group infinite medium

Configuration: `Sigma_t = 1.0`, `Sigma_s0 = 0.5`, `Sigma_s1 = 0.40`,
`Sigma_s2 = 0.45`; `Quadrature.gauss_legendre(n_ordinates=8)`; 4 cells,
all-reflective; uniform isotropic unit source; `inner_solver="krylov"`,
`inner_tol=1e-13`.  In an infinite medium the flux is isotropic, so every
`phi_l` with `l >= 1` is ZERO and `scattering_order` **cannot** matter.
Analytic answer `phi = W*Q/(Sigma_t - Sigma_s0) = 4.000000000000`.

| L | phi | rel err |
|---|---|---|
| 0 | +4.000000000000 | 1.2e-15 |
| 1 | +4.000000000000 | 4.4e-16 |
| **2** | **-3.764705882353** | **1.94** |
| **3** | **-3.764705882353** | **1.94** |

Wrong sign. 194 %.

### 1.2 On real data

`[M]` 421-group pure-Be slab (`BE009`, N = 0.1236 at/b-cm), 20 cm, 12 cells,
vacuum both faces, uniform U-235-chi source, S8 GL, Krylov:
`sum(phi)` = **2.149e3 / 1.978e3 / 2.267e4** at L = 0 / 1 / 2 — the L=2 answer
is **10.6x** the L=0 one.

### 1.3 Reachability and coverage

* **Reachable and unguarded** from the public API: `solve_sn(...,
  scattering_order=2)` and `solve_sn_fixed_source(...)`.
* `[M]` **Nothing validates `L` against the quadrature.** The sole clamp is
  `orpheus/sn/solver.py:1357-1363`, against `min(len(m.SigS) - 1)` — the DATA's
  Legendre depth, not the rule's. The `exactness` / `degree_of_exactness`
  machinery from the Q6 campaign is structurally disconnected: 0 occurrences in
  `frame.py`, `harmonic_frame.py`, `scattering.py`, `n2n.py`, `fission.py`.
* `[M]` **No witness anywhere.** The corpus uses `scattering_order` =
  0 (74 sites), 1 (46), 3 (2) — and BOTH `L=3` sites
  (`tests/sn/operators/test_scattering_operator.py:1593` and
  `tests/sn/_fixtures/wave_t_t3/_capture_pre_t3_snapshots.py:134`) use
  `Quadrature.lebedev(17)` on a 2-D mesh, where the defect does not bite.
  **P >= 2 on a 1-D chart is exercised by nothing.**

### 1.4 The scope, measured on PRODUCTION's own bound basis

`q.spherical_harmonics(L)` (i.e. whatever basis the rule actually binds), the
`l >= 1` moment of an ISOTROPIC flux, which must be 0:

| rule x L | live | rank | max offdiag | raw `l>=1` of isotropic | verdict |
|---|---|---|---|---|---|
| `gauss_legendre(8)` L=1 | 2 | 2 | 8.8e-17 | 6.9e-18 | clean |
| **`gauss_legendre(8)` L=2** | 5 | **4** | 1.155 | **5.774e-1** | RANK-DEF SPURIOUS DENSE |
| **`gauss_legendre(8)` L=3** | 9 | **7** | 1.155 | **5.774e-1** | RANK-DEF SPURIOUS DENSE |
| **`gauss_legendre(16)` L=2,3** | 5, 9 | **4, 7** | 1.155 | **5.774e-1** | same |
| `level_symmetric(4)` L=3 | 16 | 16 | 1.9e-1 | 8.8e-17 | DENSE only |
| `folded_product(2,4)` L=3 | 7 | **4** | 2.793 | 1.1e-16 | RANK-DEF only |
| the other 18 rows | — | — | — | <= 1.2e-16 | clean |

**7 of 24 rows are not clean; 4 of them are the wrong-answer defect.**  The
spurious value is exactly `1/sqrt(3) = 0.5773503` at every slab order tested
(`n_ordinates` = 2, 4, 6, 8, 16) — a structural constant, not a quadrature-
accuracy artifact.

---

## 2. ROOT CAUSE — read from the source, not inferred

`_evaluate_real_sh`, `orpheus/numerics/basis/spherical_harmonic_basis.py`.
It hard-codes `l = 1` in CARTESIAN form:

```python
Y[:, 1, 0] = mu_z ;  Y[:, 1, 1] = mu_x ;  Y[:, 1, 2] = mu_y
```

then for `l >= 2` switches to the SPHERICAL route:

```python
cos_theta = mu_x
sin_theta = sqrt(1 - mu_x**2)      # [M] 0.279 .. 0.983 on GL8 -- NOT small
on_axis   = sin_theta < 1e-15      # [M] fires 0 of 8 -- never triggers
cos_phi   = mu_y / sin_theta       # 0 / 0.279 = 0
sin_phi   = mu_z / sin_theta       # 0
phi       = arctan2(0.0, 0.0)      # = 0.0  <-- FABRICATED
```

so `cos(m*phi) = 1` and **every `m > 0` slot is non-zero**.

⭐ **`l = 1` is clean ONLY because it never takes this path** — which is exactly
the measured L boundary in §1.1. That coincidence is why the defect looks like
an "L >= 2 problem" rather than what it is.

⭐ **The `on_axis` guard proves the case was anticipated and mis-aimed.** It
fires when `sin_theta ~ 0`, i.e. `|mu| ~ 1` — a direction ALONG the polar axis.
The slab's situation is the opposite: the TRANSVERSE part is zero while
`sin_theta` is large. `[M]` over 10 shipped rules the guard fires on **2 nodes
of 1 rule** (the Lebedev poles) while the case it misses is **100 % of the
ordinates on every slab rule**.
⚠ And it is written with `np.where`, not `if` — so **no coverage tool can ever
report the branch unexercised.**

### 2.1 The convention mismatch, stated exactly

> `Quadrature` supplies `mu_y = mu_z = 0` meaning **"there is no azimuthal
> information"**. `_evaluate_real_sh` reads it as **"the azimuth is 0"**.
> Two meanings, one set of zeros.

That is AI-failure-mode #6 (convention drift), at the quadrature <-> basis
boundary, and `Quadrature.spherical_harmonics`'s own docstring
(`directional.py:531-540`) records the FIRST meaning — *"only the m = 0
harmonics carry non-zero values; the other slots are filled with zeros"* —
while the code implements the second. **The docstring states the invariant that
would make the code correct.**

### 2.2 The deeper statement

`[M]` `||Omega||` over the shipped rules:

| rule | `measure.support` | `dim` | min..max `\|\|Omega\|\|` |
|---|---|---|---|
| `gauss_legendre(2)` | `[-1,1]` | 1 | 0.5773502692 (all) |
| `gauss_legendre(8)` | `[-1,1]` | 1 | 0.1834 .. 0.9603 |
| `gauss_legendre(32)` | `[-1,1]` | 1 | 0.0483 .. 0.9973 |
| every 3-D rule | `S^2` / `S^2/sigma_y` | 3 | **1.0 exactly** |

**The slab nodes are not on the sphere the basis is defined over.** A slab
ordinate `(mu, 0, 0)` is not a direction; it is a 1-D scalar node embedded as a
3-tuple, which is only syntactically a point of `S^2`.

⟹ `[R]`, and the ontological key: **the honest image of a slab ordinate in
`S^2` is the whole latitude CIRCLE** `{(mu, sqrt(1-mu^2) cos p, sqrt(1-mu^2) sin p)}`
— an ORBIT of the axial rotation group, not a point. That is precisely why the
`phi'`-integration annihilates `m != 0` (§4.2), and why `support = "[-1,1]"` is
right: **a 1-D angular quadrature is a quadrature on the quotient `S^2/SO(2)`.**

---

## 3. THE LEAKS — non-intrinsic properties and what they signal

The user's discipline: *items carry their intrinsic properties; composability is
checked exactly where the two objects meet; properties at the meeting place are
intrinsic or DERIVED from them; a hard-code is a signal that the individual
objects are missing properties.*  Applying it:

| # | leak | evidence | what it signals |
|---|---|---|---|
| **L1** | `Quadrature.angular_frame` **fabricates** `support=SPACE_SPHERE` onto interval nodes | `directional.py:587-593`; `[M]` rule support `[-1,1]` -> frame support `S^2`, `\|\|node\|\|` 0.1834..0.9603 | `Basis` cannot declare its DOMAIN, so the frame forges a matching support to make an incompatible pairing typecheck |
| **L2** | `Basis` has **no `domain`** | `[M]` `basis/base.py` declares `evaluate / synthesize / analyze / analyze_transpose / reconstruct / reconstruct_transpose / gram_structure / mass_matrix / space` — no domain. Its own docstring: *"the nodal/domain space comes from the [measure]"* | the missing intrinsic property. Nothing at the meeting place can check anything |
| **L3** | no composability check at the meeting place | `[M]` `GalerkinFrame(basis, measure)` validates neither support nor domain | the check has no home even if L2 were fixed |
| **L4** | a 1-D rule exposes a 3-D direction tuple with fabricated zeros | `[M]` `mu_y`/`mu_z`/`axis_cosines` = **171 occurrences over 20 files** (denominator 340 `.py`) | the sweep genuinely needs `mu_x`; `mu_y = mu_z = 0` is not a property of a 1-D rule. NOT effort-excused — see §7 P4 |
| **L5** | `_harmonic_basis` selects the basis from `folded_by` alone | `directional.py:505-530`; `[M]` a slab GL rule has `folded_by is None` | a NATIVE 1-D rule was never "folded" from anything, so `folded_by` is the wrong slot to ask. The question is the measure's SYMMETRY, not its construction history |
| **L6** | `Space = str`, *"recommendations, not constraints"* | `measure.py:110-119` | even a declared domain cannot be enforced today |

⭐ **The codebase already learned L1's lesson one level down.** `measure.py:120-128`
records that `SPACE_CIRCLE` was `"[0,2pi)"` until 2026-08-02, which
*"made the tag assert a coordinate"*, and was renamed to name the **manifold**.
Writing `S^2` onto interval nodes is the same defect one level up: **the tag
asserting a membership the nodes do not have.**

---

## 4. WHAT THE LITERATURE SETTLES

Full memo: `scratch/n2n_pl_consistency_literature.md` (555 lines, page-cited).
⚠ `scratch/` is UNTRACKED — the load-bearing content is carried here so a
`git clean` cannot destroy the finding.

### 4.1 The problem has a name

**Galerkin quadrature — Morel, J. E. (1989), *NSE* 101(1), 72–87**, local at
`scratch/literature/Morel(1989)...pdf`. The moment-to-discrete / discrete-to-
moment matrices `M`, `D`.
⚠ Morel's primitive is `D` and he writes `M = D^{-1}` (Eq. 27); the modern
convention makes `M` primitive and writes `D = M^{-1}`. **Same statement** — but
quoted side by side they read as a contradiction.

### 4.2 The 1-D answer is a DERIVATION, not a convention

**Bell & Glasstone (1970), §2.6a, printed p. 103** — in the addition theorem the
`cos m(phi - phi')` terms *"vanish upon integration over `phi'`"*, leaving
Eq. (2.79): a slab transport equation whose scattering source is a **pure
Legendre sum**.
**§5.3a, printed p. 227** — spherical harmonics appear in general geometry
*"because there is no azimuthal symmetry to eliminate terms from the addition
theorem of Legendre polynomials"*.

⟹ In slab geometry the `m != 0` harmonics are not small, not degenerate:
**they are identically absent from the continuum equation**, removed before any
discretisation. Morel's 1-D `D`/`M` (Eqs. 39a/39b) carry **no `m` index at
all**; the `(l,k)` double index first appears for the CYLINDER (Eq. 45), with an
`l + k` EVEN restriction.

### 4.3 The crux: invertibility is the PREMISE, not the output

**Larsen & Morel (2010), Eqs. (1.104)/(1.106), printed p. 40-41** declare `D`
and `M` to be `N x N` **by definition**, justified by unisolvence: *"One can
uniquely define a polynomial of degree `N-1` either in terms of `N` Legendre
moments or in terms of `N` discrete function values at `N` distinct points.
Therefore `D` and `M` should be inverses of one another."*

⟹ **No pseudo-inverse can repair a rank-deficient `M`** — see §8.1, where this
was measured before the literature was consulted, and the two agree.

### 4.4 What the literature does instead

Morel 1989 §VI: 2-D `S_2` has 4 directions but a `P_1` source has 3 harmonics,
so `M` is `4 x 3`. The remedy is **not** a pseudo-inverse, **not** least-squares,
**not** regularisation: *"one must find an additional higher order spherical
harmonic function that will provide a moment-to-discrete matrix that can be
inverted."*  Corroborated by Adams & Larsen (2002) §7.

⟹ **The universal move is: CHANGE THE MOMENT SET UNTIL `M` IS SQUARE AND
INVERTIBLE.** Add harmonics when there are too few; ORPHEUS's slab has too many
AND a dependent set, so the mirror move is to remove.

### 4.5 ⚠ An honest negative

**UNFOUND**: no source in the corpus states that 3-D real spherical harmonics
become linearly dependent on a 1-D ordinate set. The literature does not pose
the question **because no established formulation ever puts them there.** The
absence is itself the finding — this is not a known-and-solved problem, it is a
configuration the canonical formulations exclude by construction.

---

## 5. THE STRUCTURAL CONDITION

Full memo: `scratch/n2n_pl_frames_attack.md` (858 lines, 7 frames + an
UNEXPLORED section). Load-bearing results carried here.

### 5.1 The exact well-posedness condition

> **`reconstruct o multiply o analyze` is well-defined iff `ker(synthesis)` is
> invariant under the multiplier** — a property of the (basis, measure,
> operator) **TRIPLE**, not of any pair.

If synthesis has a non-trivial kernel and `Lambda` does not preserve it, then two
moment vectors representing the SAME angular function produce DIFFERENT
sources: **the scattering operator is not a function of the angular flux at
all, only of the coefficient representative.** That is a well-posedness
failure, not an accuracy failure.

### 5.2 The proposed diagnostic ladder — four tiers, nested

| tier | question | object | cost |
|---|---|---|---|
| **0 categorical** | is `measure.support` inside `basis.domain`? | a `Basis.domain` property + a frame-construction check | free |
| **1 symmetry** | what is the measure's quotient group, and what is its invariant sub-basis? | `InvariantSphericalHarmonicBasis(group)`, mask DERIVED by probe | one probe at basis construction |
| **2 fidelity** | does the SHIPPED composite reproduce the multiplier? | a `Lambda`-free, graded frame verdict | one `O(K)` round trip, cached |
| **3 explanation** | on failure: is it repairable, and which modes are guilty? | `ker` invariance + the dependent circuits | one SVD, failure path only |

`[D, from the frames memo]` tier 2 implies tier 3, so the cheap gate is a
sufficient condition for the expensive one and tier 3 runs only on failure.

⚠ **Before arming tier 2 as a refusal, run its own predicate at its stated
denominator** (`plan-authoring` §10 third shape) over every shipped `(rule, L)`.
`[R]` it is expected to red on over-truncation too (`level_symmetric(4)` L=3,
`folded_product(2,4)` L=3) — a feature, but it must be MEASURED before it is
armed, not discovered as a red suite.

### 5.3 ⭐ The repair is the existing quotient concept, one symmetry up

`MirrorEvenSphericalHarmonicBasis` is documented as **"the QUOTIENT's basis"**,
and `_harmonic_basis` already carries the seam, verbatim:

> *"Only a single-mirror fold has the +/- parity split this realizes; folding by
> any other group **refuses here until a consumer exists**."*

**A consumer now exists: the slab, whose measure is a quotient by `SO(2)`.**
`SubgroupOfO3` already ships `SO2` and `Dinfh`.

`[M]` probing invariance under SO(2)-about-x, using the SAME idiom `MirrorEven`
uses (evaluate at generic directions and at their group images, classify per
slot), over **6 INCOMMENSURATE angles x 9 generic random directions** at `L = 4`:
the invariant set is **exactly** `{(0,0), (1,0), (2,0), (3,0), (4,0)}` —
**5 of 25 real slots, zero false positives, zero false negatives.**
⚠ The angles are incommensurate deliberately: a sample of right angles would
generate `C_4`, not `SO(2)` (vv-principles #13's generator-set trap, which
`MirrorEven` avoids only because a mirror has order 2).

⟹ **A zonal basis is not a hand-rolled special case.** It and `MirrorEven` are
two instances of ONE concept — *the sub-basis fixed by the group the measure
quotients by* — i.e. the trivial-representation (isotypic) subspace. And the
rule for WHEN to bind it is DERIVED from the measure's symmetry, not declared by
a chart test.

---

## 6. BLAST RADIUS — measured, and far smaller than feared

Full memo: `scratch/n2n_pl_blast_radius.md` (869 lines; Nexus-derived, every
count with its predicate and denominator).

* `[M]` **4** `FrameBase` subclasses (denominator 457 `ClassDef`s / 340 files);
  **11** production constructor call sites in 7 files. The whole angular story
  is **2 of them**, chained `Quadrature.angular_frame(L)` ->
  `HarmonicFrame.from_galerkin`, both interned per `(rule, L)`.
* `[M]` **4** production `HarmonicFrame.for_space` sites: `scattering.py:811`
  (`L = scattering_order`), `scattering.py:960`, and `fission.py:305` /
  `n2n.py:200` — the last two **hardcode `L = 0`**.
* `[M]` **One frame bypass**: `dsa.py:585` reads `angular_frame(1).table` raw.
* `[M]` **The only l-stacked per-mode apply is `material_field.py:244/270/284`**
  (`_moment_blocks` / `moment_source` / `moment_source_transpose`), cross-checked
  by an independent `conjugate`/`reconstruct_after` census (6 composition sites).
* ⭐ `[M]` **HARD NEGATIVE**: CP, MoC, MC, diffusion, homogeneous, kinetics,
  fuel, TH — **24 `.py` files, 0 angular-moment lines**. `derivations/` is out
  too (spatial Legendre; its `scattering_order` is case metadata, max 1).
* `[M]` **DSA is immune because it binds `L = 0`**, where `Lambda = sigma*I` and
  the kernel is trivially invariant. **That is why nobody noticed.**

⟹ Today this defect has **one** consumer. `[R]` It acquires more the moment
`fission.py` / `n2n.py` stop hardcoding `L = 0` — which is exactly what #426
proposes for `(n,2n)`. **That is the coupling that paused Campaign 2.**

`[R, unverified]` The frames memo additionally flags **energy condensation on
NON-NESTED grids** as a live second consumer of the same shape (`OverlapBasis`
declares `PARTITION_OF_UNITY` while its trial Gram measures DENSE,
`frame.py:206-208`). The explorer's census did not confirm an l-stacked per-mode
apply there. **Unresolved — check before relying on either statement.**

---

## 7. THE PLAN

Phases are ordered by dependency, not by size. **Effort is a measure of how many
sessions something takes, not of whether it is done** (user, 2026-08-31).

### ⚡ P0 — IMMEDIATE, land regardless of any later ruling

These are true independent of every design choice below.

| # | item | site |
|---|---|---|
| ⚡ 0.1 | **Stop fabricating the support.** `angular_frame` must not write `support=SPACE_SPHERE` onto nodes that are not on `S^2`. | `directional.py:587-593` |
| ⚡ 0.2 | **Correct `metric.py`'s "noise mode".** See §9.1 — it is a null direction, not roundoff, and the sentence is present-tense-FALSE. Re-derive the `_DENSE_METRIC_RCOND` cliff against the repaired frame. | `metric.py:105-119` |
| ⚡ 0.3 | **Correct `Quadrature.spherical_harmonics`'s docstring** — *"the other slots are filled with zeros"* is exactly the broken property. `[M]` this exact claim appears at **exactly 1 site tree-wide** (`grep -rn "filled with zeros" docs/ orpheus/`). | `directional.py:538` |
| ⚡ 0.4 | **Correct the theory page's slab-frame passage** — a DIFFERENT claim from 0.3, ⛔ corrected 2026-08-31 after I published a relayed line reference without checking it. The passage describes the slab frame as merely *"not even diagonal"* and cites P7's *"best diagonal candidate reads 1.806"*; it inherits the framing that omits RANK DEFICIENCY (§9.1). Fix alongside ⚡0.2, since both restate P7's reading. | `docs/theory/foundations/spherical_harmonics.rst:640-652` |
| ⚡ 0.5 | **A real `raise` in `_evaluate_real_sh`** refusing non-unit direction vectors, so the fabricated azimuth is unspellable. ⛔ a `raise`, NOT an `assert` — the canonical runner is `python -O`, which strips asserts (`coding-standards`). | `spherical_harmonic_basis.py` |

⚠ 0.5 will make the CURRENT slab pairing raise. It must land WITH P2, or behind
it — sequencing is P2's §6b obligation, not a reason to weaken 0.5.

### P1 — `Basis` declares its domain; the frame checks it where they meet
**Goal.** An incompatible basis–measure pairing cannot be constructed.
**Proposed means** (2026-08-31, unverified): a `domain` property on `Basis`
(closing leak L2), and a composability check at frame construction (L3), reading
the measure's own support (which requires ⚡0.1 first).
**Done when:** constructing `GalerkinFrame(SphericalHarmonicBasis(L=2), slab_measure)`
RAISES, with a test witnessing it.

### P2 — the angular basis is a consequence of the measure's symmetry
**Goal.** The basis a rule binds is derived from what the rule's measure IS, not
from its construction history.
**Proposed means** (unverified): generalise `MirrorEvenSphericalHarmonicBasis`
to `InvariantSphericalHarmonicBasis(group)` — mask DERIVED by probing, per §5.3
— with `MirrorEven` becoming an instance; and re-key `_harmonic_basis` off the
measure's symmetry rather than `folded_by` (closing L5).
**User ruling, 2026-08-31:** *"we're going to implement partially a symmetry
deriving and quotient machinery to base this upon"* — so P2 depends on that
machinery, which the user will charter separately.
**Done when:** §0's done-when 1 and 2 both hold.

⚠ **OPEN ONTOLOGICAL QUESTION, deliberately unresolved.** Two encodings:
* **(a) layout-preserving** — keep the rectangular `(L+1, 2L+1)` table and
  structurally zero the non-invariant columns, exactly as `MirrorEven` does
  (which already zeroes the `|m| > l` padding). Consumer-transparent.
* **(b) true `LegendreBasis` on `[-1,1]`** bound to the rule's OWN scalar
  measure — no 3-vector embedding, no relabel, no zero columns; the basis's
  domain IS the measure's support and nothing is fabricated anywhere. Changes
  the moment-tensor shape to `(N, L+1)`.
⟹ (b) is ontologically cleaner; (a) is precedent-shaped. **Not ruled.** If (a)
is taken, its docstring MUST say the layout is an ENCODING and name (b) as the
successor, or the encoding will later read as the ontology.

### P3 — the frame carries a well-posedness verdict
**Goal.** A pairing that cannot support a per-mode operator says so.
**Proposed means** (unverified): tiers 2/3 of §5.2.
**Deferred by user ruling** to a later session. ⚠ §5.2's arming caution applies.

### P4 — a 1-D rule stops claiming a 3-D direction
**Goal.** No object exposes a property that is not its own (leak L4).
`[M]` 171 occurrences over 20 files. **Explicitly NOT excused on effort** —
this is a multi-session item, not a rejected one.

### P5 — re-key the collateral consumers
`[M]` **12 sites repo-wide consume the slab-GL-L2 frame as a fixture**
(7 tests, 4 docs, 1 production comment). In particular
`tests/numerics/test_frame.py:906` asserts the live `l=2` Gram diagonal is
`[0.4, 0.8, 0.8]` at `rtol=1e-12` — **the two `0.8`s ARE the fabricated
`Y_2^1` / `Y_2^2` slots.** A green test PINS the defect; it must be re-pinned,
not deleted (`coding-standards`: retirement means test migration).

---

## 8. ⛔ REFUTED PREMISES — kept in place per `plan-authoring` §3

### 8.1 ⛔ "The measured Gram (pseudo-inverse) fixes it for every configuration"
**Whose:** the user's, ruled 2026-08-31 as the fix shape, with the explicit
instruction *"let's try it and see if it fixes"*.
**Refuted the same day, by measurement.** `G+` fixes every case where the Gram
is merely DENSE — including two I would have missed. It **cannot** fix rank
deficiency: `[M]` on slab GL8 L=2 the min-norm coefficients are
`c = 0.800*Y_0 + 0.200*Y_2^0 + 0.346*Y_2^2`, which reconstruct psi **exactly**
(2.2e-16) yet are not `e_0`; `Lambda` then scales the `l=0` part by `Sigma_s0`
and the `l=2` parts by `Sigma_s2`. Independently corroborated by §4.3.
**What survives:** the instinct was right that ONE general mechanism should
cover every configuration — §5.3 supplies it. And the DENSE-Gram route remains
the correct treatment for the two dense-but-full-rank rows.

### 8.2 ⛔ "rank(G) == #live is the well-posedness predicate"
**Whose:** mine.  **Refuted:** it is SUFFICIENT, not necessary — as a refusal it
over-rejects (the padding; harmless intra-block degeneracy).  The necessary and
sufficient condition is §5.1's kernel-invariance.
**But** `[M]` over the 8 shipped test frames, rank-over-LIVE-slots is a
**perfect discriminator**: slab-GL8-L2 is the only rank-deficient row, while
`product(4,4)-L2` and `LS4-L3` are DENSE for legitimate reasons and full rank.
⟹ **rank is an excellent DETECTOR and a bad REFUSAL predicate.** A guard keyed
on DENSE-ness would fire on healthy frames.

### 8.3 ⛔ "the redundancy is countable a priori: every even m>=2 slot is redundant"
**Whose:** mine. **Refuted out of sample:** `[M]` over
`n_ordinates in {2,4,8,16,32} x L in 1..6` — correct on **0 of 30 rows**.
Two errors: `(1,1)` is DEAD not live (the `l=1` Cartesian special case sets
`Y[:,1,2] = mu_y = 0`), and the true rank saturates at `2L+1` for `L >= 3`
(measured 7, 9, 11, 13 at L = 3,4,5,6).
⟹ **I derived the rule from a basis I had not READ.** Reading
`_evaluate_real_sh` gave the mechanism in one pass and explains both errors.

### 8.4 ⛔ "folded_product is spurious at l>=1"
**Whose:** mine, published in a table earlier the same day.
**Refuted:** measured with a PLAIN `SphericalHarmonicBasis`, not the basis
production binds (`MirrorEvenSphericalHarmonicBasis`). On production's own path
`folded_product` is **clean at L <= 2**. The fold machinery works.
⟹ vv-principles #28: build the operand the way production builds it.

---

## 9. COLLATERAL FINDINGS

### 9.1 ⛔⛔ Campaign 1's P7 measured this null direction and labelled it "noise"

`orpheus/numerics/metric.py:111-112` reads, verbatim:

> `[M]` 2026-08-30 on the flagship slab Gram (singular values
> `2.71 / 1.42 / 4.92e-1 / 4.74e-2` live, one `~1e-16` **noise** mode)

`[M]` my measurement of that exact frame: `2.707550, 1.419220, 4.924502e-1,
4.744684e-2, 9.483994e-18` — **bit-for-bit the same four**. And the fifth is
**not noise**: its right-singular vector `[-0.4472, 0, +0.4472, 0, +0.7746]`
= `(1/sqrt5)(-Y_0 + Y_2^0) + sqrt(3/5)*Y_2^2` evaluates to **4.98e-16 at every
node** — an exact linear dependence.

⟹ The defect had been SEEN, one day before it was found, and mislabelled. The
comment then reasons carefully about a *"cliff's ANALYTIC edge"* that treats
`sigma_3 = 4.74e-2` as the smallest GENUINE live mode, and pins
`_DENSE_METRIC_RCOND = 1e-12` to sit *"~4 orders above the noise floor"* — i.e.
**the pin was tuned to discard the signal that reports this bug.**
⭐ The transferable half: the measurement was right and the INTERPRETATION
inverted it, and the rigour of the surrounding analysis is exactly what made it
unquestionable. This is `plan-authoring` §2's `[M]`-scope defect wearing a
carefully-derived cliff.

⚠ **Consequence for P7's own conclusion:** the DENSE-metric repair makes
Parseval a theorem AND **removes the singularity from the only property that
could report it**. P7 must not be read as having addressed this defect; the two
repairs are mutually exclusive in intent (fix the metric vs fix the basis).

### 9.2 The diagnostics that exist are permissive, not inert
`[M]` `discrete_gram_structure` has exactly **one** production read
(`frame.py:286`) and it **routes a metric choice, never refuses**.
`Basis.gram_structure`'s `NotInvertible` cannot fire for SH (declared
`DIAGONAL` unconditionally, and it sits on `project`/`gram`, which the
scattering path never calls — 0 of 38 hits in `transport/operators/`).
**Nothing anywhere computes a rank or condition number of a frame Gram.**
`Basis.mass_matrix` is a genuinely inert authored diagnostic (0 production
callers).

### 9.3 A predicted axis-convention collision, one probe from settled
`[M]` `SubgroupOfO3.SO2.is_invariant(gauss_legendre(8).measure)` returns
**`False`**. `[R, frames memo]` `_embedded_nodes` uses **column 0** as the polar
axis while the `SO2`/`Dinfh` realizations use **z** (`symmetry.py:800-819` vs
`:889-899`, `:1014`), so the predicted correct answer is `True`. **The `False`
is measured; the "should be `True`" is reasoned and UNVERIFIED.** This is
directly load-bearing for P2 — the symmetry machinery must agree on an axis
before it can derive a quotient.

### 9.4 The precedent, every adjective verified by reading
`MirrorEvenSphericalHarmonicBasis` is a genuine `@dataclass(frozen=True)`
**subclass** (frozen verified by attempted mutation); field
`mirror_axis: int = 1` **has a default**; **no registry, no factory
constructor**; keeps the full `(L+1, 2L+1)` layout and zeroes columns; the mask
is a `@cached_property` **DERIVED** by evaluating `_evaluate_real_sh` at 5
hardcoded generic directions and their mirrors, raising `RuntimeError` if a slot
classifies as neither even nor odd.
`Quadrature._harmonic_basis` selects on **one field, `folded_by`** — and `[M]` a
slab GL rule has `folded_by is None`, so **the precedent's hook cannot express
the slab case as written** (leak L5).

---

## 10. WHAT THIS PLAN DOES NOT ESTABLISH

* **No production line has been changed.** Every measurement above is either
  read-only or an in-process monkeypatch with an asserted restore.
* The candidate repair is validated ONLY on §1.1's fixture:
  `[M]` in-process m=0 restriction gives L=0,1 **bit-identical** to production
  (`d = +0.000e+00`) and L=2,3 = **+4.000000000000** (6.7e-16), with both
  controls passing (patch bites on slab, `max|dY| = 0.837`; **exactly inert on a
  3-D rule**, `max|dY| = 0.0`). It has NOT been run against the suite.
* Whether encoding (a) or (b) of P2 is right is **unruled**.
* §6's energy-condensation second consumer is `[R]`, not confirmed.
* §9.3's "should be True" is `[R]`, not measured.
* The two dense-but-full-rank rows (`level_symmetric(4)` L=3,
  `folded_product(2,4)` L=3) are a SEPARATE class — *degree L asked of a rule
  too coarse to resolve it* — deferred by user ruling, with the Q6 exactness
  machinery named as its natural home.

---

## 11. ▶ RESUME SURFACE

**This plan blocks** `.claude/plans/cs4c_binding_design.md` §18.6 step 1
(#426's missing (n,2n) measurement), which is PAUSED. The coupling is §6: #426
proposes making `(n,2n)` carry `l >= 1` moments, and `n2n.py:200` currently
hardcodes `L = 0` — which is the only reason that channel is not already
affected.

**Next act, per the user (2026-08-31):** the user will charter a closely-related
plan for **symmetry derivation and quotient machinery**, which P2 is to be built
upon. Do not start P2's design before it lands.
**⚡ P0 is unblocked and should land first**, regardless of any ruling in that
plan.

**Supporting artifacts** (⚠ all in UNTRACKED `scratch/` — a `git clean`
destroys them; their load-bearing content is carried in this file):
`_pl_slab_defect_repro.py` (the reproducer), `_pl_slab_fix_probe.py` (the
in-process candidate-fix probe), `n2n_pl_diagnosis_main.md`,
`n2n_pl_consistency_literature.md`, `n2n_pl_frames_attack.md`,
`n2n_pl_blast_radius.md`.
