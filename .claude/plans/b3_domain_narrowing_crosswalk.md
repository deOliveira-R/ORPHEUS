# B3 crosswalk — narrowing the SN boundary law's domain to `Γ₊`

Mandatory pre-carve artefact per `coding-elegance` Pattern 7 / `[[lessons-L17]]`:
a carve that crosses subsystem boundaries starts with the convention crosswalk,
**written before any code**. The crosswalk IS the architecture; the code is its
transcription.

Campaign: `.claude/plans/boundary_machinery_review.md` phase **B3**.
Branch `refactor/operator-strategy-layers`. Author: main agent, 2026-07-31.

Evidence key (same as the review): `[R]` read in code · `[M]` measured by
running · `[G]` grep/graph absence · `[U]` unverified, pending.

---

## 0. The change, in one line

**[R]** The physics is `γ₋ψ = R G γ₊ψ + q` with `G : Γ₊ → Γ₊`, `R : Γ₊ → Γ₋`
(`docs/theory/foundations/boundary_conditions.rst:194-215`). SN types its
realized law `full-face → full-face` and discards the outflow rows at the
consumer; diffusion types it `Γ₊ → Γ₋` and has nothing to discard. B3 moves SN
onto the diffusion shape.

---

## 1. The convention table

| subsystem | input convention | internal convention | output convention |
|---|---|---|---|
| **SN today** `_reflect_trace` (`sn/operators/boundary.py:352-359`) | the WHOLE face slot — all `quad.N` ordinate rows (`boundary.face_view(face)`) | realized law is a full-`N` `TensorProductOperator` (`PermutationOperator(perm, axis=0) & I`, or `IncomingOrdinateMaskTensor(...) & I`) | full face slot; **caller slices** `out[sel_inflow] = full[sel_inflow]` — outflow rows discarded |
| **SN after B3** | the **outflow half** `Γ₊` — `|O_f|` rows | law is rectangular `Γ₊ → Γ₋` | the **inflow half** `Γ₋` — `|I_f|` rows; caller writes the inflow view, nothing to discard |
| **diffusion (reference, unchanged)** (`diffusion/operators.py:612-614`) | `trace.outflow_view(face)` — the outflow half | scalar `𝒜` amplitude | written to `ScalarTraceSpace.INFLOW_ROW` |
| **Bridge** | **producer-side** (`_reflect_trace` restricts before calling) | the realizer emits the narrowed operator | **producer-side** (the law's codomain IS `Γ₋`) |

The Bridge row is the Pattern-7 ruling: the narrowing lives at the **definition
site** (realizer + `_reflect_trace`), never re-applied per consumer. **[M]** B2
measured that `_reflect_trace` is the sole SN consumer of a realized law's
`apply`, so "every consumer" is one site — but the point of siting it at the
producer is that the *next* consumer inherits the honest domain for free.

---

## 2. The axis convention — the sharp edge

**[R]** The ordinate set at a face is a **three-way** partition, not two-way
(`docs/theory/foundations/boundary_conditions.rst:3130-3140`,
:eq:`ordinate-partition-inflow-outflow`):

```
{1..N}  =  I_f  ⊔  O_f  ⊔  T_f
          Ω·n < −ε   Ω·n > +ε   |Ω·n| ≤ ε   (tangential)
```

So "not inflow" ≠ "outflow". Two decisions follow, and **both are convention
choices a mis-step turns into a Mode-5 index bug**:

- **D-1. Does `Γ₊` mean `O_f` (strict) or `O_f ∪ T_f`?** **[R]** The vacuum mask
  preserves `O_f ∪ T_f` (its own equation, `vacuum-legacy-vs-trace-correct`),
  whereas `AngularTraceSpace.outflow_indices_for_face` presumably selects strict
  `O_f`. **[U]** If they differ, the narrowing silently changes which rows the
  law sees. The theory page's mitigation — *"ORPHEUS's quadrature adapters carry
  every tangential ordinate at μ = 0 so ψ = 0 on `T_f` for a properly-initialised
  flux"* — is a **claim about initialisation, not a structural guarantee**, and
  per `[[lessons-L33]]` a doc is a claim. **MUST be measured before the carve.**
- **D-2. Is `|I_f| == |O_f|`?** For an axis-aligned specular mirror the
  reflection index table is an involution pairing each inflow ordinate with
  exactly one outflow ordinate, so yes — but **[U]** unverified for the
  curvilinear / level-symmetric / product quadratures, and a face where it fails
  makes the narrowed reflective law non-square.

---

## 3. Blast radius — what the narrowing touches

| layer | site | what changes |
|---|---|---|
| trace carrier | `AngularBoundaryFlux` (`transport/fields/`) | **[G]** has only `face_view`; the scalar sibling already has `outflow_view`/`inflow_view` (`scalar_boundary_flux.py:110,119`). B3 gives the angular carrier the same two verbs — *cross-method harmonisation, one vocabulary* |
| numerics | **[G] no restriction/gather primitive exists** (class list at `numerics/operator.py`: `Permutation`, `IncomingOrdinateMaskTensor`, `PeriodicWrap`, `TensorProduct`, `Diagonal`, `RankOne`, …) | the `[G]` the review predicted — the projector the subsystem needs three times has no name |
| SN realizer | `sn/boundary/realizer.py:136-300` | each law's operator retypes; the reflective `perm` must be **re-indexed into outflow-local coordinates**; `IncomingOrdinateMaskTensor & I` for vacuum becomes the honest **zero map `Γ₊ → Γ₋`** |
| SN consumer | `sn/operators/boundary.py::_reflect_trace` | the slice-write dissolves; `apply_transpose`'s `masked[sel] = face_in[sel]` dissolves with it |
| law invariants | `assert_realizable(quad, inflow_indices=…)` | **[U]** phrased in full-`N` ordinate space — does the narrowing invalidate any? |
| **C-1 gates** | 4 outflow-must-be-zero tests + the 2-D balance gate | re-posed, never deleted — see the gate-reposing memo |

---

## 4. The algebra B3 restores — and why the primitive is a *restriction*, not a *projector*

The review asks for "one named projector". The sharper decomposition, and the
one that matches the theory page's own vocabulary, is a **restriction/injection
pair** — the trace operators `γ_±` the affine form already names:

```
γ₊ : Γ → Γ₊        restriction (gather)          ι₊ = γ₊ᵀ : Γ₊ → Γ   injection (scatter)
γ₋ : Γ → Γ₋        restriction (gather)          ι₋ = γ₋ᵀ : Γ₋ → Γ
```

Every spelling §5 of the review found is then a **composition**, not a
primitive:

| §5's finding | is |
|---|---|
| `_reflect_trace`'s slice-write `out[sel] = full[sel]` | `ι₋ ∘ γ₋` = `P_in` |
| `IncomingSourceOperator._inflow_mask`'s dense multiply | `ι₋ ∘ γ₋` = `P_in` (**[M]** bit-identical to the above, per §5) |
| `IncomingOrdinateMaskTensor` | `ι₊ ∘ γ₊` = `P_out` = `I − P_in` |
| the whole face law | `ι₋ ∘ (R G) ∘ γ₊` |

`P_in ∘ P_out = 0` stops being a curiosity and becomes `γ₋ ∘ ι₊ = 0` — two
disjoint index sets, true by construction.

**Numerically the primitive already exists as a numpy verb**: `PermutationOperator.apply`
is `np.take(x, perm, axis)`; a restriction is the SAME call with a non-square
index array. So the new type is a **generalisation of the permutation, not a
new mechanism**. It must NOT be spelled as a subclass, though — a permutation
carries invariants a restriction does not have (involution detection,
invertibility, an algebra-closed `inverse()`); a restriction is
rank-deficient and has a **scatter** transpose, not an inverse.

---

## 5. Two shapes for the carve — the decision the user makes

- **Option A — bake the narrowing into the law** (what B3 as written says). The
  realized operator genuinely IS `Γ₊ → Γ₋`. Maximum honesty; §4.2's typing
  defect is actually fixed; vacuum's realized operator becomes the zero map and
  the mask's whole "preserves the outflow" justification dissolves. Cost: the
  realizer re-indexes every law's table into outflow-local coordinates — the
  Mode-2/Mode-5 exposure sits exactly there.
- **Option B — compose the restriction around a still-full-`N` law**:
  `B_face = ι₋ ∘ γ₋ ∘ law ∘ ι₊ ∘ γ₊`. Gives the projector a name and makes the
  discard typed instead of a slice-write, but the law's *declared* domain stays
  full-face, so the review's root finding is documented rather than fixed.

**Recommendation: A.** B is a half-measure that leaves the defect in place while
adding a layer — it fails `coding-elegance`'s concept-count test (concepts go
up, not down). A is what the plan, the theory page, and the diffusion precedent
all say. The standing user guidance is explicit that effort and risk are not
deferral reasons.

**Rejected outright — Option C**, re-minting `InflowTraceSpace` /
`OutflowTraceSpace` as separate typed spaces. **[R]** That reverses the ratified
#205/#201 collapse to *one space, two selectors*
(`boundary_conditions.rst:~300`). The restriction is a **verb**, not a second
space.

---

## 6. The four opening questions — ALL ANSWERED (measured 2026-07-31)

1. **D-1 — strict `O_f`.** `[M]` And the two spellings genuinely differ:
   `outflow_indices_for_face` is strict `O_f`, the vacuum mask preserves
   `O_f ∪ T_f`. On a real CYL mesh under `product(n_mu=2, n_phi=4)`, **4 of 8
   ordinates at `xmax` are tangential**. So `I − P_in ≠ P_out`, and `P_out` is a
   **fourth map no operator type spells**. `IncomingOrdinateMaskTensor`'s
   docstring claim *"projection onto the outflow subspace"*
   (`numerics/operator.py:2319-2321`) is **measurably false**; its code
   (`:2374`) is honest.
2. **D-2 — `|I_f| == |O_f|` on every reachable fixture** `[M]`, but the
   partition is **disjoint, NOT exhaustive, NOT contiguous**: `gauss_legendre(5)`
   carries 1 tangential, `product(2,4)` carries 4 of 8, every `lebedev` carries
   4–8, and `level_symmetric(4)` is non-contiguous with 0 tangential. Only
   `gauss_legendre(4)` is the clean two-way case — **the slab is the unrepresentative
   fixture**.
3. **Bit-identity HOLDS on every production-reachable path** `[M]`. Reflective is
   exact on SLB/SPH/CYL because `perm[inflow] ⊆ outflow` — which is the **ERR-045
   invariant the law layer already certifies at realization**, so B3 *consumes* an
   existing guarantee and adds no new assumption. Vacuum trivially. White,
   albedo(α≠0) and periodic do **not** narrow neutrally — see §7.
4. **The C-1 gates stay reddenable, but their teeth change** `[M]`. The original
   bug (the law's outflow image leaking) becomes **unspellable**; only the
   write-target family survives. Two consequences: the outflow leg must **widen
   from `outflow` to `complement(inflow)`** or it is blind to the tangential rows
   (measured: `cyl_reflective` discards a **1.846** law image there), and a
   reduced `PermutationOperator` **does not validate its domain** — a full-face
   input returns wrong values with **no raise**, so the narrowed operator owes an
   explicit domain guard with its own negative test.

**Corollary that killed one of §3's proposals:** a narrowed `Γ₊` domain **cannot
be a view/slice — it must be a gather + scatter.** The diffusion precedent does
NOT transfer: `ScalarTraceSpace.OUTFLOW_ROW` is a *constant* index and a real
memory-sharing view, while the angular set is a quadrature-dependent fancy index,
so any SN `outflow_view` returns a **copy**. Same verb, opposite mutability — the
convention trap Pattern 7 exists to prevent. **The restriction is an OPERATOR,
never a `_view` accessor.**

---

## 7. ⭐ The G/R correction — USER RULING 2026-07-31

The narrowing exposed a **misassignment between the two factors B1 minted**, and
the campaign folds the correction into B3 (*clean before extending*: B3 is the
phase that makes the domains real, and B4 is about to compose `R ∘ G`).

### 7.1 The discriminator

**`G` is the composition (Koopman) operator of a measure-preserving bijection of
the boundary phase space** — `(Gψ)(x) = ψ(g⁻¹x)` for `g` acting on `∂Ω × S^d`.
Such operators are invertible, preserve the trace measure `|Ω·n| dΩ dA`, form a
group, and — the decidable test — are **multiplicative**:
`G(ψ·φ) = (Gψ)(Gφ)`. A relabeling satisfies it; **an averaging operator never
does.** Everything that fails multiplicativity is a kernel, and kernels are `R`.

In the user's framing, which this formalises: a change of direction caused by the
**geometry** is `G`; a change of direction caused by the **constitutive
assumption of the BC** is `R`.

### 7.2 The typing correction — the crossing is GEOMETRIC

**[R]** Today: `BoundaryGeometryMap` is documented `G : Γ₊ → Γ₊`, *"an
endomorphism of the outflow trace"* (`geometry/boundary/_factors.py:100-107`).

That forces the specular mirror **into `R`**, because `Ω ↦ Ω − 2(Ω·n)n` maps
`Γ₊ → Γ₋` — it is the unique ambient isometry fixing the face, it exchanges the
hemispheres, and it preserves `|Ω·n|`. The honest typing is therefore

```
G : Γ₊ → Γ₋     the deck transformation (the crossing)
R : Γ₋ → Γ₋     the constitutive kernel
```

### 7.3 Why the misassignment went unnoticed — a theorem

If `R` is rank-one, `R = u ⊗ v`, then `R ∘ G = u ⊗ (Gᵀv)`. The Lambertian's
`v = |Ω·n|` is **invariant** under both the mirror and the periodic translation,
so `R ∘ G = R`: **`G` is unobservable exactly when `R` is rank-one.** White is
precisely that case — its `G` slot was free, and the physics drifted into it.

**[R]** The code already concedes it in prose: `HemisphericalAverage`'s own
docstring says *"an all-to-all coupling, **not a relabeling**"*
(`_factors.py:251-252`) while implementing the Protocol documented as *"it only
relabels"*.

### 7.4 The corrected taxonomy

| law | `G` — deck transformation `Γ₊→Γ₋` | `R` — constitutive kernel on `Γ₋` |
|---|---|---|
| vacuum | canonical crossing (immaterial, `R=0`) | `0` |
| reflective | **mirror** | `α·I` |
| white | mirror (immaterial — `R` rank-one) | `α ·` Lambertian |
| albedo | mirror (immaterial if `R` rank-one) | `α ·` **re-emission closure** |
| periodic | **translation** | `I` |
| prescribed | — (no linear part) | `0`, plus `q ∈ Γ₋` |

Reflective and white now differ **only in `R`**; reflective and periodic differ
**only in `G`**. Orthogonal — which is the whole point.

### 7.5 What this re-diagnoses

- **Albedo's gap is in `R`, not `G`** (this corrects an earlier reading in this
  very file's §5). `ScalarResponse(α)` gives the *magnitude* but not the *angular
  distribution* of the returning flux: complete on a scalar trace (one dof),
  incomplete on an angular one. **[R]** The theory page already anticipates the
  tier — *"R is a scalar amplitude … or a full angular kernel in general weak-form
  BCs (deferred)"* (`boundary_conditions.rst:215-221`). **USER RULING: mint the
  kernel tier now.** Then `albedo(α, isotropic) ≡ white(α)` and
  `albedo(α, specular) ≡ reflective(α)` become **theorems**, not coincidences.
- **Periodic is buildable and `G` is where its cross-face-ness lives** — the
  translation maps `Γ₊(f')` to `Γ₋(f)`, and `|Ω·n|` is preserved because
  `n_f = −n_{f'}`. B1 already gave the law its `axis` field so the realizer
  derives the partner from the installation face.
- **ERR-042 needs the bijection reading.** `BoundaryGeometryMapNotMeasurePreservingError`
  must test the *multiplicativity / point-bijection* property, not flux
  conservation — the Lambertian **conserves total flux** while being no bijection
  at all, so a flux-conservation reading passes the very thing the invariant
  exists to reject.

### 7.6 The quotient picture — and which law is the orbifold

`G` is the **deck transformation of the quotient by which the physical domain is
represented**:

| BC | quotient | fixed points | what it is |
|---|---|---|---|
| periodic | `ℝᵈ/Λ` by a translation | none — **free** | a torus; a genuine **covering space**, a manifold |
| reflective | by a reflection | the mirror plane | an **orbifold** (Thurston *reflector* boundary) |
| rotational (⅛-core) | by a finite rotation | the axis | an **orbifold** (cone points) |

So the orbifold label attaches to **reflective**, not to periodic — periodic is
the free/covering-space case. And `R = I` **exactly when** the BC is a pure
symmetry statement adding no physics. Vacuum, white and albedo are not symmetry
statements at all, which is why their `G` is immaterial and all their content
sits in `R`.

---

## 8. B3 as restructured (tasks #10–#16)

| step | content |
|---|---|
| **B3.0** | the §7 G/R correction + mint the response-kernel tier |
| **B3.1** | mint `TraceRestrictionOperator` (gather; scatter transpose; **domain guard**) |
| **B3.2** | narrow the SN law's domain to `Γ₊` |
| **B3.3** | retire `IncomingOrdinateMaskTensor` (10-test suite migrates) |
| **B3.4** | build albedo (`R`-closure) and periodic (translation `G`) properly |
| **B3.5** | re-pose the C-1 gates RG-1…RG-5; **promote the mutation harness into the repo** |
| **B3.6** | theory-page corrections + the two measured live falsehoods |

---

## 9. B3.1 design — `TraceRestrictionOperator`

`[G]` Measured absent from all of `orpheus/numerics/`. The mechanism already
exists as a numpy verb: `PermutationOperator.apply` is `np.take(x, perm, axis)`,
and a restriction is the same call with a **non-square** index array. It is a
sibling, **never a subclass** — a permutation carries involution detection,
invertibility and an algebra-closed `inverse()` that a rank-deficient
restriction cannot honour.

```
γ_S : Γ → Γ_S        apply           = np.take(x, indices, axis)
ι_S : Γ_S → Γ        apply_transpose = scatter (zeros, then out[indices] = x)
```

- **Adjointable, never invertible.** `ι_S ∘ γ_S = P_S` (the projector §5 found
  spelled twice) and `γ_S ∘ ι_S = I` on the restricted space. Both are
  identities worth a test — the intrinsic-properties standard: a type carrying
  a mathematical concept ships a test of its *defining laws*, not only its
  usage.
- **Construction guards (Pattern 4).** Indices in `[0, n_total)`; **unique**
  (duplicates are not a restriction); and **sorted**.
- **⭐ Sortedness is what makes the N8 trap unspellable.** The narrowing needs
  to remap a subset of global rows into positions within the restricted space.
  `[M]` The naive `arange(sel.size)` is *exactly correct in 1-D* and wrong in
  2-D, and the existing split gate is built on a **sphere**, so it stays green
  — only end-to-end 2-D solves catch it. The correct form is
  `searchsorted(indices, sel)`, which **requires sorted indices**. So: enforce
  sortedness at construction and put the remap on the operator itself
  (`to_local(...)`), and no call site can ever spell the wrong one. That is
  Pattern 5 — build the primitive, not the product — applied to a trap the
  test-architect measured before it was written.
- **The domain guard is mandatory, with its own negative test (RG-3b).** `[M]`
  A reduced `PermutationOperator` silently accepts a full-face input and
  returns wrong values with **no raise**. `apply` validates
  `x.shape[axis] == n_total`; `apply_transpose` validates
  `x.shape[axis] == len(indices)`.
- **It makes `OperatorProduct`'s guard bite.** `[G]` Every realized boundary law
  currently advertises `domain = codomain = None`, so the composition guard
  (`operator.py:1560-1565`) is **vacuous**. B3 is the first change that gives
  those operators real spaces — B4 then composes `R ∘ G` under a guard that
  actually checks something.

### 9.1 B3.1 RESULT ✅ — and the N8 trap has **two** sites, not one

`TraceRestrictionOperator` landed with 18 tests (the first nine on its
**defining laws**: `γι = I`, `ιγ = P` idempotent+symmetric, `ι = γᵀ`
materialised, disjoint annihilation, three-way partition resolving `I`).
pyright 1 (#288 floor); numerics 956 passed.

**[M] Fitness probe on real meshes** — the narrowed composition
`g_out.apply(full_face)[g_out.to_local(perm[inflow])]` is **bit-identical** to
today's `np.take(full_face, perm, 0)[inflow]` on both fixtures:

```
SLB gauss_legendre(4): N=4 |in|=2 |out|=2 |tan|=0
   perm[inflow]=[3, 2]  outflow=[2, 3]
   local positions=[1, 0]  naive arange=[0, 1]      <-- DIFFER
   bit-identical to today: True      TensorProduct fold ok: True
CYL product(2,4):      N=8 |in|=2 |out|=2 |tan|=4
   perm[inflow]=[0, 4]  outflow=[0, 4]
   local positions=[0, 1]  naive arange=[0, 1]      <-- agree
   bit-identical to today: True      TensorProduct fold ok: True
```

⭐ **The remap appears at TWO distinct sites, and they are discriminated by
DIFFERENT fixtures — neither gate can be dropped:**

| site | quantity remapped | naive `arange` is wrong on | agrees on |
|---|---|---|---|
| **reflective narrowing** (B3.2) | `perm[inflow]` → positions in `Γ₊` | **the SLAB** — the mirror *reverses* order, so `[3,2] → [1,0]` | CYL, where `perm[inflow]` happens to be `outflow` in order |
| **schedule split** (`rows`, RG-5) | `sel ⊂ inflow` → positions in `Γ₋` | **2-D** — the lower-half rows are not a prefix of `inflow` | 1-D, where they are |

The test-architect measured the second and reported the trap as "correct in
1-D, wrong in 2-D". That is right *for the split*; for the **reflective**
remap the slab is the discriminating fixture and 1-D is exactly where it
bites. Recording both so a future reader does not conclude 1-D coverage is
sufficient — or that 2-D coverage is.

`to_local` on the operator is what makes both unspellable at the call site.

---

## 9.2 B3.2 RESULT — and three of my own claims corrected

`SNMethodSpace` gained an **`outflow_indices`** field (sibling of
`inflow_indices` — a law's domain deserves the same hand-buildability its
codomain has); `_reflect_trace` composes `ι₋ ∘ law ∘ γ₊`, transpose
`ι₊ ∘ lawᵀ ∘ γ₋`; vacuum realizes to the **zero map**, reflective to a
permutation on the **reduced** axis.

**[M] No value regressions.** The wide baseline reproduced *exactly* the 49
construction-test failures and nothing else — zero solver, adjoint-reciprocity,
DSA, Gauss-Seidel, snapshot or convergence reds. The bit-identity claim holds.

**[M] Bit-identity gate has teeth, falsified independently:** forcing the
naive `arange` remap reds **6 rows** (slab + sphere, both directions) with a
positive control proving the mutation bit (26 interceptions). `cyl` and 2-D
stayed **green** — the §9.1 complementarity, reproduced.

### Three claims of mine that were wrong

1. **"Three un-narrowed laws" — it is SIX rows across FOUR kinds.** `albedo`
   at **α=0 and α=1 too** (their fast paths return a bare
   `ZeroOperator`/`IdentityOperator`, which are *endomorphisms*), plus
   `prescribed_inflow` (now raises `ordinate axis mismatch`). Unreachable in
   production — the SN registry admits only `{vacuum, reflective}` — so the
   tree stays green, but the B3.4 todo is twice the size I stated.
2. **A shape assertion cannot detect `Γ₊→Γ₊`.** `[M]` `|Γ₊| == |Γ₋|` on
   **every quadrature × face in the tree**, so an un-narrowed endomorphism has
   exactly the right output shape. Only the anti-Mode-12 leg found the three
   above. This is the invariance-group lens biting inside B3 itself: the
   measured functional (shape) has the error class (`Γ₊→Γ₊` vs `Γ₊→Γ₋`) in its
   stabiliser.
3. **The narrowed law does NOT validate its own domain.** `[M]` Both vacuum's
   `ZeroOperator` and reflective's `TensorProductOperator` accept a **full-face
   input** and return `|Γ₋|` rows of *wrong values, no raise*. The guard I put
   on `TraceRestrictionOperator` does not travel to the operator the realizer
   *emits*. Unreachable through `_reflect_trace` (which always feeds
   `γ₊.apply(...)`, itself guarded), but reachable through `sn.bc[face].apply`.
   **Deliberately left as 4 strict-xfails for B3.4**, which restructures the
   realizer around `R ∘ G` for all seven laws and is the one place the guard
   belongs — a per-law guard added now would be seven copies of what B3.4
   collapses into one.

### Consequences to carry

- **A narrowed law can no longer compose with an un-narrowed one** —
  `0.3·specular + 0.7·white` raises on both domains. 6 mixed-BC gates are
  B3.4-blocked, which argues for **B3.4 before B3.3**.
- **`SNMethodSpace.minimal` is now a partial constructor** and after B3.4 can
  realize nothing — a retirement candidate, not a fixture.
- **16 strict xfails = the todo list**, each `--runxfail`-verified to red for
  ITS documented reason, and the B3.4 subset proven to **flip** under a landing
  simulation. (Trap met en route: a simulation built on `ScaledOperator(0.0, …)`
  **refuses a zero scalar**, so the xfail swallowed a `ValueError` and the row
  falsely read "didn't flip" — Mode 8 class 4, caught by the flip check itself.)
- ⚠ **The mutation harness is at `scratch/b3_2_mutations.py`.** `scratch/` is
  tracked-capable (44 files under it are tracked), so committing it there ends
  the "lives only in a job tmp dir" failure mode — but scratch is a holding pen
  by convention. **B3.5 owes it a discoverable home.**

---

## 9.3 ⏸ COMPACTION POINT — the B3.4 brief

**State: B3.0 `9e2139b4` · B3.1 `b39502f8` · B3.2 `7f02de15` · doc repair `b11a2ce3`.
27 commits ahead of main on `refactor/operator-strategy-layers`. NEXT = B3.4, THEN B3.3**
(user ruling 2026-07-31 — the mask retirement is cheap and independent; the un-narrowed
laws are what block the algebra). **Verify all of this against `git log`, never against a
summary.**

### What B3.4 must build — SIX rows across FOUR kinds, not the three first claimed

`[M]` Measured, and the count matters because a shape assertion cannot find them
(`|Γ₊| == |Γ₋|` everywhere):

| law | today's realization | why it is un-narrowed |
|---|---|---|
| `albedo(α=0)` | bare `ZeroOperator()` | endomorphic — no space hooks |
| `albedo(α=1)` | bare `IdentityOperator()` | endomorphism |
| `albedo(α∉{0,1})` | `α·(I & I)` | endomorphism |
| `periodic` | `PeriodicWrapOperator & I` | angular identity-with-copy |
| `white` | `AngularAverageOperator & I` | full-`N` contract-then-broadcast |
| `prescribed_inflow` | `IncomingSourceOperator` | **raises** `ordinate axis mismatch` |

All unreachable from a `BC` tag (SN admits `{vacuum, reflective}`), so the tree is green —
but the **16 strict xfails are the todo list**, `--runxfail`-verified to red for their own
reasons, and the B3.4 subset was proven to **flip** under a landing simulation.

### The three design rulings B3.4 executes

1. **Albedo's gap is an incomplete `R`, not an empty `G`.** `ScalarResponse(α)` gives the
   magnitude but not the angular distribution of the returning flux: complete on a scalar
   trace (one dof), incomplete on an angular one. Give it an explicit **re-emission
   closure**, after which `albedo(α, isotropic) ≡ white(α)` and
   `albedo(α, specular) ≡ reflective(α)` are **theorems**, not coincidences — and the
   specular closure moves its content across into `G`, exactly as the membership criterion
   predicts. **The closure belongs on the LAW**: it is the tier-1-vs-tier-2 distinction
   (§7 / `bc-method-realizability`), and a method that cannot see the choice — diffusion,
   by construction — must not be the one making it.
2. **Periodic's `G` is the translation reading the PARTNER face's `Γ₊`.** B1 already gave
   the law its `axis` field so the realizer derives the partner from the installation face.
   ⚠ **This breaks the per-face block-diagonality** that `_reflect_trace`'s docstring leans
   on to justify the Gauss-Seidel subset restriction being *exact* ("B is block-diagonal
   over faces, so the subset action is the EXACT restriction"). That justification must be
   re-derived or the schedule re-posed — **do not let it pass silently**; it is the phase's
   real work and wants a user checkpoint before a shape is committed.
3. **B3.4 owns the domain guard.** `[M]` The narrowed law does **not** validate its own
   domain — vacuum's `ZeroOperator` and reflective's `TensorProductOperator` both accept a
   full-face input and return `|Γ₋|` rows of *wrong values, no raise*. Unreachable through
   `_reflect_trace` (always fed a guarded `γ₊.apply`) but reachable through
   `sn.bc[face].apply`. Shipped as **RG-3b**, 4 strict xfails. It lands here because B3.4
   restructures the realizer around `R ∘ G` for all seven laws — one place, not seven
   copies.

### What B3.4 unblocks / retires

- **6 mixed-BC gates**: a narrowed law cannot compose with an un-narrowed one, so
  `0.3·specular + 0.7·white` raises on both domains.
- **`SNMethodSpace.minimal` becomes a partial constructor that can realize nothing** — a
  retirement candidate, not a fixture. Decide explicitly (coding-standards: never let a
  superseded path sit half-alive).
- `OperatorProduct`'s composition guard stops being vacuous once the operators carry real
  spaces — B4 then composes `R ∘ G` under a guard that actually checks something.

### Traps carried in (all measured, all cost real time)

- **`ScaledOperator(0.0, …)` REFUSES a zero scalar.** A landing simulation built on it made
  an xfail swallow a `ValueError` and falsely read "didn't flip" — Mode 8 class 4.
- **`to_local`, two sites, complementary fixtures** — §9.1. Neither covers the other.
- **The fourth search**: a signature change to a widely-consumed handle owes a **full-suite
  run BEFORE the commit**. Duck-typed surrogates are invisible to graph, grep and typecheck.
- **Doc correctness is not deferrable** (user ruling, now in
  [[feedback-articulation-lossless-disassembly]]): fix falsified prose in the SAME change,
  present-tense-false is the bug, past-tense history stays.

### Gates B3.4 must clear

`np.array_equal` bit-identity for the laws that already worked; the seven-law domain gate
with its xfails **deleted** (an XPASS(strict) failure is the point); the mutation harness at
`scratch/b3_2_mutations.py` (`ORPHEUS_B32=N1…N9|M1|M2`) still catching what it caught;
full tree `-m "not slow"` at **4335 passed / 6 skipped / 73 xfailed / 0 failed** or better,
minus the xfails B3.4 legitimately deletes; `npx pyright orpheus/` = **1** (the #288 floor);
`sphinx -E -W` exit 0.

---

## 10. Watch item — a sixth dead-capability instance

`[G]` `BoundaryGeometryMap.is_adjointable` and its concrete implementations have
**zero consumers** anywhere in `orpheus/`, `tests/` or `docs/` — measured. That
is the review's §4 pattern (declared capability, no consumer) inside the types
B1 minted two phases ago. It is kept because **B5 is its named consumer** (it
flips when white is typed as the rank-one it is), which is a concrete
architecture rather than a speculative one — but if B5 slips, this is a
retirement candidate, not a permanent fixture.

---

**Two traps carried forward into the implementation:**

- **N8** — the index remap is `searchsorted(inflow, sel)`, **not**
  `arange(sel.size)`. `[M]` The naive form is *exactly correct in 1-D* and wrong
  in 2-D, and the existing split gate (`test_psi_half_coupling.py:795`) is built
  on a **sphere**, so it stays green. Only end-to-end 2-D solves catch it.
- **The fourth search** — `[M]` six duck-typed surrogates, one a real subclass
  (`_BWithStubFace(SNBoundaryOperator)` overriding `_face_laws`) that injects a
  fake law into the production `_reflect_trace` path. Invisible to grep and to
  the graph. A full-suite run belongs **before** the commit, not after.

---

## 11. ⭐ B3.4 DESIGN RULINGS — user checkpoint 2026-08-01

Two forks were put to the user with measured options. Both answers **sharpened the
design rather than picking an option**; both are recorded verbatim below, because
each corrects something I had wrong.

### 11.1 Albedo — the specular closure belongs in `R`, not `G`

> "For albedo with specular reemission as a closure, the specular enforcement should
> probably be in R, not G (which probably has the IdentityMap, just like the
> Lambertian). So I think complete Albedo with specular closure is equivalent to
> Reflective boundary for all practical purposes maybe, but the difference is that in
> the Reflective there is geometrical inversion, and for Albedo there is a Response
> closure imposing it. I think this is part of the discipline of separating precisely
> what is Geometric and what is Response." — user, 2026-08-01

**What this corrects in B3.0.** The multiplicativity test I minted as *the* decidable
criterion is **necessary but NOT sufficient**. A specular *kernel* is a permutation,
hence multiplicative — so multiplicativity alone cannot separate "a polished wall
returning α specularly" from "a symmetry plane". The sufficient condition is the one
already written into `_factors.py`'s quotient table and then not used as a test:

> :math:`G` is the deck transformation **of an actual quotient of the physical
> domain**. A physical surface is not a quotient — the domain does *not* continue on
> the other side — so its specular pairing is **constitutive**, i.e. :math:`R`.

**The law this yields** (stronger than what B3.0 had), which is B3.0's own sentence
*"R = I exactly when the BC is a pure symmetry statement adding no physics"* finally
used as a discriminator:

| law | `G` | `R` | what it asserts |
|---|---|---|---|
| `ReflectiveBoundary(axis)` | `SpecularMirror` | `I` | symmetry plane — a **quotient**, zero physics |
| `PeriodicBoundary(axis)` | `SpatialWrap` | `I` | torus — a quotient |
| `AlbedoBoundary(α, specular)` | `IdentityMap` | `SpecularReemission(α)` | a **surface** returning α specularly |
| `AlbedoBoundary(α, isotropic)` | `IdentityMap` | `LambertianReemission(α)` | a surface returning α diffusely |
| `VacuumInflow` | `IdentityMap` | `0` | a surface returning nothing |

**⇒ EXACTLY ONE of `G`, `R` is non-trivial.** That is an illegal-states-unrepresentable
invariant, and it wants a gate.

**Consequence, flagged and NOT acted on in B3.4:** `ReflectiveBoundary(axis, albedo<1)`
is then **incoherent** — a symmetry plane cannot absorb. That object is
`AlbedoBoundary(α, specular)` wearing the geometry costume. It is unreachable from a
tag (`_law_from_tag` hard-codes `albedo=1.0`), so nothing production-facing depends on
it. Retiring the `albedo` parameter from `ReflectiveBoundary` is a **B5** item.

`SpecularReturn(axis)` and `SpecularMirror(axis)` are structurally identical and
semantically disjoint — which is the POINT, not a smell: two types make "put a
surface law in the geometry slot" unspellable.

### 11.2 Periodic — build the channel now; the quotient resolves at run time

> "Are we able to implement this now and use it to resolve option 3 at run-time? First
> of all because periodic for monte-carlo will probably need this anyway, but maybe
> there is some SCC treatment we can do in SN to use this to resolve it in a different
> way?" — user, 2026-08-01

**Ruling: build the partner-face channel (option A) now.** The quotient reading
(option 3) is then *asserted at realization* rather than baked into the mesh topology:
the identification `Γ₊(partner) ≡ Γ₋(face)` becomes a guard, not a restructure. MC will
need the same face-partner map.

`[M]` **The identification holds on every quadrature in the tree** — measured
2026-08-01, `Γ₊(partner)` and `Γ₋(face)` are EQUAL as sorted index arrays for
`gauss_legendre(8)`, `product(2,4)`, `level_symmetric(6)`, `lebedev(17)`, on both axis
pairs. So the realized periodic operator is the **identity** between two different
faces' restrictions, and the set equality is a genuine geometric invariant
(`n_f = −n_f'`) — not a shape check. (`product(2,4)` carries **4 tangential ordinates
per face**, which is also why white's strict `> 0.0` mask is wrong there.)

### 11.3 The SCC answer — the criterion already exists, unwired

`orpheus/derivations/discrete/sn/sweep_acyclicity.py` already ships
`TraceDigraph.strongly_connected_components()`, `is_acyclic`, `cyclic_components()`,
and models periodic's trace edge correctly (`outflow(f,n) → inflow(partner,n)`). Its
own docstring states the position:

> "The honest criterion is therefore not a boolean on the boundary *kind* but a
> strongly-connected-component decomposition of the trace digraph" … "nothing in the
> solver ever builds a digraph and nothing can detect a cycle."

So the SCC treatment is the **recorded algebra of record with no production consumer**.
Wiring it would replace the `permutes_ordinates` heuristic in `_reflective_faces` — the
thing that decides `B_lower`/`B_upper` — with the honest criterion, and would let
configurations that are *acyclic* (one white face, one periodic face + vacuum) sweep in
ONE pass instead of being needlessly lagged. **Its own phase; filed as an issue, NOT
built in B3.4.** B3.4 lags periodic into `B_upper`, which is what today's split already
does and is correct-but-not-minimal.

### 11.4 B3.4 split into three commits

| step | rows | design risk |
|---|---|---|
| **B3.4a** | white, prescribed_inflow | none — mechanical narrowing; dissolves the `> 0.0` twin AND the inflow mask |
| **B3.4b** | albedo ×3 | the §11.1 ruling |
| **B3.4c** | periodic | the §11.2 ruling |

Each commit carries its own doc repairs (user standing directive: a falsified doc is a
bug, fixed in the SAME change).

### 11.5 Doc falsehoods B3.4 owes (all `[M]` measured)

1. `numerics/operator.py:2622-2624` — `PeriodicWrapOperator` claims "the SN sweep
   handles the spatial wrap via its own face-pair indexing". **No such mechanism
   exists anywhere** (searched: sweep schedule, sweep graph, trace space, solver).
2. `geometry/boundary/periodic.py:42-44` — claims "the two-face plumbing is handled by
   whoever instantiates `PeriodicBoundary` and orchestrates the sweep". Nobody does.
   Mutually inconsistent with (1).
3. `numerics/operator.py:2636-2639` — xrefs `PeriodicBoundary.apply`, which **does not
   exist** (descriptors are not callable). A dangling Python-domain xref, the silent
   class `-W` never catches. Repeated at `tests/numerics/test_periodic_wrap_operator.py:12,68`.
4. `sn/operators/boundary.py:220-221` — claims periodic advertises `apply_transpose`;
   `SpatialWrap.is_adjointable` is `False`. The composite predicate reads the REALIZED
   operator (identity body ⇒ `True`), the law factor says `False` — two sources of
   truth disagreeing.
5. `sn/boundary/angular.py:55-65,71` — a **B1-era leftover** my B3.6 sweep missed: the
   note still says `response_kernel` carries "the crossing Γ₊ → Γ₋", that "that factor
   is the scalar α", and that ":math:`G_{\text{diff}}` is the geometry". All three are
   false post-B3.0 — the average IS the response.

---

## 12. B3.4b design — albedo's re-emission closure

Executing §11.1. **`AlbedoBoundary.geometry_map` is `IdentityMap()` unconditionally**; the
closure's content lives entirely in `R`.

### 12.1 The shape

```python
# _factors.py — a THIRD response kernel
@dataclass(frozen=True, slots=True)
class SpecularReemission:
    r"""R = α · (the mirror pairing) — a SURFACE returning specularly."""
    alpha: float = 1.0
    axis: str = "x"
    amplitude -> α ;  is_zero -> α == 0 ;  is_adjointable -> True (a permutation)

# _factors.py — the closure tier: amplitude-FREE shapes the law instantiates
class ReemissionClosure(Protocol):
    def kernel(self, alpha: float) -> BoundaryResponseKernel: ...

@dataclass(frozen=True, slots=True)
class SpecularReturn:
    axis: str = "x"
    def kernel(self, alpha): return SpecularReemission(alpha, self.axis)

@dataclass(frozen=True, slots=True)
class IsotropicReturn:
    axis: str = "x"; outward_sign: int = +1
    def kernel(self, alpha): return LambertianReemission(alpha, self.axis, self.outward_sign)

# albedo.py
albedo: float = 0.0
reemission: Optional[ReemissionClosure] = None
geometry_map    -> IdentityMap()                      # ALWAYS
response_kernel -> ScalarResponse(albedo) if reemission is None
                   else reemission.kernel(albedo)
```

**Why the closure is amplitude-free.** `amplitude` is a `BoundaryResponseKernel` Protocol
member (the diffusion realizer reads `law.response_kernel.amplitude`), so the kernel must
carry α. If the closure carried it too there would be two sources of one number. The
closure is therefore a *shape* the law instantiates with its own `albedo` — one α, on the
law, where the tag parser already puts it.

**Why not "put the kernel on the law directly."** That would force removing the `albedo`
field (α would live on the kernel), breaking `AlbedoBoundary(albedo=…)` at `_law_from_tag`,
the diffusion registry, and ~14 test files — for no gain, since the closure tier keeps α
single-sourced anyway. Blast radius decided this, not taste.

### 12.2 The realization — a net REDUCTION in duplication

The three albedo rows must NOT grow a second copy of the specular construction and a
second copy of the Lambertian. Extract first, then dispatch:

| realizer helper | consumed by |
|---|---|
| `_specular_operator(quad, method_space, axis, alpha, law_key)` | `ReflectiveBoundary` **and** `AlbedoBoundary + SpecularReturn` |
| `_checked_angular_average(...)` (exists, B3.4a) | `WhiteBoundary` **and** `AlbedoBoundary + IsotropicReturn` |

⇒ ONE realization body per kernel kind, reached by two law NAMES. The `≡` theorems
(`albedo(α, specular) ≡ reflective(axis, α)`, `albedo(α, isotropic) ≡ white(α)`) then hold
because the two routes literally execute the same code, not because two transcriptions
agree — and they get pinned as tests anyway, since a shared body can still be reached with
different arguments.

**This is also the B4 groundwork.** Once the bodies key on the KERNEL rather than the LAW,
B4's "realizer reads the factors instead of `isinstance`-dispatching" is a dispatch swap on
`type(law.response_kernel)`, not a rewrite.

### 12.3 `reemission=None` — the angular-resolution refusal

`ScalarResponse(α)` is **complete on a scalar trace** (one dof: `J⁻ = α J⁺`) and
**incomplete on an angular one** (`α·I` is a `Γ₊ → Γ₊` endomorphism; `IdentityMap`
supplies no crossing). So:

* **diffusion** realizes `AlbedoBoundary(α)` exactly as today — reads `amplitude`, done.
  `BC("albedo", albedo=…)` is UNCHANGED for the one method that can reach it from a tag.
* **SN** refuses it, naming the two completions. This is the first bite of the
  method-realizability taxonomy's **angular-resolution** axis (the other two — state-cone
  and spatial-topological — already bite at zero-flux and diffusion-periodic).

That is *not* the blanket refusal the user overruled: the complete law is fully built and
realizable, and only the under-determined spelling is refused.

### 12.4 Consequence to flag, NOT act on (→ B5)

`ReflectiveBoundary(axis, albedo<1)` is incoherent under §11.1's law — a symmetry plane
cannot absorb; that object is `AlbedoBoundary(α, SpecularReturn(axis))`. Unreachable from a
tag (`_law_from_tag` hard-codes `albedo=1.0`). Retiring reflective's `albedo` parameter is
**B5**, together with typing the Lambertian as the rank-one `u ⊗ v` it now visibly is.

---

## 13. ⚠ B3.4a — my equivalence measurement was WRONG, and the correction matters

The archivist's doc-repair pass refused to publish my measured-evidence block and
re-derived it. It was right; I have since reproduced the correction independently. Recording
both the fact and the mechanism, because the mechanism is the reusable part.

### 13.1 What I claimed, and why the probe could not see the truth

I claimed white was **bit-identical on `gauss_legendre(8)` and `product(2,4)`**, 1 ULP on
`lebedev(17)` / `level_symmetric(6)`. The `product` half is false, and my probe was
constructed so that it *could not* be anything else:

```python
full = np.zeros((quad.N, 3)); full[sp.outflow_indices] = x   # ← the defect
```

The pre-B3.4a operator consumed the **whole face slot**, whose tangential rows carry real
flux. My reference scattered the probe onto :math:`\Gamma_+` **only**, leaving every
mis-admitted row at ZERO — i.e. it nulled exactly the term under test. That is vv **Mode 7**
(an ansatz that cancels the math it claims to exercise) committed inside a *verification
probe* rather than an MMS. The general form: **when comparing an old path to a new one, the
probe must be built in the OLD path's domain, not the new one's** — building it in the new
domain silently assumes the very restriction being introduced.

### 13.2 The corrected measurement `[M]`

Which ordinates the retired `> 0.0` classifier admits that :math:`\Gamma_+` does not:

| quadrature | tangential / face | MIS-ADMITTED by `> 0.0` |
|---|---|---|
| `gauss_legendre(8)` | 0 | 0 |
| `product(2,4)` | 4 | 2 on `xmin`/`xmax`/`ymax`, **0 on `ymin`** |
| `product(3,4)` | 6 | 3 on `xmin`/`xmax`/`ymax`, **0 on `ymin`** |
| `level_symmetric(6)` | **0** | 0 |
| `lebedev(17)` | **12** | **0** |

So the honest statement is **the `product` family only, and there only on
`xmin`/`xmax`/`ymax`** — NOT "every production quadrature but `gauss_legendre`", which is
what `angular.py` said and what I wrote into the B3.4a briefs. `lebedev` carries the most
tangential ordinates of any production quadrature and mis-admits none of them; the
`ymin`-vs-`ymax` asymmetry is FP-sign noise about exact zero, so the defect is not even
stable under a face flip.

### 13.3 The part that changes the ARGUMENT, not just a number

On `product(2,4)` `xmax` the mis-admitted rows carry `cos_w ≈ 7.85e-17` against a norm of
`2.5650996603237282`. Adding them to that norm is below its ULP, so

> **Δnorm is EXACTLY 0.0** — the whole discrepancy sits ψ-weighted in the NUMERATOR.

which means it scales with the flux ratio between the tangential rows and :math:`\Gamma_+`,
and is **unbounded by floating point**:

| tangential : Γ₊ flux ratio | \|old − new\| | relative |
|---|---|---|
| `1e0`  | 0.0 | 0.0 |
| `1e3`  | 5.9e-14 | 5.7e-14 |
| `1e6`  | 6.0e-11 | 5.4e-11 |
| `1e12` | **6.2e-05** | **7.5e-05** |

The `1e12` figure reproduces the `6.1e-05` already recorded in `LambertianReemission`'s
docstring from a cylinder under `product(2,4)` — independently, from a different direction.

⇒ **TWO mechanically different effects both measure ≤ 1 ULP on an O(1) probe:**

1. **reduction-order** (`lebedev`, `level_symmetric` — they mis-admit NOTHING): the sum now
   runs over a restricted array instead of a zero-padded one. Genuinely FP-level and
   BOUNDED. Principled-equivalence.
2. **mis-admission** (`product` family): FP-level *only while the flux is O(1) on the
   tangential rows*, and unbounded otherwise.

**A ULP-magnitude table therefore CANNOT justify this change** — it would report ≈1e-16 for
both and read as "FP-neutral everywhere", which is exactly wrong on the one quadrature that
motivated the phase. The justification is **structural**: one classifier, not two. Keep the
magnitudes as evidence of what the second classifier *was doing*, never as the argument.

### 13.4 Repairs this correction owes

- `orpheus/sn/boundary/angular.py` — the class docstring's "That is every production
  quadrature but `gauss_legendre`" is FALSE; replace with §13.2's statement. Same line
  carries a stray `` `[M]` `` editorial marker.
- `orpheus/geometry/boundary/__init__.py` ~216-250 — **the highest-value find**: white is
  still described as `G = G_diff`, the cosine-weighted Lambertian average. That is the
  EXACT misassignment B3.0 corrected, still live in the package docstring two phases later.
  Prescribed still says `R = G = 0` and describes the retired mask.
- `orpheus/geometry/boundary/_base.py` ~334-357 —
  `assert_source_lives_on_incoming_trace` still says the guarantee "rests on the realizer
  masking the evaluation".
- `tests/geometry/test_bc_equivalence_snapshot.py` ~101-108 + its twin in `test_boundary.py`
  — `_MIXED_LAW_XFAIL`'s reason text says white is un-narrowed. The row still xfails, for a
  DIFFERENT reason than documented: vv Mode 8 class 4, a misattributed strict xfail.

---

## 14. ⭐ #325 — the ROOT cause behind B3.4a's classifier twin (user ruling 2026-08-01)

**Ruling: file it and DO it, sequenced BEFORE B4** (after B3.4b/c). The user also
reframed it correctly — *"this is probably an artifact of the equispaced quadrature in
general as implemented … it would be good to generate the nodes programmatically and
automatically for any of these quadratures (and any other family with similar problem)."*
Measurement confirmed that reframing on both counts.

### 14.1 The principle

> A node set **defined by a symmetry group** must be generated by that group's ACTION —
> sign flips, coordinate swaps, index arithmetic, all exact on floats — not by
> **evaluating a parameterization** of it. Evaluating `cos`/`sin` at symmetric angles
> destroys the symmetry to rounding, and no downstream tolerance recovers it.

`[M]` **Two families already obey it and two violate it, and the split is exactly
algebraic-vs-trig generation** — which is not a coincidence, it IS the principle:

| family | generation | exact mirror | exact tangential zeros |
|---|---|---|---|
| `lebedev` | algebraic (`O_h` orbits of ±(a,b,c)) | ✅ x,y,z | ✅ 12/axis, all exactly `0.0` |
| `level_symmetric` | algebraic (sign/permutation orbits) | ✅ x,y,z | ✅ (carries none) |
| `gauss_legendre` | `leggauss`, no azimuth | ✅ x | ✅ (carries none) |
| **`product`** | **`np.cos(2πk/n)`** | ❌ x, ❌ y | ❌ **zero** exact on `mu_x` |
| **MoC azimuthal** | **`np.cos(π(2k+1)/(2n))`** | ❌ for **EVERY** `n_azi` | ❌ when `n_azi` odd |

The MoC row is the finding the reframing bought: its `φ → π−φ` mirror — the symmetry that
pairs a track with its reflected partner — is broken for `n_azi ∈ {2,3,4,6,8,16}`, i.e.
every value measured, *including* those with no on-axis node. Looking only at `product`
would have missed it.

### 14.2 The primitive — ONE home, THREE consumers

Every circular node set in the tree sits at an exact **rational** multiple of 2π, so take
the rational directly: `circle_nodes(numerators, denominator) -> (cos, sin)` of
`2π·p/q`. Evaluate trig on the first OCTANT only; generate the other seven by sign flips
and a coordinate swap. Two fixed points need exact values, **both decided by integer
arithmetic, neither by a tolerance**: the axis (`4p ≡ 0 mod q` → `0`/`±1`) and the 45°
diagonal (`2r == q` → both components `√2/2`, since the diagonal is a fixed point of the
swap and `np.cos(π/4)`/`np.sin(π/4)` differ by 1 ULP).

| consumer | angle | call |
|---|---|---|
| `numerics/quadrature/rules_product.py:115` | `2πk/n` | `circle_nodes(k, n)` |
| `moc/quadrature.py:89` | `π(2k+1)/(2n)` | `circle_nodes(2k+1, 4n)` |
| `derivations/discrete/sn/balance.py:367` | `2πk/n` | `circle_nodes(k, n)` |

The home is a shared `orpheus/numerics/` module, **not** the quadrature package — MoC is a
consumer and MoC is not a quadrature family.

`[M]` Prototype verified: product on-axis / x-mirror / y-mirror / 45°-diagonal all exact
for `n_phi ∈ {3,4,5,6,8,10,12,16,20,24,32}`; MoC mirror exact for every `n_azi` measured
(was false for every one); unit norm unchanged at ≤ 1 ULP; values move ≤ 8.3e-16.

### 14.3 What it retires

- **`TANGENTIAL_EPS` stops deciding set membership.** It currently gates the three-way face
  partition — an epsilon deciding *which ordinates are physically grazing*, i.e. a magic
  number deciding physics rather than a convergence tolerance. It exists only because
  `product` cannot represent its own tangential set. Post-#325 it demotes to a
  provably-inert defensive guard, gated by a test that every production quadrature has
  EXACTLY zero tangential cosines.
- **`reflection_index`'s nearest-neighbour search.** Its docstring concedes the partners are
  matched *"via nearest-neighbour matching against the reflected node positions"* — which is
  what you need precisely when the reflected node is not bit-equal to any real node. With
  exact symmetry it becomes an exact lookup, and the ERR-042/044/045 involution becomes
  exact rather than 1e-16-approximate.
- **The bug CLASS B3.4a removed one instance of.** Any code classifying by the sign of a
  direction cosine is wrong on `product`; today the only defence is remembering
  `TANGENTIAL_EPS`. #325 makes the mistake unspellable instead of guarded.

### 14.4 Sequencing

**B3.4b → B3.4c → #325 → B4.** #325 before B4 because B4 composes `R ∘ G` under a guard
that actually checks spaces, and exact quadrature symmetry is worth having first. After
B3.4b/c because those are mid-campaign boundary-algebra work and their new gates are built
on `gauss_legendre` (no azimuth ⇒ untouched by #325), so the double-re-baseline risk is
negligible.

---

## 15. ⭐ B3.4b AS EXECUTED — two defects the shared-body extraction exposed

§12's design landed with **no change to its shape** — `geometry_map` is `IdentityMap()`
unconditionally, the closure tier is amplitude-free, the two `≡` theorems hold because the
routes execute one body. What §12 did *not* anticipate is that folding four call sites into
two bodies would surface two live defects. Both are fixed in the same commit; both are the
"clean before extending" rule paying out.

### 15.1 ⚠ The specular pairing's certification was welded to the WRONG law

`ReflectiveBoundary` certifies three **independent** invariants of
`quadrature.reflection_index(axis)` — measure preservation (ERR-042), involution
(ERR-044), inflow→outflow (ERR-045) — as its own methods. That was correct while reflective
was the only law standing on that table.

`AlbedoBoundary(α, SpecularReturn(axis))` realizes through **the same table**. Shipping §12
as written would have meant: a wrong `reflection_index` table is **caught on the reflective
route and silently realized on the albedo route**. That is exactly the twin-path failure
this campaign exists to remove, and it would have been introduced BY the phase correcting
the G/R conflation.

`[M]` The three are genuinely independent — each passes on a table the other two reject
(a cross-weight-class pairing survives involution + sign; a self-map survives involution +
measure; a 3-cycle survives measure + sign). So no subset substitutes.

**Fix:** the invariants moved to the *pairing*, `orpheus/geometry/boundary/_specular.py`,
with `assert_specular_pairing_{measure_preserving,involutive,maps_inflow_to_outflow}` plus
the aggregate `assert_specular_pairing_valid`. Reflective's three methods survive as
one-line delegations (the measure one MUST stay a method — it is the polymorphic hook the
base template fires among the universal five); albedo's `assert_realizable` fires the
aggregate when its closure is specular.

**Tier note, and it matters:** albedo does NOT override
`assert_geometry_map_measure_preserving`. Its `G` *is* the identity and *does* preserve the
measure, so that hook is honestly a no-op there; the pairing sits in `R`, so the check is a
law-specific extension. Reading the geometry hook as "the pairing's check" would re-commit
the very conflation the user's ruling forbids. The hook's NAME is now narrower than the
concept — noted, not renamed: it is not false, and renaming a polymorphic hook with live
overrides belongs to B4 where the realizer starts dispatching on factors.

### 15.2 ⚠ α = 0 was UNREACHABLE on the geometry-tier laws — a live crash

`ScaledOperator` refuses a zero scalar ("degenerate; use ZeroOperator explicitly"). The
pre-B3.4b reflective and white arms both ended `return float(law.albedo) * base` with only
an `α == 1.0` fast path. So `ReflectiveBoundary(axis, 0.0)` and `WhiteBoundary(..., 0.0)` —
**legal laws**, α = 0 satisfies every invariant including sub-Markov — died in the numerics
layer with a message about operator degeneracy instead of realizing the boundary they
describe. Neither is tag-reachable, so no production path; it is a crash on a legal value
that two arms carried independently and neither test suite probed.

Only visible once the four routes shared one body. **Fix:** `_attenuated_kernel_operator`
answers α = 0 with `_narrowed_zero_operator` — because a surface that returns nothing IS a
vacuum, and now says so with the same object, honest space hooks and a working transpose.
Vacuum's own arm was repointed at that helper too, so the narrowed zero map has one
construction site instead of two.

### 15.3 `[M]` Measured — the two `≡` theorems

Bit-identical (`np.array_equal`, no tolerance) across **3 quadratures × 3 amplitudes**:

| | α = 1.0 | α = 0.7 | α = 0.0 |
|---|---|---|---|
| `albedo+specular ≡ reflective` | ✅ | ✅ | ✅ |
| `albedo+isotropic ≡ white` | ✅ | ✅ | ✅ |

on `gauss_legendre(8)` (\|Γ₊\|=4), `level_symmetric(6)` (24), `lebedev(17)` (49). Bit-identity
is the right gate here and needs no justification by measurement: the two routes execute the
SAME construction, so any difference would be an argument-threading bug, not arithmetic.

### 15.4 `[M]` Blast radius — 9 tests, ZERO production

`AlbedoBoundary` is absent from `SNMesh.BOUNDARY_OPERATOR_REGISTRY`, so no `BC(...)` tag
reaches the SN realizer and the §12.3 refusal costs nothing production-facing. Diffusion's
registry keeps `"albedo"` and its whole suite is green.

Measured full-suite: `tests/geometry` + `tests/diffusion` = 452 passed / **3 failed**;
`tests/sn -m "not slow"` = 2488 passed / **6 failed**. All 9 are bare-`AlbedoBoundary`
realizations through the SN realizer:

| file | what | disposition |
|---|---|---|
| `test_sn_boundary_realizer.py::TestRealizeAlbedo` ×3 | pins Zero/Identity/Scaled-TP fast paths | the fast paths are gone — re-pose |
| `test_operator_block_role.py` `_LINEAR_LAWS` albedo ×3 | block-role over every law | migrate to completed albedos |
| `test_bc_equivalence_snapshot.py::TestAlbedo05Lebedev17Snapshot` | frozen full-face `0.5·ψ` | see §15.5 |
| `test_boundary.py::test_albedo_bc_scales_outgoing` | α scaling on the full face | re-pose |
| `test_boundary.py::test_albedo_zero_and_vacuum_agree_on_inflow_rows` | α=0 ≡ vacuum | **keep the claim** — §15.2 makes it literally true now |

### 15.5 ⚠ The three strict-xfails now xfail for the WRONG reason — created BY this phase

The `albedo_*` rows in `test_b3_domain_narrowing.py` carry `xfail(strict=True)` documenting
*"realizes FULL-FACE … echoes a Γ₊ input back"*. Post-B3.4b they still xfail — but on a
`BoundaryError` at **realize**, never reaching the assertion they document. That is
`vv-principles` Mode-8 **class 4** (the misattributed strict-xfail), manufactured by this
phase in the very file that tracks it.

They must be **migrated, not flipped**: the row's law object becomes a completed albedo
(which genuinely narrows, so the marker deletes), and whether a *bare* albedo still belongs
in that gate at all is a separate question — its honest home is the refusal negative, not a
domain row. `--runxfail` is the check that would have caught this and is mandatory before
the commit.

### 15.6 Retirement note — an axis-index twin NOT created

`_AXIS_INDEX = {"x": 0, "y": 1, "z": 2}` was a `reflective.py` module constant, and
`orpheus/sn/solver.py:1463` spells the same literal inline. The moved code uses
`AXIS_NAMES.index(axis)` — the canonical single source — so B3.4b adds no third twin. The
`solver.py` one survives; deduping it is a one-liner in a hot file with zero correctness
gain, deliberately not bundled here.

### 15.7 Dispatch reads the FACTOR, not the closure field

First draft of the albedo arm dispatched on `isinstance(law.reemission, SpecularReturn)`.
Corrected before commit: it dispatches on `law.response_kernel`
(`SpecularReemission` / `LambertianReemission` / `ScalarResponse`).

`reemission` is the field the law happens to store; `response_kernel` is the **affine
form's tier**, which is the public surface B1/B3.0 built and the thing every realizer is
supposed to consume. Reaching past it into a private field would have bypassed the
abstraction this campaign exists to establish, in the phase that extends it. Three
consequences, all good:

- α is single-sourced at the point of use (`kernel.amplitude`, not `law.albedo`).
- The certification in `assert_realizable` keys the same way, so a future closure producing
  a `SpecularReemission` inherits the table checks automatically.
- It **is** §12.2's B4 shape, already: B4 generalizes this branch's chain into the
  realizer's ONE dispatch across all laws, rather than replacing it.

`[M]` Re-verified after the switch: both `≡` theorems still bit-identical on all
3 quadratures × 3 amplitudes; diffusion still `ScaledOperator(0.4)` for all three
spellings; pyright back to the #288 floor of 1.

### 15.8 ⚠ A THIRD defect — the `≡` broke at the curvilinear ray corner

`_has_ruled_corner_action` (`orpheus/sn/operators/boundary.py`) decides whether a law's
inflow at the off-quadrature `μ = ±1` ray can be written down. It read

```python
law.geometry_map.permutes_ordinates or law.response_kernel.is_zero
```

— complete while the only specular pairing lived in `G`. The user's ruling puts one in `R`,
so `AlbedoBoundary(α, SpecularReturn(a))` has `G = IdentityMap` and was **loud-deferred at
the corner while `ReflectiveBoundary(a, α)` — the same matrix — is ruled.**

`[M]` Measured before the fix: `reflective(x, 0.7)` → `True`, `albedo(0.7, Specular(x))` →
`False`. **The `≡` theorem was false in the one consumer that reads the FACTORS rather than
the realized operator.** The realizer equivalence being bit-identical is exactly what would
have hidden this: a gate on realized output cannot see a factor-reading consumer.

**Fix:** ask the pairing of BOTH tiers,
`isinstance(law.response_kernel, SpecularReemission)` added to the disjunction. `[M]` After:
`albedo+specular` = `reflective` = True, `albedo+isotropic` = `white` = False, bare albedo
False, across all ten law spellings.

#### The tidier fix is exactly wrong — and the reason is a standing test

The obvious move is `permutes_ordinates` as a `BoundaryResponseKernel` Protocol member, so
the predicate reads both tiers uniformly with no `isinstance`. **Do not.**
`SpecularReemission` already carries `is_adjointable`; adding `permutes_ordinates` gives it
the complete member set of `BoundaryGeometryMap`, so it would satisfy that Protocol
structurally. `tests/geometry/test_boundary_factors.py:180-186` asserts the two tiers are
**disjoint** — precisely to stop a response from posing as a geometry, which is the
conflation B3.0 corrected. A convenience member is not worth disarming that guard.

**The generalizable lesson:** when a concept legitimately spans two tiers that a test keeps
structurally disjoint, it must be asked of each tier **in that tier's own vocabulary**.
Hoisting it to a shared Protocol member is how disjoint tiers silently merge.

#### Where else does a factor-reading consumer exist?

Audited: `_has_ruled_corner_action` was the one with a tier assumption. The composite
`SNBoundaryOperator.is_adjointable` reads each REALIZED law's own `is_adjointable`, so it
was already correct and only its prose enumeration ("reflective / vacuum / periodic /
albedo are adjointable") was stale — albedo now answers by its closure, and the code always
would have. Repaired in the same commit.

### 15.9 The shared predicate — user ruling 2026-08-01, and the fourth defect

§15.8 found ONE tier assumption. A systematic grep of every `.geometry_map` /
`.response_kernel` consumer found **four**, all asking the same question of the `G` tier
alone:

| site | what it decides |
|---|---|
| `sn/loss_representation/sweep_schedule.py:280` | the reflecting-face set ⟹ the `B_lower`/`B_upper` schedule split (the LAG) |
| `sn/acceleration/dsa.py:234` | DSA low-order admission |
| `sn/operators/boundary.py` `_has_ruled_corner_action` | is the off-quadrature `μ = ±1` ray expressible |
| `sn/operators/boundary.py` `_reflect_corner` | the swap itself |

All four read `sn_mesh.bc[face].law`, and `AlbedoBoundary` is absent from
`SNMesh.BOUNDARY_OPERATOR_REGISTRY` — so **none is reachable today**. Every one fires the
day #189 registers the law.

**User ruled: extract the shared predicate now** (over "fix only the corner one" and "defer
all four to B4"). Landed as `law_permutes_ordinates(law)` in
`orpheus/geometry/boundary/_base.py`, exported from the package, called from all four
sites. `[M]` after: `albedo+specular` = `reflective` = True; `albedo+isotropic` = `white` =
False; vacuum/periodic/bare-albedo unchanged.

**Why a function and not a property or a Protocol member.** A property on
`BoundaryTraceLaw` can be shadowed by a subclass, which reintroduces exactly the drift
being removed; a module-level function has one body that cannot be overridden. And the
Protocol route is affirmatively wrong — see §15.8.

**The lesson is about latency, not about the predicate.** A known-false predicate sitting
behind a registry gate is a planted landmine: the phase that makes it false is the only
phase that will ever have the context to notice, and #189 will arrive with no reason to
look here. "Not reachable" is a statement about today's registry, not about correctness.

### 15.10 Two message-contract repairs (test-architect findings, verified)

1. **The specular arm leaked an unattributed `ValueError`.** `SpecularReturn("y")` installed
   on `xmax` died inside `TraceRestrictionOperator.to_local` with a raw index complaint,
   while the diffuse sibling raised an attributed `BoundaryError` naming the mismatch. The
   asymmetry pre-dates B3.4b for `ReflectiveBoundary(axis="y")` — but nothing could *declare*
   an axis that disagreed with its installation face until the closure gave it one, so it was
   unreachable in practice. Fixed by catching and re-raising with the semantic diagnosis;
   `to_local` stays the SINGLE index authority (no second membership test to drift from it).
2. **The α = 0 refusal cited the wrong reason.** At α = 0 nothing returns, so no pairing is
   needed and the law is NOT under-determined — every closure would agree. It is refused for
   a different reason: it is a `VacuumInflow` **twin**, and admitting it would give one
   physical law two spellings with two realization paths. Refusing also keeps SN's albedo
   admission one uniform rule rather than one turning on an exact float compare (realizing at
   α = 0 while refusing at α = 1e-300 is a worse contract than refusing both). Both α = 0 and
   α ≠ 0 now carry their own accurate message. `AlbedoBoundary(0.0, SpecularReturn(a))`
   still realizes — to the narrowed zero map.

### 15.11 ⭐ `[M]` The positional pairing was a CONFIGURATION-DEPENDENT ACCIDENT

The sharpest correction of the phase, from the test-architect, reproduced independently.
Comparing the old bare-albedo positional pairing (`inflow[j] ← α·outflow[j]`) against the
specular one, on `xmax`:

| quadrature | positional == specular |
|---|---|
| `product(2, 4)` | **True** |
| `level_symmetric(6)` | **True** |
| `gauss_legendre(4)` / `(8)` | False (the slab mirror reverses order) |
| `lebedev(17)` | False |

So pre-B3.4b a bare albedo behaved **exactly like a mirror on two of the tree's
quadratures** and like nothing in particular on the others, silently. A user who validated
on `level_symmetric` and ran on `lebedev` got different physics from the same law object.

That upgrades the refusal's justification from "the pairing carries no geometry" to "the
old answer was a coincidence of index order that *looked* defensible on half the fixture
set". Written into `albedo.py` and the SN boundary theory page. It also generalizes: **a
default that is right on some fixtures and wrong on others is more dangerous than one that
is uniformly wrong**, because the fixture set decides whether anyone notices.

### 15.12 ⚠ A shared body makes the `≡` gate BLIND to bugs in the shared body

§12.2's "the theorems hold because the two routes execute the same code" is true and is
exactly why the equivalence gate cannot be the phase's verification. `[M]` (test-architect):
under a `to_local → arange` mutation the route-equivalence suite is **60 passed / 0 failed**
while the operator is wrong on 3 of 5 quadratures — both routes are wrong *identically*.

A1 (route equivalence) is **necessary and insufficient**. The catchers are the
independent-expression anchors: the specular route checked against a hand-built mirror
expectation, the diffuse route against a hand-built cosine-weighted average. Those must not
be deleted later as redundant with A1, and the plan says so in place.

Corollary, and the reason §15.9's gate matters: the factor-reading consumers do NOT share a
body with the realizer, so an `≡` gate posed over THEM has teeth that A1 structurally
cannot have.

### 15.13 The snapshot generator — user ruling 2026-08-01

`tests/geometry/_generate_bc_equivalence_snapshots.py` cannot regenerate 7 of its 8 cases.
This is **inherited breakage from B3.2/B3.4a**, not B3.4b's (which takes it 6/8 → 7/8): the
frozen `.npz` carry a **FULL-FACE pre-B3.2 `psi_in`**, which the current realizer cannot
produce, and the consuming tests deliberately restrict that full-face array to `Γ₋` — which
is exactly what made the artefacts an independent witness that the narrowing moved no values.

Three options were put to the user. **Ruled: MIGRATE the schema and re-anchor each case
against a structurally-independent expression** (over "retire the generator and freeze the
artefacts" and "reconstruct the retired full-face path inside the generator"). All 7
artefacts are invalidated by design.

**Why this is the vv-correct answer and not merely the most thorough one.** A frozen
snapshot of a retired code path is **procedurally** independent, not structurally so
(`vv-principles` L11): it certifies "the new path agrees with the old path", which is
worth exactly as much as the old path was right. The narrowing's whole premise is that the
old path was reading the wrong half-trace — so the reference being preserved is a witness
from a code path this campaign judged defective. Re-anchoring to a hand-built mirror
gather / cosine-weighted average / zero map replaces "it still does what it did" with "it
does what the math says", which is a promotion in the evidence hierarchy, not a lateral
move.

What is genuinely lost: the historical bit-identity witness across the narrowing. That
claim has already been discharged — B3.2/B3.4a each proved it at landing time, and their
proofs are recorded. A permanent artefact re-asserting a landed one-time claim is
archaeology.

**Own commit, after B3.4b** — a schema migration across 7 cases is not a tail of the
closure carve. Tracked as task #21.
