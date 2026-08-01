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
