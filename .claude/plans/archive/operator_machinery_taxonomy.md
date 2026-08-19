> **SUPERSEDED by orpheus-operator-machinery-report-v2.md (Part VI) — archived 2026-08-19 (plans triage).**

# Operator Machinery Taxonomy — the inverse/adjoint/matrix design (#226)

> **Planning artifact.** Deep design of what the SN (and cross-method) operator
> machinery should *become*, so that **every object's method surface equals its
> mathematical identity** — no vestigial promises. Inputs: the user's north-star
> principle, `neutron_transport_grand_report_v3.md` (§1/§3.1/§3.3/§5.7/§6.3/§21/§22/§38),
> a cross-domain-attacker taxonomy verdict (`operator_inverse_taxonomy_frames.md`)
> and an explorer ground-truth map (branch `refactor/inverse-as-operator @ b52bdef`).
> Code is NOT written from this yet — this is the taxonomy we agree on first.

---

## 0. The principle (the acceptance test for every object)

A `LinearOperator` is **one morphism with one promise**: `apply`. Anything it
*is* in addition — its inverse, its adjoint, its matrix — is produced by a
**factory** that returns a **distinct object**, and that object promises only
what it can deliver. An inverse-by-definition exposes the inverse action as its
`apply`; it carries **no** forward it never uses. A self-adjoint operator has
`A.H is A`. The frozenset-of-capabilities was a *string label* faking this; the
carve replaces it with **objects whose type IS the capability**.

---

## 1. THE LOAD-BEARING FACT — two layers, never conflated

The entire confusion (and the `_GaussSeidelResolvent` smell) comes from mixing
two layers that have *opposite* "one-vs-many" rules:

### Layer A — the SUBSTRATE (`LossRepresentation`)
`sn/loss_representation/__init__.py:229`. **ONE object, MANY action-views.** It
is the matrix-free per-cell coefficient kernel (WDD stencils, σ, angular
closures, the DAG/scan traversal). It is **not a morphism** — it is *data + a
traversal*. Its views:

| view | math | site |
|---|---|---|
| `loss_action(σ, ψ)` | `(L+C)ψ` (forward matvec) | :253 |
| `loss_action_transpose(σ, φ)` | `(L+C)ᵀφ` (adjoint matvec) | :294 |
| `sweep(q, σ, …)` | `(L+C)⁻¹q` (forward-substitution) | :241/:500 |
| `sweep_transpose(r, σ, …)` | `(L+C)⁻ᵀr` (adjoint solve) | **#280 — unbuilt** |

The `_OctantWalk` (`:717`) traverses ONE octant decomposition; forward and
inverse **fork only at the cell kernel + emit policy — "NEVER a boolean
`is_solve` flag"** (`:738`). So "one object, many views" is **CORRECT here** —
the substrate is exactly where it belongs. (L21 — "three actions of one object"
— is sound *for the substrate*, and was wrongly generalized to the operator.)

### Layer B — the OPERATOR (`LinearOperator`)
`numerics/operator.py:285`. **ONE morphism, ONE promise (`apply`).** `A`, `A⁻¹`,
`A.H` are **three different morphisms** that happen to *share the substrate*.
`.H` / `.inverse()` / `.as_matrix()` are **factories** returning **distinct
objects** (`_AdjointOperator`, a `SweepOperator`, an ndarray).

### The slogan correction (Q8 — your "different views" question)
> ❌ "Forward-substitution and sweep are different **views of one operator**."
> ✅ "Forward and inverse are **two operators sharing one substrate**."

The first phrasing is **flawed**, and it is *precisely* what licenses the
`_GaussSeidelResolvent` confusion (an object that owns both `apply` and `solve`).
The substrate has the views; **the operator pair does not.** Keep the word
"view" on the *substrate* and the word "operator" on the *morphism*, and the
whole design snaps into place.

---

## 2. The object graph

```
   ┌─────────────────────────────────────────────────────────────┐
   │  SUBSTRATE  (one object, many views — NOT a morphism)        │
   │  LossRepresentation:  loss_action · loss_action_transpose ·  │
   │                       sweep · sweep_transpose                │
   └───────────────▲───────────────────────────────▲─────────────┘
                   │ shared by the whole family     │
   ┌───────────────┴───────────────┐   ┌────────────┴─────────────┐
   │  FORWARD operator   A         │   │  these are DISTINCT       │
   │  A.apply(ψ) = (L+C)ψ          │   │  morphisms, produced by   │
   │                               │   │  A's factories:           │
   │  factories (return objects):  │   │                           │
   │   A.H        ─────────────────┼──►│  _AdjointOperator  (A†)    │
   │   A.inverse()─────────────────┼──►│  Sweep/Green/Dense/Resolv │
   │   A.as_matrix() ──────────────┼──►│  ndarray  (functor OUT)   │
   └───────────────────────────────┘   └───────────────────────────┘
```

Three of the four arrows are **endofunctors** (`Op → Op`): `.H`, `.inverse()`,
and composition. The fourth, `.as_matrix()`, is a **functor OUT of the category**
(`Op → Mat`, returns an ndarray) — the serialization boundary, dimension-gated.
That is why `as_matrix()` is **not** a co-equal "fourth view" (Q5 corrected): it
is a different *kind* of arrow.

---

## 3. `A.inverse()` is a STRUCTURE-KEYED FACTORY → four strategy objects

`A.inverse()` returns a `LinearOperator` whose `apply(q) = A⁻¹q`. **The concrete
subclass = the inversion ALGORITHM, chosen by A's STRUCTURE** (with an optional
`strategy=` override seam for the exceptions). All four are **sibling
`LinearOperator` subclasses** — the grand report §3.3 lists them as siblings; **no
abstract `InverseOperator` base is needed** (Q2/Q4 verdict).

| A's structure | `.inverse()` returns | algorithm | status |
|---|---|---|---|
| triangular `(L+C)` (lower-block-tri in DAG order) | **`SweepOperator`** | direct forward-substitution (LU frame) | **built** (`sweep_operator.py:44`) |
| general `(L+C−S)`, `(L+C−S−B)` (no triangular factor) | **`GreenOperator`** | Krylov / Neumann series of the `(L+C)`-preconditioned splitting | **deferred** (no literal consumer; it IS `SourceIteration`/`KrylovAcceleration`) |
| small / 0-D / **CP `[P]`** (dense by construction) | **`MatrixInverseOperator`** (§13 name; né `DenseInverseOperator`) | `lu_solve([A], ·)` | **exists un-named** (`direct_eigenvalue`, `homogeneous/solver.py:195`) |
| parametrized `(sT + A)` (α-eigenvalue, shift-invert) | **`A.resolvent(z)`** — a FACTORY name, no 4th class (§17 A3) | resolvent `(A−z)⁻¹` family; each instance = Sweep/Green/Matrix on the shift | future (α-work, §22) |

**Why structure-keyed, not free-strategy (Q4):** a CP `[P]` *cannot* sweep; a
meshed `(L+C)` *cannot* densify — the operator **owns** the selection. The
`strategy=` override is only for principled exceptions (force Krylov on a
triangular op for DSA; force Dense on a block).

---

## 4. The transformation chains — what is applied to ψ vs to q

This is the heart of your "show the transformations" request. The substrate's
coefficients are the same in every row; only the *traversal* and the *I/O role*
change.

```
FORWARD   (to ψ):   ψ ──A.apply──►  (L+C)ψ            loss_action(σ, ψ)          [matvec, forward traversal]
ADJOINT   (to φ):   φ ──A.H.apply─► (L+C)†φ           loss_action_transpose + G  [matvec, reverse coeffs + metric]

INVERSE — direct (to q):   q ──(L+C).inverse().apply──►  ψ : (L+C)ψ = q
                           sweep(q, σ): SAME coefficients, forward-SUBSTITUTION traversal (triangular solve)

INVERSE — Krylov (to q):   q ──A_loss.inverse().apply──►  ψ : (L+C−S)ψ = q
                           GreenOperator: GMRES iterating A_loss.apply(ψ_k) = (L+C−S)ψ_k,
                           PRECONDITIONED by (L+C).inverse() (the sweep). q is the rhs, ψ is the iterate.

INVERSE — dense (to q):    q ──A.inverse().apply──►  ψ : [A]ψ = q
                           DenseInverseOperator: lu_solve([A], q), [A] from A.as_matrix()

ADJOINT-INVERSE (to r):    r ──A.H.inverse().apply──►  the adjoint solve   (= #280, see §6)
```

**The clean reading:** you apply **forward operators** (`A`, `A.H`) **to fluxes
(ψ, φ)**; you apply **inverse operators** (`A.inverse()`, `A.H.inverse()`) **to
sources (q, r)**. The sweep and the matvec are the SAME substrate coefficients
walked in opposite directions — *that* is the true content of "two operators, one
substrate."

> **Note on the seed.** `SweepOperator.apply(q, *, initial_guess=ψ_prev)` carries
> an optional warm-start (the curvilinear Carlson coupled-pole seed for the
> direct sweep; the GMRES `x0` for Green). The result is `A⁻¹q` regardless of the
> seed — it is **application context (§38), not identity** — so it is a legitimate
> optional kwarg, not a leaked solve-parameter. (The cross-domain flagged the
> *appearance*; the §38 framing resolves it: context ≠ promise.)

---

## 5. The matrix-building path (the explicit loss representation)

`A.as_matrix()` materializes `[A]` by **applying A to the basis columns `eᵢ`** —
the pattern that already lives, un-promoted, as `_as_dense(op, ng)`
(`homogeneous/solver.py:109`). It is a **functor OUT** (`Op → ndarray`),
**dimension-gated** (refuse on a mesh too large to densify).

- Today only the **energy/scattering blocks** materialize (`dense_per_material`,
  `_as_dense` of `C − K_iso`). The **streaming `(L+C)` is never materialized** —
  it is matrix-free only. The future sparse-triangular path is already noted at
  `sweep_graph.py:66-71` (assemble per-octant lower-triangular `L_oct` as
  `scipy.sparse`, swap forward-substitution for `spsolve_triangular`).
- `as_matrix()` and the **`DenseInverseOperator`** are siblings of the same
  matrix view: `as_matrix()` *materializes*; `DenseInverseOperator` *inverts the
  materialization* (`lu_solve([A], ·)`). The **latent consumer that makes this
  finishable NOW** (not speculative) is `direct_eigenvalue`'s dense `A⁻¹F`
  (`homogeneous/solver.py:195`) — it already IS `A.inverse()` with the Dense
  strategy, written as free functions. **CP is the production method that earns
  the promotion** (its `[P]` is dense by construction).

---

## 6. The adjoints — a dagger category with partial inverses (Q7)

The right foreign frame is a **dagger category with restriction/partial
inverses**:
- `.H` is an **involutive endofunctor** (`A.H.H == A`), metric-aware:
  `A.H.apply(φ) = G_dom⁻¹ · (Aᵀ · (G_cod · φ))` — the Hilbert adjoint w.r.t. the
  inner products carried on the `FunctionSpace` (`G` block-diagonal on
  `FullFieldSpace`: `V` on bulk, `|Ω·n|·w` on the trace). Realized per-leaf via
  `apply_transpose` → `loss_action_transpose`.
- `.inverse()` is a **partial endofunctor** (defined only on isomorphisms);
  `is_invertible` is the **restriction idempotent** that says "defined here."
  The frozenset→predicate carve **is** that idempotent made explicit.
- **Coherence law:** `(A.H).inverse() == (A.inverse()).H` is **functoriality** —
  a *theorem*, and exactly the content of **#280** (`A.H.inverse()` = the
  transpose-solve).

**The asymmetry that explains the #280 deferral (Q2/Q6):**
- The **direct (sweep) inverse's `.H`** needs a **non-trivial reverse-DAG
  traversal** (`sweep_transpose` — reverse-order substitution on the transposed
  factor). The 2-D Cartesian `loss_action_transpose` itself *already raises*
  (`ScanMarch:1858`), so this is genuinely unbuilt. → **hard, deferred (#280).**
- The **Krylov inverse's `.H` is FREE**: `GreenOperator(A).H == GreenOperator(A.H)`
  — GMRES on the transpose-matvec. No new traversal.

So `SweepOperator.H` being deferred is **a property of which algorithm**, not an
oversight. The taxonomy makes that legible: the adjoint-inverse is cheap exactly
when the inverse is iterative, and expensive exactly when it is a direct sweep.

---

## 7. The resolvent fix — dissolving the "layer-confusion / vestigial-forward" smell

The cross-domain-attacker named a NEW smell and proved it with a test:

**Smell:** an OPERATOR-layer object (one morphism, one promise) carrying
SUBSTRATE-layer multiplexing — an `apply` AND a `solve` that are inverses of
**different** operators.

**Witness:** `_GaussSeidelResolvent` (`sn/solver.py:338`): duck-typed (not a
`LinearOperator`), `apply` = `(L+C)ψ` *documented "never called"*, `solve` =
`(L+C−B_lower)⁻¹q`. The vestigial forward exists **only** to fill a duck-typed
`CAP_APPLY` slot.

**The discriminating proof** — the round-trip `inv.apply(A.apply(x)) == x`:
- `_GaussSeidelResolvent` gives `(L+C−B)⁻¹(L+C)x ≠ x` (off by the B reflection) → **REDs.**
- A clean `(L+C−B).inverse()` (a `SweepOperator` on the boundary-folded operator) → **passes.**

**The fix (and what the resolvents BECOME — amended by the §17 re-evaluation):**
- `_GaussSeidelResolvent` → **`M.inverse()` with `M = (L+C−B_lower)` REIFIED as a
  real forward operator** (`M.apply = (L+C).apply − B_lower.apply`; `B_lower` = the
  existing boundary `LinearOperator` masked to the schedule's strictly-lower octant
  pairs — B has NO octant-diagonal (μ→−μ maps each octant elsewhere), so
  `B = B_lower + B_upper` is exact). A genuine `SweepOperator` on `M`, S-direct,
  **no vestigial forward**, round-trippable to machine precision against its OWN
  forward. The schedule stays a substrate concern; the iteration
  `ψ_{k+1} = M⁻¹(q + B_upper ψ_k)` stays in the driver. (§17 W2.)
- `_MomentWindowedResolvent` → **the typed composition `P @ A.inverse()`**
  (`OperatorProduct`, FullField→MomentField; `P` = the scattering frame's analysis
  face on the bulk ⊕ identity on the trace — a block coisometry, NOT invertible, so
  the composite makes no round-trip promise). The FUSION (never materialize ψ)
  stays a SUBSTRATE application-context: the base `SweepOperator` sweeps in
  `_SweepEmit` MOMENT mode. The `SweepOperator` itself stays ENDOMORPHIC and
  S-direct; the codomain change is owned by the explicit `P` factor, and the solver
  (`_maybe_window`) is the FACTORY that builds the product. NOT "a SweepOperator
  with a moment emit-config" — a substrate emit that silently changes the OPERATOR
  codomain would violate §1's own two-layer law. (§17 W1.)

So **neither resolvent is a method inside `SweepOperator`, nor a subclass of it**
(Q1/Q2): the G-S one is a `SweepOperator` *on a different (reified) forward
operator*; the windowed one is *a typed product with a `SweepOperator` factor*
(whose fused evaluation uses the substrate's moment emit). The "resolvent" concept
dissolves into "the inverse of the right operator, composed and fused explicitly."

---

## 8. Direct answers to the eight questions

1. **Are resolvents a method inside `SweepOperator`?** No. The G-S resolvent *is*
   a `SweepOperator` instance on the reified boundary-folded operator
   `M = (L+C−B_lower)`; the windowed one is the typed product `P @ A.inverse()`
   whose `SweepOperator` factor fuses via the substrate's moment emit (§17 W1/W2).
   They are not methods or members of `SweepOperator`.
2. **Is `SweepOperator` an abstract base the resolvents subclass?** No. `Sweep` /
   `Green` / `Dense` / `Resolvent(z)` are **sibling** `LinearOperator` subclasses
   (no `InverseOperator` ABC). The resolvents collapse into `SweepOperator` over
   the right operator / emit — not a new subclass each.
3. **Where does `GreenOperator` come from / what returns it?** `A.inverse()` when
   `A` is non-triangular (`A_loss = L+C−S`, or `L+C−S−B`). It is the **operator
   face of `SourceIteration` (Richardson) / `KrylovAcceleration` (GMRES)** and must
   *wrap* that driver, not re-implement it. Deferred until a literal
   `A_loss.inverse()` consumer exists (today the k-path is power-iteration with an
   SI inner).
4. **How is Krylov applied?** GMRES on the matrix-free `A_loss.apply` (`spla.gmres`),
   preconditioned by `(L+C).inverse()` (the sweep). `GreenOperator(A_loss).apply(q)`
   would be exactly this loop with an operator face.
5. **What happened with the loss representation? One per operator?** It is the
   **substrate, not an operator**. **One per `A`** (per streaming geometry),
   **shared by the whole family** — `A.apply`, `A.H`, `A.inverse()` all read the
   one `LossRepresentation`. (Today: `StreamingOperator` owns it; `InvertibleOperator`
   and `SweepOperator` delegate to that one instance.)
6. **Where is the explicit matrix instantiated?** `A.as_matrix()` (new) — the
   apply-to-basis pattern (`_as_dense`), a functor OUT, dimension-gated. The
   streaming `(L+C)` is matrix-free today; the sparse-triangular materialization is
   future (`sweep_graph.py:66`). The dense **inverse** of `[A]` is
   `DenseInverseOperator` (`lu_solve`), already un-named in `direct_eigenvalue`.
7. **Where are the adjoints?** `A.H = _AdjointOperator(A)` (metric `G⁻¹AᵀG`, the
   `G` on the `FunctionSpace`), per-leaf via `apply_transpose` /
   `loss_action_transpose`. The adjoint-*inverse* is **#280** — free for Krylov
   (`GreenOperator(A.H)`), hard for the sweep (reverse-DAG `sweep_transpose`,
   2-D transpose walk still unbuilt).
8. **Is "forward-substitution and sweep are different views of an operator"
   faithful or flawed?** **Flawed as stated.** Correct: forward and inverse are
   **two operators sharing one substrate**. The substrate carries the views; the
   operator pair does not. The flawed slogan is what produced the
   vestigial-forward resolvent.

---

## 9. What this means for the carve (migration implications — for the next session)

The taxonomy is **~85% already built and correct** (the substrate/operator split
exists; `InvertibleOperator` already routes apply/transpose/solve through one
`LossRepresentation` and `.inverse()` returns a `SweepOperator`). The remaining
work is to (a) make the OTHER `is_invertible` operators actually return an
inverse object, and (b) dissolve the resolvent confusion. Concretely:

- **`.inverse()` is ⅛ done.** Eight operators advertise `is_invertible=True`
  (`Multiplication`, `OperatorProduct`, `Scaled`, `Identity`, `Permutation`,
  `TensorProduct`, `Diagonal`, `Invertible`) but only `InvertibleOperator` has
  `.inverse()`. Each needs its inverse object: `Diagonal/Multiplication →`
  `DiagonalInverse` (reciprocal, the natural `DenseInverse` degenerate);
  `Product → (AB)⁻¹ = B⁻¹A⁻¹`; `Scaled → (1/α)·op⁻¹`; `Identity → Identity`;
  `Permutation → permutation inverse`.
- **The resolvents (P3) dissolve** — G-S into a `SweepOperator` on the REIFIED
  `M = (L+C−B_lower)`; windowed into the typed product `P @ A.inverse()` (fused
  via the substrate moment emit) — NOT new duck-typed classes (§17 W1/W2).
  This *replaces* the "base them on LinearOperator vs duck-typed" fork that
  started this discussion — they don't need their own class at all once they are
  the inverse of the correct forward operator. **The `SourceIteration` consumer
  then just holds `L.inverse()` and calls `.apply`** (the vestigial forward
  deletes).
- **`GreenOperator`** lands when we give `A_loss = L+C−S` (or the full within-group
  operator) an `.inverse()` that wraps the SI/Krylov driver — turning
  `SourceIteration`/`KrylovAcceleration` into the *implementation* of an operator,
  not a free-standing driver. This is the elegant endpoint: `k`-eigenvalue becomes
  `K = A_loss.inverse() @ F` with a `GreenOperator` inverse.
- **`as_matrix()` + `DenseInverseOperator`** promote `direct_eigenvalue` from a
  homogeneous free-function to production machinery — the CP method is the
  consumer that earns it.
- **`#280`** (the adjoint-inverse) is correctly deferred and now *explained*: it is
  free on the Krylov branch, hard on the sweep branch (needs `sweep_transpose` +
  the 2-D reverse walk).

### The round-trip discriminator becomes a permanent gate
`inv.apply(A.apply(x)) == x` for every `(A, A.inverse())` pair — the structural
proof that an inverse inverts *its own* forward (it REDs the confused resolvent,
passes the clean one). This is the keystone test of the whole taxonomy.

---

## 10. Cross-method pollination (this machinery is not SN-only)

- **MoC**: the ray-trace solve is **also** triangular forward-substitution (along
  characteristics) ⇒ a `TrackSweepOperator`, sibling of `SweepOperator`. This is
  the **second witness** that earns `SweepOperator` a generic `numerics`/`transport`
  home (the relocation trigger).
- **CP**: `[P]` is **dense by construction** ⇒ CP lives natively in the
  `DenseInverseOperator` / `as_matrix()` branch. CP is the consumer that promotes
  `as_matrix()` from a 0-D convenience to production.
- **Diffusion**: self-adjoint elliptic, **no triangular factor** ⇒ `.inverse()` is
  a CG/multigrid `GreenOperator` and `.H == self`. This is the **negative control**
  proving the factory keying isn't SN-overfit (it exercises Green + Dense, never
  Sweep).

---

## 11. Resolved design decisions (the open questions, pushed through)

1. **Strategy selection = type-as-structure (polymorphic `.inverse()`), NO runtime
   predicate.** Each operator's `.inverse()` override returns its natural inverse
   (`InvertibleOperator→Sweep`, `OperatorSum→Green`, `Diagonal/Mult→reciprocal`,
   `Product→B⁻¹A⁻¹`, `Scaled→(1/α)op⁻¹`, `Identity→Identity`, `Permutation→inverse
   perm`). The TYPE is the structure — adding an `is_triangular`/`structure`
   property would re-introduce stringly dispatch. The G-S variant is
   NOT structure — it is a `SweepOperator` on the schedule-masked `M` the SOLVER
   selects by geometry/performance; the windowed variant is the product
   `P @ A.inverse()` the solver composes (`_select_si_resolvent`/`_maybe_window`
   are the FACTORIES — §17 W1/W2). Two correctly-separated
   selection sites: operator = mathematical default; solver = performance strategy.
   (Realization caveat SUPERSEDED by §17 W2: `B_lower` is REIFIED — the boundary
   operator masked to the schedule's strictly-lower octant pairs — so
   `M = (L+C−B_lower)` IS a literal algebra object with its own `apply`.)
2. **`GreenOperator` WRAPS the driver; the drivers stay public primitives.**
   `A_loss.inverse()` (OperatorSum) auto-derives the split — its `InvertibleOperator`
   sub-term `(L+C)` is the preconditioner (`.inverse()` = the sweep), the rest the
   gains — and wraps a `SourceIteration` (Richardson) / `KrylovAcceleration` (GMRES)
   engine (Pattern 5; re-implementing = Smell #16). Driver-choice is a GreenOperator
   strategy. Endpoint: `K = A_loss.inverse() @ F` (the grand-report normal form).
3. **Seed kwarg STAYS: `inverse().apply(q, *, initial_guess=None)` as documented
   §38 application-context** — never changes the result (`A⁻¹q` regardless), home is
   per-apply because SI's `ψ_prev` changes each iterate. **SCOPE (amended, §17 A1):
   the ruling is NOT universal — it is the SAME claim as §13's S-direct and inherits
   its geometry split.** Slab/cylinder: TRUE (bit-exact seed-independence,
   re-verified). SPHERE pre-#282: FALSE — the lagged Carlson pole closure makes
   `.apply(q, seed)` a 2-arg relaxation step (one Picard iterate), i.e. the kwarg is
   temporarily IDENTITY-BEARING there. That is an xfail-until-#282 DEFECT note, not
   an accepted "hybrid with documented warm-start dependence".
4. **`as_matrix()` = dimension-gated `LinearOperator` method** (apply-to-basis
   default + structure-specific overrides like the streaming sparse-triangular
   assembly, `sweep_graph.py:66`). Refuses above a configurable threshold
   (`MatrixTooLarge`). `MatrixInverseOperator` (`lu_factor`/`lu_solve` on
   `as_matrix()`) is `A.inverse()` for small/dense ops; latent consumer
   `direct_eigenvalue`; CP earns the promotion. **Partiality framing (§17 A2):
   `MatrixTooLarge` is a RESOURCE effect on a TOTAL functor (every operator HAS a
   matrix; it may just not be materializable here) — NOT a structural restriction
   idempotent like `is_invertible`. Do NOT mint an `is_materializable` predicate;
   the gate is a size precheck, no structure read.**
5. **Naming = one axis (realization) + greppable role token.** Lean
   `KrylovInverseOperator` over `GreenOperator` (Green names the kernel, breaks the
   algorithm-named family — keep only with a documented axis-mix). **"Resolvent" is
   currently MISUSED** (`_GaussSeidelResolvent`/`_MomentWindowedResolvent` are sweep
   configs, not `(A−zI)⁻¹`); they dissolve, freeing `ResolventOperator(z)` for the
   genuine shift-invert α-family. (User's call on final names.)
6. **Re-scope: the resolvent dissolution AND the SI migration are ONE change.**
   Once resolvents are configured `SweepOperator`s, **`SourceIteration` holds an
   inverse operator and calls `.apply`** — never `.solve`/`.inverse()`. The solver
   builds `L_inv = (L+C).inverse()` (or the configured G-S/windowed sweep) once and
   hands it to SI (`ψ_{k+1} = L_inv.apply(rhs, initial_guess=ψ_prev)`); the
   `_has(CAP_SOLVE)`-guard + `.solve`-call + duck-typed-resolvent tangle evaporates.

## 12. Re-scoped phase order (supersedes the carve plan's P3/P4/P5)

> **🏁 CAMPAIGN MERGED TO MAIN @ `1729647` (2026-07-02, `--ff-only`;
> branch `refactor/inverse-as-operator` DELETED).** All of §12 — steps
> 1–6, 5b, and the 6-tail carve P5 — plus the vv-principles Mode-12
> uplift are on `main`. Merge gate: the FULL serial not-slow tree
> (5805 passed / 0 failed / 54 xfailed, 51:18) + ratchet EXACTLY 145 +
> Sphinx `-E -W` exit 0. `main` is ahead of `origin/main` (pushes held
> by the user); the `d82fa77` `Closes #285` trailer fires at push time.
> Everything below is the as-built record — archaeology, not active
> state. Open follow-ons live on GitHub: #280 (redesigned onto
> `A.H.inverse()`, comment 4871342146), #282, #283, #226 (burn-down
> continues).

1. **✅ DONE @ `e7115a2` (2026-07-01).** Delivered deliver-all (no retirements;
   the user-gated retire fork was not triggered). Leaves = generic
   `InverseOperator` (solve-delegation, bit-identical — NOT §9's reciprocal
   sketch; rationale in its docstring + `agent-memory/elegance-enforcer/
   issue_226_inverse_as_operator_rulings.md`); composites = structural;
   SweepOperator axis completed (solve = forward, involution by identity);
   shim forwards inverse+solve. Gates = spec PART II §12 (`tests/numerics/
   test_inverse_universal.py` + 3 extensions); 7-mutation matrix RED-verified;
   tiers 1671/0; pyright 153 unchanged. Elegance: no blockers; twin back-half
   documented with collapse trigger (mixin at 3rd sibling, steps 4–5).
   ~~Close the `.inverse()` promise-gap~~ across the remaining `is_invertible`
   advertisers — there are NINE, not eight (§17: `_BoundBoundaryOperator`,
   `geometry/boundary/_bound_compat.py:108`; decide shim-forwarding). The §17
   consumer audit found ZERO production consumers of
   `.inverse()`/`is_invertible`/`SupportsInverse` today, and 6/8 advertisers'
   `.solve` is ALSO production-unconsumed — but `is_invertible`'s docstring
   (operator.py:388) already PUBLISHES the `.inverse()` contract, so each type
   gets an explicit **deliver-vs-retire** decision: deliver the trivial natural
   inverse (Diagonal/Identity/Permutation/Scaled/Product — few lines each,
   pinned by the round-trip gate), or retire the WHOLE advertisement surface
   (`is_invertible` + CAP_SOLVE + `.solve` together, keeping the faithfulness
   keystone `is_invertible ≡ CAP_SOLVE` coherent). User call per type on the
   retire side.
2. **Generalize `SweepOperator`** with schedule/fold config + **REIFY
   `M = (L+C−B_lower)`** (the masked boundary term, §17 W2) so the G-S sweep
   round-trips against its OWN forward. The windowed path becomes the typed
   product `P @ A.inverse()` (§17 W1) — the moment emit stays substrate-level.
   **✅ DONE @ `cc293ef` (2026-07-01, `refactor/inverse-as-operator`).**
   Landed: `SNBoundaryOperator.split(schedule) → BoundarySplit(lower, upper)`
   (`SNMaskedBoundaryOperator`; split law = `SweepSchedule.lower_inflow_rows`
   — a row is lower iff its octant group is swept strictly after its face's
   reflect group); `M = (L + C) - parts.lower` fuses via
   `InvertibleOperator.__sub__` → NEW `ScheduledInvertibleOperator`
   (honest OperatorSum leaves; solve = the ONE `_solve_timed_full_field`
   body with `(schedule, reflect)` threaded through the representation's own
   `sweep` door — the solver-side `_solve_scheduled` plumbing twin RETIRED);
   `SweepOperator.inner` widened to the union (NO 3rd wrapper — the schedule
   config rides the FORWARD, not the inverse, so "generalize SweepOperator
   with schedule config" resolved to *widening*, the honest reading of §7).
   Driver congruence: `_select_si_resolvent` G-S arm returns
   `(M, (S, B_upper))` — `SourceIteration` unchanged; `_GaussSeidelResolvent`
   DELETED. W1: NEW `orpheus/sn/operators/windowing.py` —
   `BulkAnalysisOperator` (block coisometry P) `@ A.inverse()` →
   `WindowedSweep(OperatorProduct)` whose `apply` fuses via the substrate
   moment emit; the INHERITED product body is the deforestation oracle;
   public `solve_moments` RETIRED (3 sites → the one composition);
   `_MomentWindowedResolvent` survives as a thin driver adapter (dies step 3
   — **step 3 MUST delete it**: its apply→angular/solve→moments two-hat
   transiently re-introduces the codomain mismatch, enforcer note).
   In-sweep reflect = ADDITIVE row-masked `reflect_rows_inplace` (the
   inhomogeneous row `z_in = y_row + (Bz)_row`; the old whole-face overwrite
   dropped `y_row` and stamped fresh values on lagged rows). **Domain
   finding → issue #284**: the walk realizes A⁻¹ on the source subspace
   `{y.out-rows = 0}` (all production rhs; family-wide, also (L+C).solve) —
   contract for non-source consumers (M⁻¹A preconditioning) decided there.
   Gates (spec §13 + its ✅ note): W2 round-trip **5.2e-16 bulk / 4.4e-16
   trace** (was 2.667); split bit-exact; W2-FP G-S≡Jacobi (LS-S4, Mode-9-
   safe); W1 fused≡deforested 1.8e-16; mutations M-SPLIT-DIR/M-SPLIT-PART/
   M-DEFOREST in-tree and biting (spec'd M-SPLIT became UNREPRESENTABLE —
   the mask single-sources apply AND reflect). Tier 2977/1→green (one spy
   signature migrated); DD snapshot drift bit-identical pre/post; pyright
   re-baselined 153→**152** (net −1); Sphinx `-W` clean (archivist updated
   operator_algebra/discrete_ordinates/cross_section_data + Development-
   history entry); elegance review SHIP (fixes applied).
3. **Solver builds the inverse operator; SI/Krylov apply it** — dissolves the
   duck-typed resolvents AND does the consumer migration in one move (the old P3).
   **PREREQUISITE (§17 F4): land a FAST-LANE seed-threading gate first** — the
   strong het-2G sphere catcher is `@slow` (deselected by the canonical
   `-m "not slow"` run); a simulated seed-drop (|Δk| ≈ 3.5e-2 on sphere)
   otherwise reddens only a fragile 1G margin, and cylinder's seeded-value gates
   are vacuous for seed-drops (the seed telescopes away) — so a Mode-11 path-spy
   on the `initial_guess` threading is the load-bearing fast guard.
   **✅ DONE @ `1ab7429` (2026-07-02, `refactor/inverse-as-operator`).**
   Landed: PREREQ spy `tests/sn/solve/test_seed_threading_spy.py` (wraps
   `InvertibleOperator.solve` — route-invariant, green pre-AND-post rewire;
   teeth M-SEED-DROP/ZERO/STALE pre-rewire + M-PROBE post-rewire, all RED
   <1 s; the @slow sphere equivalence gate tagged
   `@catches("ERR-026","M-SEED-DROP")`). REWIRE:
   `SourceIteration(A_inv, *gains)` consumes the pre-built inverse through
   NEW `SupportsSeededApply` (canonical `apply(rhs,*,initial_guess=None)`;
   not runtime_checkable) — the `_solve_accepts_seed` inspect-probe,
   `_solve_with_seed`, and the ctor CAP_SOLVE gate DELETED (invertibility
   obligation moved to the `.inverse()` builder; apply-only step operators
   accepted BY DESIGN — the windowed product's shape);
   `_MomentWindowedResolvent` DELETED (`_maybe_window` returns
   `BulkAnalysisOperator @ sweep` directly — the factory spells `P @ A⁻¹`);
   NEW `_seeded_inverse(A)` single-sources the cast narrowing + the SCOPED
   conformance claim (reachable inverses only — composed
   `OperatorProduct.inverse()` lacks the kwarg TODAY; structural-mixin vs
   convention = **#285**, decided with steps 4–5: Green/Matrix MUST join
   the contract); `KEigenvalue(A,S,F)` guards `is_invertible` at ctor;
   Krylov fallback rewired (type:ignore deleted). Family kwarg added to
   `InverseOperator.apply` + `WindowedSweep.apply` (accepted-and-ignored,
   documented per type). **THE A-RENAME folded in** (user-caught):
   numerics/iteration.py had L = the streaming-collision COMPOSITE
   (R-1-era), colliding with the project convention (L = streaming leaf;
   invertible LHS = the A family) — all params/attrs/messages/equations
   renamed (`A_inv`/`A`); exposed + fixed two latent wrongnesses
   (KEigenvalue docstring "A_loss = L" → `A_loss = A − S`; the L0 test's
   dense `A` local = A⁻¹F → `A_inv_F`). NEW gate:
   `test_2d_windowed_product_over_gauss_seidel_M_equals_post_projection`
   (the windowed×G-S corner `P @ M⁻¹` fused ≡ deforested ≤4·N·eps,
   B_lower non-degeneracy guard — production's 2-D default composition,
   previously pinned only at ℓ=0). Consumer pins migrated per spec §3
   (single-primitive contract now asserts SweepOperator + `.inner` type
   identity; GAP-3's raise-leg retired WITH the requirement it pinned →
   `test_invertibility_obligation_lives_at_the_inverse_builder` +
   `test_keigenvalue_requires_invertible_A`). Tier 2981/0; pyright
   re-baselined 152→**148** (numerics 24→21, sn 103→102 — probe/adapter
   noise retired); Sphinx `-W` clean (6 dead adapter xrefs rewritten
   tense-triaged; NEW theory section `inverse-application-driver`;
   step-3 Development-history entry; step-2 `_GaussSeidelResolvent`/
   `solve_moments` debt audited — all historical, correct). Enforcer
   SHIP-WITH-FIXES, all conditions closed. **Step-4/5 obligations:**
   (a) Green/Matrix join `SupportsSeededApply` (accept-and-ignore, like
   InverseOperator) per #285; (b) the wrap-delegate mixin extraction at
   the 3rd sibling folds the seeded-apply back-half in; (c) OLDER doc
   sections still use L-for-composite in live prose
   (operator_algebra.rst:5424, discrete_ordinates.rst ~9331/11108) —
   one archivist tense-triage sweep, user-scheduled.
4. **`GreenOperator`** wrapping the drivers (when `A_loss.inverse()` gets a
   consumer). W3 VERIFIED VIABLE (§17): iteration.py→operator.py is one-way
   today; a `green_operator.py` leaf module + a late import inside
   `OperatorSum.inverse()` mirrors the SweepOperator↔InvertibleOperator
   precedent (streaming.py:793). MUST land with an ORDERING RULING: `L + C`
   fuses to `InvertibleOperator` (→Sweep) but `C + L` builds a plain
   `OperatorSum` (→Green) — same math, different algorithm by operand spelling;
   guard or document per the #261 canonical-ordering rule.
   **✅ DONE @ `9333305` (2026-07-02, `refactor/inverse-as-operator`;
   user unlocked the consumer gate: "let's start step 4").** Landed: (a) the **3rd-sibling MIXIN extraction fired** —
   NEW `InverseWrapMixin` (numerics/operator.py) carries the wrap-delegate
   back-half (capabilities / domain↔codomain swap / `solve→inner.apply` /
   `is_invertible` / `inverse()→inner`, `_ForwardT`-typed so each
   sibling's involution returns its narrow forward) + the **ABSTRACT
   seeded-apply signature** — **#285 resolved STRUCTURAL** (pyright
   LSP-rejects a kwarg-less override, M-MIXIN-KWARG cp-backup-mutation
   verified; ABCMeta blocks a missing `apply`); `_SolveBackedLeaf` renamed
   `_InvertibleForward` (now the mixin bound); InverseOperator +
   SweepOperator rebased (keep guard/apply/repr; no-op proven — 125
   landed gates green unchanged); residue for step 5: composed
   `OperatorProduct.inverse()` still outside the family (#285 remaining
   scope). (b) NEW **`numerics/green_operator.py`** (W3 leaf placement):
   `GreenOperator` + `ConvergenceFailure`; left-spine split derivation
   (EXACT-type walk — `InvertibleOperator` subclasses are leaves, so the
   MRO shadow survives flattening; M-GRN-FLATTEN verified), leading term
   = preconditioner via its own structure-keyed `.inverse()`, gains
   negated with the `Scaled(-1,X)→X` unwrap (proven pure deforestation by
   the no-op mutation control); wraps `SourceIteration` (§11.2 — zero
   re-implementation). **Spec §18.A integrated WITH a delta:** the
   promise is the TRUE relative residual `‖(A−B)ψ−q‖/‖q‖ < tol`, and it
   DRIVES a refinement loop (re-seed the driver until met or the TOTAL
   budget exhausts, then raise) — check-only would falsely raise for
   every ρ>1/2 since increment-stop delivers only ρ/(1−ρ)·tol (the
   spec's own het-slab gate would red); the driver stays the sole
   iteration engine. Two adversarial FP shapes found by the divergence
   tooth are caught + documented in-code: increment=nan, and the
   DENOMINATOR-overflow false convergence (`res=finite/inf=0.0` onto a
   ~1e154 iterate). (c) **`OperatorSum` contract change:**
   `is_invertible = a.is_invertible` (leading-term rule — §18.B honesty
   note in the docstring: "preconditionable at this spelling", ordered
   per the #261 `L + C` canonical-ordering rule at streaming.py:510),
   lockstep conditional `CAP_SOLVE` (faithfulness keystone),
   `solve = inverse().apply` (transitional; retires at P4),
   `inverse()` → late-import Green, base annotation
   `LinearOperator[Codomain, Domain]` (the honest factory face — the
   InvertibleOperator→Sweep override stays a legal narrowing). (d) **The
   ordering ruling landed with all four §22.1 edges GATED:** `(L+C)` →
   SweepOperator (MRO-shadow pin); `A_loss` → Green converges;
   `C + L` constructs (leading C invertible — the algebra can't read
   ρ(C⁻¹L)>1) then raises ConvergenceFailure LOUDLY end-to-end on the
   real operators; `(−S)+(L+C)` refuses at construction naming the
   canonical ordering. GATES: NEW `tests/numerics/test_green_operator.py`
   (G-I1 + dense-LU anchor + δ_j kernel fold, G-I2 object-identity
   involution, back-half-solve anchor, **G-Neumann** partial sums with
   the EXACT 4-cycle decay ratio, **G-reciprocity** via transposed
   operands (no `.H`), divergent-raise + convergent control,
   near-critical §18.A gate with the budget calibrated BETWEEN
   increment-stop (~460) and refinement-close (~920), the §23 Green
   seed spy, refusal, tol=0) + NEW
   `tests/sn/operators/test_green_operator_sn.py` (G-Neumann-L1 on a
   het-2G-vacuum slab with a **trace-consistent MANUFACTURED anchor** —
   resolves the #284 source-subspace caveat the spec's dense-LU sketch
   missed; driver BIT-IDENTITY vs hand `SourceIteration(sweep, S)`;
   dispatch/MRO pins; refusal; the C+L trap) + §24.2 pin migrations
   (survival `(L+C)−B` False→True renamed
   `test_sweep_vs_green_inverse_keyed_by_type`; predicates s_both flip +
   BOTH asymmetric ordering fixtures) + G-I3 faithfulness + universal
   dispatch row + §19.2 static pins (`assert_type` +
   `SupportsSeededApply` conformance, all three siblings). TEETH: **14
   verified** — 12 bite (SIGN/SWAP/FLATTEN/TOL/INCREMENT/SEED/ORDER/
   SUM-INV-B/CAPDRIFT/MIXIN-INV/MIXIN-SOLVE/MIXIN-KWARG), 2
   designed-green controls (UNWRAP no-op; SEED blindness proof — the
   landed §14 spy stays green while only the §23 spy catches the drop).
   Pyright ratchet EXACTLY at baseline 148 (0 new); Sphinx `-W` EXIT 0.
   Also fixed in passing: operator.py's own module docstring still
   spoke the pre-A-convention `(L−S−F)ψ=q` — now `(A−S−F)ψ=q`.
   **Consumer ruling:** the gates construct Green directly; KEigenvalue
   was examined and deliberately NOT rewired — its inner is a
   warm-started, budget-bounded RELAXATION (partial convergence by
   design in nested iteration), a DIFFERENT contract from Green's
   converged-or-raise promise; `K = A_loss.inverse() @ F` is the normal
   FORM, production power iteration is its inexact realization
   (documented in the module docstring). Production consumers arrive
   with #200 (preconditioner algebra) / future normal-form spellings /
   diffusion (.inverse() = CG Green, the negative control).
5. **`MatrixInverseOperator` + `as_matrix()`** (§13 name) — promote
   `direct_eigenvalue`; CP. The W4 threshold value + dense-vs-sparse return
   keying are decided HERE (§17 A2 already rules the partiality framing).
   **✅ DONE @ `d82fa77` (2026-07-02, `refactor/inverse-as-operator`; closes
   #285).** As built:
   - **`as_matrix()`** = universal `LinearOperator` base method (the functor
     OUT, §2's fourth arrow): the promoted `_as_dense` apply-to-basis loop,
     C-order basis + C-order output ravel, rectangular-honest (output dim
     emerges from `apply`). Resolution via the single-source
     `_resolve_basis_shape` (explicit wins → `domain.shape` → `ValueError`
     naming both remedies — the ILL-POSED arm, §27.C-discriminated from the
     resource arm). **W4 DECIDED:** `MatrixTooLarge(RuntimeError)` at
     `prod(basis_shape) > max_dimension`, default **4096** (134 MB / 4096
     applies), per-call knob; **dense-vs-sparse keying DECIDED:** dense
     ndarray is the only realization; the sparse-triangular structural
     override arrives WITH its 3-D consumer (sweep_graph.py:66), return
     type generalizes then (defer-until-consumer).
   - **`MatrixInverseOperator`** (new `numerics/matrix_inverse_operator.py`,
     leaf; NO late-import — no automatic `.inverse()` routes here, direct
     construction IS the §3 strategy-override seam): 4th `InverseWrapMixin`
     sibling; eager `lu_factor(inner.as_matrix(...))` = the guard (size /
     squareness / EXACT singularity — scipy's singular-WARNING silenced and
     replaced by a loud zero-pivot `LinAlgError` at construction;
     near-singularity priced into κ(A), not refused); `apply` = one
     backsolve, seed accept-and-ignore; `as_matrix()` override = batched
     `lu_solve(lu, I)`. **Deliberately no `inner.is_invertible` read** —
     values-not-structure; gate-witnessed: the leading-non-invertible sum
     Green REFUSES is exactly what Matrix inverts (ndarray analog per spec
     §27.B — FullField is out of as_matrix scope).
   - **§13 M-row SHARPENED (spec §27.A, supersedes the "no as_matrix()"
     parenthetical):** as_matrix is universal now, so a Green satisfies
     `G·A ≈ I` too — at DRIVER tol; the name-earner is the **machine·cond
     PRECISION GRAIN** (M-materialise both ways + M-direct residual at
     `K·ε·cond`, seed bit-identity), proven DISTINGUISHING by an in-gate
     Green contrast (its materialization residual > machine grain,
     < driver tol).
   - **Mixin bound relaxed to its true minimum:** new `_WrappedForward`
     Protocol (domain/codomain/apply — all the back-half consumes);
     `_InvertibleForward` narrows it (+is_invertible+solve) as
     `InverseOperator`'s OWN contract. Pure widening — 125+ landed family
     gates green unchanged, ratchet exactly 148.
   - **#285 product residue CLOSED:** `OperatorProduct.inverse()` →
     `InverseOperator(self)` (was the raw reversed product): bit-identical
     action (= `b.solve(a.solve(q))`), seeded kwarg now accepted (the #285
     TypeError repro flips), involution STRENGTHENS to object identity.
     Two-kinds taxonomy documented: wrap-delegate family (mixin, canonical
     seeded signature) vs ALGEBRA-CLOSED inverses (perm→perm, identity,
     scaled→scaled — first-class forwards, stay unwrapped). `Closes #285`
     goes in the commit body.
   - **Consumers:** homogeneous `_as_dense` RETIRED → 2 call sites on
     `as_matrix(basis_shape=(ng,1))`, k∞/flux byte-identical (landed SymPy
     pins untouched + a HOMOG-EQUIV local-oracle gate);
     `dense_per_material` docstrings rewritten to the honest storage-side-
     oracle role (zero production consumers — the independent side of the
     ASM-ORACLE pair); `direct_eigenvalue` stays ndarray-pure; **full
     operator spelling = step 5b below (✅ landed)**; CP the next
     production method in line.
   - **Verification:** spec PART IV (§27–§35) — 18 new gates
     (`tests/numerics/test_matrix_inverse_operator.py`) + extensions
     (universal registry row §30.10, #285 repro/value/involution/
     algebra-closed §31, faithfulness M-I3, 4th-sibling static pins);
     **14/14 mutations verified** (13 in-process + M-MINV-KWARG CLI
     `reportIncompatibleMethodOverride`); 117 combined step-4/5 gates
     green; tier 2998/0; ratchet 148. As-built deltas vs PART IV: the
     HOMOG-EQUIV gate hosted in the numerics file (the homogeneous file's
     file-level `l1+verifies` marks would write FALSE equation edges for a
     software gate); MINV-DIRECT gained a non-symmetric leg (aligns
     M-MINV-LUTRANS's activating config with its target gate — the spec's
     diag fixture is transpose-blind); the TOOLARGE per-call edges use
     tighten/lift on a 5-vector + default refusal at 4097 (same edges, no
     200 MB materialization); the non-square gate pins the domain-language
     message fragment ("SQUARE materialization") so the guard-delete
     mutation cannot hide behind scipy's own "expected square matrix";
     Green-contrast sanity at 1e3·tol (10× was tight under the G·A
     column-error cond amplification).
5b. **Homogeneous full operator spelling — MatrixInverseOperator's first
   production consumer** (task #138).
   **✅ COMPLETE — COMMITTED @ `394d8c0` (2026-07-02; 10 files,
   +860/−306; branch `refactor/inverse-as-operator`).** As-built:
   - **The carve (SC1→SC2→SC3; verification plan
     `step5b_homogeneous_k_operator_verification.md`, test-architect
     `adb6cc8773e86928e`):** `dominant_eigenpair(M, *, imag_tol)` EXTRACTED
     in `numerics/eigenvalue.py` — the ONE home of the Perron–Frobenius
     validation (argmax-real, complex-reject, real-cast, sign-normalize via
     the new `_sign_normalised`, which also folded RQI's inline flip —
     byte-identical); `direct_eigenvalue` delegates (KEPT: the (A,F)-posed
     sibling engine + RQI oracle, now zero production consumers — NOT
     superseded, and NOT under the fuller-view-oracle exception;
     test-architect AUG2 concurrence). `_assemble_loss_matrix` →
     `_assemble_loss_operator` (returns the un-materialized `C − K_iso`);
     `solve_homogeneous_infinite` spells
     `K = MatrixInverseOperator(loss, basis_shape=(ng,1)) @ production` +
     `dominant_eigenpair(K.as_matrix(basis_shape=(ng,1)))` — the explicit
     ctor IS the strategy choice (the structure-keyed `loss.inverse()`
     would return the ITERATIVE Green; the exactly-solvable 0-D problem
     earns the direct realization).
   - **Re-baseline PRINCIPLED-EQUIV** (batched `gesv` → held `lu_factor` +
     per-column `lu_solve` through the product's apply-to-basis
     `as_matrix`): measured bit-identical on this host; gated rtol=1e-12
     (κ·ULP-portable), never `==`.
   - **Gates (G1–G10, all green):** new `TestDominantEigenpair` (direct
     surface, teeth inherited from `TestGatesHaveTeeth`'s A=I bank) + the
     one-home neuter-proof PAIR (mutating the home makes the delegated
     surface's guard VANISH — relocated, not copied); the Mode-11 spy
     REWIRED `direct_eigenvalue`→`dominant_eigenpair` with the pre-impl
     skip arm converted to a hard `_require` (the SILENT-SKIP hazard: the
     old gate would have skipped green post-carve); first-consumer
     liveness spy on `MatrixInverseOperator.apply`; the K-resolvent
     intrinsic gate `K.as_matrix() ≡ np.linalg.solve(A, F)` (also the
     None-tolerant space-guard witness); the byte-equality gate retuned to
     cross-engine rtol=1e-12; the fused-oracle gate migrated to the
     operator spelling; `_solve`/`DirectEigenvalue` pre-impl scaffolding
     retired.
   - **TEETH SHARPENING (empirical T1–T5, job-tmp `step5b_teeth.py`):**
     a factor swap `F·A⁻¹` and a resolvent TRANSPOSE are **spectrally
     invisible** (similar matrices ⇒ |Δk| = 0.0 EXACTLY) — every k-level
     gate (cross-engine, SymPy anchor) is designed-green on them; the
     COMMITTED CATCHER is the matrix-level K-resolvent gate (both move the
     matrix O(1)). Dropped-factor DOES move k O(1); K_iso sign-flip reds
     the fused oracle. This CORRECTS the verification plan's G5 teeth row
     (it claimed the swap moves k O(1)). Full derivation:
     homogeneous.rst `spectral-invisibility` section.
   - **Docs (archivist `aaaa06668eb2b9713`, Sphinx `-E -W` exit 0):**
     homogeneous.rst onto the K-operator story (+ spectral-invisibility §
     with derivation + principled-equiv note + engines table),
     api/numerics.rst, operator_algebra.rst consumer-loop closed
     present-tense; SN Dev-history step-5b row (main agent); verification
     matrix regenerated (eigenvalue 30→39, homogeneous 12→14 nodes).
     Archivist flags: a vv-principles reference.md uplift candidate
     ("value gate blind to a mutation class for spectral-similarity
     reasons — pin the OBJECT, not just its spectrum") left to the user's
     curation call.
   - **Elegance (enforcer `a3dc054b61e8a7923`): SHIP — zero code defects,
     zero violations.** All three findings folded pre-commit: the
     test-module docstring swept (`DirectEigenvalue` CamelCase retired,
     `dominant_eigenpair` named); the severed-consumer parenthetical in
     `direct_eigenvalue`'s docstring repointed; the falsifiable
     per-column narration trimmed from the solver comment. Watch-item:
     when CP `[P]` (§14b) mints a meshless FunctionSpace, both
     `basis_shape=(ng,1)` args collapse — note it on the CP issue when
     it opens (no standalone issue; deliberate space-anonymous design).
   - **Verification:** trio 77/0; blast radius `tests/numerics/ +
     tests/homogeneous/` 848/0 (+1 pre-existing unrelated skip
     `test_estimators_as_functionals.py:207`); ratchet EXACTLY 145
     (numerics 21 / sn 99 / transport 25); Sphinx `-E -W` exit 0.
6. **Frozenset retirement + composition algebra + docs** (#280 explained as the
   sweep-vs-Krylov adjoint asymmetry).
   **✅ COMPLETE — W1+W2 @ `f4919b1`; W3 COMMITTED @ `cb62310`
   (2026-07-02; 12 files, +914/−280).** W3 record: the FULL §41 bank
   RED-verified under `-O` (M-ADJ-PROP / M-ADJ-FORGE / M-INV-FORGE /
   M-BRIDGE with perfect 10-row selectivity / M-PROD-REROUTE + the
   M-PROD-COMMUTING designed-green control / M-LITERAL-STRING false-red
   demo / M4′ / M-GUARD-DELETE-PYRIGHT 0→1→0 — verdict table stamped at
   spec §41); §39.3 static pins landed
   (`_typeguard_bridge_narrowing_static_pins` in
   `test_operators_apply_typed.py` + the **inverted-polarity**
   `_static_zero_inverse_unspellable.py` — green while unspellability
   holds, reds via `reportUnnecessaryTypeIgnoreComment` if `inverse()`
   ever lands on Zero/base; spec §39.3 stamp records the delta + the
   user-veto flag on the rule-scoped ignore) with both teeth
   CLI-verified (M-BRIDGE-ANNOT 3→8, M-STATIC-UNSPELLABLE 0→1);
   archivist docs pass (operator_algebra.rst three-layer rewrite +
   6 anchors, Dev-history entry, per-axis dead-xref repoints, a
   PRE-EXISTING FALSE "S/F not adjointable" claim corrected and
   empirically confirmed — all five leaves + the full
   `(L+C−S−F−B).H` construct, Sphinx `-E -W` exit 0);
   elegance-enforcer PASS (0 must-fix; POLISH-1/-2 applied — the three
   `getattr(…, "is_invertible", False)` reads flattened to direct
   predicate reads in operator.py:1151 / KEigenvalue / GreenOperator;
   the step-3 `inverter`-narrative staleness tracked as a comment on
   #226); + 2 main-agent cross-check doc catches (present-tense
   `OperatorSum.solve` in index_convention.rst:1019 and
   operator_algebra.rst:9206 — outside the archivist's
   retirement-vocabulary grep, caught by the enforcer's independent
   sweep). Full tier 3853/0 (14:23); post-polish targeted 139+43 green;
   ratchet EXACTLY 145 throughout. As-built record (W1+W2):
   - **Design "C+B" (USER-LOCKED 2026-07-02, after a full A-vs-B +
     layer-2 deliberation the user drove to the TypeGuard resolution):**
     the "per-class-needs-casts vs base-declaration-demotes-static-errors"
     compromise was a FALSE DICHOTOMY — it conflated STRUCTURAL
     non-invertibility (the TYPE has no inverse: Zero/masks/dyads/source
     leaves → `inverse()` NOT declared, misuse is a pyright error) with
     VALUE-dependent (zero-coeff mult, non-invertible-headed sum →
     declared, raises `NotInvertible` eagerly). **Design C:** per-class
     methods + PEP-647 `TypeGuard` bridges `invertible()`/`adjointable()`
     (TypeGuard NOT TypeIs — value-dependence permits only the
     one-directional promise; free functions because narrowing needs a
     call expression and a method narrows its first EXPLICIT arg, never
     self) — the runtime check IS the static permission, guards are
     type-load-bearing (guard-delete → CLI pyright RED);
     `SupportsInverse`/`SupportsAdjoint` PROMOTED to narrowing targets
     extending `LinearOperator` (casts die; the Protocols' original
     "static contract + runtime property" charter finally fulfilled by a
     CHECKED bridge). Base-hosting rule: base hosts a method IFF a
     universal realization exists (`.H` yes — one generic wrapper, now
     with the EAGER `MissingAdjoint` gate; `inverse()`/`apply_transpose`/
     `solve` no). **Design B:** `solve` pruned to the native realizations
     (the `_InvertibleForward` face = what the wrap-delegate family
     wraps): DELETED on Sum (docstring promise executed)/Identity/
     Permutation/Scaled/TP/`_BoundBoundaryOperator`; KEPT on Diagonal/
     Mult/sweep composites/mixin/Product — Product.solve RE-ROUTED
     through factor `.inverse().apply` (bit-identical per kind, total
     over solve-retired factor kinds). `MissingCapability` → two
     TypeError successors `NotInvertible`/`MissingAdjoint` (W1 bridge:
     NotInvertible born subclassing MissingCapability so the 26 landed
     raises stayed green; W2 re-parented). Composition apply-guards stay
     eager as plain TypeError; composer transpose laws guard-narrow
     (MissingAdjoint replaces raw AttributeError).
   - **Two real bugs caught by the net, fixed root-cause:**
     (1) `_seeded_inverse` crashed on ALGEBRA-CLOSED preconditioner
     heads (`(I+B).inverse().apply` — pre-existing since step 4; the §40
     sum-row fixture exposed it) → two-kinds dispatch at the ONE home
     (`_wrap_delegate_member` TypeGuard + accept-and-drop
     `_SeededExactApply`); (2) the carve deleted `RankOneOperator`'s
     per-instance caps — its ONLY adjointability advertisement — leaving
     **F† dead** behind the TP factor guard; BOTH migration agents
     converged independently, pinned as strict-xfails → the
     `is_adjointable` override landed → pins flipped → markers deleted.
   - **Verification:** keystone v2 (predicates.py rewritten; scaffold
     deleted) + §40 re-route gates (Mode-11 sentinel + dense anchor +
     5-row factor-kind matrix vs the PRE-carve npz baseline, direct rows
     array_equal / Green row driver-tol) + the full 127-read/36-file
     §37.1 migration + §37.3 completeness re-grep ZERO + the two
     S†-canary precondition rewires. Full tier `-O` serial **3853/0**;
     ratchet re-baselined **148 → 145** (sn −3; a transient +3 from a
     narrowing-losing guard swap in the dyad transpose was caught by the
     exact-floor discipline — the isinstance IS the narrowing).
   - **W3 remaining (next commit):** the §41 mutation bank (M-ADJ-PROP /
     M-ADJ-FORGE / M-INV-FORGE / M-PROD-REROUTE + M-PROD-COMMUTING
     designed-green / M-BRIDGE / M-LITERAL-STRING / M4′ /
     M-GUARD-DELETE-PYRIGHT via cp-backup), §39.3 static pins +
     `Zero().inverse()`-unspellable pyright file, Sphinx theory pass
     (archivist), elegance-enforcer review.
6-tail. **Carve P5 — composition return types + docs + #280 redesign**
   (task #135). **✅ COMPLETE — COMMITTED @ `70da74f` (2026-07-02;
   3 files, +263/−58; branch `refactor/inverse-as-operator`, now 11
   commits ahead of main).** As-built:
   - **SC1 (the one annotation):** `ScaledOperator.inverse()` →
     `"ScaledOperator[Codomain, Domain]"` (operator.py:1396) — the audit of
     ALL NINE `.inverse()` surfaces found exactly ONE parameter-dropping
     annotation. The rest are precise or documented-loose BY DESIGN:
     `SupportsInverse`/`OperatorSum` keep the loose family face
     `LinearOperator[Codomain, Domain]` (siblings, not a hierarchy);
     Product/Diagonal → the unparametrized-by-design `InverseOperator`;
     Identity `[Domain]`-precise; Permutation/TensorProduct unparametrized
     classes; the mixin returns `_ForwardT` (involution typed per sibling).
     NO further parametrization built — the bare LEAVES gate any deeper
     static gain (a separate leaf-parametrization campaign, deliberately
     not started).
   - **SC2 (the pin bank):** `_composition_algebra_return_type_static_pins`
     in `test_operators_apply_typed.py` (pyright-only, never run — the
     file's established idiom): constructors (`+`/`-`→Sum, `*`/`-A`/`/`→
     Scaled, `@`→Product with Cmid solved away, `.H`/`adjoint()`→carrier
     swap, `**`→endo), the TWO inverse kinds (algebra-closed→themselves on
     swapped carriers; wrap-delegate→InverseOperator; sum→the loose face),
     mixin involution (`green.inverse()`→OperatorSum,
     `matrix_inv.inverse()`→LinearOperator — also guards each sibling's
     `InverseWrapMixin[X]` binding). Banner records the honest boundary:
     NO `A.H.inverse()` pin ahead of #280 (`_AdjointOperator` deliberately
     lacks `inverse()` — the mixin's recorded deferral).
   - **Teeth:** M-SCALED-BARE red-verified — reverting SC1 reds the pin
     via `reportAssertTypeFailure` ("expected ScaledOperator[AngularFlux,
     ScalarFlux] but received ScaledOperator[Unknown, Unknown]").
   - **Gates:** scoped CLI pyright = the 3 pre-existing `BC.reflective`
     only; ratchet EXACTLY 145 (21/99/25) throughout; runtime 102/0
     (touched + neighbors), 17/0 re-run post-folds.
   - **Enforcer (fresh instance): SHIP, zero code defects.** All banner
     claims verified against the live tree; twin-path ACQUITTED (the
     same-valued `a_sum.inverse()` vs `two_space.inverse()` pins guard
     DIFFERENT laws — the concrete-method annotation vs the
     `SupportsInverse` bridge face; divergent red-surfaces proven);
     bare-class pins non-vacuous (widening any to `LinearOperator` reds);
     13-param signature within the file's fixture-free idiom. 1 cosmetic
     NIT folded (the docstring's "two-space throughout" over-claim
     re-worded).
   - **Docs (archivist fresh instance; Sphinx `-E -W` exit 0):** the
     stale `inverter`-parameter section (discrete_ordinates.rst
     ~:12242–12294) REWRITTEN → "Choosing the :math:`A^{-1}` realisation
     (the inverse-operator family)" + new anchor
     `choosing-inverse-realisation`; drivers take no `inverter` callable
     (pre-inverted `A_inv` through `SupportsSeededApply`;
     `KrylovAcceleration.inverter`→`preconditioner` with the
     category-mistake rationale — M approximates the FULL within-group
     inverse, not the step A⁻¹); "caller controls A⁻¹" survives as a TYPE
     choice over the four sibling KINDS; ERR-026/ERR-058 + the #195
     re-frame note preserved as history (ONE flagged deviation, correct
     per CR1: the note's false present-tense "the hook remains" clause
     re-spelled); the tracking `.. note::` DELETED (its task done); the
     Adams & Larsen present-tense straggler fixed; Dev-history P5 row at
     table top (pending-commit format; `SweepOperator.apply_transpose`
     as a code literal NOT `:meth:` — the L-002 forward-ref corollary:
     a deferred deliverable must not be a stale forward-ref).
   - **GitHub (Cardinal Rule 4):** #280 REDESIGNED onto `A.H.inverse()`
     (comment 4871342146): dead `CAP_SOLVE`/`_AdjointOperator.solve`
     spellings → the swap law `(A.H)⁻¹ = (A⁻¹).H` +
     `SweepOperator.apply_transpose` (EUCLIDEAN reverse-scan; the metric
     rides the landed `_AdjointOperator` wrapper — G_V⁺/G_W around
     apply_transpose — so NO metric code in the sweep); per-sibling
     adjoint-axis map (Green free-with-consumer, MatrixInverse
     `lu_solve(trans=1)` free-with-consumer, generic defers); the
     predicate-honesty fork ((a) construct-and-ask vs (b) land
     predicate+method together with apply_transpose — rec (b) with (a)'s
     spelling, user decides at that carve's checkpoint); gates 5.1–5.3 +
     M-ADJ-swap/M-ADJ-metric + the double-swap static pin INLINED
     self-contained (no `.claude/plans` dependency). #226 resolution
     comment posted (4871409512) — the step-3 docs-staleness tracker
     CLOSED.
   - **NOT built (deliberate):** no `_AdjointOperator.inverse()` (that IS
     #280's deliverable — implementing it now would ship a method whose
     every production path raises); no xfail gate-5.x file (an
     AttributeError-xfail verifies nothing — strict=False hides every
     failure mode; the gate specs live self-contained in the #280 comment
     instead). The sum-inverse fork 4c needed NO action — resolved at
     steps 4/5 (Sum→Green, `(L+C)` shadows→Sweep, §16 falsifier 7).

**Keystone gate threading all six:** `inv.apply(A.apply(x)) == x` for every
`(A, A.inverse())` pair (REDs the confused resolvent, passes the clean inverse).

---

## 13. Name-earning invariants — a name is a promise, backed by a test

**The bar (user, 2026-07-01):** a subclass name is earned only by a
**distinguishing invariant** — a property a *bare* inverse does not automatically
have. Round-trip alone earns only "InverseOperator"; anything fancier must be
tested to hold.

**Universal contract** (all inverse operators; necessary, NOT name-earning):
- **I1 round-trip:** `Ainv.apply(A.apply(x)) ≈ x` and `A.apply(Ainv.apply(q)) ≈ q`.
- **I2 functoriality:** `(αA)⁻¹=(1/α)A⁻¹`, `(AB)⁻¹=B⁻¹A⁻¹`, `(A⁻¹)⁻¹=A`.
- **I3 dagger law:** `A.H.inverse() == A.inverse().H` (#280).

**The family is keyed by WHICH mathematical object the inverse is — NOT by
algorithm** (the invariants prove realization is orthogonal: a dense Green and a
Krylov Green are the SAME Green operator):

| Name | Is the inverse of | Distinguishing invariant (the earning test) |
|---|---|---|
| `SweepOperator` | triangular `(L+C)` | **S-direct**: seed-INDEPENDENT result + MACHINE-precision residual (few·ULP, not iteration tol). A Krylov solve mislabeled "sweep" REDs this. |
| `GreenOperator` | full `(L+C−S)` (= `Resolvent(0)`) | **G-Neumann**: `G = Σ_k ((L+C)⁻¹S)^k (L+C)⁻¹` — the multiple-scattering expansion, which a generic `A⁻¹` has no splitting to satisfy; equals converged SI to tol + truncated partial sums. **G-reciprocity**: `⟨φ₂, Gφ₁⟩ = ⟨G.H φ₂, φ₁⟩` with `G.H == A.H.inverse()` (free NOW on the Krylov branch). **G-kernel**: `G.apply(δ_j)` = column j of the Green's matrix (unit-point-source flux). *G-positivity is CONDITIONAL* (scheme-gated: holds for step/characteristic, knowingly violated by DD in thick cells — document domain of validity, test under SC only). |
| `A.resolvent(z)` — a FACTORY name, no standalone class (§17 A3) | shifted `(A−zI)` | **R-identity** (first resolvent identity): `R(z)−R(w) = (z−w)·R(z)·R(w)` — a BINARY, family-grain law, unfakeable by a single inverse — so no single INSTANCE can earn the name; each `A.resolvent(z)` returns the structure-keyed inverse of the shift (Sweep/Green/Matrix), and the FACTORY runs the family gate. **R-continuity:** `R(0) == GreenOperator`. |
| `MatrixInverseOperator` | structureless small `A` (CP `[P]`, 0-D energy) | **M-materialise**: `Ainv.as_matrix() @ A.as_matrix() ≈ I` both ways **at MACHINE·cond precision** — ~~a matrix-free inverse CANNOT satisfy it (no `as_matrix()`)~~ *(superseded at step 5, spec §27.A: `as_matrix` became a universal base method, so an iterative Green also satisfies the identity — to DRIVER tol; the PRECISION GRAIN is the distinguisher, gate-proven by a Green contrast)*. Plus M-direct (machine·cond residual, seed-independent — bit-identical under any seed). |

**Naming rulings this forces:**
- `GreenOperator` is **EARNED** (G-Neumann + G-reciprocity are real, testable,
  un-fakeable) — and `KrylovInverseOperator` is REJECTED as a sibling: "Krylov"
  names the orthogonal realization axis, not the object. (Flips the §11.5 lean.)
- `MatrixInverseOperator` over `DenseInverseOperator` — ruled by the CP density
  probe (§14): `[P]` is dense-by-construction so "Dense" is a redundant storage
  claim, and #56's interface-current variant will be genuinely sparse.
- "Resolvent" RECLAIMED: the current `_GaussSeidelResolvent`/`_MomentWindowed-
  Resolvent` satisfy NONE of R-identity (they are sweep configs, mislabeled);
  they dissolve (§7), freeing the name for the true shift-invert α-family — at
  FACTORY grain: `A.resolvent(z) = (A−zI).inverse()`, no standalone 4th class
  (§17 A3: R-identity is BINARY/family-grain; a single instance earns only
  "inverse of the shift").
- Testability by stage: Green + Resolvent invariants testable TODAY (Krylov-branch
  adjoint is free); Sweep gets S-direct now and reciprocity after #280.

---

## 14. The two empirical rulings (2026-07-01 — literature + probes; agents:
literature-researcher ad166f2, numerics-investigators ae5dada / a5bb4b1)

### 14a. Curvilinear S-direct: OUR MAL-FORMULATION, geometry-split — fix filed as #282

The user asked: "is curvilinear seed-dependence unavoidable, or a mal-formulation
of ours?" **Answer: mal-formulation (verdict (c)), and only the SPHERE.**

- **Literature** (local: Hébert §3.9.3/3.9.4; Bailey–Morel–Chang NSE 165; Lathrop
  NSE 134; Morel NSE 101): at μ=−1 the redistribution coefficient `(1−μ²)`
  VANISHES → the starting-direction equation is CLOSED (Hébert 3.432/BMC 17), the
  full system lower-block-triangular, solved direct one-pass in every published
  formulation. Morel–Montry = weights only (τ_m, BMC 42), no lag. **No published
  formulation lags the pole seed to the previous iterate.** Full inventory:
  `agent-memory/literature-researcher/curvilinear_sweep_directness_ruling.md`.
- **Probes** (operator-level, removal-form): slab AND **cylinder** are bit-exact
  seed-independent with ~7e-16 residual (cylinder's per-level α-dome telescopes —
  the seed cancels); **sphere** has 4.6e-2 seed-sensitivity and 5.2e5 relative
  residual — exact ONLY at the SI fixed point (fixed-point probe: 1e-11). The one
  lagged quantity: the ψ½ half-angle seed (`loss_representation/__init__.py:3197`,
  default `AngularEdgeExtrapolation`). The direct solver ALREADY EXISTS
  (`carlson_inward_sweep_from_source`, `psi_half_angle_seed.py:428` = Hébert
  3.434–3.435). The MMS ladder is Mode-7 blind (all ansätze ≤ linear-in-μ, where
  the extrapolation seed is exact).
- **Taxonomy consequence:** S-direct is the RIGHT promise; the sphere is FIXED to
  meet it (route (a): carry ψ(μ=−1) per level as explicit augmented state), not
  the name weakened. The S-direct gate lands green (slab/cyl) + xfail-until-#282
  (sphere). Fix also dissolves lesson L21 (curvilinear full-angular constraint)
  and unblocks #200 (real curvilinear Krylov preconditioner). Probes preserved in
  `derivations/diagnostics/diag_curvilinear_seed_sensitivity.py` +
  `diag_sphere_fixedpoint_consistency.py` (promotion targets per #282).

### 14b. CP `[P]` density: HYPOTHESIS CONFIRMED → `MatrixInverseOperator`

`[P]` is structurally dense by construction — **0 exact zeros** in every config
(3–25 regions × thin/thick τ 0.25–25 × slab/cyl/sph), entries decay ~e^(−τ_ij)
to ~1e-12 at τ=25 but never vanish/negate; reciprocity 1e-16. Dense `np.ndarray`
all-pairs assembly (`cp/solver.py:212–397`). Interface-current sparse variant =
open #56, unbuilt. **Bonus real bug found and filed as #283:** spherical `[P]`
violates row-sum conservation for τ≳3 (interior rows saturate at exactly 1/π;
white-BC rank-1 closure amplifies to ~4.5; refinement worsens it → normalization
defect, suspected missing solid-angle factor; suite blind because all fixtures
τ≤2). Probes: `derivations/diagnostics/diag_cp_density_*.py`; memo:
`agent-memory/numerics-investigator/cp_matrix_density_and_sphere_conservation.md`.

### The meta-lesson

The name-invariant discipline CAUGHT TWO REAL BUGS before any rename happened:
demanding "SweepOperator must test S-direct" exposed the sphere seed-lag (#282);
demanding "Dense must be a true storage claim" led to the probe that exposed the
CP conservation defect (#283). Names-as-promises is not aesthetics — it is a
bug-finding instrument.

---

## 15. Ground-truth appendix (code state the design was built against)

> Pointers, not conclusions — a re-evaluation should SPOT-CHECK these against the
> live tree (or re-dispatch explorer on the 6 axes), not trust them. The Nexus
> graph predates this branch (L22 hazard): rebuild Sphinx before graph queries.

**Branch state (2026-07-01):** `refactor/inverse-as-operator` @ `b52bdef`, 4
commits over `main@f903686`: `32314cc` ratchet 420→410 · `c1126c5` P2a predicates
+ Supports{Inverse,Adjoint} · `d36507b` P2b advertisers + faithfulness keystone
(ratchet →409) · `b52bdef` P2c SweepOperator + `InvertibleOperator.inverse()`.
Tree clean (only `.claude/*` policy-uncommitted). CLI pyright oracle **409**.
Green: faithfulness gates (42), keystone equivalence (7, M-EQ teeth proven),
S†/scattering canaries (132+1 xfail), regression sweep (1919).

**The 7 promise-vs-identity mismatches** (explorer ground-truth; these are the
friction points the design must dissolve — verify each still holds):
1. **`is_invertible` on NINE ops, `.inverse()` on 1 (§17-corrected count):**
   Multiplication (multiplication_operator.py:209), OperatorProduct/Scaled/
   Identity/Permutation/TensorProduct/Diagonal (operator.py 967/1038/1071/1259/
   1536/1870), Invertible (streaming.py :682 + :781 — the only `.inverse()`),
   PLUS the shim `_BoundBoundaryOperator` (geometry/boundary/_bound_compat.py:108)
   forwarding the realized face laws' answer.
2. **THREE coexisting invertibility surfaces:** `CAP_SOLVE` string + `.solve`
   method / `is_invertible` property / `SupportsInverse`+`.inverse()`. Only
   `InvertibleOperator` spans all three (same triple on the adjoint axis).
   Endpoint = ONE surface (property + operator-returning method), carve P4.
3. **`SweepOperator.apply(rhs, *, initial_guess=None)`** vs base
   `apply(x, /)` — a keyword-only OPTIONAL extension (LSP-compatible; §11.3 §38
   ruling). Confirm pyright Protocol conformance stays clean.
4. **`InvertibleOperator` is-a `OperatorSum` by type, fused operator by
   identity** (overrides `apply` to `loss_action(σ)`). Its `.inverse()` OVERRIDE
   (→Sweep) over the future generic `OperatorSum.inverse()` (→Green) IS the
   type-as-structure dispatch (§11.1) — the design leans on this MRO fact.
5. **Resolvents' vestigial `.apply`** (§7) — solver.py:386/565, documented
   "never called."
6. **TWO energy-block materialization routes:** `dense_per_material` (iso ops)
   vs `_as_dense` apply-to-basis (homogeneous/solver.py:109) — the iso docstring
   claims the consumer the other route actually serves. `as_matrix()` (§5) is
   the dedup target.
7. **TWO sweep front doors:** module fn `transport_sweep`
   (loss_representation/__init__.py:1938, the SNSolver consumer) vs
   `InvertibleOperator.solve` — both single-source through `default_for(mesh)
   .sweep` (mild front-door twin; unification candidate at §12 step 3).

**Key pointers:** substrate `loss_representation/__init__.py` (Protocol :229;
loss_action :253 / transpose :294 / sweep :241; `_OctantWalk` :717, no-boolean-
fork :738; 2-D transpose UNBUILT ScanMarch:1858) · operators `streaming.py`
(:520 InvertibleOperator, :732 apply, :781 inverse, :797 solve; loss_rep
ownership :471/:703) · `sweep_operator.py:44` · adjoint `operator.py` (:704
_AdjointOperator, :757 G-metric apply) · drivers `iteration.py` (SI :329, step
:518-536; Krylov :576, gmres :802) · resolvents `solver.py` (:338/:486; factory
:675/:600; direct-solve site :2347) · matrix precursors `homogeneous/solver.py`
(:109 _as_dense, :129, :195 direct_eigenvalue), `sweep_graph.py:66` (sparse-
triangular future).

**Evidence artifacts:** probes `derivations/diagnostics/diag_curvilinear_seed_
sensitivity.py`, `diag_sphere_fixedpoint_consistency.py`, `diag_cp_density_*.py`
(+2 harness-caveated characterization probes) · memos `agent-memory/literature-
researcher/curvilinear_sweep_directness_ruling.md`, `agent-memory/numerics-
investigator/cp_matrix_density_and_sphere_conservation.md`, `agent-memory/
cross-domain-attacker/operator_inverse_taxonomy_frames.md` · issues **#282**
(sphere seed) / **#283** (CP conservation) · carve plan
`operator_inverse_algebra_carve.md` §0.5–0.7 (P2a-c detail; P3+ superseded) ·
verification spec `issue_226_inverse_operator_verification.md` (§2.3
faithfulness keystone + §4 P4 retirement map + literal-string-test traps REMAIN
VALID; only its P3 consumer-migration map is superseded by §12 step 3).

---

## 16. Re-evaluation guide — how to attack this plan with a clean context

**Decision authority — what is USER-LOCKED (re-open only by flagging to the user)
vs DERIVED (attack freely):**
- LOCKED: both-axes frozenset retirement via structural predicates + operator-
  returning methods (AskUserQuestion, 2026-06-29); **step-6 Design "C+B"**
  (user, 2026-07-02: TypeGuard checked bridges + per-class methods +
  structural/value split on the inverse axis; solve pruned to native
  realizations — see the §12 step-6 annotation); surgical main-agent-direct
  mode (no method-implementer); **the names-need-invariants bar** ("to be named
  GreenOperator it must be tested to have Green-operator invariants", 2026-07-01);
  `SweepOperator` in its own file; commit only when asked.
- DERIVED (fair game): the two-layers ruling (§1); type-as-structure keying
  (§11.1); resolvent dissolution (§7); GreenOperator-wraps-drivers (§11.2);
  seed-kwarg-as-context (§11.3); `as_matrix()` design (§11.4); the §12 order;
  the §13 invariant set; #282 fix route (a).

**Load-bearing claims → falsifiers** (each is checkable; a failed check breaks
the design at that pillar):
1. *Two-layer split is real in code* → find an operator whose `apply` and
   `solve` do NOT share one `LossRepresentation` instance (streaming.py:471/:703).
2. *Resolvents = same kernel, different schedule* → find a `loss_action`
   difference between `_GaussSeidelResolvent._solve_scheduled` and the operator's
   own representation (solver.py:448–462 injects the SAME `_sweep_interior`).
3. *Round-trip discriminator REDs the G-S resolvent, passes a clean inverse* →
   run `inv.apply(A.apply(x)) == x` on both.
4. *SI migration is behavior-preserving pre-#282* → the new consumer shape
   (`L_inv.apply(rhs, initial_guess=ψ_prev)`) threads the seed EXACTLY as today's
   `_solve_with_seed`, so iterates are bit-identical INCLUDING the un-fixed
   sphere. Verify via the spec §3 pinning gates (SI contract, keff curvilinear).
5. *G-Neumann is constructible/testable* → on a small case, `Σ_k ((L+C)⁻¹S)^k
   (L+C)⁻¹ q` converges to the Krylov/dense `(L+C−S)⁻¹q`.
6. *Sphere is the ONLY S-direct violator* → re-run
   `diag_curvilinear_seed_sensitivity.py` (slab/cyl bit-exact, sphere not).
7. *Type-as-structure dispatch is sound* → `InvertibleOperator.inverse()`
   (→Sweep) must shadow the generic `OperatorSum.inverse()` (→Green) by MRO, and
   no consumer depends on the sum-typed view of `(L+C)`.

**KNOWN WEAK POINTS — attacked by the 2026-07-01 re-evaluation; statuses inline
(evidence + resolutions in §17):**
- **(W1) [RESOLVED — §17] The moment-emit config changes the CODOMAIN.** A moments-emitting
  "SweepOperator" maps FullField→MomentField — strictly it realizes
  `P_frame ∘ (L+C)⁻¹` (a FUSED COMPOSITION for memory), NOT the endomorphic
  `A⁻¹`. §7/§11.1 gloss this as "an emit config." Decide honestly: either the
  windowed object IS typed as the fused product (its identity = a composition,
  its fusion = §38 context), or it is a different operator — but it must not be
  called plain `A.inverse()`. This is the same one-promise principle applied to
  our own fix.
- **(W2) [RESOLVED — §17: reify `M`] The fold-config honesty.** `(L+C−B_lower).inverse()` is the conceptual
  identity of the G-S resolvent, but `B_lower` is realized as a SCHEDULE (octant-
  group reflect), not an algebra subtraction — no forward `(L+C−B_lower).apply`
  exists to round-trip against. What does the discriminator test for THIS object?
  (Candidate: round-trip against `(L+C−B).apply` restricted to the schedule's
  fixed point; or accept the weaker "converged-SI-equivalence" pin. Unresolved.)
- **(W3) [VERIFIED VIABLE — §17] GreenOperator module placement.** `OperatorSum.inverse()` (numerics.
  operator) returning a wrapper around SourceIteration/KrylovAcceleration
  (numerics.iteration) inverts the current import direction (iteration imports
  operator). Candidate: GreenOperator lives in iteration.py; OperatorSum.inverse()
  late-imports (the SweepOperator↔InvertibleOperator precedent). Deferred with
  §12 step 4 — but the cycle must be checked before step 4 is scheduled.
- **(W4) [HALF-RESOLVED — §17 A2 rules the partiality FRAMING; threshold value +
  dense-vs-sparse keying decided at §12 step 5] `as_matrix()` policy holes**
  (§11.4).
- **(W5) [CONFIRMED + a gate-hole found — §17 F4; step-3 prerequisite added to
  §12] S-direct on sphere pre-#282**: the gate ships strict-xfail; confirm no
  step of §12 SILENTLY depends on sphere S-direct before #282 lands (step 3's
  behavior-preservation claim #4 above is the guard).

**Re-evaluation protocol for the fresh session:** treat §13/§14 as CLAIMS to
falsify, not authority; the memos/probes are evidence to RE-RUN, not citations
to trust; spot-check §15 against the live tree (or re-dispatch explorer);
rebuild Sphinx before any Nexus query (stale-graph hazard). The user's north
star is the acceptance test: every object's method surface ≡ its mathematical
identity — hunt for any object in THIS design that still promises what it
cannot do (W1 and W2 are the candidates the authors already see).

**EXECUTED 2026-07-01 — results in §17. W1/W2 resolutions are folded inline
above (§3/§7/§8/§9/§11/§12/§13/§15 carry `§17`-tagged amendments).**

---

## 17. Re-evaluation verdict (2026-07-01, clean context — the §16 protocol executed)

Four parallel independent legs (explorer / numerics-investigator / qa /
cross-domain-attacker); every §16 falsifier RUN against the live tree @
`69ed531`; memos treated as claims. The Nexus graph was confirmed STALE
(provenance `feature/sn-adjoint-transport`) and skipped — everything below is
branch-read or branch-run.

### Falsifiers — all seven CONFIRMED

| # | claim | result |
|---|---|---|
| 1 | two-layer split real | CONFIRMED — one cached `LossRepresentation` per family (`streaming.py:490 → :713 → :758/:1037`), object-identity spy-pinned (`test_one_representation_instance.py`). Honest boundary: module-fn `transport_sweep` builds FRESH per call — deliberate, operator-free door (mismatch 7). |
| 2 | resolvents share the kernel | CONFIRMED — `_solve_scheduled` injects the operator's own `_sweep_interior` (`solver.py:461`); G-S vs Jacobi differ ONLY in schedule + reflect; windowed holds NO kernel (delegates `base.solve_moments`). G-S is never constructed for 1-D (`solver.py:360–363`). |
| 3 | round-trip discriminator has teeth | CONFIRMED (run) — clean `(L+C).inverse()`: **3.8e-16**; `_GaussSeidelResolvent` round-trip vs its (vestigial) forward: **2.667, O(1) RED** (defect = `(L+C−B_lower)⁻¹B_lower x`, boundary-sourced, streamed everywhere). No pre-existing gate ran the G-S side — net-new evidence. |
| 4 | SI migration behavior-preserving | CONFIRMED-WITH-HOLE — the keystone `inverse().apply ≡ solve` pins only the WRAPPER (bit-identical by construction); the migration rewires the LOOP surface (per-iterate seed threading), and the strong het-2G sphere catcher is `@slow`-DESELECTED under the canonical gate (simulated seed-drop |Δk|=3.46e-2 reddened only a fragile 1G margin). → step-3 prerequisite added in §12. |
| 5 | G-Neumann constructible | CONFIRMED (run) — partial sums via existing `LC.inverse()` ∘ `S.apply` decay geometrically (**ratio 0.4995 ≈ c=0.5**), agree with a structurally-independent dense LU of `(L+C−S)` at **4.2e-13** by k=40. `GreenOperator` is buildable from today's machinery. |
| 6 | sphere the only S-direct violator | CONFIRMED (re-run, not cited) — slab/cyl bit-exact (8.1e-16 / 6.9e-16); sphere seedΔ **4.57e-2**, cold residual 5.2e5, fixed-point-exact 1.5e-11. NEW: production blast radius is quadrature/mesh-gated (LS-S4 / 16-cell fixed-source: SI→NaN, Krylov→negative flux; GL-S16 / 40-cell physical) — documented on #282 (comment 2026-07-01). |
| 7 | type-as-structure dispatch sound | CONFIRMED — `InvertibleOperator(OperatorSum)` MRO fact true; ZERO sum-typed consumers (no isinstance/type-is dispatch; `.terms` doesn't exist — the sum is binary `.a`/`.b`; leaves reached via typed properties). Caveat → the §12-step-4 ordering ruling (`L+C` fuses, `C+L` doesn't). |

Branch-green claims verified under the canonical gate: faithfulness **42**
(= frame 19 + predicates 12 + survival 11 — exact), keystone equivalence **7**
(exact), scattering-adjoint sample 19 incl. the 1 known xfail — 87 passed /
1 xfailed / 0 failed.

### W1/W2 resolutions (folded inline above; full frames in
`agent-memory/cross-domain-attacker/operator_inverse_w1_w2_resolutions.md`)

- **W1 RESOLVED — deforestation of a composition through a coisometry.** The
  windowed object = `P @ A.inverse()` (`P` = analysis face ⊕ Id on the trace —
  block coisometry; FullField→MomentField; NOT invertible ⇒ no round-trip
  promise). Fusion stays SUBSTRATE (`_SweepEmit` moment mode); factory = the
  solver's `_maybe_window`. The `solve_moments(rhs, frame)` cross-reach (lives
  3×) retires into the ONE product. The honest contract FACTORS:
  fusion-correctness (windowed ≡ post-projection oracle, ≤4 ULP — the existing
  gate, renamed as deforestation) ∧ S-direct on the base `A⁻¹` ∧ coisometry
  `M∘R=I` on the frame. A "moment-proxy residual" `M[(L+C)Rφ]−Mq` is REJECTED —
  category-confused (`Rφ≠ψ`; it tests the different P_N-reduced system and
  would RED a perfect fused sweep).
- **W2 RESOLVED — regular matrix splitting, M reified.** `(L+C−B) = M − N`,
  `M = (L+C−B_lower)`, `N = B_upper`; B has no octant-diagonal so the split is
  exact and M stays triangular. REIFY `M.apply = (L+C).apply − B_lower.apply`
  (B is already a `LinearOperator`; mask to the schedule's strictly-lower octant
  pairs). Gate ranking: (1) machine round-trip on M — the runnable keystone;
  (2) fixed-point equivalence vs Jacobi on a DIAGONAL cubature (Mode-9-safe);
  (3) ρ(M⁻¹N)<1 — a convergence certificate, not correctness. The
  converged-SI-equivalence fallback is REJECTED: Mode-9-DEGENERATE (the fixed
  point is splitting-invariant — it cannot even distinguish G-S from Jacobi);
  two legs (qa + cross-domain) reached this independently. "Type it as a
  preconditioner" also REJECTED — a 3rd mislabel-by-driver-role after
  "resolvent"; the object is the EXACT inverse of M.

### Three grain-sharpenings (A1–A3, folded inline)

- **A1 (§11.3):** seed-as-context ⟺ S-direct is the SAME claim — the ruling
  inherits the geometry split; sphere pre-#282 `.apply(q, seed)` is a relaxation
  step (identity-bearing kwarg), an xfail-until-#282 defect, not a documented
  hybrid.
- **A2 (§11.4):** `MatrixTooLarge` = resource effect on a TOTAL functor, not a
  restriction idempotent; no `is_materializable` predicate.
- **A3 (§13/§3):** R-identity is BINARY/family-grain ⇒ "resolvent" is a FACTORY
  name (`A.resolvent(z) = (A−zI).inverse()`); the standalone 4th class is
  dropped. §13 had conflated instance-grain and family-grain invariants.

### New findings the plan didn't contain

- **NINTH `is_invertible` advertiser:** `_BoundBoundaryOperator`
  (`geometry/boundary/_bound_compat.py:108`) — §12 step 1 decides
  shim-forwarding.
- **Consumer audit:** `.inverse()` has ZERO production callers; `SupportsInverse`
  zero consumers anywhere (deliberately not `runtime_checkable` — an
  annotation-only surface with zero annotations); `is_invertible` zero
  production readers outside composer propagation; SIX advertisers' `.solve`
  also production-unconsumed; SN Krylov bypasses the one generic
  inverse-consuming seam with an identity-lambda preconditioner (#200). BUT
  `operator.py:388` already PUBLISHES ".inverse()" inside `is_invertible`'s
  docstring — so §12 step 1 is contract-completion OR contract-retirement, per
  type (reframed above).
- **Order-dependent strategy trap** (post-step-4): `L+C` → Sweep, `C+L` → Green
  — same math, different algorithm by operand spelling; §12 step 4 carries the
  ruling requirement.
- **Fast-gate hole (F4):** the sphere seed-threading catcher is slow-deselected
  — §12 step-3 prerequisite (Mode-11 path-spy on `initial_guess` threading).
- **NEW SMELL candidate (1st sighting, cross-domain):** "codomain-changing
  emit/config" — an output-mode config presented as a config of ONE morphism
  while silently changing the codomain ⇒ it is a DIFFERENT morphism (a
  composition), not a config. TELL: a solve/apply variant taking a FOREIGN
  operator's projection as a parameter. Held for a 2nd sighting.
- **Scope note on §6's "free" Krylov adjoint:** `GreenOperator(A.H)` is free at
  the OBJECT level, but the multi-D `loss_action_transpose` still raises
  (`ScanMarch:1858`, re-confirmed) and a preconditioned adjoint solve needs
  `sweep_transpose` — so the practical adjoint-inverse is 1-D-capable today,
  multi-D gated by #280. Do not oversell "free NOW".

### Grander-scheme verdict

**The taxonomy IS the path forward.** It is the inverse-side completion of the
operator-algebra architecture the codebase has been converging on for months
(typed fields → frames → K_iso → adjoint-by-transpose → eigenvalue engines),
and the names-need-invariants discipline has now caught two real physics bugs
(#282 sphere seed-lag, #283 CP spherical conservation) BEFORE any rename — it
is a bug-finding instrument, not aesthetics. The serious rival — the PETSc/KSP
shape (operators forward-only; ALL inversion lives in solver drivers; never
reify `A⁻¹`) — was weighed and loses: it forfeits the composability the physics
is written in (`K = A_loss.inverse() @ F`, preconditioner algebra `M⁻¹A`,
dagger functoriality `(A.H).inverse() == (A.inverse()).H` = #280's content),
re-introduces twin paths at driver level, and the reified inverse ALREADY
round-trips at machine precision behind green gates. Its one valid caution is
already this plan's own rule: defer each inverse object until a consumer
exists (GreenOperator stays deferred; step 1 is now deliver-vs-retire, not
blanket delivery).
