# G2 — verification plan for `orpheus/geometry/transformation.py`

> **Status: PLAN — and G1 LANDED MID-PLAN.** The brief said
> `orpheus/geometry/transformation.py` did not exist; it now does (763 lines,
> **untracked**, written during this session). **§0 is the reconciliation and
> should be read first**: it re-scores every API finding against the shipped
> code (10 discharged, 5 open), adds three gates the landed design earns, and
> records that a **representative subset — 13 gate families, including every ⭐
> row (F3, C2, C3, C5, D2, D3, O2, O4, O7, O8, H1, A5, A6) — was smoke-run
> against the real API and every one passes**, so the matrix is writable as
> specified. §4 is written API-first and remains the specification.
>
> Campaign plan of record: `.claude/plans/geometric_transformation_consolidation.md`
> (rulings **R1, R2, R4, R5**; sequencing row **G2**).
> Prior art NOT to be re-pinned: `tests/numerics/test_symmetry.py` (105 tests) — see §2.
>
> `[M]` marks a fact MEASURED during this planning session (probe scripts in
> `/Users/rodrigo/.claude/jobs/c30e4f25/tmp/`). Everything tolerance-shaped in
> this plan is measured, not guessed.

---

## 0. ⛔ RECONCILIATION — G1 LANDED DURING THIS PLAN

**`orpheus/geometry/transformation.py` now exists** (763 lines, **untracked**,
mtime 2026-08-03 04:30 — created *while this plan was being written*). The brief
said it did not. Everything below has been re-checked against the shipped code,
and a **representative subset of §4 was smoke-run against the real API**
(probe `…/tmp/probe_real.py`, 13 gate families covering every ⭐ row): **every
one passes**, so the matrix is writable as specified. Two consequences:

* **§8 is re-scored.** Ten of the fifteen API problems are **DISCHARGED by the
  landed design** — several more elegantly than proposed. Five remain **OPEN**
  and are the actionable output.
* **Two tolerances tighten and one risk is retired** — see the "reconciled"
  column of §10.

### What the landed API actually is (and it is better than the sketch)

There is **no `Reflection` type and no `Rotation` type**. There is one frozen
`RigidMotion(linear, translation)` whose orthogonality invariant lives in
`__post_init__` (so every same-type-producing operation re-establishes it by
routing through construction — coding-elegance Pattern 4), with **named
constructors** `identity / translation_by / inversion / reflection /
rotation_from_circle_point / rotation / rotation_about_axis /
signed_permutation`, and `determinant` **computed, never declared**. That
single decision dissolves API problems A2 (where does `−I` live?) and A3 (the
normal's gauge freedom) outright: there is no type whose *label* could drift
from its *content*, and `−I` is just `inversion(d)` whose `determinant` reports
`(−1)^d` by measurement.

Three things the shipped module has that this plan did **not** propose, and
each earns a gate:

| shipped | why it matters | gate to ADD (§4-X) |
|---|---|---|
| `on_points` vs `on_directions` — two spelled actions | a direction has no position, so `t` must not act on it; applying the affine map to a direction silently denormalises unit vectors. This is exactly the BC law `ψ_in(x,Ω) = ψ_out(g⁻¹x, Q_g⁻¹Ω)` | **X1**: `on_points(x) − on_directions(x) == t` exactly; `‖on_directions(v)‖ = ‖v‖`; and `on_directions` is invariant under `seated_at(c)` for every `c` |
| `rotation_from_circle_point(plane, point)` is the **primitive**; `rotation(angle)` the convenience | this is API problem **A8** already solved — an exact `(cos, sin)` from a root-of-unity construction survives the constructor | **X2**: `rotation_from_circle_point(point=(0,1))` yields a matrix with **exact** zeros, while `rotation(angle=π/2)` does not |
| `element_order() -> int \| None` returns the **order**, and `None` for a glide/screw | the plan asked only "is the order n?"; returning the order is strictly stronger, and the `None` arm is a law no one had named — an element whose *linear* part has finite order but which translates off its fixed subspace generates an **infinite** cyclic group | **X3**: `reflection(normal=e_x) @ translation_by(e_y)` (a glide) has `element_order() is None` while its linear part has order 2 |

### §8 re-scored against the shipped code

| # | verdict |
|---|---|
| **A1** d=1 rotation | ⚠ **PARTIAL** — it *does* raise, but with the message *"the rotation plane must be given by an ORTHONORMAL pair"*, which blames the caller's plane instead of naming `SO(1) = {e}`. A gate pinning that message pins the wrong reason |
| **A2** `−I` has no home | ✅ **DISCHARGED** — one type, `inversion(d)`, `determinant` computed |
| **A3** normal gauge freedom | ✅ **DISCHARGED** — `_exact_key` hashes the realized `(Q,t)` with `+0.0` canonicalisation, routed through one function for both dunders |
| **A4** closure cap + named error | ⚠ **PARTIAL** — `_MAX_GROUP_ORDER = 2000` cap present and `[M]` both infinite cases raise; but it is a **bare `ValueError`**, so gate A6's `pytest.raises` cannot be specific |
| **A5** explicit tolerance | ✅ **DISCHARGED** — `permutes(points, *, atol)` keyword-only and **required**; `preserves(..., atol, weight_atol)` — two windows, and the `_NODE_WINDOW_FACTOR = 100` fudge did not survive |
| **A6** failure reason collapsed | ❌ **OPEN** — still `-> np.ndarray \| None`. §5.3's input-side isolation is therefore **mandatory**, not optional |
| **A7** `preserves(measure)` vs numpy-only | ✅ **DISCHARGED** — takes `(points, weights)` arrays; the layering tension is resolved in favour of R3 |
| **A8** exact rational turn | ✅ **DISCHARGED, and better** — `rotation_from_circle_point` is the primitive |
| **A9** seat-from-data | ✅ **DISCHARGED** — `seated_at(centre)`, delegating to `conjugated_by` |
| **A10** composition order | ✅ **DISCHARGED** — `a @ b` = "apply `b` first", documented; `[M]` gate O2 confirms it (0/144 violations, wrong order 102/144) |
| **A11** one action entry point | ✅ **DISCHARGED, and better** — *two* named actions (see above) |
| **A12** `det`/`is_proper`/`dim` on the element | ✅ **DISCHARGED** — plus `is_linear`, `fixed_subspace_dimension`, `element_order`, `fixes` |
| **A13** mixing dimensions raises | ✅ **DISCHARGED** — in `__matmul__` and in `close_group` |
| **A14** the embedding `ι: E(d) → E(d+1)` | ❌ **OPEN — the most consequential gap.** `[M]` `grep` for `embed`/`lift` on `RigidMotion`: **zero methods.** R2's claim that dimension genericity *dissolves* the 1-D/3-D arm split needs this operation to be **of the type**; without it, campaign step **G3** ("delete the 1-D/3-D arm split") will pad with `np.eye(3)` at the call sites and the split will have *moved*, not gone. Gate **H2** is unwritable until it exists |
| **A15** `Permutation` type | ❌ **OPEN** — bare `np.ndarray`. Lower priority than A14; the ERR-073 guard is present and `[M]` correct |

### ⭐ NEW finding — a gate must not assert TIGHTER than the type's own invariant

`__post_init__` accepts any `Q` with `max|QᵀQ − I| ≤ _ORTHOGONALITY_ATOL = 1e-10`.
So a `RigidMotion` built from a raw matrix carrying a **1e-11 shear** is a legal
value of the type. Gate **B2** as drafted asserts `atol=1e-14` — which is
**tighter than the type promises**, i.e. it would assert a property `RigidMotion`
does not guarantee, and would red on a legal instance.

Split B2 into two honest claims:

* **B2a — the type's invariant**: for *any* constructed element, including one
  built from a raw matrix and every product, `max|QᵀQ − I| ≤ _ORTHOGONALITY_ATOL`.
  Falsifier: hand a `Q` off by `2e-10` → must raise (a real production entry
  point, not a self-satisfied `raises`); hand one off by `5e-11` → must be
  accepted. This is the *only* row that pins the construction guard.
* **B2b — the named constructors are far better than the invariant**: elements
  from `reflection / rotation / rotation_about_axis / signed_permutation /
  inversion` satisfy `max|QᵀQ − I| ≤ 1e-14` (and `signed_permutation` **exactly
  0**, since its entries are integers). This is a stronger, separately-earned
  claim *about the constructors*, and it is the one that would catch a
  constructor regression.

Conflating them would leave the constructors ungated at their real quality while
falsely claiming the type guarantees 1e-14.

### Standing hazard for the implementation session

The module is **untracked**. Mutation-testing it MUST be done by `monkeypatch`
in-process or by copying to `$CLAUDE_JOB_DIR/tmp`; a `git checkout` / `restore`
/ `stash` on any path would destroy it *and* the two tracked skill files that
carry irrecoverable uncommitted state.

---

## 1. Claim layer, pillar, structural independence — the §1.5 gate

Applied before the matrix was drafted, per `vv-principles`.

**Claim layer.** Every gate in this plan is a **structural / algebraic-identity
claim** about a mathematical object — the lowest layer of the transport claim
ladder, below even a convergence-order claim (there is no discretisation here
and no `h`). None of these gates is a flux-shape or eigenvalue claim, so the
MMS-cannot-prove-eigenvalues restriction never binds and the 1-group-degeneracy
Cardinal Rule does not apply (there is no energy variable in this module).

**Pillar.** Every gate is **closed-form** — an exact identity of group theory
and linear algebra, verified either symbolically (SymPy, exact) or numerically
against an independently-constructed closed form. **No MMS row exists** (there
is no differential operator to manufacture a solution for) and **no
semi-analytical row exists** (nothing reduces to a quadrature). A plan for this
module that reached for MMS would be a category error.

**Structural independence — the chain of trust.** The whole point of ordering
G2 before G7 is that after G2 the primitive may be shared by a SUT and its
oracle without contamination, because independence lives in the INPUT ASSEMBLY.
That licence is only earned if G2 itself terminates in ground that is *outside
ORPHEUS*. The permitted references, in descending strength:

| # | reference | why it is structurally independent | covers |
|---|---|---|---|
| **P1** | **SymPy symbolic identity** under an explicit unit parameterisation (`n̂ = (sinθcosφ, sinθsinφ, cosθ)`) | exact algebra, no floats, no ORPHEUS code in the chain | `QᵀQ = I`, `det`, Cartan–Dieudonné in ℝ², Rodrigues `trace = 1 + 2cos t` |
| **P2** | **`scipy.spatial.transform.Rotation.from_rotvec`** | an externally maintained, independently tested implementation; quaternion-based, so a *different algorithm* from Rodrigues' 3-term expansion | ℝ³ rotations only |
| **P3** | **`scipy.linalg.expm(θ(vuᵀ − uvᵀ))`** | the Lie-theoretic definition `R = exp(θJ)` evaluated by scaling-and-squaring Padé — structurally different from every closed form, and **dimension-generic** | rotations at d = 2, 3, 4 (`[M]` agrees with the closed form to **1.7e-13**) |
| **P4** | **closed forms written independently in the test** — `[[c,−s],[s,c]]`; `Σ_i e_ie_iᵀ − n̂n̂ᵀ` (projector form of the Householder); `arange(n)[::-1]` | the test writes the matrix from the textbook definition, never calls the SUT to build it | reflections/rotations at every d, reversal permutations |
| **P5** | **published group-theory tables** — group orders (`\|C_n\| = n`, `\|D_nh\| = 4n`, `\|O_h\| = 48`, `\|I_h\| = 120`, `\|C_4v\| = 8`), and the **cube's Wyckoff site symmetries** (`Hamermesh 1962 §9.4`; Int. Tables Vol. A) | a number in a book | closure orders, orbit–stabiliser absolute values |
| **P6** | **exact integer arithmetic** — a permutation's own bijectivity, `π(gh) = π(g)∘π(h)`, `ord(g^k) = n/gcd(k,n)` | no tolerance can be involved | the whole `permutes` family |

**FORBIDDEN as a G2 reference** (this is the load-bearing prohibition):
`orpheus.numerics.symmetry._rotation_z` / `_reflections` / `_vertical_mirrors`
/ `_octahedral_ops` / `_icosahedral_ops` / `_rotation_about_axis` /
`_close_group` / `_orbit_closure`, and anything downstream of them
(`Quadrature`, `DiscreteMeasure`, `roots_of_unity`). G2 answers "is the new
primitive *right*", not "does it reproduce what we had". The
bit-identity-to-`symmetry.py` comparison is a **separate, later gate (G3)** with
a different purpose (no-op migration proof) and it must not be confused with,
or substituted for, a G2 row.

---

## 2. What is ALREADY pinned — do not re-pin

`tests/numerics/test_symmetry.py` (105 tests, all `@pytest.mark.foundation`).
The overlap-and-difference analysis, so G2 adds only new information:

| existing test | what it pins | why G2 still needs its own row |
|---|---|---|
| `test_realized_group_orders_are_textbook` | \|Mirror\|=2, \|C_n\|=n, \|D_nh\|=4n, \|O_h\|=48, \|I_h\|=120 via `_close_group` | pins the **3-D tag realizations**; G2 must pin the **new dimension-generic closure** at d=1,2,3 and pin the **infinite** case (which no existing test does) |
| `test_orbit_stabilizer_theorem_holds_on_every_certificate` | \|orbit\|·\|Stab\| = \|G\| on shipped quadrature rules | **self-consistency only** — \|G\| comes from the same machinery. It is **blind to a uniformly-rotated realization** (its own sibling docstring at `symmetry.py:1118` says so). G2's row pins the **absolute Wyckoff orbit types** against a published table |
| `test_duplicated_node_breaks_oh_invariance` (`catches ERR-073`) | the bijectivity guard, on `O_h` in 3-D | G2 must re-pin at the **new** call site and, crucially, **isolate which guard fired** (see §5.3) — the existing test asserts only `is None`, a 3-way OR |
| `test_vertical_mirror_planes_follow_the_dnh_setting` | the `D_nh` plane *convention* (`kπ/n` vs `kπ/n + π/2`) | a CONVENTION claim about `symmetry.py`'s realization — not a law of the element type. Stays where it is |
| `test_computed_containment_obeys_the_order_relation_laws`, `test_invariance_is_downward_closed`, the ERR-072 family, `singular_set` family | the **lattice / predicate** layer | entirely above `transformation.py`. Untouched by G2 |
| `test_a_mirror_is_improper_...`, `test_o3_contains_...` | `det` at the **tag** level | G2 pins `det` at the **element** level, dimension-generically |

**Net:** G2 is not a re-run of the symmetry suite. Seven of its gate families
(homomorphism of the induced permutation, fixed-point-set pointwise, affine
conjugation, seating single-source, Cartan–Dieudonné, order-exactly-n,
dimension genericity at d=1) have **no counterpart anywhere in the tree today**.

---

## 3. Test-file placement and tagging

* File: `tests/geometry/test_transformation.py` (new). The module lives in
  `orpheus/geometry/` per **R3**, so the test lives under `tests/geometry/`.
* Marker: **`@pytest.mark.foundation`** on every gate, **no `verifies(...)`** —
  these are software/mathematical invariants with no theory-page `:label:`
  today. This matches `tests/numerics/test_symmetry.py` exactly and the
  project's `verifies ⊥ level` doctrine.
* **Recommendation (cheap, do it in G1):** mint `:label:`s for the five
  defining equations — seating `(Q,(I−Q)c)`, composition
  `(Q₁Q₂, Q₁t₂+t₁)`, inverse `(Qᵀ, −Qᵀt)`, Householder `I − 2n̂n̂ᵀ`,
  Cartan–Dieudonné — on a theory page. Then the corresponding gates carry
  `verifies(...)` and the provenance chain (`nexus provenance_chain`) exists.
  Without labels the module is permanent `orphan_code` in the V&V audit.
* Runtime target: the whole file **< 10 s** (≈8 s of it SymPy, §6). Every gate is a small matrix
  identity; the only expensive rows are the SymPy ones (§6, budgeted
  separately) and the `O_h` orbit rows (48 ops × ≤ 12 points).
* Canonical invocation `python -O -m pytest`, **serial**. Because the file is
  *collected*, pytest's assertion rewriter protects bare `assert` (Mode 8's
  measured scope) — but §7 still requires `np.testing.assert_*` for the
  numeric rows, for the error messages, not for `-O` safety.

---

## 4. The gate matrix — 42 gates

Notation: `σ_n̂` = reflection with unit normal `n̂`; `R(P,θ)` = rotation in plane
`P` by `θ`; `g = (Q,t)` acts as `x ↦ Qx + t`; `e` = identity; `π_g` = the
permutation `g` induces on a point set. **d** is the dimension(s) the row must
be parametrized over. All rows are **closed-form pillar**, `@pytest.mark.foundation`.

> **Gate IDs vs campaign steps.** Gate IDs are `A1…A6, B1…B4, C1…C5, D1…D5,
> E1…E3, F1…F4, O1…O8, H1…H3`. The orbit/action group is lettered **O**, not G,
> so that no gate ID collides with a campaign STEP (`G1`…`G8` in
> `.claude/plans/geometric_transformation_consolidation.md` §6). Any bare
> `G<digit>` in this document is a campaign step.

### A — group axioms of the rigid-motion algebra

| ID | LAW | REFERENCE (pure math) | FALSIFIER (mutation → RED) | lvl | d |
|---|---|---|---|---|---|
| **A1** | **Associativity, two tiers.** `(g∘h)∘k = g∘(h∘k)` — **bit-identical** on the signed-permutation subgroup; `atol=1e-14` for general elements | matrix-product associativity (exact in ℝ); the two tiers are *measured*, not assumed: `[M]` bit-exact **500/500** for signed permutations, **0/500** for general rotations (‖ΔQ‖≤2.2e-16, ‖Δt‖≤1.8e-15) | compose as `(Q₁Q₂, Q₂t₁+t₂)` (operand swap in the translation rule) — breaks associativity by O(1). Nothing else in the plan pins the *order* inside the translation rule except A5 and O2 | L0 | 1,2,3 |
| **A2** | **Identity, both sides, BIT-EXACT.** `e∘g = g∘e = g` with `np.array_equal` on both `Q` and `t` | the group axiom | make `e` be `(I, ε)` for any `ε ≠ 0`, or implement `e` as `RigidMotion(np.eye(d), np.zeros(d))` with a wrong `d` | L0 | 1,2,3 |
| **A3** | **Inverse, BOTH sides.** `g∘g⁻¹ = g⁻¹∘g = e` to `atol=1e-14`; and `g⁻¹ = (Qᵀ, −Qᵀt)` compared against `np.linalg.inv` of the `(d+1)×(d+1)` homogeneous matrix | the group axiom + the closed form of the inverse of a block-triangular homogeneous matrix | drop the conjugation in the inverse (`t ↦ −t` instead of `−Qᵀt`). **`g∘g⁻¹` is BLIND to this when `Q` is symmetric** (every reflection!) — only the *general rotation* draws catch it. `[M]` both sides residual ≈1.1–1.8e-15, **0/500 bit-exact** ⟹ do NOT gate at `array_equal` | L0 | 1,2,3 |
| **A4** | **Inverse is an ANTI-homomorphism.** `(g∘h)⁻¹ = h⁻¹∘g⁻¹` | group theory | write `g⁻¹∘h⁻¹`. `[M]` the swapped order differs on **500/500** random draws — non-vacuous *provided the elements do not commute* (see §7-V3) | L0 | 2,3 |
| **A5** | **Closure of the generated group has the TEXTBOOK ORDER.** `\|⟨σ_x,σ_y,σ_z⟩\| = 8`, `\|⟨R_z(2π/4), σ_x⟩\| = 8 (C_4v)`, `\|⟨σ⟩\| = 2`, `\|⟨R(2π/n)⟩\| = n`, cube-generator set → 48 | published group orders (**P5**, Hamermesh 1962 §9.4); `[M]` all four reproduced by an independent closure written in the probe | seed the closure with `ops` but never add `ops` to the element list (the exact seeding bug this planning session hit) → orders collapse to 1. Or drop the `-0.0` canonicalisation in the key → orders inflate | L0 | 1,2,3 |
| **A6** | **A NON-finite generating set RAISES.** `⟨translation(v)⟩` and `⟨σ_{n̂,d=0}, σ_{n̂,d=1}⟩` (two PARALLEL mirrors) must raise a named error, not spin | `[M]` both exceed a 300-element cap; two parallel mirrors generate the infinite dihedral group (textbook) | remove the cap → the test hangs (so give it a `pytest.mark.timeout` or an explicit cap argument). **This row is the reason `close()` must own the cap** | L0 | 1,3 |

### B — orthogonality and determinant

| ID | LAW | REFERENCE | FALSIFIER | lvl | d |
|---|---|---|---|---|---|
| **B1** | **SYMBOLIC `QᵀQ = I` and `det`.** For `n̂ = (sinθcosφ, sinθsinφ, cosθ)`: `simplify(QᵀQ − I) == 0` and `det Q == −1`. For Rodrigues with a symbolic axis+angle: `QᵀQ − I == 0`, `det == 1`, `trace == 1 + 2cos t` | **P1** SymPy, exact — no floats anywhere in the chain | any algebraic slip in the constructor. `[M]` all five identities verified, timings 0.1 / 0.5 / 0.8 / 5.0 s | L0 | 2,3 |
| **B2a** | **THE TYPE'S INVARIANT.** *Any* constructed element — including one built from a raw matrix, and every pairwise product — satisfies `max\|QᵀQ − I\| ≤ _ORTHOGONALITY_ATOL` (shipped: `1e-10`); a matrix off by `2e-10` is **REJECTED at construction**, one off by `5e-11` is accepted | the definition of `O(n)` + the type's own declared invariant | remove the `__post_init__` guard. **The `raises` leg calls the production constructor** (§7-V8). ⚠ Do NOT assert `1e-14` here — that is tighter than the type promises (§0) | L0 | 1,2,3 |
| **B2b** | **THE CONSTRUCTORS ARE FAR BETTER THAN THE INVARIANT.** Elements from `reflection / rotation / rotation_about_axis / inversion` satisfy `max\|QᵀQ − I\| ≤ 1e-14`; `signed_permutation` satisfies it **exactly 0** (integer entries) | orthogonality; `[M]` verified on the shipped constructors | build the Householder as `I − n̂n̂ᵀ` (missing the 2), or skip the normalisation of `n̂`. This is the row that catches a constructor regression — B2a would not, because a `1e-11` shear is a *legal* element | L0 | 1,2,3 |
| **B3** | **`det` VALUES.** reflection `−1`; rotation `+1`; translation `+1`; **inversion `(−1)^d`** | closed form; `[M]` `det(−I)` = −1 / +1 / −1 at d = 1/2/3 | swap the reflection/rotation `det` claim. **d=2 is BLIND to an inversion-vs-rotation-by-π confusion** (both `det=+1`) — the row is only non-vacuous because d=1 and d=3 are present | L0 | 1,2,3 |
| **B4** | **`det` is a HOMOMORPHISM** `det(g∘h) = det(g)·det(h)`, and **seating preserves the linear part BIT-EXACTLY** (`seat(Q,c).Q is/== Q` under `array_equal`) | `det: O(n) → {±1}` is a group homomorphism | let `seat` re-derive `Q` (e.g. re-orthonormalise) → the bit-exactness leg reds. The homomorphism leg is **vacuous if only proper elements are drawn** (`1 = 1·1`) — the row MUST mix proper and improper (§7-V4) | L0 | 1,2,3 |

### C — element order, and the fixed set

| ID | LAW | REFERENCE | FALSIFIER | lvl | d |
|---|---|---|---|---|---|
| **C1** | **A reflection is an involution.** `σ² = e`, linear AND seated. Bit-exact for a COORDINATE normal, `atol=1e-15` for a general one | group theory; `[M]` `σ²==I` bit-exact for `n̂=e_x`, **0/300** for random normals ⟹ two tiers | flip a sign in the Householder. **DECLARED BLINDNESS: this gate cannot see the offset** — see C3 and §5.1 | L0 | 1,2,3 |
| **C2** | **`ord(R(2π/n)) = n` EXACTLY** — `Q^n = e` AND `Q^k ≠ e` for `0 < k < n`; and `ord(R(2πk/n)) = n/gcd(k,n)` | elementary group theory; `[M]` the separation is **6 orders of magnitude** (`‖Q^n−I‖ ≈ 1e-15` vs `min_{k<n}‖Q^k−I‖ ≥ 1.6e-1` at n=24) — the gate is NOT tolerance-fragile | generate `C_n` from `R(2·2π/n)`. `[M]` caught at **even** n (order becomes n/2) and **NOT at odd** n (order stays n) ⟹ the row MUST include even n (§7-V2) | L0 | 2,3 |
| **C3** | ⭐ **The element fixes its named fixed set POINTWISE.** `σ_{n̂,d}` fixes every point of `{x : n̂·x = d}` to `atol=1e-14`; `R` seated at `c` fixes `c`; a translation fixes nothing | the definition of the objects; `[M]` residual 0.0 for the correct `t = 2d n̂` | **the offset factor**: `t = d n̂`, `t = 4d n̂`, `t = −2d n̂`. `[M]` **all three are still involutions** (C1 green) and all three fail C3 by 0.37 / 0.74 / 1.48. **C3 is the ONLY catcher in the plan for the affine offset.** See §5.1 | L0 | 1,2,3 |
| **C4** | ⭐ **Fixed-set DIMENSION — ruling R2's structural claim.** `dim Fix(σ) = d−1`; `dim Fix(R) = d−2`; computed as `d − matrix_rank(Q − I)` | linear algebra; `[M]` verified at d=1,2,3 (reflection 0/1/2; rotation −/0/1) | parameterise the rotation by its *axis* rather than its *plane* at d=2 (there is no axis) → construction fails; or build the reflection from `n̂n̂ᵀ` instead of `I − 2n̂n̂ᵀ` → `dim Fix` becomes 0. **Use `matrix_rank` with the DEFAULT relative tol, never eigenvalue-`isclose`** — `[M]` `isclose` misreports `dim Fix = 3` for a rotation at θ = 1e-9 (§7-V6) | L0 | 1,2,3 |
| **C5** | ⭐ **Cartan–Dieudonné, sharp form.** The minimal number of reflections is `k = d − dim Fix(Q)`, and `det Q = (−1)^k` | Cartan–Dieudonné theorem; `[M]` verified on 8 cases incl. `−I` at d=1,2,3, a rotoreflection, and the identity (k=0) | any `det`/`dim Fix` inconsistency. This row is the **cross-check that ties B3 to C4** — neither can be wrong alone without it reddening | L0 | 1,2,3 |

### D — independent constructions cross-checked

| ID | LAW | REFERENCE | FALSIFIER | lvl | d |
|---|---|---|---|---|---|
| **D1** | **Householder ≡ projector construction.** `I − 2n̂n̂ᵀ` equals `Σ_{i<d} e_ie_iᵀ − n̂n̂ᵀ` built from an orthonormal basis `{e_i}` of the hyperplane (via QR/Gram-Schmidt in the test) | **P4** — the test writes the spectral form `(+1 on the hyperplane, −1 on n̂)` from the definition, never calls the SUT | a sign or factor slip in either form. Structurally independent: one is a rank-1 update, the other a spectral sum | L0 | 2,3 |
| **D2** | **Rodrigues ≡ `scipy.spatial.transform.Rotation.from_rotvec`**, `atol=1e-14` over ≥500 random (axis, angle) | **P2** — quaternion-based, an independently maintained implementation with a different algorithm | any Rodrigues term slip. `[M]` max abs diff **9.2e-16** over 2000 draws. **Gate on `atol`, NOT on ULP** — `[M]` the same data reads "5120 ULP" because entries cross zero (§7-V7) | L0 | 3 |
| **D3** | ⭐ **Planar rotation ≡ `scipy.linalg.expm(θ(vuᵀ − uvᵀ))`**, `atol=1e-12` | **P3** — the Lie-theoretic definition `R = exp(θJ)` by Padé scaling-and-squaring: a *third* structurally independent route, and the only **dimension-generic** one | any dimension-generic rotation constructor bug. `[M]` agrees with the closed form to **1.7e-13** over d=2,3,4 (600 draws) and with Rodrigues to 2.4e-15 at d=3. **The 1e-12 tolerance is measured, not guessed** — a 1e-15 gate would red on Padé noise | L0 | 2,3,(4) |
| **D4** | ⭐ **Cartan–Dieudonné composition.** `σ_α ∘ σ_β = R(2(α−β))` — reflections whose planes meet at `θ/2` compose to a rotation by `θ` | **P1** symbolic at d=2 (`[M]` residual exactly `0`, 0.3 s) + **P4** numeric at d=3 with normals OFF the coordinate planes | swap `σ_α ∘ σ_β` for `σ_β ∘ σ_α` → rotation by `−θ`. This row is the **only gate that pins the rotation's SENSE** against an independent construction | L0 | 2,3 |
| **D5** | **`−I` = product of `d` mutually orthogonal reflections**, and `−I` is NOT expressible as a single `Reflection` or `Rotation` at d=3 | closed form; `[M]` `σ_xσ_yσ_z == −I` **bit-exact**; `[M]` `−I` at d=3 has `dim Fix = 0` ⟹ neither `d−1=2` nor `d−2=1` | typing `Inversion` as a `Reflection` subclass. **This row is an API gate** — see §8-A2 | L0 | 1,2,3 |

### E — conjugation ("seat it elsewhere" ≡ "re-orient it")

| ID | LAW | REFERENCE | FALSIFIER | lvl | d |
|---|---|---|---|---|---|
| **E1** | ⭐ **`R σ_n̂ R⁻¹ = σ_{R n̂}`** | **P1** symbolic for `R = R_z(t)`, `n̂ = n(φ)` (`[M]` residual exactly `0` in **0.6 s**, and it yields the bonus identity `R_z(t)n(φ) = n(φ+t)`) + **P4** numeric for general `R`, `atol=1e-14` (`[M]` 1.2e-15 over 500 draws) | conjugate the wrong way (`R⁻¹σR`) → the normal rotates by `−t`. **The FULLY symbolic version with a general Rodrigues `R` CHOKES SymPy (>120 s, killed)** — use the specialised symbolic + general numeric pair (§6) | L0 | 2,3 |
| **E2** | ⭐⭐ **AFFINE conjugation.** `g ∘ σ_{n̂,d} ∘ g⁻¹ = σ_{Qn̂, d + (Qn̂)·t}` for any rigid motion `g = (Q,t)` | closed form derived in the test from the definition of the conjugate mirror plane; `[M]` residual **2.7e-15** over 300 draws | drop the `(Qn̂)·t` term (i.e. conjugate only the linear part) — the mirror stays at the wrong offset. **This is the law that makes R1's "seat it elsewhere" and R5's "re-orient it" ONE operation**; without it the affine part is unverified under the group action | L0 | 1,2,3 |
| **E3** | **Conjugation is a class function.** `ord(ghg⁻¹) = ord(h)`, `det(ghg⁻¹) = det(h)`, `dim Fix(ghg⁻¹) = dim Fix(h)` | group theory | any conjugation implementation that is not a group automorphism | L0 | 2,3 |

### F — the seating relation (ruling R1)

| ID | LAW | REFERENCE | FALSIFIER | lvl | d |
|---|---|---|---|---|---|
| **F1** | **`seat(Q,c) = (Q, (I−Q)c)` fixes `c`**, `atol=1e-14` | the definition; `[M]` residual 8.9e-16, **bit-exact only 58/500** ⟹ tolerance required | `t = (Q−I)c` (sign flip) → `g(c) = 2Qc − c ≠ c` | L0 | 1,2,3 |
| **F2** | ⭐ **Every row of R1's table.** reflection `(σ_n̂, 2d n̂)`; rotation about `c` `(R, (I−R)c)`; inversion through `c` `(−I, 2c)`; translation `(I, v)` — each computed as `seat(Q, c)` and compared against the row's closed form | the table is itself derived, but the test recomputes each `t` from `(I−Q)c` **symbolically** (`[M]` SymPy: `(I−σ)(dn̂) − 2dn̂ = 0`, `(I−(−I))c − 2c = 0`) | any row's `t` formula. **The inversion row `t = 2c` and the reflection row `t = 2d n̂` are the two that carry a FACTOR OF 2** — exactly the factor C1 is blind to | L0 | 1,2,3 |
| **F3** | ⭐ **SINGLE SOURCE.** `Reflection(normal=n̂, offset=d)` is BIT-IDENTICAL to `seat(householder(n̂), d·n̂)`, and `Rotation(plane=P, angle=θ, centre=c)` to `seat(R(P,θ), c)` | — | **AT RISK OF BEING TAUTOLOGICAL** — see §7-V1. The gate must compare against `t` computed **arithmetically in the test** (`2*d*n̂`), not against a second call to `seat()` | L0 | 1,2,3 |
| **F4** | **`seat(Q, 0) = (Q, 0)` bit-exactly** — the origin-seated case IS the linear case | `(I−Q)·0 = 0` | a `seat` that always allocates a fresh `t` via a floating computation that returns `-0.0` — caught by `array_equal` only if `-0.0` is canonicalised. Pin `np.signbit(t).any() == False` alongside | L0 | 1,2,3 |

### O — `permutes` / `preserves`, and the orbit action (ruling R4)

| ID | LAW | REFERENCE | FALSIFIER | lvl | d |
|---|---|---|---|---|---|
| **O1** | **POSITIVE — a genuinely invariant set yields a genuine PERMUTATION.** `permutes` returns `π` with `np.array_equal(np.sort(π), np.arange(n))`, and `‖g(x) − x[π]‖ ≤ atol` re-checked in the test | **P6** exact integers. Fixtures: cube vertices/face centres/edge midpoints under `O_h`; square vertices under `D_4`; a **centred 1-D cell lattice** under `{e,σ}` — `[M]` `π(σ) == np.arange(n)[::-1]` **exactly**, at n=6 and n=7 | return the argmin map without asserting it is a permutation (ERR-073's shape). Anti-pattern #11: **this positive row is what makes O3/O4 mean something** | L0 | 1,2,3 |
| **O2** | ⭐⭐ **THE HOMOMORPHISM LAW.** `π_{g∘h} = π_g ∘ π_h`, i.e. `pi_gh == pi_g[pi_h]`; plus `π_e = arange(n)` and `π_{g⁻¹} = argsort(π_g)` | **P6** — pure integer identity, no reference needed, no tolerance possible | **the composition order** (`pi_h[pi_g]`): `[M]` disagrees on **102 of 144** pairs on the cube under `O_h`, 0/64 violations for the correct order at d=2. Also catches a transposed action (`M @ x` vs `x @ M.T`) and an inconsistent index convention. **VACUOUS ON AN ABELIAN GROUP** — the fixture MUST be non-abelian (§7-V3) | L0 | 2,3 |
| **O3** | **NEGATIVE-1 — a genuinely non-invariant set is REFUSED.** An asymmetric point set → `None` | — | remove the distance-window guard | L0 | 1,2,3 |
| **O4** | ⭐ **NEGATIVE-2 (ERR-073) — with GUARD ISOLATION.** A bit-identical duplicated node → `None`; and the test independently proves *which* guard fired: `[M]` on the 8-vertex cube + 1 duplicate, over all 48 `O_h` ops, **guard 1 passes 48/48, guard 3 passes 48/48, guard 2 passes 0/48** | **P6**; ERR-073 | delete the `np.unique(π).size == n` check. Carry `@pytest.mark.catches("ERR-073")` **only after** re-introducing the exact bug and confirming THIS test reds (a `catches` marker is a coverage claim) | L0 | 3 |
| **O5** | ⭐ **THE WINDOW IS AN EXPLICIT PARAMETER, and it bites.** A set perturbed off-symmetry by `δ` certifies iff `δ ≲ atol`; sweep `δ ∈ {1e-15, 1e-12, 1e-9, 1e-6}` against `atol ∈ {1e-13, 1e-7}` | — | `[M]` a set off-symmetry by **1e-9** certifies under a **1e-7** window with a **perfectly injective** π ⟹ **the bijectivity guard does NOT protect against a too-wide window** — they are independent failure modes. **If `atol` is a module constant rather than a parameter this row is SIGNATURE-TAUTOLOGICAL** (§8-A5) | L0 | 1,3 |
| **O6** | **`preserves` = `permutes` + ONE weight guard — positive AND negative.** On equal weights `preserves ⟺ permutes`; on **unequal** weights arranged non-invariantly, `permutes` SUCCEEDS and `preserves` FAILS | — | `[M]` two antipodal points with weights `(1,2)` under `σ_x`: `π = [1,0]` (valid permutation) while `w[π] = (2,1) ≠ w` ⟹ isolates guard 3 exactly. **VACUOUS ON ANY EQUAL-WEIGHT FIXTURE** (§7-V5) — every shipped quadrature orbit has equal weights within an orbit, so this fixture must be hand-built | L0 | 1,3 |
| **O7** | ⭐ **ORBIT–STABILISER, ABSOLUTE.** `\|orbit\|·\|Stab\| = \|G\|` **and** the per-site absolute values against the published Wyckoff table | **P5** — cube under `O_h` (\|G\|=48): vertices `8×6`, face centres `6×8`, edge midpoints `12×4`; square under `D_4` (\|G\|=8): vertices `4×2`, edge midpoints `4×2`, centre `1×8`; 1-D lattice under `{e,σ}`: `n` even → 0 fixed points, `n` **odd** → exactly 1. `[M]` all nine verified | a uniformly rotated realization keeps the *relation* `\|orbit\|·\|Stab\|=\|G\|` true (the existing `test_orbit_stabilizer_...` is blind to it — its own sibling docstring says so). **Only the ABSOLUTE site symmetries catch it.** Falsifier: rotate the whole generator set by an arbitrary `R` while keeping the point set — the relation survives, the absolute orbit types do not | L0 | 1,2,3 |
| **O8** | ⭐⭐ **THE SEAT IS FORCED — the centroid of a `G`-preserved weighted point set is `G`-FIXED.** If `g.preserves(P, w)` for every `g ∈ G`, then `g(c) = c` for `c = Σw_i x_i / Σw_i` | **provable in three lines** from `w_{π(i)} = w_i` and linearity — the pure-math backing for **R10**. `[M]` on a cube shifted to `(2.5,−1.25,0.75)`: all **48/48** *seated* `O_h` ops preserve it with `max‖g(c)−c‖ = 4.4e-16`, while only **1/48** *unseated* ops do (the identity). 1-D: centroid `= a + nh/2` exactly, seated `σ` → `arange(n)[::-1]`, unseated `σ` → no match at all | seat at the origin instead of the centroid — `[M]` **47 of 48** operators stop preserving the set. **This is the campaign's motivating defect, measured in the abstract**, and it is the row that proves the affine part of R1 is load-bearing rather than decorative | L0 | 1,2,3 |

### X — gates the LANDED design earns (see §0)

| ID | LAW | REFERENCE | FALSIFIER | lvl | d |
|---|---|---|---|---|---|
| **X1** | ⭐ **The two actions are distinct and consistent.** `on_points(x) − on_directions(x) == t` **exactly**; `‖on_directions(v)‖ = ‖v‖` to `atol=1e-14`; and `on_directions` is **invariant under `seated_at(c)` for every `c`** | the definition of an affine map and of orthogonality | let `on_directions` apply the translation — unit vectors silently denormalise, and the BC law `ψ_in(x,Ω) = ψ_out(g⁻¹x, Q_g⁻¹Ω)` becomes wrong on the direction leg. The **seating-invariance** leg is the sharp one: it is the law that says "where the mirror sits does not change which way it points" | L0 | 1,2,3 |
| **X2** | **The circle-point constructor preserves EXACTNESS where the angle constructor cannot.** `rotation_from_circle_point(point=(0.0, 1.0))` yields a matrix with **exact** zeros; `rotation(angle=π/2)` does not | exact float arithmetic; `np.cos(np.pi/2) = 6.1e-17 ≠ 0` | make `rotation_from_circle_point` route through an angle. **This is the API-side gate for campaign step G4** — it pins that the exact path EXISTS, without importing `roots_of_unity` (forbidden as a G2 reference, §1); G4 then supplies the exact circle points | L0 | 2,3 |
| **X3** | ⭐ **A glide/screw has INFINITE order even though its linear part does not.** `reflection(normal=e_x) @ translation_by(e_y)` has `element_order() is None`, while `reflection(normal=e_x).element_order() == 2` | group theory: an element translating off its own fixed subspace generates an infinite cyclic group | report the **linear part's** order for an affine element — a lie about the element, and precisely the error a "just take the order of `Q`" implementation makes. Pairs with **A6** (infinite *groups*) as the element-level analogue | L0 | 2,3 |

### H — dimension genericity (the reason ℝ¹ is first-class)

| ID | LAW | REFERENCE | FALSIFIER | lvl | d |
|---|---|---|---|---|---|
| **H1** | ⭐ **ℝ¹ is COMPLETE and CORRECT.** `O(1) = {+1, −1}`, `SO(1) = {+1}`. `Reflection(normal=[1])` is `[[−1]]` with `det = −1`, `dim Fix = 0`; seated at `c` it is `x ↦ 2c − x`, fixing `c` **exactly**; a d=1 `Rotation` is the identity or raises (§8-A1) | closed form; `[M]` `house([1]) = [[−1]]`, `det = −1`, `2c − c == c` exactly | any construction that assumes `d ≥ 2` (a cross product, a `[2]` index, an `np.eye(3)`). **Every A–G row must carry a d=1 parametrisation where the law is defined** — this is the row that proves the 1-D/3-D arm split is deletable rather than reconcilable | L0 | 1 |
| **H2** | ⭐ **The EMBEDDING law.** `ι: E(d) → E(d+1)`, `(Q,t) ↦ (diag(Q,1), (t,0))`, satisfies `ι(g∘h) = ι(g)∘ι(h)`, preserves `det`, and acts identically on the first `d` coordinates. Composed twice: `ι³∘ι²(g) = ι_{1→3}(g)` | the definition of the block embedding | pad with `0` instead of `1` on the diagonal → `det = 0`, no longer orthogonal. **This is the operation R2 relies on to dissolve the arm split** — if it is not a NAMED method the split just moves to the call sites (§8-A14) | L0 | 1→2, 2→3, 1→3 |
| **H3** | **Mixing dimensions RAISES.** `g_1d ∘ g_3d` must raise a named error | — | silently broadcasting or padding. Must call the **production** composition, not a locally raised exception (§7-V8) | L0 | 1,3 |

---

## 5. Prose for the subtle gates

### 5.1 ⭐ The involution law is Mode-12 BLIND to the affine offset — C1 alone is a trap

This is the single most important finding in the plan, and it is a textbook
`vv-principles` **Mode 12** (invariant-functional) instance.

A seated reflection is `g = (σ_n̂, t)`. Compose it with itself:

```
g ∘ g = (σ_n̂ σ_n̂ , σ_n̂ t + t) = (I, σ_n̂ t + t)
```

and for **any** `t` in `span(n̂)` we have `σ_n̂ t = −t`, hence `σ_n̂ t + t = 0`.
So **the involution law `σ² = e` holds for every offset whatsoever**. Its
invariance group contains the entire one-parameter family the offset lives in.

Measured, on `n̂ = e_z`, `d = 0.37`:

| `t` | `σ² = e`? | max `‖g(x) − x‖` on the plane `{z = d}` |
|---|---|---|
| `2d n̂` (correct) | **True** | `0.0` |
| `d n̂` (factor-2 bug) | **True** | `3.700e-01` |
| `4d n̂` (factor-2 the other way) | **True** | `7.400e-01` |
| `−2d n̂` (sign flip) | **True** | `1.480e+00` |

⟹ **C1 must never be credited as the catcher for the affine offset.** The
only catcher is **C3** (fixes the named hyperplane pointwise). This is exactly
the family of bugs the campaign exists to fix: `Mirror('x').is_invariant(mesh)`
reads `False` for the production slab *because* the mirror plane is at `x = a`,
not `x = 0` — an offset error is the campaign's motivating defect, and the
obvious gate cannot see it.

The same blindness has a rotation twin: `R` seated at any `c` has the same
order as `R` seated at the origin, so **C2 (order) is blind to the rotation
centre** too. C3 covers both.

### 5.2 ⭐⭐ The induced-permutation homomorphism (O2) is the highest-value gate

`_orbit_closure` computes `π_M` for every operator and has never checked that
the assignment `M ↦ π_M` is a **group homomorphism**. It is the deepest law
available on this object because it cross-checks *two independent layers* — the
element algebra (matrix composition) and the action (nearest-neighbour matching
on points) — against each other, using **only exact integer arithmetic**.

`π_{g∘h} = π_g ∘ π_h`, in array form `pi_gh == pi_g[pi_h]`.

What it catches that nothing else in the plan does:

* the **composition order** inside `RigidMotion.__matmul__` (measured: the
  opposite order disagrees on **102 of 144** cube/`O_h` pairs);
* a **transposed action convention** (`nodes @ M.T` vs `M @ nodes`) — a wrong
  transpose gives `π_{M}` for `Mᵀ`, and `π_{Mᵀ}` is not a homomorphism of the
  original product;
* an **index-convention drift** in the returned permutation (`π` vs `π⁻¹`) —
  `π⁻¹` composes in the *opposite* order, so the law reddens.

**It is VACUOUS on an abelian group.** On `C_n` alone, `pi_g[pi_h]` and
`pi_h[pi_g]` are equal for every pair, so the composition-order mutation is
invisible. The fixture MUST be non-abelian: `O_h` on the cube (d=3) and `D_4`
on the square (d=2). Measured: 0 violations for the correct order at both, and
102/144 for the wrong one at d=3.

### 5.3 The three guards of `permutes` fail into ONE value — isolate them from the INPUT side

`permutes -> permutation | None` returns `None` for three different reasons
(image is not a node / the map is not injective / — for `preserves` — the
weights disagree). A negative gate asserting `is None` therefore proves only
"one of three guards fired", which is not the coverage claim the test appears
to make. Since G1's API is not written yet, the clean fix is §8-A6 (return a
typed reason). If the API keeps a bare `None`, each negative gate MUST prove
the isolation **from the input side**, in the test, using pure numpy:

* **O4 (ERR-073 duplicate).** Measured on the 8-vertex cube + one bit-identical
  duplicate, over all 48 `O_h` operators: guard 1 passes **48/48**, guard 3
  passes **48/48**, guard 2 passes **0/48**. So the duplicate is a *clean
  isolator* for the bijectivity guard, and the test can assert exactly that
  (compute the argmin map itself, show every image is a node and every weight
  matches, then show the map is not injective, then assert the SUT returns
  `None`).
* **O6 (weight guard).** Two antipodal points with weights `(1, 2)` under
  `σ_x`: measured `π = [1,0]` — a *valid* permutation — while `w[π] = (2,1)`.
  So `permutes` must SUCCEED and `preserves` must FAIL on the same input. This
  is the only construction that isolates guard 3, and it is impossible on any
  shipped quadrature (weights are equal within an orbit).
* **O3 (distance guard).** An asymmetric set whose images land far outside the
  window.

### 5.4 The window and the bijection are INDEPENDENT failure modes (O5)

It is tempting to read the ERR-073 bijectivity guard as "the check is now
sound". It is not sufficient. Measured: the two-point set
`{(1,0,0), (−1+10⁻⁹,0,0)}` is **not** `σ_x`-invariant, yet under a `1e-7`
window the nearest-neighbour map is `π = [1,0]` — a perfect bijection — and
every guard passes. The residual is `1.0e-09`.

So the window is a first-class correctness parameter, not an implementation
detail, and it must be:

1. an **explicit argument** of `permutes`/`preserves` (else O5 is
   signature-tautological — §8-A5);
2. gated by a **perturbation ladder** (`δ` vs `atol`), asserting acceptance
   below and refusal above;
3. ideally **defaulted relative to the point set's minimum pairwise
   separation**, which is a computable, set-intrinsic quantity — see §9-N3.

Note also that `symmetry.py` today carries a hidden `_NODE_WINDOW_FACTOR = 100`
multiplying `atol` for nodes while weights use `atol` raw. `permutes` has no
weights at all, so the new API should take `atol` (points) on `permutes` and
`atol, wtol` on `preserves` — the 100× factor is an artefact of one parameter
serving two questions and should not survive the carve.

### 5.5 Orbit–stabiliser must be ABSOLUTE, not merely self-consistent (O7)

`tests/numerics/test_symmetry.py::test_orbit_stabilizer_theorem_holds_on_every_certificate`
already asserts `|orbit|·|Stab| = |G|` — but `|G|` is computed by the same
machinery, so the relation is a *self-consistency* check. It survives any
mutation that rotates the entire realization by a common `R` (the relation is
conjugation-invariant), which is precisely the class the `_vertical_mirrors`
π/2 bug belonged to; that file's own docstring at `symmetry.py:1118-1121`
concedes "orthogonality, determinant, closure and group order are all preserved
by a rotated mirror set, so none of those checks can see this".

O7 therefore pins the **absolute orbit types** against the published Wyckoff
decomposition, on point sets whose site symmetries are in a book:

| point set | group | `|orbit| × |Stab|` | site symmetry |
|---|---|---|---|
| cube vertices (8) | `O_h` (48) | `8 × 6` | `C_3v` |
| cube face centres (6) | `O_h` (48) | `6 × 8` | `C_4v` |
| cube edge midpoints (12) | `O_h` (48) | `12 × 4` | `C_2v` |
| square vertices (4) | `D_4` (8) | `4 × 2` | `C_s` |
| square edge midpoints (4) | `D_4` (8) | `4 × 2` | `C_s` |
| square centre (1) | `D_4` (8) | `1 × 8` | `D_4` |
| 1-D lattice, `n` even | `{e,σ}` (2) | `2 × 1` (all) | — (no fixed point) |
| 1-D lattice, `n` **odd** | `{e,σ}` (2) | `1 × 2` at the centre | `σ` |

All nine measured and confirmed. The **odd-`n` 1-D row is the cheapest
non-trivial stabiliser in the whole plan** and it is the exact structure the
r=0 pole continuity and the sweep-reversal consumers rely on.

### 5.6 Two-tier bit-identity — earn it where it is real, do not fake it elsewhere

Three laws are **bit-exact** and should be gated with `np.array_equal` (a
strictly stronger gate, and free):

* **A2** identity: `[M]` 500/500 bit-exact.
* **B4/F4** seating preserves `Q`, and `seat(Q,0)` returns `t = 0`.
* **A1** associativity **restricted to signed permutations**: `[M]` 500/500.
* **D5** `σ_xσ_yσ_z == −I`: `[M]` bit-exact.
* **O1** `π(σ) == arange(n)[::-1]` on a centred lattice: `[M]` exact (integers).

Four are **NOT** bit-exact and a bit-exact gate there would be a false red:

* **A1** associativity for general elements: `[M]` **0/500** bit-exact.
* **A3** inverse: `[M]` **0/500**.
* **F1** seating fixes `c`: `[M]` 58/500 only.
* **C1** `σ² = e` for a **general** (non-coordinate) normal: `[M]` **0/300** —
  while it IS bit-exact for a coordinate normal. This split is exactly the
  `_reflections` (diagonal sign flip) vs `_vertical_mirrors` (Householder)
  divide the campaign is dissolving, so **both must appear in the fixture**.

And one is **bit-exact in the permutation but not in the points**: `[M]` a
seated reflection of a uniform cell-centre lattice lands its images
`5.6e-17 – 2.2e-16` away from their partners (never bit-exact, at four
different `(a, h, n)`), while the induced **permutation** is exactly
`arange(n)[::-1]`. That is why `permutes` needs a window at all, and why the
campaign's "bit-identity preservation set" (which includes
`directional.py:574`'s `arange(N)[::-1]`) is about the *integers*, not the
*coordinates*.

---

## 6. The SymPy budget, and where SymPy chokes

Per `algebra-of-record`'s choke-mode catalog, measured on this exact problem:

| symbolic claim | result | time |
|---|---|---|
| d=2 reflection `QᵀQ − I`, `det` | `0`, `−1` | **0.1 s** |
| d=3 reflection `QᵀQ − I`, `det` (spherical `n̂`) | `0`, `−1` | **0.5 s** |
| Cartan–Dieudonné at d=2, `σ_α σ_β − R(2(α−β))` | `0` | **0.3 s** |
| Rodrigues `det`, `trace` | `1`, `2cos t + 1` | **0.8 s** |
| Rodrigues `QᵀQ − I` | `0` | **5.0 s** |
| Conjugation, **specialised** `R_z(t) σ_{n(φ)} R_z(t)ᵀ − σ_{Rn̂}` | `0`, plus the bonus identity `R_z(t)n(φ) = n(φ+t)` | **0.6 s** |
| Conjugation, **general Rodrigues `R` × general Householder `S`** | ⛔ **CHOKE — killed at 120 s** | — |
| Naive "impose `Σnᵢ² = 1` by `subs`" instead of a spherical parameterisation | ⛔ **does not fire** — leaves `4n₀²(Σnᵢ²−1)` residuals | — |

Two operational rules follow, and they must be written into the test file:

1. **Parameterise the unit vector explicitly** (`n̂ = (sinθcosφ, sinθsinφ, cosθ)`);
   never try to impose `|n̂| = 1` by substitution after expansion.
2. **Specialise ONE factor of a conjugation.** The general-times-general
   conjugation is out of reach; `R_z(t)` with a symbolic angle against a
   symbolic in-plane normal closes in 0.6 s and proves the same law shape.
   The *general* conjugation is covered numerically (E1, E2) over ≥300 random
   draws.

Total symbolic budget **≈ 7 s**. Put the SymPy rows behind
`@pytest.mark.slow`? **No.** Recommendation: keep all symbolic rows in the
same file, accept ~8 s total (the Rodrigues `QᵀQ` row at 5.0 s is the outlier),
and do **not** mark them slow — a
`slow`-marked foundation gate is deselected by the project's canonical
`-m "not slow"` sweep and would stop guarding.

---
## 7. Vacuity / tautology audit — gates at risk, and how to make each reddenable

Every row below is a gate *in this plan* that could ship green-and-useless.
Each carries the fix. The precedent the campaign already paid for is the
**stable-sort tie-break gate that was VACUOUS at n ∈ {4,8,12}** because numpy's
unstable `argsort` happens to agree with stable below ~24 — a gate satisfied by
construction *for the sizes chosen*. Four rows below are that same shape.

| # | gate | risk | why | FIX |
|---|---|---|---|---|
| **V1** | **F3** (single source: `reflection(n̂,d)` ≡ `reflection(n̂).seated_at(d n̂)`) | ✅ **RISK RETIRED BY MEASUREMENT** — but the recipe still applies | the risk was "`reflection` internally calls `seated_at`, so the gate is `X == X`" (the **L34e** shape). `[M]` against the shipped code the two are **genuinely different code paths**: `reflection` computes `t = 2d n̂` directly (**bit-identical** to `2*d*n_hat` recomputed in the test) while `seated_at` goes through `(I−Q)c` — they differ by **1.11e-16** | keep BOTH legs and gate them differently: the direct leg with `assert_array_equal(g.translation, 2*d*n_hat)` (pins the FORMULA), the two-path leg with `atol=1e-15` (pins the ROUTING). Never merge them — merging is what would re-create the tautology if `reflection` is ever refactored to delegate |
| **V2** | **C2** (`ord = n` exactly) | ⚠ **SIZE-BLIND** if parametrized only over odd `n` | `[M]` the `k=2` generator mutation gives `ord = n/gcd(2,n)`: caught at n=4,6,8,12 (order halves), **NOT caught at n=3** (order stays 3) | parametrize `n ∈ {2,3,4,6,8,12}` — **even values are mandatory**. Assert `ord == n` exactly, never `Q**n == e` alone (satisfied by any element whose order divides n) |
| **V3** | **O2** (permutation homomorphism), **A4** (anti-homomorphism) | ⛔ **VACUOUS ON AN ABELIAN GROUP** | on `C_n` every pair commutes, so the composition-order mutation moves nothing. `[M]` the wrong order is caught on **102/144** `O_h` pairs and **0/n²** `C_n` pairs | fixtures MUST be non-abelian: `O_h` on the cube (d=3), `D_4` on the square (d=2). At d=1 the group `{e,σ}` is abelian and the law is *unavoidably* vacuous — declare that, do not pretend the d=1 row covers it |
| **V4** | **B4** (`det(gh) = det(g)det(h)`) | ⛔ **VACUOUS on proper-only draws** (`1 = 1·1`) | a rotation-only fixture cannot see a sign error | mix proper and improper: at least one reflection × reflection (→ +1), reflection × rotation (→ −1), and (at d=3) `−I` |
| **V5** | **O6** (`preserves` weight guard) | ⛔ **VACUOUS on any equal-weight fixture** | every shipped quadrature has equal weights *within* an orbit, so `preserves ≡ permutes` there and guard 3 never fires | hand-build the unequal-weight witness. `[M]` two antipodal points, weights `(1,2)`, under `σ_x` |
| **V6** | **C4** (fixed-set dimension) | ⚠ **NUMERICALLY FRAGILE if computed by eigenvalue-`isclose`** | `[M]` at rotation angle `1e-9`, eigenvalue-`isclose(1.0)` reports `dim Fix = 3` (wrong); `matrix_rank(Q−I)` with the DEFAULT relative tolerance reports `1` (right) down to `1e-9`, and only an *absolute* `tol=1e-10` breaks it at `1e-12` | use `d − np.linalg.matrix_rank(Q − I)` with the **default** tolerance, and parametrize angles bounded away from 0 (`θ ∈ [0.1, π]`) |
| **V7** | **D2/D3** (scipy cross-checks) | ⚠ **FALSE RED if gated on ULP** | `[M]` Rodrigues vs `from_rotvec` is `9.2e-16` in absolute terms but reads **5120 "ULP"** because matrix entries cross zero — a known campaign lesson (*ULP is the wrong metric for arrays crossing zero*) | gate `np.testing.assert_allclose(..., atol=1e-14, rtol=0)` for D2 and `atol=1e-12` for D3 (Padé noise `[M]` 1.7e-13). Never `assert_array_almost_equal_nulp` |
| **V8** | **A6** (infinite closure raises), **H3** (dimension mismatch raises) | ⛔ **SELF-SATISFIED `pytest.raises`** if the test constructs and raises the error itself | the `vv` fifth fires-but-cannot-fail class; a boundary review found **9** such legs in one file, none of which reddened under 14 guard-disabling mutations | the raising call must be the **production** `close(...)` / `g1 @ g3`. The mutation proof is "remove the cap / remove the dim check — does THIS file red?" |
| **V9** | **C1** (involution) credited for the offset | ⛔ **MODE-12 DESIGNED-GREEN** — see §5.1 | the involution functional's stabiliser contains the whole `t ∈ span(n̂)` family | never credit C1 for the offset; C3 is the catcher. Write the blindness into C1's own docstring so a future reader cannot re-mint the false claim |
| **V10** | **C1, D1, B2** on **coordinate normals only** | ⚠ **SIZE-BLIND, campaign-specific** | with `n̂ = e_x` the Householder degenerates to a diagonal sign flip — `[M]` `σ²==I` is then bit-exact (0/300 for general normals). A coordinate-only fixture tests `_reflections`'s degenerate case and is structurally blind to the general `I − 2n̂n̂ᵀ` path, which is **precisely the split this campaign exists to dissolve** | every reflection gate carries **≥1 non-coordinate normal**, e.g. `n̂ ∝ (1,2,3)` and `n̂ ∝ (1,1,0)/√2` (a `D_nh` diagonal plane) |
| **V11** | **O7** (orbit–stabiliser) if it asserts only the *relation* | ⚠ blind to a globally rotated realization | the relation `\|orbit\|·\|Stab\| = \|G\|` is conjugation-invariant — see §5.5 | assert the **absolute** `(\|orbit\|, \|Stab\|)` pairs against the Wyckoff table. Falsifier to run once: conjugate the whole generator set by an arbitrary `R`, keep the point set — the relation survives, the absolute types do not |
| **V12** | **A3** (inverse) on reflection-only draws | ⚠ blind to the missing conjugation in `t ↦ −Qᵀt` | every reflection has `Q = Qᵀ`, and for a *seated* reflection `t ∈ span(n̂)` so `−Qᵀt = t` — a `t ↦ −t` bug is invisible on the whole reflection family | include general rotations with a general `t` in the A3 fixture |

**Mutation-proof protocol for the whole file** (the standing discipline): each
gate is credited only after its named falsifier has been shown to RED **in
process** (`monkeypatch` on the module attribute, never `git checkout` — two
tracked skill files carry irrecoverable uncommitted state), under the canonical
`python -O -m pytest`, serial. Run the battery once at the end of G2 and record
the per-gate verdict; **an all-blind verdict means the harness is broken, not
that the code is perfect** — include a deliberate positive-control mutation
(e.g. `Reflection` returns `+I`) that must red *many* gates.

---

## 8. ⛔ API problems — laws the sketched G1 API CANNOT express

> **Read §0's re-score first.** This section was written against the sketched
> API in the brief. G1 landed mid-plan; **A2, A3, A5, A7, A8, A9, A10, A11,
> A12, A13 are DISCHARGED** by the shipped design (several more elegantly than
> proposed) and **A1, A4, A6, A14, A15 remain OPEN**. The reasoning below is
> retained because it is the *why* behind each requirement — a future change
> that re-opens a discharged item needs the argument, not just the verdict.

This is the section with the highest leverage: the module is not written, so
each of these is free to fix now and expensive later. Ordered by severity.

### A1 — `Rotation` at d=1 is undefined, and R2's table says `n−2 = −1`

**R2** parameterises a rotation by "a rotation plane (dim `n−2` fixed set) +
angle". At `d=1` there is no plane and `n−2 = −1`. Mathematically `SO(1) = {+1}`
is trivial: **the only d=1 rotation is the identity.** The API must decide,
explicitly:

* `Rotation(plane=…, angle=…)` requires `d ≥ 2` and raises `ValueError` at d=1; **or**
* a d=1 rotation is admitted and is always `e`.

Either is defensible; *silence* is not — a `np.cross`-flavoured implementation
will simply `IndexError`, and H1 (the row that proves ℝ¹ is first-class) will
red for the wrong reason. **Recommendation: raise, with the message naming
`SO(1) = {e}`.** ℝ¹'s non-trivial content lives entirely in `Reflection`.

### A2 — `−I` has no home: its TYPE changes with dimension

`[M]` `det(−I) = (−1)^d`: **−1 at d=1, +1 at d=2, −1 at d=3**, and at d=3
`dim Fix(−I) = 0`, which is neither `d−1 = 2` (reflection) nor `d−2 = 1`
(rotation). So:

| d | what `−I` IS |
|---|---|
| 1 | a **reflection** (through the origin point) |
| 2 | a **rotation** by π |
| 3 | **neither** — a rotoreflection `S_2`, the inversion |

R1's table lists "inversion, centre `c`" as a row, but the sketched API has only
`Reflection` and `Rotation` constructors plus seating. **`Inversion(centre=c)`
must exist as a CONSTRUCTION (a factory returning a `RigidMotion`), never as a
subclass of either**, and `RigidMotion` must be constructible from a validated
raw `(Q, t)` pair — otherwise the elements produced by `close()` (which are
products, and can be rotoreflections) cannot be re-typed at all. Gate **D5**
exists to pin this.

Corollary: the existing `symmetry.py` gap "`_inversion_op` (−I) is NOT
expressible — no `C_i`/`S_n` tag; reachable only inside `I_h`" is an *element*
gap, and this is where it closes.

### A3 — `Reflection(n̂, d)` has a GAUGE FREEDOM that breaks equality and closure

`(n̂, d)` and `(−n̂, −d)` denote the **same** mirror and `[M]` produce
**bit-identical** matrices. If `__eq__`/`__hash__` compare the stored normal,
then `⟨σ⟩` closes to **3** elements instead of 2 and A5 reds — or worse, a
closure dedup silently keeps both and every group order inflates.

**Requirement:** equality and hashing are on the realized `(Q, t)`, with a
canonical rounding (the existing `_close_group` uses `np.round(M, 9).tobytes()`
plus a `+0.0` to canonicalise `−0.0` — **keep both tricks**, and extend the key
to include `t`). Optionally canonicalise the normal at construction (first
non-zero component positive) — but the *matrix* key is what the closure needs.

### A4 — `close()` needs a cap and a NAMED error; some rigid-motion groups are infinite

R1 deliberately admits translations (periodic BC). `[M]` `⟨translation(v)⟩` is
infinite, and so is `⟨σ_{n̂,0}, σ_{n̂,1}⟩` — **two parallel mirrors generate the
infinite dihedral group**, and that is a configuration the mesh/BC layer can
easily hand you (two reflective faces of a slab!). The closure MUST carry the
`_MAX_GROUP_ORDER` cap and raise a **named** error so gate A6 can be specific
rather than catching a bare `ValueError`. `geometry/boundary/_errors.py` already
establishes the `BoundaryError(ValueError)` idiom; mint the analogue
(`NotAFinitePointGroupError` or similar) in `transformation.py`.

### A5 — `permutes` MUST take the tolerance as a PARAMETER (else O5 is signature-tautological)

If the match window is a module constant, the *production path can never supply
a varying window*, and O5's "the window bites" claim becomes a
**signature-tautological** gate in the `vv` Mode-8 sense: a hand-injected
falsifier will move it (so a naive "does it red?" check passes and gives false
confidence), but no production call can ever exercise the variation. Measured
motivation: `[M]` a 1e-9-off-symmetry set certifies under a 1e-7 window with a
perfectly injective π.

**Requirement:** `permutes(points, *, atol)` and `preserves(points, weights, *,
atol, wtol)`. Two windows, because they answer two questions — `symmetry.py`'s
hidden `_NODE_WINDOW_FACTOR = 100` exists only because one `atol` was serving
both, and it should not survive the carve.

### A6 — `permutes -> permutation | None` collapses THREE failure reasons into one value

A negative gate asserting `is None` proves "one of three guards fired". §5.3
gives an input-side workaround for the two witnesses we have, but the clean fix
is in the API: return a typed outcome carrying the reason
(`NotAPermutation.ImageNotInSet | .NotInjective | .WeightMismatch`), or expose
the guards as separately callable predicates. This is cheap now and impossible
later without touching every consumer.

### A7 — ⭐ `preserves(measure)` (R4) and "geometry is numpy-only" (R3) CANNOT BOTH HOLD

`DiscreteMeasure` lives in `orpheus/numerics/`. R3 puts the orbit primitive in
`orpheus/geometry/transformation.py`, **numpy-only, no numerics import**. So
`g.preserves(measure)` as literally spelled in R4 is **not implementable**
unless `preserves` takes bare `(points, weights)` arrays — which is
anti-pattern #13 landing on the very boundary this campaign exists to type.

Three ways out; pick one **before** G1:

1. ⭐ **Mint a geometry-side `WeightedPointSet`** (frozen dataclass: `points
   (N,d)`, `weights (N,)`, `dim`), and let `numerics.DiscreteMeasure` convert to
   it. This also gives the window a natural home (`min_separation` is a property
   of the point set — see §9-N3) and gives O5/O6/O7 a typed fixture.
2. `preserves` lives only on the **numerics wrapper**, and geometry exposes
   `permutes` alone. Costs R4's symmetry ("two predicates, not a flag").
3. Move `DiscreteMeasure` to geometry. Largest blast radius; not recommended.

### A8 — ⭐ `Rotation` must admit an EXACT RATIONAL TURN, or G4 is unreachable through the type

Campaign step **G4** is "close the checker-side `roots_of_unity` sites
(`_cyclic_ops`, `_vertical_mirrors`)", and the audit measured that
`_cyclic_ops(4)/(8)` produce **zero** exact zeros where `roots_of_unity`
produces two. If `Rotation` only accepts `angle: float`, that fix cannot be
expressed **through the new type** — the exactness would have to be re-hacked
at each call site, recreating the very duplication the campaign removes.

**Requirement:** a constructor taking the turn as an exact rational —
`Rotation.by_turn(p, q, plane=…)` or `turns: Fraction` — that routes
`(cos, sin)` through `roots_of_unity`. Then `_cyclic_ops` becomes
`[Rotation.by_turn(k, n) for k in range(n)]` and the exactness is structural.

**Note the G2/G4 boundary:** verifying the *exactness* (bit-exact mirror
partners, exact zeros) is **#325 / G4 work**, not G2. G2's obligation is only
that **the API admits the exact constructor** — a signature gate. Do not import
`roots_of_unity` into a G2 assertion (it is on the forbidden-reference list,
§1); gate the *existence and agreement* of `by_turn(p,q)` with
`Rotation(angle=2πp/q)` to `atol=1e-14`, and leave exactness to G4.

### A9 — `Reflection` needs a seat-from-data constructor, or the campaign's motivating defect stays unfixable

`[M]` (campaign audit) `Mirror('x').is_invariant(mesh)` is `False` for the
production slab and sphere **solely** because meshes start at `origin = 0.0`
while the mirror plane sits at `x = a`. R10 says the canonical seat is the
**centre of mass**. So the API needs, beside `Reflection(normal=…, offset=…)`,
a data-driven seat: `Reflection.bisecting(lo, hi, axis=…)` or
`element.seated_at(centroid)`. Without it, every consumer recomputes the offset
— Pattern 7 violated at birth.

### A10 — the composition ORDER and its spelling must be fixed and documented

`g @ h` must mean "apply `h` first" (matching `(Q₁Q₂, Q₁t₂ + t₁)`). `[M]` the
opposite reading disagrees on **102/144** cube pairs, and **gate O2 is the only
thing in the plan that pins it**. State the convention in the class docstring
and gate it with O2.

### A11 — ONE action entry point

`g.apply(points)` for `(N,d)` row-vector points, and the matrix must not be the
public way to act (`symmetry.py` today writes `nodes @ M.T` at the call site).
Otherwise the row-vs-column convention becomes a 7th spelling of `σ_a`, joining
the four BC vocabularies the campaign is consolidating.

### A12 — `det` / `is_proper` / `dim` / `order` belong ON the element

Today none of the four `σ_a` vocabularies knows its own determinant. Put
`det`, `is_proper`, `dim`, `fixed_set`, and `order()` on `RigidMotion`; the
gates then read like the math, and the four vocabularies have something to
collapse onto.

### A13 — mixing dimensions must RAISE

A d=1 element composed with a d=3 element is unspellable in `E(d)`. Numpy will
happily broadcast some of these. Gate H3.

### A14 — ⭐ the EMBEDDING `ι: E(d) → E(d+1)` must be a NAMED method

R2's claim is that dimension genericity **dissolves** the 1-D/3-D arm split.
That only happens if the embedding is an operation of the type
(`g.embed(3)` → `(diag(Q,1), (t,0))`), gated as a homomorphism (H2). If it is
ad-hoc `np.eye(3)` padding at the call sites, the arm split has *moved*, not
gone — and the campaign's own G3 deliverable ("delete the 1-D/3-D arm split")
loses its instrument.

### A15 — `permutes` should return a `Permutation` TYPE, not a bare int array

Two payoffs, both directly on this plan's gates:

* **ERR-073 becomes unrepresentable, not merely guarded.** A frozen
  `Permutation` whose `__post_init__` asserts `sort(π) == arange(n)` means that
  *returning one at all* is the bijectivity proof (Pattern 4). The guard stays
  as the `None`-returning branch; the type is belt-and-braces.
* **O2 (the homomorphism gate) reads like the math**: `pi_gh == pi_g @ pi_h`
  instead of `pi_gh == pi_g[pi_h]`, with `inverse()` and `fixed_points`
  available — and `fixed_points` is exactly `OrbitCertificate.singular_set`'s
  ingredient, so the existing certificate can be rebuilt on it.

Caveat for the falsifier: if the injectivity count-check is deleted, the
`Permutation` constructor raises instead of the function returning `None`, so
O4 reds with a *different* exception. That is still a red — but the gate should
assert `is None` and the mutation record should note the substitution, so a
future reader does not mistake the raise for the documented mechanism.

---
## 9. Laws I could NOT find a pure-math reference for

These are the honest gaps. Each is either a **convention** (no theorem exists),
a **design decision** (math is silent), or a **domain claim** (belongs to a
later step). None of them may be smuggled into G2 wearing a pure-math label.

| # | law / quantity | why there is no pure-math reference | where it belongs instead |
|---|---|---|---|
| **N1** | **The `D_nh` vertical-mirror plane placement** — planes at `kπ/n` (normals at `kπ/n + π/2`), vertex on the x-axis | this is the *standard setting* (Hamermesh; Int. Tables Vol. A), a **literature convention**, not a theorem. Every rotated mirror set is an equally valid `D_nh`; `[M]` orthogonality, determinant, closure and group order are all preserved by the rotation, so **no intrinsic law can distinguish them** | already owned by `tests/numerics/test_symmetry.py::test_vertical_mirror_planes_follow_the_dnh_setting`, whose reference is the convention + the shipped azimuthal rules. **Stays there; G2 must not claim it.** G2 gates only that *any* unit normal yields a valid reflection |
| **N2** | **Which of `(n̂, d)` / `(−n̂, −d)` is canonical** | pure gauge — `[M]` the two produce bit-identical matrices | a documented API decision (§8-A3). Gate the *consequence* (equality/hash/closure order), never the choice |
| **N3** | **The numeric VALUE of the match window `atol`** | `[M]` lower bounds exist and are measured (a seated cell-lattice reflection lands `5.6e-17 – 2.2e-16` from its partner; a general rigid-motion product drifts `1.8e-15`; Padé `expm` noise `1.7e-13`) — but **nothing bounds it from above except the point set's own minimum pairwise separation**, which is a property of the *consumer's data*, not of the algebra | **Recommendation:** default `atol` relative to `min_separation(points)` — a computable, set-intrinsic quantity — and gate *that* relation (accept below, refuse above) rather than pinning an absolute constant. This converts a modelling choice into a checkable law. It also retires `_NODE_WINDOW_FACTOR = 100` honestly |
| **N4** | **"The origin is the centre of mass" as a PHYSICAL computation** (R10, from material/nuclide densities) | the **theorem** — a `G`-preserved weighted set has a `G`-fixed centroid — IS pure math and IS gated (**O8**). The step from *material density* to *the weights of a point set* is a domain modelling claim | G7 (test migration) / the mesh consumer. G2 gates the theorem on synthetic weights only |
| **N5** | **Whether a d=1 `Rotation` raises or returns `e`** | math says only `SO(1) = {+1}`; it does not say what an API should do | a documented API decision (§8-A1) |
| **N6** | **Which failure reasons `permutes` should distinguish** | math says only "the map is / is not an equivariant bijection" | an API decision (§8-A6), gated by the *isolation* evidence in §5.3 |
| **N7** | **The choice of `π` vs `π⁻¹` as the returned permutation** | both are faithful; only the composition law (**O2**) is convention-sensitive, and it fixes the pair `(convention, composition order)` jointly, not either alone | document the convention on the return type; **O2 is what makes the pair checkable** |

---

## 10. Measured-facts appendix — every tolerance in this plan, with its provenance

Probe scripts: `/Users/rodrigo/.claude/jobs/c30e4f25/tmp/probe{1,2,3,4,7,9,10}.py`
(numpy 2.4.4, scipy 1.17.1, sympy 1.14.0, Python 3.14.3, `.venv`).

| quantity | measured | gate | consequence |
|---|---|---|---|
| Rodrigues vs `scipy…from_rotvec`, 2000 draws | max abs **9.2e-16**; "max ULP" **5120** (entries cross zero). **Re-measured against the SHIPPED `rotation_about_axis`, 300 draws: 1.11e-15** | D2 | gate `atol=1e-14`, **never ULP** |
| planar `expm` vs closed form, d=2,3,4, 600 draws | **1.7e-13**. **Re-measured against the SHIPPED `rotation(plane=…)`, 300 draws: 3.91e-14** | D3 | gate `atol=1e-12` (holds with margin) |
| `expm` vs Rodrigues at d=3, 200 draws | **2.4e-15** | D3 | — |
| `R σ_n̂ Rᵀ − σ_{Rn̂}`, 500 draws | **1.2e-15** | E1 | `atol=1e-14` |
| affine conjugation residual, 300 draws | **2.7e-15** | E2 | `atol=1e-14` |
| associativity, general elements, 500 triples | `‖ΔQ‖ ≤ 2.2e-16`, `‖Δt‖ ≤ 1.8e-15`, **0/500 bit-exact** | A1 | tolerance tier |
| associativity, signed permutations, 500 triples | **500/500 bit-exact** | A1 | `array_equal` tier |
| `e∘g = g∘e = g`, 500 draws | **500/500 bit-exact** | A2 | `array_equal` |
| `g∘g⁻¹`, `g⁻¹∘g`, 500 draws | **1.8e-15 / 1.1e-15**, 0/500 bit-exact | A3 | `atol=1e-14` |
| `(gh)⁻¹` vs `h⁻¹g⁻¹` (correct) | **2.2e-15**; swapped order differs **500/500** | A4 | non-vacuous |
| `seat(Q,c)` fixes `c`, 500 draws | **8.9e-16**, bit-exact only **58/500** | F1 | `atol=1e-14` |
| `σ² = I` bit-exact | coordinate normal **yes**; general normals **0/300** | C1 | two tiers, §7-V10 |
| `‖Q^n − I‖` vs `min_{k<n}‖Q^k − I‖` at n=24 | **1.0e-15** vs **1.6e-01** | C2 | 6 orders of margin; `atol=1e-10` is safe |
| `k=2` generator mutation | order → `n/gcd(2,n)`: caught at even n, **NOT at n=3** | C2 | even n mandatory |
| involution vs offset bug | **all** of `t ∈ {d n̂, 2d n̂, 4d n̂, −2d n̂}` are involutions; C3 residuals `0 / 0.37 / 0.74 / 1.48` | C1, C3 | §5.1 |
| `dim Fix` at rotation angle `1e-9` | eigenvalue-`isclose` → **3** (wrong); `matrix_rank` default → **1** (right) | C4 | use `matrix_rank` |
| Cartan–Dieudonné sharp `k = d − dim Fix`, `det = (−1)^k` | verified on **8** cases (refl, rot, `−I` at d=1,2,3, rotoreflection, identity) | C5 | — |
| `π(AB)` vs `π(A)[π(B)]`, cube/`O_h` | **0** violations; the swapped order violates **102/144** | O2 | non-abelian fixture mandatory |
| `π(AB)` vs `π(A)[π(B)]`, square/`D_4` | **0/64** violations | O2 | d=2 arm |
| ERR-073 duplicate, 48 `O_h` ops | guard1 **48/48** pass, guard2 **0/48** pass, guard3 **48/48** pass | O4 | clean isolation of guard 2 |
| window false-positive | a `1e-9`-off set certifies under a `1e-7` window with an **injective** π | O5 | window ⊥ bijection |
| unequal-weight witness | `π=[1,0]` valid, `w[π]=(2,1) ≠ w` | O6 | isolates guard 3 |
| Wyckoff orbit types | cube/`O_h`: `8×6`, `6×8`, `12×4`; square/`D_4`: `4×2`, `4×2`, `1×8`; 1-D: `n` odd → exactly 1 fixed point | O7 | all `= |G|` |
| centroid theorem | seated `O_h` preserves the shifted cube **48/48** (`‖g(c)−c‖ ≤ 4.4e-16`); unseated **1/48** | O8 | the affine part is load-bearing |
| seated reflection of a cell lattice | image lands **5.6e-17 – 2.2e-16** from its partner (never bit-exact, 4 configs); **permutation** exactly `arange(n)[::-1]` | O1, O5 | integers exact, coordinates not |
| `σ_xσ_yσ_z == −I` | **bit-exact** | D5 | `array_equal` |
| `det(−I)` at d = 1,2,3 | **−1, +1, −1** | B3, §8-A2 | d=2 is blind |
| `σ_n̂` vs `σ_{−n̂}` | **bit-identical**, 200 draws | §8-A3 | gauge freedom |
| closure orders | `⟨σ_x,σ_y,σ_z⟩ = 8`; `⟨R_z(π/2), σ_x⟩ = 8`; `⟨σ⟩ = 2`; translation and two-parallel-mirrors **unbounded** | A5, A6 | — |
| SymPy timings | 0.1 / 0.3 / 0.5 / 0.6 / 0.8 / **5.0** s; general×general conjugation **CHOKES (>120 s)** | B1, D4, E1, §6 | specialise one factor |

---

## 11. Acceptance criteria for G2 (what "done" means)

G2 is done when **all** of the following hold — not when the tests pass.

1. **All 42 gates exist** (A1–A6, B1, B2a, B2b, B3, B4, C1–C5, D1–D5, E1–E3,
   F1–F4, O1–O8, X1–X3, H1–H3), are `@pytest.mark.foundation`, and run under
   the canonical `python -O -m pytest`, **serial**, in the target time.
   **H2 is BLOCKED** until the embedding `ι: E(d) → E(d+1)` exists (§0/A14).
2. **Every gate has a named falsifier that has been RUN and RED**, in-process
   (`monkeypatch`), with the verdict recorded per gate. The battery includes a
   **positive control** (`Reflection` returns `+I`) that must red many gates —
   an all-blind verdict means the harness is broken (L34d).
3. **Every §7 vacuity row is discharged** — in particular: non-abelian fixtures
   for O2/A4, even `n` for C2, unequal weights for O6, a non-coordinate normal
   in every reflection row, mixed parity in B4, `matrix_rank` in C4, `atol`
   (not ULP) in D2/D3, and production entry points in every `pytest.raises`.
4. **No forbidden reference appears** — `grep` the test file for
   `numerics.symmetry`, `_rotation`, `_reflections`, `_orbit_closure`,
   `_close_group`, `Quadrature`, `DiscreteMeasure`, `roots_of_unity`: **zero
   hits**. (The G3 bit-identity gate is a *different file*.)
5. **The five OPEN §8 items have an explicit decision recorded** in the module
   docstring: **A14** (mint `RigidMotion.embedded_in(d)` — this one BLOCKS
   campaign step G3's "delete the 1-D/3-D arm split"), **A6** (typed failure
   reason for `permutes`, or accept §5.3's input-side isolation as permanent),
   **A4** (a named error class instead of the bare `ValueError` from
   `close_group`), **A1** (the d=1 rotation message should name `SO(1) = {e}`,
   not blame the caller's plane), **A15** (a `Permutation` return type).
6. **The `-rs` skip list is empty** and `--runxfail` shows every xfail (if any)
   failing for its documented reason.
7. **The measured tolerances in §10 are transcribed into the test file as named
   constants with their provenance**, so a future tightening is a deliberate act
   rather than a guess.

### The three gates to write FIRST (they shape the API most)

1. **O2** — the permutation homomorphism. It pins the composition order, the
   action convention and the index convention *simultaneously*, and it is the
   cheapest gate in the plan (pure integers).
2. **C3** — the fixed-point set, pointwise. It is the only catcher for the
   affine offset, which is the campaign's motivating defect.
3. **O8** — the centroid theorem. It proves the affine part of R1 is
   load-bearing (`[M]` 48/48 seated vs 1/48 unseated) and it is the pure-math
   backing for R10, i.e. it converts "where do we put the origin?" from a
   convention into a computed fact.

---

## 12. DELIVERED — `tests/geometry/test_transformation.py` (2026-08-03)

The gates are written and the mutation battery is run. This section is the
closeout; §4 remains the specification.

**Shipped:** 42 gates in 39 test functions, 96 parametrised cases,
`@pytest.mark.foundation` module-wide, **8.9 s** under
`.venv/bin/python -O -m pytest`, serial. `tests/geometry` = **1 failed
(the known pre-existing `TestWhiteXminPartial03GLSnapshot` red), 672 passed,
4 skipped, 1 xfailed**. pyright **0 errors**. No skips and no xfails in the
new file. Forbidden-symbol grep clean.

### Mutation verdict — 32 mutations, 32 caught, 0 blind

Method: the SOURCE of `orpheus/geometry/transformation.py` is mutated (copy to
`$CLAUDE_JOB_DIR/tmp`, mutate in place, restore from the copy — **never**
`git checkout`/`restore`/`stash`) and the gate module is run in a SUBPROCESS,
so the mutation reaches the test process. Harness:
`…/tmp/mutate.py`.

| # | mutation | gates reddened | credited gate reds? |
|---|---|---|---|
| M01 | composition: swap the operands of the translation rule | 12 | **A1** ✓ (both tiers) |
| M02 | identity carries a tiny translation | 10 | **A2** ✓ |
| M03 | inverse drops the conjugation of `t` | 4 | **A3** ✓ |
| M04 | closure never adds the generators themselves | 4 | **A5** ✓ |
| M05 | closure key drops the `-0.0` canonicalisation | 1 | **A5** ✓ (isolated) |
| M06 | infinite-group error back to a bare `ValueError` | 1 | **A6** ✓ (isolated) |
| M07 | orthogonality guard disabled | 1 | **B2a** ✓ (isolated) |
| M08 | Householder missing the factor 2 | 35 | **B2b, D1** ✓ |
| M09 | reflection skips normalising `n̂` | 14 | **B2b** ✓ |
| M10 | determinant sign flipped | 7 | **B3, B4, C5** ✓ |
| M11 | `element_order` off-by-one (starts at k=2) | 4 | **C1, C2** ✓ |
| M12 | `element_order` reports the LINEAR part's order | 2 | **X3** ✓ |
| M13 | offset factor `2δ → δ` | 5 | **C3, F2, F3** ✓ |
| M14 | offset sign flipped | 5 | **C3, F2, F3** ✓ |
| M15 | `fixed_subspace_dimension` returns the RANK | 7 | **C4, C5** ✓ |
| M16 | rotation sense reversed | 4 | **D2, D3, D4** ✓ |
| M17 | conjugation direction reversed | 8 | **E1, E2, E3** ✓ |
| M18 | seating conjugates by the NEGATIVE centre | 6 | **F1, F2, O8** ✓ |
| M19 | point action transposed (`x @ Q` for `x @ Qᵀ`) | 5 | **O2** ✓ — see below |
| M20 | `Permutation` composition order swapped | 1 | **O2** ✓ (isolated) |
| M21 | bijectivity guard deleted (ERR-073) | 1 | **O4** ✓ (isolated) |
| M22 | `Permutation` bijectivity assertion deleted | 1 | **O4** ✓ (isolated) |
| M23 | match-window guard deleted | 2 | **O5** ✓ |
| M24 | weight guard deleted | 1 | **O6** ✓ (isolated) |
| M25 | `on_directions` applies the translation | 2 | **X1** ✓ |
| M26 | circle point round-trips through an angle | 1 | **X2** ✓ (isolated) |
| M27 | embedding pads the diagonal with 0 | 2 | **H2** ✓ |
| M28 | embedding drops the translation | 1 | **H2** ✓ (isolated) |
| M29 | embedding permits lowering the dimension | 1 | **H2** ✓ (isolated) |
| M30 | composition dimension guard deleted | 1 | **H3** ✓ (isolated) |
| M31 | d=1 rotation guard deleted (the `SO(1)` message) | 1 | **H1** ✓ (isolated) |
| — | ⭐ **POSITIVE CONTROL**: `reflection` returns `+I` | **21** | reds A5, B3, B4, C1, C3, C4, D1, D4, D5, E2, F2, F3, H1, H2, O1, O3, O5, O6, O7, O8, X3 |

Twelve mutations redden **exactly one** gate — the coverage is isolating, not
diffuse, so each of those gates is that mutation's sole committed catcher.

### ⭐ The plan's §5.2 prediction, confirmed empirically

**M19 (the transposed point action) does NOT redden O1 — only O2 does.**

O1 asserts `g.on_points(points) == points[π]`, and `π` is *built by the same
mutated `on_points`*: the transposed action still maps the fixture onto
itself (`Qᵀ` is also in the group), so `permutes` still returns *a* valid
permutation — the wrong one — and O1 is self-consistent and green. The
homomorphism `π(g∘h) = π(g)∘π(h)` is what catches it, because a transposed
action is an ANTI-homomorphism.

That is exactly the claim O2's docstring makes — "it simultaneously pins the
composition order, the row-versus-column action convention, and the
``π``-versus-``π⁻¹`` choice; none of the three is checkable alone" — measured
rather than argued.

### `catches("ERR-073")` — mutation-verified, not topic-tagged

M21 re-introduces the EXACT documented bug (delete `np.unique(π).size != n`)
and reddens **O4 and only O4**. The marker is earned. Recorded substitution:
with the guard gone, the `Permutation` constructor raises rather than
`permutes` returning `None`, so O4 reds with a `ValueError` instead of an
`AssertionError` — still a red, but for a substituted mechanism, and the
docstring says so.

### Three defects the gates found — in the TEST fixtures, not production

1. **Single-pass Gram–Schmidt is not orthonormal enough for the shipped
   `1e-12` gate.** MEASURED `max|GGᵀ − I| = 2.1e-12` on near-parallel draws,
   so `rotation_from_circle_point` (correctly) refused the fixture. The
   refusal was right; the fixture was wrong. Now built by QR.
2. **An ABSOLUTE tolerance is the wrong contract for a translation residual.**
   It scales as `O(ops × ‖t‖ × eps)`, so a fixed `atol=1e-14` reds on correct
   code for large `‖t‖` draws (measured up to 2.1e-14). Re-measured
   scale-normalised over 3000–4000 draws per law at d = 1, 2, 3 — worst
   6.9e-15 — and the gates now use `_assert_close`, which scales the bound by
   `max(1, ‖desired‖_∞)`.
3. **`on_points − on_directions == t` is NOT a float theorem** (`fl(a+t) − a ≠ t`;
   measured 2.2e-16), while **`on_points == on_directions + t` IS bit-exact
   6000/6000** — it recomputes the same expression. The gate is stated in the
   direction that is true, which is also the stronger one.

All three are instances of **L35h**: bit-exactness is earned per law; measure
before choosing the assertion.

### The harness lied first — L34d, again

The battery's first run reported **32/32 BLIND** while the summary lines
plainly showed `23 failed`, `63 failed`, `38 failed`. The parser was reading
`FAILED` lines that `-q --tb=no` never emits (it needs `-rf`), and ANSI colour
codes were breaking the match. **The positive control is what exposed it**: a
mutation that makes `reflection` return the identity cannot possibly leave 42
gates green, so "0 gates" was a harness verdict, not a suite verdict. An
all-blind result remains far more likely to be a broken harness than a
perfect suite — and it costs one run to find out.

### Nothing was weakened

No tolerance was loosened to make a gate pass, no gate was deleted, and no
`# type: ignore` was added (the four pyright errors were fixed principledly:
`ArrayLike` on the helper's real contract, and `Matrix.applyfunc` instead of
`sp.simplify`, which erases the `Matrix` type). One local fixture helper was
renamed `_reflections → _reflection_family` so the forbidden-symbol grep in
§11's acceptance criteria stays a MECHANICAL check.

### Still open (unchanged by this work)

**A6** (a typed failure reason for `permutes`) is deliberately deferred by the
coordinator — no consumer, no non-identity morphism, so it fails the project's
type-minting criterion. §5.3's input-side isolators are therefore **permanent**
rather than a stopgap, and O4 carries the isolation assertion
`(guard1, guard2, guard3) == (48, 0, 48)` for that reason.
