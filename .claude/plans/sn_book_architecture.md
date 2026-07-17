# The S_N book architecture — invariant-first, broadening progression

**Status: CONFIRMED (user, 2026-07-16).** Supersedes the topic-based chapter table in
`documentation_corpus_architecture.md` §5 for the S_N book. Derived from user directives
(2026-07-16) captured verbatim in §0. This is the conceptual spine; `sn_split_catalog.md` (the
mechanical execution driver) is subordinate and re-scoped to it.

**Confirmed decisions (user, 2026-07-16):** (1) Proceed on this architecture; accept the
expanded scope (re-architecture touching `foundations/` + a connective rewrite, NOT a mechanical
split). (2) **Build `foundations/discretization.rst` FIRST** — the reusable keystone; SN chapters
cross-link it, never re-derive.

## The quality bar (governs every page) — [[feedback-articulation-lossless-disassembly]]

**Articulate = disassemble a concept so a fresh reader/LLM reassembles it LOSSLESSLY; the loss
IS the measure** (not eloquence). The docs are the LLM's **context-preservation** (insights the
user "dug from" the model, preserved so they need not be re-prompted — Cardinal Rule 3 made
literal) AND teaching material. **Target: the best neutron-transport resource ever written**,
beginner → world-class, Feynman-level articulation, polished to near-perfect. This is an
explicit licence for **maximum effort on articulation** — never trade it for brevity/speed.
Practical test per page: "if a fresh LLM read ONLY this, could it losslessly reassemble the
concept?" Get the **disassembly (outline) right before prose**; eliminate duplication/noise (that
is what raises signal); pose every concept **by its invariant** (§1).

**Methodological shift this forces:** the monolith is organized by *topic* (all balance
together, all sweep together). This design is organized by *dimension/broadening* and by
*shared-concept-then-extension*. So the remaining core content is **gathered and re-sequenced +
connective-rewritten**, NOT split by pure contiguous cuts (that was Tier-1). It also **touches
`foundations/`** (a new discretization chapter; operator intrinsic-type framing), not only SN.
The four Tier-1 chapters already extracted (`angular_quadrature`, `boundary_conditions`,
`verification`, `history`) survive untouched as shared/reference chapters.

---

## 0. The user directives (verbatim — load-bearing)

**On the progression + discretization + the invariant (2026-07-16):**
> Plan how you're going to present the move from continuous to semi-discrete to fully discrete,
> at some point you have to bring in the concept of the Step, DD and LD (i.e., upwind, box,
> linear upwind). Plan how you're going to demonstrate this for the spatial dimension, and then
> for the angular dimension (because this is essentially the cell balance). Consider making a
> chapter inside transport about discretization, with step, DD and LD as examples … it's worth
> showing it once where you don't have to rederive it again. You can cross-link that back into
> SN and show how it is used not only in space but in angle. Lean on the basic requirement for
> posing a math problem, which is finding the invariant … For the steady state transport
> equation you have sinks = sources, which is an invariant. Always find the invariant, because
> if we can nail "what is the invariant?" then we can pose the problem, and even re-pose the
> problem. For example, I believe right now our radial characteristics is particular to DD, but
> knowing what is the invariant we're representing with Radial Characteristics, we can find how
> to pose the same condition for Step or LD.

**On the operator algebra (shared vs method-extension) (2026-07-16):**
> A shared operator algebra should show the shared operators and their properties (collision as
> multiplication operator, fission as integral kernel, scattering and projection). Transport
> methods should extend the explanation when appropriate (how methods with angular flux expand
> the scattering operator and represent anisotropy). Each transport method has to show its
> invertible operator (and why THAT operator is now invertible. What does it have that makes it
> so? Then the Matrix representation which is the rawest representation of the invertible
> operator. How this invertible operator becomes an operator that encodes strategies
> (SweepOperator, etc).

**Earlier content-contract directive (2026-07-16, in `sn_split_catalog.md`):** each method
shows general→method-specific derivation + the loss-representation matrix shape (matrix =
minimal baseline) + every non-matrix rep's strategy to avoid the matrix + the operator-algebra
ops that pose the problem into the strategy-enabling shape.

---

## 1. Master principle — pose by the invariant

The organizing idea beneath everything: **a math problem is posed by identifying what it holds
invariant.** For the steady-state transport equation the invariant is **sinks = sources** (the
balance / continuity of neutrons). Naming the invariant is the act that lets us *pose* the
equation — and, crucially, *re-pose* it under a different discretization or geometry.

**Operational corollary (the generalization lever):** a construct that today is tied to one
scheme is generalized by naming the invariant it enforces, then re-posing that invariant for
the other schemes. Worked target: **radial characteristics** is currently particular to DD;
identify the invariant it represents (the curvilinear angular-redistribution balance), and the
same condition can be posed for Step and LD. Every chapter states its invariant first.

---

## 2. Shared foundations (shown ONCE, method-invariant)

### 2.1 Discretization (NEW — `foundations/discretization.rst`, or under `transport/`)

The continuous → semi-discrete → fully-discrete pipeline, derived generically:
- **continuous**: the transport operator on the continuous phase space.
- **semi-discrete**: discretize one axis (e.g. angle → ordinates), leaving others continuous.
- **fully-discrete**: close the cell balance on every axis.

The closure **schemes**, revealed under their standard numerical names (so SN never re-derives
them):
- **Step = upwind**
- **Diamond Difference (DD) = central / box**
- **Linear Discontinuous (LD) = linear upwind**

Each scheme is a way to close the **cell balance** while preserving its invariant, and the cell
balance is **dimension-agnostic** — the same closure applies whether the "cell" is a spatial
interval or an angular interval. This chapter is the single home for Step/DD/LD; SN cross-links
it (§3.3).

### 2.2 The operator algebra (`foundations/operator_algebra.rst`, EXISTS — reframe to intrinsic types)

The shared, method-invariant operators with their **intrinsic mathematical types + properties**:
- **C (collision) = a multiplication operator** (diagonal; pointwise `Σ_t·`).
- **F (fission) = an integral kernel** (`χ ⊗ νΣ_f`, rank-structured).
- **S (scattering) = a projection** onto the retained angular/energy moments.
- **L (streaming), B (boundary)** — as already in the honest algebra `A = L + C − S − B`.

Methods **extend** these when appropriate — e.g. how a method with angular flux (SN) expands
the scattering operator to represent **anisotropy** (P_N moments). The shared page states the
type + properties; the method page adds the expansion.

---

## 3. S_N realizes and extends the shared concepts

### 3.1 S_N's extension of the shared operators
How SN represents **S** with the angular flux and anisotropy (the P_N moment expansion,
`S = R∘Λ∘M` via Funk–Hecke) — the method-specific extension of the shared scattering/projection.

### 3.2 S_N's invertible operator (the load-bearing realization)
For each posing, show the chain the user specified:
1. **The invariant** (sinks = sources on the cell) → **pose** the balance.
2. **The invertible operator** `A = L + C − S` (within-group `L+C`), and **why it is
   invertible**: the **lower-triangular structure** in sweep order — *that structural property
   is what makes it invert*.
3. **The matrix representation** — the **rawest** representation of the invertible operator.
4. **How the matrix becomes a strategy-encoding operator** — the `SweepOperator` (the
   triangular structure *is* the sweep), and the other strategies (Krylov) as alternative
   encodings that avoid materializing the matrix.

### 3.3 Discretization used in BOTH space and angle (the unification)
The SN cell balance is closed by **Step/DD/LD in space** AND **in angle** (the angular
redistribution is the same closure applied to the angular cell). Cross-link to §2.1; do not
re-derive. This is the payoff of the shared discretization chapter: one concept, two dimensions.

### 3.4 The broadening progression (the book's spine)

**Part A · Cartesian S_N** — pose the whole problem, then broaden one operator at a time:
| Chapter | Broadens | Gathers (from today's topic-sections) |
|---|---|---|
| `slab_one_group` | — (the base) | 1-D slab transport eqn · invariant/balance · DD/WDD closure (via §2.1) · `(L+C)` lower-triangular · the sweep · P₀ within-group scatter · source iteration |
| `slab_multigroup` | **energy** | group-to-group + P_N (anisotropy, §3.1) scattering · fission/χ · k-eigenvalue · power iteration |
| `cartesian_multid` | **space** | streaming → gradient · 2-D/3-D balance · UBLD (LD closure, §2.1) · wavefront/octant sweep · unified dispatch |

**Part B · Curvilinear S_N** — reuse the slab, augment with the curvilinear invariant:
| Chapter | Adds | Gathers |
|---|---|---|
| `curvilinear_one_group` | curvilinear posing | curvilinear balance · **angular redistribution (re-derived from its invariant so it works for Step/DD/LD, not just DD)** · α-dome · geometry factor · radial characteristics · curvilinear closure · pole/ψ½ start · sequential sweep |
| `curvilinear_multigroup` | **energy** | curvilinear + multigroup |
| `curvilinear_numerics`? | (refinements) | Carlson coupled-pole · ERR-058 · route-(a) · the anisotropic-MMS floor |

### 3.5 Shared / reference chapters (orthogonal to the progression)
`angular_quadrature` ✅ · `boundary_conditions` ✅ · `loss_representation` (the representation
catalog the progression links to) · `verification` ✅ · `history` ✅. (✅ = Tier-1, extracted.)

---

## 4. Reconciliation + open questions

- **Supersedes** §5's topic chapter table (algebra/scattering/discrete_balance/spatial_closures/
  sweep_1d/sweep_multid/curvilinear as separate topic chapters). The *content* still moves; the
  *destination* is the progression, not topic pages.
- **Tier-1 survives** unchanged (shared/reference chapters).
- **Touches `foundations/`**: the new discretization chapter (§2.1) + the operator intrinsic-type
  reframe (§2.2). These are shared across ALL methods (CP/MoC/diffusion), not SN-only.
- **Open for the user:**
  1. Home of the discretization chapter — `foundations/discretization.rst` vs a `transport/`
     page. (User said "inside transport".)
  2. `curvilinear_numerics` as its own chapter vs folded into `curvilinear_one_group`.
  3. Multi-D curvilinear — a chapter, or out of scope for now?
  4. Chapter file names (the progression names above are provisional).
  5. Execution granularity: this is a reorganization + connective rewrite; do we drive it
     chapter-by-chapter (main-agent authored, elegance-reviewed), given it is now writing, not
     the mechanical archivist cut?
