# Boundary-condition machinery — complete review

> ## ⏸ COMPACTION POINT — cold-pickup anchor (2026-08-01, fourth rewrite)
>
> **B0, B1, B2, B3.0–B3.4c and the snapshot re-anchoring are ALL landed.
> B3.4 is COMPLETE — every law is narrowed. NEXT = #325, then B4.**
>
> ⚠ This block is rewritten at every compaction. An anchor still saying
> "NEXT = X" after X has landed **lies forward** and costs a session — it has
> now done so twice in this campaign. Trust `git log`, never this file.
>
> **Verify against the tree, in this order:**
> 1. `git log --oneline main..HEAD` on `refactor/operator-strategy-layers`. The
>    first 12 commits are the operator-strategy **P0** phase (a different
>    campaign sharing this branch —
>    `.claude/plans/operator_strategy_realization_campaign.md`); the rest are
>    this boundary campaign.
> 2. `git status --porcelain` — clean EXCEPT
>    `.claude/skills/vv-principles/{SKILL.md,error_catalog.md}`, which are
>    **forbidden to commit** by standing policy and exist ONLY in the working
>    tree. A `git checkout`/`restore`/`stash`/`clean` on those paths destroys
>    them irreversibly.
> 3. Gates: full tree `python -O -m pytest tests -m "not slow"` (**~54 min** —
>    budget for it; the ~13–16 min figure in older notes is a SUBSET);
>    `npx pyright orpheus/` → **1** (the ratchet floor, the accepted #288
>    residual — NOT a regression); `sphinx -E -W` → 0 warnings (`-E` is
>    REQUIRED on any label move — lesson L36).
>
> **Design of record for ALL of B3** is `.claude/plans/
> b3_domain_narrowing_crosswalk.md`. Its **§17** is B3.4c as executed and is
> the densest section in the file — read it before touching the boundary
> composite.
>
> | phase | commit | what landed |
> |---|---|---|
> | B3.0 | `9e2139b4` | G is the deck transformation, R the kernel; response tier minted |
> | B3.1 | `b39502f8` | `TraceRestrictionOperator` — the trace map γ± |
> | B3.2 | `7f02de15` | the SN law's domain narrowed to Γ₊ (vacuum, reflective) |
> | — | `b11a2ce3` | doc repair for every claim B3.0–B3.2 falsified |
> | B3.4a | `91f73141` | white + prescribed inflow narrow; TWO mechanisms dissolve |
> | B3.4b | `943b37c1` | albedo's re-emission closure; the specular pairing lives in R |
> | — | `87a65967` | BC-equivalence snapshots re-anchored to derived references |
> | B3.4c | see `git log` | periodic reads the PARTNER face; `B` is block-STRUCTURED |
>
> **What B3.4c changed that a fresh session will trip over:**
> * **`B` is NO LONGER block-diagonal over faces.** It is block-STRUCTURED —
>   constitutive laws on the diagonal, a quotient law (periodic) off it.
>   `SNBoundaryOperator._face_domains` is the `(row, column)` block index and is
>   certified a **permutation of the faces**. Roughly eight docstrings used to
>   assert block-diagonality; the `faces=` subset restriction never actually
>   depended on it (it filters OUTPUT faces while the whole input trace stays in
>   scope), so the conclusion survived and only the stated reason was wrong.
> * **`PeriodicWrapOperator` is RETIRED** (user ruling). Periodic realizes to a
>   bare `IdentityOperator`; the crossing lives in the CHANNEL. Do not
>   resurrect it — a shifted or rotational gluing is **#178
>   `SymmetryBoundary`**, a different object.
> * **`_face_domains` IS the face-level trace digraph**, which is what **#324**
>   (the SCC criterion) needs and never had. B3.4c unblocks it. Periodic still
>   lags into `B_upper`; `SpatialWrap.permutes_ordinates` stays `False`, and
>   that is what keeps the forward substitution triangular
>   (`tests/geometry/test_reemission_closure.py:1134-1141` is the tripwire).
> * **`periodic` is still NOT in `SNMesh.BOUNDARY_OPERATOR_REGISTRY`** (#189),
>   so `BC("periodic")` refuses at parse. The law is installed via
>   `sn_mesh.realize_boundary_law(...)`, the same production arm the tag path
>   reaches one step later. Registering it now needs only the tag-spelling
>   decision.
>
> **Then #325** (symmetry-exact node generation, crosswalk §14 — user-ruled to
> land before B4), **then B4**. Also open: **B3.3** (retire
> `IncomingOrdinateMaskTensor`) and **B3.5** (re-pose the C-1 gates, promote the
> mutation harness).
>
> **B2's result block (§B2.2) is still required reading before B4** — it records
> the two findings B4 must close (the unscaled corner swap; the `R ≠ 1` leakage
> predicate) and the family-vs-value trap that bit twice.
>
> **Rulings NOT to re-litigate** (carried forward across every compaction):
> * rank-N composition is **WIRED, not retired** (user, 2026-07-30).
> * `R` is the **response factor**; `B` is the realized **composite**. Do not let
>   the two collide again.
> * `kind`'s `"partial"` return is deliberate; B2 dissolved the question.
> * **No `as_operator` on the factors until B4 brings its consumer** — minting a
>   method with no caller and no test is the dead-capability pattern this campaign
>   exists to remove.
> * `G` carries the crossing `Γ₊ → Γ₋`; `R : Γ₋ → Γ₋` (B3.0). Membership is
>   decided by multiplicativity **plus** "is it the deck transformation of an
>   actual quotient of the physical domain?" — the second half is NOT optional
>   (user, 2026-08-01): a specular *kernel* is multiplicative too, so
>   multiplicativity alone cannot separate a polished wall from a symmetry plane.
>   ⇒ **exactly one of `G`, `R` is non-trivial.**
> * The "exactly one" law has **TWO** deliberate violators, and a gate scoped to
>   tag-reachable laws silently drops the first: `ReflectiveBoundary(axis,
>   albedo < 1)` (both factors non-trivial — a symmetry plane cannot absorb; that
>   object is `AlbedoBoundary(α, SpecularReturn)` wearing the geometry costume;
>   unreachable from a tag because `_law_from_tag` hard-codes `albedo=1.0`) and
>   bare `AlbedoBoundary(1.0)` (**ZERO** non-trivial factors — `G = I` and
>   `R = ScalarResponse(1.0) = I`). Retiring reflective's `albedo` is **B5**.
> * **`AlbedoBoundary`'s realized form is decided by its CLOSURE, not its class**
>   (B3.4b). Any predicate answering "adjointable / permutes / has a ruled
>   corner" by class NAME is wrong; read the factors. `law_permutes_ordinates`
>   is the one predicate for the permutation question and asks BOTH tiers.
> * **The geometry tier answers "whose Γ₊ does this law consume?"** (B3.4c,
>   `G.domain_face`). The response tier structurally cannot: a response is
>   constitutive — a property of the surface at THIS face — so it can never
>   reach another one. Crossing faces is an act of the deck group.
> * **The snapshot generator must never import realization code** (2026-08-01).
>   It computes references from the math; the harness imports only the case
>   registry. Both directions are AST-asserted in
>   `TestTheInversionIsStructural` — a docstring cannot hold this, because one
>   import silently turns the artefacts back into recordings while every row
>   keeps passing.

**Status:** ✅ **MAPPING COMPLETE** — all four review quadrants closed and
reconciled (diffusion arm · SN composite · law layer · test teeth). Every claim
carries an evidence tag; every **[U]** has been raised or deleted; **[X]**
corrections are recorded in §9. **B0 + B1 LANDED; B2 measured, not started.**
**Objective (user, 2026-07-30):** perfect this subsystem — make it do what it
should do *and* what it promises to do. This report is the factual substrate; it
becomes an action plan only once the mapping is complete and contradiction-free.
**Branch:** `refactor/operator-strategy-layers` @ `73627b71`.

---

## 0. Evidence key — every claim in this report carries one

| tag | meaning |
|---|---|
| **[M]** | **MEASURED** — I ran it in this session and report the number. |
| **[R]** | **READ** — quoted from source at a cited `file:line`. |
| **[G]** | **GREP** — an exhaustive-search result (absence claims live here). |
| **[U]** | **UNVERIFIED** — reported by a sub-agent, not yet independently checked. Must be resolved to [M]/[R]/[G] or deleted before this becomes a plan. |
| **[X]** | **CORRECTED** — a claim previously asserted in this session that measurement disproved. Kept deliberately, with the correction, so it is not re-derived. |

**Rule for this document:** no claim may enter the plan at **[U]**. A claim that
cannot be raised to [M]/[R]/[G] is deleted, not softened.

---

## 1. The three-layer architecture, as designed

`BoundaryTraceLaw` (descriptor, no `apply`) → `BoundaryRealizer` (the seam) →
`LinearOperator` (realized, monomorphic).

**[R]** `orpheus/geometry/boundary/_base.py:98-102` — the law has its own
docstring section titled *"No `apply`"*: *"Descriptors are **not** callable. The
§16A.3 three-layer architecture demands that operatorhood live in a separate
layer."*

**[R]** `orpheus/geometry/boundary/_realizer.py` — `realize_recursively` is
*"the ONLY function that transforms a descriptor into an operator."* The law
algebra is closed (`LawSum` / `LawScaled`, never an operator) and the walk is
structure-preserving (`LawSum→OperatorSum`, `LawScaled→ScaledOperator`).

**[G]** The realized operators ARE monomorphic — `grep "singledispatch\|@overload"`
across `orpheus/sn/operators/boundary.py` returns **zero hits**; every `apply` is a
plain `def` with one carrier argument. The four names in that file (**[G]** exactly
four, no others):

| class | line | carrier |
|---|---|---|
| `SNBoundaryOperator` | `:123` | `FullField` |
| `RadialCharacteristicBoundaryOperator` | `:478` | `RadialCharacteristicField` (a *different class*, not a dispatch arm) |
| `SNMaskedBoundaryOperator` | `:715` | `FullField` (one split half; no transpose by design, `:732-735`) |
| `BoundarySplit` | `:821` | NamedTuple `(lower, upper)` |

**[R]** Diffusion: `DiffusionBoundaryOperator` is at `diffusion/operators.py:544`
(**[X]** an earlier draft of this report said `:407` — that line is inside
`LeakageOperator.apply`; `LeakageOperator` begins at `:272`). **[G]** That file
contains exactly two operator classes plus the `_FaceClosure` dataclass at `:254`.

> **This is the property the reaction side lacks** and the reason this subsystem
> was being used as its template. It holds.

---

## 2. The tensor-network question — history and present state

### 2.1 The original vision

**[R]** `orpheus/sn/boundary/angular.py:39-41` — `AngularAverageOperator`
*"Realizes the geometric projection `G_diff` in the §15.2 sum-of-tensor-products
form of the white boundary law `R_white = G_diff ⊗ α`."*

**[R]** `orpheus/sn/boundary/realizer.py:~270` (reflective arm) — *"The grand-report
§16A.10 BC decomposition is `B = G_patch ⊗ K_omega ⊗ K_g`."*

**[R]** `orpheus/sn/boundary/angular.py` (`IncomingSourceOperator`) — *"Realises
the `q` term in the §16A.1 affine BC form `γ₋ψ = R G γ₊ψ + q` for the **rank-0
case** where `R = G = 0`."* → the **rank** vocabulary is already live.

**[R] RAISED FROM [U] — the Grand Report is IN THE REPO** and tracked:
`.claude/plans/neutron_transport_grand_report_v3.md` §16A.2 (`:2776-2801`),
verbatim:

```python
class BoundaryGeometryMap(Protocol):
    def map_outgoing_to_incoming_geometry(self, outgoing_trace_point): ...
    def as_tensor(self, trace_space): ...

class BoundaryResponseKernel(Protocol):
    def apply(self, mapped_outgoing_trace): ...
    def as_tensor(self, trace_space): ...

@dataclass(frozen=True, slots=True)
class BoundaryTraceLaw:
    geometry_map: BoundaryGeometryMap
    response: BoundaryResponseKernel
    source: BoundarySource

    def as_operator(self, trace_space):
        G = self.geometry_map.as_operator(trace_space)
        R = self.response.as_operator(trace_space)
        q = self.source.as_field(trace_space.incoming)
        return AffineBoundaryOperator(linear=R @ G, source=q)
```

**This settles three things at once:**

1. **TWO factors, `R ∘ G`.** Not a three-factor sandwich.
2. **The law owns the SPEC; the `trace_space` produces the MATRIX.** The design
   anticipated the exact spec-vs-realized-operator split, via `as_tensor` /
   `as_operator(trace_space)`. **This is the direct refutation of the
   "`geometry_map` is not populatable" error recorded in §2.3.**
3. **The three members were designed as typed dataclass FIELDS**, not as
   `-> Any` properties defaulting to `None`. The present shape is a degradation
   of the design, not an unfinished sketch of it.

**[R]** §16A.5 (`:2946-2960`) dispatches on `isinstance(law.response,
ZeroBoundaryResponse)` — on the **response object**, not the law class.

**[R]** The theory page **types** the factors (`boundary_conditions.rst:210-220`):
`G : Γ₊ → Γ₊` — *"a measure-preserving permutation, pushforward, spatial
wrap-around, or hemispheric cosine-weighted average… it carries pure geometry"*;
`R : Γ₊ → Γ₋` — *"a scalar amplitude in [0,1] for the standard sub-Markov BCs…
or a full angular kernel in general weak-form BCs (deferred)"*.

> **`G` is an endomorphism of the OUTFLOW trace and `R` is the crossing.** That
> typing independently confirms §4.2: the law's domain is `Γ₊`, the outflow trace
> — which is exactly what diffusion hands it and what SN does not.

**[G]** None of the designed types exist: `AffineBoundaryOperator`,
`ZeroBoundaryResponse`, `BoundaryResponseKernel`, `BoundaryRelation` — **0 hits in
`orpheus/`**. `BoundaryGeometryMap`'s only trace is
`BoundaryGeometryMapNotMeasurePreservingError` (`_errors.py:88`) — **an error
class named after a map that was never built.**

### 2.2 Why it stalled — and the scope of that negative result

**[R]** `docs/theory/foundations/operator_tensor_network.rst` — the **MA-Q1
master condition**:

> SOTP applies ⟺ each summand factors as `f(x₁,…,x_d) = f₁(x₁) ⊗ ⋯ ⊗ f_d(x_d)`

Three of four Wave-T operators **violate** it: scattering (the per-material per-ℓ
einsum couples group to spatial via `cells_by_material`), streaming-spatial (the
WDD recurrence sequentially couples cells), streaming-angular (the Morel–Montry
half-grid recurrence).

**[R]** `:117`, `:1002` — *"**Only T.1 (BC realizers) and T.2 (fission rank-1)**"*
succeeded. **The negative result is scoped to coupled BULK physics. Boundaries
pass MA-Q1** — a reflective face permutes ordinates identically at every surface
point and every group, so the axes are genuinely disjoint.

**[R]** `:1025` cites *"§16A.10 (BC as tensor network)"* as a north-star line.

### 2.3 What is actually built today

**[R]** The SN realizer produces:

| law | realized form |
|---|---|
| reflective | `PermutationOperator(perm, axis=0) & IdentityOperator()`, ×`α` if α≠1 |
| vacuum | `IncomingOrdinateMaskTensor(axis=0) & IdentityOperator()` |
| white | `AngularAverageOperator(...) & IdentityOperator()` |

**[X] CORRECTION (mine, this session).** I first reported the tensor-network
vision as *"materialized already."* **That is too generous and is wrong.** What
materialized is the **typing**, not the structure:

- **[G]** All five `TensorProductOperator` uses in the SN realizer are
  `X & IdentityOperator()` — a **no-op second factor**.
- **[G]** `OperatorProduct` (the `∘`) is **never used by any boundary operator**.
  The `R ∘ G` composition never happens; the realizer fuses it as `α * base`.
- **[G]** `SumOfTensorProductsOperator` has **zero production consumers
  repo-wide** (outside its own definition in `numerics/operator.py`).

**[X] CORRECTION (mine, this session).** I also claimed `geometry_map` *"is not
populatable"* because `G` needs the quadrature. **Wrong**, and the error was in
what `G` is: as a **realized operator** it needs the quadrature; as a **typed
specification** (`SpecularMirror(axis)`, `HemisphericalAverage(axis, sign)`,
`ScalarAlbedo(α)`, `ZeroResponse()`) it does not. Its `-> Any` signature
returning a realized operator is *what made it un-populatable*. That is a spec
defect, not a physics wall.

### 2.4 White is rank-one, hand-written

**[R]** `orpheus/sn/boundary/angular.py` — `AngularAverageOperator.apply`:

```python
psi_avg = (self._cos_w.reshape((-1,) + (1,)*(psi.ndim-1)) * psi).sum(axis=0) / self._norm
return np.broadcast_to(psi_avg[None, ...], psi.shape).copy()
```

Contract with a `|Ω·n|`-weighted covector, then broadcast isotropically — a
rank-one operator written inline rather than typed as one.

**[M]** `RankOneOperator` was minted **18 minutes after** the BC tensor lift, for
the same structure: `fa13e78e` *"Wave T step T.1 — lift remaining BC realizers to
TensorProductOperator"* at **16:06**; `0b2848be` *"Wave T step T.2 — fission as
rank-1 TensorProductOperator"* at **16:24**. White was walked past.

**[R]** Consequence: `AngularAverageOperator` **refuses to expose
`apply_transpose`**, on the grounds that the unweighted transpose differs from the
cosine-weighted adjoint the BC analytically possesses. Typed as a genuine
rank-one `u ⊗ v`, the transpose is `v ⊗ u` — structurally free, and the metric
question becomes explicit rather than avoided.

---

## 3. The declared-but-empty properties

**[G]** `geometry_map` and `response_kernel` are declared ONLY on the base
(`_base.py:149`, `:158` — **[X]** an earlier draft cited `:157`/`:166`, which
drifted when this session's `creates_sweep_cycle` retirement removed lines above
them), return `Any`, default `None`. **No concrete law overrides either. No
production code reads either.** The only tests assert they are `None`
(`tests/geometry/test_boundary_trace_law.py:96-102`).

**[G]** Control: `source` (the `q`, `_base.py:167`) **is** overridden
(`prescribed_inflow.py:111`) and **is** read (`_base.py:300`, `:308`;
`sn/boundary/realizer.py:345`, `:349`). Exactly one of the three "first-class
properties" is real.

### 3.1 The latent consumer — five production sites, answering in strings

**[R]** All five verified at their cited lines:

| site | the string test | the structural question it is really asking |
|---|---|---|
| `sn/loss_representation/sweep_schedule.py:267` | `bc[face] == "reflective"` | **does my `G` permute ordinates?** |
| `sn/operators/boundary.py:558` | `kind in _RULED_CORNER_KINDS` | **is my `G` adjointable?** |
| `sn/operators/boundary.py:610` | `kind in _RULED_CORNER_KINDS` | **is my `G` adjointable?** |
| `sn/solver.py:1363` | `op.kind == "vacuum"` | **is my `R` zero?** |
| `sn/acceleration/dsa.py:214` | `bc_kinds - _SUPPORTED_BC` | does this law have a proven low-order row? |

This is the stringly-typed-dispatch anti-pattern standing in for the properties,
and it is the **latent consumer** that `process-discipline.md` requires be hunted
for before declaring anything unused. The hunt returns positive.

### 3.2 The discriminator against what was retired at `ddc1ee10`

`creates_sweep_cycle` was retired **this session** with the rationale "zero
production readers". That sentence is also true of `geometry_map` /
`response_kernel`, so the discriminator must be stated explicitly:

- `creates_sweep_cycle` **could not have worked in principle** — cycle-ness
  depends on the whole face configuration, which a boolean on the boundary *kind*
  cannot express. **[M]** `reflective|vacuum` acyclic vs `reflective|reflective`
  cyclic.
- `geometry_map` / `response_kernel` **can** — both are law-intrinsic, and `R`
  is *already* in use under another name (`law.albedo`; the realizer writes
  `float(law.albedo) * base`).

Same symptom, different disease. Retiring one and populating the other is
consistent.

---

## 4. The dead-capability pattern — three instances in one subsystem

**[G]/[M]** Each is a declared capability whose stated consumer does not exist:

| declared | stated consumer | reality |
|---|---|---|
| `creates_sweep_cycle` | "§15A.2 cycle detection" | **[G]** detector never built; `assert_cycles_are_declared` unimplemented. RETIRED @ `ddc1ee10`. |
| `geometry_map` / `response_kernel` | *"Concrete BCs populate"* (`_base.py:92`) | **[G]** nothing populates |
| the vacuum mask's outflow preservation | *"downstream consumers (sensitivity adjoints, post-BC field readers) require"* | **[M]** no such caller |

> This subsystem has a habit of **declaring capability ahead of consumption**.
> That is itself a finding, and an argument for caution about materializing
> anything speculatively. `geometry_map`/`response_kernel` are the exception
> *because* §3.1 shows consumption already exists.

### 4.1 The vacuum-mask measurement

**[M]** `slab_seedless()`, S4, BCs `xmax=vacuum`, `xmin=reflective`; boundary
trace seeded strictly positive on **all** rows:

| face | law alone: outflow | composite: outflow | composite: inflow |
|---|---|---|---|
| `xmax` vacuum | **1.653** (bit-preserved ✓) | **0.000** | 0.000 |
| `xmin` reflective | 1.834 | **0.000** | 1.984 |

**[X] CORRECTION** to the reported claim *"SN vacuum's composite action is the
zero map"*: `B(ψ) ≡ 0` on a vacuum face is **correct physics**, not a defect —
`B` maps outflow trace → inflow trace and a vacuum face has zero incoming flux.

The real finding is narrower and broader at once: **the outflow preservation is
discarded for EVERY law, not just vacuum.** `_reflect_trace` writes only
`sel = inflow_indices_for_face` (`boundary.py:302`), so `B(ψ)`'s outflow rows are
zero on the reflective face too. The mask projects onto the outflow subspace; the
composite projects onto its complement. `P_in ∘ P_out = 0` by construction.

**[G]** No observer exists: the sole SN caller of a realized law's `apply` is
`_reflect_trace` (`boundary.py:301`); diffusion has its own path feeding
`trace.outflow_view(face)` (`diffusion/operators.py:613`); and the stated
beneficiaries (`outflow_indices_for_face` consumers in `streaming.py`,
`loss_representation`) read the **field's** trace, never `B(ψ)`.

### 4.2 ⭐ THE ROOT CAUSE — a domain defect, with a worked fix already in-repo

**[R]** The two methods hand their law **different domains**, and that single
difference explains the whole discard:

```python
# SN — boundary.py:295,301        the WHOLE face slot
face_in = boundary.face_view(face)
full    = law.apply(face_in)
out_boundary.face_view(face)[sel] = full[sel]     # <- discards outflow rows

# diffusion — operators.py:612-614   the OUTFLOW HALF only
out_boundary.face_view(face)[ScalarTraceSpace.INFLOW_ROW] = (
    law.apply(trace.outflow_view(face))
)                                                  # <- nothing to discard
```

**SN's law is typed `FullField → FullField` when the physics is `outflow-trace →
inflow-trace`.** The mask's outflow preservation, the discard, and the
`P_in ∘ P_out = 0` composition are all consequences of that one over-wide domain.
**Diffusion picked the honest domain and the entire problem does not arise.**

> **This converts §4.1 from "dead capability" into "a typing defect with a worked
> precedent inside this repo."** The fix is a **domain narrowing** on the SN law,
> not a new abstraction. Diffusion is the reference implementation.

**[M]** Diffusion measured end-to-end (mirror of §4.1's table): composite inflow
`== law.apply(J⁺)` **bit-identically** for albedo(0.4), zero_flux, reflective;
bulk exactly zero. The "law alone: outflow" column *cannot exist* there.

---

## 5. Structural duplication — **the original claim was WRONG; the corrected
finding is sharper**

**[X] CORRECTION.** An earlier draft carried the claim that `G_in` is implemented
"three times, three ways". **Materialized as 24×24 matrices, that is false**:

**[M]** (a) `_reflect_trace`'s slice-write `out[sel] = full[sel]` (`boundary.py:302`)
and (c) `IncomingSourceOperator._inflow_mask`'s dense multiply are **the same map,
bit-identical** — `‖a − c‖∞ = 0`. (b) `IncomingOrdinateMaskTensor` is the
**complement**: `‖b − (I − a)‖∞ = 0` and `a @ b = 0`.

So it is **two spellings of `P_in` plus one spelling of `I − P_in`** — all
diagonal, all derived from the same index set. That is still duplication worth
removing, but the honest statement is different from the original.

**[R]** And the naming is inverted: the class called
**`IncomingOrdinateMaskTensor` is documented as "projection onto the *outflow*
subspace"** (`orpheus/numerics/operator.py:~2305`). The "Incoming" class IS the
outgoing projector. **[G]** No `P_in` type exists anywhere — the projector the
subsystem actually needs three times has no name.

---

## 6. Reachability constraint

**[R] RAISED FROM [U], and it is worse than stated.**

- **[R]** SN admits `{vacuum, reflective}` (`sn/mesh/augmented_mesh.py:171-174`),
  pinned by `tests/sn/operators/test_snmesh_realizer_wiring.py:403`.
- **[R]** Diffusion admits `{vacuum, reflective, albedo, zero_flux}`
  (`diffusion/augmented_mesh.py:158-163`).
- **[M]** End-to-end on real meshes: `white`, `periodic` and `prescribed_inflow`
  raise `ValueError` on **both** methods.

> **Three of the seven laws are unreachable from a `BC` declaration under EVERY
> method.** All four alternative paths were audited and closed: **[M]**
> `realize_boundary_law` has exactly one caller (`transport/method.py:259`, gated
> by the admission table); **[G]** every `WhiteBoundary(...)` / `PeriodicBoundary(...)`
> occurrence in `orpheus/` is a docstring; `BoundaryTraceLaw.create` has no
> production caller. **They are reachable only from tests.**

**[R]** `BC.white = BC("white")` is a public, tab-completable convenience instance
(`geometry/mesh.py:126`) that **[M] no method in the codebase admits**.

**[R]** A stale advertisement of the same family: `sn/mesh/augmented_mesh.py:182-185`
lists the kinds "the realizer handles today" as `white`, `periodic`, `albedo`,
`prescribed_inflow`, **`mixed`** — and **[G]** `mixed` / `MixedBoundaryOperator`
**was deleted in Wave 11**; the only surviving references are three deletion
notices.

**[R]** `ReflectiveBoundary.kind` returns `"reflective" if self.albedo == 1.0 else
"partial"` — an exact float compare (`reflective.py:107`). **[G]** `"partial"` is
not a registry key and appears in no admission table, yet `BC.to_alpha` handles it
(`geometry/mesh.py:109`).

This is the code-side statement of issue **#189**.

---

## 7. Downstream dependency

**[R]** `.claude/plans/operator_strategy_realization_campaign.md:486-490` proposes
`SNReactionRealizer` building `j_bulk ∘ E₀ ∘ F_G ∘ M₀ ∘ π_bulk` — an
inject∘response∘project sandwich — and cites the boundary realizer as the
**landed precedent**. Per §2.3 it is not one: no boundary operator composes with
`∘`. **P2b currently rests on a precedent that does not exist.** Materializing the
boundary network makes the citation true; retiring `G`/`R` means P2b needs a
different justification.

---

## 7A. The existing issue backlog — this is a CLUSTER, not a new idea

**[M]** `gh issue list --state open --limit 100`, filtered: **16 of 100** open
issues are boundary-related. Seven of them are, individually, steps of the exact
materialization this review is converging on — filed separately, never connected:

| issue | title (abridged) | what it actually is |
|---|---|---|
| **#177** | BC: `BoundaryRelation` hierarchy for weak-form / PN / Marshak | **the typed `R` hierarchy** |
| **#178** | BC: `SymmetryBoundary` distinct from physical mirror | **a `G` distinction** (quotient gluing vs physical reflection) |
| **#183** | BC: `PeriodicWrapOperator` spatial-pushforward implementation | **the missing periodic `G`** — its apply is reportedly `x.copy()` |
| **#189** | SN: wire `WhiteBoundary`/`PeriodicBoundary`/`AlbedoBoundary` into `SNMesh.BOUNDARY_OPERATOR_REGISTRY` | **the reachability gap of §6** — already filed |
| **#219** | Architecture: grand-report §7 foundational layer — `MethodSpace` ABC + builder | **the `method_space` abstraction** — partially addressed at `1fd15f64` |
| **#265** | Data-layer invariant value-objects: Quadrature laws, **albedo α∈[0,1]** | **`ScalarAlbedo` as a value object** = the typed `R` |
| **#300** | Close the boundary SCC in CLOSED FORM (**Woodbury on the rank-1 B**) | **the payoff of typing white as rank-one** |

Plus the future realizer witnesses — **#179** (MoC), **#180** (MC), **#181** (CP)
— which are exactly the consumers that make the `BoundaryRealizer[MethodSpaceT]`
Protocol generic-ness (landed `1fd15f64`) load-bearing rather than speculative;
and **#244** (duplicate BC-registry resolution `CPMesh._resolve_bc ≈
MOCMesh._resolve_bc`), a twin-path in the admission layer.

**Consequence for the plan.** The materialization is not a speculative
architecture proposal — it is the *connective tissue* of a backlog that already
exists and is already labelled. The plan should be written to CLOSE these issues,
and any step of the plan that does not map to one should be justified as new.

**[M]** Also open and relevant to correctness, not architecture: **#313**
(reflective wall-lag rate mode), **#252** (transverse face-slope sign under
reflective BC), **#293** (analytic Marshak diffusion reference), **#306**
(Wave-O O.2/O.5 residue), **#283** (CP spherical row-sum — adjacent, not boundary).

---

## 7B. Quadrant findings — diffusion arm, SN composite, law layer

### 7B.1 Diffusion is the subsystem's only fully-factored arm

**[R]** `DiffusionBoundaryRealizer.realize` (`boundary_realizer.py:176-195`) is
**two named stages**, not an isinstance chain producing an operator:
`_partial_current_albedo(law) -> float` (`:197-249`) then
`_albedo_operator(albedo) -> LinearOperator` (`:142-157`). **The
operator-producing stage is law-blind AND method-blind.** **[M]** `method_space`
is inert there — `minimal()` and `for_face(...)` give identical results.

**[M]** `𝒜` **is** `R`. `G` is degenerate for a *dimensional* reason: a diffusion
face slot is `(2, ng)` — one outflow and one inflow DOF per group, `face_shape ==
()` — so `G : ℝ¹→ℝ¹`, `GL(1) ≅ ℝ*`, and any scalar is absorbed into `R`.
**`G = I` is forced, not chosen.**

**Proof by negation [R]:** the diffusion realizer refuses exactly two laws, and
they are exactly the two needing a non-`R` factor — `PeriodicBoundary` (*"a
trace-block permutation"*, a non-identity `G`) and `PrescribedInflow` (*"the
rank-0 AFFINE law"*, a nonzero `q`). Diffusion realizes precisely the
`{G = I, q = 0}` corner of `γ₋ψ = R G γ₊ψ + q` and refuses loudly at both walls.

**Transferable rule (narrower than "copy this"):** the stage that produces the
operator must not see the law; the stage that reads the law must produce a
**specification**, not an operator. Diffusion satisfies both because its
specification is a `float`. SN's would be a pair (`G`-spec, `R`-scalar).

**[M]** Adjoint: every realized diffusion *law* is adjointable, but
`DiffusionBoundaryOperator` and `LeakageOperator` are `is_adjointable=False` and
`.H` raises. **It would be free** — the trace metric is a single scalar per face
(Cartesian `1.0`, cylinder `2πR`, sphere `4πR²`), identical on both rows, so
`‖G⁻¹BᵀG − Bᵀ‖max = 0.0` measured in all three geometries. **Missing, not blocked.**

### 7B.2 SN composite — the adjoint promise is structurally impossible

**[M]** `_RULED_CORNER_KINDS = {"vacuum","reflective"}` (`boundary.py:100`) is
**identical to SN's admitted set** ⟹ the corner `NotImplementedError` is
**unreachable from any `BC(...)`**.

**[M]** `SNBoundaryOperator`'s docstring calls `.H` *"the one channel by which the
white-BC adjoint becomes available"* (`:152-154`). With a white face installed,
`B.is_adjointable == False` and **`B.H` raises `MissingAdjoint`** — `.adjoint()`
gates on the same predicate. **The promise is not merely unfulfilled; it is
structurally impossible as written.**

**[M]** For reachable laws the metric is inert: `‖B.H − Bᵀ‖∞ = 2.8e-17` (slab),
`4.4e-16` (2-D) — though the trace metric IS populated (`|μ_n|·w_n`).

**[M]** The `(P_sel ∘ law)ᵀ = lawᵀ ∘ P_sel` identity **is correct** — materialized,
`‖T − Fᵀ‖∞ = 0` exactly, on both fixtures. (`PermutationOperator.apply_transpose`
uses `argsort(perm)`, a genuine transpose.) This is sound and should not be
"fixed".

**[M]** `SNMaskedBoundaryOperator` **costs the adjoint**: `(L+C) − B_lower` is
`is_adjointable=False` while Jacobi `(L+C) − B` is `True`. **The boundary
Gauss-Seidel path is structurally adjoint-less** — relevant to #276.

**[M]** The split is exact: per-face `union == inflow`, `intersect == ∅`,
`‖B − (lower+upper)‖∞ = 0`, `nnz 216 = 168 + 48`.

**[M] A latent bug:** the stringly-typed `method` tag is unvalidated —
`_reflect_trace(bf, "banana")` silently performs the **transpose**.

### 7B.3 Law layer — the invariant suite is thinner than advertised

**[M]** Four of the five universal invariants have **zero statements**; only
`assert_source_lives_on_incoming_trace` has a base body. **Invariants 1 and 2 are
overridden by NOBODY** — permanent no-ops. **Four of seven laws override nothing.**

**[G]** `assert_realizable` has **exactly one** production caller
(`sn/boundary/realizer.py:194-197`). **Diffusion never calls it** — so
`WhiteBoundary(albedo=5.0)` through diffusion yields `𝒜 = 5.0` with no sub-Markov
error, directly contradicting `_base.py:328-333`'s *"every law arrives at its
primitive construction already certified."*

**[M] `law.kind` is ABSENT on four of seven laws** — measured:
`White`, `Albedo`, `Periodic`, `ZeroFlux` have no `.kind` at all. That is why
production reads it as `getattr(law, "kind", None)` (`boundary.py:558`). **The
five string-dispatch sites of §3.1 are dispatching on an attribute that does not
uniformly exist.**

**[M]** `VacuumInflow(kind="banana")` **constructs** — `kind` is a dataclass field
there but a property elsewhere. Illegal state, representable.

**[R]** `ZeroFluxBoundary` **is not an affine trace law**: `φ_Γ = 0` is
`A₋γ₋ + A₊γ₊ = 0`, a `BoundaryRelation` in the grand report's taxonomy — a tier
ABOVE `AffineTraceLaw`. **[G]** The code ships no such tier, which is *why* SN's
refusal must be a hand-written `isinstance` guard instead of a type distinction.
**This is issue #177 restated from the code side.**

**[R]** Per-law `R`/`G` derivability: **every `R` in the family is a scalar**, and
both realizers already multiply by it as a bare float
(`sn/boundary/realizer.py:288`, `diffusion/boundary_realizer.py:214`). `G` is
derivable as a **typed spec for six of seven**; the lone blocker is
`PeriodicBoundary`, which **[R] has no fields at all** — its defining parameter
(the partner face) does not exist on the law. That is a *missing field*, not a
discretization problem, and it is issue **#183**.

### 7B.4 ⭐ Rank-N composition has ZERO production consumers

**[G]/[M]** `realize_recursively`'s only non-docstring callers are its **own three
recursive calls** (`_realizer.py:307`, `:312`, `:313`). Every other occurrence is
prose or a docstring example — including `diffusion/boundary_realizer.py:93-94`,
**[M]** verified to sit inside that module's 114-line docstring.

So the entire rank-N path — `LawSum`, `LawScaled`, `LawNode`,
`realize_recursively`, and the eight algebra dunders on the ABC — is exercised
**only by tests**. And its flagship example is unreachable *in principle*:
`0.3·Reflective + 0.7·White` cannot be declared because **no method admits
`white`**.

> **This is the largest dead-capability instance in the subsystem by code volume**,
> and it must be decided rather than left half-alive: wire it (which requires
> admitting `white` first, #189) or retire it.

**[R] One gap in the closure, in two places:** the descriptor algebra has **no
product** — `R ∘ G` is not expressible at the descriptor layer, mirroring **[G]**
`OperatorProduct` being unused at the operator layer. Populating `R` and `G`
without giving the descriptor layer a compose reproduces the same half-measure one
level down.

---

## 8. NOT YET REVIEWED — the gaps this report must close

1. ✅ **CLOSED — the diffusion arm** → §7B.1. Detail `scratch/review_diffusion_arm.md`.
2. ✅ **CLOSED — the SN composite machinery** → §7B.2. Detail `scratch/review_sn_composite.md`.
3. ✅ **CLOSED — the law layer** → §7B.3, §7B.4. Detail `scratch/review_law_layer.md`.
4. ✅ **CLOSED — test coverage and its teeth** → §7C. Detail
   `scratch/review_boundary_test_teeth.md`.

---

## 7C. Test teeth — my `-O` premise was WRONG, and the suite is strong

### 7C.1 **[X] REFUTED — the Mode-8 `-O` exposure does not exist here**

I asserted twice this session that `tests/geometry/`'s bare `assert`s are stripped
under the canonical `python -O` invocation, and called it "the number I most want".

**[M] The inert-assertion fraction is `0 / 676` = 0.0 %.**

**Mechanism:** pytest's **assertion rewriter** transforms `assert X` into
`if not X: raise AssertionError(...)` in every *collected* test module at import,
before the interpreter's `-O` stripping applies. `-O` only strips asserts the
rewriter never touched. **[M]** Proven twice — a synthetic control, and falsified
copies of `test_bc_universal_invariants.py:399` and `test_law_composition.py:69`
producing **byte-identical `2 failed, 63 passed`** with and without `-O`.

**[M]** The single inert assert in the whole audit is in
`tests/geometry/_generate_bc_equivalence_snapshots.py:339` — a **non-collected**
manual generator.

> **Refined rule (supersedes the loose reading of `vv-principles` Mode 8):** bare
> asserts in **collected test modules** are safe under `-O` because of the
> rewriter. Bare asserts in **non-collected** modules — helpers, generators,
> support code, production — are stripped. The hazard is real but its scope is
> *non-collected code*, not test modules.

### 7C.2 The suite is the most mutation-live in the repo

**[M]** 31 deliberate defects injected (12 leaf-action, 14 guard-disabling, 5
adjoint), each with a bite check: **30 caught.** Every law invariant and every
realizer refusal reddens exactly its own negative test.

**The real weakness is assertion CONTENT, not execution:** **[M]** only **29.1 %**
of bare asserts pin a *value*; 44.7 % are `isinstance`/`is` structural, and
**11.1 % are `== "reflective"` tag-equality** — the test-side shadow of §3.1's
stringly-typed dispatch.

### 7C.3 Five measured coverage defects

1. **`catches("ERR-052")` is a phantom edge.** **[M]** Bite-verified: re-introducing
   the bug (renorm_calls 6→0, |φ|max 7.60→0.61) leaves the test **green**. It was
   a true catcher when written; the config drifted out of the regime (converges in
   6 outers, the bug needs 30–60). **`catches` markers DECAY.**
2. **Three permanently-inert sentinels.** `tests/sn/operators/test_bc_extraction_matvec.py:445`
   — `except Exception as exc: pytest.skip(f"…{exc}")` swallows an `IndexError`
   (1-D `Mesh1D` then `spatial_shape[1]`). The self-described "SENTINEL" has
   **never run**; a future construction bug is a green skip.
3. **`test_bc_errors.py` — 9 tautological legs**: `with pytest.raises(X): raise err`
   where `err` *is* an `X`. **[M]** Zero of 14 guard mutations touched the file.
4. **Two white tests documented "hand-computed" are blind** — **[M]** both PASS
   with the `|Ω·n|` factor dropped (hemisphere-constant input ⇒ Mode-12
   invariance). The real catchers are the cosine-weighted-average, conservation,
   and self-adjointness tests.
5. **`test_bc_equivalence_snapshot.py` is marked `l1`** while its own header records
   that the cross-implementation half was deleted — a **self-generated regression
   baseline wearing an L1 label**. Still the widest net (9/12 leaf mutations), and
   `_load_or_skip` **skips** on a missing snapshot.

### 7C.4 ⚠ TWO ANSWERS THE PLAN MUST OBEY

- **§4.2's domain narrowing is NOT free.** **[M]** The outflow discard is an
  **asserted contract**: making `B` preserve outflow rows reddens **4 tests**
  demanding they be zero (*"B emitted non-zero on the outflow rows"*) plus the 2-D
  balance gate. The narrowing is still right, but it is a **contract change
  requiring those gates be re-posed**, not a silent fix.
- **The white-transpose refusal is a TWO-part contract.** **[M]** `is_adjointable`
  is a *declared* predicate defaulting `False`; adding `apply_transpose` alone
  changes nothing. **A plan step that only adds the method is a no-op.**
- **[M]** Do not cite the keff physics gate as reflection-table coverage: it
  catches reflection→identity but is **blind to a rolled permutation** (a
  permutation preserves the balance functional).

### 8.1 Surfaced, not yet dispositioned

- **[M]** Importing `tests/geometry/test_boundary_trace_law.py` inserts
  `_stub_for_test` into the **production** `BoundaryTraceLaw.registry`, no teardown.
- **[M]** Four conflicting V&V level markers (`['foundation','l0']`) at
  `tests/diffusion/test_boundary_realizer.py:217, 228, 242, 254`.
- **[R]** `docs/theory/foundations/boundary_conditions.rst:2833-2854` shows the
  walker with the **pre-P7b 2-arg signature** and a hardcoded `SNBoundaryRealizer()`;
  `docs/theory/methods/sn/boundary_conditions.rst:240` repeats it in a user-facing
  snippet that would `TypeError`. The correct 3-arg call appears 20 lines later on
  the same page.
- **[M]** ~18 further "promises not kept" catalogued across the three quadrant
  files — five deprecated aliases advertised as *"re-exported below"* that are not
  importable; a theory page contradicting itself on the ABC's MRO; a documented
  invariant (`assert_direction_norm_preserved`) that **[G]** does not exist.

### 8.2 The original item 4, for reference — what is still being measured

   what is actually gated vs asserted-only;
   in particular the **Mode-8 `-O` exposure**: `tests/geometry/` uses bare
   `assert`, and the canonical invocation is `python -O -m pytest`, which strips
   them. If that is so across the boundary suite, an unknown fraction of this
   subsystem's "coverage" asserts nothing.

---

## 9. Contradictions resolved so far

| claim | resolution |
|---|---|
| "the tensor-network vision did not materialize" (user hypothesis) | **Partly wrong**: it materialized as *typing* (`X & I`), not as *structure*. `∘` is never used. §2.3 |
| "the vision materialized already" (mine) | **[X]** too generous; same correction. §2.3 |
| "`geometry_map` is not populatable" (mine) | **[X]** wrong — true only for a *realized operator*, false for a *typed specification*. §2.3 |
| "SN vacuum's composite action is the zero map" (sub-agent) | **[X]** factually true but mis-framed: that is correct physics. The real finding is the universal discard of outflow rows. §4.1 |
| "these have zero readers, like `creates_sweep_cycle`" | **False equivalence** — §3.2 gives the discriminator. |
| "`G_in` is implemented three times, three ways" (sub-agent) | **[X]** Materialized: it is **two spellings of `P_in` + one of `I − P_in`**. §5 |
| "the boundary test suite is inert under `-O`" (mine, twice) | **[X] REFUTED — 0/676.** pytest's assertion rewriter fires before `-O`. The hazard's real scope is *non-collected* modules. §7C.1 |
| "`BoundaryGeometryMap` has 0 hits in `orpheus/`" (sub-agent) | **Substantively right, grep imprecise**: the 3 hits are all `BoundaryGeometryMapNotMeasurePreservingError` — an error named after a map that was never built. §2.1 |
| "the outflow discard is a dead capability to be fixed" (implied, mine) | **[X]** It is an **asserted contract** — 4 tests + the 2-D balance gate demand those rows be zero. Still worth changing, but as a contract change. §7C.4 |
| "the `mixed` law is realizer-handled" (`sn/mesh/augmented_mesh.py:182-185`) | **False** — deleted in Wave 11; only deletion notices remain. §6 |

### 9.0 Two standing constraints every phase below must obey

**C-1 (from §7C.4).** The outflow discard is an **asserted contract**: 4 tests plus
the 2-D balance gate demand `B`'s outflow rows be zero. Any step that changes it
**re-poses those gates deliberately** — it never deletes them and never lands as a
silent fix.

**C-2 (from §7C.4).** `is_adjointable` is a *declared* predicate defaulting
`False`. **Adding `apply_transpose` alone changes nothing** — measured. Every
adjoint step is two-part or it is a no-op.

### 9.1 Net epistemic result

Of the concerns that opened this review, the measurements **confirmed** the
architectural ones (empty designed types, dead rank-N path, unreachable laws,
domain defect, string dispatch) and **refuted** the testing one. The subsystem's
**tests are strong** (30/31 mutations caught); its **structure** is what has
decayed away from a design that is written down and sitting in the repo.

That is the opposite of the usual failure mode, and it sets the plan's shape:
this is a **structural restoration against a known specification**, executed on a
subsystem whose gates will actually catch mistakes — with the caveat that four of
those gates encode a contract the restoration must deliberately re-pose (§7C.4).

---

# PART II — THE ACTION PLAN

**Target state = the specification already in this repo.** Grand Report v3 §16A.2:
`BoundaryTraceLaw` is a frozen dataclass of three **typed fields**
(`geometry_map`, `response`, `source`) with
`as_operator(trace_space) -> AffineBoundaryOperator(linear=R @ G, source=q)`.
This is a **structural restoration against a written spec**, not a new design.

**Technique (carried from the operator-strategy campaign):** every currently-red
gate ships as `xfail(strict=True)` naming its phase; **the strict-xfail set IS the
todo list**, and strictness means no phase can silently fix something without
deleting its marker. Verify each with `--runxfail` (Mode 8, class 4).

---

## B0 — Clean before extending *(no behaviour change)*

`coding-standards.md`: run the cleanup pass BEFORE the capability lands, so the
new structure is a no-op extension through one body rather than a third arm.

1. **Retire the false promises** (§4, §6, §8.1): the stale `mixed` advertisement
   (`sn/mesh/augmented_mesh.py:182-185`); the five "re-exported below" aliases that
   are not importable; `_zero_bulk_source`'s "single source" claim; the retired
   transient-adapter citation; the `reflect_inflow_inplace` misnomer; the
   `B_a` ray-block self-contradiction; the *"Concrete subclasses populate these"*
   line; the mask's *"sensitivity adjoints require"* justification.
2. ✅ **DONE — `law.kind` is now one derivation.** `RegistryMixin` already stored
   `key` per class, so `kind` became a single base property returning
   `type(self).key`. Before: a **mutable dataclass field** on Vacuum, a property on
   two, **absent on four**. After: uniform across all seven; `PrescribedInflow`'s
   redundant override deleted; `VacuumInflow(kind="banana")` is now a `TypeError`.
   **[M]** Gate: pyright 0; boundary set 389 unchanged; **full tree 3274 passed,
   byte-identical to baseline** (3274/5/111/57) — behaviour-neutral at scale, not
   inferred.

   **Deliberately NOT fixed — `ReflectiveBoundary.kind` → `"partial"` at α≠1.**
   It looked like a wart, but **[R]** `BC.to_alpha` handles `BC("partial",
   albedo=…)`, so the law is mirroring the *declaration* vocabulary on purpose.
   Deriving it from the registry key would have silently added partially-reflecting
   faces to the sweep schedule's reflective set — a semantic change, not cleanup.
   **B2 deletes the string dispatch that reads the tag at all**, so the question
   dissolves rather than needing an answer. The reasoning is written into the
   docstring so it is not "fixed" later by someone with less context.
3. **Coverage repairs** (§7C.3): re-pose the phantom `catches("ERR-052")` onto its
   mechanism; narrow the three skip-swallowed sentinels; replace the 9 tautological
   `pytest.raises` legs with production entry points; add the `|Ω·n|` mutation the
   two "hand-computed" white tests miss; relabel the snapshot file off `l1`; fix
   the 4 conflicting V&V markers; stop `_stub_for_test` leaking into the production
   registry.
4. **Docs**: the pre-P7b 2-arg walker snippets (`foundations:2833-2854`,
   `methods/sn:240` — the latter would `TypeError`); the MRO self-contradiction;
   the invariant table's 4 defects incl. the documented-but-absent
   `assert_direction_norm_preserved`.

### B0 — RESULT

- ✅ **B0.2 false-promise retirement.** ~18 claims corrected across 6 modules and
  2 theory pages, **[M] proven documentation-only by an AST comparison against
  HEAD with docstrings stripped** — 5 of 6 files executable-identical, the sixth
  differing only by B0.1's `kind` property (I re-ran that check independently).
  Found 2 stale walker snippets beyond the brief, and **corrected an error in this
  campaign's own `ddc1ee10`**: that commit's note said the flag was removed "from
  all six concrete laws"; measured from `ddc1ee10^`, the `ClassVar` sat on the ABC
  with explicit declarations on **three**. Sphinx baseline is now **0** warnings
  (the standing `paramref` ERROR is gone).
- ✅ **B0.3 coverage repairs**, all seven, each with a measured red. **[M] The
  auditor's 31-mutation set now catches 31/31** (was 30/31) and no mutation's red
  count fell. Skips 4 → 1 (the three inert sentinels now actually execute); the
  survivor is a legitimate 421-group-HDF5 precondition.
  - The ERR-052 finding sharpened: **the catalogued mechanism is unreachable in
    its own fixture** — with the bug re-introduced, `|φ|max` is identical at 6, 24
    and 32 outers because `F·φ/k` is scale-neutral once `k` converges. **No
    tolerance tightening could ever have rescued the old assertion.** The marker
    moved to a mechanism gate (`P(φ) = ∫νΣ_f φ dV = 1`, `rtol=1e-12`, no margin).
- ✅ **Two production defects** the doc/test phases correctly *reported instead of
  touching* (both executable, outside their briefs) — fixed by me:
  `_apply_faces`'s error said `apply` on the `apply_transpose` path
  (`sn/operators/boundary.py`), and `_base.py`'s lone `law=type(self).__name__`
  spelled the error tag as a class name where all six siblings use the registry
  key. The second is a **payoff of B0.1**: `law=self.kind` cannot drift now.

**Gate [M]:** boundary set `1229 passed, 1 skipped, 2 xfailed`; `sn/operators/` +
`geometry/` `1117 passed`; the four `match=` sites keying on the preserved
`mesh-identity invariant` substring still pass; **`npx pyright orpheus/` = 1 error
(the ratchet floor, #288)**; Sphinx `-E -W` exit 0, 0 warnings; full tree wide gate
pending. **Closes:** part of #306.

### B0 — findings handed forward (production, not fixed here)

- **[M]** `IncomingOutgoingTraceClassificationError` (**ERR-040**) has **zero raise
  sites**; its trigger is a no-op ABC default nobody overrides. **Fifth** instance
  of the §4 dead-capability pattern, this one inside the error module.
- **[M]** `orpheus/sn/solver.py:1930-1933` claims `_solve_fixed_source_si` calls
  `_reflect_outflow_into_inflow`; AST says it does not — there is exactly **one**
  production call, in the reconstruction sweep. Outside the boundary package;
  needs its own verification pass.
- **[M]** The **ERR-052 catalogue entry** needs two edits: its stated mechanism does
  not occur in that fixture, and its "Test reference" names the retired test.

> ⏸ **COMPACTION POINT** — commit, then re-anchor from this file + `git log`.

## B1 — Mint the typed specs and populate all seven laws *(pure addition)*

> ⚠ **PREREQUISITE, surfaced by B0.2.** The `R`-notation collision is unified in
> `geometry/boundary/__init__.py`, but **four sites still use `R` for the
> composite**: `sn/boundary/angular.py:41` (`R_white = G_diff ⊗ α`),
> `geometry/boundary/_composition.py:14` (`R = Σ_α c_α G_α`),
> `docs/.../foundations/boundary_conditions.rst:4876`, and
> `docs/.../methods/sn/boundary_conditions.rst:170`. These are the **§15.2
> sum-of-tensor-products framing** — a genuinely different decomposition, not a
> typo — so they need a deliberate terminology ruling rather than a sweep.
> **They collide with the `response_kernel` type this phase mints and must be
> resolved BEFORE it lands**, or the codebase will carry two meanings of `R`
> with one of them newly type-bearing.

Build what §16A.2 specifies and `orpheus/` never had.

- `BoundaryGeometryMap` / `BoundaryResponseKernel` Protocols, each with
  `as_operator(trace_space)` — **the law owns the spec, the trace space builds the
  matrix.**
- Concrete `G`: `IdentityMap`, `SpecularMirror(axis)`,
  `HemisphericalAverage(axis, outward_sign)`, `SpatialWrap(partner_face)`, `NullMap`.
- Concrete `R`: `ScalarAlbedo(α)` (with the α∈[0,1] invariant — **#265**),
  `ZeroResponse`, `UnitResponse`.
- **`PeriodicBoundary` gains its missing partner-face field** — **[R]** it has no
  fields at all today, which is why its `G` is not even expressible (**#183**).
- **`ZeroFluxBoundary` gets the `BoundaryRelation` tier it belongs to** — **[R]**
  `φ_Γ = 0` is `A₋γ₋ + A₊γ₊ = 0`, not an affine trace law. This is **#177** from
  the code side, and it is what lets SN's refusal become a type distinction
  instead of a hand-written `isinstance` guard.

**Gate:** for every law, `law.response_kernel.scalar` **is bit-identical** to the
float the realizers already multiply by (`sn/boundary/realizer.py:288`,
`diffusion/boundary_realizer.py:214`). Nothing reads the new fields yet — pure
addition, zero behaviour change. **Closes:** #177, #265, #183.

## B2 — Delete the stringly-typed dispatch

> ⚠ **B2 IS NOT MECHANICAL REPOINTING. Measured 2026-07-31, before any code was
> written — read this before starting.**

### B2.0 ⭐ THE BLOCKER — the law is DISCARDED at realization

**[R]** `orpheus/sn/mesh/augmented_mesh.py:435`:

```python
return _BoundBoundaryOperator(realized, kind=law.key)
```

**[M]** So `sn.bc[face]` is **not a law**. It is a wrapper holding the *realized
operator* (`inner`) plus a **string**. Probed on a live mesh:

```
kind             -> 'reflective'
albedo           -> <<ABSENT>>
geometry_map     -> <<ABSENT>>
response_kernel  -> <<ABSENT>>
axis             -> <<ABSENT>>
inner            -> <TensorProductOperator …>          # the REALIZED op
```

**The five dispatch sites are not being lazy — a string is the only thing that
survives realization.** They cannot be repointed at `law.geometry_map` because
at those sites *there is no law*.

**So B2 gains a first step the plan did not have:** `_BoundBoundaryOperator` must
**carry the descriptor it was realized from**, with `kind` becoming a
read-through instead of a stored copy. Small in code, but it is a production
change to a shim on the SN path (**[G]** `vacuum.py` alone has 46 external
callers), not a docstring pass. Fold in the related defect the review already
logged: **[R]** this wrapper forwards eight members to `inner` but **not**
`domain`/`codomain`, unlike every other wrapper in the codebase.

### B2.1 — The trap is RESOLVED; repointing IS behaviour-preserving

The B1 test file initially warned the opposite. **[X]** That inference was wrong
and is corrected at `668989fd`. Two measurements:

- **[M]** The wrapper stores `law.key`, **not** `law.kind`. `ReflectiveBoundary.key`
  is `"reflective"` for **every** albedo — 0.7 included — so
  `bc[face] == "reflective"` **already matches partial reflectors**, and
  `permutes_ordinates` is `True` for both. **They agree.**
- **[R]** Independently, `_law_from_tag` (`transport/method.py:299-301`)
  hard-codes `albedo=1.0` for reflective, so a partial reflector is **not
  declarable through a `BC` tag on any method** at all.

The real wrinkle the same measurement exposed, and the only one that matters:
`law.kind` and the wrapper's `key` **diverge** for a partial reflector, and only
the latter reaches production — precisely because the wrapper drops the law.
B2.0 dissolves that too.

### B2.2 — Then repoint, then retire

Six sites (§3.1 + diffusion's own) ask the structure instead:
`geometry_map.permutes_ordinates`, `geometry_map.is_adjointable`,
`response_kernel.is_zero`.

**Gate (mutation, mandatory):** change a law's spec object → each dependent
site's behaviour must change. A site that does not move is not consuming the
spec. Plus a **bit-identity** leg: with the law reachable, every realized
operator must be unchanged.

**Retire in the same change** — `_RULED_CORNER_KINDS` (**[M]** identical to SN's
admitted set, so its raise is unreachable) and the `== "reflective"`
tag-equality asserts (**[M]** 11.1 % of bare asserts, the test-side shadow of
the same defect).

### B2 — RESULT ✅ (2026-07-31)

**B2.0** `d812f7c8` · **B2.2** — the six sites now ask the factors:

| site | was | is |
|---|---|---|
| `sweep_schedule._reflective_faces` | `bc[face] == "reflective"` | `law.geometry_map.permutes_ordinates` |
| `B_b.is_adjointable` | `kind in _RULED_CORNER_KINDS` | `_has_ruled_corner_action(law)` |
| `B_b._reflect_corner` | `kind == "reflective"` | `law.geometry_map.permutes_ordinates` |
| `solver.compute_keff` | `op.kind == "vacuum"` | `op.law.response_kernel.is_zero` |
| `DSALowOrderSystem` (admission **+ row selection**) | `_SUPPORTED_BC`, `kind == "vacuum"` | the two factors |
| `DiffusionBoundaryRealizer._partial_current_albedo` | a 5-arm `isinstance` ladder | `law.response_kernel.scalar` |

**Retired:** `_RULED_CORNER_KINDS`, `_SUPPORTED_BC`, 18 test-side tag-equality
asserts (**rewired** to `isinstance(bc[f].law, …)` per coding-standards, not
deleted; two kept alongside as the deliberate legacy-surface pin on a real
mesh).

**[M] Equivalence, 7 laws × 5 sites, old expression vs new:** the retired tag
expressions are reproduced verbatim in
`tests/geometry/test_boundary_factor_consumers.py` and compared law by law.
**Every site agrees except one, deliberately:** the solver's leakage list now
includes `PrescribedInflow` (R = 0 ⟹ it leaks its whole outflow) where
`== "vacuum"` missed it. Correct, unreachable (**[M]** SN's admission table is
`{reflective, vacuum}` — five of seven laws cannot reach an SN face at all),
and given its own stating test rather than hidden.

**[M] Full gate, COMPLETE:** `4383 passed / 57 xfailed / 8 skipped / **0
failed**`, and the chunk arithmetic closes **exactly** against the 4448 tests
`tests/{geometry,sn,diffusion,numerics,transport}` collects — every one
executed, including all **60** `slow`-marked verification tests and all **3**
slow curvilinear L1 tests. (The tier is why the tree gate is slow:
`test_l1_standoff_slab_cylinder.py` alone is 10 tests in **56 min**.) pyright
`orpheus/` = 1 (the #288 floor, unmoved); sphinx `-E -W` exit 0, 0 warnings.

**[M] Mutation gate: all five sites RED.** Plus three in-process falsifications
of the gate itself (a plugin on `PYTHONPATH`, never a `git checkout`) —
wrong-factor sweep, dropped prescribed-exclusion, geometry-for-response
diffusion: 1/1/3 reds.

**The 3-search retirement audit bit on leg 3, exactly as the rule predicts.**
Graph callers and a code+tests+docs grep both looked clean; **direct callers of
the changed signature** did not. `DSALowOrderSystem._build` — a private
classmethod — has **three test call sites** passing BC *tag tuples*, and the
DSA admission test built its mesh from a `SimpleNamespace(kind="white")` **tag
surrogate** with no `law` at all. Both now speak laws: the tag→law translation
sits at the *test's* surface (`_LAW_FOR_TAG`), not pushed back into production,
and the surrogate carries real `WhiteBoundary` / `AlbedoBoundary` /
`PrescribedInflow` instances — which makes that gate strictly stronger, since a
tag surrogate could only ever exercise a string comparison the guard no longer
performs. The prescribed leg is the one that would have caught the trap below.

#### ⭐ A FOURTH SEARCH the audit needs — **duck-typed test surrogates**

Three sites in this phase, found one full-suite run at a time rather than by
audit, all the same shape: a test that stands up a **duck-typed stand-in** for
the changed handle and therefore appears in **no graph query and no symbol
grep**, because it never names the symbol.

| surrogate | where | what broke |
|---|---|---|
| `SimpleNamespace(kind="white")` | DSA admission test | `.law` absent |
| `bc={face: "vacuum"}` (bare strings) | 3-D sweep-schedule stub | `'str' has no attribute 'law'` |
| `class _NotALaw: pass` | diffusion refusal test | `.geometry_map` absent |

The graph cannot see them (no edge), a symbol grep cannot see them (no symbol),
and the type checker cannot see them (they are `SimpleNamespace` / `Any`). The
**only** thing that finds them is running the suite — so a signature change to a
widely-consumed handle owes a **full-suite run before the commit**, not after,
and the audit's greppable form is *"what do tests build in place of this
object?"*: `grep -rn 'SimpleNamespace(' tests/` near the consumer, plus the
handle's attribute names.

The third one is not a test defect at all — it is a **production regression the
surrogate caught**. Collapsing diffusion's `isinstance` ladder deleted its final
fallthrough `raise`, so a non-law that used to get a named `BoundaryError`
listing the realizable cases would have died on a bare `AttributeError` deep in
a factor read. Restored as a `response_kernel is None` guard, first in the
method, since everything after it dereferences a factor. **A collapsed dispatch
ladder loses its last arm silently — check for the fallthrough explicitly.**

**Two traps caught while writing it, both the same shape** — a *value* test
where the honest test is a *family* test. `isinstance(law.source, NoSource)`
admits `PrescribedInflow()` at its default zero source, which would have turned
a refusal into a value at BOTH the DSA admission (silently building a Marshak
row that drops `q` the day one is set) and the diffusion realizer (realizing
𝒜 = 0 on the #290 P5 path that does not exist). The disqualifying property
belongs to the FAMILY whatever `q` currently holds, so both sites keep one
`isinstance(law, PrescribedInflow)` — an *essential* type test, not the tag
smell B2 removes.

**Findings handed forward (not fixed here — both pre-existing):**
- **[M]** `_reflect_corner`'s specular swap is **unscaled** — it ignores `R`, so
  a partial reflector would re-emit its full outflow at the μ = ±1 corner. The
  tag set admitted every albedo too (`ReflectiveBoundary.key` is `"reflective"`
  regardless), so B2 changed nothing; a `.. warning::` now states it, and **B4
  closes it** by composing `R ∘ G` there.
- **[M]** `compute_keff`'s leakage list excludes partially-reflecting faces,
  which leak `(1 − α)`. Unreachable (`_law_from_tag` hard-codes `albedo=1.0`)
  and harmless (the term it skips is `net_current`, identically zero for the
  perfect reflector it means to exclude), but the honest predicate is `R ≠ 1`
  and it becomes reachable the moment **#189** admits partial reflectors.

> ⏸ **COMPACTION POINT.**
>
> **RE-ANCHOR HERE (2026-07-31).** B0, B1, B2 **and B3.0/B3.1/B3.2 are landed**, plus a
> full doc-correctness repair. **NEXT = B3.4, THEN B3.3** (user ruling: the mask
> retirement is cheap and independent; the un-narrowed laws are what block the algebra).
>
> | phase | commit |
> |---|---|
> | B3.0 — the G/R correction + the response-kernel tier | `9e2139b4` |
> | B3.1 — mint `TraceRestrictionOperator` (`γ±`) | `b39502f8` |
> | B3.2 — narrow the SN law's domain to `Γ₊` | `7f02de15` |
> | doc repair — every claim B3.0–B3.2 falsified | `b11a2ce3` |
>
> **The B3.4 brief is §9.3 of `.claude/plans/b3_domain_narrowing_crosswalk.md`** — read it
> before touching code; it carries the six-rows/four-kinds scope, the three design rulings,
> and the measured traps. **Verify merge status against `git log`, never against this
> file** (`.claude/rules/process-discipline.md`).
>
> Baselines to beat: full tree `-m "not slow"` **4335 passed / 6 skipped / 73 xfailed /
> 0 failed**; `npx pyright orpheus/` **= 1** (the #288 floor); `sphinx -E -W` exit 0.

## B3 — Narrow the SN law's domain to `Γ₊` *(the root fix — obeys C-1)*

> **DESIGN OF RECORD: `.claude/plans/b3_domain_narrowing_crosswalk.md`** — the
> mandatory pre-carve convention crosswalk (Pattern 7 / L17), carrying the
> measured answers to every opening question and the **§7 G/R correction**. Read
> it before touching code; the sub-step table is §8 there and tasks #10–#16 here.

**[R]** The theory page types `G : Γ₊ → Γ₊`, `R : Γ₊ → Γ₋`. **[R]** Diffusion
already hands its law exactly `trace.outflow_view(face)`. SN hands the whole face
slot, and every symptom in §4 follows from that.

- SN's realized law takes the **outflow trace**, returns the **inflow trace**.
- The mask's dead outflow preservation, the `P_in ∘ P_out = 0` composition, and the
  two `P_in` spellings + one `I − P_in` (§5) collapse into **one named restriction**
  — **[G]** measured absent from all of `orpheus/numerics/`.

### What the measurements changed (2026-07-31)

- **§5 needs a fourth map.** **[M]** On a real CYL mesh (`product(2,4)`) **4 of 8
  ordinates at `xmax` are tangential**, so `I − P_in ≠ P_out`; `P_out` is a map no
  operator type spells, and the mask's *"projection onto the outflow subspace"*
  docstring is **false**. The slab is the unrepresentative fixture.
- **The narrowing cannot be a view.** The diffusion precedent does not transfer —
  `OUTFLOW_ROW` is a constant index (a real view); the angular set is a
  quadrature-dependent fancy index (a copy). **The restriction is an OPERATOR.**
- **⭐ USER RULING — the G/R split is corrected inside B3.** `G` is the composition
  operator of a measure-preserving phase-space bijection (decidable test:
  **multiplicativity**); the crossing is geometric, so `G : Γ₊ → Γ₋` and
  `R : Γ₋ → Γ₋`. The Lambertian average moves **out of `G` into `R`**, and the
  **response-kernel tier is minted now**. `G` is unobservable exactly when `R` is
  rank-one — the theorem that let white's misassignment survive.
- **⭐ USER RULING — albedo and periodic are BUILT, not refused.** Albedo's gap is
  an incomplete `R` (magnitude without angular distribution), not an empty `G`;
  periodic's `G` is the translation reading the **partner** face's `Γ₊`.
- **C-1:** the gates are **2 test functions / 5 items**, and the honest blast
  radius is **16 items across 4 files** — **[M]** all signature breaks, zero value
  breaks. They are **re-posed** (RG-1…RG-5), never deleted. Two teeth findings:
  the outflow leg must **widen to `complement(inflow)`** or it is blind to the
  tangential rows, and the narrowed operator owes an explicit **domain guard**
  (a reduced permutation silently accepts a full-face input).

**Gate:** the composite's output is **bit-identical** (`np.array_equal` — no
reduction tree changes, so nulp would hide the bug class) on `{reflective,
vacuum}` at mesh level, law level for `white`/`prescribed_inflow`; the re-posed
gates state the new domain; **the mutation harness is promoted into the repo**
(**[M]** it currently lives only in a job tmp dir, so the stated gate is
unrunnable next session).

## B4 — Compose `R ∘ G` *(closes the `∘` gap at both layers)*

The realizer reads the law's two spec objects and composes them —
`OperatorProduct`, **[G]** currently used by no boundary operator — collapsing the
`isinstance` ladder into one generic body. Give the **descriptor algebra its
missing product** too (**[R]** `LawSum`/`LawScaled` have no `∘`), so `R ∘ G` is
expressible at both layers rather than one.

**Gate:** realized operators **bit-identical** to today, per law per method.
**Precedent honoured:** diffusion's two-stage shape (§7B.1) — the stage that
produces the operator must not see the law.

## B5 — White as a genuine rank-one, and its adjoint *(obeys C-2)*

**[R]** `AngularAverageOperator.apply` is already contract-then-broadcast; type it
as the rank-one it is (`RankOneOperator` was minted 18 minutes later for the same
structure and white was walked past).

**C-2: two-part or it is a no-op** — add `apply_transpose` **and** flip
`is_adjointable`. The transpose of `u ⊗ v` is `v ⊗ u`, so the cosine-weighted
adjoint the docstring says the BC "analytically possesses" becomes expressible
instead of refused.

**Gate:** `B.H` works with a white face installed (**[M]** today it raises
`MissingAdjoint` — the docstring's *"the one channel by which the white-BC adjoint
becomes available"* is currently impossible). Metric correctness gated against the
cosine-weighted inner product, not the Euclidean one. **Unblocks #300** (Woodbury
on the rank-1 `B`).

> ⏸ **COMPACTION POINT.**

## B6 — Reachability, and the rank-N verdict

- **#189**: admit `white` / `albedo` / `periodic` into SN's registry. **[M]** three
  of seven laws are unreachable under *every* method today; `BC.white` is a public
  tab-completable instance nobody admits.
- **Rank-N: ✅ USER RULING 2026-07-30 — WIRE IT, do not retire.** (§7B.4 measured
  it at zero production consumers with an in-principle-unreachable flagship.) The
  algebra is complete, closed and well-tested; the defect is that nothing reaches
  it. So B6 **admits `white` first (#189), then routes production single-BC
  realization through `realize_recursively`** so the rank-1 case is the degenerate
  path of the rank-N one — one body, not two. `0.3·Reflective + 0.7·White` becomes
  declarable, and the Marshak partial-current BC becomes a real capability rather
  than a docstring.
  **Gate:** the rank-1 path through the walker is **bit-identical** to today's
  direct realization, per law per method (the walker's leaf call IS
  `realizer.realize`, so this must hold exactly); plus a value gate on the Marshak
  mix that is not self-generated.
- **#244**: the `CPMesh._resolve_bc ≈ MOCMesh._resolve_bc` admission twin.

## B7 — Close out

Theory pages re-derived from the restored structure; **#179 / #180 / #181**
(MoC / MC / CP realizers) become mechanical — each is a new witness of
`BoundaryRealizer[MethodSpaceT]` (generic since `1fd15f64`) supplying its own
`as_operator(trace_space)`. **#219** (`MethodSpace` ABC) lands here or is closed as
superseded.

---

## Issue map — every phase closes filed work

| phase | closes / advances |
|---|---|
| B0 | #306 (partial) |
| B1 | **#177**, **#265**, **#183** |
| B2 | — (deletes an anti-pattern) |
| B3 | the §4 root defect |
| B4 | the `∘` gap, both layers |
| B5 | **#300** unblocked |
| B6 | **#189**, **#244**, rank-N verdict |
| B7 | **#179**, **#180**, **#181**, **#219** |

**Nothing in this plan is speculative architecture.** Every step either restores a
spec written in `.claude/plans/neutron_transport_grand_report_v3.md` or closes an
issue already filed and labelled.
