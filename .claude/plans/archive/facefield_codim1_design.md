# Codim-1 `FaceField` design — the ψ½ pole seed and the `BulkField`/`FaceField` duality

**Design note — 2026-07-06. DESIGN ONLY; no production code was written this
session.** This is the durable hand-off for a future implementation session, so
the reasoning is inherited rather than re-derived.

**Sequencing (user ruling, 2026-07-06):** this lands as a **codim-1 substrate
carve BEFORE Phase 3 DSA** (`#2`) in the stencil-assembly/DSA roadmap. Rationale
(clean-before-extend, `coding-standards.md`): DSA's restriction/prolongation
operate on the composite and its trace under the shared `|Ω·n|w` metric, so a
clean codim-1 substrate + the unified face measure benefits the DSA boundary
restriction directly.

**Supersession:** this design **supersedes** the "2.5f guard-retirement" framing
and **subsumes** the banked "make `FullField.__post_init__` enforce carrying ⟺
block-present" follow-up — the invariant question is answered structurally by
mesh-derived presence (§5.3), not by a `__post_init__` check.

**Provenance:** emerged from the `#280` walk-unification close-out (branch
`refactor/sn-walk-unification` @ `3f0b8c7`, campaign complete). The Gauss–Lobatto
empirical study (§3.5) is by the numerics-investigator; the tensor/bordered/
FaceField structural verdicts (§4) are by explorer + cross-domain-attacker
(agent memos + this note).

---

## ⏸ Status & compaction checkpoint (2026-07-06)

**⚠ RESHAPED 2026-07-06 (implementation session) — the metric-unification
premise (§3.2, §5.2) is REFUTED; §4 dead-end table row 3 and §5 stand corrected.**
A rigorous derivation + a numerics↔structure dialectic (both seats → the same
verdict) established:

- **ψ½'s Hilbert metric is `G_sd = V_cell`** (radial cell volume, SPD), NOT the
  through-flux zero. The shipped ghost `G_sd = 0` was a **real latent adjoint bug**
  (it severs the seed→bulk coupling from `A.H = G⁻¹AᵀG` — a wrong adjoint for any
  nonzero seed, production defect 1.3e-2; green only because every gate fed a zero
  seed). **ERR-067**, vv Mode-12. The three pole-vanishing quantities were conflated:
  **M1** moment weight (0, lives in φ) / **M2** through-flux `(1−µ²)` (0, an OPERATOR
  coefficient) / **M3** state metric (`V_cell`, the inner product). The bug installed
  M2 where M3 belongs. A state's metric is set by its operator role, not its
  integration weight.
- **`face_streaming_normal` does NOT unify the two codim-1 metrics** — it is the
  *spatial-trace* partial-current measure `|Ω·n|·w` ONLY. The pole is metrically
  BULK-like (`V_cell`); it does not route through the kernel. The through-flux equals
  the state metric only when the face's operator self-block is trivial (trace `A_tt`
  ≈ restriction map ‖2‖; pole `A_ss` ≈ banded radial transport ‖71‖).
- **Converged type design = option (c):** a `FaceField` ABC unifies **STRUCTURE ONLY**
  (FaceLayout flat buffer + presence-invariant); the **metric descends per-leaf**
  (spatial trace → `|Ω·n|·w`; pole → `V_cell`), exactly as `BulkField`'s `V·w` is
  per-leaf. The user's "trace = reduced in SPACE / pole = reduced in ANGLE, both
  grazing faces, pole is a face IN the bulk" harmonization is this design.
- **Rename ruled:** `StartingDirection*` → **`RadialCharacteristicFlux`** (retire the
  forward-sweep-role name + "ghost metric"). Lands in commit 2+, NOT commit 1.

**Execution = "correctness fix first" (user ruling):**
- **Commit 1 (THIS session) = the `G_sd = 0 → V_cell` correctness fix + Mode-12
  closure.** DONE (metric install; reciprocity/Mode-12/TestA4 gates inverted;
  diag_gsd promoted → `tests/sn/operators/test_starting_direction_metric.py`, the
  ERR-067 catcher; `face_streaming_normal` re-scoped to trace; docstrings + theory
  rewritten, anchor `sn-282-ghost-metric-face` → `sn-282-pole-state-metric`). Forward
  path bit-identical. Full wall + pyright `transport:1` + `sphinx -W` green.
- **Commit 2+ = the structural carve:** `FaceField` ABC (structure-only) + re-parent
  + the typed-key `FaceLayout` (old A(ii)) + mesh-derived presence retiring the 7
  guards + the `RadialCharacteristicFlux` rename. Then P3 DSA.

Re-anchor after `/compact` from **this note + `git log`** (trust git). The FIRST
action of the structural carve is still a **proactive `test-architect` dispatch**
(numerics spaces ↔ transport fields ↔ sn mesh — a MUST proactive trigger).

### ⏸ COMPACTION CHECKPOINT (2026-07-06) — commit 1 DONE, commit 2+ next

**Commit 1 LANDED @ `eb5943d`** (branch `refactor/sn-walk-unification`, **34
ahead of main, pushes HELD**): the `G_sd = 0 → V_cell` correctness fix + Mode-12
closure. Full wall **3211 passed**, pyright `transport:1`, `sphinx -W`=0. The
by-policy-uncommitted durable knowledge persists on disk across `/compact`:
**ERR-067** (`error_catalog.md`), the **vv-principles Mode-12 second-closure-
mechanism** uplift (`SKILL.md` — "repair the metric"), and the sub-agent memory
updates. Nothing commit-1 lives only in conversation.

**C2 LANDED (2026-07-06) — the FaceField codim-1 carve, THREE bit-identical
commits on `refactor/sn-walk-unification`:**
- **C2a `c637407`** — typed-key `FaceLayout[K]` (`FaceSlot`/`FaceLayout` →
  `Generic[K]`); ψ½ `StartingDirectionSpace` migrated onto a real
  `FaceLayout[tuple[int,int,str]]`, the hand-rolled
  `_leg_offset`/`cells_slice`/`corner_slice`/`_per_sign` retired; the 6 str
  consumers annotated `FaceLayout[str]`.
- **C2b `4081c0d`** — `FaceField(Field, Generic[K])` ABC owns the flat-buffer
  discipline ONCE; `BoundaryField(FaceField[str])` & `StartingDirectionField(
  FaceField[(level,sign,part)])` re-parented as **SIBLINGS** (the pole is NOT a
  `BoundaryField` — the `full_field.py:274` discriminator holds); `_trace_space_of`
  + `_space_of` unified to one `_face_space_of` hook; `FaceSlot.name → key`.
- **C2c `be5e7f8`** — single-walk ψ½ metric (offsets + `V_cell` weights emitted
  in ONE loop, no parallel re-derivation) + `K` bound to `Hashable` (elegance-
  review polish).

Elegance-enforcer PASS (0 must-fix). Deferred `LayoutBearingSpace` Protocol →
**#295** (`module:numerics`; retires the `FaceField.space` `# type: ignore` +
getattr-guards). Gates: transport+numerics+diffusion 1352, full tests/sn 1982,
the new `test_facefield_hierarchy.py` + `test_face_layout_typed_key.py`; pyright
`transport:1`; sphinx -W clean. **NEXT = C3** (mesh-derived presence via a
`PhaseSpaceCarrier` Protocol + retire the 7 `_require/_refuse_starting_direction`
guards) → C4 (`StartingDirection* → RadialCharacteristicFlux` rename) → C5 (docs
→ DSA).

**Commit 2+ phase plan (the structural carve — confirm ordering at plan time).**
Clean-before-extend, each step bit-identical on its own gates:

- **C2 — `FaceField(Field)` ABC + typed-key `FaceLayout`.** Introduce
  `FaceField(Field)` owning the `FaceLayout` flat-buffer discipline + the
  presence-invariant — **STRUCTURE ONLY, NO metric** (metric stays per-leaf:
  trace `|Ω·n|·w`, pole `V_cell`). Re-parent `AngularBoundaryField` /
  `ScalarBoundaryField` / the pole field under it. Generalize
  `FaceLayout`/`FaceSlot` to a typed key `[K]` so the pole space stops
  re-implementing `_leg_offset` / `cells_slice` / `corner_slice`
  (`(level,sign)` keys) and shares the ONE impl with the trace (`str` keys).
  Bit-identical offsets + metric. (This is the old A(ii) + Phase B, minus the
  refuted metric-unification.)
- **C3 — mesh-derived presence + retire the 7 guards.** A `PhaseSpaceCarrier`
  Protocol in `transport` (SNMesh satisfies it — sn→transport only) + mesh
  block-enumeration; the composite derives `starting_direction` presence from
  the mesh (never a free `Optional` arg); retire `_require_starting_direction`
  / `_refuse_starting_direction` (7 sites, `sn/loss_representation/__init__.py`
  + callers); migrate their negative tests to the construction-time invariant.
  Retirement audit = the 3-search (graph callers / text-grep code+tests+docs /
  direct constructors).
- **C4 — rename `StartingDirection*` → `RadialCharacteristicFlux`.** Broad
  mechanical rename (nexus `rename` + 3-search audit); retire the "starting
  direction" (forward-sweep-role) + residual "ghost metric" vocabulary. Its
  OWN commit (high blast radius, easy review).
- **C5 — docs + verification close-out.** Flip the theory PLANNED admonitions
  (`loss-rep-facefield-codim1`) to LANDED; full SN+numerics+transport wall;
  pyright ≤ 1; `sphinx -W`. **Then P3 DSA** (#2) — its restriction/prolongation
  ride the composite + trace metric, so the clean codim-1 substrate benefits it.

**Standing constraints (unchanged):** surgical mode (main agent writes; NO
`method-implementer`; `test-architect`/`explorer`/`qa`/`elegance-enforcer`/
`archivist` encouraged); canonical `.venv/bin/python -O -m pytest -p no:xdist
--timeout=300 -p no:cacheprovider` SERIAL; pyright oracle = `npx pyright` / the
ratchet CLI (`transport:1` baseline), streamed `<new-diagnostics>` = #226 LSP
artifact IGNORE; pushes HELD; `.claude/agent-memory/*` + `scratch/*` + the
`.claude/skills/*` accumulators NOT the main agent's to commit; never
`git checkout`/`restore` on uncommitted files (monkeypatch to revert mutations);
commits stamp **Opus 4.8** (session model); `sphinx -W` must stay clean.

**Standing constraints for the implementation session:**
- **Surgical mode** — the main agent writes the operator/type carve DIRECTLY, user
  steers step-by-step; NO `method-implementer`. `test-architect` (proactive),
  `explorer`, `qa`, `elegance-enforcer`, `archivist` allowed/encouraged.
- **Tests** — canonical `.venv/bin/python -O -m pytest -p no:xdist --timeout=300
  -p no:cacheprovider` SERIAL (xdist unstable for tests/sn + tests/numerics). Host
  env: `.venv/bin/python`.
- **pyright** — oracle is `npx pyright` / the ratchet CLI (baseline `transport:1` =
  accepted #288); streamed `<new-diagnostics>` = documented #226 LSP artifact,
  IGNORE.
- **Pushes HELD** (never `git push`); `scratch/` + `.claude/agent-memory/*` are NOT
  the main agent's to commit; never `git checkout`/`restore` on files with
  uncommitted work (mutation-revert via monkeypatch).
- **Model attribution** — commits stamp the model doing the commit, from session
  context (go-forward; this session = Opus 4.8).
- **Sphinx** must build `-W` clean; the theory pages carry the design direction
  (marked PLANNED) — flip each to current-truth as its phase lands.

**Untracked (on disk, safe across compaction; user's to keep or bin):** the
declined Gauss–Lobatto study — `derivations/diagnostics/diag_glob_0{1..5}_*.py`
(33 green diagnostics) + `scratch/experimental/glob_sphere_study/` (driver). NOT
committed (promote-only-if-adopted; the pole-node direction was DECLINED, §3.5/§4).

---

## 1. The problem

`FullField` is the **generic cross-method composite** — diffusion builds
`FullField(bulk=ScalarFlux, boundary=ScalarBoundaryFlux)` (2-block); SN
curvilinear builds 3-block. It carries an `Optional` third block,
`starting_direction: StartingDirectionField | None = None` — the SN-curvilinear
ψ½ pole seed. Presence-vs-mesh consistency is policed by **7 runtime guards**
(`sn/loss_representation/__init__.py`): `_require_starting_direction` (positive
half, "carrying ⟹ present") at 3 sites; `_refuse_starting_direction` (negative
half, "non-carrying ⟹ absent") at 4 sites.

The `Optional` + guards are the smell: **presence is a free constructor argument
that must be kept consistent with a fact only the mesh knows.** The mesh already
authoritatively computes the presence (`SNMesh.starting_direction_levels` /
`.starting_direction_space`), so there are two sources of truth for one fact, and
the guards are the machinery reconciling them (a Pattern-2 violation). Worse, the
SN-curvilinear-specific block is bolted onto the *generic* type, so the guards
also stop diffusion/Cartesian from ever seeing a block they have no concept of.

---

## 2. What ψ½ actually is (the physics, dug this session)

- **A structured per-level array, not a scalar.** `StartingDirectionSpace` (a flat
  backing buffer) holds, per carrying μ-**level** and per direction **sign**
  (μ = ∓1, both legs — R13): `cells (ng, nx)` — the angular flux ψ at the μ=±1
  pole edge at *every radial cell* — plus a `(ng,)` corner at r=R (inflow =
  given data / outflow = defect). Units `1/(cm²·s·sr)` (it IS an angular flux,
  shared with `AngularFlux`). Layout deliberately mirrors `AngularBoundaryFlux`'s
  `FaceLayout` flat-backing.
- **A genuinely solved DOF.** `FullField.to_flat()` packs
  `[bulk.values.ravel(), boundary.values, starting_direction.values]`, so ψ½ is
  in the solution vector; assembled, it is extra matrix rows *and* columns. The
  ERR-053 Krylov fix (size `n_dof = initial_guess.to_flat().size`) exists exactly
  because ψ½ grows the augmented vector.
- **Metric-null and triangular.** Its metric is the *ghost metric* — all weights
  exactly 0 (see §3.2) — so it contributes 0 to every inner product / norm /
  k-functional (it lies in the kernel of G; this is the vv Mode-12 blindness of
  the G3 reciprocity gate to seed-row errors). And the augmented matrix is
  **triangular** in ψ½ (the #284 certificate): the seed couples to the bulk in
  walk order, solved by forward substitution, not iteration.
- **Presence per level by τ_raw** (§3.4): sphere-GL carries one level;
  cylinder/slab/Cartesian carry none.

---

## 3. Physics rulings (durable — these are FACTS; they belong in the theory pages)

### 3.1 The pole is a straight characteristic

μ=±1 are the radial directions. Radial streaming has no angular redistribution —
the redistribution coefficient `(1−μ²)` vanishes at the poles (equivalently the
α-dome endpoints `α_{1/2}=α_{N+1/2}=0`). So ψ(μ=±1) satisfies a **pure 1-D
spatial transport ODE** (μ ∂ψ/∂r + σ_t ψ = Q, with μ=∓1), decoupled from the
angular recurrence — Hébert §3.9.4, "the μ=±1 rays of a sphere are straight
characteristics." This is *why* ψ½ is directly solvable (route (a), `#282`:
`carlson_inward_sweep_from_source`, the Hébert 3.434–3.435 diamond march) and why
the augmented system is triangular. The pole's special solve is **physics, not
representation** — it survives any storage choice.

### 3.2 The "ghost metric" is the entirely-grazing angular face

There is one face functional — call it `face_streaming_normal(face)` = *the normal
component of the streaming flux through the face, which vanishes when the
characteristic is tangent to the face* — and both codim-1 metrics are instances:

- **spatial face** (the boundary trace): `|Ω·n|·w` — vanishes at **grazing** Ω
  (direction tangent to the r-face);
- **angular face** (ψ½): `(1−μ²)·w` — vanishes at μ=±1 because **the pole rays
  are characteristics** (radial streaming does not cross angular cells → the flux
  is *entirely tangent* to the angular face).

So the "ghost metric" is not a special zero-measure carrier — it is *the μ=±1
angular face being entirely grazing*, the exact analog of a spatial face at
grazing incidence. Naming this one measure collapses the two current
constructions (`AngularTraceSpace` builds `|Ω·n|`; `StartingDirectionSpace`
hard-codes zeros) into one function — a falsifiable substrate win (§7 Phase A):
one `face_streaming_normal` must reproduce `G_trace == |Ω·n|·w` *and* `G_sd == 0`
at 0 ULP.

**Caveat — two different measures.** The ghost is the *through-flux /
redistribution* measure (used for the operator's angular-derivative term and the
block's inner product). It is distinct from the *scalar-flux / moment* measure
`∫ψ dμ`, for which the poles are regular endpoints of a 1-D integral. Under GL the
pole is a seed with zero moment weight (GL's interior nodes integrate the moment
to 2n−1 without it); under Lobatto the pole would carry a nonzero moment weight
(see §3.5). Do not conflate the two.

### 3.3 Circle-vs-interval — the real asymmetry

It is **neither** "sphere vs cylinder" **nor** "curvilinear vs cartesian." Two
distinct questions:

1. *Does the geometry have an angular-redistribution term* `(1−μ²)∂/∂μ`? →
   curvilinear-vs-cartesian (sphere **and cylinder** yes; cartesian no).
2. *Does the angular sweep need a separate, off-node starting DOF?* → the τ_raw
   predicate (§3.4), which is **quadrature-structural, not geometric**.

ψ½ answers question 2. The deciding fact is the **topology of the redistribution
axis**:

- **Cylinder** — redistribution is in the **azimuthal angle** ω, a **circle**
  (periodic `[0,2π)`). The production rule (`rules_product.py`) is Gauss–Legendre
  in the axial cosine × **equispaced** in φ (`linspace(0, 2π, n_φ, endpoint=False)`).
  Equispaced sampling of a periodic domain is the trapezoidal rule — *spectrally
  accurate* — **and**, for **even n_φ**, it hits φ=π (μ_x = −sinθ, μ_y = 0), the
  most-inward radial direction. So the most-inward azimuthal ordinate lands
  *exactly on the level's starting edge* → `τ_raw,0 = 0` → the seed is a bulk
  ordinate for free, at no accuracy cost. **This is #229, and it is structural,
  not a clamp** (the `[½,1]` Morel–Montry clamp is a *separate* recurrence-weight
  stabilizer; the R12a predicate reads the *un*clamped τ_raw). It is contingent on
  **even n_φ** — odd n_φ misses φ=π and the cylinder *would* carry a seed.
- **Sphere** — redistribution is in the **polar cosine** μ ∈ [−1,1], an
  **interval** whose endpoints are physical poles. The optimal rule is
  Gauss–Legendre — but GL is *open* (no endpoints), so it cannot place a node at
  μ=±1. Hence a *separate* off-node seed DOF, `τ_raw,0 ∈ (0,1)`.

So a **periodic** redistribution axis gives edge-inclusion for free; an
**interval** axis makes you pay for it. The cylinder is the existence proof that
"seed = bulk edge-ordinate" works — it works there because the axis is a circle.

### 3.4 The τ_raw presence trichotomy (R12a)

Presence per level iff the level's first-ordinate **unclamped**
`τ_raw,0 ∈ (0,1)` exclusive (`morel_montry_tau_raw_per_level`, the single source;
the clamped producer maps 0 → ½ and destroys the distinction). Bit-exact on
production rules:

- `τ_raw,0 = 0` — starting edge coincides with the first node (cylinder *product*;
  hypothetically sphere-*Lobatto*). Seed is a rank-duplicate of ψ₀ → **no block**.
- `τ_raw,0 = 1` — first node on the second edge (cylinder *level-symmetric*:
  duplicate-η nodes collapse the midpoint edge). Recurrence thread weight
  `(1−τ₀) = 0` → dead → **no block**.
- `τ_raw,0 ∈ (0,1)` — genuine independent state (sphere-GL ≈ 0.39–0.42) →
  **carries the block**.

The predicate is **per-level** and **only the mesh can compute it** (it needs the
quadrature). A geometry-name dispatch cannot express it (it would mis-handle a
Lobatto-sphere or an odd-n_φ cylinder). R12a corrects the earlier R12
"μ_start ∉ nodes" spelling, which is empirically false on level-symmetric
cylinder rules.

### 3.5 Gauss–Lobatto affordability (the empirical study)

The question: could the sphere use a quadrature that places a node at μ=±1
(Gauss–Lobatto, exact to 2n−3 vs GL's 2n−1), making the pole a bulk ordinate and
`τ_raw,0 = 0`? The numerics-investigator study (scratch, uncommitted; see §8):

- **Affordable.** At resolved order (N ≥ 8, N > L for Pₗ), GLob tracks GL at a
  bounded **~1.2× error penalty** (~1.3–1.4× at S8 → ~1.2× at S16), i.e. **~1–2
  extra ordinates** to match GL. **Not amplified by anisotropy** (P0→P5 all
  ~1.2–1.3×) and **insensitive to regime/c**. |Δk| ~30–140 pcm at S16; fine-N GL
  and GLob agree < 6 pcm. Only the under-resolved N ≲ L corner breaks (rank-
  deficient for GL too).
- **Not a drop-in.** GLob's μ=−1 node lands on the lower edge ⇒ `τ₀ = 0` ⇒ the
  production recurrence `ψ_{3/2} = (…)/τ₀` divides by zero, and the presence
  predicate keys on `τ_raw ∈ (0,1)`. Adopting a pole-node quadrature *requires*
  the seed→bulk-ordinate restructure. Handled the cylinder way (seed = the node,
  straight-char solved, recurrence starts from it).
- **Pole weighting is unbiased** (GLob converges to the same N→∞ limit; hits
  φ = Q/Σ_t exactly). The straight-char pole handling is confirmed
  redistribution-consistent (per-ordinate flat-flux residual ~1e-15).

**Ruling: affordable but architecturally DECLINED** (see §4/§5). The value of the
study is that declining the fold-in is a *principled architectural choice (keep
the bulk clean)*, not a numerical constraint.

---

## 4. Architecture exploration — dead-ends recorded (do not re-derive)

| Hypothesis | Verdict | Why |
|---|---|---|
| ψ½ is a factor/slice of an existing **tensor product** boundary | **NO** | The boundary/trace *spaces* are flat `FaceLayout` `⊕` buffers, not tensor products; the only boundary TPs are per-face *operator* laws `K_ω ⊗ I_g` (2-factor). ψ½ is off-quadrature (μ=±1 ∉ ordinates) and transverse to the trace (full-radial × angular-edge vs surface × full-ordinate). The design-intent 3-way `G_patch ⊗ K_ω ⊗ K_g` has a *spatial-patch* third factor, not ψ½. |
| A **`BorderedOperator`** (border/un-border) unifies RQI + composite + condensation | **NO** | Three distinct maths: RQI = a genuine **saddle-point/KKT** system (zero corner, λ as multiplier); the composite = a **dagger-biproduct** block algebra (blocks are traces `ι*`/`ι_*`); condensation/homogenization = **lossy Petrov–Galerkin projection**, *not* Schur. The "border ⊣ un-border" pair the intuition reached for **already exists** as the biproduct inject/project dagger pair (`BlockRole` + `FullFieldSpace` + `_AdjointOperator`). |
| The cell/face **`FaceField` duality** (mimetic/staggered/Hodge) is the codim-1 dual of `BulkField` | **YES** | `AngularBoundaryField` (spatial faces) and `StartingDirectionField` (angular edges) are both codim-1 face objects — siblings by accident, no shared parent (which is why ψ½ *hand-copies* `AngularBoundaryFlux`'s layout). Not cohomology (`ι_*∘ι* = id`, a dagger adjoint, not `d²=0`). The face measure genuinely unifies (§3.2). |
| **Dissolve the block** by folding the pole into the bulk (Lobatto or zero-weight augmentation) | **DECLINED** | Numerically affordable (§3.5) but the bulk is *cell-centered*; a pole ordinate makes it a **mixed field** — an inert (zero-weight), straight-char-solved, redistribution-special passenger that every bulk consumer (homogenization, condensation, moment extraction, `for ordinate in bulk`) must know about and skip. Cardinal Rule 2 smell (two concepts in one type, forcing a demux downstream). The zero weight means no *numerical* corruption, but the *conceptual* pollution is the cost. Reject. |

---

## 5. The design (the peak)

### 5.1 The rulings

1. **The pole ψ½ is a codim-1 `FaceField`, kept SEPARATE from the cell-centered
   bulk.** Embedding is rejected (§4 last row / the clean-bulk argument, §5.4).
2. **Presence is mesh-derived, not a free `Optional` argument.** The mesh
   enumerates its phase-space blocks (bulk always; spatial trace always; angular
   trace iff `τ_raw ∈ (0,1)` on that level); the composite's block-list mirrors
   the mesh. Mismatch becomes unconstructable → the 7 guards retire, and **no DOF
   is added where unneeded** (absent = not present, not empty-slot). Keyed on the
   §3.4 structural predicate, not geometry name.
3. **Metric and layout unified with the spatial-trace `FaceField` sibling** — one
   `face_streaming_normal` measure (§3.2); a typed-key `FaceLayout` so
   `StartingDirectionSpace` stops re-implementing the flat-buffer discipline.
4. **Bounded, sphere-only, triangular-solved** — one pole level, both signs,
   spatially resolved; it does not proliferate with N.

### 5.2 The type hierarchy

Current (`transport/fields/_bases.py`): three siblings off `Field` —
`BulkField {Angular, Scalar, Moment}`, `BoundaryField {AngularBoundary,
ScalarBoundary}`, and `StartingDirectionField` (lonely, straight off `Field`).

Target — introduce the missing codim-1 parent:

```
Field
├── BulkField   — codim-0 (cell centers);  metric = cell/volume measure
│                  Angular / Scalar / Moment
└── FaceField   — codim-1 (faces/edges);   metric = face_streaming_normal
     ├── AngularBoundaryField / ScalarBoundaryField   (spatial faces)  |Ω·n|·w
     └── StartingDirectionField (ψ½)                  (angular edges)  (1−μ²)·w ≡ 0
```

`FaceField` owns the `FaceLayout` discipline and the presence-invariant once, so
the three face families stop being siblings-by-accident.

### 5.3 The mesh-enumerated block-list (and the layering)

The composite is a direct sum whose summands are exactly the blocks the mesh
reports — *the mesh's phase-space DOF structure, reified*. Presence becomes a
*derived reading of the mesh*, not a chosen argument.

- **Layering is safe.** `FullField` (transport) consults the mesh **duck-typed**
  (the current `getattr(mesh, "starting_direction_space", None)` already does
  this); the block *types* live in `transport`/`numerics`, never `sn`. Formalize
  the duck-typing as a **`PhaseSpaceCarrier` Protocol in transport** that `SNMesh`
  (sn) satisfies — sn→transport (correct), never transport→sn.
- **Stays within "mesh provides sizes, not field init"**
  ([[feedback-mesh-owns-machinery-not-storage-init]]). Enumerating blocks returns
  *spaces* (size + layout + metric), which the field inits itself *from*
  (`role.zeros_on(space)`). The rule's test is "does the mesh initialize
  storage?" — it does not. **Watch-seam:** if the enumeration must also *name the
  role* per block, that is still "structure, not values," but verify at build.

### 5.4 Why the clean-bulk argument decides it

`BulkField` (cell centers) and `FaceField` (faces) are **genuinely different
types**. The pole is a face object; the bulk is cell-centered. Folding the pole in
(§4) crams two concepts into one type and forces every downstream bulk consumer to
demux the pole back out — the exact Cardinal-Rule-2 failure the `FaceField`
primitive exists to prevent. Homogenization/condensation must see an *untouched*
`BulkField`. So the pole stays on the face side, in its own clean codim-1
representation. **The Lobatto detour was decisive-by-elimination**: it proved the
fold-in is numerically affordable, so keeping the block separate is a chosen
architecture, not a forced one.

---

## 6. DOF accounting

Bounded and sphere-only. Per carrying level (sphere = 1), both signs (μ=∓1, needed
for the adjoint + reflective corner) plus r=R corners: **~`2·ng·(nx+1)`** DOF
(e.g. ~164 for 2g × 40 cells, vs ~640 for the S8 bulk). It is **one additional
angular level (the pole)** — it does *not* grow with N — spatially resolved (nx
per sign, because the M-M recurrence needs the seed per cell), and **solved
triangularly** (forward substitution from the source, not iterated → cheap DOF).
Not a single scalar; a radial profile bounded to one pole level.

---

## 7. Implementation phases (next session)

- **Phase A — substrate wins (independent of presence; buildable first).**
  (i) Name the one `face_streaming_normal(face)` measure; route both
  `AngularTraceSpace` and `StartingDirectionSpace` metrics through it; gate: one
  function reproduces `|Ω·n|·w` (spatial) *and* `0` (pole) at 0 ULP.
  (ii) Generalize `FaceLayout` to a typed key so `StartingDirectionSpace`
  (`(level,sign)` keys) and the trace (`FaceLabel` keys) share one flat-buffer
  impl.
- **Phase B — `FaceField` ABC + re-parent.** Introduce `FaceField(Field)`;
  re-parent `AngularBoundaryField`/`ScalarBoundaryField`/`StartingDirectionField`;
  move the `FaceLayout` + presence-invariant machinery onto it. Bit-identical.
- **Phase C — mesh-enumerated blocks + mesh-derived presence + retire guards.**
  Add the `PhaseSpaceCarrier` Protocol (transport) + the mesh block-enumeration;
  make the composite derive `starting_direction` presence from the mesh (never a
  free arg); retire the 7 `_require`/`_refuse` guards; migrate their negative
  tests to the construction-time invariant. Legitimate *dispatch* (do the seed
  step iff carrying) is ideally mesh-derived operator assembly (which arms are in
  the sum), not a runtime `if`.
- **Phase D — docs + verification.** Theory pages (§ below); the retirement audit
  (3-search: graph callers, text-grep code/tests/docs, direct constructors);
  full SN + numerics + transport wall; pyright ratchet ≤ 1.

**Verification posture:** presence-derivation is *structural* (assert the composite
block-list equals the mesh's report; a mismatch is unconstructable). Bit-identity
where the storage is merely relocated; principled re-baseline only if a metric
value genuinely changes (it should not — the face measure is the same values, one
construction).

---

## 8. Artifacts & pointers

- **Lobatto study (scratch, uncommitted — user's to keep or bin):**
  `scratch/experimental/glob_sphere_study/{driver.py, run_sweep.py}` +
  `derivations/diagnostics/diag_glob_0{1..5}_*.py` (33 green diagnostics: 01
  moment integration, 02 per-ordinate consistency, 03 end-to-end penalty, 04 pole
  τ₀=0, 05 driver faithfulness + k_inf anchor). Promotion targets (if a pole-node
  scheme is ever adopted) are recorded in each file's docstring; do NOT promote
  unless adopted.
- **Agent memos:** explorer (tensor-structure recon; cross-method map),
  cross-domain-attacker (bordered/FaceField framing), numerics-investigator
  (`glob_vs_gl_spherical_quadrature_study.md`, lesson L16). Transient — this note
  is the durable synthesis.
- **Code anchors** (re-derive by grep; structure is durable): the 7 guards
  `sn/loss_representation/__init__.py:{194,213}` (+ call sites); ψ½ storage
  `numerics/spaces/starting_direction_space.py`, `transport/fields/
  starting_direction_flux.py`; the τ_raw source `sn/spatial/
  pole_angular_closure.py:morel_montry_tau_raw_per_level`; the product quadrature
  `numerics/quadrature/rules_product.py`; `FullField.to_flat`
  `transport/full_field.py:503`; the type bases `transport/fields/_bases.py`.

---

## 9. Deferred primitives with triggers (NOT this carve)

- **`SaddlePointOperator` / bordered systems** — RQI is one genuine saddle-point
  instance (inf-sup, zero corner) but premature (unwired, `#277`). **Trigger:** the
  diffusion **mixed-form** (`#294`) is the second genuine saddle instance (flux
  bordered by the div constraint, current as multiplier); build it then, and
  retire RQI's inline `(n+1)×(n+1)` border onto it. Adopt the saddle-point/KKT
  vocabulary (Benzi–Golub–Liesen), not a generic "block matrix."
- **Structured *spatial* `FaceField`** (ergonomic cell→faces `(i±½,j)`
  addressing) — a *different* object from ψ½ (which is a µ-level-keyed angular
  trace with no spatial-face incidence). **Trigger:** `#294` mixed-form interior
  currents / CP 2-D interface currents / MoC track-face crossings — the genuine
  consumers of a persistent structured interior FaceField. `FaceLayout` today is a
  flat named-slot buffer; no structured cell→face map exists.

---

## 10. Open decisions for next session

- Phase A/B/C boundaries are gate-able as bit-identical; confirm the sequencing
  (A independent; B before C) at plan time.
- The `PhaseSpaceCarrier` Protocol surface (does it enumerate spaces only, or
  spaces + roles? — §5.3 watch-seam).
- Whether `ScalarBoundaryField`/CP traces want to move onto `FaceField` in the
  same carve (cross-method reach) or as a fast-follow.
- Confirm this carve's scope vs. the `#267` (SNMesh→MaterialMesh) and `#219`
  (MethodSpace) layout campaigns to avoid double-churn.

---

## 11. One-paragraph summary (for the fresh session)

ψ½ (the SN curvilinear starting-direction pole seed) is an irreducible codim-1
*face* DOF — the angular flux at the μ=±1 pole, which is a straight characteristic
(no redistribution → its "ghost metric" is just the entirely-grazing angular
face, the exact analog of spatial grazing). Whether a geometry needs a *separate*
seed DOF is a *quadrature-structural* fact (τ_raw ∈ (0,1)), decided by the
topology of the redistribution axis: the cylinder's azimuthal *circle* gets the
seed as a free bulk ordinate (equispaced is spectral + edge-inclusive), the
sphere's polar *interval* does not (GL is open). Making the sphere's pole a bulk
ordinate (Lobatto) is *numerically affordable* (~1.2×) but *architecturally
declined* — it would pollute the clean cell-centered bulk with a face passenger.
So the design keeps ψ½ as a separate codim-1 `FaceField` (the missing dual of
`BulkField`), present exactly where the mesh's τ_raw predicate says (retiring the
7 guards), with its metric and layout unified with the spatial trace. The
"bordering" and "structured spatial face" primitives are deferred with triggers
(diffusion mixed-form / CP-2D), because ψ½ is neither's vehicle.
