# B3.4c blast radius — PERIODIC boundary + the face-name convention

**Branch** `refactor/operator-strategy-layers` · main checkout
`/Users/rodrigo/git/nuclear/ORPHEUS` (NOT a worktree).
**HEAD at audit start** `cea73f8d`; Nexus graph built at `943b37c1` (dirty) —
code-current, only the plan-doc commit is newer.
**Line numbers are current at final read; the tree was MOVING during this audit
(see §0).** Re-derive with Nexus/grep before acting on any row.

> ### 🔴 READ §7 FIRST — THE CARVE LANDED WHILE THIS AUDIT WAS RUNNING
>
> Eight production files went from clean → modified during the audit. §7
> reconciles every finding below against the FINAL tree state: which risks the
> carve already closed, which are still open, and which findings are now
> archaeology. **§1–§6 are the audit as taken; §7 is the verdict.**

> **Scope note — worktree exclusion.** `.claude/worktrees/nexus-workspace-wiring/`
> is a *separate checkout* frozen at `f71df08` (2026-06-12). It contains its own
> copies of every file below (`orpheus/sn/boundary_realizer.py`,
> `docs/theory/boundary_conditions.rst`, …) which a naive repo-wide grep returns.
> **Every one of those hits is EXCLUDED from this audit** — they are not on this
> branch and must not be edited. All rows below are main-checkout paths.

---

## 0. ⚠ PREMISE CORRECTIONS — verify these before planning (Operating Principle 5)

Three premises in the brief have already moved. Two are in-flight *right now*.

| # | Brief said | Tree says | Evidence |
|---|---|---|---|
| **P1** | "The carve **mints** a face-name primitive in `orpheus/numerics/face_layout.py`" | **Already minted, uncommitted, and complete on disk.** `face_name()`, `face_normal()`, `face_opposite()`, `FACE_NAMES` + `_FACE_SUFFIX_SIGN`/`_FACE_SIGN_SUFFIX` are a +156-line working-tree edit to `orpheus/numerics/face_layout.py`. `git diff --stat` shows it; `git status` at my session start did **not** — it landed mid-audit. | `git diff orpheus/numerics/face_layout.py` |
| **P2** | `SNBoundaryOperator._reflect_trace` at `boundary.py:413-437` | Body is at **`orpheus/sn/operators/boundary.py:326-467`**; the per-face loop is **`:413-466`**. The `faces=`/`rows=` filters are `:401-412`. Directionally the brief is right (every law is fed THIS face's `γ₊` at `:419-422`). | Read |
| **P3** | `PeriodicBoundary` is a fieldless descriptor | **It carries `axis: str = "x"` since B1** (`orpheus/geometry/boundary/periodic.py:60`) and exposes `geometry_map → SpatialWrap(axis)` (`:63-76`) + `response_kernel → ScalarResponse(1.0)` (`:78-81`). Several `scratch/` notes and doc rows still describe the fieldless version. | `periodic.py:54-81` |

**L-007 intra-session drift is ACTIVE.** `orpheus/numerics/face_layout.py` was
clean at my first `git status` and `M` (+156) minutes later; an import probe run
between the two writes reported `face_name` / `face_normal` / `face_opposite` in
`__all__` but **missing from the module** — a transient mid-write state, since
resolved (re-probe: all 8 `__all__` names resolve; `face_opposite("xmin") ==
"xmax"`, `face_normal("zmax") == (2, 1)`). Re-run any census in this file
immediately before acting on it.

---

## 1. `PeriodicBoundary` — the complete consumption surface

### 1.1 Nexus vs grep — where they disagree

| Query | Result |
|---|---|
| `mcp__nexus__callers(py:class:orpheus.geometry.boundary.periodic.PeriodicBoundary)` | see §1.5 |
| `grep -rn "PeriodicBoundary"` (main checkout, code+tests+docs) | **43 live lines / 18 files** |

Nexus is structurally blind here for the same reason as lesson **L-009**: the
only *call* into `PeriodicBoundary` is a zero-argument constructor in tests plus
one `isinstance` arm — and `isinstance` is not a `calls` edge. **Trust grep as
primary for this symbol**; Nexus contributes the doc-node edges (§5) and the
`SpatialWrap`/`realize` call chain (§2).

### 1.2 Production consumers (code that ships)

| # | Site | Kind | What it does |
|---|---|---|---|
| P-1 | `orpheus/geometry/boundary/periodic.py:21-81` | **definition** | The law. `axis: str = "x"`; `geometry_map → SpatialWrap(axis)`; `response_kernel → ScalarResponse(1.0)`. |
| P-2 | `orpheus/geometry/boundary/__init__.py:672` | re-export | `from .periodic import PeriodicBoundary` |
| P-3 | `orpheus/geometry/boundary/__init__.py:742` | `__all__` | public API |
| P-4 | `orpheus/sn/boundary/realizer.py:144` | import | |
| P-5 | **`orpheus/sn/boundary/realizer.py:745-754`** | **the SN realization arm** | `if isinstance(law, PeriodicBoundary): return stamp_boundary_role(PeriodicWrapOperator() & IdentityOperator())` — **THE carve site.** |
| P-6 | `orpheus/sn/boundary/realizer.py:790` | error text | names periodic in the "supported laws" refusal message |
| P-7 | **`orpheus/diffusion/boundary_realizer.py:295-316`** | **the diffusion REFUSAL** | keyed on `isinstance(law.geometry_map, SpatialWrap)` — **not** on the class (see §1.4) |
| P-8 | `orpheus/numerics/operator.py:2621`, `:2638` | docstring xref | `PeriodicWrapOperator` docstring — **`:2638` xrefs `PeriodicBoundary.apply`, which does not exist** (dangling; see §5) |
| P-9 | `orpheus/geometry/boundary/_factors.py:98`, `:495` | docstring | the factor table row + the "`PeriodicBoundary` has never carried ANY field (#183)" note — **P3-stale**, it carries `axis` now |
| P-10 | `orpheus/geometry/boundary/_base.py:168` | docstring | lists periodic among laws overriding no assertion hook |
| P-11 | `orpheus/geometry/boundary/__init__.py:235-239`, `:421`, `:510`, `:522` | package docstring | the law catalogue |

**Nobody in `orpheus/` constructs a `PeriodicBoundary(...)`.** The only
production reference is the `isinstance` arm P-5 and the structural test P-7.
Confirmed: `grep -rnE "\bPeriodicBoundary\(" orpheus/` → **0 hits**.

### 1.3 Is it reachable from the BC TAG registry? — **NO, for both methods**

Measured, not inferred:

```
SNMesh.BOUNDARY_OPERATOR_REGISTRY        -> ['reflective', 'vacuum']
DiffusionMesh.BOUNDARY_OPERATOR_REGISTRY -> ['albedo', 'reflective', 'vacuum', 'zero_flux']
```

`BC("periodic")` therefore raises at `orpheus/transport/method.py:303-312`
(`law_cls is None` → `ValueError: … does not support boundary condition
'periodic' on face '…'. Supported: 'reflective', 'vacuum'.`). Admitting it is
**issue #189**.

> ### ⚠ LATENT ORIENTATION BUG — the highest-value finding in §1
>
> `_law_from_tag` (`orpheus/transport/method.py:266-347`) has explicit
> construction arms for `ReflectiveBoundary` (`:313`), `WhiteBoundary`
> (`:317`) and `AlbedoBoundary` (`:339`), then a **parameter-free
> fall-through `return law_cls()` at `:347`**.
>
> Its own docstring warning (`:294-301`) states the rule:
> *"**A law that declares its own orientation MUST be listed here.** The
> `law_cls()` fall-through hands such a law its dataclass DEFAULTS, which are
> an orientation — and one that is right for exactly one face."*
>
> **`PeriodicBoundary` gained exactly such an orientation (`axis="x"`) at B1
> and was never added to the parse.** It sits in the fall-through today. The
> day #189 admits `"periodic"`, every face — `ymin`, `ymax`, `zmax` — receives
> `SpatialWrap(axis="x")`. This is *the identical latent shape* B3.4a fixed for
> white, one law later. It is invisible now only because the registry refuses
> the tag.
>
> **Recommendation:** B3.4c should add the `PeriodicBoundary` arm
> (`axis=AXIS_NAMES[label.axis_index]`) in the same commit that makes periodic
> real — the carve is what converts this from latent to live.

### 1.4 The diffusion refusal — quoted verbatim

Brief cited `orpheus/diffusion/boundary_realizer.py:68`. That line is the
**docstring table row**:

```
``PeriodicBoundary``             —        REFUSED (needs the opposite-face wrap)
```

The **executable** refusal is at `orpheus/diffusion/boundary_realizer.py:295-316`
and it is keyed on the FACTOR, not the class:

```python
if isinstance(law.geometry_map, SpatialWrap):
    # ── REFUSAL AXIS: spatial / topological ──────────────────────
    raise BoundaryError(
        f"DiffusionBoundaryRealizer cannot realize "
        f"{type(law).__name__}: its geometry map is a SPATIAL "
        f"wrap, which couples a face to its OPPOSITE face (a "
        f"trace-block permutation) — the one geometry P1 cannot "
        f"integrate away into a per-face albedo J⁻ = 𝒜·J⁺. The "
        f"wrap lands with the boundary-operator assembly when a "
        f"diffusion consumer exists (#290 P4 seam).",
        law=type(law).key or type(law).__name__,
    )
```

The `:296-306` comment carries the structural reason in full — *"What blocks it
is that the wrap couples a face to its OPPOSITE face, and this realizer's
codomain is a per-face scalar 𝒜 — a **block-diagonal object with no slot for
cross-face coupling**."*

**This is the free oracle for the whole carve** (explorer lesson L-011): the
diffusion sibling REFUSES the identical law with the exact structural reason SN
is about to have to honour. SN's silent identity realization (P-5) is the
suspect arm; diffusion already wrote the correct diagnosis. Note that the
refusal predicate is `isinstance(law.geometry_map, SpatialWrap)` — so it fires
for **any** future law carrying a spatial wrap, not just `PeriodicBoundary`, and
needs no edit when periodic's realization changes.

### 1.5 Doc + derivation + plan consumers

Docs → §5. Derivations: **zero** hits in `orpheus/derivations/`. Plans/scratch
(`.claude/plans/`, `scratch/`) are archaeology, not consumers — but note
`scratch/b3_4_periodic_reality.md` is a prior investigation of exactly this
question and `scratch/review_law_layer.md` / `scratch/review_diffusion_arm.md` /
`scratch/boundary_chain_map.md` carry overlapping (and now partly stale) maps.

---

## 2. `SpatialWrap` and `PeriodicWrapOperator`

### 2.1 `SpatialWrap` — 15 graph edges, 13 grep lines. Tiny and clean.

| Site | Kind | Role |
|---|---|---|
| `orpheus/geometry/boundary/_factors.py:471-519` | **definition** | `axis: str = "x"`; `permutes_ordinates → False` (`:507-510`); `is_adjointable → False` (`:512-519`) |
| `orpheus/geometry/boundary/_factors.py:287` | `__all__` | |
| `orpheus/geometry/boundary/_factors.py:99`, `:340`, `:448` | docstring | factor table; the Protocol's "declared rather than derived only because `SpatialWrap` answers False"; the orbifold-vs-torus contrast |
| `orpheus/geometry/boundary/__init__.py:236`, `:624`, `:708` | re-export + `__all__` + docstring | |
| `orpheus/geometry/boundary/periodic.py:14`, `:64`, `:76` | the sole producer | `geometry_map → SpatialWrap(axis=self.axis)` |
| **`orpheus/diffusion/boundary_realizer.py:143`, `:295`** | **the sole non-definitional CONSUMER in production** | the refusal predicate |
| `tests/geometry/test_boundary_factors.py:43`, `:83` | test | the factor-pair table row |
| `docs/theory/methods/sn/boundary_conditions.rst:251`; `docs/theory/foundations/boundary_conditions.rst:442`, `:668`, `:2605`, `:2630`, `:3101` | doc | |

**Two declared-`False` predicates are the carve's flip-switches, and BOTH are
already documented as implementation gaps that B3.4 closes:**

| Predicate | Value | Comment says |
|---|---|---|
| `SpatialWrap.permutes_ordinates` (`_factors.py:507-510`) | `False` | *"Spatial, not angular: ordinate n at face A feeds ordinate n at face B."* — **stays False.** This is a statement about the ANGULAR index and remains true after B3.4c. |
| `SpatialWrap.is_adjointable` (`_factors.py:512-519`) | `False` | *"The translation IS adjointable as a map… `False` reports an IMPLEMENTATION gap… **B3.4 builds it and this flips, WITH its gate.**"* |

`permutes_ordinates == False` is what keeps periodic OUT of `_reflective_faces`
— see §3.4, where that turns out to be the single fact that saves the G-S split.

### 2.2 `PeriodicWrapOperator`

| Site | Kind | Role |
|---|---|---|
| `orpheus/numerics/operator.py:2615-2654` | **definition** | `apply`/`apply_transpose` both `np.asarray(x).copy()`; `is_adjointable → True` |
| `orpheus/numerics/operator.py:180` | `__all__` | |
| **`orpheus/sn/boundary/realizer.py:754`** | **the sole production caller** (Nexus-confirmed: `SNBoundaryRealizer.realize` is the only incoming `calls` edge from `orpheus/`) | `stamp_boundary_role(PeriodicWrapOperator() & IdentityOperator())` |
| `orpheus/sn/boundary/realizer.py:86`, `:162`, `:747-751` | import + docstring | `:747-751` is the honest in-code admission of the gap |
| `orpheus/geometry/boundary/__init__.py:239`, `:553`; `orpheus/geometry/boundary/periodic.py:38`, `:72` | docstring | |
| `tests/numerics/test_periodic_wrap_operator.py` (5 tests) | test | §6 |
| `tests/numerics/test_operator_capability_predicates.py:31`, `:67`, `:192` | test | `("PeriodicWrap", PeriodicWrapOperator(), False, True, STRUCTURAL_ABSENT)` |
| `tests/sn/operators/test_sn_boundary_realizer.py:44`, `:697-717` | test | `TestRealizePeriodic` |
| `tests/geometry/test_bc_equivalence_snapshot.py:710` | test | "*Before the re-anchoring this case pinned `PeriodicWrapOperator`'s …*" |
| `docs/theory/methods/sn/boundary_conditions.rst:490`; `docs/theory/foundations/boundary_conditions.rst:3099` | doc | |

**⚠ Three FALSITIES already live in `PeriodicWrapOperator`'s docstring** (this is
lesson L-011's delegation shape, still unrepaired at `orpheus/numerics/operator.py`):

| Line | Claim | Status |
|---|---|---|
| `:2622-2624` | *"the SN sweep handles the spatial wrap via **its own face-pair indexing**"* | **FALSE.** `grep -rn "face.pair\|partner_face" orpheus/sn/` → the mechanism does not exist. Grep the invented noun, not the symbol. |
| `:2638` | xrefs `orpheus.geometry.boundary.periodic.PeriodicBoundary.apply` | **DANGLING.** `PeriodicBoundary` has no `apply` (pure descriptor since #186/B3+β2). Nexus carries it as an `external` node of degree 1 — the graph's own signature for an unresolved xref. **Not caught by `-W`** (Python-domain roles are silent without `-n`). |
| `:2639` | *"the legacy … aliasing-safety contract (`psi_out.copy()`)"* | Refers to a method retired in Wave O step O.4a.1. |

---

## 3. The block-diagonality claim — who states it, who DEPENDS on it

### 3.1 The claim sites (`B`-over-faces only; per-ℓ / per-ordinate / metric block-diagonality is a different, unaffected claim)

| # | Site | Text | After B3.4c |
|---|---|---|---|
| BD-1 | `orpheus/sn/operators/boundary.py:43-45` | *"`B` is the block-diagonal composition over the mesh's true boundary faces … **it never mixes faces**"* | **FALSIFIED** |
| BD-2 | `orpheus/sn/operators/boundary.py:198` | class docstring: *"block-diagonal over the mesh's true boundary faces"* | **FALSIFIED** |
| BD-3 | `orpheus/sn/operators/boundary.py:220-223` | RULING P1: *"reflection is a per-system endomorphism ⇒ block-diagonal"* | **survives** — this is about the A/B **system** biproduct (trace vs ψ½ ray), not about faces. Do not repair by reflex. |
| BD-4 | `orpheus/sn/operators/boundary.py:318` | *"the block-diagonal G-adjoint metric `B.H` reads"* | **survives** — `G` is block-diagonal (bulk⊕trace), a property of the METRIC, not of `B`. |
| BD-5 | `orpheus/sn/operators/boundary.py:390-392` | *"`B` is block-diagonal over faces, so the subset action is the **EXACT restriction** (no cross-face coupling is dropped)"* | **reasoning falsified, conclusion SURVIVES** — see §3.2 |
| BD-6 | `orpheus/sn/operators/boundary.py:545-546` | `reflect_into_inflow` docstring: *"`B` is block-diagonal over faces → the subset is the exact restriction."* | same as BD-5 |
| BD-7 | `orpheus/sn/operators/boundary.py:619-620` | `split`: *"The same face set the block-diagonal law iterates"* | **cosmetic** — it just means "the face inventory" |
| BD-8 | `orpheus/sn/solver.py:1988-1989` | *"`B` is block-diagonal over faces ⟹ exact restriction."* | same as BD-5 |
| BD-9 | `orpheus/diffusion/boundary_realizer.py:304` | *"a per-face scalar 𝒜 — a block-diagonal object with **no slot for cross-face coupling**"* | **survives and is REINFORCED** — it is the refusal's whole justification |
| BD-10 | `docs/theory/foundations/boundary_conditions.rst:703` | restates BD-9 | survives |
| BD-11 | `tests/sn/operators/test_sn_boundary_operator.py:16` | module docstring: *"**block-diagonal over faces** — a single-face perturbation stays on that face"* | **FALSIFIED (and the test bites — §6)** |
| BD-12 | `tests/sn/operators/test_sn_boundary_operator.py:334`, `:352` | *"(`B` is block-diagonal over faces — it never mixes faces)"* | **FALSIFIED** |
| BD-13 | `tests/sn/operators/test_sn_boundary_operator.py:662`, `:678` | *"`B` is block-diagonal over faces, so the subset MUST be the EXACT restriction"* / *"face-block-diagonal, so xmax's output is independent"* | **FALSIFIED as stated**; the assertion may still pass (§6) |
| BD-14 | `tests/sn/operators/test_bc_extraction_matvec.py:605` | lists "block-diagonal" among the contracts a sibling file owns | pointer only |
| BD-15 | `orpheus/diffusion/operators.py:550` | *"block-diagonal over the mesh's boundary faces, `B.apply(ψ)`"* | **survives** (diffusion refuses periodic) — but becomes a claim about **diffusion's** `B`, not a universal one. Worth a scoping word. |

### 3.2 ⚠ The load-bearing distinction: the `faces=` restriction is a **ROW** restriction, and rows stay exact

This is the single most important structural finding in §3, and it is *not* what
the docstrings say.

Read `_reflect_trace` (`orpheus/sn/operators/boundary.py:401-437`): `faces` and
`rows` filter **`face_laws`** — i.e. which OUTPUT face slots get written. The
INPUT is always the whole `boundary: AngularBoundaryFlux`; nothing restricts what
a law may READ. So:

> `B.reflect_into_inflow(ψ, faces=F)` computes **the exact rows of `B` indexed by
> `F`**, whatever `B`'s off-diagonal structure is — because the full input trace
> is in scope at every face iteration.

Block-diagonality is therefore **sufficient but not necessary** for the
restriction to be exact. The docstrings gave the *sufficient* condition as if it
were the *reason*. After B3.4c the honest statement is:

> *`B` is block-**structured** over faces: `faces=` selects a set of block ROWS,
> and each block row is computed from the whole input trace, so the subset action
> is the exact restriction of `B` to those rows. (Pre-B3.4c every block row was
> supported on its own face — block-DIAGONAL — which made the same conclusion
> trivial.)*

**Verdict per call site:**

| Call site | Contract | Breaks? |
|---|---|---|
| `SNBoundaryOperator.reflect_into_inflow(faces=…)` (`boundary.py:522-548`) | "exact restriction of `B` to those faces' inflow rows" | **NO** — still exact. Docstring reason must change. |
| `SNBoundaryOperator.reflect_inflow_inplace(faces=…)` (`boundary.py:550-590`) | same + writes back | **NO, but see §3.3** — aliasing must be re-checked |
| `SNMaskedBoundaryOperator.reflect_rows_inplace` (`boundary.py:958-998`) | additive row completion | **NO for exactness; YES for freshness — see §3.4** |
| `orpheus/sn/solver.py:2002` (`_reflect_outflow_into_inflow`) | production driver; `faces=None` on the reconstruction path | **NO** |
| `orpheus/sn/loss_representation/__init__.py:4904-4910` | the inter-group reflect inside `_sweep_scheduled` | **NO** (periodic never appears in `group.reflect_faces` — §3.4) |
| `orpheus/sn/operators/scheduled_invertible.py:259-260` | binds `reflect=self.lower.reflect_rows_inplace` | **NO** (same reason) |

### 3.3 Aliasing — currently safe, and the safety is INCIDENTAL

`reflect_inflow_inplace` (`boundary.py:581-590`) reads and writes the SAME buffer:

```python
reflected = self.reflect_into_inflow(boundary_flux, faces=faces)   # full read → fresh buffer
...
for face in selected:
    boundary_flux.face_view(face)[inflow] = reflected.face_view(face)[inflow]
```

`reflect_rows_inplace` (`boundary.py:992-998`) has the identical read-then-write
shape with `+=`. **Both are safe today** because `_reflect_trace` completes every
read into a fresh `AngularBoundarySourceSink` before returning — a whole-trace
read barrier.

Two independent reasons it stays safe after B3.4c, worth pinning rather than
assuming:
1. the read barrier above (structural — survives any block structure);
2. the write touches only **inflow** rows while the periodic read consumes
   **outflow** rows, and `Γ₊ ∩ Γ₋ = ∅` within a face slot.

⚠ Reason 2 is *not* independent of the three-way partition ruling: `Γ₊ ⊔ Γ₋ ⊔
Γ₀` — the tangential set is a third bucket (explorer L-010; the campaign's own
"the face partition is THREE-way" ruling). Disjointness of inflow and outflow
holds; exhaustiveness does not. Nothing here needs exhaustiveness, but **a future
optimisation that drops the fresh-buffer allocation would silently depend on it.**
Recommend the read barrier be stated as a contract in `_reflect_trace`'s
docstring, not left implicit.

### 3.4 ✅ The G-S split is SAFE — and it is one predicate that saves it

`SweepSchedule.gauss_seidel` (`sweep_schedule.py:113-169`) only ever assigns
`reflect_faces` from `_reflective_faces(sn_mesh)` (`:122`), and

```python
# orpheus/sn/loss_representation/sweep_schedule.py:278-282
return frozenset(
    face
    for face in sn_mesh.angular_trace.layout.faces
    if law_permutes_ordinates(sn_mesh.bc[face].law)
)
```

`SpatialWrap.permutes_ordinates` is **`False`** (`_factors.py:507-510`), so a
periodic face is **never** in `reflect_faces`, never in `lower_inflow_rows`
(`sweep_schedule.py:195-215` iterates only `group.reflect_faces`), and therefore
**every periodic row lands in `B_upper`** — lagged by the SI driver as an
external gain. `lower_inflow_rows`'s docstring already says exactly this
(`:189-192`: *"a face never reflected in-sweep (vacuum, white, albedo,
**periodic** …) has ALL its rows in `B_upper`"*).

**Consequences:**

- `B_lower` remains free of cross-face coupling ⟹ the forward-substitution
  `M = (L+C−B_lower)` stays triangular-solvable. **The reified-inverse contract
  is untouched.**
- `reflect_rows_inplace` never sees a periodic face ⟹ the freshness question
  (would a periodic row read a not-yet-swept partner's outflow?) **does not
  arise**. The lagged `B_upper` route is exactly the "cyclic back-edge" the
  schedule is designed to tolerate.
- The `split` exactness argument at `boundary.py:602-604` — *"the specular map
  has no octant-diagonal, and the two row sets are complementary within each
  face's inflow by construction"* — is about the **row partition within a face**,
  which cross-face coupling does not disturb. `upper_rows` is
  `setdiff1d(inflow_indices_for_face(face), lower_rows[face])` (`:614-622`),
  a per-face set complement; still exact.

> ⚠ **The tripwire.** If B3.4c (or a follow-up) flips
> `SpatialWrap.permutes_ordinates` to `True` — an easy mistake, since periodic
> *is* an ordinate bijection (the identity) and the surrounding prose calls it a
> "pairing" — periodic enters `_reflective_faces`, gains `B_lower` rows, and the
> forward substitution acquires a cross-face read. That is a genuine
> **triangularity violation**: `M`'s inverse is computed by a sweep whose walk
> order knows nothing about a partner-face dependency. The `False` must be
> defended by a test with a comment saying WHY, not left as an unremarked
> default. (See `scratch/triangularity_question.md:205` — the same question was
> raised for `ReflectiveBoundary`/`PeriodicBoundary` earlier in the campaign.)

### 3.5 ⚠ `apply_transpose` — the one place the carve genuinely BREAKS

`_reflect_trace`'s transpose branch (`boundary.py:438-466`):

```python
out_boundary.face_view(face)[...] = gamma_out.apply_transpose(
    law.apply_transpose(restricted)
)
```

It reads face `f`'s slot and writes face `f`'s slot — hard-wired to a diagonal
block. If `B[f ← partner(f)]` is the forward block, its Euclidean transpose is
`Bᵀ[partner(f) ← f]`: **the transpose must write to the PARTNER's slot.** The
current loop cannot express that.

This is live, not hypothetical:

- `SNBoundaryOperator.is_adjointable` (`boundary.py:291-310`) reads each REALIZED
  law's own predicate. The realized periodic operator is
  `PeriodicWrapOperator() & IdentityOperator()` — a `TensorProductOperator` of
  two adjointable factors ⟹ **`True` today**, and the class docstring at `:294`
  asserts *"reflective / vacuum / periodic are"* adjointable.
- So a periodic mesh reaches `apply_transpose` **without any guard firing** and
  computes a wrong `Bᵀ`.
- Note the asymmetry the code already warns about at `:361-370`: off-diagonal
  permutation laws are *bit-identical under either spelling*, which is why the
  earlier version of this bug survived every reflective fixture and was caught
  only by het-VACUUM grid reciprocity. **The same blindness applies here.**

**Recommendation:** B3.4c must either (a) extend the transpose branch to scatter
into the partner slot, or (b) make `SpatialWrap.is_adjointable` stay `False` AND
have the realized operator report `is_adjointable = False`, so
`SNBoundaryOperator.is_adjointable` refuses the composite. Option (b) contradicts
`_factors.py:517-518` (*"B3.4 builds it and this flips, WITH its gate"*) — so
(a) is the intended path and needs a **grid-reciprocity gate on a PERIODIC
fixture**, not a reflective one.

---

## 4. The face-name convention `"{axis}{min|max}"` — full surface

### 4.1 The five brief-listed sites: ALL FIVE ARE ALREADY REPOINTED (uncommitted)

Every site the brief named has already been collapsed onto the new primitive in
the working tree. Confirmed by `git diff`.

| # | Site (brief) | Current line | Before | After |
|---|---|---|---|---|
| T-1 | `orpheus/numerics/spaces/angular_trace_space.py:181-185` (`_FACE_NORMALS`) | **`:230`** | a module-level dict comprehension over `AXIS_NAMES × (("min",-1),("max",+1))`, with its own `KeyError`→`ValueError` wrapper | `axis, sign = face_normal(face)` — the table is DELETED |
| T-2 | `orpheus/sn/boundary/angular.py:290` | **`:289`** | `f"{axis}{'max' if outward_sign == +1 else 'min'}"` | `face_name(AXIS_NAMES.index(axis), outward_sign)` |
| T-3 | `orpheus/sn/boundary/realizer.py:450` | **`:451`** | verbatim twin of T-2 | `face_name(AXIS_NAMES.index(axis), outward_sign)` |
| T-4 | `orpheus/sn/loss_representation/sweep_schedule.py:252` | **`:252`** | `f"{AXIS_NAMES[a]}{'max' if s > 0 else 'min'}"` | `face_name(a, +1 if s > 0 else -1)` |
| T-5 | `orpheus/transport/method.py:337` (reverse parse) | **`:336-338`** | `axis=AXIS_NAMES[label.axis_index]` **and** `outward_sign=+1 if label.face_name.endswith("max") else -1` — *two independent reads of one datum* | `axis_index, outward_sign = face_normal(label.face_name)` — **one parse yields both** |

### 4.2 ⚠ SITES THE BRIEF MISSED — three, and one is a genuine twin

| # | Site | Code | Status | Assessment |
|---|---|---|---|---|
| **M-1** | `orpheus/sn/loss_representation/__init__.py:730-733` (`_inflow_faces`) | `f"{AXIS_NAMES[a]}min" if s >= 0 else f"{AXIS_NAMES[a]}max"` | **NOT repointed** (file untouched by the carve) | A 6th render. Should be `face_name(a, -1 if s >= 0 else +1)`. |
| **M-2** | `orpheus/sn/loss_representation/__init__.py:739-742` (`_outflow_faces`) | `f"{AXIS_NAMES[a]}max" if s >= 0 else f"{AXIS_NAMES[a]}min"` | **NOT repointed** | A 7th render — **and it is literally `face_opposite ∘ _inflow_faces`.** Its own docstring says *"the opposite of `_inflow_faces`, axis by axis"*. This is the primitive's most on-the-nose consumer and it is currently a hand-written twin. **The strongest single argument that `face_opposite` earns its keep.** |
| **M-3** | `orpheus/sn/solver.py:1462-1463` | `areas[0] if face == "xmin" else areas[-1]` **and** `axis_index = {"x": 0, "y": 1, "z": 2}[face[0]]` | **NOT repointed** | Two separate violations: (a) a bare `face == "xmin"` sign test (the `else` silently absorbs `ymin`, `zmin`, and any unknown face → `areas[-1]`); (b) a hand-written axis literal, the exact anti-pattern `orpheus/geometry/boundary/_specular.py:106-108` names — *"`AXIS_NAMES.index` rather than a local `{"x": 0, …}` literal … a second spelling of it is a twin waiting to disagree."* Both collapse to `axis_index, sign = face_normal(face)`. |

### 4.3 Face-name PARSERS (name → axis or sign) — the complete list

| Site | Direction | Status |
|---|---|---|
| `orpheus/numerics/face_layout.py:180-206` (`face_normal`) | name → `(axis, sign)` | **THE primitive** (new) |
| `orpheus/numerics/spaces/angular_trace_space.py:230` | name → `(axis, sign)` | repointed → `face_normal` |
| `orpheus/transport/method.py:336` | name → `(axis, sign)` | repointed → `face_normal` |
| `orpheus/sn/solver.py:1462-1463` | name → sign, name → axis | **NOT repointed** (M-3) |
| `orpheus/transport/mesh/axis.py:140-147` (`FaceLabel.face_name`) | `(axis_index, endpoint)` → name | **a DIFFERENT map** — see §4.4. Deliberately not the same bijection. |
| `orpheus/transport/mesh/axis.py:398,415` (`_OUTWARD_ENDPOINTS`) | **endpoint** → sign | **structurally separate and CORRECT** — it reads `label.endpoint ∈ {"max","outer"}`, never a face name. Leave it. |
| `orpheus/geometry/boundary/_specular.py:111` | axis NAME → axis index | `AXIS_NAMES.index(axis)` — already single-sourced |
| `tests/numerics/test_angular_trace_space.py:512`; `tests/geometry/_generate_bc_equivalence_snapshots.py:263`; `tests/sn/operators/test_sn_boundary_realizer.py:527` | name → `(axis, sign)` | **test-side hand-parses — KEEP as-is.** These are mirror-not-import oracles (`vv-principles` L11): importing `face_normal` here makes the test tautological against the production derivation it verifies. Do **not** "clean these up." |

### 4.4 ⚠ The `"outer"` third endpoint — where it maps, and what a naive bijection breaks

**The mapping site is exactly one line:**
`orpheus/transport/mesh/axis.py:85` —
`_ENDPOINT_SUFFIX = {"min": "min", "max": "max", "outer": "max"}`,
consumed only at `FaceLabel.face_name` (`:140-147`), which raises `ValueError`
on any other endpoint. Measured:

```
FaceLabel(0, 'min'  ).face_name == 'xmin'
FaceLabel(0, 'max'  ).face_name == 'xmax'
FaceLabel(0, 'outer').face_name == 'xmax'      # ← the collision
```

**Would a naive `min`/`max` bijection break?** Answer in three parts:

1. **The rendering direction is SAFE.** `_ENDPOINT_SUFFIX` is many-to-one
   (`max` and `outer` both → `max`), but `face_layout.face_name(axis, sign)` is
   a bijection over `sign ∈ {−1,+1}`, a *different domain*. The two layers
   compose correctly: `FaceLabel` collapses endpoint→suffix, then the face-name
   world is a clean bijection. **They must stay separate types** — one is
   axis-geometry, the other is mesh-topology vocabulary.
2. **The inverse to `endpoint` DOES NOT EXIST, and nothing must be built that
   assumes it does.** `face_normal("xmax") == (0, +1)`, but you cannot recover
   whether the endpoint was `"max"` or `"outer"`. Any future
   `face_name → FaceLabel` helper is unsound. `orpheus/transport/method.py:329-334`
   already documents the correct discipline: read the SIGN off `face_name`, never
   off `endpoint`.
3. **⚠ THE REAL HAZARD — `face_opposite` on a curvilinear mesh returns a
   NON-EXISTENT face.** A solid radial axis has `endpoints == ("outer",)` — ONE
   face. `sn_mesh.bc` and `angular_trace.layout.faces` carry only `"xmax"`
   (documented at `orpheus/sn/mesh/augmented_mesh.py:166-169`). But
   `face_opposite("xmax") == "xmin"` — a face that **does not exist** on a
   sphere or cylinder. A periodic partner lookup written as
   `mesh.bc[face_opposite(face)]` raises `KeyError`, and
   `boundary.face_view(face_opposite(face))` raises on an unknown key.
   `face_opposite` is a **pure geometric involution and correctly knows nothing
   about mesh inventories** — so the *guard belongs at the realizer/consumer*,
   not in the primitive. B3.4c must decide and state: **periodic is a Cartesian-
   only law**, refused on any axis whose partner face is absent from the trace
   layout, with a loud `BoundaryError` naming the geometry (the diffusion
   refusal at `boundary_realizer.py:307-316` is the model).

### 4.5 A third face vocabulary exists (derivations layer) — out of scope but worth knowing

`orpheus/derivations/discrete/sn/sweep_acyclicity.py:115` uses
`_FACES = ("left", "right")`, and `orpheus/derivations/.../dsa.py` uses
`("left","right","bar")`. These are **algebra-of-record vocabularies**,
deliberately independent of the production string world (mirror-not-import).
Do NOT unify them into `face_layout`.

### 4.6 Adjacent hard-coded face consumers that a periodic pair would surprise

Not face-name *transcriptions*, but slab-shaped assumptions that a two-face
coupling touches:

| Site | Assumption |
|---|---|
| `orpheus/sn/acceleration/dsa.py:232` | `laws = (sn_mesh.bc["xmin"].law, sn_mesh.bc["xmax"].law)` — hard 1-D pair |
| `orpheus/sn/acceleration/dsa.py:687-690` | builds a `{"xmin": …, "xmax": …}` dict |
| `orpheus/transport/radial_characteristic_field.py:310-311` | `boundary_trace.face_view("xmax")` — the curvilinear outer face by name |
| `orpheus/sn/operators/boundary.py:690`, `:724`, `:797` | `sn_mesh.bc["xmax"]` — `RadialCharacteristicBoundaryOperator` (System B) is `xmax`-only |

DSA under a periodic BC is not in scope for B3.4c, but the pair at `dsa.py:232`
is where a periodic law first meets an acceleration consumer. Worth an issue.

---

## 5. Sphinx corpus

Nexus `documents` edges into `PeriodicBoundary` name **three** pages; grep of
`docs/` finds only two of them. **Nexus wins here** — the third,
`docs/api/diffusion_1d`, reaches the class through the *autodoc'd module
docstring* of `orpheus/diffusion/boundary_realizer.py`, so the string never
appears in any `.rst` source. Any repair to that module docstring silently
re-renders an API page; grep alone would have missed it.

### 5.1 Pages and claims

| Page:line | Claim | After B3.4c |
|---|---|---|
| `docs/theory/foundations/boundary_conditions.rst:96` | `PeriodicBoundary` — "was …" (rename history) | ok |
| `…:378-382` | *"The Lambertian's `v = |Ω·n̂|` is invariant under both the mirror and the periodic translation, so `R∘G = R`"* → **"G is unobservable exactly when R is rank-one"** | **survives — and it is the theorem that explains why the identity realization went unnoticed.** Periodic's `R = ScalarResponse(1.0)` IS rank-one-ish (a scalar), so `G` is invisible to any rank-one probe. Cite it in the carve's commit message. |
| `…:400-419` | the quotient table: periodic = free action = torus = covering space; reflective = orbifold | survives |
| `…:441-442` | `PeriodicBoundary(axis)` → `SpatialWrap` | **already correct** (post-B1) |
| `…:668-671` | *"`SpatialWrap` answers … Periodic is the more …"* (equivariance) | survives |
| `…:701-704` | the diffusion refusal restated + "block-diagonal object with no slot for cross-face coupling" | **survives, and is the sentence B3.4c should quote when it re-words BD-1/BD-2** |
| **`…:940-953`** | **"What remains is four rows … plus `periodic`. … periodic's `G` reads the PARTNER face's `Γ₊`, which **B3.4c** builds (#183, #189). MEASURED: they silently accept a `Γ₊` input and echo it back — i.e. `Γ₊ → Γ₊`, the wrong codomain, invisible to a shape check (vv Mode 12)."** | **This paragraph is the carve's spec, already written.** It must be REWRITTEN past-tense on landing (the "what remains" list drops to albedo-only… and per §6 albedo already landed at B3.4b, so this list is ALREADY partly stale). |
| `…:964-968` | *"`SNMethodSpace.minimal` is now a partial constructor … only `albedo` and `periodic` still realize from it, and precisely because they have not yet been narrowed."* | **FALSIFIED ON LANDING** — after B3.4c, `minimal` must no longer suffice for periodic (it cannot name a partner face). Repair in the same commit. |
| `…:2629-2630`, `:2720-2721`, `:2751`, `:2782`, `:2836-2837`, `:2930`, `:3098-3101`, `:3711`, `:4544`, `:4787` | catalogue / rename-history / "not yet narrowed" rows | `:3098` (*"`PeriodicBoundary` — **not yet narrowed**"*) and `:3101` (`PeriodicWrapOperator() & IdentityOperator()`) are the two **must-repair** rows |
| `docs/theory/methods/sn/boundary_conditions.rst:247-251` | factor table: *"spatial wrap along `axis`; **the realizer derives the partner face from the installation face**"* | **This is a FUTURE-TENSE claim stated in the present.** The realizer does no such thing today (`realizer.py:745-754` ignores the face entirely). B3.4c makes it true — until then it is falsified prose. |
| `…:373-378` | rename history (`PeriodicBoundaryOperator` →) | ok |
| `…:488-491` | *"`"periodic"` — not yet narrowed (B3.4c)"* → `PeriodicWrapOperator() & IdentityOperator()` | **must-repair** |
| `…:526` | *"…and `periodic`, whose `G` must read the PARTNER face's `Γ₊`. Those rows still emit full-N endomorphisms, are unreachable from this registry, and are pinned by strict xfails until B3.4b / B3.4c land."* | **must-repair** (B3.4b has landed) |
| `…:544` | `AXIS_NAMES[label.axis_index]` — *"so the partner is correct at any …"* | **stale after T-5**: `method.py` now reads `face_normal(label.face_name)`. Repair with the face-name carve. |
| `docs/theory/foundations/boundary_conditions.rst:5318-5344`, `:5520-5523`, `:5979-5991` | the `AXIS_NAMES` narrative + a transcribed `FaceLabel.face_name` code block | **stale after the face-name collapse** — `:5331` shows `return f"{AXIS_NAMES[self.axis_index]}{suffix}"` (still true) but `:5523` shows the retired `_law_from_tag` white spelling. Repair with T-5. |
| `docs/theory/foundations/path_integral.rst:1115` | *"named, `periodic` cyclic from a single law"* | survives |
| `docs/theory/foundations/operator_algebra.rst:1216` | *"and periodic masks, the source dyads"* | check in context; low risk |
| `docs/theory/methods/sn/cartesian_multid.rst:3894` | *"row of a never-reflected face — vacuum, white, albedo, periodic"* | **survives** (§3.4 — periodic stays in `B_upper`) |
| `docs/theory/methods/sn/slab_multigroup.rst:601`; `docs/theory/methods/sn/curvilinear_numerics.rst:101` | *"reflective, white, albedo, periodic, and mixed BCs are honoured"* | ⚠ **already over-claiming** — periodic is not reachable from either registry. Worth correcting while you are in the file. |
| `docs/theory/foundations/operator_tensor_network.rst:174` | list mention | ok |
| `docs/theory/methods/monte_carlo.rst:592-617` | MC's OWN periodic BC (`:label: periodic-bc`) | **UNRELATED subsystem** — MC has an independent periodic implementation. Do not touch; do not let a grep-and-replace reach it. |

### 5.2 Dangling xref to repair in the same commit

`orpheus/numerics/operator.py:2638` →
`orpheus.geometry.boundary.periodic.PeriodicBoundary.apply` (method does not
exist). Nexus carries it as an `external` node, degree 1. **Sphinx `-W` does NOT
catch this** (Python-domain roles are silent without `-n`), so it will survive
the build gate. See §2.2 for the two neighbouring falsities in the same
docstring.

---

## 6. Tests — full inventory, classified

**17 test files** touch periodic (`PeriodicBoundary` / `PeriodicWrapOperator` /
`SpatialWrap` / the `"periodic"` tag). Classification per
`.claude/rules/coding-standards.md` (behavioural contract / API smoke /
characterization) with the RED/XPASS prediction.

### 6.1 ⚠ WILL GO RED / XPASS — act on these

| File:line | Test | Class | Prediction |
|---|---|---|---|
| **`tests/geometry/test_bc_equivalence_snapshot.py:773-830`** | `TestPeriodicLebedev17Snapshot::test_delivers_the_partner_faces_outflow` — `@pytest.mark.xfail(strict=True, reason=_B34C_XFAIL)` | **behavioural (THE gate)** | **XPASS(strict) → FAIL.** *This is the point* — the reason string says so explicitly: *"when it lands, this test's `_apply` line hands the operator the PARTNER's Γ₊ … and this marker is DELETED — the XPASS(strict) failure is the point."* Measured today: **98 % relative error**. |
| **`tests/sn/operators/test_b3_domain_narrowing.py:312`** | `_LAWS["periodic"] = (PeriodicBoundary(), True)` — the xfail flag | **behavioural (the todo-list row)** | Flip `True → False`. `_B34C_XFAIL` (`:259-268`) is DELETED. If left, `test_realized_law_has_gamma_minus_codomain` XPASSes strict. |
| **`tests/geometry/test_boundary.py:513-530`** | `test_periodic_bc_returns_input_unchanged` — asserts `psi_in == psi_out` bit-for-bit through `_realize_for_sn(bc, quad)` = `SNMethodSpace.minimal(quad)` | **behavioural — pins the WRONG contract** | **RED.** Note `_realize_for_sn`'s own docstring (`:56-65`) already says *"**B3.4b note: PERIODIC ONLY.**"* — periodic is the last consumer of this legacy helper. When periodic narrows, **`_realize_for_sn` becomes dead and should be retired with it.** |
| **`tests/sn/operators/test_sn_boundary_realizer.py:696-720`** | `TestRealizePeriodic::test_periodic_returns_tensor_product_and_passes_through` — asserts `op.ops[0] is PeriodicWrapOperator`, `op.ops[1] is IdentityOperator`, values pass through unchanged, built from `SNMethodSpace.minimal(quad)` | **behavioural — pins the placeholder** | **RED** on at least the pass-through assertion, and on `minimal(quad)` if the new realization demands a face-aware space (which §5 `foundations:964-968` says it must). **Rewire to the successor**, do not delete. |

### 6.2 Green-but-FALSIFIED-docstring — repair prose, assertions survive

| File:line | Test | Why it stays green |
|---|---|---|
| `tests/sn/operators/test_sn_boundary_operator.py:332-357` | `test_block_diagonal_no_face_mixing` | Every fixture in this file is vacuum/reflective (`_CASES`, `:118-128`; the perturbation test hard-codes `_sn("SLB", (BC.reflective, BC.reflective))` at `:335`). **No periodic mesh ever reaches it.** The assertion is still true *for the fixtures it runs*; the module docstring (`:16`) and inline claims (`:334`, `:352`) become globally false. |
| `tests/sn/operators/test_sn_boundary_operator.py:656-730` | `TestFaceRestrictedReflect` (3 tests) | Same: reflective fixtures only. And per §3.2 the *conclusion* (exact restriction) survives even with periodic — only the stated *reason* ("block-diagonal") is wrong. Repair `:662`, `:678`. |
| `tests/sn/operators/test_bc_extraction_matvec.py:605` | comment listing sibling contracts | pointer only |
| `tests/numerics/test_periodic_wrap_operator.py:1-14` | module docstring | ⚠ **Already self-contradictory today**: `:5-6` and `:10-13` say *"returns the input by reference (not copied)"*, while `test_apply_returns_fresh_copy` (`:64-85`) asserts the opposite and its own docstring explains the Wave-7 flip. **Pre-existing falsity — fix while in the file.** |
| `tests/geometry/test_bc_errors.py:124-127` | docstring cites `diffusion/boundary_realizer.py:219` / `:229` | **stale line numbers** (actual: `:307` periodic, `:285` prescribed) |

### 6.3 Green and STAYS green (no action)

| File:line | Test | Class | Why unaffected |
|---|---|---|---|
| `tests/diffusion/test_boundary_realizer.py:300-305` | `TestRefusals::test_periodic_refused` — `match="OPPOSITE face"` | behavioural | Diffusion keeps refusing; predicate is `isinstance(law.geometry_map, SpatialWrap)`, untouched. |
| `tests/diffusion/test_boundary_realizer.py:391-396` | `TestComposition::test_composition_refusal_propagates` — `0.5*PeriodicBoundary() + 0.5*VacuumInflow()` through `realize_recursively` | behavioural | same |
| `tests/geometry/test_bc_errors.py:146-150` | `BoundaryError` base-class contract, periodic arm | behavioural | same |
| `tests/geometry/test_boundary_factors.py:83` | `(PeriodicBoundary(), SpatialWrap, ScalarResponse, 1.0, "periodic")` | behavioural (factor pair) | The FACTORS do not change at B3.4c — only the realization does. |
| `tests/geometry/test_reemission_closure.py:529` | `("periodic", PeriodicBoundary(), False)` — no re-emission closure | behavioural | |
| **`tests/geometry/test_reemission_closure.py:1134-1141`** | `test_the_bare_and_non_pairing_laws_do_not_permute` — `assert law_permutes_ordinates(PeriodicBoundary()) is False`, commented *"spatial wrap, not angular"* | **behavioural — THE TRIPWIRE GUARD** | **This is the test that defends §3.4.** It already exists and already carries the WHY. Do not weaken it; cite it in the carve. |
| `tests/geometry/test_bc_universal_invariants.py:709-721` | `NoSource` certifies masklessly for periodic | behavioural | |
| `tests/geometry/test_boundary.py:861` | `PeriodicBoundary.key == "periodic"` | **API smoke** | |
| `tests/geometry/test_boundary_trace_law.py:265-274` | unified-registry membership | API smoke | |
| `tests/sn/operators/test_operator_block_role.py:189` | `"periodic": PeriodicBoundary()` → advertises `BlockRole.BOUNDARY` | behavioural | Role is unchanged by the carve. ⚠ but see RISK R-6. |
| `tests/numerics/test_operator_capability_predicates.py:67`, `:192` | `("PeriodicWrap", PeriodicWrapOperator(), False, True, STRUCTURAL_ABSENT)` — not invertible, adjointable | behavioural | Only if the operator's predicates change (§3.5 option (b)). |
| `tests/numerics/test_periodic_wrap_operator.py` (5 tests) | identity body, fresh copy, TP composition, predicates | behavioural (primitive) | The BODY is correct — `test_the_identity_body_is_correct_on_the_partner_half_trace` proves it. **Keep all five.** |
| `tests/sn/operators/test_snmesh_realizer_wiring.py:451-468` | registry is exactly `{vacuum, reflective}`; `BC("periodic")` raises | behavioural (admission) | **Stays green unless B3.4c also does #189.** If the carve admits `"periodic"`, `:456` goes RED — a deliberate, wanted red. |
| `tests/sn/operators/test_boundary_conditions.py` (1 hit) | mention only | — | |
| `tests/geometry/_generate_bc_equivalence_snapshots.py:661-675` | the `periodic_lebedev17` snapshot GENERATOR, `partner_face="xmax"`, `reference=_periodic_image` | **generator, not a test** | Already re-anchored (`87a65967`). Do NOT regenerate — the artefact is the derived reference the gate measures against. |
| `tests/sn/sweep/test_sweep_acyclicity.py:52` | `("periodic", "vacuum", False)` | behavioural (algebra of record) | Uses `orpheus/derivations/discrete/sn/sweep_acyclicity.py`, which **already models periodic as an opposite-face wrap** (`:361-366`). See RISK R-7 — this is a ready-made oracle. |

### 6.4 Test migration summary

- **Rewire (behavioural, contract changes):** `test_boundary.py:513`,
  `test_sn_boundary_realizer.py:696`.
- **Un-xfail (the todo list):** `test_bc_equivalence_snapshot.py:773`,
  `test_b3_domain_narrowing.py:312`.
- **Retire with the last consumer:** `tests/geometry/test_boundary.py::_realize_for_sn`
  (`:56-65`, self-declared "PERIODIC ONLY").
- **Delete nothing else.** No API-smoke test here is made redundant.

---

## ⚠ RISKS

**R-1 — `apply_transpose` writes to the WRONG FACE and nothing guards it.**
`_reflect_trace`'s transpose branch (`orpheus/sn/operators/boundary.py:438-466`)
scatters into `face`, but `Bᵀ` of a cross-face block belongs in the PARTNER's
slot. `SNBoundaryOperator.is_adjointable` returns `True` for periodic (the
realized `PeriodicWrapOperator & IdentityOperator` is adjointable), so the wrong
transpose is reachable with **no raise**. The code's own §`:361-370` warning says
off-diagonal permutation laws are *bit-identical under either spelling* and were
caught only by het-VACUUM grid reciprocity. **A reflective-fixture reciprocity
gate will NOT catch this.** Needs a periodic-fixture G-adjoint/reciprocity arm.

**R-2 — The latent `_law_from_tag` orientation bug (§1.3).** `PeriodicBoundary`
carries `axis` and sits in the parameter-free `law_cls()` fall-through
(`orpheus/transport/method.py:347`), against that function's own explicit warning.
Admitting `"periodic"` (#189) without adding the arm installs `SpatialWrap("x")`
on every face. **Fix in the same commit that makes periodic real.**

**R-3 — `face_opposite` returns a non-existent face on curvilinear meshes
(§4.4).** A solid radial axis has one endpoint (`"outer"` → `"xmax"`);
`face_opposite("xmax") == "xmin"` is not in `sn_mesh.bc` or the trace layout.
Any `mesh.bc[face_opposite(f)]` raises `KeyError`. The primitive is right to
know nothing about inventories — **the guard must live at the realizer**, and
B3.4c must state whether periodic is Cartesian-only. Model the refusal on
`orpheus/diffusion/boundary_realizer.py:307-316`.

**R-4 — The `permutes_ordinates == False` tripwire (§3.4).** The ENTIRE
Gauss-Seidel split's safety rests on `SpatialWrap.permutes_ordinates → False`
(`orpheus/geometry/boundary/_factors.py:507-510`). It is easy to "fix" this to
`True` while making periodic real — periodic *is* an ordinate bijection, and the
prose around it says "pairing". Flipping it puts periodic into
`_reflective_faces` → `B_lower` → a cross-face read inside a triangular forward
substitution. `tests/geometry/test_reemission_closure.py:1134-1141` already
guards it; **do not weaken that test, and add a comment at the flip site.**

**R-5 — Two un-repointed face-name renders + one hand-parse the brief missed
(§4.2).** `orpheus/sn/loss_representation/__init__.py:730-733` and `:739-742`
(`_inflow_faces`/`_outflow_faces` — the second is *literally* `face_opposite ∘`
the first), and `orpheus/sn/solver.py:1462-1463` (`face == "xmin"` sign test +
a `{"x":0,"y":1,"z":2}` literal). Leaving them means the carve claims a
single-source-of-truth it did not achieve — a Cardinal-Rule-2 twin.

**R-6 — `SNMethodSpace.minimal` no longer suffices, and the docs say so.**
`docs/theory/foundations/boundary_conditions.rst:964-968` states minimal works
for periodic *"precisely because [it has] not yet been narrowed."* After B3.4c a
realized periodic law needs to know its installation face AND its partner's
trace. Every test building periodic via `SNMethodSpace.minimal(quad)` —
`test_sn_boundary_realizer.py:713`, `test_boundary.py:_realize_for_sn`,
`test_operator_block_role.py`, `test_b3_domain_narrowing.py` — needs a
face-aware space. **`test_operator_block_role.py:189` may go red as collateral,
not by design.**

**R-7 — A ready-made oracle is being ignored.**
`orpheus/derivations/discrete/sn/sweep_acyclicity.py:361-366` **already models
periodic correctly** as an opposite-face wrap, and `tests/sn/sweep/test_sweep_acyclicity.py:52`
already pins `("periodic","vacuum") → cyclic`. The algebra of record is ahead of
production. Use it as the structural cross-check for the carve rather than
inventing a new one — but note it uses a `("left","right")` vocabulary
(mirror-not-import), so **do not** repoint it at `face_layout`.

**R-8 — Sphinx `-W` will not catch the doc blast radius.** The dangling
`PeriodicBoundary.apply` xref (`orpheus/numerics/operator.py:2638`) and every
`:class:`/`:meth:` in §5 are Python-domain roles — silent without `-n`
(nitpicky). Nexus flags the dangler as an `external` node of degree 1; grep the
symbol names instead of trusting the build gate. Also: `docs/theory/methods/monte_carlo.rst:592-617`
documents MC's INDEPENDENT periodic BC — keep any sweep away from it.

**R-9 — Three falsities in `PeriodicWrapOperator`'s docstring (§2.2), one of
which invents a mechanism.** `orpheus/numerics/operator.py:2622-2624` claims *"the
SN sweep handles the spatial wrap via **its own face-pair indexing**"* — that
mechanism does not exist anywhere in `orpheus/sn/`. This is the exact delegation
shape that survives refactors because it is unfalsifiable by the reader. It is
the sentence B3.4c makes true; make sure it is REWRITTEN, not left standing as
an accidental prophecy.

**R-10 — The tree is moving (§0).** Six production files were modified *during*
this audit. Re-run `git status` + `git diff --stat` before acting on any line
number here.

---

## Nexus vs grep — where they disagreed

| Question | Nexus | grep | Trust |
|---|---|---|---|
| `PeriodicBoundary` consumers | 4 `calls` edges (all tests, all `PeriodicBoundary()` constructions) — **misses the `isinstance` arm in `realizer.py:745`, which is the single most important production site** | 43 lines / 18 files | **grep** — `isinstance` is not a `calls` edge (L-009 shape: the graph under-scopes a class used as a type test) |
| Doc pages documenting it | **3** — incl. `doc:api/diffusion_1d` | **2** — misses `api/diffusion_1d` | **Nexus** — the third page reaches the class through an autodoc'd *module docstring*, so the string is in no `.rst` at all |
| Dangling xref | flags `…PeriodicBoundary.apply` as an `external` node, degree 1 | finds the line but cannot say it dangles | **Nexus** |
| `PeriodicWrapOperator` production callers | `SNBoundaryRealizer.realize` only — **correct and complete** | same | either (they agree) |
| Face-name transcription sites | no help (string literals / f-strings are not nodes) | 7 sites | **grep**, exclusively |

---

## 7. 🔴 LATE UPDATE — reconciliation against the FINAL tree state

The B3.4c carve was **executed during this audit**. Files that were clean at my
first `git status` and modified by my last:

```
 M orpheus/geometry/boundary/_factors.py            +128
 M orpheus/numerics/face_layout.py                  +156
 M orpheus/numerics/spaces/angular_trace_space.py    -52/+…
 M orpheus/sn/boundary/angular.py                    +21
 M orpheus/sn/boundary/realizer.py                  +116
 M orpheus/sn/loss_representation/sweep_schedule.py    +4
 M orpheus/sn/operators/boundary.py                  +97
 M orpheus/transport/method.py                        +13
                                    511 insertions, 76 deletions
```

### 7.1 What the carve built (verified by probe, not by diff-reading)

| Mechanism | Where | Verified |
|---|---|---|
| The face-name bijection | `orpheus/numerics/face_layout.py` — `face_name` / `face_normal` / `face_opposite` / `FACE_NAMES` | `face_opposite("xmin")=="xmax"`; `face_normal("zmax")==(2,1)`; `face_normal("outer")` → `ValueError` |
| **`BoundaryGeometryMap.domain_face(face) -> str`** — the new Protocol member, *"names whose `Γ₊` the law consumes"* | `orpheus/geometry/boundary/_factors.py` | `SpatialWrap("x").domain_face("xmin") == "xmax"`; `IdentityMap`/`SpecularMirror` return `face` |
| Axis-mismatch guard | `SpatialWrap.domain_face` | `domain_face("ymin")` on an `axis="x"` wrap → `BoundaryError: SpatialWrap declares axis='x' but is installed on face 'ymin', which lies on axis 'y'…` |
| **`SNBoundaryOperator._face_domains`** — the block `(row, column)` index, **certified a PERMUTATION of the boundary faces** | `orpheus/sn/operators/boundary.py` | refuses a half-declared pair (`xmin` periodic / `xmax` vacuum) and a partner absent from the trace |
| Forward leg reads the domain face | `_reflect_trace`: `face_in = boundary.face_view(domain_face)`, `gamma_out = trace.outflow_restriction(domain_face)`, writes `out_boundary.face_view(face)` | read |
| **Transpose leg scatters into the domain face** | `_reflect_trace`: `out_boundary.face_view(domain_face)[...] = gamma_out.apply_transpose(law.apply_transpose(restricted))` | read |
| `SpatialWrap.is_adjointable` flipped `False → True` | `_factors.py` | probe: `True` |

### 7.2 Risk reconciliation

| Risk | Status | Evidence |
|---|---|---|
| **R-1** transpose writes the wrong face | ✅ **CLOSED** | transpose now writes `out_boundary.face_view(domain_face)`. Its comment states the safety argument explicitly: *"Whole-slot assignment stays safe because `_face_domains` is certified a permutation of the faces, so no two blocks scatter into one slot."* — the same permutation certificate that guards the forward leg. |
| **R-2** latent `_law_from_tag` orientation bug | 🔴 **STILL OPEN** | `orpheus/transport/method.py` has arms for `ReflectiveBoundary` (`:314`), `WhiteBoundary` (`:318`), `AlbedoBoundary` (`:344`), then `return law_cls()` at **`:352`**. **No `PeriodicBoundary` arm.** `PeriodicBoundary` carries `axis` and would be handed the `"x"` default on every face the day #189 admits `"periodic"` — against that function's own `.. warning::`. **This is the single highest-value remaining item.** |
| **R-3** `face_opposite` names a non-existent face on curvilinear | ✅ **CLOSED, and better than recommended** | Guarded at TWO levels: `SpatialWrap.domain_face` refuses an axis mismatch, and `_face_domains` refuses a partner *"not a boundary face at all — a curvilinear mesh carries `xmax` alone, so a wrap installed there names a partner the trace has no slot for."* The permutation certificate is strictly stronger than the per-face guard I proposed. |
| **R-4** `permutes_ordinates` tripwire | ✅ **HELD** | probe: `SpatialWrap.permutes_ordinates` is **still `False`**. Only `is_adjointable` flipped — exactly as `_factors.py` predicted. Periodic therefore stays out of `_reflective_faces` → out of `B_lower` → **the G-S forward substitution stays triangular.** The §3.4 analysis stands unchanged. `tests/geometry/test_reemission_closure.py:1134-1141` remains the guard. |
| **R-5** un-repointed face-name renders | 🔴 **STILL OPEN — re-verified at final read** | `orpheus/sn/loss_representation/__init__.py:731`, `:740` (`_inflow_faces`/`_outflow_faces` — the second IS `face_opposite ∘` the first) and `orpheus/sn/solver.py:1462-1463` (`face == "xmin"` + `{"x":0,"y":1,"z":2}`) are untouched. `git status` does not list either file. **The single-source-of-truth claim is not yet fully earned.** |
| **R-6** `SNMethodSpace.minimal` no longer suffices | ⚠ **NEEDS A TEST RUN** | `realizer.py` gained +116 lines; whether `minimal(quad)` still realizes periodic is now an empirical question. Run `tests/sn/operators/test_sn_boundary_realizer.py::TestRealizePeriodic` and `tests/sn/operators/test_operator_block_role.py` before assuming. |
| **R-7** the algebra-of-record oracle | ⬜ unchanged | `orpheus/derivations/discrete/sn/sweep_acyclicity.py:361-366` still the ready-made structural cross-check. |
| **R-8** Sphinx `-W` blindness | ⬜ unchanged | the doc repairs in §5 are still owed |
| **R-9** `PeriodicWrapOperator` docstring falsities | 🔴 **STILL OPEN — re-verified** | `orpheus/numerics/operator.py` is NOT in `git status`. The *"the SN sweep handles the spatial wrap via its own face-pair indexing"* claim (`:2622-2624`) is now **doubly false** — the mechanism still doesn't exist in the sweep, and the wrap now lives in `_face_domains` instead. The dangling `PeriodicBoundary.apply` xref (`:2638`) also survives. |

### 7.3 The remaining checklist

1. 🔴 **`_law_from_tag` periodic arm** (`orpheus/transport/method.py:352`) — R-2.
2. 🔴 **`orpheus/numerics/operator.py:2615-2654`** — rewrite `PeriodicWrapOperator`'s
   docstring (3 falsities + 1 dangling xref); the type's role changed and its
   prose still describes the placeholder era.
3. 🔴 **R-5's three un-repointed sites** — `loss_representation/__init__.py:731,740`
   and `solver.py:1462-1463`.
4. 🟡 **Un-xfail the two todo-list rows** — `tests/geometry/test_bc_equivalence_snapshot.py:773`
   (delete `_B34C_XFAIL`, `:678-691`) and `tests/sn/operators/test_b3_domain_narrowing.py:312`
   (`True → False`, delete `:259-268`).
5. 🟡 **Rewire the two contract-flipping tests** — `tests/geometry/test_boundary.py:513`
   (and retire `_realize_for_sn`, its last consumer) and
   `tests/sn/operators/test_sn_boundary_realizer.py:696`.
6. 🟡 **The block-diagonality prose (§3.1, 15 sites)** — BD-1/2/5/6/8/11/12/13 are
   now falsified in the tree. Replace "block-diagonal ⟹ exact restriction" with
   the honest reason: *`faces=` selects block ROWS and each row is computed from
   the whole input trace* (§3.2). BD-3/4/7/9/10/15 must NOT be touched.
7. 🟡 **The doc corpus (§5)** — the must-repairs are
   `foundations/boundary_conditions.rst:940-953`, `:964-968`, `:3098-3101` and
   `methods/sn/boundary_conditions.rst:247-251`, `:488-491`, `:526`, `:544`.
8. 🟢 **New surface deserving its own doc paragraph:** `domain_face` and
   `_face_domains` are a genuinely new concept (a boundary law's block COLUMN,
   and the permutation certificate over the face set). Neither theory page
   mentions them yet.
