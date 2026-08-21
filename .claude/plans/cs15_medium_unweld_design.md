# CS1.5 design — the Medium un-weld, grounded against the tree

**STATUS: ⏸ SUPERSEDED-AS-CHARTED (2026-08-20, same day) — DEMOTED TO
CONTESTANT.** The user re-ruled the campaign order (campaign plan §0 ruling 3's
⚠⚠): the kernel carve (CS4a) precedes the medium work, because this design's
hardest questions — the carrier surface (`mat_map`/`spatial_shape`/
`volume_measure` homes), Medium-as-carrier vs pullback-union — are interface
design for demands construction-time binding dissolves. This file remains (a)
the grounded `[M]` fact base (§2 — still authoritative, census-backed) and (b)
the INCUMBENT contestant for the joint design rounds
(`kernel_and_medium_objectives.md`). The four fork rulings of 2026-08-20 demote
with it; fork 4's regime caveat survives as a standing constraint. The
in-flight test-architect plan against this design is PARKED on the same
grounds. Original status: "⏹ GROUNDED + FORKS RULED (2026-08-20). Next:
test-architect → compact → execute." Charter = `cs1_energy_space_design.md` fork Q11 (ruled round 6,
delegated) + its §B2 (rounds 4–6). This file is the grounding pass + round-7
re-derivation Q11 said the charter lacked. Census of record:
`scratch/cs15_grounding_census.md` (601 lines — call censuses with filters +
denominators, the promise ledger, precedent adjectives). Campaign plan of
record: `space_and_kernel_binding_campaign.md` §2.5.

**Branch context:** lands on `feature/cs1-energy-space` (BRANCH-HOLD ruling,
campaign plan §2 end of COMPACTION POINT #2 — no merge to main without the
user's go).

---

## 1. The charter, restated (what Q11 ruled — none of this re-litigated)

**Goal.** The homogeneous problem's carrier stops faking a mesh: the physical
medium (what fills space) gets its own type ABOVE the discretization, the
homogeneous solver poses on it, and the fake `AxisMesh(edges=[0,1])` mint dies
by construction — in CS1.5, nowhere else.

**Ruled content** (Q11 + §B2 round 5): mint `Medium` above `MaterialMesh`
(minimal start); `MaterialMesh` keeps its name as the medium × mesh PULLBACK;
typed union `spatial_structure: MeshBackedRegions | QuotientPoint` inside the
pullback (geometry extras type-absent on the quotient member); named binding +
mesh-conformity guard; homogeneous re-homes to
`Medium.infinite_homogeneous(mix)`; `bulk_space` migrates carrier→Medium (one
planned move). Round-6 law: symmetry lives at Geometry and only shrinks down
the lattice; the quotient is licensed by the MEDIUM's surviving symmetry.

## 2. Grounding facts — `[M]` 2026-08-20 at `a5d59425` (mine + the census, convergent)

| # | fact | provenance |
|---|---|---|
| F1 | **The tree already ships the structured member's data**: `StructuredGeometry` (frozen) = `SLB\|CYL\|SPH` + `regions: tuple[Region,...]` + `bcs`; `Region = (mat_id, outer_thickness_cm)`. Docstring rules **"No infinite-medium kind"** and documents the loose-pair idiom (solvers receive `(geometry, materials dict)`, resolved at solve time). | `geometry/structured_geometry.py:146,:215`; census §6 |
| F2 | `Mesh1D.from_geometry(geom, region_meshes)` **already conforms by construction** (per-region edge accumulation + mat_id broadcast). §B2's "edges + mat_ids coupled by hand" is true of the RAW `Mesh1D` ctor only. | `geometry/mesh.py:466`; census §2d |
| F3 | `MaterialMesh` = plain class; `_init_data` = the ONE construction body; **all 10 spatial-attribute writers sit inside it** (annotation-aware filter over the 3 hierarchy files; a space-only filter missed the annotated assignments and was corrected against the known `:178` writer). `from_materials` = `cls.__new__` + fake `AxisMesh(edges=[0,1])` at `:286`. | `material_mesh.py:132-291` |
| F4 | **Complete homogeneous-path carrier surface**: `MaterialXSField` reads `{mat_map, ng, spatial_shape}` (9 sites); solver reads `bulk_space` (`:151,:205`) + hands the carrier to `MultiplicationOperator.from_mesh` (`:153`); `IntegratedReactionRate` reaches `volume_measure` via `CrossSectionField.mesh`; plus `material_xs_field()`. The quotient member must additionally supply `spatial_shape=(1,)`, `ndim=1`, `nx=1`, `volumes=[1.0]` (census §3.3 — `nx`/`volume_measure` are members the charter did not name). | this session + census §3 |
| F5 | `from_materials` call census: **1 production** (`homogeneous/solver.py:191`) + **16 test sites in 7 files** + **2 doc refs** (`spaces.rst:1030` prose already teaches the CS1.5 reading; `infinite_medium.rst:1115` is a `:meth:` xref that goes SILENTLY dead if unre-pointed). | census §1 |
| F6 | ⭐ **`mesh is None` is an ambiguous sentinel**: it holds for the quotient carrier AND genuine axis-native d≥3 meshes. `DiffusionMesh`'s typed refusal keys on exactly it (`diffusion/augmented_mesh.py:211-213`); **`SNMesh` has NO guard** — promotion of the meshless carrier dies as a bare `AssertionError` (plain) / deep `AttributeError` under `python -O` (`augmented_mesh.py:322`, `[M]` executed). The union must discriminate what the sentinel cannot. | census §1,§3 |
| F7 | **Type-absence verdict: production-clean.** No production consumer reads `.coord`/`.axis_widths`/`.areas`/`.axes`/`.reduced` on the degenerate carrier (every external reader sits behind a method layer that excludes it; `.reduced` is SN-layer machinery, not a base attribute). Test-side: `tests/transport/test_material_mesh.py:164-201` pins the SENTINEL (`mm.mesh is None`) — rewrites with the carve. Corrections to my brief: `axis_widths` has external MESHED readers (diffusion + SN radial characteristic — not internal-to-volumes as a flat claim); `mesh.reduced.coord` reads the `ReducedStreamingOperator`'s own field, not the carrier's. | census §3 |
| F8 | Fabricated readings remaining after CS1 made the SPACE honest: `coord=CARTESIAN`, `axis_widths=[1.0]`, `volume_measure` node at `0.5` (only the weight is consumed), and `.areas` raising **with a wrong message** ("2-D meshes") on the meshless path. | census §3; `material_mesh.py:495-505,:521-525` |
| F9 | **No type spells infinite/translation-invariant** (`SubgroupOfO3` = point groups; `RigidMotion` = E(d) elements; `close_group` refuses infinite groups) ⟹ the marker is minted fresh. Names `Medium`/`MeshBackedRegions`/`QuotientPoint`/`InfiniteHomogeneous` all collision-free. Precedents (adjectives verified): `AlbedoBoundary` = the registered-frozen model; the closed 2-member union's structural model is the UNREGISTERED `AxisMesh`/`RadialAxisMesh`-under-`Axis1D`-Protocol family; `EnergyGrid` = frozen `eq=False` identity-object model. | census §4,§5 |
| F10 | Layer contract: `transport` (L2) legally imports geometry + data + numerics; (input) packages cannot import transport ⟹ `orpheus/transport/medium.py` is the lowest legal home. | `tests/test_layer_imports.py:37-72` |
| F11 | The medium→medium morphism already in the tree: `Solution.homogenize` → `MaterialMesh(coarse.with_distinct_cell_ids(), homogenized)` — spelled at the PULLBACK level. Out of CS1.5 scope; the future consumer of Medium's morphism story. | census §6.4 |
| F12 | **Promise ledger** (census §7): 7 `CS1.5` comments in the tree; §T-R's D6/D7 re-point comments were promised and NOT landed (3 of 5) — D6 re-points mechanically; D7 monkeypatches `type(sn).bulk_space` (self-guarding poison, but the target moves). Discharge list = census §7's 8 items. | census §7 |

## 3. Rulings (user, 2026-08-20 — round 7 forks)

1. **Reuse `StructuredGeometry`** as Medium's structured member. Medium = the
   NAMED PAIRING `(spatial description) × (id → Mixture)` — the loose pair
   every solver already passes. Own region storage would be a Pattern-2 twin.
2. **Pullback with typed union** (charter letter): the quotient carrier stays a
   `MaterialMesh` whose `spatial_structure` is `QuotientPoint`; built via
   `MaterialMesh.from_medium`. The "Medium is its own carrier" alternative was
   set aside — structural reason: it leaks cell bookkeeping (`mat_map`,
   `spatial_shape`) up into the physics type, re-welding at the Protocol; it
   also widens `MaterialXSField.mesh` and `from_mesh(mesh=...)` annotations.
3. **Land BOTH binding arms now** (quotient + meshed-with-conformity-guard),
   with the guard's §6c witness named in the step's done-when. The meshed arm
   ships production-unreached today (`[M]` F2/F6 idiom) — accepted; Q11's own
   reason stands (CS2 binds against the complete lattice from day one).
4. **Refuse present-but-content-unequal energy grids at Medium construction**
   — ⚠ WITH THE USER'S REGIME CAVEAT, to be carried into the docstring + theory
   page: the refusal is right for COARSE multigroup data (GENDF-class — much
   structure already discarded, so "everyone must play by the same rules"), and
   **not universal**: for FINE data (ACE-class, typical for Monte Carlo)
   unequal group structures are the legitimate result of representing all
   representable structure with the fewest groups, and the modern
   reconciliation is a UNIONIZED energy grid (memory traded for speed). The
   invariant is a property of the multigroup-deterministic regime this Medium
   serves, not of media in general; a future fine-data/MC medium needs a
   different arm (or the unionization morphism as its coherence-restoring
   constructor). `eg=None` members keep mixing freely (synthetic axis).
5. Home: `orpheus/transport/medium.py` (F10; conventional default, not asked).

## 4. The design (round 7, rulings applied)

### 4.1 `Medium`

Frozen dataclass (identity semantics — `eq=False`, the `EnergyGrid` precedent;
content identity lives in the SPACES it mints, per CS1's machinery):

* `spatial: InfiniteHomogeneous | StructuredGeometry` — `InfiniteHomogeneous`
  is a small frozen marker minted fresh (F9), full-symmetry semantics in
  prose; group machinery deferred to CS2's quotient constructor (round-6 law).
  Reuse honors `StructuredGeometry`'s own "no infinite kind" ruling: the
  infinite case lives on Medium's union, never as a geometry tag.
* `materials: dict[int, Mixture]`.
* Construction invariants (parse, don't validate): assigned ids covered
  (per-region for structured; key present for infinite), uniform `ng`
  (mirrors `InconsistentMaterialsError`), **eg-coherence refusal** per ruling
  4 (each arm gets its own mutation, vv#17).
* `Medium.infinite_homogeneous(mix)` classmethod.
* `bulk_space`: the uniform formula via ONE shared helper (energy arm from
  assigned/reachable materials; NODAL spatial axis from (shape, volumes) with
  the SAME `"spatial"` label — the label is part of the S3 digest, so
  single-sourcing the helper is what makes `medium.bulk_space ==
  carrier.bulk_space` hold by the derived-name identity). Total on both
  members: infinite ⟹ Energy ⊗ quotient point (counting); structured ⟹
  per-REGION space, weights = region volumes via the existing
  `compute_volumes_1d(coord, interface_positions)` — ~5 lines riding shipped
  primitives, avoids a partial property. (If the user prefers strict
  "nothing-until-consumer", the structured arm can raise instead — flagged,
  not blocking.)

### 4.2 The pullback union

`MaterialMesh.spatial_structure: MeshBackedRegions | QuotientPoint`, members
modeled on the `AxisMesh`/`RadialAxisMesh` family (frozen dataclasses under a
small Protocol/union, unregistered — F9):

* `MeshBackedRegions`: wraps today's `axes` tuple + optional legacy `mesh`
  adapter + the derived `_volumes/_areas` block. `QuotientPoint`: no edges, no
  coord, `volumes=[1.0]` by definition, honest single-atom `volume_measure`
  (weight 1.0 — the fabricated node-at-0.5 dies with the fake axis).
* `axes/coord/axis_widths/areas/mesh/volumes/nx/ndim/spatial_shape/
  volume_measure` become forwarding properties on `MaterialMesh` — legal on
  `MeshBackedRegions`, type-absent (typed `AttributeError` with an honest
  reason) on `QuotientPoint`. `[M]` F3: all writers already funnel through
  `_init_data`, and `[M]` F7: no production reader breaks.
* The `mesh is None` sentinel ambiguity (F6) dissolves: `DiffusionMesh`'s
  refusal re-keys on `QuotientPoint` (message + witness
  `tests/diffusion/test_augmented_mesh.py:165-170` re-point); `SNMesh`
  promotion gains the typed refusal its sibling already has (closing the
  `[M]` bare-assert/-O hazard at `augmented_mesh.py:322`).
* `.areas`' wrong-message path (F8) is fixed by construction: the 2-D message
  stays for `Mesh2D`; the quotient member's absence carries its own reason.

### 4.3 The re-home and the bulk_space migration

```python
medium  = Medium.infinite_homogeneous(mix)
carrier = MaterialMesh.from_medium(medium)      # quotient pullback
mat_xs  = carrier.material_xs_field()
space   = medium.bulk_space                     # ← the migration: the QUOTIENT
# space's mint moves to Medium.  The carrier's own bulk_space (per-cell,
# meshed members) stays — a Medium cannot mint a per-cell space, so
# "migrates carrier→Medium" can only mean the quotient member's mint.
```

`medium.bulk_space == carrier.bulk_space` by the derived-name identity (same
helper, same label, same axes) — `MultiplicationOperator.from_mesh`'s chain and
the explicitly-threaded isotropic pair keep agreeing on one space; the
OperatorSum guard keeps validating. **The D5 byte gate (8/8, alive until the
merge cycle) must stay bit-identical across this rewiring — the counting
theorem measured again.**

### 4.4 The named binding, both arms (ruling 3)

* **quotient arm** `from_medium(medium)` — consumer-backed (the re-home).
* **meshed arm** `from_medium(medium, mesh)` — the conformity guard's home:
  every region interface ∈ `mesh.edges` (comparison discipline — exactness vs
  tolerance — is a test-architect question; interface positions come from
  cumulative thicknesses, edges from whatever built the mesh), and the mesh's
  material assignment consistent with (or derived from) the medium's regions.
  §6c witness named in the done-when: a hand-built non-conforming `Mesh1D` is
  REFUSED with a region-naming reason — the step's first red.

### 4.5 Fences (unchanged from the charter + census)

`Mesh1D`/`Mesh2D`/`from_geometry` untouched (the mesh keeps `mat_ids`);
method-mesh constructors and solver entries untouched (re-routing through
Medium is future work); `Solution.homogenize`/condense untouched (F11); no
Spatial/Quadrature axes, no ⊕, no S3 flip (CS2); no Optional→mandatory (CS4).

## 5. Step decomposition (hardened; §6b/§6c applied)

1. **Mint `Medium`** (+ marker + shared formula helper) — gates: construction
   refusals (per-arm mutations), the identity bridge
   `Medium.infinite_homogeneous(mix).bulk_space == <today's carrier space>`
   (licenses step 2's swap), structured-member bulk_space content. No
   production consumer yet — additive, nothing moves.
2. **Union restructure, in place** — `_init_data` populates
   `spatial_structure`; forwarding properties; `from_materials` STOPS minting
   the fake axis and builds `QuotientPoint` directly (its signature/callers
   untouched — §6b clean); quotient `volume_measure` honest; DiffusionMesh
   refusal re-keys; SNMesh promotion refusal added; intrinsic suite
   (`test_material_mesh.py`) re-pins the union member instead of the sentinel.
   Byte gate 8/8 + scoped suites green.
3. **`from_medium` (both arms) + re-home + retire `from_materials`** — ONE
   step by §6b (the solver is the sole production call site; the 16 test
   sites in 7 files migrate in the same step; `infinite_medium.rst:1115` xref
   re-points; `dead_references` sweep). `space = medium.bulk_space` lands here
   (the migration); B9/B10/B11 + D6/D7 re-point (F12); conformity-guard
   witness red-then-green.
4. **Docs + promise discharge** — `spaces.rst` chartered→landed flips
   (`:1023-1030`, `:1084-1085`); the fork-4 regime caveat into Medium's
   docstring + theory page **+ a GitHub issue** (module:data, type:feature:
   fine-data/ACE-regime medium arm — unequal grids legitimate, unionized-grid
   morphism as the coherence-restoring constructor; the docstring cites it);
   the 13-page staleness surface (census §7 item 7); campaign ledger row.
   Also: the `material_mesh.py:376-377` comment ("bulk_space moves to Medium")
   is REWRITTEN to the landed truth — under the ruled design the property
   STAYS (per-cell, meshed members); only the QUOTIENT mint moves. Corollary:
   D7's poison target `type(sn).bulk_space` does NOT move (the census flagged
   it as moving under the comment's literal reading; the design supersedes).

Every step: scoped suites + pyright ratchet (terminal 1) + copy-aside mutation
battery per CS1 practice; Sphinx `-W` clean at step 4.

## 6. Open items for the test-architect

* Battery design per step (which mutations, which gates redden, which
  anti-claims stay green — the byte gate's role in steps 2/3).
* The conformity comparison discipline (exact vs tolerance; who owns the
  interface-position arithmetic).
* The union's type-absence witnesses (per-member, per-attribute — not one
  blanket gate).
* Medium invariant arms: one mutation per arm (vv#17 granularity).
* What pins `medium.bulk_space == carrier.bulk_space` AFTER from_materials
  dies (the bridge gate's post-retirement form).
