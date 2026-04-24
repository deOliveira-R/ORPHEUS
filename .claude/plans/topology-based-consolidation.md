# Plan — topology-based consolidation of Peierls reference cases

**Author**: Claude Opus 4.7, 2026-04-23.
**Relationship to prior plan**: refines Stage 2b of
`.claude/plans/continuous-references-simplification.md` by replacing
its "unify by method module" principle with "unify by topology
class." This document supersedes the shape-keyed organization in
the earlier plan; read this FIRST if the two conflict.

## 0. The structural insight

The current Peierls reference-case organization is keyed by
**shape** (slab vs cylinder vs sphere). But the concrete machinery —
which closures apply, which accuracy ceilings hold, which V&V tests
are meaningful — is determined by **topology** (number of
boundary surfaces and their curvature class), not by shape.

Two facts make this load-bearing:

**F1**. The **F.4 scalar rank-2 per-face closure** (Stamm'ler IV
Eq. 34 = Hébert 2009 Eq. 3.323) requires a 2-surface cell. It
applies equally to slab (two parallel faces), hollow cylinder (inner
ring + outer ring), and hollow sphere (inner shell + outer shell).
Its accuracy story, its limitations (scalar-mode residual that
research-log L21 says cannot be lifted via rank-N), and its L19
stability-protocol coverage are all the same class of argument.

**F2**. On **1-surface compact cells** (solid cylinder, solid
sphere), F.4 collapses to rank-1 Mark by construction (no second-
face coupling). These shapes share the same closure story (only
rank-1 Mark is shipped), the same structural barrier at thin R
(Issue #103: 21 % err at R=1 MFP), and the same path forward (rank-N
DP\_N on the single outer face).

**Current organization obscures this.** Each of
[`peierls_cylinder.py`](orpheus/derivations/peierls_cylinder.py) and
[`peierls_sphere.py`](orpheus/derivations/peierls_sphere.py) mixes
both topology classes: the module-level `GEOMETRY` constant is the
solid (1-surface) variant, the module's `solve_peierls_*_1g` accepts
an `inner_radius` kwarg that flips to 2-surface, and the case
builders (one for solid, one for hollow F.4) sit side by side.
[`peierls_slab.py`](orpheus/derivations/peierls_slab.py) is
always 2-surface but is organizationally isolated because its
machinery (E₁ Nyström) differs from the curvilinear unified core.

**The code ALREADY does topology dispatch internally**:
`CurvilinearGeometry.n_surfaces` (line 317) returns 2 for hollow cyl
+ hollow sph + slab-polar and 1 for solid cyl + solid sph. Seven
dispatch sites in `peierls_geometry.py` branch on
`geometry.n_surfaces == 2`. What's missing is the layer that makes
this visible at the **case registration** and **documentation**
levels.

## 1. The two topology classes

### Class A — Two-surface (F.4 applies)

| Member | Implementation module | Implementation machinery |
|---|---|---|
| Slab | `peierls_slab.py` | Native E₁ Nyström |
| Hollow annular cylinder | `peierls_cylinder.py` + `peierls_geometry.py` | `CurvilinearGeometry(kind="cylinder-1d", inner_radius>0)` |
| Hollow sphere | `peierls_sphere.py` + `peierls_geometry.py` | `CurvilinearGeometry(kind="sphere-1d", inner_radius>0)` |

**What members share**:
- F.4 scalar rank-2 per-face is the production white-BC closure.
- Each has distinct inner and outer P/G primitives.
- Each has a (2×2) transmission matrix `W` with kind-specific
  reciprocity (cylinder R/r_0, sphere (R/r_0)², slab flat-equal-
  area so W_oi = W_io).
- Each has a zero-absorption conservation identity
  `W_oo + W_io = W_oi + W_ii = 1`.
- Each has the same structural residual floor class (L21: cannot be
  lifted via rank-N refinement).
- Each can be verified by the L19 two-quadrature stability
  protocol (Issue #123).

**Current shipped references**:
- `peierls_slab_2eg_2rg` (2G 2-region)
- `peierls_cyl1D_hollow_1eg_1rg_r0_{10,20,30}` (1G 1-region F.4)
- `peierls_sph1D_hollow_1eg_1rg_r0_{10,20,30}` (1G 1-region F.4)

### Class B — One-surface compact (rank-1 Mark only)

| Member | Implementation module | Implementation machinery |
|---|---|---|
| Solid cylinder | `peierls_cylinder.py` + `peierls_geometry.py` | `CurvilinearGeometry(kind="cylinder-1d", inner_radius=0)` |
| Solid sphere | `peierls_sphere.py` + `peierls_geometry.py` | `CurvilinearGeometry(kind="sphere-1d", inner_radius=0)` |

**What members share**:
- Only closure shipped is rank-1 Mark (`closure="white_rank1_mark"`
  post-Stage-5); `closure="white_f4"` silently collapses to Mark.
- Accuracy is dominated by the scalar re-entry approximation; 21 %
  err at R = 1 MFP (cylinder), similar on sphere.
- Path forward is rank-N DP\_N on the single outer face (Issue
  #103) or chord-based analytical reference (Issue #101).
- Neither has shipped continuous references in the registry today
  because the rank-1 Mark floor is too loose.

**Current shipped references**: none (by design — Issue #105).

### Deliberately NOT a third class

Multi-region solid cells (`cp_cyl1D_1eg_4rg` etc.) are Class B from
the white-BC closure standpoint — the only closure acts at the
outer white surface; internal interfaces are continuity-matched, not
BC-closed. They get grouped with Class B, not carved out as a
third class.

## 2. What this makes unambiguous

After the reorganization, the following questions resolve in O(1)
lookup without tracing code:

1. *"Which closures apply to this case?"* → look up topology class.
   Class A → {vacuum, white_rank1_mark, white_f4}. Class B →
   {vacuum, white_rank1_mark}.
2. *"What is the accuracy ceiling?"* → look up the per-class
   structural-residual table. Class A carries the F.4 residual
   table. Class B carries the rank-1 Mark table (21 % at 1 MFP …).
3. *"Does the L19 stability protocol apply?"* → Class A yes
   (multiple closures to compare); Class B vacuously (only one
   closure shipped).
4. *"What new capability would this case unlock?"* → is it Class A
   (then F.4 is already there; enrichment requires rank-N Marshak
   which L21 closed) or Class B (then rank-N DP\_N or chord-based
   analytical is open research).

## 3. Concrete restructuring — staged

### Stage T1 — Topology as a first-class concept (~1 session)

**Goal**: expose topology as a queryable property on every Peierls
geometry object.

**Scope**:

1. Add a `topology` property to `CurvilinearGeometry` (new, does
   not replace `n_surfaces`):

   ```python
   @property
   def topology(self) -> Literal["two_surface", "one_surface_compact"]:
       """Topological class for Peierls closure dispatch."""
       return "two_surface" if self.n_surfaces == 2 else "one_surface_compact"
   ```

2. Add an equivalent module-level label to `peierls_slab.py`:
   `TOPOLOGY = "two_surface"`. Slab is always 2-surface; this
   makes the claim explicit without touching its E₁ Nyström
   machinery.

3. Update the Sphinx naming glossary
   (§`theory-peierls-naming`) with the topology concept as the
   **primary organizing principle**, and update the capability
   matrix (§`theory-peierls-capabilities`) to group rows by
   topology class (Class A table first, Class B table second).

**Acceptance**:
- `cyl.topology == "two_surface"` for `inner_radius > 0`,
  `"one_surface_compact"` for `inner_radius == 0`.
- Sphinx capability matrix renders with Class A / Class B
  subheadings. A reader scanning the page sees topology before
  shape.
- Zero test changes required — the property is additive.

**Risk**: trivial. Purely additive.

**Budget**: 1 session (~2 hr).

### Stage T2 — Topology-based case builders (~1 session)

**Goal**: two shared case-builder functions replace the current
~10 per-geometry builders.

**Scope**:

1. Create `orpheus/derivations/peierls_cases.py`:

   ```python
   def build_two_surface_case(
       shape: Literal["slab", "cylinder-1d", "sphere-1d"],
       ng_key: str,
       n_regions: int,
       *,
       inner_radius: float | None = None,  # None for slab; required for cyl/sph
       closure: str = "white_f4",
       quadrature: QuadratureSpec | None = None,
   ) -> ContinuousReferenceSolution:
       """Class A: slab, hollow cylinder, hollow sphere.

       Internally dispatches the machinery call:
         - shape == "slab": peierls_slab.solve_peierls_slab_1g (E_1 Nystrom)
         - shape in {"cylinder-1d", "sphere-1d"}: solve_peierls_1g via
           CurvilinearGeometry(kind=shape, inner_radius=inner_radius)

       Normalization uses geometry.shell_volume_integral OR the slab's
       flat integral; both are unified by the CurvilinearGeometry helper
       or the slab equivalent.
       """

   def build_one_surface_compact_case(
       shape: Literal["cylinder-1d", "sphere-1d"],
       ng_key: str,
       n_regions: int,
       *,
       closure: str = "white_rank1_mark",
       quadrature: QuadratureSpec | None = None,
   ) -> ContinuousReferenceSolution:
       """Class B: solid cylinder, solid sphere."""

   def cases() -> list[ContinuousReferenceSolution]:
       """Enumerate all shipped peierls continuous references."""
       return [
           build_two_surface_case("slab", "2g", 2),
           build_two_surface_case("cylinder-1d", "1g", 1, inner_radius=0.1),
           build_two_surface_case("cylinder-1d", "1g", 1, inner_radius=0.2),
           build_two_surface_case("cylinder-1d", "1g", 1, inner_radius=0.3),
           build_two_surface_case("sphere-1d", "1g", 1, inner_radius=0.1),
           # ...
       ]
   ```

2. The per-geometry builders in `peierls_{slab,cylinder,sphere}.py`
   (`_build_peierls_slab_case`, `_build_peierls_cylinder_hollow_f4_case`,
   etc.) become **thin delegates** calling into the topology-class
   builder. They are retained as backward-compat APIs; future
   code should call the topology-class builder directly.

3. `continuous_cases()` in each per-geometry module delegates to
   `peierls_cases.cases()` filtered by shape. Auto-discovery (from
   Stage 2.4 of the prior plan) then picks up `peierls_cases` in
   its walk.

**Acceptance**:
- Same 7 shipped references produced by `peierls_cases.cases()`,
  bit-exact k_eff versus pre-migration.
- Per-geometry builders become ~30-line wrappers instead of
  ~150-line hand-codes. Net LOC reduction: ~300–400 lines.
- Regression tests from `test_peierls_rank2_bc.py` pass unchanged.

**Risk**: medium. The delegation needs to preserve slab's E₁
Nyström dispatch correctly (slab is the odd geometry out). Single-
session regression test required.

**Budget**: 1 session (~4 hr). Depends on Stage T1.

### Stage T3 — Topology-organized Sphinx (~0.5 session)

**Goal**: the Sphinx capability matrix reads by topology class
first, with shape as a sub-axis.

**Scope**:

1. Restructure §`theory-peierls-capabilities` into two sub-sections
   (Class A, Class B), each with a narrative introduction then the
   list-table.
2. Add a topology diagram (ASCII or simple SVG) showing how the
   five geometry+cavity combinations map to two classes.
3. Cross-link to the `topology` property on `CurvilinearGeometry`
   as the runtime check.

**Acceptance**:
- A fresh reader answers *"which closures apply to hollow
  cylinder?"* by reading one paragraph about Class A and skipping
  Class B entirely — without mentally filtering across shape
  boundaries.

**Risk**: trivial.

**Budget**: 0.5 session (~2 hr). Independent of Stage T2 but
natural to land together.

### Stage T4 — Test organization by topology (optional)

**Goal**: tests that exercise the same closure class share a
topology marker.

**Scope**:

1. Add `@pytest.mark.topology_class("two_surface")` to `test_hollow_*`
   tests; `@pytest.mark.topology_class("one_surface_compact")` to
   the rank-1 Mark tests on solid geometry.
2. Register the markers in `conftest.py` / `pyproject.toml`.
3. Allow running `pytest -m 'topology_class("two_surface")'` to
   exercise only F.4-class tests. Useful for bisection when a
   Class A change lands.

**Acceptance**:
- Tests can be filtered by topology class from the command line.
- The V&V audit report groups tests by topology.

**Risk**: low. Additive.

**Budget**: 0.5 session, optional.

### Stage T5 — Unified solver entry by topology (optional, future)

**Goal**: ultimate user-facing simplification — one solver function
per topology class, not per shape.

**Scope**:

```python
def solve_two_surface(shape, radii, xs, *, inner_radius=0, closure="white_f4", ...):
    """Slab / hollow cyl / hollow sph."""

def solve_one_surface_compact(shape, radii, xs, *, closure="white_rank1_mark", ...):
    """Solid cyl / solid sph."""
```

This is more aggressive than Stage 4 in the prior plan (one
`solve_peierls()` for all shapes) — it deliberately separates the
two topology classes because they have **different valid closures**.
A user calling `solve_two_surface` knows F.4 is available; a user
calling `solve_one_surface_compact` knows it isn't.

**Acceptance**:
- Per-shape wrappers (`solve_peierls_cylinder_1g`, etc.) become
  deprecated in favor of topology-class wrappers.
- The `closure` parameter's valid values are fewer per-function
  (`white_f4` is absent from `solve_one_surface_compact`'s
  signature, surfacing the F1/F2 facts structurally).

**Risk**: medium. User-facing API change across many callers.

**Budget**: 1 session, optional. Gated on Stages T1-T3.

## 4. Ordering + dependencies

```
T1 (topology property)     — landable alone, docs-only risk
  ↓
T2 (topology case builders) — depends on T1 + Stage-2 auto-discovery
  ↓
T3 (topology Sphinx)       — depends on T1; naturally lands with T2
  ↓
T4 (test markers)          — independent, landable alongside T2 or later
  ↓
T5 (topology solvers)      — depends on T2; breaking API change, opt-in
```

**Minimum-viable topology consolidation = T1 + T2 + T3** (two
sessions).

## 5. Why this is NOT a replacement for the simplification plan

The prior plan
(`.claude/plans/continuous-references-simplification.md`) handled:

- **Stage 1** (capability matrix) — still relevant; Stage T3 just
  restructures it.
- **Stage 2** (helper methods + auto-discovery) — already landed;
  Stage T2 uses auto-discovery as the pick-up mechanism.
- **Stage 5** (closure-dispatch unification) — already landed; the
  new canonical names `white_rank1_mark` / `white_f4` map naturally
  onto Class B / Class A respectively.
- **Stage 3** (`ReferenceCase` type unification) — orthogonal to
  topology; separate concern.

Topology is the **organizing principle** that the earlier plan's
Stage 2b gestured at but did not articulate. This plan supersedes
that sub-stage.

## 6. Open questions

**TQ1 — Does slab belong in Class A with the curvilinear
hollow cells?**
Yes, by the F1/F2 argument: slab has 2 surfaces, shipped closure is
F.4 (in its slab form — E_2/E_3 bilinear), residual is bounded by
scalar-mode quadrature floor. The shape differs (flat vs curved)
but the closure class is the same. **Decision: include slab in
Class A**. If Stage T2 execution reveals a hard incompatibility
(e.g., slab's E_1 Nyström doesn't factor through
`CurvilinearGeometry`), slab retains its own machinery but the
**case registration** still lists it under Class A.

**TQ2 — What about periodic / albedo / specular BCs?**
Not shipped today, but when they land: albedo (generalized white)
is naturally Class A or B depending on number of faces; specular
(with preferred direction) breaks the SO(2) symmetry on curved
surfaces and may warrant its own class. Decision deferred to that
time.

**TQ3 — Does Class B ever register references?**
Not today. For the shipped `cp_{cyl,sph}1D_*` solver tests to
acquire a Peierls reference, either rank-N DP\_N lifts the 21 % err
floor (Issue #103) or an analytical chord-based solver ships
(Issue #101). Class B's registration is structurally a placeholder
until one of those lands. Documenting it explicitly as Class B (vs
"empty registry for cylinder / sphere") makes the blocking
condition visible.

**TQ4 — Should the Sphinx page be renamed?**
`peierls_unified.rst` is currently keyed on "unified machinery." If
topology becomes the organizing principle, the page might more
accurately be named `peierls_by_topology.rst`. Decision: leave the
filename as-is to avoid cross-references breaking; clarify the
framing in the Key Facts.

## 7. What this plan does NOT propose

- **No restructuring of `peierls_geometry.py`**. The existing 4086-
  LOC unified machinery is correct and elegant. The topology
  discriminator (`n_surfaces`) is already there. This plan
  elevates an existing internal concept to the module-organization
  layer; it does not re-engineer the math.
- **No removal of `peierls_{slab,cylinder,sphere}.py`** as modules.
  They remain as backward-compat facades and as the homes of the
  slab's native E₁ Nyström. What changes is where the case
  **builders** live (move to `peierls_cases.py`) and how the
  capability matrix is organized (by topology).
- **No changes to existing test files**. Stage T4 adds markers
  optionally; pre-existing test node IDs are unchanged.
- **No deprecation of the `kind` parameter on `CurvilinearGeometry`**.
  Shape remains a useful sub-axis within a topology class (e.g.,
  cylinder vs sphere within Class A). Topology is the primary
  axis; shape is secondary.

## 8. Acceptance for the overall reorganization

After Stages T1 + T2 + T3 land, a new developer approaching the
Peierls references should be able to:

1. Open Sphinx `peierls_unified.rst`.
2. See "Class A: two-surface (F.4 applies)" as the first
   capability section.
3. Read that slab + hollow cyl + hollow sph all share this
   closure class — without needing to understand the difference
   between E_1 Nyström and CurvilinearGeometry.
4. Jump to Class B for "solid cyl / solid sph, rank-1 Mark only"
   with the Issue #103 path forward called out.

**The original session confusion (C3: "boundary-string dispatch is
non-obvious") becomes impossible** — because the reader never has
to traverse the boundary→closure→n_surfaces chain. They start at
topology, and topology tells them which closures apply directly.

## 9. Session trail

- **2026-04-23**: filed by Claude Opus 4.7 in response to the
  user's observation that "slab, annuli cylinder and hollow
  sphere share the same unified case due to shared topology."
  The shape-keyed organization in the prior plan's Stage 2b is
  superseded by this topology-first framing.
- **2026-04-23 (later, same session)**: Stages T1 + T2 + T3
  executed. See §10 below for the landed deliverables.

## 10. Landed deliverables (Stages T1 + T2 + T3, 2026-04-23)

### Stage T1 — Topology as a first-class concept

- **Property** ``CurvilinearGeometry.topology`` added
  ([peierls_geometry.py:335](orpheus/derivations/peierls_geometry.py#L335)).
  Returns ``"two_surface"`` (for hollow cyl / sph / slab-polar) or
  ``"one_surface_compact"`` (for solid cyl / sph). Semantically
  equivalent to ``n_surfaces`` (2 / 1 respectively) but named for
  the dispatch role rather than the arity.
- **Module label** ``peierls_slab.TOPOLOGY = "two_surface"``
  ([peierls_slab.py:54](orpheus/derivations/peierls_slab.py#L54)).
  Slab does not use ``CurvilinearGeometry``; this module-level
  constant carries the same label for the registration layer. Docstring
  notes the future absorption-into-curvilinear via arbitrary-precision
  quadrature (scheduled for a dedicated numerics session, not this
  plan).
- **Sphinx Key Facts** updated with the topology framing as the
  primary organizing principle; naming glossary
  (§``theory-peierls-naming``) adds a Topology block at the top
  pointing at the runtime property.

### Stage T2 — Topology-based case builders

- **New module**
  [``orpheus/derivations/peierls_cases.py``](orpheus/derivations/peierls_cases.py)
  (~200 lines). Three public entry points:
  - ``build_two_surface_case(shape, ng_key, n_regions, *, inner_radius)``
    — Class A builder; dispatches to slab's E₁ Nyström for
    ``shape="slab"``, to the curvilinear F.4 builders for
    ``"cylinder-1d"`` / ``"sphere-1d"``.
  - ``build_one_surface_compact_case(shape, ng_key, n_regions)`` —
    Class B builder; currently raises ``NotImplementedError`` with
    an explanatory message. Reserved for when Issues #103 / #101
    resolve.
  - ``continuous_cases()`` (aliased as ``cases``) — the single
    canonical registry entry point. Returns the 7 shipped Class-A
    references (1 slab + 3 cyl hollow + 3 sph hollow) plus an empty
    Class B.
- **Per-geometry ``continuous_cases()`` delegate to empty**
  ([peierls_cylinder.py:482](orpheus/derivations/peierls_cylinder.py#L482),
  [peierls_sphere.py:459](orpheus/derivations/peierls_sphere.py#L459),
  [peierls_slab.py:685](orpheus/derivations/peierls_slab.py#L685)).
  The per-geometry ``_build_*_case`` constructors remain (they are
  the actual constructor logic that ``peierls_cases`` invokes);
  only the registration hook collapses to a single source.
- Auto-discovery (shipped in Stage 2 of the prior plan) picks up
  ``peierls_cases`` automatically; the three per-geometry modules
  still pass through the walker but contribute zero references.
- **Regression**: 9/9 hollow F.4 tests pass bit-exact in 12 s
  (``test_hollow_cyl_rank2_beats_rank1_mark`` + sphere equivalents
  + conservation tests). Core F.4 math unchanged; only the
  registration layer moved.

### Stage T3 — Topology-organized Sphinx

- ``§theory-peierls-capabilities`` restructured into two
  sub-sections:
  - **Class A — two-surface (F.4 applies)**: narrative intro
    explaining slab/hollow-cyl/hollow-sph share a closure class,
    then the per-shape closure table (4 closures × 3 shapes).
  - **Class B — one-surface compact (rank-1 Mark only)**:
    narrative intro on F.4 collapse + Issues #103/#101 path
    forward, then the per-shape closure table (5 closures × 2
    shapes, with the F.4 row marked ⚠️ silent-collapse).
- The original shape-keyed capability table removed; replaced by
  the two topology-keyed tables above.
- Sphinx build clean (only pre-existing SN-MMS orphan warnings).

## 11. Acceptance test (from §8) — met

Opening ``docs/theory/peierls_unified.rst`` now shows:

- Key Facts' first bullet after the paragraph intro: **topology is
  the primary organizing principle**, with Class A and Class B
  named explicitly.
- §``theory-peierls-capabilities`` leads with "Class A — two-surface
  (F.4 applies)" before any shape-keyed content.
- The "which closures apply to hollow cylinder?" question resolves
  in one glance at the Class A table; Class B (solid geometry)
  doesn't enter the reader's consideration.

The original C3 confusion ("boundary-string dispatch is non-obvious")
becomes structurally harder to reproduce: a reader starts at topology,
and topology determines the closure set directly.

## 12. Deferred to a dedicated numerics session

As confirmed by the user 2026-04-23: **slab's native E₁ Nyström
machinery stays separate** from the curvilinear
``CurvilinearGeometry`` polar form. Absorbing it requires adaptive
arbitrary-precision quadrature on the grazing-ray stiffness in the
polar reformulation (Gauss-Laguerre on the exp-stretched coordinate
was tried and rejected because the integrand is super-exponential,
not exponential). That is its own numerics problem — not a
registration-layer simplification. Tracked as Phase G / Issue #111.

**Self-contained plan for that session**:
[`slab-into-curvilinear.md`](slab-into-curvilinear.md) (2026-04-23,
~450 lines). Contains executive summary + prerequisites-done
checklist + the mpmath.quad vs Gauss-Laguerre decision + current-
state inventory + target math + G.1–G.5 implementation plan +
regression fixture + retirement decision + acceptance criteria. A
fresh session can execute without further cross-document reading.

Slab's current position: **registered under Class A** by
``peierls_cases.build_two_surface_case("slab", ...)``, with the
implementation routing to ``peierls_slab._build_peierls_slab_case``
(the native E₁ Nyström). When Phase G lands, the routing switches
to the unified curvilinear machinery without changing the Class A
registration.

## 13. Performance note (pre-existing, not a regression)

Building the full continuous-references registry cold from a fresh
Python process takes ~15 minutes, dominated by the 7 Peierls solves
at the hollow-F.4 builder defaults (``n_panels_per_region=3``,
``p_order=5``, ``n_theta/n_rho/n_surf_quad=24``, ``dps=20``). This
is **pre-existing**: pre-refactor, ``peierls_cylinder.continuous_cases()``
ran the same 3 cyl builds, ``peierls_sphere.continuous_cases()``
ran the same 3 sph builds, and ``peierls_slab.continuous_cases()``
ran the same 1 slab build — total 7, same cost.

The registry is cached at process scope
(``reference_values._CONTINUOUS is None`` guard at
[reference_values.py:197](orpheus/derivations/reference_values.py#L197)),
so subsequent ``rv.continuous_all()`` calls in the same process
are free. Per-test overhead is not impacted.

A separate optimization opportunity — lowering the hollow-F.4
builder quadrature defaults, or deferring construction until first
access per-name — is tracked as a future perf concern. Not in
scope for the topology consolidation.
