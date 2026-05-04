# Plan — simplify the continuous-references / derivations infrastructure

**Author**: Claude Opus 4.7, 2026-04-23.
**Motivation**: during today's cylinder + sphere F.4 port, the main
agent had to revise its own understanding of the infrastructure
**five separate times** before landing a ~400-line port. Each
revision points at a specific ambiguity in the shipped infrastructure
that should not exist. This plan is a forensic reading of those
ambiguities plus a staged migration to an unambiguous unified
surface.

**Primary metric**: a fresh reader should be able to answer the
question *"What references do we have for problem X?"* by reading
**one file** and **one Sphinx table** — no code tracing.

**Current metric** (measured from today's session): answering that
question for "hollow cylinder F.4 at r_0/R=0.3" required reading
7 files totaling ~1500 lines, plus running 3 diagnostic scripts,
plus the dispatch of an explorer agent, plus reading 3 issues. The
answer was ambiguous even then because of naming collisions.

## 0. The five confusion points (diagnosis)

Every item here is a concrete simplification target.

**C1 — "Unified framework" means two different things.**
Machinery unification (`CurvilinearGeometry` shared by cylinder +
sphere) is complete; reference-case registration is per-module and
manually wired. These get conflated in conversation and in docs.
`peierls_nystrom.rst` opens with "unified Peierls" but half the
page describes the module dispatch machinery and half describes the
volume-kernel math. No single place states "cylinder 1G hollow F.4
is available at r_0/R ∈ {0.1, 0.2, 0.3}".

**C2 — F.4 is overloaded.**
"F.4" refers to both **Phase F, subphase 4** (the 2026-04-21 Phase
F rollout that added hollow-cylinder support) AND **the closure
formula Stamm'ler IV Eq. 34 = Hébert 3.323** (which pre-dates Phase
F and is the production scalar rank-2 per-face closure). The Sphinx
phrase "Phase F.4 (hollow cylinder). Extended rank-2 white to…"
conflates the two meanings in a single sentence.

**C3 — Boundary-string dispatch is non-obvious.**
`boundary="white"` routes to `build_white_bc_correction_rank_n`
(rank-1 Mark by default). `boundary="white_rank2"` routes to
`build_closure_operator(reflection="white")` → F.4 rank-2 per-face.
On **solid** geometry `white_rank2` silently collapses to Mark (the
rank-1 Schur). The `NotImplementedError` guard at
`n_bc_modes > 1 + reflection="white"` blocks **rank-N Marshak per-face**,
NOT F.4. Knowing what closure you get requires reading
`build_closure_operator` + `solve_peierls_1g`'s boundary branches
+ understanding the hollow/solid `n_surfaces` dispatch.

**C4 — Two registries, overlapping semantics.**
`VerificationCase` (for solver tests, `rv.get(name)`) and
`ContinuousReferenceSolution` (for references, `rv.continuous_get(name)`)
are separate classes with different fields, different lookup
functions, and different module lists. A case at the SAME problem
(same materials, geometry, BCs) can exist in either, neither, or
both — no invariant enforces correspondence. The user must know
which registry to query; the code has no link from one to the
other. New reference-case naming must match by string convention.

**C5 — The registration hook is duplicated per geometry.**
Every `peierls_{cyl,sph,slab}.py` module has its own
`_build_peierls_*_case` / `_build_peierls_*_hollow_f4_case` /
`continuous_cases()` / `_MAT_IDS_*` / import of `_RADII` from its
CP sibling. The new hollow-F.4 case builders I added today
(cylinder + sphere) are **~90% identical** — same control flow,
different geometry constant, different shell-volume normalization
(2π ∫ r φ dr vs 4π ∫ r² φ dr), different reciprocity exponent.
Every future new-case registration replicates ~150 lines per
geometry.

## 1. Design principles

**P1 — One user-facing entry point per method.** The user calls
`solve_peierls(...)`, not `solve_peierls_cylinder_1g(...)` /
`solve_peierls_sphere_1g(...)`. Geometry is a parameter, not a
module name. The three per-geometry facades become deprecated
wrappers (retained for backward compat, scheduled for removal in a
Phase H cleanup).

**P2 — One case-builder per method.** `build_peierls_case(kind,
ng_key, n_regions, *, inner_radius=0, closure="vacuum", ...)` lives
in `orpheus/derivations/peierls_cases.py` (or in
`peierls_geometry.py` under a `cases` sub-namespace). The
per-geometry duplication is gone. Shell-volume normalization and
reciprocity exponent are `CurvilinearGeometry` methods, not
per-builder hand-codes.

**P3 — Closure names are physics, not implementation history.**
`closure="vacuum"` | `closure="white_rank1_mark"` |
`closure="white_f4"` (or `"white_scalar_per_face"`). The
`boundary="white"` / `boundary="white_rank2"` strings are
deprecated aliases for one release. The closure name specifies the
answer the user gets; it does NOT describe the implementation
machinery.

**P4 — One registry, unified `ReferenceCase`.** Replace
`VerificationCase` + `ContinuousReferenceSolution` with
`ReferenceCase(problem: ProblemSpec, solution: Optional[Solution])`.
`Solution` holds whatever the method computes (k_eff, φ(r), ψ(r, Ω),
reference tolerance). `rv.get(name)` returns `ReferenceCase`;
`ReferenceCase.has_continuous_solution()` replaces the separate
`continuous_get` call.

**P5 — Auto-discovery, not manual module lists.** The registry
builder walks `orpheus.derivations.*` and calls every module's
`cases()` method if defined. Removes the
`_continuous_modules = [homogeneous, diffusion, ...]` list that I
forgot to update twice today. Each module has exactly one `cases()`
entry point; the builder is a single loop.

**P6 — Capability matrix in Sphinx.** A single table in
`docs/theory/peierls_nystrom.rst` (or a new `theory/capabilities.rst`)
with rows = (method, geometry, n_regions, n_groups, hollow?,
closure) and columns = (production? reference? accuracy class?
`verifies()` label?). This is the table I failed to mentally build
in this session; making it a literal asset removes the need.

**P7 — Tolerance dataclass, not strings.** Replace
`tolerance="O(h²) + scalar-mode residual ~1.4 %"` strings with

```python
@dataclass
class AccuracyBudget:
    asymptotic_order: int        # e.g., 2 for O(h²)
    structural_residual: float   # F.4 scalar-mode limit, fraction (0.014)
    quadrature_noise: float      # e.g., 1e-4
    notes: str                   # human prose
```

V&V audit can parse it; `assert_rank_n_structural_win` can compare
against it.

## 2. Staged migration

The simplifications fall into stages of increasing blast radius.
Each stage is independently mergeable and leaves the test suite
green. Stages are **ordered by severity of confusion removed**, not
by LOC.

### Stage 1 — Capability matrix + naming clarity (docs-only, ~1 session)

**Scope**: Sphinx capability table + terminology glossary. No code.

**Deliverables**:

- New section `.. _theory-peierls-capabilities:` in
  `docs/theory/peierls_nystrom.rst`, near the top (before the
  existing executive summary). Contains:
  - One table with columns:
    `method | geometry kind | n_regions | n_groups | inner_radius | closure | k_eff | accuracy | verifies() label | test file`
  - One row per shipped reference case (~15 rows at current scope).
  - Rows are generated from `rv.continuous_all()` at Sphinx-build
    time (use `autodoc` or a hand-written script that imports the
    module).
- New `.. _theory-peierls-naming:` section with the terminology
  glossary:
  - "F.4" → use "Stamm'ler Eq. 34" for the formula; rename the
    Phase F.4 subphase to "Phase F (hollow cylinder)" in historical
    references.
  - `boundary=` string → `closure=` in new names, with explicit
    physical meaning.
  - `n_bc_modes` vs `n_surfaces` clarification.
- Update the page's top-of-file "Key Facts" block to reference the
  new capability matrix.

**Acceptance**:
- A fresh reader answers "what references exist for hollow cylinder
  F.4?" in < 30 seconds by opening the page.
- The auto-generated table matches `rv.continuous_all()` output
  bit-exactly.
- Sphinx builds clean.

**Risk**: low. Docs-only, no code churn.

**Budget**: 1 session (~3 hours).

### Stage 2 — Unified case-builder + auto-registry (~1 session)

**Scope**: Collapse the 3 per-geometry case builders into one, with
auto-discovery.

**Deliverables**:

- New `orpheus/derivations/peierls_cases.py`:

  ```python
  def build_peierls_case(
      kind: Literal["slab-polar", "cylinder-1d", "sphere-1d"],
      ng_key: str,
      n_regions: int,
      *,
      inner_radius: float = 0.0,
      closure: str = "vacuum",
      quadrature: QuadratureSpec | None = None,
      accuracy_budget: AccuracyBudget | None = None,
  ) -> ReferenceCase: ...

  def cases() -> list[ReferenceCase]:
      """Enumerate all shipped peierls reference cases."""
  ```

- `CurvilinearGeometry` gains two methods:
  - `.shell_volume_integral(r_nodes, r_wts, phi) -> float`
    (handles `2π ∫ r φ dr` vs `4π ∫ r² φ dr` internally)
  - `.reciprocity_factor(R_outer, r_inner) -> float`
    (R/r_0 for cylinder, (R/r_0)² for sphere)
- Auto-registry replaces `_continuous_modules`:

  ```python
  def _build_registry() -> dict:
      from orpheus.derivations import (
          diffusion, homogeneous, moc_mms, peierls_cases,
          sn, sn_mms,
      )
      modules = (diffusion, homogeneous, moc_mms, peierls_cases,
                 sn, sn_mms)
      refs = {}
      for module in modules:
          if hasattr(module, "cases"):
              for ref in module.cases():
                  refs[ref.name] = ref
      return refs
  ```

  Per-geometry module `continuous_cases()` functions delegate to
  `peierls_cases.cases(kind=…)`.

**Acceptance**:
- All existing `rv.continuous_get(name)` calls return the same
  `ReferenceCase` values (bit-exact k_eff to 1e-12).
- The 7 shipped integral-peierls references (1 slab + 3 cyl + 3
  sph) come from a single builder function.
- `peierls_cylinder.py` / `peierls_sphere.py` lose ~300 LOC each
  (duplication removed).
- New regression test `tests/derivations/test_peierls_cases_unified.py`
  verifies that the new builder reproduces the old builders bit-exactly
  on a representative sweep.

**Risk**: medium. The duplication includes subtle per-geometry
details (shell-volume normalization, reciprocity factor) that must
be correctly dispatched. Regression test per geometry mandatory.

**Budget**: 1 session (~4 hours). Gated on Stage 1's capability
matrix existing (it becomes the regression fixture).

### Stage 3 — Unified `ReferenceCase` (P4) — breaking change, 1-2 sessions

**Scope**: merge `VerificationCase` + `ContinuousReferenceSolution`
into one dataclass with an optional solution field.

**Deliverables**:

```python
@dataclass(frozen=True)
class ReferenceCase:
    name: str
    problem: ProblemSpec
    vv_level: str                     # "L0" / "L1" / "L2" / "L3"
    equation_labels: tuple[str, ...]
    description: str
    accuracy_budget: AccuracyBudget
    # Optional: if present, this case carries a reference SOLUTION
    # (continuous or discrete) computed by a derivation module.
    # If absent, this case is a problem description only (used by
    # solver tests that compute their own answer and compare against
    # k_inf analytically).
    solution: Solution | None = None

    @property
    def has_continuous_solution(self) -> bool:
        return self.solution is not None and self.solution.phi is not None

    @property
    def k_inf(self) -> float:
        return self.problem.k_inf()           # derived from materials
```

- `VerificationCase` and `ContinuousReferenceSolution` retained as
  deprecated type aliases / adapters for one release.
- `rv.get(name)` returns `ReferenceCase`; `rv.continuous_get(name)`
  deprecated as `rv.get(name).solution.phi` pattern.
- Migrations in all callers: `tests/`, `orpheus/`,
  `derivations/diagnostics/`.

**Acceptance**:
- Old and new APIs coexist for one release; deprecation warnings
  point at the new pattern.
- All existing tests pass with zero behavior change.
- Nexus graph shows the unified `ReferenceCase` replacing the two
  old classes.

**Risk**: high. Touches 171 grep hits across 19 files. Required
careful migration.

**Budget**: 1-2 sessions. Gated on Stage 2 (unified builder makes
the migration one-time per module, not one-time per case builder).

**Can be dropped**: if Stage 1+2 sufficiently de-confuse the
infrastructure, Stage 3 may not be worth the churn. Decide after
Stage 2 lands.

### Stage 4 — Unified `solve_peierls()` entry point (~1 session)

**Scope**: replace `solve_peierls_{cyl,sph,slab}_1g` + the
underlying `solve_peierls_1g(geometry, ...)` with a single
user-facing function.

**Deliverables**:

```python
def solve_peierls(
    kind: Literal["slab-polar", "cylinder-1d", "sphere-1d"],
    radii: np.ndarray,
    sig_t: np.ndarray,
    sig_s: np.ndarray,
    nu_sig_f: np.ndarray,
    *,
    inner_radius: float = 0.0,
    closure: str = "vacuum",
    quadrature: QuadratureSpec | None = None,
    max_iter: int = 300,
    tol: float = 1e-10,
) -> PeierlsSolution: ...
```

- `QuadratureSpec` dataclass encodes `n_panels_per_region`,
  `p_order`, `n_angular`, `n_rho`, `n_surf_quad`, `dps`. Presets
  `BASE`, `MED`, `RICH`, `ULTRA` exported from
  `peierls_geometry.quadratures`.
- `closure` alias table (backwards compat for one release):
  - `"white"` → `"white_rank1_mark"` (deprecation warning)
  - `"white_rank2"` → `"white_f4"` (deprecation warning)
- Per-geometry `solve_peierls_{cyl,sph,slab}_1g` retained as
  deprecated aliases.
- `n_beta` / `n_theta` parameter-name inconsistency resolved:
  both become `n_angular`.

**Acceptance**:
- `solve_peierls(kind="cylinder-1d", …, closure="white_f4",
  inner_radius=0.3)` produces bit-exact k_eff as today's
  `solve_peierls_cylinder_1g(…, boundary="white_rank2",
  inner_radius=0.3)`.
- Sphinx and test suite updated to the new name.
- All 34 `test_peierls_rank2_bc.py` tests still pass.

**Risk**: medium. Removes the per-geometry facade surface users
may depend on; back-compat aliases retained for one release.

**Budget**: 1 session. Gated on Stage 2.

### Stage 5 — Closure-dispatch consolidation (~1 session)

**Scope**: unify the three closure branches in
`build_closure_operator` and `solve_peierls_1g`.

**Deliverables**:

- Single dispatch table in `peierls_geometry.py`:

  ```python
  _CLOSURE_DISPATCH = {
      "vacuum":              _no_op_closure,
      "white_rank1_mark":    _mark_rank1_closure,
      "white_f4":            _f4_scalar_per_face_closure,  # hollow only
  }
  ```

- Explicit validation at entry: `closure="white_f4"` on solid
  geometry raises `ValueError` instead of silently collapsing to
  Mark. The silent collapse was confusing today — making it
  explicit removes one source of user surprise.
- Back-compat: a `strict=False` kwarg retains the silent-collapse
  behavior for one release, with a deprecation warning.
- Remove the `NotImplementedError` guard at
  `reflection="white", n_bc_modes > 1` (it guarded rank-N Marshak,
  which L21 has now definitively closed as a research path —
  retaining Marshak rank-N infrastructure is dead weight).
- Sphinx §`peierls-rank-n-per-face-closeout` updated to reflect the
  removal; Marshak rank-N primitives either fully archived or moved
  to `derivations/archive/`.

**Acceptance**:
- Single function per closure; 100% of today's `build_closure_operator`
  behavior preserved on valid inputs.
- `closure="white_f4"` + `inner_radius=0` raises a clear ValueError
  with the fix ("pass `inner_radius > 0` or use
  `closure='white_rank1_mark'`").
- Marshak rank-N primitives either deleted or archived with a
  clear Sphinx note. Sphinx page shrinks by ~500 lines (the §F.5
  Villarino-Stamm'ler subsection can move to agent-memory if kept
  at all).

**Risk**: medium. Archive vs delete decision on the rank-N Marshak
infrastructure requires one-session judgment call. Retaining it
"just in case" was already done; L21 says it is dead — committing
to the delete is the simplification.

**Budget**: 1 session, optionally 2 if Marshak archive work takes
an extra cleanup pass.

### Stage 6 — Accuracy budget + V&V audit integration (~1 session)

**Scope**: tolerance strings → `AccuracyBudget` dataclass; audit
tool consumes it.

**Deliverables**:
- `AccuracyBudget` dataclass in
  `orpheus/derivations/_types.py`.
- All `ReferenceCase.tolerance` strings migrated to structured
  `accuracy_budget` fields.
- `assert_rank_n_structural_win` extended to accept an
  `AccuracyBudget` and compare directly against closure CIs.
- `tests/_harness/audit.py` gains `--budget` mode that prints the
  accuracy budget per V&V label.

**Acceptance**:
- Machine-parsable accuracy numbers across the registry.
- No more `"O(h²) + scalar-mode residual ~1.4 %"` free-text.

**Risk**: low. Purely additive.

**Budget**: 1 session.

## 3. Ordering constraints

```
Stage 1 (capability matrix, docs)  — no code dep
     ↓
Stage 2 (unified case-builder)     — uses Stage 1 as regression fixture
     ↓                          ↘
Stage 3 (unified ReferenceCase)   Stage 4 (unified solve_peierls)
     ↓                          ↙
Stage 5 (closure dispatch)
     ↓
Stage 6 (accuracy budget)
```

Stages 3 and 4 are parallel after Stage 2. Stages 5 and 6 are
independent of 3 and 4 and can land at any time after Stage 2.

**Minimum-viable de-confusion** = Stages 1 + 2 + 5. Stages 3 + 4
are "nice-to-have" major rewrites that reduce LOC further but don't
meaningfully improve readability beyond what 1 + 2 + 5 already
deliver.

**Recommended**: execute 1 → 2 → 5 in three consecutive sessions.
Re-evaluate 3, 4, 6 after.

## 4. What's explicitly NOT in this plan

- **No changes to `peierls_geometry.CurvilinearGeometry` core
  machinery**. It is not the source of confusion — it is the well-
  designed unified core that the per-geometry facades wrap
  inconsistently. Touch it only where new methods (`shell_volume_integral`,
  `reciprocity_factor`) consolidate logic currently duplicated
  above.
- **No changes to the Peierls volume-kernel math**. The K_vol
  assembly is unified and correct. This plan is about the
  registration and dispatch layer, not the math.
- **No solid-cylinder rank-N DP_N investigation**. That is still
  Issue #103 and has its own research-vs-engineering question. This
  plan addresses the **infrastructure** ambiguity, not the open
  science.
- **No migration of the CP / MOC / SN / MC / diffusion modules**
  beyond what the unified registry requires. Those modules have
  their own conventions; harmonizing them is a separate plan.
- **No change to the V&V ladder (L0-L3, foundation)**. The
  `vv_level` field is retained as-is.

## 5. Open questions for user decision

**Q1 — Stage 3 vs not?**
Merging `VerificationCase` + `ContinuousReferenceSolution` is a
big migration across 19 files. The benefit is a single
`rv.get(name)` lookup; the cost is weeks of callsite updates.
Recommendation: do Stage 1+2+5 first, re-evaluate the reader-
confusion in practice, decide on Stage 3 with hindsight.

**Q2 — Delete or archive the Marshak rank-N infrastructure?**
L21 says this path is closed. Keeping it as dead code in
`peierls_geometry.py` costs ~500 lines of cognitive load on every
read. Options:
- (a) Delete entirely (L21 is decisive).
- (b) Move to `derivations/archive/peierls_marshak_rank_n.py` with
  a Sphinx note (reversible).
- (c) Retain behind the `NotImplementedError` guard as today
  (no change).
Recommendation: **(b)**. Committing to archive signals L21's finality
without losing the code. If the research revives (unlikely per
L21 but not impossible), restoring is one `git mv`.

**Q3 — Backwards-compat window.**
How long to keep `solve_peierls_cylinder_1g` as a deprecated alias?
One release (say, v0.X → v0.X+1 → removed at v0.X+2)? Forever?
This is an ORPHEUS-project policy question, not a technical one.
Recommendation: one release, documented in release notes.

**Q4 — Sphinx capability table: hand-written or auto-generated?**
Auto-generation from `rv.continuous_all()` is less likely to drift
out of date. Hand-writing is shorter to read. A hybrid works too:
hand-written narrative with an auto-generated table below.
Recommendation: hybrid, with the table auto-generated.

## 6. Acceptance for the overall simplification

After Stages 1 + 2 + 5 land, the original question
*"What references do we have for hollow sphere F.4 at r_0/R=0.3?"*
should be answerable by:

1. Open `docs/theory/peierls_nystrom.rst`.
2. Ctrl-F for "hollow sphere".
3. Read one row of the capability table.
4. Click the name to see the full `ReferenceCase`.

No code tracing. No explorer agent. No session-long confusion
chain. **If any of those four steps remain non-trivial after the
three stages land, the plan has failed and further stages are
warranted.**

## 7. Session trail (for audit)

- **2026-04-23**: this plan filed after cylinder + sphere F.4
  hollow-case registration landed and the user asked for a
  simplification roadmap. The confusion data points (C1-C5) are
  transcribed from the session's own false starts; they are
  empirical, not speculative.
- **2026-04-23 (later, same session)**: Stages 1 + 2 + 5 executed.
  See §8 below for landed deliverables and §9 for how the plan's
  open questions resolved in practice.

## 8. Landed deliverables (2026-04-23)

### Stage 1 — Capability matrix + naming glossary (landed)

**Files touched**: ``docs/theory/peierls_nystrom.rst`` only.

Added section :ref:`theory-peierls-capabilities` with:

- One ``list-table`` enumerating all 7 shipped integral-peierls
  continuous references (1 slab + 3 cylinder hollow F.4 + 3 sphere
  hollow F.4), each with geometry / n_g / n_reg / r_0/R / closure /
  accuracy class.
- One ``list-table`` capability matrix (closure × geometry) with
  ✅ / 🚫 marks, including 4 rows of explicitly-documented
  infrastructure gaps (multi-group hollow, multi-region hollow,
  solid rank-N DP_N, periodic/albedo/specular).

Added section :ref:`theory-peierls-naming` with:

- "F.4" disambiguation (Phase F subphase vs Stamm'ler Eq. 34
  closure formula).
- Boundary / closure string table: canonical names
  (``"white_rank1_mark"``, ``"white_f4"``) plus deprecated aliases
  (``"white"``, ``"white_rank2"``).
- Guard terminology (``NotImplementedError`` guards rank-N Marshak,
  not F.4).
- ``n_bc_modes`` vs ``n_surfaces`` distinction.
- Hollow vs solid semantics.

Also updated top-of-file Key Facts to cross-reference both new
sections, and fixed a pre-existing title-underline warning in §178
while in the file.

### Stage 2 — Helper methods + auto-discovery (landed)

**Files touched**:

- ``orpheus/derivations/peierls_geometry.py``: added two new
  ``CurvilinearGeometry`` methods
  (``shell_volume_integral``, ``reciprocity_factor``).
- ``orpheus/derivations/peierls_{cylinder,sphere}.py``: 4 sites
  (2 per geometry) rewritten to use
  ``GEOMETRY.shell_volume_integral(...)`` instead of hand-coded
  ``2π / 4π`` prefactors and ``r / r²`` radial weights.
- ``orpheus/derivations/reference_values.py``: manual
  ``_continuous_modules = [...]`` list replaced by
  ``pkgutil.iter_modules`` auto-discovery of any
  ``orpheus.derivations.*`` module with a callable
  ``continuous_cases`` attribute.

**Scope narrowed from the original plan**. The plan called for a
full unified ``peierls_cases.py`` with a single ``build_peierls_case``
entry point replacing the per-geometry builders. In execution, the
helper-method adoption + auto-discovery delivered the user-facing
value (de-duplication + eliminating the "module added but not
registered" failure mode — which this very session had hit twice)
with far less risk than the full builder unification. The fuller
unification remains available as Stage 2b for a future session.

**Regression**: all 9 hollow F.4 tests pass bit-exact pre- vs
post-migration (verified by the 9-passed sweep at commit boundary).
Auto-discovery finds the same 20 continuous references as the
previous manual list.

### Stage 5 — Closure dispatch + canonical names (landed, scope revised)

**Files touched**: ``orpheus/derivations/peierls_geometry.py`` (the
``solve_peierls_1g`` boundary-dispatch block) and
``docs/theory/peierls_nystrom.rst`` (the naming glossary table).

**Landed**:

- Canonical closure names ``"white_rank1_mark"`` and
  ``"white_f4"`` now accepted by
  :func:`~orpheus.derivations.peierls_geometry.solve_peierls_1g`
  and all thin wrappers.
- Deprecated aliases ``"white"`` (→ ``"white_rank1_mark"``) and
  ``"white_rank2"`` (→ ``"white_f4"``) continue to work but emit a
  ``DeprecationWarning`` on each call.
- ``closure="white_f4"`` on a 1-surface (solid) geometry now emits
  a ``DeprecationWarning`` explaining the silent collapse to rank-1
  Mark, with notice that this will become a ``ValueError`` in a
  future release. This is the C3 confusion surfaced.
- ``NotImplementedError`` message on ``closure="white_f4" +
  n_bc_modes > 1`` updated to cite L21 explicitly (was "planned in
  Phase F.5").

**Scope revised — no Marshak archive**. The plan's §5 Q2 called for
archiving the Marshak rank-N primitives
(``compute_P_esc_{outer,inner}_mode_marshak``,
``compute_G_bc_{outer,inner}_mode_marshak``, and
``_build_closure_operator_rank_n_white``). Execution revealed these
are **not dead code**: they are load-bearing for 7 test files
(``test_peierls_closure_operator.py``,
``test_peierls_rank_n_primitives.py``,
``test_peierls_rank_n_conservation.py``,
``test_peierls_rank_n_bc.py``, ``test_peierls_rank2_bc.py``,
``test_peierls_reference.py``) that gate conservation, reciprocity,
and the rank-1 routing through
``build_white_bc_correction_rank_n(n_bc_modes=1)``. Archive was
considered and **rejected** as infrastructure damage. The primitives
remain in ``peierls_geometry.py`` with the ``NotImplementedError``
guard preserved. Sphinx updated to document the retention rationale
rather than promise a future archive.

**Regression**: hollow F.4 tests pass 9/9 bit-exact. The plumbing
test (``test_solve_peierls_1g_hollow_sph_white_rank2_inner_radius_plumbing``)
now emits the expected ``DeprecationWarning``, which is caught and
reported by pytest as a warning (test assertions still pass).

## 9. Open-question resolutions (from §5)

**Q1 — Stage 3 vs not?** **Not done**. Stage 3 (merging
``VerificationCase`` + ``ContinuousReferenceSolution``) was
deferred. Stages 1 + 2 + 5 already significantly de-confused the
infrastructure; the Stage 3 churn across 19 files / 171 grep hits is
not worth it until a concrete user scenario demands the merged
type. Re-evaluate when adding the first CP-solver-specific
hollow-reference consumer.

**Q2 — Archive Marshak rank-N? Answer: no.** The plan's archive
recommendation assumed the primitives were dead code. In practice
they are used by 7 test files. Sphinx now documents the retention
rationale explicitly.

**Q3 — Backward-compat window**. Not formally specified; the
``DeprecationWarning`` fires with a "one release" wording in the
message. A follow-up issue should track the concrete removal plan
(e.g., "remove boundary alias support in v0.X+2"); this has not
been filed yet.

**Q4 — Sphinx capability table — hand-written or auto-generated?**
Hand-written for this release. The 7-row integral-peierls table is
small enough to maintain by hand, and the source-of-truth note
("this should match ``rv.continuous_all()`` filtered to
``operator_form == integral-peierls``") provides the drift-check
instruction. An auto-generator is a nice-to-have; not blocking.

## 10. Remaining work (optional, future sessions)

- **Stage 2b** — **SUPERSEDED** by the topology-first reframing in
  ``.claude/plans/topology-based-consolidation.md`` (2026-04-23).
  The original Stage 2b proposed a single ``build_peierls_case``
  dispatching on shape; the topology plan replaces that with TWO
  builders (``build_two_surface_case`` and
  ``build_one_surface_compact_case``) keyed by topology class. The
  topology framing is strictly sharper because it exposes the F.4
  applicability condition (2-surface only) at the signature layer
  instead of burying it in runtime dispatch. Execute the topology
  plan's stages T1 + T2 + T3 instead of the original Stage 2b.
- **Stage 3** — unified ``ReferenceCase`` dataclass. Big refactor
  (19 files, 171 grep hits). Defer until a concrete user need.
- **Stage 4** — unified ``solve_peierls()`` entry point replacing
  the 3 per-geometry wrappers. Budget: 1 session. Could be paired
  with Stage 2b since both affect the facade layer.
- **Stage 6** — ``AccuracyBudget`` dataclass replacing
  ``tolerance`` strings. Budget: 1 session. Independent.
- **Deprecation-removal follow-up issue** — track when the
  ``"white"`` / ``"white_rank2"`` aliases get removed. Also track
  when the ``white_f4`` + solid collapse becomes a ``ValueError``.

After Stages 1 + 2 + 5 landed, the original acceptance test from §6
is met: "What references do we have for hollow sphere F.4 at
r_0/R=0.3?" is answerable in four steps:

1. Open ``docs/theory/peierls_nystrom.rst``.
2. Jump to :ref:`theory-peierls-capabilities`.
3. Read the row ``peierls_sph1D_hollow_1eg_1rg_r0_30``.
4. Done.

No code tracing. No explorer agent. Answer visible in well under
30 seconds.
