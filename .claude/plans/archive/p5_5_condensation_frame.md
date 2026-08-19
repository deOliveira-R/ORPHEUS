# P5.5 — Sharpen the coarsening machinery (data-native, no frame subclass)

> **STATUS: MERGED to `main` @ `68ceb9a` (ff, 2026-06-27). #274 CLOSED. Feature branch deleted. This plan is
> now ARCHAEOLOGY** — the architecture lives in the docs (`discrete_ordinates.rst` condensation section), the
> lesson in `.claude/lessons.md` L32 (cross-layer assembler placement), the durable record in MEMORY.md.
> P5+P5.5 = main range `e9a6a49..68ceb9a`: `79d0092` (P5.5a Gram-typing) · `8105da2` (keystone P5.5b–e) ·
> `364f7b8` (P5.5f docs) · `68ceb9a` (3 nice-to-haves — `OverlapBasis.from_indicator` so `as_basis()` is the
> trial foundation / vector-morphism comment / `ce`→`clamped_edges`). Gates at merge: `-m "not slow"` **2747/0**,
> pyright `orpheus/` **412**, Sphinx **-W**, elegance **PASS**. **`origin/main` is 8 behind — NOT pushed** (user
> holds pushes). Next campaign phase: **P6 (#51, adjoint-weighted homogenization)**.
>
> **REVISED 2026-06-27** after a structural (cross-domain-attacker) + elegance (elegance-enforcer)
> review of P5 (condensation draft) **and** P3 (homogenization, merged). The user chose the
> **data-native, no-frame-subclass** shape — this OVERTURNS the prior "CondensationFrame(PetrovGalerkinFrame)
> in transport/frames/" keystone. Treat P5 as a draft; reshape both P3 and P5 toward the learned shape.
> Additive on top of the 4 P5 commits (NO rebase). Branch `feature/sn-energy-condensation`, unmerged.

## 0. Starting state (git-reconciled 2026-06-27)

Branch `feature/sn-energy-condensation`, 4 P5 commits ahead of main `42427e8`, 0 behind, unmerged:
`e9a6a49` (EnergyGrid+GroupCondensation+OverlapBasis) · `f750d40` (Mixture.condense) · `ffd908d`
(Solution.condense) · `71b5768` (docs). Tracked tree clean. P5 gates: condensation 70/70, full tree
`-m "not slow"` 2744/0, pyright 412, Sphinx -W.

## 1. What the review learned (the shape)

Both reviews CONVERGED. The frame core (`FrameBase`/PG/Galerkin discipline-as-type, `conjugate`,
`OverlapBasis` one-method override, `WeightedIndicatorBasis` raising on its unbuilt synthesis) is
**elegant — leave it**. The problems are at the EDGES:

1. **No `CondensationFrame`/`HomogenizationFrame` subclass.** A verb-carrying PG subclass is *false
   symmetry* (homogenization is intrinsically 2-frame: flux + production; condensation is 1-frame + a
   bare χ-sum) and *unjustified type-minting* (the frame is just a `PetrovGalerkinFrame` with a specific
   trial/measure/test — no new frame morphism; the HarmonicFrame analogy breaks because `.condense` is a
   composite over channels, not a typed face). The "condensation"/"homogenization" identity lives in the
   **overlap factory** + the **verb**, not a frame type.

2. **Condensation is `data + numerics`, NOT transport.** `EnergyGrid`, the overlap factory, `Mixture`,
   and the PG frame are all data/numerics — condensation has ZERO transport dependency (verified: data ↛
   transport; `Mixture.condense` imports only numerics). So `GroupCondensation` dissolves into
   **`EnergyGrid.overlap_to(coarse) → OverlapBasis`** (a binary factory that STAYS in data), and nothing
   moves to `transport/frames/`.

3. **The duplicated unit is the XS channel taxonomy** (`material_xs_field.py:398-409` ≈
   `mixture.py:288-308`) — Cardinal Rule 2. The genuinely-shared part is the **Mixture assembly** (wrap
   matrices in `csr`, thread `eg`): extract a `Mixture.from_dense_channels(...)` assembler that BOTH
   `Mixture.condense` (data) and `MaterialXSField.project_through` (transport→data) call. The verbs stay
   separate (different containers, different collapsed axis, different output shape = the asymmetry law);
   only the assembler is shared.

4. **Marginalize vs average** (the χ insight): the collapse has TWO morphisms, not "weight vs no weight".
   Σ (rate-bearing) → `frame.project` = G⁻¹M (preserves a reaction RATE). χ (a probability, Σχ=1) → the
   **bare table sum** `@ T` (preserves MASS, no G⁻¹) — provably NOT the weight-1 degenerate of `project`
   (that would divide χ by the group count and break Σχ=1). The matrix `[g',g]` = source `project`
   (average) + sink `@T` (marginalize). Document this asymmetry as deliberate in both verbs.

5. **EnergyGrid/Mesh duality.** `Mesh` already yields BOTH `volume_measure` + `indicator_basis()`;
   `EnergyGrid` must gain the symmetric pair `as_measure()` (source) + `as_basis()` (target, nested
   one-hot). The fractional overlap is the BINARY `(fine, coarse)` step — `overlap_to` — NOT a unary view
   (it reads both grids' edges). `coarse.as_basis()` (nested) ⊂ `fine.overlap_to(coarse)` (fractional).

6. **Gram-structure is a basis property, not a `FrameBase` hardcode.** `FrameBase.project`/`gram` use a
   row-sum probe valid ONLY for diagonal-support or partition-of-unity trials — its own 30-line docstring
   is a precondition the type can't enforce (a tapered / GEC-rank>0 basis would silently get a wrong
   Gram). Declare it on the basis (`GramStructure` enum: DIAGONAL / PARTITION_OF_UNITY / DENSE, default
   DENSE on the ABC); `project`/`gram` RAISE on DENSE. ~Type the seam; do NOT build the dense (MR)⁻¹M
   solve (no consumer — GEC rank>0 is #275).

**Deferred (pointers, not built):**
- `OverlapBasis` for *spatial* homogenization (non-matching meshes): conceptually the same general trial
  (nested `IndicatorBasis` = one-hot degenerate), but NO consumer yet + a real hazard (n-D box-
  intersection is a **2ᵈ-dense outer-product** table, not per-axis searchsorted — a fine cell straddling a
  coarse CORNER hits 4/8 coarse cells; MUST be 2-D-tested before adoption). Trigger consumers: CP region-
  averaging / MoC FSR-merging (the spatial transpose of energy-by-lethargy straddle). Leave a doc pointer.
- The dense-Gram `(MR)⁻¹M` solve (GEC rank>0 = #275).

## 2. The reshape (each phase: green + pyright ≤ 412 + additive commit)

### P5.5a — Gram-structure typed on the basis (numerics) — touches merged `frame.py`
- `numerics/basis/base.py`: `class GramStructure(Enum){DIAGONAL, PARTITION_OF_UNITY, DENSE}`;
  `Basis.gram_structure` property defaulting to `DENSE` (the safe default — a new basis must consciously
  declare it's row-sum-collapsible).
- Override: `IndicatorBasis` → DIAGONAL (disjoint support); `OverlapBasis` → PARTITION_OF_UNITY (rows sum
  to 1); `SphericalHarmonicBasis` → DIAGONAL (orthogonal). `WeightedIndicatorBasis` keeps DENSE (it is a
  TEST basis, never `frame.basis` — `project` reads the TRIAL's structure).
- `FrameBase.gram` (and thus `project`): if `self.basis.gram_structure is DENSE` → raise
  (`MissingCapability` if it fits, else a clear `NotImplementedError` naming the unbuilt dense seam +
  #275). The row-sum probe stays for DIAGONAL / PARTITION_OF_UNITY.
- Tests: a mutation gate — a `Basis` subclass returning DENSE makes `frame.project` raise; the existing
  DIAGONAL (homogenize) / PoU (condense) paths stay green. Update the `gram`/`project` docstrings (the
  precondition is now TYPED, not just prose).

### P5.5b — EnergyGrid dual views + `overlap_to`; retire `GroupCondensation` (data) — keystone
- `EnergyGrid.as_measure() → DiscreteMeasure` (counting: nodes `arange(n)`, weights 1, support
  `"energy"` — the current `GroupCondensation.measure` body, on `self`).
- `EnergyGrid.as_basis() → IndicatorBasis` (the group-index indicator: `IndicatorBasis((arange(n+1)-0.5,))`
  — the nested/target view).
- `EnergyGrid.overlap_to(coarse, within_group=None) → OverlapBasis` (the BINARY mismatch treatment): the
  `GroupCondensation.table` body + the span/upscaling/PoU guards (`__post_init__`); returns
  `OverlapBasis(edges_per_axis=coarse.as_basis().edges_per_axis, overlap_table=table)` — so `as_basis` is
  a real consumer.
- `OverlapBasis` gains the table diagnostics (axis-neutral, it owns the table): `fractional_columns`
  (was `locally_interpolated`) + `dominant_column` (was `coarse_of_fine`).
- **Retire `GroupCondensation`** (the class + `condense_to`/`from_grids`). 3-search audit done (§4).
  `WithinGroupSpectrum`/`InverseEnergySpectrum` STAY in `energy_grid.py` (data-level within-group model).

### P5.5c — shared assembler; rewrite both verbs (data + transport)
- `Mixture.from_dense_channels(SigC, SigL, SigF, SigP, SigT, SigS_dense, Sig2_dense, chi, eg) → Mixture`
  (the shared assembler: wraps matrices in `csr_matrix`, threads `eg`). The single home for "build a
  Mixture from collapsed dense channels."
- `Mixture.energy_grid` property (`EnergyGrid(self.eg)`; eg-None guard) — the XS carries its source grid.
- Rewrite `Mixture.condense(target, spectrum, within_group=None)` (STAYS data): `trial =
  self.energy_grid.overlap_to(target, within_group)`; `frame = PetrovGalerkinFrame(trial,
  self.energy_grid.as_measure(), WeightedIndicatorBasis(trial, spectrum))`; vectors → `frame.project`,
  matrices → sink `@ trial.overlap_table` + source `frame.project`, χ → `@ trial.overlap_table`; assemble
  via `from_dense_channels`. Collapse the `within_group is None` dance to a `default`.
- Rewrite `MaterialXSField.project_through` (STAYS transport) to assemble via `Mixture.from_dense_channels`
  (kills the twin). Behaviour bit-identical.

### P5.5d — Solution wiring (sn)
- `Solution.condense`: per material, build the spectrum (unchanged), call `material.condense(coarse,
  spectrum)` (the rewritten verb does the frame-building). Minimal change.
- `Solution.homogenize`: unchanged except it benefits from the shared assembler via `project_through`.

### P5.5e — energy-grid file cleanup (data) → `data/group_structures/`
- New `data/group_structures/ornl.py`: `ORNL_421 = EnergyGrid(...)` GENERATED from the library `eg`
  (extract programmatically — `load_isotope(...).eg`, already descending; isotope-independent). Emit the
  boundary literal (no hand-transcription, L4). Self-validating test: `ORNL_421 == load_isotope(...).eg`.
- `data/wims_d_iaea_group_structures.py` → `data/group_structures/wims.py` (rename + repoint). CONTENT-
  STRIP DEFERRED (user) — keep contents; a later pass decides what to cut.
- Delete `data/EGB421.txt` (superseded by `ornl.py`). Repoint `test_condensation.py`'s EGB421 parse →
  `from ...group_structures.ornl import ORNL_421`.
- `energy_grid.py` (the TYPE + within-group spectrum) STAYS at `data/energy_grid.py` (minimize churn;
  `material_xs_field`/`mixture`/`solution` import it).

### P5.5f — tests (migration) + docs
- Migrate the 3 P5 test files to the new API: `GroupCondensation.from_grids(...)` → `fine.overlap_to(...)`
  / `EnergyGrid` dual views; `Mixture.condense` signature change; EGB421 parse → `ORNL_421`. Preserve the
  structurally-independent hand-sum oracles + the Mode-11 sentinels (now on the rewritten verb). Add: the
  Gram-structure mutation gate (P5.5a), an `ornl.py` self-validation test, the `from_dense_channels`
  assembler bit-identity.
- Docs (archivist): rewrite the condensation section (`discrete_ordinates.rst`) — data-native machinery,
  `EnergyGrid` dual views + `overlap_to`, the marginalize-vs-average χ asymmetry, the Gram-structure
  typing; repoint `:class:`GroupCondensation`` refs; the spatial-OverlapBasis deferred pointer. Sphinx -W.

## 3. Commit structure (ADDITIVE — on top of the 4 P5 commits, NO rebase)
1. `refactor(numerics): declare Gram-structure on Basis; FrameBase.project refuses a dense Gram` (P5.5a)
2. `refactor(data): EnergyGrid dual views + overlap_to; retire GroupCondensation` (P5.5b)
3. `refactor(data): Mixture.from_dense_channels assembler + rewrite condense; dedup project_through` (P5.5c)
4. `refactor(sn): Solution.condense via EnergyGrid dual views` (P5.5d) — fold into 3 if trivially green
5. `refactor(data): group_structures/ — ornl.py generated grid + wims.py; drop EGB421.txt` (P5.5e)
6. `docs(sn): data-native coarsening machinery + V&V` (P5.5f, archivist)

Keep the tree green per commit (relocate logic + migrate its tests in the same commit). elegance-enforcer
review after the keystone (P5.5b+c).

## 4. Retirement audit — `GroupCondensation` (L20, three surfaces)
- **Production:** `mixture.py` `condense` (rewrite, NOT dissolve — stays data+numerics), `solution.py`
  `condense` (reroute via `material.condense`), `energy_grid.py` (retire class + `condense_to`/`from_grids`).
- **Tests:** `tests/data/test_energy_grid.py`, `tests/data/test_mixture_condense.py`,
  `tests/sn/test_condensation.py` — migrate to `overlap_to` / dual views (NOT delete).
- **Docs:** `discrete_ordinates.rst`, `api/data.rst`, + docstring cross-refs `measure.py:296`,
  `overlap_basis.py:9`.
- **Layering reconfirmed:** data ↛ transport (grep clean); `Mixture.condense` imports only numerics ⟹
  condensation STAYS in data (refutes the prior "dissolve for layering" premise).

## 5. Verification / discipline
- `.venv/bin/python -O -m pytest tests/numerics tests/sn tests/data -m "not slow" -q -rfE --timeout=300
  -p no:xdist -p no:cacheprovider` green; condensation gate green (bit-identical where math unchanged);
  pyright `orpheus/` ≤ 412; Sphinx -E -W. CLI `npx pyright` is the oracle, NOT streamed `<new-diagnostics>`.
- Surgical, main-agent-direct, user-steered; NO `method-implementer`. elegance review after the keystone.
- Commit/push ONLY when asked; stage explicitly (NO `git add -A`); trailer `Co-Authored-By: Claude Opus
  4.8 (1M context) <noreply@anthropic.com>`; no `# type: ignore` (`cast` OK); NEVER
  `git checkout/restore/stash` on uncommitted. Branch stays unmerged until the user directs (P5 + P5.5
  merge together).

## 6. Agent memos (the review record)
- cross-domain-attacker: structural verdict (one op `Collapse(axis, table, weight, normalize?)`;
  marginalize-vs-average; the 2ᵈ-dense n-D hazard; CP/MoC pollination). agentId `a612a966e43821619`.
- elegance-enforcer: `.claude/agent-memory/elegance-enforcer/frame_projection_coarsening_shape.md`
  (6 approval conditions; false-symmetry ruling; type-the-Gram-seam). agentId `af837d908769b492b`.
