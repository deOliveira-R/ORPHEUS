# Task #55 — `transport_sweep` docs-retirement pass + two phantom verifies-markers

Branch `refactor/sn-sweep-rename`. In-tree edits only (no commit; main agent reviews/commits).
`transport_sweep` retired at coupled-block campaign **step 6, R-6.1 (2026-07-12)** — the
successor is **context-dependent** (per-site behavioral judgment; a 1:1 rename is FORBIDDEN).

## Successor surfaces grounded in LIVE code (this session)

| Context | Retired `transport_sweep` role | Current surface (verified live) |
|---|---|---|
| Production within-group solve | the operator-free entry / `A_wg.solve` | `build_within_group_system(...).resolvent` + `.solve`; the block sweep is `InvertibleOperator.solve` running its own `loss_representation` |
| Geometry/quadrature dispatch | the `reduced is not None` / curvature branch | `default_for(mesh)` selects a `LossRepresentation` (`CumprodScan` 1-D scan / DAG wavefront multi-D) via each rep's `supports` |
| 1-D array-level sweep | `_sweep_1d_*` (dissolved) | `CumprodScan.sweep` → `_OneDimScanWalk.sweep` (SI twin of the matvec `_apply_walk`) |
| 2-D array-level sweep | 2-D wavefront body | `_sweep_jacobi` (LIVE, unchanged) |
| moment-windowing 1-D guard | `transport_sweep` *raises* on 1-D+moment | `CumprodScan.sweep` raises on a `moment_frame` (kwarg is `moment_frame`, **not** `moment_projection`) |
| q½ fold consumer | fed the fused joint channel | the within-group joint march (System-B q½ source into the resolvent grid) |
| test-side | — | `sweep_once` in `tests/sn/_test_helpers` (not referenced from docs) |

`Q_aniso` is **entirely gone** from `orpheus/sn/**` (the per-ordinate source is folded into `Q`).

## PART 1 — per-site disposition (56 sites; original inventory line numbers)

Dispositions: **(a)** behavioral rewrite to current API · **(b)** past-tense/history literal ·
**(c)** delete. Every `:func:`/`:meth:` DEAD xref was literalized (L-002: they render
plain-text with no `-W` warning; a dead code-xref is a Cardinal-Rule-1 staleness bug the build
never catches).

### boundary_conditions.rst (1)
| line | disp | what I wrote |
|---|---|---|
| 2190 | b | Phase-F SI/sweep-path history; `:func:`→literal "the then-production ``transport_sweep`` entry — since retired at step 6 (R-6.1)". |

### discrete_ordinates.rst (19)
| line | disp | what I wrote |
|---|---|---|
| 71 | a | Key Facts "1-D & 2-D sweeps BARE": 1-D `:func:transport_sweep` → `:meth:CumprodScan.sweep`; kept live `_sweep_jacobi` (2-D). |
| 257 | a | **"Quadrature Dispatch" whole-framing rewrite** → `default_for`/`LossRepresentation`/`CumprodScan`; deleted the stale `_sweep_1d_*` code block; kept physics (GL fast path / Lebedev / analytical-eigenvalue); cross-linked `:doc:loss_representations`. |
| 6343 | b(banner) | **"Unified sweep dispatch" section**: added L-013 arc-head supersession `.. note::` (step 6, R-6.1, → `default_for`/`CumprodScan`, `:doc:loss_representations`); literalized dead `:func:`; past-tensed intro ("consolidated"/"branched"). |
| 6350 | b(banner) | the `def transport_sweep(...)` code block now framed historical by the banner ("read as the then-current implementation"); left as illustrative Wave-D code. |
| 7395 | b | Wave-O bare-sweep history; `:func:`→literal "(then the production sweep, since retired at step 6)". |
| 9029 | b | pre-Phase-F divergence diagnostic; `:func:`→literal "(pre-Phase-F) dispatched through the then-production ``transport_sweep`` entry". |
| 9174 | b | pre-Phase-F smoking-gun diagnostic; `:func:`→literal "went through the then-production ``transport_sweep`` entry". |
| 9966 | a | ERR-058 mirror-seed **SI twin**; `:meth:transport_sweep` → `:meth:_OneDimScanWalk.sweep` (twin of `_apply_walk`). |
| 10455 | a | coupled-pole seed row (Retraction note @10460 states this row is **"unaffected"** = current production); `:meth:transport_sweep` → `:meth:_OneDimScanWalk.sweep`. |
| 11355 | a | q½ fold-factory consumers: "operator-free `transport_sweep`" → "the within-group joint march (System-B q½ source fed to the resolvent grid)". |
| 12714 | b | BC-param-removed history; `:func:`→literal "the then-production ``transport_sweep`` entry". |
| 13190 | a | scattering-source chain: "`Q_aniso` passed to `transport_sweep`" → "added to the isotropic source … consumed by the within-group sweep (the resolvent ``solve``)". |
| 13268 | a | Normalization Chain "**Sweep** (`transport_sweep`)" → "**Sweep** (the within-group resolvent ``solve``, `:meth:InvertibleOperator.solve`)". |
| 13535/13536 | a+c | `InvertibleOperator.solve` bullet: repointed to "the operator's own `loss_representation` sweep (… selected by `default_for`)"; **(c)** deleted the "Bit-identical to a direct `transport_sweep` call" clause (meaningless once the symbol is gone). |
| 13569 | a | "**sweep operator** (`transport_sweep`, dispatching…)" → "(`:meth:InvertibleOperator.solve`, dispatching its selected representation…)". |
| 13677 | a | "`solve` operates on structured arrays matching `transport_sweep`'s contract" → "matching the within-group sweep's array contract". |
| 15551 | b | BC-param-removed history; `:func:`→literal "the then-production ``transport_sweep`` entry". |
| 20711 | b | campaign step-6 changelog literal ("… `transport_sweep` wrapper. (6a …)"); left (already history). |

### index_convention.rst (22)
| line | disp | what I wrote |
|---|---|---|
| 409 | b | PR-INDEX-1 changelog literal "Public `transport_sweep` signature unchanged"; left (already literal). |
| 436 | b | PR-INDEX-3 changelog literal "2-D wavefront `transport_sweep` body"; left. |
| 465 | b | PR-INDEX-5 changelog; `:func:`→literal "``transport_sweep`` PUBLIC contract principled". |
| 782 | b | BoundaryField row; parenthetical `:func:`→literal "(the boundary arg of the since-retired ``transport_sweep``)". |
| 801 | a | "the sweep at `transport_sweep` consumes both" → "the within-group sweep (the resolvent ``solve``) consumes both". |
| 819 | b | ScalarSourceSink *Existing counterpart*; `:func:`→literal "``Q`` arg of the retired ``transport_sweep``". |
| 828 | b | AngularSourceSink *Existing counterpart*; `:func:`→literal "``Q_aniso`` arg of the retired ``transport_sweep``". |
| 847 | a | "the sweep at `transport_sweep` avoids a wasteful N-fold broadcast" → "the within-group sweep avoids …". |
| 1018 | a | "`A_wg` `solve` routes through the fused `transport_sweep`" → "through the fused sweep (… `InvertibleOperator` … runs the sweep on its selected representation)". |
| 1118/1120/1121/1123/1124/1126/1127/1129/1130/1132 | a | **shape table** rows: `transport_sweep input/return X` → `SN sweep input/return X`; every "Defined at" `:func:transport_sweep` → `:meth:CumprodScan.sweep`; the `Q_aniso` input row → "per-ordinate source (folded into ``Q``)" (separate arg retired). Shapes unchanged (they ARE the current convention). |
| 1204 | a | `.. todo:: Archivist` stub: "worked walk-through of the `transport_sweep` AngularBoundaryFlux contract" → "the within-group sweep's AngularBoundaryFlux contract". |
| 1388/1393 | a | **"Transport-sweep typed input" whole-framing rewrite**: retitled → "Typed sweep input"; deleted the stale `transport_sweep(…)` code example; described the current typed operator surface (`InvertibleOperator.solve` on typed source) + the bare `(Q, sig_t, boundary_flux)` array contract one layer down (`CumprodScan.sweep`/`_sweep_jacobi`); added a historical note (the `Q_aniso=` deprecation alias retired with the entry at step 6, R-6.1). No inbound `:ref:` to the old title (verified). |

### loss_representations.rst (6)
| line | disp | what I wrote |
|---|---|---|
| 1503 | b | "pre-carve code scattered across `transport_sweep`" literal; left (history). |
| 1896 | b | **fixed a broken dangling "The module-level" fragment** (a partial prior edit) + added R-6.1: "The operator-free `transport_sweep` entry RETIRED at step 6 (R-6.1) with its one production caller…". |
| 2470 | b | "the retired estate (step 6): … the operator-free module-level `transport_sweep` wrapper" literal; left (history). |
| 2712 | b | carve-arc history literal ("… spelled three different ways (`transport_sweep` branching …)"); left. |
| 2731 | b | S1-phase "`transport_sweep` rewired. Bit-identical." literal; left (history). |
| 2759 | b | rejected-alternative "An enum threaded into `transport_sweep`" literal; left (history). |

### operator_algebra.rst (8)
| line | disp | what I wrote |
|---|---|---|
| 4648 | b | "at Wave T it was the free function `transport_sweep`" literal; left (history). |
| 5292 | a | dispatch-predicate "`transport_sweep` uses to select the 1-D scan vs 2-D wavefront" → "the representation dispatch (`default_for`, via each representation's `supports`) reads to select the 1-D scan (`CumprodScan`) vs the 2-D wavefront". |
| 8089 | a | "`transport_sweep` *raises* on a 1-D mesh given a moment projection" → "the 1-D scan (`CumprodScan`) *raises* on a moment frame". |
| 8398 | a | "re-run a separate full-angular `transport_sweep`" → "a separate full-angular sweep — the within-group resolvent ``solve`` (`build_within_group_system`, applying its `.resolvent`)". |
| 8409 | a | "`transport_sweep` *raises* if a moment projection reaches a 1-D mesh" → "the 1-D scan (`CumprodScan`) *raises* if a moment frame reaches it". |
| 8458 | c | **deleted** `transport_sweep` from the moment-threading entry-point list (kept live `_sweep_scheduled`/`_sweep_jacobi`); in-scope also fixed stale `moment_projection`→`moment_frame` and `_sweep_2d_scheduled`→"the scheduled 2-D body". |
| 8459 | a | parenthetical "`transport_sweep` raises on a 1-D mesh" → "the 1-D scan ``CumprodScan.sweep`` raises on a moment frame". |
| 10571 | b | campaign step-6 "… operator-free `transport_sweep` wrapper are all gone" literal; left (history). |

**All 18 (b)-history survivors are double-backtick literals** framed as past/history; **NO
site presents `transport_sweep` as the current entry point**; **NO `:func:`/`:meth:`/`:class:`
dead xref to `transport_sweep` remains** in `docs/` (grep-verified).

## PART 2 — the two phantom verifies-markers

### `inverse-as-operator` → **MINT** the equation label
- **What the tests assert** (`tests/sn/operators/test_inverse_operator_equivalence.py`,
  module-level `pytestmark=[foundation, verifies("inverse-as-operator")]`, 4 tests — build
  emitted 4 skips, not the brief's "2"): the #226 Phase-2 keystone
  `(L+C).inverse().apply(b) == (L+C).solve(b)`, plus seed-drop / returned-surface-type checks.
- **Decision**: MINT, because there is **no existing label** naming the
  `A.inverse().apply ≡ A.solve` law (the inverse family's labels are `tensor-product-inverse`,
  `streaming-inverse-direct-sum`, `product-solve-reroute`, `matrix-inverse-*` — none state this
  law). The law IS stated in prose in `operator_algebra.rst` "Design B — `solve` pruned to
  native realizations" (taxonomy §11: *solving with an operator IS applying its inverse
  object*), but was **unlabeled** — so this is the "natural anchor exists unlabeled → mint"
  branch, per L-003.
- **Placed** the new `.. math:: :label: inverse-as-operator` (`A^{-1}b = A.inverse().apply(b)
  = A.solve(b)`) right where the law is stated, with `.. vv-status: inverse-as-operator
  documented` + a rationale comment naming the foundation bit-identity gate (L-004 — it is a
  structural software invariant, NOT a solver claim). Matches the sibling `product-solve-reroute`
  pattern verbatim. Confirmed it lands in the matrix's **"Documented-only equations"** bucket
  (excluded from the orphan gate), not flagged as an unverified solver claim.

### `ld-cartesian-1d` → **REPOINT** the marker strings to `ld-ubld-d1-reduction`
- **What the tests assert** (6 marks: `test_mms_ld_slab.py` × 5 + `test_ld_slope_frame.py` × 1):
  1-D slab LD MMS convergence, Krylov==SI, scan==DAG-oracle, thick-diffusive-limit (1G+2G),
  slope-frame consistency (ERR-061). All exercise the **1-D LD operator**.
- **Decision**: REPOINT (metadata-only, no test-behavior change), because `verifies` needs an
  *equation* label and the natural 1-D LD equation is **already labeled** `ld-ubld-d1-reduction`
  (the `d=1` Kronecker reduction = the 1-D LD 2×2 operator, `discrete_ordinates.rst:2742`). One
  equation = one label, so minting a duplicate `ld-cartesian-1d` on the same equation is
  impossible, and manufacturing a NEW equation just to host the label would be backwards. This
  is the exact L-003 precedent (the sibling `ld-slab` phantom was repointed the same way; this
  session I'd previously *flagged* `ld-cartesian-1d` as out-of-scope — Task #55 is its resolution).
  The asymmetry with the real `ld-cartesian-2d` label is a genuine structural fact: the 2-D LD
  MMS earned its own dedicated verification section (legA/legB, #247/#251/#257); the 1-D case IS
  just the base operator reduction.
- **Applied**: `ld-cartesian-1d` → `ld-ubld-d1-reduction` on all 6 marks; on `test_..._converges_second_order`
  the redundant `ld-cartesian-1d` was **dropped** (that mark already carried `ld-ubld-d1-reduction`).
  Every test keeps a meaningful verifies edge; `ld-cartesian-1d` is gone from the regenerated matrix.

## GATE OUTPUTS

1. **`.venv/bin/python -m sphinx -E -W -q docs docs/_build/html`** → **EXIT=0**;
   WARNING/ERROR/CRITICAL count **0** (unchanged from the pre-edit `-E` baseline of 0);
   **0** "no matching equation node" lines for `ld-cartesian-1d` / `inverse-as-operator`.
   Skip-line diff (baseline→post): baseline had exactly 10 (6 `ld-cartesian-1d` + 4
   `inverse-as-operator`) → post **0**; **nothing else changed** (no other verifies-target orphaned).
2. **`.venv/bin/python -m tests._harness.audit`** → **EXIT=0** (69/69 ERR coverage intact;
   verifies-linkage registry consistent).
3. **`git grep -n 'transport_sweep' -- docs`** → 28 survivors, **all double-backtick literals,
   all past-tensed/history-marked (disposition b)**; none present it as current; none is a live
   code-xref. (List above in the PART-1 tables under disp=b, plus the banner-framed 6350/6364
   code-block lines.)

## SIDE-FLAG LIST (do NOT fix — reported per brief)

The brief named `_unwrap_source`, `_require_leg_pair`, the kwarg-pair protocol surface as
other step-6-retired symbols to flag if presented as current. Findings:

- **`_unwrap_source`** — absent from `docs/` entirely. Nothing to flag.
- **`_require_leg_pair` / `_refuse_radial_characteristic` / `_require_radial_characteristic`**
  (`loss_representations.rst:2462-2463`) — appear ONLY in the "retired estate (step 6)"
  section, already framed as history ("the both-or-neither pair law died with the legs it
  guarded"). **Not presented as current** → no action needed.
- **kwarg-pair protocol surface** (`radial_characteristic_source`/`_flux`,
  `seed_cot`/`seed_cot_out`) at `loss_representations.rst:2458-2469` and
  `discrete_ordinates.rst:20855-20857` — both in retired-estate / Development-history sections,
  framed as retired. **Not presented as current** → no action needed.

**Adjacent staleness I encountered but left (out of `transport_sweep` scope; L-007 flag-don't-
rewrite):**

- **The `psi_bc` boundary-dict presented as a current contract** — `discrete_ordinates.rst:13700`
  ("Persistent boundary-flux dict ``psi_bc`` carrying state between sweeps") and `:18330`
  (`psi_bc["bc_1d"]`). The `psi_bc` dict was retired for `AngularBoundaryFlux` face views
  (`index_convention.rst:781` already says so). This is the section whose `solve` "structured
  arrays" contract I repointed at 13677. Needs its own `psi_bc`→`AngularBoundaryFlux` retirement
  doc-pass (separate verify-against-live behavioral rewrite).
- **`:file:orpheus/sn/sweep.py:474`/`:634`** (`discrete_ordinates.rst:9193`) — dead file refs to
  the dissolved `sweep.py`, but inside the pre-Phase-F smoking-gun diagnostic (accurate history
  of *where the bug was*); correct as history, left as-is.
- **`transport_operator_matvec_*`** (`discrete_ordinates.rst:13527-13529`) — the "historical
  BiCGSTAB FD operator" `:func:` refs; pre-existing, out of scope, not verified this pass.

## Files touched
- `docs/theory/{boundary_conditions,discrete_ordinates,index_convention,loss_representations,operator_algebra}.rst`
- `tests/sn/verification/mms/test_mms_ld_slab.py`, `tests/transport/spatial/test_ld_slope_frame.py`
- `docs/verification/matrix.rst` — **auto-regenerated** on the build (L-008: generated artifact,
  NOT hand-edited; reflects `ld-cartesian-1d`→`ld-ubld-d1-reduction` in coverage + `inverse-as-operator`
  now in "Documented-only").

## Quality self-assessment (Directive 3)
- Derivation depth 4 · Cross-references 5 (every new/repointed ref grep-gated against live code;
  every dead xref literalized) · Numerical evidence n/a (a terminology/retirement pass — no flux
  moves; structurally absent, not a deficit) · Failed approaches 5 (the campaign's own history
  preserved via banner+literals, not rewritten) · Code traceability 5 (successors verified in
  live `coupled_system.py`/`loss_representation/__init__.py`/`streaming.py`) · Derivation source n/a.
- Weakest: numerical evidence (structural — a symbol-retirement pass has no convergence table to show).
