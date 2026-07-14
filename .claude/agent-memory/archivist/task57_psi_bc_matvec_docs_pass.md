# Task #57 — `psi_bc` + `transport_operator_matvec` family docs-retirement pass

Branch `docs/sn-matvec-psi-bc-retirement`. In-tree edits only (no commit; no
pytest; `.claude/`/`scratch/` untouched). Sibling of the #55 `transport_sweep`
pass (`task55_transport_sweep_docs_pass.md`) — same per-site behavioral
disposition, NO mechanical rename.

## Ground truths RE-VERIFIED against LIVE code this session
- **`psi_bc: dict` RETIRED** (#197). Sole production mention = the narrating
  comment `orpheus/sn/solver.py:997`. Boundary state now on typed
  `AngularBoundaryFlux` face views.
- **`Q_aniso` kwarg GONE** from `orpheus/` (folded into the one per-ordinate
  source; `git grep Q_aniso -- orpheus` → 0).
- **`transport_operator_matvec` family DELETED — zero `def`s** (`git grep -nE
  'def (transport_operator_matvec|_transport_operator_matvec)' -- orpheus tests`
  → 0). The per-geometry trio + `_unified` successor all gone.
- **Live apply surface**: `StreamingOperator.apply` (streaming.py:407→452) →
  `LossRepresentation.streaming_action`; `InvertibleOperator.apply`
  (streaming.py:768→799) → `loss_action` → `_apply_walk` (#206 Phase C, L21
  "matvec ≡ sweep"). `apply_transpose` → `loss_action_transpose` (reverse walk).
  `solve` → `loss_representation.sweep`.
- **ARCHITECTURAL FACT the old docs missed**: apply and solve now run the ONE
  loss-representation walk (`loss_representation/__init__.py:1458` "matvec ≡
  sweep, ONE discretization, L21"; :933 "one walk is a code fact"). The former
  "two distinct discretisations (FD-operator apply vs WDD sweep)" design was
  DISSOLVED — this made three big doc sections stale at the THESIS level, not
  just the symbol level (see below).

## PER-SITE DISPOSITION TABLE
Dispositions: **(a)** behavioral rewrite to live surface · **(b)** past-tense
history + campaign citation + every dead `:func:`/`:meth:` literalized ·
**(c)** delete obsolete content. Line numbers = post-edit where I note "→Nxxx".

### boundary_conditions.rst (3)
| line | disp | change |
|---|---|---|
| 246 (`.. note::` #168 Phase C) | b | dead `:func:...transport_operator_matvec_spherical` → literal "the then-production ... matvec (since deleted, #197 / #280)"; added "live forward action is `:meth:StreamingOperator.apply` through the loss-representation walk". |
| 2079→2083 (`.. list-table::` Phase D call-order, "later retired by #282") | b | dead `:func:` → literal "the then-production ...matvec — the whole per-geometry family since deleted...". |
| 2187→2193 (Phase F intro) | **a** | PRESENT-TENSE stale "the matvec **lives at** `:func:..._unified`" → "then lived at ``_transport_operator_matvec_unified`` — since deleted at the walk unification (#280); the live forward action is now `:meth:StreamingOperator.apply` through the loss-representation walk". (`InvertibleOperator` `:class:` kept — live.) |

### structured_geometry.rst (1)
| line | disp | change |
|---|---|---|
| 450→454 (follow-on issues list) | **a** | PRESENT-TENSE "StreamingOperator/InvertibleOperator consume the primitive through `:func:..._unified`" → "consume the primitive (as ``SNMesh.reduced``) through the loss-representation walk: `.apply` reads ``self.mesh.reduced.coord`` inside the walk"; parenthetical past-tenses the deleted matvec consumption. Verified live: `mesh.reduced` = `ReducedStreamingOperator`, walk reads `.reduced.coord` (`loss_representation/__init__.py:3304,3886`). |

### discrete_ordinates.rst (the big one — the brief's 13543-13545/13576/13700-13704 all sat in a THESIS-stale region)
| line | disp | change |
|---|---|---|
| 6829 (Phase C Key Facts) | b | rewrote→rewrote; dead `:func:...spherical` → literal "the then-production ...matvecs (...since deleted — #197/#280)". |
| 7104→7106 (retired-symbols bullet list) | b | dead `:func:...spherical` → literal (list already literalizes retired `boundary_face_flux.*`). |
| 8319→8321 (Phase D "corrected injection-point story") | b | possessive `:func:...spherical`'s → literal "the then-production ...matvec's (since deleted)". |
| 8988→8990 (Phase F Key Facts) | b | dead `:func:...spherical` → literal (sibling `_sweep_1d_spherical` already "(dissolved `sweep.py`)"). |
| **13531–13561** "The full operator surface" three-bullet list | **a** | `apply` bullet (dead `transport_operator_matvec`+`_spherical`+`_cylindrical` roles, "reuses FD-matvec math") → REWRITTEN to `loss_action` walk (#206 Phase C) + historical parenthetical. `apply_transpose` bullet ("dense-matrix probing") → REWRITTEN to `loss_action_transpose` reverse walk + historical note (the dense path retired with `SNStreamingOperator` D-K — this bullet had even CONTRADICTED a later paragraph @13665 in the same section). `solve` bullet already correct — kept. |
| **13575** "Apply and solve use different closures by design" | b(banner) | THESIS now false (matvec ≡ sweep = one walk). Added `.. note:: Superseded architecture (Wave D / early Wave E)` banner: the split was dissolved by #206 Phase C; today one walk; FD family deleted. Retitled "...by design **(historical)**". Past-tensed "ships→shipped", "is→was". Literalized the `:func:...matvec_*` @13588 → "the ...family... since deleted". Kept the ERR-026 closure-bias narrative under the banner. |
| **13710** "Vector layouts (apply vs solve)" | **a**(section) | THE stale-as-current psi_bc/Q_aniso bullet-list site. Whole section premise (packed-1D-vector-vs-structured-array LAYOUT DIFFERENCE, `EquationMap`, BiCGSTAB) is dissolved — both verbs take `FullField` now. REWROTE to the unified `FullField` contract: source on `rhs.interior.values` `(N,ng,nx,ny)` (Q_aniso folded in), boundary on typed `AngularBoundaryFlux` face views (replacing `psi_bc` dict). Added `.. note:: Superseded packed-vector layout` recording the old EquationMap/psi_bc/Q_aniso triple + its #197/D-H–D-J retirement. **This edit is what eliminated "Persistent boundary-flux dict psi_bc carrying state" (grep for it → 0).** |
| 18108-18109→18144-18145 (numbered bug-fix log "Results After Fix") | b | dead `:func:...spherical`,`...cylindrical` pair → literals "the then-production ... — the matvec family since deleted, #197/#280". |

### index_convention.rst (the brief's 488/1507/1524 sat in TWO thesis-stale sections + a shape-table row)
| line | disp | change |
|---|---|---|
| 483 "What stayed deliberately legacy: the FD-matvec internal contract" | b(banner) | Title now FALSE (nothing stayed — deleted). Added `.. note:: Deleted path (#197 / D-H–D-J)`: the whole packed-vector path (matvec family + `solution_to_angular_flux*` + `EquationMap`) deleted → **PR-INDEX-7 is moot**. Retitled "...**(historical)**". Past-tensed the body ("carry→carried", "enumerates→enumerated"). Literalized dead `:func:transport_operator_matvec_cylindrical` AND the in-clause `:func:solution_to_angular_flux` (inseparable co-fix — a live dead-role in a sentence I'm rewriting). |
| 1500 "PR-INDEX-7 — EquationMap packed-vector traversal flip" | b(banner) | Whole subsection = detailed spec of the now-moot migration. **PRESERVED the `.. _sn-index-convention-future-work:` label** (referenced from :280, :458) + "Future work" parent header. Retitled "...**(obsoleted)**". Added `.. note:: Obsoleted by deletion`. Past-tensed ("carry→carried", "will flip→would have"). Literalized dead `:func:...matvec_*` + in-clause `:func:solution_to_angular_flux`. |
| 1155→1168 shape-table row "FD-matvec internal `fi` (deferred to PR-INDEX-7)" + its explanatory bullet | **c** | Deleted the obsolete shape-table row (`:func:solution_to_angular_flux` role) + collapsed the "Two arrays"→"One array" note (CollisionCache is now the only exception). Replaced the deleted `fi` row/bullet with a `.. note::` tombstone recording the former shape + its #197 deletion. |

### Production docstrings/comments (5, brief's remaining sites)
| file:line | disp | change |
|---|---|---|
| `orpheus/geometry/reduced_operator.py:135` (See-also docstring) | **a** | dead `:func:..._unified` → "through the loss-representation walk (`.apply` reads ``self.mesh.reduced.coord``)"; parenthetical past-tenses the deleted matvec consumption. |
| `orpheus/sn/mesh/augmented_mesh.py:310` (`#` comment) | **a** | dead `:func:..._unified` role → named live `:meth:StreamingOperator.apply` through the walk + "(...matvec family... deleted... #197/#280)". |
| `orpheus/sn/mesh/augmented_mesh.py:1221` (method docstring) | **a** | ``transport_operator_matvec_unified`` literal presented as the CURRENT consumer of `dag_walk_cell_indices` → named the live consumer (loss-representation walk, verified `loss_representation/__init__.py:2385`) + past-tensed the `_unified` predecessor. |
| `orpheus/sn/operators/streaming.py:57` (D-J retirement record) | LEAVE | already a past-tense ``literal`` ("had retired in D-H.2-C4e.6"). Verify-only. |
| `orpheus/sn/sweep/pole_angular_closure.py:22` (module docstring) | b | dead `:func:...spherical` role → literal "the pre-Phase-B ``transport_operator_matvec_spherical`` matvec (that whole ...family since deleted — #197/#280)". |

### psi_bc verify-only sites (all confirmed history-framed, LEFT)
- `discrete_ordinates.rst:18330→18367` (`psi_bc["bc_1d"]` in "Dead End #3" debug record) — reads unambiguously as history. LEAVE.
- `discrete_ordinates.rst:6364-6368` (Wave-D code block under #55 banner) — LEAVE (banner covers).
- `discrete_ordinates.rst:21689→21726` ("psi_bc retired" changelog) — LEAVE.
- `index_convention.rst:780→793` ("the pre-#197 psi_bc dict... since-retired", my #55 pass) — LEAVE.
- `index_convention.rst:1237→1252` ("Replaces the stringly-typed psi_bc: dict") — LEAVE.
- `orpheus/sn/solver.py:997` (the one production narrating comment) — LEAVE.

## GATE OUTPUTS
1. **`.venv/bin/python -m sphinx -E -W -q docs docs/_build/html`** → **EXIT=0**;
   WARNING/ERROR/CRITICAL set **empty**, unchanged from the pre-edit `-E`
   baseline (also empty). (One transient `Title underline too short` at
   `index_convention.rst:1518` from my "(obsoleted)" retitle — fixed: underline
   66→67 code points, L-009. Rebuild → 0.)
2. **`.venv/bin/python -m tests._harness.audit`** → **EXIT=0** (69/69 ERR
   coverage intact).
3. **Survivor greps**:
   - `git grep -nE ':(func|meth):`[^`]*transport_operator_matvec' -- docs orpheus`
     → **0** (zero live matvec roles anywhere; the 3 that existed in `orpheus/`
     docstrings were the 5 sites I fixed).
   - `git grep -n 'transport_operator_matvec' -- docs orpheus` → 25 survivors,
     **all double-backtick history literals in past-tense/deleted framing**.
   - `git grep -n 'psi_bc' -- docs orpheus` → 11 survivors, **all history-framed**;
     the "Persistent boundary-flux dict psi_bc carrying state" current-contract
     bullet is GONE.
   - All new `:ref:` labels (theory-sn-index-convention,
     bc-trace-contract-respected-by-matvec, sn-index-convention-future-work,
     sn-field-vocabulary) resolve. All new code-xrefs
     (`StreamingOperator.apply`, `LossRepresentation.{loss_action,
     loss_action_transpose,sweep}`, `FullField`, `AngularBoundaryFlux`,
     `AngularFlux`) grep-gated LIVE; the three `LossRepresentation` methods
     confirmed present on the Protocol via a python `hasattr` probe.

matrix.rst NOT regenerated-with-changes (no verifies-markers touched this pass,
unlike #55).

## SIDE-FLAGS (flag only, did NOT fix — out of the matvec/psi_bc scope; L-007)
1. **The `EquationMap` / `solution_to_angular_flux*` / `build_equation_map_*` /
   `pack_with_traces` codec cluster** presented with LIVE `:func:`/`:meth:`
   roles (dead → plain-text, no `-W` warning) at:
   `discrete_ordinates.rst:6525, 6596, 7073 (:meth:EquationMap.unknowns_at_cell_for_mask),
   7113, 12614 (build_equation_map_spherical), 13490`, and the D-J retirement
   record in `streaming.py:54-56` (`:class:EquationMap`, `:func:build_equation_map_*`,
   `:func:solution_to_angular_flux*`, `:func:pack_with_traces` — renders via
   automodule). All deleted at D-J. **Needs its own codec-family retirement
   doc-pass** (sibling of #55/#57). I co-literalized only the instances that
   sat INSIDE clauses I was already rewriting (index_convention deliberately-legacy
   / PR-INDEX-7 / fi-row); the standalone ones above I left.
2. **`index_convention.rst:456-458`** — a PR-INDEX-3 changelog `.. list-table::`
   row says "`EquationMap` packed-vector traversal **deferred** to PR-INDEX-7".
   Now internally inconsistent with my "PR-INDEX-7 obsoleted" banner. It's a
   `:ref:` (resolves) in a historical changelog row → left; belongs to the
   codec-cluster pass. (The `:280` `:ref:` to the same label is FINE — it points
   at a still-live parallel-prefix-M-M research item, not PR-INDEX-7.)
3. **`index_convention.rst:1233`** — `:class:~orpheus.sn.boundary_flux.AngularBoundaryFlux`
   uses the OLD path (should be `orpheus.transport.fields.angular_boundary_flux`).
   Dead code-xref (plain-text), pre-existing, unrelated to matvec/psi_bc.
4. **`discrete_ordinates.rst:9193`** (from #55 memo) — dead `:file:orpheus/sn/sweep.py`
   refs to the dissolved `sweep.py`, correct AS history (smoking-gun diagnostic).
   Left.

## Quality self-assessment (Directive 3)
- Derivation depth 4 (retirement pass — no new math; but I dug the LIVE
  apply/solve/walk relationship to write accurate rewrites)
- Cross-references 5 (every new ref grep-gated + Protocol-hasattr-probed; every
  dead role literalized; every `:ref:` label confirmed)
- Numerical evidence n/a (structural pass — no flux moves)
- Failed approaches 5 (the campaign history preserved via banners + literals,
  NOT rewritten — Wave-D two-discretisation reasoning + ERR-026 kept under
  banners; the moot PR-INDEX-7 spec kept as a record)
- Code traceability 5 (successors verified live in streaming.py /
  loss_representation/__init__.py / augmented_mesh.py)
- Derivation source n/a
- **Weakest: numerical evidence** (structural — symbol-retirement pass).

## The distinctive lesson (→ lessons.md L-020)
The brief's per-site list undersold the job: the three "sites" 13543-13545 /
13576 / 13700-13704 (discrete_ordinates) and 488 / 1507 / 1524 (index_convention)
were not isolated dead roles — they were the visible symptoms of FIVE
doc-SECTIONS whose entire THESIS had gone stale, because the retirement I was
sent to document (the matvec family) was one facet of a deeper architectural
unification (matvec ≡ sweep = ONE loss-representation walk, #206 Phase C) that
DISSOLVED the "two distinct discretisations / packed-vector-vs-structured-array
layout difference / deliberately-legacy-pending-PR-INDEX-7" design those sections
were built to explain. A retired-symbol docs-pass must READ THE ENCLOSING
SECTION'S THESIS against live architecture, not just the line: when the symbol's
deletion is a corollary of a design unification, the section's premise is stale
too. Disposition: SUPERSESSION BANNER at the section head (state the unification
+ the current one-walk truth + cite the campaign) + past-tense the body +
literalize dead roles, PRESERVING the historical reasoning under the banner
(Wave-D two-closure narrative + ERR-026 kept). For the ONE genuinely
stale-as-current contract (the psi_bc/Q_aniso "Vector layouts" bullet list), a
full behavioral section-REWRITE to the unified `FullField` contract, with the
old triple recorded in a trailing `.. note::`. Preserve `:ref:`-target labels on
moot future-work sections (retitle "(obsoleted)", don't delete). Co-literalize
only the deleted-sibling roles (`solution_to_angular_flux`) that sit INSIDE
clauses you're already rewriting; flag the standalone codec cluster for its own
pass.
