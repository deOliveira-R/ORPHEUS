# P2-G — SN operators + mesh code-prose rebalance

Phase 2 of #231, batch P2-G. Branch `docs/sn-doc-architecture`.
Files: `orpheus/sn/operators/streaming.py`, `orpheus/sn/mesh/augmented_mesh.py`,
`orpheus/sn/operators/boundary.py`. Docstring/comment-only; ZERO code-token change.

## Headline finding (calibration vs the P2-A pilot)

Confirms the pilot's operator-file law: **the operator-algebra book is
exhaustively complete → ZERO MOVED.** Every teaching/derivation concept in these
three files has a TWIN home in the corpus (verified by grep of the landing
chapters). The rebalance is HISTORY-cut + TWIN-cut-to-pointer + trim, NOT
MOVED-writing.

**Key divergence from the P2-A pilot: these three files ARE `automodule`'d** in
`docs/api/discrete_ordinates.rst` (`streaming`/`augmented_mesh` with `:noindex:`,
`boundary` WITHOUT `:noindex:`), so their docstrings RENDER — the Sphinx build
gate is MANDATORY here (the pilot skipped it because `scattering.py` is not
autodoc'd). Good news: NO `.. math:: :label:` and NO `vv-status` in any of the
three docstrings, so cutting math blocks orphans no `:eq:` target and no
`verifies(...)` marker points at these docstrings. Baseline `-E -W` build = 0
warnings; acceptance gate = 0 post-edit.

The pre-inventory (`classify_P2-G.md`) HISTORY class was over-graded exactly as
the brief warned: many campaign-tagged docstrings are live keep-decisions,
invariants, or capability contracts. Re-adjudicated every one against the
"would a competent modifier who never leaves the file do the wrong thing
without this line?" test.

## Pointer convention

Literal greppable `docs/theory/<part>/<file>.rst §<label>`. `§` may point at an
eq-label or a section anchor; it is a human marker (files render plain-text under
`:noindex:` / by page convention), resolved by grep. Every cited label verified.
Surviving `:ref:`/`:class:`/`:meth:` roles kept as-is.

## Anchors used (all verified to resolve)

- `foundations/boundary_conditions.rst` — `:ref:`bc-extraction`` (§352) + the
  `bc-extraction-*` family (design-corrections, two-routes, reflect-trace, scope).
- `methods/sn/loss_representation.rst` — `§loss-rep-LpC` (L+C identity),
  `§loss-rep-resolution-a` (σ-free affine-in-σ), `§loss-rep-orientation-two-frames`
  (L^T adjoint two-factor / #280), `§loss-rep-inverse-adjoint-swap`,
  `§loss-rep-history` (symmetric-vs-WDD closure table).
- `methods/sn/curvilinear_numerics.rst` — `§sn-phase-d-err-026-closure-narrative`,
  `§sn-phase-d-carlson-coupled-pole-sweep`.
- `foundations/coupled_block_operator.rst` — `§coupled-block`.
- `foundations/structured_geometry.rst` — `§theory-structured-geometry`.

---

## File 1: streaming.py (adjudication)

| block | verdict | destination / record | content id |
|---|---|---|---|
| module docstring: "History" section (Wave E/D-H/D-I/D-J) | HISTORY (cut) | git log (`a614610`) | packed-vector-codec retirement narration |
| module docstring: kernel-location narration (#206/#238/S6.3/S6.9) | HISTORY-trim → CONTRACT | loss_representation.rst §loss-rep-fork-b2 | kept "kernel lives in loss_representation, ScanMarch default"; cut wave tags |
| module docstring: symmetric-closure note | CONTRACT (kept, trimmed) + TWIN | loss_representation.rst §loss-rep-history + curvilinear_numerics.rst §sn-phase-d-err-026-closure-narrative | kept the L.apply-symmetric-vs-sweep-WDD gotcha; cut the per-geometry derivation |
| module docstring: BC-handling note | CONTRACT (kept, trimmed) | :ref:`bc-extraction` | kept B-is-sibling / bare-loss_action / -B contract; cut the 1-D+multi-D dual re-derivation |
| `__all__` #206/#238 comment | COMMENT-cut | — | retired-split narration (now in module docstring) |
| `_require_typed_composite` docstring | CONTRACT (kept, trimmed) | — | the 2-guard contract; cut D-I.3d/#257/#240 step tags |
| `StreamingOperator` σ-free/Pattern-4 §§ | CONTRACT (kept, trimmed) | loss_representation.rst §loss-rep-resolution-a | wrong-simplification guard (no σ on L); cut probe/continuous teaching |
| `StreamingOperator` Notes (D-I.3d/D-J) | HISTORY (cut) | — | **STALE** ("accepts ONLY TimedFullField" — live contract is FullField) + codec-retirement |
| `StreamingOperator.block_role` comment | COMMENT (kept, trimmed) | — | kept FULL-operator + unannotated-dataclass constraint; cut apply_transpose-landed + `_eq_map`-retired |
| `is_adjointable` comment | CONTRACT (kept, trimmed) | loss_representation.rst §loss-rep-orientation-two-frames | kept two-factor + eager-raise + is_invertible; cut #280-axis teaching |
| `domain` docstring | CONTRACT (kept, trimmed) | loss_representation.rst §loss-rep-orientation-two-frames | kept composite-domain + R5 gotcha + predicate-gated; cut P4.5/#276/#280 campaign |
| `apply` docstring σ-recovery | CONTRACT (kept, trimmed) | — | cut the repeated affine-in-σ (now in class docstring) |
| `apply_transpose` summary | CONTRACT + **CORRECTNESS FIX** | — | "Hilbert transpose"→"Euclidean transpose" (body + sibling boundary.py agree; .H is Hilbert) — L-010 |
| `InvertibleOperator` (L+C)^-1≈WDD identity | CONTRACT (kept, trimmed) | operator_inverse_family / operator_algebra | cut R-1 date + taxonomy §-refs + campaign tags |
| `InvertibleOperator.apply` affine-coincidence | CONTRACT (kept, trimmed) | — | kept OWNS-the-matvec + single-source-σ; cut the multi-sentence derivation |
| `inverse()` collapse-trigger | CONTRACT (kept) | — | kept the defer-until-≥2 design guard; cut "solve retires at Phase 4" roadmap |
| `solve_moments` retired comment | COMMENT (kept, trimmed) | — | kept replacement (P@A.inverse) + guard; cut Phase-5c/#226/§17 tags |
| `_solve_timed_full_field` body comments | COMMENT (kept, trimmed) | — | kept bare-sweep + L21 constraints; cut D-H.2-C2/O.4a.2/S6.5 tags |
| `solve_transpose` wiring narration | CONTRACT (kept, trimmed) | — | kept swap-law + DD-only; cut "Phase 2.5c wires it" |

### boundary.py

| block | verdict | destination / record | content id |
|---|---|---|---|
| `_RayBoundaryFullFieldGain`-retired comment | HISTORY (cut) | git log | retired B.2b adapter |
| module docstring plan-file pointer | HISTORY-trim | :ref:`operator-algebra` + :ref:`bc-extraction` | dropped `.claude/plans/wave_o_...` path |
| `SNBoundaryOperator` O.2b-R5 parenthetical | HISTORY (cut) | — | "Before O.2b R5 B advertised the bare angular_trace" past-state |
| `_reflect_trace` adjoint-spelling gotcha | CONTRACT (kept, trimmed) + ⚠ | (module block matrix) | kept projection+transpose-spelling+the ⚠ vacuum-diagonal trap; cut legacy-sweep narration + B.2d/O.2 tags |
| `RadialCharacteristicBoundaryOperator` header | CONTRACT (kept, trimmed) | — | kept (B,B)-slot + RULING P1; cut "since B.2d / B.2b adapter retired" |
| `reflect_rows_inplace` | CONTRACT (kept, trimmed) + ⚠ | — | kept additive-forward-subst + ⚠ not-overwrite trap; cut the TWIN subspace aside (W2 gate/spec §13) |

### augmented_mesh.py

| block | verdict | destination / record | content id |
|---|---|---|---|
| TYPE_CHECKING B.5.A note | COMMENT (kept, trimmed) | — | kept "mesh provides shape only, no transport imports" rule; cut zeros_*-retirement |
| `bc` attribute doc | CONTRACT (kept, trimmed) | — | kept keyed-by-face + face-inventory-IS-BC-inventory; cut #188/C188.3 + pre-C4-retired |
| `BOUNDARY_OPERATOR_REGISTRY` comment | COMMENT (kept, trimmed) | — | kept law-classes + how-to-add-a-kind; cut pre-Wave-8-factory + #188 bypass |
| `__init__` legacy-inbound comment | COMMENT (kept, trimmed) | — | kept convert-once + self.mesh-provenance; cut C5.1/#225/retirement-path tags |
| `_init_core` method-agnostic-DATA comment | COMMENT (kept, trimmed) | — | kept _init_data-sets + materials-REQUIRED; cut "bit-identical/formerly inlined" |
| `_init_core` cell-update-strategy comment | COMMENT (kept, trimmed) | — | kept DD-default + bit-identical; cut Wave-D-#161 + "Wave C will ship" |
| `_init_core` angular-closure narrative (56 ln) | HISTORY (cut) + CONTRACT | curvilinear_numerics.rst §sn-phase-d-carlson-coupled-pole-sweep | kept default-closure + deferred-instantiation + CLASS-not-instance; cut Phase-C/D flip + PR-TYPED bullets + matvec-family-deleted |
| `_init_core` stencil-dispatch migration | COMMENT (kept, trimmed) | — | kept self.reduced-canonical; cut Wave-D/E 6-site migration narrative |
| `_init_core` boundary-trace historical-wart | COMMENT (kept, trimmed) | — | kept build-here-not-in-BC-resolve; cut C5.3/#290-wart/diffusion-twin |
| `_init_core` BC-resolution pre-C4 | COMMENT (kept, trimmed) | — | kept face-inventory-IS-BC + pole-none; cut pre-C4 bc_xmin-retired |
| backward-compat-accessors comment | COMMENT (kept, trimmed) | — | kept route-to-reduced + DeprecationWarning; cut 6-site migration narrative |
| retired-six-accessors comment | HISTORY (cut) | git log (#164) | alpha_half/redist_dAw/tau_mm retired in Wave E |

## Verdict counts

- **HISTORY (cut)**: 6 (streaming History-§ + Notes; boundary _RayBoundaryFullFieldGain; mesh angular-closure narrative + retired-6-accessors + [folded]).
- **COMMENT-cut / trimmed**: 14 (mostly campaign-tag trims keeping the constraint).
- **CONTRACT (kept, trimmed)**: 20 (every operator/mesh public + private-critical surface; embedded HISTORY/TWIN cut, latent-trap ⚠ imperatives kept inline).
- **TWIN (cut → pointer)**: folded into the CONTRACT-trims (symmetric-closure, σ-free affine, two-factor adjoint, angular-closure derivation → book pointers).
- **CORRECTNESS FIX**: 1 (streaming `apply_transpose` summary "Hilbert"→"Euclidean").
- **MOVED**: 0 (headline — the operator-algebra book is complete; grep-confirmed every teaching concept has a TWIN home).

## Before/after metrics (measured, `ast`+`tokenize`)

| file | total | docstrings | comments | code-token |
|---|---|---|---|---|
| streaming.py | 1161→1017 (−144) | 717→601 (−116) | 117→89 (−28) | IDENTICAL |
| augmented_mesh.py | 1605→1521 (−84) | 661→657 (−4) | 303→224 (−79) | IDENTICAL |
| boundary.py | 847→825 (−22) | 406→391 (−15) | 86→81 (−5) | IDENTICAL |

## Gates

- Code-token invariance (tokenize, drop COMMENT/STRING/NL/INDENT): **IDENTICAL** all 3 (1328/3034/1786).
- `import` all 3: **OK**.
- `pytest --collect-only -q`: **6652** (unchanged).
- BATCH SPECIAL `test_boundary_conditions.py`: **11 passed** (1 pre-existing warning).
- Sphinx `-E -W` build: **baseline 0 warnings → post-edit 0 warnings** (EXIT=0,
  "build succeeded"); the 3 files ARE `automodule`'d so this gate is load-bearing.
- All 5 literal `§<label>` pointers + both `:ref:` targets (`bc-extraction`,
  `operator-algebra`) verified to resolve in the theory tree.
- NO theory-page edits (`git diff --name-only docs/` empty). The other 10 changed
  `.py` files in the working tree are sibling P2 batches, not this batch.

## Judgment calls (calibrate P2-B…F operator siblings)

1. **These files ARE `automodule`'d (unlike the P2-A pilot's scattering.py) → the
   Sphinx build gate is MANDATORY.** No `:label:` in any docstring, so cutting math
   blocks is safe; but the docstrings render, so section-underline / malformed-markup
   risk is real — I kept section-title underlines verbatim (over-long is allowed) and
   only trimmed prose bodies.
2. **boundary.py is genuinely CONTRACT-dense** (a −3% cut). Over-cutting an operator's
   apply/transpose/reflect/split contract deletes load-bearing correctness content
   (Cardinal Rule 1). The HISTORY there is campaign-TAGS on live contracts, not
   standalone narration — trim the tag, keep the contract.
3. **augmented_mesh.py HISTORY lived in COMMENTS, not docstrings** (−79 comment / −4
   docstring). The docstrings are the mesh's public API surface (properties,
   classmethods) — contract. The `_init_core` comment cluster carried the wave-narrative.
4. **The pre-inventory's HISTORY class was over-graded exactly as warned.** Its
   COMMENT-cut grades for boundary.py `is_adjointable`/`domain`/`_reflect_trace`-schedule
   and mesh comments were WRONG — those state constraints (the intersection rule, the
   composite-carrier-not-bare-trace reasoning, the faces/rows subset semantics). KEPT.
5. **One CORRECTNESS FIX surfaced (L-010):** streaming `apply_transpose`'s summary said
   "Hilbert transpose" while its body + the sibling boundary.py both correctly call it
   the **Euclidean** transpose (`.H` is the metric Hilbert adjoint). Fixed to "Euclidean".
