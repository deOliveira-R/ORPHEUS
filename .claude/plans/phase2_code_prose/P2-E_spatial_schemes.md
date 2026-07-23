# P2-E — `orpheus/transport/spatial/{scheme,diamond,linear_discontinuous}.py` code-prose rebalance

Phase 2 of #231, batch P2-E (spatial discretization schemes). Branch
`docs/sn-doc-architecture`. Pre-edit totals (ast/tokenize metric):
scheme.py **1355 ln (918 doc + 47 cmt)**, diamond.py **784 ln (390 + 151)**,
linear_discontinuous.py **810 ln (393 + 99)**.

## Headline finding — ZERO MOVED (the L-033 prediction held on a NON-operator file class)

The pre-inventory (`classify_P2-E.md`, Haiku) graded **13 MOVED**; the
re-adjudication finds **ZERO**. Every one of the 13 MOVED candidates is
already TWIN in the theory corpus — I grep-confirmed each landing chapter
(per L-033: "grep the landing chapter for the concept BEFORE assuming
novelty") before crediting it. The two concepts the pre-inventory thought
needed NEW theory pages both already have homes:

- **advection–reaction interface** → `discretization.rst §discretization-closures`
  carries the blend-weight identity `discretization-blend`, the
  exactness-on-constants argument, the outflow reconstruction
  `discretization-outflow-reconstruction`, AND a dedicated `.. admonition::
  The advection–reaction reading (why these are CFD scheme names)` — even
  cross-referencing `DiscretizationSchemeBase.outgoing_face_from_average`
  and `reaction_xs`.
- **reverse-mode / transpose kernel** → `loss_representation.rst
  §loss-rep-orientation-two-frames` (the discrete Euclidean reverse-DAG
  transpose) + `§loss-rep-scanmarch-facewise` (the facewise-vs-slopewise
  #240 D5-0 Mode-2 bug).

The tensor-product / Kronecker UBLD structure the pre-inventory wanted in a
new `_ubld.rst` page is `cartesian_multid.rst §spatial-moment-space`
(`spatial-moment-space-size` = `dim = (per_axis)^d`, DD/Step=1, LD=2; the
"append iff > 1" gate; the Kronecker order) + `§ld-ubld-multidim` +
`§ld-ubld-unified-moment-matvec` + `§ld-ubld-moment-scan`. The LD Schur
elimination is `discretization.rst §discretization-ld`
(`discretization-ld-system` / `-slope-elim` / `-schur`).

**This confirms L-033 generalizes past operator files to the spatial-scheme
file class: the book is exhaustively complete → budget TWIN-cutting +
CONTRACT-trimming, not MOVED-writing.** The reason is Cardinal Rule 3 — the
theory shipped with the code — and `foundations/discretization.rst` +
`methods/sn/cartesian_multid.rst` are two very large, complete foundational
chapters for exactly these schemes.

## Pointer convention

Literal greppable `docs/theory/<part>/<file>.rst §<label>` (ratified for
this phase; the pilot form). All 9 distinct pointers verified to resolve
(file exists + `.. _label:` OR `:label: label`). Existing role-refs
(`:doc:`, `:class:`, `:meth:`, `:func:`, `:attr:`) that survive a cut are
kept unchanged. NONE of the three files is `automodule`'d
(`docs/api/transport.rst` surfaces only `transport.method` +
`transport.reaction_rate_functional`), so ALL these pointers render
plain-text greppable literals — no Sphinx link breaks, no sphinx gate.

## Correctness fix surfaced during the trim (reported — beyond a pure contraction)

**scheme.py `DiscretizationScheme.is_affine_scannable` Protocol description
was STALE/WRONG.** It stated a Linear-Discontinuous closure "does **not**"
qualify as affine-scannable — but the LIVE LD code sets
`is_affine_scannable = True` (linear_discontinuous.py:232; its slope is
Schur-eliminated to a single-upstream recurrence; DD's docstring +
`discretization.rst §discretization-spectrum` agree). The Protocol text
predated Issue #158 Increment B. Per L-001 / Cardinal Rule 1 I FIXED it
while trimming ("`True` for Diamond Difference; `True` for
Linear-Discontinuous too — its two face moments couple, but the local slope
is eliminated by the Schur complement"). This is a Cardinal-Rule-1 doc/code
alignment, not scope-creep — but flagged here explicitly per the rubric.

Also fixed one stale INTRA-DOC reference: diamond.py `update()` said "See
module docstring § 'The unified body'" — I folded that section header into
the module intro, so the named section no longer exists → repointed to "the
module docstring" generally (plain-text ref, no `-W` warning; L-002
staleness caught by the `§`-grep gate).

---

## Adjudication table — File 1: `scheme.py`

| block | verdict | destination / record | note |
|---|---|---|---|
| module: "Why this abstraction exists" (Round-1-Wave-C / #157/#158/Wave-D) | HISTORY + TWIN → cut | `discretization.rst §discretization-space-angle`; `sn/index.rst` | kept terse strategy rationale + geometry-blind pointer; cut campaign narration |
| module: "What each dataclass holds" | TWIN (intra-file dup) → cut | the UpstreamState/CellResult class docstrings (SSOT) | cut the module-level restatement of field docs that live directly below |
| module: "Geometry-as-data (Step 2.5)" | TWIN → cut to pointer | `discretization.rst §discretization-space-angle` | kept the structural-presence-not-dispatch one-liner (constraint); cut the slab/sphere/degenerate derivation |
| module: "Where downstream consumers will call this" (Wave-D dag_walk code example) | HISTORY → cut | — | cut the prospective code example; kept "consumed via SNMesh.dag_walk per-visit packets" |
| module: References / See also | CONTRACT | — | KEEP (canonical citations + navigation) |
| CellVisit intro (MoC-not-a-consumer) | CONTRACT-trimmed | — | kept neutral-defaults + method-generic + no-shared-SweepGraph constraint; cut the fiber-bundle/sheaf teaching |
| CellVisit `c_in`/`c_out`/`tau` attrs (#236 provenance) | CONTRACT-trimmed | — | kept formulas (sign/index) + neutral defaults + "consume as data, don't rebuild" + the ⚠ tau≠0.0 divide-by-zero trap; cut #236 accessor-path narration |
| CellVisit `face_area_downstream`, `cell_idx`, `streaming_terms` | CONTRACT | — | KEEP (data contracts + degenerate 0.0 signal) |
| UpstreamState / CellResult class docstrings | CONTRACT | — | KEEP (call-site input/output shapes) |
| Protocol 7-trait `Attributes` block | CONTRACT-trimmed | `loss_representation.rst §loss-rep-scanmarch-facewise`; `cartesian_multid.rst §spatial-moment-space` | kept trait semantics + DD/LD value + consumer + ⚠ traps; cut derivations; **FIXED the LD affine-scannable staleness** |
| Protocol `update` / `residual` methods | CONTRACT | — | KEEP (params/returns/round-trip invariant/matvec use-case) |
| Base class docstring: Σ-stateless + reconstruction-ops + bit-identity | CONTRACT-trimmed + TWIN | `discretization.rst §discretization-closures` | kept ABC/registry contract + Σ-stateless "no XS state" + generic-reconstruction single-source + scan-strict-vs-matvec-~1ULP note; cut the CFD teaching + exactness-on-constants derivation |
| 6 ClassVar default-docstrings (is_affine_scannable … has_transpose_kernel) | CONTRACT-trimmed | `cartesian_multid.rst §spatial-moment-space`; `loss_representation.rst §loss-rep-scanmarch-facewise` / `§loss-rep-orientation-two-frames` | kept default + opt-in + per-scheme value + ⚠ + consumer + pointer; cut the repeated derivations |
| `is_multi_moment` property | CONTRACT-trimmed | — | kept definition + SSOT role + equivalences; cut the detailed teaching |
| base `cell_kernel_batch` / `residual_kernel_batch` / `affine_scan_coefficients` / `cartesian_scan_coefficients` / `reflect_scan_coefficients` / `moment_scan_closure` (raising defaults) | CONTRACT | — | KEEP (override contracts); lightly trimmed the caller-composition narration on cartesian_scan_coefficients |
| `source_emission` / `cell_average` / `outgoing_face_from_average` staticmethods | CONTRACT-trimmed | `discretization.rst §discretization-closures` | kept formula + DD/LD special case (point-of-use); cut the `.. note::` bit-identity essay → class-docstring pointer |

## Adjudication table — File 2: `diamond.py`

| block | verdict | destination / record | note |
|---|---|---|---|
| module: "Architectural shift (Step 2.5)" | HISTORY → cut | git history | cut the cumprod-retirement / re-baseline / plan-pointer migration narration |
| module: "The unified body" (3 observations) | TWIN → cut to pointer | `discretization.rst §discretization-space-angle` | kept the 3 observations terse + closure formulas + structural-presence constraint; cut the derivation |
| module: References / See also | CONTRACT | — | KEEP |
| `_DD_W` constant comment | CONTRACT | — | KEEP (single-source-of-truth constraint) |
| class docstring + is_linear / is_positivity_preserving | CONTRACT | — | KEEP (terse; the negative-flux formula is point-of-use) |
| is_affine_scannable / transverse_coupling_is_facewise / diffusion_limit_consistent / has_transpose_kernel ClassVars | CONTRACT-trimmed | base traits + `loss_representation.rst §loss-rep-scanmarch-facewise` | kept DD value + one-line rationale + the ⚠ spatial-vs-angular-β trap; cut the derivations + campaign tags |
| supports_curvilinear ClassVar | CONTRACT | — | KEEP (terse selection guard) |
| `update()` inline comments (#236 B3 c_in/c_out, angular τ) | COMMENT-trimmed | — | kept "consume as data, don't rebuild from st.alpha_*/tau_mm" + the M-M formula; cut #236 provenance |
| `update()` spatial-closure comment | CONTRACT | — | KEEP (degenerate-None + w=½ byte-identical reconstruction) |
| `residual()` docstring + PR-TYPED-6.5/#236-B2 comment | CONTRACT-trimmed | — | kept round-trip + delegation + the "don't rebuild" constraint; cut the campaign narration |
| shared-helpers comment (#240 D5a) | COMMENT-trimmed | — | kept the single-source constraint; cut "before #240 D5a in THREE producers" |
| `_cartesian_streaming_diagonal` docstring | CONTRACT | — | KEEP — the CANONICAL operation-order bit-identity discipline (the ⚠ do-NOT-regroup) other methods point to |
| `_reflection_coeffs` docstring | CONTRACT | — | KEEP (formula + byte-identical) |
| `cell_kernel_batch` / `residual_kernel_batch` docstrings | CONTRACT-trimmed | — | kept single-source + axis-convention + math + "diamond 2 to BOTH denom+numer" bug-guard; cut the REPEATED bit-identity explanation → `_cartesian_streaming_diagonal` pointer |
| `affine_scan_coefficients` docstring | CONTRACT-trimmed | — | kept the denom/a math + slab reduction + "what NOT computed" + ordering + the ⚠ bit-identity TRAP + the "Why a SEPARATE diagonal (#242)" don't-merge constraint (condensed); cut the "before #236 inlined" HISTORY + the duplicated bullet math |
| `cartesian_scan_coefficients` docstring | CONTRACT-trimmed | — | kept the coefficient math + ×V convention; lightly trimmed the op-order note |
| `reflect_scan_coefficients` docstring | CONTRACT | — | KEEP |
| S6.4(e) storage-adapters-retired comment | HISTORY-trimmed | — | kept the THREE-capability-groups structural map (CONTRACT for extenders); cut the storage-retirement narration + retired method names |

## Adjudication table — File 3: `linear_discontinuous.py`

| block | verdict | destination / record | note |
|---|---|---|---|
| module intro + diffusion-limit `.. note::` (Increment C / ERR-061 / #247) | CONTRACT-trimmed + HISTORY | `cartesian_multid.rst §ld-ubld-sweep-global-frame` | kept the current-status contract (scattering slope threaded → recovers limit, pinned by test) + the ⚠ #247 Mode-10 gap; cut the ERR-061 frame-fix mechanism → pointer |
| module: "Why LD carries two moments" (moment defs + 2×2 matrix) | TWIN → cut to pointer | `discretization.rst §discretization-ld` | kept the ⚠ slope-row sign trap (LM-1989 memo §1.4/§6) + the SymPy-exact-on-linear validation; cut the moment integrals + the explicit 2×2 |
| module: "Schur-complement scalar contract" | TWIN → cut to pointer | `discretization.rst §discretization-ld` (`discretization-ld-schur`) | kept "Schur-eliminated → scalar, fits the existing residual contract unchanged; slope source on source[1]; the source-iteration threading is deferred #158"; cut the explicit Schur formula |
| module: "Scope (slab/Cartesian only)" | CONTRACT | — | KEEP (the curvilinear-unpublished raise-guard) |
| module: "Traits" section | TWIN (intra-file dup) → cut | the ClassVar docstrings (SSOT) | cut the module-level restatement of the trait docstrings below |
| module: References / See also | CONTRACT | — | KEEP |
| `_require_slab` / `_ld_source_moments` / `_LDCellTerms` (+ `slope()` sign site) | CONTRACT | — | KEEP (guards + the (2,ng) split + the ONE slope-row sign-convention site) |
| is_affine_scannable / spatial_basis_per_axis / diffusion_limit_consistent / has_transpose_kernel ClassVars | CONTRACT-trimmed | `cartesian_multid.rst §spatial-moment-space` | kept LD value + the ⚠ traps (flat-source-not-consistent; bilinear→DAG; #280 deferral); cut the tensor-Legendre slot detail + campaign tags |
| supports_curvilinear / theta ClassVars | CONTRACT | — | KEEP (theta = 1/3 + the ⚠ do-NOT-change-to-LLD trap) |
| `_schur_terms` (docstring + delegation comment) | CONTRACT-trimmed | — | kept the single-algebra-site + d1_closed_form delegation + ÷V/×V convention; cut "in the SAME shape ... proven ==" narration |
| `update()` / `residual()` docstrings | CONTRACT | — | KEEP |
| DAG-family batched-kernel comment (the big block) | CONTRACT-trimmed + HISTORY | `cartesian_multid.rst §ld-ubld-unified-moment-matvec` | kept the group-2 + ÷V convention + scale-free + L16-fast-path CONTRACT; cut the "former d=1 scalar arm RETIRED / flat-source artifact" HISTORY → pointer |
| `_ubld_system` docstring | CONTRACT | — | KEEP (assembly + the moment-vs-flat rank discriminator + AVERAGE_MOMENT slot) |
| `_ubld_inflow` docstring | CONTRACT (pre-inv MOVED **overturned**) | — | KEEP — the d=1 trailing-axis-append gotcha + "keys on face_moment_tail, NOT a shape probe" is exactly the "would a modifier do the wrong thing?" CONTRACT test |
| `_ubld_outgoing_faces` docstring | CONTRACT (pre-inv MOVED **overturned**) | — | KEEP — the trace-order + "out-of-cell == in-of-next-cell" consistency is the reconstruction contract |
| `cell_kernel_batch` docstring | CONTRACT | — | KEEP; trimmed a campaign tag |
| `residual_kernel_batch` docstring | CONTRACT-trimmed | — | kept the mass-diagonal-normalization (load-bearing: why ÷M for row-for-row S.apply consistency); cut the "former scalar arm retired / L21" narration |
| scan-family coefficients comment + `affine_scan_coefficients` | CONTRACT | — | KEEP (the ×V convention + the math + the slab-only guard); trimmed the delegation-comment campaign tag |
| `moment_scan_closure` docstring | CONTRACT (pre-inv MOVED **overturned**) | `cartesian_multid.rst §ld-ubld-moment-scan` | KEEP the method contract + slope-fold + slab-only; added the theory pointer for the slope-augmented-source derivation; cut the campaign tag |

---

## Verdict counts (re-adjudicated; blocks, primary verdict)

- **MOVED: 0** (headline — all 13 pre-inventory MOVED overturned to TWIN or CONTRACT).
- **TWIN (cut → `§`-pointer):** ~8 primary blocks (module geometry-as-data ×2, unified body, why-two-moments, Schur contract, the reconstruction-ops/advection-reading, the has_transpose/facewise derivations) + TWIN-cuts folded into ~20 CONTRACT-trimmed rows.
- **HISTORY (cut; record/git owns it):** ~6 (Why-this-abstraction, Where-consumers-call, Step-2.5 shift, storage-adapters-retired, the retired-scalar-matvec-arm, the ERR-061 mechanism).
- **CONTRACT (kept; many trimmed of embedded TWIN/HISTORY):** the dominant tier — every trait spec, method contract, capability-override contract, the ⚠ latent traps, the sign/index conventions, the mass-normalization, the d=1 gotchas.
- **COMMENT:** constraint-stating kept (the `_DD_W` single-source, the `_cartesian_streaming_diagonal` do-NOT-regroup, the three-capability-groups map); narration trimmed.
- **1 correctness FIX** (LD affine-scannable Protocol staleness) + **1 stale intra-doc ref FIX** (diamond `update()`).

## Before/after metrics (ast/tokenize)

| file | total | docstring | comment |
|---|---|---|---|
| scheme.py | 1355 → 1048 (−307) | 918 → 690 (−228, −25%) | 47 → 47 (0) |
| diamond.py | 784 → 680 (−104) | 390 → 331 (−59, −15%) | 151 → 123 (−28, −19%) |
| linear_discontinuous.py | 810 → 726 (−84) | 393 → 328 (−65, −17%) | 99 → 89 (−10, −10%) |
| **aggregate** | **2949 → 2454 (−495)** | **1701 → 1349 (−352, −21%)** | **297 → 259 (−38, −13%)** |

(scheme.py's comment count is 0-delta: its comments are all constraint-stating
— the Protocol scan-family NOTE, the scan-family capability comment, the
reconstruction-ops comment — CONTRACT, correctly kept.)

## Gates

- **GATE 3 (token invariance, the primary proof of zero code change):**
  non-string/non-comment token count AND sha256 sequence hash IDENTICAL to
  HEAD for all three — scheme `951 / ea8bab918c15823b`, diamond `1213 /
  dc03fa86f970d2d2`, LD `1607 / 7bb56b64d0d7733a`.
- **GATE 1 (import):** all three import OK.
- **GATE 2 (collect-only):** `6652 tests collected` — UNCHANGED.
- **GATE 4 (pointers resolve):** all 9 distinct `§label` pointers verified
  (file-ok + label-ok) + `discretization-ld-schur` eq-label + the kept
  `:doc:/theory/methods/sn/index` + the intra-doc heading "Generic affine
  reconstruction ops". No dangling `§` (all remaining `§` are literature /
  skill / issue section refs).
- **GATE 5 (no theory-page edits):** none — ZERO MOVED, so no book edits.
- **Sphinx build: NOT run, and none needed** — none of the three files is
  `automodule`'d anywhere (`docs/api/transport.rst` surfaces only
  `transport.method` + `transport.reaction_rate_functional`), so their
  docstrings are invisible to the build; the `-E -W` baseline (1 warning) is
  definitionally unchanged (L-033 point 2).

## The hardest judgment calls (calibrate the six sibling batches)

1. **ZERO MOVED held on a NON-operator file class.** The pilot (P2-A) was an
   operator file; P2-E is spatial schemes with different landing chapters
   (`discretization.rst`, `cartesian_multid.rst`, `loss_representation.rst`).
   The pre-inventory graded 13 MOVED with confidence; ALL were TWIN/CONTRACT
   after grepping. The `discretization.rst §discretization-closures` +
   `cartesian_multid.rst §spatial-moment-space` sections are so complete they
   even cross-reference the exact code symbols. **Discriminator held: grep the
   landing chapter BEFORE crediting a MOVED — Cardinal Rule 3 means the theory
   shipped with the code.**

2. **The `_ubld_inflow` / `_ubld_outgoing_faces` / `moment_scan_closure`
   "MOVED" overturn is the sharpest calibration.** The pre-inventory read the
   d-generic Kronecker reconstruction as book-teaching. But the CONTRACT test
   ("would a modifier who never leaves the file do the wrong thing without this
   line?") flips it: the d=1 trailing-axis-append, the "keys on
   `face_moment_tail`, NOT a shape probe" reasoning, and the trace-order /
   "out-of-cell == in-of-next-cell" consistency are all load-bearing gotchas a
   modifier MUST see at the method. The Kronecker *layout theory* is TWIN
   (pointed to `§spatial-moment-space`); the *reconstruction contract* is
   CONTRACT (kept). Same object, two layers — keep the contract, point the
   theory.

3. **A trim surfaced a Cardinal-Rule-1 stale claim (LD affine-scannable).**
   Verifying the Protocol trait against the LIVE LD code (L-001) caught a
   "LD does not qualify" that Issue #158 Increment B made false. A pure
   contraction would have transcribed it. Reported + fixed, not smuggled.

4. **The bit-identity discipline is a single-source (Cardinal Rule 2) even
   across docstrings.** `_cartesian_streaming_diagonal` carries the canonical
   "explicit left fold, do NOT regroup" statement; `cell_kernel_batch` /
   `residual_kernel_batch` / `cartesian_scan_coefficients` REPEATED it. Trimmed
   the repeats to a pointer at the canonical helper — a within-file dedup, not
   a TWIN-cut. The ⚠ do-NOT-regroup imperative stays at the single source.

5. **Latent-trap imperatives + open-issue keeps stayed inline (L-033 rule b/a).**
   The `tau ≠ 0.0` divide-by-zero landmine, the LD `θ ≠ LLD` trap, the
   spatial-DD-vs-angular-β conflation ⚠, the #247 Mode-10 external-slope gap,
   the #242 don't-merge-the-diagonals constraint — all kept at point of use;
   only their derivations went to pointers.

## MOVED candidates

**None.** (The book is complete for all three files.)
