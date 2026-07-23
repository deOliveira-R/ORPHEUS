# P2-F — curvilinear ψ½ operators code-prose rebalance

Phase 2 of #231, batch **P2-F** (the curvilinear-SN subsystem). Two files.
Branch `docs/sn-doc-architecture`. Docstring/comment-only; ZERO code-token
changes (proven: 2831/2831 and 2679/2679 code tokens IDENTICAL vs HEAD).

Pre-edit prose:
- `orpheus/sn/sweep/pole_angular_closure.py` — 1604 ln, 807 docstring + 283 comment.
- `orpheus/sn/operators/radial_characteristic.py` — 1456 ln, 767 docstring + 128 comment.

## Headline finding (confirms the pilot's operator-file prediction)

**ZERO MOVED. The book is complete for both files** — every derivation /
design-rationale concept has a TWIN home in the theory corpus, overwhelmingly
in the **`sn-direct-seed-*` family** (`curvilinear_one_group.rst`, the #282
route-(a) close-out chapter, lines 3115–4114) plus `balance-curvilinear`
(the curvilinear balance + α-recursion + ΔA/w proof + M-M weights) and the
τ/closure-ownership sections (`sn-tau-closure-owned`,
`sn-closure-c-constants-owned`, `sn-tau-step-c-closeout`). Verified by reading
each landing section before assuming any concept novel (pilot discriminator #1).
The rebalance is a docstring→pointer contraction, not a relocation.

## Two structural facts that changed the gate set vs the pilot

1. **`radial_characteristic.py` (file 2) IS automodule'd** at
   `docs/api/discrete_ordinates.rst:137` (`:members: :show-inheritance:
   :noindex:`). Its docstrings ARE RENDERED (unlike the pilot's
   `scattering.py`, which was not automodule'd anywhere). Consequence: I ran a
   **forced `-E -W` build gate** on file 2. `:noindex:` makes the module's own
   xref targets invisible, so its `:class:`/`:meth:` refs render plain-text by
   convention (L-002) — the literal `…rst §label` pointer form is exactly
   right, and a `.. warning::` admonition renders as a visible box.
   Baseline = **0 warnings**; post-edit = **0 warnings** (unchanged).
2. **`pole_angular_closure.py` (file 1) is NOT automodule'd anywhere** (grep
   clean; `orpheus.sn.sweep` is not package-automodule'd either) — invisible
   to the build, exactly like the pilot. Its module docstring carried **three
   `.. math:: :label:` blocks** (`pole-redistribution-balance`,
   `legacy-mm-symmetric-interpolation`, `dd-angular-recursion`). Grep-verified:
   each is defined NOWHERE ELSE and is NOT a `verifies(...)`/`catches(...)`
   target → **safe to cut** (orphans no `:eq:` target, no test edge).

## Pointer convention (per pilot, ratified for the phase)

Literal greppable `docs/theory/<part>/<file>.rst §<label>` (or `§<label>` when
the path is already named in the same sentence). Every cited label verified to
resolve as a `.. _anchor:` OR `:label:` eq-label (grep gate — 12/12 resolve;
the one `§13` a script flagged is the PRE-EXISTING "#226 taxonomy §13" prose
citation in `inverse`'s docstring, present at HEAD, kept verbatim — not a
theory pointer). Haiku's classify-pass anchors were RE-VERIFIED: several were
wrong (`sn-balance-curvilinear` is actually `balance-curvilinear`, no `sn-`
prefix; the invented "DESIGN" verdict category was folded into CONTRACT/TWIN).

## Posing check (no defect)

Both files already pose the honest ratified record `A = L + C − S − B`
(`radial_characteristic.py` module head line 24; the A_BB/A_AB/A_BA/A_AA 2×2
grid). NO fused-L or F-on-LHS defect to harmonize (contrast the P2-A pilot).

---

## Adjudication — File 1 `pole_angular_closure.py`

One row per edited block. CONTRACT-kept blocks not touched (properties, short
accessors, the three ABC abstract-method contracts, the neutral-closure
methods, `default_angular_closure_class`) are omitted.

| pre-edit block | verdict | destination / record | what was kept / cut |
|---|---|---|---|
| Module "Why this abstraction exists" (3–79) | TWIN | `curvilinear_one_group.rst §balance-curvilinear` (+ `§mm-weights`/`§wdd-face`/`§pole-mm-recurrence`), `§sn-apply-sweep-equivalence`, `§sn-direct-seed-solve`/`§sn-direct-seed-r12a` | KEPT: what a closure is + the M-M recurrence one-liner + the route-(a) seed summary. CUT: the 3 `:label:` math blocks + the legacy-interp factor-of-2 derivation → pointers |
| Module "single strategy contract" (81–127) | CONTRACT + HISTORY-cut | `§sn-pole-angular-closure-protocol` | KEPT: the ABC contract (3 methods + accessors + `cls(sn_mesh)` + one-strategy-per-level). CUT: Protocol→ABC #236/#248 retirement + LegacyTauSymmetric/BaileyFlatFlux ablation history → pointer |
| Module "Hébert citation correction" (129–148) | HISTORY | `§sn-pole-angular-closure-protocol` (carries the correction) | CUT the wrong-2009-Bailey-paper narrative → folded to a one-line "primary Hébert §3.9.4; BMC 2010 auxiliary" in the pointer |
| Module References + "See also" (150–179) | CONTRACT-kept, trimmed | — | KEPT: Hébert primary + local PDF path + BMC aux + the α-crosswalk `§normalization-alpha-crosswalk` pointer. CUT: design-memo/Phase-A-closeout plan pointers + the retired `BoundaryFaceFlux` / `DiscretizationScheme` see-also |
| Pre-class comment (Protocol retirement, 203–214) | HISTORY | `§sn-pole-angular-closure-protocol` | KEPT: "this ABC is the SOLE contract (Cardinal Rule 2)". CUT: the orphaned-and-divergent-Protocol narration → pointer |
| c_out/c_in ownership comment (341–357) | CONTRACT + TWIN-cut | `§sn-closure-c-constants-owned` | KEPT: the **α-index formulas** `c_out=α_{m+1/2}/τ`, `c_in=(1−τ)/τ·α_{m+1/2}+α_{m−1/2}` (curvilinear convention at point of use) + neutral-zeros. CUT: the Pattern-7/single-source ownership narrative → pointer |
| `tau_per_ordinate` property (432–452) | CONTRACT + TWIN/HISTORY-cut | `§sn-tau-step-c-closeout` | KEPT: `(N,)` shape, BMC-Eq-43 formula, derived c's, Cartesian 1.0, read-only cache, `CellVisit.tau` consumer. CUT: consumers list + Step-C-retirement + Leg-1 narrative → pointer |
| `_build_per_ordinate_cache` (379–400) | CONTRACT + HISTORY-cut | — | KEPT: O(1) / read-only / LAST-`__init__`-step precondition. CUT: Fix-1/B3/Step-C provenance |
| Abstract-method-contract comment (454–469) | CONTRACT + HISTORY-cut | — | KEPT: matvec typed against the ABC + loose `object` return so M-M→`_MMHalfGrid`/Identity→`None`. CUT: `DiscretizationSchemeBase` precedent + the pyright-history |
| `_MMHalfGrid` comment (536–542) + class docstring (546–597) | CONTRACT + HISTORY-cut | — | KEPT: the **`faces[g,m,i]` off-by-one index convention** + the two-consumer/accessor semantics + storage shape (curvilinear KEEP — off-by-one is the ERR-class the rule protects). CUT: Step-1.5/1.7 storage-flip + Issue-248 `__call__` history |
| τ-producer comment (663–681) | **CONSTRAINT-kept** + TWIN/HISTORY-cut | `§sn-tau-step-c-closeout` | KEPT (LATENT-TRAP, pilot #3/#4): "STRUCTURAL INDEPENDENCE (vv L11) — do NOT tidy the arithmetic into a call to `contamination.morel_montry_weights`; that would collapse the Leg-1 cross-check into a tautology". CUT: τ-ownership + Step-C teaching → pointer |
| `morel_montry_tau_raw_per_level` (687–718) | CONTRACT + TWIN-cut | `§sn-direct-seed-r12a` | KEPT: the τ_raw formula, the split rationale (clamp maps 0→½), the bit-exact trichotomy (0/1/(0,1)) + per-geometry edge conventions. CUT: the per-rule "why" (rank-duplicate / dead value / genuine state) → pointer |
| M-M kernel comment (831–840) | CONTRACT + HISTORY-cut | — | KEPT: pure-algebra, all data via args, mesh-bound composes it, tests call `compute_psi_half_per_level`. CUT: free-fn→staticmethod→C5 provenance |
| `compute_psi_half_per_level` seed param (908–914) | CONTRACT + HISTORY-cut | — | KEPT: seed shape/`None` contract. CUT: the route-(a) strategy-indirection parenthetical → one-line pointer to `precompute_psi_state` |
| Pre-`MorelMontry` comment (932–938) | CONTRACT + HISTORY-cut | — | KEPT: public surface (`precompute_psi_state`+`cell_contribution`; verification via `compute_psi_half_per_level`). CUT: #248 `__call__` retirement |
| `MorelMontry` route-(a) note (944–955) | **CONSTRAINT-kept** + TWIN-cut | `§sn-direct-seed-r12a` | KEPT: the seed dispatch (carrying→ψ½ STATE; non-carrying→`edge_extrapolated_seed`). CUT: the trichotomy detail + the retired `PsiHalfAngleSeed` zoo → pointer |
| `MorelMontryAngularSweep` class docstring (956–1021) | CONTRACT + TWIN-cut | `§pole-mm-recurrence`, `§sn-apply-sweep-equivalence` | KEPT: mesh-bound contract (coeffs precomputed, methods read `self`), the recurrence + redistribution fold as **inline one-liners**, the per-level structure, Parameters. CUT: the two display-`math` blocks + the verbose apply≡sweep teaching + the legacy-`sn_mesh=None` retirement |
| `__init__` mesh-binding + R12a comments (1053–1068) | CONTRACT + TWIN/HISTORY-cut | `§sn-direct-seed-r12a` | KEPT: mesh-binding REQUIRED + the **R12a construction-ORDER safety constraint** (predicate needs only `(quad, coord)`, both bound first). CUT: the τ_raw∈(0,1) definition + kernel-at-module-level restatement |
| `__init__` per-level-partition + Step-A comment (1075–1085) | CONTRACT + HISTORY-cut | (τ-producer note above) | KEPT: the per-level partition (M-M's concept, sphere 1 level / cylinder μ-levels). CUT: the Step-A/C τ-ownership restatement |
| `edge_extrapolated_seed` (1187–1215) | CONTRACT + TWIN-cut | `§sn-direct-seed-r12a` | KEPT: the **seed extrapolation formula** + convergence properties (O(Δμ²), linear). CUT: the per-rule trichotomy (product t=0 / level-symmetric dead weight) detail → pointer |
| `IdentityAngularClosure` comment (1434–1442) + class "Why it earns its keep" (1458–1476) | CONTRACT + HISTORY-cut | `§balance-curvilinear` | KEPT: geometry-blind zero contribution (Cardinal Rule 2), the "replaced a `pole_angular_closure is None` matvec branch (Pattern-4 leak)" one-liner. CUT: the pre-Phase-2.8 `if curvature=="cartesian"` code-block + the C4/#220 face-inventory design story → pointer |

**Kept whole (CONTRACT, no edit):** `PoleAngularClosureBase` class docstring
(the registration contract), `beta_first_order_consistent` trait docstrings ×3
(trait declarations a modifier needs — guard the opt-in-False default, pilot
#4), `_gather_per_ordinate` (+ the `np.zeros`-not-`np.empty` elegance rationale),
the three ABC abstract-method docstrings, `morel_montry_tau_per_level` (formula
+ clamp / division-by-zero conventions), `_psi_half_grid_single_level` (the
verbatim recurrence), `precompute_psi_state` (seed-dispatch + the `None`-seed
caller gotcha), `cell_contribution`, `angular_adjoint` (reverse-mode contract),
`_edge_seed_stencil`, `__repr__` test-contract, the neutral IdentityClosure
methods, `default_angular_closure_class`.

---

## Adjudication — File 2 `radial_characteristic.py`

| pre-edit block | verdict | destination / record | what was kept / cut |
|---|---|---|---|
| Module roadmap (3–25) | CONTRACT + HISTORY-cut | `§sn-direct-seed-solve` | KEPT: the four-block map (A_BB/A_AB/A_BA/Fold) + `A_AA = L+C−S−B`. CUT: the step-4c/migration/`transport→sn`-ban narration per bullet → pointer |
| Module physics/posed/solved (27–67) | CONTRACT + TWIN-cut | `§sn-direct-seed-pole-straight-characteristic`, `§sn-direct-seed-block-triangular`, `§sn-direct-seed-solve` | KEPT: the operator definition (μ∂_r+σ_t at closed rays), domain=codomain, the **two BC contract** (r=R Dirichlet / r=0 pole continuation), direct-march=resolvent + single-source. CUT: the straight-characteristic derivation + the 4e/B.2c bridge narration → pointers |
| Module "Scope"→**"Realized surface"** + "ONE solve orchestration"→**"Single source"** (69–113) | CONTRACT + HISTORY-cut | `§sn-direct-seed-solve` | KEPT: the realized-surface list (apply/apply_transpose/solve/solve_transpose/inverse) + the **single-source LATENT-TRAP** ("Do NOT re-implement the march at a call site"). CUT: step-4b/6 extract-not-twin + 4e-e2 un-weave narration → pointer. **Renamed two module section headings** (see judgment call #4) |
| `RadialCharacteristicOperator` class docstring (204–224) | CONTRACT + narration-cut + **STALE-REF FIX** | — | KEPT: endomorphic / domain=codomain / `(B,B)` slot / realized surface. CUT: B.2c re-type + 4e "unified layout retired". **FIXED** the now-stale back-ref "the module docstring's 'Scope' section" → "'Realized surface' section" |
| `…Operator.apply` (336–350) | CONTRACT + HISTORY-cut | — | KEPT: the round-trip **numeric evidence** (apply∘solve defect → 0.0 exactly; cell round-trip ~FP ULP). CUT: the step-6 "kernel's production forward" provenance |
| `…Operator.solve` (462–474) | CONTRACT + HISTORY-cut + **STALE-REF FIX** | — | KEPT: the two-leg Carlson march contract + orientation-in-data + slot-key. CUT: 4e-e2. **FIXED** "see the module docstring 'The ONE solve orchestration'" → "'Single source'" |
| `RadialCharacteristicSeeding` (A_AB) class docstring (640–738) | CONTRACT + TWIN-cut + HISTORY-cut | `§sn-direct-seed-solve` | KEPT: block-role + domain/codomain, the **cell-local-vs-A_BB distinction**, the isolation mechanism (zero-bulk fwd / discard-cotangent transpose) + the placement formula, the shared-kernel + tracked-twin-**pinned** one-liner, Parameters (the σ_t-independence). CUT: the M-M recurrence display-`math` block + the verbose tracked-transient-twin retirement narrative → pointer |
| `…Seeding.apply` (789–803) | CONTRACT + duplicate-cut | — | KEPT: the per-method contract (build seed-only grid; place −(ΔA/w)c_in ψ/V; non-carrying zero). CUT: the isolation-by-linearity teaching (now in the class docstring) |
| `RadialCharacteristicReconstruction` (Fold) class docstring (956–1033) | CONTRACT + HISTORY-cut + TWIN-cut | `§sn-direct-seed-source-fold` | KEPT: the fold-role in A_BA, `is_adjointable`/`is_invertible` predicates, the single-source constraint, the **`source_from_angular` "not a twin" distinction**, the broadcast/corners-zero convention, Parameters. CUT: the step-4c migration history + the fold display-`math` derivation → pointer |
| `RadialCharacteristicEmission` (A_BA) genericity (1207–1225) | CONTRACT + **HAZARD→warning box** | (rendered `.. warning::`) | KEPT the kernel-genericity CONTRACT; **PROMOTED the fission double-apply HAZARD to a rendered `.. warning::` admonition** (L-010 prophylactic-warning — file 2 is automodule'd so it renders as a box): "do NOT route fission's PRODUCTION seed through this operator — would DOUBLE-apply K∘integrate" |
| `…Emission` LIFT/typing/transpose (1227–1264) | CONTRACT + HISTORY-cut | `§sn-direct-seed-solve` | KEPT: within-group-gain-vs-fission-outer seam, System-A→B typing + Pattern-4 no-double-count, the `(B,A)`-slot native consumption, the transpose-role + the **nonzero-seed-cotangent-gate rationale** (a present-zero seed hides a lost pullback). CUT: THE LIFT before/after + B.2b/B.2d adapter narration |
| `…Emission.apply` / `apply_transpose` (1334–1441) | CONTRACT + tag-cut | — | KEPT: the method contracts (integrate→K→fold; the ∫dμᵀKᵀFoldᵀ pullback + weight broadcast). CUT: the residual 4e/B.2b tags |
| Trailing retired-adapter comment (1452–1456) | HISTORY | — | CUT the retired `_RayEmissionFullFieldGain` B.2b/B.2d narration → the one LIVE fact (A_BA sits at the `(B,A)` slot of `build_within_group_system`) |

**Kept whole (CONTRACT, no edit):** the `_EmissionKernel` Protocol docstring
(the by-capability typing rationale — guards the plausible wrong "annotate as
`LinearOperator`" simplification, pilot #4), the Sourcing + References module
sections, all four operators' predicate/`__init__`-attribute/`__repr__`
docstrings, `_require_member_composite`, `_check_mesh` (mesh-identity
invariant), `…Operator.apply_transpose` (the precise Euclidean-vs-Hilbert
`.H` distinction, L-010), `…Operator.inverse` (the #226-taxonomy generic
solve-backed inverse), `…Operator.solve_transpose` (the A_AB-transpose
separation), `…Reconstruction.apply`/`apply_transpose`.

---

## Verdict counts (primary verdict per edited block)

| | File 1 | File 2 |
|---|---|---|
| TWIN (teaching → pointer) | 7 | 5 |
| HISTORY (campaign narration cut) | 9 | 6 |
| CONSTRAINT-kept (LATENT-TRAP imperative preserved inline) | 2 | 1 |
| CONTRACT (kept; ~embedded TWIN/HISTORY trimmed) | rest | rest |
| POSING-HARMONIZE | 0 | 0 |
| **MOVED** | **0** | **0** |
| in-file STALE-REF fixes (section-rename fallout) | 0 | 2 |
| HAZARD → rendered `.. warning::` | 0 | 1 |

## Before/after metrics (measured; `ast` docstrings + `tokenize` comments)

| metric | file 1 before → after | file 2 before → after |
|---|---|---|
| total lines | 1604 → 1375 (−229) | 1456 → 1340 (−116) |
| docstring lines | 807 → 628 (**−179, −23%**) | 767 → 654 (**−113, −15%**) |
| comment lines | 283 → 233 (**−50, −18%**) | 128 → 125 (−3) |
| **code tokens** | **2831 → 2831 (IDENTICAL)** | **2679 → 2679 (IDENTICAL)** |

`git diff --numstat`: file 1 = 237 ins / 466 del; file 2 = 213 ins / 329 del
(insertions reflect paragraph REWRITES, not additive text). File 2's comment
count barely moves because that file is docstring-heavy (its campaign narration
lived in class docstrings, not comments).

## Gates (all pass)

1. `import` both modules → **OK**.
2. `pytest --collect-only -q` → **6652 tests collected, EXIT=0** (unchanged).
3. **Code-token invariance** — file 1 2831==2831, file 2 2679==2679, sequence
   IDENTICAL vs HEAD (non-STRING/COMMENT/layout tokens) → proves zero code change.
4. **Every cited §label resolves** — 12/12 (the `§13` a script flagged is the
   pre-existing "#226 taxonomy §13" prose citation, kept verbatim).
5. **NO theory-page edits** — `git diff` carries ZERO `docs/` files.
6. **BONUS `-E -W` build** (mandatory here — file 2 IS automodule'd/rendered):
   baseline 0 warnings → post-edit **0 warnings** (WARNING/ERROR/CRITICAL set
   unchanged), `build succeeded`, EXIT=0. The new `.. warning::` box + all
   literal pointers render clean.

## Hardest judgment calls (calibrate future curvilinear batches)

1. **File 2 is automodule'd → the build gate is MANDATORY, unlike the pilot.**
   Always grep `automodule::` for each batch file first. A rendered operator
   file (`:noindex:` or not) needs the `-E -W` build; a HAZARD/latent-trap in
   a rendered file should become a `.. warning::` box (an improvement over an
   inline comment — L-010 prophylactic-warning), which the pilot's non-rendered
   file could not do.
2. **File 1's module-docstring `:label:` blocks were cut-safe ONLY after the
   3-way grep** (not automodule'd + defined nowhere else + not a verifies/
   catches target). Never cut a docstring `:label:` without that check — a
   rendered/automodule'd file would orphan an `:eq:` target or a test edge.
3. **Curvilinear-KEEP held for every convention at point of use.** The
   α-index formulas (`c_out=α/τ`, `c_in=(1−τ)/τ·α+α`), the `faces[g,m,i]`
   off-by-one, the τ_raw formula, the seed-extrapolation `t` — all KEPT. This
   is the subsystem where sign/index conventions have bitten hardest (ERR-026
   family); when the teaching AROUND a formula is TWIN but the FORMULA is the
   convention a file-local modifier depends on, keep the formula + pointer the
   teaching.
4. **Section renames inside a rendered file leave stale back-refs — fix them.**
   Renaming file 2's module headings ("Scope of this realization" → "Realized
   surface"; "The ONE solve orchestration" → "Single source") staled two class-
   docstring back-references. Grep the file for the OLD heading text after any
   rename and repoint (found + fixed both; final grep clean).
5. **Two LATENT-TRAP imperatives are the load-bearing keeps** (pilot #3/#4):
   (a) file 1's "do NOT tidy the τ arithmetic into a call to
   `contamination.morel_montry_weights`" (would collapse the Leg-1 cross-check
   into a reference-contamination tautology — vv L11); (b) file 2's "do NOT
   re-implement the march at a call site" (single-source orchestration).
   Both are constraints a file-local modifier would violate WITHOUT the line —
   the CONTRACT test — so they stay inline even though their EXPLANATION is TWIN.

## MOVED candidates

**None.** As the pilot predicted for operator files: the `sn-direct-seed-*`
family (route-(a) close-out) + `balance-curvilinear` + the τ/closure-ownership
sections are the SSOT and are complete. Three concepts FELT possibly-novel and
each was fully TWIN after reading the landing section: the M1/M2/M3 pole-state
metric (fully at `sn-direct-seed-solve`, incl. ERR-067/Mode-12); the R14 "why
ALL Legendre moments" fold argument (fully at `sn-direct-seed-source-fold`,
incl. the isotropic-flux→linear-in-μ-source derivation); the R12a carrying-
level trichotomy (fully at `sn-direct-seed-r12a`, incl. the bit-exact table).

## Reported discrepancy (not a defect, not mine)

The working tree carries **other P2 batches' uncommitted docstring/comment
changes** (`streaming.py`, `scheme.py`, `solver.py`, `diamond.py`,
`linear_discontinuous.py`, `boundary.py`, `augmented_mesh.py`,
`iteration.py`, `operator.py`, `loss_representation/*`). My edits are isolated
to the two P2-F files (code-token-invariance proven per file). Those other
files are the parallel/prior batches the main agent stages; not adjudicated here.

Commit NOTHING — main agent reviews, certifies, commits.
