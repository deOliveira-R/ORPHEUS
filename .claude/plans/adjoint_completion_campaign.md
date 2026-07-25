# Adjoint completion campaign — #310 (residue) → #276 A4/A5/A6 → ch15

**Branch model:** per-phase `feature/sn-adjoint-*` branches cut from main; ff-merge
when green; push gated on the user's explicit ask. **Execution model:** surgical,
main-agent-direct, user-steered (NO method-implementer); AskUserQuestion checkpoint
at every phase boundary. **The binding gate spec:**
`.claude/plans/residue_verification_spec.md` (test-architect, 2026-07-23) — its
gates SHAPE the carve; read the relevant § before each commit.

**Sequencing:** #310 (capability: LD + multi-D transpose) → #276 A4 (φ* daggered
posing) → A5 (carrier) → A6 + ch15 (docs; closes #276, discharges #309's sn/adjoint
arm). Downstream consumer: #281/#51 (adjoint-weighted homogenization).

## The one cardinal rule (spec §0.4)

Object-level gates ONLY: `eig(Aᵀ) = eig(A)` makes every spectral/keff functional
DESIGNED-GREEN on the whole adjoint mutation class (Mode 12). Correctness leans on
structurally-independent OBJECT oracles — the SymPy `ld_ubld` transpose-oracle,
dense Euclidean `Mᵀ` column-probed off the FORWARD apply, assembled-`Mᵀ` LAPACK,
frozen pre-carve baselines. G-reciprocity is the composite canary, never the
keystone.

## Contract rulings (test-architect findings, spec §1.2/§13 — ratified at the
Phase-0→C1 checkpoint)

1. **The kernel-pair VJP covers the SPATIAL cell relation ONLY** —
   `(res_bar, psi_out_bar) → (psi_bar_cot, psi_in_cot [, source_cot])`. The
   Morel–Montry angular-numerator cotangent (`numer_bar`) stays on
   `PoleAngularClosure` (its `cell_contribution` transpose), matching the landed
   spatial/angular separation. A kernel folding the whole `visit` would drop
   `numer_bar` and red the sphere/cyl frozen baselines.
2. **`has_transpose_kernel` is REGISTRATION-COUPLED** — derived from "a transpose
   kernel is registered", so declared-True-with-no-kernel is unrepresentable.
3. **The LD `.H` metric MUST carry the moment mass θ** —
   `G_bulk_LD = V·w_n ⊗ diag(1, θ)` (the LD forward matvec is `M⁻¹(Aψ⃗−R)` with
   `M = ⊗diag(1,θ)`, `linear_discontinuous.py:613-620`). Without it, `.H` is a
   WRONG adjoint on the slope DOF AND reciprocity is Mode-12-blind to a slope-row
   transpose (the ERR-067 metric-repair family). A red on the C2 metric row is a
   production bug to fix, never a test to relax.

## Flip-safety (spec §10)

Flags are construction-time; walks raise apply-time. **Every flag flip is atomic
with its kernel/walk + the rewritten `test_ld_adjoint_deferral.py` pins + the
extended gate rows** — a flip before its walk is the predicate lie that file was
minted to kill.

## Stages

- **Phase 0 — setup ✅ 2026-07-24**: docs branch ff-merged → main @ `ec74be50`
  (185 commits; the doc campaign is on main, push pending user ask); stale
  pointers `docs/sn-doc-architecture` + `feature/sn-adjoint-transport` (tip
  `f9036860`, git-verified merged) deleted; **#310 filed** (the residue tracker;
  deferred-out items named there: G-S schedule reverse, d≥3 interior interleave,
  curvilinear LD, the solver.py:419 LD certificate skip); test-architect spec
  delivered (848 ln, 14 §§); `feature/sn-adjoint-residue` cut.
- **C1 — DD kernel-pair relocation** (spec §2): mint the registered pair on
  `DiscretizationSchemeBase`; relocate the hand-coded DD transpose
  (`loss_representation/__init__.py:3051-3075`) into
  `DiamondDifference` — bit-identity (frozen `walk_matvec_*.npz`, nulp=1 +
  DriftWarning escalation); Mode-11 sentinel wraps BOTH the new kernel AND the
  closure's angular transpose (count>0 sphere/cyl, ==0 slab); retirement grep
  (old hand-code deleted, zero twins). NO flag flips.
- **C2 — LD-slab adjoint** (spec §3+§4): SymPy VJP≡`AᵀM⁻¹` oracle in `ld_ubld.py`
  + dense numeric cross-check; the three NEW-algebra gates (mass-inverse order,
  slope fold, involution); LD 1-D reverse walk (moment-tailed cotangents); the
  TWO LD guards lifted together onto the trait; **flip
  `LinearDiscontinuous.has_transpose_kernel`** (= register the kernel — the C1
  derivation makes the flip structural) atomic with the deferral-file
  rewrite + LD reciprocity rows (θ-mass metric, ruling 3) + `loss_transpose_solve`
  G1/G2 + assembled-`Mᵀ` LD-slab. **FOLDS IN #311** (enforcer condition): extract
  the w-generic `outgoing_face_from_average_transpose` primitive FIRST — DD
  reroutes byte-identically (`w=½` exact), LD composes it on first landing
  (`w=1/(1+k)`), the two `psi_half_angle_seed` transposes follow; + the
  deferral-file `s/declare/derive/` NIT. ⏸ COMPACTION (R1 complete).
- **C3 — multi-D reverse `walk_full` ORACLE** (spec §5): reversed levels +
  `_CellResidualTranspose` + transposed addressing + boundary in↔out swap; the
  NEW 2-D forward-apply column-probe (the 1-D `_probe_augmented_matrix_one_group`
  does not lift); reverse-walk spy + AST tripwire (sibling of
  `test_one_dim_loop_walk`). No production flag flips (oracle arm).
- **C4 — multi-D reverse `walk_windowed` PRODUCTION**: reverse full≡window oracle
  (mirror of forward discipline); **flip `has_transpose_walk`** (ScanMarch-2D +
  wavefront) atomic per flip-safety; un-exclude cart2d in
  `test_removal_form_matvec_sweep.py`; cart2d reciprocity/coherence/phase-C rows;
  frozen baselines for the new reverse walks.
- **C5 — LD-2D** (spec §6): the moment-frame-involution reciprocity row (the
  likeliest sign-error site, ERR-066 family); LD-2D dense-`Mᵀ`; rectangular
  nx≠ny grids (Mode-2 axis-swap catchers). **Doc blast rides C3–C5**: rewrite
  the `loss_representation.rst` §loss-rep-orientation-two-frames coverage
  matrix; re-point the ~40 `#280`-citing deferral comments at #310. Closes #310.
  ⏸ COMPACTION (residue complete).
- **A4 — daggered posing + φ*** : activate `KEigenvalue(A_loss.H, S†-wrap, F.H)`
  → unchanged `power_iteration`; `solve_sn_adjoint`-shaped entry (shape steered
  in-session). Gates P1.2–P1.5 per `p6_adjoint_verification_spec.md:230-333`
  (P1.2 duality asymmetric-2G cross-group; P1.3 triple equality ∞+slab+sphere,
  ∞-only=REJECT; P1.4 4G adjoint spectrum vs `eig(Mᵀ)` — extend
  `derivations/common/eigenvalue.py:48` via `dominant_eigenpair(M.T)`; P1.5
  bi-orthogonality full cross-Gram). Landmines at spec :461-467. Promote the
  trailing-axis diagnostic → `test_scattering_adjoint.py`. ⏸ COMPACTION.
- **A5 — the φ* carrier**: type-vs-property adjudication (three shapes: property
  on `Solution` / `AdjointFlux` type / second `Solution` from `solve_sn_adjoint`);
  `Solution.homogenize`'s degenerate-φ* comment is the latent consumer. USER
  rules the shape.
- **A6 + ch15 — docs + close**: author `methods/sn/adjoint.rst` whole (exact
  discrete-transpose route rationale; μ-reversal kept only as a slab oracle;
  three-transposes landmine; posing; φ*; consumers); V&V slice rows →
  `verification/sn.rst`; wire/adjudicate #309's sn/adjoint arm; vv-status per
  wired⟹no-sentinel; update the `eigenvalue.py:26-29` seam note; close #276.
  archivist authors, qa + elegance-enforcer review. ⏸ COMPACTION (campaign done).

## Campaign-level gates (every phase)

Full-tree serial `-m "not slow"` green (pre-merge) · `-E -W` Sphinx exit 0 ·
audit exit 0 (orphans 0, phantoms 0) · pyright production floor 1 · enforcer
review · qa on math-bearing phases · every new gate class mutation-verified
(Mode-10) · k never sole evidence (Mode-12).

## Status log (append per stage)

- **A4 OPEN — branch `feature/sn-adjoint-a4` from main @ `02da1714`.**
  USER RULINGS at the A4 opening checkpoint (2026-07-24): (1) **BOTH
  entries land now** — `solve_sn_adjoint` (eigenvalue) AND
  `solve_sn_adjoint_fixed_source`, module-level siblings mirroring the
  forward family pattern; (2) **both return the unified `Solution`**
  (provisional carrier — A5 adjudicates the φ* shape; migration surface
  = the entries + gates).  Survey findings banked: `SourceIteration`
  is already composite-ready (Ravellable bridge); `FullField` carries
  Vector + `copy()`; the swap law makes `LC.H` invertible with
  `_seeded_inverse` unchanged; the daggered triple =
  `KEigenvalue(LC.H, (S+B_a).H, F.H)` off `build_within_group_system`'s
  splitting; the activation's exact new work = KEigenvalue's four
  ndarray-bound spots (asarray-guess, np.sum ×2 estimators,
  np.linalg.norm converged), the `FissionOperator.apply_transpose`
  composite variant (its docstring names this moment), and the P1.3
  sphere leg's coupled-grid dagger (F lift to the CoupledField).
- **A4-3b ✅ — THE φ* CERTIFICATION BATTERY**
  (`tests/sn/solve/test_sn_adjoint_certification.py`, 12/12, ~112 s).
  P1.3: the k triple-equality (∞ 2G+4G vs `kinf_homogeneous` — the chain
  terminates in `np.linalg.eig`) + the HET reflective slab + the SPHERE
  (the coupled daggered chain).  P1.3 mutations, each in-process
  (monkeypatch via instance-dispatch lambdas — the singledispatch
  class-descriptor cannot be setattr'd raw): F†→F shifts k @ 4G χ∦νΣf
  (precondition asserted), S†→S shifts k @ asymmetric SigS, L†→L (the
  daggered inverse riding the FORWARD sweep — "no reversal") shifts k on
  the HET leg (the ∞ legs are flat-flux-blind to it, per spec §0.6).
  P1.4: the 4G spectrum vs the corrected `(Aᵀ)⁻¹Fᵀ` reference (ψ*_cf ≠
  φ_cf materiality asserted) + THE MODE-12 PAIR EXPLICIT: F†=F leaves k
  exactly equal on ∞ (the daggered pencil (A†,F) is similar to the
  forward — designed-green) while the spectrum row reds O(1).  P1.5:
  the closed-form cross-Gram is DEGENERATE-diagonal (rank-1 F: one
  nonzero entry; both orthogonality mechanisms — Fφ_j=0 zero-right,
  χ·ψ*_i=0 zero-left — asserted; HONEST-SCOPE deviation from the
  spec's "ng distinct modes" wording recorded in the file header) +
  the right-right Gram mutation (NOT one-entry — φ*≡φ fails P1.5) +
  the SN defining-law residual row (Aᵀψ*_SN = Fᵀψ*_SN/k at the
  closed-form (A,F), atol 1e-7).  Pyright 0/0.
  **REMAINING for the A4 boundary (next session if context ends):
  the phase-end sweep — Sphinx `-E -W` + audit + matrix regen,
  elegance-enforcer review, qa review (math-bearing phase), the FULL
  serial tree (L37: sources freeze), then the boundary AskUserQuestion
  (merge / A5).  Also weigh at enforcer time: the `_g_inner`
  cross-test import in the entries file; the two `_energy_spectrum`
  spellings (entries + certification); the A5 fork is genuinely open
  (the entries return unmarked `Solution`).**
- **A4-3a ✅ — THE COUPLED (SPHERE) DAGGERED POSING.**  USER RULING:
  space-typed ZeroOperator.  Findings + carves: (1) the fission ray fold
  is ALREADY an operator — `RadialCharacteristicEmission(sn, F.kernel)`
  is kernel-generic (`A_BA_fission = Fold∘F.kernel∘integrate`, transpose
  included — the extract-shared-primitive discipline paying off); my
  first structural-zero B-row was WRONG-forward (the forward coupled
  fission source carries the q½ fold — `_radial_characteristic_fission_seed`)
  and the typed guard caught it before the k-mismatch would have; on the
  eigen-M operator the fold row BELONGS in the posing (HAZARD 5 keeps it
  out of the WITHIN-GROUP gain, not out of M).  (2) `ZeroOperator`
  already carried `codomain_zero` (the pre-#208 role-typed hook — the
  ruling was half-realized in the algebra); A4 adds the SYMMETRIC
  `transpose_zero` (the domain's dual-role zero) — first exercise of the
  zero slot's transpose, docstring era-note updated.  (3) the (B,B)
  fission slot = `ZeroOperator(codomain_zero=transpose_zero=ray
  source_zeros_on)` (w=0 closed rays never source fission).  (4) the
  ray-leg DUALITY-TYPING fix: `RadialCharacteristicOperator.
  solve_transpose` output re-classed source→FLUX members (dual-of-source
  = the adjoint ray flux — the exact sibling of A4-1's
  `InvertibleOperator.solve_transpose` fix; caught by the same
  `_check_partner` guard on the daggered sphere's first contact); the
  B.2c container pin rewritten to the duality-typed truth same-commit.
  Evidence: FORWARD coupled control `KEigenvalue(resolvent, gain,
  F_posed)` matches solve_sn sphere k @ 3.9e-10 (14 outers — validates
  the fission-fold posing pre-dagger); DAGGERED coupled k_adj vs fwd @
  5.6e-11; the ENTRY sphere dk 4.5e-10 with the ray member packaged.
  Batteries: 203/203 (psi-half + coupled + curvilinear-282 + transpose-
  solve + entries + iteration); pyright 0/0 production, legacy test file
  100 ≡ 100 HEAD.
- **A4-2 ✅ — THE ENTRIES.**  `solve_sn_adjoint` + `solve_sn_adjoint_fixed_source`
  land per the user rulings (both now; unified `Solution` return).  Shapes:
  the eigen entry mirrors `solve_sn`'s surface (no scheme param — the
  forward eigen has none; no inner_solver/schedule — ONE daggered inner
  realization), builds `_adjoint_posing_parts` (mat_xs →
  `build_within_group_system` → summed gain → `FissionOperator.
  from_solver_data` with `full_field_space` wired → the coupled F-lift
  `[[F,∅],[∅,Zero]]` on carrying meshes) and runs `KEigenvalue(resolvent.H,
  gain.H, F.H)`; the fixed-source entry takes `detector_response`
  ((ng,*spatial) Σ_d OR a FullField) with the **angle-flat dual lift** —
  NO w_n, NO 1/W: under G = V·w_n the plain broadcast is exactly the dual
  of the scalar-flux extraction ⟨1_Ω Σd, ψ⟩_G = ⟨Σd, φ⟩_V (the forward's
  /W iso lift is the dual of a DIFFERENT map — source injection; the
  asymmetry IS P1.2's content) — and drives the daggered `SourceIteration`
  via the PROMOTED `seeded_inverse` (public: 3rd consumer; green_operator
  had already been importing the underscore name cross-module).
  Carrying-mesh fixed-source = typed LOUD refusal (no gate/consumer yet —
  ships as refusal, not unexercised; eigen covers carrying via P1.3
  sphere).  Packaging `_package_adjoint_solution` mirrors the forward
  tail (cell-average view, AVERAGE_MOMENT slice, boundary from the
  converged composite, φ* = Σw_nψ*_n, honest converged flag).  Gates
  (`tests/sn/solve/test_sn_adjoint_entries.py`, 5): the entry k
  triple-equality + spectrum + ν̂Σf-degeneracy branch + ∞-ISOTROPY of ψ*
  (angle-flatness — the cheap G-consistency signal) + Solution packaging
  contract + **the cross-group cross-region DUALITY row** (fast source
  left / thermal detector right, vacuum het slab: ⟨ψ*,q⟩_G = ⟨q*,ψ⟩_G at
  1e-7 via the independent `_g_inner`, AND ⟨q*,ψ⟩_G = hand Σ V·Σd·φ at
  1e-10 pinning the lift) + both refusal rows.  Entries 5/5; pyright 0/0
  (solver.py + the test).  A4-3 hardens this file's duality row into the
  full P1.2 (mutations) and adds P1.3/P1.4/P1.5.
- **A4-1 ✅ — THE DAGGERED POSING ACTIVATION.**  `KEigenvalue(LC.H,
  (S+B_a).H, F.H)` runs through the UNCHANGED `power_iteration` on typed
  composites — k_adj == k_fwd == solve_sn @ ~5e-12 AND the ∞-medium
  adjoint SPECTRUM matches the corrected closed form digit-for-digit.
  Pieces: (a) **Carrier genericization** — `EigenvalueSolver` Protocol +
  `power_iteration` generic over the NEW UNBOUNDED `vector.Carrier`
  (numpy's stubs reject ndarray ⊨ Vector via the `__add__` bool-dtype
  self-binding overloads, so a Vector-bound TypeVar statically rejects
  every ndarray solver family the moment it is SOLVED — the bound was
  survivable before only because nothing forced solving; the missing
  bound is a stub limitation, documented at the TypeVar, NOT a wider
  contract); `KEigenvalue(Generic[V])` on the family's bounded V (body
  operator support; conforms to the unbounded Protocol); power_iteration
  stashes the previous iterate frozen-alias/ndarray-copy (the
  SourceIteration idiom).  (b) **KEigenvalue carrier-honesty** — the four
  ndarray-bound spots (`asarray`-guess, `np.sum` ×2, `np.linalg.norm`)
  rewired through the existing `_ravel`/`_l2_norm` Ravellable protocol
  (bit-identical on bare arrays).  (c) **F† composite arm**
  (`FissionOperator.apply_transpose`, the docstring's own named seam):
  the forward `(1/W)·broadcast ∘ K ∘ (w·Σ)` pulls back to
  `(w·) ∘ Kᵀ ∘ (Σ/W)` — the reduce/broadcast WEIGHTS swap sides; angular
  + scalar bulk arms, pure-bulk zero trace; 5 gates incl. the
  independent `w ⊗ hand-loop(Σχ/W)` spelling + the weight-role-swap
  discriminator (an angle-FLAT cotangent hides the swap — precondition
  pinned).  (d) **solve_transpose duality typing** — the daggered SI's
  first contact hit the `_check_partner` cross-class guard: the reverse
  sweep wrapped output in the source-sink family, but dual-of-source
  under the G-pairing IS the adjoint FLUX; `InvertibleOperator.
  solve_transpose` now wraps `AngularFlux`/`AngularBoundaryFlux` (the
  class ROLE decided at the operator boundary; rep-layer values
  untouched; 214 consumer tests green — every prior consumer was
  value-reads).  Committed gates: the honest-composite FORWARD control
  (KE's own Rayleigh on-contract for the first time — the legacy scalar
  shim kept as pre-A4 record) + the daggered smoke (k triple-equality +
  the spectrum leg + the ν̂Σf-degeneracy fail branch), both L1 in
  `test_iteration.py`.  Battery sweep 128/128; pyright: production
  floor 1 (#288 unchanged), touched files +0 (test_iteration.py 15 ≡
  HEAD 15).
- **A4-0 ✅ — the P1.4 closed-form reference + the diagnostic promotion.**
  `kinf_and_adjoint_spectrum_homogeneous` added to
  `derivations/common/eigenvalue.py` — left eigenvector via
  `dominant_eigenpair(M.T)` (the shared Perron–Frobenius primitive per
  the #276 `direct_eigenvalue` ruling; independence lives in the input
  assembly), riding the NEW `_infinite_medium_matrices` shared (A, F)
  assembly (single source with the forward sibling; forward rewired
  through it, arithmetic identical); explicit ℓ²-normalisation (the
  primitive contracts sign only — scale made comparable by
  construction).  Pin file `tests/derivations/test_adjoint_spectrum_reference.py`
  (4 foundation): the left eigen-law φ*ᵀM=kφ*ᵀ vs an INLINE-assembled
  resolvent + the output convention; k_adj==k_fwd rtol 1e-12 @ 2G+4G;
  the 4G non-degeneracy dud-guard (φ*≠φ, non-flat, χ∦νΣf precondition
  asserted).  Diagnostic `diag_276_full_scatter_kernel_ld_trailing_axis`
  PROMOTED → `test_scattering_adjoint.py::TestFullScatterKernelLDTrailingAxis`
  (LD forward-reproduction ×2 orders + the 0b3275d cells-index mutation
  teeth, monkeypatch-in-process) and the diagnostic RETIRED (git rm +
  pycache; consumers grepped: only plan/memory references).  Batteries:
  new pin 4/4; scattering-adjoint 18/18; consumer canaries 62/62;
  pyright +0 (the file's 4 = pre-existing at HEAD, proven copy-first).

- **C5 CLOSED ✅ — MERGED to main (ff, `ef71e3b5..8abd3c6c`, 5 commits:
  `15492673` C5-a flip+battery, `1f1ff09a` C5-b doc blast w/ `Closes #310`,
  `76ef830e` C5-c enforcer truth-patches, `8abd3c6c` enforcer memory) +
  phase branch deleted; user-ratified at the boundary ("Merge; compact
  before A4").  Full-tree gate: 6469 passed / 0 failed / 18 skipped /
  54 xfailed, EXIT 0, 1:45:47 (background-QoS pacing) — count EXACTLY
  accounted (6459 green-equivalent @ C4-close + 8 reverse-file items +
  2 reciprocity rows).  Enforcer verdict: code flip EXEMPLARY, every
  axis PASS (both teeth's algebra independently re-derived); its one
  blocking cluster = the L-016 recurrence ONE RING WIDER — four
  deferral-story spellings citing neither "LD-2D" nor "#310" (the
  InvertibleOperator class + method, the base trait default, the
  oracle-arm docstring) — lesson sharpened: grep the deferral story's
  PROSE spellings ("adjoint raises", "scheme-aware", "excepted"), not
  only its issue tags.  **#310 is COMPLETE — the residue is retired:
  the adjoint matvec `loss_action_transpose` covers the full registered
  scheme × representation × d grid; the deferral ledger holds ONLY R7
  (the G-S schedule-reverse, no consumer); `Closes #310` fires at
  push.**  ⏸ COMPACTION taken here (user-directed) — **A4 opens next
  session, on a fresh `feature/sn-adjoint-a4` branch from main: the
  daggered posing + φ*.**  Read FIRST: the A4 stage block above +
  `.claude/plans/p6_adjoint_verification_spec.md:230-333` (the P1.2–P1.5
  gate specs; landmine table :461-467).  The heavy lifting is banked:
  `.H` constructs on the WHOLE grid (any registered scheme ×
  representation × d — #310 C1–C5), `KEigenvalue` sits dormant at
  `orpheus/numerics/iteration.py:995` (the operator-triple posing layer,
  test-only today), and `power_iteration`'s adjoint row is the
  documented seam (`orpheus/numerics/eigenvalue.py:26-29`).  A4 = (a)
  activate `KEigenvalue(A_loss.H, S†-wrap, F.H)` → the UNCHANGED
  `power_iteration`; (b) a `solve_sn_adjoint`-shaped solver entry —
  the exact shape steered in-session with the user (surgical model);
  (c) gates P1.2 duality ⟨ψ*,q⟩=⟨ψ,Σ_d⟩ (asymmetric 2G, source/detector
  in DIFFERENT groups) · P1.3 k_adj == k_fwd == closed-form rtol≈1e-10
  on ∞ + slab + sphere (∞-only = REJECT) · P1.4 ∞-medium adjoint
  spectrum vs `eig(Mᵀ)`, 4G MANDATORY (2G Mixture A is coincidentally
  flat) — extend `derivations/common/eigenvalue.py:48` via
  `dominant_eigenpair(M.T)` · P1.5 bi-orthogonality full cross-Gram
  rtol≈1e-10; (d) promote the trailing-axis diagnostic
  (`derivations/diagnostics/diag_276_full_scatter_kernel_ld_trailing_axis.py`
  → `tests/sn/operators/test_scattering_adjoint.py`).  Mode-12: k is
  NEVER the sole evidence anywhere (eig(Kᵀ)=eig(K) — pair every k row
  with a field/pairing/object gate).  After A4: A5 (the φ* carrier —
  USER rules the shape at the checkpoint) → A6 + ch15 (closes #276).
  L37: Python sources FREEZE during long gates.**

- **C5-a ✅ @ `15492673` — THE LAST FLIP + the LD-2D gate battery (spec §6).**
  Pre-carve de-risking probe (in-process predicate narrow, no production
  edit) validated every banked claim before a line was written: forward
  LD-2D runs vacuum+reflective; pairing machine-precision (0.0 reflective);
  window ≡ full BIT-identical; dense-`Mᵀ` n=704 @ 8.9e-16 (asymmetry 5.77);
  assembled-`Mᵀ` 7.1e-15; `.H` constructs + production-metric G-reciprocity
  1.8e-16.  THE FLIP: `_multi_d_multi_moment_reverse_deferred` DELETED
  (predicate + frame-guard arm; `_DAGWavefront.has_transpose_walk` → True
  unconditional — the scheme factor stays on the registration-coupled
  `has_transpose_kernel` covering law), machine-facing contracts
  truth-patched in the SAME commit (the C4-d lesson applied UP FRONT:
  Protocol ×2 + `apply_transpose` + the `is_adjointable` two-factor
  narration — both factors' examples were pre-C2/C4 fossils).  Gates: the
  ONE dense pin extracted to `_assert_dense_mt_pins_object` (three
  duplicated probe bodies collapsed — d=2 DD, d=3, LD-2D share it);
  LD-2D dense-`Mᵀ` (FFW+window; ScanMarch construction-refusal row pins
  the parametrization's honesty); pairing vacuum+reflective; reverse
  window ≡ full at n_face_moments=2; the §6.1 assembled-`Mᵀ` KEYSTONE
  (the moment-generic `assemble_ordinate_blocks` — kernel probing + the
  `octant_moment_frame_signs` conjugation — needed ZERO changes);
  **M-R2c-MOMENTDROP with the EXACT Mode-7 pair** (anisotropic red
  4.4e-1; slope-free EXACTLY unchanged at ~1e-16 — the §6.3 blindness is
  exact, not approximate); **M-R2c-FRAMESIGN-2D as the parity split**
  (s[3]→+1 one-sided: one-backward octant rows deviate 8.66, even-parity
  rows 8.9e-16 IN ONE RUN — the involution's ∏-group theory made
  visible; fwd reference probed clean pre-patch since a both-sides frame
  error conjugates away).  Reciprocity `_BUILDERS["ld_2d_2g"]`
  (reflective nonsquare non-uniform het; the `_bulk_measure` kron loop
  was ALREADY d-generic — `[1,θ,θ,θ²]` fell out; the metric
  cross-check row pins production ⊗-mass ≡ the independent kron at d=2).
  TRACE-metric ruling recorded at the builder: the moment-resolved
  face's metric is the purely-ANGULAR Wave-O `|Ω·n|·w_n` broadcast
  uniform over the moment axis — the θ-mass is a bulk phase-space-measure
  concept (the trace metric carries no spatial factor at all); production
  (`_build_trace_metric_weights`) and `_g_inner` spell the same
  broadcast.  Deferral-file LD-2D rows → the positive surface.
  Batteries: reverse file 31/31, deferral 20/20, reciprocity 32/32,
  sweep-core+operators 1281/0; pyright CLI +0 (the reciprocity 10
  verified pre-existing via copy-first).
- **C5-b ✅ — the doc blast + the changelog truth-repair.**  The rst
  deferral ledger: LD-2D entry RETIRED into the retired-entries note
  (the ledger now holds ONLY R7); the "WHOLE ledger" statement records
  grid-completeness + the still-armed predicate machinery; orientation
  table cell → C3/C4/C5 scheme-complete; the :425-449 contract prose
  truth-patched.  The #280-citation audit (54 production sites): the
  vast majority are TRUE history (kept); genuinely stale = 3 sites, all
  fixed — the Protocol `has_transpose_walk` docstring (multi-D walks
  "are the #280 deferral"), `solver.py` certificate-skip rationale ("LD
  reverse-scan/adjoint arms are themselves #280 deferrals" — the arms
  landed C2/C5; the skip's TRUE reason, the residual-mint moment
  widening, stands alone now), `streaming.py solve_transpose` ("LD /
  multi-D reverse-scans deferred" — LD-1D landed at C2, empirically
  verified this session; multi-D = R7).  Still-true #280 deferral
  citations (sweep_operator R7 rows, generic-inverse `.H` base
  defaults, green_operator's unbuilt multi-D transpose SWEEP) kept
  unchanged — fix falsehoods, don't churn truths.  **history.rst
  truth-repair (trust-git rule)**: 22 stale *(in development)* markers
  for long-merged campaigns repaired to `merged @ <hash>` (10×
  walk-unification → per-commit + `3f0b8c74`; assembly → `b058083e`;
  k-estimator → `a4952c3e`; #226 family → `1729647`; #236 → `607b548`;
  #118 → `15185e5` — every hash git-verified an ancestor), the legend's
  marker class retired, the stale "#236 pending merge" tail note fixed,
  and the #310 campaign entry added at the table head.  NITs weighed:
  kron 2nd spelling PARKED (still exactly 1 production consumer — C5
  added none; the test-side kron is deliberately independent, anti-R1);
  `SLOPE_MOMENT = 1` PARKED-SHARPENED (at d=2 the slope is slots 1,2,3
  = ŷ,x̂,x̂ŷ — a scalar constant is 1-D-only, so the pair-completion NIT
  is even less apt post-C5); the C3 apply-frame glue trigger did NOT
  fire (LD-2D rode the existing interiors, zero new per-direction
  glue).  Sphinx `-E -W` exit 0; audit exit 0 (orphans 0/311, ERR
  69/69, collect 6730); matrix.rst regenerated.  #310 CLOSES via this
  commit's trailer (lands at push).

- **C3 CLOSED ✅ (2026-07-24)** — branch `feature/sn-adjoint-residue-c3`.
  **Full-tree gate: 6442 passed / 0 failed / 18 skipped / 54 xfailed,
  EXIT 0, 50:49** (count exactly accounted: 6427 @ C2-close + 14 C3-a +
  1 C3-b grazing tooth).  NOTE the first tree run (concurrent with the
  C3-b docstring edits) false-redded the two `inspect.getsource` AST
  tripwires with `File "<unknown>", line 1` — attributed (mechanism:
  shifted line windows under a live run), both files 33/33 on the frozen
  tree, rerun EXIT 0 → **lesson L37** (Python sources FREEZE while a
  long gate runs) @ `ce527996`.  Spec §12.2 acceptance walked: object
  gates only (no keff anywhere in the new file) · every committed gate
  class red-verified (pairing + dense under M-R2-ADDRESSING; the partial
  axis-swap; the AST tripwire's red is its identifier-scan construction)
  · Mode-12 both-legs (equivariance no-op + partial-swap red) · all
  references structurally independent (forward probe / CSR / the
  C2-proven scan reverse / hand-built labels) · regimes activate
  (rectangular, non-uniform, het, full composite) · zero flag flips ·
  out-of-scope loud (LD-2D→C5, G-S solve-reverse→R7, tail backstop).
  **C4 MERGED to main (ff, 7 commits `1efb1a61..59830618`) + phase
  branch deleted; user-ratified at the boundary checkpoint ("Merge;
  compact before C5" — the 6458-green + targeted-pin-fix evidence
  accepted as the pre-merge gate; the full-rerun option was offered
  and declined).  ⏸ COMPACTION taken here (user-directed) — **C5
  opens next session, on a fresh `feature/sn-adjoint-residue-c5`
  branch from main: the LD-2D reverse + the #310 doc blast + CLOSE
  #310.**  Read spec §6 (+§3.3(c)) and the C4-a..e status entries
  FIRST.  The heavy lifting is banked: the kernel VJP is d-generic and
  LD REGISTERS it (the walk would already run — C5 is GATES + the
  flip), and the LD-2D deferral has ONE spelling
  (`_multi_d_multi_moment_reverse_deferred`) read by both flip-safety
  faces, so the C5 flip = narrow/delete that predicate and both faces
  move together by construction.  C5 = (a) the spec §6 gates: the
  LD-2D dense-`Mᵀ` + the moment-frame-involution reciprocity row (the
  d=2 cross-moment frame sign — §3.3(c) generalization, the ERR-066
  family's likeliest sign-error site) + the §6.3 Mode-7
  mandatory-anisotropic config (random slope moments — an all-flat
  suite is blind to a mis-signed slope row) + LD-2D `_BUILDERS`
  reciprocity rows (bulk metric = the d=2 moment mass ``[1,θ,θ,θ²]``);
  note the LD-2D reverse rides the WAVEFRONT family ONLY (ScanMarch's
  facewise supports-gate rejects LD-2D in either orientation); (b) THE
  FLIP: retire the predicate + the frame guard arm + rewrite
  `test_ld_2d_wavefront_trait_stays_false` and
  `test_ld_2d_reverse_is_a_typed_deferral` to positives, same commit;
  (c) the doc blast: the `loss_representation.rst` coverage-matrix
  rewrite + re-point the ~40 `#280`-citing deferral comments at #310;
  (d) #310 CLOSES (via trailer; lands at push).  Weigh the banked NITs
  at (c): kron-loop 2nd spelling (3rd consumer?), `SLOPE_MOMENT = 1`
  pair, the C3 apply-frame glue trigger (C5's LD-2D rides the EXISTING
  interiors — NOT the trigger).  After C5: A4 (P1.2–P1.5 via the
  dormant `KEigenvalue` daggered triple) → A5 (φ* carrier — USER rules
  the shape) → A6 + ch15.  L37: Python sources FREEZE during the
  full-tree gate.**
  **C4-d ✅ @ `81aa8105` — enforcer verdict + contract truth-patch.**
  Enforcer: **code + tests ship-quality** — per-axis: twin-paths PASS
  (`_x_scan_faces_transpose` hand-verified as the genuine VJP incl. the
  shift/seed/reversal algebra AND the double seed consumption; the
  ZERO-kernel-x-out claim verified against the forward's underscore
  discards; zero duplicated reverse algebra); orientation-as-object
  PASS; scheme-door PASS (single-sourced, truthful); predicate PASS
  (grep-confirmed exactly 2 consumers); test-quality PASS (teeth
  genuine, d=3 non-vacuous, LD-2D negative kept); rst PASS.  The sole
  blocking cluster = the INVERTED L-004 blast radius: 3 machine-facing
  deferral contracts outside the diff still named the multi-D Cartesian
  reverse as deferred (the `LossRepresentation` Protocol ×2 + the
  public `StreamingOperator.apply_transpose`) — 2 MUST-FIX + 1
  SHOULD-FIX (`sweep_transpose`'s stale mirrors-the-raise analogy), all
  docstring-layer, all fixed + committed with the regenerated V&V
  matrix.  Lesson re-learned: the retirement 3-search must grep the
  DEFERRAL STORY's contract sites (Protocol docstrings, public operator
  docs), not only the rst ledger.
  **C4-e ✅ @ `a6be11c9` — full tree + the one stale pin.**  Full-tree
  serial `-m "not slow"` @ `81aa8105`: **6458 passed / 1 failed**
  (66:37, frozen sources per L37).  The 1 red =
  `test_affine_chain_transpose_single_source` — the #311 pin counting
  `outgoing_face_from_average_transpose` call sites in diamond.py at
  exactly 2; C4-b's `reflect_scan_coefficients_transpose` is a
  legitimate THIRD consumer (the ride-the-primitive discipline the pin
  itself enforces).  Count updated 2→3 (test-only edit; the
  open-coded-twin greps untouched); file 11/11.  Evidence posture: the
  6458 greens stand (no production edit after the run); the fixed file
  green targeted — a full rerun re-proves unchanged results and is
  offered at the boundary.
  **C4-a ✅ @ `1efb1a61` — the WINDOWED production reverse.**
  `MovingFrontierWindow.loss_action_transpose` = the exact windowed
  sibling of the C3 oracle: the UNCHANGED `walk_windowed` over the
  mirror graph × `_CellResidualTranspose` through the shared
  `_OctantWalk` apply-transpose frame — the mirror graph's own
  `window_plan` IS the reversed frontier, so ZERO new frontier code and
  zero forward edits; the `_DAGWavefront` base C4-deferral raise
  retired.  Gates: `test_reverse_window_equals_full` (window ≡ full at
  `array_equal`, bulk + every trace face, vacuum-het-rectangular AND
  reflective nonsquare) + the M-R2-WINDOWDRIFT tooth realized as the
  representable seed-drop (frontier-ORDER unrepresentable — window_plan
  + levels are one graph object, the M-R2-LEVELORDER finding shape).
  Green on contact (34/34); pyright 0/0.  NO flips.
  **C4-b ✅ @ `d634c1e9` — the ScanMarch-2D row-march reverse.**  The
  honest reverse-mode of the row-march's own program; the mirror label
  drives the schedule FREE (the forward's sign-reading spellings yield
  reversed rows + reversed x under mirror signs).  Per row in mirror-y
  order: (1) batched `residual_kernel_batch_transpose` with ZERO
  kernel-x-out cotangent — the forward DISCARDED its kernel `out_x`
  (the scan owns the x-chain), so that pullback vanishes STRUCTURALLY;
  `in_y_bar` = the previous physical row's `out_y` cotangent (backwards
  transverse chaining); (2) NEW `_x_scan_faces_transpose` (scan.py,
  beside its forward): v[i]=in_x_bar[i+1], v[last]=x̄_cap, ONE
  `ordinate_scan_transpose` on the SAME α sequence (transpose of a
  first-order chain = same multiplier, opposite direction); seed cot =
  direct in_x[0] + the returned ψ0 (consumed, not re-derived); (3) the
  β-pullback via NEW scheme door `reflect_scan_coefficients_transpose`
  (ABC raise-default; DD = unit application of the #311 w-generic VJP
  at w=½ — (2,−1) exact, single-sourced, no probing, no inline twin).
  `_OneDimScanWalk`'s 2-D guard reworded structural-1-D-only.  Gates:
  dense-`Mᵀ` + pairing PARAMETRIZED over all three 2-D reps;
  row-march ≈ oracle at rtol 1e-12 (bulk+trace, both configs);
  transverse-chain tooth (zeroed `in_y_bar` → red, Mode 4) + scan-seed
  tooth (dropped ψ0 → red, Mode 5); the deferral file's direct-call row
  FLIPPED positive (capability-without-predicate until C4-c).  41/41 +
  455 sweep-core + 63 operator rows; pyright 0/0 on all six files.
  **C4-c ✅ @ `28f859e8` — THE FLIPS, atomic per §10.**  USER RULING at
  the design checkpoint: scheme-aware family trait INCLUDING
  FullFieldWavefront + a d=3 gate row.  `ScanMarch.has_transpose_walk`
  → True unconditional (supports + the scheme factor narrow);
  `_DAGWavefront.has_transpose_walk` → `not
  _multi_d_multi_moment_reverse_deferred(mesh)` — the NEW module-level
  predicate is the ONE spelling of the LD-2D→C5 deferral, read by BOTH
  flip-safety faces (trait at construction, `_OctantWalk` frame guard
  at apply — R5/R6 discipline).  NEW d=3 row: pairing at machine
  precision + dense-`Mᵀ` matrix equality on a rectangular non-uniform
  2×3×4 `from_axes` spine mesh — **the mirror-octant reverse is
  verified d-generic** (the d≥3 interleave deferral on #310 is
  LD/UBLD-only; DD has no gap).  Deferral file: trait rows True
  (ScanMarch-2D/window/FFW) + NEW LD-2D negative (trait False +
  `is_adjointable` False + `.H` raises — the C5 boundary);
  `is_adjointable`/`.H` rows positive; `test_cart2d_ffw_oracle_state`
  rewritten (trait ⟺ capability — the C3 select-narrow divergence
  CLOSED).  Un-exclusions §5.6: removal-form cart2d bit-identity row;
  phase-C cart2d operator-reciprocity row; G-reciprocity `_BUILDERS`
  cart2d-DD (the 4-face `|Ω·n|·w_n` trace metric live) AND
  `_FULL_LOSS_BUILDERS` cart2d (`A = L+C−S−B` with `.H` composing
  through `OperatorSum` on the multi-D walk; het 2-material asymmetric
  transfer).  rst ledger: multi-D-production-adjoint entry RETIRED
  (recorded in the retired-entries note); 2 stale prose deferral
  mentions truth-patched; orientation-matrix cell updated.  Battery
  107+5xf; pyright CLI 0/0 on production + main gate files (the 41 in
  phase_c/reciprocity PRE-EXISTING — verified identical pre-edit);
  Sphinx `-E -W` exit 0; audit exit 0 (orphans 0, ERR 69/69).
  **C3 MERGED to main (ff, 4 commits `295053b2..d6e125f8`) + phase
  branch deleted; user-ratified at the boundary checkpoint ("Merge;
  compact before C4").  ⏸ COMPACTION taken here (user-directed) — C4
  opens next session, on a fresh `feature/sn-adjoint-residue-c4` branch
  from main: the windowed PRODUCTION reverse.  Read spec §5.1/§5.6 +
  the C3-a status entry FIRST.  The heavy lifting is banked: the mirror
  octant's `window_plan` IS the reversed frontier (built for every graph
  at construction), and `walk_windowed` already accepts the
  `_CellResidualTranspose` level op (annotation widened at C3) — C4 (a)
  adds the windowed transpose interior on `MovingFrontierWindow`
  (mirror of `_loss_action_transpose_interior`, riding
  `_inflow_to_moments`-style plumbing), (b) pins reverse
  `window ≡ full` bit-identical (`test_reverse_window_equals_full`,
  het + asymmetric + rectangular; M-R2-WINDOWDRIFT tooth), (c) the
  ScanMarch-2D reverse (`ordinate_scan_transpose` + transverse
  cotangent chaining) if in scope per spec §5, then (d) THE FLIPS,
  atomic per §10: `has_transpose_walk` True for
  ScanMarch(2-D)/MovingFrontierWindow(/FullFieldWavefront — decide the
  oracle rep's trait at the C4 design step), the deferral-file
  negatives→positives rewrite (incl. `test_cart2d_ffw_oracle_state`),
  cart2d un-exclusion in `test_removal_form_matvec_sweep`, the
  cart2d-DD `.H`-gated G-reciprocity `_BUILDERS` rows + phase-C row,
  and the cart2d baseline already frozen at C3 pins the window
  day-one.  L37 discipline: Python sources FREEZE during the full-tree
  gate.**
  **C3-a ✅ — the multi-D reverse `walk_full` ORACLE, via the MIRROR-OCTANT
  realization.**  DESIGN RULING (the C3 carve's load-bearing find): the
  reverse walk of octant ``o`` IS the forward walk over the mirror graph
  ``graph(−signs_eff)`` — ``face_in(−o) == face_out(o)`` and the mirror's
  levels ARE the forward's reversed, so the spec's four reverse ingredients
  (reversed levels + transposed addressing + boundary in↔out swap +
  reversed frontier-for-C4) all fall out of graph SELECTION;
  ``walk_full``/``walk_windowed``/``_interior_walk`` are UNTOUCHED (zero
  forward-path edits — the strongest possible "forward stays bit-identical").
  Orientation is carried by DATA (`_reverse_octant_traversal`, the multi-D
  sibling of `_reverse_traversal`: mirror AFTER the grazing map, pure-z
  self-mirror, octant order untouched — no inter-octant edge) + the THIRD
  level op `_CellResidualTranspose` (bottoming in the C2 d-generic batch
  VJP; frame conjugation with the PHYSICAL octant's signs — the mirror
  label drives addressing only).  Physics face: discrete μ→−μ is EXACT for
  the DAG addressing; the cell algebra (where μ-reversal is wrong) is the
  kernel VJP.  New surface: `_OctantWalk.loss_action_transpose` (the shared
  apply-transpose frame — trait guard, LD-2D typed deferral →C5, Pattern-4
  tail backstop, the boundary-cotangent algebra mirroring the 1-D reverse:
  out-rows seed ``streamed̄`` masked to out_idx so every forward-discarded
  path's pullback vanishes structurally, ``ψ_out† = −b̄_out``, in-rows =
  identity + walked capture) + `FullFieldWavefront.loss_action_transpose`
  / `_loss_action_transpose_interior` (the forward helpers
  `_octant_face_cochain`/`_edge_outflow` realize the transposed roles
  VERBATIM under the mirror signs).  NO flag flips (predicate False until
  C4; `_DAGWavefront` deferral rewritten to the window-only truth;
  ScanMarch's 1-D-walk guard message re-pointed at C4).
  **Gates (`tests/sn/sweep/core/test_multi_d_reverse_walk.py`, 12, all
  `-O`)**: runtime spy (3 sentinels; forward-control leg — transpose op 0×
  on forward) · AST tripwire (4 frames) · **the NEW 2-D dense-`Mᵀ`
  column-probe** (M probed off the FORWARD apply over the FULL composite
  basis incl. trace; `M_rev == M_fwdᵀ` as a MATRIX, rtol 1e-12, on
  rectangular non-uniform het + anti-vacuous asymmetry check) · Euclidean
  pairing identity ⟨Fx,w⟩=⟨x,Fᵀw⟩ at ~7e-16 (vacuum + reflective configs)
  · **d=1 cross-realization BIT-IDENTICAL** (FFW mirror-DAG reverse ≡
  CumprodScan leg reverse, DD AND LD slab — same kernel, elementwise ops,
  batching-order-free) · assembled-`Mᵀ` per-ordinate CSR blocks (LS4 4×3
  non-uniform, off-block leak 0) · M-R2-ADDRESSING tooth (traversal→
  identity: pairing red O(1); dense gate red-verified at rel 1.0 in a
  design-time one-off) · **axis-equivariance invariance gate** + partial
  axis-swap tooth · LD-2D/solve-transpose/tail-backstop deferral pins.
  **Recorded design findings**: (a) M-R2-LEVELORDER is UNREPRESENTABLE
  (levels + face roles are one graph object); (b) the TOTAL axis
  conjugation measured as an exact no-op on het σ (2e-16) — the reverse
  interior is genuinely d-generic, so the committed axis tooth is the
  PARTIAL swap (Mode-2 face-tuple cross), which cannot even shape-check on
  the rectangular primary configs (the L16 mandate = a shape-guard).
  **Baseline**: `walk_matvec_cart2d_2g.npz` captured POST-verification at
  the REP layer via `_OracleLC` (new-value re-baseline per §5.4; the same
  object C4's window pins bit-identical against); all 4 baseline rows
  green under strict DriftWarning.  Deferral file gains
  `test_cart2d_ffw_oracle_state` (trait False until C4 + direct oracle
  works — the deliberate select-narrow divergence).  Honest scope notes:
  the cart2d G-reciprocity `_BUILDERS` rows ride C4 with `.H`
  availability (C3's trace coverage = the full-composite Euclidean
  pairing + dense-`Mᵀ` trace rows); pure-z 2-D rows exist in no available
  quadrature (the branch is the trivially-self-transposed σ-diagonal,
  mirrored from the forward twin).  Targeted battery 108 passed;
  pyright CLI 0/0 on both touched production files; audit exit 0
  (orphans 0/311, ERR 69/69).
  **C3-b ✅ — enforcer verdict + doc honesty.**  Enforcer: **no
  VIOLATION, no MUST-FIX — "an exemplary carve"** (per-axis: twin-paths
  PASS with the boundary algebra hand-verified as a genuine transpose
  incl. the `+=` shed being NECESSARY — `given[in]` is consumed twice in
  the forward; illegal-states PASS incl. the grazing re-negation verified;
  reads-like-domain PASS; test-quality PASS with the d=1 `array_equal`
  judged principled-not-over-tight; streaming.py production-deferral
  claims verified still TRUE).  3 SHOULD-FIXes, all fixed same-session:
  (1) `_ApplyOperands` docstring generalized to both orientations (the
  driving-field contract); (2) `FullFieldWavefront` docstring now
  enumerates all THREE directions/kernels/level-ops; (3) the
  grazing/pure-z mirror-map corner got its pure-function tooth on
  hand-built labels (no quadrature reaches the branch — the tooth's
  first run caught a wrong EXPECTED value in the test itself; production
  logic confirmed right: mirror of (0,−1) = (−1,+1), physical recovery
  (+1,−1)).  Accepted NIT recorded: the ~15–25-line apply-frame
  scaffolding glue is per-direction BY DESIGN (leaves single-sourced;
  premature unification would parametrize over the boundary-algebra
  difference); collapse trigger = a genuinely-new third apply-direction
  frame (C5's LD-2D is NOT it — it rides the existing interior).
  Theory-page honesty patch (the C5 doc blast's down-payment): the
  deferral ledger's LD-slab entry RETIRED (dead since C2 — it denied the
  moment tail the reverse now carries), the multi-D entry rewritten to
  the oracle-landed truth, the stale "wavefront LD cell kernel is
  1-D-only" reason replaced by the honest LD-2D-reverse C5 deferral, and
  the 2×2 orientation-matrix multi-D cell updated.  Sphinx `-E -W` exit
  0 on the patch; file now 13/13 green.
  **#311 ✅ landed @ `2d226d1c`**: `outgoing_face_from_average_transpose(f̄, w)
  → (f̄/w, −((1−w)/w)·f̄)` minted as the fourth generic reconstruction
  staticmethod; DD's kernel + the two seed-march reversals rerouted
  byte-identically (frozen baselines 0-drift under strict DriftWarning; the
  seed module imports `_DD_W` — the constant's own docstring invites it);
  gates = the VJP pairing identity (Mode-10: dropW/sign/swap all red O(1))
  + the w=½ byte-identity law + the single-source retirement pin.
  **C2-a ✅ @ `d78150cc` — the batch VJP pair + two-arm reverse visit.**
  DESIGN RULING (recorded for the boundary review; forced by the standing
  lenses, mirror-the-forward + no-scheme-branches + build-the-machinery):
  the C1 cell-balance kernel CANNOT serve LD (its `denom` is a DD
  composite; LD needs raw `(g, Σ_t, θ)`), so the reverse visit now mirrors
  the forward `_apply_walk`'s TWO ARMS — Cartesian scheme-uniform through
  the NEW d-generic `residual_kernel_batch_transpose` (the exact VJP
  mirror of `residual_kernel_batch`; C3's `_CellResidualTranspose` bottoms
  in the same kernel), curvilinear keeping the C1 kernel (DD-only, M-M
  thread).  `has_transpose_kernel` re-derived as the COVERING law: batch
  VJP always ∧ (curvilinear kernel iff `supports_curvilinear`).  Evidence:
  frozen fwd+adj baselines bit-exact under strict DriftWarning through the
  arm change (power-of-2-width snapshot data); slab pairing identity
  `⟨A eᵢ, φ⟩ = (Aᵀφ)ᵢ` at 2.5e-16 (the new arm is the EXACT transpose of
  the forward kernel — tighter than the old cell-balance pairing); the
  two-arm sentinel matrix + covering-law tests + batch value tooth.
  **C2-b ✅ @ `b2ed7d9e` — the LD transpose algebra-of-record.**  SymPy:
  `derive_d1_transpose_equals_At_Minv` (VJP ≡ AᵀM⁻¹ from exact Jacobians;
  order discriminant PROVEN nonzero; inflow + face pullback rows) +
  `derive_octant_frame_sign_is_involution` (D²=I, conjugation commutes
  with transpose, d=1 + both d=2 patterns).  Numpy:
  `assemble_inflow_axis_transpose`, `D1ClosedForm.scan_reconstruct_transpose`
  (riding the new `_geom_fold` single source), `_ubld_outgoing_faces_transpose`;
  three bilinear pairing gates (Mode-10 one-offs red O(1)).  NOTE: the
  #311 "LD composes the primitive on first landing" prediction resolved
  differently — the LD batch VJP is TRACE-based (B(+1)ᵀ broadcast), not
  w-blend, so no new w-blend transpose site appeared; the primitive's
  value is the DD sites + pair completeness.  **C2-c ✅ @ `5a283f81` — the LD-slab adjoint; the atomic flip.**
  LD registers the UBLD `AᵀM⁻¹` batch VJP (generated from the SAME
  `_ubld_system` record; mass-inverse FIRST; `B(+1)ᵀ` broadcast;
  `−μB(−1)ᵀM⁻¹` inflow pullback with the d=1 axis-append mirrored) — the
  registration IS the flip via the covering derivation, atomic with: the
  `_run_transpose` LD moment arm (reframe involution →
  `scan_reconstruct_transpose` → the shared w-generic scan chain → the
  DIAGONAL self-transposes of `source_emission` + `scan_slope_face_source`
  → reframe + ×V), the R5/R6 two-guard lift unified onto the trait,
  Pattern-4 tail backstops in BOTH reverse entries (a tail-mismatched
  cotangent raises typed ValueError — the silent-broadcast hazard caught
  LIVE by the first functional probe: a tail-less composite broadcast
  through the batch VJP), **ruling 3 in production**
  (`SNMesh.full_field_space` = `V·w_n ⊗ moment_mass_diagonal` — the new
  scheme surface, = `diag(M)/V`, the same diagonal the kernel normalises
  by), the deferral-file rewrite (negatives→positives; declare→derive),
  and the gate rows: ld_slab_2g G-reciprocity + metric cross-check
  (production θ == independent raw-θ spelling) + the **ghost-metric
  stabiliser proof** (slope-flip O(1)-visible to the θ-metric, EXACTLY
  invisible to the slope-ghost metric — the Mode-12 asymmetry; note: the
  spec's "average-only blind" leg sharpened to the honest ghost/L18 form —
  a θ→1 mis-weighting is caught by the metric CROSS-CHECK, the stabiliser
  blindness needs the 0-weight ghost) + G1/G2 for `ld_slab` (moment-aware
  probe/order/mask, DD byte-identical) + the LD block-upper-triangularity
  certificate (2×2 micro-blocks; §4.3's back-substitution VALUE leg
  covered by G2's full `np.linalg.solve(Mᵀ)` — the block-backsub third
  realization left as an honest possible follow-up) + the dense `AᵀM⁻¹`
  kernel gate with the M-R1b-MASSORDER tooth.  Batteries: 97 + 207 green
  (frozen baselines strict-DriftWarning; MMS-LD forward; C1 sentinels);
  pairing `⟨A eᵢ, φ⟩ = (Aᵀφ)ᵢ` ≈1e-12 on the non-uniform het LD slab;
  pyright CLI 0/0; `-E -W` Sphinx exit 0; audit exit 0 (orphans 0/311,
  ERR 69/69); matrix collect 6668→6688 (+20 foundation).  #311 closes
  with C2-c's commit trailer.
  **C2 CLOSED ✅ (2026-07-24)** — **Enforcer verdict: every code axis
  PASS ("the engineering is exemplary"** — hand-verified VJPs; the
  two-arm reverse RETIRES a genuine pre-C2 twin, the slab reverse's
  cell-balance spelling vs the forward's batch kernel; covering law
  sound; Mode-12 gates with proven teeth); 3 MUST-FIX + 1 SHOULD-FIX =
  docstrings the flip made FALSE (the pre-flip LD-deferral contract), all
  four rewritten @ `26b28e85`.  Enforcer NITs recorded: (a)
  `moment_mass_diagonal`'s kron loop is an accepted 2nd spelling —
  collapse into `_ubld.py` on a THIRD consumer; (b) `SLOPE_MOMENT = 1`
  constant would complete the `AVERAGE_MOMENT` pair (forward + reverse
  together, pre-existing).  **Full-tree gate: 6427 passed / 18 skipped /
  54 xfailed, EXIT 0, 1:12:21** (nohup E-core pacing).  **C2 MERGED to
  main (ff, 8 commits `2d226d1c..26b28e85`) + phase branch deleted;
  user-ratified at the boundary checkpoint.  ⏸ COMPACTION taken here
  (user-directed) — C3 opens next session: read spec §5 (R2a) +
  `sweep_graph.py` (`walk_full` — the full-cochain oracle walk;
  `_CellSolve`/`_CellResidual` level ops) FIRST, on a fresh
  `feature/sn-adjoint-residue-c3` branch from main.  C3's kernel is
  ALREADY BUILT: the d-generic `residual_kernel_batch_transpose` (DD +
  LD) landed at C2 — C3 mints `_CellResidualTranspose` bottoming in it +
  reversed levels + transposed addressing (gather at face_out,
  scatter-accumulate at face_in, boundary in↔out swap) + the NEW 2-D
  forward-apply column-probe oracle + the reverse-walk spy/AST tripwire.
  NO production flag flips (the oracle arm).**

- **C1 ✅ 2026-07-24** — the DD kernel-pair relocation: `streaming_cell_transpose`
  minted on `DiscretizationSchemeBase` (STORAGE-FREE spatial VJP
  `(res_bar, psi_out_bar, denom, abs_mu_A_total, volume) → (psi_bar_cot,
  psi_in_bar)`; default typed raise); DD's override = the relocated hand-code
  VERBATIM (op order preserved; the single scatter-add proven bit-safe by the
  one-visit-per-slot property); **`has_transpose_kernel` now DERIVED in
  `__init_subclass__` from the override** (ruling 2 — declared-True-with-no-kernel
  unrepresentable, pinned by the `_DeclaredLiar` law test; DD/LD declarations
  DELETED, values reproduced by derivation); the walk's reverse visit calls the
  registered kernel, the `numer_bar` accumulation stays with the walk (ruling 1);
  the R5 guard message → derived-trait + #310 C2. **Degenerate-ordinate branch KEPT
  slot-local** (face-free volumetric relation; its |μ|≈6e-17 rows would make a
  kernel routing compute-then-discard a fake face cotangent) — the retirement pin
  allows exactly that one `denom * ob` survivor. NEW gates
  (`tests/sn/operators/test_streaming_cell_transpose_relocation.py`, 9 foundation):
  Mode-11 sentinel (count>0 adjoint / ==0 forward control) · angular-thread
  asymmetry (MorelMontry fires sphere/cyl, silent slab) · M-R1a-VALUE tooth (1.9
  perturbation → O(1) red + bit-identical recompute control) · retirement grep ·
  the derivation law. Gates: targeted battery 109 passed + 1 xfail (frozen
  `walk_matvec_*` baselines bit-exact through the relocation); new file 9/9;
  pyright production floor 1 (CLI; the streamed scheme.py complaint was LSP
  noise); `-E -W` exit 0 zero warnings; matrix regen honest (+9 foundation,
  collect 6659→6668). Out-of-scope survivors noted: `psi_half_angle_seed.py`'s
  own DD-chain transposes (System-B seed machinery, spec-scoped out).
  **Enforcer verdict: COMMIT-READY, 0 MUST-FIX / 0 SHOULD-FIX** — the
  reassociation-safety proof independently re-derived (zero-init +
  one-visit-per-slot ⟹ `fl(a+b)` both ways); the derivation mechanics verified
  against `RegistryMixin` (key= forwarding, dataclass re-run idempotent,
  MRO-correct); pyright scheme.py CLI 0/0; doc-honesty grep clean tree-wide.
  ONE architectural opportunity → **#311 FILED** (the w-generic affine-outflow
  TRANSPOSE primitive; the forward is single-sourced, the adjoint carries 3
  hardcoded `w=½` spellings; C2's `w=1/(1+k)` is the L-002 collapse trigger —
  folded into C2) + 1 NIT (deferral-file "declare"→"derive" messages, rides
  C2's rewrite). **Full-tree gate: 6405 passed / 18 skipped / 54 xfailed; the
  only 2 reds were the STALE `TestKernelSourceOfRecord` sha-pins** — broken by
  the #231 P2-E docstring trim (execution-only pin, first full run since that
  merge), attribution PROVEN by 3-revision hash comparison (pre-docs-merge
  `126c4e40` matches old pins; C1 hashes identical to pre-C1 `ec74be50`);
  re-hashed per the pin's own procedure @ `8b5d34d8` (anchors were green in
  the same run). Run note: a nohup-detached run lands in macOS background QoS
  (E-cores) — ~1.6 h instead of ~52 min; fine for a gate, don't mistake it
  for a hang. **C1 MERGED to main (ff) + phase branch deleted; user-ratified
  at the boundary checkpoint. ⏸ COMPACTION taken here (user-directed) — C2
  opens next session: read the spec §3+§4 + #311 FIRST (the w-generic
  transpose primitive lands before the LD kernel), on a fresh
  `feature/sn-adjoint-residue-c2` branch from main.**
