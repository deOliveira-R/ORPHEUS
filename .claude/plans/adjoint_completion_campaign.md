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
