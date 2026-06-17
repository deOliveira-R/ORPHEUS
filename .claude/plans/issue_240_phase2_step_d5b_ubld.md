# #240 Phase 2 Step D5b — N-dim Linear-Discontinuous on the DAG wavefront (the UBLD)

> **Durable in-repo recovery anchor** (project rule: plans live in ORPHEUS/.claude/).
> Parent: `.claude/plans/issue_240_phase2_step_d_homing.md` (§D5). Subsumes **#158 Increment D / #38**.
> Branch `feature/sn-space-angle-tier2` (off `main`@`cba6d2f`), NOT pushed/merged.
> **STATUS (2026-06-17): D5b-S3 DONE + COMMITTED (`e74eafb` feat / `3e8f101` chore). The unified all-d LD moment matvec + the diffusion-limit closure (ERR-061) + OWED-2 scan + the pure-z gate (ERR-062) are landed; #37 (Inc C) + #38 (Inc D) CLOSED (folded). ⭐ RESUME AT D5b-S4 (#62, strengthened 2-D MMS + the `ld-cartesian-2d` verifies label) and/or D6 (#55, docs expansion — the archivist stubs `two-moment-axes`/`spatial-moment-space`/`ld-ubld-unified-moment-residual`/`ld-ubld-moment-scan`/`ld-ubld-pure-z-collision` await rich-narrative expansion). Branch NOT pushed. ⚠ error_catalog.md (ERR-060/061/062) restored in the working tree but UNCOMMITTED — lands via the instruction-architecture flow; the `catches` markers ARE committed.**
>
> ⭐ **STATUS UPDATE (2026-06-17): D5b-S4 DONE + COMMITTED (`882341d` test / chore next).**
> The strengthened 2-D Cartesian LD stress MMS (`SN2DCartesianLDStressMMSCase` + Branch-1
> SymPy algebra-of-record + 8 `@foundation` derive-tests + the L1 D5b.2/.3/.4 gates, all
> `@verifies("ld-cartesian-2d")`) is landed; the `ld-cartesian-2d` label STUB is minted
> (archivist expands in D6). Flow: method-implementer build → elegance PASS-WITH-NITS + qa
> SUPPORTED-WITH-CONCERNS → 4 nits resolved (parametrize Branch2≡Branch1 over BOTH groups +
> 2 cells; extract `_build_per_cell_hetero_materials` killing the `build_materials` byte-twin
> with the DD het case; sharpen the honest-scope note to vv Mode 10; rename the
> over-promising struct-indep test). Main-agent re-verify GREEN (15 S4 gates + 632
> no-regression + 50 MMS-dir + Sphinx clean; Nexus graph refreshed). ⭐ **HONEST SCOPE (do
> NOT overclaim):** verifies the slope-UNKNOWN half of the LM-1989 trap; the slope-SOURCE
> half is UNVERIFIED — the external Q̂ is not consumed (`_lift_external_source_to_moments`
> zeros the slope rows; `solve_sn_fixed_source` rejects a moment-resolved external source)
> AND the scattering channel `Σ_s·φ̂` EXERCISES-but-does-not-CONSTRAIN the slope-source sign
> (vv **Mode 10** activated-but-unconstrained, qa-discovered: a slope-source sign flip leaves
> the O(h²) order + value band unchanged because `Σ_s·φ̂` is O(h)-small, error enters above
> O(h²)). DEFERRED → **#247** (moment-resolved external source + boundary trace + a
> slope-source-sign-sensitive gate; candidate D5b-S5). vv Mode 10 + qa L-034 added to the
> skills (instruction-architecture flow). ⭐ **RESUME AT D6 (#55):** archivist rich-narrative
> expansion of the committed doc stubs (`ld-cartesian-2d` + S3's `two-moment-axes` /
> `spatial-moment-space` / `ld-ubld-unified-moment-residual` / `ld-ubld-moment-scan` /
> `ld-ubld-pure-z-collision`) + mint the orphan UBLD eq labels + ff-merge readiness for
> `feature/sn-space-angle-tier2`→`main` + apply the Phase-0 hand-off.
>
> ⭐⭐⭐ **THE D5b-S3 ARCHITECTURE DECISION (2026-06-17, settled with the user — do not lose):**
> The S3-A consumer dispatch surfaced a genuine fork: the d=1 LD **matvec** is Schur-reduced to a
> SCALAR residual, but `A=(L+C)−S` (OperatorSum) subtracts a now-`(2,)`-moment `S·ψ` element-wise →
> a φ̂-carrying iterate cannot route through the scalar d=1 matvec. RESOLUTION (the user's framing +
> my analysis): **the scalar d=1 matvec was a FLAT-SOURCE artifact, NOT a d=1 Schur degeneration.**
> The Schur reduction is a *solve* technique (it degenerates naturally on the SWEEP and already
> carries the slope via `schur_xV`'s `s_hat`); a *matvec is a forward apply* — applying the per-cell
> `2^d×2^d` operator to the moment vector is intrinsically moment-valued. Pre-Inc-C, Q̂≡0 left the
> slope ψ̂ globally-uncoupled, so the Krylov unknown was legitimately scalar. **Inc C makes `Σ_s·φ̂`
> couple the slope GLOBALLY in EVERY dimension** → the slope is a genuine global DOF for all d.
> ⇒ **UNIFY the matvec across all d: `(L+C)·ψ⃗` returns the `spatial_basis_per_axis^d`-moment
> residual for ALL d** (DD/Step=1=scalar=byte-identical; LD-1D=2; LD-2D=4). `A=(L+C)−S_full` is then
> ONE honest operator in every dimension. This **FOLDS the former S3-B (d≥2 matvec raise) INTO S3**
> — land complete `(L+C−S_full)` for LD in all d at once. (User Q2 chose "fold all-d".)
> **BRANCH-FREENESS is the acceptance bar (Cardinal Rule 2, user-verified achievable):** the carve
> NET-REMOVES branches. (1) Drop the `len(s_axes) > 1` conjunct in the moment gates
> (`sweep_graph.py:883`/`:929`) → the pure scheme trait `spatial_basis_per_axis > 1` alone (admits
> d=1 LD). (2) Collapse LD's `if len(s_axes)==1/!=1` kernel forks (`linear_discontinuous.py:468/618/656`)
> to ONE d-generic moment path (S1 proved `assemble_ubld`'s d=1 reduction == the 1-D algebra); RETIRE
> the scalar `kernel_rhs` Q̂=0 matvec arm (flat-source artifact). (3) Close the d≥2 `_CellResidual.cell`
> raise (widen `_ApplyOperands.probe`; InvertibleOperator flat-ravel is shape-agnostic). (4) SHAPE-AGNOSTIC
> carriers via the `face_moment_tail`/`spatial_moment_tail` formula (`per_axis==1 → tail==() → scalar →
> byte-identical`), NEVER `if LD-d≥2: widen`. (5) d=1 SWEEP stays CumprodScan (perf via STRATEGY
> selection, not an `if d==1`); `scan_xV` gains `s_hat` single-sourced with `schur_xV` through
> `D1ClosedForm`. (6) ZERO `isinstance(scheme)` (audit confirms none today); curvilinear LD stays a
> DECLARED structural exclusion. Audit confirmed (grep, 2026-06-17): zero scheme-type branches today;
> the only dimension branches are the `len(s_axes)>1` moment conjunct (collapses) + the LD kernel forks
> (collapse) + legitimate geometry-keyed `supports()`/`Compatibility` (stay) + the curvilinear exclusion
> (stays). GATES now ALL in-scope: 4 (DD bit-id 513/1/4) / 1 (both tripwires flip, 2G-het) / 3 (genuine
> Mode-9 SI≡Krylov on `(L+C−S_full)`, aniso+diagonal-cubature+nonzero-c) / 2 (MMS O(h²)) / 5 (1G guard).
> Dispatch brief is the authoritative spec; closeout → `issue_240_d5b_s3_unified_matvec_closeout.md`.
>
> ⭐ **IMPL STATUS (2026-06-17, NOT committed — working tree holds 11 production + 2 test files + doc stubs):**
> The unified-moment-matvec ARCHITECTURE LANDED + branch-free (grep clean: no `len(s_axes)>1`, no
> scheme-isinstance). **GATE 4 byte-id 513/1/4 (main-agent re-run); GATE 3 SI≡Krylov 4.99e-11; matvec≡sweep
> round-trip 1e-16; GATE 2 2-D 3-pass.** ⚠⚠ **BUT the DIFFUSION LIMIT is NOT recovered** — thick-cell LD
> converges to DD only with refinement (main-agent CONFIRMED: nx=4 rel=0.389 / nx=16 0.079 / nx=64 0.008 =
> the flat-source signature Inc C must KILL). The operator is internally self-consistent (matvec≡sweep,
> SI≡Krylov) but is NOT the diffusion-consistent operator. The 1-D tripwire CORRECTLY LEFT xfailed (no false
> green). **numerics-investigator DISPATCHED (background) for OWED-1.** PRIME SUSPECT = the moment-source
> weighting convention: S3-A0 chose `S_full=Σ_s⊗I` (raw) + the agent M-normalized the matvec (`M⁻¹·(L+C)`),
> algebraically self-consistent but maybe NOT diffusion-consistent vs the principled `S_full=Σ_s⊗M` +
> natural residual (the slope rows `M_ii=θ=1/3` may be under-driven 3× while `M_00=1` average agrees → the
> φ̂-TINY smoking gun, ~−0.02 vs φ̄~1.5). SECOND SUSPECT = the d=1 dense `assemble_ubld` slope row (NEVER
> exercised pre-carve; S1's bilinear-exactness oracle NULLS the diffusion regime — Mode-7 at the primitive).
> Investigator EMPOWERED to revisit the Σ_s⊗I decision (DD-no-op at M=1). OWED-2 (the `scan_xV` `s_hat` for
> the d=1 CumprodScan SI path) is coupled, flagged. Archivist DEFERRED until the convention settles (L10).
> Recovery: `issue_240_d5b_s3_unified_matvec_closeout.md` (the M-norm derivation) + the investigator closeout
> `issue_240_d5b_s3_diffusion_limit.md`.
>
> ⭐ **UPDATE 2 (2026-06-17) — OWED-1 FIXED (ERR-061); OWED-2 is a COMMIT BLOCKER (3 red tests), being fixed.**
> numerics-investigator found + fixed OWED-1: **ERR-061** = the per-ordinate LD slope `ψ̂_n` lived in the
> SWEEP frame but `φ̂ = Σ_n w_n ψ̂_n` (the scattering slope source) assumed the GLOBAL-x frame → backward
> ordinates (μ<0) CANCELLED the forward slopes → φ̂ ~6× under-driven (sign-flip #1 + convention-drift #6).
> FIX = a single-sourced `2^d` moment-frame involution `octant_moment_frame_signs(octant_signs, per_axis)`
> via `_reframe` at the `_CellSolve`/`_CellResidual` octant boundary (DD/Step per_axis=1 → None → byte-id).
> Cracked by a from-scratch LM-1989 solver reproducing the WRONG value bit-for-bit (sweep-frame) then the
> RIGHT value (global-frame). **Diffusion limit RECOVERED (main-agent re-run): nx=4 rel 0.389→0.041 / 16
> 0.079→0.002 / 64 0.008→0.000; tripwire flips to genuine PASS.** ⚠⚠ **OWED-2 is NOT a deferrable follow-on
> — it REDS 3 tests** (main-agent re-run, `test_mms_ld_slab.py`): `test_sn_1d_slab_ld_mms_converges_second_order`,
> `test_sn_1d_slab_ld_mms_krylov_matches_si`, `test_ld_two_paths_scan_equals_dag_oracle` — the d=1 LD
> CumprodScan SI/scan path CRASHES on the moment source (`broadcast (N,1,nx,2) vs (1,1,nx)`; a regression,
> the scan was valid pre-S3). FIX (investigator-specified) = `scan_xV` gains `s_hat` single-sourced with
> `schur_xV` through `D1ClosedForm` + apply the SAME `octant_moment_frame_signs((dir_sign,),2)` involution
> on the scan (backward: `-Σ_s·φ̂/W` IN, sign-flip ψ̂ OUT) → scan≡DAG. **method-implementer DISPATCHED for
> OWED-2 (background); NOT committable until the 3 reds are green.** Then: full re-verify (GATE 4 513/1/4 +
> the 3 reds green + diffusion repro + GATE 3) → elegance + qa on the COMPLETE S3 → commit. #246 filed
> (typed SpatialMomentSpace shape-probe predicate, S4 deadline). ERR-061 written to error_catalog.md
> (forbidden-to-commit; land via instruction-architecture flow); `catches("ERR-061")` markers ARE committed.
>
> ⭐ **UPDATE 3 (2026-06-17) — OWED-2 DONE + full re-verify GREEN; reviews IN; one qa BLOCKER being fixed.**
> OWED-2 single-sourced (`D1ClosedForm._slope_fold` shared by matvec + scan; scan consumes the same
> `octant_moment_frame_signs`/`_reframe` involution — Pattern 2). **Main-agent FULL re-verify ALL GREEN:**
> imports OK; branch-grep clean; diffusion nx=4 rel 0.041 BOTH solvers; GATE 4 513/1/4; LD slab MMS 7 passed;
> spatial+operators 575P/7F (the 7 git-stash-confirmed pre-existing curvilinear-SPH). **REVIEWS: elegance
> PASS-WITH-NITS** (zero BLOCK; single-sourcing provable; optional nit = hoist `frame_signs_for` free func).
> **qa SUPPORTED-WITH-CONCERNS** — headline correctness GENUINE (the from-scratch LM-1989 solver reproduces the
> WRONG value bit-for-bit sweep-frame then the RIGHT global-frame, VALUE-anchored to analytical diffusion 2.36;
> `catches("ERR-061")` mutation-verified; `@foundation` independent ground stays green). ⚠ **qa CONCERN A
> (COMMIT BLOCKER, being fixed):** the d≥2 LD **matvec** `pure_z` arm (`loss_representation.py:742`) lacks the
> **sweep**'s moment-broadcast guard (`:654-655`) → 2-D LD Krylov on a PURE-Z-bearing quad (MMS N=110 has 2;
> `level_symmetric` none = why it hid) CRASHES `broadcast (2,6,6) vs (1,2,6,6,4)`. L21/L14 twin-path asymmetry;
> the matvec twin was never gated vs the sweep on a pure-z quad (3rd recurrence of "matvec needs a committed
> gate"). Loud crash, NO false-green. **method-implementer DISPATCHED (background): single-source/port the
> pure_z guard + add a COMMITTED 2-D LD Krylov≡SI gate on a pure-z quad + bundle the elegance hoist.** Concern B
> (non-blocking): 4 runtime-safe pyright nits; the `verifies("ld-cartesian-1d"/"ld-slab")` labels lack math
> blocks (pre-existing, archivist/D6). NOT committable until Concern A green; then re-verify → commit.
> S1 (the unified per-cell UBLD primitive) is complete on both branches:
> **S1 commit chain** (branch `feature/sn-space-angle-tier2`, NOT pushed/merged):
> `7abba5d` (design plan, pre-session) → `cb84b7b` (Branch-1 SymPy algebra-of-record + oracles)
> / `9382ace` (chore) → `69b19c9` (Branch-2 production numpy primitive `orpheus/sn/spatial/_ubld.py`
> + the LD Rule-of-Three collapse: `_schur_terms`/`_kernel_terms`/`affine_scan_coefficients`
> single-source through `d1_closed_form`) / `effdfc1` (chore) → `c8489e4` (status) →
> `3567567` (hindsight `_xV` extraction) / `e080b10` (audit chore).
> **Verification:** bit-identical DD negative control (513 strict-gate pass, no golden `.npy` moved);
> LD principled ~1-ULP re-baseline (no golden moved); the d-generic dense `assemble_ubld`/
> `per_cell_solve` primitive is built + d=2-exact-on-bilinear verified (the S2 d≥2 path) but NOT yet
> wired (the `linear_discontinuous.py` `len(s_axes)!=1` raise STAYS). ERR-060 caught (dropped
> |μ_axis| inflow factor). elegance PASS + qa SUPPORTED both branches.
> **Hindsight audit (post-S1, independent elegance-enforcer): PASS — the code landed in the right
> place.** `d1_closed_form`/`D1ClosedForm` correctly live in `_ubld.py` (co-located with the dense
> reduction they're proven `==` to — NOT moved into `LinearDiscontinuous`); the dense primitive
> ahead of its S2 consumer is the justified reference-oracle exception (it caught ERR-060); the
> sympy↔numpy split is the justified algebra-of-record (do NOT single-source the 1-D factors —
> defeats structural independence; `balance.py` confirmed clean of a 3rd copy). ONE real duplicated
> expression FIXED (`3567567`: `schur_xV`/`scan_xV` shared `g·V`/`V·eff_denom` → `D1ClosedForm._xV`).
> ⚠ The ERR-060 `error_catalog.md` entry is written in the working tree but UNCOMMITTED
> (`.claude/skills/` is in the standing forbidden set — land via the instruction-architecture flow;
> the `catches("ERR-060")` marker IS committed in the test).
>
> ✅ **S2 (#60) DONE + committed (`495af60`; chore records in the follow-up `chore(claude)`).**
> The d≥2 UBLD LD kernel runs on the DAG wavefront: `cell_kernel_batch`/`residual_kernel_batch`
> fork on `len(s_axes)` — d=1 stays on `d1_closed_form`→CumprodScan (no perf regression), d≥2
> assembles `A=G+F_out+Σ_t·M` (Kronecker, `_ubld_system`) + the per-axis upwind inflow
> (`_ubld_inflow`→`assemble_inflow_axis`) and solves the batched `2^d` system (`per_cell_solve`);
> outgoing faces = the tensor-Legendre downstream-node trace (sum of `o_a∈{0,1}` blocks,
> `_ubld_outgoing_faces`). The contract is indexed by ONE trait
> `spatial_basis_per_axis` (DD/Step=1, LD=2; per-cell=`per_axis**d`, per-face=`per_axis**(d-1)`);
> the `is_multi_moment` (per_axis>1) gate widens the interior face cochain
> (`_MovingFrontier._win`, `FullFieldWavefront._octant_face_cochain`) + the moment-reducing emit
> (`_CellSolve`/`_CellResidual`) — DD/Step **byte-identical** (no length-1 axis appended). `Q_cells`
> scalar source is lifted to slot 0 (`AVERAGE_MOMENT`) inside `_ubld_system` (no `_SolveOperands.Q`
> carrier change this step). Single sources minted in `_ubld.py`: `AVERAGE_MOMENT` (slot-0 layout)
> + `face_moment_tail` (the two storage policies). The D5b.0 pin INVERTED
> (`test_cell_kernel_batch_admits_multi_d`).
> **Gates green:** DD bit-id strict (sweep/core+solve+cartesian_2d) **562 / 2 skip / 4 xfail** no
> golden moved; D5b.1 round-trip (full `2^d`, non-flat inflow, het); D5b.4 kernel matvec twin (faces
> both directions); D5b.3 two-paths FFW≡MFW (NON-SQUARE, forced reps); D5b.5 DD≠LD routing-flip; the
> carry-forward d=2 numpy↔symbolic A==A CELL pin + `test_d2_exact_on_bilinear` (the genuine ERR-060
> inflow catcher); 2-D LD MMS O(h²) smoke (vacuum, `@l1` NO `verifies` — Mode-7 honest); LD
> spatial+primitive+symbolic 40; Sphinx clean (matrix regen). **Reviews:** elegance PASS-WITH-NITS
> (orphaned `is_affine_scannable` docstring; `AVERAGE_MOMENT`/`face_moment_tail` single-source; gate
> invariant comment) ALL FIXED; qa SUPPORTED-WITH-CONCERNS (dropped a misattributed
> `catches("ERR-060")` off the cell-A==A pin; `_ubld_inflow` `|None` seed) ALL FIXED.
> ⚠ The method-implementer DIED on an API-overload AFTER completing the impl+tests+crosswalk+docs
> (NOT mid-logic); the main agent triaged the tree (imports+gates green), ran the reviews, fixed all
> nits, and committed. The crosswalk `.claude/plans/issue_240_d5b_s2_crosswalk.md` is the design
> record (moment ordering `[bar,ŷ,x̂,x̂y]`, the Q lift, the domain-edge minimal touch).
>
> ⭐ **S3 = fold Increment C (the φ̂ spatial-moment iterate; closes the diffusion limit d=1 #37 AND
> d≥2). PROACTIVE test-architect + explorer dispatches DONE.** Authoritative design record =
> `.claude/plans/issue_240_d5b_s3_crosswalk.md` (READ FIRST). ⭐ **FP RESOLUTION (load-bearing):**
> Inc C is a PHYSICS-COMPLETION, NOT iteration-only — S2 converges to `(L+C−S_flat)` (diffusion-limit-
> INCONSISTENT), S3 to `(L+C−S_full)` (CONSISTENT, `S_full=Σ_s⊗I_spatial`). The FP CHANGES (that IS
> the point). **DO NOT gate S3-FP==S2-FP (a Mode-9 mis-application — REJECTED).** The genuine Mode-9
> invariant = SI-with-lag ≡ Krylov on the SAME `(L+C−S_full)`. S3 SPLIT along sweep/matvec (L14):
>
> **S3-A0 (#64) ✅ DONE + committed (`d313d16` feat / `96dfc96` chore) — the typed-field-space
> foundation.** User chose option (b): a first-class typed space factor. Minted
> `SpatialMomentSpace(FunctionSpace)` (peer of
> `SphericalHarmonicSpace`; the within-cell tensor-Legendre DG basis, size `per_axis^d`); the flux +
> source-sink factories OPTIONALLY compose it via `*` (the `_compose_spatial_moments` single-source,
> "append iff >1"), DEFAULT-OFF — CONSTRUCT GENERAL, no production field carries the axis yet →
> byte-identical (DD strict 513/1/4, mutation-verified teeth both ways). Also closed #207
> (`TensorProductSpace.find_factor`, now generic). The Σ_s⊗I scattering einsum lift
> (`material_xs_field.py`, the prior dispatch's landed half) folded into this commit. elegance
> PASS-WITH-NITS + qa SUPPORTED, all nits fixed (find_factor generic; from_mesh_and_L routed through
> the single-source). ⚠ #245 filed (relocate `AVERAGE_MOMENT`/`face_moment_tail` `_ubld`→`numerics` —
> the numerics→sn UP-import smell; deferred per Pattern 6).
>
> **S3-A (#61) RESUME HERE — the φ̂-iterate consumer (forward/SI).** The SpatialMomentSpace primitive
> is BUILT (S3-A0); S3-A SELECTS it. READ the S3-A0 API + the prior in-flight analysis:
> `.claude/agent-memory/method-implementer/issue_240_d5b_s3_a0_spatial_moment_space_closeout.md`
> (the `spatial_moments` factory param + `_compose_spatial_moments` + `find_factor`) +
> `.claude/agent-memory/method-implementer/issue_240_d5b_s3_a_inc_c_closeout.md` (the §BLOCKER now
> resolved by S3-A0). SELECT the SpatialMomentSpace factor for
> the LD iterate: `_CellSolve.cell` accumulate `φ̂` (stop dropping at `psi_avg[..., AVERAGE_MOMENT]`);
> the between-sweep iterate carriers (un-windowed `AngularFlux` + windowed `HarmonicMomentField`) gain
> the axis (via the S3-A0 `spatial_moments` factory param, selected from
> `scheme.spatial_basis_per_axis**ndim`); the two source seams (d≥2 `_ubld_system` genuine `(2^d,ng)`
> Q; d=1 scan `D1ClosedForm.kernel_rhs`/`schur_xV` Q̂ — currently hard-coded 0); `ScatteringOperator.apply`
> producer is ALREADY moment-agnostic (S3-A0). Flips BOTH thick-diffusive tripwires (1-D #37 krylov-slab;
> 2-D via SI). Gates GATE1 (both tripwires, 2G-het MANDATORY)/GATE2 (transport MMS O(h²))/GATE4 (DD
> bit-id)/GATE5 (1G guard). **S3-B (#63)** = close the d≥2 Krylov matvec raise (`_CellResidual.cell` —
> the `2^d` probe + flat ravel); GATE3 (SI≡Krylov on `(L+C−S_full)`, aniso+diagonal-cubature+nonzero-c).
> **S4 (#62)** = the strengthened Mode-7 MMS + the non-vanishing-boundary `BoundaryFlux` trace +
> `@verifies("ld-cartesian-2d")`.
>
> **DOCS (user-requested) ✅ DONE (`a2dbb39`):** the angular-vs-spatial moments discussion
> (`two-moment-axes`) + both D5b-S3 stubs (`spatial-moment-space`, `ld-ubld-scattering-moment-lift`)
> expanded to full narrative; Sphinx clean. **S2 carry-forwards STILL OPEN for D6:** `ld_ubld.py` oracle
> line-number→symbol citations; mint/link the UBLD orphan equation family.
> The design pass (plan mode) resolved the two open architecture questions WITH THE USER and
> EXPANDED the scope — see "⭐⭐ THE TWO DECISIONS" below. The approved implementation plan is
> `.claude/plans/mellow-swinging-breeze.md` (the S1–S5 sub-step sequence + the unified
> `n_spatial_moments` lattice). Tasks: #58 (umbrella, in_progress) → #59 (S1, in_progress) →
> #60 (S2) → #61 (S3, the Inc C fold) → #62 (S4, MMS) → #55 (D6). D5b now subsumes BOTH
> #37 (Inc C, FOLDED IN) and #38 (Inc D, = this).
> D5-0 + D5a DONE (the routing substrate is honest). The D5 design pass (test-architect spec +
> literature formulation + cross-domain-attacker frame check) is DONE — the three memos below
> are the authoritative inputs.

## ⭐⭐ THE TWO DECISIONS (design pass 2026-06-16 — do not lose across compaction)

**Decision 1 — FOLD Increment C (#37) into D5b.** D5b ships the d-generic *global
moment-source iterate* (the scattering-slope feedback `φ̂`), closing the diffusion limit for
BOTH d=1 (#37) AND d ≥ 2 in ONE arc. (User chose "Fold Inc C into D5b" over "defer".)

**Decision 2 — UNIFY the per-cell algebra WITHOUT a 1-D perf regression** (user: "why can't we
have unification while keeping performance?"). Single-source the UBLD per-cell algebra as ONE
d-generic **primitive** (Kronecker-assembled M/G/F + moment elimination); *derive* the d=1
`affine_scan_coefficients` from it (the closed-form `2×2` Schur = the affine scan recurrence) so
d=1 keeps riding **CumprodScan** (cumprod-vectorized; zero 1-D regression; L16-safe); d ≥ 2
rides the wavefront via a **batched `2^d×2^d` solve** (`np.linalg.solve` over the anti-diagonal
cell stack — vectorized across cells, NOT a per-cell Python loop). Do NOT retire the scan path —
re-derive it from the primitive.

**The unifying architecture — ONE `n_spatial_moments` reduction lattice.** A class-level scheme
trait `n_spatial_moments: ClassVar[int]` (DD/Step = 1, LD-d = `2^d`) indexes the whole contract:
cell unknown (`ψ̄` → `2^d`), face payload (scalar → trailing `2^{d-1}` moment axis), iterate flux
(`φ̄` → `2^d` spatial moments), scattering (`Σ_s` → `Σ_s ⊗ I_{spatial-moment}`), source build
(`(ng,)` → `(2^d, ng)`). DD/Step at `n=1` are the trivial reduction — the carve MUST keep them
**bit-identical** by construction (the backward-compat invariant gating every sub-step).

**Two structural findings that de-risked the carve (the two explorer passes):**
- `WavefrontFlux`/`InteriorFaceSpace` were RETIRED at S6.4(f). The interior face cochain is now
  raw per-axis ndarray tuples, and the DAG walk (`SweepDependencyGraph.walk_full`/`walk_windowed`)
  is ALREADY moment-agnostic and ALREADY exercises the per-axis multi-face collection (DD 2-D
  drives it). ⇒ NO new face type — append a trailing moment axis to the per-axis ndarrays; the
  walk plumbing carries it untouched; the KERNEL owns the `2^d×2^d` solve.
- Today ONLY `φ̄` exists anywhere (iterate, `ScatteringOperator`, production cell kernel are all
  single-spatial-moment; the LD `(2,ng)` two-moment path is reachable only from a UNIT TEST). ⇒
  Inc C *introduces* a spatial-slope contract; it does not extend one that already flows. Minimal
  sites: scalar-flux reconstruction (`_CellSolve.cell`, accumulate `φ̂`); the iterate flux field
  (new `(n_spatial_moments, ng, *spatial)` slot); `ScatteringOperator.apply` (lift over the
  spatial-moment axis); `Q_cells` (grow a slope-source row).

## ⭐⭐⭐ THE LOAD-BEARING FINDING (do not lose this across compaction)

The plan's original framing — "LD's multi-D closure is bilinear = average + **one slope per axis**
(`1+d` moments)" — is **WRONG for Cartesian cells**. The literature pass (Adams 2001, NSE 137;
Maginot-Ragusa-Morel 2016 = **MRM-2016**, NSE 185; Börgers-Larsen-Adams 1992, JCP 98) established a
**design fork**:

* **Simplex-P1 LD** (`{1, x, y}`, `1+d` moments) — the literal "one slope per axis." On a
  **quadrilateral/Cartesian cell it FAILS the thick-diffusion limit** (Adams 2001). It is correct
  ONLY on unstructured triangle/tetra meshes (WMMP-2001).
* **Bilinear/trilinear DG-P1 = UBLD** (`{1, x, y, xy}` in 2-D, `{1,x,y,z,xy,xz,yz,xyz}` in 3-D —
  the **tensor product of 1-D linear Lagrange functions**, `2^d` moments) — **PRESERVES the
  thick-diffusion limit on quadrilaterals**. The `xy` **cross moment** is exactly what simplex-P1
  lacks and is the reason bilinear survives.

**⇒ D5b builds the tensor-product UBLD (`2^d` moments), NOT the simplex `1+d`.** Building the naive
`1+d` object would ship a SILENT physics bug (a scheme that fails the diffusion limit but passes
every non-diffusion test). The `ld-thick-diffusive` tripwire is the catcher: it flips xfail→PASS on
UBLD, stays FAIL on simplex (correct physics — document the basis-dependence).

**Per (group, ordinate): 4 moments in 2-D, 8 in 3-D.** Moment basis = tensor-product Legendre
`{ψ̄, ψ̂_x, ψ̂_y, ψ̂_xy}` (2-D). The current 1-D LD `{ψ̄, ψ̂}` (θ=1/3 Legendre) is the d=1 reduction.

## Authoritative inputs (READ FIRST — these three memos are the contract)

1. **Formulation** — `.claude/agent-memory/literature-researcher/multi_d_ld_closure.md` (the per-cell
   UBLD system, the Kronecker-assembly recommendation, the face-contract widening, the lumping +
   diffusion-limit verdicts, the citations). **THE implementer's contract for the math.**
2. **Verification** — `.claude/agent-memory/test-architect/d5_nd_polymorphism_verification.md` §"D5b"
   (the gates D5b.0–D5b.5) + §3 (the MMS design) + §1 (the coverage matrix) + §4 (the sequence).
3. **Frame check** — `.claude/agent-memory/cross-domain-attacker/d5_trait_and_mms_frames.md` Frame 2
   (the MMS activation analysis + **the slope-driver strengthening** — see §MMS below).

## The math — the per-cell UBLD system (MRM-2016 §2, Eqs. 1–12)

Galerkin weak form of `Ω·∇ψ + Σ_t ψ = S` with upwind flux on the INCOMING faces, per cell K:

```
(12)   A_UBLD · ψ⃗ = R⃗            A = G + F_out + Σ_t·M       (2^d × 2^d, non-symmetric, DENSE)
```
- `M_ij = ∫_K B_i B_j` — mass (collision), tensor-product Legendre → diagonal weights `θ^{#active axes}`
  (`1` for ψ̄, `θ` for ψ̂_x/ψ̂_y, `θ²` for ψ̂_xy).
- `G_ij = ∫_K B_i (Ω·∇B_j)` — gradient/streaming (couples all `2^d` moments).
- `F` — surface: OUTFLOW (`Ω·n>0`, implicit, into `A`) + INFLOW (`Ω·n<0`, **upwind**: upstream
  neighbour's outflow trace → into `R`).
- `R⃗ = M·S⃗_moments + F_in·ψ_in_traces` (the source moments + the upwinded inflow-face traces).

Moment-balance reading (ORPHEUS notation, `s_a=|μ_a|/Δ_a` raw `g`):
- **0th (ψ̄)** = cell balance — the multi-axis generalization of the 1-D `|μ|(ψ_out−ψ_in)+Σ_t h ψ̄=Q̄h`.
- **1st (ψ̂_x, ψ̂_y)** = per-axis slope rows — each the 1-D slope row ON THAT AXIS,
  **PLUS coupling to ψ̂_xy through the transverse face integral** (the multi-D content).
- **cross (ψ̂_xy)** = the bilinear row — NO 1-D analog; makes the system 4×4 not 3×3; the
  diffusion-limit-load-bearing term.

## ⭐ The architecture — the elegant build (Kronecker assembly) + the contract widening

**The build (single-source, Pattern 5 "build the primitive not the product"):** assemble `M, G, F`
as **Kronecker products of the verified 1-D LD operators** (which ORPHEUS already has). The 2-D mass
matrix is `M_1D ⊗ M_1D`, etc. This (a) single-sources the 1-D math, (b) makes the **d=1 reduction a
Kronecker-with-empty-axes identity** (the existing `linear_discontinuous.py` IS the d=1 case),
(c) the diffusion-limit-load-bearing `xy` coupling falls out of the linear algebra rather than
hand-transcription. **Do NOT hand-transcribe the 4×4 / 8×8 entries.** SymPy-verify the assembly
against (i) the 1-D per-axis reduction and (ii) the "exact on a bilinear flux `ψ=a+bx+cy+dxy`" oracle
(the multi-D analog of the 1-D "exact on linear-in-x").

**⚠ THE CONTRACT WIDENING (the architecture-review item — Cardinal Rule 2):**
- The per-cell unknown is a **`2^d`-vector**, not a scalar `cell_average_flux`. `CellResult` carries
  `cell_average_flux` + ONE `outgoing_spatial_flux` today — no slot for the `2^d−1` slope/cross moments.
- Each downstream face is a **`2^{d-1}`-moment object** (average + transverse slope), NOT a scalar.
  In 2-D the downstream-x face carries a linear-in-y profile
  `ψ_out^x(y) = (ψ̄+ψ̂_x) + (ψ̂_y+ψ̂_xy)·P₁(y)`. The interface-flux cochain (`WavefrontFlux` /
  the moving-frontier face storage) must carry per face a `2^{d-1}`-moment object.
- This is a genuine widening over the DD face (scalar) AND the current 1-D LD face (scalar, because
  the 1-D downstream face is 0-dimensional). **The face-cochain type is the load-bearing extension —
  it is the main design-pass decision.** The 1-D Schur-scalar collapse likely does NOT survive to
  multi-D.

**Where it plugs in (LD rides the wavefront, NOT the scan-march):**
- Close the `linear_discontinuous.py:430-436` raise (`len(s_axes)!=1`) by supplying the multi-axis
  bilinear kernel in `cell_kernel_batch`/`residual_kernel_batch`.
- LD's 2-D path is `FullFieldWavefront`/`MovingFrontierWindow` — ALREADY pure scheme-delegators
  (`sweep_graph.py:849/890` → `scheme.cell_kernel_batch`/`residual_kernel_batch`). So closing the
  kernel raise + widening the face cochain makes LD run N-D on the wavefront with (ideally) no
  loss-rep change.
- LD stays the **structural EXCL on the scan-march** (`transverse_coupling_is_facewise=False`,
  landed in D5-0): the `affine_scan_coefficients` raise on non-neutral curvature is NOT inverted —
  only the DAG-kernel `len(s_axes)!=1` raise inverts.

**Lumping:** ship **UBLD (unlumped)** first (most accurate; the object the diffusion-limit proofs are
stated for; clean linear `Aψ=R` solve). FLBLD/SCB (Wareing 2001 / Adams) = a later robustness option
for void/grazing (never default — large numerical diffusion). **`is_positivity_preserving=False`
STAYS** — NO linear bilinear DFE (UBLD or FLBLD) is strictly positive; strict positivity needs a
NON-linear scheme (defer, like the 1-D positivity fix-up).

## Current code surface (file:line — verify at pickup; D-series moved things)

- `orpheus/sn/spatial/linear_discontinuous.py`: `_kernel_terms` (`:411`, **the raise `:430-436`**),
  `cell_kernel_batch` (`:450`), `residual_kernel_batch` (`:475`), `_ld_source_moments` (`:200`,
  requires `(2,ng)` → generalize to `(2^d, ng)`), `_LDCellTerms` (`:223`), `_schur_terms` (`:295`),
  `theta=1/3` (`:286`). `is_affine_scannable=True` (`:278`), `transverse_coupling_is_facewise`
  inherits `False` (the D5-0 EXCL).
- `orpheus/sn/spatial/scheme.py`: the `cell_kernel_batch`/`residual_kernel_batch` contract; `CellResult`;
  the reconstruction staticmethods (`source_emission`/`cell_average`/`outgoing_face_from_average`).
- `orpheus/sn/sweep_graph.py:849/890`: the scheme-delegator call sites (the walk passes the cell to
  `scheme.cell_kernel_batch`/`residual_kernel_batch`); the `WavefrontFlux` / `_MovingFrontier` face
  storage = the contract-widening target.
- `orpheus/sn/spatial/diamond.py:293-401`: DD's d-generic kernel (the explicit `zip(s_axes,psi_in)`
  left-fold) — the STRUCTURAL shape model, though LD's is a coupled `2^d` Schur solve, not a left-fold.
- The MMS home: `derivations/continuous/mms/sn.py` (`SN2DCartesianMMSCase` `:642` is the mirror;
  `_default_hetero_2d_xs_functions` `:950`); `tests/sn/verification/mms/test_mms_ld_slab.py` (1-D, the
  pattern); NEW `tests/sn/verification/mms/test_mms_ld_2d.py`.

## Verification — the D5b gates (test-architect spec §2-D5b; new capability, MMS-gated, NOT bit-id)

- **D5b.0 — INVERT the raise pin.** `test_linear_discontinuous.py::TestLDKernel::
  test_cell_kernel_batch_rejects_multi_d` (`:361-367`) asserts the `d=1` raise → REPLACE with
  `test_cell_kernel_batch_admits_multi_d` (the `(1+2)×…`/`2^d` system runs, returns
  `(psi_avg, (out_x, out_y))`, `len(faces)==2`). A false gate if left.
- **D5b.1 — round-trip (foundation).** `residual_kernel_batch` at the `psi_bar` `cell_kernel_batch`
  solves for vanishes. NON-flat per-axis `psi_in` (the bilinear slope per axis ACTIVE); feed back the
  FULL solved `2^d` state (partial feed passes spuriously). `n_groups∈{1,2}` het. `atol≈1e-12` (the
  `2^d` solve is a few division ULP deeper than DD's `1e-13` — document the bound). Catches #4
  slope-elimination Schur index drift, #2 x-slope↔y-slope swap.
- **D5b.2 — multi-D LD MMS O(h²) (L1, the structurally-independent ref).** NEW
  `test_mms_ld_2d.py`. Driver `solve_sn_fixed_source(..., Mesh2D(...), scheme=LinearDiscontinuous())`
  (`scheme=` threads end-to-end, `solver.py:1972`→`_as_sn_mesh`); the 2-D LD mesh routes to
  `FullFieldWavefront`. Ladder `nx=ny∈{20,40,80}` (or non-square) → `orders[-1]>1.95 ∧ all>1.85` +
  the VALUE band vs `phi_exact` (rate≠correctness). `@l1 @slow @verifies("ld-cartesian-2d",
  "transport-cartesian-2d")` — **`ld-cartesian-2d` is a NEW label D6 mints.** The MMS = §MMS below.
- **D5b.3 — two-paths FFW≡MFW (foundation, the headline invariant).** Same LD via two DAG schedules
  (both pure delegators → same `cell_kernel_batch`; agreement proves the walk storage doesn't perturb
  LD's cell math). Mode-9 config (the §MMS stress ansatz: het, 2G asym, μ-non-trivial, NON-SQUARE
  `nx≠ny`). `rtol≈1e-9` (SAFETY×conv_tol) + a single-sweep `assert_array_almost_equal_nulp`. ⚠ VERIFY
  each leg ran its rep (`default_for` is the expected wavefront; force MFW for the oracle leg) — a
  two-paths gate where both legs secretly ran the same rep is a silent false green.
- **D5b.4 — Krylov≡SI matvec twin (foundation).** LD forward matvec (`residual_kernel_batch`) ≡ LD SI
  sweep (`cell_kernel_batch`) on the 2-D wavefront. ⚠ **HAZARD:** the matvec reconstructs the outgoing
  face — the multi-D LD must supply a WORKING `outgoing_face_from_average` in BOTH directions per axis.
  **A sweep without a verified matvec is half a carve (L14 leg-3 standoff).** `rtol=1e-9,atol=1e-11`.
- **D5b.5 — DD≠LD routing-flip (foundation, the cross-thread guard).** Drive `solve_sn_fixed_source`
  twice (same mesh/materials/source), `scheme=DiamondDifference()` vs `LinearDiscontinuous()`; assert
  the converged fluxes DIFFER at coarse mesh (`not np.allclose(phi_dd, phi_ld, rtol=1e-3)`). **This is
  the catcher that LD is GENUINELY computing LD, not DD** — if the D5-0 misroute regressed, this sees
  DD≡LD and FAILS. (No `verifies` — a "these are different schemes" contract.)
- **Coverage-matrix delta (definition-of-done):** LD `FullFieldWavefront d=2` + `MovingFrontierWindow
  d=2` flip `NOT-YET → T`. **d=3 LD stays `NOT-YET`** (out of D5b scope — flag/track, #227-adjacent).
- The `ld-thick-diffusive` tripwire (1-D today, `test_mms_ld_slab.py:227`): a multi-D analog flips
  xfail→PASS on UBLD (stays FAIL on simplex — the basis-choice catcher).

## The multi-D LD MMS design (vv Mode-7 override MANDATORY — the hardest gate to get right)

**FORBID the existing isotropic 2-D MMS** (`SN2DCartesianMMSCase`, `phi=sin·sin`): it is angularly
flat → NULLS the per-axis slope rows → tests NOTHING about the bilinear coupling (the exact thing D5b
introduces). AI has no SymPy derivation cost → pick the stress ansatz.

**Trial (per ordinate n, group g):** `ψ_{n,g}(x,y) = [ A_g + μ_{x,n}·B_g + μ_{y,n}·C_g ] / W`
with mixed-scale spatial drivers. The test-architect base ansatz:
```
A_g = a0 + a1·sin(πx/Lx)sin(πy/Ly) + a2·cos(2πx/Lx)cos(3πy/Ly)
B_g = b0 + b1·sin(2πx/Lx)sin(πy/Ly)        # x-slope driver
C_g = c0 + c1·sin(πx/Lx)sin(2πy/Ly)        # y-slope driver
```
**⭐ cross-domain-attacker STRENGTHENING (Frame 2 — MANDATORY, the base ansatz is INSUFFICIENT):**
`B` and `C` as written are **x↔y reflections of each other** → a same-sign slope-row transcription
bug (the LIKELY bug, shared code path) **cancels in the measured flux → false green.** Fix: add
x↔y-**broken** cross-harmonics to the SLOPE drivers `B` and `C` THEMSELVES (not only `A`), with `B`
and `C` **NOT reflection-related**, e.g. `B += b2·cos(πx)cos(2πy)`, `C += c2·sin(3πx)cos(πy)`. KEEP
`a0>0` (boundary non-vanishing — the BC closure must satisfy non-trivial interior), 2G asymmetric
downscatter, NON-SQUARE domain. (A first test: weak vs strengthened ansatz against a deliberately
same-sign-flipped LD slope row — weak still converges O(h²) = false green; strengthened breaks it.)

**Activated/nulled (declare in the test, Mode 7):** ACTIVATES streaming(both axes) + collision +
both per-axis slopes + the cross moment + BC + group coupling; NULLS nothing.

**Derivation (Branch 1, L11):** a NEW `SN2DCartesianLDStressMMSCase` in `derivations/continuous/mms/
sn.py` — substitute into `μ_x ∂_x ψ + μ_y ∂_y ψ + Σ_t ψ = (1/W)(Σ_s^T φ + Q^ext)`, solve
SYMBOLICALLY (SymPy) for `Q^ext_n`. `φ_g = ∫ψ dμ = A_g` (the μ_x·B+μ_y·C terms integrate to 0 over a
symmetric quadrature); the streaming derivative carries the full μ-weighted ansatz. Pin the symbolic
source with a `@foundation` derive-test BEFORE consuming it. Structurally independent of the LD kernel.

**⚠ The Q̂≠0 slope-SOURCE posture depends on whether D5b threads the MOMENT SOURCE** (see open question
below): if D5b is flat-source (Q̂=0, like 1-D Increment B), the MMS verifies the slope-UNKNOWN sign
only (document "slope-SOURCE deferred to the moment-source increment"); if D5b threads the moment
source, the MMS MUST supply Q̂≠0 per axis (the strengthened B/C make this natural — the MMS is
slope-source-READY). **The method-implementer declares the posture; the MMS Q̂ follows it.**

## Implementation sequence (test-architect §4, adapted — D5-0/D5a already DONE)

`[Phase 0 lit-formulation ✅ DONE] ∥ [Phase 1 cross-domain MMS frame ✅ DONE]` → D5-0 ✅ → D5a ✅ →
**D5b (here)** → D6. The two proactive dispatches that FEED D5b are already done (the memos above).
At D5b pickup:

0. **DESIGN PASS FIRST (the contract widening).** This is the architecture decision the user must
   shape: the `CellResult`/face-cochain (`WavefrontFlux`) type extension (scalar → `2^d` cell /
   `2^{d-1}` face moments); whether the 1-D scalar path unifies with or stays parallel to the multi-D
   moment path; where the `2^d×2^d` per-cell solve lives. Likely a `cross-domain-attacker` or
   `Plan`-agent pass + user sign-off. ENTER PLAN MODE.
1. **The Kronecker M/G/F assembly + the per-cell `2^d` solve** (method-implementer, against the
   literature contract; SymPy-verify the d=1 reduction + the bilinear-exactness oracle).
2. **Widen the face cochain + wire `cell_kernel_batch`/`residual_kernel_batch` for `len(s_axes)>1`;
   close the raise; INVERT the D5b.0 pin.**
3. **The MMS module** (SymPy `SN2DCartesianLDStressMMSCase`, strengthened ansatz) + the D5b.1–D5b.5
   gates. elegance-enforcer + qa per sub-step.
4. **(deferred)** FLBLD/SCB robustness option (a follow-up issue); d=3 LD (a follow-up issue).
5. **D6** (archivist): mint `ld-cartesian-2d` + link the orphan labels (`loss-rep-facewise-separable`,
   `loss-rep-scanmarch-solve-affine`, `loss-rep-scanmarch-apply-residual`, `ld-cartesian-1d`,
   `ld-slab`); expand the LD theory `.. todo:` stub (`discrete_ordinates.rst:1498`) into the rich
   UBLD derivation + the Step/DD/LD↔advection narrative.

## DESIGN QUESTIONS — RESOLVED (design pass 2026-06-16)

1. **The face-cochain type — DISSOLVED by the code's current state.** `WavefrontFlux`/
   `InteriorFaceSpace` are RETIRED; the face cochain is raw per-axis ndarrays and the walk is
   already moment-agnostic. RESOLUTION: append a trailing `2^{d-1}` moment axis to the per-axis
   face payload; the kernel owns the `2^d` solve; the 1-D scalar path is the d=1 (trailing-dim-1)
   reduction. NO new type. (See "⭐⭐ THE TWO DECISIONS" + the `n_spatial_moments` lattice above.)
2. **D5b vs Increment C — USER CHOSE TO FOLD INC C IN.** D5b ships the d-generic moment-source
   iterate (`φ̂`); the multi-D thick-diffusive tripwire PASSES within D5b, and the same contract
   closes 1-D #37. Sequenced as S2 (flat scattering, external/MMS moment source threaded) → S3
   (the scattering-slope feedback `Σ_s·φ̂`). So #37 does NOT block S1/S2; it IS S3.
3. **UBLD-only vs UBLD+FLBLD:** ship unlumped UBLD first (decided); FLBLD/SCB = a tracked follow-up.
4. **d=2 only vs d=2+d=3:** D5b targets d=2 (the matrix is 4×4, the MMS is 2-D). d=3 (8×8, trilinear)
   is deferred (no 3-D quadrature production path; #227-adjacent). The Kronecker assembly is d-generic
   by construction, so d=3 is "widen the assembly + a 3-D MMS" later.

## Citations (algebra-of-record)
- **MRM-2016** Maginot, Ragusa & Morel (2016), NSE 185(1):17-42, `10.13182/nse16-38` — the UBLD/FLBLD
  derivation, Eqs. (1)-(12), basis (8a-d). OSTI preprint `osti.gov/pages/servlets/purl/1343843`.
- **Adams-2001** Adams (2001), NSE 137(3):298-333, `10.13182/nse00-41` — the thick-diffusion-limit
  verdict (simplex FAILS / bilinear PASSES on quadrilaterals). THE `ld-thick-diffusive` reference.
- **BLA-1992** Börgers, Larsen & Adams (1992), JCP 98(2):285-300, `10.1016/0021-9991(92)90143-m` —
  the original 2-D LD asymptotic diffusion-limit analysis.
- **WMMP-2001** Wareing, McGhee, Morel & Pautz (2001), NSE 138(3):256-268, `10.13182/nse138-256` —
  the 3-D LLD/FLBLD on tetra (the simplex object; the lumping recipe).
- Local (d=1 reduction oracles): Larsen-Morel 1989 (JCP 83), Larsen-Morel-Miller 1987 (JCP 69) in
  `scratch/literature/`.

## Recovery
Read this file + the three memos (top of "Authoritative inputs") + the Step-D homing plan §D5. The
D5-0/D5a commit chain: `4b465b7`(D5-0 fix)→`f7e0324`/`9ddd296`(D5-0 chore/docs)→`66dbd9a`(D5a)→
`377ccb3`/`6197351`(D5a docs/chore). Related: #38 (= this), #37 (Inc C, the moment-source question),
#242 (the deferred DD curvilinear-diagonal merge), #227 (d≥3 kernels). EXIT after D5b+D6: the
Spatial×Angular tensor-product campaign (the NEXT campaign — extract the curvilinear angular
redistribution into a distinct AngularScheme) + ff-merge `feature/sn-space-angle-tier2`→`main`.
