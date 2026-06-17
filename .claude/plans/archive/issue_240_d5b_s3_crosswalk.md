# D5b-S3 crosswalk — fold Increment C (the spatial-moment iterate φ̂ + scattering Σ_s·φ̂)

Pre-implementation design record (Pattern 7 / L17), from the proactive test-architect
(`.claude/agent-memory/test-architect/d5b_s3_inc_c_moment_iterate_verification.md`) + explorer
passes on `feature/sn-space-angle-tier2` @ `495af60` (S2 done). READ FIRST alongside the durable
plan `.claude/plans/issue_240_phase2_step_d5b_ubld.md`.

## ⭐ THE FP RESOLUTION (load-bearing — do NOT gate on S2-FP==S3-FP)

Inc C is a **PHYSICS-COMPLETION, not an iteration-only change** (test-architect confirmed against the
live code + literature). S2 converges to `(L+C − S_flat)ψ = Q_ext`, `S_flat = Σ_s ⊗ e₀e₀ᵀ` (scatter
ONLY the spatial-average moment; the slope rows get ZERO scattering source) — O(h²) but
**diffusion-limit-INCONSISTENT**. S3 converges to `(L+C − S_full)ψ = Q_ext`, `S_full = Σ_s ⊗ I_{spatial-moment}`
(scatter every spatial moment; slope rows carry `Σ_s·φ̂`) — **diffusion-limit-CONSISTENT**. The
converged FP CHANGES (that is the point — the thick-diffusive tripwire flips xfail→PASS *because* the
FP becomes correct). `S_full ≡ S_flat` at DD/Step n=1 (`Σ_s ⊗ I_1`) — the negative-control bit-identity.
The plan's "FP-invariance vs the flat-source FP" wording is a **Mode-9 mis-application** — REJECTED.
The genuine Mode-9 invariant is: SI-with-lagged-moment-iterate ≡ direct/Krylov solve of the SAME
`(L+C − S_full)` operator (the within-group analog of D5b.4, lifted to the full operator).

## ⭐⭐ RESOLUTION (2026-06-17) — A/B split SUPERSEDED by the unified moment matvec

The A/B split below was drawn on the assumption that "the d=1 matvec just works" (S3-A flips the 1-D
tripwire via the existing d=1 krylov-slab matvec). The S3-A consumer dispatch REFUTED that: the d=1 LD
matvec is **Schur-reduced to a SCALAR** residual, and `A=(L+C)−S` (OperatorSum) subtracts `S·ψ`
ELEMENT-WISE — so a φ̂-carrying `(2,)`-moment `S·ψ` cannot subtract from the scalar `(L+C)·ψ`. The user +
main agent resolved it: **the scalar d=1 matvec was a FLAT-SOURCE artifact (Q̂≡0 ⇒ slope globally
uncoupled ⇒ scalar Krylov unknown), NOT a d=1 Schur degeneration.** A matvec is a forward APPLY (nothing
to eliminate); Inc C makes `Σ_s·φ̂` couple the slope GLOBALLY in every dimension ⇒ the matvec is
moment-valued for ALL d. **DECISION: unify the matvec — `(L+C)·ψ⃗` returns the `spatial_basis_per_axis^d`
moment residual for all d (DD/Step=1=byte-identical; LD-1D=2; LD-2D=4); fold the former S3-B (d≥2 matvec
raise) into ONE all-d step.** Branch-freeness is the acceptance bar (drop the `len(s_axes)>1` moment
conjunct → the scheme trait alone; collapse LD's kernel forks to one d-generic path; retire the scalar
`kernel_rhs` matvec arm; shape-agnostic carriers via `face_moment_tail`; CumprodScan SWEEP untouched,
`scan_xV` gains `s_hat`). Full record: the durable plan's "⭐⭐⭐ THE D5b-S3 ARCHITECTURE DECISION" block.
The "S3-A/S3-B" labels below are now ONE step "S3 unified moment matvec"; GATES 1–5 all in-scope.

## The A/B split (along the sweep vs matvec L14 boundary) — HISTORICAL (see RESOLUTION above)

> ⭐ **UPDATE (2026-06-17): S3-A split again — a typed-field-space prerequisite (S3-A0) emerged + is
> DONE.** The first S3-A dispatch surfaced that carrying φ̂ in the iterate needs the typed-field
> space-shape contract widened first (the `Field.__post_init__` shape gate has no slot for a moment
> axis). That became **S3-A0 (#64) ✅ DONE + committed (`d313d16`/`96dfc96`):** minted
> `SpatialMomentSpace` (user chose the typed-factor design) + the optional `spatial_moments` factory
> param (DEFAULT-OFF, byte-identical) + closed #207 (`find_factor`) + folded in the `Σ_s ⊗ I`
> scattering einsum lift. So the scattering-lift + the typed-space half of "S3-A" below are ALREADY
> DONE; the REMAINING S3-A (#61) is the CONSUMER that SELECTS the now-built SpatialMomentSpace factor.
> API to consume: `.claude/agent-memory/method-implementer/issue_240_d5b_s3_a0_spatial_moment_space_closeout.md`.

- **S3-A (#61) — the forward / SI path (CONSUMER, the remaining work).** SELECT the SpatialMomentSpace
  factor for the φ̂ iterate carrier (both windowed + un-windowed, via the S3-A0 `spatial_moments`
  factory param) + the φ̂ cell-emit accumulation + the source seams. (The `Σ_s ⊗ I` scattering lift is
  DONE in S3-A0.) Flips BOTH thick-diffusive tripwires: the 1-D (#37) via the existing d=1 krylov-slab
  matvec (no raise there — `S` enters the matvec via `S.apply`'s slope lift); the 2-D via SI (the
  sweep path). Gates 1/2/4/5.
- **S3-B (#63) — the d≥2 Krylov matvec.** Close the `_CellResidual.cell:929` raise + widen the apply
  probe + the flat ravel → 2-D Krylov works. Gate 3 (Krylov≡SI on `(L+C − S_full)`).

## The seams (explorer map, file:line @ `495af60` — verify at pickup)

**S2 STALE-PLAN CORRECTIONS:** the plan's `linear_discontinuous.py:430-436 raise` is GONE (S2 landed
the d≥2 dispatch); the S3 raise to close is `_CellResidual.cell` `sweep_graph.py:929-936`. The named
traits `n_spatial_moments`/`is_multi_moment` DO NOT EXIST — only `spatial_basis_per_axis: ClassVar[int]`
(`scheme.py`, DD=1/LD=2); the gate is the inline `len(s_axes) > 1 and spatial_basis_per_axis > 1`
(+ `_n_face_moments` property `loss_representation.py:284` + `face_moment_tail` `_ubld.py`).

1. **The iterate — TWO carriers (the scalar flux is NEVER the between-sweep iterate):**
   - **2-D Cartesian windowed (production):** `_MomentWindowedResolvent` (`solver.py:458`) →
     `InvertibleOperator.solve_moments`; iterate = `TimedFullField(bulk=HarmonicMomentField
     (L+1,2L+1,ng,*spatial))`. The in-sweep accumulator is `moment_buf` (`loss_representation.py:459`).
     φ̂ rides as a trailing `2^d` spatial-moment axis on `moment_buf`/the moment field.
   - **1-D / curvilinear / un-windowed:** iterate = `TimedFullField(bulk=AngularFlux (N,ng,*spatial))`.
     The slab LD cell solve produces a `2^d`-moment `psi_avg` but `_CellSolve.cell:883-884` DROPS to
     slot 0 — so the slope is NOT carried between sweeps today. S3-A must carry it (a spatial-moment
     slot on the un-windowed iterate, gated on `spatial_basis_per_axis > 1`).
   - Construction: `_within_group_si` (`solver.py:672`); windowing `_maybe_window` (`:575`, gated
     `is_cartesian and ndim==2`); cold-start `_windowed_cold_start` (`:607`). SI loop:
     `SourceIteration.solve` (`numerics/iteration.py:448`), `rhs = q_ext + Σ gains.apply(psi)` (`:503`).
2. **`ScatteringOperator.apply`** (`scattering.py:926`, `@singledispatchmethod`): arms for
   `HarmonicMomentField` (`:1093`, windowed) + `AngularFlux` (`:1067`, un-windowed) +
   `TimedFullField` (`:1002`). The `Σ_s ⊗ I` lift = make the per-ordinate einsums
   (`_assemble_per_ordinate_source`, `add_iso_source`, `_aniso_source_from_moment_values`)
   spatial-moment-axis-agnostic: Σ_s carries NO spatial-moment index, so a trailing `2^d` axis rides
   through every einsum as a spectator broadcast. DD n=1 → no trailing axis (`face_moment_tail(1)==()`)
   → byte-identical. ⚠ VERIFY `MaterialXSField.apply_p0_in_scatter` broadcasts over the trailing axis
   (test-architect's flagged assumption).
3. **The source seams (TWO, distinct):**
   - **d≥2 wavefront:** `_ubld_system` (`linear_discontinuous.py:531-536`) lifts a scalar `Q_cells`
     into slot 0 today. S3 threads a genuine `(2^d, ng)` moment source via `_SolveOperands.Q`
     `(N,ng,*spatial)` → trailing `2^d` axis → `Q_cells (N_oct,ng,n_diag,2^d)` (the slice already
     preserves a trailing axis). `_ubld_system`'s `R_source = M·S_moments` einsum already accepts a
     full `2^d` vector — NO `_ubld.py` primitive change. ⚠ `_SolveOperands.Q` is the ONE carrier
     shared by DD+LD — gate the widening on LD-d≥2 (or keep the dual-rank `Q_cells` contract).
   - **d=1 scan (the #37 tripwire path):** `D1ClosedForm.kernel_rhs` (`_ubld.py:324`) HARD-CODES
     Q̂=0; `schur_xV` (`:343-381`) already has an `s_hat` arg but `_kernel_terms`
     (`linear_discontinuous.py:486`) calls `kernel_rhs` FLAT. S3-A threads the scattering-slope into
     the d=1 scan source. (`_ld_source_moments` `:207` is a THIRD legacy ×V per-cell path with a
     `(2,ng)` guard — generalize only if that path is exercised.)
4. **φ̂ accumulation:** `_CellSolve.cell` (`sweep_graph.py:858`), the moment-mode emit `:892-895`
   (`einsum("nlm,ngd,n->lmgd", Y_octant, psi_avg, weights_octant)`). Keep the full `2^d` axis (don't
   drop at `:883-884`) and accumulate `φ̂_ℓ^m,p = einsum("nlm,ngdp,n->lmgdp", ...)`. The spatial-moment
   axis (p) is ORTHOGONAL to the harmonic (ℓ,m) axis. Stay gated for DD (n=1 → no p axis → original).
5. **The d≥2 matvec raise (S3-B):** `_CellResidual.cell` (`sweep_graph.py:929-936`). Entry:
   `StreamingOperator.apply` → `_TwoDimWavefront.loss_action` (`loss_representation.py:966`) →
   `walk_windowed(level_op=_CellResidual(...))`. Closing needs: `_ApplyOperands.probe` carries the
   `2^d` probe; drop the raise + pass full `psi_bar`; `InvertibleOperator` flat ravel
   (`n_unknowns_flat`/`to_flat`) accounts for the widened probe. The kernel `residual_kernel_batch`
   d≥2 arm is ALREADY implemented + verified (D5b.4).

## The gates (test-architect spec — full detail in its memo)

- **GATE 4 (land FIRST) — DD/Step bit-identity (negative control).** `python -O -m pytest
  tests/sn/sweep/core tests/sn/solve -W "error::tests.sn.regression._regression_assert.DriftWarning"`
  stays at the S2 baseline (re-confirm 562/2 skip/4 xfail). The `Σ_s⊗I` lift + φ̂ read are NO-OPs at
  `spatial_basis_per_axis==1`.
- **GATE 1 (PRIMARY) — thick-diffusion limit.** Ref = the continuous diffusion solution (scaled-transport
  ε→0; Adams-2001 / BLA-1992 / LMM-1987 — structurally independent of the LD kernel). Leg 1a: the 1-D
  tripwire `test_mms_ld_slab.py::test_ld_thick_diffusive_limit_xfail` (`:235`, σ_t·h≈10, c=0.99) flips
  xfail→PASS — REMOVE the xfail, rename, **Mode-8-migrate the bare `assert:270`** to
  `np.testing`/`pytest.fail` (inert under `-O`). Leg 1b: a NEW 2-D Cartesian analog (`test_mms_ld_2d.py`)
  via SI, thick box, ref `LD-2D ≈ DD-2D` (diffusion proxy) + continuous band. **2G-het companion
  MANDATORY** (the slope source `Σ_s^T·φ̂` is group-coupled, Mode-6 — a 1G-only tripwire is a degeneracy
  guard failure, Gate 5).
- **GATE 2 — transport MMS still O(h²).** The 1-D LD slab MMS + the 2-D smoke STAY GREEN (their flux
  has spatial slopes → `Σ_s·φ̂` now active; the MMS Q_ext absorbs it). Drift = bug, NOT relax.
- **GATE 3 (S3-B) — the genuine Mode-9.** SI-with-lag ≡ Krylov on the SAME `(L+C−S_full)`. Config:
  ANISOTROPIC vacuum + DIAGONAL `level_symmetric` cubature + NON-ZERO self-scatter c≈0.5-0.8 + 2G-asym
  + NON-SQUARE (non-zero scatter MANDATORY — c=0 makes S_full≡S_flat, the Mode-9 degeneracy). The SUT
  for the closed `:929` raise. Pair with GATE 1 (GATE 3 alone is necessary-not-sufficient — both could
  solve the same wrong operator).
- **GATE 5 — 1G degeneracy guard.** GATES 1+3 carry a 2G-asym het config. Reject a 1G-only scatter gate.

## Claims S3 unlocks vs S4
- **S3 unlocks:** the d≥2 Krylov 2-D LD matvec (S3-B); the thick-diffusion limit (both d); the genuine
  `(2^d, ng)` moment-source carrier.
- **Stays S4:** the `@verifies("ld-cartesian-2d")` flux-SHAPE claim (needs the strengthened Mode-7
  stress ansatz `SN2DCartesianLDStressMMSCase` + the non-vanishing-boundary `BoundaryFlux` moment
  trace). S3's MMS gates stay a convergence smoke + the thick-diffusion VALUE anchor.

## Hazards
- **Mode-9 mis-application** (the headline — do NOT gate FP-invariance vs S2).
- **Mode-8:** the flipped tripwire's bare `assert` must become a function call (fires under `-O`).
- **Mode-6:** the slope source `Σ_s^T·φ̂` is group-coupled → 2G-asym het mandatory.
- **DD fragility:** `_SolveOperands.Q` is the shared carrier — gate its widening on LD-d≥2.

## Citations
Adams-2001 (NSE 137, thick-diffusion-limit verdict), BLA-1992 (JCP 98, 2-D LD asymptotic),
LMM-1987 (JCP 69, the diffusion-limit asymptotic), MRM-2016 (NSE 185, the UBLD system).
