# #240 Phase 2 Step D — home every method into its native base (DiscretizationScheme vs LossRepresentation)

> **Durable in-repo recovery anchor** (project rule: plans live in ORPHEUS/.claude/, not ~/.claude).
> Parent: `.claude/plans/next_principled_polymorphism.md` / `issue_240_phase2_layer_separation.md`.
> Branch `feature/sn-space-angle-tier2`. Approved 2026-06-16.
> **STATUS: D1–D4 + D5-0 DONE + committed. D5 DESIGN PASS DONE** (test-architect spec +
> literature-researcher + cross-domain-attacker — memos under `.claude/agent-memory/`).
> **D5-0 (routing honesty) committed `4b465b7`:** scheme-named ClassVar trait
> `transverse_coupling_is_facewise` (DD `True`, LD default `False`) closes the LIVE
> 2-D-LD→ScanMarch silent-DD misroute (a 2-D LD mesh was silently computing DD); the
> footgun of a presence-only `@runtime_checkable` isinstance is closed with the
> `TestCapabilityTraitsAreGenuineBools` registry-wide bool teeth. **⭐ D5b RESHAPED by the
> literature (Adams 2001 / Maginot-Ragusa-Morel 2016 / Börgers-Larsen-Adams 1992):** the
> Cartesian N-D LD must be the **tensor-product UBLD (`2^d` moments per cell — `{1,x,y,xy}`),
> NOT the simplex `1+d` (`{1,x,y}`)** — the `xy` cross moment is diffusion-limit-load-bearing;
> the simplex object FAILS the thick-diffusion limit on quadrilaterals (would have shipped a
> silent physics bug). The cell/face contract WIDENS (cell unknown `2^d`, each face a
> `2^{d-1}`-moment object) → needs an architecture pass. NEXT = **D5a** (DD scan-march fold),
> then **D5b** (UBLD, needs its own design pass), then **D6**. See §D5 + [[project-issue-158-ld-dag]].
> Commits: A+B `f0d68c3`/`4937c3a` · D1 `8bc1a49` · D2 `784edeb` · D3 `4f04126` · D4 `c40a341`
> (+ chore records each). D4 finding: scheme was ALREADY Σ-stateless → D4 = the
> diffusion-readiness contract gate + the Base interface note (no code change). The deferred
> model-agnostic param rename (cross-domain-attacker) is filed as **#241**
> (`total_xs`→`reaction_xs`, `cell_average_weight`→`face_blend_weight`; `streaming` KEPT).
> All bit-identical/principled; gates green throughout (strict 505/1/4, full 1083→1093).

## Context

Steps A + B are committed. Step B established the loss-rep reads Σ from the operator's C
(single-sourced). Step C (affine-closure → methods) was INVERTED to run after D, then the
inversion DISSOLVED once a literature evaluation (Larsen-Morel 1989, Larsen-Morel-Miller
1987, Bailey-Morel-Chang 2010, Hébert + CFD box/κ-scheme refs) **confirmed the user's
structural claim**:

* **Step ↔ first-order upwind** (donor cell, `ψ̄=ψ_out`, w=1).
* **DD ↔ central — precisely the Keller box scheme** (central advection + trapezoidal
  reaction; box ≡ diamond proven). DD's negative-flux pathology IS the central `Pe>2` wiggle.
* **LD ↔ DG-P1 with upwind flux** (linear/2nd-order upwind; θ=⅓ = P1 Legendre normalization).
* The coefficient triple `(a, 1/denom, w)` is a **generic advection–reaction discretization**
  (swap `|μ|→a`, `Σ_t→r`; formulas unchanged). `w(Σ)` ½→1 is the generic **Péclet/κ-scheme
  blend** (SymPy-verified `w→½` as Σ→0), NOT an SN artifact.
* Genuinely transport-specific = the **ANGULAR axis**: quadrature ordinates + the curvilinear
  angular redistribution (`(1−μ²)/r·∂_μψ`, α-recursion, Morel-Montry τ). Literature splits
  the diffusion-limit analysis into a *spatial* paper (LMM-1987) and an *angular* paper
  (BMC-2010) — space ⊗ angle are factorizable axes.

**Decisive consequence (refutes the original Step C/D framing in
`issue_240_phase2_layer_separation.md`):** the assembly + closure ops are NOT SN-loss-rep
realization that leaked into the scheme — they are a **genuinely generic advection–reaction
*spatial* discretization** that correctly belongs ON the scheme. The assembly does NOT
relocate to the loss-rep; the closure ops stay in the spatial layer → no circular import, no
C-before-D inversion needed.

**User directive (2026-06-16):** task unchanged (home the dangling functions/leaked methods
into their native base) — only the *how*. Each member → **DiscretizationScheme Base**
(generic advection–reaction spatial) or **LossRepresentation Base** (SN angular/sweep); the
cut is the proposal below (up for discussion). The generic advection–reaction form **must be
REALIZED** (Σ-stateless, parameterized by wave-speed + reaction-rate) so SN consumes it
cleanly AND **in preparation for diffusion** — the CONFIRMED next model (standalone diffusion
solver AND consistent-DSA preconditioner, #2). The full **Spatial × Angular tensor product**
is the NEXT campaign, from this clean base.

## The homing (proposed — "what belongs where")

Native-base test: *function of (wave-speed, reaction-rate, source, basis) → generic
advection–reaction → DiscretizationScheme Base; angular/ordinate/sweep-structure →
SN-specific → LossRepresentation Base.*

**→ DiscretizationScheme Base (generic advection–reaction SPATIAL — diffusion-consumable):**
* Reconstruction ops — `cell_average` (`(1−w)·in+w·out`), `source_emission` (`QV·inv/w`),
  NEW `outgoing_face_from_average` (`(ψ̄−(1−w)·in)/w`). Currently DANGLING free functions in
  `spatial/affine_closure.py` → become Base methods/staticmethods.
* Coefficient assembly — `affine_scan_coefficients`, `cell_kernel_batch`,
  `residual_kernel_batch`, the LD Schur (`_schur_terms`/`_kernel_terms`). Already on the
  scheme (`diamond.py`/`linear_discontinuous.py`) — CONFIRMED correctly placed; Σ-stateless (D4).
* Basis — `θ` (LD), the `w`-rule, the moment layout `_ld_source_moments`, traits.

**→ LossRepresentation Base (SN-specific angular + sweep):**
* Sweep/matvec walks (octant/scan/wavefront) — already there.
* The **angular redistribution** (Morel-Montry α-recursion/τ) — FUSED in `diamond.py`'s
  curvilinear `update` + loss-rep `_apply_walk` curvilinear thread. *Cartesian* inlined
  `2ψ̄−in` homes NOW via `outgoing_face_from_average`; the deep **curvilinear angular
  extraction is the NEXT campaign** (the spatial⊗angular fusion already guarded with
  `NotImplementedError`).
* Ordinate/octant structure, per-cell sweep vocab (`CellVisit`/`UpstreamState`/`CellResult`).

## Sequenced sub-steps (each verifiable + committable; route around the 7 reds; NEVER all `tests/sn` — #212)

**D1 — `outgoing_face_from_average` + collapse the inlined closures.** Add the generic
inverse `(ψ̄−(1−w)·in)/w`; route the inlined sites through it: DD `diamond.py:169`/`:346`/
`:384`, LD `linear_discontinuous.py:459`/`:484`, `scan.py:344` + `loss_representation.py:1431`/
`1452`. BIT-IDENTICAL at each scheme's w. Kills the realization-leak duplication. [old-Step-C core]

**D2 — home the reconstruction ops onto the DiscretizationScheme Base.** Move
`cell_average`/`source_emission`/`outgoing_face_from_average` from `affine_closure.py` free
functions → Base methods/staticmethods. Migrate callers (loss-rep scan-solve
`loss_representation.py:2632/2647/2800/2818` + D1 sites). No circular import (spatial layer
owns them; loss-rep already imports the scheme).

**D3 — the rename.** `CellUpdateBase→DiscretizationSchemeBase` FIRST (substring collision),
then `CellUpdate→DiscretizationScheme`, then `cell_update→scheme`. Nexus `rename`, MAIN tree
only (NEVER the worktree). Blast radius: ~29 prod+16 test (`*Base`), ~26 prod+12 test
(`CellUpdate`), ~78 prod+65 test (`cell_update` attr), + docs. Registry `key=` strings
UNCHANGED. Mechanical, bit-identical. (Module-file rename OPTIONAL.)

**D4 — realize the Σ-stateless generic advection–reaction form.** Coefficient methods take
reaction-rate (Σ) as an EXPLICIT param sourced from C (extends Step B); scheme holds NO Σ
state. Interface → `(wave-speed, reaction-rate, source, geometry) → coefficients/cell-update`,
diffusion-ready. `cell_balance.total_xs` already supports this.

**D5 — complete spatial-closure polymorphism in N-D (#239 + N-dim LD).** ⭐ RE-SCOPED by the
user (2026-06-16): *"LD, like any other spatial closure, needs to be able to run on all sweep
strategies."* Complete-polymorphism principle ([[feedback-principled-over-bit-identical]],
"equivalent capabilities across polymorphism"): every spatial closure (DD/Step/LD) must have a
VALID path on every sweep strategy in every dimension; the only permitted exclusions are
STRUCTURAL (declared via `supports()`/`Compatibility(reason)`), never arbitrary code gaps.
**CURRENT GAP: LD is 1-D-ONLY** — it rides `CumprodScan` (1-D scan) but its multi-D kernel
(`_kernel_terms`/`cell_kernel_batch`) RAISES `NotImplementedError` for `len(s_axes) != 1`. DD runs
1-D + N-D (scan + wavefront); LD's N-D is the hole. D5 is now a CAMPAIGN (own design pass + a
proactive test-architect FIRST), two threads:

* **D5a — #239: the N-D DD scan-march onto the coefficient model.** The 2-D `ScanMarch`
  (`loss_representation.py` `_sweep_interior`/`_loss_action_interior`) inlines DD `2g`/`alpha`/`beta`;
  fold onto the coefficient model via `outgoing_face_from_average` (D1) + an **N-D coefficient cache**
  (`CollisionCache`/`from_geometry` is 1-D-only) + the transverse `s_y·ψ_y` folded into
  `source_emission`. Principled ~1-ULP re-baseline (#158-B1; two-paths oracle + cart2d apply snapshot,
  Mode-9 het/non-square/diagonal cubature).
* **D5b — N-dim LD on the DAG wavefront — SUBSUMES #158 Increment D (#38).** ⭐ **RESHAPED by the
  D5 literature pass** (`.claude/agent-memory/literature-researcher/multi_d_ld_closure.md`): the
  "BILINEAR = one slope per axis (`1+d` moments)" intuition is **WRONG for Cartesian cells**. Two
  distinct objects: the **simplex-P1 LD** (`{1,x,y}`, `1+d` moments) **FAILS the thick-diffusion
  limit on quadrilaterals** (Adams 2001, NSE 137); the **tensor-product UBLD** (`{1,x,y,xy}`, **`2^d`
  moments** — 4 in 2-D, 8 in 3-D) **preserves it** (the `xy` cross moment is the load-bearing term;
  Maginot-Ragusa-Morel 2016, NSE 185, UBLD Eqs. 1-12; Börgers-Larsen-Adams 1992). **So D5b builds the
  UBLD, not the simplex** — building the naive `1+d` object would ship a silent diffusion-limit
  physics bug. The per-cell system is a **`2^d×2^d` dense Galerkin solve** (`A = G + F_out + σ_t·M`);
  the elegant build = **assemble M/G/F as Kronecker products of the verified 1-D LD operators**
  (d=1 reduces to identity, single-sources the 1-D math, oracle = "exact on a bilinear flux"). LD
  still rides the **DAG wavefront** (`FullFieldWavefront`/`MovingFrontierWindow`, pure scheme-delegators
  via `sweep_graph.py:849/890`); close the `linear_discontinuous.py:430` `len(s_axes)!=1` raise with
  the `2^d` kernel. ⚠ **Cell/face contract WIDENS**: the cell unknown is `2^d`-long and each outgoing
  face is a `2^{d-1}`-moment object (avg + transverse slope), NOT scalar — exceeds the current
  `CellResult`/`WavefrontFlux` scalar slots → **architecture pass before coding** (the 1-D Schur-scalar
  collapse likely does NOT survive to multi-D). Lumping (FLBLD/LLD) deferred — `is_positivity_preserving`
  stays `False` (no bilinear DFE is positivity-preserving). The scan-march stays a DD/Step (slopeless)
  optimization whose `supports()` already (D5-0) excludes LD via `transverse_coupling_is_facewise=False`
  — the STRUCTURAL exclusion. Verify (test-architect spec §2-D5b + cross-domain-attacker MMS frame):
  each scheme ≥1 valid strategy per dim; routing-flip DD≠LD; multi-D UBLD MMS O(h²) with a `Q̂≠0`
  slope source AND **x↔y-broken cross-harmonics in the slope drivers B,C** (the cross-domain-attacker
  strengthening — a reflection-symmetric same-sign slope-row bug would else cancel in the flux); the
  `ld-thick-diffusive` tripwire flips xfail→PASS on UBLD (would stay FAIL on the simplex — correct
  physics, document the basis-dependence). INVERT `test_cell_kernel_batch_rejects_multi_d`. D5b needs
  its OWN design pass (the contract widening + the Kronecker assembly).

**Connections:** the `SweepStrategy` first-class abstraction (`.claude/plans/sn_sweep_strategy.md` —
protocol `sweep`+`residual`+`supports`; `Compatibility(ok,reason)` selection on `is_cartesian`+`ndim`)
is the natural home for the `supports()` exclusions; #158 Increment D (#38) is SUBSUMED. `is_affine_scannable`
likely needs to be dimension-aware (LD: True 1-D / False multi-D) OR `supports()` carries it.

**D6 — docs + next-campaign issue.** Archivist: add Step/DD/LD ↔ advection-scheme
correspondence to `docs/theory/discrete_ordinates.rst` + scheme docstrings (upwind/box-central/
DG-P1-upwind, box≡diamond nuance, `w(Σ)` Péclet coupling, spatial⊗angular factorization,
lumped-step-vs-step-characteristic). File a NEW issue for the full **Spatial × Angular tensor
product** (extract curvilinear angular redistribution → distinct AngularScheme).

## Verification (test-architect spec)

* **D1/D2/D3: BIT-IDENTICAL** — strict DriftWarning gate (CORRECT path
  `tests.sn.regression._regression_assert.DriftWarning`):
  `python -O -m pytest tests/sn/sweep/core tests/sn/solve -W "error::tests.sn.regression._regression_assert.DriftWarning"`
  (baseline 505/1/4). Rename: migrate name/attr-pinned tests; oracle/solver tests survive.
* **D4: bit-identical** prod + a NEW capability gate (explicit reaction-rate ≠ mesh Σ → matches
  direct build), mirroring Step B's removal-form teeth.
* **D5a (#239 DD scan-march): ~1-ULP** — two-paths oracle `test_scanmarch_sweep_equals_oracle`/
  `_residual_equals_oracle` + a cart2d apply arm on `TestT4bPreT4RegressionSnapshot` (Mode-9:
  het/2G-asym/non-square, diagonal cubature).
* **D5b (N-dim LD): NEW capability** (not bit-id) — multi-D LD MMS O(h²) on the wavefront (the
  LM-1989 slope-sign trap; MUST include a non-vanishing slope source `Q̂≠0`, the recurring coverage
  gap); a per-(scheme×strategy×dim) coverage matrix (every scheme has ≥1 valid strategy per dim;
  `supports()` rejects LD-on-scan-march-in-N-D with the structural reason, asserted as a NEGATIVE
  gate); routing-flip DD≠LD on the same N-D config. ⚠ proactive **test-architect FIRST** — this is a
  cross-strategy/cross-dim carve (the L17 trigger) AND a new L0/L1 reference (multi-D LD MMS).
* **Basis/contract gate (vv L11):** thin scheme exposes basis/traits; recon ops reachable as
  Base methods; Mode-8 negative arms use `pytest.fail`.
* **Routing-flip (#158-B3 analogue):** DD vs LD same 1-D het 2G → DIFFERENT + each scheme's
  absolute anchor (LD `TestLDLinearExactness`/`test_mms_ld_slab` O(h²); DD k∞).
* **Route-arounds:** `-k "not (vacuum_bulk_bit_identical_1d and SPH) and not (sphere_1g_apply_bit_identical or sphere_2g_apply_bit_identical) and not test_2d_mesh_resolution and not two_d_cartesian_loss_action"`.
  Dirs: `tests/sn/operators tests/sn/spatial tests/sn/sweep/core tests/sn/sweep/cartesian_2d tests/sn/solve`.
* Per sub-step: elegance-enforcer + qa; Sphinx clean; `tests._harness.audit` exit 0.

## Execution discipline
method-implementer for D1/D2/D4/D5 against the spec; elegance + qa before each commit;
archivist for D6. One commit per sub-step + a `chore(claude)` records commit; trailer
`Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`; explicit paths only
— NEVER the 3 forbidden untracked / parallel-session set / `docs/_build/`. `.claude/agent-memory/*`
→ chore commit. Branch NOT pushed/merged. **Step C (#48) is SUBSUMED by D1+D2.**

## EXIT
The full Spatial × Angular tensor-product campaign (the D6 issue). Then resume
`issue_158_spatial_cellupdate_carve.md` §"Increment C" (#37) → Step-4 (#36) + Inc D (#38);
ff-merge → `main` + apply the Phase-0 hand-off.
