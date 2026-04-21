# Next-session plan — rank-N hollow-sphere closure: Hébert extraction and beyond

**Branch entering this work**: `investigate/peierls-solver-bugs` (at commit `4a169ea`).
**Author of this plan**: Claude Opus 4.7, 2026-04-21.
**Audience**: fresh Claude Code session with zero context. Read §§0-4 in full before touching any code.

---

## 0. One-paragraph summary (read first)

After ~150k tokens of investigation across multiple agent sessions, we've established that ORPHEUS's Phase F.4 scalar rank-2 white-BC closure on hollow sphere gives **0.077 % k_eff residual at σ_t·R = 5, r_0/R = 0.3** and IS literally the textbook closure (Stamm'ler 1983 Ch. IV, Ligou 1982 Ch. 8, Stacey 2007 Ch. 9 — three independent derivations). The "canonical" rank-N Legendre-moment extension of Sanchez-McCormick 1982 §III.F.1 that we spent most of the investigation trying to implement **plateaus at 1.42 % and is uncross-validated by any successor textbook** — it's a theoretical construction that didn't cross over into the canonical literature. The structural reason is proven: at σ_t = 0, `W_oo[0,0] + W_io[0,0] = 1` exactly (F.4 mode-0 conservation), but at modes n ≥ 1 the naive conservation `W_oo[n,n] + W_io[n,n]` falls to 0.28, 0.13, 0.09 because the inner-surface angle remapping `c_in = √(1 − (R/r_0)²(1 − µ²_emit))` makes `P̃_n(c_in) ≠ P̃_n(µ_emit)` — which the Legendre-moment basis cannot represent. The only remaining authoritative reference not yet checked is **Hébert 2009 *Applied Reactor Physics* Chapter 3**, which the user has just dropped at `/workspaces/ORPHEUS/Hebert(2009)Chapter3.pdf` (122 pages, Adobe Paper Capture scanned OCR). **Hébert is the last reference that might contain a correctly-reducing Legendre rank-N ladder with c_in-aware machinery — check it first, decide from there.**

---

## 1. Context state — everything that has landed

### 1.1 Commits on `investigate/peierls-solver-bugs` (most recent first)

- `4a169ea` — docs: F.4 quadrature-floor data added to close-out.
- `a640a83` — docs: **four-reference synthesis** (Ligou, Sanchez 2002, Stamm'ler IV, Stacey 9). **THIS is the key commit — the headline finding.**
- `a2e2205` — diag: closure characterization + cross-σ_t·R parameter scan.
- `0b0533b` — diag: Sanchez-McCormick §III.F.1 recipe investigation (60+ variants, plateaus at 1.43%).
- `53fae60` — fix: `solve_peierls_1g` forwards `inner_radius` to `composite_gl_r` (one-line pre-existing bug fix) + 33 MC cross-check tests of the W matrix (all pass — W is CORRECT).
- `ca9d68f` — feat: Marshak partial-current per-face primitives (infrastructure only, dead code behind guard).
- `cf6ab48` — docs: earlier Marshak rank-N plan (superseded by four-reference synthesis).
- `b9bc3df` — diag: measure-mismatch diagnostics (earlier hypothesis, now partially refined).
- `d890a1e` — feat: rank-N per-face infrastructure (Lambert primitives + (2N×2N) W).
- Phase F.4 and earlier: slab rank-2 + hollow cyl rank-2 + hollow sph rank-2 scalar closure.

### 1.2 What's VERIFIED (tested and green)

- Phase F.4 slab rank-2 scalar: machine-precision bit-exact vs legacy `E_2`/`E_3` bilinear.
- Phase F.4 hollow cyl/sph rank-2 scalar: 1 %-13 % residual across r_0/R (per-face closure works).
- `compute_hollow_sph_transmission_rank_n` at N=1 reduces to scalar W bit-exact.
- Sanchez-McCormick reciprocity `A_k · W_{jk}^{mn} = A_j · W_{kj}^{nm}` at 1e-13.
- 4M-sample MC cross-check of W at (m,n) up to (2,2) — analytical formula is faithful integral.
- σ_t = 0 closed-form values for W_oo, W_io match to 1e-10+.
- F.4 mode-0 conservation `W_oo[0,0] + W_io[0,0] = 1` holds exactly at σ_t=0.
- Inner_radius bug in `solve_peierls_1g` fixed + regression test pins it.

### 1.3 What's CHARACTERIZED (numerically documented, not bit-exact)

- F.4 N=1 residual across σ_t·R ∈ {1, 2.5, 5, 10, 20}: 3.27%, 0.55%, **0.077%**, 0.26%, 0.45%.
- F.4 quadrature floor at σ_t·R = 5: ~0.04-0.1 % (Mark DP_0 truncation, NOT pure quadrature).
- Sanchez µ-ortho recipe N=2 through N=4: all 1.42 % at σ_t·R = 5 (literal plateau).
- Sanchez conservation `W_oo[n,n] + W_io[n,n]` at σ_t=0: n=0 gives 1.000 (F.4 identity), n≥1 gives 0.28, 0.13, 0.09 (structural failure).
- 60+ recipe variants tested across Lambert/Marshak P/G, with/without (ρ/R)² Jacobian, (2n+1) Gelbard factor placements, B^µ Gram conversions: **best is 1.42 %**, none below 1 %.

### 1.4 What's IMPLEMENTED but DEAD CODE behind the guard

- `compute_P_esc_outer_mode_marshak` (lines 2195-2290 of `peierls_geometry.py`).
- `compute_P_esc_inner_mode_marshak`, `compute_G_bc_outer_mode_marshak`, `compute_G_bc_inner_mode_marshak`.
- `_build_closure_operator_rank_n_white` — rank-N assembly helper using Marshak primitives.
- `compute_hollow_sph_transmission_rank_n` — (2N × 2N) transmission matrix.

All reachable ONLY if the `NotImplementedError` guard in `build_closure_operator` is lifted (which we WILL NOT do until a correctly-reducing recipe is found).

### 1.5 Literature state

Extracted canonical formulations (in `.claude/agent-memory/literature-researcher/`):

- `sanchez_mccormick_rank_n_per_face.md` — the 1982 §III.F.1 recipe (uncross-validated).
- `rank_n_closure_four_references_synthesis.md` — **the critical 4-way comparison**. Ligou, Sanchez 2002, Stamm'ler Ch. IV, Stacey Ch. 9 all scalar/DP-0; Sanchez 2002 is piecewise-constant sectors, not Legendre.
- `stammler_1983_ch6_interface_currents.md` — Ch. 6 negative finding (Ch. 6 is SN, Ch. IV is CP).
- `rank_n_interface_current_canonical.md` — tagged Sanchez-McCormick 1982 §III.F.1 as not cross-validated.
- `cp_moment_integrals.md` — CP integration formulas (Ch. IV specifically).

**NOT YET EXTRACTED**: Hébert 2009 *Applied Reactor Physics* Ch. 3.

---

## 2. The structural diagnosis (proven, no need to re-test)

At σ_t = 0 on hollow sphere, the SHIPPED `compute_hollow_sph_transmission_rank_n` (MC-verified) gives:

```
Mode 0:   W_oo[0,0] + W_io[0,0] = 1.000     ✓ F.4 identity (conservation)
Mode 1:   W_oo[1,1] + W_io[1,1] = 0.281     ✗ naive conservation fails
Mode 2:   W_oo[2,2] + W_io[2,2] = 0.134     ✗
Mode 3:   W_oo[3,3] + W_io[3,3] = 0.092     ✗
```

**Why mode 0 works**: P̃_0 = 1 is angle-independent. The inner-surface cosine remapping `c_in(µ_emit) = √(1 − (R/r_0)²(1 − µ²_emit))` is INVISIBLE at mode 0 because P̃_0(c_in) = P̃_0(µ_emit) = 1 trivially.

**Why modes n ≥ 1 fail**: `P̃_n(c_in) ≠ P̃_n(µ_emit)`. The W_io integrand couples the emission basis P̃_n(µ_emit) with the arrival basis P̃_m(c_in), but the two µ variables are related by the non-trivial `c_in` mapping. The Legendre-moment basis doesn't diagonalize this mapping.

**For a true rank-N closure to work**, either:
1. The basis must diagonalize the `c_in → µ_emit` map (e.g., a geometry-adapted basis), OR
2. A completely different angular representation must be used (e.g., Sanchez 2002 piecewise-constant sectors, which handle c_in as a direct ray-by-ray reprojection within the sector framework).

---

## 3. The plan for the next session

### 3.1 PRIORITY 1 — Dispatch literature-researcher on Hébert 2009 Ch. 3

**Why Hébert is the last-best-hope reference**:
- Hébert is a MODERN graduate textbook (2009, now Applied Reactor Physics 3rd ed 2020).
- Ch. 3 covers "collision probability and interface-current methods" based on our prior bibliographic research.
- It's the only reference in the user's priority list (per `AGENT.md`) that might rigorously derive a Legendre rank-N IC method for curvilinear cells.
- If anyone has documented the c_in-aware rank-N closure, it's Hébert.

**Dispatch prompt (pre-staged, use directly)**:

```
Task: Extract rank-N per-face interface-current closure from Hébert 2009
*Applied Reactor Physics* (2nd ed) Chapter 3. PDF at
/workspaces/ORPHEUS/Hebert(2009)Chapter3.pdf (122 pages, Adobe Paper
Capture scanned — use Read with pages=X-Y to OCR via vision).

CRITICAL context to read first:
1. /workspaces/ORPHEUS/.claude/agent-memory/literature-researcher/rank_n_closure_four_references_synthesis.md
   (the 4-way comparison showing all known references use scalar/DP-0, NOT Legendre rank-N).
2. /workspaces/ORPHEUS/.claude/plans/next-session-rank-n-hebert-and-beyond.md
   (this plan, §2 structural diagnosis).

Focus: identify whether Hébert Ch. 3 contains a Legendre-moment rank-N
IC closure for CURVILINEAR (hollow cylindrical or spherical) cells
that CORRECTLY handles the inner-surface angle remapping
c_in = √(1 − (R/r_0)²(1 − µ²_emit)) ≠ µ_emit.

Specific questions:
1. Does Hébert derive an IC method with per-face N-mode basis?
2. What basis is used (P̃_n? µ-weighted orthonormal? geometry-adapted?)?
3. How does the transmission matrix W handle outer→inner ray transitions
   where the arrival cosine differs from emission cosine?
4. Is there a worked hollow-cylinder or hollow-sphere example demonstrating
   k_eff → k_inf as N → ∞?
5. Does the method reduce to F.4's scalar closure (Stamm'ler Eq. 34) at N=1?
6. Is there a discussion of the c_in remapping or any related angular
   mapping? Any proofs of convergence for curvilinear geometry?

Binary-search page probe: 122 pages is a large chapter. Start with
TOC (pages 1-3), then jump to where §3.N treats IC/CP per-face modes
(typical section names: "interface current method", "multi-mode
collision probability", "anisotropic boundary conditions").

Output: /workspaces/ORPHEUS/.claude/agent-memory/literature-researcher/hebert_2009_ch3_interface_currents.md
(full extraction with equation numbers).

Reply ≤ 500 words focused on:
- Does Hébert solve the c_in remapping problem?
- If YES: what's the recipe? Does it reduce to F.4 at N=1?
- If NO: what does Hébert use instead (scalar? DP-0? PCA sectors?)?

No need to re-extract the already-covered references. Focus entirely
on Hébert 2009 Ch. 3.
```

### 3.2 PRIORITY 2 — Decision matrix based on Hébert finding

Three branches:

**Branch A — Hébert has a correctly-reducing Legendre rank-N ladder**:
1. Transcribe the key equations carefully (equation numbers + page references).
2. Dispatch numerics-investigator to implement the Hébert recipe.
3. Verify: at N=1 it matches F.4 bit-exactly; at N=2+ the Sanchez 1.42 % plateau breaks.
4. If it closes below 0.1 %, land it:
   - Rewrite `_build_closure_operator_rank_n_white` with the new primitives.
   - Lift the `NotImplementedError` guard in `build_closure_operator` (line 3481 area).
   - Lift the sibling guard in `solve_peierls_1g(boundary="white_rank2")` (line ~4004).
   - Add acceptance tests in `tests/derivations/test_peierls_rank2_bc.py`.
   - Update Sphinx at `docs/theory/peierls_unified.rst` (the Phase F.5 subsection at line 2989 already has a "Phase F.5 open" section — replace with a "Phase F.5 landed" description of the correct closure).
   - Retire the `NotImplementedError` regression test in `test_peierls_rank2_bc.py`.

**Branch B — Hébert does NOT solve the c_in problem but describes Sanchez 2002-style piecewise-constant angular sectors**:
1. Decision: implement the PCA sector closure (significant infrastructure lift).
2. Prep work:
   - New data structure for angular sector partitioning (hemisphere split into Nθ × Nφ cones).
   - New P, G, W primitives computed with sector-averaged integrands.
   - Sanchez 2002 Eq. 37: basis functions = characteristic functions on sectors.
   - Orthogonality: `∫ f^ρ · f^ν · (Ω·n) dΩ dA = c_{a,r} δ` with `c_{a,r} = ∫_a dS ∫_r |Ω·n| dΩ`.
   - At N² = 1 (single sector) this reduces to F.4 scalar — good reduction target.
3. Land if it closes below 0.1 %.

**Branch C — Hébert also just does scalar/DP-0 (most likely outcome)**:
1. **Close Issue #119 formally** with F.4 as the production closure.
2. Deprecate the Marshak primitive infrastructure OR keep it in case someone picks up Sanchez 2002 PCA later. Recommend KEEP as it's well-tested and didn't introduce regressions.
3. Update documentation:
   - `docs/theory/peierls_unified.rst`: rewrite the Phase F.5 subsection to explain the scalar closure's textbook status (cite all 4 references), present the cross-σ_t·R residual data, and describe the angle-remapping structural diagnosis.
   - Update Sphinx Key Facts to reflect "F.4 rank-2 scalar IS the method" rather than "F.4 is a temporary rank-0 limit awaiting rank-N".
4. Close Issue #119 with a referenced summary comment: link commits, the four-reference synthesis memo, and the conservation-probe diagnostic.
5. File a research-tag follow-up issue for Sanchez 2002 PCA if future work desires.

### 3.3 PRIORITY 3 — If the session has bandwidth after Branch closure

A promising novel research direction NOT required for closure:

**Derive a c_in-aware Legendre rank-N closure from first principles.**

The idea: parameterize the angular distribution at the INNER surface in terms of `c_in` directly, not `µ_emit`. Then P̃_n(c_in) acts as the natural basis at inner, and the outer→inner transmission involves a Jacobian `dµ_emit / dc_in` that can be computed analytically.

Sympy exploration:
- Define the outer→inner map `c_in(µ_emit) = √(1 − (R/r_0)²(1 − µ²))` symbolically.
- Compute the Jacobian `dc_in / dµ_emit = (R/r_0)² · µ_emit / c_in`.
- Build a TRANSFORMED basis pair: {P̃_n(µ_emit) at outer, P̃_m(c_in) at inner}.
- Work out the rank-N transmission matrix in this basis.
- Check if the closure `(I - W)⁻¹ P` now reduces to F.4 at N=1 AND converges to k_inf at N→∞.

This is genuinely novel work (not in the four references) and is where the user's "exercise math/physics" spirit applies. If it works, it's publishable. If it doesn't, it's still great diagnostic evidence for Branch C.

---

## 4. Tests and diagnostics to run after any implementation

### 4.1 Regression pins (MUST stay green after any changes)

```bash
python -m pytest tests/derivations/test_peierls_rank2_bc.py -v
```

All 30+ Phase F tests. Specifically:
- `test_solve_peierls_1g_hollow_sph_white_rank2_inner_radius_plumbing` — inner_radius fix.
- `test_sphere_marshak_G_equals_4x_P_sigt_zero` — Marshak primitive identity.
- `test_rank_n_white_closure_raises_pending_normalisation` — STAYS until recipe works.
- `test_hollow_sph_rank2_beats_rank1_mark` — F.4 residual at r_0/R ∈ {0.1, 0.2, 0.3}.

### 4.2 Diagnostic probes

All in `derivations/diagnostics/`:

- `diag_rank_n_W_mc_crosscheck.py` — W is correct (33 tests).
- `diag_rank_n_sph_marshak_primitives_sigt_zero.py` — Marshak primitive σ_t=0 verification.
- `diag_rank_n_15_N1_reduction_model_split.py` — F.4 Model A/B split.
- `diag_rank_n_16_mode1_scale_sensitivity.py` — optimal-c scaling scan.
- `diag_rank_n_sph_keff_probe.py` — k_eff probe across 7 conventions.
- `diag_sanchez_fractional_scan.py` — 16 α-weight combinations.
- `diag_sanchez_N_convergence.py` — Sanchez N=1..4 plateau proof.
- `diag_hybrid_f4_plus_sanchez.py` — F.4 mode-0 + Sanchez n≥1 hybrid.
- `diag_rank_n_sanchez_conservation_probe.py` — conservation structural diagnosis (new this session).
- `diag_rank_n_closure_characterization.py` — F.4 vs Sanchez σ_t·R scan.

### 4.3 Target acceptance at hollow sphere, σ_t·R = 5, r_0/R = 0.3

| Closure | N=1 err | N=2 err | Target |
|---------|---------|---------|--------|
| F.4 (scalar) | 0.077 % | N/A | baseline |
| Any rank-N that closes | 0.077 % (reduce-target) | ≤ 0.1 % (gate) | ≤ 0.01 % at N ≥ 4 |

Regress if ANY of:
- N=1 bit-exactness with F.4 breaks (>1e-10 relative deviation).
- N=2 residual >0.1 % (the gate).
- Non-monotone convergence under N refinement (N=3 worse than N=2 by >10 %).

---

## 5. Context artefacts (reading order for the next session)

### Essential (read first)

1. **This plan** (`.claude/plans/next-session-rank-n-hebert-and-beyond.md`) in full.
2. **Four-reference synthesis memo** (`.claude/agent-memory/literature-researcher/rank_n_closure_four_references_synthesis.md`) — the critical cross-reference.
3. **Previous close-out** (`.claude/plans/post-four-reference-synthesis-close-out.md`) — cross-parameter data and structural diagnosis.

### Code

4. `orpheus/derivations/peierls_geometry.py` — full 4000-line module. Key sections:
   - Lines ~1500-1900: scalar P_esc / G_bc primitives (Lambert, Model A/B).
   - Lines ~2200-2500: Marshak per-face primitives (dead code behind guard).
   - Lines ~2920-3050: `compute_hollow_sph_transmission_rank_n` (W matrix, MC-verified).
   - Lines ~3400-3500: scalar `compute_hollow_sph_transmission` (F.4 N=1 W).
   - Lines ~3480-3530: `build_closure_operator` with guard.
   - Lines ~3650-3700: `_build_closure_operator_rank2_white` (F.4 scalar assembly).
   - Lines ~3800-3900: `_build_closure_operator_rank_n_white` (Marshak rank-N dead code).

### Diagnostic / investigation

5. `derivations/diagnostics/diag_rank_n_sanchez_conservation_probe.py` — the structural diagnosis (mode-0 conservation holds, n≥1 fails).
6. `derivations/diagnostics/diag_rank_n_closure_characterization.py` — cross-σ_t·R scan + quadrature refinement.
7. `derivations/diagnostics/derive_mu_weighted_basis.py` — sympy µ-ortho basis derivation (Jacobi P^{(0,1)}_n).

### Agent memory

8. `.claude/agent-memory/literature-researcher/sanchez_mccormick_rank_n_per_face.md` — Sanchez 1982 §III.F.1 full extraction.
9. `.claude/agent-memory/numerics-investigator/peierls_rank_n_sanchez_closure_failed.md` — 60+ recipe empirical scan.
10. `.claude/agent-memory/numerics-investigator/peierls_rank_n_W_mixed_basis.md` — earlier MC cross-check findings.

---

## 6. If the Hébert extraction surprises us (expect Branch C but be open)

Flags that would pivot to Branch A (Hébert has the closure):
- Hébert §3.N presents the IC method with a per-face mode expansion on curvilinear cells.
- Worked examples show k_eff converges to k_inf under N refinement.
- Explicit treatment of the `c_in` mapping with an "angular projection operator" or similar.
- Any discussion of why Sanchez-McCormick's §III.F.1 is incomplete for curvilinear.

Flags that would confirm Branch C (Hébert also scalar/DP-0):
- Hébert's IC section equals Stamm'ler Eq. 34 with modern notation.
- No Legendre/Marshak expansion beyond DP-0.
- Comment like "for homogeneous cells the scalar closure is exact up to quadrature."
- Any explicit statement that rank-N Legendre IC is theoretical / not used in practice.

---

## 7. Budget estimate

- **Hébert extraction (literature-researcher dispatch)**: ~30-60 min. The PDF is 122 pages scanned — OCR time + analysis.
- **Branch A implementation (if applicable)**: ~3-5 hours. Mode primitives + closure assembly + tests + Sphinx.
- **Branch B implementation (Sanchez 2002 PCA)**: ~1-2 days. Major new infrastructure.
- **Branch C close-out**: ~30-60 min. Update docs, close issue, file follow-up.
- **Direction C (novel derivation)**: open-ended. 4+ hours of careful sympy + analysis.

Realistic 1-session target: complete Hébert extraction + either implement Branch A/B or execute Branch C close-out.

---

## 8. If you're fresh Claude and unsure where to start

1. Read this plan top to bottom.
2. Read `rank_n_closure_four_references_synthesis.md` — you need the 4-way comparison context.
3. Read `post-four-reference-synthesis-close-out.md` — you need the quantitative data.
4. Dispatch the Hébert extraction with the pre-staged prompt in §3.1.
5. Based on the result, branch to §3.2 A, B, or C.

---

## End of plan

This is the 8th commit / 3rd session / ~150k tokens of investigation on Issue #119. The structural picture is clear. The remaining unknown is whether Hébert 2009 Ch. 3 presents a correctly-reducing Legendre rank-N ladder. If it does, great — implement it. If it doesn't, close Issue #119 with F.4 as production, and the rank-N problem becomes open research (either PCA sectors or the c_in-aware Legendre derivation).
