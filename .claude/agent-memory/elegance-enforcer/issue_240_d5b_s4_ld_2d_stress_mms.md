---
name: issue-240-d5b-s4-ld-2d-stress-mms
description: PASS-WITH-NITS verdict on #240 D5b-S4 — the 2-D Cartesian LD bilinear-UBLD stress MMS (new L1 ref + 8 foundation gates + theory stub); Branch-1/Branch-2 split, prescribed_inflow reuse, cross-check coverage gap.
metadata:
  type: project
---

# #240 D5b-S4 — 2-D Cartesian LD STRESS MMS (elegance review)

Verdict **PASS-WITH-NITS**. Branch `feature/sn-space-angle-tier2`, NOT committed.
A NEW L1 MMS verifying the S3-landed multi-D bilinear UBLD LD closure (μ-bilinear
ansatz `ψ=(A+μ_x B+μ_y C)/W`, φ=A, non-vanishing trace a0>0, non-square het-2G).
8 foundation gates pass (ran -O, 8/8). Files: `mms/sn.py` +592, `test_mms_ld_2d.py`
+177, new `test_sn_mms_ld_2d_stress_symbolic.py`, `discrete_ordinates.rst` +73 stub.

## What is EXEMPLARY (reinforce)
- **`_LD2D_STRESS_COEFFS` as `(num,den)` pairs single-sourcing the amplitudes**
  across Branch-1 (`sp.Rational(n,d)`) and Branch-2 (`n/d` float). This is the
  RIGHT pattern: the `Branch2==Branch1` cross-check then pins the two EVALUATORS
  agree, not just the symbolic identity. The amplitudes cannot drift. Textbook
  Pattern-2 (single source of truth) + Pattern-7 (normalise at definition).
- **`prescribed_inflow` reuses the Phase-4.6 `BoundarySourceSink.prescribed_inflow`
  primitive CLEANLY** + `build_nonvacuum_fixed_source` (generic over the
  `(external_source(mesh), prescribed_inflow(sn_mesh))` protocol — ONE def shared
  by every non-vacuum case). No boundary-trace re-implementation. The 2-D face slot
  IS `(N,ng,ny)` (face_layout.py:84), so the case's `(N,ng,n_t)` matches — the
  primitive's docstring saying `(N,ng)` is a 1-D-centric DOC in the PRIMITIVE (a
  pre-existing doc-staleness there, NOT this carve's bug).
- **Named A/B/C drivers + their derivatives** (`_drivers` returns the 9-tuple),
  `phi_exact = A`, named `streaming`/`removal`/`in_scatter` intermediates in
  `external_source`. Reads like the PDE residual. Pattern-3 honoured.
- **Honest-scope discipline is best-in-class.** The slope-SOURCE half of the
  LM-1989 trap + transverse face-slope inflow moment are DEFERRED and documented
  in THREE places (case docstring, theory `.. note::`, gate docstrings) with the
  exact production reason (`_lift_external_source_to_moments` zeros slope rows;
  `mesh.trace` carries no `2^{d-1}` transverse-moment axis). Frame-2 §232 honoured.
  This is the anti-pattern-#11 EXCEPTION done right (deferred scope tracked + named
  as a candidate D5b-S5).
- Quadrature-exactness claims VERIFIED (`level_symmetric(4)`: N=24, `<μ_x²>=W/3`
  exact, `<μ_x μ_y>=0`, no pure-z ordinate). Docstring honest.

## NITS (none blocking; main agent resolves)
- **NIT-1 (CONCERN — cross-check coverage gap, the dispatch's Q1).** The Branch-2
  `external_source` is a HAND-TRANSCRIBED twin of Branch-1 `Q_closed` (it does NOT
  lambdify `Q_closed`; it re-types `streaming+removal−in_scatter` term-by-term in
  numpy). The `_LD2D_STRESS_COEFFS` single-sources the AMPLITUDES but NOT the source
  ASSEMBLY FORMULA. That hand-transcription is JUSTIFIED algebra-of-record
  independence (Branch-2 must not import the symbolic path — L11) — NOT a smell to
  collapse. BUT the cross-check `test_ld2d_stress_branch2_matches_branch1_source`
  pins it only at **group 0** (`c_0=1`, scaling=identity), **one cell** `(3,2)`,
  with 2G in-scatter folded into an effective `sig_s_eff`. So the Branch-2 path for
  **group 1** (the `c_g=0.4` per-group scaling + the genuine 2G downscatter coupling
  `Σ_s[g',g]·A_{g'}`) is pinned ONLY by the end-to-end O(h²) headline gate, NOT by a
  foundation-tier cross-check. The FD residual test operates on Branch-1 SYMBOLIC
  only — it never touches Branch-2 numpy. Bug habitat: a `c_g` mis-scale or a
  downscatter transpose (`Σ_s[g,g']` vs `Σ_s[g',g]` — ERR-002/ERR-009 class) in the
  group-1 assembly would NOT fire the foundation gate; only the slow L1 headline.
  Elegant fix: parametrize the cross-check over `g∈{0,1}` and ≥2 cells (the
  symbolic layer scales by `c_g` mechanically — multiply `Q_closed` by `c_g` and
  set the per-group `sig_t`/`sig_s_eff`). Cheap, closes the group-axis blind spot.
- **NIT-2 (the dispatch's Q5 — weak structural-independence guard).** ACCEPTABLE
  as-is. `test_v_ld2d_stress_is_structurally_independent_of_ld_kernel` checks only
  finiteness + case type, NOT the import graph. It does NOT assert the source-build
  imports no `_LDCellTerms`/`_ubld`. Verdict: leave it (an AST/import-graph assert
  is brittle + the real independence is structural — Branch-1 builds from the
  continuous PDE, Branch-2 from `_drivers`, neither imports the kernel; that's
  visible by inspection). BUT the test NAME over-promises ("is_structurally_
  independent") for what it checks. Either rename to `..._builds_without_ld_scheme`
  or make it assert the import (`importlib`-introspect the module's references). NIT.
- **NIT-3 (`build_materials` verbatim duplication — the dispatch's Q3).**
  `SN2DCartesianLDStressMMSCase.build_materials` is BYTE-DUPLICATE of
  `SN2DCartesian2GHeterogeneousMMSCase.build_materials` (sn.py:863-896 vs the new
  copy) — the per-cell `Mixture` assembly loop (`SigC=sig_a, SigL=0, ... SigS=[csr],
  Sig2, chi=0`). This is MMS-INFRASTRUCTURE boilerplate, NOT a reference-independence
  requirement (the materials build is shared mechanism; the MMS independence lives in
  `external_source`/`phi_exact`, which ARE distinct). Per Cardinal-Rule-2 +
  unify-after-TWO: there are now exactly TWO identical copies → the trigger fires.
  Elegant fix: extract `_build_per_cell_hetero_materials(mesh, sigma_t_fn, sigma_s_fn,
  n_groups) -> dict[int,Mixture]` (a free helper next to `_default_hetero_2d_xs_
  functions`); both cases call it. Bug habitat: a future `Mixture` field addition (a
  P1 `SigS=[S0,S1]` row, or a fission `SigF`) lands on one copy, the twin keeps the
  old shape → silent divergence between the DD-MMS and LD-MMS references. `build_mesh`
  is NOT a twin (the LD case adds non-square `ny` default — a genuine difference).

## Recurring tells confirmed this review (for AGENT.md institutional knowledge)
- The MMS-case family has a STANDING twin habitat in `build_materials` (the per-cell
  Mixture loop). Every new 2-D MMS case copies it. Two copies now = extract.
- algebra-of-record Branch-2 cross-checks tend to pin only the IDENTITY-SCALING group
  (`c_0=1`) + one cell — the group-axis + cell-axis coverage is the blind spot. Demand
  the cross-check parametrize over groups when `c_g≠1` scaling is in the numpy path.
