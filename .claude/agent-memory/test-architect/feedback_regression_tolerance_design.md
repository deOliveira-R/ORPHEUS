---
name: regression-tolerance-design
description: Principled regression-snapshot tolerance — SAFETY×conv_tol gate + DriftWarning tripwire + -O-safe helper; read conv_tol off the run config (SoT), never hardcode
metadata:
  type: feedback
---

A frozen-snapshot regression gate MUST NOT carry hand-picked magic
tolerance floors (`rtol=1e-12`, `rtol=5e-6`). Replace them with a
single-source-of-truth helper (`assert_regression`) carrying two layers:

1. **Correctness gate (hard fail)** tied to what the computation
   PROMISES. **Iterative** results (k_eff / flux from power-iteration or
   source-iteration) → `SAFETY × conv_tol`, where `conv_tol` is the
   solver's OWN convergence stopping criterion for the pinned quantity —
   `keff_tol` for k_eff, `flux_tol` for the eigen flux, `inner_tol` for
   the fixed-source flux. `SAFETY=10` is the iteration-map amplification
   headroom `ρ/(1−ρ)` for `ρ≲0.9` (a principled bound, NOT a fudge).
   **Direct** results (single sweep, no outer iteration) →
   `assert_array_almost_equal_nulp(nulp=reduction_depth)`.
2. **Drift tripwire (informational)**: `DriftWarning(UserWarning)` when a
   value clears the gate but moved beyond bit-identity (ULP>0).
   Escalatable via `-W error::<module>.DriftWarning` = the strict
   bit-identity gate for a pure-refactor PR. `logging.debug` carries the
   per-element ULP forensic breakdown.

**Why:** ORPHEUS SN regression redesign 2026-06-01. The magic floors
encoded no claim about what the solver promised; the slab floor was too
tight for FP-non-associativity (~1e-11) and the curvilinear band papered
over a since-fixed bug. The principled gate IS the claim. Implemented at
`tests/sn/regression/_regression_assert.py`.

**How to apply (load-bearing details):**
- Read `conv_tol` off the run config (the dict that drove the solve),
  NEVER hardcode it in the test. Make `run_config(cfg)` (defaults +
  per-case override) the SINGLE source of truth shared by the generator
  and the test, so they cannot drift on what "converged" means
  (coding-elegance Pattern 7). Map quantity→config-key explicitly
  (k_eff→keff_tol, eigen-flux→flux_tol, fs-flux→inner_tol).
- The helper MUST be `-O`-safe: `np.testing.*` (raises unconditionally)
  + explicit `raise`, NEVER a bare `assert` ([[default-test-mode-is-O]]
  / vv-principles Mode 8 — `-O` strips bare asserts in non-rewritten
  helper modules). DEMONSTRATE it: a deliberately-wrong expected value
  must still FAIL under `-O`; a within-tol value must PASS.
- Register the DriftWarning filter as `always` in the local conftest so
  pytest surfaces every occurrence in the summary (not once-per-location).
- DEMONSTRATE both tripwire modes: informational (test passes, warning
  in summary on a +3 ULP perturbation) AND escalated (`-W error::...`
  → hard FAIL on the same perturbation). The demo perturbation is how you
  prove the tripwire works.

**⭐ CARVE-VERIFICATION COROLLARY (#236 §2, 2026-06-13): an ITERATED
end-to-end snapshot CANNOT serve as the bit-identity gate for a
"zero-numerical-change" refactor — descend to a single-step DIRECT
snapshot.** MEASURED on clean `main` (no carve): the committed iterated DD
snapshots (`test_dd_regression.py`) ALREADY drift under `-W
error::DriftWarning` (e.g. `cyl_1g_homogeneous` 272005 ULP / 3.56e-11 rel;
`sphere_2g_homog` 30580 ULP) — cross-run/cross-machine iterative FP jitter
against the frozen `.npz`. So the prior "DriftWarning is dormant in the
real suite" claim is FALSE for iterated cases (it holds only same-run
same-machine). CONSEQUENCE for a carve that must prove DD didn't move: the
`SAFETY×conv_tol` iterated gate stays green either way (too loose to
witness bit-identity), and `-W error::DriftWarning` is noisy. The
bit-identity PROOF must live where there is NO outer iteration — a
single-sweep / single-matvec snapshot on a fixed-seed RANDOM ψ (het, ≥2g,
non-zero inflow → activates the curvilinear redistribution; flat ψ nulls
it), captured pre-carve via the root-conftest `--capture-baseline` flag
(`tests/conftest.py:42`, the Wave-O `test_bc_extraction_matvec.py`
mechanism), asserted `np.testing.assert_array_equal` (rtol=0) post-carve;
the only legitimate drift is `reduction_depth × ULP` (→ `kind="direct"`,
`nulp=nx`, WITH the 3-criteria note if the fold order shifts). Pair the
snapshot ("did not move a bit") with a closed-form `(diagΣ_t−Σ_s0ᵀ)⁻¹Q`
anchor ("value is correct") — vv §bit-identity criterion 2.

See [[eigen-on-nonfissile-mixture-is-malformed]] for the snapshot-regen
guard that surfaced alongside this; [[snapshot-migration-when-production-goes-bare]]
for the broader "schema = persisted∩compared + vacuum-bit-id correctness
gate" recipe.
