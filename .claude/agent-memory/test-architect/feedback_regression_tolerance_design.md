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
  → hard FAIL on the same perturbation). A freshly-regenerated snapshot
  is bit-identical to the current solver, so DriftWarning is dormant in
  the real suite — the drift only appears cross-run / cross-machine; the
  demo perturbation is how you prove the tripwire works.

See [[eigen-on-nonfissile-mixture-is-malformed]] for the snapshot-regen
guard that surfaced alongside this.
