---
name: issue-240-phase2-step-d4-closeout
description: #240 Phase 2 Step D4 — the diffusion-readiness contract gate locking the DiscretizationScheme as a Σ-stateless generic advection–reaction interface + the Base interface note (gate + docstring only; bit-identical)
metadata:
  type: project
---

# #240 Phase 2 Step D4 — Σ-stateless advection–reaction contract lock

Branch `feature/sn-space-angle-tier2` (host env, `.venv/bin/python`, `-O`).
NOT committed (user commits after qa). Plan: `.claude/plans/issue_240_phase2_step_d_homing.md` §D4.

**Why:** D4 had no production-numerics work to do — the `DiscretizationScheme`
was ALREADY Σ-stateless (verified: `DiamondDifference()` and
`LinearDiscontinuous()` are stateless dataclasses, `vars()` empty; every
coefficient method takes Σ as an explicit param `total_xs`/`sigt_cells`/`sig_t`;
`theta` is a `ClassVar`, not Σ-state). So "realize the generic advection–reaction
form" was structurally DONE by Steps A–D3. D4 = LOCK it in (a contract gate) +
make the interface explicit (a Base docstring note). NO numerical change.

## Deliverables

1. **Contract gate** — `tests/sn/spatial/test_scheme_reaction_rate_contract.py`
   (NEW, 10 tests, all `@pytest.mark.foundation`, NO `verifies` — a
   software-contract gate, not equation verification). Proves:
   - **POSITIVE (closed-form at arbitrary reaction-rate):** DD + LD
     `affine_scan_coefficients` called with an arbitrary `sig_t` (a removal-like
     `Σ_r = 0.35`, a pure-advection `Σ → 1e-13`) match the closed form AT THAT Σ.
     DD → `w=½`, `denom = Σ·V + 2|μ|·A_down`, `a = 2|μ|·A_total/denom − 1`.
     LD → `w = 1/(1+k)`, `k = (|μ|/θ)/(Σ·V + |μ|/θ)`, `inverse_denom = 1/S`,
     `a = m(1+k)²/S − k` (slab-neutral geom: A_down=1, A_total=2, dA_w=0, c_out=0).
     rtol=1e-14 (discriminates a 1e-12 perturbation — verified).
   - **Péclet limits:** LD `w → ½` as `Σ → 0` (central/advection limit),
     `w → 1` as `Σ → 1e9` (full-upwind limit), strictly monotone between (the
     κ-scheme blend). DD `w = ½` Σ-INDEPENDENT (central / Keller box, verified
     at Σ ∈ {1e-13, 0.35, 1e9}).
   - **NEGATIVE / statelessness teeth (the lock-in):** `vars(scheme)` empty
     (no `self.<sigma>` smuggle — verified the teeth WOULD fire if one were
     stored via `object.__setattr__`); a reused instance called twice with
     DIFFERENT Σ returns coefficients that (a) DIFFER in `inverse_denom`
     (Σ-responsive) and (b) each match the closed form at their OWN Σ (no hidden
     Σ memory); re-call at Σ_a after Σ_b reproduces bit-for-bit
     (`assert_array_equal`).
   - **Mode-8 safe:** every assertion is `np.testing.assert_*` / `pytest.fail`
     (no bare `assert` — fires under `-O`). 10/10 pass under `-O`.

2. **Interface note** — `orpheus/sn/spatial/scheme.py`, `DiscretizationSchemeBase`
   docstring, NEW "Generic advection–reaction interface (Σ-stateless — #240
   Phase 2 Step D4)" section (~26 lines incl. math block): coefficient methods
   parameterized by (streaming/wave-speed, reaction-rate Σ, source, geometry),
   hold NO Σ state; the model-agnostic advection–reaction discretization the
   diffusion solver consumes (standalone + consistent-DSA #2); Step↔upwind /
   DD↔central-box / LD↔DG-P1, `w(Σ)` Péclet/κ blend; pointer to D6 narrative on
   `:doc:/theory/discrete_ordinates`; deferred #241 rename note (`total_xs`→
   `reaction_xs`, `cell_average_weight`→`face_blend_weight`; `streaming` KEPT).

## Verification (bit-identical — gate + docstring only)

- New gate: **10 passed** under `-O`.
- Strict DriftWarning gate `tests/sn/sweep/core tests/sn/solve
  -W error::...DriftWarning` → **505 passed / 1 skipped / 4 xfailed** (= baseline 505/1/4).
- Route-around `tests/sn/operators spatial sweep/core sweep/cartesian_2d solve`
  with the standard `-k` route-around → **1093 passed** (= 1083 baseline + 10 new),
  6 skipped / 7 deselected / 5 xfailed. No DriftWarnings, no regressions.
- Sphinx `-W --keep-going` → **exit 0**, ZERO warnings from `scheme.py` /
  `DiscretizationScheme` (the `:meth:`/`:mod:`/`:doc:`/math note is clean). The
  only build items are pre-existing (`mesh.py` `paramref` ERROR; unrelated
  `SyntaxWarning`s in `test_projection_operators.py` / `test_fission_operator.py`).
- `tests._harness.audit` → **exit 0**.

## Working-tree note

`docs/verification/matrix.rst` shows a build-regenerated diff (foundation
3076→3086, new `spatial/test_scheme_reaction_rate_contract, 10` row) — the
CORRECT auto-regeneration reflecting the 10 new foundation tests, NOT a
hand-edit (file header: regenerated on every build). Left in place as a coherent
unit with the new gate. `.claude/agent-memory/cross-domain-attacker/*` changes
are pre-existing (the naming-signal memo), not mine.

## Remaining D-series

D5 (#239 2-D Cartesian ScanMarch onto the coefficient model — needs a 2-D
coefficient cache + transverse `s_y·ψy` fold; ~1-ULP re-baseline; MAY split off),
D6 (archivist: Step/DD/LD ↔ advection-scheme rich narrative + the Spatial ×
Angular tensor-product next-campaign issue). #241 holds the deferred
model-agnostic param rename. D6 DISPATCH_REQUEST to archivist owed once D5 lands
(or now — the Base docstring note is the stub D6 expands).
