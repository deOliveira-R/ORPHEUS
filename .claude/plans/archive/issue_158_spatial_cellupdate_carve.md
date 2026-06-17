# Tier 2a — LinearDiscontinuous as a FULL selectable discretization protocol (#158)

> **STATUS (2026-06-15):** Step 1 + **Increment A committed LOCAL** (`2b56348`/`30aadb9`/`3f0cbc2`)
> + **Increment B (the COEFFICIENT MODEL) CODE DONE + VERIFIED + REVIEWED** (awaiting user "commit").
> on `feature/sn-space-angle-tier2`. **Authoritative B plan + prototype = `mellow-swinging-breeze.md`**
> (plan-mode file). B = `affine_scan_coefficients`→`(a, inverse_denom, w)`; closure methods retired →
> generic `affine_closure.py`; matvec APPLY via the `matvec_via_kernel` capability (LD→`residual_kernel_batch`,
> DD→`cell_balance` byte-id); LD `is_affine_scannable=True` + slab guard → LD rides `CumprodScan`.
> DD BYTE-IDENTICAL (505/1/4); two-paths + group3≡group2 + O(h²) + matvec≡sweep green; elegance PASS
> (×2) + qa SUPPORTED; #239 filed (2-D ScanMarch lift). NEXT after B-commit: Increment C (diffusion
> limit) or Step-4 close-out (#36: theory page + SymPy derivation into `derivations/`).
> Parent: `.claude/plans/sn_space_angle_discretization_plan.md`.
>
> **⭐ INCREMENT A LANDED (the corrected architecture — user directive: NO twin DAG):**
> LD runs on `FullFieldWavefront` (the polymorphic any-d DAG oracle, incl. 1-D slab) via the
> group-2 kernel (`cell_kernel_batch`/`residual_kernel_batch`, ÷V Schur, 1-D, flat Q̂=0) — the
> SAME contract DD fills. NO new strategy (the `OneDimPerCellWalk` 1-D-DAG attempt was a twin
> path → REVERTED). `default_for` routes 1-D LD→FullFieldWavefront (DD→CumprodScan unchanged).
> `cell_update` threaded through `solve_sn_fixed_source`→`_as_sn_mesh`→`SNMesh`. The
> `SNSolver.__init__` scan-cache now built only for `is_affine_scannable` (DD; LD skips — the
> scan cache feeds the DAG-FREE scans only). Files: `linear_discontinuous.py` (group 2 +
> `_kernel_terms`), `loss_representation.py` (revert only), `solver.py` (cache guard + kwarg).
> Gates: `test_linear_discontinuous.py` (gate-1 occupant + `TestLDKernel` group-2 round-trip +
> group1≡group2) + `test_mms_ld_slab.py` (routing + O(h²) MMS + Krylov matvec≡sweep +
> `test_ld_thick_diffusive_limit_xfail` forward tripwire). **526 passed/1 skip/5 xfail** strict
> (505 DD-bit-id baseline + 21 LD). elegance PASS-w-NITS-fixed; qa SUPPORTED-w-fixes-applied.
> ⚠ Flat-source LD (Q̂=0) is O(h²) but LOSES the diffusion limit (~58% deficit thick-diffusive) —
> HONESTLY documented (class docstring `.. warning::`) + GATED (the xfail tripwire flips at Inc C).
> FOLLOW-UPS: stress-ansatz gate-3 (qa: the sin/1G/homog MMS under-stresses — the test-architect's
> `SNSlabLDStressMMSCase` is the strong gate); NIT-2 named per-axis carrier (Inc D).

## Why this REVISION (the original per-cell+flat-source plan was refuted)

Two pieces of evidence gathered this session + the user's full-protocol directive corrected the plan:

1. **LD IS affine-scannable** (decisive): `ψ_out = a·ψ_in + b`, `a` source-independent to machine
   ε; `ψ_out` exactly linear in `ψ_in`. The literature memo's "not affine-scannable" conflated
   *two cell unknowns* with *non-affine chain*. CATCH: LD's `ψ̄` is NOT a face-pair function
   (affine in `ψ_in` + the source) → DD's `cell_average_from_faces` doesn't fit; the scan needs
   a generalized cell-average emission (DD's `½(in+out)` the special case).
2. **Flat-source LD (Q̂=0) loses the diffusion limit** (decisive): O(h²) at c=0 (even k=7 stress),
   but thick diffusive (c=0.999) flat error 17–57% vs canonical 1–2%. The scattering-source slope
   `Σ_s·φ̂` (→ `φ̂` in the iterate, the "global moment-contract") is REQUIRED for the #233/#236
   payoff.

**Evidence scripts** (`$CLAUDE_JOB_DIR/tmp/`, promote during implementation): `ld_slab_derivation.py`
(SymPy 2×2 + linear-exactness oracle), `ld_slab_order_experiment.py` (full vs flat O(h²) c=0),
`ld_probe_affine_and_diffusion.py` (affine proof + diffusion-limit discriminator).

## The full protocol (the user's directive)

LD = a full selectable scheme like DD, working on ALL strategies. The `CellUpdate` contract
(`diamond.py:506-516`) has THREE capability groups:

| Group | Methods | Family | LD status |
|---|---|---|---|
| 1 per-cell ref | `update`/`residual` | per-cell DAG (1-D) — **basic qualification** | ✅ done (Step 1) |
| 2 batched kernel | `cell_kernel_batch`/`residual_kernel_batch` | multi-D wavefront DAG — **qualification (multi-D)** | future (2-D bilinear) |
| 3 affine-scan triple | `affine_scan_coefficients` + cell-avg emission | DAG-free scan — **performance bonus** (affine) | this plan (1-D) |

"DAG families = basic qualification; affine + DAG-free scan = performance bonus." **Do NOT retire
the per-cell/DAG path** (the wavefront needs it). `is_affine_scannable=True` (set in B).

## Staged increments (each landable: green + bit-id where applicable)

- **A — per-cell DAG qualification (1-D) + canonical fixed-source O(h²) [IN PROGRESS].**
  `OneDimPerCellWalk(_LossRepresentation)` (supports=`is_1d AND NOT is_affine_scannable`); per-cell
  `cell_update.update` walk (generalize the degenerate-cyl loop `_OneDimScanWalk._run:2681-2705`
  for non-affine schemes); thread the (2,ng) moment source; register in `LOSS_REPRESENTATIONS`;
  thread `cell_update=` through `solve_fixed_source` (`solver.py:1843/1971`). Gates: gate1 ✅,
  gate2 seam-control (DD-via-PerCell≡DD-scan nulp=nx), gate3 slab MMS O(h²) canonical Q̂,
  gate5 DD bit-id. LD stays `is_affine_scannable=False` in A. LOW-risk (no DD scan touch).
- **B — affine-scan bonus.** Generalize group 3 (scheme-provided source emission `b` + cell-avg
  emission affine-in-inflow; DD bit-id special case). `is_affine_scannable=True`; LD rides
  CumprodScan. Gates: CumprodScan-LD≡OneDimPerCellWalk-LD (nulp), gate3 via scan, DD strict bit-id.
- **C — diffusion limit.** `φ̂` in the within-group iterate → scattering source `(Q̄, Q̂=Σ_s·φ̂)`.
  Diffusion-limit gate (thick slab c→1; promote the probe Part 3). #233 pole-cell lift.
- **D — multi-D group 2 (future).** `cell_kernel_batch` for 2-D bilinear LD → wavefront.

## Verification
- Strict DD bit-id (every increment): `python -O -m pytest tests/sn/sweep/core tests/sn/solve
  -W "error::tests.sn.regression._regression_assert.DriftWarning"` (505/1skip/4xfail baseline).
- LD: gate1 ✅ (`tests/sn/spatial/test_linear_discontinuous.py`), gate3 slab MMS O(h²), two-paths
  gate (B), diffusion-limit gate (C), #233 char gate run-with-LD.
- `python -m tests._harness.audit` exit 0; Sphinx clean; elegance-enforcer + qa per increment.
- test-architect extends its spec with the two-paths + diffusion-limit gates before B/C.
- Route-around reds: #232, #195, #212, #237.

## Critical files
- `orpheus/sn/spatial/linear_discontinuous.py` (occupant; extend group 3 [B], group 2 [D]).
- `orpheus/sn/loss_representation.py` (`OneDimPerCellWalk` [A]; scan source-emission+cell-avg [B];
  `_OneDimScanWalk._run` per-cell [A]; `LOSS_REPRESENTATIONS`+`default_for`).
- `orpheus/sn/spatial/{diamond.py (ref), scan.py, cell_balance.py}` + `CollisionCache` [B].
- `orpheus/sn/solver.py:1843/1971` (thread `cell_update=` [A]).
- `docs/theory/discrete_ordinates.rst` (LD + 3-group protocol + affine/diffusion findings).

## DO NOT
- Don't retire the per-cell/DAG path; don't ship flat-source as the END state; don't break DD
  bit-id in the group-3 generalization (B); no method-implementer (surgical inline carve).
