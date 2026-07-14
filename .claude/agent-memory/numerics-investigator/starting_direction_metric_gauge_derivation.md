---
name: starting-direction-metric-gauge-derivation
description: Derivation-of-record + reusable recipe — the SN curvilinear starting-direction (ψ½) block Hilbert metric G_sd is GAUGE-FREE (any SPD), because apply_transpose is the EXACT Euclidean transpose (T=Aᵀ); the ghost G_sd=0 is the one forbidden value (breaks the adjoint for nonzero seed). Recommend V_cell. ψ½ is BULK-like (transport self-block), NOT trace-like; the plan's face_streaming_normal is an OPERATOR coefficient mis-installed as the STATE metric. Backs L17.
metadata:
  type: reference
---

# Deriving the SN starting-direction (ψ½) block metric G_sd (#282/#280 2.5d)

**Question posed:** derive, from the adjoint constraint `Aᵀ G = G A†`, the correct
Hilbert state metric `G_sd` for the curvilinear starting-direction block of the augmented
loss composite `A` on `FullField = bulk(AngularFlux) ⊕ trace(AngularBoundaryFlux) ⊕
SD(StartingDirectionFlux)`, so the G-adjoint reciprocity gate stops being Mode-12 blind to
seed-row errors. Ghost metric `G_sd≡0` (`starting_direction_space.py:for_levels` builds
`np.zeros`) justified as `(1−µ²)|_{µ=±1}=0`.

## Verdict (measured, not guessed)

- **`apply_transpose` T == plain Euclidean transpose Aᵀ EXACTLY** (‖T−Aᵀ‖=3.6e-16 global,
  every block incl. `T_sb==A_bsᵀ`). `A.H=G⁺·apply_transpose·G` (operator.py:1146) is therefore
  the honest metric-adjoint `G⁻¹AᵀG` for ANY invertible G.
- ⟹ **reciprocity `⟨Aψ,φ⟩_G=⟨ψ,A.Hφ⟩_G` is GAUGE-FREE in G_sd**: it holds for EVERY SPD
  `G_sd`. The determining equation pins G_sd ONLY up to "non-degenerate" (production diagonal
  architecture: any strictly-positive diagonal). NOT uniquely determined.
- **`G_sd=0` is the single forbidden value** — not SPD, puts the seed in null(G), makes A.H
  SEVER the seed block (`H[seed,:]=H[:,seed]=0`). Sharper than "blind": with nonzero seed data
  the shipped A.H BREAKS reciprocity at 1.3e-2 (production path) — a WRONG adjoint the instant
  the seed carries data (the `A_bs` seed→bulk coupling is unmatched). Green only because the
  gate feeds a present-but-ZERO seed.
- **Recommended fixing (physical, gauge): `G_sd = V_cell`** — per (level,sign): cells (ng,nx)
  = radial cell volume `V_cell` group-broadcast; corner (ng,) = any positive boundary measure
  (outer-cell V, or the r=R shell area — gauge). Matches the bulk convention `G_bulk=V_cell·w_n`
  (measured g_bulk∈[1.46,101]; `full_field_space.py:170`, `_bulk_measure`=`w·V`). An angular
  weight `w` may be folded in (V·w) but is neither required nor canonical (no natural quadrature
  weight for the single µ=±1 ray).

## The tension resolved

`(1−µ²)|_{µ=±1}=0` is the angular-INTEGRATION weight of the pole ray (correctly zero — the seed
is absent from `φ=Σ_n w_n ψ_n`). It is IRRELEVANT to the STATE metric. ψ½ is a first-class radial
state field (nonzero self-block `A_ss`‖4.0‖, nonzero coupling `A_bs`‖6.0‖), whose Hilbert metric
is set by its role in the operator algebra (must be SPD so the adjoint is non-degenerate). The
radial-field VOLUME makes G_sd nonzero.

## Block structure (measured, order seed≺bulk≺trace, A block-lower-triangular)

`A_bs≠0` (seed→bulk M-M recurrence), `A_sb=A_st=0` (seed rows self-contained), `A_ss≠0`.
`A.H=G⁻¹AᵀG` is block-UPPER-triangular, seed at the TOP ⟹ bulk⊕trace rows of A.H are BITWISE
gauge-invariant (Δ=0.0 exact under identity/V_cell/10·V_cell); only φ†_seed moves. No physical
observable (bulk φ†, adjoint reaction rates) reads G_sd — the gauge is genuinely free.

## Reusable recipe (any augmented-composite block Hilbert metric)

1. Assemble A (`A.apply`), T (`A.apply_transpose`), and — cross-check — A.H (`A.H.apply`) as
   dense matrices by unit-vector probing over `FullField.to_flat` (layout `[bulk.ravel, trace,
   seed]`). Partition indices from the block sizes.
2. **Check `T==Aᵀ` against numpy's `Afwd.T`** (structurally-INDEPENDENT ground — NOT the
   operator's own metric machinery; avoids reference contamination). Exact ⟹ metric gauge-free
   (need SPD). Weighted ⟹ pinned by `g_s[i]=g_b[coupled bulk j]` on the coupling blocks.
3. Faithfulness: confirm `G⁺TG (dense) == A.H (production)` so the dense analysis IS production.
4. Sweep candidate `g_s ∈ {0, 1, V, V·w}`; measure reciprocity defect for (a) ZERO block data
   [the current gate regime — all pass, blind], (b) RANDOM block data [0 BREAKS, SPD holds],
   (c) RANDOM + a self-block SIGN FLIP [SPD REDS = Mode-12 closes, 0 stays blind].
5. Forward stays bit-identical under the install (metric read only by A.H + inner_product;
   monkeypatch `space.inner_product_weights`, assert `A.apply` unchanged — #208 trace precedent).

## Install caveats for the main agent (surgical)

- Replace the `np.zeros` ghost in `StartingDirectionSpace.for_levels` with V_cell (SPD).
- The Mode-12/reciprocity gates MUST also switch to NONZERO seed probes (`_random_composite`
  `seed_block="random"`; the Mode-12 test's `_zero_seed_composite`) — a zero-block probe can't
  see the fix. `test_mode12_..._blind_to_a_seed_row_flip` INVERTS: assert the flip now REDS.
- Rewrite the docstrings that document the blindness as correct: `starting_direction_space.py`
  "ghost metric — ALL weights zero (structural)" module block + the A4 honest-scope note, and
  `test_g_adjoint_reciprocity.py::_random_composite` / the Mode-12 docstring.

**PROMOTED + fix LANDED (2026-07-06).** `G_sd=V_cell` is live (`starting_direction_space.for_levels`
builds `np.tile(cell_volumes,ng)`, strict-positive guard; `augmented_mesh:858` passes
`cell_volumes=self.volumes`). The 5 `diag_gsd_*` diagnostics are DELETED, consolidated into the
permanent gate **`tests/sn/operators/test_starting_direction_metric.py`** (`pytestmark=foundation`,
19 gates green under `-O`; the Mode-12 gate `test_derive_gsd_and_close_mode12` carries
`@catches("ERR-067")`). ERR-067 is in the catalog (§4589). **Coverage split (mutation-verified,
non-obvious):** reverting production `G_sd→0` REDs the three PRODUCTION-PATH value-catchers
(`test_production_path_vcell_honest_adjoint`, `test_shipped_metric_block_values`,
`test_trace_and_pole_faces_both_hold_reciprocity`) but LEAVES the dense `test_derive_gsd_and_close_mode12`
GREEN — the dense gate is structurally self-contained (builds its own {0,id,V_cell,V·w} candidates,
faithfulness robust), so it catches a broken CLOSURE MECHANISM (SPD-flip-reds + ghost-blind contrast +
the `<1e-12` control leg), NOT the shipped value; the production-path trio catches the value. Together
they cover ERR-067. R4 must use a freshly-computed `2·V_cell` perturbation (NOT `2×`shipped) or it
spuriously reds on the ghost via the non-vacuous-sanity check while its real forward-blindness claim
holds. Artefacts (retired): `diag_gsd_0{1..5}`. Backs [[lessons-L17]].

## FaceField-plan dialectic (2026-07-06) — two NEW probes, the metric is a per-leaf property

`facefield_codim1_design.md` parents ψ½ under a `FaceField` ABC as a sibling of the
spatial boundary trace, unifying BOTH metrics into one `face_streaming_normal`
measure: spatial `|Ω·n|·w` / pole `(1−µ²)|_pole·w ≡ 0` (§3.2, §7 Phase A gate demands
`G_sd==0 at 0 ULP`). Two probes REFUTE the metric-unification (structure-sharing
survives):

- **diag_gsd_04 (trace-vs-pole asymmetry):** under the SHIPPED metric, activating the
  SPATIAL face keeps reciprocity at 1.6e-16; activating the POLE face breaks it to
  9.4e-1 — a **5.8e15 asymmetry from ONE measure**. Shipped blocks: bulk V·w [1.46,101],
  trace |Ω·n|w [0.22,0.30] NONZERO, seed [0,0]. `face_streaming_normal` is CORRECT for
  the spatial face, WRONG for the pole.
- **diag_gsd_05 (seed self-block is transport):** A block-norms (nx=6): bulk A_bb offdiag
  89.3, trace A_tt offdiag **2.0** (diag∈[−1,1], a reflection/restriction map), seed A_ss
  offdiag **71.2** (diag∈[−1,3.65], a BANDED radial streaming+collision operator). Ratio
  A_ss/A_tt = 35.6. ψ½ has genuine INTERIOR radial dynamics (bulk-like); the trace is a
  pure restriction (no interior dynamics). A_sb=A_st=0 (triangular), A_bs≠0 (feeds bulk).

**The category error (crisp):** `face_streaming_normal` = the streaming-NORMAL / angular-
redistribution coefficient that legitimately lives INSIDE the operator A. It equals the
STATE metric G ONLY when the face's self-block A_ff is trivial (trace: A_tt=reflection ⟹
through-flux ≡ state-metric ≡ partial current). ψ½'s A_ss is a full transport operator ⟹
through-flux (0) ≠ state-metric (V_cell). The plan installs an operator coefficient as the
Hilbert metric. **The user's "reduced metric in ANGLE" harmonization, read as
MARGINALIZE-the-angle (keep the radial volume), GIVES V_cell and refutes the plan; read as
MEASURE-the-angular-face (through-flux) it gives 0. The plan conflated the two readings.**
Verdict = FaceField ABC owns STRUCTURE only (FaceLayout + presence-invariant); metric stays
per-leaf (spatial→partial current, pole→V_cell), same as the bulk's metric is per-leaf V·w.
Diagnostics `diag_gsd_0{4,5}_*.py` were the diagnostics-of-record; the SPD fix LANDED and both
are now PROMOTED/inverted into `tests/sn/operators/test_starting_direction_metric.py`
(`test_seed_selfblock_is_transport_not_reflection` as-is; the trace-vs-pole probe INVERTED to
`test_trace_and_pole_faces_both_hold_reciprocity` + `test_shipped_metric_block_values`). All 5 diag
files deleted. **Doc-drift flagged (out of numerics-investigator scope):** the UPPER module docstring of
`starting_direction_space.py` still says "ghost metric — ALL weights zero", contradicting the LOWER
section (V_cell, Mode-12 CLOSED) and the live `for_levels` — the archivist should reconcile.
