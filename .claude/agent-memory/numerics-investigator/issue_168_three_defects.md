---
name: Issue #168 SN curvilinear MMS — three independent boundary defects
description: Issue #168 attributed the ~1.26 order MMS gap to a single bug (cell-center used as outer-face flux); empirical investigation found THREE independent O(1) boundary truncation defects in the symmetric-closure FD operator. Option A as the issue proposes it does NOT close ERR-026; Option C (boundary-face-flux Protocol + sphere-pole correction) is the only mechanically correct path.
type: project
---

# Issue #168 — three independent boundary truncation defects

## Why this memo exists

Issue #168 (open as of 2026-05-09) frames the curvilinear MMS
convergence gap as a **single bug**: cell-center used as face-flux at
`i = nx-1` for outgoing μ. Empirical investigation refutes the
single-bug framing — there are at least **three independent O(1)
boundary truncation defects** in the symmetric-closure
`transport_operator_matvec_spherical` / `transport_operator_matvec_cylindrical`.

**Why:** future implementers (and this same investigator on a
fresh-context return) need to know that Options A and B from the
issue are insufficient — fixing only Defect 1 produces an apparent
improvement at fixed nc (1.91 → 1.42 L2 error at nc=10) but the order
**degrades** with refinement (1.23 → 1.11 → 1.05 across nc 10→20→40→80),
because Defects 2 and 3 dominate as h→0.

**How to apply:** before any Issue #168 fix dispatch, read the
design memo at `.claude/plans/issue_168_design.md` and the
diagnostic scripts under `scratch/derivations/diagnostics/diag_issue168_*.py`.
DO NOT ship Option A alone — it would close GH #168 with a fix that
regresses the very property it claims to deliver.

## The three defects (empirical signatures from
`diag_issue168_05_apply_vs_solve.py` truncation-residual table)

For ψ_exact = sin(πr/R)/W (isotropic ansatz), nc-refinement of
``||L · ψ_exact - Σ_s/W·φ - Q/W||`` per cell:

| nc | inner i=0 | outer i=N-1 | next i=N-2 | mid-domain |
|----|-----------|-------------|-------------|------------|
| 10 | 0.117     | 0.152       | (varies)    | 7.4e-2     |
| 80 | 0.135     | 0.137       | 6.8e-2      | 1.4e-4 (O(h²) ✓) |

- **Defect 1 (outer cell, μ > 0)**: `psi_right = fi[N-1]` is
  cell-center, not face value at r=R. O(h) truncation. Fixable by DD
  extrapolation: `psi_right = 1.5·fi[N-1] - 0.5·fi[N-2]`.
- **Defect 2 (next-to-outer cell, μ < 0)**: BC fill at i=N-1
  overwrites `fi[N-1, n_inward, 0]` with the BC face value (=0 for
  vacuum). Matvec at i=N-2 then uses that 0 as a "cell-center" in
  the symmetric face average `0.5*(fi[N-2] + fi[N-1])` → corrupts
  the interior face flux at the boundary-adjacent cell. **Root cause
  is structural**: the array `fi[..., -1, ...]` conflates cell-center
  storage with BC face-value storage. Pre-Round 3 hardcoded reflective
  BC accidentally wrote a faithful cell-center value back; vacuum BC
  exposed the conflation.
- **Defect 3 (sphere pole, i=0)**: Bailey 2009 ΔA/w redistribution
  `-μ ΔA[0] · ψ_0/V[0] = -(3μ/h)·ψ_0` cancels the streaming for FLAT
  ψ but overcorrects for ψ varying linearly near r=0, giving the
  discrete operator HALF the correct streaming at i=0. Empirically:
  `(3μ/2)·ψ'(h/2)` discrete vs `3μ·ψ'(0)` continuum.

## Cartesian path: not affected

`transport_operator_matvec` uses upwind FD (different stencil), not
symmetric closure. The bug is curvilinear-only. Fix should NOT
generalize to Cartesian.

## What the issue's three options actually fix

| Option | Defect 1 | Defect 2 | Defect 3 | Empirically gives O(h²)? |
|--------|----------|----------|----------|---------------------------|
| A — DD extrap at outer face | ✓ | ✗ | ✗ | NO (orders 1.23→1.05) |
| B — ghost cell             | ✓ | ✗ | ✗ | NO (same as A) |
| C — DD throughout, with separate cell/face storage | ✓ | ✓ | partial | YES if Phase B sphere-pole correction added |

## Recommended implementation path

**Phase A** (~250 LOC):
- Add `BoundaryFaceFlux` Protocol (parallel to Wave C `CellUpdate`).
- `DDExtrapolation` strategy: `psi_face_out = 2·psi_cell - psi_face_in`.
- Wire through `transport_operator_matvec_spherical/_cylindrical`.
- Separate cell-center from face-value storage in
  `solution_to_angular_flux_*` — DO NOT overwrite `fi[N-1]` with
  BC values; thread the BC face-value as a separate argument to the
  matvec.
- Regenerate 6 curvilinear regression snapshots (sphere/cyl
  homogeneous + 3-region + P1 aniso). 5 Cartesian snapshots stay
  bit-identical.
- Phase A alone gives orders ~1.5-1.7 (Defect 3 still limiting).

**Phase B** (~100 LOC, needs literature):
- Sphere-pole correction at i=0. Lewis & Miller §4.5; Bailey 2009;
  Carlson starting-direction; Lathrop pole stencil.
- After Phase B: orders > 1.9, all 4 xfail markers off, ERR-026 fully
  closed.

## Risk: Defect 3 may not have a clean symmetric-closure fix

The sphere pole's factor-of-2 mismatch in Bailey ΔA/w may be a
fundamental limitation of the symmetric-closure FD approach. If
Phase B literature search reveals no clean fix, escalate: either ship
Phase A with relaxed xfail (orders > 1.5), or rebuild the apply
method to mirror the WDD sweep math entirely (most expensive option).

## Diagnostic scripts and design memo

- `.claude/plans/issue_168_design.md` — full design memo.
- `scratch/derivations/diagnostics/diag_issue168_01_characterize.py` — reproduces ~1.26 order.
- `scratch/derivations/diagnostics/diag_issue168_02_option_a_dd_extrap.py` — falsifies single-bug framing.
- `scratch/derivations/diagnostics/diag_issue168_05_apply_vs_solve.py` — truncation residual scaling.
- `scratch/derivations/diagnostics/diag_issue168_08_unit_residual.py` — per-cell structure of all three defects.

## Related work

- ERR-026 in `.claude/skills/vv-principles/error_catalog.md` —
  status currently PARTIAL CLOSURE; Issue #168 fix would change to
  CLOSED (Phase B) or remain PARTIAL with updated narrative
  (Phase A only).
- Wave E Round 3 closeout at
  `.claude/agent-memory/method-implementer/sn_reshape_wave_e_round_3_closeout.md`.
- The Wave E orchestrator's framing — that the BC plumbing fix
  alone would close ERR-026 on MMS — was empirically falsified
  during Round 3. Issue #168 inherited that falsified framing.
