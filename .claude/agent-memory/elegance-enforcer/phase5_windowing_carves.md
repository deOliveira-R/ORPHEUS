---
name: phase5-windowing-carves
description: Phase 5 SN angular/face/moment windowing carves (5a/5b/5c) + the WavefrontFlux mint — durable facts on the rolling-window cochain, the in-sweep moment carve (first 5-series twin-collapse), and the reference-oracle retirement-exception. All LANDED.
metadata:
  type: project
---

Durable facts from the Phase 5 SN memory-windowing carves (5a angular → 5b face → 5c moment)
plus the `WavefrontFlux`/`InteriorFaceSpace` Phase-1 mint they build on. All LANDED on
`refactor/field-role-typing` (`93807aa`..`63719a2`). Verdicts were PASS-WITH-NITS throughout;
the nits are landed history — these are the reusable structural lessons.

## The arc (what each carve did, structurally)
- **5a** holds the 2-D Cartesian within-group SI iterate as `HarmonicMomentField` (N→(L+1)²)
  instead of full per-ordinate `AngularFlux`, via `_MomentWindowedResolvent` (gated
  `sn_mesh.reduced is None`; curvilinear + Krylov stay full-angular). Commit-1 win:
  `_aniso_source_from_moment_values` single-sourced the `R·Λ` moment→source assembly, retiring
  the per-ℓ `_PerLegendreOrderScattering` kernel.
- **5b** replaces full per-axis interior face backing `O(N·ng·nx·ny)` with a rolling 2-diagonal
  `_MovingFrontier` cochain `O(N·ng·nx)` for BOTH the sweep and the matvec twin. ⭐ design-space
  win: the window's contiguous slots → basic-slice zero-copy VIEWS where the full-field grid
  diagonal forced fancy copies ⟹ a SPEEDUP (~0.77×), not a memory-vs-time trade. The cochain API
  (`seed_x`/`seed_y`/`incoming`/`emit`) reads like trace algebra, abstraction ~free.
- **5c** moves the angular→moment projection FROM a post-sweep flat `MomentProjection.apply` INTO
  the windowed anti-diagonal walk (`moment_buf += einsum` per level); `_MomentWindowedResolvent.solve`
  COLLAPSES to a 2-line forward to `base.solve_moments`. Measured 3.06× peak-memory win.

## Durable lesson 1 — 5c is the FIRST 5-series carve that RETIRES a twin (not adds one)
Pre-5c the resolvent did `base.solve` + flat `MomentProjection.apply` → moment arithmetic in TWO
places (sweep scalar fold + post-projection). 5c folds the post-projection into the sweep; concept
count DROPPED (the right-abstraction-shrinks-the-count test). Contrast 5a/5b which were optimizations
that kept a fuller view as an oracle. When reviewing a windowing carve, ask which kind it is:
twin-collapse (concept count drops) or fuller-view-relinquish (needs an oracle pin, below).

## Durable lesson 2 — the reference-oracle retirement-exception, done right (promoted AGENT.md #3)
5b keeps BOTH the full-field walk (`graph.apply`/`residual`, oracle) and the window walk
(`apply_windowed`/`residual_windowed`, production); 5c keeps the flat `MomentProjection.apply` as a
fuller-view oracle. LEGITIMATE because: ONE shared cell kernel (`DiamondDifference.cell_kernel_batch`/
`residual_kernel_batch` — Pattern 2, math cannot drift) AND a foundation `window≡full` /
`in-sweep≡post-projection` equivalence test pins it. The oracle test must PROBE THE CORNERS: the
full moment tensor incl. ℓ≥1 (the scalar ℓ=0 cross-check is blind to ℓ/m drift), the cross-octant
shed-capture across octants × het meshes. Principled-equiv (cross-octant `+=` reorders the ordinate
sum) is fine with a documented `4·N·eps` bound + a structural anchor (SI≡Krylov≡k_inf), and the
guard must `raise AssertionError` (not bare `assert`, which `-O` no-ops). Without such a pin, "keep
the reference impl" is the superseded-code-obscures-signal anti-pattern.

## Durable lesson 3 — DI-over-flag + the polymorphic-sink REJECTION
5c's moment mode is an optional OBJECT (`moment_projection=None` vs given, mirroring the existing
`reflect: Callable|None`); the public surface is NAMED methods `solve`/`solve_moments`, no boolean
flag. The author considered + REJECTED a uniform sink/accumulator Protocol: the angular path is
octant-local-buffer-then-scatter, the moment path is global cross-octant `+=` — a REAL asymmetry
(it IS what makes moments principled-equiv-not-bit-id). A uniform sink would add Protocol+2 classes
for a 2-line output-branch delta. Pattern-6 ENDORSED: two instances DELIBERATELY not aligned. Use
this as the template for "when NOT to extract a Protocol from 2 sites."

## Durable lesson 4 — the FaceField-ABC unify-after-two deferral (the mint)
`BoundaryField` (instance 1) and `WavefrontFlux` (instance 2) overlap heavily (flat buffer,
FaceLayout-on-space, `_check_partner` mesh guard, `zeros_on`/`from_mesh`) BUT the face KEY differs
structurally — boundary = face-name STRING (`"xmin"`), interior = axis-index INT (`face(0)`). Lifting
`FaceField(Field)` now would abstract over {string-keyed vs int-keyed} = abstracting over the
DIFFERENCE = the unify-after-two trap. Correct to defer until both concretes + the sweep/matvec
consumers show the real shared surface. Bound the duplication; point each copy at the deferred lift.
(Note: WavefrontFlux briefly appeared orphaned after 5b's window bypassed it; 5b commit 7 `60e6eeb`
recovered the full-field oracle to consume it — verified still in production use: `sweep.py`,
`operator.py`, `sweep_graph.py`. The orphan was resolved, not retired.)

## Recurring nits caught (all promoted to AGENT.md #7)
- Source-arg asymmetry: angular arm passes a named typed `AngularSourceSink.from_mesh(...)`,
  moment arm passes RAW `.values` — same payload today, but the moment arm bypasses the typed
  envelope (a future `/W` re-home would hit one arm only). Route both through the same named source.
- Aliased return slot: `(moment_buf, moment_buf[0,0])` where the 2nd slot is a LIVE VIEW (angular
  mode returns an independent scalar) → return `(moment_buf, None)` for the discarded slot.
- 2-site duplication of the windowed-gate + moment-cold-start across eigenvalue + fixed-source SI
  → factored to `_maybe_window` / `_windowed_cold_start` (the `windowed` bool must escape).
- Reaching past a typed primitive to a raw ndarray: `l_block(0)[0]` instead of `.scalar_flux()`
  (which carries the Y₀⁰=1 identity in one place).
- seed/absorb byte-identical-but-direction twin in the mint → factor `_edge_slot(name)` so the
  edge-locating slice math is single-sourced.

Related: [[wave-o-operator-algebra]] (the typed source/sink these emit),
[[field-role-and-block-role-attrs]].
