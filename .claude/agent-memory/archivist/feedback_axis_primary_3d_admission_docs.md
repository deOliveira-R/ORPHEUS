---
name: axis-primary inversion + 3-D admission docs (C5 / Issue #225)
description: Documenting the SN N-D campaign's C5 — a constructor data-flow INVERSION (axes→mesh→axes round-trip retired, axes tuple now PRIMARY) plus axis-native 3-D admission with NO Mesh3D dataclass. The DIRECT sequel to the C4 face-name-carve docs ([[feedback_face_name_carve_docs]]). NEW C5 section in boundary_conditions.rst (the C4 predecessor's home), 6-substep arc (C5.1 inversion / C5.2 phantom retire / C5.3 geometry-blind trace / C5.4 vv-Mode-9 gate retarget / C5.5 admission+value-gates / C5.6 pin flips), Mode-9 gate-retarget framing, lossy-round-trip ASCII data-flow, d=3 value-gate evidence table mapped to V&V pillars, method-rename stale-ref sweep (from_mesh_and_quadrature→from_quadrature_and_layout) with LIVE-vs-HISTORICAL triage, and a #223-pointer note on a worked-example carved OUT of scope.
type: feedback
---

# Axis-primary inversion + 3-D admission docs (C5 / Issue #225)

Sequel to [[feedback_face_name_carve_docs]] (C4): C4 made the *boundary
keying* dimension-agnostic; C5 makes the *whole mesh* dim-agnostic and
admits the first 3-axis Cartesian SNMesh. The decisive design fact:
**axis-native, NO Mesh3D** (a Mesh3D would have exactly one consumer →
"unify after two instances" forbids it). The axes tuple, made PRIMARY by
C5.1, IS the N-D entry. The campaign is sequenced *clean before extend*:
C5.1–C5.4 invert+de-phantom the mesh layer, C5.5 admits d=3 as a one-line
gate removal.

## Where it landed (home selection)

The C5 section goes in `boundary_conditions.rst`, immediately AFTER the
C4 `bc-face-name-carve` arc (the C4 section already anticipated C5: "when
`Mesh3D` lands in C5" — which C5 *contradicts* by going axis-native, so
that forward-ref had to be corrected, not just appended-to). The SNMesh
*construction narrative* ALSO lives in `discrete_ordinates.rst`'s
"Two-Layer Mesh Pattern" — its "`SNMesh` **wraps** the base mesh" claim
is exactly what C5.1 inverts (post-C5 `self.mesh` is inbound-provenance,
None at d≥3), so that section needs a focused axis-primary correction +
a cross-ref to the new C5 anchor. Two homes, one canonical section
(`boundary_conditions`) + one construction-narrative touch
(`discrete_ordinates`).

## The decisive doc moves (reusable for a constructor-inversion + N-D-admission carve)

1. **Lossy-round-trip as an ASCII data-flow diagram.** The pre-C5
   pathology is `axes → mesh → axes` (from_axes synthesizes a legacy
   mesh; __init__ discards the tuple and re-derives). Two `.. code-block::
   text` diagrams — the lossy BEFORE (`from_axes → legacy_mesh_from_axes
   → Mesh → axes=axes_from(mesh) ← re-derived`) and the one-body AFTER
   (both surfaces funnel into `_init_core`, conversion is one-directional
   `mesh → axes` only on the legacy surface). The diagram makes "lossy in
   two ways" (custom labels reset + d=3 has nowhere to round-trip
   through) concrete. This is the load-bearing WHY — d=3 *appeared* to
   need a "third arm" ONLY because the inverted flow mandated a legacy
   mesh at every dimension.

2. **The bitwise-identity equation gets an eq-label + the mesh.py line
   citations.** `axis_widths[a] = np.diff(axes[a].edges)` is the single
   identity that dissolves the isinstance metadata branch; tag it
   `:label: sn-axis-widths` and cite that it's byte-identical to
   `Mesh1D.widths` / `Mesh2D.dx` / `.dy` (mesh.py:287/:567/:572 — same
   np.diff upstream). This is the "structural carve, no numbers change"
   proof at the metadata level (the sha256 affine goldens are the proof
   at the output level).

3. **The Mode-9 gate retarget is a 2-row before/after `.. list-table::`.**
   C5.4 is the highest-risk edit (vv-principles Mode 9 — a splitting/opt
   verified only where the wrong gate is accidentally satisfied). The
   pathology: `reduced is None` is a COINCIDENCE proxy true for EVERY
   multi-D Cartesian incl. d=3, but the in-sweep moment kernel is 2-D
   only → at d=3 it would have SILENTLY moment-windowed (corrupted
   iterate, not a principled refusal). Table cols = Gate / pre-C5.4
   coincidence-proxy / post-C5.4 genuine-condition (`is_cartesian and
   ndim==2` for windowing; `is_cartesian and not is_1d` for G-S).
   Explicitly call out the ONE branch that STAYS `reduced is not None`
   (the 1-D sweep-cache — there it keys on data AVAILABILITY, the genuine
   guard, not a dimension proxy). Frame the G-S "2-D ONLY" docstring as
   STALE Phase-3 NARRATION (the schedule+scheduled-sweep are d-generic
   since C3) — corrected, with the d=3 FP-invariance VALUE-GATED by the
   C5.5 Mode-9 box (never gate a splitting's FP-invariance on a
   degenerate box — vv Mode-9 made operative).

4. **The d=3 value-gate table maps 1:1 to V&V pillars (vv-principles
   discipline).** 4 gates: k_inf 3-D≡2-D≡1-D (L1 closed-form eigenvalue
   — `λ_max(A⁻¹F)`, A=diag(Σt)−Σs0ᵀ, NEVER the sweep; 2g+4g NEVER 1g —
   1g is the degenerate νΣf/Σa ratio; atol 1e-8, d=3 box solved
   1.8750000050 vs 1.875); per-ordinate ψ=Q/(WΣt) (L1 closed-form flux,
   rtol 1e-10, sharpest Mode-1/3/4 probe, DD flat-flux exact); scattering
   φ=(diag(Σt)−Σs0ᵀ)⁻¹Q (L1, Mode-6 group-coupling catcher because
   mixture C's SigS is ASYMMETRIC so SigS vs SigSᵀ observable; measured
   2.6e-9 = SI-convergence-limited NOT discretization); Mode-9 G-S≡Jacobi
   FP-invariance (L2, on the degenerate-BREAKING box: mixed BCs +
   nx≠ny≠nz + het + DIAGONAL cubature/ERR-056 shared-face). Explicitly
   note the eigenvalue reference is STRUCTURALLY INDEPENDENT of the sweep
   (matrix eigensolve, not transport solve) — satisfies the
   vv-principles "eigenvalue needs closed-form, not MMS" boundary.

5. **The two-default-BC-conventions note.** A `.. note::` capturing the
   genuine subtlety: the SOLVER entry defaults un-declared faces to
   VACUUM (fixed-source: un-specified boundary leaks); a fresh SNMesh
   with no BC declarations defaults to REFLECTIVE (infinite-lattice
   eigenvalue convention). Both preserved per-surface at d=3; the gates
   exercise reflective (eigenvalue k_inf) + mixed (Mode-9).

6. **"What runs the d=3 path, what's deferred" = the honest-scope
   subsection.** d=3 runs on the d-generic FullFieldWavefront ORACLE
   spine (correct from day one) NOT the optimized kernels; ScanMarch
   d≥3 + MovingFrontierWindow d≥3 deferred to #227, MEASUREMENT-gated
   (C3.6 "construct general, select narrow, specialize on measured
   cost"); multi-D adjoint = pre-existing orthogonal deferral (G-adjoint
   territory). Naming the deferrals + their issue keeps the capability
   map complete and pre-empts "did you forget the d=3 kernels?".

## Stale-ref sweep — LIVE-vs-HISTORICAL method-rename triage

C5.3 renamed `InflowTraceSpace.from_mesh_and_quadrature` →
`TraceSpace.from_quadrature_and_layout` (geometry-blind, layout-driven,
mesh param retired). The ~6 elegance-flagged refs split into:
- **LIVE-API (WRONG → fix to the new method):** the "construction goes
  through" claims in boundary_conditions + operator_algebra. Rewrite to
  `from_quadrature_and_layout` + drop the 2-D-cylindrical-NotImplementedError
  "live" claim (the refusal is UNREACHABLE — such a Mesh2D can't become
  an SNMesh; it lives at the construction surface, not the trace factory).
- **HISTORICAL (#188-era):** the "Pre-cleanup history" / "Why the split
  existed" blocks. KEEP as history, render the old method name as a
  parenthetical literal + name the live successor: "(then named
  `InflowTraceSpace.from_mesh_and_quadrature`, since C5.3 the
  geometry-blind `TraceSpace.from_quadrature_and_layout`)".
- **PERVASIVE class-name drift (`InflowTraceSpace`/`OutflowTraceSpace`,
  ~12 sites) = NOT C5's scope.** That class NEVER existed in code (live
  class is always `TraceSpace`); the names describe the pre-#188
  Γ₋/Γ₊ split-trace CONCEPT. This is squarely #223 ("still narrates the
  pre-#188 split trace"). DEFER — touching only the C5-changed sites.
- **The worked-example CODE BLOCK (`InflowTraceSpace.from_mesh_and_quadrature(mesh,
  quad, faces=...)`, ~886-909) is #223's, NOT trivial to half-fix.** The
  signature itself changed (mesh/faces args gone), so renaming JUST the
  method leaves a self-inconsistent block — STRICTLY WORSE. The right
  C5-scoped move: add a `.. note:: Stale narration — pre-#188 split trace
  (tracked by #223)` ABOVE the block naming the live unified API, so a
  reader knows it's historical illustration not current API — WITHOUT
  doing #223's holistic rewrite. (Cardinal Rule 3: never leave a reader
  thinking stale code is live.)

## index_convention.rst — the rank-generic note (NOT a rewrite)

The `(N, ng, nx, ny)` page (Issue #196 storage layout) is NOT directly
invalidated — d≤2 layout is byte-identical. But it's now a d≤2
SPECIALIZATION of the rank-generic `(N, ng, *spatial_shape)`. ONE
`.. note::` in Key Facts: at d=3 → `(N, ng, nx, ny, nz)`; energy-first /
ordinate-first hold at every rank, only the spatial-tail LENGTH changes;
every field/xs/scattering read since C5.2 keys on rank-generic
`spatial_shape` NOT a hard-coded `(nx, ny)` pair (an `(nx,ny)`-keyed read
SILENTLY TRUNCATES a 3-D tensor — the live d=3 landmine C5.2 retired).
Cross-ref the C5 anchors. Minimal touch — do NOT rewrite the whole
2-D-shaped page.

## Cross-ref reality on this carve (matches C4 memo)

- `orpheus.numerics.spaces.trace_space` is NOT automodule'd → my
  `:meth:`TraceSpace.from_quadrature_and_layout`` / `:class:`TraceSpace``
  render PLAIN TEXT (no -W warning). BUT the original `InflowTraceSpace`
  refs were ALSO plain text (that class never existed in code) — so the
  fix keeps the same render with the CORRECT name. Same for `sn.axis`,
  `face_layout`, `material_xs_field`, `scattering`, `method_space`,
  `pole_angular_closure`, `sweep_schedule` (all NOT automodule'd → plain
  text by existing-page convention; surfacing them is a separate
  architectural-docs task — DEFER).
- RESOLVE: `orpheus.sn.geometry` (SNMesh.*), `orpheus.sn.solver`,
  `orpheus.sn.loss_representation` (FullFieldWavefront's owning module),
  `orpheus.numerics.measure.DiscreteMeasure` ARE automodule'd. A bare
  `:class:`DiscreteMeasure`` (no module) does NOT resolve — use the FQN.
- All 8 new C5 `:ref:` anchors (`sn-axis-primary-c5`, `sn-c5-*`) defined
  in the new section; the cross-doc `:ref:` from index_convention →
  boundary_conditions RESOLVES (verified live `href` in the rendered
  HTML). grep the C5 anchor set against the cited set BEFORE building.

## Build gate (this worktree — DIFFERENT from C4 memo's 11)

Task brief said "zero-warning baseline". RECONCILED: the INCREMENTAL
build (env populated — `_generated/*.rst` materialized by a prior `-E`)
is genuinely ZERO Sphinx WARNING/ERROR/CRITICAL (the only noise = 6
python `SyntaxWarning`s from autodoc importing 2 test files w/ bad escape
seqs — NOT doc warnings, untouched files). The forced `-E` COLD build is
"build succeeded, 11 warnings" — all env-artifact + pre-existing
(homogeneous plotting ×2 / verification.rst include CRITICAL ×8 / mesh.py
paramref ERROR ×1; the `_generated`/`.h5` lines are the InputError/
FileNotFound DETAIL of those same entries — raw grep double-counts to 21,
summary says 11). GATE: incremental clean + `-E` warning-SET byte-identical
to the pre-edit `-E` baseline (proved by `diff` of line-normalized sort
-u sets = empty = PASS) + grep `boundary_conditions|discrete_ordinates|
operator_algebra|index_convention|sn-c5|Title underline|Inconsistent` in
the log EMPTY. EXIT=0 authoritative.

## Verification matrix (auto-regenerates)

The 4 `test_d3_admission.py` `verifies()` tags (matrix-eigenvalue,
multigroup, transport-cartesian) fold in AUTOMATICALLY on build —
`solve/test_d3_admission` appears in matrix.rst (`0,5,0,0,2,0`) and the
3 labels' counts include its edges. Confirmed in rendered HTML. No manual
matrix edit needed (registry is single source of truth — see
[[feedback_autogen_tables]]).

## Quality scores (this task)

| Dimension | Score | Note |
|---|---|---|
| Derivation depth | 5 | lossy-round-trip data-flow (before/after) + axis-widths eq-label w/ mesh.py byte-identity citations + Mode-9 gate-retarget table + 4-gate value-evidence table w/ exact numbers |
| Cross-references | 4 | every SNMesh/solver/FullFieldWavefront/DiscreteMeasure symbol linked; trace_space/sn.axis/face_layout/method_space refs plain-text (unsurfaced modules, deferred, MATCHES the original behavior) |
| Numerical evidence | 5 | k_inf 1.8750000050 vs 1.875 atol 1e-8, per-ordinate rtol 1e-10, scattering 2.6e-9 measured, G-S≡Jacobi atol 1e-8; each mapped to its V&V level + failure mode + structural-independence note |
| Failed approaches | 5 | the lossy round-trip + WHY it mandated a legacy mesh at every d (the d=3 "third arm" illusion); the Mode-9 coincidence-proxy that would have silently corrupted a d=3 iterate; the unreachable-2D-cyl-refusal |
| Code traceability | 5 | every claim → :meth:/:attr:/:func: + the 2 dr-consumer + field/xs/scattering migrations named; code-fact grep-verified (from_axes/axes/coord/volume_measure/coord_system/_as_sn_mesh/mat_map/_maybe_window/_select_si_resolvent all exist as cited) |
| Derivation source | n/a | structural/architecture carve, no SymPy derivation needed (value gates are the evidence, not a manufactured solution) |
