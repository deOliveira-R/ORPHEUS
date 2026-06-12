---
name: c5-4-5-6-d3-admission
description: C5.4/5.5/5.6 (#225) 3-D Cartesian admission review — PASS-with-nits; the headline FINDING = the C5.3-carry fail-loud guard is DEAD (mu_z is always-present @property zero-fallback, `is None` never fires)
metadata:
  type: project
---

C5.4 `580bc9b` / C5.5 `1da1e2f` (headline) / C5.6 `294f9c0` (#225,
branch worktree-sn-nd-layout) — the 3-D Cartesian admission. Reviewed
PASS-with-nits. Siblings [[c5-1-axis-primary-snmesh]] [[c5-2-phantom-shim-retirement]]
[[c5-3-geometry-blind-trace]]. d≥3 now constructs+SOLVES through the
d-generic FullFieldWavefront spine, mesh-LESS from birth (`mesh is None`),
NO Mesh3D dataclass. Entry = the axes tuple via ONE `_as_sn_mesh` seam.

**★ THE HEADLINE FINDING (CONCERN, latent-not-active) — the C5.3 carry is
NOT closed; the fail-loud guard is DEAD CODE.** C5.5 added to
`trace_space._build_omega_dot_n` (trace_space.py:236):
`if getattr(quadrature, f"mu_{AXIS_NAMES[axis]}", None) is None: raise`.
But `mu_x`/`mu_y`/`mu_z` are CLASS-LEVEL @property on `Quadrature`
(directional.py:293-317) delegating to `axis_cosines(i)`, which RETURNS
`np.zeros(n_points)` for `axis_index >= nodes.shape[1]` — NEVER `None`.
So `getattr(q,"mu_z",None) is None` is ALWAYS False. VERIFIED by driving
a genuine 2-cosine `DiscreteMeasure` quad through `_build_omega_dot_n`
with a `("...","zmin","zmax")` faces tuple: NO RAISE, zmin row all-zero
(all-tangential) — EXACTLY the rank-mismatch the comment names ("a z face
on a 2-D quadrature"). The guard cannot catch the case it documents.
PROOF-OF-PATTERN: `test_angular_average_operator.py:167-174` ALREADY
documents that a sibling "requires mu_z" early-return does NOT fire on
Gauss-Legendre "because 1-D GL carries mu_z=zeros(N)" — the team already
hit this; a downstream "no outgoing ordinates" guard catches it there.
NO test pins the new guard (grep: zero refs to `_build_omega_dot_n` /
"rank-mismatch" in tests) → invisible, no false-green but no protection.
Not-active TODAY: production d=3 uses `level_symmetric` = genuine
3-cosine, `mu_z` non-zero, table correct (VERIFIED zmin==−mu_z). FIX
destination: discriminate "all-zero cosines on a face the layout names as
NORMAL" (`np.allclose(_quadrature_axis(q,axis), 0.0)`), NOT `is None`;
plus a test driving the 2-cosine-quad-with-z-face raise.

**C5.4 `580bc9b` — PASS clean.** Both SI gates retargeted off the
`reduced is None` coincidence-proxy: `_maybe_window`→`is_cartesian and
ndim==2` (single site, Pattern 7); `_select_si_resolvent` G-S→
`is_cartesian and not is_1d` (multi-D; `_sweep_scheduled` d-generic since
C3 — the "2-D ONLY" was stale Phase-3 narration, correctly corrected, NOT
a behavior change). 4 stale proxy-narration docstrings swept. The 1-D
sweep-cache `reduced is not None` KEPT is CORRECT discrimination — it keys
on AVAILABILITY of the reduced data it READS, not a dim proxy (data-keyed
≠ dim-keyed; honest). Mode-9 risk (silent d=3 moment-window of the SI
iterate) correctly identified as the highest-risk C5 edit and value-gated.

**C5.5 `1da1e2f` — PASS-with-nits.**
- volumes `reduce(np.multiply.outer, axis_widths)` inline = FINE (named
  LHS `_volumes`, dimensional comment `V[i,j,k]=Δx_i·Δy_j·Δz_k`; Pattern-3
  satisfied). `volume_measure` meshgrid arm REUSES `self.volumes.ravel()`
  for weights (single source for magnitudes; only nodes recomputed) → NO
  twin volume formula. `from functools import reduce` function-local =
  cosmetic nit only.
- `_apply_default_bcs` axis branch: reader uses polymorphic `ax.bc`
  (right); WRITER reaches behind via `isinstance(ax,RadialAxisMesh)` →
  `bc_outer` vs `bc_low/high`. CONCERN-not-VIOLATION: SINGLE production
  writer site (solver.py:93; adapters branch on AxisCoord not type, and
  CONSTRUCT fresh not fill-default), only 2 axis types → unify-after-two
  defers. Destination when a 3rd axis type lands: `Axis1D.with_default_bc(bc)`
  polymorphic writer on the Protocol (it already has polymorphic `endpoints`
  + `bc` reader, lacks the writer).
- `mat_map` kwarg + legacy-mesh RAISE (solver.py:125): runtime check of an
  invariant the type system could carry (illegal `mesh=legacy & mat_map=X`
  combo IS representable; anti-pattern 8). CONCERN-acceptable = LEAST-BAD:
  no `Mesh3D` type to hang `mat_map` on, alternative blocks heterogeneous
  3-D entirely; fail-loud, dissolves when a `Geometry` sum-type per-variant
  material channel lands. Crosswalk reasoning is sound (legacy carries
  `mat_ids` field; axis tuple has no material channel).
- Tests EXEMPLARY: k_inf 3-D≡2-D≡1-D headline (matrix-eigenvalue ref never
  touches sweep = structurally independent; ≥2g only, 1g excluded as
  degenerate); pure-absorber per-ordinate ψ=Q/(WΣ_t) rtol 1e-10 (DD
  flat-flux exact, c=0 tolerance-tight story HONEST) + scattering companion
  φ=(diag Σ_t−Σ_s0ᵀ)⁻¹Q (Mode-6 asymmetric-S catcher, convergence-limited
  rtol 1e-7 justified); Mode-9 G-S≡Jacobi on the degenerate-BREAKING box
  (mixed x-refl/y-vac/z-refl, nx≠ny≠nz=5,3,4, het 2G split, DIAGONAL
  level-symmetric cubature ERR-056, non-flat guard max/min>1.2). np.testing
  only (Mode-8). C5.1 refusal pin DELETED per its own removal trigger
  (retirement=test-migration honored).

**C5.6 `294f9c0` — PASS-with-CONCERN (stale-doc cluster).** G-c3
`TestD3SupportsMatrix` FLIPPED to LIVE `default_for(real_from_axes_d3)`
landing on `FullFieldWavefront`; uses `type(selected) is not
FullFieldWavefront` EXACT-type (CORRECT — `FullFieldWavefront` and
`MovingFrontierWindow` are sibling `_DAGWavefront` subclasses; `isinstance`
would falsely pass for the window). order-INVARIANT/order-RELATIVE split is
PRINCIPLED: invariant facts (8 unmerged groups, reflected-set==reflective-
faces, each shared face assigned EXACTLY once ERR-056, reflect⊆outgoing)
go LIVE on the real mesh; order-RELATIVE ("which group is last") stays
synthetic because it depends on the cubature's octant enumeration — a real
discriminator, right call.

**THE C5.6 CONCERN — stale "d=3 unconstructible / when Mesh3D lands" doc
cluster the admission newly invalidated, NOT swept (de-roling discipline
violated):**
1. `test_sweep_schedule_nd.py:13-17` MODULE HEADER still asserts "a d=3
   SNMesh is unconstructible today (no Mesh3D; legacy_mesh_from_axes
   raises) ... never a d=3 flux VALUE — the value claim lands with Mesh3D"
   — DIRECTLY refuted by the SAME COMMIT's new `test_gs_d3_schedule_from_real_mesh`
   which constructs a real d=3 mesh in that very file. Sharpest instance.
2. `loss_representation.py:1475-1476` `default_for` docstring Raises:
   "a d≥3 Cartesian mesh — WHEN Mesh3D LANDS — → FullFieldWavefront" — the
   function C5.6's headline test now exercises live still says its d=3 path
   is contingent on Mesh3D.
3. `loss_representation.py:1194-1195` `ScanMarch` docstring "when one
   becomes constructible (Mesh3D, C5)". Routing CLAIM still correct
   (d≥3 falls through ScanMarch.supports → FullFieldWavefront); the
   contingency narration stale + miscredits Mesh3D for the mesh-less
   mechanism.
4. `axis.py:683-688` `legacy_mesh_from_axes` d≥3 NotImplementedError —
   BEHAVIOR still correct (legacy adapter genuinely d≤2; C5.5 routes
   around with `mesh=None`), MESSAGE rationale stale ("the 3-D admission
   gate exercises pure shape functions WITHOUT constructing an SNMesh" —
   C5.5 DOES construct one).
DISCRIMINATOR for this cluster (consistent with [[c5-3-geometry-blind-trace]]
+ [[c3-6-honest-d-dispatch]] doc rules): a claim "d=3 needs Mesh3D / is
unconstructible" = STALE→fix (mechanism is mesh-less `from_axes` now);
a routing claim "d≥3 falls through to FullFieldWavefront" = STILL CORRECT,
leave. Behavior-correct-message-stale raises (axis.py:685) = message-fix.

Gates VERIFIED GREEN: d3_admission+si_gate+axis_native 26✓; schedule_nd+
unified_dispatch 21✓. The dead guard does not false-green (no test pins it).
