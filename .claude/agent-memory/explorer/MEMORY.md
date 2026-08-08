# Explorer memory index

Slim index. Behavioral lessons live in `lessons.md` (read FIRST each dispatch).
This index holds only (2) git-true active-state and (3) durable reference. The
per-campaign `file:line` carve/audit maps are archaeology — they go stale in days
and are re-derivable in seconds via Nexus. Keep notes here DURABLE, not transient.

## 1. Lessons (read first)

- [lessons.md](lessons.md) — exploration lessons. The spine (blast-radius =
  graph+grep+constructors+doc-nodes; verify-premise-first; durable-shape vs
  line-map; git-is-authoritative-for-merge-status) is now PROMOTED to AGENT.md
  Operating Principles 4–7 — L1/L2/L3/L5 remain as forensic war-stories with
  "→ now in AGENT.md" pointers. Lesson-only (narrower question shapes): carve
  verdicts name both retire AND keep-as-anchor with the discriminator (L4);
  boolean-presence ≠ integer-width before a typed probe-swap (L6); **the BRIEF's own
  timeline is a claim (verify commit dates vs the target doc's mtime), and a prior
  audit's "X cannot be expressed" EXPIRES when a SIBLING campaign lands a substrate —
  read the new module, don't just re-run the emptiness greps (L20)**; **"what breaks
  if this numeric primitive changes?" → SWAP IT AND RUN; the dangerous assertion
  class is a FROZEN rhs, not an exact comparison (L13)**; **a "fold this redundant
  DOF" hypothesis is refuted (or not) by the FUNCTIONALS' parity, not the
  algorithm — and fold-the-algorithm ≠ fold-the-state (L15)**; **a stored numeric
  TAG (exactness/order/rank) is a claim — sweep it, and first measure what the
  declared SYMMETRY gives for free; a tag-pinning `assert x == <tag>` is evidence
  the property is UNTESTED (L16)**; **a static table behind a computed fast path has
  DEAD ROWS — measure which are consulted before costing a change; and a tag routing
  through two dispatch branches is invisible on the fixture both accept, so find the
  input that FAILS (L18)**; **a brief's TYPE table is a claim about MATERIALIZED objects
  (grep the accessor's return type before mapping construction sites), and an all-green
  "what breaks?" run may have measured INERT rather than SAFE — find the gate and check
  the path routes through it (L21)**; **a "what REMAINS of issue X" recon greps the
  campaign's own MID-FLIGHT PROSE ("today still…", "remaining half…", cautions) — the
  interim honesty notes are the falsified class; and an acceptance-gate check reads the
  FIXTURE LIST, not the assertion (L22)**.

## 2. Active / in-flight state

> **2026-07-13 hygiene pass:** the #280 walk-unification → coupled-block operator
> campaign is CLOSED AND MERGED (step-7 "campaign CLOSES" @ `b23e972e`; #34 ray-leg
> retirement @ `015dcc73`/`03e275e8`; all git-verified ancestors of main `6e3ebad0`).
> The six transient maps it carried (step-6 audit, 4e un-weave, B.2 blast radius,
> 2.5 P0-A, 2.5d carrier, product-cyl seed) hit their delete-conditions and were
> DELETED. Post-campaign durable facts: A_BB = `RadialCharacteristicOperator`
> (`orpheus/sn/operators/radial_characteristic.py`) now WRAPS the ψ½ march — its
> `.solve` is the sole production caller of `carlson_inward_sweep_from_source`;
> the walk executors (`_OneDimScanWalk`/`_loop_walk`/`_dag_legs`) live in
> `sn/loss_representation/__init__.py`.

> **2026-07-21:** task #54 (`sn/spatial` → `sn/sweep`) is MERGED (`588f2429`
> git-verified ancestor); its transient audit memory hit the delete-condition and
> was DELETED.

The durable SHAPE of the SN operator-algebra subsystem lives in
`.claude/agents/explorer/AGENT.md` ("SN operator-algebra subsystem — durable
shape") — read that, not a frozen file-list here. Every SN campaign this agent
audited — the operator-algebra unification, the Wave-O / role-typing / g-adjoint
work, the matvec carve onto `_OneDimScanWalk`, LD-on-the-DAG, the foundation
cleanup (moment-resolved source + trace widening + predicate scoping), and the
field-typed algebra map — is **merged to main** (git-verified 2026-06-22; #236
merged too, per 2026-07-03 recheck — the only surviving local SN branch is
`feature/sn-adjoint-transport`, the paused #276 campaign).

> Merge-status in memory goes STALE (lessons L5). ALWAYS reconcile any "resume X"
> against `git merge-base --is-ancestor <hash> HEAD` before acting; never trust a
> frozen "NOT pushed". Landed milestones live in the SN theory page's development
> history (`docs/theory/discrete_ordinates.rst`), not here.

## 3. Durable reference (survives code churn)

These are convention/units facts that a line-number drift cannot invalidate —
they pin WHY a quantity carries the units it does. Keep them; everything else
in the old index was a frozen carve-map and is proposed for retirement.

- [Spatial transform category — durable](spatial_transform_category_durable.md) —
  the spatial layer's mirrors are **E(d) = O(d) ⋉ ℝ^d**, not O(3): `_orbit_closure`
  on a CENTRED cell lattice already returns the sweep DAG's `arange(n)[::-1]`, but
  production meshes start at `origin=0` so `Mirror.is_invariant` is False; the gap
  is the TRANSLATION (no affine type exists). Two genuine group objects:
  `octant_moment_frame_signs` = the character rep of (Z₂)^d; `reflection_index("x")`
  = the r=0 quotient's deck transformation. Chart map never written (only its
  integrated Jacobian). Sweep reversal spelled 11×; done right once (the adjoint).
- [Angular layer — hidden geometric transformations](angular_layer_hidden_transformations.md)
  — where the angular/quadrature/moment layer rotates/reflects WITHOUT naming a group
  element: the SH basis's polar axis is `μ_x` ⟹ every `Y_ℓ^m` carries an unnamed 120°
  `O_h` rotation about `(1,1,1)` that `SubgroupOfO3` cannot tag; "octants" is really the
  `(Z₂)³=D_2h` orbit stratification (26 parts on Lebedev-17, not 8); `_orbit_closure` is
  one of SIX partner-map engines; the curvilinear fiber circle is only ever an ordered
  interval. Two measured defect leads (slab `ℓ≥2` SH docstring is false; the `arctan2`
  round trip destroys the roots-of-unity exactness) + the six open #325 sites.
- [Quadrature landscape — durable shape](quadrature_landscape_durable.md) — which
  of {range, spacing, rule-on-circle, rule-on-interval, exactness-space,
  node-generation} has ≥2 realizations; MoC's unnamed `[0,π)` quotient + Σω=1 vs
  SN's `[0,2π)` + Σw=4π; `level_symmetric_sn` is EQUAL-weight and measures degree
  **3 for every N** (tag says N−1) ⟹ discrete SH Gram 25–45% off at L≥2 while the
  frame assumes the analytic Gram; `min(2n_mu−1, n_phi−1)` is CORRECT and sharp.
  **+ the generation-KIND census (2026-08-02):** GROUP-ACTION is exactly ONE
  (`roots_of_unity`, still 0 consumers/unexported); `level_symmetric` is HYBRID not
  exact (sign subgroup exact, orbit rep by formula); the CHECKER's own ops split the
  same way (O_h exact, `C_n`/`D_nh` by cos/sin); `is_invariant`'s only production
  caller is unreachable ⟹ every shipped `invariance_group` is a DECLARED tag; and
  the partner maps are 3 SEARCH (no distance guard) vs 1 FORMULA.
- [SN multigroup axis structure](sn_multigroup_axis_structure.md) — three-tier
  group-blindness (structure / group-diagonal data / S+F coupling); NO group loop
  ("within-group" = fission-external, full multigroup S inside); τ/c are ANGULAR
  closure weights not optical thicknesses; AGENT.md drift: `_within_group_triple`
  → `build_within_group_system` (coupled_system.py).
- [HarmonicMomentField UNITS convention](harmonic_moment_field_units_convention.md)
  — why a stored SH moment carries SCALAR-flux units (no-prefactor SH, Y₀⁰=1,
  weights sum to 4π → sr cancels); R≠M*; the ERR-039/ERR-051 history behind it.
- [Harmonic frame + folded-quadrature plumbing](harmonic_frame_folded_quadrature_plumbing.md)
  — Q5.6: NO computed Gram on the kernel path (raw Y^TW analysis; R carries the
  CONTINUUM (2l+1); exactness ASSUMED ⟹ the folded ξ-odd garbage); 3 shape-contract
  tiers (table-derived adapts, L-derived + fixed-slot hardcode rectangular ⟹
  zeroed-columns sub-basis, not flat); σ_y parity per slot (cos branch odd-m ∪ sin
  branch even-|m|); cylinder-P1 gate is route-equivalent = folded-blind.
- [Cylindrical SN level-order sensitivity (#326)](cylindrical_sn_level_order_sensitivity.md)
  — α IS `−W·ξ` at a half-angle boundary (Hébert 3.399, exact via Dirichlet kernel); the
  recursion is a cumulative integral in ω ⟹ the level must be a HALF range (Hébert
  `0≤ω≤π`; BMC Eq. 52 `Σw̄=2sinθ`) but ORPHEUS spans `[0,2π)`; every existing α gate is
  telescoping-blind; the ξ-mirror invariant (not the MMS) is the adjudicator.
- [Non-SN geometric-transform census](nonsn_geometric_transform_census.md) — zero
  hand-built rotation/reflection matrices outside `numerics/symmetry.py`; trig is 5
  lines; MoC owns the only 2 hand-rolled `_orbit_closure` clones (both guard-free);
  MC has no reflection (periodic only); CP's images are an ADDITION; 4 spellings of
  one `(I−TR)⁻¹` orbit sum. Two durable identities: MoC's azimuthal set IS the upper
  half of `periodic_trapezoid(2n, STAGGERED)`, and its mirror partner is `k↦n−1−k`.
- [Phase 5 µ-resolved primitive inventory](phase5_mu_resolved_primitive_inventory.md)
  — µ-resolved vs µ-integrated primitives in `peierls_geometry.py` for the
  continuous-µ specular multibounce closure.
- [Pyright ignored-package measurement](reference_pyright_ignored_package_measurement.md)
  — true error count for a [tool.pyright]-ignored package: CLI file args +
  /tmp ignore-free config; discount the editable-install import artifacts.
- [Cylindrical SN level-order sensitivity](cylindrical_sn_level_order_sensitivity.md)
  — #326: the 1-D cyl sweep reads ONLY (η, w); ξ enters solely via the source; a
  mirror-pair tie splits τ into {1, 0.5} by 1 ULP of cos φ. Measured: α bit-id,
  τ/c permute, ψ a pure permutation, φ 1-2 ULP isotropic but ~20% ξ-dependent.
  **+ the HALF-RANGE fix map**: the fold belongs in `LevelStructure` (the quotient
  S²/Z₂), NEVER in the `DiscreteMeasure` (a full-S² cubature 2-D Cartesian
  consumes); (A) fold-existing-nodes vs (B) Hébert-midpoint are separated by the
  **R12a predicate**, not by α (both give the same `2κ`) — (B) flips every cyl
  level to ψ½-CARRYING; the ONE real break is the **ξ-odd SH moment**, and it
  vanishes if you fold the ALGORITHM (+lift) rather than the STATE; the
  `circle-vs-interval` theory section's periodicity principle is FALSIFIED.
