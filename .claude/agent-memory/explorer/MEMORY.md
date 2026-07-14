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
  boolean-presence ≠ integer-width before a typed probe-swap (L6).

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

- [task #54 sn/spatial rename audit](task54_sn_spatial_rename_audit.md) — TRANSIENT @
  `6e3ebad0` (2026-07-13): full 3-search blast radius for renaming `orpheus/sn/spatial/`
  (4 production import files / 7 statements, all direct-module-path — the package
  `__init__` has ZERO live importers; 20 test import statements / 15 files; 147 doc
  role occurrences, 10 already dangling; zero string-literal/config refs; mutation
  TOMLs + diag_276 pre-broken since #272; name-candidate evidence: walk lives in
  loss_representation, tests/sn/sweep is broader than the package). **Delete when #54 lands.**

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

- [HarmonicMomentField UNITS convention](harmonic_moment_field_units_convention.md)
  — why a stored SH moment carries SCALAR-flux units (no-prefactor SH, Y₀⁰=1,
  weights sum to 4π → sr cancels); R≠M*; the ERR-039/ERR-051 history behind it.
- [Phase 5 µ-resolved primitive inventory](phase5_mu_resolved_primitive_inventory.md)
  — µ-resolved vs µ-integrated primitives in `peierls_geometry.py` for the
  continuous-µ specular multibounce closure.
- [Pyright ignored-package measurement](reference_pyright_ignored_package_measurement.md)
  — true error count for a [tool.pyright]-ignored package: CLI file args +
  /tmp ignore-free config; discount the editable-install import artifacts.
