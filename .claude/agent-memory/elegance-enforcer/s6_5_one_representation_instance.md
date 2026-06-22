---
name: s6-5-one-representation-instance
description: S6.5 (#222 capstone) — the "one representation instance" carve closed the operator's FOUR construction doors to ONE; PASS clean. The twin-DELIVERY plumbing pattern finally collapsed at the INSTANCE level (not just operator level).
metadata:
  type: project
---

S6.5 of #222 — "one representation instance" unification. Reviewed on
`refactor/field-role-typing`-lineage worktree `sn-nd-layout`, working-tree diff. **PASS clean
(no conditions).** This is the carve that finally collapses the standing
twin-DELIVERY plumbing smell (AGENT.md Institutional Pattern #1) at the INSTANCE
level — the prior carves single-sourced the OPERATOR; S6.5 single-sources the
representation INSTANCE the operator's two actions consume.

**Why:** `LossRepresentation` is a STATELESS frozen dataclass (mesh-only field),
so "matvec ≡ sweep — two actions of ONE operator" (L21) was only a *coincidence
of FOUR `default_for` calls agreeing*, not a type fact. S6.5 makes it a type fact
by construction.

**How to apply (verification recipe for any "collapse N doors to 1" carve):**
1. Grep the construction-site inventory: `MovingFrontierWindow(`, `CumprodScan(`,
   `ScanMarch(`, `FullFieldWavefront(`, `default_for(` across `orpheus/`. After
   the carve the ONLY live production call sites must be (a) the operator's
   cached_property door and (b) any DELIBERATELY-retained operator-free
   functional entry. S6.5 result: exactly 2 (`operator.py:1549` cached_property,
   `loss_representation.py:1613` inside `transport_sweep`). Everything else is a
   comment or class/def.
2. The bit-identity-by-construction claim "dropped the `AngularSourceSink` wrap"
   was VERIFIED at the source: `from_mesh` stores `values=values` by reference
   (`_bases.py:212`, no copy/asarray) and `_unwrap_source` returns `.values`, so
   the unwrapped object IS `rhs.bulk.values` — passing it directly is the same
   object. `Q` is read-only in sweep bodies. Always trace the wrap/unwrap pair to
   the storage primitive before accepting "bit-identical by construction".
3. Scope-boundary honesty: the retained `transport_sweep` has exactly ONE live
   caller (`solver.py:1634`, the `solve_sn` post-convergence reconstruction). It
   genuinely has no `InvertibleOperator` in scope (`SNSolver` is the outer
   eigenvalue driver, NOT the `(L+C)` operator) AND it needs the typed
   `from_isotropic` factory (`/W` projection) — so the typed boundary is
   load-bearing for THIS caller, not dead ceremony. Retention is principled.

**Delegation SSOT (axis 2):** `InvertibleOperator.loss_representation` is a PLAIN
`property` forwarding to `self.streaming.loss_representation` (the leaf's
`cached_property`). Correct: caching the forward too would mint a SECOND handle =
the very coincidence being eliminated. Mesh-identity invariant enforced at
`InvertibleOperator.__init__` (streaming.sn_mesh is diagonal.sn_mesh). The
plain/cached asymmetry is the right call, NOT a staleness risk.

**Test discrimination (axis 4) — the S6.4 one-walk SPY pattern reused:** spies
capture the call-time `self` of `MovingFrontierWindow.loss_action`/`sweep`/
`_sweep_interior` and assert object identity with `A.loss_representation`. Proves
the SOLVE PATH *runs on* the object (behavioral), not that two property handles
compare equal (spelling) — the docstring's non-tautology argument is exact.
`pytest.fail` only (no bare asserts) → `-O`-safe. The 3 re-pointed spies in
`test_invertible_operator.py` moved from `patch("...transport_sweep")` to
`patch.object(CumprodScan, "sweep")` — the FALSE-GREEN-AFTER-RELOCATION class:
solve no longer routes through the free function, so the old patch would pass
vacuously on zero calls. Class-level patch == call-time instance-method
resolution site (1-D meshes → `CumprodScan`). Correct spy-target migration.

**Test-helper promotion:** `cart2d_2g_nonsquare`/`het_operands` moved from
`test_one_octant_walk.py` to `tests/sn/_test_helpers.py` VERBATIM (byte-identical
bodies incl. rng seed 20260611) when the S6.5 tests became the 2nd consumer —
textbook unify-after-two-instances.

Cross-ref [[sweep_strategy_carve]] (S1/S2 of #222), [[s6_4_dag_ownership_move]],
[[s6_4_e_walk_levelop_collapse]], [[s6_4_f_module_geography_wavefrontflux_retire]]
— S6.5 is the FINAL door-closing capstone of the #222 sequence after the (f)
module-geography landing.
