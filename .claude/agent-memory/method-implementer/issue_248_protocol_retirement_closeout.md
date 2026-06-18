# Issue #248 — retire the orphaned `PoleAngularClosure` Protocol + dead legacy `__call__`/`tau_mm` surface

**Branch** `feature/sn-spatial-angular-product` · 2026-06-18 · NOT committed (main agent commits after review).
**Type** RETIREMENT + TEST-MIGRATION carve (no new numerics; deletion + re-pin of dead-path tests onto the live path).
**Owed by** #236 Phase 2 (B2 retyped consumers Protocol→ABC + made strategy methods abstract on the ABC; B3/Step-C closeouts both named #248 as OWED).

## What was retired (the architecture: Cardinal Rule 2 — one contract per concept)

After #236 Phase 2 B2, the `@runtime_checkable PoleAngularClosure` **Protocol** was orphaned AND
divergent: every production consumer (matvec / sweep / geometry / scheme / cell-balance) was retyped
onto the `PoleAngularClosureBase` **ABC**, and the three strategy methods
(`precompute_psi_state` / `cell_contribution` / `angular_adjoint`) were made `@abstractmethod` on the
ABC. The Protocol carried `is_linear` / the `c_*`/`tau` accessors / `__call__` but NOT the three
strategy methods. The legacy bundle `__call__` (with its `tau_mm` argument) had ZERO production
callers — the live matvec/sweep use `precompute_psi_state` + `cell_contribution`, never `__call__`.

DELETED (all in `orpheus/sn/spatial/pole_angular_closure.py`):
- `class PoleAngularClosure(Protocol)` + `@runtime_checkable` decorator (the whole block).
- `MorelMontryAngularSweep.__call__` (the legacy bundle) + the now-orphaned staticmethod
  `_weighted_angular_recurrence_single_level` + the dead instance method `_redistribution_for_level`
  (zero call sites — only `__call__` and the deleted `_weighted_...` reached it).
- `IdentityAngularClosure.__call__` (slab zeros).
- The ABC `__call__` `@abstractmethod` stub.
- `Protocol` / `runtime_checkable` from the `typing` import (now unused as code symbols).
- `PoleAngularClosure` from this module's `__all__`, from `orpheus/sn/spatial/__init__.py` import + `__all__`.

KEPT (re-verified they feed the LIVE path, grepped callers before touching):
- `_psi_half_grid_single_level` (staticmethod, the recurrence kernel) — consumed by
  `_psi_half_grid_for_level` (live) AND the public `compute_psi_half_per_level`.
- `_psi_half_grid_for_level` (instance) — consumed by `precompute_psi_state` (the live matvec/sweep path).
- `compute_psi_half_per_level` (public PR-TYPED-6b method) — NOT in the deletion set; calls the kept kernel.
- `morel_montry_tau_per_level` (the τ producer, #236 Step-A/C) — untouched.

## Test migration (the load-bearing piece — vv Mode 11: the dead `__call__` gate NEVER executed the live path)

The hand-calc tests pinned the verbatim Hébert §3.9.4 M-M weighted-DD redistribution algebra
(sphere `R_0 = 2/√3`, `R_1 = -2/√3`; cylinder per-level `R = ±1.0`/`±6.3`; linearity) — but through
the DEAD `MorelMontryAngularSweep()(psi, α, dAw, τ, V)` bundle. Re-pinned onto the LIVE surface:

`tests/sn/sweep/curvilinear/test_pole_angular_closure.py`:
- Module docstring/import: dropped bare `PoleAngularClosure`; kept `MorelMontryAngularSweep`,
  `PoleAngularClosureBase`. `TestProtocolConformance.test_morel_montry_satisfies_protocol`:
  `isinstance(MorelMontryAngularSweep(), PoleAngularClosure)` → `PoleAngularClosureBase`; docstring
  reworded off "structural typing / no inheritance" → ABC conformance.
- NEW module helper `_redistribution_via_live_path(closure, psi_level, α, dAw, τ, V)` — reconstructs the
  SAME `R_m = (ΔA/w)_{i,m}/V_i·(α_{m+1/2}ψ_{m+1/2} − α_{m-1/2}ψ_{m-1/2})` the dead bundle returned, but
  the half-angle ψ-thread comes from the LIVE `closure.compute_psi_half_per_level(psi_level, τ)` (the
  same `_psi_half_grid_single_level` kernel the matvec's `precompute_psi_state` consumes) and the test
  owns the explicit α/ΔA/w/V geometry fold. **WHICH live method pins which hand-calc quantity:**
  `compute_psi_half_per_level` pins the half-angle recurrence (ψ_{3/2}=2, ψ_{5/2}=4 on the sphere
  fixture); the explicit fold pins the α-weighted redistribution R_n. Works on an UNBOUND
  `MorelMontryAngularSweep()` (the kernel is a staticmethod / arg-passing public method — no
  `_require_mesh_bound`).
- `TestMorelMontryHandCalc` (sphere 2-ord), `TestCylindricalLevelDispatch` (the 2-level cylinder is
  reconstructed per-level: one `_redistribution_via_live_path` call per level on the level's ordinate
  slice, scattered back to global slots — mirrors how `precompute_psi_state` loops `level_indices`
  calling `_psi_half_grid_for_level`), and `test_strategy_is_linear_in_psi` all driven through the
  helper. SAME expected redistribution values, rtol 1e-14.

`tests/sn/sweep/curvilinear/test_compute_psi_half_per_level.py`:
- MIGRATED `test_recurrence_output_random_seed` off `_weighted_angular_recurrence_single_level`: now
  pins the recurrence formula DIRECTLY against the verbatim Hébert ψ_{m+1/2}=(ψ_m−(1−τ_m)ψ_{m-1/2})/τ_m
  via the live `compute_psi_half_per_level` grid (plus the geometry fold reconstruction).
- DELETED `test_redistribution_from_psi_half_matches_call` (a `__call__`-vs-helper cross-check —
  vacuous once `_weighted_...` is gone: it verified two dead/helper paths agree). Its companion
  `test_method_delegates_to_psi_half_grid_helper` STAYS (pins `compute_psi_half_per_level ==
  _psi_half_grid_single_level`, both KEEPs).
- Dropped the now-unused `_mm_weighted_angular_recurrence_single_level` import alias.

**MUTATION-CHECK (the migration genuinely asserts the algebra, not vacuously green):**
- Sphere hand-calc: perturbed `expected` `2/√3 → 2/√3·1.0000001` → REDS (`DESIRED 1.1547` vs produced
  `1.1547001`); reverted by re-edit (NOT git checkout, L28) → green.
- Cylinder hand-calc: perturbed level-1 `6.3 → 6.30001` → REDS; reverted → green.

## Docstring / cross-ref repointing (NO behaviour change)

Repointed 8 `:class:`/`:attr:`/`:meth:`PoleAngularClosure...` cross-refs → `PoleAngularClosureBase`:
`operator.py:508`, `solver.py:2056`, `geometry.py:1315`, `scheme.py:244/258/267`, `cell_balance.py:143`,
`psi_half_angle_seed.py:174` (the last reworded "Protocol's call signature" → "ABC's strategy
contract"). `operator.py:146` TYPE_CHECKING import of `PoleAngularClosure` was docstring-only (no
annotation usage) → DROPPED. In-file: rewrote the module-docstring "Architectural mirror" section to
the current ABC-only reality (+ #248 retirement-history bullet); updated the ABC docstring's "Subclasses
MUST declare" + the one-line-registration example off `__call__` → the three strategy methods; repointed
the `_MMHalfGrid` "redistribution body" docstring + the M-M class comment + the `_require_mesh_bound`
error message + the `__init__`/staticmethod docstrings off the deleted `__call__`/helper. Repointed the
2 dangling `:meth:_weighted_angular_recurrence_single_level` Sphinx refs in `psi_half_angle_seed.py`
(module docstring + `PsiHalfAngleSeed` Protocol docstring) → the surviving `_psi_half_grid_single_level`.

## VERIFICATION (literal stdout)

Migrated test files (`.venv/bin/python -O -m pytest <two files> -p no:cacheprovider -q`):
`41 passed, 1 warning in 0.05s` (the `-O` bare-assert warning is the standard ORPHEUS canonical-invocation
notice; all asserts are `np.testing.*`). Identical 41-count before AND after the production deletion →
the migration is bit-stable across the carve. (43 → 41 from the single deleted test +
parametrization; reconciles exactly.)

Live curvilinear MMS admission gate
(`tests/sn/verification/mms/test_curvilinear_operator_admits_mms.py`): `2 passed`.
Broader live-path sanity: `tests/sn/spatial/` + both migrated files `117 passed`; live curvilinear
sweep/cache gates (`test_tau_producer_equivalence` + sph/cyl sweep regressions +
`test_streaming_equilibrium_curvilinear` + `test_w1_clamp_silent_on_flat` (ERR-026 catcher) +
`test_coupled_pole_mu_level_invariant` + `test_psi_half_angle_seed`) `94 passed`. The production
matvec/sweep path (consumes `precompute_psi_state` + `cell_contribution` + the `c_*`/`tau` accessors,
never `__call__`) is fully intact.

Import check: `MorelMontryAngularSweep()` is NO LONGER callable (`callable(mms) == False`; calling it
raises `TypeError: ... object is not callable`; `'__call__' not in MorelMontryAngularSweep.__dict__`).
`PoleAngularClosure` gone from `orpheus.sn.spatial`; `PoleAngularClosureBase` present.

pyright `pole_angular_closure.py`: **26 errors** (down from the ~31/54 #226 baseline the B1/B3
closeouts recorded — the deletion removed `|None`-carrying code). ALL 26 are pre-existing #226-class
noise (`reportArgumentType`/`reportOptionalSubscript`/`reportOptionalMemberAccess`/
`reportAttributeAccessIssue` on the `tuple|None` per-level annotations subscripted after
`_require_mesh_bound()` pyright can't narrow). ZERO reference any deleted/repointed symbol. ZERO new.

FINAL `grep -rn "PoleAngularClosure\b" orpheus/ tests/ | grep -v PoleAngularClosureBase`: 7 remnants,
ALL intentional narrative — 4 in `pole_angular_closure.py` (#248 retirement comments at :102/:105/:198 +
the `:doc:` theory-section-name ref at :149, flagged below) and 3 in the migrated test file (docstrings
describing the retired Protocol). ZERO live-symbol references.

## Sphinx — EXPLICITLY OUT OF SCOPE per the brief ("I handle Sphinx; FLAG for archivist")

No `docs/theory/*` edits, no theory-label edits, NO docs stub written (the brief overrides the default
method-implementer Sphinx-stub manifest item, exactly as the #236 Step-C brief did). FLAGGED for the
archivist (dangling `:class:`/`:meth:`/`:attr:`PoleAngularClosure cross-refs to the now-deleted
Protocol — these are SOURCE `.rst`, will break the Sphinx build until repointed to
`PoleAngularClosureBase`):
- `docs/theory/discrete_ordinates.rst` lines 211, 964, 1035, 1037, 1098, 1100, 1143, 1145, 1245, 1258,
  1260, 1424, 1433, 5371, 5382, 5476, 5571, 5682, 5693, 5830, 7208 (the `:attr:`c_*/tau accessor refs +
  the section heading "PoleAngularClosure (Issue #168 Phase B)" that the in-file `:doc:` at
  `pole_angular_closure.py:149` points at).
- `docs/theory/operator_algebra.rst` lines 1077, 1306, 1339, 1373, 1633 (`:meth:`cell_contribution +
  `:class:` refs).
- The theory label `sn-pole-angular-closure-protocol` (if present) — brief said do not touch; archivist
  to rename/repoint when expanding.

## DISPATCH

archivist DISPATCH_REQUEST emitted (followup:false) for the `docs/theory/*` `PoleAngularClosure`→
`PoleAngularClosureBase` cross-ref cleanup + the "Protocol retired by #248" narrative.

## FOLLOW-UP (qa-caught missed call site, 2026-06-18 — `a51ae0288e1960cdb`)

The L12-grep MISSED ONE `__call__` site because it was VARIABLE-BOUND (`pac(...)`, `pac =
sn_mesh.pole_angular_closure`), not `MorelMontryAngularSweep()(...)` nor a `PoleAngularClosure`-name
hit. Site: `tests/sn/sweep/curvilinear/test_unified_matvec_cylinder.py:170` inside the L0 hand-reference
helper `_hand_reference_cyl_matvec`. It SLIPPED VERIFICATION because its sole consumer
`test_unified_cylinder_matches_hand_reference` is `@xfail(strict=False, reason="#206")` — the post-deletion
`TypeError` (pac no longer callable) was SWALLOWED as an xfail (27 xfailed in 0.61s = errored at the
`pac(...)` line, never ran the 27-case matvec). A dead test masquerading as a tracked xfail, xfailing for
the WRONG reason.

FIX (same migrate-then-re-pin pattern): replaced `redist_full = pac(psi_g_first, ..., level_indices=...,
carlson_context=...)` with a per-level loop over `level_indices` calling the SHARED helper
`redistribution_via_live_path` (HOISTED `_redistribution_via_live_path` from `test_pole_angular_closure.py`
into `tests/sn/_test_helpers.py`, +`carlson_context` kwarg; both files now import it — Cardinal Rule 2, no
2nd private copy). Each level passes its `reduced.alpha_per_level[p]` / `reduced.redist_dAw_per_level[p]` /
`tau_mm_per_level[p]` / V + the per-level `carlson_ctx_per_level[p]` the helper ALREADY built — mirroring how
production `precompute_psi_state` loops levels with each level's Carlson coupled-pole seed. Removed the dead
`if pac is None: pac = MorelMontryAngularSweep()` fallback (a cylindrical SNMesh ALWAYS binds
`MorelMontryAngularSweep(self)`) + the now-unused `MorelMontryAngularSweep` import.

PROOF the reconstruction is faithful + load-bearing: (1) in-process `np.array_equal` of the new
reconstruction vs a faithful replica of the DELETED `_weighted_angular_recurrence_single_level` kernel
(`_psi_half_grid_single_level(seed=Carlson) ∘ α/ΔA/w/V fold`) = `True`, max-diff EXACTLY 0.0 on the
product(2×4) cylinder — BYTE-IDENTICAL to what `__call__` produced. (2) `--runxfail` now surfaces 27
AssertionErrors (the documented #206 `assert_allclose` divergence, values like 3.87e+02) — NOT a TypeError;
the matvec comparison actually runs. Normal run = 27 clean `xfailed` (right reason now). (3) Mutation: scaling
the helper output `×1.0000001` shifts the hand-reference `out[]` by 1.07e-06 > 0 → `redist_full` is genuinely
consumed in `streaming+redistribution+collision`, not vacuous. GATES: cylinder file 2P/27xf/2desel; closure
file 13P (hand-calc sphere/cyl values unchanged); full curvilinear non-slow 166P/27xf/0F; `tests/sn/spatial/`
76P. Final grep `pac(`/`pole_angular_closure(`/`MorelMontryAngularSweep()(` = ZERO. pyright: zero NEW error
classes (the 2 flagged lines are `reduced.alpha_per_level[p]` = the SAME pre-existing #226 `reduced | None`
noise the HEAD `pac(... reduced.alpha_per_level ...)` call already carried).

## LESSON

A pure-RETIREMENT carve of a dead `__call__`-style bundle interface has a test-migration blast radius
the audit's named-test list under-counts in TWO directions: (1) every dangling DOCSTRING
`:meth:`/`:class:` cross-ref to a deleted symbol becomes a Sphinx-build break, and these hide in
NEIGHBOUR modules (here `psi_half_angle_seed.py`'s module docstring + Protocol docstring pointed at the
deleted `_weighted_angular_recurrence_single_level`) — L12-grep `_weighted_angular...`/`_redistribution_for_level`
across the WHOLE prod tree (not just the closure file) before deleting, and repoint each to a surviving
sibling. (2) ⭐ a VARIABLE-BOUND call site (`pac(...)` where `pac = obj.attr`) is INVISIBLE to both the
`ClassName()(...)` grep AND the deleted-symbol-NAME grep — and when that site sits inside an `xfail`'d test,
the post-deletion `TypeError` is SWALLOWED (the test xfails for the WRONG reason, a dead test masquerading as
a tracked one), so the carve's own gate run shows nothing. To catch BOTH: after deleting a `__call__`/dunder,
ALSO grep variable-bound invocations `[^a-zA-Z_.]\bvarname\(` for the common closure handle names
(`pac|closure|sweep|op|...`) repo-wide, AND treat every `xfail(strict=False)` test that TOUCHES the deleted
surface as suspect — run it under `--runxfail` and confirm it reaches its REAL assertion (the documented
reason), not an upstream `TypeError`/collection error; a sub-second wall-clock on a multi-case matvec xfail is
the tell it errored early. The migration win is vv-Mode-11-shaped: the hand-calc algebra was real and valuable
but pinned the DEAD `__call__` that the rewired production never executes — re-pin it through the SAME
recurrence kernel the live `precompute_psi_state` consumes (`compute_psi_half_per_level` →
`_psi_half_grid_single_level`), keeping the verbatim expected values + mutation-verifying RED (and, for a
migrate-then-delete oracle relocation, in-process `array_equal` vs a faithful replica of the deleted kernel),
so the algebra is now verified through the path production actually runs.
