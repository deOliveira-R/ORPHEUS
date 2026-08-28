# P4.9 — pre-design measurements, 2026-08-28 (taken at P4.3's branch, gate running)

> ✅ **STATUS 2026-08-28 (post-P4.9a `7a0f434c`): rows 1–3 + the
> `mm_a_in_coeff` handing are LANDED** (the split ruling assigned them to
> P4.9a; this memo's row-1..3 radii and the cache-constant fork are
> EXECUTED history). **Rows 4–5 are P4.9b's §6b starting point** — the
> ~150 `StreamingOperator(sn_mesh)` ctor sites and ~60 `SNMesh.scheme` /
> `pole_angular_closure` reads, with the closure-construction question
> ("where do closures come from after row 5") ruled toward the PRODUCTION
> posing-head factory. ⚠ Re-verify the row-4/5 counts at P4.9b design
> time (`plan-authoring` §7) — they were measured pre-P4.9a.

Written while P4.3's exit gate ran. **Nothing here is committed to the plan
yet.** Re-verify per `plan-authoring` §7 before designing. Every enumeration
below ran WHOLE (no `| head`), AST where membership matters.

## ⭐⭐ THE SIZING FINDING — the charter's edit list is 5 rows; the tree says two
## of them are migrations, not edits

| edit-list row | `[M]` measured blast radius |
|---|---|
| 1. DD sheds the M-M block | small — the block + `cell_balance_terms` route at `diamond.py:194` |
| 2. `cell_balance_terms` retires | 1 production consumer; **4 test files, ~7 call sites** — and the twin-equivalence gates lose their subject (see below) |
| 3. protocol sheds angular members | **5 production files** read them (`diamond`, `cell_balance`, `linear_discontinuous` read `upstream_state.angular_upstream`; `loss_representation` reads `result.outgoing_angular_state`; cache/loss read `geom.c_*` = the CACHE's own fields, not the protocol's) + test surface not yet enumerated |
| 4. `StreamingOperator` gains the 4-arg ctor | ⛔ **~150 construction sites: 1 production (`coupled_system.py:441`) + ~148 in ~40 test files**, all spelling `StreamingOperator(sn_mesh)` |
| 5. `SNMesh` sheds `scheme` + `pole_angular_closure` | ⛔ **~60 production attr-reads in 11 files** — `loss_representation/__init__.py` alone ~45 `mesh.scheme` + 6 `pole_angular_closure`; `solver.py` 8; `streaming.py` 2; `dsa.py`, `mms/sn.py`, `loss_kernel_gauge.py` (4), `assembly.py`, `radial_characteristic.py` (2), `cache.py`, and ⚠ `transport/radial_characteristic_field.py:360` — an **L2 module reading the closure off the mesh** (the reach-through, again) |

⟹ rows 1–3 are one coherent small step (the unweld proper — kills the twin,
closes #407); rows 4–5 are a **re-plumbing migration** an order of magnitude
larger. **Candidate decomposition: P4.9a (rows 1–3 + the is-identity gate) then
P4.9b (rows 4–5).** ⚠ §6b check for that split: rows 1–3 do NOT change the
`StreamingOperator` signature or any `SNMesh` field, so the boundary cuts no
call-site set — the interval state is the current construction path with a
twin-free DD. Verify at design time.

## The done-when tell, run at design time (§10 discipline)

`[M]` `1.0 - tau` (all spellings) over `orpheus/`:

* `orpheus/transport/`: **exactly 1** — `diamond.py:229`, the twin. The
  done-when "`grep -c "1.0 - tau" orpheus/transport/` → 0" is precise and
  reachable.
* owner spellings (stay): `closure.py:1329` (the march), `:1593` (`c_in`),
  `:1993` (adjoint), `:936` (prose); L0 ladder: `angular_differencing.py:383,
  :504`, `balance.py:252,:267` (legitimate — the ladder's own math).
* ⚠ `sn/sweep/cache.py:377`: `mm_a_in_coeff = (1.0 - tau) / tau` — the cache
  builder derives a closure-owned CONSTANT itself rather than being handed it.
  Not the twin (a coefficient, not the update relation), and not covered by the
  done-when's `transport/` scope — but it is the same smell one notch down.
  Decide at design time whether the closure hands the cache its constants.

## The is-identity witness (§6c) — baseline confirmed

`[M]` test files naming `precompute_psi_state`: **8**; naming
`outgoing_angular_state`: **3**; **intersection EMPTY** — nothing gates the two
M-M spellings against each other, and the planned `is`-identity gate cannot
pass today (the degenerate branch reaches diamond's inline copy, which is no
closure object at all). The gate's red-then-green is the phase's witness.

## `angular_adjoint` — P4.9 does NOT touch it (that is §5b/O-3's step)

`[M]` family: `closure.py:279` (property), `:557` (ABC), `:1915` (M-M),
`:2128` (Identity). Production callers: `loss_representation/__init__.py:3769`
(reverse walk), `radial_characteristic.py:962`. ⚠ Two spy/count gates pin call
counts (`test_psi_half_coupling.py:2633`, `test_streaming_cell_transpose_relocation.py:161`)
— they survive P4.9 untouched but will need re-derivation when §5b retires the
hand-written reversal.

## The twin-equivalence gates lose their subject when `cell_balance_terms` dies

`[M]` `tests/sn/sweep/core/test_cell_balance_for_streaming.py:23` ("matches
`cell_balance_terms` exactly"), `test_cache.py:474`
(`test_cache_populator_matches_cell_balance_terms`), `test_diamond.py:1400,
:1444`, `test_ordinate_scan.py:64` — all compare the two spellings. Retiring
one twin makes each of these a `coding-standards` rewire case: re-derive what
each still tests (most become hand-written pins of the SURVIVING helper, which
is a promotion — the values stop being cross-checked and start being anchored).
Do NOT delete-and-forget; the migration is part of row 2's step.

## Where the closures COME FROM after row 5 (the open design question)

Today `augmented_mesh.py:268` builds `scheme` and `:422` builds
`pole_angular_closure` — the mesh is the composition site. After rows 4–5 the
posing head owns both choices. ⚠ The ~150 test constructions all spell
`StreamingOperator(sn_mesh)`; a test-side builder (in `tests/_test_helpers`)
that poses `(domain, codomain, spatial_closure, angular_closure)` from a mesh
is probably the migration's lever — decide whether it may live in production
(a posing-head convenience) or tests only. ⚠ The direction-supplier ruling
(P4.9's charter, 2026-08-28) constrains the spatial closure's signature:
assembled face coefficients, no bare μ.

## Row 3's test surface (enumerated whole; triage by meaning at design time)

`[M]` **12 test files** touch the protocol's angular members by attr or kwarg —
`tests/sn/sweep/core/` (cache, cell_balance_for_streaming, cell_visit_c_stamp,
diamond, discretization_scheme_protocol, ordinate_scan),
`tests/sn/sweep/slab/test_dd_recurrence.py`,
`tests/sn/verification/mms/test_mms_ld_slab.py`, and 4 under
`tests/transport/spatial/` (ld_ubld_primitive, ld_ubld_symbolic,
linear_discontinuous, scheme_reaction_rate_contract). ⚠ the kwarg filter is
deliberately wide (`c_out=` also names LD-balance kwargs) — each hit needs its
meaning read before it counts as a §6b site.

⚠ `sweep_graph.py`'s six `self.scheme` reads are its OWN field (the walk holds
a scheme reference) — downstream of `mesh.scheme` at construction, so row 5
re-plumbs where that field is FILLED, not these reads.
