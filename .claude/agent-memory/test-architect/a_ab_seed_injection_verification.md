---
name: a-ab-seed-injection-verification
description: The reusable recipe for gating A_AB — the ray→bulk seed injection (Coupled Block Operator campaign Step 3), a CELL-LOCAL RECTANGULAR coupling operator (domain=ray cells(p,-1), codomain=bulk AngularSourceSink) extracted from the in-sweep Morel-Montry seed injection. Bit-identity inheritance from (L+C).apply/(L+C).apply_transpose, the Euclidean forward↔transpose adjoint-consistency as the CORRECTNESS cross-check, the shared-closure-method blindness of the equivalence gates, the "reads-only-cells(p,-1)" asymmetry pin, the defer-forward fallback.
metadata:
  type: reference
---

# Verifying A_AB — the ray→bulk seed injection (Step 3)

**When this recipe applies.** A per-cell coupling term is welded into an
in-sweep matvec; a refactor poses it as a NAMED RECTANGULAR operator
(different domain and codomain) whose forward AND reverse each already exist
as production closure methods. The un-weld is *bit-identity-preserving by
intent*. Distinct from A_BB (an endomorphic resolvent, `test_ray_operator`
sibling) and A_BA (`RadialCharacteristicReconstruction`, a moment→ray fold):
A_AB is `ray-cells → bulk-residual`, σ-INDEPENDENT, cell-local (no spatial
march). Sibling of `aba_schur_fold_unweld_verification.md`.

## The object (worked instance)

`A_AB` (suggested class `RadialCharacteristicSeedInjection`, home
`orpheus/sn/operators/radial_characteristic.py` beside A_BB) injects the ψ½
starting-direction seed into the bulk angular recurrence (Morel-Montry,
Hébert §3.9.4). It is the `A_bs` block of `(L+C)` — pinned at `A_bs=7.505` by
`TestRegressionFloor` (A_AB is that block).

- **domain** = ray space (`sn_mesh.radial_characteristic_space`); reads ONLY
  the inward leg `radial_characteristic.cells(p, -1)`.
- **codomain** = bulk residual (`AngularSourceSink`, `(N,ng,nx)`).
- **σ-INDEPENDENT** — with bulk=0 the collision term `denom·ψ = 0` drops out;
  `A_AB.apply(s) = −(ΔA/w)·c_in·ψ_recur(s)/V` is pure geometry+angular.
  Constructor takes `sn_mesh` ONLY (contrast A_BB, which needs σ_t).
- `.apply(seed) → bulk`: `closure.precompute_psi_state(np.zeros((N,ng,nx)),
  radial_characteristic=seed)` (bulk-zeroed → isolates the seed by
  LINEARITY of the M-M recurrence) + per (carrying level, cell)
  `closure.cell_contribution(...)` → place `−upstream_numer/V[i]` into the
  bulk. **may be DEFERRED to step 4** (user-pending, like A_BB.apply).
- `.apply_transpose(v) → ray`: `numer_bar[p][:,within,i] = −v[:,ords,i]/V[i]`,
  `closure.angular_adjoint(tuple(numer_bar)) → (_, seed_cells_bar)`, place
  `seed_cells_bar[p]` onto ray `cells(p,-1)` (other ray slots zero). = the
  in-sweep `cb_minus += seed_cells_bar[p]` (`loss_representation:2846`).
  Euclidean transpose (NOT the metric `.H` — L19). `is_adjointable=True`.

## The structural-independence spine (the two-anchor keystone)

Both equivalence gates are BIT-IDENTITY INHERITANCE, not fresh value claims
— A_AB WRAPS the SAME `precompute_psi_state`/`cell_contribution`/
`angular_adjoint` the production `_OneDimScanWalk` calls. So:

- **They are BLIND to a bug in the shared closure methods** (mutate
  `cell_contribution` → BOTH A_AB and the `(L+C)` reference flip → gate 1
  stays GREEN). The equivalence gates pin A_AB's OWN glue (the loop, the
  bulk-zeroing, the seed-threading, the numer_bar build, the cells(p,-1)
  slotting) + ROUTING (Mode-11) — NOT the closure correctness.
- The closure correctness is INHERITED: `TestRegressionFloor`
  (`A_bs` pin + `welded=exact-inverse 3.5e-16` + `extract=dense-LU 5.5e-16`)
  and `test_radial_characteristic_metric.py::test_r1_transpose_is_euclidean_across_configs`
  (`‖T − Aᵀ‖ < 1e-12` including `T_sb = A_bsᵀ` = A_AB's transpose, cross-
  checked vs numpy `Afwd.T` — a structurally-INDEPENDENT ground, ≥2G+3G).
- **The CORRECTNESS cross-check A_AB adds is gate 3** (Euclidean forward↔
  transpose adjoint consistency): forward uses precompute+cell_contribution,
  transpose uses angular_adjoint — SEPARATELY implemented closure methods, so
  a mismatch between them (a sign, a tau, a missing /V in the reverse
  recurrence) reds. This is the one gate whose teeth are demonstrable by
  in-process monkeypatch of a shared method (the flip lands on ONE side only).

## The gate suite `TestA_AB_SeedInjection` (foundation; sphere-GL S4 only)

Home: `tests/sn/operators/test_psi_half_coupling.py` (`pytest.mark.foundation`).
Reuse `_sphere`/`_graded_sphere`/`_loss`/`_template`/`_dense`/`_blocks`/
`_ray_source`/`_ray_cotangent`/`_install_engine_spy` pattern. ≥2G every value
row. `-O`-safe (`pytest.fail`/`np.testing`, never bare `assert`).

| # | gate | law | oracle | tooth (in-process) | apply-dep |
|---|------|-----|--------|--------------------|-----------|
| 1 | forward-bit-id | `A_AB.apply(s).values ≡ (L+C).apply(FF(bulk=0,bnd=0,ray=s)).bulk.values` | production `(L+C).apply` (structurally indep of the A_AB CLASS; bulk=0,bnd=0 isolate A_bs by linearity) | Mode-11 counter (precompute==1 w/ `psi_view` ALL-ZERO + `radial_characteristic is s`; cell_contribution ≥ nx); `array_equal` placement | YES |
| 2 | transpose-bit-id | `A_AB.apply_transpose(v).cells(p,-1) ≡ (L+C).apply_transpose(FF(bulk=v,bnd=0,ray=0)).radial_characteristic.cells(p,-1)`; other ray slots (cells(p,+1)/corners) EXACTLY 0 | production `(L+C).apply_transpose` (ray=0 nulls the A_BB self-block transpose → isolates seed_cells_bar) | Mode-11 counter (angular_adjoint==1, `numer_bar == −v/V`); `array_equal` placement + zero-slot | NO |
| 3 | adjoint-consistency (Euclidean) | `⟨A_AB.apply(u), v⟩ = ⟨u, A_AB.apply_transpose(v)⟩` < 1e-11, ≥2 draws (PLAIN dot — NOT the V_cell metric) | forward (precompute+cell_contribution) ↔ transpose (angular_adjoint) — separately-implemented duals | patch `cell_contribution` sign → reds (apply flips, transpose doesn't); patch `angular_adjoint` seed sign → reds (transpose flips, apply doesn't) | YES (+fallback) |
| 4 | A_sb=0 triangular | seed-row/bulk-col of `(L+C)` = 0 (#284) | cross-ref `TestRegressionFloor.test_loss_operator_is_block_triangular_in_the_ray` | (in floor) | NO |
| 5 | seed-consumed asymmetry | (5a fwd) `A_AB.apply(full random seed) ≡ A_AB.apply(cells(p,-1)-only seed)` — reads ONLY the inward leg; +1/corner ignored; AND non-vacuous (`max|out|>1e-6`). (5b transpose) writes ONLY cells(p,-1) | two-seed identity + non-vacuity guard | patch `precompute` to read `cells(p,+1)` → 5a reds (only-`−1` seed has +1=0 ⇒ differs) | 5a YES(+fallback), 5b NO |
| 6 | non-carrying control | ctor over a seedless mesh (slab/cyl) RAISES, `match=` the `radial_characteristic_space is None` message | `pytest.raises` + positive control (sphere constructs) | (the raise IS the tooth; mirror A_BB ctor guard) | NO |

**Defer-forward variant** (user may realize ONLY `.apply_transpose`, mirroring
A_BB). Gates 1, 3-primary, 5a-primary are `.apply`-dependent → wrap the
`A_AB.apply` spelling in `xfail(strict=False, reason="A_AB.apply deferred to
step 4")` so they flip to xpass when the forward lands (L4; NEVER strict=True).
The transpose stays FULLY verified without the forward:
- gate 2 (bit-id vs production transpose) — always runs.
- **gate 3-FALLBACK**: replace `A_AB.apply(u)` with the production forward
  surrogate `(L+C).apply(FF(0,0,u)).bulk` → `⟨(L+C).apply(0,0,u).bulk, v⟩ =
  ⟨u, A_AB.apply_transpose(v)⟩`. Structurally sound (the surrogate uses the
  forward closure, A_AB.apply_transpose uses angular_adjoint) and pins the
  transpose against the production forward. Its `cell_contribution`-flip tooth
  still reds (surrogate flips, transpose doesn't).
- **gate 5a-FALLBACK**: same surrogate — `(L+C).apply(0,0,full)` ≡
  `(L+C).apply(0,0,cells(-1)-only)`.
- gate 5b + the metric-file `test_r1` cross-ref carry the transpose asymmetry.

## Refutations pre-empted (vv-principles)

1. **The sphere carries ONE level (R12a)** → the per-carrying-level loop is
   length 1; broadcast≡per-level is a COINCIDENCE. A multi-level indexing bug
   is UNTESTABLE (cylinder is non-carrying → A_AB doesn't exist there). NOTE
   the blind spot; the loop structure is inherited single-source, verified
   end-to-end on the sphere. (Same as A_BA refutation #1.)
2. **≥2G mandatory** — the seed varies per group; a group-axis bug is nulled
   at 1G (anti-#3). Use `_sphere(ng=2)`; a 4G row strengthens asymmetry.
3. **Non-trivial seed** — fixed-seed RANDOM `cells(p,-1)`, varying per cell AND
   group, so the M-M recurrence's per-cell propagation is exercised (a
   constant/zero seed under-drives it).
4. **GRADED mesh** for the value rows (`_graded_sphere`) — uniform V/ΔA nulls a
   cell-index drift (Mode-5); the grading breaks it (§0.6). A_AB has no spatial
   march, but the cell LOOP + ÷V[i] index still benefit.
5. **Structural independence** — the forward reference `(L+C).apply` does NOT
   route through the A_AB CLASS (production uses precompute+cell_contribution
   directly). It SHARES those closure methods → bit-id inheritance (necessary-
   NOT-sufficient, L13); the value is anchored by the floor + metric `test_r1`.
6. **Mode-11** — the WRAP sentinels prove A_AB EXECUTES the closure methods
   (precompute==1 with the seed threaded + bulk zeroed; angular_adjoint==1 with
   `numer_bar==−v/V`), NOT merely that the gate is green.
7. **Mode-12 SAFE** — gate 3 is EUCLIDEAN (plain dot, no null space). Do NOT
   use `_recip_defect`/`_v_cell_seed` (the V_cell metric `A.H` = the COMPOSITE
   Hilbert adjoint, verified separately in `test_radial_characteristic_metric.py`;
   A_AB's own transpose is the Euclidean `Aᵀ`, L19). The metric reciprocity's
   Mode-12 closure (ERR-067, G_sd=V_cell) is what makes the composite gate
   catch A_bs=A_AB errors — cross-ref, don't re-derive.
8. **σ-independence** — verify at the prompt: with bulk=0, `denom·ψ=0` and
   `spatial_upstream=|μ|A·0=0`, so `(L+C).apply(0,0,s).bulk = −angular_numer/V`
   (σ cancels, `0+x=x`, `0·d=0` exact) ⇒ gate 1 is `array_equal` (0-ULP), not
   principled-equiv. A_AB's ctor needs NO σ_t.

## Mutation ledger (target file:line → expected RED)

- `pole_angular_closure.py:1325` `cell_contribution` `upstream_numer_term` sign
  flip → **gate 3 REDS** (apply flips; transpose via angular_adjoint doesn't).
  gate 1 GREEN (shared method — documents the bit-id-inheritance blindness).
- `pole_angular_closure.py:1407` `angular_adjoint` `seed_bar` sign flip →
  **gate 3 REDS** (transpose flips; apply doesn't). gate 2 GREEN (shared).
- `pole_angular_closure.py:1280` `precompute_psi_state` read `cells(p,+1)` not
  `cells(p,-1)` → **gate 5a REDS** (only-`−1` seed has +1=0 ⇒ differs). gate 1
  GREEN (shared).
- A_AB forgets to zero `psi_view` (passes bulk) → **gate 1 Mode-11 sentinel
  REDS** (`psi_view` not all-zero).
- A_AB passes `radial_characteristic=None` (drops the seed) → precompute seeds
  carrying levels at ZERO → **gate 1 REDS** (A_AB.apply=0 vs non-zero ref) +
  sentinel (`radial_characteristic is None`).
- A_AB slots `seed_cells_bar[p]` on `cells(p,+1)` → **gate 2 REDS** (`array_equal`
  vs production which slots `−1`). Divergent reimplementation bypassing the
  closure → **gate 1/2 Mode-11 counter 0**.
- construct on a slab/cyl `SNMesh` → **gate 6 RAISES** the `None`-space message.

## Result contract

Land the transpose+structure gates (2, 4, 5b, 6) + gate 3-fallback LIVE now;
the `.apply`-spelling gates (1, 3-primary, 5a-primary) `xfail(strict=False)`
until the forward lands. Every demonstrable tooth mutation-verified in-process
(monkeypatch, NEVER `git checkout` an uncommitted file). Mirror A_BB's foundation
suite shape; the convergence-ORDER/physics claim (if any A_AB has) belongs in an
L1 sibling, not here (L9 — don't stack `verifies` on foundation).
