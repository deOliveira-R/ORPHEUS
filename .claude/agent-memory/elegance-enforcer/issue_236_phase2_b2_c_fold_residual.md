---
name: issue-236-phase2-b2-c-fold-residual
description: #236 Ph2 STEP B2 elegance verdict — fold P2 (DD.residual) c_in/c_out onto CellVisit + #226 Protocol→ABC matvec binding fix
metadata:
  type: project
---

# #236 Phase 2 STEP B2 — fold P2 (`DD.residual`) c onto `CellVisit` + ABC-expansion binding fix

**Verdict: PASS-WITH-NITS** (none blocking; uncommitted, branch `feature/sn-spatial-angular-product`). BIT-IDENTICAL 0-ULP.

Folds the 2nd of 4 c_in/c_out rebuild sites (crosswalk `.claude/plans/issue_236_phase2_tau_carve_crosswalk.md` §7). B1 folded P3 (sweep_cache CollisionCache); B2 folds P2 = `DiamondDifference.residual` (matvec twin @ n_mask=1). P1 (`cell_balance.py:313-314` = DD.update solve path) is the TRACKED B3 deferral; P0/producer = M-M `__init__:989-990`.

## The mechanism (RIGHT)
DD is deliberately STATELESS (reads only CellVisit/UpstreamState, never the closure) → c travels as DATA on the frozen+slotted `CellVisit` (new `c_in/c_out: float=0.0`, `scheme.py:274-277`). NEW single stamp `SNMesh._make_cell_visit` (`geometry.py:1294`) that ALL 4 dag_walk yield sites funnel through (Pattern 2 — one c-lookup, not N), reading `closure.c_{in,out}_per_ordinate[global_ordinate]` (the B1 accessor). Global-ordinate index = `direction_idx` (slab/sphere) / `level_indices[p][m]`=`global_n` (cylinder) — the SAME index `streaming_terms` resolves (verified geometry.py:1418 + _make_cell_visit use global_n). DD.residual reads `visit.c_out/c_in`, dropped dead `tau`/`alpha` locals; assembly byte-unchanged.

## ⭐ KEY RULINGS
1. **ABC-expansion = RIGHT, RESOLVES the B1-flagged asymmetry, NOT scope creep.** This is the exact ARCH-OPP I filed in B1 (the two mint-on-both bindings: `transverse_coupling_is_facewise` consumer→ABC vs `pole_angular_closure` consumer→Protocol). B2 retyped `SNMesh.pole_angular_closure` Protocol→ABC at 4 annotation sites → forced declaring `precompute_psi_state`/`cell_contribution`/`angular_adjoint` `@abstractmethod` on the ABC (they lived ONLY on concretes). Matches the cited `DiscretizationSchemeBase` precedent (declares update/residual abstract). Both concretes implement all 3 → `__abstractmethods__` empty → instantiable (verified). `psi_state: Any` for contravariance is correct + minimal.

2. **SSOT verified at the ONE source.** matvec route (`MorelMontryAngularSweep.cell_contribution:1236-1238` reads `self._c_{in,out}_per_level[p]`) and sweep route (`DD.residual`←`visit.c`←`_make_cell_visit`←`c_*_per_ordinate`←`_gather_per_ordinate(self._c_*_per_level)`) BOTH single-source the same `_c_*_per_level` instance arrays. Cross-route equivalence tests green (11 matvec≡sweep). 460 sweep-core green, 4 xfail (#195), 62 diamond+cell_balance green.

3. **Defaulted-0.0 = sound.** IdentityAngularClosure yields neutral zeros → default 0.0 IS the correct slab value (free correctness); ~25 direct test CellVisit(...) constructions still build; round-trip honesty forced the test surrogates (see NIT-2).

## NITS (all follow-up, none block commit)
- **NIT-1 (the ASSEMBLY-DUPLICATION the brief asked about) — ACCEPTABLE for B2, flag as follow-up.** The 2-line `dAw·c_out` / `dAw·c_in·psi_ang` assembly exists TWICE: `DD.residual` inline (diamond.py:310-321, n_mask=1) and `cell_contribution` (pole_angular_closure.py:1239-1244, vectorized n_mask). The c-FORMULA is now single-sourced (B2's win); only the thin (ΔA/w)-scaled assembly differs in granularity (scalar-vs-vectorized). This is the `diamond.py:289-303` TODO (route DD.residual through `closure.cell_contribution`), explicitly a Phase-6 cleanup target. Pattern-2 inst-#2 (twin matvec/fold, verified-identical, acceptable-for-now): leaves are single-sourced, byte-verified by the 11 equivalence tests. STOPS being acceptable if a 3rd assembly site appears OR an edit lands on one-not-other. Latent hazard: a future metric-weighted (ΔA/w) reform landing in one assembly only. Collapse destination = DD.residual routes through cell_contribution.
- **NIT-2 (test-surrogate twin).** `_visit_c` (test_cell_balance) + `_c_from_streaming_terms` (test_diamond) BOTH recompute `c_out=α_out/τ`, `c_in=(1-τ)/τ·α_out+α_in` — the very formula B2 retired from prod — to stamp test visits (closure not in fixture scope). Two copies of the retired formula now in test code. JUSTIFIED (L11: the surrogate IS the independent reference the round-trip pins prod against; importing the prod closure would make the round-trip vacuous). But 2 copies = unify-after-two trigger → extract ONE `_visit_c` test helper (module-shared). Follow-up.
- **NIT-3 (Protocol now orphaned — ARCH-OPP, file/track).** `PoleAngularClosure` (@runtime_checkable Protocol, pole_angular_closure.py:191) is now used in ZERO production type-annotations/isinstance — only docstrings + `__init__` export. The ABC won; the parallel Protocol contract is dead weight (two contracts for one concept, must be kept in sync by hand — divergence habitat). Retire the Protocol OR document why both survive. This is the natural close of the B1 ARCH-OPP. Out of B2 scope; track for Step-C tail.

## Crosswalk status (anti-#11 OK)
P1/cell_balance.py = B3 (tracked, deferred, DD.update:230 5th-τ-consumer). No untracked leak. `sweep_cache.py:323 mm_a_in_coeff=(1-τ)/τ` is the WDD half-angle `a_in`, NOT a c-constant — separate concern, not a leak.
