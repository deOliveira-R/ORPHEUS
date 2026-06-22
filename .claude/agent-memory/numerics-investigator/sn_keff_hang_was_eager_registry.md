---
name: sn-keff-hang-was-eager-registry
description: OPEN (#212). A reported "SN heterogeneous keff non-convergence hang" was NOT a solver bug — continuous_get() eagerly walks the whole reference registry into O(minutes) Peierls mpmath solves. Solver converges in 0.3s. Bug STILL LATENT on main; lazy-by-name fix prototyped, NOT landed.
metadata:
  type: project
---

# SN keff "hang" = eager continuous-reference registry walk — OPEN (#212)

**STILL LATENT on origin/main (verified 2026-06-08):**
`orpheus/derivations/reference_values.py` still builds the registry eagerly via
`pkgutil.walk_packages` calling every module's `continuous_cases()` (lines ~187-196);
there is NO `continuous_case_builders` lazy-by-name contract. Tracked by **OPEN
Issue #212** ("Eager continuous-reference registry: continuous_get() pays
O(minutes) Peierls-Nyström cost to fetch ANY reference"). The prototyped fix was
never committed. This is the diagnostic of record for #212 — preserve as an open
breadcrumb.

## What the "hang" actually is
A reported "SN heterogeneous-2G/1G k-eigenvalue NON-CONVERGENCE hang"
(`refactor/field-role-typing`) is NOT a solver bug. The SN power/source iteration
converges to the Case reference in ~0.3 s at the exact failing config (n_per=320,
S8 GL, keff_tol=inner_tol=1e-12, max_outer=max_inner=500), |Δ|=3.4e-8 — three
orders below the test's 1e-5 gate.

Both failing tests open with `ref = continuous_get("sn_slab_1eg_2rg_S8")`.
`continuous_get` builds the WHOLE continuous-reference registry on first access via
`pkgutil.walk_packages` over `orpheus.derivations.*`, calling every module's
`continuous_cases()` EAGERLY. One hook —
`orpheus.derivations.continuous.peierls_nystrom.cases.continuous_cases`
(`_class_a_cases`) — eagerly builds 13 hollow-cyl/sph/slab references via adaptive
`mpmath.quad` eigenvalue solves (its own `capability_rows` docstring calls this
"O(minutes)"). The SN test pays that O(minutes) cost to fetch ONE unrelated SN ref.

**First-bad commit:** `cfe37ec` ("feat(derivations): topology-organised Peierls
infrastructure"), child of GOOD `d1021df`. Before it, peierls registered cheaply.
(cfe37ec confirmed present on main 2026-06-08.) Any suggested-good predating cfe37ec
by weeks (e.g. B.5.A `82e82c9`) is itself BAD — the hang predates B.5; B.5.2 is
exonerated.

## METHODOLOGY LESSON (the durable, reusable part — [[lessons-L1]] cascade order)
**A "non-convergence hang" in a test that consumes `continuous_get(...)` / any
registry-walk fixture: bound the SOLVER directly FIRST, before any bisect.**
- Tiny max_outer/max_inner, OR bypass the fixture by calling the reference's
  producing module's builder directly. If the solver converges fast, the hang is in
  the FIXTURE, not the solver. Do NOT bisect the solver — the cascade's Step-3
  fixed-source/solver isolation immediately exonerates it.
- Reproduce the registry-walk hang on a fast proxy:
  `timeout 60 python -c "from orpheus.derivations.reference_values import
  continuous_get; continuous_get('sn_slab_1eg_2rg_S8')"` (exit 124 = bad). Bisect on
  THIS proxy, not the full test.
- `continuous_all_names()` is cheap (names only); `continuous_all()` /
  `continuous_by_operator_form` / building any peierls ref is O(minutes).
- Generalises: distinguish fixture-cost from solver-cost before attributing a
  timeout to numerics. A timeout is not a convergence failure until the solver is
  bounded independently of its data fixtures.

## Proposed fix (prototyped + proven, NOT committed — for #212)
Lazy-by-name registry. Walker collects `{name: thunk}` via a new opt-in
`continuous_case_builders()` contract (peierls exposes it; names derived from
(shape,ng,nr,r0) WITHOUT building); `continuous_get(name)` materialises only the
requested ref (memoised). Eager `continuous_cases` retained for audit consumers.
Verified at prototype time: SN tests pass in ~2 s (was hang), all 28 ref names
enumerable cheaply, peierls refs still build on demand.

Promotable build-cost gates (xfail-strict that flip to xpass under the fix →
self-retiring): solver-fast-on-failing-config; SN producer fast + peierls producer
xfail-strict blocks build; slab-2G-2rg peierls case xfail-strict slow.
