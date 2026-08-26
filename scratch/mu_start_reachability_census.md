# `mu_start` edge-extrapolation branch — reachability census (#361)

**Status: COMPLETE**
Date: 2026-08-26 · Branch `refactor/unweld-phase-b` · Host `.venv` (Py 3.14)

## 1. The object under measurement — `file:line`

The **edge-extrapolation branch** is the `else` arm of the seed dispatch in
`PoleAngularClosure.precompute_psi_state`:

- `orpheus/sn/sweep/pole_angular_closure.py:1748-1755`
  ```python
  if p in self._carrying_levels:
      ...                                     # first-class psi-half STATE
  else:
      psi_half_seed_arr = self.edge_extrapolated_seed(psi_level, p)   # ← THE BRANCH
  ```
- Body: `edge_extrapolated_seed` at `pole_angular_closure.py:1641-1670`
- Stencil (carries the **#361 `argsort` hazard**): `_edge_seed_stencil` at
  `pole_angular_closure.py:1671-1706`; the bare `order = np.argsort(mu)` is
  **line 1696**.
- A second, mirror consumer of `_carrying_levels` lives in the adjoint at
  `pole_angular_closure.py:1972` (checked separately below).

### The predicate that decides "carrying"

`PoleAngularClosure.__init__` (`pole_angular_closure.py:1536-1541`):
```python
self._carrying_levels = frozenset(
    p for p, s in enumerate(march_start_structure_per_level(quad, coord))
    if s.consumes_independent_seed
)
```
`consumes_independent_seed` (`pole_angular_closure.py:811-822`) is
`not (on_edge_node or degenerate)`, and the two conjuncts are produced by
`march_start_structure_per_level` (`pole_angular_closure.py:823-874`),
**bit-exactly, no tolerance**:

| coord | `on_edge_node` | `degenerate` |
|---|---|---|
| SPHERICAL | `mu_x[0] == -1.0` | `len(mu)>1 and mu_x[0] == mu_x[1]` |
| CYLINDRICAL (per level) | `xi[level_idx[0]] == 0.0` | `len>1 and eta[level_idx][0] == eta[level_idx][1]` |

`assert_carrying_quadrature` (`pole_angular_closure.py:893-946`) raises unless
`non_carrying_levels(starts) == ()`, i.e. unless **every** level is carrying.

## 2. ⚠ FIRST STRUCTURAL FINDING — the gate is CYLINDER-ONLY

`grep -rn assert_carrying_quadrature orpheus/` → **exactly one production call
site**: `orpheus/sn/mesh/augmented_mesh.py:347`, inside
`case CoordSystem.CYLINDRICAL:` of the coord `match`.

The `case CoordSystem.SPHERICAL:` arm (`augmented_mesh.py:350-354`) does **not**
call it. So "the population `assert_carrying_quadrature` admits" is, read
literally, a statement about the **cylinder arm only** — and taking it as the
whole population would silently drop the sphere arm, which is exactly where the
`argsort` comment says the sort is load-bearing. **Both arms are censused below.**

(more sections appended as measured)

---

## 3. THE VERDICT

> **0 of 92 admitted rules reach the edge-extrapolation branch** — where the
> 92 are every (shipped-factory rule × chart) pair that the **production
> `SNMesh` constructor** accepts, drawn from a 548-row population
> (274 constructed rule instances × 2 curvilinear charts).
>
> Restated against the literal question asked (`assert_carrying_quadrature`'s
> own admitted set, which is CYLINDER-ONLY): **0 of 88**.

⚠ **But the branch is NOT dead, and the correct disposition is NOT retirement.**
A hand-built **Gauss–Lobatto** polar rule — legal, `assert_carrying_quadrature`
never runs on the sphere arm — **builds a production `SNMesh(SPHERICAL)` and
ENTERS the branch** at `n = 6, 8, 9, 10, 11, 17` (`[M]`, 6 of 11 orders probed).
See §7.

⭐ And on the #361 `argsort` hazard specifically: **0 of 75** reachable
non-carrying levels have the tie-break live (§8).

## 4. Why the naive check is TAUTOLOGICAL, and what was measured instead

`assert_carrying_quadrature` raises iff `non_carrying_levels(starts) != ()`,
and `MorelMontryAngularSweep.__init__` computes `_carrying_levels` from the
**same producer, same `(quad, coord)` pair**. So "an admitted rule has a
non-carrying level" is a contradiction *in terms*: it can never be observed,
and a probe reporting `0` for it carries **zero information** (`[M]` I ran it
anyway: `0` on both charts, exactly as the algebra forces).

What actually needed measuring, and was:

1. **Is the gate on the constructor path at all, unconditionally?**
   `_init_core` is *"The ONE construction body both surfaces funnel into"*
   (`augmented_mesh.py:248`); all three entry points (`SNMesh(...)`,
   `from_axes`, `from_material_mesh`) route through it. ✅
2. **Is the gate's `quad` the closure's `quad`?** `assert_carrying_quadrature(
   quadrature, self.coord)` at `:347` and `angular_redistribution(
   angular_measure, coord)` (which stores `quadrature=quadrature`) receive the
   same object; the closure reads `angular.quadrature` / `angular.coord`. ✅
3. **Does the branch stay silent when the code actually RUNS?** ✅ §5.
4. **Which OTHER gates bound the population** — since the sphere has no
   carrying gate at all? §6.
5. **Would the instrument have fired?** ✅ positive controls, §9.

## 5. Execution evidence — the branch never runs on any admitted rule

Probe: `scratch/probe_361_census_02_production.py`. For every (rule, chart),
build the production `SNMesh`, then drive `precompute_psi_state` on a random
ψ with `_edge_seed_stencil` **and** `edge_extrapolated_seed` wrapped by a call
counter installed on `MorelMontryAngularSweep` itself.

| chart | rules probed | `SNMesh` BUILT | forward exec errors | levels exercised | branch entries |
|---|---:|---:|---:|---:|---:|
| SPHERICAL | 274 | **15** | 0 | 15 | **0** |
| CYLINDRICAL | 274 | **77** | 0 | 511 | **0** |
| **total** | 548 | **92** | 0 | **526** | **0** |

Every one of the 92 executed cleanly (no error silently masking a skipped
branch), and `edge_extrapolated_seed` was called **0** times in total.

## 6. The DENOMINATOR — how M was enumerated, and validated

### 6.1 Enumeration is STRUCTURAL, never a name pattern

```python
FACTORIES = {name: getattr(Quadrature, name)
             for name, val in vars(Quadrature).items()
             if isinstance(val, (classmethod, staticmethod))}
# -> ['folded_product', 'gauss_legendre', 'lebedev', 'level_symmetric', 'product']
```
`vars()` + `isinstance` cannot silently drop a member the way a regex can
(the `[a-z0-9-]+` failure mode). **5 shipped families, complete.**

### 6.2 Filter validation — a known member of every SHAPE

Asserted in the probe (it `assert`s, it does not print):

| shape expected to exist | caught by | witness |
|---|---|---|
| 1-D polar, NO level structure | `gauss_legendre` | `(4)` → N=4, `level_structure=None` |
| sphere cubature, NO level structure | `lebedev` | `(5)` → N=14, `level_structure=None` |
| level-structured, ABS_MU_Z levels | `level_symmetric` | `(4)` → N=24 |
| level-structured, full-circle product | `product` | `(4,8)` → N=32 |
| level-structured, σ_y **QUOTIENT** | `folded_product` | `(4,8)` → N=16, `folded_by=Mirror('y')` |

### 6.3 Parameter domains discovered by EXECUTION

| family | admissible parameters | kind |
|---|---|---|
| `gauss_legendre(n)` | every `n ≥ 1` | unbounded |
| `lebedev(order)` | **{3,5,…,29,31,35,41,47} — 18 members, closed roster** | finite |
| `level_symmetric(sn_order)` | **{2,4,…,18} — 9 members, closed roster** (S20/S22 refuse: no positive weights) | finite |
| `product(n_mu, n_phi)` | every `n_mu ≥ 1, n_phi ≥ 1` | unbounded |
| `folded_product(n_mu, n_phi)` | `n_mu ≥ 1`, `n_phi` **even** ≥ 2 | unbounded |

### 6.4 The ladder — and where it breaks its own pattern

- `lebedev`: **ALL 18** shipped orders (vv #31 finite-roster corollary — a
  ladder over an enumerable set is a *sample* of it).
- `level_symmetric`: **ALL 9** shipped orders.
- `gauss_legendre` n ∈ {1,2,3,4,5,6,7,8,9,11,12,13,16,17,24,32} — **odd**
  (1,3,5,7,9,11,13,17), **prime** (2,3,5,7,11,13,17), **non-power-of-two**
  (3,5,6,7,9,11,12,13,17,24).
- `product` n_mu ∈ {1,…,9,12,16} × n_phi ∈ {1,…,10,12,14,16} = 143 —
  **odd n_phi** (1,3,5,7,9) is essential here: the `degenerate` conjunct is a
  mirror η-tie on a full circle, and `on_edge_node` (ξ=0 at φ∈{0,π}) is
  parity-gated.
- `folded_product` n_mu ∈ {1,…,9,12,16} × n_phi ∈ {2,…,16 even} = 88 —
  n_phi=2 included deliberately (the degenerate all-tangential edge).

**274 rule instances × 2 charts = 548 classified rows, 0 unclassified.**

### 6.5 What bounds the population — the gate is not alone

| chart | gate | refuses |
|---|---|---|
| SPHERICAL | `_assert_tau_within_unit_interval` (#336) | 192 |
| SPHERICAL | specular-pairing (`BoundaryError`, axis 'x') | 44 |
| SPHERICAL | trace-rank "genuine mu_x cosines" | 12 |
| SPHERICAL | `_assert_alpha_dome_closes` | 11 |
| SPHERICAL | **`assert_carrying_quadrature`** | **NEVER CALLED** |
| CYLINDRICAL | **`assert_carrying_quadrature`** | **141** |
| CYLINDRICAL | `cylindrical_streaming` structure-less | 34 |
| CYLINDRICAL | `_assert_alpha_dome_closes` | 11 |
| CYLINDRICAL | trace-rank "genuine mu_x cosines" | 11 |

⟹ **the whole admitted population is two families:**

- **SPHERE: `gauss_legendre(n)`, `n ≥ 2` only.** Its carrying-ness is a
  **theorem, not a 15-point sample**: Gauss–Legendre nodes are the roots of
  `P_n` — strictly interior to `(−1,1)` and **simple**. So `mu[0] == -1.0` is
  impossible (`on_edge_node` false) and `mu[0] == mu[1]` is impossible
  (`degenerate` false), at **every** order, forever.
- **CYLINDER: `folded_product(n_mu, n_phi)`, `n_phi ≥ 4` even.** Carrying is
  the factory's advertised and gated consequence (T25: the σ_y fold is free
  exactly at the STAGGERED shift with even `n_phi`).

## 7. ⭐ The branch IS reachable — Gauss–Lobatto on the sphere

`assert_carrying_quadrature` never runs on the sphere, and admission is by
**structure, not provenance** (its own docstring), so a hand-built
`Quadrature(measure=DiscreteMeasure(...))` is a first-class member of the
admitted set. Probe `scratch/probe_361_census_04_second_gate.py`:

| rule (hand-built) | `mu[0]` | `on_edge_node` | carrying | production `SNMesh(SPHERICAL)` | branch |
|---|---|---|---|---|---|
| `lobatto(4)` | −1 | True | **False** | REFUSED (`tau_raw[3]=1.0000000000000007`) | – |
| `lobatto(5)` | −1 | True | **False** | REFUSED (`tau_raw[4]=1.0000000000000022`) | – |
| `lobatto(6)` | −1 | True | **False** | **BUILT** | **FIRES (1)** |
| `lobatto(7)` | −1 | True | **False** | REFUSED (τ) | – |
| `lobatto(8)` | −1 | True | **False** | **BUILT** | **FIRES (1)** |
| `lobatto(9)` | −1 | True | **False** | **BUILT** | **FIRES (1)** |
| `lobatto(10)` | −1 | True | **False** | **BUILT** | **FIRES (1)** |
| `lobatto(11)` | −1 | True | **False** | **BUILT** | **FIRES (1)** |
| `lobatto(12)` | −1 | True | **False** | REFUSED (τ) | – |
| `lobatto(16)` | −1 | True | **False** | REFUSED (`tau_raw[15]=1.0000000000000266`) | – |
| `lobatto(17)` | −1 | True | **False** | **BUILT** | **FIRES (1)** |

⚠ The refusals are **round-off-gated, not structural**: `tau_raw[N-1]` lands at
`1 + O(1e-15)` and the `[0,1]` membership guard is exact. So the 5 refusals are
an accident of the last ULP, not a design refusal — a re-derivation of τ, or a
different Lobatto node routine, flips them.

⭐ This matters because Gauss–Lobatto on the sphere is **live project work**,
not a hypothetical: the agent-memory topic file
`glob_vs_gl_spherical_quadrature_study.md` records it as an affordable
alternative whose blocker is exactly this — *"μ=−1 node ⟹ τ_0 = 0 ⟹ must be
straight-char"* — and `rules_sphere.py:876` names *"a future Gauss–Lobatto
sphere"* in `MarchStart.on_edge_node`'s own docstring. **The branch is the
mechanism that future rule needs.**

## 8. The #361 `argsort` hazard — measured where the branch DOES fire

`_edge_seed_stencil` (`pole_angular_closure.py:1671-1706`) returns `(m0, m1, t)`
and the seed is `(1−t)·ψ[m0] + t·ψ[m1]`. The comment's worry is that `argsort`
permutes within tied blocks, so `m1` is an arbitrary member of the next
η-group. **Two conditions must BOTH hold for that to bite: a tie, and `t ≠ 0`.**

Over every configuration that reaches the branch when admission is bypassed
(`scratch/probe_361_census_06_all13.py`) — **13 rules, 75 non-carrying levels**:

| statistic | count |
|---|---|
| levels where the tie-break is LIVE (`ties > 0` **and** `t ≠ 0`) | **0 of 75** |
| levels with `t == 0` exactly (⟹ `m1` **annihilated** from the answer) | **74 of 75** |
| levels where `m0 != ` stored slot 0 (the comment's `[M] order[0]==0`) | **0 of 75** |

The single `t ≠ 0` level is `gauss_legendre(2)` forced onto the CYLINDRICAL
chart (`t = −3.66e-01`) — a slab rule with no `LevelStructure`, refused in
production by `cylindrical_streaming`'s structure-less guard. It has **0 ties**.

**Direct mutation of the `m1` slot on the one witness that passes production
admission** (`lobatto(8)`, sphere): `(m0, m1, t) = (0, 1, +0)`;
`seed == psi[m0]` **bit-exactly**; perturbing slot `m1` by `1e3` moves the seed
by **`0.000000e+00`**; perturbing slot `m0` by `1e3` moves it by
**`1.000000e+03`** (positive control, same instrument).

Mechanism, and it is a theorem rather than a coincidence: on the sphere
`mu_start ≡ −1.0` (Hébert §3.9.4) and `on_edge_node` *means* `mu[m0] == −1.0`,
so `t = (mu_start − mu[m0])/(mu[m1] − mu[m0]) = 0/Δ = 0` **exactly**. The very
degeneracy that makes a level non-carrying is the one that annihilates `m1`.
The same identity holds on a NODE_ALIGNED cylinder level: the most-inward
azimuth is `φ = π`, where `η = −sinθ_p = mu_start`.

⟹ **the `argsort` in `_edge_seed_stencil` is un-exercised, and where the branch
IS reachable its output is annihilated.** #361 is neither a live repair nor a
clean retirement — see §10.

## 9. Positive controls (vv #17 — an all-zero verdict needs one)

`scratch/probe_361_census_03_controls.py`:

| leg | configuration | stencil calls | state |
|---|---|---:|---|
| **C3** negative | `folded_product(4,8)` CYL, as shipped | **0** | `1.458557510e+02` |
| **C2** positive | same, `_carrying_levels := frozenset()` | **4** | `8.391913149e+01` (**rel Δ 4.246e-01**) |
| **C3** negative | `gauss_legendre(8)` SPH, as shipped | **0** | `8.751460996e+01` |
| **C2** positive | same, forced non-carrying | **1** | `5.640624122e+01` (**rel Δ 3.555e-01**) |
| **C4** adjoint | `folded_product(4,8)` CYL as shipped | **0** | `seed_cells_bar = {0,1,2,3}` |
| **C4** adjoint | same, forced non-carrying | **4** | `seed_cells_bar = {}`, `\|psi_bar\|` 3.80e+02 → 1.80e+02 |

The counter is wired to the **live dispatch** (both consumers:
`precompute_psi_state:1748` and `angular_adjoint:1972`), fires when the
predicate says it should, and moves the numbers. The `0` of §5 is a
measurement, not a dead probe.

**Second, independent control** (`probe_361_census_04_second_gate.py` leg A) —
closure built DIRECTLY, `SNMesh` and the carrying gate bypassed entirely:

| class | closure BUILT | closure REFUSED | entered the branch |
|---|---:|---:|---:|
| carrying | 105 | 237 | **0** |
| NON-carrying | 13 | 193 | **13 (all of them)** |

Perfect separation: every non-carrying rule that can build a closure fires the
branch; no carrying rule ever does.

## 10. Disposition — #361 is NEITHER a repair NOR a retirement

- **Not a live repair.** No shipped-factory rule reaches the branch through
  production admission (0 of 92), and in the 75 levels reachable by bypass the
  tie-break is dead (0 of 75) because `t == 0` annihilates `m1`.
- **Not a clean retirement.** The branch is the seed mechanism the *sphere*
  arm needs for exactly the rule class the project has an open study on
  (Gauss–Lobatto, `on_edge_node`), and a hand-built Lobatto sphere rule
  **passes production admission today and fires it**. Retiring it removes the
  only seed path a `μ = −1`-noded sphere rule has.
- **The honest third state** (vv Mode 10 — *exercised-but-unconstrained*, here
  one rung lower: *unexercised*): keep the branch, and record in its docstring
  the measured fraction (`live on 0 of 92 production-admitted rules; the only
  witness is a hand-built μ=−1-noded sphere rule`) so a future audit cannot
  read it as covered. The `argsort` itself can be replaced by the producer's
  own ordering **or** left with the `t == 0` theorem stated beside it — but
  either change lands with **no witness** unless a Lobatto row is added in the
  same step (`plan-authoring` §6c).

## 11. Un-constructible members — holes in M, named

| candidate | why it is a hole |
|---|---|
| **Gauss–Lobatto polar** | no shipped factory: `grep -rin "lobatto" orpheus/` hits only PROSE (`generating_measure.py:133,145` — *"Node constraints (Gauss / Radau / Lobatto) are not families"*); `GeneratingMeasure` exposes `.gauss()` and `.on()` only. **Hand-built in §7 — it is the sharpest member of M and it FIRES.** |
| **Gauss–Radau polar** | same: named in prose, no constructor. A Radau rule with the fixed node at `μ = −1` is the same `on_edge_node` witness. NOT probed. |
| **Gauss–Chebyshev polar** | `gauss_chebyshev(n)` exists as a `DiscreteMeasure` (`numerics/measure.py:1401`) but is NOT a `Quadrature` classmethod and is NOT registry-registered (a known open item). Its nodes are interior + simple, so it is carrying by the same theorem as GL — probing it would not change the verdict. NOT probed. |
| **`Quadrature.quotient(g)` for `g ≠ Mirror('y')`** | reachable API (`quotient` is a public instance method), but the fold must be *free* on the node set; only the σ_y fold of a STAGGERED product is shipped as free. NOT swept — an open sub-population. |
| **arbitrary hand-built `Quadrature(measure=…)`** | admission is by STRUCTURE not provenance, so this population is genuinely **OPEN**, not enumerable. M is closed over *shipped factories*; it cannot be closed over all legal quadratures. §7 is the proof this matters. |
| `level_symmetric(20)`, `(22)` | refused by the factory itself (no positive weights) — not members. |

## 12. Refuted candidates (structural reason each)

1. **"`assert_carrying_quadrature` ADMIT ⟹ branch dead" as a MEASUREMENT** —
   refuted as vacuous. Both sides read `march_start_structure_per_level(quad,
   coord)`; the implication is definitional. Replaced by (a) verifying the
   linkage (same object, same coord, unconditional call site) and (b) execution.
2. **"The gate bounds the whole population"** — refuted: `grep` gives exactly
   one production call site, inside `case CoordSystem.CYLINDRICAL`. The sphere
   arm is bounded by four *other* gates and admits `gauss_legendre` only.
3. **"The closure constructor is a second, independent gate ⟹ dead twice
   over"** — refuted by measurement: **13** non-carrying rules build a closure
   directly and **all 13** fire the branch. The constructor's τ-monotonicity /
   `[0,1]` guards refuse most non-carrying rules but not all.
4. **"`product(4,8)` NODE_ALIGNED is the canonical live witness"** (implied by
   `edge_extrapolated_seed`'s own docstring: *"NODE_ALIGNED full-circle product
   rules hit its t = 0 degenerate"*) — refuted: it cannot even build a closure
   (`angular_cell_edges_per_level: level 0 is not a monotone…`), so it is not a
   witness for anything. The docstring describes a *retired* code path's
   behaviour.
5. **"a MIXED rule (some levels carrying, some not) could slip through"** —
   refuted structurally: `non_carrying_levels` returns *every* offender and the
   gate raises on a non-empty tuple, so mixed is refused a fortiori. (The
   `non_carrying_levels` docstring already says mixed is unreachable through
   every shipped factory; measured here, `0` mixed rows in 548.)
6. **`psi_view` shape `(N, ng, nx, 1)`** (from `precompute_psi_state`'s own
   docstring, *"canonical layout"*) — refuted by the production caller:
   `loss_representation/__init__.py:3589` passes `np.zeros((N, ng, nx))`, 3-D.
   The 4-D form raises `too many values to unpack` at `:1752`. **The docstring
   is present-tense wrong** — worth a separate one-line fix.

## 13. Probes

- `scratch/probe_361_census_01_enumerate.py` — population + classification (548 rows) → `scratch/_probe_361_rows.json`
- `scratch/probe_361_census_02_production.py` — production `SNMesh` admission + execution → `scratch/_probe_361_prod.json`
- `scratch/probe_361_census_03_controls.py` — positive controls C1–C4
- `scratch/probe_361_census_04_second_gate.py` — closure-direct sweep + hand-built Lobatto / duplicate-node witnesses
- `scratch/probe_361_census_05_hazard.py` — `(m0, m1, t)` + the `t == 0` annihilation mutation
- `scratch/probe_361_census_06_all13.py` — the 13 firing rules / 75 levels, hazard statistics

Run with `PYTHONPATH=. .venv/bin/python scratch/probe_361_census_0N_*.py`
(the `tests.sn._test_helpers.placeholder_materials` import needs the repo root
on the path). Single-process, no suite invocation, no file under `orpheus/` or
`tests/` touched.
