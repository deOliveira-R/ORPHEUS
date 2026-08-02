# Q3 — Gauss-rule retirement blast radius (audit)

**Branch:** `refactor/operator-strategy-layers`
**Tree state at audit open:** clean for `orpheus/` and `tests/` — the only
working-tree modifications are `.claude/skills/vv-principles/{SKILL,error_catalog}.md`
plus untracked `scratch/_from_instruction_arch_worktree/`. **No in-flight carve**
touches any audited file, so every `file:line` below is at HEAD.

**Method:** four searches per target (graph callers via Nexus, text grep across
`orpheus/` + `tests/` + `docs/`, `__all__`/re-export scan, doc-node scan). Line
numbers are current at this HEAD; re-derive via Nexus if the tree has moved.

---

## 0. NAMING HAZARD — read this before any grep-based edit

There are **two distinct symbols spelled `gauss_legendre`** in the production tree
and they must never be conflated:

| Spelling | Kind | Retirement target? |
|---|---|---|
| `orpheus.numerics.measure.gauss_legendre(n)` | module-level function, `measure.py:1064` | **YES — target 1** |
| `Quadrature.gauss_legendre(n_ordinates=16)` | **classmethod**, `quadrature/directional.py:530` | **NO** |

A bare `grep -rn "gauss_legendre" tests/` returns ~570 lines, of which **450 are
`Quadrature.gauss_legendre(...)`** — the classmethod, which is the SN slab
quadrature entry point and is emphatically NOT retiring. The true `tests/` blast
radius of target 1 is **two files**. Any grep in the retirement's execution must
be anchored (`^|[^.]gauss_legendre\(` or an import-line grep), never bare.

A third, unrelated spelling lives in the oracle tree:
`orpheus.derivations.common.quadrature.gauss_legendre(a, b, n, *, dps=53)`
(`orpheus/derivations/common/quadrature.py:249`) — mpmath-backed, different
signature, explicitly out of scope.

---

## 1. `orpheus/numerics/measure.py::gauss_legendre(n)` — line 1064

### (a) GRAPH CALLERS

`mcp__nexus__callers` → **29 callers, ALL in `tests/`. ZERO production callers.**

- **Production (`orpheus/`): NONE.**
- **Tests:** 28 in `tests/numerics/test_measure.py`, 1 in
  `tests/numerics/test_measure_phase.py:85`.
- **Docs:** the symbol is referenced from `docs/theory/foundations/discrete_measures.rst`
  (see (b)); no `implements`/`cites` edge.

### (b) TEXT GREP

**`orpheus/` (production, excl. `derivations/`)** — every hit is the definition
itself or a re-export; there is **no call site**:

| file:line | nature |
|---|---|
| `orpheus/numerics/measure.py:1064` | definition |
| `orpheus/numerics/measure.py:1097` | its own `ValueError` message |
| `orpheus/numerics/measure.py:1133` | docstring xref **from `gauss_chebyshev`** (`:func:`gauss_legendre``) — dangles if target 1 retires and target 2 does not |
| `orpheus/numerics/measure.py:1229` | `__all__` |
| `orpheus/numerics/__init__.py:14`, `:90` | package re-export + `__all__` |

**`tests/`** — 2 files:

- `tests/numerics/test_measure.py` — import at `:30`; ~40 call sites (`:52, :117,
  :140, :148, :151(×2), :159, :172, :173, :192, :193, :215, :216, :238, :255,
  :266, :313, :325, :334, :346, :355, :387, :411, :442, :507, :520, :533, :543,
  :553, :563, :573, :579, :587, :596`).
- `tests/numerics/test_measure_phase.py` — import at `:22`; used at `:84` inside a
  `parametrize` lambda (`lambda: gauss_legendre(5)`).

**`docs/` (excl. `_build`)** — 2 live references to the module-level symbol:

- `docs/theory/foundations/discrete_measures.rst:345` —
  `:func:`~orpheus.numerics.measure.gauss_legendre`` in the "1-D primitive
  constructors" list; `:347` states its `degree_of_exactness = 2n-1` / weight sum 2.
- `docs/theory/foundations/discrete_measures.rst:357` — a **prose code sample**:
  ``μ_S2 = gauss_legendre(n_mu) * equispaced(0.0, 2*np.pi, n_phi)``. Plain literal,
  no role; no build check of any kind.

> **Sphinx-gate note.** `:345`/`:348`/`:352` are **Python-domain** roles
> (`:func:`). Per the project's measured behaviour these render as plain text with
> **no `-W` warning** unless the build runs `-n`. The doc blast radius of this
> retirement is therefore **NOT caught by the Sphinx gate** — it must be repointed
> by hand.

### (c) RE-EXPORTS — yes, two levels

- `orpheus/numerics/measure.py:1229` — `__all__` entry `"gauss_legendre"`.
- `orpheus/numerics/__init__.py:14` — `from orpheus.numerics.measure import (… gauss_legendre …)`.
- `orpheus/numerics/__init__.py:90` — `__all__` entry `"gauss_legendre"`.

So the symbol is part of the **public `orpheus.numerics` surface**.

### (d) TEST CLASSIFICATION

`tests/numerics/test_measure.py` splits sharply into two groups. Only ~6 tests are
*about* `gauss_legendre`; the other ~30 use it as a **convenient fixture** for
testing `DiscreteMeasure` algebra (`*`, `|`, `pushforward`, `restrict`,
`with_metadata`, `__iter__`, `__getitem__`, `repr`, `BundleMeasure`).

| file:line | test | class | justification |
|---|---|---|---|
| `tests/numerics/test_measure.py:43` | `test_gauss_legendre_polynomial_exactness` | **BEHAVIORAL** | asserts the 2n−1 exactness contract against `sympy`-free monomial integrals — this is the correctness contract of the rule; rewire to the successor. |
| `tests/numerics/test_measure.py:113` | `test_gauss_legendre_weight_sum_is_two` | **BEHAVIORAL** | pins Σw = 2, i.e. the *generating measure* (Lebesgue on [−1,1]) — exactly the property the Golub-Welsch parameterization must reproduce; rewire. |
| `tests/numerics/test_measure.py:504` | `test_constructor_invalid_args_raise` (`:507` GL arm, `:509` GC arm) | **BEHAVIORAL** | pins the `n < 1` domain guard; the successor must keep an equivalent guard — rewire (or fold into the generic constructor's guard test). |
| `tests/numerics/test_measure.py:138, 146, 156, 168, 187, 212, 236, 245, 310, 322, 331, 342, 353, 384, 400, 436, 518, 527, 541, 551, 561, 572, 578, 586, 594` | `test_n_points_property`, `test_dim_property_1d_and_2d`, `test_integrate_vector_valued_integrand`, `test_tensor_product_*`, `test_direct_sum_*`, `test_pushforward_*`, `test_restrict_*`, `test_with_metadata_*`, `test_bundle_measure_*`, `test_call_*`, `test_integrate_*`, `test_iter_*`, `test_len_*`, `test_getitem_*`, `test_repr_*` | **fixture-only** (neither of the three buckets) | these test `DiscreteMeasure`, not the rule. They need a **one-line substitution** of the constructor, not migration or deletion. Do not count them in the retirement's "test migration" scope beyond a mechanical rename. |
| `tests/numerics/test_measure.py:255, 266` | `test_direct_sum_composite_quadrature` (uses `gauss_legendre(4).pushforward(...)`) | **fixture-only** | same. |
| `tests/numerics/test_measure_phase.py:84` | `test_untagged_generic_rule_phase_raises` | **BEHAVIORAL (about a different SUT)** | parametrized over `[lambda: gauss_legendre(5), lambda: gauss_chebyshev(4)]` — the contract under test is "an *untagged* rule raises on `.phase`". The two constructors are chosen **because they carry no `invariance_group`**. If the successor tags rules, this test's premise changes; flag it. |

**No API-SMOKE and no CHARACTERIZATION tests exist for target 1.** There is no
frozen `.npz`/hardcoded-node baseline anywhere for `measure.gauss_legendre`.

---

## 2. `orpheus/numerics/measure.py::gauss_chebyshev(n)` — line 1107

### (a) GRAPH CALLERS

`mcp__nexus__callers` → **4 callers, ALL in `tests/`. ZERO production callers.**

- `tests/numerics/test_measure.py:68` `test_gauss_chebyshev_polynomial_exactness`
- `tests/numerics/test_measure.py:123` `test_gauss_chebyshev_weight_sum_is_pi`
- `tests/numerics/test_measure.py:504` `test_constructor_invalid_args_raise`
- `tests/numerics/test_measure_phase.py:85` `test_untagged_generic_rule_phase_raises`

### (b) TEXT GREP

**`orpheus/`** — definition (`measure.py:1107`), its own error string (`:1149`),
the equation label `discrete-measure-gauss-chebyshev` (`:1114`), `__all__`
(`:1228`), and the package re-export (`orpheus/numerics/__init__.py:13, :89`).
**No call site anywhere in production.**

**`tests/`** — only the two files above (`test_measure.py:29` import, `:86`,
`:128`, `:509`; `test_measure_phase.py:21` import, `:84`).

**`docs/`** — one reference:
`docs/theory/foundations/discrete_measures.rst:348` (`:func:` role) + `:350`
(the prose describing its weight function and `degree_of_exactness`). Same
silent-Python-domain-role hazard as target 1.

> **Equation-label hazard.** `measure.py:1114` carries
> `:label: discrete-measure-gauss-chebyshev`. That is a **`math:equation` node**.
> If any page `:eq:`-references it, retiring the docstring breaks the reference —
> and unlike `:func:`, an `:eq:` failure IS a build warning. Grep confirmed **no
> `:eq:` consumer** of that label, so it retires cleanly, but the successor's
> docstring should carry the derivation forward (Cardinal Rule 3).

### (c) RE-EXPORTS — yes, two levels

- `orpheus/numerics/measure.py:1228` — `__all__`.
- `orpheus/numerics/__init__.py:13`, `:89` — import + `__all__`.

### (d) TEST CLASSIFICATION

| file:line | test | class | justification |
|---|---|---|---|
| `tests/numerics/test_measure.py:68` | `test_gauss_chebyshev_polynomial_exactness` | **BEHAVIORAL** (`@l1`, `@verifies("discrete-measure-integrate")`) | asserts the **weighted** exactness against the closed-form `∫x^k(1−x²)^{−1/2}dx = π(k−1)!!/k!!`. This is the ONE place in the tree that pins the weighted semantics — it is the load-bearing evidence for section (e). Rewire to the successor; it becomes the Chebyshev-measure arm of the generic construction. **Note `:88` also asserts the tag (`rule.degree_of_exactness == 2*n-1`) — a tag-pin, not a property (lessons L-016); the property itself IS tested here, so the tag-pin is redundant, not vacuous.** |
| `tests/numerics/test_measure.py:123` | `test_gauss_chebyshev_weight_sum_is_pi` | **BEHAVIORAL** | Σw = π *is* the total mass of the generating measure `(1−x²)^{−1/2}` on [−1,1]. In the Golub-Welsch re-pose this becomes `μ_0` — the first moment of the weight — so it is exactly the assertion the new parameterization must reproduce. Rewire. |
| `tests/numerics/test_measure.py:504` (`:509` arm) | `test_constructor_invalid_args_raise` | **BEHAVIORAL** (shared with target 1) | domain guard `n < 1`. |
| `tests/numerics/test_measure_phase.py:84` | `test_untagged_generic_rule_phase_raises` | **BEHAVIORAL (different SUT)** | GC is used here as an *untagged* rule fixture; the contract under test belongs to `.phase`, not to GC. See target-1 note. |

---

## 3. `orpheus/numerics/measure.py::equispaced(a, b, n)` — line 1161

### (a) GRAPH CALLERS

`mcp__nexus__callers` → **8 callers, ALL in `tests/numerics/test_measure.py`.
ZERO production callers.**

`:224, :236, :287, :400, :436, :474, :484, :495`.

### (b) TEXT GREP

**`orpheus/`** — definition (`:1161`), two error strings (`:1198`, `:1200`),
`__all__` (`:1227`), package re-export (`orpheus/numerics/__init__.py:12, :88`).
**No call site.** Every other `orpheus/` hit for the word "equispaced" is
**English prose in a docstring** describing the φ-grid, never the symbol:
`quadrature/rules_product.py:2, :25, :73, :100, :162`;
`quadrature/directional.py:208, :598`; `quadrature/__init__.py:32`;
`quadrature/registry.py:475`; `sn/loss_representation/__init__.py:2893, :3672`.
(This is the L-009 substring trap in a mild form — "equispaced" as a common
adjective vs. `equispaced(` as a symbol. Anchor the retirement grep on
`equispaced(`.)

**`tests/`** — `tests/numerics/test_measure.py` only (import `:28`, calls at the
8 sites above plus `:226, :227, :239, :296, :414, :447`).
All other `tests/` hits are prose (`test_rules_product.py:3`,
`test_symmetry.py:845`, `test_quadrature_directional.py:77`,
`test_spherical_harmonic_space.py:286`, `test_g_adjoint_reciprocity.py:135`,
`sn/primitives/test_quadrature.py:8`, `sn/sweep/core/test_sweep_regression.py:245`).

**`docs/`** — one `:func:` reference at
`docs/theory/foundations/discrete_measures.rst:352` (+ `:354` prose,
+ `:357` the code sample, + `:95` prose). All other doc hits are prose.

### (c) RE-EXPORTS — yes, two levels

- `orpheus/numerics/measure.py:1227` — `__all__`.
- `orpheus/numerics/__init__.py:12`, `:88` — import + `__all__`.

### (d) TEST CLASSIFICATION

| file:line | test | class | justification |
|---|---|---|---|
| `tests/numerics/test_measure.py:474` | `test_equispaced_midpoint_integrates_constants_exactly` | **BEHAVIORAL** | pins `∫1 = b−a` for the midpoint rule — the degree-0 half of the `degree_of_exactness = 1` claim. |
| `tests/numerics/test_measure.py:484` | `test_equispaced_midpoint_integrates_linears_exactly` | **BEHAVIORAL** | pins `∫x = (b²−a²)/2` — the degree-1 half. Together these two ARE the exactness contract. |
| `tests/numerics/test_measure.py:495` | `test_equispaced_invalid_args_raise` | **BEHAVIORAL** | pins BOTH guards (`n >= 1` and `a < b`); note `a < b` is a guard the other two primitives do not have. |
| `tests/numerics/test_measure.py:224, 236, 287, 400, 436` | `test_direct_sum_concatenates_on_shared_space`, `test_direct_sum_rejects_mismatched_spaces`, `test_pushforward_invertible_map_change_of_variables`, `test_bundle_measure_separable_smoke`, `test_bundle_measure_varying_fiber` | **fixture-only** | `equispaced` is the convenient *arbitrary-interval* rule (it is the only primitive that takes `(a, b)`). `:287` in particular uses `equispaced(0,1,200)` specifically because the midpoint rule converges slowly enough to make the change-of-variables identity a real test. A successor must still offer an `(a,b)`-parameterized rule or these need re-fixturing. |

**No API-SMOKE and no CHARACTERIZATION tests.**

---

## 4. `orpheus/numerics/quadrature/rules_1d.py::gauss_legendre_on_mu(n)` — line 33

**This is the only one of the four with LIVE PRODUCTION CONSUMERS.**

### (a) GRAPH CALLERS

`mcp__nexus__callers` → **7 callers: 1 production + 6 tests.**

- **Production:** `orpheus/numerics/quadrature/directional.py:530`
  `Quadrature.gauss_legendre` classmethod (the call is at `:541`).
- **Tests:** all 6 in `tests/numerics/test_rules_1d.py` (`:32, :48, :55, :65,
  :89, :114`).

> ### ⚠ THE GRAPH MISSES THE SECOND PRODUCTION CONSUMER
>
> `orpheus/numerics/quadrature/registry.py:517` holds
> `factory=gauss_legendre_on_mu` inside the `QuadratureSpec(name="GaussLegendre1D", …)`
> entry of `quadrature_registry` (`registry.py:514-525`). This is a **captured
> function object in a dataclass field**, not a call — so it produces **no
> `calls` edge** and `mcp__nexus__callers` reports it nowhere. It is nonetheless a
> fully live production path: `select_quadrature("slab", target_degree=…)` reaches
> the rule through `spec.factory(**params)`.
>
> Consequence for the carve: **anything that swaps this symbol must patch the
> registry's captured reference too** (the spec is a frozen dataclass — the known
> technique is `object.__setattr__`). Text-grep, not the graph, is the primary
> evidence for this symbol.

### (b) TEXT GREP

**`orpheus/`** — 6 sites:

| file:line | nature |
|---|---|
| `orpheus/numerics/quadrature/rules_1d.py:33` | definition |
| `orpheus/numerics/quadrature/rules_1d.py:78` | its own error string |
| `orpheus/numerics/quadrature/rules_1d.py:98` | `__all__` |
| `orpheus/numerics/quadrature/directional.py:87` | `from .rules_1d import gauss_legendre_on_mu` |
| `orpheus/numerics/quadrature/directional.py:541` | **call site** (inside `Quadrature.gauss_legendre`) |
| `orpheus/numerics/quadrature/__init__.py:80`, `:96` | re-export + `__all__` |
| `orpheus/numerics/quadrature/registry.py:226` | import |
| `orpheus/numerics/quadrature/registry.py:517` | **captured `factory=` reference (graph-invisible)** |
| `orpheus/numerics/__init__.py:35`, `:91` | package re-export + `__all__` |

**`tests/`** — 3 files import it:
`tests/numerics/test_rules_1d.py:20`, `tests/numerics/test_registry.py`,
`tests/numerics/test_symmetry.py` (`:656`, `:927`, `:938` — the registry-spec /
declared-symmetry rows).

*Indirect* test consumers via `Quadrature.gauss_legendre(...)` number **450 lines
across `tests/`**, plus `tests/numerics/test_symmetry.py:307` and `:321`. These do
not name the symbol but every one of them executes it.

**`docs/`** — one `:func:` reference:
`docs/theory/foundations/discrete_measures.rst:776` (+ `:778` describing its
`SO(2)` tag and `degree_of_exactness = 2n − 1`).

### (c) RE-EXPORTS — yes, THREE levels

- `orpheus/numerics/quadrature/rules_1d.py:98` — `__all__`.
- `orpheus/numerics/quadrature/__init__.py:80` (import), `:96` (`__all__`).
- `orpheus/numerics/__init__.py:35` (import), `:91` (`__all__`).

### (d) TEST CLASSIFICATION — `tests/numerics/test_rules_1d.py`

| file:line | test | class | justification |
|---|---|---|---|
| `:32` | `test_gauss_legendre_returns_discrete_measure` | **API-SMOKE + TAG-PIN** | asserts `isinstance`, `n_points`, `nodes.shape`, `weights.shape`, `support`, then `invariance_group == SO2` and `degree_of_exactness == 2n−1`. The shape/type half is smoke → delete. The two metadata asserts are **tag-pins, not property tests** (lessons L-016) — they assert what the constructor wrote, so they cannot fail unless the constructor changes. Keep only if the successor is expected to carry the same tags, and then as a **metadata-contract** test, not evidence of exactness. |
| `:48` | `test_gauss_legendre_weight_sum_is_two` | **BEHAVIORAL** | `Σw = 2 = ∫_{-1}^{1}dx` — the total mass of the generating measure (`μ_0` in Golub-Welsch). Rewire. |
| `:55` | `test_gauss_legendre_rejects_zero` | **BEHAVIORAL** | domain guard for `n ≤ 0`. Rewire. |
| `:65` | `test_gauss_legendre_bit_identical_to_legacy_adapter` | **nominally CHARACTERIZATION — but see the warning below; in fact it is ROUTE-EQUIVALENCE and will not bite** | see ⚠ below. |
| `:89` | `test_gauss_legendre_polynomial_exactness` | **BEHAVIORAL** (`@l1`) | closed-form `∫x^k dx` over `k ≤ 2n−1`. The real exactness contract. Rewire — this is the test that must remain green through the Golub-Welsch swap. |
| `:114` | `test_gauss_legendre_fails_above_degree_of_exactness` | **BEHAVIORAL (anti-vacuity)** | explicitly written to prove `:89` is not vacuous by showing failure at `k = 2n`. This is the sharpness gate; keep and rewire — it is exactly the guard that would catch an over-claimed `degree_of_exactness` in the successor (cf. the `level_symmetric_sn` over-claim, #327). |

> ### ⚠ `test_gauss_legendre_bit_identical_to_legacy_adapter` (`:65`) DOES NOT BITE
>
> It self-describes as "the **load-bearing contract** for the Issue 4 refactor: if
> the rule's nodes/weights drift even at the last bit, the SN regression snapshots
> will silently mis-compare", and uses `np.array_equal` (exact bits).
>
> **Measured structure:** the LHS is `gauss_legendre_on_mu(n)`; the RHS is
> `Quadrature.gauss_legendre(n)`, which at `directional.py:541` **calls
> `gauss_legendre_on_mu(n_ordinates)`**. Both sides are computed in the same
> process from the same function. This is a **route-equivalence assertion, not a
> frozen right-hand side** — under a Golub-Welsch swap both sides move together
> and the test stays GREEN while the nodes change. It is immune to the very drift
> it claims to guard.
>
> The **actual** frozen baselines that would move are the `.npz`/`.npy` snapshots
> that embed GL nodes: `tests/sn/sweep_ref_2g.npy`,
> `tests/numerics/data/step6_product_reroute_baseline.npz`,
> `tests/geometry/snapshots/bc_equivalence_white_xmin_partial_03_GL.npz`, and the
> `tests/sn/_data/bc_extraction_baseline/*_SLB_*.npy` family. Those are the real
> characterization surface for target 4 — **and none of them name the symbol**, so
> a grep-driven audit alone would miss them.
>
> Verdict: `:65` is an **API-SMOKE-equivalent** in its current form. Either delete
> it with the symbol, or (better) convert it to a genuine characterization test
> against a stored node/weight baseline before the swap.

---

## (e) THE `degree_of_exactness` SEMANTIC SPLIT

### The split, MEASURED (not inferred)

Probe run against the current tree (`.venv/bin/python`, `orpheus.numerics.measure`).
"unweighted-exact to k" = largest `k` such that every monomial `x^j`, `j ≤ k`,
reproduces `∫_{−1}^{1} x^j dx`; "weighted-exact to k" = same against
`∫_{−1}^{1} x^j (1−x²)^{−1/2} dx`. `−1` means **not even exact for constants**.

| n | TAG (both rules) | GL unweighted | GL weighted | GC unweighted | GC weighted |
|---|---|---|---|---|---|
| 3 | 5 | **5** | −1 | **−1** | **5** |
| 4 | 7 | **7** | −1 | **−1** | **7** |
| 5 | 9 | **9** | −1 | **−1** | **9** |

The two rules carry the **same integer** and are exact in **disjoint** senses. GC
is not merely "less accurate" unweighted — it fails at degree 0 (Σw = π ≠ 2). The
field as it stands is a **number without a stated measure**.

### Consumers that READ `.degree_of_exactness`

**Production — 5 read sites, all in 2 files:**

| file:line | consumer | affected by a weight-relative reinterpretation? |
|---|---|---|
| `orpheus/numerics/measure.py:467-470` | `__mul__` (tensor product) → `new_dx = min(self.dx, other.dx)` | **YES — this is the primary risk site.** `min()` over two degrees is sound only when both operands share a weight function. It performs **no weight reconciliation**. |
| `orpheus/numerics/measure.py:528-531` | `__add__` (direct sum) → `new_dx = min(self.dx, other.dx)` | **YES — same defect, and it is already demonstrably wrong.** |
| `orpheus/numerics/measure.py:693` | `consolidate` → `degree_of_exactness=self.degree_of_exactness` (preserved) | **NO** — pure pass-through, weight-agnostic. |
| `orpheus/numerics/measure.py:899` | `with_metadata` → pass-through/override | **NO** — pure pass-through. |
| `orpheus/numerics/quadrature/registry.py:872` | `params = spec.degree_of_exactness_for(target_degree)` | **INDIRECTLY YES.** This reads the **spec's inverter callable**, not a measure field. But the selector's whole contract is "give me a rule of degree ≥ d", and today all four registry entries are **unweighted-on-their-domain** rules, so the contract is implicitly unweighted. A weighted rule entering the registry would need a weight in the query, or the selector silently returns an incompatible rule. |

> #### MEASURED latent defect — the propagation rule is already unsound
>
> ```
> GL4 + GC4  (__add__):  degree_of_exactness = 7   (min(7,7))
>                        ∫1 dμ  = 5.141592653589793   (should be 2)
> GL4 * GC4  (__mul__):  degree_of_exactness = 7   support = '[-1,1] × [-1,1]'
>                        ∫1 dμ  = 6.283185307179586   (should be 4)
> ```
>
> Both composites **claim degree 7 while failing at degree 0**. Nothing in the
> code prevents this composition and no test covers it. This is not a
> consequence of the planned redesign — it exists today, and it is the concrete
> argument that the field MUST become weight-relative. **The `min()` rule at
> `measure.py:470` and `:531` needs a "same generating measure" precondition
> regardless of how the Gauss constructors are re-posed.**
>
> No production code path currently constructs such a mixture (GC has zero
> production callers — see target 2), so this is **latent, not live**.

**Tests — 18 read sites.** Classified by whether a weight-relative reinterpretation
changes their meaning:

| file:line | reads | affected? |
|---|---|---|
| `tests/numerics/test_measure.py:55` | `rule.degree_of_exactness == 2n−1` (GL) | NO — GL is unweighted; the number is unchanged. |
| `tests/numerics/test_measure.py:87` | `rule.degree_of_exactness == 2n−1` (GC) | **YES — this is the assertion whose MEANING changes.** The surrounding test already verifies the *weighted* integral, so the code is right and only the field's contract needs to say so. |
| `tests/numerics/test_measure.py:219` | `p.degree_of_exactness == 5` after `gauss_legendre(3) * gauss_legendre(2)` | NO — both operands unweighted, `min(5,3)`… **verify: asserts 5, so it is `min` over the operands as constructed.** Same-weight composition, sound. |
| `tests/numerics/test_measure.py:315, 318, 349` | pushforward/restrict drop the field to `None` | NO — drop semantics are weight-independent. |
| `tests/numerics/test_measure.py:667` | `out.degree_of_exactness == 3` (consolidate, hand-built measure) | NO. |
| `tests/numerics/test_measure_partition.py:237, 243` | field present on parent / `None` on partitions | NO. |
| `tests/numerics/test_rules_1d.py:43` | `m.degree_of_exactness == 2n−1` | NO — GL, unweighted. |
| `tests/numerics/test_rules_product.py:43` | `== min(2n_mu−1, n_phi−1)` | NO — product rule is unweighted on `S²`. |
| `tests/numerics/test_rules_sphere.py:39, 111` | `== order` / `== sn_order−1` | NO (unweighted). *(Separately: `:111` is the tag-pin behind the #327 over-claim.)* |
| `tests/numerics/test_registry.py:220, 243` | `measure.degree_of_exactness == 15` after `select_quadrature(..., target_degree=15)` | **YES, indirectly** — these pin the selector's degree contract. If the field becomes weight-relative, the selector query must state a weight (or the registry must declare all entries unweighted) for these to keep meaning what they mean now. |
| `tests/numerics/test_registry.py:77`, `test_symmetry.py:938` | `spec.degree_of_exactness_for(...)` | **YES, indirectly** — same reason. |
| `tests/numerics/test_spherical_harmonic_basis.py:120` | `deg = measure.degree_of_exactness; L = deg // 2` | **YES — the sharpest semantic consumer.** It derives the safe SH order from the field: `ℓ+ℓ' ≤ deg`. That inference is **only valid for an unweighted cubature on `S²`**. Currently only reached with `lebedev_sphere`, so it is safe today — but it is the consumer that would silently change meaning if a weighted rule ever flowed in. |

**Docs — 1 semantic consumer:**
`docs/theory/foundations/spherical_harmonics.rst:293` — "*For a discrete angular
cubature … whose `degree_of_exactness` is at least `2L`, the no-prefactor real
`Y_ℓ^m` satisfy the discrete orthogonality …*". This states an **unweighted**
premise without saying so. It is the doc-side twin of
`test_spherical_harmonic_basis.py:120` and should gain the weight qualifier.

`docs/theory/foundations/discrete_measures.rst:347/350/354` already state the
three primitives' degrees, and `:350` **already says "(in the weighted sense)"**
for `gauss_chebyshev` — so the *documentation* is ahead of the *type*. That
parenthetical is the existing evidence that the split was known and left
unencoded.

### Verdict on (e)

- **No production consumer would break**, because no production code reads a
  `degree_of_exactness` produced by a weighted rule today (GC has zero production
  callers).
- **Two production sites are unsound by construction** and would be *fixed*, not
  broken, by making the field weight-relative: `measure.py:470` and `:531`.
- **The one consumer that would silently change meaning is
  `tests/numerics/test_spherical_harmonic_basis.py:120`** (and its doc twin
  `spherical_harmonics.rst:293`), which derives `L = deg // 2` on the unstated
  assumption that `deg` is unweighted.
- **The registry selector (`registry.py:872`) is the place where a stated weight
  must be threaded**, since "give me degree ≥ d" is a weight-relative request.

---

## (f) `equispaced` — IS IT USED BY THE PRODUCT RULE'S φ AXIS?

**No. `equispaced` has ZERO callers in the entire production tree, and the product
rule does not use it.** Its own docstring is the only place that connects it to
the φ axis, and it connects them **contrastively** ("*the project's existing code
uses left-endpoints, but this primitive offers midpoints*").

### The evidence

`product_mu_phi` builds its φ axis inline at
`orpheus/numerics/quadrature/rules_product.py:127-130`:

```python
# Equispaced φ in [0, 2π) — left-endpoints, matching the legacy
# ProductQuadrature contract pinned by the regression snapshots.
phi_pts = np.linspace(0, 2 * np.pi, n_phi, endpoint=False)
w_phi = 2.0 * np.pi / n_phi
```

Measured, `n_phi = 4`:

```
equispaced(0, 2π, 4).nodes           = [0.785, 2.356, 3.927, 5.498]   (MIDPOINTS)
np.linspace(0, 2π, 4, endpoint=False)= [0.000, 1.571, 3.142, 4.712]   (LEFT-endpoints)
identical? False
```

So the docstring's claim is **accurate**: the two conventions differ, and the
production path uses the left-endpoint one, which `rules_product.py:128` pins as a
regression-snapshot contract.

### Three further truths worth recording

1. **`equispaced` is used by NOTHING but its own tests.** 8 graph callers, all
   `tests/numerics/test_measure.py`. It is a **speculative primitive that never
   acquired a consumer** — the strongest retire-cleanly case of the four.

2. **The product rule does not even build a tensor product of measures.**
   `product_mu_phi` returns nodes of shape `(n_mu·n_phi, 3)` — **Cartesian
   direction cosines on `S²`** (`support='S^2'`) — assembled by an explicit double
   loop (`rules_product.py:139-156`), not by `DiscreteMeasure.__mul__`. So even a
   left-endpoint `equispaced` would not slot in without a change of representation.

3. **⚠ A LIVE DOC FALSITY.** `docs/theory/foundations/discrete_measures.rst:357`
   states:

   > ``μ_S2 = gauss_legendre(n_mu) * equispaced(0.0, 2*np.pi, n_phi)`` gives a
   > product quadrature on `μ ∈ [-1,1] × φ ∈ [0,2π)` — the standard
   > polar-azimuthal split for transport on the unit sphere.

   Measured, that expression yields `nodes (8, 2)`, `support='[-1,1] × [0.0,6.28…]'`
   — a 2-D measure on a **rectangle in `(μ, φ)`**, whereas the production `S²` rule
   yields `nodes (8, 3)`, `support='S^2'`. They are **different objects**, and the
   sentence's name `μ_S2` asserts they are the same. This is a doc bug **today**,
   independent of the retirement; the retirement makes it un-runnable as well as
   untrue.

4. **Pre-existing dangling xrefs, retired for free.** `orpheus/sn/quadrature`
   **does not exist** in this tree (verified: no such path; no `class
   GaussLegendre1D` or `class ProductQuadrature` anywhere in `orpheus/` or
   `tests/`). The following are therefore **already-dangling** Python-domain
   xrefs — silent, per the no-`-n`-build rule:
   - `orpheus/numerics/measure.py:1092` — `:class:`orpheus.sn.quadrature.GaussLegendre1D``
   - `orpheus/numerics/measure.py:1178` — `:class:`orpheus.sn.quadrature.ProductQuadrature``
   - `orpheus/numerics/quadrature/rules_1d.py:73` — same `GaussLegendre1D`
   - `orpheus/numerics/quadrature/rules_product.py:116` — same `ProductQuadrature`
   - `tests/numerics/test_rules_1d.py:6` (module docstring)

   Retiring targets 1, 3 and 4 removes four of these five for free. The
   `rules_product.py:116` one survives and needs separate repair.

   `measure.py:1093-1094` additionally says "*The Issue 4 adapter will rebuild it
   on top of this primitive in Wave B*" — a **future-tense promise that never
   happened** and is the historical reason `measure.gauss_legendre` has no
   production caller.

---

## (g) OTHER IN-TREE GAUSS/CHEBYSHEV NODE COMPUTATIONS (candidate twin paths)

### Production tree

| file:line | expression | verdict |
|---|---|---|
| `orpheus/numerics/measure.py:1098` | `np.polynomial.legendre.leggauss(n)` | **TWIN PATH #1** — target 1. |
| `orpheus/numerics/quadrature/rules_1d.py:79` | `np.polynomial.legendre.leggauss(n)` | **TWIN PATH #2** — target 4. **Byte-for-byte the same statement as #1**, differing only in the `DiscreteMeasure` metadata attached (`invariance_group=SO2` here, absent there). |
| `orpheus/numerics/quadrature/rules_product.py:125` | `mu_gl, w_gl = np.polynomial.legendre.leggauss(n_mu)` | **TWIN PATH #3 — production, and NOT on the retirement list.** The cylindrical/2-D product rule computes its own GL polar axis rather than calling `gauss_legendre_on_mu(n_mu)`. Three independent `leggauss` call sites in one package. **Flag for the campaign: retiring #1 and #2 into one Golub-Welsch construction leaves #3 as a surviving twin unless it is re-routed.** |
| `orpheus/numerics/measure.py:1150-1152` | `np.cos((2i−1)π/2n)` / `np.full(n, π/n)` | **TWIN PATH #4** — target 2; the closed-form Chebyshev rule, trig-evaluated. (Related to the #325 trig-vs-algebraic node-generation thread.) |
| `orpheus/numerics/measure.py:1202-1203` | `a + (arange(n)+0.5)·h` / `np.full(n, h)` | target 3; not a Gauss rule. |
| `orpheus/numerics/quadrature/rules_product.py:129` | `np.linspace(0, 2π, n_phi, endpoint=False)` | the φ axis — the left-endpoint convention discussed in (f). Not a Gauss rule; a **twin of `equispaced` in concept but not in convention**. |
| `orpheus/transport/radial_characteristic_field.py:299` | `np.polynomial.legendre.legvander(mu_p, ords.size − 1)` | **NOT a twin path.** `legvander` builds a Legendre **Vandermonde matrix** (basis evaluation for moment projection); it generates no nodes or weights. Leave alone. |

**No `scipy.special.roots_*`, no `roots_chebyt/chebyu`, no `eigh_tridiagonal`, and
no existing Golub-Welsch machinery anywhere in the production tree.** The
Golub-Welsch construction would be new code, not a consolidation onto an existing
implementation.

### Oracle / derivations tree (`orpheus/derivations/**`) — DO NOT TOUCH

Reported per brief, for awareness only. These are the deliberately
structurally-independent oracles.

| file:line | expression |
|---|---|
| `orpheus/derivations/common/quadrature.py:227-242` | `_leggauss(n, dps)` — dispatches to `np.polynomial.legendre.leggauss` at `dps == 53` (`:239`) and to **`mpmath.gauss_quadrature(n, "legendre")`** above it (`:242`). This is the arbitrary-precision oracle. |
| `orpheus/derivations/common/quadrature.py:249` | `gauss_legendre(a, b, n, *, dps=53)` — the oracle's GL on `[a,b]`. |
| `orpheus/derivations/common/quadrature.py:268` | `gauss_legendre_visibility_cone(...)` |
| `orpheus/derivations/common/quadrature.py:332` | `composite_gauss_legendre(...)` |
| `orpheus/derivations/common/quadrature.py:527` | `np.polynomial.laguerre.laggauss(n)` inside `gauss_laguerre` |
| `orpheus/derivations/common/shifted_legendre.py:77` | `np.polynomial.legendre.leg2poly(...)` — basis conversion, not nodes |
| `orpheus/derivations/discrete/sn/balance.py:360` | `mu_gl, w_gl = np.polynomial.legendre.leggauss(N)` — a **fifth** independent GL computation, in the oracle tree. Legitimate: an oracle that called the production rule would lose its independence. |

### ⚠ ONE PRODUCTION → DERIVATIONS IMPORT EXISTS

The brief asked to report any production import from the oracle tree. **There is
one, and it is a Gauss-Legendre import:**

- `orpheus/cp/solver.py:44` —
  `from orpheus.derivations.common.quadrature import composite_gauss_legendre`
- used at `orpheus/cp/solver.py:170` —
  `q = composite_gauss_legendre(self.mesh.edges.tolist(), p.n_quad_y)`
- with explanatory comments at `orpheus/cp/solver.py:100` and `:166`.

So the **production CP solver consumes the oracle tree's GL construction**. That
is not one of the four targets and nothing in this retirement touches it, but it
means the derivations GL path is **not purely an oracle** — it is load-bearing for
`orpheus/cp`. Any future "make derivations strictly independent" pass must
reckon with this edge. (Corroborating breadcrumb:
`orpheus/derivations/common/quadrature.py:347` says "*`_composite_gauss_legendre`
in the CP solver — both are now [unified]*", i.e. this consolidation was
deliberate.)

Additionally, `orpheus/derivations/discrete/sn/contamination.py:214` and
`orpheus/derivations/continuous/mms/sn.py` (8 sites: `:229, :591, :1923, :2007,
:3169, :3403, :3478, :3674`) call **`Quadrature.gauss_legendre(...)`** — the
production classmethod — i.e. the dependency also runs derivations → production.
Since that classmethod wraps target 4, **the MMS reference solvers execute
`gauss_legendre_on_mu` transitively**. They are not textual consumers of the
symbol but they are runtime consumers of its numbers.

---

## SUMMARY

| target | production callers | test files | doc refs | re-exported | retire cleanly? |
|---|---|---|---|---|---|
| 1 `measure.gauss_legendre` | **0** | 2 (`test_measure.py`, `test_measure_phase.py`) | 2 (`discrete_measures.rst:345, :357`) | `measure.__all__`, `numerics/__init__` | **YES** — no production reachability |
| 2 `measure.gauss_chebyshev` | **0** | 2 (same) | 2 (`discrete_measures.rst:348, :350`) | same | **YES** |
| 3 `measure.equispaced` | **0** | 1 (`test_measure.py`) | 2 (`discrete_measures.rst:352, :357`) | same | **YES** — never acquired a consumer |
| 4 `rules_1d.gauss_legendre_on_mu` | **2** (`directional.py:541` call; `registry.py:517` captured factory — graph-invisible) | 3 (`test_rules_1d.py`, `test_registry.py`, `test_symmetry.py`) + 450 indirect lines via `Quadrature.gauss_legendre` | 1 (`discrete_measures.rst:776`) | 3 levels of `__all__` | **NO** — live path; needs the registry capture patched and the SN snapshot surface checked |

### Items the four-search audit found that a `callers()` query alone would have missed

1. `registry.py:517` `factory=gauss_legendre_on_mu` — a **captured function object**, zero graph edges, fully live.
2. `rules_product.py:125` — a **third production `leggauss`** that is not on the retirement list and survives the consolidation as a twin.
3. `docs/theory/foundations/discrete_measures.rst:357` — a **false code sample** that no build checks.
4. Four **already-dangling** `orpheus.sn.quadrature.*` xrefs in the audited docstrings.
5. `orpheus/cp/solver.py:44` — a **production → derivations** GL import.
6. `test_rules_1d.py:65` — a self-declared "load-bearing bit-identity contract" that is **route-equivalence and will not bite**.
7. The **measured** `GL + GC` / `GL * GC` unsoundness at `measure.py:470` / `:531`.

---

## RECONCILIATION — the carve began DURING this audit

At audit open the tree was clean for `orpheus/`. At audit close it was not. Two
changes landed under me:

| path | state |
|---|---|
| `orpheus/numerics/generating_measure.py` | **NEW, UNTRACKED** (`?? `), 531 lines — not at HEAD |
| `orpheus/numerics/measure.py` | **MODIFIED**, `+41 / −4`, uncommitted |

### What landed

`orpheus/numerics/generating_measure.py` (untracked, in-flight) is the
Golub-Welsch construction: `class GeneratingMeasure` (`:175`) with `.mass` (`:213`),
**`.gauss(n) -> DiscreteMeasure`** (`:224`) and `.on(a, b)` (`:279`); recurrence
coefficient generators `_legendre_recurrence` (`:340`), `_chebyshev_t_recurrence`
(`:350`), `_chebyshev_u_recurrence` (`:361`), `_hermite_recurrence` (`:368`);
family constructors `jacobi(a, b)` (`:409`) and `laguerre(a=0.0)` (`:491`).

`orpheus/numerics/measure.py` gained:
- `SPACE_HALF_LINE = "[0,inf)"` (`:108`)
- a `TYPE_CHECKING` forward-ref import of `GeneratingMeasure` (`:91`)
- a new `DiscreteMeasure` field **`generating_measure: "GeneratingMeasure | None" = None`** (`:239`)
- a rewritten `degree_of_exactness` docstring that now **states the claim is relative to `generating_measure`**.

### Effect on this audit — item by item

- **Targets 1–4 are ALL STILL PRESENT.** Line numbers for the three `measure.py`
  targets have **shifted by +37**: `gauss_legendre` 1064 → **1101**,
  `gauss_chebyshev` 1107 → **1144**, `equispaced` 1161 → **1198**.
  `gauss_legendre_on_mu` is unmoved at `rules_1d.py:33`. Every other `measure.py`
  line number in sections 1–3 above shifts by the same +37; `__all__` is now at
  `:1254-1266`. **All non-`measure.py` findings are unaffected.**
- **Section (e) is now HALF-CLOSED.** The semantic split I measured is *already
  acknowledged in the new docstring*, which cites the same phenomenon with an
  independently-measured number ("*misses the latter by 0.696 at n=4, q=x⁶*").
  The field is now documented as weight-relative.
- **But the (e) DEFECT I measured is STILL OPEN — and the carve did not touch it.**
  `__mul__` and `__add__` still do a bare `new_dx = min(...)` at
  **`measure.py:507` and `:568`** (was `:470`/`:531`), with **no reconciliation of
  `generating_measure`**. Verified: the diff contains no `new_dx`, no `__mul__`,
  no `__add__` hunk. So `GL4 + GC4` still returns `degree_of_exactness = 7` while
  integrating the constant 1 to 5.14. **Now that `generating_measure` exists as a
  field, the fix is available and cheap** — the composition rule can require
  matching generating measures (or drop the degree to `None`), which is precisely
  the guard the field was minted to make expressible. This is the highest-value
  still-open item in the audit.
- **Nothing consumes `generating_measure` yet.** Grep across `orpheus/` + `tests/`
  finds it only inside `generating_measure.py` itself and the `measure.py`
  docstring/field declaration. No rule sets it; no test exercises it; it is not
  re-exported from `orpheus/numerics/__init__.py`. The four targets have **not**
  been re-posed onto `GeneratingMeasure.gauss` — the retirement is still ahead,
  and this audit's blast radius stands.
- **`SPACE_HALF_LINE` + `laguerre()` signal scope beyond the brief's three
  families.** The construction already admits Jacobi (⊇ Legendre, Chebyshev-T,
  Chebyshev-U), Hermite and Laguerre. Note `orpheus/derivations/common/quadrature.py:500`
  ships an independent `gauss_laguerre` (oracle tree, `np.polynomial.laguerre.laggauss`
  at `:527`) — after the carve these become a **production/oracle pair**, which is
  the correct relationship, not a twin path, provided the oracle stays
  independent (it uses `laggauss`/mpmath, the new one uses Golub-Welsch — genuinely
  independent constructions).

### Confidence statement

All findings in sections 0–4 and (e)–(g) were taken against the tree as it stood
at HEAD `b0a003b4` with a clean `orpheus/`, and were re-verified against the final
tree except where noted above. The only invalidated content is the `measure.py`
line numbering (+37 below line 105) and the "already-dangling xref" line numbers
in (f) for `measure.py` (`1092` → `1129`, `1178` → `1215`). The **structural**
claims — zero production callers for targets 1–3, the graph-invisible registry
capture for target 4, the third `leggauss` in `rules_product.py`, the false doc
sample at `discrete_measures.rst:357`, the non-biting bit-identity test, and the
open `min()` unsoundness — are unchanged.

