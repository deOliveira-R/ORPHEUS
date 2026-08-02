# Issue #325 — verification plan (symmetry-exact circle node generation)

**Author:** test-architect · **Date:** 2026-08-01 · **Branch:** `refactor/operator-strategy-layers`
**SUT:** `orpheus/numerics/roots_of_unity.py` — `roots_of_unity(numerators, denominator) -> (cos, sin)`
**Status:** PRE-carve. Primitive exists on disk; NOT exported from `orpheus/numerics/__init__.py`;
the three consumers are NOT repointed.

Every `[M]` claim was measured in this session. Scripts:
`/Users/rodrigo/.claude/jobs/c30e4f25/tmp/{proto,real_mutants,m01…m17}.py`.

---

## 0. HEADLINE — read this before anything else

### 0.1 ⛔ #325 is BLOCKED by a live latent defect it makes visible

`[M]` **The per-level ordinate ordering in `rules_product.py:139` is under-determined by its
sort key, and the ambiguity is a ~1 % angular-flux error that does NOT converge away under
angular refinement.**

Four-leg end-to-end, cylindrical fixed-source, 2G, **heterogeneous** (3/4/3 cells, mixtures
A/B), random per-ordinate source, `inner_tol=1e-13`, `n_phi=32`:

```
A baseline (today: linspace+cos, argsort default)
B exact node values + LEGACY ordering forced back
C exact node values + default argsort ordering
D exact node values + kind="stable" ordering

max|A - B| = 1.054712e-14   bit-identical=False   rel=1.063e-14      <- node-VALUE drift alone
max|A - C| = 9.996633e-03   bit-identical=False   rel=1.008e-02
max|A - D| = 1.866051e-02   bit-identical=False   rel=1.881e-02
max|B - C| = 9.996633e-03   bit-identical=False   rel=1.008e-02      <- ORDERING alone
max|B - D| = 1.866051e-02   bit-identical=False   rel=1.881e-02      <- ORDERING alone
max|C - D| = 1.762316e-02   bit-identical=False   rel=1.775e-02      <- argsort KIND alone
```

`B` and `C` carry **bit-identical node values**; they differ only in `level_indices`. The
solver moves **1.0 %**. `C` and `D` differ only in `np.argsort(kind=)` — a numpy
implementation detail — and the solver moves **1.8 %**.

And it does not wash out:

```
 n_phi   max|C-D| (order ambiguity)          rel   max|A-B| (value drift)          rel
     8                 0.000000e+00    0.000e+00             1.221245e-15    1.205e-15
    16                 0.000000e+00    0.000e+00             3.108624e-15    3.088e-15
    32                 1.762316e-02    1.777e-02             1.054712e-14    1.063e-14
    64                 1.028758e-02    1.019e-02             7.438494e-15    7.365e-15
   128                 1.034785e-02    1.041e-02             1.831868e-14    1.843e-14
```

**The issue body's blast-radius sentence — "Values shift ≤1 ULP, so anything pinned
bit-exactly re-baselines" — is TRUE for the node values and FALSE for the solver results.**
The value drift is `~1e-14`; the ordering-induced drift is `~1e-2`, **twelve orders larger**.

Root cause, and it predates #325: the sort key `η = sinθ·cos φ` is **2-to-1** over
`φ ∈ [0, 2π)` — `φ` and `2π−φ` share `η`, differing only in `ξ = μ_y`. Today those two floats
differ in the last bits so the sort is accidentally total; #325 makes them bit-equal, which is
the whole point of the change, and the accident evaporates. `[M]` distinct-η census per level:

```
 n_phi  #distinct eta OLD  #distinct eta NEW  #tied pairs NEW  #tied pairs OLD
     3                  3                  2                1                0
     4                  4                  3                1                0
     6                  6                  4                2                0
     8                  8                  5                3                0
    16                 16                  9                7                0
    20                 19                 11                9                1     <- ties ALREADY today
    32                 30                 17               15                2     <- ties ALREADY today
    64                 60                 33               31                4     <- ties ALREADY today
```

`[M]` The per-level ordering changes for **every** `(n_mu, n_phi)` in
`{2,4,6,8} × {3,4,5,6,8,12,16,24,32}` — 36 of 36 — including `n_phi ∈ {3,5,6}` which have no
axis node at all (there the tie is the *mirror pair*, not the axis).

⇒ **Deliverable (7) — the ordering ruling — must land BEFORE the node-value change, as its own
commit with its own reference.** §7.

### 0.2 The other four headline findings

| # | Finding | Where |
|---|---|---|
| **F1** | `[M]` The acceptance criterion **"agreement with `np.cos`/`np.sin` to ≤1 ULP" is FALSE as written, and epistemically backwards.** Max \|shipped − `np.cos(2πp/q)`\| = **8.33e-16 = 3.75 ulp(1.0)**. `np.cos(2πp/q)` is not the true value. Against 100-digit mpmath: shipped ≤ **1.257e-16 = 0.57 ulp(1.0)**; legacy ≤ **8.256e-16 = 3.72 ulp(1.0)**. | §1 L6 |
| **F2** | `[M]` The gate named in `angular_trace_space.py` as `test_eps_below_min_genuine_cosine` **does not exist** (dangling Python-domain xref at lines 90 and 159 — silent, no `-W` warning). The real gate is `test_eps_sits_in_the_round_off_to_genuine_gap`, and **`product` is NOT in its fixture list and would FAIL today.** | §3 |
| **F3** | `[M]` The nearest-neighbour reflection search is **NOT mis-pairing today** for even `n_phi` — perm ∧ involution hold up to `n_phi=1024`, margin `5.0e-3` vs a `1e-16` perturbation. The live-bug hypothesis is **REFUTED**. But for **odd `n_phi` the `axis=0` map is garbage** (residual up to **0.94**) while still being a valid permutation AND a valid involution. | §4 |
| **F4** | `[M]` `tests/numerics/test_rules_product.py::test_product_bit_identical_to_legacy_adapter` — whose docstring claims *"Pins the cylindrical regression snapshots (`cyl_*_product_*.npz`)"* — is a **TAUTOLOGY**. `Quadrature.product` calls `product_mu_phi`; the named "legacy adapter" `orpheus.sn.quadrature.ProductQuadrature` **does not exist**. Proven: replacing `product_mu_phi` with random nodes leaves it **green**. | §5.4 |
| **F5** | `[M]` **Unit norm is designed-green on 11 of 11 mutants**, including the constant map `(1,0)`. **The x-mirror law is also satisfied by the constant map.** Both of the user's suspicions **CONFIRMED by measurement**. | §2 |

---

## 1. Deliverable (1) — the primitive's intrinsic-property test

**Home:** `tests/numerics/test_roots_of_unity.py`.
**Level:** `@pytest.mark.foundation` — these are software/mathematical invariants of a
numerics primitive with no theory-page `:label:`; per `feedback_vv_tagging` a foundation
module carries **no** `verifies(...)`. (If the archivist later mints a
`:label: roots-of-unity-quadrant` theory anchor — the module already declares
`:label: roots-of-unity-quarter-split` and `:label: roots-of-unity-quadrant` in its own
docstring — the L1/L4 rows may additionally carry `verifies(...)`. That is a docs decision,
not a test decision.)

### 1.1 The thirteen laws

`[M]` All thirteen hold on the shipped body over
`q ∈ {1,2,3,4,5,6,7,8,9,10,12,16,20,24,32,40,48,64,96,128}`.

| ID | Law | Exact statement | Comparator |
|----|-----|-----------------|------------|
| **L1** | on-axis exact | `4p ≡ 0 (mod q)` ⟹ `(cos,sin) ∈ {(1,0),(0,1),(−1,0),(0,−1)}` | `np.array_equal` |
| **L2** | x-mirror, **as an index permutation** | `cos == cos[(q−k)%q]` ∧ `sin == −sin[(q−k)%q]` | `np.array_equal` |
| **L3** | y-mirror, **as an index permutation** | `cos == −cos[(q/2−k)%q]` ∧ `sin == sin[…]` (even `q`) | `np.array_equal` |
| **L4** | 45° diagonal, **VALUE-PINNED** | `2r == q` ⟹ `\|cos\| == \|sin\| == np.sqrt(0.5)` | `np.array_equal` |
| **L5** | unit norm | `\|cos²+sin²−1\| ≤ 2 ulp` | `np.testing` |
| **L6** | **≤1 ulp of TRUTH** | `\|component − mpmath(dps=50)\| ≤ np.spacing(1.0)` | `np.testing` |
| **L7** | dihedral closure | node SET invariant under `p→p+1` and `p→−p` | set equality on floats |
| **L8** | injectivity | exactly `q` distinct `(cos,sin)` pairs | `len(set(...)) == q` |
| **L9** | rotation is one `q`-cycle | `p→p+1` permutes the set as a single orbit | orbit walk |
| **L10** | **no negative zero** | `not any(signbit(x) & (x == 0))` | `np.signbit` |
| **L11** | periodicity mod `q` | `f(p+q) == f(p) == f(p−3q)` | `np.array_equal` |
| **L12** | `ζ⁰ = (1,0)` | `f(0,q) == (1.0, 0.0)` | `==` |
| **L13** | addition theorem | `ζ^p·ζ¹ == ζ^(p+1)` (complex product) ≤ 8 ulp | `np.testing` |

Three **guard** rows (`[M]` all confirmed raising):

| ID | Input | Behaviour |
|----|-------|-----------|
| **G1** | `np.array([0.0, 1.0])`, or a Python `float` | `ValueError: … requires integer numerators …` |
| **G2** | `denominator ∈ {0, −4}` | `ValueError: … requires denominator >= 1 …` |
| **G3** | `np.array([True, False])` | `ValueError` — **`np.issubdtype(np.bool_, np.integer)` is `False`.** Pin it whichever way the user rules; today it refuses. |

One **shape** row (`[M]` holds for scalar / 0-d / 1-d / 2-d):

```
  python int 3   in.shape=()       out.shape=()       MATCH=True
  np.int64(3)    in.shape=()       out.shape=()       MATCH=True
  0-d array      in.shape=()       out.shape=()       MATCH=True
  1-d (4,)       in.shape=(4,)     out.shape=(4,)     MATCH=True
  2-d (2,3)      in.shape=(2, 3)   out.shape=(2, 3)   MATCH=True
```

### 1.2 The `q` sweep — every member is load-bearing

| class | members | ACTIVATES | NULLS |
|---|---|---|---|
| `4 ∣ q`, `8 ∤ q` | 4, 12, 20, 28 | axis (L1) + octant swap | the 45° branch |
| `8 ∣ q` | 8, 16, 24, 32, 40, 48, 64, 96, 128 | **the 45° fixed point** (L4) | — |
| odd `q` | 3, 5, 7, 9 | generic octant fold | axis beyond `p=0`; L3 vacuous |
| `q ≡ 2 (mod 4)` | 6, 10 | y-mirror with no x-axis node | 45°; x-axis |
| degenerate | **1, 2** | `r == 0` for every `p`; the `quad==2` select lane | fold and swap entirely |

**`q ∈ {1,2}` is mandatory.** `[M]` `q=2` is the only sweep member that enters the `quad==2`
lane where `−sin_quad` produces the `−0.0` that L10 pins. A sweep starting at `q=3` never
reaches it. This is the config-blindness discipline applied to a pure-numerics primitive: the
convenient sweep (`q ≥ 3`, "the interesting cases") nulls the exact branch L10 exists to test.

### 1.3 Per-law mutation matrix — which edit reds which law

`[M]` Mutants are single-edit perturbations of the **shipped body** (`real_mutants.py`).
`.` = law REDS (has teeth). `T` = law HOLDS (blind).

```
                               M0    M1    M2    M3    M4    M5    M6    M7    M8    M9   M10   M11
L1  on-axis exact               T     T     .     T     .     .     T     T     T     .     .     .
L2  x-mirror (index perm)       T     .     .     .     .     T     .     T     .     .     .     T
L3  y-mirror (index perm)       T     .     .     .     .     T     .     T     .     .     .     .
L4  45deg diagonal              T     .     T     T     T     T     T     T     .     T     .     .
L5  unit norm                   T     T     T     T     T     T     T     T     T     T     T     T
L6  <=1 ulp of TRUTH            T     T     .     .     .     .     T     T     T     T     .     .
L7  dihedral closure            T     .     .     .     .     T     .     T     .     .     .     T
L8  q distinct nodes            T     T     T     .     T     T     T     T     T     T     T     .
L9  rotation = q-cycle          T     T     T     .     T     T     T     T     T     T     T     .
L10 no negative zero            T     T     T     T     T     T     T     .     T     T     T     T
L11 periodic mod q              T     T     T     T     T     T     T     T     T     T     .     T
L12 zeta^0 == (1,0)             T     T     .     T     .     T     T     T     T     .     T     T
L13 addition theorem            T     T     .     .     .     T     T     T     T     T     T     T

  M0  control (shipped)        M6  no octant fold
  M1  no 45deg patch           M7  no signed-0 canonicalisation
  M2  octant swap inverted     M8  diag from cos/sin(pi/4) instead of sqrt(0.5)
  M3  quad-2 sin sign flip     M9  theta perturbed (kills the axis)
  M4  global cos/sin swap      M10 legacy linspace+cos
  M5  conjugate p -> -p        M11 constant map (1, 0)
```

Committed catcher per mutation:

| Mutation | Caught by | Sole catcher |
|---|---|---|
| M1 drop the 45° patch | L2, L3, L4, L7 | — |
| M2 invert the octant swap | L1, L2, L3, L6, L7, L12, L13 | — |
| M3 quadrant sign flip | L2, L3, L6, L8, L9, L13 | — |
| M4 global `cos`↔`sin` swap | L1, L2, L3, L6, L7, L12, L13 | — |
| M5 conjugate `p→−p` | **L1, L6 only** | near-sole |
| M6 drop the octant fold | L2, L3, L7 | — |
| M7 drop `+ 0.0` canonicalisation | **L10** | **YES** |
| M8 diag from `cos/sin(π/4)` | L2, L3, L4, L7 | — |
| M9 perturb `theta` | L1, L2, L3, L7, L12 | — |
| M10 legacy `linspace`+`cos` | L1, L2, L3, L4, L6, L7, L11 | — |
| M11 constant map `(1,0)` | L1, L3, L4, L6, L8, L9 | — |

`{L2, L6, L10}` is a size-3 covering set. **Ship all thirteen anyway** — the project standard
is one test per defining law, and the matrix is the evidence for which are load-bearing.

### 1.4 The gate module is WRITTEN and every gate is PROVEN able to red

`tests/numerics/test_roots_of_unity.py` — **251 passed** under both `python -O -m pytest` and
plain `pytest`, 0.5 s, `@pytest.mark.foundation`, zero bare `assert` (all `np.testing` /
`pytest.fail` / `pytest.raises`).

**Method note (vv Mode 8 METHOD WARNING).** My first teeth harness reported *"18 gates never
red"*. That was an artefact: it patched a privately-loaded copy of the test module while
`pytest.main` imported its own, so **the mutation never bound**. The corrected harness patches
`orpheus.numerics.roots_of_unity.roots_of_unity`, purges the test module from `sys.modules`,
and carries a **positive-control probe** reporting whether the freshly-imported test module
actually bound the mutant. Every row below reads `[BIT]`. Do not trust a "gate is blind"
verdict without that probe.

`[M]` Control clean, all eleven mutations bound:

```
CONTROL (unmutated): 17 green, 0 red -> none
  positive-control probe: test module bound a mutant? False  (expect False)
  M1 no 45deg patch                reds  6 gates   [BIT]
  M2 octant swap inverted          reds  9 gates   [BIT]
  M3 quad-2 sin sign               reds  9 gates   [BIT]
  M4 global cos/sin swap           reds  7 gates   [BIT]
  M5 conjugate p->-p               reds  3 gates   [BIT]
  M6 no octant fold                reds  5 gates   [BIT]
  M7 no signed-0 canon             reds  3 gates   [BIT]
  M8 diag = cos/sin(pi/4)          reds  6 gates   [BIT]
  M9 theta perturbed               reds  7 gates   [BIT]
  M10 legacy linspace+cos          reds  9 gates   [BIT]
  M11 constant (1,0)               reds  8 gates   [BIT]
```

`[M]` Per-gate teeth:

```
                                                                 M1   M2   M3   M4   M5   M6   M7   M8   M9  M10  M11
test_addition_theorem                                             .    R    R    R    .    .    .    .    .    .    .
test_agrees_with_arbitrary_precision_to_one_ulp                   .    R    R    R    R    .    .    .    .    R    R
test_diagonal_components_are_the_correctly_rounded_sqrt_half      R    .    .    .    .    .    .    R    .    R    R
test_no_negative_zero                                             .    .    .    .    .    .    R    .    .    .    .
test_node_set_is_closed_under_the_dihedral_action                 R    R    R    R    .    R    .    R    R    R    .
test_nodes_are_distinct                                           .    .    R    .    .    .    .    .    .    .    R
test_non_integer_numerators_are_refused                           R    R    R    .    R    R    R    R    R    R    R
test_non_positive_denominator_is_refused                          R    R    R    .    .    R    R    R    R    R    R
test_on_axis_nodes_are_exactly_zero_and_unit                      .    R    .    R    R    .    .    .    R    R    R
test_output_shape_matches_input_shape                             .    .    .    .    .    .    .    .    .    .    .
test_periodic_in_the_numerator_mod_q                              .    .    .    .    .    .    .    .    .    R    .
test_rotation_acts_as_a_single_q_cycle                            .    .    R    .    .    .    .    .    .    .    R
test_sqrt_half_constant_is_correctly_rounded                      .    .    .    .    .    .    .    .    .    .    .
test_unit_norm                                                    .    .    .    .    .    .    .    .    .    .    .
test_x_mirror_is_an_exact_index_permutation                       R    R    R    R    .    R    .    R    R    R    .
test_y_mirror_is_an_exact_index_permutation                       R    R    R    R    .    R    .    R    R    R    R
test_zeroth_root_is_exactly_one                                   .    R    .    R    .    .    .    .    R    .    .
```

The three all-`.` rows have threat models the body mutants do not touch. `[M]` targeted
mutations prove each one reddens:

```
CONTROL
  unmutated                                  reds: NONE

--- T1: gross magnitude error in the diagonal constant (unit-norm's tooth) ---
  _SQRT_HALF = 0.707 (3 s.f.)                reds: [... 'test_unit_norm']
  _SQRT_HALF = 1/sqrt(2)  (-1 ULP)           reds: [...]   <- unit_norm does NOT red: its documented blind spot

--- T2: respell the module CONSTANT (sqrt-half gate's tooth) ---
  _SQRT_HALF := np.sin(pi/4)  (-1 ULP)       reds: ['test_diagonal...', 'test_sqrt_half_constant_is_correctly_rounded']
  _SQRT_HALF := 1/np.sqrt(2)  (-1 ULP)       reds: ['test_diagonal...', 'test_sqrt_half_constant_is_correctly_rounded']

--- T3: shape promotion (shape gate's tooth) ---
  np.atleast_1d on the return (0-d -> (1,))  reds: ['test_output_shape_matches_input_shape']
  .ravel() on the return (2-d -> 1-d)        reds: ['test_output_shape_matches_input_shape']
```

⇒ **17 of 17 gates are provably able to red.**

### 1.5 A tautology the harness caught in my OWN gate — recorded, not hidden

The first draft carried `test_y_mirror_is_not_claimed_for_odd_q`, written as an
activation/companion guard asserting that odd `q` has no node at `π`
(`not np.any(2 * p == q)`). The teeth harness flagged it as never-reddening, and `[M]` it is a
**theorem about parity**, not a property of the code:

```
for odd q, can 2*p == q ever hold for integer p?
  q=3: any(2p == q) = False  (2p is EVEN, q is ODD)
  q=101: any(2p == q) = False  (2p is EVEN, q is ODD)
=> the predicate is a THEOREM about parity, not a property of the code.
=> it executes under every runtime mode and can NEVER red: vv Mode 8 class 2
```

Removed, with a comment in its place explaining why. This is the exact class
`vv-principles` Mode 8 names — *"audit COMPANION/activation guards for reddenability: ask
what input makes this assertion fail; if no input can, the guard is dead."* Worth noting that
it survived writing and a full green run, and only the mutation harness exposed it.

---

## 2. Deliverable (2) — Mode-12 stabiliser analysis (DESIGN TIME, pre-mutation)

For each gate: the measured functional, its invariance group, and the intersection with the
threat model. **Both of the user's suspicions are CONFIRMED, and I found four more.**

### 2.1 CONFIRMED — (a) unit norm cannot see a transposed octant

- **Functional:** `N(c,s) = c² + s² − 1`.
- **Invariance group:** `O(2)` acting on the pair — every rotation and reflection of
  `(cos, sin)`, hence the swap `c↔s`, every sign flip, conjugation, and **any relabelling of
  which root goes where**.
- **`[M]` Verdict: all 11 mutants pass L5, including M11 the constant map `(1,0)`.**
- **NOT a dead gate — its honest scope.** `[M]` L5 *does* red on a magnitude error:

```
mutant                                        unit norm   x-mirror
M0 control (shipped)                               True       True
M12 _SQRT_HALF off by 1 ULP                        True       True
M13 _SQRT_HALF = 0.707 (3sf)                      False       True     <- L5's only tooth
M14 _SQRT_HALF = 1/sqrt(2) (-1 ulp)                True       True
M15 global scale (1+2e-16)                         True       True
M16 cos<->sin swap (control)                       True      False
```

  ⇒ **Ship L5, but NEVER credit it against any symmetry/orientation/permutation mutation.**
  Its docstring must say: *"invariance group O(2); catches gross magnitude error only
  (`[M]` reds at `_SQRT_HALF = 0.707`, green at ±1 ULP)."*

### 2.2 CONFIRMED — (b) a mirror law is satisfied by a degenerate map

- **Functional:** `cos == cos[mirror]` ∧ `sin == −sin[mirror]`.
- **Invariance group:** contains (i) the **constant map** — `[M]` M11 `(1,0)` passes L2 and
  L7; (ii) **global conjugation** `p → −p` — `[M]` M5 passes L2, L3 **and** L7, because
  conjugation is an automorphism of the dihedral group; (iii) **signed zero**, since
  `-0.0 == 0.0` — `[M]` M7 passes L2 and L3.
- This is structurally the `face_opposite` involution/no-fixed-point pair from B3.4c: an
  involution law alone is blind because **the identity is an involution**.
- **Closure (three legs, all required):**
  1. **L8 injectivity** — `q` distinct nodes. Kills the constant map. (`[M]` L8 reds on M11.)
  2. **L1 / L6 pointwise anchor** — the ONLY catchers of M5 conjugation. A mirror law can
     never see a global orientation flip; only a pointwise reference can.
  3. **L10 signbit** — the ONLY catcher of M7.

### 2.3 NEW — the `|cos| == |sin|` form of the 45° law is blind to a wrong constant

`[M]`:

```
  sqrt(0.5) [correct]    |cos|==|sin| law: True   value-pinned law: True
  1/sqrt(2) [-1ulp]      |cos|==|sin| law: True   value-pinned law: False
  nextafter up [+1ulp]   |cos|==|sin| law: True   value-pinned law: False
```

- **Invariance group of `|cos| − |sin|` at the diagonal:** any common value `v` — a
  one-parameter family. `[M]` `np.cos(np.pi/4)` **is** bit-equal to `np.sqrt(0.5)`; it is
  `np.sin(np.pi/4)` and `1/np.sqrt(2)` that are 1 ULP low, and `np.sqrt(2)/2` that is equal.
- ⇒ **L4 MUST pin the value `np.sqrt(0.5)`, not merely the equality.** Add a companion
  `test_sqrt_half_constant_is_correctly_rounded` asserting `_SQRT_HALF == np.sqrt(0.5)` with
  the four alternative spellings enumerated in its body — that pins the choice against a
  future "equivalent" respelling.

### 2.4 NEW — the ERR-044 involution certification is blind to a WRONG partner map

`[M]` The odd-`n_phi` `axis=0` reflection map has residual **0.94** and is nevertheless a
**valid permutation** and a **valid involution**:

```
  product(4,3)  max|Rx - x_ref| = 9.404e-01
      ERR-042 measure-preserving: RAISE(BoundaryGeometryMapNotMeasurePreservingError)
      ERR-044 involutive: PASS
      ERR-045 inflow->outflow: RAISE(ReflectionDidNotMapInflowToOutflowError)
```

- **Invariance group of "is an involution":** the whole centraliser of order-2 elements — any
  self-inverse pairing, however geometrically wrong.
- ⇒ **the only functional outside the stabiliser is the RESIDUAL** `max‖R·x_n − x_{ref[n]}‖`.
  That IS deliverable (4)'s gate. §4.

### 2.5 NEW — "sortedness by η" is blind to tie permutation (deliverable 7's trap)

- **Functional:** `np.all(np.diff(eta[level_indices]) >= 0)`.
- **Invariance group:** the **symmetric group on each block of equal η**. `[M]` legs C and D
  differ by exactly such a permutation and BOTH are sorted; they differ by **1.8 %** in the
  converged flux.
- ⇒ **A "the level is sorted by η" assertion is exactly the wrong gate.** The functional
  outside the stabiliser is the **full index tuple** `tuple(level_indices[p])` compared against
  an independently-constructed expected tuple. §7.

### 2.6 NEW — the canonical curvilinear diagnostic is blind to the ordering

`[M]` The flat-flux fixed-source `Q/Σ_t` leg (homogeneous, uniform source, reflective,
`n_phi=4`) gives `phi.ravel()[:3] = [0.6631456 0.6631456 0.6631456]` on **all four** legs, and

```
max|A - B| = 1.332268e-15   max|B - C| = 1.332268e-15   max|C - D| = 0.000000e+00
```

The single most powerful curvilinear diagnostic in the playbook **cannot discriminate the
orderings** — the homogeneous flat-flux config nulls exactly the term the ordering perturbs.
This is `AGENT.md` §0.6 config-blindness in its sharpest form: the *recommended minimum catch*
is here a designed-green gate. **The heterogeneous anisotropic-source config at `n_phi ≥ 32` is
mandatory** for any ordering claim.

### 2.7 NEW — `np.array_equal` is blind to signed zero across the whole suite

`-0.0 == 0.0` is `True`, so every `array_equal` / `allclose` / `assert_array_equal` comparison
in the tree is blind to the `+ 0.0` canonicalisation. `[M]` M7 passes L1, L2, L3, L4, L7, L11,
L12. Only `np.signbit` sees it. Two consequences:
- L10 is the **sole** catcher and must be written on `np.signbit`.
- **A byte-level snapshot hash WOULD see it** (different bit pattern) while `array_equal`
  would not — so a `.npz` regenerated with and without the canonicalisation compares equal by
  value and unequal by `sha256`. Note this in the re-baseline record (§5).

### 2.8 Summary table — invariance group ∩ threat model

| Gate | Measured functional | Invariance group | Designed-green on |
|---|---|---|---|
| L5 unit norm | `c²+s²−1` | `O(2)` on the pair | **all 11 mutants** |
| L2/L3 mirror | mirror index permutation | constant map; global conjugation; signed zero | M5, M7, M11 |
| L7 closure | node SET under `⟨rot, refl⟩` | same as L2/L3 | M5, M7, M11 |
| L4 as `\|c\|==\|s\|` | equality of magnitudes | any common value `v` | wrong `_SQRT_HALF` |
| L4 **value-pinned** | `\|c\| == sqrt(0.5)` | — | M2, M3, M5 |
| ERR-044 involution | `ref[ref] == id` | any self-inverse pairing | **residual-0.94 garbage map** |
| "sorted by η" | `diff(eta) >= 0` | `S_k` on each equal-η block | **the 1.8 % tie permutation** |
| flat-flux `Q/Σ_t` | homogeneous converged flux | the whole ordering group | **the ordering entirely** |
| any `array_equal` | value equality | signed zero | M7 |

---

## 3. Deliverable (3) — the `TANGENTIAL_EPS` demotion gate

### 3.1 What exists today, stated plainly

- `[M]` **The gate named in the production docstring does not exist.**
  `orpheus/numerics/spaces/angular_trace_space.py:90` and `:159` both cite
  `` :func:`test_eps_below_min_genuine_cosine` ``. `grep -rn` across `orpheus/ docs/ tests/`
  returns **only those two production lines** — no such test. The real gate is
  `tests/numerics/test_angular_trace_space.py:325::test_eps_sits_in_the_round_off_to_genuine_gap`.
  It is a Python-domain xref, so Sphinx `-W` never warned (`coding-standards.md`: the silent
  retirement class). **Fix the two docstring lines in the same commit as the gate.**
- `[M]` **`product` is NOT in the fixture list** — it is
  `gauss_legendre(2,3,4,5,7,8,16,32,64)` + `lebedev(3,5,7,11,17,29,53)`. Neither `product`
  nor `level_symmetric` appears.
- `[M]` **The module docstring's claim "nominally-tangential cosines are exactly `0.0`" is
  FALSE for `product` today**, and the existing gate would FAIL if `product` were added:

```
--- families ALREADY in the gate's fixture list ---
PASS  gauss_legendre(8)        mu_x: exact0=  0 in-band=  0 max_band=0.000e+00 …
PASS  lebedev(17)              mu_x: exact0= 12 in-band=  0 max_band=0.000e+00 …

--- families NOT in the gate's fixture list ---
FAIL  product(n_mu=2, n_phi=4)   mu_x: exact0=  0 in-band=  4 max_band=1.500e-16  mu_y: exact0=  2 in-band=  2 max_band=9.999e-17
FAIL  product(n_mu=4, n_phi=8)   mu_x: exact0=  0 in-band=  8 max_band=1.728e-16  mu_y: exact0=  4 in-band=  4 max_band=1.152e-16
PASS  product(n_mu=2, n_phi=3)   mu_x: exact0=  0 in-band=  0 max_band=0.000e+00  mu_y: exact0=  2 in-band=  0 max_band=0.000e+00
FAIL  product(n_mu=4, n_phi=6)   mu_x: exact0=  0 in-band=  0 max_band=0.000e+00  mu_y: exact0=  4 in-band=  4 max_band=1.152e-16
FAIL  product(n_mu=8, n_phi=32)  mu_x: exact0=  0 in-band= 16 max_band=1.806e-16  mu_y: exact0=  8 in-band=  8 max_band=1.204e-16
PASS  level_symmetric(4)        mu_x: exact0=  0 in-band=  0 max_band=0.000e+00 …

existing fixture list passes : True
product would pass           : False
level_symmetric would pass   : True
```

  Refinement the issue body does not carry: `product(2,3)` **passes**, because with odd
  `n_phi` no angle is a quadrant multiple so no cosine is nominally zero at all. The violation
  is structured: `μ_x` violates iff `4 ∣ n_phi`; `μ_y` violates iff `n_phi` is even (the
  `φ = π` node, `sin(π) = 1.22e-16`). `[M]` `μ_y` has exactly `n_mu` exact zeros in every case
  — the `φ = 0` nodes, where `sin(0)` **is** exactly `0.0`. That asymmetry is the issue body's
  "not stable under a face flip" claim, and it is confirmed.

- **TANGENTIAL_EPS is genuinely load-bearing today.** `[M]` every in-band value is
  `≤ 1.81e-16 < TANGENTIAL_EPS = 8.881784197001252e-16`, so those ordinates *are* currently
  classified tangential. Removing the eps today would break the partition. The demotion is
  earned only after §6 step 4.

### 3.2 The gate — a REGISTRY WALK, not a fixture list

The enumeration surface is `orpheus/numerics/quadrature/registry.py:434`
`quadrature_registry: list[QuadratureSpec]` — four specs today (`GaussLegendre1D`,
`LebedevSphere`, `LevelSymmetricSN`, `ProductQuadrature`), each carrying `factory`,
`parameters`, and `degree_of_exactness_for` (a target-degree → parameter-dict inverter).

**Gate T1 — `test_every_registered_quadrature_has_exactly_zero_tangential_cosines`**
(`tests/numerics/test_angular_trace_space.py`, `@pytest.mark.foundation`):

1. `for spec in quadrature_registry:` — **walk it, never a hand-listed table.** This is
   directly the C5.3/#225 lesson: a hand-listed 4-face transcription silently lacked the z
   faces.
2. For each spec, invert `spec.degree_of_exactness_for(d)` over
   `d ∈ {1,3,5,7,9,11,15,19}`, drop `None`, dedupe the parameter dicts, and build the measure
   via `spec.factory(**params)`.
3. For each of the three cosine columns, assert
   `not np.any((|μ| > 0.0) & (|μ| <= TANGENTIAL_EPS))` via `np.testing`/explicit `raise`
   (never a bare `assert` — this file is collected so the rewriter protects it, but the helper
   that walks the registry may not be; §6 checks that).
4. **Coverage self-check (the anti-vacuity leg):** assert
   `len(specs_walked) == len(quadrature_registry)` and `len(param_sets) >= 2` per spec, and
   assert the walked name set **equals** `{s.name for s in quadrature_registry}`. A future
   family added to the registry is then automatically in the gate; a future family added
   *outside* the registry is a separate (documented) hole — see T3.

**Gate T2 — the demotion statement itself.**
`test_tangential_eps_is_provably_inert`: for every registry-walked measure, assert the
inflow/outflow/tangential partition computed with `TANGENTIAL_EPS` is **bit-identical** to the
partition computed with `eps = 0.0`. That is the operational meaning of "demoted from
classifier to defensive guard", and it is the sentence the docstring should then make.

**Gate T3 — the enumeration hole, made explicit.** `Quadrature` also has factory
classmethods (`gauss_legendre`, `lebedev`, `level_symmetric`, `product`). Assert the set of
public `Quadrature` classmethod factories is a subset of the registry's `factory` targets —
so a fifth factory added to `Quadrature` without a registry spec **reds T3**. Without this,
"registry walk" is only as complete as the registry.

**Mutation proof for T1/T2:** monkeypatch `roots_of_unity` back to `linspace`+`cos` inside
`rules_product` and confirm T1 and T2 both red on the `ProductQuadrature` spec. `[M]` the
predicate is already known to fail there (§3.1), so this reddening is guaranteed — but run it
so the record shows the gate fired, not the prediction.

### 3.3 One more flag

`registry.py:473` tags `ProductQuadrature` with `axis_aligned=False`, commented
*"depends on phi grid; conservatively False"*. After #325 the rule **is** exactly axis-aligned
whenever `4 ∣ n_phi`. Either tighten the flag (making it `n_phi`-dependent) or leave it and
record why. It is a structural-flag claim that the selector reads; leaving a now-false
conservative tag is the kind of decayed claim §2 is about.

---

## 4. Deliverable (4) — the `reflection_index` exactness gate

### 4.1 Measured state of the search TODAY — the live-bug hypothesis is REFUTED

`[M]` For even `n_phi` the NN search is a valid permutation and a valid involution at every
configuration measured, and the mis-pair margin is enormous:

```
=== stress: LARGE n_phi (near-degenerate spacing) ===
  n_phi=   64: min node separation=8.013e-02  non-perm/non-invol per axis = {0:(F,F),1:(F,F),2:(F,F)}
  n_phi=  128: min node separation=4.008e-02  non-perm/non-invol per axis = {0:(F,F),1:(F,F),2:(F,F)}
  n_phi=  256: min node separation=2.004e-02  non-perm/non-invol per axis = {0:(F,F),1:(F,F),2:(F,F)}
  n_phi=  512: min node separation=1.002e-02  non-perm/non-invol per axis = {0:(F,F),1:(F,F),2:(F,F)}
  n_phi= 1024: min node separation=5.010e-03  non-perm/non-invol per axis = {0:(F,F),1:(F,F),2:(F,F)}
```

Node separation shrinks as `2π/n_phi`; the perturbation is `~1e-16`. A mis-pair needs
`2π/n_phi ~ 1e-16`, i.e. `n_phi ~ 6e15`. **Not reachable. Report this to the user as a
refutation:** the NN search on `product` is an inelegance and an `O(N²)` cost
(`directional.py:9` concedes it), **not** a live mis-pairing bug.

`[M]` But the match is *not exact* today, and becomes exact after the change:

```
  product(4,8)  OLD  axis=0  exact=False  max|Rx - x_ref|=3.455e-16
  product(4,8)  NEW  axis=0  exact=True   max|Rx - x_ref|=0.000e+00
  product(6,16) NEW  axis=0/1/2  exact=True  max|Rx - x_ref|=0.000e+00
  lebedev(17)   axis=0/1/2  exact=True   max|Rx - x_ref|=0.000e+00   (control: already algebraic)
  level_symmetric(8) axis=0/1/2 exact=True max|Rx - x_ref|=0.000e+00  (control)
```

### 4.2 The genuine live defect — odd `n_phi`

`[M]` For **odd** `n_phi` the `axis=0` map is garbage and `reflection_index(0)` returns it
without complaint:

```
          product(4,3)     0   perm=True invol=True   max|Rx-x_ref|=9.404e-01  mean=6.597e-01
          product(4,5)     0   perm=True invol=True   max|Rx-x_ref|=5.812e-01  mean=4.047e-01
          product(4,7)     0   perm=True invol=True   max|Rx-x_ref|=4.185e-01  mean=2.909e-01
          product(4,9)     0   perm=True invol=True   max|Rx-x_ref|=3.266e-01  mean=2.268e-01

product(4,5).reflection_index(0) returned WITHOUT error: [ 3  2  1  0  4  8  7  6  5  9 13 12 11 10 14 18 17 16 15 19]
```

The azimuthal set `{2πk/n}` for odd `n` is **not invariant** under `φ → π−φ`; there is no
x-mirror partner and `argmin` returns the nearest node however far. `[M]` **ERR-042 and
ERR-045 DO raise** on it (so a reflective/albedo BC construction is protected), **ERR-044 does
not** (§2.4). `[M]` odd `n_phi` is essentially unreachable in production — one grep hit,
`tests/numerics/test_registry.py:165`, and that is a comment on `_product_invert(4)` returning
`{"n_mu": 3, "n_phi": 5}`. **So the selector CAN mint an odd-`n_phi` product rule.** That is a
latent landmine independent of #325.

### 4.3 Recommendation: **KEEP the search, ADD the exactness gate.** Do NOT replace it.

Reasons, in priority order:

1. **The search is the only thing that works for `lebedev`.** Lebedev orbits are not indexed
   by an azimuthal integer; there is no closed-form partner index. Replacing the search with a
   `product`-specific index formula would mint a **twin path** (coding-elegance anti-pattern
   #1): one partner-map derivation for `product`, another for `lebedev`, both claiming to
   compute the same object. That is exactly the architecture Cardinal Rule 2 forbids.
2. **Exactness is a property of the INPUT, not of the algorithm.** With an exactly-symmetric
   node set the NN search *is* an exact lookup — `[M]` residual `0.000e+00` on all three axes.
   The right move is to assert the property the input now has, not to change the algorithm
   that consumes it.
3. **Replacing it would delete the only detector of the odd-`n_phi` defect.** A hard-coded
   index formula would return a "partner" for odd `n_phi` with no residual to inspect.

**Gate R1 — `test_reflection_partner_match_is_bit_exact`** (`tests/numerics/`,
`@pytest.mark.foundation`): registry-walk (as §3.2), and for every axis assert
`np.array_equal(R @ nodes, nodes[ref])` — **`array_equal`, residual exactly zero**, not a
tolerance. **Skip/xfail the `axis=0` row for odd `n_phi`** with an explicit
`reason="no x-mirror exists for odd n_phi — see R2"`, `strict=False`.

**Gate R2 — `test_reflection_partner_map_refuses_a_nonexistent_mirror`**: a NEW production
guard inside `_find_reflections` (or `_compute_sphere_reflection_partners`) that **raises**
when `max‖R·x_n − x_{ref[n]}‖` exceeds a structural threshold. Threshold: `0.0` is correct
after #325 for every shipped family (`[M]` all three algebraic families give exactly `0.0`),
so the guard is `if resid != 0.0: raise`. That converts §4.2 from a silent landmine into a
construction-time error, and it is the same "make it unspellable rather than guarded" move the
issue is about. Negative test: `Quadrature.product(4, 5)` must raise.

**Gate R3 — mutation proof:** monkeypatch `roots_of_unity` back to `linspace`+`cos` and
confirm R1 reds with residual `~3.5e-16`. `[M]` that value is already measured.

**Do NOT** add an "is an involution" or "is a permutation" row as the exactness evidence —
`[M]` §2.4 proves both pass on a residual-`0.94` garbage map. They are legitimate as
*separate* ERR-044 regression rows; they are **not** catchers for this claim.

---

## 5. Deliverable (5) — the re-baseline discipline

### 5.1 The three `vv-principles` criteria, concretely

**Criterion 1 — the new formulation is principled at every step.**
Every intermediate is named and integer-exact until the single transcendental call:
`quad, r = divmod(4p, q)` (quarter turn + residual), `lo = min(r, q−r)` (octant fold),
`theta` (first-octant angle), `swap`/`diagonal` (integer predicates), the quadrant sign table.
Nothing is "whatever the reduction order produced". Contrast the status quo:
`np.linspace(0, 2π, n, endpoint=False)` produces an *unnamed* angle array whose error is
argument-reduction noise nobody tracks. **PASS.**

**Criterion 2 — verified against a structurally-independent reference.**
The reference is **`mpmath` at 50–100 dps**, which shares nothing with either implementation
above the trusted-library line (`algebra-of-record` §"structural independence applies above the
trusted-library line": `mpmath.cos` is trusted upstream; the *reduction* is the load-bearing
part and here there is none — the reference evaluates `cos(2π p/q)` directly at 100 digits).
`[M]`:

```
  max |circle_nodes - exact| = 1.2573e-16  = 0.57 ulp(1.0)
  max |legacy       - exact| = 8.2562e-16  = 3.72 ulp(1.0)
  components compared           : 1040
  circle_nodes strictly closer  : 592
  legacy       strictly closer  : 30
  exact tie                     : 418
  worst legacy-advantage: 1.1102e-16 absolute (new err 1.2573e-16 vs legacy 1.4709e-17) at n=128, k=17
```

**The honest claim is not "always more accurate".** It is: *the new value is within
`0.57 ulp(1.0)` of the true value everywhere; the legacy is within `3.72 ulp(1.0)`; on the
`30/1040 = 2.9 %` of components where the legacy happens to be closer, the new value is still
within `1.26e-16`.* **PASS**, stated that way.

> ⚠ Two claims I must correct rather than inherit.
> - **The issue body's "agreement with `np.cos`/`np.sin` to ≤1 ULP" is FALSE.** `[M]` max
>   `|new − np.cos(2πp/q)| = 8.33e-16 = 3.75 ulp(1.0)`. It is also the wrong reference —
>   `np.cos(2πp/q)` is not the true value. Restate the criterion as *"≤1 ulp of the
>   arbitrary-precision value"*, which `[M]` holds at `0.57 ulp`.
> - **The main agent's "max node drift 9.99e-16 at `n_phi=48`" does not reproduce.** `[M]` over
>   `n_phi ∈ {3…128}` the max is **8.326673e-16 at `n_phi=24`** (also `8.33e-16` at 48 and
>   100). `9.99e-16` appears only at `q=360`, outside any plausible `n_phi`. The issue body's
>   "8.3e-16 at n=24" is **correct**; use it. Likewise the "312 vs 8" split is sweep-dependent
>   — over `q ∈ {3…128}` I measure **592 / 30 / 418 tie**. Cite the ratio, not the counts.

**Criterion 3 — the drift is FP-non-associativity, dimensionally explainable.**
Node values: `[M]` `≤ 8.33e-16 ≈ 3.75 ulp(1.0)` — a single-step computation with reduction
depth ~3 (`4p`, `divmod`, one multiply). **PASS.**
Solver values: `[M]` `A vs B` (value drift, ordering held) `= 1.05e-14` relative `1.06e-14`,
against `inner_tol = 1e-13`. Bounded by `iteration_count × condition × ULP` and **well inside
the `SAFETY(10) × conv_tol` regression gate**. **PASS.**
⛔ **The ORDERING drift `1.0e-2` FAILS criterion 3 outright** — it is not FP noise, it is an
algorithmic change. It must not be smuggled into the same re-baseline. §7 / §6.

### 5.2 Snapshot store inventory

| Store | Files | Regenerated by | Compared by | Touches `product`? |
|---|---|---|---|---|
| `tests/sn/regression/snapshots/` | 20 `.npz` | `_generate_snapshots.py` | `_regression_assert.py` (`SAFETY×conv_tol` iterative / `nulp(reduction_depth)` direct + `DriftWarning`) | **YES** — `cyl_1g_homogeneous_product_dd_n20.npz`; generator lines 157 & 186 `Quadrature.product(n_mu=2, n_phi=4)` |
| `tests/sn/regression/snapshots/walk_matvec_*.npz` | 4 | `_generate_walk_baselines.py` | `test_walk_matvec_baselines.py` | `walk_matvec_cyl_2g.npz` — verify the quadrature |
| `tests/sn/regression/snapshots/2d_octant_*.npz` | 6 | `_generate_2d_octant_snapshots.py` | `test_dd_regression.py` | LS4/LS6/Lebedev5 — **no** |
| `tests/geometry/snapshots/` | 7 `.npz` | `_generate_bc_equivalence_snapshots.py` | `test_bc_equivalence_snapshot.py` | LS4/LS6/GL/Lebedev17 — **no** (B3.4b re-anchored these to derived expressions) |
| `tests/numerics/data/step6_product_reroute_baseline.npz` | 1 | **hand-pinned** (no generator; referenced only at `tests/numerics/test_operator.py:726`) | `test_operator.py` | **name says `product` — audit it** |
| `tests/sn/_data/bc_extraction_baseline/*.npy` | many | `test_bc_extraction_matvec.py` / `_2d.py` write them | same files | audit |
| `tests/sn/_fixtures/wave_t_t{3,4}/` | `_capture_pre_t{3,4}_snapshots.py` | those scripts | wave-T tests | audit |
| `tests/sn/sweep_ref_2g.npy` | 1 | unknown — **hand-pinned** | `test_sweep_regression.py` | audit |
| `tests/_harness/pyright_baseline.json` | 1 | pyright | harness | no |

### 5.3 The re-baseline mechanics

1. **Use the existing tripwire as the discovery instrument.** `_regression_assert.py` already
   emits `DriftWarning` for any sub-tolerance bit movement and documents
   `-W error::DriftWarning` as the strict bit-identity contract. Run the affected suites with
   `-W error::DriftWarning` **before** touching any `.npz`: the failures ARE the exact
   re-baseline list, derived rather than guessed. No new machinery is needed.
2. **Regenerate only what that run names.** Everything else must stay bit-identical, and
   proving that is half the evidence.
3. **Per-snapshot justification record.** Append to
   `tests/sn/regression/README.md` (and the geometry twin) one row per regenerated file:

   | field | content |
   |---|---|
   | file | `cyl_1g_homogeneous_product_dd_n20.npz` |
   | old → new | `keff`/`flux` max ULP distance and absolute delta |
   | why not bit-identical | "node values moved ≤ 8.33e-16 (`[M]`, `roots_of_unity` vs `linspace`+`cos`)" |
   | criterion 1 | the named-intermediate chain (§5.1) |
   | criterion 2 | the **independent** anchor used (see §5.4 — NOT the old snapshot) |
   | criterion 3 | measured drift vs `SAFETY × conv_tol` |
   | ordering | **`unchanged` / `changed` — and if changed it belongs to the §7 commit, not this one** |

### 5.4 ⚠ CONTAMINATED RE-BASELINE — one confirmed, one to audit

**CONFIRMED — `tests/numerics/test_rules_product.py::test_product_bit_identical_to_legacy_adapter`
is a tautology.** Its docstring says *"Bit-identical match against
`orpheus.sn.quadrature.ProductQuadrature.create`. Pins the cylindrical regression snapshots
(`cyl_*_product_*.npz`)"* — but its body compares `product_mu_phi(n,m)` against
`Quadrature.product(n_mu=n, n_phi=m)`, and `directional.py:613` is literally
`measure, structure = product_mu_phi(n_mu, n_phi)`. `[M]`:

```
=== does ProductQuadrature (the named 'legacy adapter') exist? ===
  NO -> ModuleNotFoundError: No module named 'orpheus.sn.quadrature'

=== CONTROL: unmutated ===
  test passes: True

=== MUTATION: product_mu_phi replaced by RANDOM nodes ===
  test STILL passes: True
  -> TAUTOLOGY CONFIRMED: the gate cannot red.

=== positive control: the mutation really did bite ===
  mutated Quadrature.product mu_x[:3] = [0.18881712 0.16021416 0.7415052 ]
  restored                mu_x[:3] = [ 5.08374127e-01  3.11289374e-17 -5.08374127e-01]
```

This is the B3.4b/§16 pattern exactly: a gate whose "independent" reference is the code under
test. Consequences:
- **The claim "pins the cylindrical regression snapshots" is FALSE** and must be deleted.
- **The snapshot has no independent anchor today.** Regenerating it and comparing only against
  the old value would be a contaminated re-baseline.
- **Required replacement (Gate S1):** rewrite the test to compare `product_mu_phi` against an
  **expression derived from the definition**, in the test body, from
  `roots_of_unity` + `leggauss` — `μ_x = √(1−μ²)·cos`, `μ_y = √(1−μ²)·sin`, `μ_z = μ`,
  `w = w_GL · 2π/n_φ` — plus a `Σw == 4π` closed-form leg and the mpmath anchor for the node
  values at one `(n_mu, n_phi)`. That is the same inversion B3.4b applied to the bc_equivalence
  generator: **derive from the equation, do not re-read the producer.**

**TO AUDIT — `tests/numerics/data/step6_product_reroute_baseline.npz`.** Named for `product`,
hand-pinned, no generator script, consumed at `tests/numerics/test_operator.py:726`. Before
regenerating, establish what it is a baseline *of* and whether its regeneration path reads the
same code. If it does, it needs an independent expression before it moves.

**Also audit** every `_generate_*.py`: a generator that imports the module under change and
records its output is a recorder, not a reference. `_generate_snapshots.py` builds
`Quadrature.product(...)` and runs the **solver** — that is legitimate (the snapshot is a
solver-output regression baseline, and its independent anchor is the closed-form
`k_∞ = νΣ_f/Σ_a` / `Q/Σ_t` leg elsewhere in the suite), **provided** those independent legs
exist and are not themselves ordering-blind (§2.6 says the `Q/Σ_t` leg IS ordering-blind).

---

## 6. Deliverable (6) — ORDER OF OPERATIONS (green at every commit boundary)

Ten steps. Steps 1–3 are pure additions (nothing can red). Step 4 is the **ordering ruling**
and must precede any node-value change. Steps 5–8 carry the value change.

| # | Commit | Content | Gate at the boundary | Can anything red? |
|---|---|---|---|---|
| **1** | `test(numerics): intrinsic-property gates for roots_of_unity` | **ALREADY WRITTEN** — `tests/numerics/test_roots_of_unity.py`, §1.4. 17 gates, 251 params, `[M]` green under `-O` and plain, all 17 mutation-proven. NO production edit. | `[M]` `251 passed in 0.52s`; teeth matrix §1.4 in the commit body. | No — pure addition against existing code. |
| **2** | `docs(numerics): repoint the dangling eps xref` | Fix `angular_trace_space.py:90` + `:159` → `test_eps_sits_in_the_round_off_to_genuine_gap`. Docstring only. | `sphinx-build -W` clean; `pytest tests/numerics -q`. | No. |
| **3** | `test(numerics): registry-walk gates (currently xfail for product)` | §3.2 T1/T3 + §4.1 R1, written as registry walks. `product` rows land as **`xfail(strict=True)`** with `reason="#325 step 6 — product still uses linspace+cos"`. | Suite green with the xfails; `--runxfail` and READ each message to confirm it reds for ITS documented reason (vv Mode 8 class 4). | No — the xfail set IS the todo list. |
| **4** | ⛔ **`fix(numerics): make the cylindrical per-level ordering total`** | §7. Ruling + injective sort key + `test_level_ordering_is_deterministic`. **This step changes results by ~1 % and re-baselines `cyl_*_product_*` on its own.** | §7 gates. Independent anchor required (§7.4). | **YES, deliberately** — this is the step that moves numbers, and it must move them alone. |
| **5** | `test(numerics): replace the tautological product pin` | §5.4 Gate S1 — derive the expected product nodes from the definition. | Mutation: replace `product_mu_phi` with random nodes → **S1 must now RED**. Record it. | No (still on legacy nodes). |
| **6** | `feat(numerics): export roots_of_unity; repoint the three consumers` | `orpheus/numerics/__init__.py` export; `rules_product.py:115`, `moc/quadrature.py:89`, `derivations/discrete/sn/balance.py:367`. | Run affected suites with **`-W error::DriftWarning`** — the failures are the re-baseline list (§5.3). | **YES** — snapshots move. |
| **7** | `test(sn,geometry): re-baseline the product-derived snapshots` | Regenerate exactly what step 6 named; write the §5.3 justification rows. | `pytest tests/sn tests/geometry tests/numerics -q` green (serial). | No, by construction. |
| **8** | `test(numerics): flip the step-3 xfails` | Delete the `xfail(strict=True)` markers — they XPASS now. Add §3.2 T2 (the demotion statement) and §4.3 R2 (the odd-`n_phi` production guard) + its negative test. | T2 mutation: monkeypatch `roots_of_unity` → `linspace`+`cos`, confirm T1/T2/R1 red. | No (strict-xfail deletion is forced by XPASS). |
| **9** | `docs(numerics): TANGENTIAL_EPS is now a defensive guard` | Rewrite the `angular_trace_space.py` docstring + the constant's comment to state the demotion, citing T2. Fix the `axis_aligned` flag question (§3.3) or record why not. | `sphinx-build -W`. | No. |
| **10** | `docs(theory): the PRINCIPLE — generate by the group action` | The issue's last acceptance criterion. Name `lebedev` / `level_symmetric` as the families that already obeyed it. **Dispatch the archivist.** | `sphinx-build -W`; Nexus `staleness` clean. | No. |

**Full-tree gate:** once, after step 8, `python -O -m pytest -m "not slow"` **serial**
(≈52 min per `reference_test_execution_env`). Not after every step.

**Why step 4 precedes step 6.** If they land together, the `cyl_*_product_*` snapshots move by
`~1e-2` and the commit message says "node values moved ≤1 ULP". That is a false justification
record — and per §5.1 criterion 3 the `1e-2` drift is not FP noise at all. Separating them
gives each its own honest story and lets the DriftWarning run in step 6 be interpretable.

---

## 7. Deliverable (7) — the ordering gate

### 7.1 The finding, restated as a claim

`[M]` **`rules_product.py:139` `order = np.argsort(mu_x[level_arr])` produces an ordering that
is (a) not determined by the physics, (b) not determined by the sort key, and (c) currently
determined by float noise that #325 removes.** The consequence is a **1.0–1.8 % angular-flux
change that does not converge away** (§0.1).

`[M]` It is already latent today: `n_phi ∈ {20, 32, 64}` carry 1/2/4 tied pairs under the
legacy nodes. #325 takes `n_phi=64` from 4 tied pairs to **31**.

### 7.2 The sort key is not injective — and that is the real defect

`η = sinθ·cos φ` is 2-to-1 over `φ ∈ [0, 2π)`: `φ` and `2π−φ` share `η` and differ only in
`ξ = μ_y`. The Bailey Eq. 50 / Hébert §3.9.4 recursion
`α_{m+1/2} = α_{m−1/2} − w_m η_m` (`orpheus/geometry/reduced_operator.py:778-786`) is a
**sequential running sum** over that order, and `pole_angular_closure` then hands ordinate `m`
the interval `(α_{m−1/2}, α_{m+1/2})`. Two ordinates with equal `η` and equal `w` leave the
α **array** unchanged under a swap but receive **different intervals** — hence the 1 %.

**⚠ This is a physics ruling for the user, not a test-design call.** Three candidate keys:

| key | deterministic | injective | physical reading |
|---|---|---|---|
| (a) `kind="stable"` on `η` | yes | no | "ties resolve to increasing `φ` index" — an *implicit* convention |
| (b) lexicographic `(η, ξ)` | yes | **yes** | mirror partners ordered by azimuthal sign — explicit |
| (c) sort by `φ` directly | yes | **yes** | the α recursion IS an azimuthal sweep; `η` is a proxy for it |

**Recommendation: (c), with (b) as the fallback.** Reason: the recursion's own derivation
(`balance.py::verify_alpha_closure`, "symmetry of `η = sinθ cos φ` over equally-spaced `φ`")
describes an **azimuthal** sweep. Sorting by `η` is sorting by a monotone-only-on-a-half-turn
function of the real variable, and its 2-to-1-ness is exactly why the order is ambiguous.
Sorting by `φ` makes the key injective by construction and makes the closure
`α_{1/2} = α_{N+1/2} = 0` read directly. **(a) is NOT recommended** — it makes the answer
reproducible without making it *determined*, and it leaves the next reader with the same
under-determination, now silent.

**Whichever is chosen, the closure `α_{N+1/2} = 0` must be re-verified** — `[M]`
`balance.py::verify_alpha_closure` currently asserts `|Σ w η| < 1e-14`, which is
order-independent, so it is **Mode-12 blind to the whole question** and cannot be cited as
evidence.

### 7.3 The gates

**Gate O1 — `test_level_ordering_is_a_determined_permutation`** (`tests/numerics/`,
`@pytest.mark.foundation`). The measured functional is the **full index tuple**, not
sortedness (§2.5):
- For every `(n_mu, n_phi)` in `{2,4,6,8} × {3,4,5,6,8,12,16,20,24,32,64}`, assert
  `tuple(level_indices[p])` equals a tuple **independently constructed in the test** from the
  chosen key — e.g. for (c), `np.argsort(phi_index)` restricted to the level.
- **NOT** `np.all(np.diff(eta) >= 0)`. `[M]` legs C and D are both sorted and differ by 1.8 %.

**Gate O2 — `test_level_ordering_is_sort_kind_invariant`.** For each configuration, build the
ordering with `kind ∈ {"quicksort", "stable", "heapsort", "mergesort"}` and assert all four are
**bit-identical index tuples**. `[M]` today this REDS at `n_phi ∈ {32, 64}`:

```
  n_phi= 16: all kinds agree = True
  n_phi= 32: all kinds agree = False
      quicksort: (16, 15, 17, 18, 14, 13, 19, 20, 12, 21, 11, 22, 10, 9, 23, 8, …)
         stable: (16, 15, 17, 14, 18, 13, 19, 12, 20, 11, 21, 10, 22, 9, 23, 8, …)
  n_phi= 64: all kinds agree = False
```
This is the gate that proves the key is injective — a non-injective key cannot pass it.

**Gate O3 — `test_level_ordering_key_is_injective`.** Directly: the sort key array has `q`
distinct values per level. The structural statement behind O2.

**Gate O4 — the α-closure re-verification.** `α_{N+1/2} == 0` to machine precision under the
CHOSEN ordering, for every configuration, **and** a companion asserting the α **array** is not
identically zero (non-vacuity — a level whose α dome is flat proves nothing).

**Gate O5 — the ordering-sensitivity regression, with an independent anchor.** The
`Q/Σ_t` flat-flux leg is designed-green here (§2.6). The anchor must be one of:
- **angular convergence to a common limit** — solve at `n_phi ∈ {16,32,64,128}` and assert the
  heterogeneous flux converges; `[M]` today the C-vs-D gap is **flat at ~1e-2** across
  `n_phi ∈ {32,64,128}`, so a converging sequence under the chosen key IS the evidence the
  ambiguity is gone;
- **cross-family agreement**: `level_symmetric` carries no η ties (`[M]` `level_symmetric(2…8)`
  all pass the exact-zero probe and the reflection exactness), so a matched-accuracy
  `product` vs `level_symmetric` comparison is a genuinely independent cross-check — L4-class,
  so it corroborates, it does not verify.

**Mandatory config for O5** (config-blindness, §2.6): **heterogeneous** (≥2 regions),
**≥2 groups**, **anisotropic per-ordinate source**, **`n_phi ≥ 32`**. `[M]` at `n_phi ∈ {8,16}`
the sort kinds happen to agree and the gate reads `0.000000e+00` — a designed-green row.

### 7.4 Where it lands

**Step 4 of §6, before the node-value change, as its own commit**, because:
- it changes results by `~1e-2` **independently** of the `≤1 ULP` node shift;
- its justification record is a *physics ruling*, not an FP-drift argument, and the two must
  not share a commit message;
- landing it first means step 6's `DriftWarning` run measures **only** the node-value drift,
  which is what makes that run interpretable.

---

## 8. Open items requiring a USER ruling

| # | Question | Why it is not mine to decide |
|---|---|---|
| **U1** | **The ordering key**: `(c)` sort by `φ`, `(b)` lexicographic `(η, ξ)`, or `(a)` `kind="stable"`? | It changes the cylindrical α-recursion's physical reading, and the converged answer by ~1 %. |
| **U2** | Should `_find_reflections` **raise** on a non-exact match (§4.3 R2)? It would make `Quadrature.product(4, 5)` an error at construction. | A construction-time refusal narrows the API; `[M]` `_product_invert(4)` can mint `n_phi=5`, so the selector would need a matching constraint. |
| **U3** | Is refusing **`bool` numerators** intended (`np.issubdtype(np.bool_, np.integer)` is `False`)? | Pin it either way; today it raises. |
| **U4** | Tighten `registry.py` `ProductQuadrature.axis_aligned` from the conservative `False` to `4 ∣ n_phi`? | The selector reads it; a now-false conservative tag is a decayed claim. |
| **U5** | Does the theory page (issue acceptance criterion 6) land in `docs/theory/methods/sn/` or a new `docs/theory/foundations/` anchor? | Archivist + corpus-consistency call. |

---

## 9. Corrections to the record

| Source | Claim | Verdict |
|---|---|---|
| Issue body, acceptance criteria | "agreement with `np.cos`/`np.sin` to ≤1 ULP" | **WRONG.** `[M]` `3.75 ulp(1.0)` against numpy; `0.57 ulp(1.0)` against mpmath. Restate against the true value. |
| Issue body, "Verified (prototype)" | "Node values move by ≤ 1 ULP (max 8.3e-16 at n=24)" | **Number CORRECT** (`[M]` `8.326673e-16 at n_phi=24`); **"≤1 ULP" mislabels it** — `8.33e-16 = 3.75 ulp(1.0)`. |
| Issue body, Blast radius | "Values shift ≤1 ULP, so anything pinned bit-exactly re-baselines" | **INCOMPLETE and misleading.** `[M]` solver results shift `~1e-2` via the level ordering — twelve orders larger. §0.1. |
| Issue body, cost item 2 | "`reflection_index` … the specular involution is, on `product`, approximate" | **Values approximate (`[M]` `3.5e-16`), but the resulting MAP is exact** — `[M]` perm ∧ involution hold to `n_phi=1024`, margin `5.0e-3`. Not a live mis-pairing bug. §4.1. |
| Crosswalk §14.3 | "`reflection_index`'s nearest-neighbour search … becomes an exact lookup" — implying replacement | **Keep the search.** Replacing it mints a twin path (lebedev has no index formula) and deletes the only detector of the odd-`n_phi` defect. §4.3. |
| Main agent, round 2 | "max node drift 9.99e-16 at `n_phi=48`" | **Does not reproduce.** `[M]` `8.326673e-16` at `n_phi=48`; the max over `{3…128}` is `8.326673e-16` at `n_phi=24`. `9.99e-16` occurs at `q=360`. |
| Main agent, round 2 | "closer on 312 components, further on 8" | **Direction confirmed, counts sweep-dependent.** `[M]` over `q ∈ {3…128}`: **592 closer / 30 further / 418 tie**. Cite the ratio. |
| `angular_trace_space.py:90`, `:159` | `` :func:`test_eps_below_min_genuine_cosine` `` | **Dangling xref** — no such test. §3.1. |
| `angular_trace_space.py` docstring | "nominally-tangential cosines are **exactly** `0.0`" | **FALSE for `product` today.** True after step 6. §3.1. |
| `tests/numerics/test_rules_product.py` | "Bit-identical match against `ProductQuadrature.create`. Pins the cylindrical regression snapshots" | **TAUTOLOGY + false coverage claim.** `ProductQuadrature` does not exist. §5.4. |
| `derivations/.../balance.py::verify_alpha_closure` | cited as the α-closure algebra of record | **Mode-12 blind to the ordering question** — `|Σ w η| < 1e-14` is order-independent. §7.2. |
