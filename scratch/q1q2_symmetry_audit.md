# `orpheus/numerics/symmetry.py` — correctness audit

**Scope**: pre-carve defect hunt (the carve widens `_orbit_closure` to return its
orbit permutation and adds an orbifold singular-set concept).
**Method**: every claim carries a numerical demonstration run under
`.venv/bin/python -O`. Probes: `/Users/rodrigo/.claude/jobs/c30e4f25/tmp/p{1..12}*.py`.
**Legend**: CONFIRMED = demonstrated numerically. SUSPECTED = reading only.
CLEAN = probed, found correct.

**Live-severity qualifier (measured).** `is_invariant` has **zero production
consumers** — `grep -rn "\.is_invariant(" --include=*.py .` returns only 7 sites,
all in `tests/numerics/test_symmetry.py`. The quadrature-selection layer
(`quadrature/registry.py:677`) reads the **declared** `spec.invariance_group` tag
and calls `is_subgroup_of`, never `is_invariant`. So the false-certification
defects below are **latent today and live the moment the carve wires
`_orbit_closure`'s permutation into anything**. The containment defects (S-7,
S-9, S-10) ARE on the live selection path.

---

## Severity index

| # | Sev | Site | Defect |
|---|-----|------|--------|
| S-1 | **CRITICAL** | `:769-779`, `:665-675` | `_so2_representatives()` samples only `C_4` ⟹ a non-`SO(2)` rule certifies `SO(2)`-invariant |
| S-5 | **CRITICAL** | `:904-954` | `_orbit_closure` never verifies the match is a **bijection** ⟹ a non-invariant measure certifies invariant |
| S-2 | HIGH | `:582-616` | 1-D path: `SO3`/`Cn` return `True` unconditionally; every tag returns `True` for a symmetric 1-D rule |
| S-3 | HIGH | `:683-697` | `SO3` and `O3` use the *same* 120-element `I_h` op set ⟹ `SO3.is_invariant ≡ O3.is_invariant`, and neither is `SO(3)` |
| S-7 | HIGH | `:515-518` | `D_nh ⊆ O(2)` asserted **TRUE**, false under both embeddings of `O(2)`; pinned by a test |
| S-9 | MOD | `:204` | `Z2 ⊆ SO3` asserted TRUE while the module's own `Z2` realisation is **improper** (det = −1) |
| S-10 | MOD | `:486`, `:511-513` | `Cn(1)` ≠ `Trivial` in the lattice ⟹ 1 antisymmetry + 3 transitivity violations |
| S-8 | MOD | `:359-364` | Docstring: "neither `C_n` nor `D_nh` sits inside `O_h`, `I_h`" — false, and self-contradictory in the same sentence; 10 true relations missing |
| S-6 | MOD | `:947-951` | Bucket-hit branch can select a **farther** `j` than the brute-force branch ⟹ false negative |
| S-11 | MOD | `:945`,`:950`,`:952` | Node window `100·atol`, weight window `1·atol`; both absolute, neither scales |
| S-4 | MOD | `:619-648` | `_is_reflection_invariant_1d` also non-bijective; node window `10·atol` (a third, different constant) |
| S-12 | LOW | `:428-438` | `is_invariant` docstring describes fingerprint strategies the code does **not** implement (and the fingerprint is provably necessary-only) |
| S-13 | LOW | `:558-562` | Stringly-typed `support` dispatch overrides the geometric `dim` |
| S-14 | LOW | `:926`, `:400` | Public `atol` unvalidated; `atol ≤ 1e-19` silently saturates `int64` |
| S-15 | LOW | `:732-747` | `_rotation_x` / `_rotation_y` are dead code |
| — | **CLEAN** | `:750-896` | **All operator generators verified correct** — see §CLEAN |

---

## S-1 (CRITICAL) — `_so2_representatives()` falsely certifies `SO(2)` invariance

`symmetry.py:769-779` (generator) → `:665-670` (`SO2` dispatch) → `:672-675`
(`O2` dispatch inherits it).

### What is wrong

`_so2_representatives()` returns **four** rotations about z: `0°, 90°, 180°, 270°`.
The first is the **identity** (a vacuous check), so the set has exactly **three**
non-trivial elements, all in `C_4`.

`_orbit_closure` is a *generator-set* check: closure under every listed matrix
implies closure under every **product**, i.e. under the group they generate. That
makes the strategy sound for `C_n`, `D_nh`, `O_h`, `I_h` (whose listed sets do
generate the intended group) and **unsound for `SO2`/`O2`**, whose listed set
generates only `C_4`.

Consequence: **any** `C_4`-invariant rule certifies as `SO(2)`-invariant. A
product quadrature with `n_phi ∈ {4, 8, 12, 16, …}` is exactly such a rule.

The `is_invariant` docstring's framing ("*necessary but not sufficient* … but
**sufficient by construction** for the rules ORPHEUS ships") is **false for the
product family**, one of the four rules ORPHEUS ships. It is not a conservative
upper bound — it is an unconditional `True` on `n_phi ≡ 0 (mod 4)`.

### Numerical demonstration (`p1_so2.py`)

```
=== the representative set ===
n ops: 4
  op[0]: rotation about z by    0.000 deg   det=+1.000000000000000  is_identity=True
  op[1]: rotation about z by   90.000 deg   det=+1.000000000000000  is_identity=False
  op[2]: rotation about z by  180.000 deg   det=+1.000000000000000  is_identity=False
  op[3]: rotation about z by  270.000 deg   det=+1.000000000000000  is_identity=False

=== product quadratures: DECLARED SO2 vs is_invariant(SO2) vs TRUTH ===
                        rule     N   declared  is_inv(SO2)   C_nphi only?
     product(n_mu=4,n_phi=4)    16        SO2         True           True
     product(n_mu=4,n_phi=8)    32        SO2         True           True
    product(n_mu=4,n_phi=12)    48        SO2         True           True
    product(n_mu=4,n_phi=16)    64        SO2         True           True
     product(n_mu=4,n_phi=3)    12        SO2        False           True
     product(n_mu=4,n_phi=5)    20        SO2        False           True
     product(n_mu=4,n_phi=6)    24        SO2        False           True
     product(n_mu=4,n_phi=7)    28        SO2        False           True
     product(n_mu=2,n_phi=2)     4        SO2        False           True

=== the decisive falsification ===
rule            : Quadrature.product(n_mu=4, n_phi=8), N = 32
is_invariant(SO2)          : True
closure under rot_z(90 deg): True
closure under rot_z(45 deg): True   <- 2*pi/n_phi = 45 deg, still a symmetry
closure under rot_z(22.5 d): False   <- HALF the azimuthal spacing: NOT a symmetry
closure under rot_z(1 deg) : False

=> is_invariant(SO2) says True; the rule is only C_8-invariant.

=== does it also falsely certify C_n for n it does NOT have? ===
  is_invariant(Cn( 2)) =  True   (truth: n | 8  ->  True)
  is_invariant(Cn( 3)) = False   (truth: n | 8  ->  False)
  is_invariant(Cn( 4)) =  True   (truth: n | 8  ->  True)
  is_invariant(Cn( 5)) = False   (truth: n | 8  ->  False)
  is_invariant(Cn( 6)) = False   (truth: n | 8  ->  False)
  is_invariant(Cn( 8)) =  True   (truth: n | 8  ->  True)
  is_invariant(Cn(16)) = False   (truth: n | 8  ->  False)
```

The `C_n` sibling is **sound** (its generator set really does generate `C_n`),
which makes the `SO2` answer self-inconsistent with it: `Cn(16)` is correctly
rejected while `SO2 ⊋ C_16` is accepted. See S-16 for the full monotonicity tally.

### Coverage: the path has ZERO tests

`p11_bite.py`, instrumenting `_so2_representatives` while running the real test:

```
=== MUT-SO2: why the counter read 0 ===
  test_gauss_legendre_1d_so2_invariant does: SO2.is_invariant(1-D measure) -> True  ; _so2_representatives calls: 0
  an S^2 measure would reach it:            SO2.is_invariant(S^2 measure) -> True  ; _so2_representatives calls: 1
  => the ONLY SO2 test in the suite takes the 1-D short-circuit; the 3-D SO(2)
     path (and with it the C_4-sample defect S-1) has ZERO test coverage.
```

Whole-suite mutation (widen the sample `C_4 → C_8` via an in-process plugin):

```
### BASELINE (no mutation)
182 passed, 1 warning in 1.78s
### MUTATION: mut_so2
[MUT-SO2] _so2_representatives called 0 time(s) during the run
182 passed, 1 warning in 0.96s
```

0 calls across 182 tests — the function is **never executed** by the suite.

### Minimal correct fix

A **finite** node set on `S²` is `SO(2)`-invariant iff every node lies on the
rotation axis (the `SO(2)` orbit of an off-axis point is a circle — uncountable):

```python
def _so2_invariant(nodes, weights, atol) -> bool:
    return not (np.hypot(nodes[:, 0], nodes[:, 1]) > atol * 100).any()
```

The rules that *declare* `invariance_group=SubgroupOfO3.SO2`
(`rules_1d.py:84`, `rules_product.py:147`, `registry.py:438`, `:474`) assert a
property of the **continuum measure they discretise**, not of the discrete
measure. That is a legitimate but *different* claim and must not share a
predicate with the discrete one. The honest discrete answer for a product rule is
`Cn(n_phi)`, which the existing `Cn` path already computes correctly.

If instead the intent is a deliberately-conservative *screen*, rename it
(`may_be_invariant`) and document the `False`-only-is-meaningful semantics —
the current name is consumed as a certificate.

---

## S-5 (CRITICAL) — `_orbit_closure` never verifies the match is a bijection

`symmetry.py:904-954`.

### What is wrong

The docstring claims: *"find a permutation π such that M(nodes)_i = nodes_{π(i)}"*.
The body computes, for each `i`, **some** `j` with `‖nodes[j] − M·nodes[i]‖ ≤
100·atol` and `|w_j − w_i| ≤ atol`. It never checks that `i ↦ j` is **injective**
(hence, on an equal-size set, bijective). A many-to-one match satisfies every
assertion in the loop.

Invariance is `M#µ = µ` as *measures*. A many-to-one match certifies a set whose
pushforward has different **multiplicities** — i.e. a genuinely non-invariant
measure. No tolerance games are needed: an **exact duplicate node** suffices.

This is the defect that directly blocks the carve. If `_orbit_closure` is widened
to *return* the orbit permutation, the returned array is **not a permutation** on
exactly the inputs this check wrongly accepts.

### Numerical demonstration (`p5_orbit_closure.py`)

```
=== S-5: EXACT-DUPLICATE node makes the 'permutation' many-to-one ===
  N: 24 -> 25   total mass: 12.566370614359172 -> 13.089969389957471
  positive control: node 0 and node 24 identical?  True
  is_invariant(O_h)                    : True
  instrumented verdict                 : True
  ops whose match map is NON-INJECTIVE : 48 of 48
     e.g. op[0]: the map i->j hits only 24 distinct targets out of 25 sources

  TRUTH: this measure is NOT O_h-invariant.  Under any M with M.p0 != p0 the
  pushforward puts mass 2*w0 at M.p0 where the measure itself has only w0.
     p0        = [-0.408248 -0.816497 -0.408248]
     M.p0      = [0.408248 0.408248 0.816497]   (M.p0 == p0 ? False)
     mass at p0 in mu     = 1.047197551196598
     mass at p0 in M#mu   = 0.523598775598299
  => M#mu != mu, yet is_invariant returns True
```

The instrumented copy is byte-for-byte the production body plus a record of
`i → j`; it reproduces the production verdict exactly on the unperturbed rule
(`instrumented verdict: True   production verdict: True`, `anomalies logged: none`).

### The fix is free for every shipped rule (`p11_bite.py`)

```
=== does the FIX (require a bijection) break any shipped rule? ===
  lebedev(7)   vs O_h    production= True   with-bijection= True
  lebedev(17)  vs O_h    production= True   with-bijection= True
  lebedev(23)  vs O_h    production= True   with-bijection= True
  LS(2)        vs O_h    production= True   with-bijection= True
  LS(4)        vs O_h    production= True   with-bijection= True
  LS(6)        vs O_h    production= True   with-bijection= True
  LS(8)        vs O_h    production= True   with-bijection= True
  product(4,8) vs C_8    production= True   with-bijection= True
  => the bijection requirement is FREE for every shipped rule (no false negatives).
```

### Minimal correct fix

Collect `j` into an array and assert `len(set(jj)) == n` per operator (equivalently,
mark targets consumed and reject a second claim). This is the natural shape for the
carve anyway: build the permutation, then *validate* it is one.

---

## S-2 (HIGH) — the 1-D path certifies every group for a symmetric 1-D rule

`symmetry.py:582-616`.

`_check_invariance_1d` splits tags into `rotational_only` (`SO2`, `SO3`, `Cn`) →
unconditional `True`, and the rest → `μ → −μ` closure. Both halves are wrong.

- **`SO3` → `True` unconditionally is unsound.** `SO(3)` contains the π-rotation
  about x, whose action on `μ = cos θ` is `μ → −μ`. So `SO(3)`-invariance of a
  μ-measure *requires at minimum* reflection closure — and in fact requires the
  measure to be uniform (the only `SO(3)`-invariant measure on `S²`). Returning
  `True` for an asymmetric μ-rule is a false certificate, not a conservative one.
- **`Oh` / `Ih` → reflection-only is unsound.** A μ-measure carries no azimuthal
  information; `O_h`/`I_h` invariance is simply not decidable from it. The code
  answers `True` whenever the μ nodes are `±`-paired.

### Numerical demonstration (`p4_1d.py`)

```
legend: T = is_invariant returns True, . = False
  measure                                   Trivial     Z2    SO2     O2     Oh     Ih    SO3     O3  Cn(3) Dnh(3)
  gauss_legendre_on_mu(8)   [symmetric]           T      T      T      T      T      T      T      T      T      T
  asymmetric [0.1,0.5,0.9], w=[1,1,1]             T      .      T      .      .      .      T      .      T      .
  nodes +/-0.5 but weights [1,2]                  T      .      T      .      .      .      T      .      T      .

=== S-2  monotonicity violation on the 1-D path ===
Lattice asserts Z2 subset-of SO3   : True
On the ASYMMETRIC measure:
   Z2.is_invariant  = False  (correct: no x->-x partner)
   SO3.is_invariant = True  (FALSE POSITIVE)
   O3.is_invariant  = False
   Cn(3).is_invariant = True  (FALSE POSITIVE for a mu-measure? see below)

SO(3) contains the pi-rotation about x, whose action on mu = cos(theta) is mu -> -mu:
   R_x(pi) =
 [[ 1.  0.  0.]
    [ 0. -1. -0.]
    [ 0.  0. -1.]]
   => (x,y,z) -> (x,-y,-z);  mu = z  ->  -mu.
   So SO(3)-invariance of a mu-measure REQUIRES mu -> -mu closure, which the
   asymmetric measure fails.  Returning True is unsound, not conservative.
```

The first row is the headline: for the shipped `gauss_legendre_on_mu(8)` rule,
**every one of the ten tags returns `True`** — including `I_h`, which the suite's
own negative test (`test_lebedev_is_NOT_icosahedral_invariant`) exists to prove
is *not* satisfied by an `O_h` rule.

### S-3 corollary, and the contradiction with the containment half

```
=== S-3  1-D O_h / I_h certification vs the containment lattice ===
is_invariant(O_h) on the shipped 1-D GL rule  : True
is_invariant(I_h) on the shipped 1-D GL rule  : True
registry declares GaussLegendre1D invariance_group = SO2; and
  Oh.is_subgroup_of(SO2) = False  <- the SELECTION layer says a cartesian2d(O_h) geometry may NOT use GL1D
  but is_invariant says the SAME rule IS O_h-invariant.  The two public
  predicates of this module contradict each other on the same object.

Honest test: lift the SAME GL mu-rule to S^2 (axisymmetric, n_phi=8) and run the
real 48-op O_h closure on the lift:
   _orbit_closure(lift, _octahedral_ops()) = False  <- NOT O_h-invariant
   is_invariant(O_h) on the 1-D form       = True  <- but the 1-D form certifies True
```

### Minimal correct fix

Decide the semantics of a 1-D measure explicitly (it is the μ-marginal of an
axisymmetric `S²` measure) and answer only for groups that *act* on that marginal:

- `Trivial`, `Cn`, `SO2` (any rotation about z fixes μ) → `True`;
- `Z2` (as `σ_h`), `O2`, `Dnh` → `μ → −μ` closure;
- `SO3`, `O3`, `Oh`, `Ih` → **not decidable from the marginal**: either `False`
  (conservative and honest) or raise `NotImplementedError`. Do **not** return
  `True`.

---

## S-3 (HIGH) — `SO3` and `O3` are the same check, and neither is `SO(3)`

`symmetry.py:683-697`.

`SO3` dispatches to `_orbit_closure(..., _icosahedral_ops(), ...)`.
`O3` dispatches to `_icosahedral_ops() + [_inversion_op()]`. But `−I` is
**already** in `_icosahedral_ops()` (the function returns `group + [inv @ M for M
in group]`, and `inv @ I = −I`). So the two op *sets* are identical and

```
  |SO3 ops| = 120   |O3 ops| = 121   |unique(O3 ops)| = 120
  -I already present in _icosahedral_ops(): True
  => SO3.is_invariant(m) == O3.is_invariant(m) for EVERY m.
  BUT the SO3 op set contains 60 IMPROPER (det=-1) matrices, which are NOT in SO(3).
     So the 'SO3' check is really the I_h check: it REJECTS a measure that is
     SO(3)-symmetric-by-sampling but not inversion-symmetric.
```

Two distinct defects: (a) `SO3` and `O3` are indistinguishable although
`SO(3) ⊊ O(3)`; (b) the `SO3` check is *stricter* than `SO(3)` — it demands
closure under 60 improper elements that are not in `SO(3)` at all.

**Fix**: `SO3` must use the proper subgroup `I` (the 60 `det = +1` elements of
`_icosahedral_ops()`); `O3` keeps all 120. The docstring's "necessary condition"
caveat then becomes true rather than merely stated.

---

## S-7 (HIGH) — `D_nh ⊆ O(2)` is asserted TRUE and is false

`symmetry.py:515-518` — `_contains(outer=_NamedSubgroup, inner=Dnh)` returns
`True` for `outer ∈ (O2, O3)`. Pinned by `tests/numerics/test_symmetry.py:216`
(`assert SubgroupOfO3.O2.contains(SubgroupOfO3.Dnh(n))`).

### Numerical demonstration (`p9_monotonicity.py`)

```
=== S-7: the asserted relation  D_nh <= O2  is FALSE ===
  measure: product(n_mu=4, n_phi=4) with all phi offset by 22.5 deg, N = 16
  O2.contains(D4h)  asserted : True
  O2.is_invariant   : True
  D4h.is_invariant  : False
  -> O2-invariant does NOT imply D4h-invariant, so D4h is NOT a subgroup of this O2.

  WHY: the module realises O2 as SO(2)-about-z + sigma_h  (i.e. C_inf_h).
       D_nh additionally contains n VERTICAL mirrors / C_2' axes, which are
       not in C_inf_h.  Under the other common embedding O(2) -> C_inf_v the
       relation fails too (D_nh contains sigma_h, C_inf_v does not).
       D_nh sits inside D_inf_h only.  Direct check of the individual ops:
        closure under sigma_h (z->-z)          : True
        closure under the 4 z-rotations of C_4 : True
        closure under the 4 VERTICAL mirrors   : False  <- the failing generator
```

The module realises `O2` as `_so2_representatives() + _reflections("z")` = the
`C_∞h` embedding. `D_nh ⊄ C_∞h` (vertical mirrors). Under the other standard
embedding `O(2) ↪ O(3)` as `diag(A, 1)` = `C_∞v`, `D_nh ⊄ C_∞v` either (`σ_h`).
`D_nh` sits inside `D_∞h` only.

**Fix**: either drop the `(O2, Dnh)` relation and its test, or mint the group the
relation actually needs (`D_∞h` = `O2 × {e, σ_h}` under the `C_∞v` reading) and
route `D_nh` under it. The module docstring's `O2` gloss
("`SO(2) ⋊ Z_2` — axial rotations + reflection", `:111`) is ambiguous between the
two embeddings and should name which one.

---

## S-16 — the law that ties the two halves together, and its 68 violations

For every asserted edge `A ⊆ B` and every measure `m`, `B.is_invariant(m)` must
imply `A.is_invariant(m)`. Measured across 11 measures × 19 groups
(`p9_monotonicity.py`):

```
=== MONOTONICITY violations:  A <= B asserted, B-invariant TRUE, A-invariant FALSE ===
  68 violation(s) across 11 measures x 19 groups
         O2.contains(C3     ) asserted TRUE, but O2-inv & NOT C3-inv on: ['lebedev(7)', 'lebedev(17)', 'level_symmetric(4)', 'level_symmetric(8)', 'product(4,8)', 'product(4,4)', 'product(4,8)+22.5deg', 'product(4,4)+22.5deg']
        SO2.contains(C3     ) asserted TRUE, but SO2-inv & NOT C3-inv on: ['lebedev(7)', 'lebedev(17)', 'level_symmetric(4)', 'level_symmetric(8)', 'product(4,8)', 'product(4,4)', 'product(4,8)+22.5deg', 'product(4,4)+22.5deg']
         O2.contains(C5     ) asserted TRUE, but O2-inv & NOT C5-inv on: [... same 8 ...]
        SO2.contains(C5     ) asserted TRUE, but SO2-inv & NOT C5-inv on: [... same 8 ...]
         O2.contains(C6     ) asserted TRUE, but O2-inv & NOT C6-inv on: [... same 8 ...]
        SO2.contains(C6     ) asserted TRUE, but SO2-inv & NOT C6-inv on: [... same 8 ...]
         O2.contains(D1h    ) asserted TRUE, but O2-inv & NOT D1h-inv on: ['product(4,4)+22.5deg']
         O2.contains(D2h    ) asserted TRUE, but O2-inv & NOT D2h-inv on: ['product(4,4)+22.5deg']
         O2.contains(D3h    ) asserted TRUE, but O2-inv & NOT D3h-inv on: [... same 8 ...]
         O2.contains(D4h    ) asserted TRUE, but O2-inv & NOT D4h-inv on: ['product(4,4)+22.5deg']
         O2.contains(D6h    ) asserted TRUE, but O2-inv & NOT D6h-inv on: [... same 8 ...]
        SO3.contains(Z2     ) asserted TRUE, but SO3-inv & NOT Z2-inv on: ['asym 1-D [.1,.5,.9]']
```

Root causes: the `C3`/`C5`/`C6` rows are **S-1** (`SO2 ⊇ C_n` is a *true*
relation; the `SO2` *checker* is wrong). The `D_nh` rows are **S-7** (a false
relation). The `Z2` row is **S-2 + S-9**.

This monotonicity law is the single strongest regression gate the module could
carry and there is no test for it today.

---

## S-9 (MOD) — `Z2 ⊆ SO3` asserted TRUE while `Z2`'s realisation is improper

`symmetry.py:204` stores `(Z2, SO3)`. `symmetry.py:659-663` realises `Z2` as
`_reflections("z")`:

```
  lattice asserts Z2 <= SO3 : True
  the 3-D check realises Z2 as _reflections('z') =
[[ 1.  0.  0.]
     [ 0.  1.  0.]
     [ 0.  0. -1.]]
  det = -1.0  -> IMPROPER.  SO(3) contains only det=+1 elements.
  => under the module's OWN realisation of Z2, Z2 is NOT a subgroup of SO(3).
```

The enum comment `Z2 = "Z2"  # {e, σ} — single reflection / 180° rotation`
(`:109`) is the ambiguity: the lattice reads `Z2` as the *proper* `C_2`, the
invariance check realises it as the *improper* `σ_z`. **Fix**: pick one. `Cn(2)`
already exists for the proper reading, so `Z2` should be the reflection group and
the `(Z2, SO3)` edge must go.

---

## S-10 (MOD) — `Cn(1)` and `Trivial` are the same group but not in the lattice

`symmetry.py:486` (`Cn` outer accepts only `Trivial`) and `:511-513` (`Trivial`
outer accepts `Cn(1)`). The two are mutually contained, i.e. equal — but the
named entries disagree (`p8_lattice.py`):

```
=== antisymmetry: contains both ways but NOT equal ===
  violations: [('Trivial', 'C1')]

=== transitivity: A<=B and B<=C but NOT A<=C ===
  3 violation(s):
     C1 <= Trivial <= Z2  but NOT C1 <= Z2
     C1 <= Trivial <= Oh  but NOT C1 <= Oh
     C1 <= Trivial <= Ih  but NOT C1 <= Ih

=== Cn(1) is the trivial group — is the lattice consistent about that? ===
  Trivial.contains(Cn(1)) = True    Cn(1).contains(Trivial) = True    => they are mutually contained, i.e. the SAME group
   Trivial.contains(Trivial) =  True   .contains(Cn(1)) =  True
        Z2.contains(Trivial) =  True   .contains(Cn(1)) = False  <-- INCONSISTENT
       SO2.contains(Trivial) =  True   .contains(Cn(1)) =  True
        O2.contains(Trivial) =  True   .contains(Cn(1)) =  True
        Oh.contains(Trivial) =  True   .contains(Cn(1)) = False  <-- INCONSISTENT
        Ih.contains(Trivial) =  True   .contains(Cn(1)) = False  <-- INCONSISTENT
       SO3.contains(Trivial) =  True   .contains(Cn(1)) =  True
        O3.contains(Trivial) =  True   .contains(Cn(1)) =  True
```

**Fix**: normalise `Cn(1) → Trivial` and `Dnh(1) → ` its named equivalent at
construction (`SubgroupOfO3.Cn` classmethod), so the tag space has no aliases.
That is the illegal-states-unrepresentable spelling and removes all four
violations at once.

The named table itself is fine — it *is* transitively closed:

```
=== the named table is transitively closed? ===
  stored edges: 19   transitive closure: 19
  edges implied by transitivity but MISSING from the table: none
```

Reflexivity holds for all 19 groups probed.

---

## S-8 (MOD) — the `contains` docstring's `O_h`/`I_h` claim is false and self-contradictory

`symmetry.py:359-364` states:

> `C_n` and `D_nh` both sit below `O(3)`; **neither sits inside** the named
> compact subgroups `O_h`, `I_h` (those have specific finite cyclic subgroups —
> `C_4 ⊂ O_h`, `C_5 ⊂ I_h`, etc. — but exhaustive enumeration … is out of scope)

The first clause contradicts the parenthetical in the same sentence. Computing
the truth **axis-fixed**, from the module's own operator sets (`p10_axis_truth.py`):

```
Cn and Dnh are realised ABOUT THE Z AXIS (_rotation_z / _vertical_mirrors).
Truth is therefore computed axis-fixed, directly from the module's own O_h / I_h sets.

  group    in O_h (z-fixed)   contains() says    in I_h (z-fixed)   contains() says
     C1                True             False                True             False
     C2                True             False                True             False
     C3               False             False               False             False
     C4                True             False               False             False
     C5               False             False               False             False
     C6               False             False               False             False
    D1h                True             False                True             False
    D2h                True             False                True             False
    D3h               False             False               False             False
    D4h                True             False               False             False
    D6h               False             False               False             False

=> AXIS-FIXED, the true relations MISSING from the lattice are:
      C1 <= Oh      C1 <= Ih      C2 <= Oh      C2 <= Ih      C4 <= Oh
      D1h <= Oh     D1h <= Ih     D2h <= Oh     D2h <= Ih     D4h <= Oh
```

These are **false negatives** (conservative — they reject a legal pairing rather
than admit an illegal one), so they are not a correctness hazard for the
selection layer today. They are a documentation defect plus 10 missing relations.

**The deeper issue the docstring never resolves**: is containment *up to
conjugation* or *axis-fixed*? The two answers differ:

```
   rot_z(2pi/5) in _icosahedral_ops() : False
   (I_h in this construction has its 5-fold axes along the icosahedron
    vertices, none of which is +z.)
```

`C_5 ⊂ I_h` is TRUE abstractly and FALSE for the z-axis realisation the module
actually uses. `contains` and `is_invariant` must agree on which; today the
docstring speaks abstractly and the code acts axis-fixed. **Fix**: state the
convention in the class docstring, then fill in the 10 axis-fixed relations
(they are computable in one line from the operator sets, so "exhaustive
enumeration is out of scope" is no longer a reason).

---

## S-6 (MOD) — bucket-hit and brute-force paths select different `j`

`symmetry.py:938-951`. On a bucket **miss** the code takes the global
`argmin`; on a bucket **hit** it takes the min over the bucket's members only.
When the target sits near a quantisation boundary, the true nearest node can fall
into the adjacent bucket while a farther decoy stays inside the target's bucket.

### Numerical demonstration (`p6_paths_tol.py`)

```
=== S-6: bucket-HIT selects a FARTHER j than brute-force would ===
  target z  = 0.5000000000000491   bucket 5000000000000
  true  z   = 0.500000000000051   bucket 5000000000001   |dz| = 0.019 * atol
  decoy z   = 0.49999999999996   bucket 5000000000000   |dz| = 0.891 * atol
  -> decoy shares the target's bucket (True); true partner does NOT (False)
  -> bucket path returns the DECOY (0.89*atol away); brute force would return the TRUE partner (0.02*atol away)
  nodes   =
 [[ 0.5                 0.25               -0.5000000000000491]
    [ 0.5                 0.25                0.500000000000051 ]
    [ 0.5                 0.25                0.49999999999996  ]]
  weights = [1. 1. 9.]
  _orbit_closure(sigma_z) = False
  POSITIVE CONTROL (decoy removed): True
  => the decoy, though 44x farther, decides the weight comparison and the verdict.
```

This answers the brief's question directly: **yes**, the two paths can select
different `j`, and the consequence here is a **false negative**. The precondition
is near-duplicate nodes (< `atol` apart), so it is not a live hazard for the four
shipped rules — but it *is* a correctness hazard for the carve, because the
permutation returned would be wrong.

Related: the dict index is **not** inert on the true cases (contrary to what a
quick read suggests). Branch statistics:

```
=== bucket-HIT vs BRUTE-FORCE branch statistics on the shipped rules ===
  lebedev(order=7)  vs O_h     bucket-hit    1248  brute-force       0  (  0.0% fall through)
  lebedev(order=17) vs O_h     bucket-hit    5280  brute-force       0  (  0.0% fall through)
  level_symmetric(4) vs O_h    bucket-hit    1152  brute-force       0  (  0.0% fall through)
  level_symmetric(8) vs O_h    bucket-hit    3840  brute-force       0  (  0.0% fall through)
  product(4,8) vs C_4-as-SO2   bucket-hit     512  brute-force    1024  ( 66.7% fall through)
```

For genuinely-invariant rules the index hits 100 %; the fall-throughs in the last
row are genuine non-membership (the image node really is absent), where the
brute-force scan then correctly rejects. **No performance defect.**

**Fix**: probe the 3³ = 27 neighbouring buckets on a hit, or (simpler and exactly
as fast at these sizes) drop the dict and use `scipy.spatial.cKDTree` /
`np.argmin`, which has no boundary pathology at all.

---

## S-11 (MOD) — the 100× node / 1× weight tolerance asymmetry is accidental

`symmetry.py:945` and `:950` use `atol * 100` for the node match; `:952` uses
bare `atol` for the weight match. Both are **absolute** and neither scales with
the data.

```
=== tolerance asymmetry: node window = 100*atol, weight window = 1*atol ===
  level_symmetric(4): weights are all 5.235988e-01
  lebedev(17): weight range [4.810747e-02, 1.249451e-01]
  relative weight resolution of the check at default atol=1e-13:
     level_symmetric(4): 1.910e-13   lebedev(17): 2.079e-12

  node[0] displaced by 5e-12 (< 100*atol): is_invariant(O_h) = True
  weight[0] scaled by (1+1e-9) => abs delta 4.811e-11: is_invariant(O_h) = False
```

A **5e-12 geometric error** (≈ 2×10⁴ ULP on the unit sphere — far beyond any
legitimate FP noise) is accepted, while a **1e-9 relative weight error** is
rejected. There is no principle that makes 100× the right ratio; `_orbit_closure`
uses 100×, `_is_reflection_invariant_1d:644` uses **10×** for the same kind of
node match, and neither is documented. The `atol` docstring (`:447-450`) claims it
"matches the floating-point noise on a 1-D Gauss-Legendre weight computation" —
which describes the *weight* window only and says nothing about the node windows
that are 10×/100× larger.

**Fix**: one named node tolerance and one named weight tolerance, both
parameters, weight relative (`rtol`) rather than absolute, and the constants
documented. `np.isclose`-style `atol + rtol*|w|` is the standard spelling.

Mutation-verified as ungated (in-process plugin widening the weight window
1× → 100×): `[MUT-WTOL] mutated _orbit_closure called 9 time(s)` and
`182 passed`. Positive control that the mutation bites:

```
=== MUT-WTOL bite check ===
  lebedev(17) with weights[0] += 5e-12 (between 1*atol and 100*atol)
    production (weight window = 1*atol)   : False
    MUTATED    (weight window = 100*atol) : True
  => the mutation DOES change the answer; the green suite is blindness.
```

---

## S-4 (MOD) — `_is_reflection_invariant_1d` is also non-bijective

`symmetry.py:619-648`. The logic **does** match its docstring ("closure of a 1-D
measure under `x → −x`") in shape, but carries the same missing-bijection defect
as S-5, plus a third distinct node tolerance (`atol * 10`, `:644`).

```
=== S-4  _is_reflection_invariant_1d: non-bijective match ===
  nodes   = [-0.5  0.5  0.5]   (3 nodes; -0.5 has TWO near-partners at +0.5)
  weights = [1. 1. 1.]
  _is_reflection_invariant_1d -> True
  TRUTH: the reflected set is {+0.5, -0.5, -0.5-3e-13}, which is NOT the input set
         (the input has two nodes near +0.5 and one near -0.5; the image has the
          opposite multiplicity).  A bijection does not exist.  Reported True.
```

(The third node is `0.5 + 3e-13`, inside the `10·atol = 1e-12` window; numpy
prints both as `0.5`.)

`np.searchsorted` index handling (`:634-643`) is otherwise correct — `j == n` and
`j == 0` are both handled, and `candidates` cannot be empty for `n ≥ 1`.

---

## S-12 (LOW) — the `is_invariant` docstring describes strategies the code does not use

This answers the addendum directly. **Two separate results.**

### (a) The docstring is wrong (doc/code divergence, in the SAFE direction)

`symmetry.py:428-431` claims the `OctahedralOh` path decides by *"the multiset of
`(|x|, |y|, |z|)` triples sorted, weighted"*, and `:432-438` claims the
`IcosahedralIh` path uses *"a 12-element representative orbit"*. Neither is what
runs (`p2_strategy_table.py`):

```
=== docstring claim vs code: OctahedralOh ===
docstring (symmetry.py:428-431): 'The fingerprint is the multiset of (|x|,|y|,|z|) triples sorted, weighted.'
code actually calls _orbit_closure with n_ops = [48] (= |O_h| = 48)
=> the docstring describes a strategy the code does NOT implement.

=== docstring claim vs code: IcosahedralIh ===
docstring (symmetry.py:432-438): 'checked via a 12-element representative orbit on a regular icosahedron vertex set'
code actually calls _orbit_closure with n_ops = [120]
```

The module docstring (`:44-54`) likewise calls the whole strategy a "fingerprint"
and describes "nearest-neighbour matching" — the second half is accurate, the
"fingerprint" framing is not.

### (b) The fingerprint the docstring names is provably **necessary-only**

Had the code implemented it, it would be the severe defect the addendum
suspected. Demonstrated with two fingerprint-preserving, closure-breaking
mutations (`p3_oh_fingerprint.py`):

```
base rule: level_symmetric(sn_order=4), N = 24
  is_invariant(O_h)             : True
  direct 48-op _orbit_closure   : True

MUTATION: flip sign of x on node 0 ONLY (|x|,|y|,|z| and w unchanged)
  node 0 before: [-0.40824829 -0.81649658 -0.40824829]  after: [ 0.40824829 -0.81649658 -0.40824829]
  positive control  np.array_equal(base, mutated) = False  (must be False)
  positive control  #rows differing              = 1
  fingerprint(base) == fingerprint(mutated)      = True  <- the DOCSTRING'S fingerprint cannot see the mutation
  is_invariant(O_h) on mutated                   = False  <- the CODE rejects it
  direct 48-op _orbit_closure on mutated         = False

MUTATION 2: randomise every node's octant (all |x|,|y|,|z| and w preserved)
  fingerprint(base) == fingerprint(mut2)         = True
  #rows differing from base                      = 21
  is_invariant(O_h) on mut2                      = False
```

`I_h` positive + negative controls (`p3`):

```
=== I_h ===
  12 icosahedron vertices, equal w:  is_invariant(I_h) = True
                                     is_invariant(O_h) = False
  POSITIVE CONTROL w[0] *= 1.001:    is_invariant(I_h) = False (must be False)
  POSITIVE CONTROL node[0] += 1e-6:  is_invariant(I_h) = False (must be False)
```

**Verdict**: the `O_h` and `I_h` paths are **sound** (modulo S-5's missing
bijection, which applies to every path). The defect is the docstring. Fix it —
otherwise a future reader "optimises" the code to match the documented
fingerprint and lands the severe bug.

---

## THE STRATEGY TABLE (load-bearing for the carve)

Measured by wrapping `_orbit_closure` and recording the op count per call
(`p2_strategy_table.py`).

```
STRATEGY TABLE — n_ops passed to _orbit_closure (S^2 measure: lebedev order=7, N=26)
             tag                           strategy   n_ops  verdict
         Trivial  SHORT-CIRCUIT (no _orbit_closure)       -     True
              Z2                     _orbit_closure     [1]     True
             SO2                     _orbit_closure     [4]     True
              O2                     _orbit_closure     [5]     True
    OctahedralOh                     _orbit_closure    [48]     True
   IcosahedralIh                     _orbit_closure   [120]    False
             SO3                     _orbit_closure   [120]    False
              O3                     _orbit_closure   [121]    False
           Cn(4)                     _orbit_closure     [4]     True
          Dnh(4)                     _orbit_closure     [9]     True

STRATEGY TABLE — 1-D measure (gauss_legendre_on_mu(8), support='[-1,1]', dim=1)
             tag                           strategy   n_ops  verdict
         Trivial       1-D path (NO _orbit_closure)       -     True
              Z2       1-D path (NO _orbit_closure)       -     True
             SO2       1-D path (NO _orbit_closure)       -     True
              O2       1-D path (NO _orbit_closure)       -     True
    OctahedralOh       1-D path (NO _orbit_closure)       -     True
   IcosahedralIh       1-D path (NO _orbit_closure)       -     True
             SO3       1-D path (NO _orbit_closure)       -     True
              O3       1-D path (NO _orbit_closure)       -     True
           Cn(4)       1-D path (NO _orbit_closure)       -     True
          Dnh(4)       1-D path (NO _orbit_closure)       -     True
```

**For the carve**: on an `S²` (or any `d ≥ 2`) measure, **every** non-`Trivial`
tag goes through `_orbit_closure`, so a permutation is available for all of them.
The tags with **no** permutation available are:

1. `Trivial` (and `Cn(1)`) on any measure — short-circuits at `:549-552`. The
   permutation is trivially the identity; return it explicitly.
2. **Every tag on a 1-D measure** — `_check_invariance_1d` never calls
   `_orbit_closure`. `Z2`/`O2`/`O3`/`Oh`/`Ih`/`Dnh` route to
   `_is_reflection_invariant_1d`, which *computes* a match (`best`) but discards
   it; `SO2`/`SO3`/`Cn` return `True` with no computation at all.

The clean carve therefore has two jobs, not one: widen `_orbit_closure`, **and**
give the 1-D path a permutation-returning kernel (which for the reflection case
is a 5-line change to `_is_reflection_invariant_1d`, since it already computes
`best`). If the S-2 fix collapses the 1-D path onto a 1-D `_orbit_closure` over
`{+1, −1}` matrices, both jobs merge into one — the preferable shape.

---

## S-13 (LOW) — stringly-typed `support` dispatch overrides the geometric `dim`

`symmetry.py:558-562`: `is_1d = measure.dim == 1 or measure.support in
(SPACE_INTERVAL_M11, "[0,1]", "R")`. The `or` means the **string tag alone**
routes a genuinely 3-D node array to the 1-D branch, where `:621` silently
reduces `(N, 3)` nodes to their **x-column**.

```
lebedev(7): nodes.shape = (26, 3)  dim = 3  support = S^2

SAME 3-D nodes re-tagged support='[-1,1]':  dim = 3
  is_invariant(Oh) on the correctly-tagged measure : True
  is_invariant(Oh) on the mis-tagged measure       : False
  -> the support STRING alone flips the dispatch to the 1-D branch, which then
     silently reduces the (N,3) node array to its x-COLUMN (symmetry.py:621)
     and answers from that.  No shape/dim guard.
```

Also: `"[0,1]"` is spelled as a literal although `SPACE_INTERVAL_01` exists
(`measure.py:105`), so a rename of the constant would silently drop the branch.
**Fix**: dispatch on `measure.dim` alone (geometry is a property of the array,
not of a label), or guard with `if measure.dim != 1: raise`.

---

## S-14 (LOW) — public `atol` is unvalidated; `≤1e-19` silently saturates `int64`

```
=== int64 quantisation safety of np.round(nodes/atol).astype(np.int64) ===
  int64 max = 9.223372e+18
  atol=1e-13:  1.0/atol = 1.000e+13   quantised int =        10000000000000   round-trips? True
  atol=1e-15:  1.0/atol = 1.000e+15   quantised int =      1000000000000000   round-trips? True
  atol=1e-16:  1.0/atol = 1.000e+16   quantised int =     10000000000000000   round-trips? True
  atol=1e-18:  1.0/atol = 1.000e+18   quantised int =    999999999999999872   round-trips? True
  atol=1e-19:  1.0/atol = 1.000e+19   quantised int =   9223372036854775807   round-trips? False
  atol=1e-20:  1.0/atol = 1.000e+20   quantised int =   9223372036854775807   round-trips? False
```

Answering the brief's specific question: **at the default `atol = 1e-13` with
coordinates in `[-1, 1]`, the `int64` conversion is safe** (`1e13 ≪ 2⁵³`, so the
quantised integer is exactly representable and the bucket key is deterministic).
It stays safe down to `atol = 1e-18`. At `1e-19` every coordinate saturates to
`int64` max — **all nodes collapse into one bucket** with no exception and no
warning. `atol` is a public keyword with no validation (`:396-401`). **Fix**:
validate `atol > 0` and reject/clamp below the safe floor, or key the bucket on
the rounded `float` directly and skip the `int64` round-trip entirely.

---

## S-15 (LOW) — dead code

`_rotation_x` (`:732`) and `_rotation_y` (`:741`) have **zero call sites** in the
repository (`grep -rn "_rotation_x\|_rotation_y" --include=*.py --include=*.rst .`
returns only their definitions and the worktree copy). Both are correct
(verified below); they are simply unused. Per the retire-as-you-go rule, delete
them or wire them into the `S-3` fix (a proper-`SO(3)` sample would naturally use
off-z axes).

---

## CLEAN — the operator generators

All eight generators were probed for orthogonality, determinant, uniqueness,
generated-group order, and closure under multiplication (`p7_generators.py`).
**No defects found.**

```
generator                       n  uniq    orth err                  dets  closed?  group order
_reflections('x')               1     1    0.00e+00                [-1.0]    False       2 (exp 2)
_reflections('y')               1     1    0.00e+00                [-1.0]    False       2 (exp 2)
_reflections('z')               1     1    0.00e+00                [-1.0]    False       2 (exp 2)
[_inversion_op()]               1     1    0.00e+00                [-1.0]    False       2 (exp 2)
_so2_representatives()          4     4    0.00e+00                 [1.0]     True       4 (exp 4)
_cyclic_ops(1)                  1     1    0.00e+00                 [1.0]     True       1 (exp 1)
_cyclic_ops(2)                  2     2    0.00e+00                 [1.0]     True       2 (exp 2)
_cyclic_ops(3)                  3     3    0.00e+00                 [1.0]     True       3 (exp 3)
_cyclic_ops(4)                  4     4    0.00e+00                 [1.0]     True       4 (exp 4)
_cyclic_ops(6)                  6     6    0.00e+00                 [1.0]     True       6 (exp 6)
_vertical_mirrors(1)            1     1    0.00e+00                [-1.0]    False       2 (exp 2)
_vertical_mirrors(2)            2     2    0.00e+00                [-1.0]    False       4 (exp 4)
_vertical_mirrors(3)            3     3    2.22e-16                [-1.0]    False       6 (exp 6)
_vertical_mirrors(4)            4     4    0.00e+00                [-1.0]    False       8 (exp 8)
_octahedral_ops()              48    48    0.00e+00           [-1.0, 1.0]     True      48 (exp 48)
_icosahedral_ops()            120   120    1.11e-15           [-1.0, 1.0]     True     120 (exp 120)
```

(`closed? = False` for the single-matrix generators is expected — a lone
reflection is a *generator*, not a group; the "group order" column confirms each
generates the right group. `_orbit_closure` checking generators rather than full
groups is **mathematically sound**: closure under each generator implies closure
under every product, and weight preservation composes. This is precisely why S-1
is a defect — the `SO2` generator set generates `C_4`, not `SO(2)`.)

### `_octahedral_ops()` — CLEAN BILL

```
  count               : 48  unique: 48
  all orthogonal      : True
  determinants        : [-1.0, 1.0]   #(+1): 24  #(-1): 24
  entries in {-1,0,1} : True
  closed under mult   : True
  contains identity   : True
  contains inversion  : True
```

Exactly the 48 signed permutation matrices = the full hyperoctahedral group
`B_3 ≅ O_h`, with the correct 24 proper / 24 improper split. The docstring's
count claim (`8 × 6 = 48`) is accurate.

### `_icosahedral_ops()` — CLEAN BILL

```
  returned count      : 120  unique: 120
  #proper (det +1)    : 60  (I has 60)
  #improper (det -1)  : 60  (I_h \ I has 60)
  orthogonality err   : 1.11e-15
  closed under mult   : True
  element orders      : {1: 1, 2: 31, 3: 20, 5: 24, 6: 20, 10: 24}
```

The order census matches `I_h` exactly: 1 identity; 31 order-2 (15 `C_2` + 15 `σ`
+ inversion); 20 `C_3`; 24 `C_5`; 20 `S_6`; 24 `S_10`. Total 120. The BFS
(`iteration_cap = 200`, dedup at `atol = 1e-10`) terminates correctly and the
chosen face `(raw[0] + raw[2] + raw[8])/3` is a genuine face (all three pairwise
edge lengths equal).

### `_rotation_x/y/z` — CLEAN BILL

```
  _rotation_x(pi/2): e1 -> [0. 0. 1.]   right-handed? True   det = +1.000000000000000
  _rotation_y(pi/2): e0 -> [ 0.  0. -1.]                     det = +1.000000000000000
  _rotation_z(pi/2): e0 -> [0. 1. 0.]   right-handed? True   det = +1.000000000000000
```

`_rotation_y` maps `e_x → −e_z` and `e_z → +e_x`, which **is** the standard
right-handed `R_y = [[c,0,s],[0,1,0],[−s,0,c]]` (the `y`-rotation cycle is
`z → x → −z`, not `x → z`). The probe's initial "right-handed? False" was a
probe-expectation error, not a code defect — verified by `R_y(π/2) · e_z = e_x`.

### `_rotation_about_axis` — CLEAN (Rodrigues, `orth err ≤ 1.11e-15` via `_icosahedral_ops`)
### `_vertical_mirrors(n)` — CLEAN (n distinct mirrors, normals at `kπ/n`, generating `C_nv` of order `2n`)

---

## Existing test coverage

`tests/numerics/test_symmetry.py` — 71 tests, all `@pytest.mark.foundation`.
Baseline: `182 passed` across the five numerics files.

### GATED today

| Behaviour | Test |
|---|---|
| Reflexivity of the 8 named entries | `test_named_reflexivity` |
| `Trivial ⊆` every named entry | `test_trivial_inside_every_named_group` |
| `O(3) ⊇` every named entry | `test_o3_contains_every_proper_subgroup` |
| `SO2 ⊂ O2 ⊂ O3` | `test_so2_chain` |
| `Z2 ⊂ Oh ⊂ O3` | `test_z2_chain_via_octahedral` |
| `Oh ⊄ SO3`, `Ih ⊄ SO3` | `test_octahedral_not_in_so3`, `test_icosahedral_not_in_so3` |
| `C_n ⊆ C_m ⟺ n∣m` | `test_cyclic_containment_by_divisibility` |
| `C_n ⊆ SO2`, `C_n ⊆ D_nh`, `Z2 ⊆ D_nh` | 3 tests |
| `is_subgroup_of` reverses `contains` | `test_is_subgroup_of_reverses_contains` |
| `O_h` invariance of Lebedev (11/17/23) and LS (2/4/6/8) | 7 params |
| `I_h` **negative** on Lebedev-17 | `test_lebedev_is_NOT_icosahedral_invariant` |
| `O_h` **negative** on a 2-point asymmetric set | `test_asymmetric_measure_not_octahedral_invariant` |
| `Trivial` short-circuit | `test_trivial_invariance_is_universal` |
| 1-D GL `SO2`/`Z2` (both via the 1-D short-circuit) | 6 params |

### UNGATED

- **`SO2` / `O2` invariance on any `S²` measure** — the S-1 defect. `_so2_representatives` is called **0 times** in the whole 182-test run.
- **`Cn` and `Dnh` invariance** — no `is_invariant` test at all (containment only).
- **`SO3` / `O3` invariance** — no test at all; S-3 is entirely unobserved.
- **`I_h` positive** — only the Lebedev negative exists; nothing pins that `I_h` *accepts* an icosahedral set.
- **Bijectivity of the orbit match** — S-5; no test builds a duplicated-node measure.
- **Monotonicity** (`A ⊆ B` ∧ `B`-invariant ⟹ `A`-invariant) — the single law tying the two halves together; 68 violations, 0 tests.
- **Transitivity / antisymmetry of `contains`** over the parameterised families — S-10; 4 violations, 0 tests.
- **The tolerance constants** — S-11; both mutations pass green.
- **Non-`O_h` groups' operator sets** — see the `MUT-OH` result below.

### One committed test pins a FALSE relation

`tests/numerics/test_symmetry.py:216`
```python
assert SubgroupOfO3.O2.contains(SubgroupOfO3.Dnh(n))
```
This is S-7. Fixing the lattice requires deleting this assertion.

### Mutation results (Enforcement #11 — every gate cited was made RED or proven blind)

In-process pytest plugins under `PYTHONPATH`; **no tracked file was edited**.

```
### BASELINE (no mutation)
182 passed, 1 warning in 1.78s
### MUTATION: mut_so2      (SO(2) sample widened C_4 -> C_8)
[MUT-SO2] _so2_representatives called 0 time(s) during the run
182 passed, 1 warning in 0.96s
### MUTATION: mut_oh       (O_h crippled to its 8 diagonal sign-flips = D_2h)
[MUT-OH] _octahedral_ops called 8 time(s) during the run
182 passed, 1 warning in 1.66s
### MUTATION: mut_wtol     (weight window 1*atol -> 100*atol)
[MUT-WTOL] mutated _orbit_closure called 9 time(s) during the run
182 passed, 1 warning in 1.32s
```

`MUT-OH` is the strongest coverage finding: replacing the 48-element `O_h` with
its 8-element diagonal subgroup leaves the whole suite green, so **no test can
distinguish `O_h` from `D_2h`** — the coordinate-permutation half of `O_h` is
unconstrained. Positive control that the mutation genuinely changes an answer:

```
=== MUT-OH bite check ===
  measure: all sign combinations of (0.6, 0.8, 0)/|.|, N = 4
  closure under the 8 SIGN FLIPS only (the MUTATED O_h) : True
  closure under the real 48-element O_h                 : False
  is_invariant(O_h) (production)                        : False
  => the mutation DOES change the answer (True vs False) on this input,
     so the green 182-test suite under MUT-OH is genuine blindness.
     x<->y swap of (0.6,0.8,0) gives (0.8,0.6,0), absent from the set:
       min distance from (0.8,0.6,0) to the node set = 0.282843
```

(`MUT-SO2`'s 0-call count is **not** an inert mutation — it is the finding: the
mutated code is never reached. Verified separately in `p11_bite.py`, where an
`S²` measure does reach it, 1 call.)

---

## Is `is_invariant`'s public contract testable in its current `-> bool` shape?

**No, not adequately — and this is exactly what the carve fixes.**

A `bool` collapses three distinct claims into one bit:

1. **which** group elements were checked (the S-1 defect is invisible: a `True`
   from a `C_4` sample is indistinguishable from a `True` from the real `SO(2)`);
2. **whether the match was a bijection** (the S-5 defect is invisible: `True` from
   a genuine permutation is indistinguishable from `True` from a many-to-one map);
3. **which strategy decided it** (the 1-D-vs-3-D-vs-short-circuit split of S-2 is
   invisible: the caller cannot tell that `Oh.is_invariant(gl_1d)` answered from a
   μ-reflection rather than from 48 operators).

Everything a test can assert against `-> bool` is a *single sample of the truth
table*, which is why the suite's 7 `is_invariant` calls leave S-1, S-2, S-3 and
S-5 wholly unobserved, and why crippling `O_h` to `D_2h` is green.

Returning the **orbit permutation** (the carve's stated goal) makes all three
testable directly:

- the permutation's existence and **bijectivity** become assertable
  (`sorted(perm) == range(n)`) — closing S-5 by construction rather than by an
  extra check;
- `perm` per operator exposes **which operators ran** (`len(perms) == 48`)
  — closing the `MUT-OH` blindness and making S-1 assertable
  (`len(perms) == 4` for `SO2` is immediately visibly wrong);
- the permutation is the object the **orbifold singular set** needs anyway (a
  node is singular iff it is a **fixed point** of some non-identity element, i.e.
  `perm[i] == i` for some `M ≠ I`) — so the two carve halves are one change, and
  the singular set is a *derived* quantity rather than a second code path.

Recommended return shape: `dict[int, np.ndarray]` (op index → permutation array)
or a small frozen `OrbitAction` dataclass carrying `ops`, `perms`, and a
`stabiliser` accessor; keep `is_invariant -> bool` as a thin
`... is not None` façade so existing call sites are untouched.

---

## Recommended fix order (for the carve)

1. **S-5** (bijection) — free for all shipped rules, and the carve's return value
   is wrong without it. Do it *inside* the widening.
2. **S-1** (`SO2`/`O2`) — the only defect that lets a non-invariant rule certify
   through a *sound-looking* path. Requires deciding the discrete-vs-continuum
   semantics of the `SO2` tag, which also touches the four `invariance_group=SO2`
   declarations.
3. **S-2 / S-3** (1-D path, `SO3` vs `O3`) — same decision as (2), applied to the
   marginal; collapsing the 1-D path onto `_orbit_closure` over `{±1}` merges the
   carve's two jobs.
4. **S-7 / S-9 / S-10** (lattice) — on the live selection path; each is a
   few-line change plus one test edit (`test_symmetry.py:216`).
5. **S-11 / S-6 / S-4 / S-13 / S-14** — tolerance and dispatch hygiene; naturally
   folded into the rewrite.
6. **S-8 / S-12 / S-15** — documentation truth + dead-code retirement.
7. **New gate**: the monotonicity law of S-16, parametrised over every
   (edge, measure) pair. It is one loop, it catches S-1/S-2/S-7/S-9 simultaneously,
   and it is the natural home for the `catches(...)` markers once these are
   catalogued.
