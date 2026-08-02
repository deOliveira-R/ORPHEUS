# #325 — main-agent measurements

Everything here is `[M]` measured by the main agent, not inherited from the issue
body. Two of my own round-1 probes were WRONG and are recorded as such, because
the wrong ones are the more instructive half.

Scripts: `$CLAUDE_JOB_DIR/tmp/verify_circle_nodes{,2}.py`, `verify_argsort_ties.py`.

---

## 0. Two probes of mine that were WRONG — do not inherit their conclusions

| probe | what I wrote | why it measured nothing |
|---|---|---|
| the 45° diagonal mask | `(8p) % q == np.where(2*((4p)%q)==q, 1, 0)` | garbage algebra — it asks whether `(8p)%q` equals 0 or 1. The correct condition is `2*((4p) mod q) == q`. |
| the x-mirror on the status quo | `np.cos(-phi)` vs `np.cos(phi)` | **bit-equal for ANY `phi`** — libm's cosine is even by construction. It tests the library, not the node set. |

The x-mirror one produced a confident *false negative* ("the legacy mirror is
already fine"), which would have killed the issue's premise. The honest test
applies the mirror as an **INDEX PERMUTATION on the generated set**
(`k → (n−k) mod n`), because that is what production consumes.

**Rule this earns:** when testing a symmetry of a node SET, permute the set's
indices. Never negate the parameter — the parameterizing function's own parity
will satisfy you regardless of what the set does.

---

## 1. The primitive is exact where the status quo is not `[M]`

`roots_of_unity` vs `np.linspace(0, 2π, n, endpoint=False)` + `np.cos`/`np.sin`,
each mirror applied as an index permutation:

```
     n |  NEW x  OLD x |  NEW y  OLD y | NEW #exact0 OLD #exact0
     3 |   True  False |   None   None |           1           1
     4 |   True  False |   True  False |           4           1
     6 |   True  False |   True  False |           2           1
     8 |   True  False |   True  False |           4           1
    12 |   True  False |   True  False |           4           1
    16 |   True  False |   True  False |           4           1
    24 |   True  False |   True  False |           4           1
    32 |   True  False |   True  False |           4           1
    48 |   True  False |   True  False |           4           1
    64 |   True  False |   True  False |           4           1
```

x-mirror: exact for **every** `n` tested {3,4,5,6,8,10,12,16,20,24,32,48,64};
legacy exact for **none**. y-mirror likewise on every even `n`.

45° diagonal (`2·((4p) mod q) == q`), `|cos| == |sin|` bitwise:

```
  q=  8: 4 diagonal nodes | NEW: True | LEGACY: False
  q= 16: 4 diagonal nodes | NEW: True | LEGACY: False
  q= 24: 4 diagonal nodes | NEW: True | LEGACY: False
  q= 32: 4 diagonal nodes | NEW: True | LEGACY: False
  q= 48: 4 diagonal nodes | NEW: True | LEGACY: False
  q= 64: 4 diagonal nodes | NEW: True | LEGACY: False
```

`np.cos(pi/4) - np.sin(pi/4) = 1.110e-16` — the 1-ULP gap the diagonal override
closes. Without the override the construction would break its own swap symmetry
at precisely the fixed point that defines it.

## 2. It is also strictly MORE ACCURATE — the criterion-(2) reference `[M]`

Against a **100-digit mpmath** reference (structurally independent of both
implementations), per-component, over the same `q` sweep:

```
  TOTAL components where NEW is worse: 8;  where OLD is worse: 312
```

So this is not a lateral move justified only by structure — the values become
more correct. That is `vv-principles` "Bit-identity vs principled-equivalence"
criterion (2) satisfied by a genuine independent ground, not by an old-vs-new
ULP distance.

Node drift new-vs-old, i.e. the re-baseline magnitude:

```
  n_phi=  3: 4.996e-16     n_phi= 16: 6.106e-16
  n_phi=  4: 1.837e-16     n_phi= 20: 2.776e-16
  n_phi=  5: 1.665e-16     n_phi= 24: 8.327e-16
  n_phi=  6: 6.106e-16     n_phi= 32: 6.106e-16
  n_phi=  8: 2.220e-16     n_phi= 48: 9.992e-16   <- worst
  n_phi= 12: 6.106e-16     n_phi= 64: 6.106e-16
```

⚠ **Correction to the issue body**: it claims max drift `8.3e-16 at n=24`. The
true worst over the wider sweep is **9.992e-16 at n_phi=48**. Cite this number.

Unit norm `max |c²+s²−1| = 2.220e-16` (1 ULP, unchanged). Negative zeros after
canonicalization: **0**.

---

## 3. ⚠⚠ THE FINDING THE ISSUE BODY MISSED — exact symmetry breaks a production sort

`rules_product.py:139` sorts each μ-level by direction cosine:

```python
order = np.argsort(mu_x[level_arr])
```

Today the azimuthal mirror pair `(k, n−k)` has `cos φ` values differing in the
last bits, so `mu_x` has a strict total order. Under #325 they become
**BIT-EQUAL** — a genuine tie — and the sort's resolution among tied entries is
no longer the same. `n_phi=8`, one μ level (idx 1,3,5,7 are all nominally
±√2/2):

```
    idx  OLD mu_x                     NEW mu_x
    1   7.07106781186547572737e-01   7.07106781186547572737e-01
    3  -7.07106781186547461715e-01  -7.07106781186547572737e-01
    5  -7.07106781186547683760e-01  -7.07106781186547572737e-01
    7   7.07106781186547350693e-01   7.07106781186547572737e-01
```

Four distinct floats collapse to two. Measured over
`n_mu ∈ {2,4,8} × n_phi ∈ {3,4,6,8,12,16,24,32}`:
**24 of 24 configurations have their per-level ordinate order permuted.**

Why it is not cosmetic: `LevelStructure.level_indices` is consumed in production
at `orpheus/transport/radial_characteristic_field.py:289` (the curvilinear pole
angular closure), and the Bailey Eq. 50 α-recursion is a **sequential running
sum over the level** — a permutation changes the recurrence order, which is not
a ≤1-ULP perturbation of anything.

This is crosswalk §16.5 ("⭐⭐ THE RECURRING HAZARD — index-order coincidences on
symmetric quadratures") firing again, one phase later.

**Two questions it raises, both out for measurement:**

1. Is the α-recursion actually order-sensitive end-to-end? (If yes, #325 is
   reclassified from "≤1 ULP re-baseline" to "changes a production recurrence
   order" — a different conversation.)
2. **Is the sort KEY itself defective?** `η = sinθ·cos φ` is **2-to-1** over
   φ ∈ [0, 2π) — φ and 2π−φ share η. So sorting a full azimuthal circle by η
   never had a total order; the old code was ordering the two half-circles by
   rounding noise. If the sweep needs ordering by azimuthal ANGLE, the key is
   wrong TODAY and #325 merely makes the wrongness visible.

A gate asserting "the level is sorted by η" is **Mode-12 blind** to all of
this — sortedness is invariant under permuting equal elements.

---

## 3b. The ordering change is TIGHTLY BOUNDED — and it collapses to ONE decision `[M]`

The reordering is **confined to permuting entries with BIT-EQUAL η and w**
(azimuthal weights within a level are all equal by construction). Measured over
`(n_mu, n_phi) ∈ {(2,4),(4,8),(4,12),(8,16),(6,6),(8,32),(3,5),(4,24)}`:

```
 n_mu  n_phi | eta seq bit-identical?  w*eta partial sums identical?
    2      4 |                   True                           True
    4      8 |                   True                           True
    4     12 |                   True                           True
    8     16 |                   True                           True
    6      6 |                   True                           True
    8     32 |                   True                           True
    3      5 |                   True                           True
    4     24 |                   True                           True
```

So the **α-recursion's running sum is bit-identical** old vs new — the Bailey
Eq. 50 partial sums at each POSITION do not move at all. What changes is only
*which ordinate index sits at each position*, and the swapped ones carry
identical η and w.

That leaves exactly two order-sensitive consumers, and only one of them is real:

**(i) The per-level moment integration** — `radial_characteristic_field.py:300`,
`np.einsum("n,nl,ngx->lgx", w_p, legendre, q_p)`. Same SET of summands, different
enumeration order ⟹ pure FP re-association, bounded by (reduction depth) × ULP.
That is `vv-principles` criterion (3), satisfied.

**(ii) `most_inward` — `radial_characteristic_field.py:314` — a DISCRETE flip.**

```python
most_inward = ords[int(np.argmin(mu[ords]))]
```

`np.argmin` returns the FIRST minimum. The minimum η within a level is
`−sinθ·max(cos φ ... )`; for **even** `n_phi`, `φ = π` is a node so the minimum
is unique, but for **ODD** `n_phi` it is achieved by a mirror PAIR. Measured:

```
 n_mu  n_phi | min eta TIED? | multiplicity | most_inward flips?
    2      3 |          True |            2 |               True
    2      4 |         False |            1 |              False
    4      5 |          True |            2 |               True
    4      6 |         False |            1 |              False
    4      7 |          True |            2 |               True
    4      8 |         False |            1 |              False
    4      9 |          True |            2 |               True
    6     11 |          True |            2 |               True
    4     12 |         False |            1 |              False
    8     16 |         False |            1 |              False
```

⚠ Note what this says about TODAY: for odd `n_phi` the two candidates already
differ only by rounding noise, so **the current selection of "the most inward
ordinate" is already arbitrary** — decided by which way the last bit fell.
#325 does not introduce an ambiguity; it makes an existing one *visible and
deterministic*.

### …and line 314 is a Pattern-2 duplicate of the sort contract `[M]`

`level_indices[p]` is sorted **ascending by η** — so the most-inward ordinate is
`ords[0]` by construction, and the `argmin` re-derives it. Verified over
`n_phi ∈ {3,4,5,7,8,9,11,12,16,32}`:

```
  every case: argmin == ords[0]  True;  level really sorted ascending  True
```

**Consequence:** the whole ordering question collapses to a SINGLE decision —
how the level sort breaks ties among equal-η ordinates — after which
`most_inward = ords[0]` follows for free and the second order-sensitive site
stops existing. Two arbitrary choices become one documented one.

### The sort key is not injective — that is the root

`η = sinθ·cos φ` is **2-to-1** over φ ∈ [0, 2π). Sorting a full azimuthal circle
by η therefore never had a total order; the pre-#325 code was ordering the two
half-circles by rounding noise. Three candidate resolutions:

| option | what it does | cost |
|---|---|---|
| (a) `kind="stable"` | ties resolve to construction order (increasing φ index) | one keyword |
| (b) `np.lexsort((mu_y, mu_x))` | makes the key **injective** — ξ = μ_y separates the mirror pair | one line, fully principled |
| (c) represent the degeneracy | pair mirror partners explicitly in `LevelStructure` | a type change |

⚠ `np.argsort`'s default `kind="quicksort"` (introsort) is **not** a stable sort,
and its tie resolution is an implementation detail: measured, at
`n_mu=8, n_phi=32`, `quicksort != stable` on the tied arrays. So under #325 the
intra-level order is **not well-defined** unless the tie-break is pinned. A gate
asserting "the level is sorted by η" is **Mode-12 blind** to all of this —
sortedness is invariant under permuting equal elements.

## 3c. ⛔ RECONCILING THE TWO AGENTS — the ordering effect is REAL, ~1 % `[M]`

The two sub-agents returned **directly contradictory** verdicts on the load-bearing
fact:

| agent | verdict |
|---|---|
| test-architect | ordering alone moves the solve **~1e-2**, flat across `n_phi ∈ {32,64,128}` |
| explorer | ordering effect is **FP-only** (≤50 ULP), with a structural proof |

I reconciled it by my own measurement rather than adjudicating the arguments.

### My hypothesis — and it was REFUTED

The test-architect's harness uses a **random per-ordinate source**
(`m16_e2e_het.py:78-79`). I hypothesised that this breaks the mirror symmetry
between the tied partners, so the 1 % was an artefact of an unphysical source and
would collapse for a symmetric one. **Wrong.** Same harness, source varied:

```
SOURCE = random
  max|B - C| = 9.996633e-03   rel = 1.008e-02    <- ORDERING alone
  max|C - D| = 1.762316e-02   rel = 1.777e-02    <- argsort KIND alone
SOURCE = isotropic
  max|B - C| = 1.267787e-02   rel = 6.387e-03    <- ORDERING alone
  max|C - D| = 1.362045e-02   rel = 6.862e-03    <- argsort KIND alone
SOURCE = mirror_symmetric
  max|B - C| = 9.134825e-03   rel = 7.238e-03    <- ORDERING alone
  max|C - D| = 1.054731e-02   rel = 8.357e-03    <- argsort KIND alone

CONTROL: are the mirror partners' SOURCES actually equal?
  random            : max |src[a] - src[partner(a)]| = 9.554e-01
  mirror_symmetric  : max |src[a] - src[partner(a)]| = 5.551e-16
```

The control proves the symmetric source really is symmetric across partners, and
the effect **survives at ~7e-3**. The test-architect's number stands.

### Where the explorer's proof goes wrong — it is true, but about the wrong object

The proof: swapping a tied pair changes α only by `η(w_a − w_b)`, and weights are
constant within a μ-level, so the change is zero.

That is correct **about the α SEQUENCE at the half-positions** — `α_{m+½}` really
is invariant. It does NOT cover the **per-ordinate α ASSIGNMENT**. Ordinate `a`
moving from position `m` to `m+1` takes its α-pair from
`(α_{m−½}, α_{m+½})` to `(α_{m+½}, α_{m+3/2})`, and those differ by `w·η ≠ 0`.
The sequence is invariant; who receives which element of it is not. The solver
sees the second.

So both agents measured honestly and one drew an unsupported conclusion from a
true lemma. The `≤50 ULP` figure is consistent with a configuration where the two
orderings coincide (`quicksort == stable` for `n_phi ≤ 16` — measured in §3b) or
with an eigenvalue readout, which is Mode-12 blind here.

### What this means

**#325 is BLOCKED, and by a defect that predates it.** The per-level ordering is
under-determined by its own sort key, and the cylindrical angular flux depends on
the resolution at ~1 % — a discretization-INDEPENDENT error (flat across
`n_phi ∈ {32,64,128}`). Today that resolution is decided by rounding noise: the
legacy order is visibly scrambled relative to the symmetric interleave —

```
  D (stable):  [16, 15, 17, 14, 18, 13, 19, 12, 20, 11, 21, ...]   clean interleave
  B (legacy):  [16, 15, 17, 18, 14, 13, 19, 20, 12, 21, 11, ...]   noise-scrambled
```

⚠ **The open question is no longer "is it deterministic" but "which order is
CORRECT".** The user's ruling (lexsort on `(η, ξ)`) makes it deterministic and
geometry-determined, which is necessary — but at 1 % the choice also has to be
*right*, and no gate in the tree currently adjudicates it. That needs an
independent reference (the curvilinear MMS / L1 suite run under each candidate
ordering), not a convention.

## 3d. ⛔ WHY 3024 TESTS SAW IT ONCE — the coverage is in the blind cell `[M]`

The explorer swapped in the algebraic primitive and ran the whole consuming
surface: **3024 tests, exactly 1 failure**, and that failure is a guard whose
fixture depends on the defect existing. One snapshot moved —
`cyl_1g_homogeneous_product_dd_n20.npz` — by **2 ULP**.

A 1 % effect that breaks nothing is a coverage statement, not a magnitude
statement. Varying heterogeneity × groups × `n_phi`, isotropic source,
ORDERING-alone (legacy vs default argsort, node values bit-identical):

```
                            config  n_phi |  ORDERING-alone rel
------------------------------------------------------------------
1G homogeneous (the snapshot's regime)  4 |           6.627e-15
1G homogeneous (the snapshot's regime)  8 |           8.720e-15
1G homogeneous (the snapshot's regime) 16 |           3.366e-14
1G homogeneous (the snapshot's regime) 32 |           3.749e-14

                  1G heterogeneous      4 |           3.954e-02
                  1G heterogeneous      8 |           7.211e-02
                  1G heterogeneous     16 |           4.152e-02
                  1G heterogeneous     32 |           2.607e-02

                    2G homogeneous      4 |           9.132e-16
                    2G homogeneous      8 |           1.370e-15
                    2G homogeneous     16 |           1.979e-15
                    2G homogeneous     32 |           1.522e-15

                  2G heterogeneous      4 |           9.294e-03
                  2G heterogeneous      8 |           1.795e-02
                  2G heterogeneous     16 |           1.055e-02
                  2G heterogeneous     32 |           6.387e-03
```

**The discriminator is HETEROGENEITY, not group count.** Homogeneous nulls the
effect to FP noise in BOTH 1G and 2G; heterogeneous exposes it at **0.6 % – 7.2 %**,
at every `n_phi`, with no sign of converging away.

That is `vv-principles` anti-pattern #4 verbatim — *"NEVER accept a
homogeneous-only verification — flat flux nulls every redistribution and
weight-cancellation term"* — and Mode 9 (a defect verified only in a degenerate
regime). The single cylindrical `product` snapshot in the tree sits in the
**1G homogeneous** cell, which is precisely the cell where the error is
identically zero. It could never have caught this.

⚠ Note this also revises my §3c framing: I wrote that source symmetry was not the
discriminator, which is right, but the discriminator is *spatial heterogeneity*
— and homogeneous is exactly what the existing regression pins.

## 3e. ⭐ THE PROVABLE CRITERION — α has a closed form, and it picks the ordering `[M]`

The α-recursion is not bookkeeping. It is a **cumulative quadrature in the
azimuthal angle**, and that partial sum has an exact closed form — so the
correct ordering is decidable by algebra, with no MMS and no convention.

### The derivation

    α_{m+½} = α_{m−½} − w_m η_m,  α_½ = 0   ⟹   α_{m+½} = −Σ_{m'≤m} w_{m'} η_{m'}

With `η = sinθ·cos ω` and equal azimuthal weights `2π/n`, apply the Dirichlet
identity `2 sin(t₀/2)·Σ_{m'=0}^{m} cos(m' t₀) = sin((m+½)t₀) + sin(t₀/2)` with
`t₀ = 2π/n`:

> **α_{m+½} = −sinθ · (π/n)/sin(π/n) · [ sin(ω_{m+½}) + sin(π/n) ]**
>
> where `ω_{m+½} = ω_m + π/n` is the HALF-ANGLE boundary.

So α is an exact **affine function of `sin` at the half-angle boundary** — it is
`−ξ` there, up to the discrete quadrature factor `(π/n)/sin(π/n) → 1` and the
starting constant. This holds **iff the ordinates are enumerated in increasing
ω**, because the recursion is a cumulative integral in that variable.

### The measurement — `max |recursion − closed form|`

```
     n |  omega-order (W)  eta-order legacy (L)   lexsort (X)
     4 |        2.319e-16             2.954e+00     2.954e+00
     6 |        3.331e-16             2.954e+00     2.954e+00
     8 |        7.003e-16             3.044e+00     3.044e+00
    12 |        7.366e-16             3.003e+00     3.003e+00
    16 |        4.466e-16             2.970e+00     2.970e+00
    24 |        1.137e-15             2.928e+00     2.928e+00
    32 |        8.604e-16             2.905e+00     2.905e+00
```

**ω-ordering reproduces the closed form to machine precision. Both η-based
orderings — including the lexsort ruling — are off by O(1), ≈3.0, at every n.**

And the continuum limit confirms the identification (ω-order):

```
  n=   4: max|alpha - (-xi_half)| = 8.122e-01   factor (pi/n)/sin(pi/n) = 1.1107207345
  n=   8: max|alpha - (-xi_half)| = 3.920e-01                            1.0261721530
  n=  16: max|alpha - (-xi_half)| = 1.906e-01                            1.0064545428
  n=  32: max|alpha - (-xi_half)| = 9.383e-02                            1.0016081891
  n=  64: max|alpha - (-xi_half)| = 4.654e-02                            1.0004017082
  n= 128: max|alpha - (-xi_half)| = 2.318e-02                            1.0001004059
```

Halving with n — `α → −ξ` at O(1/n), exactly as the closed form predicts. Under
an η-ordering the "half-angle boundary" has no geometric meaning at all, because
consecutive entries in the enumeration are not adjacent in ω.

### Why no existing gate could have caught this — closure is Mode-12 blind

The one invariant the code does assert is closure, `α_{M+½} = 0`. Measured under
four orderings including a **random shuffle**:

```
  n=  5         W: -1.110e-16    L: -0.000e+00    X: -0.000e+00   shuffled: +1.110e-16
  n=  8         W: -1.110e-16    L: +1.110e-16    X: +1.110e-16   shuffled: -0.000e+00
  n= 16         W: -0.000e+00    L: +5.551e-17    X: +5.551e-17   shuffled: -0.000e+00
  n= 32         W: +8.327e-17    L: -0.000e+00    X: -0.000e+00   shuffled: +8.327e-17
```

Closure is a **telescoping sum**: it holds for EVERY permutation, including a
random one. It is `vv-principles` anti-pattern #8 / Mode 12 exactly — the
functional's invariance group contains the entire error class, so no tolerance,
refinement, or regime change could ever expose the ordering through it.

### ⚠ Status of the conclusion

The algebra above is unconditional. The step that makes it a verdict on
PRODUCTION is the premise that ORPHEUS's α is this cumulative sum taken over
`level_indices` order — **not yet verified against the sweep code**; an explorer
is establishing the recursion of record and whether the level spans the full
circle or a half-range. On a HALF range `ω ∈ [0,π]`, `η` IS monotone in `ω`, so
η-ordering and ω-ordering coincide up to direction and there is no defect at all.
**Do not act on this section until that premise is confirmed.**

## 4. MoC's reflection is a closed form wearing a search `[M]`

`orpheus/moc/geometry.py:222` `_reflected_azi_index` finds the reflected
azimuthal index by nearest-neighbour search in ANGLE space:

```python
phi_refl = np.pi - phi[azi_index]
...
return int(np.argmin(np.abs(phi - phi_refl)))
```

But the MoC node set is `φ_k = π(2k+1)/(2n)`, and

    π − φ_k = π(2n − 2k − 1)/(2n) = π(2(n−1−k) + 1)/(2n) = φ_{n−1−k}

exactly, as **integer index arithmetic**. Measured — the search returns exactly
`n−1−k` for every `n_azi` tested:

```
 n_azi                            search result                        n-1-k  agree?
     2                                   [1, 0]                       [1, 0]  True
     3                                [2, 1, 0]                    [2, 1, 0]  True
     4                             [3, 2, 1, 0]                 [3, 2, 1, 0]  True
     5                          [4, 3, 2, 1, 0]              [4, 3, 2, 1, 0]  True
     6                       [5, 4, 3, 2, 1, 0]           [5, 4, 3, 2, 1, 0]  True
     8                 [7, 6, 5, 4, 3, 2, 1, 0]     [7, 6, 5, 4, 3, 2, 1, 0]  True
    16 / 32 / 64                                                              True
```

And the reason the search exists is visible: the angle identity is NOT bit-exact
in float, which is what a nearest-neighbour match is for.

```
  n=  2: max|pi-phi_k - phi_(n-1-k)| = 0.000e+00   bit-exact: True
  n=  3: max|pi-phi_k - phi_(n-1-k)| = 4.441e-16   bit-exact: False
  n=  4: max|pi-phi_k - phi_(n-1-k)| = 0.000e+00   bit-exact: True
  n=  8: max|pi-phi_k - phi_(n-1-k)| = 2.220e-16   bit-exact: False
  n= 16: max|pi-phi_k - phi_(n-1-k)| = 4.441e-16   bit-exact: False
```

Not a live bug — the noise is ~4e-16 against an angular spacing of π/n, so the
argmin is never in doubt. But it is an O(n_azi) search per call, called twice
per track in `_build_reflection_links` (`geometry.py:385,393`), computing
something that is exactly `n−1−k`.

**Both wrap guards are dead code.** `φ_k ∈ (0, π)` for every `k`, so
`π − φ_k ∈ (0, π)` always and neither `if phi_refl < 0` nor `if phi_refl >= np.pi`
can ever fire.

This replacement is exact, O(1), and **independent of the node-value change** —
it stands on its own whether or not the rest of #325 lands.

---

## 5. The primitive as landed

`orpheus/numerics/roots_of_unity.py` — new file. Doctests clean;
`npx pyright orpheus/numerics/roots_of_unity.py` → `0 errors, 0 warnings`.

Named `roots_of_unity`, not the issue's placeholder `circle_nodes`: these ARE
the q-th roots of unity (`ζ_q^p = e^{2πip/q}`), and that is the object's precise
native term.

Three deltas from the issue body's prototype:

- **The axis needs no branch.** `r == 0` ⟹ `theta == 0.0` ⟹ IEEE gives exactly
  `1.0` / `0.0`. It falls out of the construction; the issue implied an explicit
  case.
- **Signed zero canonicalised** via `+ 0.0` (exact identity elsewhere). Without
  it the sign flips emit `-0.0`, 1–2 per q. Note any gate written with `==` is
  blind to this, since `-0.0 == 0.0`.
- **Integer dtype required**, float numerators refused rather than truncated.

Home is `orpheus/numerics/`, deliberately NOT the quadrature package — MoC is a
consumer and MoC is not a quadrature family. The module docstring names
`orpheus/numerics/symmetry.py` as the *checking* face of the same concept
(`is_invariant` verifies what this module generates).
