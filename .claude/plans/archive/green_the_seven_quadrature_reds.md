# The full tree is GREEN and the 162-commit branch is merged

> **Goal.** `main` carries the operator-strategy campaign. Today the branch is
> **162 commits ahead of `main`, 0 behind, 7 red**, so none of it has shipped —
> boundary machinery B0–B3.4c, geometric consolidation G1–G6.1, quadrature
> Q0–Q5.E, and the affine-boundary campaign P1–P6 are all stranded.
>
> **Done when** `tests/ -m "not slow"` is 0 failed and the branch has ff-merged.
>
> **⛔⛔ THE DEFERRAL IS REFUTED — 2026-08-06, same day, by measurement.**
> See "§ REFUTED" below. The rationale was that 4 of the 7 are #327-fork
> dependent; measured, **none of the seven is**. Greening is UNBLOCKED and was
> the right call all along.
>
> ~~**DEFERRED 2026-08-06 — gated on #327, and the reason is structural, not
> scheduling.** Adjudicate all seven in ONE pass **after** the `level_symmetric`
> fork is resolved.~~ Verification criterion (1) is already discharged (below),
> so the remaining work is the re-capture, not the investigation.
>
> **Proposed means** (*hypothesis until each red is adjudicated*): re-baseline
> the quadrature-induced reds citing the L0/L1/L2 evidence, and re-pose the
> gates whose instrument has outlived its claim (**#333**).

⚠ **NOT the same seven as `.claude/plans/archive/green_the_seven_baseline_reds.md`.**
That plan greened a DIFFERENT seven on 2026-06-26 (`574cff8`), merged, and is
archaeology. Its **discipline** is the playbook here; its inventory is not.

---

## ⚠️ THE NON-NEGOTIABLE DISCIPLINE (inherited verbatim; `vv-principles` L11 + #6/#10)

**NEVER regenerate a bit-identity snapshot just to make it green.** A snapshot
re-captured from a buggy path masks the bug forever. Per red:

1. **Verify the CURRENT production output is correct** against a
   **structurally-independent** reference.
2. Only if (1) passes, re-capture, **citing the verification in the commit**.
3. If (1) FAILS it is a **real regression** — STOP, do NOT re-capture, dispatch
   `numerics-investigator`.

---

## `[M]` Root cause — one cause, seven symptoms

All six non-#33 reds trace to **`81689a58`** ("re-capture the SN baselines the
quadrature consolidation moved"), which rebuilt the Gauss construction to impose
the measure's reflection symmetry rather than inherit it, moving the nodes by a
few ULP.

**L0 check against `numpy.polynomial.legendre.leggauss`** — a structurally
independent implementation:

| n | max node Δ | ULP | antisymmetry (ours / numpy) |
|---|---|---|---|
| 4 | 1.110e-16 | 1.0 | 0.0 / 0.0 |
| 8 | 2.220e-16 | 3.0 | 0.0 / 0.0 |
| 16 | 3.331e-16 | 3.0 | 0.0 / 0.0 |

confirming the commit's "1–4 ULP in the nodes" claim.

**L0 exactness against the closed form** `∫₋₁¹ xᵏ dx` (the rule's *defining*
property, degree ≤ 2n−1):

| n | ORPHEUS worst rel. residual | numpy | |
|---|---|---|---|
| 4 | 2.220e-16 | 5.551e-16 | ORPHEUS tighter |
| 8 | 1.110e-16 | 4.441e-16 | ORPHEUS tighter |
| 16 | 4.996e-16 | 3.331e-16 | **numpy tighter** |

⛔ **State this honestly: NEITHER DOMINATES.** Both sit at 2–5e-16 on the
defining property, i.e. both are correct to floating point and the difference is
root-finder noise, not a quality gap. The re-baseline is justified by "the new
value is *correct*", **not** by "the new value is *better*" — and the n=16 row is
kept precisely because it refutes the tidier claim.

---

## The seven — inventory, instrument, and reading

| # | test | instrument | `[M]` reading |
|---|---|---|---|
| 1 | `test_streaming_operator::…cart2d_1g_vacuum_apply_principled_equiv` | ULP budget 256 | **1152 ULP** |
| 2 | `…cart2d_2g_specular_apply_principled_equiv` | ULP budget 256 | **296 ULP** |
| 3 | `test_affine_carve_bit_identity[si_2d_p1_aniso_het]` | **sha256** | hash ≠ hash |
| 4 | `…[krylov_2d_p1_aniso_het]` | **sha256** | hash ≠ hash |
| 5 | `…[si_slab_2g_het]` | **sha256** | hash ≠ hash |
| 6 | `test_diamond::…test_spherical_inward_bit_identical` | `np.array_equal` | values print identically; differ in last bits |
| 7 | `TestWhiteXminPartial03GLSnapshot` (task #33) | `rtol = 4ε` | 8/60 elements, **max abs 2.22e-16** = 1 ULP at 0.2 |

Reds 1/2/6/7 are FP-scale. **Reds 3/4/5 have no reading at all**, which is the
finding below.

---

## ⭐⭐ Finding — a `sha256` gate is unfalsifiable-in-magnitude by construction

`test_affine_carve_bit_identity` stores **only hex digests** (`GOLDEN` is a dict
of `sha256`s frozen at pre-carve commit `63719a2`). The pre-carve VALUES were
never kept. So when it reddens:

* the magnitude of the change **cannot be computed**;
* the new value **cannot be compared** with the old at all;
* 1 ULP and a catastrophic error are **the same red**.

⚠ **The choice was defensible when made** and the docstring says why: the DD
regression snapshots compare at ~1e-11 and already pre-drift ~6920 ULP, so they
*cannot* verify a ZERO-change claim; a `sha256` was the sharpest instrument
available for it. The defect is not "sha was wrong".

**The defect is that a zero-change claim pinned against a FROZEN PAST has a
shelf life.** It is falsified by the first legitimate change anywhere upstream,
and at that moment the instrument gives no way to decide whether the new value
is fine. Worse, the gate's *claim* ("the #208 carve changed nothing") and its
*instrument* ("bytes equal a golden from `63719a2`") are the same statement only
while nothing else changes — the quadrature rebuild broke that coupling, so the
gate now reddens for a true statement about a **different** change.

⛔ **It cannot be re-posed as a same-run A/B** (which is what the claim actually
is): the pre-carve untyped path was retired, so there is no second side left to
compare against. This is the `coding-standards` "a retirement's rewire can
demote a gate" family, one step further along — the other side is gone entirely.

**Disposition (proposed, pending the L1/L2 evidence):** re-pose from a
zero-change-vs-frozen-past claim to a **drift tripwire with stored VALUES and a
documented ULP budget**, and record that #208's zero-change claim is now
*historical* — verified at `63719a2`, not re-verifiable, because only a hash
survives. A future red then reports a magnitude instead of a hash mismatch.

---

## Verification status

- [x] L0 — quadrature nodes vs `numpy.leggauss`: 1–3 ULP.
- [x] L0 — exactness vs the closed form: both at ε; neither dominates.
- [x] **L1/L2 — the SN physics suite**: `[M]` 2026-08-06,
      `pytest tests/sn -m "(l1 or l2) and not slow"` → **462 passed, 0 failed,
      7 xfailed** in 6:38. MMS convergence, analytical references and eigenvalue
      benchmarks are all structurally independent of any bit-identity snapshot,
      so this is criterion (1) discharged for every red below: **the current
      production values are correct.**
      ⚠ Re-measured rather than inherited — `81689a58` claimed "347 L1/L2 SN
      tests pass UNCHANGED" and the count today is 462. The claim holds; the
      number in it was already stale, which is exactly why the discipline says
      verify rather than cite.
- [ ] per-red adjudication + re-capture — **BLOCKED on #327**
- [ ] full tree `-m "not slow"` = 0 failed
- [ ] ff-merge

---

## ⛔ THE DEFERRAL DECISION (2026-08-06) — read this before restarting

### `[M]` Which quadrature family each red rides on

| red | fixture | quadrature |
|---|---|---|
| 1, 2 | `_cart2d_mesh` | **`level_symmetric(sn_order)`** |
| 3, 4 | `*_2d_p1_aniso_het` | **`level_symmetric(sn_order=4)`** |
| 5 | `si_slab_2g_het` | `gauss_legendre(8)` |
| 6 | spherical diamond | `gauss_legendre(4)` |
| 7 | `white_xmin_partial_03_GL` | `gauss_legendre` |

**Four of seven ride on `level_symmetric`.**

### Why that blocks the re-capture

**#327's acceptance criteria contain an UNDECIDED FORK:**

> `level_symmetric_sn` **either implements the Carlson–Lathrop moment-matched
> weights it claims, or is renamed** to what it is (an equal-weight `O_h` orbit
> rule) and its tag corrected to 3.

* **Fork A (implement Carlson–Lathrop):** replaces `w_octant = 4π/(8·n_octant)`
  — `[M]` **one distinct weight in the whole rule at every order** — with
  per-level moment-matched weights. Every weighted sum moves. Those four
  baselines move **materially, not by ULP**.
* **Fork B (rename + correct the tag):** metadata only; they do not move at all.

⟹ Re-capturing now freezes a value with a **known possible expiry**, and — worse
than wasted work — blesses it with this plan's documented verification ritual,
so the later #327 change reads as a **regression from a baseline we ourselves
created and certified**. That inverts the meaning of the gate.

### Why not green the three `gauss_legendre` ones now

**Reds 3, 4, 5 are ONE gate** — one module, one instrument, one `parametrize` —
and two of its three arms are `level_symmetric`. Re-posing one arm while two stay
red fragments the record and spends the adjudication twice. That leaves reds 6
and 7, which cannot unblock the merge alone. Splitting buys nothing and costs the
single-pass argument.

### Why the merge can wait (the fact that tips it)

`[M]` **`main` is 0 commits BEHIND** — 162 ahead, nothing landing on it. There is
no divergence pressure and the `--ff-only` merge stays trivially available, so
this is *not* the ordinary long-lived-branch risk where waiting compounds. Had
anything been landing on main, the balance would go the other way and greening
first would be correct.

⚠ **Standing condition:** that argument expires the moment a commit lands on
`main` from anywhere. Re-check `git rev-list --count HEAD..main` before relying
on it again.

### The reds are not hiding a live error

#327 measures ℓ=0 and ℓ=1 **exact** ("which is why isotropic and P1 transport
look healthy"), and reds 3/4 are **P1**-anisotropic. So the deferred configs are
correct today; what is uncertain is only whether their VALUES will move.

### ⛔⛔ REFUTED 2026-08-06 — the fork-dependence claim was FALSE

**The claim above** — "4 of 7 ride on `level_symmetric`, and #327 Fork A moves
every weighted sum materially, not by ULP" — **is measured false. It is left in
place because the way it was wrong is the lesson.**

`[M]` **`O_h` orbit count of the node set** (an orbit = one multiset
`{|x|,|y|,|z|}`; a single orbit + `Σw = 4π` *uniquely forces* equal weights, so
there is no freedom for a moment-matched solve to exploit):

| order | nodes | orbits | consequence |
|---|---|---|---|
| **S2** | 8 | **1** | **equal weights FORCED** |
| **S4** | 24 | **1** | **equal weights FORCED** |
| S6 | 48 | 2 | 2 free weights |
| S8 | 80 | 3 | 3 free weights |
| S12 | 168 | 5 | 5 free weights |
| S16 | 288 | 8 | 8 free weights |

`[M]` and **every one of the four `level_symmetric`-riding reds uses
`sn_order=4`** — `_capture_pre_t4_snapshots.py:175` (default `sn_order=4`),
`test_streaming_operator.py:101,480`, `test_affine_carve_bit_identity.py:158`.

⟹ At S4 a genuine Carlson–Lathrop solve returns **the weights already there,
bit-for-bit**. Fork A cannot move those four. Reds 5/6/7 are `gauss_legendre`
and were never fork-dependent either. **None of the seven is.**

⭐ **How the error was made, so it is not repeated: the blast radius was
inferred from the ISSUE'S PROSE ("implement the moment-matched weights") rather
than measured from the REALIZATION (count the orbits).** One `np.unique` over
sorted `|node|` triples — six lines — refutes a committed rationale. This is the
"a scope read from PROSE gets refuted by the REALIZATION" failure already on
record from the G6.3 campaign, recurring in a new tier.

⭐ **The deeper correction to #327 itself:** the defect is NOT "the rule is wrong
at every order". At S2/S4 the equal-weight implementation is **already correct
and provably optimal for its node set** — degree 3 is the ceiling there, and the
advertised 3 at S4 is RIGHT. The defect is that **S6+ leaves 2/3/5/8 free
weights unsolved**. That is a much narrower statement than the issue's title.

### What survives, and what does not

* **Does NOT survive:** the fork-dependence rationale, refuted candidate **R5**
  (whose premise was the same false claim), and the "adjudicate after #327"
  gating.
* **Survives unchanged:** every verification measurement below (L0 nodes, L0
  exactness, L1/L2 462-green), the finding that a `sha256` gate is
  unfalsifiable-in-magnitude (**#333**), and the observation that reds 3/4/5 are
  one gate that should move together.
* **New coupling, in the other direction:** Fork A moves exactly **one** GREEN
  baseline — `2d_octant_equivalence_04_vacuum_2g_het_gradientQ_LS6.npz` (S6, so
  2 free weights). So #327 and the baseline set ARE coupled, just not through
  these seven.

### Correction to this plan's own earlier framing

This file opened by calling greening the top priority. On the evidence it is
not — **#327 is**, and greening is downstream of it. The earlier framing is left
in place per `plan-authoring` §3.

---

## Refuted candidates

| # | candidate | structural reason it fails |
|---|---|---|
| **R1** | Bump the ULP budgets (256 → 2048) and move on | Hides the question. A budget raised to fit the observed drift asserts nothing; the next legitimate change raises it again. The budget must follow from the reduction depth, not from today's reading. |
| **R2** | Re-freeze the `sha256`s | Each re-freeze silently weakens the claim to "nothing changed since I last re-froze", and leaves the *next* red as undiagnosable as this one. Now **#333**. |
| **R5** | Green the 4 `level_symmetric` reds now and re-do them after #327 | The re-capture is cheap; the *certification* is not. A baseline blessed by this plan's verification ritual makes the subsequent #327 move look like a regression against a trusted reference — the ritual is what must not be spent twice, not the keystrokes. |
| **R6** | Green only the 3 `gauss_legendre` reds now | Reds 3/4/5 are one gate with two `level_symmetric` arms; a partial re-pose fragments it. And 3 of 7 does not unblock the merge, which is the only thing greening is *for*. |
| **R3** | Retire reds 3/4/5 outright | The configs (2-D 2G P1-aniso het, SI and Krylov; slab 2G het) are good regression coverage. Retiring loses it; the instrument is what needs replacing, not the fixture. |
| **R4** | Green by pinning the OLD quadrature construction | Inverts the campaign — the reconstruction is the deliverable, and the L0 evidence says the new rule is correct. |
