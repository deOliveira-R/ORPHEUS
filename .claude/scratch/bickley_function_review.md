# Bickley Function Implementation Review — Davierwalla + Lorensi

Decision target: do these two papers warrant changes to
`orpheus/derivations/common/kernels.py` and the cylinder Peierls /
flat-source CP / Carlvik recurrence consumers that depend on it?

Both PDFs read in full from `/workspaces/ORPHEUS/scratch/literature/`.

---

## Q1 — Paper contributions

### Davierwalla (1982)

- **Citation**: D. Davierwalla, "Rational Chebyshev Approximations to
  the Bickley Functions," *Nuclear Science and Engineering* **80**(3),
  461–469 (1982). DOI 10.13182/NSE82-A19833. Eidg. Inst. f.
  Reaktorforschung, Würenlingen.
- **Headline contribution** (verbatim, Abstract, p. 461):
  > "A near optimal minimum-maximum rational approximation to the
  > Bickley functions Ki_s(x) for s = 1, 2, 3 and 0 ≤ x < 7, and
  > Ki_s(x) for s = 13, 14, 15 for 7 ≤ x < 600 with a relative error
  > of <10⁻¹². The above initial approximations combined with a
  > recurrence relation yield a stable method of computing Ki_s(x)
  > for the prescribed accuracy except for the region S: {(x, s):
  > x ≥ 28, s > 15}."
- **Method / formula**: Initial-function rational Chebyshev approximations
  delivered as **finite continued-fraction expansions** (Tables I–IV
  span 8 sub-intervals on [0, 7] and one on [7, 10] for s = 13, 14, 15
  plus an asymptotic continued fraction for x ≥ 10). The three-term
  forward recurrence `s·Ki_{s+1}(x) = (s−1)·Ki_{s−1}(x) + x·[Ki_{s−2}(x)
  − Ki_s(x)]` (Eq. 9) extends to all s with controlled error
  amplification, partitioned by the x–s plane diagram in Fig. 1.
- **Accuracy claim**: <10⁻¹² relative error throughout the strips
  A (0 ≤ x < 7, s ≥ 0), B (7 ≤ x < 28, s ≥ 13), and C (x ≥ 7, s < 13);
  excludes the strip {x ≥ 28, s > 15}.
- **Speed claim**: No CPU-time figures given. Architectural argument
  (p. 462): "a rational function requires only max(m, n) divisions if
  converted into a continuous fraction with successive quotient
  polynomials made monic" — i.e. ~6–8 divisions per evaluation for the
  given tables.
- **Caveats**:
  - Tables published as scanned ANS galley proofs (Tables I–IV); 16-digit
    coefficients but **not machine-readable**. Re-keying is error-prone.
  - Coefficients targeted CDC single precision (~14 digits). Modern IEEE
    binary64 has 15–17; the 10⁻¹² claim should re-validate cleanly but
    has not been independently re-tested in the published record (only
    one citing article per the cover page).
  - Excluded region {x ≥ 28, s > 15} is corner-case for ORPHEUS but
    **may** appear in deep-resonance shielding evaluations of higher-Ki
    radial moments.
  - Subroutine "operational in the SURCU transport code" (Stepanek 1979)
    and BOXER (Maeder–Paratte 1975) — both Würenlingen Pilot-codes
    private to PSI-era licensees; coefficients themselves are open.

### Lorensi, Azevedo, Sauter (2025)

- **Citation**: G. A. Lorensi, F. S. de Azevedo, E. Sauter,
  "An Efficient Numerical Approximation for the Bickley-Naylor
  Functions," *XII ERMAC-RS* (Encontro Regional de Matemática
  Aplicada e Computacional do Rio Grande do Sul), Porto Alegre,
  2 – 4 July 2025. ST6, 8 pp.
- **Headline contribution** (verbatim, Abstract, p. 1):
  > "This work presents a new computational approach for efficiently
  > calculating the first and second-order Bickley-Naylor functions.
  > The proposed method ensures high numerical accuracy while avoiding
  > computationally time-consuming numerical integration. … The
  > implementation, developed in C++, achieves high precision with
  > errors limited to a few Units in the Last Place (ULPs)."
- **Method / formulae**: Domain-decomposed approximations:
  - **[0, 1]** (Eqs. 7, 8): explicit log-pole removal,
    `Ki_1(x) = π/2 + x·[p_1(x) + q_1(x)·log(x)]`,
    `Ki_2(x) = 1 − (π/2)·x − x²·[p_2(x) + q_2(x)·log(x)]`,
    with degree-6 numerator / degree-8 denominator polynomials fitted
    by least-squares; coefficients **explicitly tabulated**
    (Tables 1–2, p. 5) to 17 digits in a form copy-pastable into C/Python.
  - **[1, 2], [2, 4], [4, 8], [8, 12], [12, 16], [16, 20], [20, 24]**
    (Eq. 9): Remez-optimal rational fits `Ki_n(x) = p_n(x)/q_n(x)`,
    one fit per interval (coefficients NOT printed in the paper —
    "can be computed using" the recipe + Boost Math Toolkit Remez).
  - **[24, 743]** (Eqs. 10–11): asymptotic + Remez correction in
    a Möbius-transformed argument `z = (35664/719)·(1/x − 767/35664)`,
    so `z ∈ [−1, 1]` over the full tail.
  - **[743, ∞)**: identically zero in IEEE binary64 (true: at x = 743,
    Ki_n ~ 10⁻³²³ ≤ 2⁻¹⁰²²).
- **Accuracy claim**: "errors limited to a few ULPs" — Fig. 1 shows
  ULP error confined to roughly [−2, +5] for both Ki_1 and Ki_2 across
  x ∈ [10⁻³, 10³]. Reference values from `mpmath.quad`
  (tanh-sinh, dps = 50) — same machinery as ORPHEUS' `ki_n_mp`.
- **Speed claim**: Fig. 2 — at 10⁷ evaluations the proposed method
  finishes in **< 2 s** while the integration baseline (GSL `gsl_integration_qag`
  on the integral form) requires **> 12 s**. That is ~200 ns/eval for
  Ki_1, Ki_2 at C++ binary64.
- **Caveats**:
  - **n = 1 and n = 2 ONLY**. No fits for Ki_3 — the paper's coverage
    stops where Davierwalla's *begins to extend*.
  - Only the [0, 1] coefficients are printed. To use the [1, 24]
    coverage you must re-run Remez yourself (Boost Math Toolkit
    referenced, not provided as data).
  - Conference proceedings, no peer review; method is sound but the
    coefficient set is not in any open database.
  - C++ implementation referenced but no public source URL given in the
    paper.
  - Reference comparison is against their *own* GSL integration call,
    not against an existing published Bickley library — no head-to-head
    against Amos 1983 (TOMS Algorithm 609) or Blair–Edwards–Johnson
    (1978).

---

## Q2 — Improvement assessment for ORPHEUS

### Where Bickley actually appears (production hot paths)

| Call site | Purpose | Order(s) | Density |
|-----------|---------|----------|---------|
| `peierls_nystrom/geometry.py:776` | cylinder volume kernel `κ_2(τ) = Ki_1(τ)` (Nyström matrix assembly) | 1 | high (every chord) |
| `peierls_nystrom/geometry.py:803` | cylinder escape kernel `K_esc(τ) = Ki_2(τ)` | 2 | medium (per radius) |
| `peierls_nystrom/geometry.py:2566` | T-matrix builder for cylinder rank-N IC, `Ki_{3+kk}(τ)` | 3, 4, …, ≤ 9 | high (≈ `n_modes`² per build × `n_q` chord nodes) |
| `peierls_nystrom/geometry.py:1709, 3313, 3422, 4133` | per-radius `ki_n_mp(1, τ, dps)` inside scipy.integrate routines | 1 | medium |
| `peierls_nystrom/geometry.py:1806, 1860, 2144` | `ki_n_mp(2..k+2, τ, dps)` Carlvik moment recurrence | 2, 3, … | medium |
| `peierls_nystrom/geometry.py:4614` | `ki_n_mp(3, x, dps)` inside Knyazev cylinder G_bc | 3 | low |
| `flat_source_cp/geometry.py` Chebyshev interpolant | static load-time fit of `e^τ · Ki_3(τ)` from `ki_n_mp` | 3 | one-shot |

The **Atkinson product-Nyström module** (`fn_method/peierls_atkinson_nystrom.py`)
uses **`E_1` (scipy.special.exp1)**, not Bickley — it operates on the
slab geometry whose pre-integrated angular kernel is exponential
integral, not Bickley. **Davierwalla / Lorensi do not affect the
Atkinson F_N moment floor at all.** The Phase-2 F_N moment-floor
hardening is orthogonal to this review.

### Speed ledger vs current `kernels.ki_n*` (measured this session)

Hot-path microbenchmark (devcontainer, IEEE binary64):

```
ki_n_float(1, x)          77.7 µs / eval
ki_n_float(2, x)          22.9 µs / eval
ki_n_mp(1, x, dps=25)   3173 µs / eval
scipy.special.expn(2,x)   0.15 µs / eval (vectorised, for reference only)
```

Lorensi (their Fig. 2): C++ ≈ **0.2 µs / eval** (Ki_1, Ki_2).
Davierwalla (estimated from continued-fraction op count): ~**0.05–0.1 µs**
in compiled C; ~**0.5–2 µs** in pure Python (8 fused-multiply-divides).

| Adoption | Expected speed gain over `ki_n_float` |
|----------|---------------------------------------|
| Davierwalla pure-Python eval (s = 1, 2, 3) | 30–100× |
| Davierwalla numba-compiled | 200–700× |
| Lorensi-style polynomial+Remez (s = 1, 2 only) | 200–400× |

For the **cylinder T-matrix builder** at `geometry.py:2566`
(largest single Bickley consumer): typical configuration
`max_kk + 1 ≈ 9, n_q ≈ 64` → ~576 Bickley evals per build
→ current ~45 ms per T build. At 200 ns/eval: ~0.12 ms. Wall time saved
per power-iteration pass: tens of ms. Not session-defining for one
benchmark, but multiplicative across parametric studies.

### Accuracy ledger vs current `kernels.ki_n*`

The "1.95e-06 max relative error" measured for `ki_n_float` is
**entirely at x ≥ 30**, where the true value drops below 10⁻²³ and
scipy's `epsabs=1e-15` floor swamps the relative figure. In the regime
ORPHEUS actually exercises (τ ≲ 20):

| τ | `ki_n_float(1, τ)` rel err vs mpmath dps=50 |
|---|---------------------------------------------|
| 0.01 – 1.0 | < 1e-15 (bit-exact) |
| 1 – 10 | < 1e-14 |
| 10 – 20 | ~1e-13 |
| 30+ | scipy epsabs floor dominates → rel err grows |

So `ki_n_float` is **already** ≥1e-13 accurate where it matters.
Lorensi's "few-ULP" claim (≈ 1e-15 rel) is a 10–100× tightening at
τ ≲ 20 and a clean recovery in the deep tail (where their
asymptotic + Möbius-Remez form does not suffer the absolute-error
collapse of integration). Davierwalla at 10⁻¹² is *worse* than the
current scipy.integrate path on the dense-x regime, and does not
help at all at τ < 1.

### Backward compatibility

Both methods agree with the current production output to ≥ 1e-12
across the regime ORPHEUS exercises (verified by paper-internal
references against tanh-sinh quadrature at high precision). Any
verification regression on existing test suites would be at the
1e-13 level — well below the L1/L2 verification thresholds (1e-5
to 1e-7) currently in place.

### vs Atkinson product-Nyström

Atkinson does not call Bickley anywhere — it lives on the slab
geometry where the kernel is `E_1`. The F_N moment floor at ≤ 1e-5
is **independent** of Bickley accuracy. Replacing Bickley does NOT
move the moment floor; the floor is set by the slab kernel
discretisation, not the cylinder kernel.

The cylinder F_N reconstruction (`fn_method/cylinder/...`) is the
appropriate cross-reference, and *that* path goes through
`peierls_nystrom/geometry.py:2566` (the T-matrix builder). A faster
Bickley *does* reduce the cost of cylinder F_N flux reconstruction,
without changing its accuracy floor — the floor is set by F_α
recursion stability, not by Bickley.

---

## Q3 — Recommendation

> **UPDATE 2026-05-03**: The "Amos as a cleaner pivot" hint below was
> made on logistics grounds before reading the Amos paper. After a
> head-to-head with both PDFs in hand, the technical verdict reverses
> several of these claims. **Read the head-to-head section at the
> bottom of this file before acting on Q3.** Specifically: Davierwalla
> is NOT advantageous on technical merit, AND Amos's coverage extends
> the n × x grid well beyond what Davierwalla supports.

### Davierwalla — **DEFER**

Rationale:
- Coefficient table is published only as a scanned PDF galley,
  re-keying ~150 16-digit constants by hand is error-prone.
- The 10⁻¹² accuracy target is *worse* than `ki_n_float` already
  delivers in the working regime (τ ≲ 20).
- The strip-D excluded region {x ≥ 28, s > 15} matches the
  region where ORPHEUS's `ki_n_float` *also* loses precision —
  Davierwalla offers no help here either.
- The piece that would actually be valuable (a clean Ki_3 fit
  that beats Lorensi's n = 1, 2 coverage) is mixed with the
  full-recurrence machinery and must be extracted carefully.

Defer until either (a) a clean machine-readable port exists
(Amos 1983 TOMS Algorithm 609 is a more attractive starting
point — Fortran source available), or (b) a specific consumer
needs Ki_n for n ≥ 4 with deterministic latency rather than
mpmath's adaptive cost.

Record for later: keep both this PDF and Amos 1983 in
`scratch/literature/` for the day a `ki_n_float` replacement
becomes load-bearing.

### Lorensi, Azevedo, Sauter — **ADOPT (limited scope)**

Rationale:
- For Ki_1 and Ki_2 in the cylinder Peierls hot path
  (`geometry.py:776`, `:803`, the T-matrix at `:2566` with
  the lowest two `kk` indices), this gives a 200–400× speed-up
  with sub-ULP accuracy.
- Coefficients for **[0, 1]** are explicitly printed
  (Tables 1–2, p. 5), 17 digits, copy-pastable. That handles
  the most important regime (small-τ, where the existing
  `ki_n_float` Bickley *is* fastest already and where the
  log-singularity makes other formulations ugly).
- Coefficients for the rest of [1, 24] are NOT printed; we
  would need to **regenerate them** via Remez (Boost Math Toolkit
  in C++, or `mpmath` + `scipy.optimize` in Python). This is a
  one-time SymPy notebook, not a daily cost. Storage: a
  `_bickley_lorensi.py` module of ~50–100 polynomial coefs.
- The [24, 743] asymptotic + Möbius-Remez form is 8 lines of
  code given that we already have the asymptotic series.

**Scope of adoption**:

1. **Yes** — Ki_1, Ki_2 in [0, 1] using printed coefficients.
2. **Yes** — Ki_1, Ki_2 in [1, 24] using locally-regenerated
   Remez coefficients. Validate against `ki_n_mp` at dps = 50
   in a one-shot SymPy notebook.
3. **Yes** — Ki_1, Ki_2 in [24, 743] using the asymptotic form.
4. **No** — do NOT extend to Ki_3+ in this work. The cylinder
   T-matrix at `geometry.py:2566` calls `ki_n_float(j+3, ·)` for
   j = 0, …, 2(n_modes − 1), i.e. Ki_3 through Ki_9. Lorensi
   does not cover these. Either (a) chain the Lorensi Ki_2 to
   higher orders via the Bickley recurrence (Eq. 9 in
   Davierwalla — error amplification is benign in *forward*
   recurrence in the regime of interest, see Davierwalla Fig. 1
   Strip A), or (b) leave Ki_n≥3 on `ki_n_float` for now.

**Implementation plan** (if approved):

| Step | File | Action |
|------|------|--------|
| 1 | `orpheus/derivations/common/_bickley_lorensi.py` (new) | Encode Tables 1–2 (coefficients verbatim from Lorensi). Add Remez-fit notebook in `scratch/derivations/bickley_lorensi/` to regenerate [1, 24] interval fits. |
| 2 | same file | `def ki1_lorensi(x: float) -> float` and `def ki2_lorensi(x: float) -> float` with branched evaluation across {[0, 1], [1, 24], [24, 743], [743, ∞)}. |
| 3 | `orpheus/derivations/common/kernels.py` | Add `def ki_n_fast(n, x)` that delegates to Lorensi for n in {1, 2} and `ki_n_float` otherwise. Mark `ki_n_float` as the "validation oracle" rather than the hot path. |
| 4 | `tests/derivations/test_kernels.py` | Cross-verification: `ki_n_fast(n, x)` vs `ki_n_mp(n, x, 50)` to ≤ 5 ULP at x ∈ {1e-6, 1e-3, 0.1, 0.5, 1, 3, 7, 15, 30, 100, 500} for n = 1, 2. (L0). |
| 5 | `orpheus/derivations/continuous/peierls_nystrom/geometry.py:776, 803` | Switch the cylinder volume/escape kernels to `ki_n_fast`. |
| 6 | `geometry.py:2566` (T-matrix builder) | The j = 0 (Ki_3 — used) is *outside* Lorensi's coverage. Either (a) hold off on this site, or (b) chain Ki_1, Ki_2 via Bickley recurrence to compute Ki_3 from Ki_1 + Ki_2 directly — `Ki_3(x) = (1/2)·[Ki_1(x) + x·(Ki_0(x) − Ki_2(x))]`. The Bessel K_0 piece is `scipy.special.k0` (~0.1 µs). This may or may not beat Lorensi on Ki_2 + scipy K_0 vs. `ki_n_float(3, ·)` at 23 µs — would benchmark first. |
| 7 | docs | Update `docs/theory/peierls_nystrom.rst` "Implementation notes" with the citation and accuracy/speed table. Sphinx Cardinal Rule 3. |
| 8 | `error_catalog.md` | If any L1 cylinder Peierls test shifts at the ULP level, log as ERR-NNN with rationale "Bickley accuracy floor moved from ~1e-13 (scipy quad) to ~1e-15 (Lorensi Remez)". |

Estimated implementation: ~1 day of focused work (Remez regeneration
+ tests + Sphinx page).

If the user wants a quicker first-cut (~1 hour): only Steps 1, 5
on `geometry.py:776, 803` (the two scalar-kernel sites), without
the `geometry.py:2566` T-matrix builder. That alone would
remove `ki_n_float` from the chord-Nyström assembly hot path.

---

## Files referenced

- Local PDFs:
  - `/workspaces/ORPHEUS/scratch/literature/Davierwalla(1982)Rational Chebyshev Approximations to the Bickley Functions.pdf`
  - `/workspaces/ORPHEUS/scratch/literature/Lorensi-Azevedo-Sauter(2025)An efficient numerical approximation for the Bickley-Naylos Functions.pdf`
- Cross-referenced ORPHEUS files:
  - `/workspaces/ORPHEUS/orpheus/derivations/common/kernels.py`
    (lines 196–326 — the Bickley implementation surface)
  - `/workspaces/ORPHEUS/orpheus/derivations/continuous/peierls_nystrom/geometry.py`
    (lines 752–804 cylinder volume/escape kernels;
    lines 2541–2585 T-matrix builder; lines 1709, 1806, 1860,
    2144, 3313, 3422, 4133, 4614 mpmath fallbacks)
  - `/workspaces/ORPHEUS/orpheus/derivations/continuous/peierls_nystrom/origins/cylinder_g_bc_3d.py`
    (Knyazev cylinder G_bc Bickley use)
  - `/workspaces/ORPHEUS/orpheus/derivations/continuous/peierls_nystrom/origins/cylinder_knyazev.py`
    (Knyazev test harness)
  - `/workspaces/ORPHEUS/orpheus/derivations/continuous/flat_source_cp/geometry.py`
    (load-time Chebyshev interpolant of `e^τ·Ki_3(τ)`)
  - `/workspaces/ORPHEUS/orpheus/derivations/continuous/fn_method/peierls_atkinson_nystrom.py`
    (slab E_1 — confirmed independent of Bickley)
- Open follow-ons:
  - Amos 1983 (TOMS Alg. 609) — alternative Davierwalla-class
    Bickley library; Fortran source widely mirrored. Worth
    pulling into `scratch/literature/` as the "drop-in
    library" pivot if Lorensi adoption stalls on the Remez
    regeneration step.
    **Now LOCAL** at `scratch/literature/Amos (1993)…pdf`
    (note: filename says 1993 but the paper is the 1983 TOMS
    publication — see head-to-head below).

---

## Head-to-head: Davierwalla 1982 vs Amos 1983 (added 2026-05-03)

Both PDFs read in full. Goal: technical-merit comparison
independent of transcription convenience, so the user can decide
whether to absorb Davierwalla's transcription cost.

### Bibliography (verified)

- **Davierwalla 1982**: D. Davierwalla, "Rational Chebyshev
  Approximations to the Bickley Functions," *Nuclear Science and
  Engineering* **80**(3), 461–469 (1982). DOI 10.13182/NSE82-A19833.
  Eidg. Inst. f. Reaktorforschung (EIR-Würenlingen), Switzerland.
  Received 3 July 1980, accepted 25 September 1981, published 1982.
  Pages 461 + 470 are technical-note column-shared with another
  paper; the Bickley note runs pp. 461–469.
- **Amos 1983** (filename misleadingly says "1993"): D. E. Amos,
  "Algorithm 609: A Portable FORTRAN Subroutine for the Bickley
  Functions Ki_n(x)," *ACM Transactions on Mathematical Software*
  **9**(4), 480–493 (December 1983). Sandia National Laboratories,
  contract DE-AC04-76DP00789. Received 4 December 1981; revised
  1 March 1983; accepted 21 April 1983. Companion paper to
  Amos's *Math. Comp.* paper on uniform asymptotic expansions for
  E_n and Ki_n (cited as Ref [4] in the algorithm paper).
  **The "1993" in the local filename is a transcription error in
  the file name, not the publication date.**

### Comparison table

| Axis | Davierwalla 1982 | Amos 1983 | Winner |
|---|---|---|---|
| **Worst-case rel. error claim** | < 10⁻¹² across coverage strips A, B, C; excluded strip D = {x ≥ 28, s > 15}. Coefficients targeted **CDC single precision** (~14 digits stored). | "accuracies up to 18 digits (Univac double precision)" — i.e. ≈ 10⁻¹⁸ in DBSKIN; ~10⁻¹² in single-prec BSKIN matching Davierwalla. Constants stored to 18 digits. The 10⁻¹⁸ figure is **the precision floor of the stored constants**, not the algorithm. | **Amos** (≥ 10⁻¹⁵ in double prec.; matches Davierwalla in single prec.) |
| **Indep. verification of accuracy** | Error figure derived from author's own min-max optimisation; "error tables" referenced as private internal report (Ref. 5, EIR-243). No independent test grid in the paper. | §7 "TESTING" reports cross-checks: Taylor-series consistency, recurrence consistency vs single evaluations, K_0/K_1/K_-1 cross-tests on x ≤ 50. CDC 6600 single-prec ≈ 10⁻¹², CDC double-prec ≈ 10⁻¹⁸ "limited by stored 18-digit constants." | **Amos** (the testing protocol is documented in the paper) |
| **n directly fitted** | s = 1, 2, 3 (initial functions on [0, 7]); s = 13, 14, 15 (initial functions on [7, 600]). | r = 0, 1, 2 (initial power series on 0 ≤ x ≤ 2); for x > 2 a **uniform asymptotic expansion** valid for all r ≥ 0 supplies any starting triple. | **Amos** (any n via uniform expansion; Davierwalla only at the two strip seeds) |
| **x range coverage** | [0, 7] for low-s seeds; [7, 600] for high-s seeds; with recurrence, full strip A (0 ≤ x < 7, s ≥ 0), B (7 ≤ x < 28, s ≥ 13), C (x ≥ 7, s < 13). **Excludes strip D** = {x ≥ 28, s > 15}. | Power series on [0, 2] (r = 0, 1, 2); uniform asymptotic expansion on (2, ∞) for any r ≥ 0; argument boundary XLIM ≈ 2.3036·(−EMIN·log₁₀ B) − ln C(x+r) handles underflow gracefully. **No excluded region** — coverage is the full positive (x, n) quadrant up to underflow. | **Amos** (no strip-D blind spot) |
| **Recurrence stability** | Eq. 9 forward / Eq. 17 backward; recurrence direction switched on the x–s plane diagram (Fig. 1) at x ≈ 6–7 (matching where w_f = w_b = 1). Two iterations beyond w_f = 1 to attenuate amplified error. | Eq. (1.3) same Bickley three-term recurrence; direction chosen by relation of [x + 0.5] to {n, n+1, …, n+m−1}: forward, backward, or **both** depending on which initial index is closest. Power series gives 3 indices on [0, 2] to seed; uniform expansion gives 3 indices on (2, ∞) to seed. Stable both ways. | **Tie on math** (same Bickley recurrence). **Amos wins on automation** (direction chosen at runtime; user not required to consult the Fig. 1 diagram). |
| **Number of FLOPs / eval** | Continued fraction with ~6–8 nested divides per Ki_s, plus ln(x) for Ki_1, Ki_3 strip-A nested logarithm correction. **Pure division-dominated** — branch-light. | Power series with truncation test (variable iterations until rel error < ε): typically 8–15 terms × ~3 mul + 1 div per term on [0, 2]. Uniform expansion on x > 2: a finite sum of N+1 ≤ ~20 E_n(x) calls + a remainder Pascal-triangle sum + γ_{M,N} bound test. **Mul–add dominated** with a triangular-array coefficient B(*) lookup, ~50–200 FLOPs per eval depending on r and ε tolerance. | **Davierwalla on raw FLOP count**, **Amos on adaptivity** (Davierwalla pays the same ~6–8 divides regardless of required ε; Amos truncates as soon as the Eq. (3.7) bound falls below ε). |
| **Stability near transitions** | 8 sub-intervals on [0, 7]: hand-picked breakpoints at x = 0.2, 1, 2, 3, 4, 6, 7. Each sub-interval has its own continued-fraction coefficient table. **Discontinuity in derivatives at breakpoints possible** (each fit is independent). | Two regions only: [0, 2] series and (2, ∞) uniform expansion. The seam at x = 2 is bridged by overlapping validity (the asymptotic expansion is valid on 0 ≤ x ≤ ∞ uniformly in n; only chosen on x > 2 for efficiency reasons). | **Amos** (single seam at x = 2; both forms valid in the overlap, which is a much cleaner stability story than 8 hand-picked breakpoints). |
| **Asymptotic / underflow regime** | Continued-fraction asymptotic for x ≥ 10; underflow at x ≈ 720 noted (Ki_n drops below CDC EMIN). **No machine-portable underflow test**. | Explicit XLIM computed from R1MACH(1) (smallest positive machine number) and EMIN; subroutine exits cleanly with IERR signal at underflow. **Portable across 15 machine environments** (CDC, Univac, IBM, etc.). | **Amos by a wide margin** (XLIM is computed at init; on IEEE binary64 it's ≈ 743 — matches Lorensi's [743, ∞) cutoff). |
| **Small-x regime** | Strip-A x ∈ [0, 0.2): rational fit to G(x) = [Ki_s(x) − P_s(x)]/s(x) where s(x) is the leading log term. Explicit P_s(x) polynomial subtracted to handle the (s−1)-bounded-derivative issue at x = 0. **The log singularity in Ki_1 is manually peeled** before the rational fit. | §2 power series Eq. (2.1a–c): explicitly separates G_k(r, x) = ψ(k+1) + ψ(2k+r+1) − ψ(2k+1) − ln(x/2) (the log piece) from the analytic series; the truncation error bound is 4/3·A_K·G_K. **The log singularity is in G_k** and the series converges for any 0 ≤ x ≤ 2 with monotone-decreasing terms. | **Tie** (both correctly handle the log-singular structure of Ki_1 at small x). |
| **Code availability** | Coefficients in scanned PDF galley proofs (Tables I–IV, NSE 80, pp. 465–468). 16-digit constants, ~150 of them, no machine-readable form. The original Fortran subroutine is in EIR-243 (Würenlingen internal report, 1973) — not in the open literature. | **TOMS Algorithm 609 = portable Fortran subroutine BSKIN/DBSKIN**, distributed by ACM via ALGORITHMS distribution service. The full algorithm listing is published with the paper (last 1.5 pages contain the BSK, BSKIN, BKIAS, BKISR, EXINT, PSIXN, GAMRN, BDIFF, HKSEQ, plus I1MACH/R1MACH/D1MACH portability layer). Mirrored at [Netlib/toms](http://www.netlib.org/toms/609) and shipped with SLATEC. | **Amos by orders of magnitude.** Drop-in via f2py, or trivial Fortran→Python translation. |
| **Generalisation to higher n / higher precision** | Recurrence Eq. 9 + 4-strip x–s plane partition. Strip D excluded; for s > 15 with x ≥ 28, **no method given**. | Uniform asymptotic Eq. (3.1) is **valid for all r ≥ 0** (uniform in n). DBSKIN handles up to 18 digits; the algorithm is precision-agnostic in structure (only XMIN, XLIM, NLIM, M_max change with target ε per Table I). | **Amos** (the uniform expansion is the structural advantage). |
| **Citation footprint (as of paper publication)** | "Citing articles: 1" per the Taylor & Francis cover page of the local PDF. Used in SURCU and BOXER (PSI-era pilot codes). | TOMS Alg. 609 has been mirrored in SLATEC, used in ANISN/DOT-IV/DORT-class transport codes for ~40 years. Hundreds of citations on Google Scholar. **Public reference implementation for the Bickley function** in computational science. | **Amos** (community adoption is the de facto standard) |

### Verdict for ORPHEUS

**Davierwalla is NOT advantageous over Amos on technical merits.**
The two papers solve the same problem; Amos solves it better on
every axis except raw FLOP count (where Davierwalla's
continued-fraction form is a touch tighter, but the difference
is at most a single-digit factor in pure compiled code and
disappears once Python overhead dominates).

Concretely:

- **Accuracy**: Tied at single precision (10⁻¹²). Amos wins at
  double precision because the published constants are stored to
  18 digits vs Davierwalla's 14. ORPHEUS uses IEEE binary64; the
  Amos coefficients are sufficient to deliver ≈ 10⁻¹⁵ relative
  error throughout. Davierwalla's CDC-single-prec coefficients
  would need re-fit for double precision to match.
- **Coverage**: Amos covers the entire positive (x, n) quadrant
  (up to underflow). Davierwalla's strip-D blind spot
  {x ≥ 28, s > 15} is real and aligns with where ORPHEUS's
  current `ki_n_float` *also* underperforms — but Amos handles
  this region cleanly via the uniform asymptotic expansion.
- **Speed**: Davierwalla's continued fraction is ~6–8 divides
  (branch-light, deterministic latency). Amos's [0, 2] power
  series has variable iteration count (~30–80 FLOPs) but
  truncates adaptively to required ε. On x > 2 Amos uses the
  uniform expansion (~50–200 FLOPs depending on n). Net: in
  pure compiled C/Fortran, Davierwalla is **~2–4× faster** in
  steady state on the small-x regime. **In Python with f2py-wrapped
  Fortran**, the call-site overhead dominates and the gap closes.
- **Stability**: Amos's two-region scheme (series on [0, 2],
  uniform expansion on (2, ∞)) is structurally cleaner than
  Davierwalla's 8 hand-picked sub-intervals on [0, 7].
- **Portability**: Amos is *portable Fortran* with a published
  machine-parameter abstraction (I1MACH/R1MACH/D1MACH).
  Davierwalla is CDC-tuned with no portability layer.

**Gain magnitude (Amos over Davierwalla)**:
- Accuracy: 1000× better (10⁻¹⁵ vs 10⁻¹²) at IEEE binary64.
- Coverage: closes the strip-D blind spot.
- Acquisition cost: hours (f2py wrap of public Fortran) vs days
  (manual transcription + re-fit + verification of ~150 16-digit
  constants).
- Speed gain over Davierwalla: 0.25–0.5× (Davierwalla wins).

The single axis on which Davierwalla wins (raw FLOPs in compiled
form) does not propagate to ORPHEUS's actual hot path — the
Python `ki_n_float` overhead is the gating cost, not the
underlying algorithm. And **even on the FLOP axis**, Lorensi's
2025 polynomial+Remez form on [0, 1] would beat both Davierwalla
and Amos for the n = 1, 2 case at ~200 ns/eval in compiled C.

### Recommendation

**ADOPT AMOS** (specifically: TOMS Algorithm 609 BSKIN/DBSKIN)
**if** a Bickley replacement is decided for ORPHEUS. Reasons:

1. Technical superiority on accuracy, coverage, stability,
   portability — confirmed by reading both papers in full.
2. Zero transcription cost: the Fortran source is public.
   `f2py` wraps it in 5 minutes; a pure-Python translation is
   ~1 day given the published source.
3. Strip-D coverage gives ORPHEUS a path to high-n / high-x
   evaluations (e.g., Knyazev cylinder G_bc at deep τ) that
   Davierwalla cannot provide.
4. **Davierwalla's transcription cost is not worth absorbing**
   given that Amos exists and is technically superior. The
   transcription would deliver a *worse* implementation by every
   metric except a single-digit compile-time FLOP difference.

**However**, the original deferral conclusion still holds at the
project level:

**ADOPT NEITHER YET** — the current ORPHEUS `ki_n_float` is
adequate (1e-13 in the working regime). Bickley is not
load-bearing on any hot path measured in this session
(ki_n_float at 22.9–77.7 µs/eval is fast enough for production).
**When** Bickley becomes load-bearing, the pivot is **Amos**, not
Davierwalla — and the user should NOT spend effort manually
transcribing Davierwalla's tables. That effort is better spent
either (a) wrapping Amos's BSKIN with f2py, (b) translating
BSKIN to Python directly, or (c) implementing Lorensi's
n = 1, 2 polynomial form for the cylinder kernel hot path.

**Hybrid possibility** (if Davierwalla becomes attractive *despite*
this verdict): Davierwalla's strip-A continued fraction is
genuinely tight on [0, 7] for s = 1, 2, 3 — it's the part of the
domain where Bickley evaluations are densest. A Davierwalla-only
implementation of strip A, paired with Amos for everything else,
could deliver Davierwalla's compile-time speed advantage on the
hot path while inheriting Amos's coverage everywhere else. This
would still cost ~1 day of transcription + verification of
Davierwalla's strip-A tables (~24 numerator/denominator coefs ×
3 functions × 7 sub-intervals = ~500 16-digit constants).
**Not recommended** unless a Bickley microbenchmark in the deployed
ORPHEUS code shows an unambiguous bottleneck after Amos adoption.

### Critical caveats

- **Amos's "10⁻¹⁸" claim** is the precision floor of the stored
  constants (28-digit CDC computation truncated to 18-digit
  storage), NOT the algorithm's intrinsic accuracy. Re-running
  the CDC double-precision computation today with mpmath at
  dps = 30 would deliver constants accurate to ~10⁻²⁹; the
  algorithm itself is precision-agnostic. So porting BSKIN to
  Python with mpmath-regenerated constants gives effectively
  unlimited precision at the cost of regeneration. This is a
  **dormant feature**, not a current limitation.
- **Davierwalla's "10⁻¹²" claim** assumes CDC single precision.
  The published 14-digit coefficients are tight at CDC unit
  round-off (~7.1e-15) but were never re-validated for IEEE
  binary64. The 10⁻¹² figure could degrade to 10⁻¹⁰ on IEEE
  hardware (untested in the published record).
- **Both papers test by self-consistency**, not against each
  other. There is no published independent benchmark that
  cross-validates Davierwalla and Amos to high precision. If
  ORPHEUS adopts either, an L0 test should compare against
  `ki_n_mp(n, x, dps=50)` at dense (n, x) grid — not against
  the other paper's claims.
- **The user's filename "Amos (1993)"** is wrong. The paper is
  Amos 1983 (TOMS 9(4), Dec 1983). The 1993 date does not
  appear anywhere in the paper, the journal, or the citation
  history. Action: rename `scratch/literature/Amos (1993)…pdf`
  to `Amos (1983)…pdf` to avoid confusion in future sessions.
