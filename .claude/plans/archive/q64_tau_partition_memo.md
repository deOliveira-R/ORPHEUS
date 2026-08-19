# Q5.6.4 — the cylinder τ partition: what was tested, what was assumed, what it showed

**Status: LANDED ON THE BRANCH, ONE TEST RED, THE DESIGN'S STRUCTURAL PREMISE
REFUTED.** Written 2026-08-11 as a compaction-proof handoff. Self-contained on
purpose: assume the reader has NO memory of the session that produced it.

Companion documents (all committed):
* `.claude/plans/archive/quadrature_machinery_campaign.md` — the campaign plan; its
  Q5.6.4 item carries the same measurements inline with the refutations.
* `scratch/q64_red_triage.md` — the per-test red triage, §D is the open red.
* `scratch/q64_tau_edge_convention_literature.md` — the literature deliverable.
  **Line 943 is the sentence that refutes the design.**
* `scratch/q64_absorber_retirement_blast_radius.md` — the retirement audit.
* `orpheus/derivations/discrete/sn/angular_differencing.py` — the analysis
  module: the P0–P4 predicate ladder, ν-closure, both β's, the τ/β nomenclature.

---

> ✅✅ **CLOSED 2026-08-11 ON THE PRIMARY SOURCE — READ §9bis.12 FIRST.**
> Morel & Montry 1984's cylinder appendix (A2)–(A4) **IS** what the carve
> implements: `[M]` production reproduces all four properties their appendix
> states, including `τ₁ → ¼`, `τ_M → ¾`, sinθ-independence, and
> `(τ₁−½)/w = −M/4`. The endpoint behaviour earlier layers diagnosed as the
> defect is the primary's own designed value, and the BMC-43-vs-R&L-16 tension
> is the trade M&M **announce** (printed 616–617): R&L buy second-order angle
> by moving the ordinates; M&M buy diffusion-limit consistency on an arbitrary
> quadrature and pay in order. **The design question is settled; only the SEED
> remains** (§9bis.12's last section quantifies the lever).
>
> ⛔⛔⛔ **THIS MEMO HAS FOUR LAYERS. §9bis.10 was the verdict; §9bis.12 is the
> primary-source confirmation AND one correction to §9bis.10's own reasoning.**
> * **§1–§9** = attempt 1 (2026-08-11, pre-compaction). Its diagnosis is refuted.
> * **§9bis.1–.9** = attempt 2, which refuted attempt 1 and then proposed
>   re-posing τ into ω. **That proposal is ALSO refuted.**
> * **§9bis.10** = the verdict: **the landed carve is CORRECT and stays.** It is
>   the first time the cylinder satisfies BMC Eq. (43) exactly. The accuracy
>   regression is the recurrence's transient error amplification, which the
>   retired `[½,1]` absorber had been masking by forfeiting Eq. (43).
>
> Every **NUMBER** in §4 and §9bis stands. What moved, twice, is the
> interpretation. Refuted claims are preserved in place with markers (§2's table
> and each section's banner are the map) — per `plan-authoring` §3, because the
> two wrong frames are what a fresh reader would otherwise re-derive.
>
> ⛔ **Settled, do not reopen:** candidate (3) is not a definition; τ is NOT
> re-posed into ω; the absorber is NOT restored; Hébert is NOT cited against BMC.

## 0a. ▶▶ WHERE TO RESUME — the SEED, and nothing about τ

**Everything about the τ partition is finished.** Landed, literature-ratified,
gated, citations repaired. `[M]` commits on `refactor/operator-strategy-layers`:
`3dda18ca` (carve) · `d5067c4d` (analysis module) · `39b46a31` (snapshots) ·
`c9bb61b4`+`c33178ef` (corpus + gates) · memo layers `1a0bad08` · `8db88596` ·
`4e1c6090` · `3ae4ac70` · `b77ce7bd` · `6cd1d20e` · pointers `9412ee81`.
`main` is untouched. Verify any of these with
`git merge-base --is-ancestor <hash> HEAD`.

### The open problem, stated in the domain's terms

A curvilinear level's angular march needs a starting value at the most-inward
direction (`μ = −1` sphere, `ω = π` per cylinder level). **That seed's error
does not decay** — it is carried through the whole level by the M-M recurrence,
and it is the last unexamined term in the cylinder's accuracy.

### Why it is the seed and not the closure — the evidence, in one place

* **All four accuracy rankings track the recurrence's transient amplification
  and are INVERTED against closure fidelity** (§9bis.8c). The closure is right;
  something feeding it is not.
* `[M]` **the τ → ¼ endpoint behaviour is Morel & Montry's own designed value**
  (§9bis.12) — `τ₁ → ¼`, `τ_M → ¾`, `(τ₁−½)/w = −M/4`, all four reproduced.
  Not a defect to fix.
* `[M]` **ψ̂'s sign is a property of SEED CONSISTENCY** (§9bis.5 correction):
  production seed ⟹ ψ̂ positive; zero seed ⟹ negative, bounded by `A(M)`.
  The seed is already demonstrably the controlling term.

### ⭐ THE LEVER IS QUANTIFIED AND TWO-SIDED — this is the whole design constraint

`[M]` M&M Eq. (14) admits an **arbitrary edge-march seed `μ_s`**, and for
Gauss-S₂ diamond **`β(μ_s) = μ_s + 2/√3`** — linear, **unit gain, no damping**,
and **sign-determining**:

| seed | β | consequence |
|---|---|---|
| **under**estimate | `β < 0` | **the flux dip** |
| exact | `β = 0` | the M-M property |
| **over**estimate | `β > 0` | no dip, **but the diffusion limit is LOST** |

⟹ **there is a two-sided target, not a "make it more accurate" direction.** Any
sweep of the seed must report β's SIGN, not just an error norm. Their own
diagnostic is the *effective starting cosine* (their Eq. 10).

### The concrete candidate, from the primary

`[M]` Reed & Lathrop printed p. 239: **α is known at BOTH ends** —
`α_{1/2} = 0` when neutrons flow out of but not into the first cell, and
`α_{M+1/2} = 0` at the last — *"either `α_{m+1/2}` or `α_{m−1/2}` is known,
depending on the direction in which the equations are being solved."* R&L march
from either end. **Production marches ONCE, one-sided**
(`orpheus/sn/sweep/psi_half_angle_seed.py`,
`MorelMontryAngularSweep.psi_half_seed`, route (a) since #282).
⟹ **the experiment is a two-ended march against the one-sided one.**
⚠ M&M say **nothing** about march direction or conditioning; that claim is
R&L's alone, and "ill-conditioned" is *their* framing, not a measurement of
ours. Measure it before repeating it.

### The instruments that now exist and WORK (three of the old five do not)

| use this | for | where |
|---|---|---|
| **R&L Eqs. 15/16** — `(τ_m−½)/w_m` bounded ⟺ node is the μ-midpoint to `O(w²)` | second-order angular accuracy. **POINTWISE**, so the σ_y fold cannot annihilate it | §9bis.11 |
| **`A(M) = max_m Π(1−τ)/τ`** (`2.41 … 9.44`, M = 2…32) | the transient amplification bound | `test_psi_half_positivity.py` |
| **the seed-regime pair** (marched vs zero) | ψ̂ sign | same module |
| MMS aniso-cyl at **`nx = 320`** | truncation order in `n_φ`, range to `9e-6` | §4.5 |

⛔ **DEAD instruments — do not reach for them:** P1 `c = Σwη²`, BMC β and
Lathrop β are **τ-blind** (bit-identical under garbage τ, §9bis.4); ν-closure
only reports which CHART you ran it in (§9bis.3); and the trajectory-resolvent
cross-check is **REFERENCE-limited** at `≈3e-2` — refining the SN side moves it
the WRONG way (§9bis.8b).

### ⛔ Settled — re-opening any of these is re-doing measured work

1. ⛔⛔ **REFUTED 2026-08-12 — was:** *"τ is **not** re-posed into ω
   (that is BMC's *diamond*, leading-order only)."* The ruling conflated
   two orderings: BMC's "leading-order" is about the **diffusion limit**,
   not about **truncation order**. `[M]` Lathrop (2000) NSE 134 Eq. (30),
   LOCAL and scan-verified — with δ = 2τ−1 the truncation is
   `O(δΔμ + Δμ²)` and *only* δ = 0 (the midpoint node) gives `O(Δμ²)`;
   every weighted diamond is FIRST order. `[M]` our ω-edges are exactly
   equispaced (4.4e-16) and each node is exactly its cell's ω-midpoint
   (1.6e-15), so τ-in-ω would be exactly ½ — yet τ-in-η reads [0.25, 0.75]
   and `(τ−½)/w` is unbounded (= M/4). `[M]` at spatially-converged
   nx=320 the measured angular orders are **1.83 (shipped) vs 2.88
   (τ≡½)**, an 8× accuracy gap at n_φ=64. §9bis.3 had already named the
   mechanism ("the partition is read in ω, the WEIGHT in η") and then
   ruled the fix out. ⛔⛔ **AND THAT REFUTATION IS ITSELF REFUTED, same day.** The
   accuracy inference rested on `build_cylindrical_anisotropic_mms_case`,
   whose exact solution `A(r) + B(r)η` (sn.py:3803) lies **in the
   closure's kernel** — τ is the barycentric coordinate, so the closure
   is exact on `span{1, η}` by construction. On a fixture that actually
   EXCITES the closure (`CylHarmonicMMS`, h₁+h₂), `[M]` the SHIPPED τ
   **wins** by 34 % / 12 % at n_φ = 16 / 32 and ties at 64, spatially
   converged (nx = 80/160/320 agree to 3 digits). **The original ruling
   — do not re-pose τ into ω — STANDS, now on measurement as well as on
   Morel & Montry's appendix.** Lathrop Eq. (30) remains true about the
   EQUATION's truncation order; it is not a ranking of solution accuracy.
   Full chain + both fixtures:
   `.claude/plans/issue_235_angular_accuracy_campaign.md` §4bis. Full chain:
   `.claude/plans/issue_235_angular_accuracy_campaign.md` §4bis.
2. Candidate (3), "τ on α's recursion-defined edges", is **not a definition**.
3. The `[½,1]` absorber is **not** restored — 5 of 5 primaries prescribe no
   limiter, and the interval is Grant's on the **spatial** weight.
4. Hébert is **never** cited against BMC, nor as the source of any τ.
5. The `tol_per_cell` on the resolvent gate is **not** relaxed, and that gate is
   **not** re-baselined.

### Open, tracked elsewhere

* **#352** — the resolvent fixture is a starved solve (`max_inner=300`, needs
  ≈528); the re-baseline must ride whatever moves that snapshot.
* The resolvent L1 gate is reference-limited and its docstring does not say so.
* `bailey-dome-recursion` still encodes a refuted attribution (2 definitions,
  6 `:eq:` consumers, a `vv-status`, a V&V-matrix row) — its own pass.
* Not local, not read: **Grant 1968** (need the TITLE for a bib entry),
  **Alcouffe & O'Dell 1986** (user has a scan awaiting OCR), **Lathrop &
  Carlson 1966**, and the lead **Dudziak/O'Dell/Alcouffe LA-7911-PR (1979)**.

---

## 0. ⭐ READ THIS FIRST — the one-paragraph state

The cylinder's Morel–Montry angular closure weight τ was re-posed from a
**chord-midpoint-in-η** cell partition to a **midpoint-in-ω** one, and the
`[½,1]` "absorber" (a clamp) was retired. The change is landed on branch
`refactor/operator-strategy-layers` across five commits. It is **structurally
tidy and empirically worse**: `[M]` an independent-method flux-shape
cross-check moved from `6.593e-02` to `1.268e-01` against a `1e-01` tolerance
(1.92× worse, now RED), and the anisotropic-cylinder MMS floor is ~1.8–2×
worse at `n_φ ≥ 16`. Worse, **the structural argument that justified the change
is false**: it assumed α lives at the *geometric* arc half-angle boundary, and
α actually lives at a **recursion-defined** edge that no source identifies with
the geometric one. A third, unmeasured candidate — τ built on α's own
recursion-defined edges — is the partition the campaign's stated goal (α and τ
share ONE object) actually points at. **Do not close Q5.6.4 and do not relax
the tolerance until that third candidate is measured.**

`main` is unaffected. All five commits are branch-only, so revert is cheap.

---

## 1. What Q5.6.4 was originally supposed to be, and why that died

The campaign plan declared 6.4's acceptance test as:

> *"the #229 azimuthal floor — today flat at ≈1.9e-2 on the anisotropic
> curvilinear MMS with no convergence order. It must fall and recover an order.
> … This is the acceptance test."*

⛔ **REFUTED before any code was written**, on two independent grounds:

1. The number was stale. `1.9e-2` was the pre-6.3 fixture; the gate's own
   docstring already carried `3.538e-3 → 6.782e-4` post-fold.
2. The claim was false. `[M]` retiring the absorber **alone** makes the floor
   **1.8–3.4× WORSE** at every rung (§4.1). No τ choice makes "the floor falls
   and an order comes back" true, because (§4.3) the "floor" measured at
   `nx=80` is substantially the **spatial** error, not an angular one.

So 6.4 was re-scoped mid-session to a structural goal: *a polar level's
azimuthal cells are ONE partition, and every coefficient that references a cell
boundary reads it from there.* **That goal is still right.** What is wrong is
the partition that was chosen to realise it.

---

## 2. THE PREMISE CHAIN — every claim, and its status

This is the most important section. The design was reached by a chain; one link
is false and it is not the last one.

| # | premise | status |
|---|---|---|
| **P-A** | α and τ both reference "the boundary between azimuthal cell m and m+1", derive it independently, and DISAGREE | ⛔ **REFUTED as a DUPLICATION claim, 2026-08-11 (§9bis.2)** — α derives no boundary at all: it derives edge-INDEXED amplitudes from `(w, η)`. There was no twin derivation. (Was: ✅ TRUE, a genuine Cardinal-Rule-2 defect, `[M]` §4.2.) |
| **P-B** | α lives at the **geometric** arc half-angle boundary `ω_{k−1/2}` | ⛔ **the REFUTATION is itself refuted, 2026-08-11** — `[M]` §9bis.1: α ∝ `ξ(e_arc)` to `1e-16`; κ is an AMPLITUDE, not a displacement. TRUE in shape. (Was: ⛔ FALSE, "the fatal link" — see §3, now superseded.) |
| **P-C** | ∴ τ must be built on the geometric arc partition | ✅ **REINSTATED 2026-08-11** on a different ground: not "α lives at the arc edge" but "α's *shape* is the arc-edge dome". The PARTITION was right. (Was: ⛔ INVALID.) |
| **P-D** | the `[½,1]` absorber exists to compensate the chord partition's stretched end cells | ⛔ **REFUTED 2026-08-11 (§9bis.5)** — it is a **non-amplification guard**: `τ ≥ ½ ⟺ \|(1−τ)/τ\| ≤ 1`. `[M]` worst running product `1.000000` clamped vs `40.7` (chord) / `9.44` (arc) at `n_φ = 64`. |
| **P-E** | ν-closure `= 1.000000` validates the arc partition | ⛔ **INVALID AS VALIDATION**, and worse than §5.3 said — `[M]` §9bis.3 it is **CHART-RELATIVE**: it reports which variable you ran it in, nothing else. |
| **P-F** | `τ ≡ ½` is the principled M-M value expressed in ω | ⛔ **the REFUTATION is itself refuted, 2026-08-11** — `[M]` §9bis.3: `τ ≡ ½` closes ν-closure **exactly** in ω at every order, and it is P2 in the MARCH variable (§9bis.7), not a foreign method. (Was: ⛔ REFUTED — Hébert's diamond, fails ν-closure by 16.5 %; that 16.5 % was the η reading.) |
| **P-G** | BMC's β is the diffusion-limit oracle that adjudicates a partition | ⛔ **STRONGER THAN §5.2 SAID** — `[M]` §9bis.4 BMC β takes τ **nowhere**; it is a function of `(quadrature, partition)` only, so it cannot grade a τ on ANY fixture, folded or not. Same for P1 `c` and Lathrop β. |
| **P-H** | the literature's cumulative-weight edges can be transplanted to the cylinder | ⛔ **REFUTED** — violates P3, worsening with refinement. §4.4. *(unchanged)* |

⛔⛔ **THE 2026-08-11 RE-READING SUPERSEDES THE PARAGRAPH BELOW.** The defect is
NOT "α and τ disagree about the partition" — they agree (P-C). It is that the
**partition is read in ω and the WEIGHT in η**: two charts on one level
(§9bis.3). Kept for the record because the wrong frame is what produced §8:

> ⭐ **P-A survived and P-C did not.** The defect is real; the fix was aimed at
> the wrong target. A fresh session should NOT conclude "the chord partition was
> fine all along" — it should conclude "α and τ must share a partition, and we
> have not yet identified which one α actually uses."

---

## 3. ⛔ HOW P-B IS FALSE — the sentence that refutes the design
## ⛔⛔ **THIS WHOLE SECTION IS REFUTED — 2026-08-11, `[M]` §9bis.1/§9bis.2.**
## κ is an AMPLITUDE defect of the quadrature, not a displacement of the edge:
## α ∝ ξ(e_arc) to `1e-16`, so α DOES reference the arc partition. The section
## is kept verbatim (plan-authoring §3) because it is the argument a fresh
## reader would otherwise re-derive, and because its own ⚠ process note below
## is the lesson that repeated itself.

The literature deliverable, `scratch/q64_tau_edge_convention_literature.md:943`:

> *"…the **recursion-defined edge** — which is what α *is*, by definition, in
> both sources — and the geometric arc-half-angle edge, **which no source
> uses**."*

The reasoning that makes this decisive:

* The campaign's own theorem **T3** (gated, `tests/sn/sweep/curvilinear/test_alpha_closed_form.py`)
  states the production α has an exact closed form
  `α_k = −w_gl · κ · [ξ(ω_{k−1/2}) − ξ(ω_{−1/2})]`, with
  `κ = Δω / (2 sin(Δω/2))`, `[M]` κ − 1 = 2.6 % at `n_φ = 8`.
* The design read this as *"α equals ξ at the arc boundaries (up to a small
  factor), therefore α's partition IS the arc partition."*
* **That is backwards.** κ ≠ 1 is exactly the statement that α is **not** at the
  geometric arc boundary. α is defined by its own recursion
  (`α_{m+1/2} = α_{m−1/2} − w_m η_m`, Hébert 3.397–3.399 after Alcouffe &
  O'Dell); the edge that recursion implies is κ-scaled relative to the
  geometric arc edge. T3's closed form is a *theorem about our rule* relating
  the two, not an identification of them.

⟹ **α's partition is the recursion-defined one.** Unifying τ onto the geometric
arc partition unified it onto a partition α does not use, i.e. it did not
achieve P-A's goal at all — it replaced one mismatch with a different one.

⭐ **And this re-reads P-D.** If the chord partition happens to sit closer to
α's recursion-defined edges than the geometric arc does, then the chord
partition was approximately *right*, and the `[½,1]` absorber was patching its
**endpoint pathology** (`τ_raw,0 = 0` bit-exactly when a node sits on the arc
start, which divides by zero in the recurrence) rather than compensating a wrong
partition. That is the *opposite* of what the session concluded. **It is not yet
measured either way.**

⚠ **Process note worth carrying forward:** the refuting sentence was in a
sub-agent report that was read, quoted, and recorded — and the design it
invalidated was built anyway, because the correction ("κ is not an error in α")
was filed as a detail rather than propagated to the conclusion resting on it.

---

## 4. WHAT WAS MEASURED — every instrument, with configuration

All probes are in `$CLAUDE_JOB_DIR/tmp/` (**ephemeral — the job dir is cleaned
up; re-create from the descriptions below if needed**). Host env,
`.venv/bin/python`. ⚠ Do NOT run probes under `python -O` if they use bare
`assert` — `-O` strips them (lesson L42).

### 4.1 The MMS floor, absorber live vs retired — `nx = 80`

Anisotropic cylindrical MMS, volume-weighted L2 of the scalar flux,
`build_cylindrical_anisotropic_mms_case(n_phi=·)`, `nx=80`, `max_inner=500`,
`inner_tol=1e-13`. The absorber-live row **reproduces the shipped gate's own
docstring to every printed digit**, so the instrument is sound.

| τ convention | `n_φ`=8 | 16 | 32 |
|---|---|---|---|
| chord-midpoint + `[½,1]` absorber (pre-carve production) | `3.5384e-3` | `6.7824e-4` | `2.4837e-4` |
| chord-midpoint, unclamped (the naive retirement) | `6.2244e-3` | `2.3020e-3` | `6.0065e-4` |
| geometric arc ω half-angle, unclamped | `3.1503e-3` | `1.1326e-3` | `3.2611e-4` |
| cumulative-WEIGHT edges (the sphere's) | `1.9252e-2` | **diverges (nan)** | **diverges (nan)** |
| `τ ≡ ½` (angular diamond) | `3.4055e-3` | `3.7443e-4` | `1.3279e-4` |

### 4.2 The α/τ partition disagreement — the defect P-A names

`folded_product(n_mu=4, n_φ)`, level 0. "ω-width spread" = the spread of the
ω-widths implied by pushing the chord partition's η-edges back through `arccos`,
against a quadrature whose weights are **bit-exactly equal**.

| `n_φ` | max rel η disagreement | implied ω-width spread | quadrature's own weight spread |
|---|---|---|---|
| 8 | `5.3825e-2` | 18.71 % | `0.00e+00` |
| 16 | `1.7752e-2` | 17.59 % | `0.00e+00` |
| 32 | `4.7227e-3` | 17.48 % | `0.00e+00` |
| 64 | `1.1987e-3` | **17.46 %** | `0.00e+00` |

Algebra, verified to `3.3e-16`: for INTERIOR edges
`chord = cos(Δω/2) × arc` (from `½[cos ω_a + cos ω_b] = cos(½(ω_a+ω_b))·cos(½Δω)`),
while the ENDPOINTS are pinned at `∓sinθ` **unscaled** — so the outermost cells
stretch to absorb the shrink. The η error → 0 as `Δω → 0`; the ω-width spread
does **not**, converging to ≈17.45 %.

⚠ **This measurement is still valid. It is the interpretation that changed:**
a mismatch between τ's cells and the *quadrature's* cells is only a defect if
α's cells are the quadrature's — and per §3 they are not.

### 4.3 ⚠ THE CONFOUND — `nx = 80` does not isolate the angular floor

HEAD convention, sweeping the MESH at fixed `n_φ`:

| `n_φ` | nx=40 | 80 | 160 | 320 | spatial orders |
|---|---|---|---|---|---|
| 32 | `6.1936e-4` | `2.4837e-4` | `1.6358e-4` | `1.4449e-4` | 1.32 · 0.60 · 0.18 |
| 64 | `5.3570e-4` | `1.5503e-4` | `6.1967e-5` | `4.0491e-5` | 1.79 · 1.32 · 0.61 |
| 128 | `5.1562e-4` | `1.3397e-4` | `3.8772e-5` | `1.5488e-5` | 1.94 · 1.79 · 1.32 |

At `n_φ=128`, refining `nx` 80→320 still drops the error **8.6×** ⟹ the
`≈1.3e-4` that every convention "saturated" to at `nx=80` is the **MESH**.
The shipped gate's `n_φ` 8→16 leg IS angular-dominated (`3.5e-3 ≫ 1.3e-4`) and
remains sound; anything read at `n_φ ≥ 32, nx = 80` is mixed.

⭐ Consequence: **"the cylinder has no O(h²) window at any practical
quadrature" is over-general.** At `nx=320` the angular error converges cleanly
at ~O(n_φ⁻²) with **no flat floor in range** (§4.5).

### 4.4 P3 admissibility — why the sphere's convention cannot transplant

BMC Eq. 52's cumulative-weight edges (correctly renormalised to
`Σ w̄ = 2 sinθ`), `folded_product(n_mu=4)`, level 0. Lathrop's admissibility
condition is that the ordinate lie INSIDE its own cell:

| `n_φ` | M | min `τ` | max `τ` | ordinates OUTSIDE their own cell |
|---|---|---|---|---|
| 8 | 4 | 0.1522 | 0.8478 | 0 — admissible |
| 16 | 8 | −0.3259 | 1.3259 | **4/8** |
| 32 | 16 | −1.1841 | 2.1841 | **12/16** |
| 64 | 32 | −2.8552 | 3.8552 | **28/32** |

Mechanism: an arc cell's η-measure is
`2 sinθ · sin ω_m · sin(Δω/2) ∝ sin ω_m`, **not constant**, while a trapezoid
weight is — and the mismatch **WIDENS** with refinement (cell η-measure ÷ its
mean): `[0.586,1.414]` / `[0.305,1.531]` / `[0.154,1.561]` / `[0.077,1.568]` at
`n_φ = 8/16/32/64`. ⟹ BMC Eq. 52 is not a law; it is the statement that *in
their* quadrature the weight equals the cell's η-measure. This is the
`derivations/.../angular_differencing.py` inherited comment, sharpened:
*"weights are uniform in φ-space, not η-space."*

### 4.5 The ANGULAR-ISOLATED comparison — `nx = 320`

Spatial contribution `≤ 1.5e-5`, so `1e-4 … 3e-3` is angular.

| τ convention | 8 | 16 | 32 | 64 | orders in `n_φ` |
|---|---|---|---|---|---|
| chord + absorber (pre-carve) | `3.5111e-3` | `5.8890e-4` | `1.4449e-4` | `4.0491e-5` | 2.58 · 2.03 · 1.84 |
| chord, unclamped | `6.2063e-3` | `2.2824e-3` | `5.7024e-4` | `1.4175e-4` | 1.44 · 2.00 · 2.01 |
| geometric arc (LANDED) | `3.1281e-3` | `1.1078e-3` | `2.8285e-4` | `7.1658e-5` | 1.50 · 1.97 · 1.98 |
| `τ ≡ ½` | `3.4258e-3` | `3.4907e-4` | `3.9321e-5` | `9.1485e-6` | 3.29 · 3.15 · 2.10 |

⚠ `τ ≡ ½` wins here **and is still wrong** (P-F): it fails ν-closure by 16.5 %
and is a different published method. This table measures **truncation order**,
which is precisely what `τ=½` optimises — a Mode-12 blindness of the
instrument, not a verdict.

### 4.6 ⭐⭐ THE DECISIVE MEASUREMENT — an INDEPENDENT method

`tests/sn/verification/analytical/test_phase_c_crosscheck.py::test_phase_e_trajectory_resolvent_flux_shape_crosscheck[cyl_2g_3reg_folded_4x8_dd_n40]`.
Compares the SN cylinder flux SHAPE (L∞-normalised per group) against the
**trajectory-resolvent** reference — a structurally different method. BOTH
snapshots computed against the SAME reference:

| τ convention | per-group max `|Δφ_norm|` | max | vs `tol = 1e-01` |
|---|---|---|---|
| OLD chord + `[½,1]` absorber | `[0.06593, 0.01335]` | **`6.593e-02`** | PASS, comfortable |
| NEW geometric arc, no clamp | `[0.12676, 0.01341]` | **`1.268e-01`** | **RED** |

⟹ **1.92× WORSE.** The old margin was 6.6 %, not marginal. And note group 2 is
essentially unchanged (`0.01335 → 0.01341`) — the whole regression is in
group 1.

⚠ Caveat, so this is not over-read: the resolvent reference carries its own
discretisation error, this is a SHAPE comparison after L∞ normalisation on a
heterogeneous 3-region closed cylinder, and the bar is a loose `1e-1`. It is
not a gold standard. **But it is independent, and 1.92× is not noise.**

### 4.7 Snapshot / k_eff movement

`[M]` of the FOUR cylinder artifacts only TWO moved, and only ONE is a genuine
angular-closure catcher:

| artifact | moved? | τ-sensitivity (flux, under a deliberate `τ := 0.7`) |
|---|---|---|
| `cyl_2g_3reg_folded_4x8_dd_n40` | **yes** (k_eff `1.2302082296342958 → 1.2310212585879858`) | `8.78e-02` |
| `walk_matvec_cyl_2g` | **yes** | — |
| `cyl_1g_homogeneous_folded_4x8_dd_n20` | no | `1.14e-10` |
| `cyl_1g_homogeneous_folded_2x4_dd_n20` | no (τ unchanged, M=2) | `1.07e-11` |

Two independent blindnesses, both worth knowing:
* **M = 2 is degenerate**: `[M]` the ω-midpoint and chord partitions are
  **BIT-IDENTICAL** at M = 2 (the single interior chord midpoint
  `(η₀+η₁)/2 = 0` IS the arc edge at `ω = π/2`), diverging only from M = 3
  (`3.17e-2` / `4.46e-2` / `6.82e-2` at M = 3/4/6). **No `folded_product(·,4)`
  fixture can ever see a partition change.**
* **1-group homogeneous is degenerate**: its `k_eff = 1.5` **is**
  `k_∞ = νΣ_f/Σ_a`, flux-shape independent ⟹ τ-independent by construction
  (`vv-principles` anti-pattern #3).

---

## 5. THE INSTRUMENTS — what each can and cannot see

Carry this table forward; three of these five mislead.

| instrument | what it measures | ⚠ blind spot |
|---|---|---|
| MMS aniso-cylinder L2 | truncation order | **contaminated by the MESH at `nx=80`** (§4.3); and optimised by `τ=½`, so it cannot adjudicate a closure |
| **BMC β** (`contamination_beta`) | the diffusion-limit contamination | ⛔ **IDENTICALLY ZERO on a σ_y-folded arc for ANY antisymmetric edge set** — §5.2 |
| **ν-closure** (`nu_closure_residual`) | "did this τ come from a contiguous partition?" | ⛔ **cannot say WHICH partition** — §5.3 |
| trajectory-resolvent shape | agreement with an INDEPENDENT method | own discretisation error; loose bar; shape-only |
| P3 admissibility | ordinate inside its own cell | a *necessary* condition; several partitions satisfy it |

### 5.2 ⛔ β is a symmetry identity on the folded arc — verified with garbage

β = 0 is the entire reason the Morel–Montry τ exists, so it looked like the
right adjudicator. It is not, here. `folded_product(4,16)` level 0:

| edge set fed to β | β |
|---|---|
| production | `+6.94e-18` |
| GARBAGE — edges scaled 0.5× | `+3.47e-18` |
| GARBAGE — edges **CUBED** | `+1.73e-18` |
| GARBAGE — **random**, antisymmetrised | `−3.47e-18` |
| one edge nudged (breaks antisymmetry) | **`−3.53e-03`** |

Proof: the fold makes nodes antisymmetric (`max|η + η[::-1]| = 0.000e+00`) and
the α dome symmetric (`2.78e-17`), so
`term_{M−1−m} = (−η_m)(α_m(−e_m) − α_{m+1}(−e_{m+1})) = −term_m` — the sum
cancels pairwise for ANY antisymmetric edges. β sees only antisymmetry.
**It also certified the cumulative-weight convention that DIVERGES the solve**
(§4.4), which is how the blindness was caught.

⚠ Two β's exist and they are near-opposites — see the nomenclature table in
`angular_differencing.py`. BMC's β is a SCALAR, zero iff τ is Morel–Montry.
Lathrop's β is a SEQUENCE (α's pointwise defect), zero iff `τ ≡ ½`.

### 5.3 ⛔ WHY ν-CLOSURE DID NOT VALIDATE THE DESIGN

ν-closure runs BMC Eq. 43 forward: `ν_{1/2} = η_start`,
`ν_{m+1/2} = (η_m − (1−τ_m) ν_{m−1/2}) / τ_m`, and asks whether it lands on the
level's far endpoint. `[M]` `ν/sinθ` at close:

| `n_φ` | geometric arc | chord | `[½,1]` clamp | `τ ≡ ½` |
|---|---|---|---|---|
| 8 | `1.000000` | `1.000000` | `1.016389` | `1.164784` |
| 16 | `1.000000` | `1.000000` | `1.001930` | `1.039182` |
| 64 | `1.000000` | `1.000000` | `1.000030` | `1.002412` |

**Read this correctly.** It proves the clamp and `τ≡½` correspond to **no**
partition — a genuine and useful result, and the honest principled condemnation
of the absorber. It does **NOT** distinguish arc from chord: *both* close
exactly, because both ARE contiguous partitions. Treating `1.000000` as
validation of the arc choice was the error; it is a necessary condition that
the rejected alternative also satisfies.

---

## 6. THE THREE CANDIDATE PARTITIONS — defined for re-implementation

For a polar level with radial cosine `η`, azimuthal cosine `ξ`,
`sinθ = √(1−μ_z²)`, `M` ordinates stored **η-ascending** (⟹ ω-DESCENDING from
near π to near 0), `ω_m = atan2(ξ_m, η_m)`:

**(1) CHORD (pre-carve production; the pre-6.4 shipped convention)**
```
e[0] = -sinθ ;  e[M] = +sinθ ;  e[m] = ½(η_{m-1} + η_m)   for m = 1..M-1
```
Endpoint pathology: a node ON the arc start (`ω=π`, i.e. on Σ) gives
`τ_0 = 0` bit-exactly ⟹ division by zero in the recurrence. **This is what the
`[½,1]` absorber was blocking.**

**(2) GEOMETRIC ARC ω-MIDPOINT (LANDED — the one under suspicion)**
```
ω_edge[0] = π ;  ω_edge[M] = 0 ;  ω_edge[m] = ½(ω_{m-1} + ω_m)
e[m] = sinθ · cos(ω_edge[m])
```
On an equispaced-ω rule this equals `ω_m ± Δω/2`, and P2 gives the closed form
`τ_m = ½ + ½ cot(ω_m) tan(Δω/4)` — `[M]` reproducing production to
`1.1e-16 / 2.2e-16 / 7.8e-16 / 7.4e-15 / 2.3e-14` at `n_φ = 4/8/16/32/64`
(⚠ the agreement DEGRADES two orders; the gate uses `atol=1e-13`).

**(3) RECURSION-DEFINED —** ⛔ **ILL-POSED. DO NOT ATTEMPT. `[M]` §9bis.1:**
exactly one edge per level (the mid-level one, where `κξ > sinθ`) has no real
solution, at every quadrature order. α is the arc-edge dome scaled by κ — an
amplitude, not a displaced position. The text below is preserved because it is
what a reader would otherwise re-derive; its premise is the refuted P-B.
⟹ **the live candidate is (C): P2 in the MARCH variable — see §9bis.7.**

*(superseded text follows)* **NOT YET MEASURED. This is the next step.**
The edge α's own recursion implies. Two equivalent routes, and they should be
checked against each other:
* invert α: from `α_{m+1/2} = α_{m−1/2} − w_m η_m` with `α_{1/2} = 0`, the edge
  is whatever value makes α equal its tangential target — per T3 this is
  `κ`-scaled relative to (2), with `κ = Δω/(2 sin(Δω/2))`;
* or run BMC Eq. 43's ν recursion with the τ that zeroes BMC's β, and read off
  the ν ladder as the partition.

⚠ Expect (3) to sit BETWEEN (1) and (2), since `κ = 1 + Δω²/24 + …` and
`cos(Δω/2) = 1 − Δω²/8 + …` bracket the geometric edge from opposite sides.
**If (1) is closer to (3) than (2) is, that explains §4.6 completely** and the
correct outcome is a partition change *plus* a principled endpoint guard, not a
revert to the clamp.

Also worth measuring for completeness: **(4) cumulative-weight** — already
refuted (§4.4), keep it in the table as the negative control.

---

## 7. WHAT IS LANDED — the tree's actual state

Branch `refactor/operator-strategy-layers`. Verify each with
`git merge-base --is-ancestor <hash> HEAD`. **None is on `main`.**

| commit | content |
|---|---|
| `3dda18ca` | the carve: `angular_cell_edges_per_level` minted as ONE partition producer; the two τ producers collapsed into one geometry-free body (P2 carries no geometry); `morel_montry_tau_raw_per_level` RETIRED; the `[½,1]` absorber RETIRED; an explicit refusal for non-arc levels; 17 gate re-poses |
| `d5067c4d` | `orpheus/derivations/discrete/sn/angular_differencing.py` — P0–P4 ladder, ν-closure, both β's, the τ/β nomenclature; **retires `contamination.py`**; fixes its sphere-arm API bit-rot. ⚠ `morel_montry_weights` now DELEGATES to production, so it is no longer an independent reference — the τ legs of `_c_surrogate`-based gates are TAUTOLOGICAL and say so |
| `39b46a31` | the two moved cylinder snapshots re-baselined, with the falsified "all three will move" prediction corrected in place |
| `c9bb61b4` | corpus pass: `morel-montry-clamp` → `morel-montry-closure`, a NEW `angular-cell-partition` equation carrying the per-geometry `cases`, a new `sn-tau-absorber-retirement` section. `sphinx -E -W` exit 0, 0 warnings; `DEAD TARGETS: 0` |
| `c97128d7` | the resolvent verdict + the refuted premise, recorded in the triage |

**Also changed and NOT to be silently reverted with the carve:**
* the SPHERE is untouched (cumulative-weight edges, literature-confirmed
  verbatim). `[M]` 4/8 of its τ are below ½ at S₈, 8/16 at S₁₆ — which is the
  literature's own point that `[½,1]` was never the admissible range in EITHER
  arm. Any revert must keep this.
* ⭐ **cylinder-P3 became a THEOREM** under (2): on a monotone arc the
  ω-midpoint edges bracket their own node, so `τ ∈ (0,1)` is forced (`[M]` 4000
  random arcs: `min(τ) = 4.7e-07`, `min(1−τ) = 7.6e-10`), and its only equality
  case is a node ON Σ ⟹ cylinder-P3 reduces to the fold criterion `Σ = ∅`. P3
  keeps its teeth on the SPHERE. **If the partition changes again, re-check
  whether this theorem survives.**
* several gates were re-posed with their superseded claims kept in place (τ
  bound `1/5 → 1/4`; the reversal identity from bit-exact to 64 ULP; the τ
  trichotomy narrowed; two retire-with-tombstones). A revert must re-pose them
  BACK, not delete them.

### Current red state

**ONE red attributable to this work:** the §4.6 resolvent row. Everything else
green: `tests/sn/sweep` + MMS + admission = 800 passed / 8 pre-existing #51
ledger reds; the slow curvilinear tier 310 passed / 27 xfailed / 0 failed;
`tests/sn/regression` 17 passed; the other 7 rows of `test_phase_c_crosscheck`
pass.

⚠ **Separately found, PRE-EXISTING, not caused by this work:**
`test_xi_odd_companion_sees_the_tie_break_but_does_not_adjudicate_it` had been
RED **in the slow tier since the 6.3 flip** (its full-circle fixture is
inadmissible), uncounted by the `-m "not slow"` ledger. Retired-with-tombstone
in `3dda18ca`. **The "#51 = 9 reds" ledger was measured with slow DESELECTED —
do not treat it as the whole tree.**

---

## 8. THE DECISION, AND WHAT TO DO NEXT

⛔⛔ **§8 IS ANSWERED AND ITS INSTRUCTION IS VOID — 2026-08-11, see §9bis.**
The answer to its own opening question is the SECOND branch: *the arc partition
is right and something else must move with it* — the CHART the weight is read
in. Candidate (3), which the rest of this section instructs you to build, is
ill-posed (§9bis.1). **Go to §9bis.7 for the live candidate.** Preserved below
because §9bis's whole argument is a correction OF this text, and because the
probe shape in the second paragraph is still the right shape.

*(superseded text follows)*

**The open question:** is the geometric arc partition wrong, or is it right and
something else must move with it?

**Do this before deciding:** implement candidate (3) from §6 and measure it on
BOTH instruments that discriminate — §4.6 (the resolvent shape, the decisive
one) and §4.5 (the `nx=320` angular-isolated ladder, the supporting one) —
alongside (1) and (2) on the same fixtures. Four rows, two instruments. That is
the whole experiment, and it is cheap: the MMS solves are ~0.3 s each and the
resolvent reference need be computed ONCE and reused (see `q64_phaseE_old_vs_new.py`'s
shape — it calls `T._run_cyl_2g_3reg_full()` once, then diffs each snapshot's
`scalar_flux` against it via `T._interpolate_per_region`).

Then one of:

* **(a) (3) wins** ⟹ re-pose the partition onto the recursion-defined edges,
  keep everything else (the single producer, the collapsed τ body, the
  nomenclature, the analysis module, the corpus work). Add a principled
  endpoint guard for the `τ_0 = 0` pathology instead of the retired clamp. This
  is the outcome the structure and the numbers would BOTH want.
* **(b) (1) wins and (3) is no better** ⟹ the chord partition was
  approximately right. Revert the partition change, KEEP the absorber
  retirement only if the endpoint pathology is otherwise guarded, and keep
  every non-partition deliverable. Record that P-A's defect is real but its fix
  is unknown, and re-open it as its own issue.
* **(c) nothing beats (1)+absorber** ⟹ revert the carve, keep the analysis
  module and the nomenclature and the corpus corrections, and file the α/τ
  partition mismatch as an open architectural defect with this memo as its
  evidence base.

⛔ **DO NOT:** relax `tol_per_cell`; close Q5.6.4; re-baseline the resolvent
comparison; or cite `τ ≡ ½` as principled. Any of those erases the signal.

---

## 9bis. ⭐⭐ ATTEMPT 2 (2026-08-11, after the compaction) — §8's experiment
## could not be run, and the reason rewrites §2 and §3

Raw probe output: `scratch/q64_attempt2_probe_outputs.md` (committed).
Probes: `$CLAUDE_JOB_DIR/tmp/q64_probe{A,B,C,D,E}_*.py` (**ephemeral**; each
one's purpose and shape is described below well enough to rebuild).

### 9bis.1 ⛔ CANDIDATE (3) IS ILL-POSED — there is no "recursion-defined edge"

§6's candidate (3) — "the edge α's own recursion implies" — presumes α has an
edge POSITION. It does not. `[M]` probe B1, `folded_product(4, n_φ)`, level 0:

| `n_φ` | rel spread of `α_prod / ξ(e_arc)` over interior edges | fitted constant |
|---|---|---|
| 8 | `1.555e-16` | `0.713918 = κ · 0.695709690275` |
| 16 | `1.110e-15` | `0.700200 = κ · 0.695709690275` |
| 32 | `1.912e-15` | `0.696829 = κ · 0.695709690275` |
| 64 | `2.313e-14` | `0.695989 = κ · 0.695709690275` |

α's **shape is exactly ξ at the geometric arc edges** and its **amplitude is
κ times the exact dome** — the fitted `w_gl` is the same number to 12 digits at
every order, which is the level's polar weight. An amplitude is not a position.
Asking "at which edge does the exact dome take the value α has?" means solving
`ξ(ẽ) = κ · ξ(e_arc)`, and `[M]` probe A2: **exactly one edge per level has no
real solution — the mid-level one, where ξ is maximal and `κξ > sinθ`** (min RHS
`−1.371e-02 / −3.347e-03 / −8.319e-04 / −2.077e-04` at `n_φ = 8/16/32/64`). It
fails at every quadrature order, and refining does not rescue it.

⟹ **§8's stated next step cannot be taken.** Not "was measured and lost" —
*is not a well-formed question*.

### 9bis.2 ⛔ §3 IS ITSELF REFUTED: P-B is TRUE in shape

§3 concluded from T3's `κ ≠ 1` that "α is **not** at the geometric arc
boundary". B1 measures the opposite where it matters: α is proportional to
`ξ(e_arc)` to `1e-16`, i.e. **α references exactly the arc partition**, and κ is
an amplitude defect the *quadrature* carries (midpoint-rule vs exact integral of
`cos ω` over the cell), not a displacement of the edge. The literature sentence
at `scratch/q64_tau_edge_convention_literature.md:943` distinguishes a
recursion-defined edge from the geometric one *in general*; on OUR rule the two
coincide up to that scalar, so the distinction has no numerical content here.

⟹ **the carve's PARTITION choice (2) was correct**, and §2's P-C is reinstated
on a different ground than P-B: not "α lives at the arc edge" but "α's shape is
the arc-edge dome". P-A, however, is refuted as a *duplication* claim: α never
derives a boundary. It derives edge-INDEXED amplitudes from `(w, η)` by the
recursion `α_{m+1/2} = α_{m−1/2} − w_m η_m`. There was never a second derivation
of the same object, so there was no Cardinal-Rule-2 twin to collapse — though
naming the partition once was still right.

### 9bis.3 ⭐⭐ THE REAL DEFECT: the partition is read in ω, the WEIGHT in η

The landed convention takes the cell boundaries in ω (correct, §9bis.2) and then
takes P2 — the barycentric coordinate — **in η**. Those are two different charts
on the same level, and the mismatch is the whole regression.

`[M]` probe A3, ν-closure run in BOTH variables (`0` / `1.000000` = exact close):

| `n_φ` | convention | in η | in ω |
|---|---|---|---|
| 8 | chord | `1.000000` | `0.159076` |
| 8 | chord + `[½,1]` absorber | `1.016389` | `0.079538` |
| 8 | arc / ω-mid (LANDED) | `1.000000` | `0.148847` |
| 8 | `τ ≡ ½` | `1.164784` | **`0.000000`** |
| 64 | chord | `1.000000` | `0.020411` |
| 64 | arc / ω-mid (LANDED) | `1.000000` | `0.018719` |
| 64 | `τ ≡ ½` | `1.002412` | **`0.000000`** |

⟹ **§5.3's reading of ν-closure was chart-relative.** It says "τ came from a
contiguous partition **in the variable you ran it in**" — nothing more. §2's P-F
condemned `τ ≡ ½` for failing it by 16.5 %; run in the march variable, `τ ≡ ½`
closes **exactly** at every order and the other three fail. The instrument
cannot adjudicate; it only reports which chart you chose.

### 9bis.4 ⭐⭐ THREE OF THE FIVE INSTRUMENTS ARE τ-BLIND — measured with garbage

§5's table lists blind spots per instrument. The sharper fact, `[M]` probe D on
`folded_product(4,16)` level 0 with the quadrature and partition held FIXED and
τ replaced by garbage:

| instrument | production (arc) | `τ ≡ ½` | GARBAGE `τ ≡ 0.7` | GARBAGE random |
|---|---|---|---|---|
| P1 `c = Σwη²` | `0.282432589924` | `0.282432589924` | `0.282432589924` | `0.282432589924` |
| **BMC β** (Eq. 41/75) | `0` | `0` | `0` | `0` |
| **Lathrop β** (α defect, Eq. 25) | `0.741555747146` | `0.741555747146` | `0.741555747146` | `0.741555747146` |
| ν-closure (in η) | `1` | `1.03918231642` | `1.01319215547` | `0.972815198734` |
| ⭐ amplification (new) | `4.72887003107` | `1` | `0.428571428571` | `2.18309886184` |
| ⭐ P1 closure defect (new) | `0.139106457701` | `0.00777890941286` | `0.0211225798756` | `0.274253487045` |

⟹ **P1, BMC β and Lathrop β are functions of `(quadrature, partition)` ONLY.**
They take τ nowhere. No fixture, tolerance, or refinement can make them grade a
τ — they certify a random τ as readily as the production one. That is why 6.4
could reach a worse answer with every structural instrument green: **the whole
published suite was τ-blind, and the one instrument that sees τ (ν-closure)
was read in the chart the chosen convention had just made exact.** Mode 12,
applied to the adjudication apparatus rather than to a gate.

### 9bis.5 ⭐ WHAT THE `[½,1]` ABSORBER ACTUALLY WAS: a non-amplification guard

Not "compensating the chord partition's stretched end cells" (§2's P-D). The
angular recurrence is `ψ̂_{m+1/2} = (ψ_m − (1−τ_m)ψ̂_{m−1/2})/τ_m`, so an
upstream error reaches the next face multiplied by `−(1−τ_m)/τ_m`:

> **`τ_m ≥ ½ ⟺ |(1−τ_m)/τ_m| ≤ 1 ⟺ the recurrence does not amplify.`**

`[M]` probe B2, worst running product `Π|(1−τ)/τ|` over the level:

| `n_φ` | chord | chord + absorber | arc (LANDED) | `τ ≡ ½` | #τ<½ (chord / arc) |
|---|---|---|---|---|---|
| 8 | `5.03` | **`1.000000`** | `3.36` | **`1.000000`** | 2 / 2 of 4 |
| 16 | `1.02e+01` | **`1.000000`** | `4.73` | **`1.000000`** | 4 / 4 of 8 |
| 32 | `2.04e+01` | **`1.000000`** | `6.68` | **`1.000000`** | 8 / 8 of 16 |
| 64 | `4.07e+01` | **`1.000000`** | `9.44` | **`1.000000`** | 16 / 16 of 32 |

min τ → `1/5` (chord) and `1/4` (arc) as `n_φ → ∞`; the amplification grows
**∝ M** for both. This explains, with one mechanism, every number in §4:
retiring the absorber alone made the MMS floor 1.8–3.4× worse (§4.1 — the
amplification was unleashed); the arc is worse than chord+absorber (§4.5/§4.6 —
`9.44` vs `1.0`); and `τ ≡ ½` sits exactly at the non-amplifying boundary.

⚠ Honest caveat: `τ ≡ ½` is *marginally* non-amplifying (product exactly 1 —
errors neither grow nor decay), and `ψ̂_{m+1/2} = 2ψ_m − ψ̂_{m−1/2}` is linear
extrapolation, so it has the diamond scheme's usual positivity exposure. `τ > ½`
damps. That is a real trade to state, not a refutation.

### 9bis.6 ⭐ THE HONEST τ INSTRUMENT: the P1 closure defect at the arc faces

Near the diffusion limit `ψ` is affine in the direction cosines
(`ψ = a + bη + cξ`), so feed `η` and `ξ` through the closure and compare the
produced face values against their true values at the arc edges — weighted by
`ξ(e_arc)`, the coefficient α actually multiplies in the balance. `[M]` probe B3:

| `n_φ` | convention | max\|Δη\| | max\|Δξ\| | ξ-weighted |
|---|---|---|---|---|
| 64 | chord | `6.09e-04` | `6.64e-01` | `3.37e-01` |
| 64 | chord + absorber | `1.22e-03` | `1.63e-02` | `5.79e-04` |
| 64 | arc (LANDED) | **`6.94e-16`** | `1.42e-01` | `7.16e-02` |
| 64 | `τ ≡ ½` | `1.23e-03` | `6.13e-04` | **`4.79e-04`** |

⭐ Read the LANDED row: it is exact on η **by construction** — it is P2-in-η on
its own partition — and pays with a ξ defect that barely converges
(`3.89e-1 → 1.42e-1` over an 8× refinement, ≈ O(M^−½)). **The convention made
one basis function exact and was then certified by the instrument measuring that
same basis function.** `τ ≡ ½` is the balanced choice, both defects O(M⁻²).

No single τ can be exact for both: requiring it gives
`sin p + sin q = sin(p+q)` for the two half-cell angles, whose only solutions
are `p = 0` or `q = 0` (the node ON an edge). So a τ closure is intrinsically
first-order-in-one-chart, and the choice of chart IS the design decision.

### 9bis.7 THE RE-POSE — one rule, both arms, no clamp

Not "hardcode ½ on the cylinder". The rule is:

> **τ is the barycentric coordinate of the ordinate between its own cell's two
> edges, measured in the variable the level's angular march RUNS IN.**
> Sphere: the march is in μ. Cylinder: the march is in ω, arc by arc (T22b).

`[M]` probe E:
* **E1 — the sphere does not move: `0.000e+00` at N = 4/8/16/32/64.** BMC Eq. 12
  stays verbatim, because there the march variable *is* the radial cosine.
* **E2 — the cylinder gives `max|P2_in_ω − ½| = 0.0 … 9.0e-15`** on every
  shipped `folded_product(n_mu ∈ {2,4}, n_φ ∈ {4..64})`. So ½ is a
  **consequence** on equispaced arcs, never an assertion — and on a deliberately
  non-equispaced monotone arc the body returns a non-trivial
  `[0.408, 0.700, 0.250, 0.581, 0.481]`, all inside `[0,1]`.
* **E3 — R12a is untouched**: the march-start facts are integer facts read off
  `ξ = 0` and `η_0 = η_1` bit-exactly (Q5.4/T26) and never consulted τ.
  `_assert_tau_within_unit_interval` passes trivially, and the chord partition's
  `τ_0 = 0` endpoint pathology (§6, the thing the absorber was blocking) becomes
  **unspellable**: τ is ½, never 0. No clamp, no endpoint guard.

⟹ this **keeps** the carve (one named partition producer, the collapsed τ body,
the nomenclature, the analysis module, the corpus work) and changes the CHART
the partition and the weight are read in. `angular_cell_edges_per_level` should
return the partition **in the level's march variable** (μ sphere, ω cylinder)
instead of projecting the cylinder's to `sinθ cos ω_edge`; consumers that need
the radial cosine (BMC's β functionals) convert there, where the cosine is the
functional's own variable. That also makes the carve's existing claim — *"P2
carries no geometry, so ONE body serves both arms"* — true, which today it is
not: one body reads both arms in the radial cosine, and only the sphere marches
in it.

### 9bis.8 THE DECISIVE ROW — the independent method

Probe C re-solves `cyl_2g_3reg_folded_4x8_dd_n40` under each convention and
compares the flux SHAPE against the trajectory-resolvent reference (computed
once, reused), exactly as
`test_phase_e_trajectory_resolvent_flux_shape_crosscheck` does but with a live
solve in place of the frozen snapshot. Gate tolerance `1.2e-1`.

**⭐ The instrument is VALIDATED**: it reproduces both of §4.6's anchors to
every printed digit *and* both k_eff values to all 16 — `6.5934e-02` /
`1.2302082296342958` for chord+absorber and `1.2676e-01` /
`1.2310212585879858` for the landed arc. Whatever it says about the other rows
can be trusted to the same standard.

`[M]` 2026-08-11, `folded_product(4, 8)`, `n = 40`, per-group max
`|Δφ_norm|` after L∞ normalisation, tol `1.2e-1`:

| convention | τ (level 0) | per group | max | verdict |
|---|---|---|---|---|
| **(1c) chord + `[½,1]` absorber** | `0.500, 0.500, 0.586, 0.780` | `[0.06593, 0.01335]` | **`6.593e-02`** | PASS ✅ **best** |
| **(C) `τ ≡ ½`** | `0.500, 0.500, 0.500, 0.500` | `[0.10181, 0.01507]` | **`1.018e-01`** | PASS |
| (2) arc, P2 in η (LANDED) | `0.260, 0.459, 0.541, 0.740` | `[0.12676, 0.01341]` | `1.268e-01` | **RED** |
| (1) chord, unclamped | `0.220, 0.414, 0.586, 0.780` | `[0.14409, 0.01305]` | `1.441e-01` | **RED** |

⚠⚠ **THE RANKING DISAGREES WITH THE OTHER TWO INSTRUMENTS, and this is the
open question of attempt 2.** `τ ≡ ½` passes but is 1.54× worse here, while
§4.5's MMS ladder puts it 4.4× *better* and §9bis.6's closure defect 1.2×
better. All three agree on the bottom two, and in the same order. So the live
disagreement is exactly `τ ≡ ½` vs `chord + absorber`, on a **single coarse
rung** (`n_φ = 8`, where `κ − 1 = 2.6 %` is largest).

⛔⛔ **AND THE ROWS ABOVE ARE STARVED SOLVES.** `[M]` every one of the four
emitted
`ConvergenceWarning: inner(source-iteration) hit max_inner=300 without reaching
tol=1e-10 (last residual ≈1.9e-06) … needs about 528`. The residuals are
near-identical across rows (`1.913 / 1.886 / 1.906 / 1.893 e-06`) and are five
orders below the `1e-1` shape differences being graded, so the ranking is
probably robust — but it is **not** established until re-run converged.

⟹ ⭐ **A finding in its own right, independent of τ: the shipped fixture
`cyl_2g_3reg_folded_4x8_dd_n40` is a starved solve at the snapshot defaults**
(`max_inner = 300`, needs ≈528). That is a live #340-class defect on a
regression artifact AND on the L1 gate that consumes it — file it separately;
it is not caused by, and does not depend on, anything in Q5.6.4.

### 9bis.8b ⛔⛔ THE RESOLVENT CROSS-CHECK IS REFERENCE-LIMITED — §4.6 IS NOT DECISIVE

`[M]` probe C2, re-run at `max_inner = 4000` with **every row reporting
`fully_converged`**, as a ladder in `n_φ`. One cached reference grades every
rung (the trajectory-resolvent method is independent of the SN quadrature):

| convention | `n_φ`=8 | 16 | 32 |
|---|---|---|---|
| (1) chord, unclamped | `1.4409e-01` | `4.6313e-02` | `2.9175e-02` |
| (1c) chord + `[½,1]` absorber | `6.5934e-02` | `3.7739e-02` | `3.0323e-02` |
| (2) arc, P2 in η **(LANDED)** | `1.2676e-01` | **`2.2008e-02`** | `2.7606e-02` |
| (C) arc, P2 in ω (march variable) | `1.0181e-01` | `3.1985e-02` | `2.7626e-02` |

⚠ The **converged** `n_φ=8` column is bit-identical to §9bis.8's starved one to
every printed digit — so the starvation was NOT driving that ranking. And (2) is
**non-monotone** (`2.20e-2 → 2.76e-2`): it undershoots at 16 by error
cancellation against the reference, then settles.

⭐ **All four collapse to ≈2.8–3.0e-2 by `n_φ = 32`, within 10 % of each other.**
`[M]` probe F decides WHOSE error that is — hold `n_φ = 32`, sweep `nx`:

| convention | `nx`=40 | 80 | 160 |
|---|---|---|---|
| (1c) chord + absorber | `3.0323e-02` | `3.2314e-02` | `3.3682e-02` |
| (2) arc, P2 in η | `2.7606e-02` | `2.9204e-02` | `2.9234e-02` |
| (C) arc, P2 in ω | `2.7626e-02` | `2.9655e-02` | `3.1045e-02` |

⟹ **the floor does not fall — it RISES, monotonically, in every column.**
Refining the SN mesh 4× makes the agreement *worse*. That is the signature of a
comparison in which **the SN side is converging past the reference**: the floor
is the **trajectory-resolvent reference's OWN discretisation error**
(`n_r=24, n_mu_axial=16, n_phi_az=32, n_traj_quad=64`), not the SN's spatial
error, and **no refinement on the SN side can push past it.**

⛔ **Consequences, and the first retires the memo's own §4.6 framing:**

1. **§4.6 is NOT "THE DECISIVE MEASUREMENT".** It is one coarse rung
   (`n_φ = 8`) read against a reference whose own error is `≈3e-2` — **20–40 %
   of the difference it is used to grade.** Its dynamic range for an angular
   claim is exhausted by `n_φ = 16`.
2. It does carry real signal AT `n_φ = 8` (SN error `6.6e-2 … 1.44e-1` exceeds
   the floor), and there it ranks `(1c) ≺ (C) ≺ (2) ≺ (1)`.
3. ⚠ **A finding about the L1 gate itself:**
   `test_phase_e_trajectory_resolvent_flux_shape_crosscheck[cyl_2g_3reg_folded_4x8_dd_n40]`
   carries `tol_per_cell = 1.2e-1` against a reference floor of `≈3e-2`. The
   gate **cannot be tightened by improving the SN** — only by refining the
   *reference*. Its docstring credits it with pinning "the resulting shape
   agreement"; that claim is bounded by the reference and the bound is nowhere
   stated. It is also pinned to the single coarsest cylinder rule, which is the
   only rule at which it can still see a closure at all.
4. ⟹ **§5's instrument table needs a demotion**: the resolvent's blind spot is
   not merely "own discretisation error; loose bar; shape-only" — it is a **hard
   floor that SN-side refinement moves the WRONG WAY**.

### 9bis.8c THE VERDICT ACROSS ALL INSTRUMENTS

| instrument | can it grade a τ? | ranking |
|---|---|---|
| P1 `c = Σwη²`, BMC β, Lathrop β | ⛔ **τ-blind** — garbage-identical (§9bis.4) | — |
| ν-closure | ⛔ chart-relative label only (§9bis.3) | — |
| resolvent shape | ⚠ only at `n_φ=8`, its own error 20–40 % of signal | `(1c) ≺ (C) ≺ (2) ≺ (1)` |
| MMS aniso-cyl, `nx=320` (§4.5) | ✅ range to `9e-6`, 3 orders below the resolvent floor | `(C) ≺ (1c) ≺ (2) ≺ (1)`; (C) best order (3.29/3.15/2.10) |
| P1 closure defect (§9bis.6) | ✅ structural, no fixture | `(C) ≺ (1c) ≺ (2) ≺ (1)` |
| amplification (§9bis.5) | ✅ the `τ ≥ ½` half | `(1c) = (C) ≺ (2) ≺ (1)` |

⭐ **Every discriminating instrument ranks the LANDED convention (2) below both
(1c) and (C), so the landed state must change regardless.** Three of the four
rank **(C) first**; the fourth is reference-limited and speaks only at the
coarsest rule. With the literature (§9bis.9 — Hébert's own rule *and* his own
closure ARE (C), both arms verbatim, Lathrop's objection measured absent), the
decision is **(C)**.

⚠ The one dissenting number, stated plainly: at `n_φ = 8` chord+absorber really
is 1.5× closer to the resolvent reference than (C). On the MMS at that same
`n_φ = 8` the two sit within 2.5 % (`3.4258e-3` vs `3.5111e-3`), and (C) pulls
4.4× ahead by `n_φ = 64`. Two different fixtures — and the coarsest rung is
where an asymptotic argument has least force.

### 9bis.10 ⛔⛔⛔ THE VERDICT — (C) IS REFUTED AND THE LANDED CARVE IS CORRECT

**Read this before §9bis.3 through §9bis.9.** An adversarial review
(`scratch/q64_attempt2_qa_review.md`) broke the chain at C6, and re-reading BMC
directly confirms it. Every number in §9bis stands; two of my inferences do not.

#### The three sentences that decide it, `[M]` read myself in the BMC sidecar

`scratch/literature_ocr/Bailey-Morel-Chang(2010)…md`:

* **line 366, Eq. (41)**: `Σ_m μ_m[α_{m+1/2}μ_{m+1/2} − α_{m−1/2}μ_{m−1/2}] = 0`
  — the first-order condition. **Genuinely τ-free**, so §9bis.4's *measurement*
  is right.
* **line 372/376, Eqs. (42)/(43)**:
  `τ_m = (μ_m − μ_{m−1/2})/(μ_{m+1/2} − μ_{m−1/2})`, i.e.
  `μ_m = τ_m μ_{m+1/2} + (1−τ_m)μ_{m−1/2}`.
* **line 657**: *"…the same factor that Morel and Montry forced to be zero to
  preserve the Galerkin diffusion approximation. **Forcing this β factor to be
  zero determines the Morel and Montry weights.**"*

⟹ ⛔ **§9bis.4's INFERENCE is wrong.** β is τ-free *by construction*: β = 0 is a
condition on `(α, edges)`, and **τ is then DEFINED by Eq. 42 as the barycentric
coordinate of those edges.** τ is *downstream* of β, not an input to it. "τ-free
functional" therefore does NOT mean "cannot grade a τ" — it means **τ is
determined, not chosen.**

⟹ ⛔⛔ **§9bis.3's whole framing collapses.** There is no "chart freedom".
**Eq. 43 interpolates μ — the radial cosine — full stop.** So P2-in-η IS BMC
verbatim, and the LANDED convention is the literature's own closure on the
correct partition. My "the partition is read in ω and the weight in η, two charts
on one level" was not a defect; it was the published requirement.

#### And BMC refute (C) explicitly, in their own CYLINDER section

* `[M]` **line ~382, Eq. (53)** is the cylinder closure, and the next sentence
  reads: *"As in the spherical-geometry case, `τ = 1` gives the step scheme and
  **`τ = ½` gives the diamond scheme**."* So (C) is, in BMC's own cylinder
  section, named the diamond.
* `[M]` **line 172**: *"…only the **weighted** diamond difference discretization
  developed by Morel and Montry produces a good diffusion approximation for the
  **first-order** scalar flux… preserving first-order consistency eliminates the
  flux dip **in general**, while preserving leading-order consistency only
  eliminates the flux dip in highly diffusive problems."*
* `[M]` **lines 408 / 637**: *"as the number of quadrature points is increased
  with these methods [step and diamond], their respective solutions approach the
  correct diffusive behavior through first order."*

⟹ **§9bis.9's literature defence is REFUTED.** Hébert 2009 presents the plain
diamond; **BMC 2010 exists precisely to show the diamond is insufficient in
curvilinear geometry.** Citing the textbook against the paper that improved on it
is backwards. `beta_first_order_consistent = True` on the M-M closure is exactly
BMC's first-order property, and (C) would forfeit it.

⭐ **Lines 408/637 also explain the MMS ladder that had favoured (C).** The
diamond's contamination → 0 as M grows, so at fine `n_φ` the distinction BMC care
about *vanishes* — which is exactly the regime §4.5 measures (`nx=320`,
`n_φ` up to 64). The MMS L2 ladder is Mode-12 blind to first-order
consistency, and BMC state the mechanism themselves.

#### My C6 instrument was wrong in its BASIS, and corrected it inverts

A 1-D cylindrical flux is **even in ξ** — that is exactly what the σ_y fold
quotients — so its ω-expansion is **cosine-only**: `ψ ~ Σ a_m cos mω`, with
`η = sinθ cos ω` the **m = 1** harmonic and `ξ = sinθ sin ω` **not a realisable
mode at all**. `[M]` probe G0: a folded rule stores only `ξ ≥ 0`, so
`Σwξ = +6.703` folded vs `0.000` unfolded — the signed ξ is not represented.
§9bis.6 fed ξ through the closure and *weighted by it*: it graded a mode the
solution cannot have. (The review also tested five weightings — the ranking was
stable, so the weighting was not the defect.)

`[M]` probe G1, re-run on `{cos mω}`, max over m = 0…4:

| `n_φ` | (1) chord | (1c) chord+abs | **(2) arc, P2 in η (LANDED)** | (C) arc, P2 in ω |
|---|---|---|---|---|
| 16 | `5.161e-01` | `8.284e-01` | **`4.430e-01`** | `8.284e-01` |
| 32 | `1.064e-01` | `1.648e-01` | **`8.584e-02`** | `1.648e-01` |
| 64 | `2.517e-02` | `3.918e-02` | **`1.994e-02`** | `3.918e-02` |

and on the **m = 1 column, which IS BMC Eq. (43)**, (2) is EXACT —
`4.2e-16 / 2.8e-16 / 1.1e-15 / 2.1e-15` at `n_φ = 8/16/32/64` — while chord
reads `5.4e-2 → 1.2e-3`, the absorber `1.4e-1 → 2.4e-3` and (C)
`1.6e-1 → 2.4e-3`. ⟹ ⛔ **§9bis.6's "the convention made one basis function
exact and was then certified by the instrument measuring it" is WRONG.** That
basis function is `μ`, and reproducing `μ` at the faces is the *published
requirement*. (2) meets it to machine precision; nothing else comes within 12
orders. (⚠ at `n_φ = 8`, M = 4, the m ≥ 3 columns exceed the rule's resolving
power and are meaningless — read only m ≤ 2 there.)

#### My other two overreaches, corrected

* ⛔ **§9bis.1's "at every quadrature order" is FALSE — an even-M parity
  artefact.** `[M]` probe G2: `0 of M+1` failing edges at **odd** M
  (`n_φ = 6/10/14/18`, all constructible), `1 of M+1` at even M. At even M an
  edge sits exactly at `ω = π/2` where ξ is maximal; at odd M the edges straddle
  it. The honest claim: **on every rule with even M — which is every shipped
  fixture (`n_φ = 4/8/16/32/64 ⟹ M = 2/4/8/16/32`) — one edge has no real
  recursion-defined position.** Candidate (3) stays disqualified (a definition
  that exists for some legal rules and not others is not a definition), but my
  quantifier was wrong, in a memo about exactly that. → `plan-authoring` §2.
* ⛔ **§9bis.5 is mis-scoped in BOTH directions.** The END-TO-END product
  `Π|(1−τ)/τ|` is exactly `1.000000` for chord, arc **and** `τ ≡ ½` — because
  `τ(π−ω) = 1−τ(ω)` telescopes. What my table reports is the worst **transient
  partial** product (still real: every intermediate face is consumed by its
  neighbours' redistribution). The only genuinely **dissipative** convention is
  chord+absorber (`2.45e-02`, 40×). And positivity is **inverted** from what I
  implied: `τ ≡ ½` is the *safest* derived candidate (`min ψ̂ = −24.2` vs arc
  `−77.2` vs chord `−230`). ⚠ **No gate covers ψ̂ positivity on either arm** —
  the two curvilinear positivity gates are both on the sphere's converged
  *scalar* flux.

  ⛔⛔ **THE `−77` NUMBERS ARE MIS-SCOPED — corrected 2026-08-11, and this is
  MY error, propagated into two sub-agent briefs before it was caught.** They
  come from a probe driving the recurrence with an **inconsistent seed** (a
  positive constant seed against an `exp(−6cos ω)` field), not from the
  production value path. `[M]` I then re-measured and reproduced `−0.149` with
  120/400 faces negative on the shipped `cyl_2g_3reg_folded_4x8_dd_n40`
  fixture — **but that run passed `radial_characteristic=None`, which the
  method's own docstring says seeds at ZERO and is "legitimate only for the
  ψ-independent COEFFICIENT use".** So my re-measurement is the same class of
  artefact and settles nothing either.

  ⚠ **What IS established:** the exposure is `(1−τ)/τ`, a property of the
  angular-diamond FAMILY, and `[M]` **the retired absorber was never a
  positivity guarantee** — under the inconsistent seed it cuts the magnitude
  ~10× (`−23.3` vs arc `−77.2`) and still leaves 6/17 faces negative, and
  `τ ≡ ½` behaves the same (`−24.2`). That comparison is valid because all
  four conventions were driven identically.

  ✅ **SETTLED 2026-08-12, and it settles AGAINST the `−77`.** `[M]` on the
  production value path — converged flux plus the marched route-(a) ψ½ state —
  **ψ̂ is POSITIVE**: `+0.1337 / +0.1286 / +0.1287` at `n_φ = 6/8/16`, within
  12 % of `min ψ`. `−76.9` reproduces **only** with a random ψ AND a zero seed.
  ⟹ **the sign is a property of SEED CONSISTENCY, not of the closure.** Both
  regimes are now pinned separately, with the amplification `A(M)` (`2.41 …
  9.44` at `M = 2…32`) as the stated bound on the zero-seed one:
  `tests/sn/sweep/curvilinear/test_psi_half_positivity.py`.

  ⚠ **Standing rule from this episode: never cite a ψ̂ positivity number
  without stating which seed produced it.** Two sub-agent briefs and two
  user-facing reports carried the `−77` as a production property before it was
  caught.

#### ⟹ THE DECISION: KEEP THE LANDED CARVE

**Q5.6.4's partition change is correct and is the first time this codebase
satisfies BMC Eq. (43) exactly on the cylinder.** The chord+absorber it replaced
does not (ν-closure `1.016`), and (C) does not (`1.039`).

**The measured accuracy regression is NOT the closure being wrong.** The review
found the four candidates' accuracy ranking in *perfect rank correlation with the
recurrence's transient error amplification* and *inverted against closure
fidelity*. The retired `[½,1]` absorber bought damping by **forfeiting Eq. 43** —
it traded the published correctness property for stability, silently. So the
carve did not introduce a defect; **it removed a mask and exposed a pre-existing
amplification problem.**

▶ **The real remaining work is therefore the amplification / the march SEED, not
τ.** The review's unexamined third option — hold the arc closure and sweep
`MorelMontryAngularSweep.psi_half_seed` — is the one candidate that could recover
the accuracy *without* giving up Eq. 43. That is the next experiment.

⛔ **DO NOT** re-pose τ to ω; **DO NOT** restore the `[½,1]` absorber; **DO NOT**
cite Hébert against BMC. All three are now measured errors, not open options.

### 9bis.12 ⭐⭐⭐ MOREL & MONTRY 1984 — THE PRIMARY SOURCE. Q5.6.4 IS CLOSED.

J. E. Morel & G. R. Montry, *Analysis and Elimination of the Discrete-Ordinates
Flux Dip*, TTSP 13(5):615-633 (1984). Local since 2026-08-11; findings
`scratch/q64_morel_montry_findings.md`. Page map: **printed = PDF + 613**.

#### ⭐⭐⭐ THE LANDED CARVE IS MOREL & MONTRY'S OWN CYLINDER APPENDIX

Their appendix states the structural fact outright (printed 632): for the
cylinder the cell-edge **cosines** are NOT weight partial sums — the cell-edge
**AZIMUTHS** are. (A4) partitions in φ by cumulative weight from
`φ_{1/2} = −π`; (A3) maps `μ = sinθ_ℓ cos φ`; (A2) is barycentric **in the
radial cosine**. On an equal-weight level (A4) collapses to **equispaced-ω
midpoints**, and (A2)–(A4) give exactly
`τ_m = ½ ∓ ½ cot(φ_m) tan(Δφ/4)` (sign = march orientation).

`[M]` production against the four properties M&M's appendix states:

| `n_φ` | M | `τ₁` (→¼) | `τ_M` (→¾) | sinθ-independent? | `(τ₁−½)/w` | `−M/4` |
|---|---|---|---|---|---|---|
| 8 | 4 | `0.259892` | `0.740108` | spread `1.67e-16` | `−0.9604` | `−1.00` |
| 16 | 8 | `0.252425` | `0.747575` | `5.55e-16` | `−1.9806` | `−2.00` |
| 32 | 16 | `0.250603` | `0.749397` | `7.16e-15` | `−3.9903` | `−4.00` |
| 64 | 32 | `0.250151` | `0.749849` | `2.11e-14` | `−7.9952` | `−8.00` |
| 128 | 64 | `0.250038` | `0.749962` | `8.45e-14` | `−15.9976` | `−16.00` |

⟹ **all four reproduced.** ⭐ **The `τ → ¼` endpoint value that §9bis.11
diagnosed as the accuracy defect is Morel & Montry's own designed value, and
`(τ₁−½)/w = −M/4` is their formula's own O(M) divergence.** The carve did not
introduce it and no τ choice inside this method removes it.

#### ⭐⭐ THE BMC-43 vs R&L-16 CONFLICT IS THE PRIMARY'S ANNOUNCED TRADE

Not an open question. M&M cite **Reed & Lathrop as their reference 1**, reject
them at **printed 616** — R&L's *induced* quadratures cannot integrate degree
> 3, so conservation breaks under anisotropic scattering — and concede in the
next sentence, **printed 617** (verified on the scan), that unlike R&L theirs
works with **any** S_N quadrature set but is **"only first-order accurate."**

> **R&L buy second-order angle by MOVING THE ORDINATES. M&M buy
> diffusion-limit consistency on an ARBITRARY quadrature and pay in order.**

⟹ §9bis.11's "genuine conflict at the arc endpoints" is real but **already
adjudicated by the primary, deliberately**. Clamping τ toward ½ would trade
away β = 0 — the paper's entire result — for an accuracy order the paper never
claims. **Q5.6.4's design question is CLOSED.**

#### ⛔ A correction to §9bis.10's own reasoning (the decision is unchanged)

§9bis.10 argued "β = 0 *determines* τ, so τ is downstream, not chosen",
quoting BMC's line 657. `[M]` **that sentence of BMC's is false by dimension**:
β is ONE scalar over an (N−1)-dimensional solution family — three τ-vectors were
found with `β = O(1e-16)` at `‖τ − τ_MM‖_∞ = 0.238 / 0.308 / 0.242`.

What actually determines τ in the primary is **per-ordinate coincidence of the
scheme-implied edges (4a) with the weight-defined edges (15) — N equations**,
inverted to give (16b); **β = 0 then FOLLOWS by parity** (17a)–(19). It is a
consequence, not the defining condition, and the derivation is a one-term
ansatz `ψ_m = φ + 3Jμ_m` under the three standard quadrature requirements —
**not Galerkin**.

⟹ the conclusion (barycentric-in-μ is right) stands and is now better founded:
N per-ordinate equations, not one scalar condition. And §9bis.4's measurement is
fully explained — β cannot grade a τ because it is one scalar over a family,
which is a stronger statement than "τ-free".

#### The clamp question is now CLOSED at 5 of 5

`[M]` **NOT ADDRESSED** in M&M — τ is introduced bare at (16b)/(A2) and never
bounded, sign-conditioned, or special-cased; `τ ∈ (0,1)` follows automatically
from a monotone edge march plus node-in-cell. **Five of five primaries (M&M,
R&L, BMC, Hébert, Lathrop) prescribe no limiter.** The absorber retirement is
settled on the whole literature, and its origin is the Grant spatial-weight
transplant (§9bis.11).

#### Positivity: no authority either way — our gate stays a CHARACTERISATION

`[M]` **NOT ADDRESSED.** The only negativity M&M analyse is of the **diffusion
coefficient** via `β < 0`: `r_D = −2β/σ_t` is where D blows up and the slope
reverses. **That IS the dip — a scalar-flux local minimum, not a negative
angular flux.** So our `min ψ̂ ≈ −77` has no literature bearing in either
direction, and the ψ̂ gate must pin behaviour, not certify a contract.

#### ⭐ THE SEED LEVER, QUANTIFIED — this is the remaining work

M&M analyse the starting direction in depth (Miller–Alcouffe, Eqs. 10–14).
`[M]` **Eq. (14) admits an arbitrary edge-march seed `μ_s`, and for Gauss-S₂
diamond `β(μ_s) = μ_s + 2/√3` — LINEAR, UNIT GAIN, no damping, and
SIGN-DETERMINING:**

* underestimate ⟹ `β < 0` ⟹ **the flux dip**;
* overestimate ⟹ `β > 0` ⟹ no dip, **but the diffusion limit is lost**.

Their diagnostic is the *effective starting cosine* (10). Also printed 629: the
spatial truncation is `Δr²/r²` for the weighted fluxes vs `Δr²` for the
starting flux. ⟹ **the seed controls β directly and linearly; there is a
two-sided target, not a "make it more accurate" direction.** Both terminal
conditions are printed on one line (`α_{1/2} = 0` and `α_{N+1/2} = 0`) but
M&M say **nothing** about march direction or conditioning — that part remains
R&L's (§9bis.11).

#### ⚠ Four hazards to carry

1. **Notation collision across the three primaries.** M&M cylinder: `μ` radial,
   `φ` azimuth. BMC cylinder: `η` radial, `ω` azimuth. Ours: `η` radial, `ω`
   azimuth. Reading an M&M equation with BMC letters silently swaps the roles.
2. `[M]` **M&M normalise `Σ W = 1`, so their weight is HALF the cell μ-measure**
   (hence the `2W_m` in their Eq. 15 and the `2/r` in Eq. 1). Any factor-of-2
   comparison against our weights must account for it.
3. **"weight = cell μ-measure" is SPHERE-ONLY** — the primary says so directly
   (printed 632). §4.4 measured this from the other side.
4. ⛔ **New published typo (L-010 class): M&M Eq. (1)'s two ψ half-indices are
   CROSSED** (`α_{m+1/2}` printed against `ψ_{m−1/2}`). Their own (4)/(6a) are
   uncrossed, and only the uncrossed form telescopes under the zeroth moment —
   which is what makes (5) the exact continuity equation. Do not "fix" our code
   to match Eq. (1) as printed.

#### Lead, UNREAD

M&M's in-text "Miller–Alcouffe" credit resolves to their reference 2 =
**Dudziak, O'Dell & Alcouffe, LA-7911-PR (July 1979)** — no Miller in the author
list. Strong candidate for Hébert's unresolved ref. [36]. ⛔ The OSTI client
hit `ConnectionResetError` twice; this is a LEAD, not a source.

### 9bis.11 ⭐⭐⭐ REED & LATHROP 1970 — the absorber's origin, and the criterion we lacked

Source: W. H. Reed and K. D. Lathrop, *"Truncation Error Analysis of Finite
Difference Approximations to the Transport Equation,"* NSE **41**, 237 (1970).
Local since 2026-08-11; sidecar
`scratch/literature_ocr/07-Truncation Error Analysis…md`. Findings:
`scratch/q64_reed_lathrop_findings.md`. ⚠ the sidecar misreads a half-integer
subscript `1` as `3` on PDF p. 3 — the print is right.

#### ⭐⭐⭐ THE ABSORBER'S ORIGIN — a SPATIAL interval transplanted onto the ANGULAR weight

`[M]` **read on the rendered page** (PDF p. 4 = printed p. 239), footnote 8:

> *"Page 393 of Ref. 6 Grant lets the **a** weights (which he calls θ) depend on
> the sign of **a**, but this is necessary only to keep **a** between ½ and 1.
> Grant does not determine angular weights and limits his degrees of freedom by
> assuming an explicit form for α."*

On that same page, R&L's notation is unambiguous: **`a_{i+1/2}` is the SPATIAL
weighted-diamond parameter** (Eqs. 5a/5b expand `ψ_{i±1}` in `Δr` with it, and
Eq. 11 gives `a_{i+1/2} = ½ + ⅙(r_{i+1}−r_i)/(r_{i+1}+r_i)`), while **`τ_m` is
the ANGULAR weight** (Eqs. 6a/6b, in `w_m ∂ψ/∂μ`) and `α_{m±1/2}` is the
redistribution coefficient. The footnote's closing clause says outright that
Grant *does not determine angular weights at all*.

⟹ **ORPHEUS's retired `[½,1]` absorber was Grant's SPATIAL weighted-diamond
interval applied to the ANGULAR closure weight.** The retirement is therefore
not merely "uncited" — it corrected a transplant between two different schemes.
That is the origin story, sourced, and it closes P-D for good.

#### The R&L / Morel–Montry fork — the same closure, one extra equation

Both use the identical closure and the identical barycentric node — `[M]` R&L
**Eq. 13c** `μ_m = τ_m μ_{m+1/2} + (1−τ_m)μ_{m−1/2}` **is BMC Eq. 43, forty
years earlier** — both crediting Grant 1968. The fork is a THIRD equation R&L
impose and M–M do not: **Eq. 13b** `(1−τ_m)α_{m+1/2} + τ_m α_{m−1/2} =
(1−μ_m²)/2`. With it the triple becomes a quadratic in the ORDINATE, marched
from `α_{1/2}=0` — R&L take the μ mesh as given and solve for `(α, μ_m, τ)`,
i.e. **edges in, ordinates out**. That is Lathrop-2000's "fixes a quadrature",
and R&L admit the cost themselves (printed p. 248): 4th and higher even moments
are not integrated correctly, *"<3 % for the fourth moment and S₈ Gauss
weights"*. **We take the quadrature as given and solve for τ — the M–M branch.**

#### ⭐⭐ THE INSTRUMENT WE LACKED, and it diagnoses the arc

`[M]` **Eq. 15**: `τ_m = ½ − (1/2w_m)(μ_{m+1/2} + μ_{m−1/2} − 2μ_m)`, and
**Eq. 16**: the angular truncation error is second order **iff**
`μ_m = ½(μ_{m+1/2} + μ_{m−1/2}) + O(w_m²)` — the ordinate must be the **μ-midpoint
of its own cell**, equivalently `τ = ½ + O(w)`. Otherwise the angular term's
coefficient `[(1−τ)²α_{m+1/2} − τ²α_{m−1/2}]/w_m` is `O(1/w)` and the scheme
degrades to first order.

⭐ Unlike BMC's β this is **POINTWISE**, so the σ_y fold does not annihilate it.
It is solve-free. `[M]` on `folded_product(4, n_φ)`, level 0:

| `n_φ` | END cell μ-width ÷ Δω² | INTERIOR ÷ Δω | min τ | **max\|τ−½\|/w** |
|---|---|---|---|---|
| 8 | `0.2414` | `0.4577` | `0.259892` | `4.394e-01` |
| 16 | `0.2509` | `0.4954` | `0.252425` | `9.062e-01` |
| 32 | `0.2534` | `0.5051` | `0.250603` | `1.826e+00` |
| 64 | `0.2540` | `0.5076` | `0.250151` | `3.658e+00` |
| 128 | `0.2541` | `0.5082` | `0.250038` | `7.319e+00` |

⟹ **the endpoint cells are QUADRATICALLY narrow in μ (`∝ Δω²`) while carrying
LINEAR weight (`∝ Δω`)**, so the node cannot be their μ-midpoint to `O(w²)`:
`τ → ¼` there, refinement-INDEPENDENT, and `|τ−½|/w` **doubles every
refinement** — `O(M)` divergence.

⭐⭐ **This is the accuracy loss, named and sourced: it is GEOMETRIC, it lives at
the GRAZING ordinates, and no τ choice removes it** — the same `∝ sin ω`
cell-measure fact as §4.4, now with a primary-source criterion and a stated
consequence (second order → first). BMC Eq. 43 and R&L Eq. 16 **genuinely
conflict at the arc endpoints**; the absorber "resolved" it by forcing τ=½ there
and silently forfeiting Eq. 43. **That conflict is the real content of Q5.6.4,
and it is a design question, not a defect.**

#### For the seed experiment, and two citation upgrades

* `[M]` printed p. 239: *"If the M directions are ordered so that neutrons flow
  out of, but not into, the first μ cell, then `α_{1/2} = 0`. If neutrons flow
  into, but not out of, the last cell, then `α_{M+1/2} = 0`."* ⟹ **α is known
  at BOTH ends and R&L march from either**; a one-sided march is the
  ill-conditioned choice. Production marches once. **This is the concrete,
  sourced candidate for the seed work.**
* `[M]` R&L call the α recursion (their Eq. 9 = 13a, `α_{m+1/2} − α_{m−1/2} =
  −μ_m w_m`) *"a requirement commonly invoked"*, citing **K. Lathrop and
  B. Carlson, J. Comp. Phys. **1**, 173 (1966)** — a more specific primary
  source for our recursion than Alcouffe & O'Dell. Neither is local.
* ⭐ `[M]` **Eq. 28b** (their 3-D CYLINDER section — R&L DO treat cylinders,
  contra Lathrop 2000) keeps the weight `w` and the cell width `Δμ` **separate**:
  `[(1−τ)α_{m+1/2} + τα_{m−1/2}]Δμ = η²w`. **Primary-source licence for
  weight ≠ cell-measure**, which is §4.4's finding from the other side.
* ⚠ **Grant 1968, JCP **2**(4):381–402, DOI `10.1016/0021-9991(68)90044-2`** is
  the origin of BOTH the weighted-diamond ansatz AND the `[½,1]` interval. NOT
  local, closed access. Acquiring it is the user's call.

⚠ **What cuts against us, recorded honestly:** R&L's second-order result is a
**pro-midpoint** argument, and their measured accuracy is dominated by the
quadrature's second moment (5.3× at S₈) with τ worth ~11 %. But `[M]` their own
τ never leaves `[0.44, 0.56]`, capping their transient amplification at ~1.8 —
so their configuration *cannot* exhibit our 9× transient, and their silence
about it is not evidence of absence.

### 9bis.9 ⭐⭐ THE LITERATURE — Hébert's side, and a production citation that is false three ways

⛔⛔ **THIS SECTION'S CONCLUSION IS REFUTED BY §9bis.10** — it reads Hébert
without BMC's refutation of exactly Hébert's scheme. Its *factual* content
stands and matters (the section numbers, the equation identities, the
quadrature-reduction measurement, and above all the false production citation);
its inference — "⟹ `τ ≡ ½` on the cylinder is not a foreign method" — does not.

Deliverable: `scratch/q64_cylinder_closure_literature.md`. All three sources were
LOCAL; nothing acquired. **Every claim below I re-verified myself** in
`scratch/literature_ocr/Hebert(2009)Chapter3.md` (line numbers are that file's).

⛔ **Two premises of my own brief were wrong, and correcting them is the answer:**

1. **Hébert's cylinder is §3.9.3, not §3.9.4.** `[M]` sidecar line 2796 =
   *"3.9.3 The difference relations in 1D **cylindrical** geometry"*, line 2946 =
   *"3.9.4 … 1D **spherical**"*. The whole 3.418–3.439 range is **spherical**.
2. **Eqs. 3.437/3.439 are NOT a weighted recurrence.** `[M]` line 3053 reads
   `φ_{n+1/2,i} = 2φ_{n,i} − φ_{n−1/2,i}` — that IS `τ ≡ ½`, being (3.431)
   rearranged for the sweep. The cylinder's azimuthal counterpart is (3.412)/
   **(3.414)** at line 2919: `φ_{p,q+1/2,i} = 2φ_{p,q,i} − φ_{p,q−1/2,i}`.
   **Hébert defines no τ anywhere in chapter 3, in either geometry.**

`[M]` the two defining diamonds, verbatim from the sidecar:
* line 2876, **§3.9.3 CYLINDER, Eq. (3.406)**:
  `φ_{p,q,i} = ½(φ_{p,q,i−1/2} + φ_{p,q,i+1/2}) = ½(φ_{p,q−1/2,i} + φ_{p,q+1/2,i})`
  — the second equality is the diamond **in the azimuthal index q**.
* line 3010, **§3.9.4 SPHERE, Eq. (3.431)**: the identical statement in the
  polar index n.
* line 2839, **§3.9.3, Eq. (3.399)**: `α_{p,q+1/2} = α_{p,q−1/2} + W_{p,q} μ_{p,q}`
  — the α recursion, cylinder, indices (level, azimuthal). Sign opposite to ours;
  the crosswalk is `conventions/normalization.rst §normalization-alpha-crosswalk`.

#### ⛔⛔ A PRODUCTION CITATION THAT IS FALSE THREE WAYS — and it is the ROOT CAUSE

`orpheus/sn/sweep/pole_angular_closure.py` module docstring, lines 10–14, claims:

> *"The production closure is the per-cell Morel--Montry **weighted**-DD angular
> recurrence of **Hébert Eqs. 3.437 / 3.439**,
> `φ_{n+1/2,i} = (φ_{n,i} − (1−τ_n)φ_{n−1/2,i})/τ_n`"*

1. 3.437/3.439 are **§3.9.4 SPHERICAL**, cited as the source for a module whose
   whole point is that one strategy covers **both** curvilinear arms.
2. They are **not weighted** — they are `2φ_n − φ_{n−1/2}`, i.e. `τ ≡ ½` exactly.
3. Hébert writes **no τ**, so the weighted formula is attributed to an author who
   does not state it. Its real sources are BMC Eq. 43 / Lathrop Eq. 23 (and
   before them Morel & Montry 1984, **not in the local library**).

⭐ **This is the causal account of the whole affair.** Someone read "Hébert 3.439
is the angular recurrence" (true), transcribed it as "the *weighted* recurrence
with a τ" (false), and then needed a τ — which arrived from BMC's
**µ**-barycentric formula, written for a **different quadrature**. The citation
defect and the design defect are one defect. Fixing the citation is mandatory
under EVERY outcome, because "weighted-DD … of Hébert 3.437/3.439" is false
whether we keep a weighted τ or not.

#### The reconciliation: each source's recipe is the PREDICATE on ITS OWN rule

| source | geometry | its quadrature | its τ |
|---|---|---|---|
| **Hébert** §3.9.3 (3.406/3.412/3.414) | cylinder | **equispaced ω, weight = the cell's ω-measure** — `[M]` (3.370)/(3.375) reduce exactly to `ω_q = π(q−½)/M` (≤4.4e-16, N=4…32) | **`τ ≡ ½`**, his default and only scheme |
| **Hébert** §3.9.4 (3.431/3.437/3.439) | sphere | his own | `τ ≡ ½` |
| **BMC** Eq. 43 (sphere) / **Eq. 74** (cylinder, printed p. 160, *"analogous to Eq. (42)"*) | both | cumulative-**weight** edges, Eq. 52, per-level renormalised to `2√(1−ξ²)` | barycentric in the **radial cosine µ**; ω never appears in (53)/(74)/(75) |
| **Lathrop** Eq. 23 | ⚠ **sphere only** — `[M]` "cylind" count = **0** | cumulative-weight in µ | barycentric in µ; `τ ≡ ½` is a **deficiency** |

⭐⭐ **The decisive fact: Hébert's own recommended cylindrical quadrature places
every ordinate at the exact ω-midpoint of an equal-ω cell weighted by that
cell's ω-measure.** That is *our* `folded_product` rule. So on the rule we ship,
`τ ≡ ½` **IS** BMC-43 / Lathrop-23 — the barycentric predicate — merely read in
the variable the cells are equal in. And the mirror: `[M]` BMC's Eq.-52 cosine
partition on that same rule drives τ to `[−0.33, 1.33] → [−1.18, 2.18] →
[−2.86, 3.86]` at M = 8/16/32, i.e. **outside the admissible `[0,1]`, worse with
refinement** — an independent reproduction of §4.4, with the mechanism named:
an equal-ω cell's µ-measure ∝ sin ω, spreading `5.0× / 10.2× / 20.4×` (≈ 2M/π).

⟹ **`τ ≡ ½` on the cylinder is not a foreign method** (P-F's charge). It is this
project's own predicate realized on the quadrature its own primary source
recommends. This is `lessons-L48` — *take the PREDICATE, not the recipe* —
applied to the exact case that produced the campaign's own ruling.

#### ⭐ LATHROP'S DEFICIENCY IS BOUNDED, AND IT IS ABSENT HERE — measured

Lathrop condemns `τ ≡ ½` twice (printed pp. 245 and 249–250), and the argument
runs through the **NODES**: `δ = 0` *implies* the ordinates are cell midpoints,
and midpoint nodes give **Eq. (63)** `c = ⅔(1 − 1/N²) ≠ ⅔`, missing the
diffusion condition. It is an argument against *choosing your nodes* to make
τ ≡ ½ on a cumulative-weight-in-µ partition.

Our cylinder nodes are fixed by the quadrature and are midpoints **in ω**, not in
the cosine — and `[M]` 2026-08-11 they satisfy the condition **exactly**:

| rule | rel err, 2nd moment vs `w_gl sin²θ · π/2` | 0th moment |
|---|---|---|
| `folded_product(2, 4…64)` | `2.1e-16 … 0.0e+00` | `0.000e+00` |
| `folded_product(4, 4…64)` | `1.2e-16 … 2.0e-16` | `0.000e+00` |

The mechanism: `cos²ω = ½(1 + cos 2ω)`, and the ω-midpoint rule annihilates the
`cos 2ω` term exactly on an equispaced grid over the half-circle. **Lathrop's
`⅔(1 − 1/N²)` deficit is a property of midpoints in the COSINE; ours are
midpoints in the ANGLE, and they are exact.** This is also why probe D found P1
τ-blind: it grades the nodes, and the nodes are already right.

#### What the literature does NOT say — record it, do not paper over it

* **No source prescribes any clamp or limiter on τ.** Range is `[0,1]`
  everywhere; Lathrop *derives* it from node ∈ cell. Positivity is analysed only
  for the **spatial** diamond (Hébert 3.387–3.389). ⟹ the absorber retirement
  stands on the literature. ⚠ **And §9bis.5's `τ ≥ ½ ⟺ non-amplifying` argument
  is in NO source** — BMC even print the amplifying recurrence (their Eq. 54)
  and say nothing about it. It is an ORPHEUS observation, measured and true;
  label it as ours, never cite it as literature.
* `[M]` **BMC Eq. (47) CONFIRMED**: S₂ Gauss–Legendre gives
  `τ₁ = 1 − 1/√3 = 0.42265 < ½`, and exactly **half** the set is sub-½ at every
  order — so a universal `τ ≥ ½` was never the admissible range in either arm.
* **Neither source discusses the ω-vs-µ choice explicitly.** The answer is
  reconstructed from which quadrature each one assumes, as above.
* ⚠ **Three published typos**, verified on the rendered pages: BMC Eqs. (50) and
  (52) both print self-referential RHS subscripts (`m+1/2` must be `m−1/2`), and
  BMC printed p. 156's coordinate line prints `η = sinθ cos ω` where it must be
  `sin ω` (else `η ≡ µ`). Do not "correct" our code to match those as printed.
* ▶ **Not local, and it is the one place a clamp/positivity argument could
  live: Morel & Montry (1984), TTSP 13(5):615** — the primary τ source that both
  BMC and Lathrop only cite. Also absent: Alcouffe & O'Dell (Hébert ref. [36]),
  Reed & Lathrop (Lathrop ref. 7). **Acquiring these is the user's call.**

---

## 9. TRAPS FOR THE NEXT SESSION

1. **`τ` has three meanings and `β` two** — closure weight / optical depth /
   critical half-thickness in mfp; BMC's scalar contamination / Lathrop's α
   defect. Nomenclature table: `angular_differencing.py` module docstring.
2. **A math symbol has three spellings plus its number** — `tau_raw`, `τ_raw`,
   `\tau_{\rm raw}`. A two-spelling grep reported clean and missed a live
   falsehood; it was found by grepping `tfrac15|tfrac45`. And `absorber` is
   also a *material*, `clamp` also a GMRES restart clamp — triage concept-grep
   hits by MEANING (now in `coding-standards`).
3. **`folded_product(·, 4)` is M = 2 and cannot see a partition change.** Use
   `n_φ ≥ 6`.
4. **1-group homogeneous fixtures cannot see τ** (`k_eff = k_∞`).
5. **`nx = 80` mixes spatial and angular error.** Use `nx = 320` for any
   angular claim.
6. **`python -O` strips a standalone probe's asserts** (L42). Use prints, or
   run inside pytest.
7. **NEVER `git checkout <path>` / `restore` / `stash` / `clean`** on any
   tracked path — the tree carries irrecoverable uncommitted-by-policy state.
   Mutation-test by in-process monkeypatch; compare against history with
   `git show <hash>:<path>` (that is how §4.6's OLD snapshot was obtained).
8. **`pole_angular_closure` has no `automodule`** anywhere in `docs/api/`, so no
   Sphinx severity can see its Python-domain refs. `tools/check_docstring_xrefs.py`
   is the only gate.
9. The `[½,1]` absorber's retirement is independently justified by the
   literature (no source prescribes any limiter; BMC's own S₂ gives
   `τ₁ ≈ 0.4226 < ½`) **and by ν-closure**. That part of the carve does not
   depend on the partition choice — but the `τ_0 = 0` endpoint pathology it was
   blocking is real, and under partition (1) it needs a guard.
