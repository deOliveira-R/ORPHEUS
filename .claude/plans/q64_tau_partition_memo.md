# Q5.6.4 — the cylinder τ partition: what was tested, what was assumed, what it showed

**Status: LANDED ON THE BRANCH, ONE TEST RED, THE DESIGN'S STRUCTURAL PREMISE
REFUTED.** Written 2026-08-11 as a compaction-proof handoff. Self-contained on
purpose: assume the reader has NO memory of the session that produced it.

Companion documents (all committed):
* `.claude/plans/quadrature_machinery_campaign.md` — the campaign plan; its
  Q5.6.4 item carries the same measurements inline with the refutations.
* `scratch/q64_red_triage.md` — the per-test red triage, §D is the open red.
* `scratch/q64_tau_edge_convention_literature.md` — the literature deliverable.
  **Line 943 is the sentence that refutes the design.**
* `scratch/q64_absorber_retirement_blast_radius.md` — the retirement audit.
* `orpheus/derivations/discrete/sn/angular_differencing.py` — the analysis
  module: the P0–P4 predicate ladder, ν-closure, both β's, the τ/β nomenclature.

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
| **P-A** | α and τ both reference "the boundary between azimuthal cell m and m+1", derive it independently, and DISAGREE | ✅ **TRUE** — and it is a genuine Cardinal-Rule-2 defect. `[M]` §4.2 |
| **P-B** | α lives at the **geometric** arc half-angle boundary `ω_{k−1/2}` | ⛔ **FALSE** — see §3. This is the fatal link. |
| **P-C** | ∴ τ must be built on the geometric arc partition | ⛔ **INVALID** — follows only from P-B |
| **P-D** | the `[½,1]` absorber exists to compensate the chord partition's stretched end cells | ⚠ **UNPROVEN, possibly inverted** — see §3 |
| **P-E** | ν-closure `= 1.000000` validates the arc partition | ⛔ **INVALID AS VALIDATION** — see §5.3. It is necessary, not sufficient. |
| **P-F** | `τ ≡ ½` is the principled M-M value expressed in ω | ⛔ **REFUTED** — it is Hébert's *diamond* scheme (Eqs. 3.406/3.431), a different method, and `[M]` it fails ν-closure by 16.5 % |
| **P-G** | BMC's β is the diffusion-limit oracle that adjudicates a partition | ⚠ **TRUE in general, BLIND here** — see §5.2 |
| **P-H** | the literature's cumulative-weight edges can be transplanted to the cylinder | ⛔ **REFUTED** — violates P3, worsening with refinement. §4.4 |

⭐ **P-A survived and P-C did not.** The defect is real; the fix was aimed at the
wrong target. A fresh session should NOT conclude "the chord partition was fine
all along" — it should conclude "α and τ must share a partition, and we have not
yet identified which one α actually uses."

---

## 3. ⛔ HOW P-B IS FALSE — the sentence that refutes the design

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

**(3) RECURSION-DEFINED — NOT YET MEASURED. This is the next step.**
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
