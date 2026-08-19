# Campaign — the cylinder's angular error converges at its design order,
# and the scheme that achieves it can be SHOWN to be the better one

**Issue of record: [#235].** Opened 2026-08-12 out of Q5.6.5
(`.claude/plans/archive/q64_tau_partition_memo.md` closed; evidence
`scratch/q65_endpoint_defect_findings.md`). Branch
`refactor/operator-strategy-layers` @ `bea6a367`.

> ⚠ **#235's TITLE names a MECHANISM** ("a 2-D (η,φ) angular closure").
> That is the leading hypothesis, not the goal — see §1. Its body is
> bannered STALE; the 2026-08-12 comment is the live state.

---

## 1. The GOAL, separately from any proposed MEANS

**Goal (outcome, durable).** A curvilinear SN solve's angular error falls
at the scheme's design order instead of sitting on a floor — **and** the
project can demonstrate, with an instrument it trusts, that the shipped
angular differencing is better than the alternatives.

⭐ **The second clause is currently the harder one, and it is the reason
this is a campaign rather than a patch.** Six plausible instruments have
been measured and all six are dead (§3). We can propose schemes; we
cannot currently rank them.

**Proposed means** (as of 2026-08-12, all UNVERIFIED hypotheses):
* (M1) a 2-D (η,φ) angular closure — #235's own proposal;
* (M2) a different 1-D scheme (characteristic / exponential / step in ω)
  rather than a weighted diamond;
* (M3) a Galerkin/spectral treatment in ω on the folded arc's realisable
  `{cos mω}` basis;
* (M4) impose BOTH angular endpoints (a BVP in angle) rather than shooting
  from one.

**Done when:** an angular-refinement ladder on a *discriminating* fixture
shows the shipped scheme at its design order with no floor in the
n_φ = 16–64 band, and a documented instrument suite ranks it above the
alternatives — **including plain diamond**, which currently beats it.

---

## 2. THE PREMISE TABLE — every claim, with status

All `[M]` rows measured 2026-08-12 unless dated otherwise. Configurations
are stated because a number without its fixture is not reusable.

| # | claim | status |
|---|---|---|
| P1 | The cylinder's aniso MMS has no clean O(h²) spatial window; the error is an angular floor | `[M]` #229 (CLOSED — a measurement RECORD, not open work) |
| P2 | The floor scales with **azimuthal** n_φ, flat in polar n_mu | `[M]` #229 |
| P3 | Post-fold + post-unclamp the floor is `3.1503e-3 → 1.1326e-3 → 3.2611e-4` at n_φ = 8/16/32, nx=80 (ratios **2.78× / 3.47×**) — ~10× below #235's body ladder, and scaling BETTER than the pre-6.3 2.58×/2.38× | `[M]` probe I |
| P4 | ⭐ **Plain diamond τ≡½ BEATS the shipped M-M τ by 3.0× / 2.5× at n_φ = 16 / 32** (shipped wins only at n_φ=8). Each variant re-solved to its own fixed point; all 12 solves converged | `[M]` probe I |
| P5 | That reproduces the cost the Q5.6.4 carve explicitly RATIFIED ("~1.8–2× worse at n_φ=16/32/64; principled ≠ more accurate") | `[M]` — `1.1326e-3` = 1.67× the gate docstring's pre-6.4 `6.782e-4` |
| P6 | The n_φ=32 shipped-vs-**reversed-garbage** gap is only **1.11×** ⟹ this fixture STOPS DISCRIMINATING at the orders we care about | `[M]` probe I |
| P7 | The level's angular march is affine with linear part **exactly `(−1)^M·I`** (gain to 2e-12, pointwise in (g,i)); `Π(1−τ)/τ = 1` exactly | `[M]` probe B — a CONSEQUENCE of an already-gated identity, confirmed through the production kernel |
| P8 | ⟹ `A(M) = 2.41…9.44` is purely **transient**; end-to-end gain is 1 | `[M]` |
| P9 | Both angular endpoints are independently solvable (α=0 at each) and production **marches both**; the recurrence's trailing face is consumed by nobody | `[M]` probe A + explorer map |
| P10 | At the pole the recurrence MEETS the isotropy requirement: `r/seed` = 0.0088 (even M) / 2.0088 (odd M), improving with N | `[M]` probe H |
| P11 | At a **vacuum** surface `L∞(D)` diverges spatially (1.52e-1 → 4.01e-1 over nx 6→96, always the outermost cell); **reflective** converges (ratios 0.95/0.98) and the hot cell migrates inward ⟹ a grazing-discontinuity angular boundary layer | `[M]` probes D/E |
| P12 | BMC 2010: the weighted τ is first-order **diffusion-limit** consistent, plain diamond leading-order only | *literature, unverified numerically here* — **this is what Phase 0 tests** |
| P13 | The reconciliation of P4 with P12: the aniso-cyl MMS is **not in the diffusion limit** | **HYPOTHESIS**, 2026-08-12 — Phase 0 |

### ⛔ Refuted premises, kept in place (plan-authoring §3)

* ⛔ #235 body: "the most-inward ordinate sits on −sinθ exactly (τ_raw=0 →
  a structural ÷0 that keeps the cylinder clamp)". **False since the σ_y
  fold** — `[M]` nodes are strictly-interior arc midpoints
  (`ω/π = [0.875, 0.625, 0.375, 0.125]` at n_φ=8).
* ⛔ #235 body: "unclamp → NaN". **False** — the `[½,1]` absorber was
  RETIRED at Q5.6.4; production runs unclamped on both arms.
* ⛔ "τ is the last unexamined term" (Q5.6.5's founding premise) — the
  seed and τ are both settled; see the memo's §0a and F1/F4.

---

## 3. ⛔⛔ THE INSTRUMENT GRAVEYARD — six dead. Do NOT re-derive.

| instrument | why it is dead |
|---|---|
| P1's `c = Σ w η²` | **τ-blind** — bit-identical under garbage τ |
| BMC's contamination β | **τ-blind** (τ-free by construction) |
| Lathrop's β | **τ-blind** |
| ν-closure residual | reports only which **chart** it was run in |
| trajectory-resolvent cross-check | **reference-limited** at ≈3e-2 — refining the SUT makes agreement WORSE ([[lessons-L49]]) |
| **the endpoint defect `D`** | reference-free, pointwise, genuinely τ-loaded (ranks garbage 2.6–45× above shipped) — **but uncorrelated with accuracy**: Pearson r on log = **+0.75 → +0.26 → +0.06** at n_φ = 8/16/32 vs the analytic MMS, 0/4 ranks agreeing at the two finer orders |

⭐ **The pattern, and the campaign's central methodological problem:** every
one was plausible, and every one measured **something other than what it
was credited with**. Five were caught by a blindness probe; `D` was caught
only by correlating against an independent reference. **No instrument
ranks a scheme in this campaign until it has passed that correlation.**

**What survives:** the manufactured solutions
(`orpheus.derivations.continuous.mms.sn`) — analytic reference, hence
structurally independent (`vv-principles` L11). Limitations, stated: one
fixture family, floor-not-rate on the cylinder, `slow`, and per P6 it
stops discriminating by n_φ=32.

---

## 4. Phases

### ⏱ Phase 0 — GET AN INSTRUMENT (in flight, 2026-08-12)

**Goal.** The campaign can tell a better angular scheme from a worse one.
Until this closes, every design phase is working blind.

Four perspectives dispatched in parallel:

| perspective | question | deliverable |
|---|---|---|
| `numerics-investigator` | **#319's flux-dip experiment**: does the centre-dip DECAY RATE vs optical thickness split M-M from plain diamond? Resolves P4-vs-P12 | `scratch/q68_flux_dip_discriminator.md` |
| `literature-researcher` | What curvilinear angular differencing exists beyond weighted diamond? What do production codes ship? Is the floor a NAMED phenomenon? Is the P13 trade documented? | `scratch/q68_curvilinear_angular_differencing_survey.md` |
| `cross-domain-attacker` | Is the 1-D η-march the right STRUCTURE? (connection/parallel transport, characteristics-in-angle, spectral in ω, 2-D cell complex, BVP-not-IVP) | `scratch/q68_angular_frame_attack.md` |
| `test-architect` | The instrument SUITE + the mechanical anti-blindness acceptance test + a fixture that still discriminates at n_φ=16–64 | `scratch/q68_angular_instrument_design.md` |

**Done when:** at least one instrument has passed a documented
correlation-against-independent-accuracy test over **more than four**
deliberately-varied schemes, and the P4/P12 tension is resolved either way.

⚠ A **refutation is a success here.** If the flux-dip decay does NOT split
the two schemes, that refutes the angular-consistency reading the landed
carve rests on — and that is more valuable than a confirmation.

### Phase 1 — decide the STRUCTURE (blocked on Phase 0)
Choose among M1–M4 (or a synthesis) on the evidence. Not before.

### Phase 2 — build it
Surgical: main agent writes, user steers (`delegation.md` — no
`method-implementer` for operator-algebra carves).

### Phase 3 — the retirement + the corpus
Whatever is replaced gets the 3-search audit; theory pages follow.

---

## 4bis. ⭐⭐⭐ PHASE 0 RESULT (2026-08-12) — the closure measures
## midpoint-ness in the WRONG VARIABLE, and it costs exactly one order

Two perspectives converged on this independently (literature + structural
frame attack); every link below is verified by the main agent directly.

**The chain, each link measured:**

1. `[M]` Production's angular cell edges are **exactly equispaced in ω**
   — `max|e_ω − linspace(π,0)| = 4.4e-16` at M = 4/8/16.
2. `[M]` Each ordinate is **exactly the midpoint of its own cell in ω** —
   `max|ω_m − ½(e_m+e_{m+1})| = 1.6e-15`. Equivalently
   `ω_k = (k+½)π/M` to 2e-15 (⟹ the nodes are the roots of `T_M`;
   Chebyshev–Gauss).
3. ⟹ **τ measured in ω would be exactly ½.**
4. `[M]` But the shipped τ is the barycentric coordinate in
   **η = sinθ·cos ω** — a *nonlinear* reparameterisation — so it reads
   `[0.2506, 0.7494]` at M=16, deviation → **¼**, not 0.
5. `[M]` Hence `(τ−½)/w` is **UNBOUNDED**: `0.96 → 1.98 → 3.99 → 8.00` at
   M = 4/8/16/32 (exactly `M/4`), and on the sphere `0.43 → 0.88 → 1.76
   → 3.53` at N = 8/16/32/64. `max|τ−½|` converges to ≈0.25 / ≈0.11 —
   it does **not** go to zero.
6. **Lathrop (2000) NSE 134:239–264 Eq. (30)** — LOCAL, verified against
   the scan — with `δ = 2τ−1` the truncation is `O(δΔμ + Δμ²)` and
   *"only with μ_m = μ̄ (δ = 0) is the truncation order O(Δμ²)"*; §IV:
   *"Even if Reed's 'optimum' weighted diamond difference relations are
   used, the truncation errors are O(δΔμ + Δμ²)."*
   ⟹ **every weighted diamond is FIRST-order in angle.** This is the
   same criterion as R&L Eqs. 15/16, which §0a already listed as a LIVE
   instrument — it was never run against the shipped τ.
7. `[M]` **The predicted consequence, measured at spatially-converged
   nx = 320** (aniso-cyl MMS, each τ re-solved to its own fixed point,
   all converged):

| τ | n_φ=8 | 16 | 32 | 64 | local orders | global |
|---|---|---|---|---|---|---|
| shipped (M-M weighted) | 3.1281e-3 | 1.1078e-3 | 2.8285e-4 | 7.1658e-5 | 1.50, 1.97, 1.98 | **1.83** |
| plain diamond τ ≡ ½ | 3.4258e-3 | 3.4907e-4 | 3.9321e-5 | **9.1485e-6** | 3.29, 3.15, 2.10 | **2.88** |

⟹ a **~one-order gap**, and diamond is **8×** better at n_φ=64.

### ⛔ CORRECTION (same day) — that 8× is FUNCTIONAL-specific; the honest number is 3.2×

The `test-architect` pass found that the **graded functional** matters more
than the fixture, and it applies to the table above. Two verifications by
the main agent:

* ✅ **The closure is exact on `span{1, η}` BY CONSTRUCTION** — pure algebra,
  no run: τ *is* the barycentric coordinate, so
  `τη₊ + (1−τ)η₋ ≡ η_m` and `τ·1 + (1−τ)·1 ≡ 1`. **A seventh blindness
  mechanism**, and the one that explains the graveyard: any instrument whose
  signal lives in `span{1, η}` is structurally τ-blind.
* ✅ **The flagship fixture's exact solution IS in that kernel** — verified
  from the source, not inferred: `sn.py:3803`
  `psi_exact(r, eta_n) = self.A(r) + self.B(r) * eta_n`, i.e. exactly affine
  in η (`[M]` residual ≤ 4.4e-16 at r = 0.3/0.7/1.3/1.9).

`[M]` Re-grading the SAME solves on the **angular** flux instead of the
scalar flux, nx=320:

| functional | shipped order | diamond order | ratio s/d at n_φ = 8/16/32/64 |
|---|---|---|---|
| scalar flux (as first reported) | 1.83 | 2.88 | 0.91 / 3.17 / 7.19 / **7.83** |
| angular flux | 1.85 | 2.44 | 0.94 / 1.69 / 2.40 / **3.24** |

⟹ **The DIRECTION survives** (diamond converges faster in both; the
advantage grows monotonically with n_φ). **The MAGNITUDE was inflated
≈2.4×** by the scalar-flux functional, and the order gap is **≈0.6, not
≈1.05**. Mechanism: `Σ_n w_n cos(m ω_n) = 0` — the quadrature identity that
makes `φ_exact = A(r)` closed-form is the same identity that annihilates
the angular defect in φ. Every gate in `tests/sn/verification/mms/` grades
the scalar flux.

### ⛔⛔⛔ RESOLVED, AND IT INVERTS THE CONCLUSION — the fixture was the difference, not nx

`[M]` 2026-08-12, main agent, running the **test-architect's** harmonic
fixture (`CylHarmonicMMS`, `h₁+h₂`, σ_t=5, c=0.5) across nx, ψ-L2, ratio
`err_shipped / err_diamond` (**< 1 means SHIPPED wins**), all converged:

| n_φ | nx=80 | nx=160 | nx=320 |
|---|---|---|---|
| 16 | 0.746 | 0.743 | **0.742** |
| 32 | 0.910 | 0.896 | **0.892** |
| 64 | 1.052 | 1.031 | **1.022** |

**Spatially converged to three digits — nx is NOT the confound**, and these
are the exact reciprocals of the test-architect's `1.34 / 1.10 / 0.95`, so
two independent measurements agree. The difference is **the FIXTURE**:

| fixture | angular content | verdict at n_φ = 16 / 32 / 64 |
|---|---|---|
| shipped `build_cylindrical_anisotropic_mms_case` | `A + Bη` — **IN the closure's kernel** | diamond better, 1.7× / 2.4× / **3.2×**, growing |
| `CylHarmonicMMS` `h₁+h₂` | **EXCITES the closure** | **SHIPPED better by 34 % / 12 %**, equal at 64 |

⟹ ⛔ **MY CONCLUSION OF THIS MORNING IS REFUTED.** "Plain diamond converges
faster / the shipped scheme is one order down / re-posing τ into ω is a
live candidate on accuracy grounds" was a property of **a fixture that is
blind to the operator under test**. On a fixture that actually excites the
angular closure, **the shipped Morel–Montry τ WINS at the orders that
matter** and is level at n_φ=64.

⭐ **What survives, and what does not:**

* ✅ Lathrop Eq. (30) is still a true statement about the *equation's*
  truncation order. ✅ `(τ−½)/w` is still unbounded (`= M/4`). ✅ The
  affine-fixture ladder (1.83 vs 2.88 scalar / 1.85 vs 2.44 angular) is
  still a correct measurement **of that fixture**.
* ⛔ The **inference** from those to "the shipped scheme is worse" does
  **not** hold. A truncation-order statement about the equation is not a
  ranking of solution accuracy, and the fixture that appeared to confirm
  it cannot see the closure at all.
* ⭐ **Q5.6.4's landed carve is now vindicated by MEASUREMENT**, not only
  by Morel & Montry's appendix — on the only fixture built to see it.

⭐⭐ **The campaign's real prize, and it is a methodological one:** the
manufactured fixture that has graded this campaign for months lies in the
**kernel of the operator it grades**. Every accuracy conclusion ever drawn
from it *about the angular closure* is void — including one I drew, with
confidence, three hours before this measurement. `Σ_n w_n cos(mω_n) = 0`
and `psi_exact = A(r) + B(r)η` are the two identities that did it.

### ⚠ ORIGINAL OPEN-CONTRADICTION NOTE (superseded by the block above; kept per plan-authoring §3)

The `test-architect`'s **purpose-built** fixture (harmonics
`h_m = Re[(η+iξ)^m]`, graded on ψ) reportedly does **NOT** resolve shipped
vs diamond (`1.34 / 1.10 / 0.95×`, sign-flipping). The table above, on the
**shipped** fixture graded on ψ, gives `3.24×` monotone. Both cannot be a
property of "the schemes"; at least one is a property of its fixture.
**Unresolved as of 2026-08-12.** Read `scratch/q68_angular_instrument_design.md`
§7 and reconcile the configurations before either number is used.

### ⭐ The synthesis that dissolves the apparent conflict with #319

`[M]`/derived: the diffusion limit's angular content **is** `span{1, μ}`
(`ψ → φ/4π + 3μJ/4π`), which is exactly the closure's kernel. So the
barycentric τ is exact there and τ≡½ is not ⟹ **#319's flux-dip test is
structurally biased toward the incumbent.** That does not make it
worthless — it makes it a **CONSTRAINT** (does the scheme preserve the
diffusion limit? pass/fail), never a **RANKER**.

⟹ The two instruments answer different questions, and the trade is real:
**M-M passes the constraint and loses the convergence ranking; τ≡½ may win
the ranking and fail the constraint.** That vindicates Q5.6.4's
literature route as sound rather than a fallback — while showing the cost
it ratified is larger and better-characterised than recorded.

⛔ **AND IT REFUTES #235's OWN PREMISE.** Neither scheme flatlines, so
there is **no 2-D (η,φ) obstruction** visible: the angular error
converges. The "floor" the campaign has been chasing was measured at
**nx = 80**, where `[M]` both schemes read ~1.5 because the SPATIAL error
dominates — diamond's fine-end order collapses to **0.06** there, purely
from hitting the spatial floor. **The measurement the whole issue rests
on was spatially contaminated.**

⚠ **What this does NOT settle, and why nothing may be switched yet.**
M&M/BMC chose the η-weighting to buy **first-order diffusion-limit
consistency**. The trade on the table is therefore *one order of angular
accuracy away from the diffusion limit, in exchange for correct behaviour
in it*. #319's flux-dip experiment (running) tests whether that purchase
is real. **Do not recommend τ ≡ ½ until it lands.**

Note also the arms differ: `[M]` GL nodes on the **sphere** are NOT their
cells' midpoints (`1.09e-2 / 2.98e-3 / 7.73e-4` at N = 8/16/32, ~1/N²), so
δ ≠ 0 is *forced* there and weighting is unavoidable. The cylinder's arc
rule is the equal-interval case where δ = 0 is available for free. **A
per-arm answer is likely the right one.**

## 4ter. ⏹ PHASE 0 CLOSED (2026-08-12) — the shipped scheme is VINDICATED,
## and the campaign's problem was never the instrument's reference

### #319's prediction is REFUTED on the axis it named

`[M]` `numerics-investigator`, 251 solves, **all `converged` AND
`fully_converged`, zero warnings, `max|balance_defect| = 0`**, every variant
re-solved to its own fixed point. Record:
`scratch/q68_flux_dip_discriminator.md`, raw
`$CLAUDE_JOB_DIR/tmp/q68_results.jsonl`.

* At constant cells-per-mfp the angular-consistency defect is **constant to
  four significant figures** over `Σ_t·R = 5…50` — for the shipped τ, for
  `τ ≡ ½`, **and for the garbage control**. Fitted decay rate `0.000` for
  all three. The dip proper decays for (A) but **changes sign** mid-range
  for (B)/(C), with a *steeper* pre-crossing slope — the opposite of
  "persists".
* ⛔ Morel & Montry themselves (p. 13) call the diamond dip *"much more
  severe for this optically thin problem"* — **#319's framing contradicted
  its own primary source.**
* ⛔ **#319's literal design self-destructs**: at fixed `c`, `Σ_a·R` grows
  with thickness, the interior becomes a source plateau, and at
  `c = 0.99, Σ_t·R = 100` all three schemes agree to three significant
  figures. A naive run reports "equal ⟹ REFUTED" for a reason unrelated to
  angular differencing. Correct scaling: `Σ_t = 1/ε, Σ_a = ε, Q = ε`.

### …but the angular-consistency reading SURVIVES — the axis was `h → 0`

* `[M]` the shipped τ's defect vanishes at **exactly first order in h**
  while `τ ≡ ½` **saturates**: ratio `3.2× → 8.8× → 21× → 46× → 98× →
  204×` over 2→64 cells/mfp.
* `[M]` a λ-continuum `τ(λ) = λ τ_MM + (1−λ)/2` has its optimum converge to
  **λ = 1 exactly** (`0.748→0.993` and `1.300→1.001`, **two independent
  instruments**), and the anomaly is **linear in M&M's β**, whose zero is
  at the shipped τ by their own identity.
* This is verbatim M&M p. 16: consistency is recovered *"as the spatial
  mesh is refined"*. Their Fig. 4/6 reproduced to 2.5 %
  (`c_s(r_0) = −1.3157` vs their ≈ −1.35).

### The #235 reconciliation, and it is two-part

1. ✅ **Verified by the main agent from the case's own attributes**: the
   aniso-cyl fixture is `sigma_t = 1.0`, `sigma_s = 0.5`, `radius = 5.0`
   ⟹ **`c = 0.5`, `Σ_a·R = 2.5`**. Half of every collision is an
   absorption — **it is nowhere near the diffusion limit**, which is where
   the weighted τ's advantage lives.
2. `[M]` Even *in* the diffusion limit the benefit **decays with angular
   order and inverts**: sphere `14× (S2) → 1.4× (S4) → 0.9× (S8/S16)`;
   cylinder `5.4× (n_φ=8) → 2.3× (n_φ=16)`. **That is #235's own
   crossover**, reproduced on a different fixture, observable and geometry.

⟹ Combined with §4bis's kernel finding, the flagship fixture failed the
campaign **twice over, independently**: its exact solution is in the
closure's kernel, *and* its material is far from the regime where the
scheme's advantage exists.

### A cheap reference-free τ instrument DOES exist — on the sphere only

`[M]` M&M's β evaluated with the **τ-implied** edge cosines: solve-free,
**exactly zero** for the shipped τ at S2/4/8/16, and it catches `τ→1−τ`
(the Mode-12 reflection the membership / fold-box / reversal gates are all
blind to). ⛔ **Identically zero for both schemes on the folded cylinder at
every n_φ — refuted there.** This does not contradict the graveyard's "BMC
β is τ-blind": a β built from the *standard* partition edges is blind **by
construction** — that substitution *is* M&M's β = 0 proof.

### Main-agent verification status — stated because it is uneven

* ✅ **Verified directly**: the fixture's `c = 0.5 / Σ_a·R = 2.5` (read from
  source); the committed gate runs — `12 passed in 28.6 s` under
  `python -O -m pytest`.
* ⚠ **NOT independently reproduced**: the `h → 0` separation and the λ→1
  optimum. My own quick check was **starved** (`ρ = 0.999553`, needed
  ≈61 694 iterations, I budgeted 6 000) and used a wrong analytic reference
  (`φ_∞` does not hold at `Σ_a·R = 0.08`, which is leakage-dominated).
  Inconclusive — **not** a refutation. Those claims rest on the agent's
  251-solve record plus the gate. Note the #340 warning machinery caught my
  error instantly and printed the required iteration count: N4.7 paying off.

### ✅ Deliverable PROMOTED 2026-08-12 (the diagnostic is RETIRED — no twin)

`derivations/diagnostics/diag_q68_angular_diffusion_limit_consistency.py`
(12 gates, 29 s) has been split and deleted. Successors — **17 gates**:

* the shared instrument →
  :func:`orpheus.derivations.discrete.sn.angular_differencing.morel_montry_beta`,
  beside its τ-BLIND sibling `contamination_beta`, with the edge-set
  distinction documented on it (that is the whole difference);
* `tests/sn/sweep/curvilinear/test_angular_beta_identity.py` — solve-free,
  `foundation`, 14 rows incl. a NEW `τ → 1−τ` reflection leg;
* `tests/sn/verification/analytical/test_angular_diffusion_limit_consistency.py`
  — `l1`/`slow`, 3 rows — a solve-free β identity with a negative leg
(→ `tests/sn/sweep/`), the M&M Fig. 4 literature gate, the mesh-refinement
discriminator, and `test_optical_thickness_is_NOT_the_discriminating_axis`,
which **pins this report's own refutation** so the decay-rate framing
cannot be silently re-adopted (→ `tests/sn/verification/analytical/`, slow).

### ⭐⭐ PHASE 0's ANSWER, and it is not the one the phase was set up to find

The campaign's problem was **never** that its instruments lacked a
trustworthy reference. The MMS reference was analytic and sound the whole
time. The problem is that **nobody asked whether the fixture's content
intersects the operator's kernel, or whether its material sits in the
regime where the scheme's advantage exists.** Both checks are cheap,
design-time, and would have been decisive years of campaign-time ago.

⟹ **Phase 1 does not open on "replace the closure".** The shipped scheme
is the measured λ-optimum. What opens instead is: promote the gate, and
make the two entry checks mandatory for every future fixture.

## 5. Settled — re-opening any of these is re-doing measured work

1. ⛔⛔ **REFUTED 2026-08-12 — was: "τ is *not* re-posed into ω (that is
   BMC's *diamond*, leading-order only)."** The ruling conflated two
   different orderings: BMC's "leading-order" is about the **diffusion
   limit**, not about **truncation order**, and Lathrop Eq. (30) says the
   midpoint diamond is the only *second*-order option while every
   weighted diamond is first. §4bis measures the consequence: 1.83 vs
   2.88. The memo's §9bis.3 had already named the mechanism — *"the
   partition is read in ω, the WEIGHT in η"* — and then ruled the fix out.
   Re-posing τ into ω is now a **live candidate**, gated on the
   diffusion-limit trade.
2. The `[½,1]` absorber is **not** restored (no primary prescribes a
   limiter; the interval is Grant's on the **spatial** weight).
3. Hébert is **never** cited as the source of any τ — he defines none.
4. The seed is **not** the limiter (P7/P10).
5. `D` may be gated as a consistency residual; it may **never** vote on a
   scheme (§3).
6. The α-dome closure contract is enforced, single-sourced, at
   `bea6a367`. Not to be re-litigated.

## 6. Adjacent issues — read before designing

**#319** (the flux-dip / DSA angular-consistency reading — Phase 0's
critical path) · **#326** (the fold; largely landed) · **#233** (the r→0
spatial O(h), WONTFIX for DD — a *different* defect, do not conflate) ·
**#237** (degenerate-azimuthal stressor) · **#1** (Gauss-type azimuthal
quadrature) / **#3** (φ-based cell edges) — the *quadrature* side, where
this is the *closure* side · **#67** (task) the endpoint defect's fate.

## ⏸ COMPACTION POINT — insert after Phase 0 closes
Carry: the phase→commit table, the Phase-0 verdict on P13, whichever
instrument survived and its correlation evidence, and any new graveyard row.
