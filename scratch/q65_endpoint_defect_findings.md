# Q5.6.5 — the discarded endpoint condition, and the campaign's first
# reference-free τ-loaded instrument

`[M]` 2026-08-12, branch `refactor/operator-strategy-layers` @ `adb73fd5`.
Entry point was task #66 / memo `q64_tau_partition_memo.md` §0a.

**Configuration for every number below** (plan-authoring §4 — a number
without its fixture is not reusable): `_heterogeneous_2g_cylinder` from
`tests/sn/sweep/curvilinear/test_psi_half_positivity.py` — 2-region
(fuel r<1 / moderator r>1), 2-group, `Quadrature.folded_product(n_mu=4,
n_phi=·)`, `R_out = 2.0`, **vacuum** outer face unless a row says
reflective, `nx = 12` unless a row says otherwise, fixed-source
`external_source = 1`, `inner_tol = 1e-12`. **Every solve quoted here
reported `Solution.converged` True** — the #340 starvation check, run
because probe D's L∞ growth has exactly the signature a starved solve
would fake. Probes: `$CLAUDE_JOB_DIR/tmp/q65_probe_{a..g}*.py`.

---

## 0. The one-paragraph state

A curvilinear level's angular march has **two** independently solvable
endpoints (α = 0 at each), and production **computes both** — yet it
imposes only one and discards the other. The discarded condition's
violation `D` is computable from arrays that already exist at the end of
every curvilinear solve, and **nothing in the tree compares them**. `D`
turns out to be the campaign's first instrument that is simultaneously
reference-free, pointwise, and τ-loaded — and, measured confound-free, it
ranks the landed Q5.6.4 partition **first by 2.6–20×**.

⛔ **Task #66's founding premise is refuted in the form it was written.**
"The seed is the last unexamined term" — `D`'s sign does not flip with
M parity while the seed's gain does, so a seed offset is *not* what the
endpoint inconsistency measures (§F4). And R&L's "ill-conditioned
one-sided march", which §0a already flagged as *their* framing rather
than our measurement, **does not apply to this scheme**: the end-to-end
gain is exactly 1 (§F1). The two-ended march has no conditioning
motivation. Something better was found in its place.

---

## F1. The seed → trailing-face map is exactly `(−1)^M · I`

`[M]` `Π_m (1−τ_m)/τ_m = 1.000000000000000` (≤ 5e-15 over M = 3/4/8/16),
and a finite seed perturbation through the **production kernel**
(`_psi_half_grid_for_level`) gives gain

| M | measured gain | predicted | spread over (g,i) |
|---|---|---|---|
| 3 | `−1.000000000000` | `−1` | 1.55e-12 |
| 4 | `+1.000000000000` | `+1` | — |
| 8 | `+1.000000000000` | `+1` | 1.65e-12 |
| 16 | `+1.000000000000` | `+1` | 2.11e-12 |

The recurrence never mixes cells or groups, so the map is diagonal: an
r-dependent seed error `ε(g,i)` lands on the trailing face as
`(−1)^M ε(g,i)`, **pointwise and unamplified**.

⚠ **This is a CONFIRMATION, not a discovery, and the distinction matters
for what may be claimed.** `Π(1−τ)/τ = 1` is *already gated* for both
arms at `atol=1e-12` —
`test_psi_half_positivity.py::test_cylinder_recurrence_amplification_matches_the_closed_form`
and `::test_sphere_recurrence_is_neutrally_stable_across_the_level`. The
gain follows from it by one line of algebra. What the measurement adds is
that the *realized code* matches the algebra end-to-end.

⟹ **`A(M) = 2.41…9.44` is purely TRANSIENT** (interior of the level).
The endpoint amplification is exactly unit. The memo suspected this;
it is now measured through production.

## F2. One of two available boundary conditions is discarded

At **both** angular endpoints of a level the redistribution coefficient
vanishes (`α_{1/2} = 0`, `α_{M+1/2} = 0`), so the balance decouples into
a plain radial DD ODE. Production solves **both**, with one engine:

* `carlson_inward_sweep_from_source` marches the inward leg (ω = π,
  η = −sinθ_p) from the r = R inflow corner to the pole →
  `cells(p, −1)`; **this is the recurrence's seed**;
* the same engine on reversed data marches the outward leg (ω = 0,
  η = +sinθ_p), pole-continued, → `cells(p, +1)`
  (`orpheus/sn/operators/radial_characteristic.py:528-538`).

The M-M angular recurrence, marched from the seed over all M ordinates,
**also** produces the ω = 0 edge: `faces[:, M, :]`. That array is
consumed by nobody — `_MMHalfGrid.upstream_per_ordinate` drops it
explicitly (`pole_angular_closure.py:569`: *"Excludes the trailing face
… not consumed as anyone's upstream"*).

⟹ **`D := faces[:, M, :] − cells(p, +1)`** is the discarded condition's
residual: two structurally different discretizations of one point of
phase space, currently uncompared.

## F3. `D` is a pure diagnostic — there is no parameter to absorb it

Because the map's linear part is `(−1)^M I` and *both* boundary values
are supplied by physics, imposing both gives
`cells(p,+1) = (−1)^M cells(p,−1) + b(ψ)` — an **over-determination**,
a constraint on the interior solution, not an equation for a free seed.

⟹ ⛔ **`D` must NOT be used to "correct" the seed.** There is no free
parameter, and zeroing `D` would only force the marched endpoint onto the
directly-marched one with no evidence the latter is the better of the
two — the L49 reference-limited trap, one campaign step later. `D`'s role
is **a-posteriori error indication**, and that is a structural conclusion,
not a cautious one.

## F4. The sign is systematic — and it is NOT a seed offset

`[M]` signed mean of `D` is **negative at every** n_φ ∈ {6,8,12,16,24,32}
**and every level** (79–100 % of (g,i) entries negative):

| n_φ | M | signed mean (p=0) | max\|D\| | rel |
|---|---|---|---|---|
| 6 | 3 | `−3.935e-2` | 2.110e-1 | 2.94e-2 |
| 8 | 4 | `−2.941e-2` | 2.194e-1 | 3.09e-2 |
| 16 | 8 | `−5.381e-3` | 5.248e-2 | 7.38e-3 |
| 32 | 16 | `−1.753e-3` | 1.206e-2 | 1.70e-3 |

The seed's gain is `(−1)^M` (F1) and flips with parity; `D`'s sign does
not. **A seed offset of either sign cannot produce a parity-independent
`D`**, so the endpoint inconsistency is dominated by something else.
This is the measurement that refutes #66's framing.

## F5. In L2, `D` is spatially converged and second-order in n_φ

`[M]` at fixed n_φ, refining `nx` 6 → 96 leaves L2(`D`) converged to 4–5
digits (n_φ=8: `8.9242e-2 → 9.0753e-2`; n_φ=16: `1.8725e-2 → 1.9256e-2`)
and the signed mean converged to 4–5 digits. **The L2 defect is purely
angular.** Over n_φ (at nx=12, all levels): `8.2414e-2 → 1.8539e-2 →
4.8138e-3`, observed order **2.15 then 1.95**.

## F6. In L∞ it does not converge — a vacuum-surface angular boundary layer

`[M]` max\|D\| sits at the **outermost cell** at every nx and grows:

| nx | vacuum max\|D\| | growth | hot cell | reflective max\|D\| | growth | hot cell |
|---|---|---|---|---|---|---|
| 6 | 1.5236e-1 | — | 5/5 | 9.2485e-2 | — | 5/5 |
| 12 | 2.1943e-1 | 1.440× | 11/11 | 1.2299e-1 | 1.330× | 7/11 |
| 24 | 2.9625e-1 | 1.350× | 23/23 | 1.7089e-1 | 1.390× | 13/23 |
| 48 | 3.5895e-1 | 1.212× | 47/47 | 1.6247e-1 | **0.951×** | 26/47 |
| 96 | 4.0124e-1 | 1.118× | 95/95 | 1.5937e-1 | **0.981×** | 52/95 |

**Mechanism, stated as a falsifiable prediction and then tested** — at a
vacuum surface the exact angular flux is *discontinuous across grazing*
(`ψ(R, η<0) = 0`, `ψ(R, η>0) > 0`); the angular chain marches straight
through that jump, so the endpoint inconsistency at r = R is O(1) per Δω
and **concentrates** as the last cell shrinks onto the surface. A
**reflective** face has no such jump, so the growth must stop.

⟹ It does. Under reflection the growth ratio turns over to 0.95/0.98,
max\|D\| converges to ≈1.6e-1, and the hot cell **migrates off the
surface** to a fixed physical location (cell fraction 0.83 → 0.64 → 0.57
→ 0.55 → 0.55, i.e. r ≈ 1.1 — just outside the material interface).
The L∞ divergence is a property of the vacuum surface, not of the scheme.

## F7. ⭐⭐ `D` SEES τ — and ranks the landed carve first, confound-free

Three of the campaign's five instruments were bit-identical under garbage
τ (memo §9bis.4). `D` is not.

**The confound, and why the first reading was inadmissible.** Probe F
held the converged flux fixed while swapping τ — but that flux was
converged *under the shipped τ*, so it is self-consistent with it by
construction and the shipped τ wins for free. **Probe G re-solves the
whole fixed-source problem under each τ**, so each variant is measured
against its own converged fixed point.

`[M]` L2(`D`), all levels, nx=12, every solve `converged=True`:

| τ variant | n_φ=8 | n_φ=16 | n_φ=32 | order 8→16 | 16→32 |
|---|---|---|---|---|---|
| **shipped (Q5.6.4)** | **8.2414e-2** | **1.8539e-2** | **4.8138e-3** | 2.15 | 1.95 |
| diamond (τ ≡ ½) | 2.1756e-1 (2.6×) | 7.5306e-2 (4.1×) | 1.6950e-2 (3.5×) | 1.53 | 2.15 |
| reversed | 8.0132e-1 (9.7×) | 3.4405e-1 (18.6×) | 9.9433e-2 (20.7×) | 1.22 | 1.79 |
| shuffled | 1.1407e0 (13.8×) | 8.3325e-1 (45.0×) | 1.8310e-1 (38.0×) | 0.45 | 2.19 |

⭐ **Every variant preserves `Π(1−τ)/τ = 1` exactly** (permutations of the
shipped multiset; τ ≡ ½ gives ratio 1 per step), so the tree's
neutral-stability gate is **GREEN for all four**. `D` separates them by
up to 45×. It sees precisely what that gate structurally cannot: the
**per-step partition**, as against the product identity.

⛔ **A two-point slope said this was an ORDER result. Three points refute
that** — recorded here rather than dropped, per plan-authoring §3. From
{8,16} alone the shipped τ read 2.15 against the diamond's 1.53, which
looks like a recovered order; with n_φ=32 the diamond reaches 2.15 and
shuffled 2.19. **All four are asymptotically ~2nd order.** What the
shipped partition buys is (a) reaching that rate at the *coarsest rung
tested* while the others are pre-asymptotic, and (b) an error constant
2.6–20× smaller. That is still a positive reference-free result for the
landed carve — it is a **constant** result, not an order one.

## F8. ⛔ THE SPHERE DOES NOT BEHAVE LIKE THE CYLINDER — and the difference
## is a DIFFERENT MECHANISM, not a different norm

A parallel `explorer` map (`scratch/q65_alpha_ends_map.md`) found the same
discarded-endpoint gap independently and reported the **sphere**'s gap
*growing* with N — against my cylinder L2 *falling* at order ≈2. I
assumed a norm artefact (its metric max-like = my F6 L∞, which does
grow). **That assumption is wrong.** `[M]` measured directly, same
2-region 2-group fixture, `CoordSystem.SPHERICAL`, GL-N, nx=12, vacuum,
all `converged=True`:

| N | L2(D) | ratio | max\|D\| | rel | hot cell |
|---|---|---|---|---|---|
| 4 | 1.0603e-1 | — | 1.9765e-1 | 3.84e-2 | 11/11 |
| 8 | 4.3486e-2 | 2.44× | 5.9822e-2 | 1.16e-2 | **1**/11 |
| 16 | 4.8633e-2 | **0.89×** | 7.3638e-2 | 1.43e-2 | **1**/11 |
| 32 | 6.2456e-2 | **0.78×** | 7.9358e-2 | 1.54e-2 | **1**/11 |
| 64 | 7.2456e-2 | **0.86×** | 8.3420e-2 | 1.62e-2 | **2**/11 |

**Both norms agree and both RISE monotonically after N=8.** The sphere's
endpoint defect does **not** converge in angular order — it plateaus at
L2 ≈ 7e-2 / rel ≈ 1.6e-2.

⭐ And the hot cell is at the **opposite end from the cylinder's**:
`i = 1` (r ≈ 0.17), i.e. **at the pole**, not at the outer vacuum
surface. The sphere's two endpoints are the two senses of one diameter,
stitched by the pole continuation `ψ½⁺(0) = ψ½⁻(0)` — a constraint the
*radial* march imposes and the *angular* recurrence knows nothing about.
So the sphere's defect looks like a **pole-consistency** defect and the
cylinder's like a **vacuum-surface grazing-discontinuity** defect: two
mechanisms, not one. ⚠ The pole reading is one probe and is NOT
established — the hot-cell location is measured, the mechanism is a
hypothesis.

⟹ **Scope correction to F5 and F7: both are CYLINDER results.** The
"second-order reference-free indicator" and the τ ranking were measured
on the folded cylinder only. Whether `D` ranks τ correctly on the sphere
is **unmeasured**.

### ⛔⛔ F8 IS ITSELF REFUTED — `nx = 12` IS TOO COARSE FOR THE SPHERE

Kept in place per plan-authoring §3, because the wrong frame is the one a
fresh reader would otherwise re-derive — and because **both** the
explorer's "grows with N" and my own confirmation of it above were the
same measurement error, made independently, twice.

`[M]` Re-measured at `nx = 96`, same fixture, all `converged=True`:

| N | L2(D) @ nx=48 | ratio | L2(D) @ nx=96 | ratio | hot r/R @96 |
|---|---|---|---|---|---|
| 4 | 1.0671e-1 | — | 1.0678e-1 | — | 0.9948 |
| 8 | 3.1050e-2 | 3.44× | 3.0681e-2 | 3.48× | 0.9948 |
| 16 | 1.4319e-2 | 2.17× | 1.2277e-2 | 2.50× | 0.8490 |
| 32 | 1.1659e-2 | 1.23× | 5.6881e-3 | 2.16× | 0.0260 |
| 64 | 1.4812e-2 | **0.79×** | 6.0136e-3 | **0.95×** | 0.0365 |

⭐ **The turnover point MOVES OUT WITH nx** — `nx=48` turns at N=32,
`nx=96` at N=64. That is the textbook signature of the SPATIAL error
becoming the floor: once the angular error falls below it, the total stops
falling. **The sphere's endpoint defect does converge in angular order**;
at `nx = 12` the spatial floor dominates from N=8 on, which is what
manufactured the "rises with N" reading.

⟹ The two arms are **not** qualitatively different after all. What IS
measured and asymmetric: the cylinder's L2(D) is essentially
**nx-independent** (probe C: `8.92e-2 → 9.08e-2` over nx 6→96 at fixed
n_φ), so its n_φ sweep at nx=12 sits nowhere near a spatial floor and
**F5/F7 stand unaffected**. The sphere's L2(D) is strongly
nx-dependent (`4.86e-2 @ nx=12` vs `1.43e-2 @ nx=48` at N=16). Why the
spatial component differs so much between the arms is **not explained**
here.

### F8b — the pole MEETS the isotropy requirement; that hypothesis is dead too

F8's mechanism guess was a *pole-consistency* defect: at r = 0 the exact
flux is isotropic, so the two angular endpoints must coincide, which the
radial march imposes (`ψ½⁺(0) = ψ½⁻(0)`) and the angular recurrence knows
nothing about. Since the recurrence's linear part is exactly `(−1)^M`
(F1), isotropy at the pole requires its inhomogeneous part to satisfy
`b(pole) = (1 − (−1)^M)·ψ_{1/2}(pole)` — **0 for even M, 2ψ_{1/2} for odd
M**. A parity-dependent prediction no pure-truncation story makes.

`[M]` `r ≡ marched − (−1)^M·seed` at the innermost cell, nx=24, as a
fraction of `seed`:

| N (=M) | 4 | 5 | 8 | 9 | 16 | 17 | 32 |
|---|---|---|---|---|---|---|---|
| `r/seed` | 0.0102 | **2.0096** | 0.0093 | **2.0091** | 0.0090 | **2.0088** | 0.0088 |

⟹ **The prediction holds, and therefore the hypothesis is refuted.** The
alternation is exact, and the residual off it is `0.88 %` and *shrinking*
with N (1.02 % → 0.88 %). The pole is where the recurrence does **best**,
not worst. Its ~0.9 % is the smallest term in sight, not the limiter.

(Also `[M]`, and consistent: `|direct⁺ − direct⁻|` at the innermost CELL
is `6.17e-2` and **N-independent** to 4 digits over N=4…64 — the pole
continuation is imposed at the r=0 FACE, so the two legs' cell-CENTRE
values differ by the O(Δr) march across the first cell. Spatial, hence
flat in N.)

## F9. ⛔ The α-closure contract is UNENFORCED in the canonical run

`α_{M+1/2} = 0` is not an axiom of the recursion (which is strictly
one-sided, `α_{m+1/2} = α_{m−1/2} − w_m μ_m` from `α_{1/2} = 0`) — it is a
**consequence of the measure's antisymmetry**, i.e. a genuine contract on
every admitted quadrature. Its enforcement:

* **sphere** — `orpheus/geometry/reduced_operator.py:784`, a **bare
  `assert`**: `assert abs(alpha[N]) < 1e-12`. Not a type-narrowing
  assert; a numerical contract check.
* **cylinder** — `:885-893`, **no check at all**.

`[M]` The project's canonical runner is `python -O -m pytest`
(`.claude/rules/vv-testing.md:12`), and `-O` sets `__debug__ = False` and
strips every `assert`. Demonstrated on the verbatim recursion: a measure
with `alpha[N] = +0.2000` is **REFUSED** under plain `python` and
**ACCEPTED** under `python -O`. The one production guard on this contract
is inert exactly where the project runs its tests.

This is independent of everything else here and is small to fix.

### F9 — ✅ FIXED, 2026-08-12 (uncommitted at time of writing)

The fix was **not** "add a check to the cylinder arm", because that would
have been adding a guard to a *duplicate*. `[M]` the recursion had **three
spellings** — sphere factory, cylinder factory, and
`derivations/discrete/sn/angular_differencing.py::alpha_dome` — which is
exactly why the contract could be enforced on one arm only. Cardinal
Rule 2 applies before the bug fix does.

What landed:

* **`orpheus/geometry/reduced_operator.py::alpha_dome(cosines, weights)`**
  — ONE body, now public (`orpheus.geometry.alpha_dome`). Both curvilinear
  factories call it; the derivations name **delegates** to it.
* **`_assert_alpha_dome_closes(alpha, *, coord, level=None)`** — a real
  `raise ValueError`, per level on the cylinder so the failure is
  locatable. Shaped after its sibling
  `pole_angular_closure._assert_tau_within_unit_interval`.
* ⭐ **The recursion and the contract are deliberately SEPARATE
  functions.** The derivations P0/P4 predicate ladder exists to
  characterise measures whose dome does *not* close; a guard welded into
  the shared recursion would make that analysis unspellable. Two gates pin
  that split as load-bearing rather than stylistic.
* `_ALPHA_CLOSURE_ATOL = 1e-12`, justified in a comment against measured
  data: worst shipped residual `2.78e-16` (`folded_product(4,32)`), GL
  `5.6e-17`…`8.2e-17` over N=4…128 — ≈1 ULP of the dome peak, not
  drifting with N, so 1e-12 clears every shipped rule by ~4 orders.

**Verification.** `[M]` α is **bit-identical** to the retired bodies on
every shipped rule (`np.array_equal`: sphere GL 4/8/16/32/64/128; folded
(4,8), (4,16), (4,32), (2,6) — all levels; all three spellings agree).
The guard REFUSES `alpha[M+1/2] = 0.2` on both arms **under `python -O`**,
names the offending level on the cylinder, and admits `5e-13` (the
positive control — `vv-principles` #11). 18 new gates in
`tests/geometry/test_reduced_operator.py::TestAlphaDomeClosureContract`
(40 → 58 passed in that module); the module's 3 failed / 4 errors are
**pre-existing at HEAD** (verified in a detached worktree at `adb73fd5`,
byte-identical failure set) and belong to the #51 ledger.
`tests/geometry + tests/sn/primitives + tests/sn/mesh`: **1216 passed**,
only those pre-existing reds.

---

## F10. ⛔⛔ THE DECIDING EXPERIMENT: `D` IS **NOT** A PROXY FOR ACCURACY

F7 established `D` as a **discriminator** (it sees τ; it ranks garbage
2.6–45× above the shipped partition). It did **not** establish `D` as a
**selection criterion** — that lower `D` means a more accurate solution.
Probe I tests exactly that, against a **structurally independent**
reference: the anisotropic cylindrical MMS
(`build_cylindrical_anisotropic_mms_case`), whose reference is **analytic**,
not another ORPHEUS solver — so it is not the L49 reference-limited trap.

`[M]` nx = 80, `max_inner=500`, `inner_tol=1e-13`, all 12 solves
`converged=True`. Harness validated: the fixture is a **folded** rule (all
ξ > 0, 4 levels × M = n_φ/2, **all 4 carrying**), so `D` is the same object
probes A–G measured; and the shipped MMS numbers reproduce the tree's own
record — `1.1326e-3` at n_φ=16 is **1.67×** the gate docstring's pre-6.4
`6.782e-4`, matching the retirement note's stated "1.8–2× worse".

| n_φ | ranks agreeing | Pearson r on log(metric) |
|---|---|---|
| 8 | 2 / 4 | **+0.7515** |
| 16 | 0 / 4 | **+0.2608** |
| 32 | 0 / 4 | **+0.0630** |

⟹ **The correlation degrades monotonically with angular refinement, to
essentially zero.** At coarse angle both metrics are dominated by the same
gross error; as that resolves they decouple completely.

Per-variant MMS L2 (the independent metric), rank in brackets:

| τ variant | n_φ=8 | n_φ=16 | n_φ=32 |
|---|---|---|---|
| shipped (Q5.6.4) | **3.1503e-3** [1] | 1.1326e-3 [2] | 3.2611e-4 [3] |
| diamond τ ≡ ½ | 3.4055e-3 [2] | **3.7443e-4** [1] | **1.3279e-4** [1] |
| reversed | 7.2701e-3 [3] | 1.2548e-3 [3] | 2.9387e-4 [2] |
| shuffled | 1.7521e-2 [4] | 3.4577e-3 [4] | 8.5636e-4 [4] |

**⟹ `D` MAY NOT VOTE ON τ. Decided, negatively, in one 20-minute
experiment — before any machinery was built on it.**

### Reading the τ column honestly

* The **diamond beats shipped by 3.0× / 2.5× at n_φ = 16 / 32**. This is
  **not news**: it is the *already-ratified* honest cost of the Q5.6.4
  carve, recorded in `pole_angular_closure.py`'s own retirement comment
  ("~1.8–2× WORSE at n_φ = 16/32/64 … Principled ≠ more accurate; the
  scheme that satisfies P2/P3 wins over one with a smaller number on a
  single manufactured fixture"). Probe I reproduces it against a
  *different* alternative and at a different nx, which strengthens the
  record rather than disturbing it.
* ⚠ **The n_φ=32 shipped-vs-reversed inversion is 1.11×** (`3.26e-4` vs
  `2.94e-4`) and **no claim is built on it.** An 11 % gap on one
  manufactured fixture is not evidence that a permuted τ is a better
  scheme; it is evidence that at n_φ=32 this fixture no longer separates
  these two.

### What this does and does not do to `D`

`D` remains what it was measured to be: a **cheap, reference-free,
pointwise consistency residual** for a boundary condition production
computes and discards. What it is NOT is an accuracy indicator — and that
distinction is now a *measurement*, not a caution.

⭐ **This is the L49 failure mode caught prospectively.** Had step 2 of the
proposed plan run first, `D` would have "confirmed" the landed carve
(F7's 2.6–45× ranking) using an instrument that does not track accuracy,
and the confirmation would have entered the record as evidence. The
campaign's τ-blind trio and the reference-limited resolvent both got that
far. This one was stopped by asking what the instrument measures *before*
letting it decide anything.

⛔ **Standing consequence: the campaign still has NO reference-free
instrument that can rank τ.** §9bis.8c's situation is unchanged by this
session. Any future `D`-based τ argument must cite this section first.

---

## What this changes for the campaign

The memo's §9bis.8c recorded that **all four** surviving accuracy
rankings were *inverted* against closure fidelity — every instrument made
the landed convention look worse, which is why Q5.6.4 closed on
literature ratification (M&M's own cylinder appendix) rather than on a
measurement. `D` is the first instrument that is **not** inverted, and it
is the first that needs **no reference at all** (contrast the resolvent
cross-check, reference-limited at ≈3e-2 — [[lessons-L49]]).

⛔ **REFUTED as a campaign advance — see F10 (2026-08-12).** The paragraph
above is true as written and irrelevant as an argument: `D` being
"not inverted" is worthless once measured against an independent
reference, because `[M]` `D` does not track accuracy at all
(Pearson r on log = **+0.06** at n_φ=32, 0/4 ranks agreeing). An
instrument that ranks the landed carve first while being uncorrelated
with accuracy is not evidence for the carve. **§9bis.8c stands
unchanged: there is still no reference-free instrument that can rank τ.**
What `D` genuinely is: a consistency residual for a discarded boundary
condition, worth gating on its own terms and forbidden from voting on a
scheme.

## Open / not done

* `D` is not computed, named, or gated anywhere in production. Whether it
  should be a first-class named quantity is the open design question.
* The vacuum-surface angular boundary layer (F6) is a real, unquantified
  accuracy limiter at the outermost cell. Plausibly related to the #229
  azimuthal-floor story; **not** established — no measurement links them.
* Only the σ_y-folded **cylinder** was measured. The sphere arm has the
  same two-endpoint structure and was not probed.
* `_carrying_levels`, `_tau_per_level` are private; a production
  realization needs a public surface.
