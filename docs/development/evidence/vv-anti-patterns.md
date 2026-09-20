# V&V anti-patterns — evidence

The thirty-six anti-patterns of the `vv-principles` skill — each one's NEVER
sentence, mechanism, CHECK sub-clauses and founding measurement — moved VERBATIM
from `.claude/skills/vv-principles/SKILL.md` at commit `207d2b07` (2026-09-20).
The review-time instrument is [the core](../skills/vv-principles.md); every
`[case]` link in its anti-pattern list lands on a heading of this page. Text here
is history and keeps its original glyphs and its original list numbering.

## AP1 l4 alone

1. **NEVER** claim verification on the basis of L4 agreement alone —
   **instead** require an L0–L2 evidence chain pointing at a
   structurally-independent reference. Two ORPHEUS solvers agreeing is
   _cross-implementation agreement_, NOT _correctness evidence_.

## AP2 allclose solver

2. **NEVER** assert `np.allclose` against another solver in this
   codebase — **instead** match the claim to a reference at the right
   level (analytical for L1, MMS for spatial convergence, MC tally for
   L4 cross-check only after MC itself is verified).

## AP3 one group

3. **NEVER** accept a 1-group eigenvalue test as evidence of solver
   correctness — **instead** demand ≥2 groups. k = νΣ_f/Σ_a is
   flux-shape independent; 1G is degenerate. (An instance of the
   Mode-12 invariance-group lens — test-design table below.)

## AP4 homogeneous only

4. **NEVER** accept a homogeneous-only verification — **instead**
   demand at least one heterogeneous, mesh-refined, multi-group case.
   Flat flux nulls every redistribution and weight-cancellation term.

## AP5 rate only

5. **NEVER** read "convergence rate is correct" as "result is correct"
   — **instead** verify the converged-to value. O(h²) to the wrong
   limit is still O(h²).

## AP6 untraced reference

6. **NEVER** trust a reference that has not been traced back to a
   structurally-independent analytical or symbolic ground —
   **instead** treat it as **reference contamination** until the
   trace is shown. The most seductive failure mode: MC vs MC,
   CP vs unverified MC, method-of-images converged to the wrong BC.

## AP7 derivations agree

7. **NEVER** treat "two derivations agree" as proof — **instead**
   check whether they are _structurally_ independent. ERR-032 (two
   antiderivatives both using `∫E_2 = 1 − E_3` instead of `½ − E_3`)
   agreed at 1e-39 because they shared the upstream identity, not
   because either was right.

## AP8 particle balance

8. **NEVER** accept "particle balance holds" as L0 evidence —
   **instead** require per-ordinate flat-flux residual. Telescoping
   sums hold by construction even with wrong per-ordinate balance.
   (A Mode-12 instance: the balance functional annihilates per-ordinate
   errors that cancel in the sum.)

## AP9 validation verification

9. **NEVER** conflate validation with verification — **instead** state
   which screw is being turned. If the equation itself is wrong,
   verification can pass cleanly; only L3 catches it.

## AP10 reasonable numbers

10. **NEVER** accept "it produces reasonable numbers" — **instead**
    enumerate every term, isolate it, and verify sign AND magnitude.
    Sign-flipped small terms look reasonable.

## AP11 negative legs

11. **NEVER** test a contract-validation method (`assert_X`, `check_X`,
    `verify_X`) ONLY against a deliberately-broken instance —
    **instead** require AT LEAST one positive test (correct instance,
    MUST NOT raise) AND AT LEAST one negative test (broken instance,
    MUST raise). Negative-only testing validates the *raising
    behaviour* but NOT the *invariant claim* — the test cannot tell
    you the method's claim is correct, only that the method raises
    when told to. ERR-051: `assert_galerkin_idempotency` asserted
    `Π R = I` instead of `Π R = 4π · I` under the no-prefactor SH
    convention; the bug hid for an entire merge cycle because the
    sole test fed it a deliberately-non-orthogonal Y so the wrong
    invariant produced the expected failure. The test was
    self-referential: the broken Y was constructed precisely to make
    the wrong assertion succeed at raising. The
    structural-independence requirement ([L11](lessons.md#l11-structural-independence)) applies to ALL test
    design, not just numerical cross-checks.

## AP12 neutral proxy

12. **NEVER** credit a "behavior-neutral field-zeroing / relabel /
    no-op retype" claim on the basis of a fast proxy (snapshots didn't
    move, no guard errors raised, type-check passes) — **instead**
    re-prove neutrality for EVERY consumer with a direct old-vs-new
    VALUE comparison (`np.array_equal` / `assert_array_almost_equal_nulp`
    on the consumer's actual output). A neutrality claim holds only for
    the ONE fission/emission/operator contract it was proven against;
    proxies are blind to a per-consumer divergence whose precondition
    the zeroing itself breaks. ERR-063: "zeroing χ on non-fissile
    regions is inert" was TRUE for the SN/`compute_macro_xs` contract
    (χ gated by the SAME region's νΣf) and FALSE for `solve_peierls_mg`
    (source-region νΣf weighted by SINK-region χ) — the same-χ snapshots
    masked it because χ_i ≡ χ_j everywhere until the zeroing broke that
    equality. A green snapshot pinned the masked regime, not the claim.

## AP13 finite sample

13. **NEVER** accept a finite "representative sample" of a group /
    parameter family / operator set as a check for the WHOLE thing —
    **instead** compute the object the sample actually generates and
    compare it to the claimed one. A generator-set check is sound
    exactly when the listed elements GENERATE the claimed group (then
    closure under each generator implies closure under every product);
    it is a false certificate the moment they generate a proper
    subgroup, and the failure is *designed-green* in the Mode-12 sense —
    no tolerance, refinement, or regime change can expose it. ERR-072:
    `SubgroupOfO3.SO2.is_invariant` sampled `{0°, 90°, 180°, 270°}`
    about z — four rotations that generate `C_4`, not `SO(2)` — so
    every product quadrature with `n_phi ≡ 0 (mod 4)` certified as
    `SO(2)`-invariant while being invariant only under `C_8`. Two
    review tells travel with this pattern: (a) a docstring that
    **pre-authorises the gap** ("necessary but not sufficient in
    general, but sufficient by construction for the rules we ship") —
    a named risk reads as an assessed risk, so CHECK the enumerated
    "rules we ship" against the sample, one by one; (b) for a
    CONTINUOUS group, the honest discrete predicate is usually a
    *different* question entirely (a finite node set is
    `SO(2)`-invariant iff every node is ON the axis), which means the
    tag being asserted describes the CONTINUUM object being
    discretised, not the discrete one — two different claims must
    never share one predicate name.
    ⭐ **The same defect wears a REFINEMENT LADDER as its disguise**, and there
    it looks like the opposite of a sample. `n = 8, 16, 32, 64` reads as "four
    orders, an 8× range, every order" — but it is a **single congruence class**,
    and any failure mode gated on a divisibility property of `n` is either
    present in all four or absent from all four. (`[M]` 2026-08-11, SN Q5.6.4:
    "exactly 1 of `M+1` edges has no real solution, **every quadrature order**"
    was measured at `n_φ = 8/16/32/64`, i.e. `M = 4/8/16/32`, all EVEN. The
    failing edge is the one at `ω = π/2`, which is an edge only when `M` is even;
    at odd `M` — `n_φ = 6, 10, 14, 18, 26, 34, 66`, every one a legal shipped
    rule — **0 of `M+1`** fail, and the conclusion "ill-posed at every order"
    that a design was built on is a parity artefact.) Review rule: a refinement
    ladder must break the arithmetic pattern of its own step — include one
    non-power-of-two, one odd, one prime — before any universal is claimed over
    it.
    ⭐ **And the third disguise, which needs no arithmetic pattern at all: the
    sample drawn from the STALE CLAIM'S OWN NEIGHBOURHOOD.** When you spot-check
    an inherited formula, the orders/parameters that come to mind are the ones
    the claim itself made salient — and those are exactly where a formula that
    was *once* right still agrees. The bias is in the sampling, not in the
    arithmetic, so "I checked three values and it held" carries no information
    about the values you did not think of. ⟹ **check a claim at a point its own
    text does NOT name**, and prefer the point where the superseding change was
    largest. (`[M]` 2026-08-14, quadrature Q6-A: `degree_of_exactness` for
    level-symmetric `S_N` was documented as `N-1`. The realized degrees are
    `S2→3, S4→5, S6→7, S8→9, S10→11, S12→11, S14→15, S16→15, S18→17` — the
    formula is right at **S2/S12/S16/S18** and under-claims by 2 everywhere
    else. `S_12` is the *retired frontier* the stale prose kept naming, so the
    order a checker most naturally reaches for is one of the four that confirm
    it. Same file, same day, same cause: the frontier itself still read `S_12`
    while the tree served `S_18`.) This is the dual of the ladder rule above —
    that one says *break the pattern of your own steps*; this one says *do not
    let the claim choose your steps for you*.
    ⭐ **And the FOURTH disguise, where the sample looks continuous and is
    not: a TOLERANCE sweep on an iterative solver.** A tolerance is a knob
    with a *discrete* effect — it acts ONLY through the iteration count it
    induces — so four decades of it can land in ONE equivalence class, and
    the honest reading *"this tolerance does not move the error"* is really
    *"none of my values changed the iteration count."* The failure is
    flattering in the specific way that matters when a BAND is being
    derived: it attributes the whole error floor to the one tolerance that
    happened to bite, and the band then omits a term that is real and
    binding at the shipped configuration. ⟹ **report the iteration count
    beside every row of a tolerance sweep; two rows with the same count are
    one measurement, not two**, and a sweep whose count never moves has not
    tested that tolerance at all. (`[M]` 2026-09-06, ORPHEUS #448: a
    pre-carve memo derived the finalize band from a sweep reading *"four
    decades of `flux_tol` move it not at all; the empirical driver is
    `inner_tol` alone"* — `1e-6/1e-7/1e-8/1e-9` all gave `n_outer = 10`. At
    `flux_tol = 1e-11` the count goes to 12 and the reconstruction's
    deviation from the converged iterate falls `3.43e-11 → 6.96e-13`, **49×**
    — so the outer term is not merely bounded but *dominant* at the gate's
    own `inner_tol = 1e-11`, which is the opposite of what the sweep was
    read to say. The band itself was correct; its stated mechanism was not,
    and a later session tightening `inner_tol` alone would have found the
    floor immovable and gone looking for a bug.)

## AP14 not bijection

14. **NEVER** read "every element found a matching partner" as "the map
    is a bijection" — **instead** assert the structure the docstring
    names. A nearest-neighbour / lookup loop that finds, for each `i`,
    *some* `j` within tolerance computes a **relation**, not a
    permutation; many-to-one maps satisfy every assertion in the body.
    ERR-073: `_orbit_closure` documented "find a permutation π such
    that `M(nodes)_i = nodes_{π(i)}`" but never checked injectivity, so
    duplicating one node of an `O_h`-invariant rule (bit-identical
    duplicate — no tolerance games) produced a measure with
    `M#µ ≠ µ` (mass `1.047` vs `0.524` at the same point) that
    certified invariant, with the match map non-injective for 48 of 48
    operators. Generalisation: whenever a docstring names a STRUCTURE
    (permutation, bijection, isomorphism, partition, basis) that the
    body only *implies*, either assert it or weaken the docstring —
    and prefer **returning the structure** to returning a `bool` about
    it, because a returned permutation makes its own bijectivity
    assertable while a `bool` makes it unfalsifiable.

## AP15 order predicate

15. **NEVER** ship a module that exposes BOTH an order relation
    (`contains` / `is_subgroup_of` / `refines` / `⊆`) AND a predicate
    that must respect it (`is_invariant` / `satisfies` / `admits`)
    without a test of the **compatibility law** — **instead** gate
    `A ⊆ B  ∧  P(B, x)  ⟹  P(A, x)` over every (edge × fixture) pair.
    The law is one loop, needs no external reference, and cross-checks
    the two halves against each other: neither half can be wrong alone
    without the law reddening. Measured on `numerics/symmetry.py`
    (2026-08-02): 68 violations over 11 measures × 19 groups, isolating
    a false lattice edge (`D_nh ⊆ O(2)`, itself pinned by a committed
    test), a sampled-group checker (ERR-072), and a realisation
    mismatch (`Z2 ⊆ SO(3)` asserted while `Z2` is realised as an
    improper reflection) — three independent defect classes surfaced
    by a single invariant that no per-relation and no per-predicate
    test could see.

## AP16 tighter invariant

16. **NEVER** assert a property TIGHTER than the type's own
    construction invariant — **instead** split it into two gates: one
    on the invariant the type actually promises (plus the threshold at
    which it rejects), one on the constructors' far better realised
    quality. A gate demanding `‖QᵀQ − I‖ ≤ 1e-14` of an *arbitrary*
    element of a type whose `__post_init__` admits `1e-12` is asserting
    something the type does not guarantee: the `1e-13` shear is a
    **legal value**, so the gate is a latent false red that a future
    legitimate input will trip, and it silently mis-states the
    contract to every reader. (2026-08-03, `geometry/transformation.py`:
    the constructors are exact — `signed_permutation` measures `0`,
    Householder a few ULP — so the two claims differ by four orders and
    conflating them buys nothing.) The general shape: **a gate on a
    type's invariant must quote the type's own threshold**, and any
    stronger claim belongs to the *producer* that achieves it, not to
    the type that permits it. Corollary for review: when a tolerance in
    a test is tighter than the one in the production guard it is
    testing, one of the two numbers is wrong — find out which before
    relaxing either.

## AP17 positive control

17. **NEVER** run a mutation battery without a **positive control** —
    **instead** include one deliberate mutation that MUST redden many
    gates, and treat an all-blind verdict as *the harness is broken*
    until that control proves otherwise. The harness lies before the
    code does, and it lies in the SAFE-LOOKING direction: a parser that
    fails open reports "0 caught" — which reads as "write more tests"
    rather than "your instrument is dead". (2026-08-03: a battery
    reported **32/32 BLIND** while its own captured summaries plainly
    read `23 failed` / `63 failed` — it scanned for `FAILED` lines that
    `-q --tb=no` never emits, and ANSI codes broke the match. A control
    mutation making `reflection` return `+I` cannot leave 42 gates
    green; that contradiction is what exposed it. Cost one run.) Same
    family as the earlier subprocess-mutation failure (monkeypatching
    the PARENT while pytest re-imports a clean module in a CHILD reads
    GREEN for every mutation): in both, the *evidence pipeline* failed
    while the code under test was fine. **Verify the instrument on a
    known-positive before trusting any negative it reports.**
    ⭐ **And the AIM trap, which produces a false RED verdict about an
    arm: a mutation that reddens MANY rows while the row it was designed
    for stays GREEN hit a different object than the one it was aimed at.**
    Read an arm's red SET against its TARGET row, never its count — collateral
    reds certify the instrument, not the aim. (`[M]` 2026-09-14, the C3b-2
    battery: the "no refusal" arm patched `solve_sn`, which is NOT on the
    multiplying-source refusal path (the driver reads
    `power_iteration(SNSolver(...))`); 25 rows reddened, the refusal leg
    stayed green. Re-armed on the NAMED predicate `_admits_multiplying_source`
    — a name the concept had lacked, which is why the arm had nothing to
    mutate — it reads exactly 1 red = the refusal leg.)
    ⭐ **And the granularity trap, which produces a false GREEN verdict
    about coverage: mutating a MULTI-ARM guard as a unit certifies only
    its first reachable arm.** A guard written as a sequence of
    early-return checks (`if x is None: reject` / `if wrong_kind: reject`
    / `if too_small: reject`) reads as one thing and is N independent
    claims. Break the whole guard and *something* reddens — the arm whose
    input the suite happens to supply — and the run reports "the guard is
    gated". The arms nothing exercises are then indistinguishable from
    the ones that bite, and they are precisely the arms most likely to be
    inert, because no shipped input reaches them.
    ⟹ **mutate each arm separately and record which gate each one
    reddens.** The verdict is a table, not a boolean. An arm that reddens
    *nothing* is a guard with no witness: either construct one, or say in
    its docstring that it is unfalsifiable and why.
    ⭐ **And the way an arm SILENTLY loses its witness later: someone adds a
    correct guard EARLIER on the same path.** The clauses above are about a
    guard that never had a witness. This is about one that had exactly one
    and stopped — and it is invisible, because the new guard makes the old
    gate go RED (a wrong-message mismatch), which reads as "update the
    expected string" rather than "that gate was the only pin on a different
    guard". Relax the match and the displaced guard is now unfalsifiable,
    with its docstring still advertising the protection.
    ⟹ **when a new guard preempts an existing one, ask what was pinning the
    OLD one — and if the answer is "the test you are about to re-point",
    write it a direct witness in the SAME commit.** Prefer pinning the
    displaced guard at *its own predicate*: a guard keyed on a value signal
    is reachable by calling it directly, which no earlier guard can preempt,
    and is better evidence than an end-to-end path that merely happened to
    reach it.
    > `[M]` 2026-08-26, un-weld P1 item 9. A new moment-mass guard (a
    > curvilinear multi-moment metric cannot be spelled) fired one frame
    > before the #158 LD scan-closure guard, so
    > `test_ld_curvilinear_scan_rejected` failed on
    > `match="slab/Cartesian"`. `[M]` `grep "slab/Cartesian" tests/` returned
    > **exactly one hit** — that test was the scan guard's *only* witness, and
    > the one-line fix (re-point the regex) would have orphaned it while
    > leaving the whole suite green. The scan guard keys on a VALUE signal
    > (`dA_w`/`c_out` non-neutral, not a chart tag), so a direct
    > positive+negative witness with no mesh at all was both possible and
    > stronger; both landed in the same commit.
    (`[M]` 2026-08-14, the quadrature selector's V conjunct: three arms —
    no claim / wrong reference / degree short. Mutated separately, two
    reddened exactly one gate each and the **`claim is None` arm left all
    54 gates green**, because no registered rule and no test-local spec
    had `exactness=None`. A whole-guard mutation would have reddened on
    the reference arm and certified all three. Fixed by constructing a
    claimless spec; the real future occupant — `folded_product` — cannot
    be registered yet for unrelated structural reasons, which is exactly
    why the arm had no witness.)
    ⭐ **A third pipeline failure, and the cheapest to introduce: a
    mutation that makes production RAISE kills COLLECTION, and pytest
    then reports `FAILED = 0`.** Any production call in a `parametrize`
    ARGUMENT LIST (or at test-module scope) runs at import, i.e. AFTER
    the mutation plugin installs — so `Interrupted: 1 error during
    collection` / `rc=2` is what the battery records, and a summary
    scanner counting `^FAILED` reads it as "0 caught", again in the
    safe-looking direction. (`[M]` 2026-08-11, the SN angular-closure
    seam: 6 of 13 mutations — *including the positive control* — came
    back as clean zeros in 2 s each, because the module's own
    orientation gate built its cases in the `parametrize` list and the
    P3 guard raises one frame downstream of the mutated producer.) Two
    fixes, both one line: **never call production in a `parametrize`
    argument list** — parametrize by a LABEL and build inside the body;
    and run the battery with `--continue-on-collection-errors`, counting
    `^ERROR` separately from `^FAILED` so a collection kill can never be
    read as a green.
    ⭐ **And the shape the granularity trap takes after a Pattern-2 hoist: a
    guard moved to ONE shared home has as many arms as it has CALL SITES.**
    Single-sourcing the guard *body* is right and does not single-source the
    *wiring* — each site passes its own operands and its own owner label, and
    those expressions can differ from one another. So the elegance move that
    removes the duplication **creates** the blind spot, and a mutation of the
    shared body (which reddens *something*) certifies all of them. ⟹
    **enumerate the call sites, then enumerate the distinct expressions they
    pass, and mutate per SITE.** The verdict is a table with one row per call
    site. A site that reddens nothing has no witness, and the site most
    likely to be miswired is the one whose operand expression differs from
    its siblings'.
    `[M]` 2026-08-21, ORPHEUS CS4a-R `assert_energy_extent_conforms`: one
    body, four call sites (`fission.py:201`, `multiplication_operator.py`,
    `isotropic_scattering.py` ×2). Disabled per site over 655 rows — **F 1
    red, C 0, IsoS+IsoN2N 0** — and `grep -rn "energy extent" tests/`
    returned exactly **one** assertion tree-wide. The C site was the only one
    keying on `self.coefficient.values.shape[0]` rather than
    `self.mat_xs.ng`; it was also the one with no witness. The module's
    docstring carefully wrote the *inertness* denominator ("live on 4 of 13
    production bindings") and never the *witness* denominator. Repaired the
    same day: one per-site wrong-ng row each, message asserting the
    constructing operator's name.
    ⭐ **And the verdict that must be read by IDENTITY rather than by SIZE:
    when a mutation's red set is EXACTLY the tests that NAME the mutated
    symbol, the symbol has no consumer and its pins are a MIRROR, not a
    gate.** Every clause above reads a red count; this reads the red set's
    *composition*, and it is the one case where a healthy-looking count is
    the finding. Note #18 does **not** catch it — asking *"by what mechanism
    does THIS gate see THIS property?"* returns a perfectly good answer
    ("it reads the field directly"), because the defect is not that the pin
    is blind but that **nothing downstream of the value exists**. Such a pin
    asserts only that a producer wrote what a producer wrote; it is green
    forever, it names a real-sounding contract, and a coverage audit counts
    it. ⟹ after any mutation, diff the red set against
    `grep -rln "<symbol>" tests/`. Equality ⟹ retire the symbol (or wire the
    consumer the contract implies); a red OUTSIDE the naming set is what
    proves a consumer exists. ⚠ Two mechanics this depends on: the mutation
    must be patched at **every rebinding site** (a name re-exported through a
    package `__init__` keeps the ORIGINAL — that under-reddened by 50 % on the
    first attempt below), and it needs a call counter, or a small red set is
    ambiguous between *no consumer* and *no bite* (#17's positive control).
    `[M]` 2026-08-26, ORPHEUS `ReducedStreamingOperator`: FLIPPING both
    `requires_upstream_angular_state` and `angular_marching_axis` on **997**
    constructed operators, over `tests/{sn/sweep,geometry,sn/primitives,
    transport}` at `-m "not slow"` — **2585 → 6 failed / 2579 passed** — and
    all six reds are the six assertions that name the fields (3 in
    `test_reduced_operator.py::TestProperties`, 3 in
    `test_snmesh_consumes_reduced.py`). Zero production readers; the chart
    predicate the fields claim to carry is live under two OTHER spellings
    (`upstream_state.angular_upstream is None`, `SNMesh.is_cartesian`).
    A grep had said the same, weakly; the battery also rules out dynamic and
    inherited readers.
    ⭐ **And the POSITIVE CONTROL can itself be Mode-12 blind: a mutation
    inside the SUT's own invariance group is a null control.** #17 says
    include a mutation that MUST redden many gates. The one that comes to
    mind is the *simplest* corruption of the object — negate it, flip it —
    and the simplest corruptions are exactly the symmetries a well-built
    fixture is invariant under. `[M]` 2026-09-02, #429 tracker 2.3: the
    control for a chart battery negated μ in the Archimedes chart
    (`out[:, a] = -mu`); Gauss–Legendre nodes are symmetric about 0 and
    the product rule carries σ_h, so the node SET was unchanged and only
    the **9** order-sensitive pins moved, while two ordinary arms reddened
    **58** and **57**. Read naïvely, the battery says the suite is blind to
    its control. ⟹ **choose the control OUTSIDE the fixture's symmetry
    group** — break a codomain, swap an axis, scale a weight — and, if an
    arm you meant as the control reddens fewer gates than an ordinary arm,
    name the ordinary arm as the effective control rather than lowering
    the bar. ⭐ **And the NEGATIVE control of an invariance test has the
    same trap one level up: the SUT's stabiliser can be strictly larger
    than the group it is DECLARED under, so "an element outside H" can
    still be a null control.** `[M]` 2026-09-02, #429 tracker 3.1
    (archivist): the quotient map π_a of S²/SO(2)_a is bit-exactly
    invariant under the mirrors σ_b, b ≠ a — O(2)_a and SO(2)_a induce the
    same orbit partition, so a quotient map fixes the PARTITION and the
    partition does not fix the group. A negative leg spelled as "a mirror,
    which is not in SO(2)" reads as a control and moves nothing. ⟹ before
    choosing the outsider, measure the SUT's actual stabiliser (or reason
    it out: what else preserves the partition?), and pick the outsider
    outside THAT — here a rotation about another axis. Same session, same battery, the harness's other lie: the
    driver's inline `$(...)` capture returned EMPTY summaries for the four
    widest arms; read arm output from a FILE.

    ⭐ **And the FILTER form of the positive-control rule: a two-stage census
    (a name-net THEN a literal pattern) needs a control per STAGE — a synthetic
    fixture you author passes stage 1 by its own NAME and certifies only
    stage 2.** When a census first selects candidate functions by a token net
    and then scans their bodies, a control whose function is named with one of
    the net's tokens clears stage 1 for free, so every arm validates the
    literal patterns and none validates the net. ⟹ name the control fixture
    with the spelling you are LEAST sure the net catches, or add a member that
    is only reachable through a token you did not think of.

    > `[M]` 2026-09-03, ORPHEUS #428 census. `test_n2n_multiplicity_census`'s
    > net `("n2n", "sig2", "sig_2n", "_2n")` misses `sig_2` — the spelling
    > `derivations/` uses — so its claim *"a thirteenth literal home is
    > unspellable"* was false: widening the net by that ONE token finds 2
    > literals (`derivations/common/eigenvalue.py:61, :290`, correct by design
    > there — the reference tree must not read the SUT's constant — but
    > escaping by a filter gap, not a named exclusion). Its control names its
    > synthetic function `n2n_source_assembly`, so it never exercised the net.

### 2026-09-20 the prescribed repair

**Surprise.** A mutation battery on `tools/docs/generate_harness.py --check` (arm 4b, `scratch/_harness_eval/review/qa_generator_mutation.md`): deleting the role-block markers from an `AGENT.md` made `--check` red with a DRIFT message prescribing a regeneration; applying that regeneration inserted a SECOND role block, and `--check` then read green. The detector's own remedy converged on a wrong state it blessed.

**Clause.** #17 (i): a gate whose message names a fix owes an arm that APPLIES the fix and re-runs the gate. Rider on #17 (a), same battery: the budget key was mutated on a rule entry and caught, and on an agent entry and ignored — a config option is mutated once per CONSUMER kind.

## AP18 breaks structure

18. **NEVER credit a mutation's reds as coverage of a property when the
    mutation also breaks a STRUCTURAL law the object obeys** (linearity,
    symmetry, positivity, conservation, a shape/type contract) —
    **instead** mutate INSIDE the object's algebraic class, so the only
    thing that changed is the property under test. This is the exact dual
    of #17 and it fails in the *opposite*, more flattering direction:
    #17's broken harness reports a false "0 caught" (which reads as
    *write more tests*); an over-powered mutation reports a false "richly
    caught" (which reads as *nothing to do here*) — and a coverage audit
    closes on it. The tell is a red count wildly out of scale with the
    property's reach, and reds concentrated in gates that have no
    plausible view of the property (end-to-end eigenvalue / convergence
    gates reddening for a boundary-slot bookkeeping claim). Ask of every
    red: *by what mechanism does THIS gate see THIS property?* — if the
    honest answer is "it doesn't; it sees the law I broke", the gate is
    not a catcher. (2026-08-03, the SN `(L+C)` matvec's tangential trace
    slots: writing a CONSTANT sentinel into those rows made the operator
    **affine**, and 60 gates reddened — every one of them a Krylov/SI
    solve that diverges when its operator stops being linear. The
    realistic bug is linear (`out[tan] = ±ψ[tan]` — what you get by
    initialising the output block from the input, or by the "not inflow
    == outflow" trap), and re-run in that class it reddened **exactly
    1 of 5076** tests, with 94 148 rows mutated. The affine verdict
    over-stated coverage by 60×.) Corollary — the honest way to size a
    property's reachable audience BEFORE mutating: check whether the
    measured quantity's **metric/functional annihilates it** (the Mode-12
    stabiliser check). Here the trace metric `G = |Ω·n|·w_n` is *exactly*
    zero on tangential rows — a `1e6` perturbation there moves
    `⟨·,·⟩_G` by `0.0`, bit-identical — and the rows are decoupled from
    the bulk, so no solver-level, norm-level or reciprocity gate can EVER
    be a catcher and a direct array assertion is the only instrument that
    can exist. Knowing that first would have flagged the 60 reds as
    impossible on sight.

## AP19 positive reading

19. **NEVER** cite a gate's POSITIVE reading as evidence that the gate is
    *loaded* on the structure it is credited with (metric-loaded,
    weight-loaded, transpose-loaded) — **instead** cite the reading under
    the DELIBERATELY-WRONG structure. A tiny residual is exactly what a
    *blind* gate produces too, so the positive leg cannot discriminate
    loaded from blind; only the negative leg can. This is #11's
    positive+negative pairing applied to the *structure* a gate rides on
    rather than to a contract-validation method, and it is the operational
    form of the Mode-12 stabiliser question: "is the thing I claim this
    gate is sensitive to actually outside the measured functional's
    stabiliser?" The tell in review is a comment of the shape "X is the
    metric-loaded partner: `[M]` residual 1.8e-15" — the number quoted is
    the one measurement that carries **zero** information about loading.
    (2026-08-06, the SN affine-boundary P5 reciprocity rows: two new
    `_BUILDERS` cases were added with a prose argument that a partner face
    is *mandatory* because the zero morphism is metric-blind (`0ᵀ = 0`
    under every metric) — and then neither case was added to the
    committed wrong-metric control, which stayed `["slab", "sphere"]`.
    Measured by dropping `|Ω·n|` from the trace metric while `A.H` stays
    built for the true one: `slab_declared_prescribed_2g` reads
    `1.98e-16` true / **`2.410e-01` wrong**, `..._white_2g` reads
    `1.68e-16` / **`1.351e-01`** — against the already-listed `sphere`'s
    `1.05e-16` / `1.213e-03`. The two ungated cases fire **100–200×
    harder** than the one that IS gated, and both clear the control's own
    anti-dud precondition (`|Ω·n|` spread `0.5212 > 0.1`). The argument
    for needing the partner was right; the evidence offered for it was the
    wrong measurement, and the control that would have supplied the right
    one was one list entry away.) Review rule: for every "this fixture is
    loaded on S" claim, demand the S-broken reading, and if a negative
    control for S already exists in the module, the new fixture belongs in
    its parametrize list — a partner argued for in prose and ungated in
    code is an unverified coverage claim.

## AP20 rows cases

20. **NEVER** count the ROWS a new case multiplies into as new coverage —
    **instead** count the CASES, and for each row the case reaches, name
    the line in that row's BODY that reads the thing the case varies. A
    shared `parametrize` dict (`_BUILDERS`, `_CASES`, a fixture registry)
    is consumed by every test in the module, so adding one case adds one
    row per consumer — and a row whose body never touches the varied
    thing is *structurally* incapable of reddening for it, in the exact
    Mode-12 sense (not sub-floor, not under-tested: annihilated). The
    inflation is invisible in a diff, which shows `+3 lines` and a test
    count that jumps by 6, and it lands in the closeout as "+6 rows" of
    coverage. (2026-08-06, SN affine-boundary P5: 3 new cases → 6 new
    rows, of which **3 are provable non-catchers**. Two are
    `test_full_field_space_metric_matches_independent_reference[…]`,
    whose body builds the metric from `volumes` / `weights` /
    `omega_dot_n` and never reads `sn.bc` — measured **bit-identical** to
    the pre-existing `slab_2g` row, `g_inner = inner_product =
    -0.6830574021861343` for all three cases. The third is
    `test_full_loss_reciprocity_per_group_one_hot[…]`, whose one-hot
    composite zeroes the WHOLE trace block, which is exactly where `B`'s
    range and co-range live, so `⟨Bψ,φ_g⟩_G ≡ 0 ≡ ⟨ψ,B.Hφ_g⟩_G`
    identically. All three stayed green under **four** mutations —
    including a positive control that reddened 17 sibling reciprocity
    rows.) The rows are free and harmless; the *claim* is the defect.
    Write the closeout as "3 cases → 3 mutation-verified rows (+3 that
    ride along and are blind to the varied thing)", and if a row is blind
    for every case, say so once rather than re-counting it per case.

## AP21 line grep

21. **NEVER** audit a negated claim with a LINE-based grep — **instead**
    search a multi-line WINDOW (subject within ±2 lines of the negation),
    because prose wraps and the subject and its negation routinely land on
    different lines. This is the missing mechanic in the
    "grep the CONCEPT, not only the symbol" retirement-audit rule
    (`.claude/rules/coding-standards.md`): that rule tells you to widen
    the *vocabulary* you search for; this one tells you to widen the
    *window* you search in. A correction pass that greps
    `white.*not adjointable` finds every instance the formatter happened
    to keep on one line and silently reports the rest as clean.
    (2026-08-06, SN affine-boundary P5: a pass corrected the two sites
    where the subject and the negation shared a line —
    `SNBoundaryOperator.apply_transpose`'s docstring ("the white BC has no
    Euclidean transpose") and a test-module header ("white would drop
    it") — and missed the `SNBoundaryOperator` class docstring's
    re-emission-closure paragraph, where "white" and "(not adjointable)"
    sit on ADJACENT lines. Measured: `WhiteBoundary()` →
    `is_adjointable = True`, `AlbedoBoundary(0.7, IsotropicReturn(...))`
    → `True`, i.e. the surviving sentence is present-tense FALSE. It sits
    in the SAME class-docstring family, ~470 lines from a site the same
    pass corrected to say the opposite. ⭐ Sharpest detail: that paragraph
    was ITSELF a correction pass — its own closing sentence reads "it is
    the enumeration in prose that had to stop naming classes" — and the
    fix was applied to the *subjects* and not to the parenthetical
    *verdicts* beside them. A correction is not evidence its own paragraph
    is now clean.) **The aggravator is the reason
    this ranks as an anti-pattern rather than a grep tip:** a half-done
    correction leaves the stale claim and its correction coexisting in ONE
    FILE, which is strictly worse than either alone — a reader who lands
    on the stale one gets no signal it was superseded, and the file now
    contradicts itself, so *whichever* sentence a future contributor
    trusts, they can cite the file for it. Review rule: after any
    claim-correction pass, re-run the audit as a windowed search over the
    whole tree and reconcile every hit BY TENSE (past-tense history stays;
    a present-tense claim is a MUST-FIX) — and check the corrected file
    itself first, since it is the likeliest place for a survivor.

    ⚠ **The file-granularity twin: a finding that names ONE `file:line` will
    be repaired at that `file:line` only.** A review report's site list
    silently becomes the repair's denominator, so a finding about a CLAIM
    owes the claim's full site census (every spelling, every file), not the
    site where it was noticed — or the same sentence survives one module
    over, present-tense false, and the report that "closed" it reads as
    evidence it is gone.

    > `[M]` 2026-08-31, ORPHEUS CS4c step 0/3. F9 reported the multi-model
    > sharing sentence on `ScatteringOperator.isotropic_kernel`; commit
    > `0ed9dca2` corrected exactly that docstring (its message even quotes
    > the measured denominator), step 3 then retired `isotropic_kernel`
    > outright — and the identical claim in the iso energy module's header
    > (*"CP / MoC / diffusion feed raw scalar-flux arrays"*) survived,
    > present-tense-false: `[M]` cp/moc/mc reference 0 of the 10 roster
    > classes and 0 of the array verbs; 1 of the 3 named models fed it.
    > Landed as a rule 2026-09-04 (CS4c step 5), when the SAME sentence was
    > about to be re-published in the module's rewritten header.

## AP22 shared object

22. **NEVER** read "neither side calls the other" as "the two sides are
    independent" when the test CONSTRUCTS one input object and hands the SAME
    object to both — **instead** ask independence **per axis**: the
    *derivation* axis (how each side computes its answer) AND the *input
    resolution* axis (how each side decided WHICH question to ask). A
    single-sourcing refactor closes the second axis silently, because the
    diff only shows a call being re-pointed and both derivations plainly
    remain distinct. This is the *shared-input* sibling of the
    `coding-standards` rewire-demotion clause, which covers only the
    *caller* case (the survivor calls the other side); here neither calls
    the other and the gate is still demoted. Tell in review: the test body
    builds a domain object (`SelfPairedDeck.mirror(axis).motion`, a
    `RigidMotion`, a config) and passes it to the SUT while ALSO using it
    for the reference — previously the SUT resolved that object from its own
    tag (`"x"` → an index → its own construction). The refactor is usually
    right (one source is the point); what must not survive is the
    docstring's unqualified "genuinely independent routes / neither consults
    the other". (2026-08-07, G6.3 §7d.3.
    `tests/geometry/test_specular_response_pins_to_geometry.py` compares a
    geometric permutation against `quad.ordinate_permutation(motion)`; before
    the retirement the response side was `quad.reflection_index(axis)`,
    which resolved the LETTER through the quadrature tier's own
    `_resolve_axis_to_index` → `RigidMotion.reflection(normal=eye[axis])`, so
    the gate cross-checked the axis-letter convention across two tiers. `[M]`
    mutating `_mirror_motion` to swap the x and y letters leaves that file
    **15/15 GREEN** while reddening **78** sibling gates whose reference
    resolves the letter through its own literal. The class was preserved
    tree-wide only because the new reference helper deliberately keeps a
    LOCAL `_AXIS_INDEX` map instead of importing the production
    `AXIS_NAMES`.) Review rule: for every "independent routes" docstring,
    enumerate the shared objects flowing INTO both sides and name, in the
    docstring, which convention the gate can therefore no longer see — and
    check that some other gate still sees it, by mutating the shared
    resolution and confirming *something* reds.

    ⭐ **The shared object need not be BUILT by the test — a cached
    production PROPERTY of the SUT, reached by both sides, is the same
    defect.** Everything above assumes the test constructs one input and
    hands it to both sides. When the reference is a production composite
    (a frame form, a conjugated product) whose middle factor is the SUT's
    own cached `isotropic_energy`, both routes read the SAME instance and
    "neither side calls the other, and the test builds nothing shared" does
    NOT establish independence. ⟹ enumerate the objects each side READS,
    cached properties included, and expect the answer to differ per
    parametrised row — two ids that read as equally strong may not be.

    > `[M]` 2026-09-05, ORPHEUS CS4c step 5 (qa review F-3).
    > `test_the_l0_conjugation_identity_is_bit_exact_on_the_base[F]`
    > compares `F.apply(ψ)` against `F.full_fission_kernel.apply(ψ)/W`, and
    > `full_fission_kernel = conjugate(FissionMomentOperator(self.isotropic_energy))`
    > — the SAME `IsotropicFission` the fast path applies. Perturbing
    > `IsotropicFission.apply` by 1e-7: the `[S_L0]` row reds (its reference
    > builds `LegendreMomentTransfer` from the DATUM), the `[F]` row stays
    > green. The row still pins `M₀` against `∫ψ dΩ` and `R₀` against the
    > `/W` broadcast — a real claim — and F's energy binding IS pinned by the
    > sibling reciprocity row; the defect was the PLAN's pillar ("built by
    > hand in the test, `frame.conjugate` is NOT called"), which over-claimed
    > independence the shipped reference does not have.

## AP23 invariance sizing

23. **NEVER** size an A-vs-B **invariance** gate ("the answer is bit-identical
    when knob K is varied") by the breadth of what it computes — **instead**
    enumerate the production lines that **READ K**, and prove the fixture makes
    **each one** discriminating. An invariance gate compares two runs of the
    *same* code, so by construction it is blind to every mutation that is not
    K-dependent: its entire reachable coverage is the set of K-readers, which is
    usually one or two lines and is *enumerable by grep*. Two traps travel with
    it, and both point the flattering way. (a) **The ordinary positive control
    does not work.** A deliberately catastrophic mutation (replace the kernel
    with the identity) leaves the gate fully GREEN — correctly — so vv #17's
    "one mutation must redden many gates" reads as *the harness is dead* when it
    is merely the wrong control; the control for an invariance gate must itself
    be **K-dependent** (canonically: neuter the knob so A and B become the same
    object, and require the *activation* leg to red). (b) **The fixture can
    annihilate the only K-reader**, and then the gate is green forever while
    reading like a strong bit-identity acceptance criterion — Mode 12, at the
    fixture rather than the functional. (2026-08-07, SN G6.5's
    `TestFacePackingOrderIsBookkeeping`: K = the `FaceLayout` face packing
    order; the **only** line in the whole realization path that reads a flat
    offset is `slot.slice_view(metric_flat)` in `AngularTraceSpace._face_spaces`.
    The fixture's two faces are `("xmin","xmax")` — **the same axis** — so
    `|Ω·n| = |μ_x|` and the two slots' metric slices are bit-identical.
    `[M]` a wrong-slot read (offset 0 instead of the face's own) executed on
    every face at every offset and produced `changed=False` every time: **0 of
    10 rows red**. The same mutation with a y-face in the layout moves
    `Γ₊(xmax)`'s weights by **max |Δw| = 0.963**. A genuinely K-dependent
    in-class mutation elsewhere — shifting the deck permutation iff the face's
    slot is not at offset 0 — did red 4 of 10, so the gate has teeth, just not
    over its own subject.) Review rule: for every invariance claim, demand the
    K-reader list in the docstring, one measured mutation per reader, and an
    explicit note naming the rows that **structurally cannot** see K (here
    `vacuum` and `lambertian` never reach the deck kernel) so their green is not
    silently counted as coverage — the #20 row-inflation rule, applied to an
    invariance gate's parametrize grid.

## AP24 metric adjudicates

24. **NEVER** let a metric ADJUDICATE between design candidates until you have
    validated the metric against the mechanism it claims to measure —
    **instead** run three checks first: the **BASIS** check, the
    **RANK-CORRELATION** check, and the **cost-against-alternatives** check.
    #1–#23 all concern *gates* (a pass/fail coverage claim about one
    implementation). An **adjudicating instrument** is a different object: it
    *ranks* two or more candidate formulations, and a design is then built on
    its verdict. It fails in ways a gate cannot, because there is nothing to
    mutate — every candidate is "correct code", and the instrument is green for
    all of them by construction.

    (a) ⭐ **The BASIS check — probe your test functions against the problem's
    SYMMETRY GROUP and against what the DISCRETIZATION can represent.** A
    ranking instrument feeds trial modes through the candidates and scores the
    residual. If one trial mode is *forbidden by the problem's symmetry* — or,
    worse, is a mode the discretization cannot even represent — the instrument
    ranks candidates by their behaviour on a mode the solver will never see, and
    the ranking can invert against the physically-realisable modes. **The trap
    is that the natural audit is of the WEIGHTING, and the weighting can be
    provably robust while the basis is wrong.** Checking five weightings and
    finding the ranking stable feels like a validation; it validates nothing
    about the basis. The cheap discriminator is a **discrete moment of the trial
    mode on the actual rule**: if `Σ w · f(Ω_n) ≠ ∫ f dΩ` by an O(1) amount, the
    rule does not represent `f`. (`[M]` 2026-08-11, SN Q5.6.4: an instrument fed
    `η` and `ξ` through the cylindrical angular closure, justified as "a
    P1/diffusion-limit flux is affine in the direction cosines". True of P1 in
    3-D; false here — the `ξ→−ξ` reflection makes `J_φ ≡ 0`, so the diffusion
    limit at a level is affine in `η` ALONE, and BMC 2010's own Eq. (1) is
    `φ/4π + 3J_r μ/4π`, one cosine. And on a σ_y-**folded** rule every node has
    `ξ > 0`, so `quad.mu_y` samples `|ξ|` — measured `Σwξ = +6.703 ≠ 0` folded
    vs `0.000` unfolded. Re-run on the realisable basis `{cos mω}` = Chebyshev
    in `η`, the ranking **INVERTED**: the candidate the instrument ranked
    worst-but-one became best, by ≈2× on *every* harmonic `m = 1…4`.)

    (b) ⭐⭐ **The RANK-CORRELATION check — with ≥3 candidates, ask which
    mechanism the metric is actually ordered by.** Tabulate the metric beside
    one column per candidate *mechanism*. A metric in perfect rank correlation
    with mechanism A and **anti**-correlated with mechanism B **cannot
    adjudicate a question about B**, no matter how authoritative it looks or how
    tight its tolerance. This is Mode 12's stabiliser question asked
    *empirically* instead of algebraically, and it is available whenever you have
    three candidates — which a design debate always does. (`[M]` same case: the
    shipped `1.2e-1` flux-shape cross-check ordered the four candidates
    `6.59e-02 < 1.02e-01 < 1.27e-01 < 1.44e-01`, which is **exactly** the order
    of their recurrence error-amplification `{1.00, 1.00} < 9.44 < 40.7` and the
    **reverse** of their closure accuracy `1.99e-02 … 5.64e-02`. So the number
    the whole campaign was steering by measured the *recurrence*, while the
    campaign was arguing about the *closure*.) Corollary: when the ranking is
    explained by a mechanism nobody was debating, the debate was mis-framed —
    report the mechanism, do not pick a side.

    (c) **The cost-against-alternatives check — an "honest cost" must be
    measured against the candidates on the table, not stated in the abstract.**
    A caveat of the form "candidate X has the usual exposure to Y" is a claim
    about a *comparison*; if every candidate has Y and X has the LEAST of it,
    the caveat is inverted and it argues *for* X. (`[M]` same case: "`τ ≡ ½` has
    the diamond scheme's usual positivity exposure" — the destabilising
    coefficient is `(1−τ)/τ`, so `τ ≡ ½` is the *least* exposed of the four;
    measured `min ψ̂ = −24.2` vs `−77.2` and `−230` for the alternatives on the
    same profile.)

    (d) ⭐⭐ **The ZERO-SET check — solve `instrument(candidate) = 0` for the
    candidate. If the solution IS the incumbent, the instrument measures
    distance-to-the-incumbent and cannot adjudicate.** This is the cheapest of
    the four (one line of algebra, no run) and it catches the most dangerous
    shape: an instrument that is reference-free, tight, pointwise, honestly
    parameter-loaded, *and* confirms whatever is already shipped — at every
    order, on every fixture, for every material. It is (a)'s BASIS failure seen
    from the other side: (a) asks whether your trial functions are realisable,
    (d) asks whether they lie in the SUT's own kernel. The tell in a plan is an
    instrument whose scoring rule restates the design's *defining property*
    ("τ is the barycentric coordinate", "the limiter keeps `a ∈ [½,1]`") rather
    than a consequence of it. (`[M]` 2026-08-12, SN #235: the curvilinear
    angular closure `ψ_m = τψ̂_+ + (1−τ)ψ̂_−` is EXACT on `span{1, μ}` — `4.4e-16`
    cylinder, `8.9e-16` sphere, every order — *because* τ is defined as the
    barycentric coordinate in the radial cosine. Three proposed instruments have
    that same zero set: the diffusion-limit test (the diffusion limit's angular
    content IS `span{1, μ}`), the shipped anisotropic-MMS ansatz `A(r) + B(r)η`
    (affine in η at every `r` ⟹ closure residual `6e-17` ⟹ the fixture is in the
    incumbent's kernel and cannot rank it), and the η-weighted closure defect. A
    fourth, Reed & Lathrop's `|τ−½|/w`, has zero set `τ ≡ ½` — i.e. it is one
    *candidate* wearing a criterion's clothes.)
    ⭐ **Corollary, measured, on the choice of the graded FUNCTIONAL: an
    INTEGRATED functional admits signed cancellation and can rank garbage above
    production.** When the error under test enters the measured quantity through
    a sum (`φ = Σ_n w_n ψ_n`), a candidate can score well by CANCELLING against
    the error floor it does not control, and the sum's own moment identities can
    annihilate whole error classes exactly. Grade the un-integrated field when
    one exists. (`[M]` same case, same solves, `n_φ=64`: the scalar-flux `L2`
    ranks two garbage τ permutations **1.6×/2.0× BETTER** than production and is
    blind to a 2 % τ jitter (`1.04×`); the angular-flux `L2` ranks the same
    permutations `17.0×`/`8.0×` WORSE and resolves the jitter at `3.96×` —
    dynamic range `2.1×` vs `40×`. The mechanism is exact: the closure defect is
    `∝ cos(mω_m)`, whose discrete azimuthal moment vanishes, so the *identity
    that makes the manufactured reference closed-form is the identity that
    annihilates the defect in the graded quantity*.)

    (e) ⭐ **The REGIME check — is the fixture's PHYSICS in the regime where
    the differentiating mechanism is even active?** Distinct from (d), and the
    two are easy to conflate: (d) asks whether the fixture's *content* lies in
    the scheme's kernel; (e) asks whether the fixture's *material and geometry*
    put the scheme's claimed advantage in play at all. A fixture can pass (d)
    outright — rich angular content, large residual, healthy dynamic range —
    and still rank the candidates by a mechanism nobody cares about, because
    the one under test is dormant at that `c`, that optical thickness, that
    cells-per-mfp. The check is two lines of arithmetic on the fixture's
    cross-sections, not a run. (`[M]` 2026-08-12, SN #235: the flagship
    anisotropic-cylinder MMS is `σ_t = 1.0`, `σ_s = 0.5`, `R = 5` ⟹ `c = 0.5`
    and `Σ_a·R = 2.5` — **half of every collision is an absorption**, about as
    far from the diffusion limit as a scattering problem gets. The property
    under test, first-order diffusion-limit consistency, is *by definition*
    invisible there. That fixture had therefore failed the campaign **twice
    over, independently**: (d) because its exact solution `A(r) + B(r)η` sits
    in the closure's kernel, and (e) because its material sits outside the
    regime — two unrelated blindnesses in one artefact, either sufficient to
    void every accuracy conclusion drawn from it about the closure. And they
    have OPPOSITE fixes, so diagnosing only one leaves you confidently wrong:
    (d) wants richer angular content, (e) wants a different material.)

    Review rule: an instrument that decides a design owes the same evidence a
    gate owes — but the evidence is *these five checks*, not a mutation. And
    when a plan says "instrument I ranks the candidates as …", the first
    questions are **"what basis, what is the ranking correlated with, what
    scheme makes this instrument read zero, and is the mechanism under test
    even awake in this fixture?"**

## AP25 bundled mechanisms

25. **NEVER accept a "this artefact was unaffected because X" explanation for
    a null result when the change RETIRED MORE THAN ONE MECHANISM — instead
    enumerate every mechanism the commit bundled and check the artefact
    against each.** A re-baseline is the one place a reviewer *welcomes* a
    null result, so its stated reason gets the least scrutiny — and the reason
    is usually written about the mechanism that MOTIVATED the change, not the
    one that also rode along. The output is worse than a wrong number: it is a
    durable **certificate of blindness**, and it points the flattering way
    (future work skips the fixture the certificate exonerates). It survives
    even when the null result itself is correct, because a right conclusion
    reached through a void argument reads as verified.
    ⭐ The grep-decidable tell is a commit that both (a) reports a per-artefact
    null and (b) names ONE cause, while its own diff/subject line contains a
    conjunction — *"the partition is taken in ω **and** the absorber was
    compensating for it"*, *"rename X **and** drop the shim"*. One conjunction
    in the subject = two mechanisms = two checks owed per artefact.
    (`[M]` 2026-08-12, SN Q5.6.4 cylinder re-baseline `39b46a31`. Its case
    list, landed in `tests/sn/regression/_generate_snapshots.py` and carrying
    an `[M]` marker, reads: *"folded_2x4 has M = 2 ordinates per level, and
    `[M]` the new ω-midpoint partition is BIT-IDENTICAL to the retired
    η-midpoint one at M = 2 … **So this case's tau did not change at all**,
    and no M = 2 fixture can ever see a partition change."* The partition half
    is TRUE and measured — the interior edge is `5.0e-17 ≈ 0`. But the same
    commit retired **two** things, and at `M = 2` the *absorber* was the
    binding one: raw `τ = 0.2929 < ½`, so the retired `max(0.5, min(1.0, τ))`
    had been clamping it. Measured on `folded_product(2, 4)`:
    `τ: 0.292893 → 0.5`, i.e. **`Δτ = 2.071e-01` on every level** — τ changed
    a great deal. The artefact genuinely did not move, but for the *other*
    reason the same commit gave for a *different* case (a homogeneous medium's
    near-flat ψ nulls the redistribution). Cost: the sentence teaches every
    future session *"M = 2 ⟹ τ-insensitive"*, and a whole-suite differential
    refutes it directly — the `n_φ = 4` (`M = 2`) row of
    `test_cyl_tau_equals_the_ANALYTIC_closed_form_not_the_chord_convention`
    is among the 32 gates that redden when the old τ is restored. A later
    closure change touching only a limiter would be waved past every `M = 2`
    fixture on the strength of an `[M]` marker.)
    This is the test-design face of `plan-authoring` §2's sharpening — the
    marker certifies the half that was measured and lends its authority to the
    half that was not — and it is why the **audit denominator** matters twice
    over: `39b46a31` also sha256-swept *"all 23 snapshots"* and concluded
    *"these are the only two that changed"*, having looked in ONE directory
    while 7 further frozen references in three other directories had moved.
    ⟹ Two review rules, both cheap: **(i)** for every bundled mechanism,
    one null-check per artefact, and write the mechanism's name beside the
    number; **(ii)** a re-baseline's radius is the set of frozen REFERENCES
    reachable by the changed code — run the whole module tree, and grep the
    non-`.npy` carriers too (`sha256`/`hexdigest` literals, expected-value
    formulas inside test bodies, prose naming the fixture's quadrature).

## AP26 output path

26. **NEVER gate a claim about the PATH by asserting the OUTPUT — instead
    instrument the path.** A short-circuit, a cache, a no-op fast lane, a
    "costs nothing when there is nothing to do" promise: each says the code
    *did not do something*, and a correct implementation and one that does the
    work and then throws it away are **indistinguishable in the return value**.
    So the gate is green under both, its name says the promise is covered, and
    the promise is free to rot. This is #19's discrimination rule moved from a
    gate's *sensitivity* to a function's *route*, and it has the same tell: the
    reading you naturally take is the one that carries no information.
    ⟹ The witness must observe the route, and it is usually one counter —
    monkeypatch the expensive call in the module's own namespace and assert the
    call list is empty. Choose the *cheapest observable on the skipped path*
    (the parse, the subprocess, the `stat`), not a timing measurement, which is
    a flaky proxy for the same question.
    ⚠ The companion trap: assert **identity**, not equality, whenever the
    promise is "unchanged". `to_json(json.loads(x)) == x` for anything `to_json`
    produced, so an equality assertion cannot see a needless round-trip; `is`
    can.
    (`[M]` 2026-08-16, nexus Track 1.1. A staleness pass at the MCP tool
    boundary promised to be *absent to the byte* on a fresh graph. Two gates
    asserted it — one by identity, one by the flag's absence — and the mutation
    battery came back **GREEN** on both arms it was meant to cover: deleting the
    early return still returned the same object, because the walk marks nothing
    on a fresh graph; and forcing a re-serialisation was unreachable, because
    the early return preempts it. Two gates were then written that DO
    discriminate — a parse counter for the route, and an identity assertion in
    the *dirty-tree-but-unaffected* case, which is the only configuration that
    actually reaches the walk. Both reddened their arm. Note the second one is
    the general lesson in miniature: **the case that exercises a fast lane's
    fallback is not the case the promise is about**, so a claim of this shape
    needs a gate on each side of the branch, never one on the happy side.)

## AP27 workaround idiom

27. **NEVER treat a retired type's leftover WORKAROUND IDIOM as stale prose to
    re-word — instead ask what error class the detour's functional
    ANNIHILATES, because the idiom is a coverage claim, not a style one.**
    When a type restricts an operation (`flux + flux` refused, a unit forbidden,
    a nullable barred), every gate that needed that operation was written
    through a **detour** — an algebraically equivalent spelling that the type
    *did* permit. Retire the type and the restriction vanishes; the detour does
    not, because it is still green and still "correct". So a retirement sweep
    reads it as vocabulary and fixes the comment. That is the trap: **the
    detour is a DIFFERENT functional from the direct statement, and a weaker
    one** — it was chosen for expressibility, never for sensitivity, so nobody
    ever computed its invariance group. This is Mode 12 asked of the **idiom**
    rather than of the fixture or the metric, and it is the retirement-audit
    intersection the `coding-standards` three-search rule cannot reach: grep
    finds the words, and the words are not the defect.
    ⟹ **The check is design-time and costs ten lines: model both functionals in
    pure arithmetic and evaluate them on the error class the gate is credited
    with.** No SUT import, no fixture, no mutation — which also makes it the
    right instrument when the tree is under concurrent edit and mutating a
    production file is unsafe (`process-discipline`'s mutation-revert hazard).
    ⚠ And when a sweep claims to have fixed such a site, **verify it moved the
    ASSERTION, not the prose** (`git show <old>:<file>` vs the working tree). A
    prose-only fix leaves the blindness wearing a corrected comment — strictly
    worse than before, because the file now reads as audited.
    > `[M]` 2026-08-19, ORPHEUS campaign-1 CS3-R. The retired `FluxDisplacement`
    > torsor made `ψ₁ + ψ₂` a `TypeError`, so **five** operator-linearity gates
    > spelled additivity as the affine detour
    > `op(ψ₁ + λ(ψ₂−ψ₁)) = (1−λ)op(ψ₁) + λop(ψ₂)`. Affine maps *preserve affine
    > combinations*, so the detour is **exactly** blind to an affine regression
    > `A(x) = Lx + q`. Pure-numpy probe (`n = 6`, random `L`, `λ = 0.7`, no SUT
    > import): retired form reads `4.440892e-16` at `q ≠ 0` — **bit-identical to
    > its own `q = 0` control** — while the direct `A(ψ₁ + ψ₂)` reads
    > `1.288361e+00`. No tolerance, refinement, or fixture change could ever
    > have exposed it. Corroborated independently by the codebase's own battery
    > (`tests/sn/operators/test_declared_law_is_linear.py`): under an `affine`
    > mutation, **19 of 69** rows reddened and *"neither base-point-independence
    > row reddened"*, so re-spelling to the direct form **upgraded** that
    > battery by two rows. ⭐ The sharpest residue: one of the five had its
    > assertion upgraded and its body comment rewritten while its **docstring**
    > still described the old, blind contract — the file now misdescribes a
    > correct gate in the direction that invites someone to "restore" the
    > weaker form (anti-pattern #21's self-contradicting-file aggravator,
    > created by the repair itself).

## AP28 optional metadata

28. **NEVER design a guard against an operand's OPTIONAL METADATA without
    probing that field on a PRODUCTION instance — instead build the object the
    production path builds and read the attribute.** A guard of the form
    *"refuse unless `operand.<record>.<property>` agrees with mine"* type-checks,
    reviews as structural, and reads as the strongest possible construction-time
    refusal. It is inert wherever `<record>` is `None` — and an optional
    metadata slot is exactly the kind of field a *convenience factory* populates
    and a *composite factory* forgets, so the inert region is systematically the
    COMPOSITE, i.e. production. This is Mode 8's SIGNATURE-tautological class one
    layer in: there the producer's *signature* cannot admit the knob and
    `inspect.signature` answers it; here the signature is correct, the annotation
    is correct, and only a **runtime read of a production-built instance**
    answers it.
    ⚠ The aggravator, and the reason it is its own item: the guard is loaded
    exactly where it is cheapest to test and inert exactly where the traffic is.
    A hand-built fixture reaches for the *simple* constructor — the axis-built /
    `of_axes` / `from_parts` one, which populates the record — so every gate the
    author writes goes red on demand and the guard ships certified. Nothing in
    the test suite ever holds the object production holds.
    ⟹ Two checks, both one line. **(a)** Construct the operand the way the
    production call site constructs it and print the field
    (`sn.full_field_space.axes`), rather than reading the class. **(b)** Count the
    production construction sites by which factory they use, and write the
    fraction into the guard's docstring — a guard live on 5 of 13 bindings is a
    different object from one live on 13, and only the written fraction stops the
    next audit reading it as universal. If the inert fraction is non-zero, the
    honest guard is the one keyed on what the object *always* carries (a shape
    slot), with the richer check named as the successor and its enabling phase
    cited — otherwise the shape arm becomes a silent twin the day the record is
    populated.
    > `[M]` 2026-08-20, ORPHEUS campaign-1 CS4a round-1 review. Two of three
    > independent design assemblies specified the kernel↔space conformity refusal
    > as *"kernel `ng` == the space's `EnergyAxis` shape"*. `FunctionSpace.axes`
    > is `Optional[tuple[Axis, ...]]` defaulting to `None`
    > (`numerics/space.py:196`); `of_axes` populates it, and
    > `FullFieldSpace.from_blocks` — the ONLY producer of the SN/diffusion
    > composite — passes `name`/`shape`/`interior_space`/`trace_space` and **not
    > `axes`** (`numerics/spaces/full_field_space.py:238-243`), while
    > `SNMesh.full_field_space` builds its interior block as a bare
    > `FunctionSpace(name="sn_bulk", shape=…, inner_product_weights=…)`, also
    > axes-less. Measured on a live production instance: the quotient
    > `bulk_space` carries `axes = (EnergyAxis(...), Axis('spatial', (1,)))`; the
    > SN `full_field_space` carries `None`. `[M]` **7 of 13** production
    > constructions of the four rebound operators thread `full_field_space`
    > (⛔ corrected same day from "8 of 13" by the sibling review's runtime
    > probe: `from_mesh`'s chain NAMES `full_field_space` but `MaterialMesh`
    > has none, so that site resolves `bulk_space` — the rule's own
    > build-the-operand directive is what caught its own founding number), so
    > the "construction refusal" would have been inert on the majority path,
    > green, and unfalsifiable — while every axis-built test fixture reddened it
    > on demand. Caught by building the operand rather than by reading either
    > design; a third assembly had reasoned to the same conclusion from the class
    > definition and was the only one to state it.

    ⭐ **The TEMPORAL twin, and it is the one a RETIREMENT creates: a guard whose
    predicate reads an attribute through a DEFAULTED `getattr` goes silently
    inert the day that attribute retires.** #28 above is about a guard that is
    born inert on the majority path; this is about a guard that is *live today*
    and is killed by an unrelated refactor — and the defensive spelling is
    exactly what makes the death silent instead of loud.
    `getattr(x, "mesh", None)` was written so a duck-typed carrier without the
    attribute would SKIP the check; after the attribute retires, *every* operand
    skips it, the branch is unreachable, and nothing fails. An `AttributeError`
    would have been a loud, one-line fix.
    ⟹ **The retirement audit's fourth search** (beside graph callers, text grep,
    direct constructors — `coding-standards`): grep the retiring name **inside a
    defaulted `getattr`/`hasattr`**, and re-key every hit IN THE SAME STEP with a
    red witness (`plan-authoring` §6c). And the design-time countermeasure:
    prefer a guard that reads the attribute DIRECTLY when the operand's type
    guarantees it — a defaulted `getattr` in a *condition* is a coverage claim
    with a hidden expiry date.
    > `[M]` 2026-08-21, ORPHEUS campaign-1 CS4b verification design (pre-carve).
    > The design record named ONE such site (`transport/full_field.py:265-274`,
    > the composite cross-slot mesh-identity gate). `grep -rn
    > "getattr([^,]*, *['\"]mesh['\"]"` over `orpheus/` + `tests/` returns
    > **four**, and the fourth is a second live guard nobody had looked at:
    > `sn/solver.py:338-345` refuses a bare System-A residual on a
    > starting-direction-carrying mesh, and its own comment says it exists to
    > prevent the Mode-12(b) blindness its removal re-opens. `[M]` it has **no
    > test witness** (three distinctive message fragments, 0 hits in `tests/`),
    > so its silent death would also be an invisible one. Same session, the
    > related measurement that makes the rule worth carrying: of the **22**
    > mesh-identity guard call sites the phase re-keys, **8 redden nothing**
    > across a 3936-row measured denominator — and the two sharpest are the
    > `apply_transpose` and second-`solve`-arm TWINS of witnessed forward arms,
    > i.e. exactly #17's "the site most likely to be miswired is the one whose
    > operand expression differs from its siblings'".

## AP29 traffic census

29. **NEVER accept a design that replaces runtime dispatch with a
    construction-time KEY on the strength of a class-level inventory of ARMS —
    instead run a per-INSTANCE traffic census: instrument the boundary and log
    `(bound key, observed operand type)` over one real workload per family.**
    An inventory of a dispatch table enumerates what the class *can* receive;
    a design that collapses the table asserts what each *instance* actually
    *does* receive, and those are different populations. The inventory is
    cheap, static, greppable and reads like evidence; the census needs a run.
    So the design ships with the arms counted and the traffic unmeasured, and
    the failure is silent in the flattering direction: the selected arm is
    *an* arm the class supports, so nothing type-errors — production simply
    stops reaching the body it needs.
    ⚠ Three distinct ways the key fails, all invisible statically: (a) **wrong
    arm** — the instance is bound to key K₁ and fed only carrier C₂;
    (b) **non-determination** — one instance, one key, two carrier families in
    one solve; (c) **asymmetric arrow** — the arm accepts a typed carrier and
    returns a bare one, so no single key names its domain AND codomain.
    ⟹ The census is ~15 lines and needs a **positive control** (#17): wrap
    `cls.__dict__["apply"]` through the descriptor protocol so
    `singledispatchmethod` still dispatches, log `type(operand).__name__` per
    `id(instance)`, and confirm the workload's headline number is bit-identical
    with and without the wrapper — an instrumentation that perturbs the run
    measures its own perturbation.
    > `[M]` 2026-08-20, ORPHEUS campaign-1 CS4a round-2 review. All three
    > independent assemblies proposed *"select the one apply body from the
    > bound space at construction; apply-time dispatch retires"*, each citing
    > the same static inventory of arms. Census over one solve per family
    > (SN k-eigen S4 2-region, 1-D diffusion, homogeneous k∞; control keff
    > `0.18764940308862563` bit-identical instrumented vs not): **6 of 12
    > production instances refute it.** `SNSolver.fission_op` is bound to
    > `sn_mesh.full_field_space` and receives **`ndarray` ×17 and nothing
    > else** (`sn/solver.py:1339` → `:1439`, reached from
    > `numerics/eigenvalue.py:420` — the k-eigenvalue outer iteration), so
    > selection from the space picks the composite arm and orphans the only
    > arm production uses. Diffusion's `IsotropicScattering`/`IsotropicN2N`,
    > both bound to one `FullFieldSpace`, receive `ndarray` ×27 **and**
    > `FullField` ×25. The SN iso pair, bound to no space at all, is fed
    > `ScalarFlux` ×225 while returning bare ndarray
    > (`isotropic_scattering.py:96-98`). One assembly had listed the census as
    > its own strongest self-attack and deferred it as "a plan, not a fact" —
    > run, it was the refutation.

    ⭐ **The FOURTH way, and the one a census will MIS-REPORT rather than miss:
    (d) NO arm — the instance is never applied at all, because a FUSED PARENT
    overrides the sum's body.** (a)–(c) all assume an arm runs and the key names
    the wrong one. Here the answer to *"which arm?"* is *"none, ever"*, and the
    key selects a body production never executes. The mechanism is invisible to
    both instruments: statically the operand is a legitimate member of an
    `OperatorSum`, so `apply` looks reachable; dynamically the parent's override
    never calls `b.apply`, it reads `b`'s **data**. ⟹ two consequences a census
    must state separately, because they point opposite ways:
    * **for the DESIGN** — there is no action body to select at construction, so
      the collapse buys nothing at that binding; the honest question becomes
      whether the operand is an operator at all or a data field wearing one;
    * ⛔ **for the RETIREMENT** — *zero applies is NOT zero consumers.* The
      object is load-bearing through an attribute read, so a "dead, retire it"
      verdict inferred from the traffic is exactly backwards. This is #17's
      red-set-by-IDENTITY clause with the polarity flipped: there, pins that
      only NAME a symbol prove no consumer exists; here, no `apply` traffic
      coexists with a live consumer one frame up.

    ⭐⭐ **And the measurement discipline the same census forces: count BODIES
    EXECUTED, not ARMS DISPATCHED.** An arm may itself be a *re-dispatcher* — it
    reads a sub-carrier off the composite and re-enters the same dispatcher — so
    one call is two counted arm entries and construction-time selection
    *relocates* the branch one frame in rather than removing it. The tell in a
    census is two arms with **exactly equal counts on every scenario**; the tell
    in the source is a `self.apply(...)` inside a registered arm.

    > `[M]` 2026-08-30, ORPHEUS CS4c step 0 (HEAD `2f44ed4e`; 11 production
    > entries; 23-verb denominator; every arm fired by an activation control;
    > all 11 headline numbers bit-identical instrumented vs not).
    > **(d):** `MultiplicationOperator` at `sn/coupled_system.py:446` is minted
    > **22 / 22 / 24 / 25 / 20** times per k-solve (once per outer) on
    > slab-SI / slab-Krylov / sphere / cylinder / 2-D — and **every instance is
    > silent, in all 9 SN scenarios, under BOTH inner solvers.**
    > `StreamingCollisionOperator` subclasses `OperatorSum` holding it as `b`,
    > and **overrides** `apply` (`sn/operators/streaming.py:723`) to call
    > `loss_action(self.sigma, psi)` with
    > `sigma = self.diagonal.coefficient.values` (`:712-719`);
    > `apply_transpose` is overridden the same way (`:744`). I had predicted the
    > silence was source-iteration owning the sweep and that Krylov would
    > re-enter `OperatorSum.apply` — **refuted**: the override is on the
    > composite, not on the solve strategy.
    > **Re-dispatch:** `ScatteringOperator`'s `FullField` arm runs
    > `self.apply(cast(...), psi.interior)` (`scattering.py:1189`), so every
    > composite apply scores twice — `FullField ×N` **and**
    > `AngularFlux|HarmonicMomentFlux ×N`, equal counts on all 6 scenarios.
    > Reading "4 arms" as four alternatives over-counts the bodies and
    > under-counts the branching.

    ⭐ **A census's NOT-RUN row is a measurement; its EXPLANATION is a separate,
    usually UNMEASURED claim — and the two candidate explanations need OPPOSITE
    repairs.** Never attach a reason to a zero-traffic arm without discriminating
    *"no consumer exists"* from *"this workload never reaches the consumer"* —
    instead enumerate the production BRANCHES on the path to that arm and check
    whether the workload took **both** sides of each. The two readings demand
    opposite work — a design decision (retire / declare a future consumer) versus
    one more scenario — so guessing produces a durable **false-dead certificate**,
    and it points the flattering way: "no consumer" reads as a finished
    investigation, and nothing prompts a re-check. ⟹ the discriminator is cheap
    and mechanical: for each zero row, name the `if`/`is None` that gates its
    consumer, and add ONE scenario on the unexercised side.

    > `[M]` 2026-08-31, ORPHEUS CS4c step 5 (re-run of the 2026-08-30 step-0
    > census). Step 0 recorded `IsotropicScattering.apply_transpose` and
    > `IsotropicN2N.apply_transpose` as NOT-RUN with the reason *"declared future
    > consumer = the adjoint diffusion chain (#281); no such entry exists at
    > HEAD"*. The consumer existed at HEAD and was neither diffusion nor future:
    > the **ray-system adjoint**, `sn/operators/radial_characteristic.py:1536`,
    > fed `S.isotropic_energy + N2N.energy` at `sn/coupled_system.py:591`.
    > `solve_sn_adjoint` branches on `radial_characteristic_field_space is None`
    > (`sn/solver.py:2919`) and step 0's workload had **only slab adjoints**,
    > which take the early return. Adding ONE curvilinear adjoint scenario
    > (0.9 s) moved both rows from 0 to **985 calls each** — **2 of step 0's 5
    > dead-arm rows were configuration artefacts**. ⭐ The aggravator: step 0's
    > own caveat said *"the residual risk is not another solver family, it is
    > another configuration of a driven family"* — the memo stated the right
    > hazard and then attributed 2 of its own rows the wrong way, because a
    > *stated* hazard reads as an *assessed* one (#13's tell, at census scale).

    ⭐ **A verb reads dead when a consumer is fed the operator's KERNEL instead
    of the operator — a fifth way, distinct from (d)'s parent-override.** (d)
    covers a fused parent that overrides `apply` and reads the child's *data*.
    This is the sibling case: **one slot fed two different LEVELS of the same
    abstraction**, so one sibling's verb is the hottest in the census and the
    other's reads dead — with no override anywhere. ⟹ before crediting a zero
    on a verb whose body is a one-line delegation
    (`return self.kernel.apply_transpose(...)`), grep the **delegate's** call
    sites too: a caller holding the kernel bypasses the verb while computing
    identical arithmetic. Retiring on that zero removes a real spelling of a
    live path; recording "no consumer" is equally wrong.

    > `[M]` 2026-08-31, ORPHEUS CS4c step 5. `RadialCharacteristicEmission` is
    > constructed twice in production and its `emission_kernel` slot is fed
    > `S.isotropic_energy + N2N.energy` — the **operators** — at
    > `sn/coupled_system.py:591`, but **`F.kernel`** — the
    > `TensorProductOperator` — at `sn/solver.py:2947`. So
    > `radial_characteristic.py:1536`'s `emission_kernel.apply_transpose(...)`
    > routes into each iso operator's own verb on the loss side (**985 calls
    > each**) and steps OVER `IsotropicFission.apply_transpose` on the fission
    > side (**0**), whose whole body is
    > `return self.kernel.apply_transpose(_values_of(chi))`. One
    > abstraction-level inconsistency, two opposite census verdicts.
    > ✅ REMEDIED 2026-09-04 by CS4c step 5: the posing feeds
    > `F.isotropic_energy` (the operator) and `FissionMomentOperator` calls
    > the operator's verbs — one level, both routes.

    ⭐ **An ARM counter is one frame too coarse: two operators can be
    ARM-IDENTICAL and BODY-DIFFERENT, and the identity reads as
    interchangeability.** Every clause above assumes the ARM is the unit of
    observation — (a) wrong arm, (b) non-determination, (c) asymmetric
    arrow, (d) no arm at all, the re-dispatch and the kernel-bypass. This is
    the remaining case: **the right arm fires, the counter is correct, and
    the body takes an early return the counter cannot see.** ⟹ **NEVER**
    report two bindings as "fed alike" from equal arm counts — **instead**
    instrument the BRANCHES inside the arm (the `is_X` predicate, the
    `if … is None`, the short-circuit) and report, per role, how many of
    those calls reached the body the arm is named for. The failure points
    the flattering way: equal counts read as *"one operator, two roles,
    nothing to decide"*, which is exactly the conclusion a role-unification
    carve invites. And the discriminating datum is usually a property of the
    bound **DATUM**, not of the space, the role, or the arm — so no amount
    of construction-time space analysis can recover it.

    > `[M]` 2026-09-04, ORPHEUS CS4c step 5b (HEAD `f90f7914`, 16 scenarios,
    > all 16 headline numbers bit-identical instrumented vs control). After
    > #426 step 2 made `ScatteringOperator` and `N2NOperator` thin roles of
    > one `TransferOperator`, the census read `apply[FullField] ×4687`,
    > `apply[AngularFlux] ×4422`, `apply[HarmonicMomentFlux] ×265` —
    > **byte-for-byte equal on both roles, in all 9 SN forward scenarios**.
    > Instrumenting the arm's own branch: on **12 of the 13** legacy
    > scenarios `N2NOperator.is_isotropic` is **True on 100 %** of those
    > calls and `build_aniso_source` returns `None`, while
    > `ScatteringOperator`'s returns a source every time. The cause is the
    > PAD: the fixture's `Sig2` is a length-1 stack, `TransferKernel.at_order`
    > appends exact zeros to reach the solve's `L`, and `is_isotropic` reads
    > the padded VALUES — so the (n,2n) role is *shaped* anisotropic and
    > *valued* isotropic. A positive control (a synthetic ℓ=1 `Sig2`) flips
    > it and moves k by **−36.9 Δk·10⁵**; the 421-group Be-reflected
    > production fixture flips it too (`build_aniso_source = SOURCE ×748`).
    > ⭐ The aggravator: without that ONE production-data scenario, every
    > (n,2n) anisotropic body in the census reads 0 **three days after #426
    > landed the anisotropy**, and a design round reading a correct arm table
    > would have concluded the path has no traffic.

## AP30 codomain constructor

30. **NEVER credit an "X is not data of this operation" claim from the
    ARITHMETIC — check the CODOMAIN constructor.** Purity / locality /
    diagonality tell you X is not needed to COMPUTE the action; they say
    nothing about producing the RESULT, and in a typed-carrier codebase the
    result's constructor is where the dependency actually lives.
    > `[M]` 2026-08-20, CS4a physics assembly: all four ORPHEUS SN energy
    > operators are spatially diagonal AND all read `mesh` off the operand to
    > stamp `…SourceSink.from_mesh(v, mesh)` / `zeros_on(mesh)` (≥11 sites),
    > with the bound space carrying no mesh on any block — while a production
    > guard asserts the thesis's opposite (`sn/operators/streaming.py:589`:
    > "its mesh is carried by its CrossSectionField coefficient"). The
    > locality argument was TRUE and the "a mesh is never data of the
    > interaction, only of the pullback" corollary FALSE at the output mint.

## AP31 single draw

31. **NEVER pin ``np.array_equal`` — or publish "bit-exact" — from a single
    draw's green reading** — **instead** sweep seeds (or prove the float
    re-association exact) before claiming the bit tier; a one-draw exact
    reading is compatible with a law that fails ``array_equal`` on 40 % of
    inputs, and the resulting gate is seed-fragile: green today, red on any
    innocent fixture edit, and its "BIT-EXACT" docstring reads as a stronger
    claim than the sweep supports. "Bit-exact" is a property of the DRAW
    until a sweep makes it a property of the fixture.
    > `[M]` 2026-08-24, CS4b S6 gates: ``R∘E = id`` pinned ``array_equal``
    > on one seed and documented "BIT-EXACT". Seed sweep on the same
    > fixture: **844 of 2000** seeds fail (worst rel 1.5e-16 — Σ w_n(φ/Σw)
    > re-association); idempotence fails 57/200. The shipped SN carrier IS
    > exact (200/200 — Σw = 2 exactly, symmetric weights), which is why the
    > production-facing rows keep ``array_equal`` honestly. Found by an
    > archivist audit, re-pinned ``nulp=1`` same day.

    ⭐ **The finite-roster corollary — for a SHIPPED finite family, probe
    every member; a ladder is for unbounded families.** The #13 ladder
    rule (break your steps' arithmetic pattern) is the right discipline
    when the family is infinite; when the population is an enumerable
    shipped set, any ladder is a SAMPLE of it, and the member you skip is
    where the counterexample lives — while "n ∈ {…}, k of k" reads as
    exhaustive.
    > `[M]` 2026-08-24, the section-divisor probe: gram-einsum vs
    > ``weights.sum()`` probed at n ∈ {2,4,5,6,16,64} — 8 of 8 equal,
    > pattern dutifully broken (odds, non-powers). The shipped family's
    > ONE divergent member is **GL8** (1 ULP), inside the probed range and
    > skipped; the published consequence ("the array_equal licence
    > survives") was false. Caught by the docs audit's independent census;
    > the correction shipped a falsifiable GL8 gate row pinning the bound.

## AP32 consistency leg

32. **NEVER rank scheme CANDIDATES by positivity properties alone —
    **instead** add a **CONSISTENCY leg**, because sign-preservation and
    monotonicity are BOTH blind to the perturbation that breaks it.  This
    is anti-pattern #5 ("convergence rate correct ≠ result correct")
    relocated from a *solver result* to a *design candidate*, and it is
    nastier there: at the design stage nothing is converging yet, so the
    reviewer has only structural properties to rank by — and the two
    structural properties that a positivity review naturally reaches for
    are exactly the two that cannot see a wrong limit.
    The failing candidate does not look sick.  It converges **cleanly**,
    at the right rate, monotonically, positively — to the wrong answer.
    The check is one line on the cell transmission `a(τ_opt)`:
    ``sp.series(a - sp.exp(-t), t, 0, 2)`` — or equivalently
    ``a'(0) == -1``, since `a` must reproduce `exp(-τ_opt)` to first order
    for the scheme to be consistent at all.
    ⟹ every candidate table ranking closures by "positive?" /
    "monotone?" owes a third column, and a candidate that passes the
    first two and fails the third must be struck, not ranked lower.
    > `[M]` 2026-08-26, the LD lumping family (#158 / #408). A derivation
    > memo proposed the `(λ, ν) = (0, ½)` member as "monotone, at first
    > order", and it was carried into an issue as a recommendation. Its
    > transmission is `2/((1+τ)(2+τ))`: `a(0) = 1` ✓, strictly positive
    > for all `τ_opt` ✓, `A⁻¹ ≥ 0` ✓ — and `a'(0) = **−3/2**`, not `−1`.
    > Refined over a fixed thickness (`Σ_t X/|μ| = 1`) at 10 / 100 / 1000
    > / 10000 cells it converges to `0.2367 → 0.2245 → 0.2233 → 0.2231`,
    > i.e. cleanly to **`e^{−3/2}`** instead of `e^{−1}`. Both properties
    > the memo checked are TRUE of it; consistency is a third property
    > neither sees. Caught by an archivist re-deriving the family while
    > documenting it — not by the review that proposed it. The corrected
    > statement: `a'(0) = −1` forces **ν = 1 − λ** (a ONE-parameter
    > family, not two), monotonicity forces `λ ≤ 0`, and the nearest
    > monotone consistent member is `(0, 1)`, `a = 1/(1 + τ_opt/2)²`
    > (independently re-derived: `a'(0) = −1`, and `0.3769 → 0.3688 →
    > 0.3680 → 0.3679 = e^{−1}`).
    ⭐ The review tell, and it is cheap: **a positivity claim about a
    scheme is a claim about its NUMERATOR's roots; a consistency claim is
    a claim about its first DERIVATIVE at zero.** Two different objects —
    so no amount of care about the first can substitute for looking at
    the second.

## AP33 fact twice

33. **NEVER let a fact recorded for ONE job be spent on ANOTHER** — a
    symmetry (or tolerance, or tag) written down because predicate A needs
    it must not be read by predicate B asking a different question —
    **instead** give each job its own field, and name the field by the job.
    The failure is silent in both directions: B is answered from a value
    that was never a claim about B's question, and it stays green because
    the two questions happen to agree on every shipped input.
    The review tell is one grep: **a field read by two predicates whose
    docstrings ask different questions.**  If the second reader cannot say
    why the first job's answer is also its answer, it needs its own entry.
    > `[M]` 2026-09-03, #434 D1 (`qa`, the symmetry-machine review). The
    > quadrature registry recorded, per geometry, the closure a reflecting
    > FACE needs (`D_2h`: the ordinates must be closed under the coordinate
    > mirrors so a face reflection is an exact permutation). Stage 0 then
    > read the SAME field as a FOLD licence — "a rule folded by a subgroup
    > of what the geometry owes is admissible" — so `folded_product(4,8)`,
    > a σ_y fold, was admitted for `cartesian2d`, whose z-uniform solution
    > is even in μ_z ONLY: 2 of the 4 `(sign μ_x, sign μ_y)` sweep
    > quadrants were empty. The two questions ("must the nodes be closed
    > under σ?" and "is ψ even under σ?") agree on the cylinder and differ
    > on the plane, and every registered rule lives on the cylinder side.
    > Fixed at R3 by a THIRD field, `unspent` (the finite symmetry the
    > solution still has in the local frame), and a total coverage test
    > `H ⊆ unspent · spent` on the group the rule was folded by.

## AP34 brute force

34. **NEVER credit a "brute-force control" — or any second implementation
    offered as an independent oracle — on the strength of its NAME**
    — **instead** check that its BODY is structurally independent of the
    thing it checks, and do it mechanically: compare the two ASTs after
    α-normalisation (rename every local to a canonical placeholder). Two
    functions that are α-equivalent are ONE implementation wearing two
    names, and their agreement is a tautology, not evidence — anti-pattern
    #7's shared-upstream-identity failure at the level of code rather than
    of derivation.  A reviewer reading the two bodies separately does not
    see it: each reads correctly, and the copy is what *made* them read
    the same.
    > `[M]` 2026-09-03, #434 C1 (`qa`). The symmetry module shipped an
    > invariance predicate and a certificate builder whose three
    > docstrings claimed "ONE closure"; the predicate's body INLINED a
    > character-for-character copy of the closure's lambda. Nothing
    > declared the copy, and the two agreed on every input by
    > construction. The α-normalised AST check found it in one pass; R2
    > of #434 then made the copy unspellable (the closure's `images_of`
    > became REQUIRED and the second body was deleted), and the structural
    > claim is now GATED — `_orbit_closure` has exactly one call site,
    > asserted by AST in `tests/numerics/test_invariance.py`.

## AP35 overloaded unit

35. **NEVER report a derived COMPARISON quantity by its unit name alone when
    that name is overloaded — instead write its definition beside it, and check
    whether the choice survives your own fixture set.** `pcm` = 10⁻⁵ and says
    nothing about what was divided by what: `Δk·10⁵`, `Δk/k₀·10⁵` and
    `Δρ·10⁵ = (1/k₀ − 1/k)·10⁵` are three different numbers differing by `k₀`.
    Below `k ≈ 1.01` they are indistinguishable, which is why the habit never
    forms; above it they diverge without limit. ⟹ the tell that this is
    load-bearing rather than pedantic is a **fixture set spanning a range of
    the normalising quantity** — when one exists, emit all three columns. Same
    family as #24's ranking instruments: a number adopted to *compare* owes more
    than a number reported to *record*.

    > `[M]` 2026-09-03, ORPHEUS #426. A study reported `−413.55 pcm` (absolute)
    > for a `k₀ = 1.0953` fixture and `−529.26 pcm` for a `k₀ = 1.5262` one. In
    > reactivity those are `−346.01` and `−228.00` — so *"the thicker reflector
    > makes the truncation worse"* is **TRUE** in one convention and **FALSE**
    > in the other, a 2.3× spread, with no statement anywhere of which was
    > meant. Caught by the `qa` reproduction, whose own convention (relative)
    > disagreed with the claim's on the SAME k pair — the disagreement was the
    > finding. ⭐ The companion lesson from the same review: **two probes over
    > the same production code cannot see a SHARED convention** — a
    > reproduction certifies the arithmetic, not the premises both instruments
    > inherit. Test those against PHYSICS instead: here, strict
    > upper-triangularity of an energy-losing transfer matrix (8195 of 8195
    > entries) and the hard entrywise bound `|Σ_ℓ|/Σ_0 = |⟨P_ℓ(μ)⟩| ≤ 1` (a
    > stray `(2ℓ+1)` on ℓ = 1 would have read 2.9), each one `.todense()` and a
    > comparison.

## AP36 deselected catcher

36. **NEVER read a `catches("ERR-NNN")` / `verifies(...)` marker as coverage
    without reading the test's OTHER markers — a catcher deselected by the
    canonical invocation is a gate that cannot RUN.** Mode 8's nine classes are
    all about a gate that cannot *fail*; this is its dual. The catalogue counts
    the ERR as caught, the test genuinely reds when the defect is re-introduced,
    and `-m "not slow"` — the project's canonical gate — deselects it, so the
    regression lands green in every run that decides a merge, and `nexus
    errors` reports the entry covered. ⟹ when auditing an ERR's coverage, read
    the catcher's marker SET, and say plainly when an ERR's not-slow coverage is
    **zero**; a `slow`-only catcher is an ERR that needs a second, fast witness
    or an honest catalogue note.

    > `[M]` 2026-09-05, ORPHEUS CS4c step 5 (qa review F-5) — the SAME class
    > with a TYPE CHECKER as the absent enforcer: the C6 gate's static half
    > (`_c6_static_typing_pins`, `assert_type` rows "pyright-checked, never
    > run") is enforced by nothing, because the only pyright gate runs
    > `pyright orpheus/` and the pins live under `tests/`. `[M]` running
    > pyright on the file by hand: 1 pre-existing error in a pin nobody had
    > run (#452). A gate that never RUNS is not only a `slow`-marker shape.

    > `[M]` 2026-09-03, ORPHEUS #428 census. ERR-023's ONLY catcher is
    > `tests/mc/test_gaps.py:718` (`slow` + `catches("ERR-023")`). Under the
    > mutation ν₂ₙ: 2 → 1 the MC tree reads **39 passed / 0 red** at
    > `-m "not slow"`; the same test **FAILS in 84 s** run alone. Real teeth,
    > never engaged by the gate that matters.
