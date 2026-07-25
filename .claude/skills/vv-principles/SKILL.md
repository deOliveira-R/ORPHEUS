---
name: vv-principles
description: PROACTIVELY use when reviewing claims of correctness, designing verification plans, or evaluating whether evidence supports a claim. Provides the V&V hierarchy (L0–L3 + foundation), the 6 AI failure modes catalogue, the reference hierarchy by structural independence, anti-patterns, and the hierarchical claim taxonomy. Preloaded by qa, test-architect, numerics-investigator, and archivist.
---

# V&V Principles — claim taxonomy, evidence hierarchy, anti-patterns

This skill is the **decision instrument**. The pedagogy lives in
[reference.md](reference.md). Open this file during reviews,
verification-plan design, and bug triage. Open `reference.md` when
you need the philosophy or the worked case studies.

The corpus (Sphinx) home of this doctrine is
`docs/theory/verification/principles.rst` — the normative ladder /
pillar / claim-layer definitions render there. This skill is the
agent-side operational instrument: new failure modes and
anti-patterns land HERE first; the corpus page carries the doctrine
and its rationale.

---

## CRITICAL: Anti-patterns to flag immediately

Each line below is a redirect: **NEVER** do X — **instead** do Y. If
you see the left-hand pattern in a PR, claim, or doc, raise it before
any other review work.

1. **NEVER** claim verification on the basis of L4 agreement alone —
   **instead** require an L0–L2 evidence chain pointing at a
   structurally-independent reference. Two ORPHEUS solvers agreeing is
   _cross-implementation agreement_, NOT _correctness evidence_.
2. **NEVER** assert `np.allclose` against another solver in this
   codebase — **instead** match the claim to a reference at the right
   level (analytical for L1, MMS for spatial convergence, MC tally for
   L4 cross-check only after MC itself is verified).
3. **NEVER** accept a 1-group eigenvalue test as evidence of solver
   correctness — **instead** demand ≥2 groups. k = νΣ_f/Σ_a is
   flux-shape independent; 1G is degenerate. (An instance of the
   Mode-12 invariance-group lens — test-design table below.)
4. **NEVER** accept a homogeneous-only verification — **instead**
   demand at least one heterogeneous, mesh-refined, multi-group case.
   Flat flux nulls every redistribution and weight-cancellation term.
5. **NEVER** read "convergence rate is correct" as "result is correct"
   — **instead** verify the converged-to value. O(h²) to the wrong
   limit is still O(h²).
6. **NEVER** trust a reference that has not been traced back to a
   structurally-independent analytical or symbolic ground —
   **instead** treat it as **reference contamination** until the
   trace is shown. The most seductive failure mode: MC vs MC,
   CP vs unverified MC, method-of-images converged to the wrong BC.
7. **NEVER** treat "two derivations agree" as proof — **instead**
   check whether they are _structurally_ independent. ERR-032 (two
   antiderivatives both using `∫E_2 = 1 − E_3` instead of `½ − E_3`)
   agreed at 1e-39 because they shared the upstream identity, not
   because either was right.
8. **NEVER** accept "particle balance holds" as L0 evidence —
   **instead** require per-ordinate flat-flux residual. Telescoping
   sums hold by construction even with wrong per-ordinate balance.
   (A Mode-12 instance: the balance functional annihilates per-ordinate
   errors that cancel in the sum.)
9. **NEVER** conflate validation with verification — **instead** state
   which screw is being turned. If the equation itself is wrong,
   verification can pass cleanly; only L3 catches it.
10. **NEVER** accept "it produces reasonable numbers" — **instead**
    enumerate every term, isolate it, and verify sign AND magnitude.
    Sign-flipped small terms look reasonable.
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
    structural-independence requirement (L11) applies to ALL test
    design, not just numerical cross-checks.
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

---

## The 6 AI failure modes — mechanism and detection

These failure modes are mechanically explainable, NOT arbitrary —
they are the observable signature of sub-word tokenizer co-location.
**L0 verification is the only defense.** See reference.md §2 for the
mechanism (tokenization grounding, AI-targeted but not AI-exclusive).

| #   | Mode             | Example                                | Detection (L0 strategy)                                   |
| --- | ---------------- | -------------------------------------- | --------------------------------------------------------- |
| 1   | Sign flip        | `(a − b)` vs `(b − a)`                 | Heterogeneous eigenvalue diverges under refinement        |
| 2   | Variable swap    | `mu_x` vs `mu_y`; `SigS` vs `SigS^T`   | Per-ordinate flat-flux residual; asymmetric 2G inputs     |
| 3   | Missing factor   | Missing `ΔA/w`, `2π`, volume           | Fixed-source flux spike at r=0 vs `Q/Σ_t` analytic        |
| 4   | Wrong recursion  | `α_{m+1/2}` index drift                | Per-ordinate flat-flux residual                           |
| 5   | Index error      | `face[i]` vs `face[i+1]`               | Non-uniform mesh produces detectably different keff       |
| 6   | Convention drift | Definition site vs usage site disagree | 2G heterogeneous with asymmetric SigS — wrong group ratio |

The catalogued instances live in `error_catalog.md` (ERR-NNN entries).

---

## Test-design failure modes — when the test cannot see the solver bug

Modes 1–6 above are **solver bugs**: the code is wrong, the test must
catch them. The modes below are **test-design failures**: the solver
bug is real, but the test is structured so it cannot mathematically
observe the bug. These are mechanically distinct from 1–6 and require
a different defense (test review, not L0 verification).

| #   | Mode             | Example                                | Detection (test-review strategy)                          |
| --- | ---------------- | -------------------------------------- | --------------------------------------------------------- |
| 7   | MMS simplification bias | Curvilinear ψ_chosen = sin(πr/R)/W chosen isotropic-in-μ "to isolate the radial closure" — the angular redistribution term IS the sweep's hardest math, but the test cannot see it because it cancels by ansatz design. ERR-026 (curvilinear sweep WDD) is invisible to this MMS. | Every multi-dim test must declare which terms its ansatz **activates** AND which it **nulls**. If the nulled set includes a term covered by an active ERR-NNN, redesign the ansatz. Add an angularly-non-trivial companion case (e.g., ψ = (A(r) + B(r)μ)/W) so the redistribution path is exercised. **NEVER** ship only the simpler case. |
| 8   | Compiled-out assertion (runtime-mode strip) | A test asserts via a bare `assert` statement, but the suite runs under `python -O`, which strips `assert` to a NO-OP. The test collects, passes, and reports green — while asserting **nothing**. Bites hardest for always-on canary/sentinel gates: a tripwire that cannot trip is a *false green* worse than no gate. (ORPHEUS canonical invocation is `-O`; the SN sentinel set mixed bare `assert` with `np.testing.assert_*` — the bare-assert sentinels were inert under `-O`.) | Per gate, decide the runtime mode **explicitly**. A bare-`assert`-bearing gate MUST run **without `-O`** (or be rewritten to use `np.testing.assert_*` / `pytest.fail`, which are function calls and fire under `-O`). Review: grep the gate's tests for `^\s*assert ` and confirm the invocation does not pass `-O`. **NEVER** assume an assert fired just because the test passed under an unknown optimisation level. |
| 9   | Splitting / acceleration verified only in a degenerate (FP-coincident) regime | An iteration-only change — a splitting (Gauss-Seidel, σ_r-removal), a synthetic accelerator (DSA), a preconditioner — MUST NOT change the converged fixed point, only the rate. But the FP-invariance is often verified ONLY on a regime where the wrong formulation is *accidentally exact*, so the gate is blind to the real bug. Two ORPHEUS instances: (a) the σ_r-fold (#215) — `A_wg.solve(S_residual)` with a σ_r-sweep equals the true solve only for ISOTROPIC flux (`Σ_s0·I` vs `Σ_s0·P_iso`); exact on a fully-reflective uniform box, **46–56 % wrong** on vacuum / heterogeneous. (b) the octant-group G-S shared-face bug (ERR-056) — correct on an AXIS-ALIGNED quadrature (each face one octant), **wrong fixed point** on a diagonal/spherical cubature (shared faces). Both pass the degenerate gate and ship silent errors. | The FP-invariance gate MUST run on a config that BREAKS the degenerate coincidence: an ANISOTROPIC flux (vacuum / heterogeneous / streaming — not the fully-reflective isotropic box) AND, for angular-schedule changes, a DIAGONAL cubature (`lebedev` / `level_symmetric` — not an axis-aligned `product`). Assert the converged flux equals the UN-accelerated (Jacobi / plain-SI) fixed point to solver tolerance, separately from the rate claim. A synthetic accelerator is correctness-safe BY CONSTRUCTION only if its correction → 0 at convergence (DSA); a *splitting* is correct only if every consistent split shares ψ\* — verify it, don't assume it. **NEVER** gate a splitting/acceleration FP-invariance on the isotropic-reflective box or an axis-aligned quad alone. |
| 10  | Activated-but-unconstrained term (the term runs but the MMS is blind to its sign) | An MMS ansatz's Mode-7 declaration marks a term **activated** (its code path IS exercised — the rows are populated and consumed) — yet the test is still blind to a sign/magnitude error in that term, because the term enters the measured quantity as a HIGHER-ORDER-small forcing that gets absorbed below the convergence floor. Mechanically distinct from Mode 7: Mode 7 is *nulled* (the term cancels by ansatz design, code path NOT run); Mode 10 is *run-but-not-constrained* (the code path executes, but flipping its sign does not move the converged value above the O(hᵖ) floor / value band). ORPHEUS instance (#240 D5b-S4): the 2-D LD stress MMS feeds a non-zero slope-moment SCATTERING source `Σ_s·φ̂` through the LD slope-source rows (instrumented: slope moments 0.26/0.13/0.07, scattered `Σ_s⊗I` and consumed) — so the slope-source code path is genuinely exercised, BUT a sign flip on those rows (and even a 3× magnitude error) leaves the convergence order at ~1.97 and the value band passing, because `Σ_s·φ̂` is an O(h)-small DG-internal forcing whose error enters above O(h²). The "activated" declaration was true; the term was still unverified for its sign. | The Mode-7 activated/nulled declaration is NECESSARY but NOT SUFFICIENT — for every term the declaration marks **activated** AND that carries a sign/convention trap (a slope-row sign, a transpose, a recursion direction), MUTATION-verify the term is also **constrained**: re-introduce the exact sign/factor error and confirm the gate goes RED. If the mutation passes (order/value-band unmoved), the term is *exercised-but-unconstrained* — declare it so in the honest-scope note (NOT "verified", NOT "nulled" — a third state) and, if the trap matters, add a companion gate that isolates the term so its error is O(1) in the measured quantity (e.g. a fixed-source problem where that term is the DOMINANT forcing, not a higher-order perturbation). **When no such isolating regime EXISTS** — the term is localized and never the dominant forcing in ANY configuration (#251 Leg B: a boundary-trace transverse-slope sits below the bulk O(h²) discretization floor everywhere, so there is no "dominant-forcing" regime and an *improves-on-flat* leg is unachievable — a correctly-consumed slope can even make the converged value slightly WORSE) — the companion gate is UNAVAILABLE, and the complete resolution is the STRUCTURAL pair alone: (a) assert the producer threads the projected moment through at MACHINE PRECISION (the stamp/threading proof, with a leggauss-only / structurally-independent reference), AND (b) mutation-verify a CONSUMED sign flip moves the converged value O(1) above the solver tolerance (≫ a named `_CONSUMPTION_TOL`, NOT the sub-floor value band), paired with the no-op control leg (a scalar / zeroed input → byte-identical) that pins the asymmetry. There is then NO value-improvement leg to add — do not manufacture one (it would falsely RED a correctly-consumed term). **NEVER** read "the code path runs" or "the ansatz activates the term" as "a sign error in the term is caught" — only a red mutation proves that. |
| 11  | Gate-never-executes-the-rewired-path (the named "twin" is green AND its asserts fire, but it never calls the changed production code) | A refactor moves logic onto a NEW production reader (a helper, an accessor, a data field on a packet) and a closeout names an EXISTING gate — typically a slow apply/matvec/round-trip twin — as the bit-identity evidence for the new path. The gate is green and its assertions DO fire (not Mode 8), and the term it would test IS reachable in some configuration (not Mode 10) — but the gate's actual execution path NEVER calls the rewired reader, because the production consumer routes around it (reads the source array directly, uses a batched kernel, or the per-element method has zero non-docstring callers). Distinct from Mode 8 (assertion compiled out) and Mode 10 (path runs but error is sub-floor): here the rewired production line is simply *not on the gate's call graph at all*. The gate's green proves the UNCHANGED siblings are unchanged, not that the new path is correct. ORPHEUS instance (#236 Phase 2 B2): the c-fold moved `c_out=α_out/τ` onto `CellVisit.c_in/c_out` (stamped by new `SNMesh._make_cell_visit`), read ONLY by `DiamondDifference.residual` (diamond.py:308-309). The closeout named the 640s matvec-twin as proof the global-ordinate mapping is byte-correct — but a file-write sentinel in `DD.residual` proved it was **never called** across the entire twin (curvilinear matvec reads `closure.cell_contribution`→`_c_per_level` directly; the per-visit `DD.residual` has ZERO production callers — its lone `scheme.py` "caller" is a docstring). A c_in↔c_out swap in the stamp left the twin AND the full `sweep/core` suite green; only an in-process probe walking real `dag_walk` visits caught it. | For any gate named as evidence that a NEW production line/reader is correct, SENTINEL-INSTRUMENT that exact line (a file-write or counter — **NOT** a bare `assert`, which `-O` strips, and NOT a print that scrolls past) and confirm the gate's run actually EXECUTES it before crediting the gate. If the sentinel never fires, the gate is vacuous FOR THAT CLAIM regardless of its green/assert-firing status — find (or write) a gate whose call graph reaches the new reader, and mutation-verify it reddens. Separately, when the only catchers are tests that build the input packet DIRECTLY with a surrogate that recomputes the production formula, those tests pin the CONSUMER's threading (which field → which slot — mutate the consumer, they red) but are structurally blind to the PRODUCER/stamp (mutate the stamp, they stay green — the surrogate carries the same wrong value on both sides). **NEVER** read "the named twin is green" as "the named twin exercises the rewired code" — only a fired sentinel + a red stamp-mutation proves the new path has a committed catcher. **Sharpening (NEW private adapter/reader):** when the rewired line is a fresh private helper/accessor with no public surface, the gold-standard "the gate executes the rewired line" proof is a **pytest-plugin sentinel that WRAPS the internal call** (an in-process autouse fixture / `monkeypatch` that increments a counter or appends to a list each time the production reader is entered), asserting the counter > 0 at gate end. This is strictly stronger than a file-write probe: it runs IN the test process on the SAME object the production path constructs, so a green twin that routes around the new line (batched kernel, direct source read, zero-caller per-element method) leaves the counter at 0 and reddens the gate — the routed-around path cannot fake the wrap. |
| 12  | Invariant-functional gate (the measured quantity's invariance group contains the error class) | A gate measures a DERIVED functional `f(K)` of a constructed object `K` — an eigenvalue, a spectrum, a balance sum, a normalised shape — and the mutation class it is credited against lies inside `f`'s **invariance group**, so the gate is blind *exactly*: not sub-floor (Mode 10) but identically-zero error in the measured quantity, at every tolerance, in every regime, under every refinement. ORPHEUS instance (#226 taxonomy step 5b, the homogeneous ``K = A⁻¹F`` carve): the verification plan's teeth row claimed a factor-swap mutation (``F·A⁻¹`` for ``A⁻¹F``) would move k∞ O(1) and red the value gates — but ``A·(A⁻¹F)·A⁻¹ = F·A⁻¹`` (the swapped product is SIMILAR to the true one) and ``eig(Mᵀ) = eig(M)``, so both the swap AND the resolvent-transpose mutations are spectrally invisible: measured ‖Δk‖ = 0.0 EXACTLY while the matrix itself moves O(1) (‖ΔK‖ ≈ 1.46 swap / 1.43 transpose). Every k-level gate — the cross-engine rtol=1e-12 consistency gate AND the structurally-independent SymPy closed-form anchor — was DESIGNED-GREEN on the whole class; tightness and reference-independence are irrelevant when the functional annihilates the error algebraically. The committed catcher is the matrix-level intrinsic gate ``K.as_matrix() ≡ np.linalg.solve(A, F)``. Anti-patterns #3 (1G k is flux-shape independent — a degenerate functional) and #8 (particle balance holds by telescoping — the balance functional is invariant under per-ordinate errors that cancel in the sum) are prior instances of the same lens. | At GATE-DESIGN time — **before any mutation is run** — enumerate the measured functional's invariance group (spectra: similarity conjugation + transpose; balance/telescoping sums: any per-term error cancelling in the sum; normalised shapes: global scaling; trace/determinant: similarity) and intersect it with the threat model's mutation classes. A mutation inside the stabiliser is DESIGNED-GREEN — no tolerance tightening, mesh refinement, or regime change can ever expose it through that functional (contrast Mode 10, where an isolating regime MAY exist); the remedy is a gate on a functional OUTSIDE the stabiliser — canonically the constructed OBJECT itself (a matrix/operator-level intrinsic gate against an independently-posed reference: **pin the OBJECT, not just its spectrum**) — then mutation-verify the object-level gate reds (the Mode-10 discipline; the analytic stabiliser check is what the empirical mutation cannot give you when the plan's EXPECTED outcome is itself wrong, which is precisely how the step-5b overclaim arose). Live application (#276 A4, as MEASURED at the phase sweep): the daggered eigenvalue has ``eig(Kᵀ) = eig(K)`` BY CONSTRUCTION, so "k* matches k" gates the posing identity while carrying ZERO vector information — and it is EXACTLY blind to the factor-ORDER/similarity family (``eig(Mᵀ) = eig(M)``: A4's own P1.4 reference encoded precisely this wrong law, its rank-1 dominant eigenvector degenerating to ν̂Σf with zero A-physics, and every k row was designed-green on it; the structurally-independent SN daggered solve's VECTOR row caught it). Do NOT overstate the stabiliser: leaf-transpose DROPS (F†→F, S†→S, L†→L) are NOT inside it — transposing ONE factor is not a similarity of the pencil, and k measurably moves (1.488→0.153 under F†=F on the 4G ∞ fixture) — so k-equality rows ARE legitimate mutation teeth for drops in regimes with the asserted visibility preconditions (χ∦νΣf, asymmetric SigS, spatial structure); the committed catchers for factor-order and flux-shape remain eigenvector/bilinear functionals (the adjoint spectrum row, biorthogonality, duality pairings), never the shared spectrum alone. **NEVER** credit a value gate as a mutation class's catcher without the stabiliser check or a red mutation — a green value gate is not a caught mutation. **A SECOND closure mechanism (ERR-067): repair the METRIC, not (only) the gate.** When the invariance is an *artefact of a wrong metric* — a zero-weight / degenerate state block places the error class *inside* the measured functional's stabiliser — the remedy can be to make the metric non-degenerate (SPD) so the error class EXITS the stabiliser and the SAME functional catches it. Available exactly when the metric itself was the bug, and then **the correctness fix and the Mode-12 closure are one and the same** (ERR-067: the SN ψ½ block Hilbert metric ``G_sd ≡ 0`` put the seed rows in ``ker G``, so G-adjoint reciprocity ``⟨Aψ,φ⟩_G = ⟨ψ,A†φ⟩_G`` was identically blind to any seed-row error — and worse, ``A.H`` was a *wrong adjoint* for any nonzero seed; installing the SPD ``G_sd = V_cell`` both fixes the adjoint AND makes reciprocity a real catcher). Closing this way carries its OWN trap: the blindness is exact precisely on the input the old gate fed (here a *zero* seed), so the closure gate MUST (a) exercise the previously-nulled input (a *nonzero* seed) AND (b) carry a **control leg** — the unmutated honest baseline still holds ``< tol`` — else a still-broken baseline (also off-tolerance on the new input) *mimics* "caught" and the closure is itself false. Full derivation: `docs/theory/foundations/infinite_medium.rst` (``spectral-invisibility`` anchor); worked case: [issue226_spectral_invisibility.md](scripts/issue226_spectral_invisibility.md); metric-repair case: ERR-067 (`error_catalog.md`) + `tests/sn/operators/test_starting_direction_metric.py::test_derive_gsd_and_close_mode12`. |

**The mechanism is non-tokenizer.** Modes 1–6 are observable signatures
of sub-word tokenizer co-location (see `reference.md` §2). Mode 8 is a
toolchain/runtime-mode failure: the assertion is real in source but
compiled out by the interpreter flag, so the bug is unobservable at run
time regardless of how good the assertion is. Mode 11 is a
call-graph/coverage failure: the assertion is real AND fires, but the
gate's execution simply never reaches the rewired production line, so
the gate measures the unchanged siblings rather than the change (the
defense is sentinel-instrumenting the named line, not trusting that a
green twin exercised it). Mode 12 is an algebraic-invariance failure:
the gate executes everything and asserts on the intended quantity —
but that quantity is a functional whose invariance group contains the
error class, so the error is annihilated before any assertion sees it
(the defense is analytic and design-time: enumerate the functional's
stabiliser, then gate the OBJECT itself). Mode 7 is
human cognitive bias: the simplest trial function that satisfies the
BCs is also the most error-resistant to derive, and so wins by default
even when stronger trials would stress more of the solver. AI agents
using SymPy have no derivation cost, so this defense is no longer
needed — and yet the bias survives because the existing canonical MMS
examples are isotropic-by-construction (Lewis & Miller §6.4 ansatz set,
NIST MMS reference set). **Always** pair an isotropic ansatz with an
angularly non-trivial companion in curvilinear / Pℓ contexts.

This mode does not get its own ERR-NNN entry until a real solver bug
is shown to have hidden behind an MMS ansatz in production. The
abstract risk is documented here (skill table); a concrete instance
becomes an ERR entry per the "Log every caught bug" directive below.

**Documentation-layer companion to the Mode-10 sub-floor defense.**
When a term is *exercised-but-unconstrained* (Mode 10 — no isolating
regime exists, so the verification is a STRUCTURAL pair, not a
value-improvement leg), the honest-scope note has a second home beyond
the test: a prophylactic `.. warning::` block IN the theory/RST page
itself, pre-empting the future over-claim a fresh reader would
naturally make from the code (e.g. "do NOT read this as 'recovers 2nd
order at the boundary' — the boundary face-slope sits below the bulk
O(h²) floor and is verified only structurally, NOT by an
error-improvement leg"). The warning is a doc-authoring move, not a
test: it inoculates the *next* session's claim taxonomy at the exact
page where the over-claim would otherwise be minted. **Always** pair
the Mode-10 honest-scope note (test side) with a prophylactic
`.. warning::` (doc side) when the verification is structural-only —
the test pins the math, the warning pins the language.

---

## Hierarchical claim taxonomy — verify the lower layers first

Claims are layered. Each layer adds dependencies. Verify lower layers
before higher ones, and match evidence to the _claim's_ layer.

```
              ┌────────────────────────────────┐
              │  Eigenvalue claim              │  depends on eigenvalue solver
              │  (k_eff, k_inf)                │  + flux shape + discretisation
              └────────────────────────────────┘
                            ↑ depends on
              ┌────────────────────────────────┐
              │  Flux-shape claim              │  depends on the discrete model
              │  (ψ(r,μ,E), φ(r))              │  + boundary conditions
              └────────────────────────────────┘
                            ↑ depends on
              ┌────────────────────────────────┐
              │  Convergence-order claim       │  pure math; lowest dependency
              │  (O(h^p), MMS slope)           │  verifies parts AND whole
              └────────────────────────────────┘
```

Layer reclassifications to apply when reading a claim:

- **Convergence-order results are _math claims_, NOT _solver claims_.**
  They prove the discretisation is consistent — nothing about the
  solved value being correct. MMS lives at this layer.
- **Flux-shape results are _model claims_, NOT _eigenvalue claims_.**
  They depend on the equation and the BC, not on the eigenvalue
  iteration. MMS reaches this layer when the source is structurally
  independent of the code's primitives.
- **Eigenvalue results are _solver claims_.** They bring the iteration
  scheme, normalisation, and convergence test into consideration. MMS
  does NOT directly reach this layer — k-eigenvalue verification needs
  an analytical eigenvalue (homogeneous infinite medium, transfer
  matrix) or a structurally-independent semi-analytical reference.

---

## CRITICAL: The three pillars of verification

Every verification reference is one of three kinds. Each kind proves a
different thing. **NEVER** name a reference vaguely as "analytical" —
**instead** identify which pillar it belongs to, because each pillar
has a different evidence boundary.

### The duality at the centre

Two questions reveal the pillar split:

- **"Given an equation, find the solution"** → **closed-form** analytical solutions
- **"Given a solution, find the equation source"** → **MMS** (Method of Manufactured Solutions)

When neither question closes algebraically:

- **"Reduce the equation to a single integral, evaluate to arbitrary precision"** → **semi-analytical**

Closed-form and MMS are both *analytical* (exact by construction).
Semi-analytical is *exact via arbitrary-precision numerics*. The
distinction matters when judging what claims a pillar can support.

### What each pillar proves

| Pillar          | Convergence-order                  | Flux-shape            | Eigenvalue            | When it applies                                         |
| --------------- | :--------------------------------: | :-------------------: | :-------------------: | ------------------------------------------------------- |
| Closed-form     | ✓ (against exact)                  | ✓ (under assumptions) | ✓ (exact)             | Limited regimes (homogeneous, simple geometry)          |
| **MMS**         | ✓ (great flexibility)              | ✓ (any imposed shape) | **✗** (source-driven) | Any operator that admits a non-vanishing trial solution |
| Semi-analytical | ✓ (against arb-precision integral) | ✓                     | ✓                     | Hard cases with no closed form                          |

**MMS does NOT prove eigenvalues.** This is mechanical, not a
limitation. By construction MMS is a *source-driven* problem — you
imposed the solution, derived the source that makes it true, and the
eigenvalue is whatever k you started with. There is no eigenvalue
information in MMS to verify against. **NEVER** make eigenvalue claims
on the basis of MMS evidence — **instead** match the eigenvalue claim
to a closed-form or semi-analytical reference.

### MMS operational rules

- **Trial solution MUST NOT vanish under derivatives.** Trigonometric
  and exponential functions are the canonical candidates. Polynomials
  vanish at finite derivative order and produce trivial residuals.
- **Trial solution MUST be non-trivial at boundaries** to verify
  boundary-condition handling. A solution that vanishes at the
  boundary by construction tests nothing about the BC.
- **Trial solution MUST stress-test the numerical method, NOT
  minimise source complexity.** Human MMS designs trend toward simple
  sources because hand-derivation of Q^ext is error-prone. AI agents
  using SymPy have no such constraint — the source is derived
  programmatically. **NEVER** pick "the simplest trig that satisfies
  the BCs" when stronger trial functions exist — **instead** pick
  ψ_chosen for stress-test value: high-frequency oscillation, mixed
  scales, near-singular boundary behaviour, non-trivial group-coupling
  for multi-group transport. The simplification heuristic that
  protects humans from arithmetic errors does not serve verification.
  See reference.md §4.3 for the mechanism.
- **Manufactured source MUST be structurally independent of the
  code's primitives.** If the source is generated by the same
  numerical primitives the code uses, MMS becomes a tautology.

### Semi-analytical correctness ladder

Semi-analytical correctness rests on a two-step chain:

1. **Integrator correctness.** For `scipy.integrate`, `mpmath.quad`,
   etc., correctness is commonly assumed (well-tested upstream). For
   custom integrators, integrator correctness is itself a
   verification requirement before this pillar applies.
2. **Reduction correctness.** The reduction from equation to single
   integral is the pillar's load-bearing math. If the reduction is
   wrong, the integral is exact for the wrong equation — a reference
   contamination instance (see anti-patterns).

If both steps hold, the integral evaluation gives the solution to
arbitrary precision. The Peierls reference solver in
`orpheus.derivations.continuous.peierls_nystrom` is the canonical
ORPHEUS instance.

### Structural independence — applies across all three pillars

Whichever pillar you use, the chain of trust **MUST** terminate in a
structurally-independent ground. **Procedurally-independent ≠
structurally-independent.** Two derivations that use different code
paths but exercise the same integrand or identity are *procedurally*
independent only. When shipping a new reference, force the cross-check
to come from a different *structural* angle — a kernel check (row-sum,
particle balance) AND a closed-form check (eigenvalue, asymptotic
limit) — **NEVER** two derivations of the same closed form.

### Ancillary references — NEVER pillars

These are NOT pillars; they are ancillary uses of references that
already exist:

- **Independent re-derivation** — a different mathematical path to the
  same closed form. Strong cross-check if the paths are *structurally*
  independent (different identity / different integrand). Weak if only
  procedurally independent.
- **Code-to-code (L4)** — Reserve **exclusively** for cross-implementation
  agreement. **NEVER** proves correctness — both implementations could
  be wrong. Every L4 claim **MUST** name its L0–L3 backing.
- **Monte Carlo** — itself a numerical method that needs verification
  (geometry tracker, free-flight sampler, collision physics, tally
  estimators). Useful as a *consumer* of references; **NEVER** a
  *source* of them. Comparing CP-vs-MC is L4 benchmarking, not
  verification, until MC itself has been verified against an
  analytical / probability-chain reference.

---

## V&V level taxonomy — the ladder

```
VERIFICATION — "Are we solving the equations right?"
  L0  Term verification        hand calc vs code, per term
  L1  Equation verification    analytical solutions, MMS, convergence order
  L2  Integration testing      multi-group + heterogeneous, self-convergence

VALIDATION — "Are we solving the right equations?"
  L3  Validation               experimental data (ICSBEP, IRPhE, SINBAD)

INFORMATIONAL — parallel to the ladder
  L4  Benchmarking             code-to-code — produces zero correctness info

ORTHOGONAL TO THE LADDER
  foundation                   software invariants — no theory-page :label:
                               (data structures, factory outputs, algebraic
                               reduction invariants). Use @pytest.mark.foundation;
                               NEVER carry verifies(...).
```

- **L4 is parallel to the correctness ladder, not part of it.**
  L4 produces information about whether two implementations agree —
  it produces zero information about whether either is correct.
  Every L4 claim **MUST** name its L0–L3 backing.
- **L3 is sequenced, not aspirational.** ICSBEP / IRPhE / SINBAD data
  exists; L3 starts after L1 maturity (when the verification matrix
  has populated, verified entries below it). L3 without L2 is
  accidental agreement.
- **Necessity chain.** L1 without L0 = compensating errors. L2 without
  L1 = masked components. L3 without L2 = accidental agreement. L4
  without L0–L2 = proves nothing.

---

## CRITICAL: Bit-identity vs principled-equivalence

**Bit-identity is an implementation property, not a math property.** A
regression contract that demands `np.array_equal` on numerical outputs
is a strong gate when the implementation is unchanged — you get free
verification by inheritance from a previously-verified reference.
That same gate becomes the WRONG gate when a refactor deliberately
changes the floating-point reduction tree (a wiring through a new
primitive, a vectorization, a measure-based integration replacing a
broadcast-multiply-then-flat-sum). The two implementations compute
the same value in real arithmetic and disagree at IEEE-754 ULP because
addition is not associative.

**MUST** accept a non-bit-exact change ONLY when ALL THREE of the
following hold. Reject if any fails.

1. **The new formulation is principled at every step**, meaning each
   intermediate is a named, inspectable quantity — not "whatever the
   reduction order happened to produce". Per-group integrated
   reaction rate `r_g = ∫ Σ_g φ_g dV` is principled (a reactor-physics
   quantity); the per-cell-per-group product field `V_i Σ_(i,g)
   φ_(i,g)` summed across all axes by `np.sum` is unprincipled (the
   intermediate is a `(N, ng)` array no consumer ever names). Refactors
   that move from unnamed-intermediate to named-intermediate are
   principled even if they cost bit-identity.
2. **The new value is verified against a structurally-independent
   reference.** Old-vs-new ULP-distance is necessary but **NEVER**
   sufficient — proving "the new value is close to the old value"
   does not prove the new value is correct (both could be wrong by
   the same systematic offset). The reference must come from a
   different structural angle: closed-form analytical (e.g. `k_∞ =
   νΣ_f / Σ_a` for homogeneous reflective), higher-precision
   recomputation (`mpmath`, `float128`), MMS, or any of the three
   pillars. If no structurally-independent reference is reachable,
   the change is REJECTED.
3. **The drift is FP-non-associativity, dimensionally explainable.**
   For an iterative solver: drift bounded by `(iteration count) ×
   (condition number) × ULP`. For a single-step computation: drift
   bounded by `(reduction depth) × ULP`. Drift that exceeds these
   bounds signals an algorithmic change masquerading as FP noise —
   investigate.

When all three hold, **MUST** narrow the regression contract for the
specific touched primitive (e.g. relax `np.array_equal` →
`assert_array_almost_equal_nulp(nulp=K)` for the affected outputs);
preserve bit-identity elsewhere. The contract narrows in scope, gains
a documented relaxation justified by the three criteria above, and
stays principled. **NEVER** silently relax the contract without
documenting all three.

**Worked example (issue #169)**: `compute_keff` rewired from
`np.sum(Σ_p · φ · V[:, None])` (single flat reduction over `(N, ng)`,
unnamed intermediate) to `compute_group_production_rate(φ).sum()`
(per-group rate vector intermediate, then sum over groups). The
intermediate IS the per-group production rate — a reactor-physics
diagnostic quantity. Verified against `k_∞ = νΣ_f/Σ_a` for the
homogeneous reflective snapshots (analytical limit), bit-identical
agreement at the cell-averaged-flux test. Drift on heterogeneous
snapshots: ≤ `iteration_count × ULP`, well under the existing
`rtol=1e-12` regression tolerance — no contract relaxation needed in
that case. The principled refactor passed all three criteria.

**Anti-pattern to flag**: an API method whose only purpose is to
reproduce a specific legacy FP reduction tree (e.g. a `mu.total(M)`
verb that exists because `mu(M).sum()` doesn't bit-match
`np.sum(M * V[:, None])`). The legacy FP order is an arbitrary
historical choice; encoding it in the API is reverse-engineering the
abstraction to fit the implementation. Prefer composing the
principled chain `mu(M).sum()` and accepting the FP order it produces.

**Anti-pattern to flag — an OFFLINE-isolated error is NOT
automatically "the floor."** An error measured in isolation (a
component's residual, a per-kernel discrepancy, a matvec
self-consistency round-trip ≈ 0, an offline reconstruction's
truncation) does NOT, by that fact alone, earn the label "dominant
error floor" OR "the improvement this change buys." Internal
self-consistency is necessary but **NEVER** sufficient — a
matvec≡sweep round-trip at 1e-16 proves SI and Krylov solve the SAME
operator, NOT that the operator's fixed point is correct (ERR-061:
every component individually correct, the bug was the FRAME-CONSISTENCY
between two correct components — "O(h²) to the wrong limit is still
O(h²)"). Before crediting an isolated error as the floor or as an
improvement, it MUST survive THREE end-to-end checks: (1) an
**end-to-end swap** — wire the isolated piece into the full solver and
confirm the claimed effect persists in the converged answer (not just
the offline residual); (2) a **term-silent control** — zero / scalar-
ize the term and confirm the converged answer is byte-identical where
the term should not matter (pins the asymmetry); (3) **amplification —
the sharpest disproof** — grow the term (3×, 10×) and confirm the
converged answer gets strictly WORSE against a structurally-independent
reference. If amplifying the term does NOT degrade the end-to-end
result, the term is not the floor and the "improvement" is offline-only
— the claim is REJECTED. AMPLIFY is the strongest single test because a
genuinely-dominant error term cannot be scaled up without the
converged value moving; a sub-floor / inert / compensated term stays
silent under amplification, exposing the false credit.

---

## CRITICAL: 1-group degeneracy — canonical statement

**k = νΣ_f / Σ_a is flux-shape independent.** A 1-group eigenvalue
test cannot detect any error in the spatial, angular, or scattering
operators — the result is a material-property ratio, computable
without solving the transport equation. **Multi-group (≥2G) is MUST
for any verification claim.** This section is the **canonical home**
of the 1-group-degeneracy rule — historically shorthanded "Cardinal
Rule 6" across the codebase, a citation retired 2026-06-21 (CLAUDE.md
has only Cardinal Rules 1–5; this rule lives here, in `vv-principles`,
with anti-pattern #3 as its operational form). `qa/AGENT.md` and
`test-architect/AGENT.md` cite this section.

---

## CRITICAL: Log every caught bug

This skill owns the L0 error catalog (`error_catalog.md` in this
skill directory). Every agent that loads `vv-principles` is bound by
the following directive.

**MUST** log every bug caught during development → `error_catalog.md`
with:

- **ERR-NNN** (next sequential ID)
- **Failure mode** (1–6 from the AI failure modes table)
- **How it hid** — what evidence-class fooled the previous tests
- **Which test catches it** — linked via `@pytest.mark.catches("ERR-NNN")`
- **Lesson** — one sentence

The catalog is a QA publication artifact and the skill's primary
self-improvement vehicle. **NEVER** close a numerical-bug
investigation without an ERR entry. The catalog grows the skill;
gaps in the catalog mean lessons did not propagate.

**A `catches("ERR-NNN")` marker is a COVERAGE CLAIM, not a topic
tag.** **NEVER** attach `catches("ERR-NNN")` to a test on the basis
that it lives in the same area / same module / same equation family
as the bug — **instead** mutation-verify that THIS specific test goes
red when the EXACT documented bug is re-introduced. A test can carry
the marker while being structurally blind to that bug: it may catch a
DIFFERENT failure class in the same code region (e.g. a cell-matrix
assembly pin `A==A` that is blind to a dropped *inflow*-assembly
factor, because the inflow term is not in the matrices it checks).
The blind marker inflates the catalog's per-ERR coverage count with a
non-catcher and creates a false sense that the regression is pinned.
This is L7's level-conflation argument applied to `catches`: the
marker writes a coverage edge the audit trusts, so an unverified
marker is a phantom. **Verification recipe**: for every NEW
`catches(ERR-NNN)`, re-drop the exact bug the ERR entry documents
(the entry names the file + the dropped/flipped factor) and confirm
THIS test — not merely *some* test in the run — fails under the
canonical `-O` invocation. If a different test catches it and this one
stays green, the marker belongs on the OTHER test; drop it here.
(#240 D5b-S2: `test_d2_assembled_matrices_match_symbolic` carried
`catches("ERR-060")` but ERR-060 was the dropped `|μ_axis|` factor in
`assemble_inflow_axis`; the pin checks `assemble_ubld`'s A/M/G/F_out
and PASSED under the |μ_axis| drop — only `test_d2_exact_on_bilinear`
caught it.)

---

## Sign-pattern + magnitude fingerprint diagnostic

Sign-pattern + magnitude scaling form a 2-D fingerprint that pins bug
class before debugger steps. The full fingerprint catalog lives in the
adjacent skill — see
[../numerical-bug-signatures/SKILL.md](../numerical-bug-signatures/SKILL.md).
**Read fingerprints before opening mpmath.**

---

## Pointers

- **Catalogued bugs (ERR-NNN):** `error_catalog.md` in this skill
  directory. Every L0-caught bug carries: failure mode (1–6), how it
  hid, which test catches it, lesson.
- **Worked case studies:** `scripts/` in this skill directory, one
  file per epistemic-failure case. See `scripts/_template.md` to
  add a new one.
- **Adjacent skills:** [`numerical-bug-signatures`](../numerical-bug-signatures/SKILL.md)
  (recognition catalog), [`probe-cascade`](../probe-cascade/SKILL.md)
  (factor isolation), [`nexus-verification`](../nexus-verification/SKILL.md)
  (graph-based coverage audit — invoke its tools during a V&V review).

For the philosophy (structural independence, Oberkampf–Roache frame,
tokenization grounding, reference contamination), read
[reference.md](reference.md).
