# Phase 1 self-review — moment-space + layering plan (ERR-039 fix)

Date: 2026-05-26
Branch: refactor/moment-space-and-layering (commits 1ab6233..f54295c)
Reviewer: main-agent self-review (paired with `qa` agent independent review)

This is the implementer's perspective on what Phase 1 shipped, what
could have been done better, and what the plan-as-executed reveals
about the original plan-as-drafted. Companion to the QA agent's
independent review at `.claude/agent-memory/qa/phase1_moment_space_review.md`.

---

## What Phase 1 actually shipped

8 commits, 17 files, +2517/-563 LoC:

| Commit | Subject | Net effect |
|---|---|---|
| 1ab6233 | docs(plan): P0 reconnaissance landed | Plan refined with 8 reconnaissance findings + resolved 4 open decisions |
| 0eb9cf3 | feat: SphericalHarmonicBasis | New numerics/basis/ package; relocated evaluate_real_sh body |
| f27ba47 | feat: SphericalHarmonicSpace | New numerics/spaces/ package; canonical g_C metric carrier |
| 6ba9102 | refactor: MomentProjection + ReconstructionOperator | Renamed HarmonicMomentProjection → MomentProjection; new ABC sibling |
| 36159f5 | refactor: .T = w_n · S_0; .H = g_C · S_0 | The ERR-039 root fix; broadcast helper generalisation |
| 28d9275 | test: V&V suite (13 tests) | 5 plan identities + 4 API tests + 4 retired/replaced |
| be92b92 | docs: collapse warnings + retire assert_galerkin_idempotency | Sphinx labels added; ABC docstring fixed; dead method gone |
| f54295c | refactor: retire dead _build_rhs_* | 176 lines of dead noise deleted; inline (2ℓ+1) duplicate gone |

ERR-039 is closed across the codebase. The (2ℓ+1) literal lives in
ONE production place. The four operators (S_0, Π^T, Π*, R) are
separately typed with executable identities.

---

## What surprised me

1. **`_build_rhs_*` were dead code already.** The test
   `TestNoLegacyMachineryInCallPath` had already retired them in
   practice; my P1.7 work was just making it official. Expected
   complex migration, got clean deletion.
2. **LS quadratures are NOT exact for L≥2 SH integration.** I
   assumed S_8 would integrate degree-4 SH products. LS_8 gives
   24% diagonal error. They are optimised for moment integration
   in the SN transport equation, not arbitrary `⟨Y_l, Y_l'⟩`
   products. Saved as memory `feedback_no_tolerance_relaxation`.
3. **`numerics/iteration.py` already L1-clean** (CC.7 from
   reconnaissance). P3.4 Problem/Solver split is greenfield, not
   a carve — much smaller scope than I'd budgeted.
4. **The original ERR-039 fix was wrong too.** Round 1 (2026-05-10)
   replaced `R.apply` with bare `S_0` inside `apply_transpose` — but
   still mislabeled it the Hilbert adjoint. Round 2 (this work) is
   the actual structural fix: `.T = w_n · S_0` (representation
   transpose), `.H = g_C · S_0` (Hilbert adjoint via generic
   metric-aware machinery).

---

## What I'd organize differently in hindsight

### O1. The `range` / `codomain` dual-name was a transition cost

P3.5 renames `range → codomain` framework-wide. During Phase 1, I
needed `MomentProjection.codomain` for the new test pin AND `range`
for the existing framework (`_AdjointOperator.apply` reads `range`).
I added both, with `range` aliasing `codomain`. **Cost:** every
operator now has to think about both names during the P3.5 transition.

Better: would have been to do P3.5 first (mechanical rename across
operator.py + all consumers), then build Phase 1 on top of a single
canonical attribute name. Cost: P3.5 has to land first (sequencing
flip). Estimated additional commit: 1-2 medium.

Recommendation: **defer the `range` alias removal to P3.5**; document
in the plan that `range` is the legacy framework name and any new
operator should expose `codomain` first with `range` as alias until
P3.5 completes.

### O2. `MomentProjection.domain` should be `@cached_property`

Current implementation creates a new `FunctionSpace` on every property
access. Cost: small ~µs per access, accumulates in hot Krylov paths.
Fix: add `@cached_property` to `domain`, `codomain`, `range`.

This is a smell `coding-elegance` Pattern 3 (named intermediates)
catches: the FunctionSpace is a named intermediate that should be
materialised once. Easy P1.8 fix.

### O3. `HarmonicMomentReconstruction.from_Y` could retire

After P1.3, `from_Y(Y)` is a back-compat shim that delegates to
`from_spherical_harmonic_space(SphericalHarmonicSpace.from_L(L), Y)`.
The only remaining production caller is
`sn/scattering.py:build_aniso_source` (lines 627-650). If that
caller migrates to `from_spherical_harmonic_space`, `from_Y` can
retire.

Per `feedback_aggressive_retirement`, the shim should retire one
merge cycle after Phase 1 lands. Easy P1.8 fix; one consumer to
rewire.

### O4. The `synthesize(c, directions)` signature re-evaluates Y

Each call to `SphericalHarmonicBasis.synthesize(c, nodes)` builds
the full `Y(N, L+1, 2L+1)` table via `self.evaluate(nodes)`. For
callers that already hold Y (the production case in
`build_aniso_source`), this is wasted work.

Better signature:
```python
def synthesize(self, coefficients, *, Y=None, directions=None) -> NDArray:
    if Y is None:
        if directions is None:
            raise ValueError("either Y or directions required")
        Y = self.evaluate(directions)
    return np.einsum("nlm,lm...->n...", Y, coefficients)
```

This is `coding-elegance` Pattern 5 (build the right primitive):
the consumer can bring its precomputed Y. P1.8 fix.

### O5. `ReconstructionOperator` ABC has only one subclass

Per `feedback_unify_after_two_instances`, the ABC was premature.
But the Grand Report §5.7 sibling-of-ProjectionOperator structure
justifies it as architectural anchor for the future PN solver's
moment-space reconstruction. The cost is small (1 abstract class,
1 abstract method); the benefit is the type-level signal that
`HarmonicMomentReconstruction` is a `ReconstructionOperator`.

**Verdict: keep**, but note in the plan: "ABC ships eagerly because
§5.7 explicitly anticipates ≥2 instances (addition-theorem
reconstruction + future PN method-space reconstruction)."

### O6. 5 API/type tests dropped `@verifies` markers in P1.6

`test_moment_projection_codomain_is_spherical_harmonic_space`,
`test_from_spherical_harmonic_space_roundtrip`,
`test_spherical_harmonic_space_equality_by_name_shape`, and the
mass-matrix parametrised tests don't reference Sphinx equation
labels (they're API/type contracts, not math identities).

Per `vv-principles` §"V&V level taxonomy", API contract tests
should be tagged `@pytest.mark.foundation` (software invariants
without an equation label). Currently they carry only `@l0`/`@l1`.

**Verdict: should add `@pytest.mark.foundation` to the 5 API tests**
in a follow-up. Easy P1.8 fix.

### O7. P1.7 retirement scope expanded beyond strict P1.7

P1.7 plan: retire the inline `(2*l+1)` literal in
`_build_rhs_cartesian`. I deleted ALL THREE
`_build_rhs_{cartesian,spherical,cylindrical}` functions, citing
`feedback_aggressive_retirement` (superseded code is noise) AND the
test `TestNoLegacyMachineryInCallPath` having already declared them
dead.

**Was this scope creep?** Strictly: yes. The plan said cartesian
only. Defensible: the other two were equally dead, the test pinned
their absence, retirement-as-a-group is cleaner than one-at-a-time.

**Verdict: principled, document in plan as a deliberate scope
expansion** — the cost of half-retiring is a future PR that does
the other half with the same justification.

---

## What the original plan got right

1. **P1.0 reconnaissance demand.** The explorer's consumer inventory
   (CC.1-CC.8) shaped Phase 1 in load-bearing ways: `GroupFlux`
   dropped, `IsotropicSource`/`PerOrdinateSource` added to P3.3
   list, `KrylovAcceleration` flagged, `numerics/iteration.py` L1
   status confirmed. Without it, multiple Phase 1 commits would have
   chased dead ends.
2. **Anti-recommendation 1-6.** Each one prevented a structural
   mistake during implementation — the explicit "do NOT" list is
   `lessons-L13` (briefs name existing types) applied to a plan
   level.
3. **The §1 four-operator decomposition.** The mathematical insight
   that `R / Π* / Π^T / S_0` are all `S_0` post-multiplied by a
   diagonal was the conceptual key. With that framing, the
   implementation is mechanical; without it, every commit would have
   been an ad-hoc patch.
4. **Sequencing P1.1 → P1.2 → P1.3 → P1.4 → P1.5 → P1.6 → P1.7.**
   Each step ships independently testable behaviour; tests at every
   step gate the next; no big-bang.

---

## What the plan needs updated for Phase 2/3

### Plan updates from Phase 1 execution

1. **P1.0 — mark COMPLETE** (currently says ✅ — confirm).
2. **All P1.x — mark COMPLETE** with commit SHA references.
3. **Phase 2 — explicit deferral note**, citing
   `feedback_unify_after_two_instances`:
   - P2.1 (DualSpace, Space.dual()) — would generalise the
     `_broadcast_for_leading_axes` helper I added inline in P1.4
     into a Space method. Defer until a second use of leading-axis
     metric broadcasting arises.
   - P2.2 (TensorProduct / DirectSum) — defer per plan's own note
     ("sequence when consumers are scheduled"); issues #172, #173
     are the consumer triggers.
4. **Phase 3 — re-scope P3.4**: per CC.7, this is greenfield (no
   import-graph carve), so much smaller than originally framed.
5. **Phase 3 — add P3.0.5**: the `range → codomain` rename (P3.5)
   should sequence EARLIER, perhaps after P3.1 (linter) before
   P3.2 (numerics reorg). Reason: the rename touches every operator;
   doing it first means the reorg moves files that already use the
   new name.

### Issues to file (or close)

1. **ERR-039 catalog entry** needs update: Round 1 fix (2026-05-10)
   was incomplete; Round 2 (this work) is the structural fix. The
   test references need updating to the new class names. The entry
   should add a "Round 2 fix" section.
2. **No new GitHub issues** needed from Phase 1 execution — all
   surprises were absorbed into the plan or feedback memories. P1.8
   close-out items (O1-O7 above) could be a single GitHub issue OR
   bundled into the Phase 3 work.
3. **Sphinx orphan-label markers** (5 in `test_spherical_harmonic_space.py`
   after my edit, plus a few pre-existing in other test files) —
   could file a tracked issue "audit @verifies markers vs Sphinx
   labels" but it's not Phase 1 blocker.

---

## Verdict

Phase 1 ships clean. ERR-039 is closed across the codebase. The
seven Os above are polish, not corrections — none of them invalidate
the Phase 1 contract.

**Recommendation for the user's "review + organize" question:**

- Mark Phase 1 COMPLETE in the plan with commit SHA references.
- Add a P1.8 close-out commit for O1 (`@cached_property` on
  domain/codomain), O2 (retire `from_Y`), O3 (improve `synthesize`
  signature), O4 (`@pytest.mark.foundation` on API tests). Estimated:
  1 small commit.
- Update ERR-039 catalog entry (Round 2 fix + new test refs).
- Defer Phase 2 entirely (per `feedback_unify_after_two_instances`
  + plan's own deferral language).
- Re-scope Phase 3 per CC.7 (greenfield P3.4) and re-sequence P3.5
  before P3.2.
- Compact context and start a fresh session for Phase 3 (it's a
  many-file architectural overhaul that benefits from clean
  context).
