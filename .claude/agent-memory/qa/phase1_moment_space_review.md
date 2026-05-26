---
name: phase1-moment-space-review
description: QA review of Phase 1 (P1.1–P1.7) moment-space + layering work — ERR-039 fix landing on refactor/moment-space-and-layering
metadata:
  type: project
---

# Phase 1 review — moment-space + layering (ERR-039 endpoint)

Review scope: `git diff refactor/sn-operator-algebra..HEAD` on
`refactor/moment-space-and-layering` (17 files, +2517/−563 LoC).
Reviewer applies `vv-principles` + `coding-elegance` skills.

Severity tags: 🔴 BLOCKER · 🟡 SMELL · 🟢 NOTE.

---

## Category A — Architectural / organizational concerns

### A1 🟢 `.codomain` and `.range` are intentionally redundant — but the redundancy is honest

Files: `orpheus/numerics/projection.py:393–426`.

`MomentProjection.codomain` and `.range` both return
`SphericalHarmonicSpace.from_L(self.L)`. The plan
(`moment_space_and_layering_plan.md:559–561`) schedules
`P3.5 — range → codomain` as a framework-wide rename across
`operator.py` and every operator implementation. Until then, the
generic `_AdjointOperator.apply` machinery in `operator.py:490, 494,
511–516` reads `inner.range` to find the codomain metric; for the
new `MomentProjection.H` chain to work TODAY without modifying the
framework, `MomentProjection` must expose a `.range` alias.

Verdict: this is principled transitional state, NOT
`coding-elegance` anti-pattern 1 (two ways to spell the same thing).
The two names exist because they belong to different abstractions
(operator algebra book uses `codomain`; ORPHEUS framework still
reads `range`); P3.5 retires `.range` codebase-wide. **The dual
exposure is correctly scoped to the transition window.**

The one improvement worth flagging: the `.range` docstring at line
417 already names P3.5 and explains the alias. Good — but the
`__all__`/Sphinx surface lists neither attribute as deprecated. Add
`.. deprecated:: P3.5 use .codomain` to the `.range` docstring at
line 418 so users see the gate explicitly. Without it, a contributor
landing P3.5 may forget to delete the alias here. Minor; not
load-bearing.

### A2 🟡 `from_Y` shim is still load-bearing — schedule its retirement explicitly

Files: `orpheus/numerics/projection.py:616–632`,
`orpheus/sn/scattering.py:639`,
`tests/numerics/test_projection_operators.py:151,161,171,201,284,415,469`,
`tests/numerics/test_spherical_harmonic_space.py:380`.

Production callers of `HarmonicMomentReconstruction.from_Y`:
- `orpheus/sn/scattering.py:639` — `build_aniso_source` builds R from
  the cached Y table.

Test callers: 9 in `test_projection_operators.py` (these are the
legacy V&V suite — many tests construct R directly via `from_Y(Y)`
without first building a `SphericalHarmonicSpace`).

Per `feedback_aggressive_retirement` in the user's memory:
"superseded code = noise that obscures signal. Every refactor must
retire its predecessor pattern; deprecation shims only one merge
cycle." The shim is currently TWO merge cycles into its life (Wave 0
introduction → Phase 1 fix). The plan does not explicitly schedule
its retirement; CC.3 retires the inline `(2*l+1)` literal in
`_build_rhs_*`, but **`from_Y` itself remains a way to build R
without going through SphericalHarmonicSpace** — preserving the
exact "two ways to spell the same thing" anti-pattern at the
construction surface.

Recommendation: file a P3.x follow-up to retire `from_Y`:
1. Rewire `scattering.py:639` to
   `HarmonicMomentReconstruction.from_spherical_harmonic_space(
   self._sh_space, Y)` where `self._sh_space` caches
   `SphericalHarmonicSpace.from_L(L)` once at `__post_init__`. This
   is also the cure for A3 below.
2. Rewire the 9 test callers to construct via `from_spherical_harmonic_space`.
3. Delete `from_Y`.

If P3.2/P3.5 already encompass this, the plan does not say so —
the `from_Y` retirement is currently invisible in the roadmap.
Either explicitly add it to a phase, or accept the indefinite-shim
debt with a tracked issue.

### A3 🔴 `MomentProjection.domain` / `.codomain` / `.range` rebuild a fresh space + a fresh `FunctionSpace` on every property access — and the hot adjoint path hits them per call

Files: `orpheus/numerics/projection.py:393, 416, 428`,
`orpheus/numerics/operator.py:511–516`.

The `_AdjointOperator.apply` body reads:
```python
inner_range = getattr(self.inner, "range", None)   # line 511
inner_domain = getattr(self.inner, "domain", None) # line 512
```
on **every adjoint application**. For `M.H.apply(c)`:
- `M.range` → builds a fresh `SphericalHarmonicSpace` (allocates a
  new `(L+1, 2L+1)` weight array via `_padded_metric_tensor`).
- `M.domain` → builds a fresh `FunctionSpace(name="angular_ordinate",
  shape=(N,), inner_product_weights=self.weights)`.

For a Krylov iteration that calls `M.H.apply` once per matvec ×
~50 iterations × ~100 source iterations, that's ~5000 allocations of
the padded metric (and ~5000 `FunctionSpace` constructions, which
involve `__post_init__` validation on `SphericalHarmonicSpace`).
Each allocation is `O(L²)` memory. For `L=7` that's 64 floats per
space × 5000 → 320 KB of throwaway allocations per outer iteration
just for the space objects, plus the GC pressure.

**Why I'm calling this BLOCKER**: the `coding-elegance` skill
Pattern 7 (normalise at the definition site) demands that the
SphericalHarmonicSpace metric be computed once, not per access.
There's no algorithmic justification for rebuilding it — `L` is
frozen at `__post_init__` time. The same applies to `.domain`:
`self.weights` never changes, so the `FunctionSpace` is also
constant.

**Fix (one line each)**:
```python
@cached_property
def codomain(self) -> "SphericalHarmonicSpace":
    from orpheus.numerics.spaces.spherical_harmonic_space import SphericalHarmonicSpace
    return SphericalHarmonicSpace.from_L(self.L)

@cached_property
def range(self) -> "SphericalHarmonicSpace":
    return self.codomain

@cached_property
def domain(self) -> FunctionSpace:
    return FunctionSpace(
        name="angular_ordinate",
        shape=(self.weights.shape[0],),
        inner_product_weights=self.weights,
    )
```

`@cached_property` IS available on frozen dataclasses (it stores on
`__dict__`; frozen blocks `__setattr__` on declared fields only,
but `@cached_property` writes through `__dict__` directly — see
CPython issue 87634 / fix in 3.12+). If the project still supports
3.11, the pattern needs `functools.cached_property` plus a
non-frozen dataclass (which is a larger surgery) — OR a manual
`__post_init__` that uses `object.__setattr__` to install the
spaces once. Either way: **do not rebuild per call**.

Same fix applies to `SphericalHarmonicSpace.basis` (`spaces/spherical_harmonic_space.py:200–206`), which builds a fresh `SphericalHarmonicBasis(L=self.L)` on every access. Since `SphericalHarmonicBasis` is also frozen and the per-ell arrays are `@cached_property` on it, the construction itself is cheap — but it's still a Python-level allocation in a hot path.

This issue ALONE warrants a P1.8 amendment commit before Phase 3 begins. Phase 3 layers more callers onto `M.H` / `M.codomain`; landing the per-access rebuild as the "interface" calcifies it.

### A4 🟡 `_broadcast_for_leading_axes` is a module-level free function consumed only by `_AdjointOperator.apply`

File: `orpheus/numerics/operator.py:197–237`.

Per `coding-elegance` Pattern 6 (defer abstraction until two
instances), and Pattern 5 (build the right primitive), the
function's API ("reshape an inner-product weight tensor to broadcast
on a data tensor with metric-axes-leading") is generic enough that
it has a future home — but TODAY there is exactly one caller.

Three reasonable placements:

1. **Method on `_AdjointOperator`** (private) — keeps the helper
   where it's used; future callers can promote it.
2. **Method on `FunctionSpace`** — `space.broadcast_weights_to(ndim)`
   reads as "ask the space for its weights pre-shaped for this
   tensor", which IS a natural domain operation. Co-locates with
   `inner_product` / `norm` already on `FunctionSpace`.
3. **Stay where it is** — module-level free function.

I lean toward (2) because the operation is conceptually "the space
broadcasts its metric to act on a tensor whose leading axes are
this space's shape". The `FunctionSpace.inner_product` body at
`space.py:165–167` already does an implicit broadcast (the
right-aligned numpy default); the explicit leading-axis pad is a
legitimate `FunctionSpace`-level concern.

But this is plainly a "right answer is opinion" choice — none of
the placements is wrong. **Flag as smell; not a blocker; address if
the helper acquires a second caller**.

### A5 🟢 `SphericalHarmonicBasis.synthesize(coefficients, directions)` recomputes Y from directions internally

File: `orpheus/numerics/basis/spherical_harmonic_basis.py:265–305`.

The signature takes `(coefficients, directions)`. Inside, line 304
calls `Y = self.evaluate(directions)` — recomputing the Y table.
This duplicates work when the caller already has Y in hand (which
is the case in the SH-space tests at `test_R_equals_2l_plus_1_times_S0`
line 181: `basis.synthesize(c * factor[:, None], measure.nodes)` —
the test already has Y at line 172).

For test code this is fine. For production code (currently zero
production callers — `synthesize` is documented as "where the
ERR-039 naked S₀ lives") it would matter. Two options:

1. **Add an overload** `synthesize_from_Y(coefficients, Y)` that
   skips the recomputation. Production paths use the Y-overload;
   convenience tests use `synthesize`.
2. **Accept the recomputation** — `_evaluate_real_sh` is cheap
   (~O(N·L²)) and the naked-S₀ operation is itself O(N·L²); the
   recomputation is a constant-factor doubling, not asymptotic.

Recommendation: defer until a production caller appears. The
current callers are V&V tests where clarity beats performance.

### A6 🟡 `ReconstructionOperator` ABC has only one concrete subclass — but the deferral rule has a load-bearing exception

File: `orpheus/numerics/projection.py:207–231`.

Per `feedback_unify_after_two_instances`: "Build each new
geometry/topology standalone first; only unify after ≥2 working
instances."

`ReconstructionOperator` has exactly one subclass:
`HarmonicMomentReconstruction`. The plan justifies the early ABC
via Grand Report v3 §5.7 ("sibling-of-ProjectionOperator
structure"). Is this enough?

Counter-argument **for** the ABC: `ProjectionOperator` itself
also had only one concrete subclass when it shipped, but the
`(Galerkin, Petrov-Galerkin)` discipline split makes the ABC family
load-bearing for type-tagging. The ReconstructionOperator ABC
mirrors that structure precisely.

Counter-counter-argument: there is NO concrete
`Petrov-Galerkin Reconstruction` yet (the §5.7 structure pairs each
projection discipline with a reconstruction). The current
`HarmonicMomentReconstruction(ReconstructionOperator)` inheritance
adds zero invariant beyond the `capabilities = {CAP_APPLY}` default
that `LinearOperatorMixin` would have given it.

**Verdict**: the ABC is on the edge of acceptable. It has zero
behavioural payload (no methods beyond the abstract `apply`), it
mirrors a planned-but-not-shipped split, and the user's memory
explicitly says "accuracy is non-negotiable" with respect to the
unify-after-two rule. Strictly applied, this ABC should have been
deferred until a `PetrovGalerkinReconstruction` was needed.

In practice: the cost is small (one class with one abstract method),
the benefit (type-tagging readability) is real, and the §5.7
reference grounds it. **Smell, not blocker**. Worth recording the
rule-departure explicitly in a CC.x entry in the plan so future
reviewers know the tradeoff was made consciously.

### A7 🟢 P1.7 retirement scope is principled, NOT overreach

Files: `orpheus/sn/solver.py:838–862` (the retirement note),
`tests/sn/test_fixed_source_g1.py:307–360`.

The user asked: was deleting `_build_rhs_{spherical,cylindrical}`
overreach beyond the moment-space scope? P1.7 carves all three.

Strict-moment-space scope: only `_build_rhs_cartesian` carried the
inline `(2*l+1)` literal. Therefore strictly, only Cartesian needed
to retire under the §P1.3 "exactly one place" claim.

Aggressive-retirement scope: all three were dead code (zero
production callers, verified by
`test_no_legacy_eq_map_or_decoder_in_g1_path`). Per
`feedback_aggressive_retirement`: "superseded code = noise that
obscures signal."

The Phase 1 fix correctly applied the looser rule. The dead-code
test sentinel at lines 333–342 (`assert not hasattr`) is the
strongest form of "never called" — better than the runtime spy.
**Principled per the user's own rule. Note, not flag.**

One thing to verify: the cleanup note at
`solver.py:843–862` includes a multi-line comment block explaining
the retirement. Make sure the equivalent docs page (the index theory
page if there is one) doesn't still reference these functions.
Spot-check: `docs/theory/spherical_harmonics.rst` was edited in
P1.6 — confirm via Sphinx build, but the docs probably never named
the private `_build_rhs_*` helpers.

### A8 🟡 `MomentProjection.from_measure` and `from_spherical_harmonic_space` are two construction paths that could collapse to one

Files: `orpheus/numerics/projection.py:333–391`.

`from_measure(measure, L)` and
`from_spherical_harmonic_space(space, *, weights, Y)` are
distinct classmethods. The first is the "I have a measure" path;
the second is the "I already built the Y table" path.

The second exists because the SN scattering operator caches Y on the
`ScatteringOperator` instance and wants to reuse it (avoid double
SH evaluation in the hot scattering source builder). That's a
legitimate need.

But the design is mildly awkward: the natural cross-product is
`(measure | space) × (precomputed Y | recompute Y)`, four
combinations. Two ship; two don't (`from_measure_with_Y` and
`from_space_recomputing_Y`). The latter pair is implementable as
trivial dispatches:

```python
@classmethod
def from_space(cls, space, *, measure):
    Y = space.basis.evaluate(measure.nodes)
    return cls(weights=measure.weights, Y=Y, L=space.L)
```

If the hot-cache requirement is the only justification for the
two-path API, factor it into a single constructor that takes
`measure` and an optional `Y` (cached or `None`), defaulting to
"recompute":

```python
@classmethod
def from_measure(cls, measure, L, *, Y=None):
    Y = Y if Y is not None else SphericalHarmonicBasis(L=L).evaluate(measure.nodes)
    return cls(weights=measure.weights, Y=Y, L=L)
```

One constructor, one path. **Smell** — the surface is
multiplying without payoff. Address as P1.8 amendment if `from_Y`
also retires (A2); otherwise it's defensible.

---

## Category B — V&V coverage concerns

### B1 🔴 `test_from_spherical_harmonic_space_roundtrip` is an L4-style code-to-code check, NOT verification

File: `tests/numerics/test_spherical_harmonic_space.py:357–387`.

The test asserts
`np.array_equal(M_new.apply(psi), M_legacy.apply(psi))` — two
ORPHEUS code paths agreeing on a numerical output. Per
`vv-principles` anti-pattern 1: "NEVER claim verification on the
basis of L4 agreement alone — instead require an L0–L2 evidence
chain pointing at a structurally-independent reference."

**What the test actually proves**: the new API constructor path
produces a `MomentProjection` whose `(weights, Y, L)` fields are
bit-identical to the legacy `(weights, Y, L)` constructor. So the
two constructors install the same internal state, and `apply` is the
same einsum on the same fields.

**What it does NOT prove**: that either path is correct. Both
`M_new.apply` and `M_legacy.apply` could be wrong by the same
factor — the test cannot distinguish.

Mitigations already in the file:

- The four mathematical-identity tests (`test_basis_mass_matrix_against_lebedev`,
  `test_R_equals_2l_plus_1_times_S0`,
  `test_pi_R_is_4pi_identity_on_band_limited`,
  `test_H_equals_g_C_times_S0`) DO pin the production semantics
  against structurally independent references. So the L0–L2 backing
  exists in the file.

So the roundtrip test is not load-bearing for ERR-039 verification;
it's a regression guard against the API migration. **The
classification needs to be honest** — the test is not L0/L1, it's a
foundation-level software invariant ("the new constructor preserves
field assignment"). Drop the `@pytest.mark.l0` marker (it currently
carries `@pytest.mark.l0`) and use `@pytest.mark.foundation` per
`CLAUDE.md`'s tagging convention. The test stays; the LABEL is
wrong.

This is BLOCKER-severity because it actively misleads the V&V
matrix `python -m tests._harness.audit` will print — "ERR-039 has
12 L0/L1 tests" will include this one, inflating the count.

### B2 🟢 The unit-vector cross-check in `test_R_equals_2l_plus_1_times_S0` IS structurally independent

File: `tests/numerics/test_spherical_harmonic_space.py:191–202`.

The previous (un-shipped) per-ℓ Python-loop check accumulated
`(2*l+1) * Y[:, l, m_off] * c[l, m_off]` and compared to
`R.apply(c)` for random `c`. That had FP drift because the
einsum's reduction order differs from the Python loop.

The new check sets `c = e_{l₀, m_off₀}` (single non-zero entry).
Then `R.apply(e) = Y[:, l₀, m_off₀] * (2*l₀+1)` exactly — no
accumulation, no FP-non-associativity, bit-identical (`np.array_equal`).

Per `lessons-L11` (structural-independence-via-elimination-of-FP):
the unit vector eliminates the reduction step entirely. The einsum
collapses to a single multiplication per ordinate. **The check is
bit-exact by structure, not by tolerance**. This is the right
fingerprint — the test is structurally independent of the
particular reduction tree.

(Plus, the loop over `(ell_0, m_off_0)` covers every `(ℓ, m)`
combination, so any sign/index error in the `(2ℓ+1)` weighting
would be caught at exactly one of the (L+1)² entries.)

**Excellent test**. Note it; don't change it.

### B3 🟢 The 5 API/type tests dropping `@pytest.mark.verifies` is the correct call

File: `tests/numerics/test_spherical_harmonic_space.py:329–354,
357–387, 443–472`.

Tests that check "the codomain is a SphericalHarmonicSpace", "the
roundtrip is bit-identical", "equality is by (name, shape)" — these
are software-construction invariants, not mathematical identities.
They have no equation label in `docs/theory/`.

The proper marker is `@pytest.mark.foundation` per `CLAUDE.md`:
"software invariants — no theory-page :label: (data structures,
factory outputs, algebraic reduction invariants)". Currently they
carry `@pytest.mark.l0` (or `l1` for the equality test).

Recommendation: rewire to `@pytest.mark.foundation`. They will
still run with `pytest -m foundation` and won't pollute the L0/L1
verification count.

This is the same root issue as B1 — the file conflates
"verification level" with "test severity". The foundation marker
exists precisely for these.

### B4 🟢 The `assert not hasattr` retirement test is acceptable

File: `tests/sn/test_fixed_source_g1.py:333–342`.

Concern raised: "what if a future contributor accidentally
re-introduces `_build_rhs_*` with a different name?"

Answer: they can't, by structure. The retirement test pins the
attribute name `_build_rhs_cartesian` etc.; if someone re-introduces
the function under the original name, the test fires. If they
re-introduce it under a DIFFERENT name, the test won't fire — but
**that's fine**, because the failure mode being guarded is
"re-introduction of the inline `(2*l+1)` literal via the same
code path that ERR-039 traveled". A new function with a new name
would face fresh review.

The test is doing what it should: pin the specific retired symbol
by name. The shape of the test ("symbol absent") is the strongest
form of the retirement contract. **No change needed.**

(One micro-improvement: the helpful error message at line
337–340 names "the inline (2*l+1) duplicate the moment-space plan
§P1.3 ('exactly one place') retires". This is the right level of
context for a future contributor.)

### B5 🟢 P1.7 retirement has no principled-equivalence debt

Per `vv-principles` §"Bit-identity vs principled-equivalence", the
three criteria apply when a refactor moves from one principled
implementation to another. P1.7 is a **deletion** of dead code —
there is no new implementation to satisfy criterion 1 (principled)
or criterion 2 (verified against an independent reference). The
"no callers" sentinel (test_no_legacy_eq_map_or_decoder_in_g1_path)
is itself the verification: the production path doesn't traverse
this code.

The only residual concern is whether someone DOES build a Cartesian
Krylov path in the future, expecting `_build_rhs_cartesian` to
exist. If so, they'd write a new RHS builder against the typed
`AngularFlux + (L+C).apply` algebra (per the comment block at
solver.py:1294–1300), not resurrect a packed-1D helper. **The
retirement is principled because the new path is architecturally
incompatible with the retired helper**, not because the new path
reproduces the old result.

### B6 🟡 The retired `assert_galerkin_idempotency` was the bug — but is its retirement an ERR entry?

The user's question (C2 below) — the broken `assert_galerkin_idempotency`
asserted `Π R c = c` instead of the correct `Π R c = 4π · c` under
the no-prefactor SH convention. The plan §1.5 CC.5 records: "method
shipped broken — sole caller fed a non-orthogonal Y so the wrong
invariant FAILS as a contract test."

This IS a Mode #6 (convention drift) failure: the **convention** of
the SH basis (no-4π/(2ℓ+1) prefactor) and the **method's identity
claim** (`Π R c = c`) drifted at the method-definition site. The
genuine identity is `Π R c = 4π · c`, established at
`test_pi_R_is_4pi_identity_on_band_limited`.

Logging it as a new ERR-NNN is the right call per the `vv-principles`
directive "Log every caught bug". See C2 for the ERR draft.

---

## Category C — Lessons / catalog updates needed

### C1 🟡 ERR-039 should be re-statused — fix shipped, but tag remains valid

File: `.claude/skills/vv-principles/error_catalog.md:2872–2944`.

Status currently: **"CAUGHT IN QA REVIEW 2026-05-10, fixed in same
commit as introduction (Wave 0 / C0.5 of SN performance plan)."**

The status is true for the Wave 0 fix (which made
`apply_transpose` return the bare `S_0(c)` — the "naked synthesis").
But Phase 1 has REWORKED that fix: post-P1.4, `apply_transpose`
carries `w_n` (it's now the **representation transpose** Π^T, not
the naked synthesis); the Hilbert adjoint Π* = g_C·S₀ is reached
via `.H` through the generic adjoint machinery.

The original ERR-039 test (`TestApplyTransposeIsWWeightedAdjoint`)
referenced in the catalog at line 2943 has been RENAMED to
`TestApplyTransposeIsRepresentationTranspose` (in
`tests/numerics/test_projection_operators.py:236`) and joined by
`TestHilbertAdjointViaGenericMachinery` plus the new full file
`tests/numerics/test_spherical_harmonic_space.py`.

**Action**: append a Phase-1 follow-up to the ERR-039 entry:

> **2026-05-26 Phase 1 update (refactor/moment-space-and-layering):**
> The Wave 0 fix made apply_transpose return the naked S₀(c) — which
> Phase 1 audit revealed was ALSO mislabeled (it's actually the
> representation transpose Π^T, which carries w_n, not the naked
> S₀). Phase 1 P1.4 split the three operators:
>
> - `M.apply_transpose` → Π^T = w_n · S₀ (representation transpose)
> - `M.H` → Π* = g_C · S₀ (Hilbert adjoint, via metric-aware adjoint machinery)
> - `R.apply` → R = (2ℓ+1) · S₀ (addition-theorem reconstruction)
>
> All three differ by diagonal multiplications; conflating any pair
> hides bugs. The endpoint is `test_T_carries_w_n_and_H_carries_g_C`
> at `tests/numerics/test_spherical_harmonic_space.py:393–435`,
> which pins both adjoint identities against the production code.

Then leave the status as "FIXED". Do NOT open as a new ERR — this
is the deepening of the same convention-drift failure mode.

### C2 🔴 The broken `assert_galerkin_idempotency` warrants a new ERR-NNN entry

Per `vv-principles` directive "Log every caught bug": the
`assert_galerkin_idempotency` method shipped at Wave 0 with the
wrong invariant. It claimed to check Galerkin idempotency
`Π R c = c`, but under the no-prefactor SH convention the genuine
identity is `Π R c = 4π · c`. The method would have falsely-failed
on every correct production pairing.

The plan §1.5 CC.5 acknowledges this and notes "sole caller fed a
non-orthogonal Y so the wrong invariant FAILS as a contract test."
That description is itself the smoking gun: the method was so
clearly wrong that even its only test was constructed to make the
wrong assertion produce the expected failure. It was effectively
self-canceling, then deleted.

**Draft ERR entry** (the QA agent should land this in
`error_catalog.md` before completing the review per the
self-improvement trigger):

```
## ERR-048 — assert_galerkin_idempotency asserted Π R = I instead of Π R = 4π · I (no-prefactor SH convention)

**Status:** CAUGHT IN PHASE 1 P1.6 (2026-05-XX, moment-space + layering plan).
Method deleted; sole caller deleted alongside.

**Failure mode:** #6 (convention drift) — the method's docstring
named the discipline ("Galerkin idempotency Π R c = c") but the
SH basis convention installs a 4π factor that the addition-theorem
reconstruction carries (R = (2ℓ+1) · S₀, and the discrete Gram
diagonal is 4π/(2ℓ+1), so the composition Π R = 4π · I, not I).

**Date discovered:** 2026-05-XX Phase 1 audit during ERR-039 endpoint
work.
**Module:** `orpheus.numerics.projection.GalerkinProjection.assert_galerkin_idempotency`.

**Mechanism:** The method was added during Wave 0 alongside the
ProjectionOperator ABC. The implementer drafted the assertion from
the abstract Galerkin idempotency statement "Π R = I on the
coefficient space", forgetting that the convention-bearing SH basis
carries an implicit metric (the 4π/(2ℓ+1) Gram diagonal) which
modifies the discrete identity by a factor of 4π.

**How it hid:** The method was never called against the production
(Π, R) pair. The sole call site (`tests/numerics/test_projection_operators.py:368–381`,
since deleted) deliberately built a NON-orthogonal Y matrix so the
method would raise — testing that the method **detects violations**,
not that it **describes the correct invariant**. Net effect: the
method shipped with a wrong-identity error and the test silently
agreed because no one verified what "correct" should look like.

**How caught:** Phase 1 P1.5 designed a test suite for the ERR-039
endpoint, including `test_pi_R_is_4pi_identity_on_band_limited`.
The test discovered that Π R = 4π · I, not I — directly contradicting
the assert_galerkin_idempotency docstring. Audit traced back to the
sole caller's deliberately-broken Y construction; the method's wrong
identity was confirmed.

**Fix:** Delete the method and its sole caller (P1.6). The genuine
identity is now pinned by `test_pi_R_is_4pi_identity_on_band_limited`
(tests/numerics/test_spherical_harmonic_space.py:208–232) and the
sibling test in test_projection_operators.py:194–231.

**Lesson:** When a verification method ships with NO production caller,
its test cannot validate the method against the production semantics
— at best it validates the method against an artificially-broken case.
The "method that asserts a property" pattern requires AT LEAST one
test that exercises the method against a CORRECT instance (where it
must NOT raise) and AT LEAST one that exercises against an INCORRECT
instance (where it MUST raise). Single-test-against-broken-instance
is a discipline failure that hides convention drift indefinitely.
This generalizes: every contract-validation method needs both a
positive and a negative test, or its claim cannot be trusted.

**Test reference:** `tests/numerics/test_spherical_harmonic_space.py::test_pi_R_is_4pi_identity_on_band_limited`.
Tagged `@pytest.mark.l1` and `@pytest.mark.catches("ERR-048")`.
```

The catalog entry numbering: I see ERR-040 through ERR-047 are
already taken (Wave 3 BC trace-law); ERR-048 is the next free slot.
Confirm against the live catalog before assigning.

### C3 🟢 Project-lessons addition — "every contract-validation method needs both a positive and a negative test"

This is the load-bearing lesson from C2. It generalizes ERR-035
(specialised-before-general) and ERR-048 (validation-method-without-
correct-instance): both are failures of the contract surface to
encode what "correct" looks like.

If the user agrees, this is worth a new `feedback_validation_needs_positive_test.md`
entry in `.claude/agent-memory/qa/lessons.md` (or the project-level
equivalent). One sentence:

> **A validation method (`assert_X`, `check_X`, `verify_X`) MUST be
> tested against at least one CORRECT instance where it must NOT
> raise. Tests-against-broken-only is a discipline failure that
> hides identity drift indefinitely.**

Not load-bearing for Phase 1 closeout, but the lesson is genuinely
new. Add or skip per discretion.

### C4 🟢 The "no tolerance relaxation" feedback is already user-memory; project-lessons need is unclear

The user has the feedback in user-memory already
("feedback_unify_after_two_instances" / "feedback_aggressive_retirement"
imply it but don't say it directly). The Phase 1 work USED a unit
vector cross-check specifically to avoid relaxing tolerance — that's
exactly the spirit of the rule.

If a project-level `.claude/agent-memory/qa/lessons.md` does not
have it, add a brief entry: "When a comparison drifts at FP-non-
associativity (e.g. einsum vs Python loop), DO NOT relax tolerance.
Instead, restructure the test to eliminate the FP gap (unit vector,
single-term input, bit-exact construction). See
`test_R_equals_2l_plus_1_times_S0` for the canonical pattern."

This is a sharper version of `lessons-L11` and worth duplicating
into qa/lessons.md if it's not there. Worth checking; not blocking.

---

## Verdict

**Phase 1 is NOT ready to land as-is.** Two BLOCKER-severity items
must be addressed before Phase 3:

1. **A3** — `MomentProjection.{domain, codomain, range}` rebuild
   `SphericalHarmonicSpace` / `FunctionSpace` per access. The hot
   adjoint path (`M.H.apply`) hits both per call. Fix via
   `@cached_property` (or `__post_init__` + `object.__setattr__` for
   pre-3.12 compat). One commit; ~5 minutes of work.

2. **B1** — `test_from_spherical_harmonic_space_roundtrip` carries
   `@pytest.mark.l0` but it's a code-to-code check (L4-style),
   not L0 verification. Rewire to `@pytest.mark.foundation` per
   `CLAUDE.md`. Same patch can do B3 (the 4 other API/type tests in
   the file that should also be foundation, not l0/l1).

These two land as a single **P1.8 amendment commit** ("perf+test-tag
fixes; phase 1 closeout"). Estimated cost: 20 minutes.

Two SMELL-severity items are worth landing in P1.8 if cheap, else
filing as Phase 3 follow-ups:

3. **A2** — `from_Y` shim retirement schedule (file or close in P3.x).
4. **A8** — `from_measure`/`from_spherical_harmonic_space`
   collapse to one constructor.

One CRITICAL catalog item:

5. **C2** — Add ERR-048 entry to `error_catalog.md` (the broken
   `assert_galerkin_idempotency`). This is REQUIRED per
   `vv-principles` "Log every caught bug" — Phase 1 caught it, the
   catalog has not absorbed the lesson yet. **Do this before
   Phase 3 even with P1.8 deferred** — the catalog update is independent of code.

Lower priority (notes only):

- A4 (`_broadcast_for_leading_axes` placement) — opinion call;
  defer until second caller.
- A5 (synthesize duplicates Y compute) — fine for test-only consumer.
- A6 (`ReconstructionOperator` ABC has one subclass) — defensible
  per Grand Report v3 §5.7; record the rule-departure if the
  Phase 1 closeout edits the plan.
- A7 (P1.7 retirement scope) — principled, no action.
- B2, B4, B5, B6 — note for the V&V record; no action.
- C1, C3, C4 — catalog/lessons housekeeping; do at convenience.

**Phase 1's mathematical content (the ERR-039 endpoint) is solid.**
The five identity tests at
`tests/numerics/test_spherical_harmonic_space.py:84–270` pin the
production semantics against structurally-independent references
(Lebedev semi-analytical Gram diagonal; closed-form algebraic
identities; the unit-vector bit-exact cross-check). The R, Π^T, Π*
trichotomy is now correctly separated; the metric lives in exactly
one place (`SphericalHarmonicBasis.metric_per_ell`); the
addition-theorem factor lives in exactly one place
(`SphericalHarmonicBasis.addition_theorem_factor`); and the
two prose warnings collapsed to one Sphinx note. This is the
right shape.

The issues above are about CONSTRUCTION cost (A3 — perf) and TEST
LABELING discipline (B1, B3) — not about the math. Phase 1 should
land with the P1.8 amendment; the foundation is sound.
