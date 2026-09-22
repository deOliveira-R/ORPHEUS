---
harness:
  kind: rule
  budget_tokens: 1600
  paths:
    - "tests/**"
    - "tests/_harness/**"
---

# V&V test harness & test-execution standards

The *why* behind these — the V&V
hierarchy, the six AI failure modes, structural independence — is the
`vv-principles` skill.

## Canonical test invocation: `python -O -m pytest`

ORPHEUS treats **`python -O -m pytest`** as the canonical invocation: the
production path strips `assert`. Plain `pytest` (`__debug__ == True`) is
reserved for tests that exercise the Layer-3 defensive dimensional checks
living under `assert` (`Field._check_partner`).

- check: scripts, hooks and CI default to `python -O -m pytest …`, and list the
  `-O` job FIRST.
- check: a Layer-3 catch test carries `@pytest.mark.skipif(not __debug__,
  reason=…)` with a companion `@pytest.mark.skipif(__debug__, …)`
  "strip-verified" test; any new test depending on a production-code `assert`
  carries this pair.
- The scope of `-O` (a collected test module keeps its `assert`s; a helper, a
  fixture, a `conftest`, a generator and production code lose them) is
  `coding-standards` § "A bare `assert`"; outside a collected module assert
  with `np.testing.assert_*` or a `raise`.

## Never relax a tolerance to fit an inexact method

A test tolerance is a **contract**. If a test fails because the chosen
quadrature or cubature is not exact at the required degree, do NOT loosen the
assertion; implement or substitute one that IS exact.

1. check: whether another existing factory IS exact (SH integration at L=2:
   `Lebedev` order ≥ 2L and `product_mu_phi(n_μ≥L+1, n_φ≥2L+1)` are exact; LS
   rules are NOT, being optimised for moment integration in transport, not
   arbitrary SH products).
2. If none exists, implement what is needed rather than weakening the
   assertion.
3. Investigate WHY a failure exceeds the FP-non-associativity bound; usually it
   is a mathematical issue, not a numerical one. Cross-checks must be
   structurally independent (`vv-principles` #7 and its structural-independence
   section): prefer a bit-identical unit-vector cross-check over a merely
   procedurally-independent loop-vs-einsum comparison.
4. A KNOWN approximation gap (a method inexact by design) is no reason to
   loosen either: document the gap in the docstring AND the test, pin the
   residual against a structurally-independent reference proving it is
   approximation and not bug, and file the closure issue (ERR-036 and ERR-038
   both shipped `tol=5e-2`).

## Tagging & linking (the `tests/_harness/` registry)

`tests/_harness/` carries the verification metadata; the architecture is
`docs/theory/verification/harness.rst`.

- **Tag a level**, one way (all feed the registry): `@pytest.mark.l0`/`l1`/
  `l2`/`l3`/`foundation`; or `class TestL0Foo:`; or a file-level `pytestmark =
  [...]`; or inherited via `case_name=` from a `VerificationCase`. Precedence:
  explicit > class decorator > class name > case inheritance.
- **Link to an equation:** `@pytest.mark.verifies("label")` (a `.. math::
  :label:` in `docs/theory/`); Nexus writes a `tests` edge from the test node
  to the equation node.
- **Link to a caught bug:** `@pytest.mark.catches("ERR-NNN")` for every
  `error_catalog.rst` entry. Either marker is a coverage CLAIM with a shelf
  life, adjudicated per `vv-principles` § "Log every caught bug".

## A type that embodies a mathematical concept ships the test of its defining laws

When a type or an operator embodies a mathematical concept (a cone, a
probability simplex, a multiplier algebra, an affine torsor, a frame), its
defining properties are enumerated and each is one gate that lands WITH the
type, at L0 or L1: closure under the operations, the identities, the
invariants, and one failing-input negative. A usage test or a round-trip
asserts nothing about the concept: an untested law is a claim the code makes
about the mathematics that nothing verifies, and the concept drifts silently.
`[R]` the user, 2026-06-19, on `SpectrumField`'s `Σχ = 1` (issue #257).

- check: for a new math-bearing type, list its laws before its first consumer
  and write one gate per law (a cone: closed under `+` and under `λ ≥ 0`, with
  `σ = 0` its origin; a simplex: `Σχ = 1` at construction plus the refusal; a
  multiplier algebra: `M_f M_g = M_{fg}`, `M_1 = I`, `M_0 = 0`, `M_f.H = M_f`).
- tell: a type whose only tests are its consumers'.

## A hand-built mixture is gated on its own consistency identity

A hand-built verification fixture prints and asserts the mixture's identity
`σ_t == σ_c + σ_f + Σ_to SigS[0][g,:]` before any reference value is trusted:
an inconsistent mixture gives the transport balance and the
production/absorption balance two DIFFERENT answers, and two correct solvers
then report different ones with no bug in either (`[M]` 2026-08, #340 N5: a
brief's "benign pole" was 30 % off because `sig_s` was written `[to, from]`
where the constructor reads `[from, to]`).

- check: assert the identity in the fixture, and its companions in the same
  line: a group with `φ ≡ 0` is a 1-group problem wearing a 2-group shape
  (`vv-principles` #3); a negative `σ_c` is unphysical.
- tell: a brief's reference value on a hand-built mixture, quoted but never
  re-derived.

## Trivial execution & audit

- `pytest -m l0` — term verification; `pytest -m "l1 and not slow"`;
  `pytest -m "verifies('matrix-eigenvalue')"`.
- `python -m tests._harness.audit` prints the V&V matrix (level × module ×
  equation), the orphan equations and the ERR-NNN coverage; Sphinx regenerates
  `docs/theory/verification/matrix.rst` from the same registry on every build.
