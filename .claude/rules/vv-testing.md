---
paths:
  - "tests/**"
  - "tests/_harness/**"
---

# V&V test harness & test-execution standards

(Applies when working under `tests/`. The *why* behind these — the V&V hierarchy, the 6 AI
failure modes, structural independence — lives in the `vv-principles` skill.)

## Canonical test invocation: `python -O -m pytest`

ORPHEUS treats **`python -O -m pytest`** as the canonical/default invocation — the
production path strips `assert`. Plain `pytest` (debug, `__debug__ == True`) is reserved for
tests that exercise the Layer-3 defensive dimensional checks living under `assert`
(`Field._check_partner`).

- Default in scripts/hooks/CI: `python -O -m pytest …`; list the `-O` job FIRST.
- Layer-3 catch tests: `@pytest.mark.skipif(not __debug__, reason=…)`, with a companion
  `@pytest.mark.skipif(__debug__, …)` "strip-verified" test. Any new test depending on a
  production-code `assert` MUST carry this pair.
- Note: pytest's assertion rewriting still works under `-O` for `assert` *inside test files*
  (AST-rewritten); only production-code (`orpheus/`) asserts get stripped.

## Never relax a tolerance to fit an inexact method

A test tolerance is a **contract**. If a test fails because the chosen quadrature/cubature
isn't exact at the required degree, do NOT loosen the assertion — implement or substitute
one that IS exact.

1. First check whether another existing factory IS exact (e.g. SH integration at L=2:
   `Lebedev` order ≥ 2L, `product_mu_phi(n_μ≥L+1, n_φ≥2L+1)` are exact; LS rules are NOT —
   they're optimised for moment integration in transport, not arbitrary SH products).
2. If none exists, implement what's needed rather than weakening the assertion.
3. Investigate WHY a failure exceeds the FP-non-associativity bound — usually it's a
   mathematical issue, not a numerical one. (`vv-principles` L11: cross-checks must be
   structurally independent; prefer a bit-identical unit-vector cross-check over a merely
   procedurally-independent loop-vs-einsum comparison.)

## Tagging & linking (the `tests/_harness/` registry)

`tests/_harness/` carries verification metadata; architecture: `docs/theory/verification/harness.rst`.

- **Tag a level** (pick one — all feed the registry): `@pytest.mark.l0`/`l1`/`l2`/`l3`/
  `foundation`; or `class TestL0Foo:`; or file-level `pytestmark = [...]`; or inherited via
  `case_name=` from a `VerificationCase`. Precedence: explicit > class decorator > class name
  > case inheritance.
- **Link to an equation:** `@pytest.mark.verifies("label")` (matches a `.. math:: :label:`
  in `docs/theory/`); Nexus writes a `tests` edge from the test node to the equation node.
- **Link to a caught bug:** `@pytest.mark.catches("ERR-NNN")` for every `error_catalog.md` entry.

## Trivial execution & audit
- `pytest -m l0` — term-verification; `pytest -m "l1 and not slow"`;
  `pytest -m "verifies('matrix-eigenvalue')"`.
- `python -m tests._harness.audit` — prints the V&V matrix (level × module × equation),
  orphan equations, ERR-NNN coverage. Sphinx regenerates
  `docs/theory/verification/matrix.rst` from the same registry on every build.
