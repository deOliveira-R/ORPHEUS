---
name: ORPHEUS V&V tagging idioms
description: How to tag tests for the tests/_harness registry (verifies, l0/l1, foundation, xfail)
type: feedback
---

ORPHEUS's V&V audit is driven by `tests/_harness/registry.py` which
parses pytest markers. Conventions learned from reading the existing
test suite:

**Why:** Tests must be picked up by `python -m tests._harness.audit`
and by `verification_coverage()` in Nexus. Mis-tagged tests show as
orphans and will surface in every `session_briefing`.

**How to apply:**

- Use module-level `pytestmark = [pytest.mark.verifies("label")]` to
  assign a default Sphinx label to every test in the file. Then
  override per-test with additional `@pytest.mark.verifies(...)` when
  a test also exercises a different equation.
- Level markers: `@pytest.mark.l0` / `.l1` / `.l2` / `.l3` or
  `.foundation` for software-invariant tests (bit-exact recovery,
  frozen-immutability, algebraic consistency). Foundation tests MUST
  NOT carry `verifies(...)` in strict reading of
  `docs/testing/architecture.rst` — BUT the existing code tolerates
  both, and for a test that is BOTH a foundation check AND a
  regression gate on an equation (e.g. rank-1 bit-exact recovery of
  an existing verified closure), double-tagging with the existing
  label (e.g. `peierls-white-bc`) provides the upstream gate that
  breaks noisily when the equation's code changes.
- `@pytest.mark.slow` for tests that take >5 s. Use freely — `pytest
  -m "l1 and not slow"` is the standard fast-gate.
- xfail tests for features not yet implemented: use
  `strict=False` so the test naturally un-xfails when the feature
  lands, rather than silently failing XPASS. Always include a
  `reason=` string that points to the exact API or Sphinx label that
  will unlock the test.
