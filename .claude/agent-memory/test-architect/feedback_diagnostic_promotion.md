---
name: Diagnostic-to-permanent-test promotion patterns
description: Patterns and pitfalls when graduating derivations/diagnostics/diag_*.py scripts into tests/derivations/ permanent V&V harness
type: feedback
---

When promoting `derivations/diagnostics/diag_*.py` scripts into
`tests/derivations/`, follow these patterns. They are derived from
the 2026-04-30 four-file batch promotion.

**Why:** Diagnostics carry investigation context (helper imports,
print narratives, hard-coded constants) that must be sanitised
before promotion — otherwise the test inherits brittleness from a
scratchpad. Three failure modes occurred in this batch.

**How to apply:**

1. **Verify the diagnostic actually runs as-is** before promotion.
   The synthesis-style diagnostics frequently import from sibling
   diag files that were deleted in earlier triages — `pytest
   --collect-only` on the source surfaces this immediately. If the
   source imports a dead helper, prefer the public solver API
   (`solve_peierls_1g`, `solve_peierls_mg`) and reproduce the
   numbers locally rather than vendoring orphaned helpers.

2. **Reproduce pin numbers via the public API before writing the
   test.** A 1-line python invocation that prints the regime values
   gives you confidence the test will pass AND lets you choose
   tolerances that absorb quadrature jitter without being so loose
   they hide regressions. For NEGATIVE-regression gates (pin a
   known-bad result), use a ±5% range around the empirical value.

3. **Three classes of foundation promote** based on what the test
   gates:
   - **Software invariant** (per-face vs aggregate parity, branch
     reduction): `@pytest.mark.foundation` + `@pytest.mark.catches("ERR-NNN")`
     where the audit logged an ERR.
   - **Negative-regression gate** (pins a documented WONTFIX
     plateau): `@pytest.mark.foundation`. The test asserts the gap
     stays in the published range; if a future closure improvement
     lands, the test fails — the FAILURE is the signal to update
     the gate.
   - **Math-origin promote** (general identity unrelated to a
     theory equation): `@pytest.mark.foundation`. No
     `verifies(...)` because there's no `:label:` in `docs/theory/`
     to point at.

4. **Suppress documented warnings**, don't ignore them. Phase-4
   sphere/cyl multibounce emits a `UserWarning` at N≥4 by design
   (overshoot pathology). When the test pins the overshoot
   *behaviour* itself, wrap the call in
   `warnings.catch_warnings(); warnings.simplefilter("ignore", UserWarning)`.
   When the test pins the *absence* of warnings (slab geometric
   immunity), use `simplefilter("error", UserWarning)`.

5. **Source delete order**: only `rm` the diag AFTER the new test
   passes in isolation. Never delete before verifying — the
   diagnostic is the only place those pin numbers existed.
