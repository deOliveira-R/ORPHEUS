---
name: pyright-ratchet-harness-invariant
description: PASS-with-nits ruling on the #226 pyright error ratchet (ede47dc); reusable verdicts for harness software-invariant reviews
metadata:
  type: project
---

# Pyright error ratchet (#226, commit `ede47dc`, branch `test/pyright-ratchet`)

A foundation+slow harness test that holds per-module pyright error counts monotone
(must only go down) until the #226 burn-down hits zero. No CI yet → the harness IS
the gate. VERDICT: PASS-with-nits (landable as-is; 2 optional CONCERNs).

**Why:** establishes the reusable "what good looks like" for a harness software-invariant
in ORPHEUS, distinct from the SN operator-algebra carves that dominate the rest of memory.

**How to apply** — the rulings that generalize to the NEXT harness-invariant review:

1. **Dual-direction failure = illegal-slack-state unrepresentable (Pattern 4).** A
   ratchet that fails ONLY on increase admits the slack state (burn down, forget to
   re-baseline, contract silently re-permits regression to the old number). Failing on
   DECREASE too makes `live == baseline` the only green config. This is the load-bearing
   design move; endorse it. Verify both leak-corners: a module burned to 0 (drops out of
   the live Counter → still caught via `live.get(m,0)==0 < b`) and a brand-new module
   (caught via implicit `baseline.get(m,0)==0`).

2. **Version-hint-NOT-skip is the correct vv-Mode-8 call.** Skip-on-tool-version-mismatch
   is the inert-tripwire false-green Mode 8 warns against (newer-pyright machines silently
   stop ratcheting). Surface drift as an appended triage string on a REAL failure. Same
   ruling applies to any external-tool-versioned gate.

3. **`pytest.fail()` not bare `assert`** — function call, fires under canonical `python -O`.
   ⚠ **[CORRECTED 2026-09-20, `[M]`]** the stated REASON is over-broad: a bare `assert` in a
   COLLECTED test module does fire under `-O`, because pytest's rewriter replaces the `Assert`
   AST node with an explicit `raise` before compile, so `-O` has nothing to strip (measured:
   `assert 1 == 0` -> `1 failed` under both `python` and `python -O`; the
   `_warn_about_missing_assertion` warning `-O` prints concerns NON-rewritten modules). The
   preference still holds where it bites — a contract asserted in a helper, a conftest, a
   plugin or in `orpheus/` itself is stripped. ⟹ do NOT flag a bare `assert` in a test module
   on this ground; DO flag one in production or shared-helper code (`coding-standards`).
   `foundation` marker, NEVER `verifies(...)` (software invariant ≠ physics equation). Direct
   sibling precedent = `tests/test_layer_imports.py` (import-linter, same artifact-invariant
   class). Both `foundation`+`slow` registered in pyproject.toml.

4. **Single-source the counting, not just the data.** Test imports `collect_module_counts`
   / `read_baseline` from the harness module; `--update` CLI writes the SAME
   `collect_module_counts` output the test reads back → baseline definition can't drift from
   live definition. `write_baseline` normalizes (sorted keys + `total=sum`) so every baseline
   diff is deterministic. This is the producer-side normalization done right.

**The 2 optional CONCERNs (local code-shape, NOT architecture):**

- **C1 — `regressions`/`improvements` are two spellings of one comparison** keyed on
  opposite sides (`live.items()` vs `baseline.items()`) so each catches the absent-on-
  -other-side corner. Correct today but load-bearing-by-accident: a future THIRD category
  added to one comprehension not the other reintroduces the keying asymmetry. Twin-fold smell
  (institutional pattern #2) at small scale. Collapse to ONE keyed pass over
  `baseline.keys() | live.keys()`, partition by `sign(n-b)`.

- **C2 — `--update` command string transcribed at 3 live sites** (write_baseline comment,
  improvements failure msg, + docstrings). Pattern-7 convention-re-applied. A module rename
  stalens the printed regenerate-command. Single-source the 2 NON-doc copies via a module
  constant; leave the `.rst` as an independent reader restatement.

**Retirement trigger to track (non-blocking, recommend a #226 note):** when `total==0` the
dual-direction ratchet degenerates to a plain `live=={}` assertion and the `--update` /
baseline-file scaffolding becomes dead weight (institutional anti-pattern #11). The endgame
elegant move is flip to a zero-error assertion + delete baseline+CLI. Scaffolding is JUSTIFIED
now; just shouldn't outlive the burn-down silently.

**`_ROOT_BUCKET="(root)"`** is a CORRECT unrepresentable-collision sentinel (parens aren't
legal in a package name) — not a latent collision, good instinct. Bucketing predicate
`len(rel.parts) > 1` is right (a file directly in `orpheus/` has exactly one part).
