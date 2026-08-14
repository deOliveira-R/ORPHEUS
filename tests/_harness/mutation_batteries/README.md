# Hand-written mutation batteries

A **battery** is a small set of hand-authored, *named*, hazard-specific
mutations — each one the plausible transcription of a real bug — installed by
monkeypatching production **in process** and run against a chosen slice of the
suite. It answers one question: **do these gates actually have teeth, and which
gate catches which hazard?**

## ⛔ This is NOT `tests/_mutation/`

The two directories exist for different questions and their recipes are not
interchangeable.

| | `tests/_mutation/` | here |
|---|---|---|
| mechanism | **cosmic-ray**: automatic AST mutation of a whole module | hand-authored named mutations |
| question | *what fraction* of mutants does the sentinel set kill? | *does gate X redden for hazard N?* |
| output | a mutation **score** (a percentage) | a per-hazard, per-gate colour |
| reverts by | `git checkout -- <module>` (its README, step 5) | **nothing** — it never writes to disk |

⚠ **The `git checkout` in the cosmic-ray recipe is exactly the practice these
batteries exist to avoid.** A `git checkout` / `restore` / `stash` / `clean` on
a path you hold uncommitted edits in destroys them irrecoverably
(`.claude/rules/process-discipline.md`). Every battery here monkeypatches in
process precisely so that reverting is *not an operation* — when the process
exits, the tree is untouched because it was never touched.

## Why these live in the repo at all

`[M]` The failure this directory prevents has happened **twice**, and the
evidence is still on disk:

- The B3 plan's acceptance gate names a **31-mutation harness** ("the sweep
  still catches 30/31"). It lived only in a job-scratch directory and
  **evaporated with the session**. `b3_2_boundary.py` is its reconstruction,
  and it carries 10 mutations, not 31.
- `tests/sn/operators/__pycache__/` still contains
  `conftest_mutate_kernel.*.pyo` and `conftest_mutate_pr.*.pyo` — but
  `git log --all -- 'tests/sn/operators/conftest_mutate*'` is **empty**. Two
  more hand-written plugins were written into `tests/`, used, and lost,
  leaving only bytecode.

A battery that cannot be re-run is not evidence. It is a claim.

## Running one

Batteries are pytest **plugins**, loaded with `-p`, and selected by an
environment variable. Canonical invocation is `-O` and **SERIAL** (xdist is
unstable on this interpreter):

```sh
# CONTROL leg — nothing patched. Everything must be GREEN.
.venv/bin/python -O -m pytest tests/sn/operators \
    -p no:randomly -p tests._harness.mutation_batteries.b3_2_boundary -q

# One mutation. Expect RED, and read WHICH gates went red.
ORPHEUS_B32=N1 .venv/bin/python -O -m pytest tests/sn/operators \
    -p no:randomly -p tests._harness.mutation_batteries.b3_2_boundary -q -rf
```

**Always run the CONTROL leg first.** A battery whose control is red is
measuring its own breakage, and `vv-principles` #17 is explicit that the
harness lies before the code does — and lies in the *safe-looking* direction.

## The inventory

| module | what it mutates | mutations | env var |
|---|---|---|---|
| `b3_2_boundary.py` | `SNBoundaryOperator._reflect_trace`, the SN boundary realizer, `TraceRestrictionOperator.to_local` | 10 (`N1`–`N5`, `N7`–`N9`, `M1`, `M2`) | `ORPHEUS_B32` |

`[M]` 2026-08-14 at `0f5ca91c`, `b3_2_boundary.py` against `tests/sn/operators`:
CONTROL → **1167 passed, 1 skipped, 5 xfailed** in 38.9 s;
`ORPHEUS_B32=N1` → **66 failed**, 1101 passed in 138.7 s.

⚠ A full 10-leg sweep costs **~20–25 min**: a mutated leg runs ~140 s rather
than ~39 s, because a broken boundary makes downstream fixed-source solves
iterate to their budget caps. It also runs all 1174 tests in
`tests/sn/operators`, so its reds are *not* scoped to the gates being
justified — use `-rf` and read the names.

## Not yet migrated

These belong here and are still in `scratch/`. Each needs its **own CONTROL-leg
verification at HEAD** before it moves — a battery that no longer runs should be
repaired or retired, not relocated:

- `scratch/b3_4a_mutations.py`
- `scratch/mutate_angular_closure_seam.py`
- `scratch/mutate_convergence_contract.py`

## Wanted: make these run in CI

Today a battery is a manual run, so nothing stops it rotting. The honest fix
lifts each gate's assertion into a callable helper, so a fast meta-test can
monkeypatch and call it directly — seconds rather than 140 s a leg.

⛔ **Do not implement that as pytest-inside-pytest.** `vv-principles` #17
records the measured failure: monkeypatching the parent while a child process
re-imports a clean module reads **GREEN for every mutation**, which is a
mutation harness failing in the direction that looks like success.
