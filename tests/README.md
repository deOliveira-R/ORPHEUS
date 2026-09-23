# tests/ — every check of ORPHEUS, one folder per regimen

Everything that tests the code lives under `tests/`, organised by *regimen*: how a suite runs, how its cases depend on each other, and what a result means. What a test is *about* (its V&V level, foundation, the equation it verifies, the bug it catches) is a marker on the test, never a folder.

The canonical statement is the verification principles page, section "Where a case lives — one folder per regimen" (`docs/theory/verification/principles.rst`, label `vv-test-suite-layout`): the regimen table, why the folder is the regimen and the kind a marker, why only the gates stack, why there is no code-to-code (L4) suite, and where a new case goes. Read it before adding a case; this file does not restate it.

- `gates/`: the pass/fail checks, run by pytest on every commit (`testpaths = ["tests/gates"]`; the canonical invocation is `.venv/bin/python -O -m pytest`). Below its top level it mirrors the package: `gates/sn/` gates `orpheus/sn/`.
- `performance/` and `validation/`: the measured regimens. A folder is created when its first case lands, so an absent one means no case yet.
- `_harness/` and `conftest.py`: shared by every regimen (the V&V registry and audit, reference helpers, the pyright ratchet; the collection hooks). The harness contract is `docs/theory/verification/harness.rst`.
