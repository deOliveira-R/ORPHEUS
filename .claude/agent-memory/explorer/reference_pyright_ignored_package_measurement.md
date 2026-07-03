---
name: pyright-ignored-package-measurement
description: How to get a TRUE pyright error count for a package listed in pyproject [tool.pyright] ignore — and the 2-artifact discount rule for the /tmp-config bypass
metadata:
  type: reference
---

Measuring a pyright-IGNORED package (the `[tool.pyright] ignore` list in
`pyproject.toml` suppresses the files' OWN diagnostics, so a plain
`npx pyright orpheus/<pkg>` reports 0 regardless of truth):

- **Bypass that works:** pass the package's files as CLI args together with
  `--project /tmp/<dir>/pyrightconfig.json` where the config has NO ignore
  list, absolute `venvPath` = repo root, `venv: ".venv"`, `pythonVersion`
  matching pyproject. CLI file args escape the config's include; an
  ignore-free config suppresses nothing.
- **Discount rule:** this bypass ALWAYS adds spurious `reportMissingImports`
  for `orpheus.*` imports — pyright cannot follow the PEP-660 editable
  install from an external project root (an absolute-path `include` entry is
  rejected outright: "not relative"; `executionEnvironments.extraPaths`
  doesn't apply to files outside the config root either). Those imports
  resolve fine in repo-root runs. Subtract them; only non-import diagnostics
  are the package's true count.
- Verified 2026-07-03 on `orpheus/diffusion`: bypass reports 3, true count
  is 1 (`solver.py` scipy `LinearOperator(matvec=...)` reportCallIssue —
  scipy stubs don't model the ctor kwargs).

How to apply: the pyright ratchet campaign's Level 2 lifts ignore entries
one at a time — any "how dirty is ignored package X?" pre-lift scoping
question uses this method + discount rule instead of trusting the plain
run's 0 or misreporting artifacts as real errors. (The alternative —
deleting the ignore entry from pyproject.toml — mutates a tracked file;
explorer is read-only, so the /tmp bypass is the explorer-safe route.)
