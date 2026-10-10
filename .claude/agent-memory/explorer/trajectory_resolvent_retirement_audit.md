---
name: trajectory-resolvent-retirement-audit
description: P1 step (e) audit of the old trajectory-resolvent family. Origins stay; the gates that string-pin it; the ERR arm that reds; the missing successors
metadata:
  type: project
---

Stamped at HEAD `62468273` (2026-10-10). The full table is `scratch/characteristic_architecture/p1_step_e/explorer/dependency_audit.md`. The write-scope hook refuses `.claude/plans/`, so the orchestrator copies it.

**Shape:**
- **`origins/specular/` is NOT retired.** The spec's B7 and its 61 KEEP SymPy rows keep 18 of the 23 `derive_*` functions; 5 retire with the 26 duplicate rows. This is a clean AST split.
- **0 production importers.** The family owns no sibling module outside its package.

**The gates that red on their own when files are deleted, and that a symbol grep misses:**
- `withdrawal_506_placement.txt`, together with `test_withdrawal.py`'s `assert len(files) == 27`;
- `test_error_catalogue_reconciles.py`: arm 2 (an entry with no catcher) and arm 5 (the catalogue's test-path citations);
- `test_content_identity.py:64`, which imports `ROSTER` from a test file;
- `cross_method` tolerance keys keyed on adapter NAMES, and a case whose only tolerance is the old adapter.

**Why:** step (e) retires the family in one commit. These are the surfaces it must carry.

**How to apply:**
- A test-file deletion owes a grep of the file NAME across `tests/*.txt`, `tests/_harness/` and the error catalogue.
- Read every literal count of files or ids in a placement or withdrawal gate.

Related: [[characteristic-reference-sn-consumers]], [[trajectory-resolvent-family-inventory]].
