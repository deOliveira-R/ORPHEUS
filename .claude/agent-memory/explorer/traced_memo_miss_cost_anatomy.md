---
name: traced-memo-miss-cost-anatomy
description: Where a cheap traced-memo miss spends 2.7 s (2026-10-15, b5b5166a) and how to probe the child without editing orpheus
metadata:
  type: project
---

Measured 2026-10-15 at `b5b5166a` (#602). Record of the tree at that date; re-verify before repeating.

- A cheap miss takes 2.70 s. The parts:
  - the manifest, 1.79 s. It splits into `importlib.metadata.version` called once per Distribution code object
    (1245 calls for 6 names, 0.65 s), `_source_text` parsing each of 129 files twice and digesting all 1948 defs
    (0.79 s), and `packages_distributions` (0.22 s);
  - the child's teardown, 0.35 s; `os._exit` measured 4 ms;
  - the imports, 0.47 s;
  - the interpreter start, 0.036 s.
- The recorder is already `sys.monitoring` PY_START. Its hook costs 0 on the compute and 37 ms on the imports.
- A cold validation (the first in a process) costs about 0.8 s, about 6 ms per distinct pinned file. A warm one
  costs about 9 ms for 448 pins.
- The probe technique: point `tm._BOOT` at a scratch wrapper that re-spells the boot at MODULE LEVEL, so no
  wrapper def is recorded, and extend `tm._CHILD_ENVIRONMENT` to pass the probe's variables.
- A forked recording template produced manifests equal to a fresh child's, apart from `_generate_here`.
- Probes, logs and the report are in `scratch/characteristic_architecture/p3/m602/`.

**Why:** #602 cuts the miss cost; the next profile should start from this anatomy, not re-derive it.
**How to apply:** re-run the probes in `scratch/characteristic_architecture/p3/m602/` after any change to the
miss path. See [[traced-memo-client-census]].
