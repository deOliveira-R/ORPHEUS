---
name: characteristic-memo-stage-census
description: #405 P3 memo boundaries in the characteristic reference — the dry-run stage census technique, and why the group block is the seam (keyability, PointValue re-assembly, skeleton pins)
metadata:
  type: project
---

Stamped 2026-10-14 at `ad573eb1`. The full report is `scratch/characteristic_architecture/p3/speed/memo_boundaries.md`; the tools are in `.../speed/memo/`. These are claims about the tree on that date: re-verify before repeating them.

**The technique (reusable).** A dry-run census answers "which stage invocations recur with identical arguments, per test" in about 40 min instead of 2.6 h.
- Run `python -O -m pytest --assert=plain -p <plugin>`: the test's own asserts are stripped.
- The plugin stubs the expensive stage with a right-shaped synthetic and records a digest of exactly the content the stage reads.
- It no-ops `numpy.testing`, makes the pencil answers synthetic, and patches the ladder judges (`geometric_error` and kin) per test module in `pytest_collection_modifyitems`.
- `Δ = warm1 − dry` per test splits stage time from test-side (mpmath) time.

**Why:** v1 ended 88 of 109 system tests at the first synthetic `NoFundamentalMode`, and 3 of 6 ladder tests at the ladder judge. Each stub layer was found by the census's own failures.

**How to apply:** list the outcomes per file after every census run. A non-passed test is a lower bound on its invocations.

**The shape found:**
- The memo's value is cross-RUN. 0 of the top 30 tests share a block, and within one run the tests' `lru_cache`s already share everything else (175 repeats, about 69 s).
- `evaluate(PointValue)` entries pin `LineRule.transport`: each memoised point read re-assembles every group block in its fresh process. `[[traced-memo-client-census]]` has the earlier family.
- `PanelBasis`, `ConcentricPartition` and `RigidMotion` are `eq=False` and not `ContentIdentity`, so `encode_exact` refuses them as keys. `Walls` and the arrays key fine; the `GroupTransport` payload round-trips exactly.
- A solve entry pins 129 module skeletons, all 13 characteristic modules among them, because the package `__init__` imports everything. Only body edits are stage-local.
- The gates behind `bypass()` hold 2726 s of the 7800 s stage time. `bypass()` is all-or-nothing.
