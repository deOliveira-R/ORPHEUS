---
name: sood_registry Cache (Phase B4)
description: Persistent solver-output cache for Sood-family benchmark cases — feature/peierls-greens-cylinder Phase B4, parallel with B1/B2/B3.
type: project
---

# sood_registry Phase B4 — solver-output cache

**Branch:** `feature/peierls-greens-cylinder`
**Parent SHA:** `f391f6f` (Phase A closeout)
**Commits (this phase):**

- `562f37b feat(sood_registry): persistent solver-output cache for Sood cases`
- `e75a36a test(sood_registry): foundation gates for solver-output cache`
- `9e8fbb9 docs(sood_registry): cache subsection in Sphinx stub`

**Files created:**

- `orpheus/derivations/continuous/sood_registry/cache.py` (NEW, ~350 NCSL incl. docstrings)
- `tests/derivations/test_sood_registry_cache.py` (NEW, 15 foundation tests, 322 lines)

**Files modified:**

- `orpheus/derivations/continuous/sood_registry/__init__.py` — re-exports `SoodResultCache`, `sood_cache`, `clear_cache`, `cache_info`
- `.gitignore` — added `.cache/` exclusion
- `docs/theory/sood_registry.rst` — new "Solver-output caching" subsection + Phase B preview updated to point at it

**Test gates:** 15/15 cache foundation tests pass in 0.22s; 80/80 prior sood_registry + fn_method tests still green at the same tolerances. Sphinx -W exits 0.

---

## Design choices

### API choice: BOTH decorator and class form

**Rationale.** The two APIs serve different audiences:

- **Decorator** (`@sood_cache()`) is the recommended path for tests that wrap a solver call once and want the cleanest call site. It composes naturally with existing solver functions (`solve_fn_sphere_bare_critical`, etc.) — wrap a new function, call with kwargs, done.
- **Class** (`SoodResultCache`) is for tests that need programmatic invalidation, introspection of stored entries, or fine-grained control (e.g., the cache test suite itself, which uses `cache.put()` and `cache.get()` directly to model corruption / version-mismatch scenarios). It also makes the underlying state visible — `cache.info()` returns a list of every entry with timestamps and paths.

The decorator delegates to the class internally, so they share all storage / hashing / version logic. No code duplication.

### Cache directory: `.cache/sood_registry/` at project root

**Rationale.** Project-root keeps cache state out of the source tree (`orpheus/`) and out of the test tree (`tests/`). The `.cache/` namespace is conventional (npm, cargo, etc. all use it). Project-root resolution is done by walking upward from the cache module to find `pyproject.toml` — independent of `cwd`, which sub-agents and pytest both reset between calls.

### Hash function: SHA-256 (16-char prefix) over canonicalized JSON

- **JSON canonicalization** with `sort_keys=True`, `separators=(",", ":")` for deterministic byte sequence.
- **Custom canonicalizer** for numpy: `np.float64(1.30)` → `1.30`, `np.ndarray` → `{"__nd__": True, "shape": [...], "dtype": "...", "data": [...]}`. This makes solver kwargs hash-stable regardless of how the caller materializes scalar values (Python float vs numpy float64 are the most common collision case in tests).
- **16-char prefix of SHA-256** → 2^64 keyspace, far above what a few thousand cache entries could collide on. Filename-safe.

### Version: auto-detect via `git rev-parse --short HEAD`

- Default `version="auto"` → invokes `subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=project_root)`.
- Falls back to `"nogit"` outside a git checkout (offline/tarball runs).
- Callers may pin an explicit string for reproducibility.
- A version mismatch on read = miss; the entry is overwritten on the next `put()`. This auto-invalidates entries on solver-implementation changes that bump the SHA.

### Storage format: pickle, one file per entry

- **Pickle over JSON** because solver result types are dataclasses with numpy arrays. JSON would need a serialization shim per result type (and would silently lose precision on long-decimal floats).
- **One file per entry** (`<solver-slug>_<hash>.pkl`) means a corrupt write only poisons one entry — the rest of the cache is untouched. Compare to a single shared JSON: a single bad write poisons everything.
- **Atomic writes** via tempfile + `os.replace()` — protects against half-written files if the process is killed mid-write.

---

## Speedup measurements (no solver was migrated this phase)

Per the optional-low-priority directive, NO existing slow test was migrated to use the cache. The sphere F_N test at N=15 takes ~1s (already fast enough) and the Variant α benchmarks have their own scope held by Phase B1/B2/B3 — touching them would have crossed scope boundaries.

The L1 smoke test in the cache suite wraps the transfer-matrix `kinf_homogeneous` solver (microseconds; the cache lookup adds overhead, not subtracts). It exists to confirm end-to-end wiring with a real production solver, not as a speedup demonstration.

**Where the cache will pay off (priority order):**

1. **Sphere F_N at high N.** `solve_fn_sphere_bare_critical(c=1.30, n_modes=N)` runtime grows polynomially with N. The convergence-with-N test sweeps N ∈ {5, 8, 10, 12, 15} on every run. Cache adoption would reduce this from ~5 s (sum of computations) to ~few ms (5 disk reads) on warm-cache runs.
2. **Variant α at fine quadrature.** Sphere/cylinder/slab benchmarks at `n_r=24, n_mu=128` are the slowest existing tests (a couple of seconds each). The Sood-2003 cylinder cross-check at `(24, 20, 64, 96)` takes ~36 s for two test functions. A cache wrapper at the test level would compress repeat runs to <1 s.
3. **Future Westfall-Metcalf cylinder F_N.** When Phase B1 lands, the cylinder F_N solver will likely be in the same speed class as sphere F_N — high-N convergence sweeps will benefit equivalently.

The wiring is now in place; adoption is a per-test decision.

---

## Honest verdict

### Cache utility — high-value primitive, but adoption-driven

The cache's value is realized only when slow solvers are wrapped. The cache itself adds zero cost to existing tests (none import it). Tests that DO opt in get a 1-2 order-of-magnitude wall-clock speedup on the second run forward, with auto-invalidation on solver-implementation changes.

### Cache poisoning — real but bounded

A buggy solver's wrong answer survives across runs until either (a) the SHA changes (auto-invalidation) or (b) a user runs `clear_cache()`. **Mitigation in current usage:** every Sood-case test asserts the solver result against a published reference value, so a poisoned entry STILL fails the test — the assertion is the protection, not the cache. The cache cannot mask a bug that the test would otherwise catch; it only speeds up correct solver runs.

**Hazard scenario:** a developer adds a cached test against a solver that does NOT have an external-reference assertion (e.g., a self-consistency check that compares two solver invocations). In that case, both invocations could hit the same poisoned entry and agree trivially. **Recommendation:** when adopting the cache in a new test, ensure the assertion is against a registry truth value or a structurally-independent reference — never against another cached call.

### Architectural seams

1. **The graph-rebuild pattern.** Nexus-style "rebuild only on graph change" workflows fit this cache's shape closely. The cache key + version semantics could be reused for archivist-rebuild-on-graph-change patterns: hash the relevant document/code state, store the rebuild output, invalidate on hash bump.
2. **The Variant α MR matrix cache.** The 2024 Variant α multi-region work has its own implicit memoization concerns (large rank-2 matrices recomputed across tests). The pattern here generalizes: any deterministic numerical procedure that returns a picklable artifact and is called repeatedly with the same kwargs is a cache candidate.
3. **CI cache hits.** A future CI improvement: persist `.cache/sood_registry/` between CI runs (keyed by branch SHA) to amortize the slow-test cost across PR iterations. The current design supports this directly — the cache directory is already a self-contained artifact.

### Concerns / known limitations

- **No file locking.** Concurrent pytest workers (xdist) writing to the same key produce a "last writer wins" outcome. Harmless for deterministic solvers (identical bytes), but not transactional. The module docstring documents this.
- **Pickle is Python-version-sensitive.** A cache populated by Python 3.12 may not be readable by Python 3.10 if a stored type's `__module__` path differs. This is a soft-stale: caught by `pickle.UnpicklingError` and evicted gracefully.
- **No size cap or LRU eviction.** Long-lived caches grow unboundedly. Acceptable for ORPHEUS's scale (a few thousand entries at most); a future enhancement could add an LRU cap.

---

## Manifest cross-check

- [x] Branch-1 SymPy module: **N/A** — this is infrastructure, not a verified mathematical method. The cache stores other solvers' outputs; it does not derive equations.
- [x] Foundation-tagged test gate: `tests/derivations/test_sood_registry_cache.py` (15 tests, all `@pytest.mark.foundation`).
- [x] Branch-2 production solver: `orpheus/derivations/continuous/sood_registry/cache.py`.
- [x] L1 cross-check: the smoke test wraps `kinf_homogeneous` (transfer-matrix reference) and verifies byte-equal payload + value matches Sood truth at 1e-5. (NOT a structural-independence cross-check in the rigorous V&V sense — it's a smoke of the cache's plumbing against a known-good reference.)
- [x] Sphinx stub: `docs/theory/sood_registry.rst` — new "Solver-output caching" section with TODO marker for archivist expansion.
- [x] Closeout memo: this file.
- [ ] DISPATCH_REQUEST to archivist: per scope, deferred — the cache TODO marker will be expanded into rich narrative when the archivist next sweeps the sood_registry page (likely after Phase B converges and the wide-enumeration content lands too).

---

## Self-improvement notes

No new anti-patterns surfaced. No skill edits proposed. The discipline mapped cleanly:

- `algebra-of-record` — N/A for infrastructure; the State-1A/B/C taxonomy is for mathematical references.
- `vv-principles` — the L1 smoke is honestly described as plumbing-verification, not a structural-independence claim. No false elevation of the cache's verification status.
- `numerical-bug-signatures` — N/A; no numerical method involved.
- `cross-domain-frames` — checked: caching infrastructure is its own domain, not a foreign-frame match for a transport problem.
- `subagent-handoff-protocol` — followed: archivist dispatch deferred to a later sweep, not emitted here per scope.
