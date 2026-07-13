# SN mutation-testing harness (sensitivity metric for the sentinel set)

This directory holds **cosmic-ray** configs that measure how close the
SN sentinel set (`pytest -m sentinel`) comes to "for any regression,
the set changes". A regression ≈ a behaviour-changing code edit ≈ a
*mutant*; "the set changes" = a sentinel *kills* it. The fraction of
killed mutants is the **mutation score** — the honest sensitivity
number. No small set reaches 100 %; the score quantifies the gap.

See `.claude/plans/sn_sentinel_harness.md` for the full design.

## Tool verdict (Phase S0, 2026-06-01, Python 3.14.3)

**cosmic-ray 8.4.6 — WORKS. Chosen.** mutmut 3.5.0 installs and its
CLI runs, but cosmic-ray won is the right tool here because:

| Criterion | cosmic-ray 8.4.6 | why it matters |
|---|---|---|
| Runs on Py 3.14.3 | yes (CLI + `local` distributor) | xdist already broke on 3.14; compat was a real risk |
| Scope to ONE module | `module-path = "orpheus/transport/spatial/diamond.py"` | the taxonomy makes mutation tractable: mutate one tier's module, run only that tier |
| No xdist needed | `distributor.name = "local"` | xdist 3.8.0 DEADLOCKS on Py 3.14.3 (see `sn_taxonomy_reorg_mapping.md`) |
| Per-mutant speed | ~1.0 s/mutant (≈ baseline cost) | 374 mutants for diamond.py ≈ 6 min single-process |
| Survivor localization | `dump` carries `definition_name` + `start_pos` per mutant | a survivor maps to the enclosing function → the capability node |
| Does NOT mutate strings | confirmed (zero docstring-noise mutants) | survivors are signal, not docstring churn |

**diamond.py full score: 373/374 killed = 99.7 %** (update 100 %,
update_batch 100 %, residual 99.3 %). The 1 survivor is a
`FloorDiv_Sub` in `residual()` — see "Known survivors" below.

## Per-capability scoping recipe

The capability taxonomy (`.claude/plans/sn_test_taxonomy.md`) makes
mutation tractable: mutate the module(s) of ONE capability node, run
ONLY that node's tier tests (or the node's sentinels). Survivors that
map (via `definition_name`) to a function tested by a DIFFERENT tier
are a *scoping artifact*, not a coverage hole — include that tier's
test file in the `test-command`.

Worked example from S0: scoping diamond.py against only
`test_diamond.py` + `test_cell_balance_for_streaming.py` left
`update_batch` (the 2-D Cartesian path) at 0 % — because that path is
tested by `test_cell_update_batch.py`. Adding the third file → 100 %.

```toml
[cosmic-ray]
module-path = "orpheus/sn/sweep/<module>.py"   # one module / one tier
timeout = 30.0                                     # bound mutant infinite loops
test-command = "env PYTHONPATH=<worktree> <venv>/bin/python -O -m pytest <that tier's tests> -q -p no:cacheprovider --timeout=30"
[cosmic-ray.distributor]
name = "local"                                     # NEVER xdist on Py 3.14.3
```

```sh
# 1. baseline: unmutated suite MUST pass (cosmic-ray raises otherwise)
python -m cosmic_ray.cli baseline <cfg>.toml
# 2. init: generate the work order (mutant DB)
python -m cosmic_ray.cli init <cfg>.toml <cfg>.sqlite
# 3. exec: run every mutant (single-process, local distributor)
python -m cosmic_ray.cli exec <cfg>.toml <cfg>.sqlite
# 4. score: dump JSON, classify by definition_name → capability node
python -m cosmic_ray.cli dump <cfg>.sqlite | <score script>
# 5. ALWAYS restore the module — see gotcha below
git checkout -- orpheus/sn/sweep/<module>.py
```

## GOTCHA — cosmic-ray leaves the LAST mutant on disk if killed

The `local` distributor applies a mutation, runs tests, then reverts.
If `exec` is interrupted (Ctrl-C, timeout-kill), the in-flight
mutation stays applied to the source file. **Always
`git checkout -- <module-under-test>` after a run** and before
trusting any subsequent test result. S0 caught this: diamond.py was
left with `* cell_avg // (...)` (the surviving mutant) on disk.

## Sentinel-set sensitivity on diamond.py (Phase S2)

`tests/_mutation/diamond_sentinels.toml` runs `pytest -m sentinel` (no
`-O`) against every diamond.py mutant. The sentinel subset is a FAST
TRIPWIRE, not the deep net, so its score is lower than the full suite's
99.7 % — that is expected and honest.

| sentinel-set version | update_batch | update | residual | TOTAL |
|---|---|---|---|---|
| S1 (15 node sentinels) | 100 % | 33 % | 0 % | 42.5 % |
| S2 (+spatial-closure +linearity) | 100 % | 96 % | 57 % | 81.3 % |
| **S2 final (+round-trip 2G/4G het)** | **100 %** | **96.1 %** | **97.9 %** | **96.8 %** |

S2 final = **96.8 %** (362/374); **98.1 % excluding the 5 equivalent
mutants** (`ReplaceTrueWith/FalseWith` on the `is_linear` /
`is_positivity_preserving` ClassVar metadata + `RemoveDecorator` on
`@dataclass`). Exceeds the ≥90 % S2 target for the core sweep
capability. The S2 gap-closers cost ~0 wall-clock (the gate stays at
~4.1 s): they promote already-cheap `DiamondDifference` unit tests
(`test_diamond.py::TestResidual`, `TestBitIdenticalCurvilinear`,
`test_dd_recurrence`) to sentinels.

The 12 remaining survivors are all equivalent / near-equivalent: 5
ClassVar-flag + decorator mutants (no numeric effect); 4 in `update`
are `Gt`→`GtE`/`NotEq` flips on the `face_area_downstream > 0.0` /
`angular_upstream is not None` **structural-presence** guards (flipping
`>` to `>=` at a geometric `0.0` boundary is equivalent for
non-degenerate cells) + `Sub_Mod`/`NumberReplacer` on cold branches; 3
`NumberReplacer` constant tweaks in `residual` the round-trip
tolerates. The full per-tier `test_diamond.py` is the deep net (S0:
99.7 %); the sentinel set is the fast tripwire.

## Known survivors (diamond.py)

- `residual()` line 269, `-` → `//`
  (`denom[:,0]*cell_avg - (...)` → `... // (...)`). Survived because
  `residual()` is exercised mostly at the converged cell-average where
  the true residual ≈ 0; `x // y` for `x ≈ y > 0` gives `1.0` while
  `x - y` gives `≈ 0`, and no current test asserts the residual
  *magnitude* against a strongly non-converged `cell_avg` input. This
  is a (minor) sensitivity gap in the apply-direction round-trip, NOT
  a production bug. Candidate for a future cheap negative test
  (assert `residual(cell_avg + δ) ≈ denom·δ`). Logged, not papered
  over.
