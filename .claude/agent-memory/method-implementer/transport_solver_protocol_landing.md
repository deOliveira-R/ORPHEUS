---
name: TransportSolver Protocol landing
description: 5-commit landing on feat/transport-solver-protocol — 38 foundation tests, structural unification of math-heart classes + discrete adapters
type: project
---

# TransportSolver Protocol — landing closeout

`feat/transport-solver-protocol` 2026-05-04. Spec at
`.claude/plans/transport_solver_protocol_spec.md` from test-architect
agent `aeeadd1991ba4d1ca`.

## Commit sequence (verbatim from spec)

| Commit | Hash | What |
|--------|------|------|
| 1 | `221b5f9` | Protocol + DiscretizationSpec + File 1 conformance tests (8 tests) |
| 2 | `014bc6f` | Billiard dual-factory + Protocol surface + File 2 (11 bit-equality tests) |
| 3 | `65c1ce5` | MomentSpace properties + File 3 (6 bit-equality tests) |
| 4 | `09eaf44` | File 4 polymorphism regression net (6 tests) + test_eigenvalue.py docstring |
| 5 | `1e0b44b` | Discrete-adapter scaffold + File 6 (7 tests, 1 skipped) + Sphinx page |

Plus 2 commits from parallel agents that landed on this branch
inadvertently: `35d6baa` (Spectrum, between commits 1 and 2) and
`be4a524` (BasisSpace memo, between commits 4 and 5).

## Verification status

- **Foundation tests added: 38 + 1 skipped TODO**
  - `tests/derivations/test_transport_solver_protocol.py`: 8/8.
  - `tests/derivations/test_billiard_from_problem_unified.py`: 11/11
    (bit-equality with asymmetric Sig 3/4/5 catchers).
  - `tests/derivations/test_moment_space_from_problem_unified.py`: 6/6.
  - `tests/cross_method/test_polymorphism.py`: 6/6.
  - `tests/cross_method/test_discrete_adapters_smoke.py`: 6/6 + 1 skip.
- **Existing suites preserved**:
  - `tests/derivations/test_trajectory_resolvent_billiard.py`: 15/15.
  - `tests/derivations/test_fn_method_moment_space.py`: 14/14.
  - `tests/cross_method/test_eigenvalue.py`: 84 tests collected
    (IDs preserved); 53 foundation pass.
- **Sphinx -W**: clean (build succeeded; nexus "skipping" lines are
  pre-existing label-equation gaps, not Sphinx warnings).

## Architectural decisions

### Dual factory on Billiard

Billiard's existing `materials: dict[str, Any]` field clashed with
the Protocol's `materials: dict[int, Mixture]`. Resolution: rename
the legacy field to `xs_payload` (raw XS dict consumed by the
dispatchers), add new fields `materials: dict[int, Mixture]` and
`geometry_spec: GeometrySpec` for the Protocol surface. The
`from_problem` factory accepts both signatures and synthesizes the
absent representation. Bit-equality preserved because the
dispatchers read `xs_payload`, never `materials`. Verified by
`np.array_equal(...)` IEEE-754-exact on every geometry × group
combination.

### geometry_spec property on MomentSpace

MomentSpace's existing `geometry: GeometrySpec` field is heavily
used internally (8 self.geometry refs) and externally (test
asserts `ms.geometry is geom`). Protocol requires
`geometry_spec`. Resolution: add `geometry_spec` as a
`@property` aliased onto `self.geometry`. Property identity is
stable: `ms.geometry_spec is ms.geometry`. Zero breakage.

### test_eigenvalue.py NOT refactored

The spec called for "84 tests preserved + bodies collapsed to
parametrize". Strict reading: this would change pytest collection
IDs (e.g., `test_fn_slab_matches_truth[Ua-1-0-SL-c1.30]` →
`test_critical_dimension_slab[fn_slab-Ua-1-0-SL-c1.30]`) which
contradicts "Refactored file MUST yield identical pytest collection
IDs (CI / pytest-xdist preserved)". I resolved by interpreting
the load-bearing constraint as **identical IDs** and shipping File
4 (`test_polymorphism.py`) as a foundation regression net
demonstrating polymorphism works, while leaving `test_eigenvalue.py`
test bodies unchanged (only docstring updated to point at
`test_polymorphism.py`).

### sphere_mr / hollow geometries on the new factory path

GeometrySpec's schema doesn't carry R_in (only
`critical_dimension_cm`). The new (Protocol) factory therefore
raises `ValueError("hollow_sphere ... use legacy signature")` for
hollow / multi-region geometries. The legacy factory still
supports them. File 2 documents this — `hollow_sphere` /
`annulus` tests use the legacy path only; multi-region sphere_mr
test confirms `materials.keys() == {0, 1}` (per-region Mixture
synthesis from the legacy raw XS dict).

## Branch hygiene incidents

The R0.5 + R2 + MomentSpace round had multiple branch-hygiene
incidents from parallel agents. R3 had zero by being defensive.
This landing had **frequent** trampling — 6+ branch swaps,
working-tree files vanishing, my Commit 2 even landing on the
wrong branch (`feat/basis-space-galerkin-spectral` instead of
`feat/transport-solver-protocol`).

### Recovery procedures used

1. **Cherry-pick from wrong branch**: when Commit 2 landed on
   basis-space, I cherry-picked the commit hash to my branch
   (`git cherry-pick 15fb087`) and reset basis-space back to
   its proper tip with `git branch -f`.
2. **Re-apply working-tree changes**: when checkout swept
   uncommitted edits, I re-read the affected files and re-applied
   the Edit calls from scratch. Lost roughly 30 minutes per
   incident.
3. **Untracked-file rescue**: my Commit 5 untracked files
   (`discrete_adapters.py`, `test_discrete_adapters_smoke.py`,
   `transport_solver_protocol.rst`) vanished entirely once when
   I checked out from a parallel agent's branch. I recreated
   them from memory + my current conversation transcript.

### Prevention recommendations

- **Commit fast** — every working unit gets staged + committed
  immediately. Don't accumulate 3 files of changes before
  committing.
- **Untracked files are vulnerable** — they don't survive
  checkout to a different branch's tip. Stage them as soon as
  they're created (`git add -N` to mark as intent-to-add even
  before committing).
- **Verify branch BEFORE every Bash command that mutates state**
  — `git branch --show-current` is cheap; running tests on the
  wrong branch wastes minutes.
- The R3 lesson "be defensive about branch state" applies in
  triplicate when 4 parallel agents are running. Sometimes
  parallelism isn't worth the coordination cost.

## Spec deviations / non-issues

1. **File 5 NOT collapsed to parametrize** — see "test_eigenvalue.py
   NOT refactored" above. Pragmatic interpretation; the
   regression net ships in File 4 instead.
2. **No L4 cross-check between discrete and reference** — spec
   § "Deliberately NOT verified" item explicitly defers this.
   Placeholder skipped test in File 6 documents the gap.
3. **No archivist dispatch** — per the user-control directive
   from prior memory entries; the Sphinx stub carries TODO
   markers for archivist expansion.

## Foundation-tier discipline preserved

- All new test files carry `pytestmark = [pytest.mark.foundation]`.
- NO `@pytest.mark.verifies(...)` decorators — these tests
  verify software invariants (Protocol conformance, bit-equality
  preservation, structural typing), not equation labels.
- Bit-equality via `np.array_equal()` and `float.hex()`; NEVER
  `np.allclose`.

## Next dispatch (optional)

Archivist expansion of `docs/theory/transport_solver_protocol.rst`
— the page carries one TODO marker requesting a worked
cross-method substitution example (e.g., `MomentSpace` swapped in
for `Billiard` on closed-sphere k_inf). Not blocking.
