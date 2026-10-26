# The moment-basis carve: the dependency audit of `orpheus/numerics/moment_layout.py` (S4, #423)

This is the `retirement-audit` G.24 table for the retirement in step S4 of `moment_basis.md`. Its rows are the module's eight public names. Its columns are the three surfaces of `retirement-audit` A, and the name's new home. The retirement is one commit, done leaves first: the rank predicates and the counts first, then the tail rule, then the constants, then the module.

Counts are `[M]` 2026-10-25 on the pre-S4 tree (`a838c31e`). The import census comes from the S0 spec (`scratch/layer_dag/moment_basis/gate_spec.md` §6): 23 `import` statements, 14 of them in production across 11 files and 9 in tests across 9 files.

| name | production uses (outside the module) | test uses | docs / harness mentions | new home |
|---|---|---|---|---|
| `is_moment_valued_by_flat_rank` | 3: `sn/loss_representation/__init__.py` ×2, `sn/solver.py` ×1 | 0 | `verification/sn.rst` ×1 | `orpheus/transport/fields/_bases.py` (ruling 4) |
| `is_moment_valued_by_rank` | 2: `sn/loss_representation/__init__.py`, `sweep_graph.py` | 0 | `verification/sn.rst` ×1 | `orpheus/transport/fields/_bases.py` (ruling 4) |
| `cell_moment_count` | 14 sites, 7 files | 2 (`test_spatial_moment_tail_is_the_schemes_axis.py`) | `cartesian_multid.rst` ×1 | `TensorProductBasis.dimension`, read as `scheme.cell_basis(d)` or `spatial_moment_basis(p, d)` (the user, 2026-10-25) |
| `face_moment_count` | 5 sites, 5 files | 0 | `verification/sn.rst` ×1; `.claude/inventories/implements_declaration_inventory.md:574` | `TensorProductBasis.face_dimension(a)` (the user, 2026-10-25) |
| `face_moment_tail` | 22 sites, 8 files | 0 | `cartesian_multid.rst`, `verification/sn.rst`, `history.rst`; the inventory | `orpheus.numerics.space.factor_tail` (ruling 4; named by the user, 2026-10-25) |
| `spatial_moment_tail` | 2 sites (`fields/_bases.py`) | 1 | `cartesian_multid.rst` ×3 | `factor_tail`: it was an alias of `face_moment_tail`, so two names for one rule become one |
| `AVERAGE_MOMENT` | 10 sites, 6 files (2 of them late imports) | 4 files | `cartesian_multid.rst` ×3, `history.rst` ×1 | `orpheus.numerics.basis.interval_legendre_basis.AVERAGE_MOMENT` (the user, 2026-10-25) |
| `SPATIAL_MOMENT_AXIS_LABEL` | 1 (`fields/_bases.py`); since S3 the label is the cell basis's layout label | 5 files | `cartesian_multid.rst` ×1 | `orpheus.numerics.basis.interval_legendre_basis.SPATIAL_MOMENT_AXIS_LABEL` (defined there at S1, so two definitions stood from S1 to S4) |

## The hand spellings re-pointed with the retirement

S3 re-pointed four sites that spelled the basis without importing the module:
- the float root at `fields/_bases.py:533`;
- the Kronecker decode at `_ubld.py:144`;
- the face count at `scheme.py:518`;
- the literal `2` in the LD statics.

S4 found two more hand spellings of the tail rule with `git grep -nE "> 1 else \(\)|== 1 else \(|if [a-z_]+ > 1 else \("` over `orpheus/`:
- `sn/loss_representation/assembly.py:214`, `(cm,) if cm > 1 else ()`;
- `transport/spatial/scheme.py:523`, `(1, 1, 1, m) if m > 1 else (1, 1, 1)`.

Both now call `factor_tail`. The grep's positive control is the rule's own bodies: it found `moment_layout.py:127` and `space.py:1335` through the `== 1 else (` alternative. Its 3 other hits are unrelated conditionals: `walls.py:175`, `face_transmission.py:1043` and `_traced_memo_boot.py:628`.

The literal `2` has no reliable AST shape, so it was audited by grep, as the S0 spec asked. After S3, `cell_moment_count(2` and `face_moment_count(2` return 0 hits in `orpheus/transport/spatial/`. The control is `main`, where the same pattern (`(cell|face)_moment_count\(2`) finds 6 lines in `linear_discontinuous.py` and 1 in `_ubld.py`. `[M]` 2026-10-25, `git show main:<file> | grep -cE`.

## Witnesses

`tests/gates/numerics/test_moment_layout_retired.py` (6 rows) holds the retirement's witnesses:
- the module is gone;
- no import of it remains, by an AST census over `orpheus/`, `tests/` and `tools/` with a positive control;
- no float root, Kronecker decode or θ literal survives;
- the tail rule lives in `orpheus.numerics.space`.

It was red on the pre-S4 tree, 6 of 7 rows by the S0 spec. The stale-rank row moved to `tests/gates/numerics/test_concept_order.py::test_the_concept_order_names_only_modules_that_exist`.

`[M]` 2026-10-25: an AST census of the four retired function names over 1224 files under `orpheus/`, `tests/` and `tools/` finds 0 uses. Its positive control is a synthetic `face_moment_tail(1)`, which it finds.

## Path pins

- `tests/gates/numerics/test_concept_order.py:43`: the tier-0 entry is deleted.
- The test `test_width_one_appends_nothing_and_the_policy_lives_in_the_layout_module` is renamed `..._and_the_tail_rule_is_space_composition`. The location check now lives in the retirement gate's tail row.
- `.claude/hooks/session-health.sh:36,49`: the pyright canary is now `orpheus/numerics/coefficient_layout.py`. The old canary imported nothing first-party, so it could not detect a first-party resolution failure (S0 spec §6). The hook reads `PYRIGHT: OK` on the new canary `[M]`.
- `.claude/inventories/implements_declaration_inventory.md:574`: the implementers are re-pointed to the new homes.
