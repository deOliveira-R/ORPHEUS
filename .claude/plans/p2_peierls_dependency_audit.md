# P2 (#578) — the Peierls ray layer's retirement: dependency audit

The table `retirement-audit` G.24 asks for, for the symbols the P2 carve
retired from `orpheus/derivations/continuous/peierls_nystrom/geometry.py`
(branch `refactor/kernel-p2-reference`). It was written after the deletion
rather than before it. The order it records is the order the carve followed:
leaves first, then the helpers that only the retired set called, then the
methods themselves. The counts are a measurement of the tree after the carve.

**Instrument** `[M]` 2026-10-11. The script is
`scratch/characteristic_architecture/p2/residual.py`, written to scratchpad
and copied there. It uses Python `re`, with the symbol bounded on both sides
by "not an identifier character". It reads `.py`, `.rst`, `.md`, `.txt` and
`.toml` files under `orpheus/`, `tests/`, `docs/` (excluding `docs/_build/`),
`.claude/` and `tools/`.

**Positive controls.** The same filter finds the live successors: `rays_from`
has 35 hits in `orpheus/`, 29 in `tests/` and 22 in `docs/`; `lengths_between`
has 5, 20 and 6.

**What the import column rests on.** The full non-slow suite ran on the carved
tree (`scratch/characteristic_architecture/p2/full_nonslow_2.log`). A live
import of a deleted symbol would fail at collection.

| retired symbol | role | `orpheus/` hits | `tests/` hits | `docs/` hits | `.claude/` hits | successor |
|---|---|---:|---:|---:|---:|---|
| `_chord_tau_mu_sphere` | leaf; dead (0 callers at `a7aeebb5`) | 0 | 0 | 0 | 7, all archaeology | none |
| `ObserverRays`, `Run` | helpers, called only by the retired methods | 0 | 1, a history docstring | 0 | 0 | `orpheus.geometry.chord.Rays` |
| `source_position` | method | 0 | 2, history docstrings | 5, all past tense or the "Replaced" column | 4, archaeology | `Rays.orbit_coordinate_at` |
| `which_annulus` | method; outer-biased locator | 0 | 4, history docstrings | 3, past tense or "Replaced" | 5, archaeology | `CurvilinearGeometry.region_of` (inner-owns) |
| `optical_depth_along_ray` | method; shell and segment loops plus the inline locator | 0 | 8, history docstrings | 6, past tense or "Replaced" | 11, archaeology | `Chord.line_integral` over `Rays.ahead`, `Rays.to_wall` or `Rays.lengths_to` |
| `rho_inner_intersections` | method | 0 | 4, history docstrings | 2, "Replaced" | 2, archaeology | `Rays.wall`, `Rays.to_wall` |
| `rho_max` | method | 1 (a comment naming the quantity ρ_max, read through `rho_max_val`) | 5: one is the quantity in a message string; 4 are history | 6, past tense or "Replaced" | 12, archaeology | `Rays.ahead.sum(-1)` |

**Sorted by tense** (`retirement-audit` item 19). There are 0 present-tense
claims and 0 imperatives naming a retired method, of the 18 non-archaeology
hits read one by one. Each remaining hit is one of three kinds:

- a history docstring, "until #578" or "retired in #578";
- a row of the "Replaced" column in the `peierls-ray-layer` table of
  `docs/theory/references/peierls_nystrom.rst`;
- the mathematical quantity ρ_max, which survives the method.

**The docstring surface.** `dead_references` was run after the Sphinx rebuild;
its result is recorded in the plan's landing note.
