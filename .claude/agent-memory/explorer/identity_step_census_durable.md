---
name: identity-step-census-durable
description: Durable findings from the #459 identity-step census (2026-09-12) — where scattering_order is clamped vs padded per entry, which generating data have NO content equality, the one flipping pin, the one hash-keyed SNMesh consumer, the closure-class contest. Read before any Problem-identity / scattering_order carve.
metadata:
  type: project
---

# Identity step (R-cc3/R-cc4, #459) — the durable half of `scratch/_consumers/explorer_identity_step_census.md`

Snapshot 2026-09-12 @ `4207c1d6` (no `orpheus/`/`tests/`/`docs/` change during the census). The
`file:line` map is in the scratch memo and rots; these are the SHAPE facts.

**Why:** the consumers campaign's first step makes a Problem's identity the CONTENT of its
generating data and moves `scattering_order` onto the Problem; the carve crosses `sn/mesh/`,
`homogeneous/`, `data/macro_xs/`, `sn/solver.py` and every order consumer.

**How to apply — the facts that shape the design, not the counts:**

- **The retained order has THREE spellings, and the entries disagree.** `SNSolver.__init__`
  clamps `min(L, min over materials of len(SigS)−1)` (ruling O-1: the scattering stack alone);
  the two ADJOINT entries bypass `SNSolver` (`_adjoint_posing_parts` → `build_within_group_system`
  → `TransferKernel.at_order`), which PADS zero moments above the stored order — no clamp, no
  record; `solve_sn_fixed_source` hands DSA the RAW request before the solver clamps, and DSA
  re-derives "is ℓ≥1 live" from `len(mats) > 1` per material. ⟹ a hub-held CLAMPED order is what
  makes "P0 forward == P3 adjoint" true; the key must hold the clamped value.
- **No result object records `L`** — `Solution`/`SolutionBase` have 0 occurrences; the 36
  regression `.npz` carry no order key; only the cfg table and the moment field's SPACE carry it.
- **Only ONE truncation order needs a home.** The spatial moment order is a `ClassVar` on the
  scheme TYPE (DD 1 / LD 2 — covered by `type(scheme)`); the frame degree is minted from
  `scattering_order`, a field's own order, or fission's hard 0 (`frame.rst` defers a
  `max(L_flux, L_scatter)` rule until a detector/PN consumer exists).
- **Every array-bearing generating datum RAISES on `==` and `hash` today**: `Mixture`,
  `Quadrature`/`DiscreteMeasure`/`LevelStructure`, `AxisMesh`/`Mesh1D`/`Mesh2D`; `BC` compares
  but is unhashable (dict params); `Materials` is `eq=False` by ruling; `LinearDiscontinuous`
  is a plain class (identity), `DiamondDifference` a frozen 0-field dataclass. `Quadrature`
  carries NO rule name — content = nodes/weights bytes (+ hashable `support`/`invariance_group`/
  `exactness`/`folded_by`). `PrescribedInflow`/`AlbedoBoundary.reemission` hold Protocol
  callables — no content identity is possible there. `Mixture` is mutable and 33 test sites
  mutate it post-construction (0 in production).
- **The precedent is `Axis._identity_key`**: per-class tuple of primitives + `tobytes()` for
  arrays, generator/provenance EXCLUDED, `NotImplemented` across classes, `hash((cls, key))`;
  a persisted hash needs a digest (`hash()` of bytes is process-salted).
- **`is_same_phase_space` has exactly ONE message pin and it FLIPS**
  (`test_compare_cross_mesh_rejected` builds two same-data meshes and expects a refusal); the
  discriminator fixture (`_quad8_mesh`) already ships in the same file. The #281 collapse arms
  (`homogenize`/`condense` with `adjoint=`) are all over shared constituents and stay green.
- **The one hash-keyed `SNMesh` consumer is `_GEOM_CACHE_INTERN` (a `WeakKeyDictionary`)**
  whose entry carries a closure-INSTANCE identity check — two content-equal live meshes would
  alias the key and fight over the entry. A production `is not` guard on `SNMesh` also lives in
  `ScheduledInvertibleOperator.__init__` (0 pins). `Solution.__eq__` is dataclass-generated and
  includes `mesh` (benign, but the hub's `__eq__` must never raise).
- **R-cc3's "closure class" row is CONTESTED by the tree**: the predicate's docstring and
  `sn/index.rst` say "Do not strengthen the predicate by adding the closure"; the first attack's
  row 4 says the opposite. A user fork, with two doc sites + one docstring stating the opposite.
- **Homonym**: `La13511Case.scattering_order` (54 sites, a derivations record FIELD) is not the
  solver kwarg — exclude it from any name-keyed sweep.

Related: [[lessons]] L-041, L-042, L-043; the prior census `problem_solution_split_census.md`.
