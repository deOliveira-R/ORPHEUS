# P3 step 3f: the dependency audit for product integration under route B (retirement-audit G.24)

Destination: `.claude/plans/p3_product_integration_dependency_audit.md` (the explorer's write scope is `scratch/`; the
orchestrator moves it). Written 2026-10-13 by the explorer on `main` at `114141d1`, clean tree under `orpheus/`, `tests/`,
`docs/`. Plan: `.claude/plans/characteristic_reference_architecture.md`, "### 3f" sections to the end. Prototype read:
`scratch/characteristic_architecture/p3/product/{product2,product3,law,rule_stieltjes}.py`, `gate_spec.md`.

**Instruments.** Nexus (`callers`, `neighbors`, `dead_references`): the graph was built at `4ebfa7a4` (dirty, branch
`refactor/characteristic-axial-lift`), 6 commits behind HEAD with 16 indexed files changed, the whole package among them;
every node it returned carries the stale warning. Its caller sets agreed with grep on every symbol checked (5 of 5:
`exponential_ends`, `ExponentialGrading.of`, `graded_ends`, `halvings`, the `grading` module's importers). Grep is
`git grep` over `orpheus tests tools docs .claude` excluding `docs/_build`; methods and constructors by grep and a Python
AST pass (`ast_pin.py` in the session scratchpad: test function to symbol). `dead_references` baseline: 0 dead of 65
checked (stale graph). Every `docs/api` automodule of the package renders `:members:` with no `:noindex:`
(`docs/api/derivations.rst:239-271`), so a dead `:class:`/`:func:`/`:meth:` to a PUBLIC name is caught by `sphinx -n`
(set diff); a private name (`_PieceSet`, `_attenuated`, `_Toward`) is never rendered, so grep is its only gate.

Line numbers are from `grep -n` of the anchor named beside them.

## 1. The table: symbols that place nodes or grade pieces along the attenuation

"prod" = `orpheus/`; "tests" = `tests/`; counts are call or construction sites of the symbol, excluding its definition.
Leaves first within each fate group (G.24 order).

### 1a. Retire under route B

| symbol (def, anchor) | callers | tests that pin it | docs + docstrings | pins | fate, reason |
|---|---|---|---|---|---|
| `grading._WIDEST_FIRST` (`grading.py:38`, `_WIDEST_FIRST = 1.0`) | prod 1: `ExponentialGrading.of` (`grading.py:93`, `min(resolved, _WIDEST_FIRST)`) | `test_characteristic_lift.py::test_the_first_graded_piece_is_capped_at_one_mean_free_path` (docstring `:461`, behaviour via `.of`) | none outside `grading.py` | `grading.py` skeleton (module constant) | **retire** (leaf). The measured cap on the first exponential piece exists only for `.of`. |
| `ExponentialGrading.of` (`grading.py:86`, `def of(cls, points`) | prod 2: `_PieceSet.of` (`transport.py:280`, `grading = ExponentialGrading.of(rule.points)`), `TraversalRule._attenuated` (`transport.py:689`, `ExponentialGrading.of(self.inner_points)`) | lift `:439` `test_the_derived_grading_integrates_the_attenuation_to_double_precision` (`l0`, `verifies("characteristic-traversal-integrals")`); `:459` the cap test; `:505` `test_the_exponential_ends_follow_the_gradings_first_depth_and_ratio` (param `which`); `:521` `test_a_familys_slots_are_graded_from_its_slowest_and_fastest_members` (asserts `gradings == {ExponentialGrading.of(8), ExponentialGrading.of(12)}`, `:542`) | `characteristic.rst:1975` (`:meth:` xref, `sphinx -n` sees it), `:2039`, `:2107-2109` (`.. implements:: characteristic-exponential-grading :by: ...ExponentialGrading.of`), `:9292`; `matrix.rst:1082` lists the equation; docstrings `transport.py:43,246,678`, `lines.py:72` (`:class:` to the class, survives) | `grading.py` skeleton (classmethod signature) + DefPin `ExponentialGrading.of` | **retire.** The node-count law it derives is replaced by `law.node_counts` (m from the contour bound, independent of x). The `implements` directive and the equation `characteristic-exponential-grading` (no `verifies` marker, gate spec item 7) lose their implementer. |
| `transport._Toward` (`transport.py:238`, `class _Toward(Enum)`) | prod 3: `_toward_exits` / `_toward_entries` / `_toward_both` (`transport.py:595,600,605`) | none by name (0 of 6 characteristic test files; control: `_Toward` greps 8 hits elsewhere) | `characteristic.rst:1837` (literal); docstrings `transport.py:45,262` | `transport.py` skeleton (class) | **retire.** The three piece sets collapse into one geometric set (plan 3f "Pieces"). |
| `TraversalRule._attenuations` (`transport.py:617`, `def _attenuations`) | prod 2: `_PieceSet.of` (`:279`), `_attenuated` (`:689`) | none by name | docstring `:621` cites `exponential_ends` | DefPin | **retire.** Feeds only the along-line exponential gradings. |
| `TraversalRule._attenuated` (`transport.py:665`, `def _attenuated(self, slot`) | prod 2: `_to_end` (`:703`, def), `volterra` (`:802`, `own = self._attenuated(`) | none by name; behaviour pinned through `volterra`/`at` rows (T5, T7, AS*) | `characteristic.rst:2431`; `error_catalog.rst:9944` (ERR-100 **Fix**, present tense "One body, ``TraversalRule._attenuated``, integrates"); docstrings `transport.py:68,110` | DefPin | **retire.** The two-level inner rule becomes the triangle's Gauss rule of e^{-x tau} (`product3.WeightGaussSet.triangle`). ERR-100's Fix prose is a present-tense must-fix (E.19). |
| `_attenuated_length` (`transport.py:104`) | prod 3: `_PieceSet._integral` (`:333`), `_attenuated` (`:696`), `volterra` (`:799`) | none by name | `characteristic.rst:2403`; docstring `transport.py:62` | `transport.py` skeleton (signature) + DefPin | **retire.** Under route B the weight IS the attenuation; `volterra`'s carried-in term becomes a contraction with `to_start` (`product3` `volterra`). |
| `_PieceSet._integral`, `_PieceSet.nodes` (`transport.py:324`, `def nodes`; `:330`, `def _integral`) | `nodes`: prod 3 (`to_end`, `to_start`, `volterra` `:793`); `_integral`: prod 2 | none by name | `characteristic.rst:2421` (`_PieceSet.to_end`/`to_start`), `:2504` (`on_columns`); docstrings `transport.py:65,109,830` | DefPin | **replace.** GL nodes per piece shared over the family become offsets per (piece, speed) (`rule_stieltjes.weight_rule_offsets_batched`). `to_end`/`to_start` keep their names in the prototype (subclass). |
| `_PieceSet.of` (`transport.py:276`) | prod 3 (the `_toward_*` properties) | indirect: lift `:521` spies `transport_module.exponential_ends` while building `TraversalRule.of(...).outflow()` | `transport.py:262-269` | DefPin | **replace** by the geometric piece set (slot ends, hp walk; the slot-middle cut is rejected, plan 3f law). |
| `TraversalRule._branch_edges` (`transport.py:566`) | prod 1: `_PieceSet.of` (`:289`) | `test_characteristic_transport.py::test_the_turning_grading_halves_toward_the_closest_approach_at_seven_points` (behaviour; Nexus lists it under the grading query) | `characteristic.rst:2261`; `error_catalog.rst:9822` (ERR-099) | DefPin | **replace** by the closed-form hp walk of tightness c (`product2.graded_ends`, QA4 fix: D^2 = (b/|POmega|)^2 + (e + a)^2, bounded step count). It is the only along-line caller of `grading.halvings`; the walk does not call `halvings`. |
| `TraversalRule._nodes`, `TraversalRule._orbit_coordinate` (`transport.py:652`, `:638`) | `_nodes`: prod 2 (`_PieceSet.nodes`, `_attenuated`); `_orbit_coordinate`: prod 1 (`_nodes`) | none by name | `transport.py:325` | DefPin | **replace.** Node placement through the parameter (`chord.orbit_coordinate_at(start + distance)`) is the route RB3 refutes (93 eps held in s; 8e13 to 4.5e15 eps round trip on the slab face). Successor: the offset form, `Rays.orbit_coordinate_at` (section 3). |
| `TraversalRule.inner_points` field (`transport.py:428`) and `TransportResolution.inner_points` (`assembly.py:109`) | **signature census** (F, members spelled without the symbol): `TraversalRule.of`/`.on_chord` prod 2 (`lines.py:481`, `reading.py:283`), tests 12 (lift `:85,105,352,538`; transport `:145,159,684,1047,1050,1054,1067,1072`); `Lines.chunks(points, inner_points)` prod 2 (`assembly.py:242`, `reading.py:247`), tests 3 (assembly `:321,494`; lift `:689`); `LineRule.transport(points, inner_points, support)` prod 1 (`system.py:203-204`, multi-line), tests 22 (assembly 19, reading 3); `TransportResolution(...)` prod 0 (control: tests 26 sites in 7 files: `_characteristic_ladders.py` 3, convergence 3, reading 5, reference 13, system 2, `numerics/test_traced_memo_clients.py:149` 1; docs 3 at `characteristic.rst:4138,6604,6959`) | `test_characteristic_reference.py:1196-1201` (the content-identity legs, `"inner_points": (leg("more inner points", ...))`); `:461` `test_a_malformed_resolution_is_refused`; the ladders' `_RUNG_AXES` (`_characteristic_ladders.py:94`) | `characteristic.rst:1693,1978,2430-2432,3482,4591,5407,8953`; `assembly.py:103` docstring; `.claude/agent-memory/elegance-enforcer/lessons.md` (memory) | `assembly.py` skeleton (dataclass fields are class attributes) AND the memo KEY (the key encodes `Resolution`, whose `transport` is a `TransportResolution`) | **retire** (`inner_points`: the triangle has no inner rule of its own, plan 3f). `points` leaves the resolution too under ruling (a) "once the derivation is measured": same census, and `TransportResolution` then holds `line_points` alone. Every memo key moves once. |

### 1b. Survive (rules over lines, the basis, the hp halving outside the line)

| symbol | callers | tests that pin it | docs + docstrings | pins | fate, reason |
|---|---|---|---|---|---|
| `ExponentialGrading` (class, `grading.py:42`) | prod 1 construction: `DYADIC = ExponentialGrading(1.0, 2.0)` (`grading.py:98`); `exponential_ends` reads `.first`/`.ratio` | lift `:459` (`assert DYADIC == ExponentialGrading(1.0, 2.0)`, `:472`) | `characteristic.rst:3506` (`:class:`), `transport.py:43,246,678`, `lines.py:72` | `grading.py` skeleton | **survive as the dyadic grading's (first, ratio) holder.** Its docstring (the A&S 25.4.30 bound, the table at `:65-73`, "Each piece of the grading is integrated by an n-point Gauss-Legendre rule") is present-tense about along-line pieces that no longer exist: must-fix (E.19), or fold into `DYADIC`. |
| `exponential_ends` (`grading.py:116`) | prod 4 today: `_PieceSet.of` x2 (`transport.py:291,293`) and `_attenuated` (`:689`) **retire**; `lines._grazing_ends` (`lines.py:86`) and `impact_rule` (`lines.py:160`) **stay**, both with `DYADIC` | lift `:427` `_graded_integral`, `:505`, `:277` (docstring); **COUNTING-SPY** lift `:521` `monkeypatch.setattr(transport_module, "exponential_ends", recording)` (`:536`): raises once `transport` stops importing the name (setattr's default `raising=True`), so the row is deleted or re-pointed, not left; lift `:556` spies `lines_module` (stays) | `characteristic.rst:1875,2110,2114` (`.. implements:: characteristic-exponential-grading :by: ...exponential_ends`), `:3499,9304`; docstrings `transport.py:41,265,621`; qa memory `lessons_archive.md:6121` (past tense) | `grading.py` skeleton unchanged unless its docstring-free signature changes | **survive, DYADIC callers only.** Its docstring's "the grading's first depth and ratio (ExponentialGrading, from the integrating rule's node count)" becomes false: must-fix. |
| `DYADIC` (`grading.py:98`) | prod 2: `lines.py:86,160` | lift `:459,505,521` (docstring), `:556` (`assert seen and all(g is DYADIC ...)`, `:579,582`) | `characteristic.rst:2121,2138,2191,3514,3581,6660,9296` | skeleton | **survive.** Polar and cosine coordinates, the impact rule's rim: the exponent is not linear in them (plan 3f "What stays"). |
| `VANISHING_DEPTH` (`grading.py:30`) | prod 2: `exponential_ends` (internal), `lines._grazing_ends` (`lines.py:84`) | lift `:427,439,475,505`; assembly `:1229` (docstring) | `characteristic.rst:1880,3503,4282` | skeleton | **survive.** |
| `graded_ends` (`grading.py:101`) | prod 3 sites in 2 symbols: `PanelBasis.of` (`basis.py:127`), `lines._grazing_ends` (`lines.py:85`), `impact_rule` (`lines.py:158,159`) | assembly `:260` `test_the_impact_rule_is_graded_at_the_slowest_polar_node`, `:1108` `_slab_rule_with_fixed_halvings`, `:1219` `test_every_block_entry_of_a_thin_panel_slab_is_resolved_at_eight_points` (`l1`, `verifies("characteristic-galerkin-assembly")`) | `characteristic.rst:1361,3494` | skeleton | **survive untouched.** No along-line caller. |
| `halvings` (`grading.py:150`) | prod 4: `_branch_edges` (`transport.py:584`) **replaced**; `lines.py:84,158,159` **stay** | indirect through the slab/impact rows above (assembly `:1125` `verifies(...)`, `catches("ERR-101")`; `:1166,1195,1219`) | `characteristic.rst:4282`; grading module docstring bullet "hp" names "the branch points of a line's orbit coordinate (the traversal rule)": must-fix | skeleton | **survive** for the impact rule and the grazing direction; loses its along-line caller to the closed-form walk. |
| `_Packing`, `_Slots`, `_Pieces`, `_carried`, `at`, `outflow`, `entry_response`, `transit_depth`, `per_transit`, `on_columns`, `piece_counts`, `extent` (`transport.py`) | unchanged | assembly `:309` (`extent`, budget) | `transport.py` module docstring "The pieces" (`:33-55`) and "One integrand, two levels" (`:58-70`): present tense, must-fix | unchanged unless signatures move | **survive** (the prototype reuses each). `extent` keeps counting pieces while node arrays grow to F x J x m: the chunk budgets (`assembly.LineRule.of` `chunk=512, budget=1024`; `reading.PointRule.of` `chunk=4096, budget=2**18`) are re-measured (plan 3f open item 3). |
| `PanelBasis.even` (`basis.py:158`) | gains a consumer: the law's even-panel count (p + 1 / 2p + 1), `product2._slot_geometry` | | | | **survive**, new reader. |

## 2. Where the source is evaluated today, and how the speeds share it

- `_PieceSet.nodes` (`transport.py:324`) -> `TraversalRule._nodes` (`:652`) -> `_orbit_coordinate` (`:638`) ->
  `basis.values(orbit, panel)`: shape `(L, J, n)` and `(L, J, n, p + 1)`, **no family axis**. Built once per piece set,
  three times per rule (exits, entries, both).
- `_PieceSet._integral` (`:330`) broadcasts `_attenuated_length(self.rule._speed[..., None, None], ...)` over F against those
  shared values.
- `TraversalRule._attenuated` (`:665`): the inner rule's nodes and `basis.values` (`:694`) are formed once for the whole
  family; only the graded ends depend on the family, through its steepest and gentlest member (`_attenuations`).
- `volterra` (`:781`) reads `both.nodes` for the test functions and `_attenuated` for the triangle; `at` (`:825`) reads
  `exits.to_end`, evaluating nothing new.
- The family is `TraversalRule.speed` `(F,)` / `_speed` `(F, L)`, set by `Lines.speeds` from `axial_lift` (`lines.py:341`):
  the cylinder only; the sphere and the slab have F = 1.

Under route B the pieces and their counts stay speed-free (m does not depend on x), but the nodes are offsets per
(piece, speed): `product3.WeightGaussSet.values_at_offset` evaluates `basis.values` on `(F, L, J, m)` and the triangle's
`values_at` on `(F, L, J, m, m)`. Sharing across speeds is lost: F x m evaluations per piece instead of n.

## 3. How a `Rays` is minted today

- **Constructors:** 2 of 2 `Rays(` expressions tree-wide (`orpheus tests tools scratch/characteristic_architecture`,
  PCRE `(?<![A-Za-z_.])Rays\(`), both inside `ConcentricPartition.rays_from` (`chord.py:194`, `def rays_from`): the slab
  at `:214`, the radial charts at `:218`. No direct construction elsewhere; 0 in the prototype (it reads nodes through
  `rule._orbit_coordinate(t)`, the parameter route, `product3.py:67,78`).
- **Call sites in this package:** 1 of 1, `reading.angular_flux` (`reading.py:282`, `rays = basis.partition.refined_at(c).rays_from(c, omega)`).
  Elsewhere: `fn_method/sphere/flux_reconstruction.py:444`, `peierls_nystrom/geometry.py` (via `CurvilinearGeometry.rays_from`, `:521`).
- **What `Rays.orbit_coordinate_at(distance)` reads** (`chord.py:675`):
  - radial: `image.level`, `image.level_half_chord`, `image.speed` and `self.side`; it does **not** read
    `self.orbit_coordinate` or `self.start`. It is the point's coordinate exactly only when the point sits at the
    image's level (`s = side * y*`). `rays_from` arranges that by chording at the level `(c, c|Omega_x|/|POmega|)`.
  - slab (`AxialImage`): `self.orbit_coordinate + image.rate * ell`; independent of any level.
- **At a piece end or a crossing:**
  - slab: nothing in `Rays` prevents a direct mint `Rays(chord, start=<crossing parameter>, side, orbit_coordinate=<breakpoint>)`,
    and its verb would be exact at the crossing; no site does it.
  - radial: a production chord carries ONE level per line, the impact rule's `(top, half_chord)` (`lines.py:479-480`,
    `self.partition.chord(..., level)`). A `Rays` at another crossing r_j needs the image re-levelled:
    `RadialImage.at_level(radius, half_chord)` exists (guarded to `_LEVEL_AGREEMENT = 64` ulp, `chart.py:599`) and
    `half_chord_at(r_j)` forms y_j cancellation-free for r_j at or above the level, but `at_level` has 1 caller
    (`chord.py:241`, inside `chord`), so nothing mints a re-levelled chord or `Rays` today. A piece end inside a slot (an
    hp-walk end) is on no level at all.
  - `distance` is used linearly with no sign check; `tests/gates/geometry/test_rays.py:192`
    (`test_the_orbit_coordinate_at_the_point_is_the_points`) samples only `ell` in [0, ahead], so a negative offset
    (backward from a piece's end) is unpinned.
- Homonym: `RadialImage.orbit_coordinate_at(t)` / `AxialImage.orbit_coordinate_at(t)` / `Chord.orbit_coordinate_at(t)`
  take a PARAMETER (`chart.py:726,786`, `chord.py:418`); only `Rays.orbit_coordinate_at` takes a distance.

## 4. Traced-memo pins

The two memo clients are `CharacteristicDerivation.solve` and `.evaluate` (`reference.py:348,392`, `@traced_memo`; 2 of 2
in `orpheus/`). A generation pins every first-party MODULE that ran by its skeleton (imports, constants, decorators,
signatures, class attributes) and every DEF that ran by its normalised body (`traced_memo.py:14-17`, `ModulePin`,
`DefPin`). The carve stales every entry through:
- skeleton changes: `grading.py` (`_WIDEST_FIRST`, `of`), `transport.py` (imports, `_Toward`, signatures, fields),
  `assembly.py` (`TransportResolution` fields), `lines.py` (the `from .grading import` line; `chunks`/`_rule` signatures);
- def bodies: `reading.PointRule.row`, `reading.angular_flux`, `system.GalerkinSystem.groups` (`system.py:198`);
- new pins: the rule and law modules, and the generated table if read from a file (`const_rules.npz`: a file pin);
- the key itself, when `TransportResolution`'s fields change.
`tests/gates/numerics/test_traced_memo_clients.py:149` and `test_characteristic_reference.py:1267` (`c3`, memo) construct
`TransportResolution` positionally.

## 5. Agent surfaces (B.4)

`.claude/agents/` and `.claude/skills/`: 0 hits for the retiring names (control: `ExponentialGrading` in
`.claude/plans/characteristic_reference_architecture.md`). Memory: `qa/lessons_archive.md:6121` (`exponential_ends`,
`_PieceSet.of`, past tense, stays), `elegance-enforcer/lessons.md` (`inner_points`, to triage by tense).
