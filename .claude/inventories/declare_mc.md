# `implements` declaration candidates — `docs/theory/methods/monte_carlo.rst`

Derived from the tree at HEAD `58e46c6f` (2026-08-18), NOT from
`.claude/inventories/implements_declaration_inventory.md` (refuted on this
population; ignored entirely as instructed).

⚠ Drift check run at close: HEAD moved `a1c90aac → 58e46c6f` during this pass.
`git diff a1c90aac..HEAD` touched `docs/theory/methods/monte_carlo.rst` in
exactly **two** hunks (`7c250811`), both re-pointing the ERR-023/ERR-024
catalogue path from `.claude/skills/…/error_catalog.md` to
`docs/theory/verification/error_catalog.rst`. No equation, label, prose pointer
or rationale comment changed, and no file under `orpheus/mc/`, `tests/mc/` or
`orpheus/geometry/structured_geometry.py` changed at all. All reads below are
against the working tree at `58e46c6f` (clean).

⚠ The graph (`.nexus/graph.db`) was built at 18:10 while HEAD is 19:57 — the
zero-implements measurement therefore predates `bb075c93`
("eleven equations declare that NOTHING implements them"). That commit did not
touch this page; re-confirmed independently against the working tree: the MC
page carries **no** declaration directive on any of the 22 (its only
`.. vv-status:` lines are the four `documented` markers listed below).

## Facts established before the per-equation work

- ⭐ **All 22 equations currently have ZERO `implements` edges** — not one
  declared, not even a name-token guess. Measured against `.nexus/graph.db`:
  every `math:equation:<label>` node in the list has `n_impl = 0`, while the
  corpus carries 13 206 `implements` edges overall (13 121 inferred, 85
  directive). ⟹ **the "partial answer is worse than none" hazard does NOT
  bite here**: there is no inference to stand down, so every declaration is
  pure gain and a missed implementer costs a missing edge, not a lost guess.
  (Enumeration is still complete below.)
- **The MC implementation is confined to ONE module**, `orpheus/mc/solver.py`
  (640 lines). Measured: `grep -rln "monte_carlo|MonteCarlo" orpheus/` returns
  only `orpheus/mc/__init__.py` (a pure re-export list) and `orpheus/mc/solver.py`;
  the free-flight sampler, roulette, splitting, majorant and delta-tracking
  vocabulary appear nowhere else in `orpheus/`. The only equation whose
  implementer lives outside the package is `ws-pitch` (geometry) and the only
  one with no production site at all is `hetero-tolerance` (a test tolerance).
- **Only 4 `.. (vv-status rationale)` comments exist on this page**, and *none
  of them belongs to one of the 22 equations in scope* — they annotate
  `majorant-no-collision`, `majorant-real-collision`,
  `delta-tracking-collision-pdf` and `mc-analog-multiplication`, all four of
  which are outside this list. So for these 22 the answer came from the
  equation + its surrounding prose + the code block the page itself quotes
  (the page inlines the implementing snippet under almost every equation —
  that quoted snippet IS the authored pointer here, in place of a rationale
  comment).

## ⚠ Finding the parent must weigh BEFORE landing these declarations

**The L0 property tests REPLICATE the solver logic instead of calling it.**
`tests/mc/test_properties.py` carries verbatim comments *"Replicate solver's
chi sampling"* (`:397`), *"Apply roulette (replicate solver logic)"* (`:460`),
*"Replicate solver splitting logic"* (`:524`), *"Replicate solver's batch
statistics"* (`:620`), and `test_periodic_bc_wrapping` (`:420-439`) asserts
`3.7 % 3.6 == 0.1` on bare Python floats — it never runs `_random_walk`.
`[M]` `grep -rn "_random_walk|_russian_roulette|_split_heavy|_precompute_xs" tests/`
returns **0 hits across the whole test tree**: none of the four private kernels
is ever imported or called by any test. `tests/mc/test_gaps.py::test_splitting_copy_count`
(`:640-646`) replicates the splitting rule the same way.

⟹ once these declarations land, a contexts-carrying coverage run will
**REFUTE** the `verifies` claims of the L0 replicating tests against the
declared implementers, unless the test happens to also drive
`solve_monte_carlo` end-to-end (the `tests/mc/test_gaps.py` and
`test_monte_carlo.py` module-level `pytestmark` claims mostly DO, because those
tests call `solve_monte_carlo`). **That refutation is a true finding, not a
declaration error** — the equation really is implemented where I say, and the
L0 test really does not execute it. Do not weaken a declaration to make a
replicating test go green.

---

## `keff-mean`  (22 claims)
- **verdict**: DECLARABLE
- **rationale comment on the page**: none (no `.. (vv-status rationale)` block on
  this equation). The page's own pointer is the prose at
  `monte_carlo.rst:961-962`: *"The final k-effective is the **cumulative mean**
  of M active cycle values"*, plus the named gate
  `test_mc_properties.py::test_batch_statistics_formula`.
- **what the equation says**: the reported keff is the running arithmetic mean of
  the per-cycle keff over the M active cycles.
- **implementers** (complete):
  - `orpheus.mc.solver.solve_monte_carlo` — `orpheus/mc/solver.py:614`
    (`keff_history[ia] = keff_active[:i_active].mean()`; the returned scalar is
    `keff=keff_history[-1]` at `:626`) — the only site in the tree that forms a
    cumulative mean of cycle eigenvalues. `grep` for any shared batch-statistics
    helper in `orpheus/` returns nothing.
- **confidence**: high. `MCResult` is a pure container (no computation); nothing
  else in `orpheus/` averages a keff history. Would change only if a statistics
  helper were extracted out of the orchestrator.

## `sigma-keff`  (22 claims)
- **verdict**: DECLARABLE
- **rationale comment on the page**: none. Pointer is the prose at
  `monte_carlo.rst:969 (prose) / :978-984 (CLT + gates)` (*"The standard deviation of the mean"*, CLT
  `σ_M ~ 1/√M`) and the two named gates
  (`test_mc_convergence.py::test_sigma_scales_with_sqrt_n`,
  `test_mc_properties.py::test_batch_statistics_formula`).
- **what the equation says**: the reported uncertainty is the sample standard
  deviation OF THE MEAN over the M active cycles (the `1/(M(M−1))` normalisation).
- **implementers** (complete):
  - `orpheus.mc.solver.solve_monte_carlo` — `orpheus/mc/solver.py:616-619`
    (`sigma_history[ia] = np.sqrt(((keff_active[:i_active] - keff_history[ia])**2).sum() / (i_active - 1) / i_active)`;
    returned as `sigma=sigma_history[-1]` at `:627`) — the `/(M−1)/M` split is
    literally the equation's `1/(M(M−1))`.
- **confidence**: high. Same argument as `keff-mean`; the two are computed three
  lines apart in one function and share the same sole site.

## `free-flight`  (20 claims)
- **verdict**: DECLARABLE
- **rationale comment on the page**: none on this label. The page inlines the
  implementing line itself (`monte_carlo.rst:410 (the quoted `free_path = ...` line)`: *"In the code:: free_path =
  -np.log(rng.random()) / sig_t_max[ig]"*) — that quoted snippet is the authored
  pointer.
- **what the equation says**: the flight distance is an inverse-CDF exponential
  sample at the MAJORANT rate (not the local total) — `s = −ln ξ / Σ_maj,g`.
- **implementers** (complete):
  - `orpheus.mc.solver._random_walk` — `orpheus/mc/solver.py:407`
    (`free_path = -np.log(rng.random()) / xs.sig_t_max[ig]`) — the sole
    exponential sampler in the tree (`grep -rn "np.log(rng.random())" orpheus/`
    returns this one line).
- **confidence**: high.

## `decompose`  (17 claims)
- **verdict**: DECLARABLE
- **rationale comment on the page**: none on this label (the neighbouring
  derivation steps `majorant-no-collision` / `majorant-real-collision` /
  `delta-tracking-collision-pdf` DO carry rationale comments; `decompose` does not).
- **what the equation says**: at every point the majorant splits into the real
  total plus a fictitious virtual cross section, `Σ_maj = Σ_t(r) + Σ_virtual(r)`.
- **implementers** (complete):
  - `orpheus.mc.solver._random_walk` — `orpheus/mc/solver.py:429-430`. Line 429
    assembles the REAL total from its three parts
    (`sig_t = sig_a + sig_s_sum + sig_2n_sum`) and line 430 solves the
    decomposition for the virtual part (`sig_v = xs.sig_t_max[ig] - sig_t`).
    This is the equation rearranged, and it is the only place `Σ_virtual` exists.
- **confidence**: high. The decomposition is not a standalone helper — it is these
  two lines. Nothing else in the tree names a virtual cross section.

## `scattering-cdf`  (17 claims)
- **verdict**: DECLARABLE
- **rationale comment on the page**: none. Pointer is the inlined snippet at
  `monte_carlo.rst:711-713` plus the `.. note::` on the `Σ_s[from, to]`
  row-vs-column convention (anti-ERR-002).
- **what the equation says**: the exit group is drawn by inverting the cumulative
  distribution of the scattering ROW (from-group fixed), normalised by the row sum.
- **implementers** (complete):
  - `orpheus.mc.solver._random_walk` — `orpheus/mc/solver.py:440-442`
    (`cum_s = np.cumsum(sig_s_row)` ; `ig = np.searchsorted(cum_s, rng.random() * sig_s_sum)`
    ; `ig = min(ig, xs.ng - 1)`), with the row selected at `:426`
    (`sig_s_row = xs.sig_s_dense[mat_id][ig, :]` — the `[from, to]` convention the
    note is about).
  - `orpheus.mc.solver._precompute_xs` — `orpheus/mc/solver.py:369-371`
    (`sig_s_dense[mat_id] = np.array(mix.SigS[0].todense())`) — SUPPORTING: it
    materialises the `[from, to]` P0 matrix the CDF is taken over, and the P0
    (`SigS[0]`) selection is part of what the page documents under this equation.
    Declare only if you want the convention's source in the chain; the CDF itself
    is entirely in `_random_walk`.
- **confidence**: high for `_random_walk`; medium on whether `_precompute_xs`
  should carry the edge (it supplies the matrix, it does not form the CDF).

## `chi-sampling`  (17 claims)
- **verdict**: DECLARABLE
- **rationale comment on the page**: none on this label. Nearby authored pointer:
  the `**Limitation:**` paragraph at `monte_carlo.rst:814-818` — *"The fission
  spectrum χ is taken from a single material (`_any_mat`), not the material at the
  absorption site"* — which names the implementing variable directly, and
  `_any_mat` exists only inside `_precompute_xs`.
- **what the equation says**: the rebirth group is `searchsorted(cumsum(χ), ξ)`,
  i.e. inverse-CDF sampling of the fission spectrum.
- **implementers** (complete — THREE sites, all needed):
  - `orpheus.mc.solver._precompute_xs` — `orpheus/mc/solver.py:372`
    (`chi_cum = np.cumsum(_any_mat.chi)`) — the `cumsum(χ)` half of the equation,
    and the site the page's Limitation paragraph is about.
  - `orpheus.mc.solver._random_walk` — `orpheus/mc/solver.py:456-457`
    (`ig = np.searchsorted(xs.chi_cum, rng.random())` ; `ig = min(ig, xs.ng - 1)`)
    — the `searchsorted` half, applied at fission rebirth after absorption.
  - `orpheus.mc.solver.NeutronBank.initialize` — `orpheus/mc/solver.py:292-296`
    (`group = np.array([np.searchsorted(chi_cum, rng.random()) for _ in range(max_n)])`
    ; `group = np.clip(group, 0, ng - 1)`) — the SAME sampling rule applied to the
    initial source population. ⚠ This is the site a `_random_walk`-only
    declaration would miss.
- **confidence**: high. The `min(..., ng-1)` / `clip(..., ng-1)` guard is the same
  clamp in both consumers, which is what identifies them as one rule with two
  call sites.

## `direction-sampling`  (12 claims)
- **verdict**: DECLARABLE
- **rationale comment on the page**: none. The page inlines the four implementing
  lines at `monte_carlo.rst:536-539` and then spends a whole subsection
  (`mc-direction-sampling-err018`) proving this sampler is NOT isotropic and
  explaining why it is retained (MATLAB compatibility) — i.e. the equation is
  deliberately a description of the code, bias included.
- **what the equation says**: `θ = πξ₁`, `φ = 2πξ₂`, then `Ω_x = sinθ cosφ`,
  `Ω_y = sinθ sinφ` — UNIFORM in θ (the documented ERR-018 non-isotropy), 2-D
  projection only (no `Ω_z`).
- **implementers** (complete):
  - `orpheus.mc.solver._random_walk` — `orpheus/mc/solver.py:410-413`, guarded by
    `if not virtual_collision:` at `:409` (direction is resampled only on a real
    collision — the "virtual collisions preserve direction" property the page
    calls essential).
- **confidence**: high. Sole site; no other direction sampler exists in `orpheus/mc/`.

## `roulette-prob`  (12 claims)
- **verdict**: DECLARABLE
- **rationale comment on the page**: none. Inlined snippet at
  `monte_carlo.rst:847-850`.
- **what the equation says**: the kill probability is `1 − w/w₀`, the shortfall of
  the post-walk weight against the start-of-cycle weight.
- **implementers** (complete):
  - `orpheus.mc.solver._russian_roulette` — `orpheus/mc/solver.py:478-480`
    (`terminate_p = 1.0 - bank.weight[i_n] / weight0[i_n]`, with the `w₀ = 0`
    branch setting `terminate_p = 1.0`) — sole site.
  - `orpheus.mc.solver.NeutronBank.save_start_weights` — `orpheus/mc/solver.py:306`
    (`return self.weight[:self.n].copy()`) — SUPPORTING: it is what makes `w₀`
    the *start-of-cycle* weight (the copy is taken at `solve_monte_carlo:591`,
    before transport). Declare only if you want the `w₀` provenance in the chain.
- **confidence**: high for `_russian_roulette`.

## `roulette-conservation`  (12 claims)
- **verdict**: DECLARABLE (with a stated alternative — read the confidence line)
- **rationale comment on the page**: none. The page frames it as a
  *"Weight Conservation Proof"* subsection and names two statistical gates
  (`test_roulette_weight_conservation`, `test_roulette_restore_weight`).
- **what the equation says**: the roulette is unbiased — `E[w_after] = w_before`,
  because the survivor is restored to `w₀` with probability exactly `w/w₀`.
- **implementers** (complete):
  - `orpheus.mc.solver._russian_roulette` — `orpheus/mc/solver.py:477-484`
    (the whole body: `terminate_p` from the ratio, kill at `:481-482`, restore to
    `weight0[i_n]` at `:483-484`). The identity is the expectation OF THIS
    PROCEDURE and of nothing else; the restore-to-`w₀` (rather than keep-at-`w`)
    is precisely the choice that makes the expectation come out to `w`.
- **NOT an implementer, checked explicitly**: no guard/assert enforces the
  conservation. `grep -n "^\s*assert " orpheus/mc/solver.py` returns **0 hits** —
  the module contains no assertions at all, so there is no runtime enforcement
  site to declare (and under `python -O` there could not be one anyway).
- **confidence**: medium-high. If the project's convention is that an
  *expectation identity* is only "implemented" by something that COMPUTES it,
  then this is `NOTHING:identity` (nothing computes `E[w_after]`). I judge
  DECLARABLE because the equation is the unbiasedness contract of a specific
  procedure and `_russian_roulette` is the unique realisation of that procedure —
  the same relation `splitting-weight-conservation` has to `_split_heavy`. What
  would change it: a project ruling that consequence-of-an-algorithm equations
  take no `implements` edge.

## `keff-cycle`  (12 claims)
- **verdict**: DECLARABLE
- **rationale comment on the page**: none. Inlined snippet at
  `monte_carlo.rst:944`, plus the **Weight normalisation** paragraph immediately
  below (`:948-956`) which inlines the normalisation the denominator depends on.
- **what the equation says**: the per-cycle eigenvalue is the ratio of summed
  end-of-cycle weights to summed NORMALISED start weights.
- **implementers** (complete):
  - `orpheus.mc.solver.solve_monte_carlo` — `orpheus/mc/solver.py:603`
    (`keff_cycle = bank.weight[:bank.n].sum() / weight0.sum()`) — the ratio itself.
    ⚠ Note the numerator sums over `bank.n` AFTER roulette+compaction+splitting,
    which is what the page's *"weights after the random walk, roulette, and
    splitting"* means.
  - `orpheus.mc.solver.NeutronBank.normalize_weights` — `orpheus/mc/solver.py:299-302`
    — SUPPORTING but load-bearing: it is what makes `w⁰` the *normalised* start
    weights the equation stipulates (called at `:590`, one line before the copy).
  - `orpheus.mc.solver.NeutronBank.save_start_weights` — `orpheus/mc/solver.py:306`
    — SUPPORTING: produces the `w⁰` array (`:591`).
- **confidence**: high for `solve_monte_carlo`; medium on including the two
  `NeutronBank` methods (they are genuine parts of the equation's premise —
  breaking either falsifies it — but neither forms the ratio).

## `fission-weight`  (12 claims)
- **verdict**: DECLARABLE
- **rationale comment on the page**: none on this label; but the adjacent
  `mc-analog-multiplication` DOES carry one and it names this label explicitly:
  *"the fission-weight factor it derives is verified by
  `tests/mc/test_properties.py::test_fission_weight_adjustment` (the wired
  `fission-weight` label)"*. The page also inlines the three implementing lines at
  `monte_carlo.rst:748-750`.
- **what the equation says**: analog absorption multiplies the weight by the
  expected fission neutrons per absorption, `w ← w · νΣ_f/Σ_a`.
- **implementers** (complete):
  - `orpheus.mc.solver._random_walk` — `orpheus/mc/solver.py:452-455`
    (`if sig_a > 0: w *= sig_p / sig_a` `else: w = 0.0`), with the numerator
    `sig_p = mat.SigP[ig]` at `:425` and the denominator
    `sig_a = mat.SigF[ig] + mat.SigC[ig] + mat.SigL[ig]` at `:423`. The `else`
    branch is the page's documented **non-fissile** case (`w ← 0`).
- **confidence**: high. Sole site.

## `ws-pitch`  (5 claims)
- **verdict**: DECLARABLE
- **rationale comment on the page**: none, but the page carries an EXPLICIT prose
  pointer immediately under the equation (`monte_carlo.rst:327-328`): *"This is
  the convention used by
  :meth:`~orpheus.geometry.structured_geometry.StructuredGeometry.wigner_seitz_pin_cell`
  (`r_cell = pitch / sqrt(pi)`)."* — the page names its own implementer.
- **what the equation says**: the square unit cell and the Wigner-Seitz cylinder
  must have equal area, so `p = R_cell·√π` (equivalently `R_cell = p/√π`).
- **implementers** (complete, production):
  - `orpheus.geometry.structured_geometry.StructuredGeometry.wigner_seitz_pin_cell`
    — `orpheus/geometry/structured_geometry.py:438` (`r_cell = pitch / np.sqrt(np.pi)`)
    — the site the page names.
  - `orpheus.mc.solver.ConcentricPinCell.default_pwr` — `orpheus/mc/solver.py:92`
    (`radii=[0.9, 1.1, pitch / np.sqrt(np.pi)]`) — the MC package's own application
    of the same conversion: the outermost annulus radius IS the WS radius. ⚠ A
    geometry-only declaration would miss this, and it is the one inside `orpheus/mc/`.
- **⛔ Deliberately NOT declared (cross-page collision — the parent should rule)**:
  the identical math is separately LABELLED on two other pages, so declaring these
  against `ws-pitch` would cross-claim another page's equation:
  - `orpheus.moc.geometry.MOCMesh.__init__` — `orpheus/moc/geometry.py:272`
    (`self.pitch = mesh.edges[-1] * np.sqrt(np.pi)`) → belongs to
    `:eq:`pitch-recovery`` (`method_of_characteristics.rst:159`).
  - `orpheus.derivations.continuous.mms.moc.build_moc_mesh` —
    `orpheus/derivations/continuous/mms/moc.py:409` (`ws_r = P / np.sqrt(np.pi)`)
    → same math as `moc-wigner-seitz` (`method_of_characteristics.rst:145`) /
    `wigner-seitz` (`collision_probability.rst:1449`).
  If the project wants ONE canonical equation for the equal-area conversion, that
  is a page-level decision, not a declaration decision — flagging, not deciding.
- **Test-side sites** (not declared; listed for completeness because they compute
  the equation in the *exact* direction the MC page states it, `p = R√π`):
  `tests/mc/test_monte_carlo.py:128`, `tests/mc/test_cross_verification.py:79`,
  `tests/mc/test_convergence.py:85` and `:122`, `tests/mc/test_gaps.py:678` — five
  copies of `pitch = r_cell * np.sqrt(np.pi)`, i.e. the ERR-017 fix, duplicated per
  test module rather than routed through a helper.
- **confidence**: high for the two production sites. Medium on scope: whether the
  MoC/CP sites should also carry this label depends on a page-ownership ruling.

## `periodic-bc`  (5 claims)
- **verdict**: DECLARABLE
- **rationale comment on the page**: none. Inlined snippet at
  `monte_carlo.rst:612-614` plus the prose *"Python's modulo returns non-negative
  for positive pitch"*.
- **what the equation says**: a particle leaving the square cell is wrapped by
  modulo in both coordinates, `x ← x mod p`, `y ← y mod p`.
- **implementers** (complete):
  - `orpheus.mc.solver._random_walk` — `orpheus/mc/solver.py:417-418`
    (`nx_ = nx_ % pitch` ; `ny_ = ny_ % pitch`) — the only wrapping arithmetic in
    the tree; applied every flight, before the material lookup.
  - `orpheus.mc.solver._mc_bc_periodic` — `orpheus/mc/solver.py:200-206` —
    SUPPORTING / declaration site only. ⚠ It computes NOTHING: it returns the
    string `"periodic"` and its docstring restates the equation
    (*"x' = x mod pitch"*). It is the `BC_REGISTRY` entry that ADMITS the BC
    (`MCMesh.__init__`, `orpheus/mc/solver.py:177-184`, rejects any other kind), not the site that applies
    it. Declare it only if the ledger wants the admission point; the arithmetic is
    exclusively in `_random_walk`.
- **confidence**: high for `_random_walk`.

## `roulette-restore`  (2 claims)
- **verdict**: DECLARABLE
- **rationale comment on the page**: none. Inlined snippet at
  `monte_carlo.rst:847-850`, immediately after the `roulette-prob` block.
- **what the equation says**: the roulette outcome is two-valued — `0` with
  probability `P_kill`, else restored to `w₀` (NOT left at `w`).
- **implementers** (complete):
  - `orpheus.mc.solver._russian_roulette` — `orpheus/mc/solver.py:481-484`
    (`if terminate_p >= rng.random(): bank.weight[i_n] = 0.0` /
    `elif terminate_p > 0: bank.weight[i_n] = weight0[i_n]`). The `elif
    terminate_p > 0` third arm is the documented supercritical case (`w > w₀` ⟹
    weight left UNCHANGED, neither killed nor restored) — the page discusses it at
    `:874-882`, so the same function carries all three outcomes.
- **confidence**: high. Same body as `roulette-prob` / `roulette-conservation`;
  the three labels partition one function (probability / outcome / expectation).

## `collision-estimator`  (1 claim)
- **verdict**: DECLARABLE
- **rationale comment on the page**: none. Inlined snippet at
  `monte_carlo.rst:1008-1010` with the in-code marker
  *"`tally[ig] += w / sig_t   # <-- collision estimator, ERR-024 fix`"* — the page
  quotes the exact line, comment included.
- **what the equation says**: the group flux is accumulated as `Σ w/Σ_t` over
  every REAL collision (scatter, (n,2n) or absorption alike), before the
  three-way branch.
- **implementers** (complete):
  - `orpheus.mc.solver._random_walk` — `orpheus/mc/solver.py:437`
    (`tally[ig] += w / sig_t`), sited inside the `else:` (real-collision) arm at
    `:435` and BEFORE the branch sample at `:438` — the placement is the whole
    point of the equation (ERR-024).
- **⚠ Partial implementation, worth recording**: the code implements only the
  SUM. The `1/(N_act·V)` normalisation in the equation is applied **nowhere** —
  `MCResult.tally` is documented as *"raw scattering tally"* and returned unscaled
  (`solve_monte_carlo:635`); `flux_per_lethargy` divides by `|Δu|` only. So the
  declared implementer realises `Σ w/Σ_t`, not `φ_g`. This is a genuine
  code-vs-equation gap, not a declaration ambiguity.
- **confidence**: high on the site; the gap above is a finding for the parent, not
  a reason to withhold the declaration.

## `mc-lethargy-width-sign`  (1 claim)
- **verdict**: DECLARABLE
- **rationale comment on the page**: none (it sits inside a `.. note::`). The note
  itself is the pointer: *"fixed by taking |Δu_g| at the definition site (the code
  computes `flux_per_lethargy = tally / np.abs(xs.du)`)"* and names the regression
  gate `tests/mc/test_gaps.py::test_flux_per_lethargy_nonnegative`.
- **what the equation says**: with ORPHEUS's fastest-first (descending `eg`) group
  convention, `Δu_g = ln(eg[g+1]/eg[g])` is NEGATIVE — hence the `abs`.
- **implementers** (complete):
  - `orpheus.mc.solver._precompute_xs` — `orpheus/mc/solver.py:360`
    (`du = np.log(eg[1:ng + 1] / eg[:ng])`) — computes exactly the signed quantity
    the equation defines, in exactly the equation's index order.
  - `orpheus.mc.solver.solve_monte_carlo` — `orpheus/mc/solver.py:628`
    (`flux_per_lethargy = tally / np.abs(xs.du) if xs.du is not None else None`)
    — the `|Δu|` the note calls the fix; the equation exists to justify this `abs`.
- **⛔ NOT an implementer — checked and rejected**:
  `orpheus.data.energy_grid.EnergyGrid.lethargy_widths`
  (`orpheus/data/energy_grid.py:190`) computes `np.log(edges[:-1] / edges[1:])`,
  which is the NEGATION of this equation — its docstring says *"Positive by
  construction (edges are descending)"*. Same physics, opposite sign convention;
  declaring it here would assert the equation of a symbol that contradicts it.
- **confidence**: high.

## `splitting-weight-conservation`  (1 claim)
- **verdict**: DECLARABLE (with a stated alternative — read the confidence line)
- **rationale comment on the page**: none. The page frames it as a *"Splitting
  Weight Conservation Proof"* subsection and names two gates
  (`test_mc_properties.py::test_splitting_weight_conservation`,
  `test_mc_gaps.py::test_splitting_copy_count`).
- **what the equation says**: the stochastic rounding is unbiased —
  `E[N] = n(1−f) + (n+1)f = n + f = w`, so the split population carries the
  original weight in expectation (and exactly, since each copy gets `w/N`).
- **implementers** (complete):
  - `orpheus.mc.solver._split_heavy` — `orpheus/mc/solver.py:492-495`
    (`N = int(np.floor(bank.weight[i_n]))` ; `if bank.weight[i_n] - N > rng.random(): N += 1`
    ; `new_w = bank.weight[i_n] / N`). Line 493-494 IS the Bernoulli(`f`) draw whose
    expectation the equation evaluates; line 495 is the exact-conservation half
    (`N·(w/N) = w`).
- **NOT an implementer, checked**: no assert/guard enforces it — `orpheus/mc/solver.py`
  contains **zero** `assert` statements.
- **confidence**: medium-high, for exactly the reason given under
  `roulette-conservation`: nothing COMPUTES `E[N]`; `_split_heavy` is the unique
  procedure the expectation is of. Alternative verdict if the project rules that
  consequence-equations take no edge: `NOTHING:identity`. The two conservation
  labels must be ruled the same way — they are structurally identical.

## `hetero-tolerance`  (1 claim)
- **verdict**: DECLARABLE — but **test-only**, and the parent should read the
  caveat before landing it
- **rationale comment on the page**: none. The prose above it
  (`monte_carlo.rst:345-348`) is the pointer: *"Even with matched areas, the square
  corners contain extra moderator not present in the circle, introducing a ~3–6%
  systematic bias. The heterogeneous test tolerance accounts for this."*
- **what the equation says**: the MC-vs-CP acceptance criterion is
  `|k_MC − k_ref| < 5σ_MC + 0.06·k_ref` — statistical band plus a 6 % geometric
  systematic allowance.
- **implementers** (complete — there is NO production site; grep for `0.06` in
  `orpheus/` finds nothing related):
  - `tests.mc.test_monte_carlo.test_mc_heterogeneous` —
    `tests/mc/test_monte_carlo.py:159`
    (`tol = 5.0 * max(result.sigma, 1e-10) + 0.06 * case.k_inf`)
  - `tests.mc.test_monte_carlo.test_mc_heterogeneous_extended` —
    `tests/mc/test_monte_carlo.py:190` (the same expression, duplicated)
  Both resolve as `py:function:tests.mc.test_monte_carlo.*` nodes.
- **⚠ Caveat the parent must weigh**: `test_mc_heterogeneous` is ALSO the sole
  `verifies("hetero-tolerance")` claimant. Declaring it as its own implementer
  makes that claim adjudicate to *corroborated* trivially ("did the test execute
  itself?"), which is a vacuous green in the ledger. The alternative verdict is
  `NOTHING:definition` — the equation DEFINES an acceptance tolerance whose
  declaration site is the test's own prose. I lean DECLARABLE because the two
  test functions genuinely compute the expression and no other site does, but
  this is the one row in the 22 where I would take a ruling rather than assume.
- **confidence**: medium. What would change it: a project rule on whether a test
  function may be an `implements` source at all.

## `majorant`  (1 claim)
- **verdict**: DECLARABLE
- **rationale comment on the page**: none. Inlined snippet at
  `monte_carlo.rst:387-389 (the quoted majorant loop)` (verbatim the three implementing lines) plus a
  hand-checked value (*"for 2G materials A and B the majorant is [0.60, 2.00]"*)
  and the gate `test_mc_properties.py::test_majorant_computation`.
- **what the equation says**: the Woodcock majorant is the per-group maximum of
  the total cross section over ALL materials, `Σ_maj,g = max_m Σ_t,m,g`.
- **implementers** (complete):
  - `orpheus.mc.solver._precompute_xs` — `orpheus/mc/solver.py:362-364`
    (`sig_t_max = np.zeros(ng)` ; `for mix in materials.values():` ;
    `sig_t_max = np.maximum(sig_t_max, mix.SigT)`) — the running max over the
    material dict, i.e. the equation. Sole site (`grep -rn "majorant" orpheus/`
    returns only this module).
- **confidence**: high. Optional secondary: `orpheus.mc.solver._PrecomputedXS`
  (`orpheus/mc/solver.py:331-347`), the dataclass whose `sig_t_max` field is
  annotated *"(ng,) majorant per group"* — declare only if the ledger wants the
  carrier type as well as the computation.

## `virtual-collision-probability`  (1 claim)
- **verdict**: DECLARABLE
- **rationale comment on the page**: none on this label (its complement
  `majorant-real-collision` has one, and that comment points HERE: *"the
  complement of the wired :eq:`virtual-collision-probability` (P_real = 1 −
  P_virtual)"*). The page inlines the implementing branch at
  `monte_carlo.rst:492-497`.
- **what the equation says**: at a sampled collision site the rejection
  probability is `(Σ_maj,g − Σ_t,g)/Σ_maj,g`.
- **implementers** (complete):
  - `orpheus.mc.solver._random_walk` — `orpheus/mc/solver.py:430-434`
    (`sig_v = xs.sig_t_max[ig] - sig_t` then
    `if sig_v / xs.sig_t_max[ig] >= rng.random(): virtual_collision = True`) —
    the ratio IS the equation, and comparing it to a uniform draw is the
    Bernoulli realisation of it. The `else` arm setting `virtual_collision =
    False` is what makes the direction-preservation property (page: *"Key
    property"*) hold at `:409`.
- **confidence**: high. Note this shares lines 429-430 with `decompose` — the
  decomposition and the probability are one arithmetic step apart, which is
  correct: both labels belong on `_random_walk`.

## `splitting`  (1 claim)
- **verdict**: DECLARABLE
- **rationale comment on the page**: none. The page states the two-case
  distribution and immediately proves its expectation
  (`splitting-weight-conservation`).
- **what the equation says**: the copy count is `⌊w⌋+1` with probability
  `w−⌊w⌋`, else `⌊w⌋` — stochastic rounding of the weight to an integer.
- **implementers** (complete):
  - `orpheus.mc.solver._split_heavy` — `orpheus/mc/solver.py:490-494`
    (`if bank.weight[i_n] > 1.0:` gate at `:491`, `N = int(np.floor(...))` at
    `:492`, `if bank.weight[i_n] - N > rng.random(): N += 1` at `:493-494`) —
    sole site; the `> 1.0` gate is the page's *"Neutrons with weight w > 1"*.
- **confidence**: high.

## `branching`  (1 claim)
- **verdict**: DECLARABLE
- **rationale comment on the page**: none. Two inlined snippets — the Σ assembly
  at `monte_carlo.rst:639-645` and the three-way sample at `:659-670` — plus the
  `Σ_T` / `Σ_2n` counting-convention paragraph.
- **what the equation says**: at a real collision the reaction is drawn from
  `{Σ_s, Σ_2n, Σ_a}/Σ_t`, with `Σ_a = Σ_f+Σ_c+Σ_L`, `Σ_2n = Σ_g' Σ_2n(g→g')`, and
  `Σ_t` the sum of the three.
- **implementers** (complete):
  - `orpheus.mc.solver._random_walk` — `orpheus/mc/solver.py:423-429` (the three
    partial sums and `sig_t = sig_a + sig_s_sum + sig_2n_sum`) and `:438-452`
    (the single scaled uniform `r = rng.random() * sig_t` with the
    `if r < sig_s_sum` / `elif r < sig_s_sum + sig_2n_sum` / `else` cascade — the
    inverse-CDF realisation of the three probabilities).
  - `orpheus.mc.solver._precompute_xs` — `orpheus/mc/solver.py:369-371`
    — SUPPORTING: builds `sig_s_dense` / `sig_2n_dense`, the dense rows whose sums
    are `Σ_s,g` and `Σ_2n,g`. Same call as under `scattering-cdf`; declare
    consistently with that decision.
- **confidence**: high for `_random_walk`.

---

## Machine-readable summary (primary declarations only)

All 14 node ids below were resolved against `.nexus/graph.db`; **0 unresolved**.
Every one is a `function` / `method` / `class` node — no module appears.

| equation | primary implementer(s) |
|---|---|
| `keff-mean` | `orpheus.mc.solver.solve_monte_carlo` |
| `sigma-keff` | `orpheus.mc.solver.solve_monte_carlo` |
| `free-flight` | `orpheus.mc.solver._random_walk` |
| `decompose` | `orpheus.mc.solver._random_walk` |
| `scattering-cdf` | `orpheus.mc.solver._random_walk` (+`_precompute_xs`, supporting) |
| `chi-sampling` | `orpheus.mc.solver._precompute_xs`, `orpheus.mc.solver._random_walk`, `orpheus.mc.solver.NeutronBank.initialize` |
| `direction-sampling` | `orpheus.mc.solver._random_walk` |
| `roulette-prob` | `orpheus.mc.solver._russian_roulette` (+`NeutronBank.save_start_weights`, supporting) |
| `roulette-conservation` | `orpheus.mc.solver._russian_roulette` |
| `keff-cycle` | `orpheus.mc.solver.solve_monte_carlo` (+`NeutronBank.normalize_weights`, `NeutronBank.save_start_weights`, supporting) |
| `fission-weight` | `orpheus.mc.solver._random_walk` |
| `ws-pitch` | `orpheus.geometry.structured_geometry.StructuredGeometry.wigner_seitz_pin_cell`, `orpheus.mc.solver.ConcentricPinCell.default_pwr` |
| `periodic-bc` | `orpheus.mc.solver._random_walk` (+`_mc_bc_periodic`, admission only) |
| `roulette-restore` | `orpheus.mc.solver._russian_roulette` |
| `collision-estimator` | `orpheus.mc.solver._random_walk` |
| `mc-lethargy-width-sign` | `orpheus.mc.solver._precompute_xs`, `orpheus.mc.solver.solve_monte_carlo` |
| `splitting-weight-conservation` | `orpheus.mc.solver._split_heavy` |
| `hetero-tolerance` | `tests.mc.test_monte_carlo.test_mc_heterogeneous`, `tests.mc.test_monte_carlo.test_mc_heterogeneous_extended` ⚠ test-only |
| `majorant` | `orpheus.mc.solver._precompute_xs` |
| `virtual-collision-probability` | `orpheus.mc.solver._random_walk` |
| `splitting` | `orpheus.mc.solver._split_heavy` |
| `branching` | `orpheus.mc.solver._random_walk` (+`_precompute_xs`, supporting) |

## Summary

All 22 equations are **DECLARABLE**; none is a `law` / `identity` /
`canonical-form` / `definition` with no computational site. The MC chapter is
almost entirely sampling rules and estimators, and every one of them has a
concrete site — 20 of the 22 inside the single module `orpheus/mc/solver.py`,
which is why the enumeration closes cleanly: five functions
(`_precompute_xs`, `_random_walk`, `_russian_roulette`, `_split_heavy`,
`solve_monte_carlo`) plus three `NeutronBank` methods carry the whole chapter.
The three rows I could not settle by evidence alone, and which need a project
ruling rather than more searching, are: (1) **`roulette-conservation`** and
(2) **`splitting-weight-conservation`** — both are *expectation identities of a
procedure*; nothing computes `E[w_after]` or `E[N]`, and `orpheus/mc/solver.py`
contains zero `assert`s so no guard enforces them either, so they are DECLARABLE
against `_russian_roulette` / `_split_heavy` only if the ledger accepts "the
unique procedure whose expectation this is" as an implementer (they must be
ruled the same way — they are structurally identical); and (3)
**`hetero-tolerance`**, whose only computational sites are the two test
functions that also *claim* it, so declaring them makes that claim
self-adjudicating — the alternative is `NOTHING:definition`. One scope question
sits outside those three: `ws-pitch` states the same equal-area conversion that
`moc-wigner-seitz`, `pitch-recovery` and `wigner-seitz` state on the MoC and CP
pages, so I declared only the two sites the MC page itself names or owns and
left `MOCMesh.__init__` / `build_moc_mesh` to their own pages. Finally, the
finding that most affects what these declarations will *do*: the L0 gates in
`tests/mc/test_properties.py` (and `test_gaps.py::test_splitting_copy_count`)
**replicate** the solver logic inline rather than importing it, so their
`verifies` claims will adjudicate as REFUTED against the declared implementers —
correctly, and that is a coverage gap to file, not a declaration to soften.
