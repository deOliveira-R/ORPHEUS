# `implements` declaration candidates — 17 equations across 9 theory pages

Derived from the tree at HEAD `a1c90aac` (2026-08-18), branch `main`.
**Does NOT use** `.claude/inventories/implements_declaration_inventory.md` (refuted on
this population). Every symbol below was resolved against `.nexus/graph.db` with the
node-type filter `('function','method','class','data')` — a module is never proposed.

> Note on the count: the brief says 18 equations; the table lists **17** labels, and
> `grep -rn ":label: <name>$" docs/theory/` finds all 17, each exactly once. No label
> is missing; the 18 is an off-by-one in the brief's prose.

---

## `bare-slab-eigenfunction`  (7 claims) — `docs/theory/methods/diffusion_1d.rst:624`
- **verdict**: DECLARABLE
- **rationale comment on the page**: none. (The page's only `.. (vv-status rationale)`
  block is at `:574`, and it is scoped to the **#290 operator-family** labels
  — `diffusion-operator-family`, `diffusion-removal-xs`, `diffusion-scalar-composite`,
  `diffusion-partial-current-dictionary`, `diffusion-albedo-law`,
  `diffusion-interior-conductance`, `diffusion-boundary-closure`. It closes with
  "These sentinels are co-located here, with their #290-family siblings" and does
  **not** name any `bare-slab-*` label. The bare-slab section that follows carries no
  sentinel of its own.)
- **what the equation says**: the fundamental spatial mode of a zero-flux bare slab is
  `φ(x) = sin(πx/L)` — group-independent, so multigroup keeps the same shape.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.derivations.continuous.cases.diffusion.derive_1rg_continuous` —
    `orpheus/derivations/continuous/cases/diffusion.py:721` — builds and **returns the
    callable** `phi(x, g) = phi_spectrum[g] * np.sin(np.pi * x / L)` (`:750`); its own
    docstring says "See :eq:`bare-slab-eigenfunction`, :eq:`bare-slab-buckling`, and
    :eq:`bare-slab-critical-equation`" and its `equation_labels` tuple already carries
    the label (`:787`). This is the sine, materialised.
  - `orpheus.derivations.continuous.cases.diffusion.derive_1rg` —
    `orpheus/derivations/continuous/cases/diffusion.py:99` — the `VerificationCase`
    sibling; carries `"bare-slab-eigenfunction"` in `equation_labels` (`:138`) and
    relies on the sine mode to reduce the PDE to the 2x2 algebraic problem it solves.
  - `orpheus.derivations.continuous.cases.diffusion._bare_slab_spectrum` —
    `orpheus/derivations/continuous/cases/diffusion.py:684` — the kernel both of the
    above call. Its docstring states the sine mode as the *premise* of the reduction
    ("The fundamental spatial mode is `sin(πx/L)`, so the PDE reduces to the 2x2
    algebraic eigenvalue problem") and it produces `phi_spectrum`, the `c_g` the sine
    is scaled by.
- **confidence**: high. The label is already spelled in two `equation_labels` tuples
  and in a docstring `:eq:` reference; `np.sin(np.pi * x / L)` appears exactly once in
  `orpheus/derivations/continuous/cases/diffusion.py`, and the other tree-wide
  `np.sin(np.pi * x / L)` hits are all SN/MoC **MMS ansätze** (different equations).
  What would change it: if `_bare_slab_spectrum` is judged too indirect (it consumes
  the mode rather than evaluating it), drop it and keep the two public functions.

---

## `bare-slab-buckling`  (7 claims) — `docs/theory/methods/diffusion_1d.rst:631`
- **verdict**: DECLARABLE
- **rationale comment on the page**: none (see the note under
  `bare-slab-eigenfunction` — the `:574` sentinel block is #290-family scoped).
- **what the equation says**: the geometric buckling of the fundamental sine mode on a
  slab of thickness `L` is `B² = (π/L)²`.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.derivations.continuous.cases.diffusion._bare_slab_spectrum` —
    `orpheus/derivations/continuous/cases/diffusion.py:706` — literally
    `B2 = (np.pi / L) ** 2`, then used as `D * B2` on the diagonal of `A`.
  - `orpheus.derivations.continuous.cases.diffusion.derive_1rg` —
    `orpheus/derivations/continuous/cases/diffusion.py:103` — its **own** second copy,
    `B2 = (np.pi / fuel_height) ** 2`, and it carries `"bare-slab-buckling"` in
    `equation_labels` (`:137`). (Independent evaluation, not a call into
    `_bare_slab_spectrum` — worth knowing when declaring: two sites compute it.)
  - `orpheus.derivations.continuous.cases.diffusion.derive_1rg_continuous` —
    `orpheus/derivations/continuous/cases/diffusion.py:721` — carries the label at
    `:786` and reaches the value through `_bare_slab_spectrum`; its `Provenance`
    `derivation_notes` state `B^2 = (pi/L)^2` verbatim.
- **confidence**: high. `(np.pi / <length>) ** 2` occurs exactly twice in `orpheus/`,
  both in this file, both listed.

---

## `bare-slab-critical-equation`  (7 claims) — `docs/theory/methods/diffusion_1d.rst:639`
- **verdict**: DECLARABLE
- **rationale comment on the page**: none (see note above).
- **what the equation says**: substituting the sine mode into the diffusion equation
  gives the eigenvalue condition `D·B² + Σ_r = (1/k)·νΣ_f` — i.e. the buckled removal
  operator balances the scaled production operator.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.derivations.continuous.cases.diffusion._bare_slab_spectrum` —
    `orpheus/derivations/continuous/cases/diffusion.py:707-709` — assembles the
    multigroup matrix form of exactly this condition:
    `A = np.diag(D * B2 + absorption + scat_out) - _downscatter_matrix(...)` (that
    diagonal **is** `D·B² + Σ_r`), `F = np.outer(chi, production)` (`νΣ_f` with its
    emission spectrum), then `M = np.linalg.solve(A, F)` — i.e. `A φ = (1/k) F φ`.
  - `orpheus.derivations.continuous.cases.diffusion.derive_1rg` —
    `orpheus/derivations/continuous/cases/diffusion.py:105-107` — the same assembly,
    independently written; carries `"bare-slab-critical-equation"` at `:139`.
  - `orpheus.derivations.continuous.cases.diffusion.derive_1rg_continuous` —
    `orpheus/derivations/continuous/cases/diffusion.py:721` — carries the label at
    `:788`; docstring `:eq:`-cites it.
  - `orpheus.derivations.continuous.cases.diffusion._diffusion_coeffs` —
    `orpheus/derivations/continuous/cases/diffusion.py:94` — supplies the `D` the
    condition is stated in (`D = 1/(3Σ_tr)`). **Weakest of the four**: it implements
    `diffusion-coefficient`, which is a *different* label on the same page; include it
    only if you want the `D` factor attributed. I would NOT declare it here.
- **confidence**: high for the first three; deliberately flagged the fourth as a
  judgement call rather than silently including it.

---

## `bare-slab-keff`  (7 claims) — `docs/theory/methods/diffusion_1d.rst:646`
- **verdict**: DECLARABLE
- **rationale comment on the page**: none on the page. **But there is an authored
  rationale for this exact label in the TEST**, and it is decisive —
  `tests/diffusion/test_continuous_reference.py:61-65`:
  > "The 1G bare-slab continuous reference (derive_1rg_continuous) solves the
  > bare-slab ODE and produces its closed-form keff; its derivation siblings
  > (buckling / eigenfunction / critical-equation) are pinned below, so the governing
  > ODE and its terminal keff are pinned by the same reference"
- **what the equation says**: solving the critical equation for `k` gives the closed
  form `k = νΣ_f / (D·B² + Σ_r)`.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.derivations.continuous.cases.diffusion._bare_slab_spectrum` —
    `orpheus/derivations/continuous/cases/diffusion.py:711-717` — the multigroup
    generalisation of the quotient: `k` is the dominant real eigenvalue of
    `M = A⁻¹F`, which at `ng = 1` **is** `νΣ_f/(D B² + Σ_r)` exactly. Returns
    `(k_val, phi_normalised)`.
  - `orpheus.derivations.continuous.cases.diffusion.derive_1rg_continuous` —
    `orpheus/derivations/continuous/cases/diffusion.py:721` — the function the test's
    rationale names by name; sets `k_eff=k_val` on the returned
    `ContinuousReferenceSolution`.
  - `orpheus.derivations.continuous.cases.diffusion.derive_1rg` —
    `orpheus/derivations/continuous/cases/diffusion.py:99` — computes
    `k_val = float(np.max(np.real(np.linalg.eigvals(M))))` (`:108`) and emits it as
    the `VerificationCase.k_inf` and into the rendered LaTeX.
- **confidence**: high. ⚠ One asymmetry worth landing knowingly: unlike its three
  siblings, `bare-slab-keff` is **absent from every `equation_labels` tuple in the
  code** — the three functions declare buckling / eigenfunction / critical-equation but
  not keff. It is claimed only from the test's `pytestmark`. That is an authoring gap,
  not evidence against the verdict: all three functions demonstrably produce `k`.

---

## `sn-homogenization-bilinear`  (11 claims) — `docs/theory/foundations/frame.rst:1006`
- **verdict**: DECLARABLE
- **rationale comment on the page** (verbatim, `frame.rst:1014-1024` — it is NOT a
  `vv-status` sentinel; the author deliberately wrote a *negative* sentinel):
  > `.. (Wired P6, #281 — no vv-status sentinel.) This bilinear identity —`
  > `   the eigenvalue-consistent (adjoint-weighted) effective cross section`
  > `   that keeps k_eff first-order stationary — is now a VERIFIED solver`
  > `   claim, not documented-only. Solution.homogenize / Solution.condense`
  > `   build the collapse under the ``adjoint=`` parameter, and the`
  > `   full-taxonomy discriminator gates C1 (tests.sn.test_homogenization)`
  > `   and C4 (tests.sn.test_condensation) stack`
  > `   verifies("sn-homogenization-bilinear") against structurally-`
  > `   independent per-region hand rules. The label is covered by tests, so`
  > `   it carries no ``documented`` sentinel.`

  ⭐ This names the two entry points outright. It is the whole answer for the API
  layer; the kernel that actually evaluates the ratio is one hop below.
- **what the equation says**: the eigenvalue-consistent effective cross section is the
  **adjoint-weighted (Petrov-Galerkin) ratio** `Σ_R = ∫_R φ* Σ φ dV / ∫_R φ* φ dV`;
  the functional preserved is the bilinear form `⟨φ*, Σφ⟩`, not a linear rate.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.sn.solution.Solution.homogenize` — `orpheus/sn/solution.py:715` — the
    **spatial** entry named by the rationale. Its `adjoint=` arm type-guards the
    importance (`AdjointSolution`, refusing a forward `Solution`), builds
    `phi / phi_star / rho`, and dispatches to `project_through_bilinear`.
  - `orpheus.sn.solution.Solution.condense` — `orpheus/sn/solution.py:926` — the
    **energy** entry named by the rationale; its `adjoint=` arm builds the
    representative importance spectrum and passes it as
    `Mixture.condense(..., adjoint_spectrum=...)`.
  - `orpheus.transport.mesh.material_xs_field.MaterialXSField.project_through_bilinear`
    — `orpheus/transport/mesh/material_xs_field.py:407` — **the production kernel**.
    Constructs the `PetrovGalerkinFrame`s whose test bases carry `phi_star * phi`
    (pair frame) and `rho` (collision frame), and projects every channel. This is
    `Σ_R = ⟨φ*, Σφ⟩/⟨φ*, φ⟩` materialised.
  - `orpheus.data.macro_xs.mixture.Mixture.condense` —
    `orpheus/data/macro_xs/mixture.py:313` — the energy-axis kernel; its
    `adjoint_spectrum=` branch (`:453`) is the bilinear collapse on the group axis
    (data-native, no transport dependency).
  - `orpheus.derivations.common.homogenization.vector_bilinear_rule` —
    `orpheus/derivations/common/homogenization.py:249` — **the algebra of record**.
    Its docstring is the equation, verbatim: "T1: `Σ_{R,g} = Σ_i V_i φ*_{i,g} Σ_{i,g}
    φ_{i,g} / Σ_i V_i φ*_{i,g} φ_{i,g}` — the pair weight φ*⊙φ". The single source
    every theorem and the production rule are welded to.
  - `orpheus.derivations.common.homogenization.collapse_rules` —
    `orpheus/derivations/common/homogenization.py:327` — the `Sum`-form display of the
    same rules that the theory page *renders* (`"vector"` key is this equation), proof-
    welded to the builders by `_display_matches_builder`.
  - `orpheus.numerics.frame.PetrovGalerkinFrame` — `orpheus/numerics/frame.py` — the
    frame TYPE the page argues the bilinear form forces ("ORPHEUS therefore builds it
    as a `PetrovGalerkinFrame` with an explicit flux-weighted test basis, *not* a
    `GalerkinFrame` with a flux-weighted measure"). Include if you want the structure
    attributed as well as the formula; the page's own prose makes the case.
- ⚠ **Adjacent, deliberately NOT included** — the other four channel rules
  (`angular_sigma_t_rule` T1b, `matrix_per_pair_rule` T2,
  `fission_nsf_mixed_fold_rule` / `fission_chi_canonical_rule` T3) are *siblings* of
  T1 under a different label. The page cites `sn-homogenization-adjoint-weighted` for
  the taxonomy (the C1/C4 test classes stack **both** markers). If you are also
  declaring that label, those four belong there, not here.
- **confidence**: high. The rationale names two of the seven by name; the other five
  are reached by reading the named ones' dispatch one hop.

---

## `sn-homogenization-balance-preservation`  (1 claim) — `docs/theory/foundations/frame.rst:779`
- **verdict**: DECLARABLE (but read the caveat — this is the one I would most want a
  second opinion on)
- **rationale comment on the page**: none on this label. The nearest one is at `:755`
  and belongs to the **sibling** label `sn-homogenization-balance` (`:747`):
  > `.. (vv-status rationale) Literature-transcribed definition: the total-XS`
  > `   balance identity every Mixture carries (the same identity as`
  > `   :eq:`sigT-computed`). Its preservation under homogenization is gated`
  > `   by ``test_homogenized_materials_balance`` (which calls`
  > `   ``Mixture.assert_balanced``); the identity itself is a data-layer`
  > `   definition, not an SN solver claim.`

  ⭐ Note the middle clause: the author explicitly says the **preservation** (i.e. THIS
  label) is gated by `test_homogenized_materials_balance` calling
  `Mixture.assert_balanced`. That is the authored pointer, sitting on the neighbour.
- **what the equation says**: after the forward collapse, the coarse cell's balance
  residual equals the `w_{i,g}`-weighted average of the fine residuals — hence zero.
  No "rebalance the homogenized total" step is needed; preservation is automatic
  *because every removal channel shares the one weight* `w_{i,g} = V_i φ_{i,g}`.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.data.macro_xs.mixture.Mixture.balance_residual` —
    `orpheus/data/macro_xs/mixture.py:194` — computes **the LHS of this equation
    verbatim**: `|Σ_t − (Σ_c + Σ_L + Σ_f + rowsum(Σ_s0) + rowsum(Σ_2n))|` per group.
    Evaluated on the homogenized `Mixture`, it *is* the quantity the equation asserts
    to be 0.
  - `orpheus.data.macro_xs.mixture.Mixture.assert_balanced` —
    `orpheus/data/macro_xs/mixture.py:223` — the `raise`-backed enforcement routine
    the page's own rationale names. This is exactly the "ships an `assert_`/`check_`
    routine enforcing it" case the brief says WOULD be an implementer.
  - `orpheus.transport.mesh.material_xs_field.MaterialXSField.project_through` —
    `orpheus/transport/mesh/material_xs_field.py:343` — the forward (two-morphism)
    collapse that MAKES the property true: it is the code that collapses every removal
    channel with the *same* per-`(R, g)` weight, which is the entire content of the
    theorem. `Solution.homogenize`'s docstring states the preservation as a property
    of this path.
  - `orpheus.derivations.common.homogenization.derive_balance_tradeoff` —
    `orpheus/derivations/common/homogenization.py:832` — the **theorem of record** for
    this residual. It symbolically forms `B_{g'} = Σ_{t,R,g'} − Σ_{a,R,g'} −
    Σ_g Σ_{s,R}[g'→g]` and proves it vanishes identically at flat `φ*` (= the forward
    case, which is precisely this label) while being generically non-zero under the
    bilinear rule. Its docstring: "At flat `φ*` the residual vanishes identically —
    the forward collapse enjoys both properties".
- **caveat / what would change my confidence**: an alternative reading is
  `NOTHING:identity` — the equation is a *theorem about* the collapse map rather than a
  formula anything evaluates. I reject that reading because (a) the LHS is literally
  `Mixture.balance_residual`, a shipped computed property; (b) a `raise`-backed checker
  exists and is invoked on the homogenized product; (c) a symbolic derivation of the
  exact residual ships in the algebra of record. If you want the narrowest defensible
  set, take `balance_residual` + `assert_balanced` and drop the other two.
- **confidence**: medium-high on the verdict, high on the four symbols being the right
  candidate pool.

---

## `moc-mms-psi-ref`  (3 claims) — `docs/theory/verification/method_of_characteristics.rst:306`
- **verdict**: DECLARABLE
- **rationale comment on the page**: not on this label directly, but the neighbouring
  sentinel at `:371-378` (owned by `moc-mms-reference-equilibrium`) names it
  explicitly, verbatim:
  > `.. (vv-status rationale) reference-derivation context: the analytically-`
  > `.. computed FSR-average equilibrium of the manufactured solution, the datum`
  > `.. of the MMS reconstruction :eq:`moc-scalar-flux-reconstruction`. It enters`
  > `.. the reconstruction AND the reference, so it cancels in the convergence`
  > `.. error metric — the MMS gate is blind to its sign. The MOC operators the`
  > `.. MMS convergence test verifies are ``moc-mms-psi-ref`` / ``moc-mms-qext``
  > `.. (wired on ``tests/moc/test_mms.py``).`

  ⭐ The last sentence is the author distinguishing the *context* labels (equilibrium,
  reconstruction — sentinel'd `documented`) from the two that are genuinely verified.
  These two are the real MMS operators.
- **what the equation says**: the manufactured scalar-flux ansatz on the pin cell is
  the smooth radial `φ_ref(r) = 1 + A cos(π r²/R²)`, `r` measured from the cell centre,
  `R = P/2`.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.derivations.continuous.mms.moc.MOCPinCellMMSCase.phi_ref` —
    `orpheus/derivations/continuous/mms/moc.py:165` — **the equation, one line**:
    `1.0 + self.amplitude * np.cos(np.pi * self._r2(x, y) / self._R ** 2)`. Its
    docstring is the formula.
  - `orpheus.derivations.continuous.mms.moc.MOCPinCellMMSCase` —
    `orpheus/derivations/continuous/mms/moc.py:107` — the owning case class. It carries
    `equation_labels = ("characteristic-ode", "bar-psi", "boyd-eq-45",
    "moc-mms-psi-ref", "moc-mms-qext")` as a class-level default (`:142-148`) — i.e.
    the tree already declares this label on this class — and it holds `amplitude`,
    `pitch` and the derived `_R`/`_cx`/`_cy`/`_r2` that parameterise the ansatz.
  - `orpheus.derivations.continuous.mms.moc.MOCPinCellMMSCase._r2` —
    `orpheus/derivations/continuous/mms/moc.py:162` — the `r² = (x−P/2)² + (y−P/2)²`
    half of the equation (the choice the page argues for at length: "the `r²` argument
    ensures `C^∞` smoothness at `r = 0`"). Optional if you prefer public symbols only.
  - `orpheus.derivations.continuous.mms.moc._build_moc_mms_continuous_reference` —
    `orpheus/derivations/continuous/mms/moc.py:424` — wraps `phi_ref` into the
    `ContinuousReferenceSolution.phi` callable (a `y = c_y` cut) and forwards
    `equation_labels=mms_case.equation_labels` (`:472`), so the label reaches the
    reference registry through it.
- **confidence**: high. ⚠ Worth knowing before landing: the **module docstring of
  `orpheus/derivations/continuous/mms/moc.py` re-declares `.. math:: :label:
  moc-mms-psi-ref` at `:17`** (and `moc-mms-qext` at `:36`) — a second definition of
  the same label in a docstring. It raises no Sphinx warning today only because that
  module is not `automodule`'d; if the declaration campaign ever renders it, expect a
  duplicate-label collision.

---

## `moc-mms-qext`  (3 claims) — `docs/theory/verification/method_of_characteristics.rst:329`
- **verdict**: DECLARABLE
- **rationale comment on the page**: same as `moc-mms-psi-ref` above — the `:371-378`
  sentinel names `moc-mms-qext` as one of the two "MOC operators the MMS convergence
  test verifies".
- **what the equation says**: substituting the isotropic `ψ_ref = φ_ref/4π` into the
  characteristic ODE gives the **angle-dependent** manufactured source
  `Q_mms(x,y,φ_a,θ_p) = (1/4π)[sinθ_p(cosφ_a ∂_xφ_ref + sinφ_a ∂_yφ_ref) + Σ_t φ_ref]`.
  The streaming half vanishes on a `4π` average, so an isotropic per-FSR source cannot
  carry it — the source must be injected per segment.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.derivations.continuous.mms.moc.mms_sweep` —
    `orpheus/derivations/continuous/mms/moc.py:249` — **evaluates the source per
    segment**, at the segment midpoint, on both the forward and backward passes:
    `streaming = cos_phi*dphi_dx + sin_phi*dphi_dy`;
    `q_seg = inv_4pi * (phi_val + sin_p * streaming / sig_t)`
    — i.e. `Q_mms/Σ_t`, the equilibrium form the exponential closure consumes. This is
    the only site in the tree that assembles the equation.
  - `orpheus.derivations.continuous.mms.moc.MOCPinCellMMSCase.dphi_dx` —
    `orpheus/derivations/continuous/mms/moc.py:169` — the `∂_x φ_ref` factor:
    `−A·(2π/R²)·(x−c_x)·sin(π r²/R²)`, which the page states as its own displayed
    (unlabelled) equation immediately under `moc-mms-qext`.
  - `orpheus.derivations.continuous.mms.moc.MOCPinCellMMSCase.dphi_dy` —
    `orpheus/derivations/continuous/mms/moc.py:177` — the `∂_y φ_ref` factor.
  - `orpheus.derivations.continuous.mms.moc.MOCPinCellMMSCase` —
    `orpheus/derivations/continuous/mms/moc.py:107` — carries `"moc-mms-qext"` in its
    `equation_labels` default (`:147`) and supplies `sigma_t`, the `Σ_t φ_ref` factor.
- **confidence**: high. `mms_sweep` is the unambiguous computing site; the two
  derivative methods are the equation's own named partials.

---

## `cp-escape-from-p-cell`  (3 claims) — `docs/theory/references/peierls_nystrom.rst:6556`
- **verdict**: DECLARABLE
- **rationale comment on the page**: no `.. (vv-status rationale)` sentinel on this
  label. But the prose **immediately below the equation names the code line verbatim**,
  which is stronger than a sentinel:
  > "which is the code line ``P_out = np.maximum(1.0 - P_cell.sum(axis=1), 0.0)`` in
  > all three CP derivation modules."
- **what the equation says**: the flat-source escape probability is recovered from the
  CP matrix itself as the row-sum defect,
  `P_esc,i^CP = 1 − Σ_j P_ij^cell = 1 − (1/(Σ_t,i V_i)) Σ_j rcp_ij^cell`.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.derivations.continuous.flat_source_cp.geometry.build_cp_matrix` —
    `orpheus/derivations/continuous/flat_source_cp/geometry.py:237` — carries **both
    halves** of the equation: the normalisation
    `P_cell[i,:] = rcp[i,:] / (sig_t_g[i] * volumes[i])` (`:362`) and then the named
    line `P_out = np.maximum(1.0 - P_cell.sum(axis=1), 0.0)` (`:365`). This is also the
    SUT of `tests/derivations/test_cp_geometry.py::TestEscapeFromPCell` — a test class
    named after the equation.
  - `orpheus.cp.solver.CPMesh._apply_white_bc` — `orpheus/cp/solver.py:401` — the
    **production** copy of the same line (`:414`), inside the shipped CP solver.
  - `orpheus.cp.solver.CPMesh._normalize_rcp` — `orpheus/cp/solver.py:388` — the
    production half that forms `P_cell = rcp/(Σ_t V)`, i.e. the equation's second
    equality. Include if you want both equalities attributed on the production side;
    `_apply_white_bc` alone covers the headline `1 − Σ_j P_ij`.
  - `orpheus.derivations.continuous.flat_source_cp.slab._slab_cp_matrix` —
    `orpheus/derivations/continuous/flat_source_cp/slab.py:32`
  - `orpheus.derivations.continuous.flat_source_cp.cylinder._cylinder_cp_matrix` —
    `orpheus/derivations/continuous/flat_source_cp/cylinder.py:46`
  - `orpheus.derivations.continuous.flat_source_cp.sphere._sphere_cp_matrix` —
    `orpheus/derivations/continuous/flat_source_cp/sphere.py:45`
    → these three are the "all three CP derivation modules" the page names. ⚠ **They no
    longer contain the line.** Each is now a 5-line pre-bound delegation to
    `build_cp_matrix` (measured: `grep -n "sum(axis=1)"` over
    `orpheus/derivations/continuous/flat_source_cp/` returns exactly **one** hit,
    `geometry.py:365`). I list them because they are the geometry-named public faces the
    page and the reference registry point at, but the honest declaration is
    `build_cp_matrix`; declaring the three wrappers as well is defensible-but-redundant.
- ⛔ **Doc-drift found** (report, do not fix here): the page's clause "in all **three**
  CP derivation modules" is **present-tense false** — the three modules were
  consolidated onto `flat_source_cp/geometry.py`'s single `build_cp_matrix`. The
  sentence should read "in `flat_source_cp/geometry.py` (shared by all three
  geometries) and in `orpheus/cp/solver.py`". Also `:6548-6552` still cites
  `_cylinder_cp_matrix` / `_sphere_cp_matrix` as the white-BC closure sites; those
  symbols do still exist, so no dead reference — only the "three copies" claim is wrong.
- **confidence**: high on `build_cp_matrix` + `CPMesh._apply_white_bc`; medium on
  whether to also declare the three wrappers (a scope preference, not a fact question).

---

## `hebert-3-323`  (2 claims) — `docs/theory/references/peierls_nystrom.rst:3233`
- **verdict**: DECLARABLE
- **rationale comment on the page**: none — this label sits in the "Status /
  keystone fact" close-out section of the rank-N per-face investigation, which carries
  no `vv-status` sentinels at all. The *prose* is the rationale and it is explicit:
  "The same equation is the production F.4 closure on BOTH 1D curvilinear geometries
  shipped in ORPHEUS, with the geometry-specific transmission matrix `W` supplied by
  `compute_hollow_cyl_transmission` … or `compute_hollow_sph_transmission`", and
  "rank-0 scalar `P_SS` … with the scalar geometric series `(1 − P_SS)^{-1}` is exactly
  the `(I − W)^{-1}` we ship at `N = 1`".
- **what the equation says**: the white-BC interface-current closure adds a rank-1
  correction to the collision-probability matrix,
  `P̃ = P + [β⁺/(1 − β⁺ P_SS)] · P_iS · p_Sj^T`; at `β⁺ = 1` this is Stamm'ler Eq. 34
  and, in the rank-0 limit, ORPHEUS's `(I − W)^{-1}` F.4 closure.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.derivations.continuous.peierls_nystrom.geometry._build_white_hebert_op` —
    `orpheus/derivations/continuous/peierls_nystrom/geometry.py:5777` — its own
    docstring's first line is "**Hébert (2009) Eq. (3.323) rank-1 white-BC closure
    operator**". Builds `R = [[1]]`, `T = [[P_ss]]` so that `(I − TR)^{-1}_{00} =
    1/(1 − P_ss)` — the equation's scalar geometric factor, exactly.
  - `orpheus.derivations.continuous.peierls_nystrom.geometry._build_closure_operator_rank2_white`
    — `orpheus/derivations/continuous/peierls_nystrom/geometry.py:5093` — the **F.4
    rank-2 per-face** assembly: `R = (I − W)^{-1} ∈ ℝ^{2×2}`, the matrix form the page
    identifies with 3.323 at `N = 1`.
  - `orpheus.derivations.continuous.peierls_nystrom.geometry._build_white_f4_op` —
    `orpheus/derivations/continuous/peierls_nystrom/geometry.py:5749` — the F.4
    entry point (`build_closure_operator(..., reflection="white")`), which the page
    calls "**ORPHEUS F.4**".
  - `orpheus.derivations.continuous.peierls_nystrom.geometry.build_closure_operator` —
    `orpheus/derivations/continuous/peierls_nystrom/geometry.py:4926` — the dispatcher
    the page names by name (it is the one carrying the `NotImplementedError` guard the
    Status section discusses).
  - `orpheus.derivations.continuous.peierls_nystrom.geometry.compute_hollow_cyl_transmission`
    — `.../geometry.py:4537` — supplies `W` for the cylinder (`Ki₃` Bickley fold,
    `W_oi = (R/r_0) W_io`). Named by the page as the geometry-specific half.
  - `orpheus.derivations.continuous.peierls_nystrom.geometry.compute_hollow_sph_transmission`
    — `.../geometry.py:4663` — supplies `W` for the sphere (bare `exp(−τ)`,
    `W_oi = (R/r_0)² W_io`). Also named by the page.
  - `orpheus.derivations.continuous.peierls_nystrom.cylinder._build_peierls_cylinder_hollow_f4_case`
    — `orpheus/derivations/continuous/peierls_nystrom/cylinder.py:188` — **already
    declares the label**: `equation_labels = ("hebert-3-323", …)` at `:343`, and its
    docstring says "Uses Stamm'ler IV Eq. 34 = Hébert 2009 Eq. 3.323 (see
    :math:numref:`hebert-3-323`)".
  - `orpheus.derivations.continuous.peierls_nystrom.sphere._build_peierls_sphere_hollow_f4_case`
    — `orpheus/derivations/continuous/peierls_nystrom/sphere.py:201` — same, label
    declared at `:339`.
- ⚠ **Considered and NOT included**: `compute_P_ss_cylinder` / `compute_P_ss_sphere`
  (`.../geometry.py:2078`, `:2166`) compute the `P_SS` that appears *inside* the
  equation. They implement `P_SS`, which the page treats as its own quantity with its
  own derivation section; including them would attribute a factor rather than the
  closure. Add only if you want maximal breadth. `BoundaryClosureOperator`
  (`.../geometry.py:4192`) is the carrier TYPE of the `(G, R, P, T)` tuple — the same
  judgement call as `PetrovGalerkinFrame` above.
- **confidence**: high. Two `equation_labels` tuples already name the label, and the
  page names four of the eight symbols in prose.

---

## `real-sh-discrete-orthogonality`  (2 claims) — `docs/theory/foundations/spherical_harmonics.rst:302`
- **verdict**: DECLARABLE — ⭐ **this is the case the brief warns about: it LOOKS like a
  pure orthogonality identity, and ORPHEUS does ship the routine that computes it.**
- **rationale comment on the page**: none on this label (the page's `vv-status`
  sentinels stop at `real-sh-l2plus`, `:288-290`; the orthogonality section that follows
  carries none, correctly — it is test-covered, not `documented`).
- **what the equation says**: on a cubature exact to degree `≥ 2L`, the no-prefactor
  real `Y_ℓ^m` are discretely orthogonal:
  `Σ_n w_n Y_ℓ^m(Ω_n) Y_ℓ'^m'(Ω_n) = (4π/(2ℓ+1)) δ_ℓℓ' δ_mm'`.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis.mass_matrix`
    — `orpheus/numerics/basis/spherical_harmonic_basis.py:239` — **computes the LHS,
    literally**: `Y = self.evaluate(measure.nodes);
    np.einsum("n,nlm,nLM->lmLM", measure.weights, Y, Y)`. Its docstring calls the
    residual "a quadrature-exactness diagnostic". This is the SUT of both verifying
    tests (`test_basis_mass_matrix_against_lebedev`,
    `test_mass_matrix_under_multiple_quadratures`).
  - `orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis.metric_per_ell`
    — `orpheus/numerics/basis/spherical_harmonic_basis.py:172` — **the RHS**:
    `4π/(2ℓ+1)`. Its docstring states exactly this equation's content: "This is the
    THEORETICAL … metric … The DISCRETE counterpart against a quadrature is the diagonal
    of `mass_matrix`; the two agree iff the quadrature is exact for `Y_ℓ^m Y_ℓ'^m'` of
    degree `ℓ + ℓ' ≤ 2L`." That sentence *is* the equation, in prose.
  - `orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis.evaluate` —
    `orpheus/numerics/basis/spherical_harmonic_basis.py:186` — tabulates the
    `Y_ℓ^m(Ω_n)` the sum runs over (the no-prefactor convention the identity is stated
    in). Second-tier: it implements `real-sh-l0/l1/l2plus`, so include only for breadth.
  - `orpheus.numerics.spaces.spherical_harmonic_space.SphericalHarmonicSpace` —
    `orpheus/numerics/spaces/spherical_harmonic_space.py` — carries the `4π/(2ℓ+1)`
    Gram as its `inner_product_weights`; built by
    `SphericalHarmonicSpace.from_L` (`:166`) via `_padded_metric_tensor` (both resolve).
    The page calls this metric "the single source of truth for the SH convention".
- **confidence**: high. ⛔ Do **not** classify this as `NOTHING:identity` — the discrete
  Gram is a shipped computation, not a statement about the algebra.

---

## `pi-r-equals-4pi-i`  (1 claim) — `docs/theory/foundations/spherical_harmonics.rst:318`
- **verdict**: DECLARABLE
- **rationale comment on the page**: none on the label. The prose immediately below is
  the rationale and names every object: "`Π` is the spherical-harmonic frame's
  **analysis face** (``frame.analysis``, `M = Y*W`), `R` is its **reconstruction face**
  (``frame.reconstruction``, with the `(2ℓ+1)` factor), and the `4π` factor comes from
  the no-prefactor convention summing the `4π/(2ℓ+1)` orthogonality with the `(2ℓ+1)`
  reconstruction weight — i.e. the frame is **4π-tight** (frame operator
  `S = T*T = 4π I`)."
- **what the equation says**: analysis-after-reconstruction on the band-limited
  coefficient space is `4π` times the identity — the SH Galerkin frame is `4π`-tight
  (NOT 1-tight; the retired `assert_galerkin_idempotency` checked the wrong `Π R = I`).
- **implementers** (complete list, each verified to resolve):
  - `orpheus.numerics.frame.FrameBase.gram` — `orpheus/numerics/frame.py:254` — **the
    one production site that actually FORMS `M R`**: it "takes a **single**
    `analysis ∘ reconstruction` probe of the all-ones coefficient vector — the **row
    sum** of `M R` — and installs it as the coefficient space's metric". For the SH
    frame (a disjoint-support/orthogonal trial, `gram_structure = DIAGONAL`) that row
    sum *is* the diagonal, i.e. `4π`. This is the equation, computed.
  - `orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis.analyze` —
    `.../spherical_harmonic_basis.py:306` — the `Π = M = Y*W` factor
    (`einsum("n,nlm,n...->lm...", weights, table, values)`).
  - `orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis.reconstruct`
    — `.../spherical_harmonic_basis.py:329` — the `R` factor, carrying the `(2ℓ+1)`
    (`einsum("nlm,l,lm...->n...", table, addition_theorem_factor, coefficients)`).
  - `orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis.addition_theorem_factor`
    — `.../spherical_harmonic_basis.py:162` — **the `4π` itself, factored**: its
    docstring is "Used by the addition-theorem reconstruction `R = (2ℓ+1) Y` and equal
    to `4π · g_C^{-1}`". The `4π` on the RHS of this equation is exactly
    `metric_per_ell × addition_theorem_factor`, and this property is one of the two.
  - `orpheus.numerics.frame._FrameAnalysis` / `orpheus.numerics.frame._FrameReconstruction`
    — `orpheus/numerics/frame.py` — the two operator faces the equation is stated
    between (`frame.analysis`, `frame.reconstruction`), both named by the page. Include
    if you want the arrows attributed as well as the contraction.
  - `orpheus.numerics.frame.GalerkinFrame` — `orpheus/numerics/frame.py` — the frame
    class the verifying test constructs; the page's "the frame is 4π-tight" is a
    property of this object at an SH basis.
- **confidence**: medium-high. High that `FrameBase.gram` + `analyze` + `reconstruct` +
  `addition_theorem_factor` are the right pool; medium on whether `FrameBase.gram`
  should headline, because `gram` is generic over any basis and its SH specialisation is
  what equals `4π`. What would change it: if you would rather attribute only
  SH-specific symbols, drop `FrameBase.gram` and `GalerkinFrame` and keep the three
  `SphericalHarmonicBasis` members.

---

## `hilbert-adjoint-equals-metric-times-S0`  (2 claims) — `docs/theory/foundations/spherical_harmonics.rst:371`
- **verdict**: DECLARABLE
- **rationale comment on the page**: no `.. (vv-status rationale)` block, but there IS
  a **`vv-status` sentinel** at `:407` reading
  `.. vv-status: hilbert-adjoint-equals-metric-times-S0 documented`.
  ⚠ **That sentinel is stale/contradicted**: two L1 tests carry
  `@pytest.mark.verifies("hilbert-adjoint-equals-metric-times-S0")`
  (`tests/numerics/test_spherical_harmonic_space.py:225` and `:348`), so the label is
  *test-covered*, not `documented`. Same authoring situation the `sn-homogenization-
  bilinear` author handled by replacing the sentinel with a "no vv-status sentinel"
  note. Worth flagging to whoever lands the declaration.
  The page's surrounding prose is the real rationale and names the objects:
  "`Π* = g_C · S_0` — the W-weighted Hilbert adjoint, exposed by ``frame.analysis.H``
  and computed generically by the metric-aware ``_AdjointOperator`` wrapper. The
  `g_C^ℓ = 4π/(2ℓ+1)` factor is the SH Gram-matrix diagonal carried on the analysis
  face's codomain … (which IS the canonical `SphericalHarmonicSpace`)."
- **what the equation says**: the Hilbert adjoint of the SH analysis face is the naked
  synthesis pre-scaled by the SH metric:
  `(Π* c)_n = Σ_ℓ (4π/(2ℓ+1)) Σ_m Y_ℓ^m(Ω_n) c_ℓ^m` — i.e. `Π* = g_C · S_0`, distinct
  from the representation transpose `Π^T = w_n · S_0`.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.numerics.operator._AdjointOperator` — `orpheus/numerics/operator.py:1239`
    — the generic wrapper the page names; it composes
    `(1/w_V)·Π^T(w_W·c)` from the domain/codomain metrics, which for the SH frame
    evaluates to `g_C · S_0`. Its `apply` (`orpheus.numerics.operator._AdjointOperator.apply`,
    resolves) is the site if you prefer method granularity.
  - `orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis.analyze_transpose`
    — `.../spherical_harmonic_basis.py:317` — the `Π^T = w_n · S_0` the adjoint is
    built from; its own docstring states this equation's conclusion: "The metric-aware
    ``_AdjointOperator`` combines it with the domain/codomain Gram to give
    `M* = g_C · S_0`, so the Frame's analysis face gets ``.H`` for free."
  - `orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis.synthesize`
    — `.../spherical_harmonic_basis.py:282` — the `S_0` factor
    (`einsum("nlm,lm...->n...", table, coefficients)`), explicitly identified by the
    page as "the bare synthesis (the frame-theory synthesis operator `T*`)".
  - `orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis.metric_per_ell`
    — `.../spherical_harmonic_basis.py:172` — the `g_C^ℓ = 4π/(2ℓ+1)` factor.
  - `orpheus.numerics.spaces.spherical_harmonic_space.SphericalHarmonicSpace` —
    `orpheus/numerics/spaces/spherical_harmonic_space.py` — the codomain that CARRIES
    `g_C` as `inner_product_weights`; the page says it "IS the canonical
    `SphericalHarmonicSpace`", and `_AdjointOperator` reads the metric off it.
  - `orpheus.numerics.frame._FrameAnalysis` — `orpheus/numerics/frame.py` — owns the
    `.H` property (`frame.analysis.H`) that is the equation's LHS as a shipped API.
- **confidence**: high on the symbol pool; the only judgement is granularity
  (`_AdjointOperator` the class vs `.apply` the method).

---

## `bc-single-delivery`  (2 claims) — `docs/theory/foundations/boundary_conditions.rst:2045`
- **verdict**: DECLARABLE
- **rationale comment on the page**: no `.. (vv-status rationale)` block — and the page
  says so on purpose, in a `.. note::` at `:2163-2177`:
  > "`:eq:`bc-single-delivery`` carries **no** ``vv-status`` sentinel because it needs
  > none: it is a genuine L1 equation claim with a committed gate.
  > ``tests/sn/solve/test_declared_inflow_reaches_the_rhs.py``'s
  > ``test_the_declared_boundary_law_holds_on_the_answer`` carries
  > ``@pytest.mark.verifies("bc-single-delivery")`` on both inner-solver parameters,
  > plus ``@pytest.mark.catches("ERR-075")`` …"

  ⭐ An explicit authored statement that this is a real claim, not documentation.
- **what the equation says**: for prescribed inflow (`L = 0`), the converged inflow
  trace equals the declared boundary source **exactly and exactly once**:
  `γ_-ψ|_f = q_f`. Three defects are distinguishable in one reading — doubled delivery
  `2q`, lost channel `0`, wrong linear factor `q + Lγ_+ψ`.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.transport.source_sinks.angular_boundary_source_sink.AngularBoundarySourceSink.from_mesh_laws`
    — `orpheus/transport/source_sinks/angular_boundary_source_sink.py:426` — **the top
    rung**: reads every face's declared `mesh.bc[face].law` and materialises `q` from
    the `PrescribedInflow` ones (all other laws contribute zero, "their entire content
    is the linear factor `L`"). This is the `q_f` of the equation, produced from the
    declaration. Its docstring calls itself "⭐ The DECLARED boundary conditions' `q` —
    the user's path".
  - `orpheus.transport.source_sinks.angular_boundary_source_sink.AngularBoundarySourceSink.from_specs`
    — `.../angular_boundary_source_sink.py:324` — rung 2, evaluates a per-face
    `InflowSourceSpec` at the face's inflow shape.
  - `orpheus.transport.source_sinks.angular_boundary_source_sink.AngularBoundarySourceSink.prescribed_inflow`
    — `.../angular_boundary_source_sink.py:212` — rung 3, **the packing rule itself**:
    writes the inflow ordinate slots and leaves every other slot zero. The page: "the
    packing rule is stated exactly once" — here.
  - `orpheus.transport.source_sinks.angular_boundary_source_sink.AngularBoundarySourceSink`
    — `.../angular_boundary_source_sink.py:166` — the carrier type of `q_f`, if you
    want the object as well as the constructors.
  - `orpheus.sn.solver._build_fixed_source_rhs` — `orpheus/sn/solver.py:3190` — **the
    "once"**. Its docstring: "This is the ONE object that represents a source everywhere
    in the solve; this helper is its single construction point (Cardinal Rule 2 — the SI
    and Krylov inner paths both consume what it returns, rather than each re-deriving
    it)." The single-delivery property is a consequence of this being the single
    construction point, and it is the gate's direct SUT (`from orpheus.sn.solver import
    _build_fixed_source_rhs`, test `:46`).
  - `orpheus.geometry.boundary.prescribed_inflow.PrescribedInflow` —
    `orpheus/geometry/boundary/prescribed_inflow.py:67` — the law the equation is the
    `L = 0` specialisation OF ("the rank-0 case `R = 0` of the affine boundary law").
  - `orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.inflow_indices_for_face`
    — `orpheus/numerics/spaces/angular_trace_space.py:474` — **`γ_-` itself**, the
    inflow-half gather the equation's LHS is. (Its SN-side delegator is
    `SNMethodSpace.inflow_indices_for_face`, `orpheus/sn/mesh/method_space.py:180`;
    note the class is `SNMethodSpace`, not the `SNAngularTraceSpace` a naming guess
    would suggest.) Second-tier — include only if you want the trace selector
    attributed as well as the delivery.
- **confidence**: high for the four `AngularBoundarySourceSink` entries +
  `_build_fixed_source_rhs`; medium on `PrescribedInflow` and the gather (both are
  legitimately implementers of *neighbouring* labels on the same page —
  `bc-affine-law` / the `γ_±` notation label — so declaring them here risks
  double-attribution).

---

## `en-kernel-derivative`  (1 claim) — `docs/theory/verification/reference_solutions.rst:250`
- **verdict**: DECLARABLE
- **rationale comment on the page**: none on this label, but the sentinel on its
  neighbour `en-definition` (`:237-241`) names the evaluator and this identity's role:
  > `.. (vv-status rationale) definition: the canonical A&S 5.1.4 defining integral`
  > `.. of the exponential integral. Definitional — the implementing evaluator`
  > `.. ``kernels.e_n`` is pinned by the derived identities (special values,`
  > `.. derivative, full-line integral) in ``tests.derivations.test_kernels``.`

  ⭐ The author sentinels `en-definition` as `documented` and points at the *derived*
  identities — of which this is one — as the things that carry the verification. So the
  neighbour is documentation and THIS is a claim. The page also names the tests
  outright: "Both engines are exercised against the three identities above in …
  `test_en_derivative_identity` …".
- **what the equation says**: `E_n'(x) = −E_{n−1}(x)` — the exponential integral's
  order-lowering derivative recurrence.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.derivations.common.kernels.e_n_derivative` —
    `orpheus/derivations/common/kernels.py:173` — **the equation, as a function.** Its
    first docstring line is "Analytical derivative :math:`E_n'(x) = -E_{n-1}(x)`", and
    its body is `return float(-mpmath.expint(n - 1, x))`, plus the documented `n == 1`
    fallback `E_1'(x) = −e^{−x}/x` (because `E_0` diverges at 0). Nothing else in the
    tree spells this identity.
  - `orpheus.derivations.common.kernels.e_n` —
    `orpheus/derivations/common/kernels.py:140` — the `E_n` the identity is *about* (and
    the LHS of the test's finite-difference side). Second-tier: it primarily implements
    `en-definition`, which is sentinel'd `documented`. Include only for breadth; I would
    NOT declare it here.
- ⚠ **`e_n_derivative` has NO production consumer** — `grep -rn "e_n_derivative"` over
  `orpheus/` + `tests/` returns the definition, the module docstring's cross-reference,
  and the one test import. Its own docstring is honest about this: "This function exists
  primarily as the RHS of the kernel-identity tests." That does **not** disqualify it —
  it is still the tree's implementation of the identity, and it is a legal `function`
  node — but the declaration will link a label to a symbol whose only caller is a test.
  Flagging so you land it knowingly.
- ⚠ **Claim-count discrepancy**: the brief says 1 claim; `docs/theory/verification/matrix.rst:600`
  reads `` `en-kernel-derivative`, 20 ``. Not resolved here — just noting the two
  sources disagree by 20×, which may matter for prioritisation.
- **confidence**: high.

---

## `e1-decomposition`  (1 claim) — `docs/theory/references/peierls.rst:549`
- **verdict**: DECLARABLE
- **rationale comment on the page**: none on this label. ⚠ And the page's *prose around
  it is a trap*: it introduces the decomposition, then says it "motivates the classical
  **singularity-subtraction** approach (**used in the original implementation**)" and,
  four paragraphs later, "**Why the classical singularity-subtraction scheme failed**"
  + "The unified ``_pg.solve_peierls_*`` adaptive-quadrature K-matrix assembly …
  **is the production path**". Read alone, that reads like "retired — nothing implements
  it". **It is not retired**, it just moved: the decomposition is live in a *different*
  module (the Atkinson product-Nyström `fn_method` path), which this page never
  mentions. The verifying test's docstring is what names the live consumer:
  > "This identity is the foundation of the Atkinson fix … the load-bearing
  > notation-bridge identity used by ``…peierls_atkinson_nystrom`` to split the kernel
  > into log-singular + smooth parts."
  (`tests/derivations/test_path_ai_legacy_plain_gl_signature.py:15-19, 53-70`)
- **what the equation says**: the `E_1` kernel's log singularity is separable —
  `E_1(z) = [−ln z − γ] + R(z)` with `R` analytic and `R(0) = 0` — which is what lets a
  Nyström assembly handle the singular part by product-integration weights and the
  remainder by ordinary quadrature.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.derivations.continuous.fn_method.peierls_atkinson_nystrom.E1_smooth_remainder`
    — `orpheus/derivations/continuous/fn_method/peierls_atkinson_nystrom.py:190` —
    **the `R(z)` half, computed**: `return float(exp1(tau)) + math.log(tau)`, with the
    A&S 5.1.11 Taylor branch `−γ_E + τ` below `τ = 1e-15`.
  - `orpheus.derivations.continuous.fn_method.peierls_atkinson_nystrom.product_simpson_log_weights`
    — `.../peierls_atkinson_nystrom.py:140` — **the `−ln z` half**: the
    product-integration weights against the log-singular piece.
  - `orpheus.derivations.continuous.fn_method.peierls_atkinson_nystrom._simpson_panel_smooth_weights`
    — `.../peierls_atkinson_nystrom.py:213` — the quadrature applied to the `R` half
    ("Standard Simpson weights … against the smooth remainder kernel `(c/2) R(|t − s|)`").
  - `orpheus.derivations.continuous.fn_method.peierls_atkinson_nystrom.build_peierls_operator`
    — `.../peierls_atkinson_nystrom.py:250` — the assembly that **composes the two
    halves** into `K`; the split only exists because this reassembles it.
  - `orpheus.derivations.continuous.fn_method.peierls_atkinson_nystrom._F_k_log_primitives`
    — `.../peierls_atkinson_nystrom.py:90` — the closed-form `∫ s^k log|t−s| ds`
    antiderivatives the product weights are built from. Second-tier.
  - `orpheus.derivations.continuous.fn_method.peierls_atkinson_nystrom._GAMMA_EULER`
    — `.../peierls_atkinson_nystrom.py:87` — resolves as a **`py:data:` node** (legal
    under the widened ontology): the `γ` of the equation, as a module-level constant.
- ⚠ **Convention note, worth knowing before landing**: the page writes
  `R(z) ≡ E_1(z) + ln z + γ` with `R(0) = 0`; the code's `E1_smooth_remainder` returns
  `E_1(τ) + ln τ`, i.e. the page's `R` **minus γ**, with `R_code(0) = −γ_E` (its own
  docstring says so). Both are correct decompositions — the code folds `γ` into the
  smooth part instead of the singular part — but a reader diffing them will see a
  mismatch. The verifying test uses the *page's* split (`R = E1 − (−γ_E − ln τ)`), so
  the two conventions coexist in the tree today.
- ⚠ **Also live, and NOT the same equation**:
  `orpheus/derivations/continuous/peierls_nystrom/slab.py:25-26` restates this
  decomposition **verbatim in its module docstring** and then says "The Nyström method
  uses **singularity subtraction** …". Per the theory page, that description is the
  RETIRED scheme — the module's production path is adaptive `mpmath.quad`
  (`_basis_kernel_weights`), and `slab.py` imports only `e_n_mp`, no log split. So the
  `slab.py` module docstring is **present-tense-false stale prose**, not an implementer.
  Do not declare it. (Report to the archivist — a `coding-standards` present-tense
  falsehood.)
- **confidence**: medium-high. High that the `fn_method` module is where the
  decomposition lives; medium on how many of its six symbols you want, since the
  equation is a *split* and the split is realised across a small family rather than in
  one function. Minimum defensible set: `E1_smooth_remainder` +
  `product_simpson_log_weights`.

---

## `sigT-computed`  (1 claim) — `docs/theory/foundations/cross_section_data.rst:704`
- **verdict**: DECLARABLE
- **rationale comment on the page**: none on this label. Two *other* pages carry the
  rationale for it, both verbatim:
  - `docs/theory/foundations/frame.rst:755-760` — "the total-XS balance identity every
    Mixture carries (**the same identity as** `:eq:`sigT-computed``) … the identity
    itself is a data-layer definition, not an SN solver claim."
  - `docs/theory/foundations/path_integral.rst:771` — "(`:eq:`sigT-computed``); gated by
    ``Mixture.assert_balanced``. A data-layer …"

  The page's own prose gives the reason it is COMPUTED rather than read: "MF=3 MT=1 does
  **not** include upscattering" + "Computing from components ensures self-consistency".
- **what the equation says**: the total cross section is *derived* from its components,
  `σ_t,g(σ_0) = σ_c + σ_f + σ_α + Σ_g' σ_s0,g→g' + Σ_g' σ_2n,g→g'` — never read from
  ENDF MF=3 MT=1.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.data.micro_xs.gendf._build_isotope` —
    `orpheus/data/micro_xs/gendf.py:244` — **the equation exactly as written on the
    page** (microscopic, `σ_0`-resolved). Under the comment
    `# --- Total cross section (computed from components) ---` (`:315`):
    `sigT[i_sig0] = sigC[i_sig0] + sigF[i_sig0] + sigL[i_sig0] + row_sums`, then
    `+= rowsum(sig2)` when `sig2.nnz > 0` (`:319-321`), looped over every `σ_0`
    dilution. This is the site the page's whole section is about, and it is the one a
    name-token heuristic would never find.
  - `orpheus.data.macro_xs.mixture.compute_macro_xs` —
    `orpheus/data/macro_xs/mixture.py:578` — the **macroscopic** twin, one line
    (`:658`): `SigT = SigC + SigL + SigF + rowsum(SigS[0]) + rowsum(Sig2)`. The page's
    `σ_α` is ORPHEUS's `SigL`. `Mixture.balance_residual`'s docstring calls this "**VERBATIM
    the line that derives ``SigT`` in `compute_macro_xs`**".
  - `orpheus.data.macro_xs.mixture.Mixture.balance_residual` —
    `orpheus/data/macro_xs/mixture.py:194` — computes the identity's **defect**
    `|Σ_t − (Σ_c + Σ_L + Σ_f + rowsum(Σ_s0) + rowsum(Σ_2n))|`; its docstring displays
    the equation.
  - `orpheus.data.macro_xs.mixture.Mixture.assert_balanced` —
    `orpheus/data/macro_xs/mixture.py:223` — the `raise`-backed enforcement ("the
    canonical whole-identity check"), which `compute_macro_xs` itself invokes as a
    "free regression guard against a future derivation bug".
- ⛔ **The brief's own worked example is ILLEGAL under the stated ontology.**
  `orpheus.data.macro_xs.mixture.Mixture.SigT` resolves as
  `py:attribute:orpheus.data.macro_xs.mixture.Mixture.SigT` — node type **`attribute`**,
  not `data`. `SigT` is a plain dataclass FIELD (`mixture.py:61`), not a property and
  not a module-level constant, so it is refused exactly like a module would be. The
  ontology widening to `data` admits `py:data:` nodes (module-level constants — e.g.
  `_GAMMA_EULER` above resolves fine); it does not admit `py:attribute:`. If you want
  the carrier attributed, use the CLASS
  `orpheus.data.macro_xs.mixture.Mixture` (resolves, type `class`).
- **confidence**: high. The two computing sites are unambiguous and are the micro/macro
  pair; the two `Mixture` methods are the enforcement half.

---

# Summary

**All 17 equations came back DECLARABLE. Zero are `NOTHING:<kind>`.** That is the
headline, and it is the opposite of the refuted inventory's verdict on this population
— which is consistent with the brief's warning that "nothing implements it" fails
flatteringly. Every one of the 72 proposed symbols was re-verified against
`.nexus/graph.db` at close: **72/72 resolve** with node type in
`('function','method','class','data')`.

**The method the brief prescribed did most of the work, but with a twist worth
recording: for 8 of the 17, the authored rationale that named the answer was on a
NEIGHBOURING label, not on the equation itself.** `bare-slab-keff`'s rationale is in
the *test's* `pytestmark` comment; `moc-mms-psi-ref` / `moc-mms-qext` are named inside
`moc-mms-reference-equilibrium`'s sentinel; `sn-homogenization-balance-preservation`'s
is inside `sn-homogenization-balance`'s; `en-kernel-derivative`'s is inside
`en-definition`'s; `sigT-computed`'s lives on two *other pages*. So the search radius
has to be the surrounding section, not the block above and below the `.. math::`. Three
labels (`sn-homogenization-bilinear`, `bc-single-delivery`,
`hilbert-adjoint-equals-metric-times-S0`) carry an authored note that explicitly
discusses their `vv-status` state — those are the strongest evidence on the page and
they read as prose, not as sentinels.

**Confidence, and where I am least sure.** High on 14. Medium-high on three, and the
uncertainty is about *scope*, never about existence:
- `sn-homogenization-balance-preservation` — an alternative reading is
  `NOTHING:identity` (a theorem about the collapse map). I rejected it because the
  equation's LHS *is* the shipped `Mixture.balance_residual`, a `raise`-backed checker
  is invoked on the homogenized product, and the symbolic residual is derived in
  `derive_balance_tradeoff`. Narrowest defensible set stated in the entry.
- `pi-r-equals-4pi-i` — `FrameBase.gram` is the only site that forms `M R`, but it is
  generic over any basis; whether it should headline over the three SH-specific members
  is a preference I flagged rather than silently resolved.
- `e1-decomposition` — the split is realised across a family of six rather than in one
  function; minimum set given.

**Three findings the declaration campaign should carry forward, none of which I fixed:**
1. ⛔ **The brief's own worked example is illegal under the stated ontology.**
   `Mixture.SigT` resolves as **`py:attribute:`**, not `py:data:` — it is a dataclass
   field. The `data` widening admits module-level constants (verified:
   `_GAMMA_EULER` resolves as `py:data:`) and not class attributes. Any declaration list
   built on "a field is a `data` node" will fail the same way the 11 module rows did.
2. ⛔ **Two present-tense-false doc claims found while deriving.**
   (a) `peierls_nystrom.rst:6560` says the escape line lives "in all three CP derivation
   modules" — measured, it lives in **one** (`flat_source_cp/geometry.py:365`); the
   three geometry modules are now 5-line delegations.
   (b) `peierls_nystrom/slab.py:25-31`'s module docstring still describes
   singularity-subtraction as the method it uses; the theory page records that scheme as
   failed (#113) and the module's production path as adaptive `mpmath.quad`.
   Also: `spherical_harmonics.rst:407` carries
   `.. vv-status: hilbert-adjoint-equals-metric-times-S0 documented` while **two L1
   tests** `verifies` that label — the sentinel contradicts the tree.
3. ⚠ **A duplicate equation label exists.**
   `orpheus/derivations/continuous/mms/moc.py:17` and `:36` re-declare
   `.. math:: :label: moc-mms-psi-ref` / `moc-mms-qext` in the **module docstring**.
   Harmless today only because that module is not `automodule`'d; it will collide the
   moment it is rendered.

**Count note.** The brief says 18 equations; the table lists 17 labels and
`grep -rn ":label: <name>$" docs/theory/` finds all 17, each exactly once. Nothing is
missing. Separately, `en-kernel-derivative` is billed at 1 claim here and at **20** in
`docs/theory/verification/matrix.rst:600` — worth reconciling before prioritising.
