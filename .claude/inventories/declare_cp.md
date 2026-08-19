# `implements` declarations — `docs/theory/methods/collision_probability.rst`

Scope: the 15 equations named in the dispatch brief (207 coverage claims between
them). Derived FRESH from the tree at HEAD `a1c90aac` (2026-08-18); the refuted
`.claude/inventories/implements_declaration_inventory.md` was not consulted.

**Node-resolution status.** Every dotted path listed under **implementers** below
was resolved against `.nexus/graph.db` with the loader snippet from the brief and
returned exactly one node of type `function` / `method` / `class`. No module nodes
are proposed.

**Durable structural fact that governs 13 of the 15 verdicts** (record it once,
here, rather than repeating it per equation):

> The CP method has **exactly two** independent implementations of the flat-source
> collision-probability construction, and they are structurally parallel line for
> line:
>
> 1. **Production** — `orpheus.cp.solver.CPMesh`, split across
>    `_compute_slab_rcp` (slab / E₃, scalar) and `_compute_radial_rcp`
>    (cylinder + sphere, y-quadrature), then `_normalize_rcp`, then the
>    `BC_REGISTRY`-selected closure `_apply_white_bc`. `compute_pinf_group` is the
>    orchestrator that runs those four in order.
> 2. **Derivations reference** — `orpheus.derivations.continuous.flat_source_cp.geometry.build_cp_matrix`,
>    ONE geometry-polymorphic body covering all three geometries (the slab is the
>    degenerate 1-point y-rule at `y=0, w=0.25`), with the rcp assembly, the
>    normalisation, and the white-BC closure all **inlined in that single
>    function**.
>
> ⚠ `slab._slab_cp_matrix`, `cylinder._cylinder_cp_matrix`, `sphere._sphere_cp_matrix`
> are **thin facades** post-Phase-B.2: each pre-binds a `FlatSourceCPGeometry`
> singleton and returns `build_cp_matrix(...)`. They contain **no arithmetic**. I do
> NOT list them as implementers (see "Refuted candidates" at the end for the
> reasoning) — but note the theory page still cites them as if they carried the
> formula, which is stale prose (documented at the end).

---

## `optical-path`  (23 claims)  — page line 168

- **verdict**: DECLARABLE
- **rationale comment on the page**: none (no `.. (vv-status rationale)` block within
  ±30 lines; the nearest one on the page is at line 1108 and belongs to
  `self-collision-probability-slab`). The authored evidence is instead the module-level
  `pytestmark = pytest.mark.verifies(...)` in `tests/cp/test_verification.py:108`, whose
  own comment says the label is "exercised by the per-geometry CP eigenvalue tests" as
  part of "the integral-transport derivation chain underlying every CP matrix
  computation".
- **what the equation says**: the optical path between two points is the line integral
  of Σ_t along the connecting ray — i.e. `τ = ∫ Σ_t ds`, geometric length converted to
  mean free paths.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.cp.solver.CPMesh._compute_slab_rcp` — `orpheus/cp/solver.py:266,269` —
    the slab realization: `tau = sig_t_g * t` (piecewise-constant Σ_t × region width)
    followed by `bnd_pos = np.concatenate(([0.0], np.cumsum(tau)))`, which is the
    discrete `∫ Σ_t ds` from the reflective centre.
  - `orpheus.cp.solver.CPMesh._compute_radial_rcp` — `orpheus/cp/solver.py:338,341-343` —
    the curvilinear realization on the y-chord family: `tau = sig_t_g[:, None] * chords`,
    then the cumulative `bnd_pos` recurrence. This is the same integral, one per impact
    parameter.
  - `orpheus.derivations.continuous.flat_source_cp.geometry.build_cp_matrix` —
    `orpheus/derivations/continuous/flat_source_cp/geometry.py:311,313-315` — the
    independent derivations realization of the same two lines, for all three geometries.
  - `orpheus.derivations.common.kernels.chord_half_lengths` —
    `orpheus/derivations/common/kernels.py:52` — **the `ds` half of the integral**: it
    computes the per-region geometric chord length `ℓ_k(y)` that Σ_t multiplies. It does
    not itself apply Σ_t, so declaring it is a judgement call; I include it because
    without it the equation's `∫ ... ds` measure has no computational site at all in
    curvilinear geometry, and it is a genuine primitive rather than a wrapper.
    ⚠ If you prefer a strict reading ("the implementer is the site that produces τ"),
    drop this one and keep the three above — that reading is also defensible and loses
    nothing, because all three τ sites already resolve.
- **confidence**: high for the three τ-producing sites; medium for `chord_half_lengths`.
  What would change it: a ruling on whether a geometric-measure primitive counts as
  implementing an equation whose integrand it supplies but does not evaluate.
- **scoping note** (checked, not assumed): `optical-path` is **not** a shared label.
  `grep -rn "optical-path" docs/theory/` outside this page returns only prose in
  `references/peierls_nystrom.rst` and the generated `verification/matrix.rst`. MoC's
  own optical path is a **separate label**, `optical-thickness`
  (`docs/theory/methods/method_of_characteristics.rst:238`), so MoC/SN τ sites belong to
  their own equations and are correctly excluded here.

---

## `complementarity`  (4 claims)  — page line 221

- **verdict**: DECLARABLE
- **rationale comment on the page**: none as a `.. (vv-status rationale)` block, but the
  page names the code inline and unambiguously: *"In the code:
  ``P_out = 1 - P_cell.sum(axis=1)`` (:meth:`CPMesh._apply_white_bc`). Verified by
  ``tests/cp/test_properties.py::test_row_sums`` for all three coordinate systems."*
  The test module's own `pytestmark` comment (`tests/cp/test_properties.py:34`) repeats
  it: `"complementarity",  # test_row_sums: sum_j P_ij + P_out = 1`.
- **what the equation says**: within one cell, the within-cell collision probabilities
  out of region `i` plus the escape probability sum to one — so the escape probability
  is *defined* by the row deficit.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.cp.solver.CPMesh._apply_white_bc` — `orpheus/cp/solver.py:414` —
    `P_out = np.maximum(1.0 - P_cell.sum(axis=1), 0.0)`. This is the equation solved for
    `P_i,out`; it is the only site in production that produces `P_out`.
  - `orpheus.derivations.continuous.flat_source_cp.geometry.build_cp_matrix` —
    `orpheus/derivations/continuous/flat_source_cp/geometry.py:365` — the identical line
    in the independent derivations path.
- **confidence**: high. The equation is not merely an identity the algebra satisfies —
  it is the *constitutive definition* of `P_out` in both implementations, so there is a
  literal computational site. What would change it: nothing short of `P_out` being
  computed some other way (it is not; `grep "P_out = " orpheus/` returns exactly these
  two lines plus one unrelated `P_out_n` in the Peierls Nyström geometry).

---

## `reciprocity`  (4 claims)  — page line 237

- **verdict**: DECLARABLE
- **rationale comment on the page**: none as a `.. (vv-status rationale)` block. The
  authored knowledge is (a) the page's own "Practical consequence" paragraph naming
  `:meth:`CPMesh._normalize_rcp``, and (b) the `pytestmark` comment at
  `tests/cp/test_properties.py:35`:
  `"reciprocity",      # test_reciprocity: Sigma_t V P_ij = Sigma_t V P_ji`.
- **what the equation says**: detailed balance — the reduced quantity
  `r_ij ≡ Σ_t,i V_i P_ij^cell` is **symmetric**, because a chord traversed i→j is the
  same chord traversed j→i.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.cp.solver.CPMesh._compute_slab_rcp` — `orpheus/cp/solver.py:259-306` — this
    method's whole product is the array `rcp[i,j] = Σ_t,i V_i P_ij`, i.e. **the left-hand
    side of the equation as a first-class object**. It is symmetric *by construction*:
    `gap_d` swaps `(i,j)→(j,i)` consistently (`:287` vs `:289`), `gap_c = bnd_pos[i] +
    bnd_pos[j]` (`:301`) is symmetric, and the four-term difference `dd`/`dc` is symmetric
    in `(tau_i, tau_j)`. A sign or index error here is what breaks reciprocity.
  - `orpheus.cp.solver.CPMesh._compute_radial_rcp` — `orpheus/cp/solver.py:309-386` — the
    same symmetric `rcp` assembly for cylinder/sphere (`:364/:366` and `:379`).
  - `orpheus.cp.solver.CPMesh._normalize_rcp` — `orpheus/cp/solver.py:388-399` — the map
    `P_cell[i,:] = rcp[i,:] / (Σ_t,i V_i)` that converts the symmetric reduced object into
    the asymmetric probability matrix. This is the site the page names and the site
    `test_reciprocity` pins.
  - `orpheus.derivations.continuous.flat_source_cp.geometry.build_cp_matrix` —
    `orpheus/derivations/continuous/flat_source_cp/geometry.py:316-363` — the independent
    derivations realization of both halves (symmetric rcp assembly + row normalisation).
- **confidence**: high that these are the sites; medium on the *framing*. Positive
  evidence that this is not a bare law: the symmetric reduced object `rcp` is
  **materialized as an array** in three functions, so the equation's left-hand side has a
  literal computational site rather than being only an emergent property. What would
  change it: a project ruling that "an identity satisfied by construction, with no
  `assert_`/`check_` guard in `orpheus/`, is a `law`". ⚠ I checked for such a guard:
  `grep -rni "reciprocit" orpheus/ --include=*.py` returns **zero** hits in `orpheus/cp/`
  and zero CP-relevant guards anywhere (the hits are Peierls-Nyström
  `CurvilinearGeometry.reciprocity_factor`, an unrelated surface-area ratio, plus SN
  adjoint prose). Enforcement lives only in `tests/cp/test_properties.py::test_reciprocity`.

---

## `reciprocity-lower-triangle`  (4 claims)  — page line 254

- **verdict**: DECLARABLE
- **rationale comment on the page**: none as a `.. (vv-status rationale)` block. The
  page says: *"In the code, :meth:`CPMesh._normalize_rcp` divides the reduced collision
  probability by Σ_t,i V_i for each row."* And the test module states it as a pinning
  claim (`tests/cp/test_properties.py:36-38`): *"test_reciprocity also pins the
  rearranged lower-triangle form the code uses to fill P_ji (**a factor error in
  `_normalize_rcp` breaks the identity**)"*.
- **what the equation says**: the rearranged reciprocity — `P_ji = P_ij · (Σ_t,i V_i)/(Σ_t,j V_j)`
  — so the lower triangle is determined by the upper one.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.cp.solver.CPMesh._normalize_rcp` — `orpheus/cp/solver.py:388-399` — the
    per-row division `rcp[i,:] / (Σ_t,i V_i)` applied to a symmetric `rcp` is *exactly*
    this ratio map; it is the only site whose factor can break the identity, and it is
    the one both the page and the test name.
  - `orpheus.derivations.continuous.flat_source_cp.geometry.build_cp_matrix` —
    `orpheus/derivations/continuous/flat_source_cp/geometry.py:360-362` — the inlined
    equivalent in the derivations path.
- ⛔ **Doc-drift finding, please fix while landing this**: the page's own sentence
  *"The CP matrix need only be computed for j ≥ i; the lower triangle follows from
  [this equation]. **This halves the computation cost.**"* is **present-tense false**.
  Both `_compute_slab_rcp` (`orpheus/cp/solver.py:273,282` — `for i in range(N): ... for
  j in range(N):`) and `_compute_radial_rcp` (`:347,359`) and `build_cp_matrix`
  (`geometry.py:319,331`) run the **full N² double loop** and compute both triangles
  explicitly. No lower-triangle fill-in exists anywhere in the tree. The equation is
  therefore realized as an **identity the produced matrix satisfies**, not as an
  implemented shortcut — which is precisely why `_normalize_rcp` (the factor) and not a
  triangle-fill routine is the implementer.
- **confidence**: high on the implementer; the doc-drift above is the thing that would
  otherwise have made me hunt for a nonexistent "fill lower triangle" function.

---

## `surface-to-region`  (23 claims)  — page line 277

- **verdict**: DECLARABLE
- **rationale comment on the page**: none as a `.. (vv-status rationale)` block, but the
  page quotes the code verbatim as a literal block: *"ORPHEUS uses this convention in the
  white-BC transform::  `# White-BC closure (geometry-agnostic)` / `P_in = sig_t * V *
  P_out / S_cell`"*, and names the surface-area accessor
  (``mesh.areas[-1]``, :func:`~geometry.coord.compute_areas_1d`). The
  `tests/cp/test_verification.py` `pytestmark` comment adds that this label is one of the
  "B.2 additions (issue #87)" whose correctness is implied by a converged `k` matching the
  analytical CP eigenvalue to 1e-8.
- **what the equation says**: the probability that a neutron entering the cell isotropically
  through the surface has its first collision in region `j`, obtained by reciprocity against
  the volume source: `P_in,j = Σ_t,j V_j P_j,out / S`.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.cp.solver.CPMesh._apply_white_bc` — `orpheus/cp/solver.py:415` —
    `P_in = sig_t_g * V * P_out / S_cell`, with `S_cell = self.mesh.areas[-1]` at `:412`.
    Character-for-character the equation.
  - `orpheus.derivations.continuous.flat_source_cp.geometry.build_cp_matrix` —
    `orpheus/derivations/continuous/flat_source_cp/geometry.py:366-367` — the identical
    line, with `S_cell = geom.surface_area(R_cell)`.
  - `orpheus.derivations.continuous.flat_source_cp.geometry.FlatSourceCPGeometry.surface_area` —
    `orpheus/derivations/continuous/flat_source_cp/geometry.py:202-215` — supplies the `S`
    of the equation per geometry (`1`, `2πR`, `4πR²`), i.e. the three-row table the page
    prints directly under the equation. It is the only site in `derivations/` that decides
    `S`, so a wrong `S` breaks this equation and nothing else catches it.
  - `orpheus.geometry.coord.compute_areas_1d` — `orpheus/geometry/coord.py:72` —
    **verified to resolve** (`py:function:orpheus.geometry.coord.compute_areas_1d`). This
    is the production counterpart of `surface_area`: the page names it explicitly, and
    `mesh.areas[-1]` (read at `orpheus/cp/solver.py:412`) is where production gets `S`.
    ⚠ Scoping caveat, flagged not decided: it is a *general-geometry* primitive also
    consumed by SN/MoC/diffusion, so declaring it against a CP-page equation widens the
    equation's reach beyond CP. If you prefer CP-local declarations only, the three above
    are complete for `orpheus/cp/` + `orpheus/derivations/` — but then production's `S`
    has no declared site while the derivations' `S` does, which is an asymmetry worth
    making deliberately rather than by omission.
- **confidence**: high for the three listed; the fourth is deliberately left as a decision.

---

## `surface-to-surface`  (23 claims)  — page line 328

- **verdict**: DECLARABLE
- **rationale comment on the page**: none as a `.. (vv-status rationale)` block. The page's
  own pointer is the paragraph immediately below: *"The same formula appears in all three
  derivation scripts (e.g., `orpheus/derivations/continuous/flat_source_cp/slab.py`, line
  `P_in = sig_t_g * t_arr * P_out` … and `…/cylinder.py`, line `S_cell = 2.0 * np.pi *
  r_cell`)"*. ⛔ **That citation is stale** — see the doc-drift section at the end; those
  lines moved into `geometry.build_cp_matrix` at Phase B.2 and no longer exist in
  `slab.py`/`cylinder.py`.
- **what the equation says**: whatever fraction of the isotropic inward surface source does
  not collide anywhere inside the cell leaves again — `P_in,out = 1 − Σ_j P_in,j`.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.cp.solver.CPMesh._apply_white_bc` — `orpheus/cp/solver.py:416` —
    `P_inout = max(1.0 - P_in.sum(), 0.0)`.
  - `orpheus.derivations.continuous.flat_source_cp.geometry.build_cp_matrix` —
    `orpheus/derivations/continuous/flat_source_cp/geometry.py:368` — identical.
- **confidence**: high. The quantity is materialized and named in both paths and exists
  nowhere else: `grep -rn "P_inout" orpheus/ --include=*.py` returns **8** lines total —
  per file, one definition, one `if P_inout < 1.0` guard, one use in the `p-inf` closure,
  plus (production only) two docstring restatements at `solver.py:406,431`.

---

## `p-inf`  (25 claims)  — page line 346

- **verdict**: DECLARABLE
- **rationale comment on the page**: none as a `.. (vv-status rationale)` block. Authored
  pointer in the prose immediately below the equation: *"This formula is **identical for
  all three geometries** when expressed in terms of V_i and S. It is implemented in the
  white-BC transform selected from :attr:`CPMesh.BC_REGISTRY` and independently in all
  three derivation scripts."* Plus `tests/cp/test_verification.py:86` and
  `tests/cp/test_slab.py:16` and `tests/cp/test_peierls_flux.py:31` carry it in their
  `verifies(...)` lists, and `orpheus/derivations/continuous/flat_source_cp/slab.py:131`
  names it in the verification case's own `labels` list.
- **what the equation says**: the infinite-lattice CP matrix — the within-cell matrix plus
  the rank-one escape/re-entry outer product, resummed over infinitely many lattice
  traversals.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.cp.solver.CPMesh._apply_white_bc` — `orpheus/cp/solver.py:401-424`, the
    formula at `:418-420` (`P_inf = P_cell.copy()` … `P_inf += np.outer(P_out, P_in) /
    (1.0 - P_inout)`) — the production implementation.
  - `orpheus.derivations.continuous.flat_source_cp.geometry.build_cp_matrix` —
    `orpheus/derivations/continuous/flat_source_cp/geometry.py:370-373` — the independent
    derivations implementation (this is the line the theory page quotes, though it
    attributes it to `slab.py`).
  - `orpheus.cp.solver.CPMesh.compute_pinf_group` — `orpheus/cp/solver.py:237-257` — the
    orchestrator that IS `P^∞` as a produced object: it dispatches the rcp assembly by
    coordinate system, applies `_normalize_rcp`, and returns `self._bc_transform(P_cell,
    sig_t_g)`. Include it: it is the public name the whole architecture section of the page
    uses for this equation, and a dispatch error in it (slab routed through the radial
    kernel) produces a wrong `P^∞` with every sub-step individually correct.
- **deliberately EXCLUDED, with reason** (so a later audit does not re-derive it):
  - `orpheus.cp.solver._cp_bc_white` — `orpheus/cp/solver.py:427-435` — a **registry
    factory**; its body is `return cp_mesh._apply_white_bc`. It restates the formula in its
    docstring but computes nothing. Declaring it would attribute the equation to a
    dispatch table entry.
  - `orpheus.cp.solver._cp_bc_vacuum` — `:437-443` — implements `P_inf = P_cell`, i.e. the
    *absence* of the white-BC closure; that is the complementary boundary condition, not
    this equation.
  - `_slab_cp_matrix` / `_cylinder_cp_matrix` / `_sphere_cp_matrix` — thin facades, no
    arithmetic (see the header block).
- **confidence**: high.

---

## `cp-infinite-lattice-sum`  (4 claims)  — page line 512

- **verdict**: DECLARABLE — **same implementers as `p-inf`**, and this is the point of the
  equation rather than a defect: it is the *derivation* of `p-inf` (the geometric series
  `1 + P_in,out + P_in,out² + ⋯`), and its second line **is** `p-inf`.
- **rationale comment on the page**: none as a `.. (vv-status rationale)` block. The
  authored pointer is the `pytestmark` comment at `tests/cp/test_properties.py:39-41`:
  *"test_row_sums builds P_inf via the white-BC geometric series and asserts row sums = 1;
  **a sign flip in the series makes the row sum != 1**"*.
- **what the equation says**: summing the escape→re-enter→escape geometric series in
  closed form gives the `1/(1 − P_in,out)` resummation factor; the series converges because
  `P_in,out < 1`.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.cp.solver.CPMesh._apply_white_bc` — `orpheus/cp/solver.py:418-420` — the
    `/(1.0 - P_inout)` factor **is** the summed series, and the guard `if P_inout < 1.0`
    at `:419` is the convergence condition the equation's closing sentence states.
  - `orpheus.derivations.continuous.flat_source_cp.geometry.build_cp_matrix` —
    `orpheus/derivations/continuous/flat_source_cp/geometry.py:370-372` — same, with the
    same `if P_inout < 1.0` guard at `:371`.
- ⚠ **Judgement note for the declarer**: `p-inf` and `cp-infinite-lattice-sum` will end up
  with an overlapping implementer set. That is correct — they are the same closed form,
  one stated and one derived — and it is why `cp-infinite-lattice-sum` carries only 4
  claims while `p-inf` carries 25. Do not "de-duplicate" by declaring only one; the
  4 claims on `cp-infinite-lattice-sum` are made by `tests/cp/test_properties.py`, a
  different file from `p-inf`'s 25.
- **confidence**: high.

---

## `e3-def`  (25 claims)  — page line 654

- **verdict**: DECLARABLE
- **rationale comment on the page**: none as a `.. (vv-status rationale)` block, but two
  authored pointers, both explicit about the symbols: (a) directly under the equation —
  *"Computed analytically via :func:`scipy.special.expn` (**wrapped as :func:`_e3` in the
  solver and :func:`~orpheus.derivations.common.kernels.e3` in the derivations**)"*; (b)
  the "Derivation source" paragraph at page line ~630 — *"The kernels are implemented in
  `orpheus/derivations/common/kernels.py` (slab, arbitrary-precision) and
  `orpheus/derivations/continuous/flat_source_cp/geometry.py` … Slab E₃:
  :func:`~orpheus.derivations.common.kernels.e3` and
  :func:`~orpheus.derivations.common.kernels.e3_vec` (wrappers over
  :func:`scipy.special.expn`)"*. This is the whole answer, authored.
- **what the equation says**: the slab angular-averaging kernel `E₃(τ) = ∫₀¹ μ e^{−τ/μ} dμ`
  — the third-order exponential integral, with `E₃(0) = 1/2`.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.cp.solver._e3` — `orpheus/cp/solver.py:118-120` —
    `return expn(3, np.maximum(x, 0.0))`; the production slab kernel, the only kernel call
    in `_compute_slab_rcp`.
  - `orpheus.derivations.common.kernels.e3` — `orpheus/derivations/common/kernels.py:110-112` —
    scalar `float(expn(3, max(x, 0.0)))`.
  - `orpheus.derivations.common.kernels.e3_vec` — `orpheus/derivations/common/kernels.py:115-117` —
    vectorised; this is the one `build_cp_matrix` actually reaches.
  - `orpheus.derivations.continuous.flat_source_cp.geometry.FlatSourceCPGeometry.kernel_F3` —
    `orpheus/derivations/continuous/flat_source_cp/geometry.py:165-172` — the `kind ==
    "slab"` branch **is** `E₃`; this is the polymorphic seam through which the derivations
    CP builder selects it, and a wrong branch here silently swaps the kernel.
  - `orpheus.derivations.continuous.flat_source_cp.geometry.FlatSourceCPGeometry.kernel_F3_at_zero` —
    `orpheus/derivations/continuous/flat_source_cp/geometry.py:174-183` — returns `0.5` for
    the slab, i.e. **`E₃(0)`**, the special value the page states in the sentence right
    after the equation and the literal `0.5` in the production self-collision line
    (`orpheus/cp/solver.py:280`). ⚠ Borderline by construction (it is a returned constant,
    not an evaluation of the integral); include it if you want the special value covered,
    drop it if you read `implements` strictly as "evaluates the integrand". I lean include:
    a wrong constant here breaks the equation's stated `E₃(0) = 1/2` with nothing else to
    catch it.
- **deliberately EXCLUDED, with reason**:
  - `orpheus.derivations.common.kernels.e_n` / `e_n_mp` / `e_n_derivative` — these are the
    **general** `E_n` at arbitrary precision (mpmath), and the general definition has its
    **own equation label**: `en-definition` at
    `docs/theory/verification/reference_solutions.rst:232` (with siblings
    `en-kernel-special-values`, `en-kernel-derivative`, `en-kernel-integral`). Declaring
    them here would attach a Peierls-Nyström reference primitive to a CP-page equation and
    duplicate a label that already exists. Their consumers are
    `orpheus/derivations/continuous/peierls_nystrom/slab.py`, not CP.
  - `_ki3_mp` / `_ki3_scaled_cheb` / `_exp_kernel` — the cylinder/sphere kernels; they
    belong to `ki3-def` and the sphere's `e^{−τ}`, not to `e3-def`.
- **confidence**: high.

---

## `s-integral`  (23 claims)  — page line 1007

- **verdict**: DECLARABLE — but read the mechanism note; this is the one where the
  authored text tells you the array is *never materialized*.
- **rationale comment on the page**: none as a `.. (vv-status rationale)` block. The
  authored pointer is an explicit **Implementation** paragraph 20 lines below the
  equation: *"**Implementation.** ORPHEUS does **not** pre-compute the S array. Instead,
  :meth:`CPMesh._compute_radial_rcp` evaluates the second-difference at each y-quadrature
  point and integrates numerically. This avoids storing S and allows the same code path
  for all three geometries by parameterising the kernel function. **The derivation scripts
  use the same approach.**"* That names both implementers.
- **what the equation says**: the impact-parameter integral of the kernel difference
  between the through-centre and same-side optical paths, `S(i,j,g) = ∫₀^{R_i}[Ki₃(τ_p) −
  Ki₃(τ_m)] dy` — the boundary-indexed antiderivative whose four-corner difference
  (:eq:`pcell-from-smat`) gives `P_ij^cell`.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.cp.solver.CPMesh._compute_radial_rcp` — `orpheus/cp/solver.py:309-386`, the
    integration at `:357` and `:384` (`np.dot(y_wts, ...)`) over the four-corner differences
    built at `:371-382` — the named production implementation.
  - `orpheus.derivations.continuous.flat_source_cp.geometry.build_cp_matrix` —
    `orpheus/derivations/continuous/flat_source_cp/geometry.py:331-357` — "the derivation
    scripts use the same approach", now consolidated into this one body.
  - `orpheus.derivations.continuous.flat_source_cp.geometry._second_difference` —
    `orpheus/derivations/continuous/flat_source_cp/geometry.py:71-89` — the free function
    that computes the four-term combination `F(g) − F(g+τ_i) − F(g+τ_j) + F(g+τ_i+τ_j)`
    which, integrated in `y`, **is** the `S`-difference of :eq:`pcell-from-smat`. Include:
    it is the only place the combination exists as a named object rather than an inline
    expression.
    ⚠ Note its **native** equation is a different one: its docstring already says
    *"See :eq:`cp-second-difference-operator` for the derivation"*, and that label is real
    and `tested` (`docs/theory/references/peierls_nystrom.rst:6071`). Declaring it against
    `s-integral` as well is legitimate (a symbol may implement several equations) but is
    an addition to an existing home, not a rehoming — say so in the directive if the
    project cares about primary ownership.
- **mechanism note (important, keep it with the declaration)**: the `S` array itself is
  never formed. The code takes the second difference **first** and the `y`-integral
  **outermost**; the equation states the `y`-integral first and the difference outermost.
  These are equal by linearity of the `y`-integral, so the declaration is honest — but a
  reader hunting for a literal `S = np.zeros(...)` will not find one, and that is *by
  design*, stated by the page.
- **confidence**: medium-high. High that these are the only candidate sites; medium on the
  framing, because a strict reader could call `s-integral` a `canonical-form` (a form
  exhibited to show structure that no production path takes). ⛔ I do **not** recommend
  that verdict: the page names a method for it in the imperative present tense, and the
  functional it defines is computed — only its evaluation order differs. Calling it
  `canonical-form` would suppress 23 claims that a real code path does carry.

---

## `self-double-integral`  (23 claims)  — page line 1054

- **verdict**: DECLARABLE
- **rationale comment on the page**: the nearest `.. (vv-status rationale)` block is at
  page line 1108 and belongs to the *downstream* label
  `self-collision-probability-slab`, not to this one; it reads verbatim:
  > `.. (vv-status rationale: definitional normalisation P_ii = r_ii/tau_i`
  > `   of the wired reduced form self-slab (a verifies-target of`
  > `   tests/cp/test_slab.py); the optical limits below are pinned by`
  > `   tests/cp/test_verification.py::TestOpticalLimits.)`
  For `self-double-integral` itself the authored pointer is in the derivation prose: the
  evaluation ends *"which is exactly the reduced form :eq:`self-slab` and,
  **character-for-character, the solver line** ``rcp[i, i] += sti * t[i] - (0.5 -
  _e3(tau_i))``"*, and the page continues *"For **cylindrical and spherical** geometries,
  the same structure holds … The term ``kernel_zero - kernel(tau_i)`` is F(0) − F(τ_i)"*
  (page lines ~1128-1140).
- **what the equation says**: the self-collision reduced probability is the double integral
  of the point-to-point kernel over the region's own optical thickness,
  `r_ii = ∫₀^{τ_i}∫₀^{τ_i} F₁(|s−t|) dt ds` — which collapses to `τ_i − (F(0) − F(τ_i))`.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.cp.solver.CPMesh._compute_slab_rcp` — `orpheus/cp/solver.py:280` —
    `rcp[i, i] += sti * t[i] - (0.5 - _e3(tau_i))`; the page calls this
    character-for-character the evaluated equation.
  - `orpheus.cp.solver.CPMesh._compute_radial_rcp` — `orpheus/cp/solver.py:354-357` —
    `self_same = 2.0 * chords[i, :] - (2.0 / sti) * (kernel_zero - kernel(tau_i))`, then
    `rcp[i, i] += 2.0 * sti * np.dot(y_wts, self_same)`; the curvilinear form of the same
    double integral, per chord.
  - `orpheus.derivations.continuous.flat_source_cp.geometry.build_cp_matrix` —
    `orpheus/derivations/continuous/flat_source_cp/geometry.py:325-329` — the derivations
    realization, covering all three geometries through `geom.kernel_F3`.
- **confidence**: high. The self term is the one contribution that cannot come from the
  gap formula, so it is a distinct, separately-written expression in every implementation —
  there is no risk that this equation has no site.

---

## `tau-m`  (1 claim)  — page line 1286

- **verdict**: DECLARABLE
- **rationale comment on the page**: none for `tau-m` itself. The *neighbouring* label
  `cumulative-optical-path-radial` (page line 1309) carries one, and it is worth quoting
  because it calibrates how this page's author thinks about this cluster:
  > `.. (vv-status rationale) definition: a cumulative-sum construction of the`
  > `.. boundary-position array from the per-region optical thicknesses. It is`
  > `.. computed identically in the solver and all three derivation scripts, so`
  > `.. it is a shared definitional primitive (not a structurally-independent`
  > `.. pinning target); the downstream radial rcp it feeds is verified by the`
  > `.. per-geometry CP eigenvalue / pinf-vs-derivation gates.`
  Two lines below the `tau-m`/`tau-p` pair the page also states the code form outright:
  *"The gaps: ``gap_d = bnd_pos[j] - bnd_pos[i+1]`` (same-side) and
  ``gap_c = bnd_pos[i] + bnd_pos[j]`` (through-centre)."*
- **what the equation says**: the same-side optical gap between region `i`'s outer boundary
  and region `j`'s inner boundary along the chord at height `y`, `τ_m(y) = Σ_{k=i+1}^{j} τ_k(y)`.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.cp.solver.CPMesh._compute_radial_rcp` — `orpheus/cp/solver.py:364,366,371` —
    `gap_d = bnd_pos[j, :] - bnd_pos[i + 1, :]` (and the mirrored `j < i` branch), clamped
    at `:371`. The telescoped form of the equation's sum.
  - `orpheus.derivations.continuous.flat_source_cp.geometry.build_cp_matrix` —
    `orpheus/derivations/continuous/flat_source_cp/geometry.py:336,338` — identical, with
    the clamp fused into the same expression.
- **deliberately EXCLUDED, with reason**: `CPMesh._compute_slab_rcp` computes a
  structurally identical `gap_d` at `orpheus/cp/solver.py:287,289` — but the **slab has its
  own labels** for this construction (`cumulative-optical-path-slab` p.1385, `dd-slab`
  p.1402, `dc-slab` p.1411), and `tau-m`/`tau-p` live in the chapter *"Optical Path
  Construction Along a Chord"* whose own opening sentence scopes it: *"the most
  geometrically intricate part of the **cylindrical and spherical** CP implementations"*.
  Declaring the slab site here would double-cover a construction the page already
  partitions by geometry. ⚠ `build_cp_matrix` is listed anyway because it is ONE body that
  serves all three geometries — there is no slab-only site in `derivations/` to exclude.
- **confidence**: high on the sites; medium on the slab exclusion, which is a page-scoping
  judgement rather than a fact about the code. What would change it: a ruling that
  geometry-partitioned labels should each claim every structurally-equivalent site.

---

## `tau-p`  (1 claim)  — page line 1297

- **verdict**: DECLARABLE
- **rationale comment on the page**: none for `tau-p` itself; same neighbouring
  `cumulative-optical-path-radial` rationale quoted under `tau-m` above, plus the same
  *"The gaps: … ``gap_c = bnd_pos[i] + bnd_pos[j]`` (through-centre)"* sentence.
- **what the equation says**: the through-centre optical path — the same-side gap plus
  twice the optical thickness of everything interior to region `i`,
  `τ_p(y) = τ_m(y) + 2 Σ_{k=1}^{i} τ_k(y)`.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.cp.solver.CPMesh._compute_radial_rcp` — `orpheus/cp/solver.py:379` —
    `gap_c = bnd_pos[i, :] + bnd_pos[j, :]`. This **is** the equation, telescoped:
    `bnd_pos[i] + bnd_pos[j] = Σ₁^i τ_k + Σ₁^j τ_k = τ_m + 2 Σ₁^i τ_k` exactly.
  - `orpheus.derivations.continuous.flat_source_cp.geometry.build_cp_matrix` —
    `orpheus/derivations/continuous/flat_source_cp/geometry.py:351` — identical line.
- **deliberately EXCLUDED, with reason**: `CPMesh._compute_slab_rcp:301` computes the same
  `gap_c` — excluded for the same page-scoping reason given under `tau-m` (the slab's
  through-centre path is the page's own `dc-slab`, p.1411).
- **confidence**: high on the sites; medium on the slab exclusion (same caveat as `tau-m` —
  the two should be decided together, never one each way).

---

## `neutron-balance`  (23 claims)  — page line 1703

- **verdict**: DECLARABLE
- **rationale comment on the page**: none as a `.. (vv-status rationale)` block within the
  section. The authored pointers are two, both explicit: (a) the "Matrix Form" paragraph —
  *":eq:`neutron-balance` in matrix form: `AΦ = (1/k)BΦ` … The analytical verification
  eigenvalue is `λ_max(A⁻¹B)`, **computed by
  :func:`~orpheus.derivations.common.eigenvalue.kinf_from_cp`**, which builds the full
  NG×NG matrices and uses ``numpy.linalg.eigvals``. **The solver does NOT form these
  matrices** — see :ref:`why-not-full-matrices`."*; (b) the "Power Iteration" numbered
  list, which maps each term of the balance onto a `CPSolver` method by name
  (`compute_fission_source`, `solve_fixed_source`, `compute_keff`).
- **what the equation says**: the multi-group CP balance — the collision rate in region `i`
  group `g` equals the transported total source (fission scaled by `1/k`, plus scattering,
  plus twice the (n,2n) emission) collected over all birth regions through `P^∞`.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.cp.solver.CPSolver._solve_fixed_source_jacobi` — `orpheus/cp/solver.py:557-584` —
    the production realization solved for `φ`: builds `Q = fission + Σ_sᵀφ + 2Σ_2ᵀφ`
    (`:563-568`), applies `P^T (V Q)` per group (`:572-574`), divides by `Σ_t V` (`:576-580`).
    Every term of the equation appears in that block.
  - `orpheus.cp.solver.CPSolver._solve_fixed_source_gs` — `orpheus/cp/solver.py:586-672` —
    the Gauss-Seidel arm; the same balance, swept group-by-group with the latest `φ`
    (`:617-641`). It is a **second, independently-written** realization of the same
    equation (a row-extraction of the same operators), not a wrapper — it must be declared
    or one of the two shipped solver modes has no declared implementer.
  - `orpheus.cp.solver.CPSolver.solve_fixed_source` — `orpheus/cp/solver.py:550-555` — the
    two-arm dispatcher on `self.solver_mode`, and the method the page names for step 2 with
    the equation `φ_ig^new = Σ_j P_ji,g V_j Q_jg / (Σ_t,ig V_i)` written out. ⚠ It carries
    no arithmetic; include it only if you want the page's own named symbol covered. I lean
    include (the page states the formula against this name), but it is the weakest of the
    six and dropping it costs nothing since both arms are declared.
  - `orpheus.cp.solver.CPSolver.compute_fission_source` — `orpheus/cp/solver.py:544-548` —
    the `χ_ig Σ_g' νΣ_f,ig' φ_ig' / k` term of the bracket, as a separate method (the page's
    step 1).
  - `orpheus.cp.solver.CPSolver._compute_balance_residual` — `orpheus/cp/solver.py:675-704` —
    ⭐ **the most direct implementer on the page**: it forms both sides of this equation
    (`collision = Σ_t V φ` at `:688`; `transported = P^T (V Q)` at `:700-702`) and returns
    their norm. It exists for no other purpose than to evaluate `neutron-balance`.
  - `orpheus.derivations.common.eigenvalue.kinf_from_cp` —
    `orpheus/derivations/common/eigenvalue.py:237-300` — the **independent derivations
    realization**: assembles the same balance as the dense `A`/`B` pencil and solves it
    with `numpy.linalg.eigvals`. This is the other half of the pair whose agreement to 1e-8
    is what makes the 23 claims on this label meaningful.
    ⚠ Its more precise homes are `matrix-A-def` (p.1730) and `matrix-B-def` (p.1736), which
    are outside this brief's 15 but on the same page — the page itself introduces them as
    ":eq:`neutron-balance` in matrix form". Declaring it here is an addition, not a
    rehoming; if you also declare A/B-def in a later pass, keep all three.
- **deliberately EXCLUDED, with reason**:
  - `orpheus.cp.solver.solve_cp` — `orpheus/cp/solver.py:881+` — the canonical entry point,
    but its body is mesh defaulting, `P_inf` assembly per group, and a `power_iteration`
    call. No term of the balance is written in it.
  - `orpheus.cp.solver.CPSolver.compute_keff` — implements `cp-keff-update` (p.1762), a
    different labelled equation on this page with its own (n,2n) net-removal subtlety.
- **confidence**: high.

---

## `p-transpose-flux-balance`  (1 claim)  — page line 2070

- **verdict**: DECLARABLE
- **rationale comment on the page**: none as a `.. (vv-status rationale)` block. The
  authored pointer sits immediately below the equation and is exact: *"In code:
  ``phi[:, g] = P_inf[:, :, g].T @ source`` (:meth:`CPSolver._solve_fixed_source_jacobi`)."*
  The section also carries the ERR-009 history (`P @ source` instead of `P.T @ source`:
  correct for homogeneous problems, 8 % wrong for the 1G 2-region slab), which is exactly
  the defect a declared `implements` edge makes traceable.
- **what the equation says**: because `P` is indexed `[birth, collision]`, the balance sums
  over the **first** index — a column sum — so the flux update contracts with `Pᵀ`, not `P`.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.cp.solver.CPSolver._solve_fixed_source_jacobi` — `orpheus/cp/solver.py:574` —
    `phi[:, g] = self.P_inf[:, :, g].T @ source`; the line the page quotes.
  - `orpheus.cp.solver.CPSolver._solve_fixed_source_gs` — `orpheus/cp/solver.py:633` —
    `transported_g = self.P_inf[:, :, g].T @ source_g`; the same contraction in the
    Gauss-Seidel arm. **The page does not name it** — it is exactly the kind of second site
    a page-driven declaration misses, and it is a live production path (`solver_mode =
    "gauss_seidel"`).
  - `orpheus.cp.solver.CPSolver._compute_balance_residual` — `orpheus/cp/solver.py:702` —
    `transported[:, g] = self.P_inf[:, :, g].T @ (V * Q[:, g])`; the third `.T` site, and
    the one whose disagreement with the other two would be silent (it would make the
    residual diagnostic lie rather than the flux wrong).
- **completeness evidence**: `grep -n "P_inf\[:, :, g\]" orpheus/cp/solver.py` returns
  **4** lines — `:574`, `:633`, `:702` (the three `.T` contractions above) and `:953`,
  which is not a contraction at all but the *assembly* assignment
  `P_inf[:, :, g] = cp_mesh.compute_pinf_group(xs.sig_t[:, g])` inside `solve_cp`. So:
  **3 of 4 hits are contraction sites, and all 3 are transposed.** No fourth contraction
  exists in `orpheus/`.
- **confidence**: high.

---

# Cross-cutting findings

## A. Measured graph state — these 15 have **no** guesses to stand down

`[M]` 2026-08-18, against `/Users/rodrigo/git/nuclear/ORPHEUS/.nexus/graph.db`
(13 206 `implements` edges total in the graph): **all 15 equations have ZERO incoming
`implements` edges — not one, inferred or declared.** Verified per label by walking
`g.in_edges(<eq node>)` and filtering `type == 'implements'`; e.g. `math:equation:p-inf`
has degree 29 made of `tests: 25`, `contains: 2`, `equation_ref: 2` and **no**
`implements`.

⟹ **The brief's stand-down hazard does not bite on this population.** The name-token
heuristic never fired for these labels (their spellings — `p-inf`, `tau-m`,
`self-double-integral` — share no token with any symbol name), so every declaration below
is pure gain and there is no guessed edge that a partial answer could destroy. That
*lowers the cost of an omission* but does not remove it: an omitted implementer stays
unlinked, and the 207 claims stay partly unadjudicable.

## B. Doc-drift found while deriving (all present-tense-false; not fixed — `docs/` is
out of my scope)

1. ⛔ **`p-inf` / `surface-to-surface` cite lines that no longer exist.** The page says
   *"implemented … independently in all three derivation scripts (e.g.,
   `orpheus/derivations/continuous/flat_source_cp/slab.py`:
   `P_inf_g[:,:,g] = P_cell + np.outer(P_out, P_in) / (1.0 - P_inout)`)"* and, under
   `surface-to-surface`, *"the same formula appears in all three derivation scripts (e.g.
   `…/slab.py`, line `P_in = sig_t_g * t_arr * P_out` … and `…/cylinder.py`, line
   `S_cell = 2.0 * np.pi * r_cell`)"*.
   `[M]` `grep -n "P_in = \|S_cell = \|P_inf_g\[:, *:, *g\]" …/flat_source_cp/{slab,cylinder,sphere}.py`
   returns **0 hits in all three files**. Phase B.2 collapsed the three scripts into the
   single `geometry.build_cp_matrix`; the three modules are now facades with no arithmetic.
   The page still describes the pre-B.2 layout, and it says "all three … independently",
   which is the substantive half of the falsehood — there is **one** derivations
   implementation, not three, so the independence the sentence advertises is gone.
2. ⛔ **`reciprocity-lower-triangle`'s "This halves the computation cost."** No
   lower-triangle fill-in exists; all three implementations run the full N² loop
   (`orpheus/cp/solver.py:273,282` and `:347,359`;
   `…/flat_source_cp/geometry.py:319,331`). Detail under that equation above.
3. ⚠ **`ki3-def`'s naming caveat is already documented on the page** (the A&S-convention
   note at page line ~712 explaining that the page's `Ki₃` is canonical `Ki₂^{A&S}`), so it
   is *not* drift — but any declarer touching the kernel labels should read it first, since
   the label `ki3-def` is deliberately mis-named for cross-consistency with existing
   `verifies("ki3-def")` decorators.

## C. Refuted candidates — symbols I considered and rejected, with the structural reason

| candidate | why NOT an implementer of any of the 15 |
|---|---|
| `…flat_source_cp.slab._slab_cp_matrix` | Post-B.2 thin facade: pre-binds `SLAB` and returns `build_cp_matrix(...)`. Zero arithmetic. (Same for `cylinder._cylinder_cp_matrix`, `sphere._sphere_cp_matrix`.) Declaring them would make three names claim one body's equation. |
| `orpheus.cp.solver._cp_bc_white` | `BC_REGISTRY` factory; body is `return cp_mesh._apply_white_bc`. Restates `p-inf` in its docstring, computes nothing. |
| `orpheus.cp.solver._cp_bc_vacuum` | Implements `P_inf = P_cell` — the *absence* of the white-BC closure, i.e. the complementary BC, not `p-inf`. |
| `orpheus.cp.solver.solve_cp` | Canonical entry point, but its body is mesh defaulting + per-group `compute_pinf_group` assembly + `power_iteration`. No equation term is written in it. |
| `orpheus.cp.solver.CPSolver.compute_keff` | Implements `cp-keff-update` (page line 1762) — a different labelled equation with its own (n,2n) net-removal derivation (ERR-015). |
| `orpheus.derivations.common.kernels.e_n` / `e_n_mp` / `e_n_derivative` | The **general** `E_n`; the general definition has its own labels on another page (`en-definition`, `en-kernel-special-values`, `en-kernel-derivative`, `en-kernel-integral` — `docs/theory/verification/reference_solutions.rst:232-255`). Their consumers are Peierls-Nyström, not CP. |
| `…flat_source_cp.geometry._ki3_mp`, `_ki3_scaled_cheb`, `_exp_kernel` | The cylinder / sphere kernels — they implement `ki3-def` and the sphere's `e^{−τ}`, not `e3-def`. |
| `…peierls_nystrom.geometry.CurvilinearGeometry.reciprocity_factor` | Named "reciprocity" but is an outer/inner **surface-area ratio** for the F.4 Nyström closure — a different relation on a different page. This is the only `reciprocit*` hit in `orpheus/` that could be mistaken for a CP guard; it is not one. |
| `orpheus.derivations.common.quadrature_recipes.chord_quadrature` | Supplies the `y`-quadrature *rule* for the outer integral. It is measure machinery, not the CP equation; if the project wants it declared it belongs to the chord/quadrature labels, not to `s-integral`. |
| `CPMesh._compute_slab_rcp` (for `tau-m` / `tau-p` only) | Computes a structurally identical `gap_d`/`gap_c`, but the page partitions this construction by geometry and gives the slab its **own** labels (`cumulative-optical-path-slab`, `dd-slab`, `dc-slab`). See the note under `tau-m`. |

## D. What I could NOT settle

Nothing was left unresolved as a **verdict** — all 15 are DECLARABLE and none is a
`NOTHING:<kind>`. Three items are judgement calls I have flagged rather than decided,
because each turns on a project convention rather than on a fact about the code, and each
is one sentence of ruling away:

1. **`optical-path` — does `chord_half_lengths` count?** It supplies the `ds` measure but
   never multiplies by `Σ_t`. Include ⟹ the curvilinear `ds` has a site; exclude ⟹ the
   equation is covered by the three `τ = Σ_t · ℓ` sites alone.
2. **`surface-to-region` — does the general-geometry `compute_areas_1d` count?** It is
   where *production* gets `S`, and the page names it; but it also serves SN/MoC/diffusion,
   so declaring it widens a CP-page equation past CP. Note the asymmetry if you exclude it:
   the derivations' `S` (`FlatSourceCPGeometry.surface_area`) would be declared while
   production's would not.
3. **`tau-m` / `tau-p` — slab in or out?** Decided OUT above on the page's own geometry
   partition. The two must be decided together, never one each way.

One equation deserves a second reader's eye on the *framing* rather than the site list:
**`s-integral`**, whose `S` array is provably never materialized (the page says so
outright). I ruled it DECLARABLE because the functional it defines *is* computed — only
the order of the difference and the `y`-integral is swapped, which is an equality by
linearity — and because calling it a `canonical-form` would suppress 23 live claims. If
the project's convention for `implements` is "the code forms the object the equation
names", then `s-integral` (and only `s-integral`, of these 15) would flip to
`NOTHING:canonical-form`, and `pcell-from-smat` — outside this brief but adjacent — would
need the same ruling in the same pass.

---

# Appendix — ready-to-land declaration list

Directive form taken from the shipped precedent
`docs/theory/methods/sn/loss_representation.rst:219` (read, not recalled): one
`.. implements:: <label>` block per implementer, `:by:` carrying the dotted path, body
prose saying *why this IS the implementation*. The `:by:` values below are exactly the
21 distinct symbols verified above; **every one resolves to exactly one node of a legal
type** (`function` / `method` / `class`), re-checked at close.

**43 core declarations + 4 flagged-optional = 47.**

| label | `:by:` (core) | flagged-optional |
|---|---|---|
| `optical-path` | `orpheus.cp.solver.CPMesh._compute_slab_rcp`<br>`orpheus.cp.solver.CPMesh._compute_radial_rcp`<br>`orpheus.derivations.continuous.flat_source_cp.geometry.build_cp_matrix` | `orpheus.derivations.common.kernels.chord_half_lengths` |
| `complementarity` | `orpheus.cp.solver.CPMesh._apply_white_bc`<br>`orpheus.derivations.continuous.flat_source_cp.geometry.build_cp_matrix` | — |
| `reciprocity` | `orpheus.cp.solver.CPMesh._compute_slab_rcp`<br>`orpheus.cp.solver.CPMesh._compute_radial_rcp`<br>`orpheus.cp.solver.CPMesh._normalize_rcp`<br>`orpheus.derivations.continuous.flat_source_cp.geometry.build_cp_matrix` | — |
| `reciprocity-lower-triangle` | `orpheus.cp.solver.CPMesh._normalize_rcp`<br>`orpheus.derivations.continuous.flat_source_cp.geometry.build_cp_matrix` | — |
| `surface-to-region` | `orpheus.cp.solver.CPMesh._apply_white_bc`<br>`orpheus.derivations.continuous.flat_source_cp.geometry.build_cp_matrix`<br>`orpheus.derivations.continuous.flat_source_cp.geometry.FlatSourceCPGeometry.surface_area` | `orpheus.geometry.coord.compute_areas_1d` |
| `surface-to-surface` | `orpheus.cp.solver.CPMesh._apply_white_bc`<br>`orpheus.derivations.continuous.flat_source_cp.geometry.build_cp_matrix` | — |
| `p-inf` | `orpheus.cp.solver.CPMesh._apply_white_bc`<br>`orpheus.derivations.continuous.flat_source_cp.geometry.build_cp_matrix`<br>`orpheus.cp.solver.CPMesh.compute_pinf_group` | — |
| `cp-infinite-lattice-sum` | `orpheus.cp.solver.CPMesh._apply_white_bc`<br>`orpheus.derivations.continuous.flat_source_cp.geometry.build_cp_matrix` | — |
| `e3-def` | `orpheus.cp.solver._e3`<br>`orpheus.derivations.common.kernels.e3`<br>`orpheus.derivations.common.kernels.e3_vec`<br>`orpheus.derivations.continuous.flat_source_cp.geometry.FlatSourceCPGeometry.kernel_F3` | `…FlatSourceCPGeometry.kernel_F3_at_zero` |
| `s-integral` | `orpheus.cp.solver.CPMesh._compute_radial_rcp`<br>`orpheus.derivations.continuous.flat_source_cp.geometry.build_cp_matrix`<br>`orpheus.derivations.continuous.flat_source_cp.geometry._second_difference` | — |
| `self-double-integral` | `orpheus.cp.solver.CPMesh._compute_slab_rcp`<br>`orpheus.cp.solver.CPMesh._compute_radial_rcp`<br>`orpheus.derivations.continuous.flat_source_cp.geometry.build_cp_matrix` | — |
| `tau-m` | `orpheus.cp.solver.CPMesh._compute_radial_rcp`<br>`orpheus.derivations.continuous.flat_source_cp.geometry.build_cp_matrix` | — |
| `tau-p` | `orpheus.cp.solver.CPMesh._compute_radial_rcp`<br>`orpheus.derivations.continuous.flat_source_cp.geometry.build_cp_matrix` | — |
| `neutron-balance` | `orpheus.cp.solver.CPSolver._solve_fixed_source_jacobi`<br>`orpheus.cp.solver.CPSolver._solve_fixed_source_gs`<br>`orpheus.cp.solver.CPSolver.compute_fission_source`<br>`orpheus.cp.solver.CPSolver._compute_balance_residual`<br>`orpheus.derivations.common.eigenvalue.kinf_from_cp` | `orpheus.cp.solver.CPSolver.solve_fixed_source` |
| `p-transpose-flux-balance` | `orpheus.cp.solver.CPSolver._solve_fixed_source_jacobi`<br>`orpheus.cp.solver.CPSolver._solve_fixed_source_gs`<br>`orpheus.cp.solver.CPSolver._compute_balance_residual` | — |

## Drift check at close (L-012)

- `git status --short -- orpheus/ docs/theory/methods/collision_probability.rst` → **empty**
  at close; nothing I audited carries uncommitted state.
- HEAD moved `a1c90aac → 58e46c6f` during this dispatch (≈20 commits, almost all
  `.claude/` memory + plans). `git log a1c90aac..HEAD -- docs/theory/methods/collision_probability.rst
  orpheus/cp/ orpheus/derivations/` names two commits touching the CP page, but the page is
  **byte-stable across my read window**: re-running the label grep at close returns the
  identical 15 line numbers and the identical 6 `vv-status rationale` positions, and
  `wc -l` is 2296 both times.
- ⭐ **Checked for conflict with the concurrent NONE campaign**: `bb075c93` ("eleven
  equations declare that NOTHING implements them") landed `.. no-implementation::`
  directives for `operator-solve`, `solve-does-not-distribute`,
  `apply-solve-denominator-inequality`, `apply-solve-neumann-series`,
  `apply-solve-neumann-expansion`, `apply-solve-parallel-identity`, `carrier-grid-cell`,
  `streaming-as-tensor-product-sum`, `loss-rep-leaf-sum`, `loss-rep-removal-sigma`,
  `loss-rep-facewise-separable` — **all on `foundations/operator_algebra.rst` and
  `methods/sn/loss_representation.rst`. None of my 15 is among them**, so nothing here
  contradicts a landed declaration.
