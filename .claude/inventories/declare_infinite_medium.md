# `implements` declarations for `docs/theory/foundations/infinite_medium.rst`

Derived from the tree at HEAD `a1c90aac` (2026-08-18), branch `main`.
Scope: the 12 equations named in the brief (159 `tests` claims between them).

## Two facts that govern how to read every verdict below

1. **⚠ The page carries NO `.. (vv-status rationale)` comment on any of these
   12 equations.** The comment form the brief describes DOES exist on this page
   (`.. vv-status: <label> documented` + free-comment prose), but only on the
   ten OTHER labels: `boltzmann`, `group-flux`, `group-xs`, `sigs-convention`,
   `sigs-in-scatter-transpose`, `one-over-E`, `flux-per-lethargy-plateau`,
   `maxwellian`, `resolvent-similarity`, `pin-cell-volume-fractions`
   (`grep -n "vv-status" docs/theory/foundations/infinite_medium.rst`).
   So the "read the authored rationale first" step returns *nothing* here and
   every verdict below is derived from the equation + prose + tree.
   **The authored knowledge that DOES exist, and that nothing currently reads,
   is somewhere else** — see fact 2.

2. ⭐ **`[M]` NONE of the 12 has any `implements` edge today — not one, not even
   an inferred one.** Measured against `.nexus/graph.db`: in-edge types are
   only `contains` / `equation_ref` / `tests`.

   | label | contains | equation_ref | tests | implements |
   |---|---|---|---|---|
   | `mg-balance` | 2 | 5 | 70 | **0** |
   | `inf-hom-balance` | 2 | 1 | 12 | **0** |
   | `two-group-A` | 2 | 1 | 12 | **0** |
   | `two-group-F` | 2 | 0 | 12 | **0** |
   | `two-group-Ainv` | 2 | 0 | 12 | **0** |
   | `two-group-M` | 2 | 0 | 12 | **0** |
   | `two-group-charpoly` | 2 | 1 | 8 | **0** |
   | `two-group-roots` | 2 | 1 | 8 | **0** |
   | `keff-update` | 2 | 1 | 8 | **0** |
   | `number-density` | 2 | 0 | 3 | **0** |
   | `normalisation` | 2 | 0 | 1 | **0** |
   | `resolvent-object-gate` | 2 | 0 | 1 | **0** |

   ⟹ **the brief's "a partial answer is worse than none" hazard does not bind
   on this population**: there is no name-token guess here to stand down. A
   declaration can only add information. (It still binds *forward* — a wrong
   declaration is a wrong edge — so each row below carries its own confidence.)

## ⭐ The authored knowledge that answers most of this, and that nothing reads

`orpheus/derivations/continuous/analytical/homogeneous.py` already carries a
hand-written `equation_labels=(...)` tuple on each `derive_*` case — the page
author's own statement of which equations that reference realizes. It is the
`.. (vv-status rationale)` of this page, living in the code instead of the
page:

- `derive_2g` / `derive_2g_continuous` → `matrix-eigenvalue`, `removal-matrix`,
  `fission-matrix`, `mg-balance`, and **all six** `two-group-*`
  (`homogeneous.py:269-280`, `:594-598`).
- `derive_1g` / `derive_1g_continuous` → `one-group-kinf`, `inf-hom-balance`
  (`:201`, `:550`).
- `derive_4g` / `derive_4g_continuous` / `derive_2g_n2n` /
  `derive_2g_n2n_continuous` → `matrix-eigenvalue`, `removal-matrix`,
  `fission-matrix`, `mg-balance`.

⚠ **But the tuple is a CLAIM, not a measurement, and two of its six
`two-group-*` entries are not honoured by the function's own body** — see
`two-group-charpoly` / `two-group-roots` below. It is the best available lead
and it is not self-verifying.

---

## `mg-balance`  (70 claims)

- **verdict**: DECLARABLE
- **rationale comment on the page**: none (see fact 1). The nearest authored
  statement is `equation_labels` on six `derive_*` cases in
  `orpheus/derivations/continuous/analytical/homogeneous.py` (`:273`, `:346`,
  `:418`, `:596`, `:642`, `:694`), plus the sibling label
  `sigs-in-scatter-transpose`'s own rationale, which says its "terminal use is
  the removal matrix" — i.e. the in-scatter term of THIS equation
  (`infinite_medium.rst:361-366`).
- **what the equation says**: per energy group `g`, the collision rate
  `Σ_t,g φ_g` equals the in-scatter from every group `Σ_g' Σ_s(g'→g) φ_g'`
  plus the fission source `(χ_g/k)·Σ_g' νΣ_f,g' φ_g'`. It is the group-
  discretised form of `inf-hom-balance`, one step before the matrix collection
  (`matrix-eigenvalue` / `removal-matrix` / `fission-matrix`).

### ⛔ Read this before declaring — the 70 claims are NOT one population

`[M]` the 70 `tests` edges split across **five solver families**, not one:

| claiming family | claims | what it actually runs |
|---|---|---|
| SN (`tests/sn/*`) | 26 | `orpheus/transport/operators/scattering.py` moment operators |
| CP (`tests/cp/*`) | 26 | `orpheus/cp/solver.py` — its own `SigS[0].T @ φ` transcription |
| homogeneous (`tests/homogeneous/*`) | 12 | `orpheus/homogeneous/solver.py` + `orpheus/transport/operators/isotropic_scattering.py` |
| MC (`tests/mc/*`) | 5 | `orpheus/mc/` group-transfer sampling |
| MoC (`tests/moc/*`) | 1 | `orpheus/moc/core.py:183` — its own `sig_s0[i].T @ phi` |

⚠ `[M]` **there is no shared multigroup-balance implementer across the five.**
`orpheus/transport/operators/scattering.py` is imported by `orpheus/sn/` only
(`solver.py`, `coupled_system.py`); CP and MoC each hand-roll the in-scatter
from `mix.SigS[0]` (`cp/solver.py:511`, `moc/core.py:93`). That is a genuine
4-way twin path of one equation (Cardinal Rule 2 observation, out of scope
here) — and operationally it means **a Tier-1-only declaration will leave 58
of the 70 claims refuted on execution evidence**, because those tests never
enter the homogeneous chain. Declare Tier 2 as well, or accept that outcome
knowingly.

- **implementers — TIER 1** (the page's own 0-D form + the model-shared terms;
  all verified to resolve):
  - `orpheus.homogeneous.solver.solve_homogeneous_infinite` — `orpheus/homogeneous/solver.py:150` — poses and solves the collected balance `Aφ = (1/k)Fφ` for the infinite medium; the page's `direct-eigensolve` section is a walkthrough of this function.
  - `orpheus.homogeneous.solver._assemble_loss_operator` — `orpheus/homogeneous/solver.py:115` — builds the balance's LHS-minus-in-scatter, `A = C − K_iso = diag(Σ_t) − Σ_s0ᵀ − 2Σ_2ᵀ`.
  - `orpheus.transport.operators.isotropic_scattering.IsotropicScattering.apply` — `orpheus/transport/operators/isotropic_scattering.py:263` — computes the in-scatter term verbatim; its own docstring is `(Σ_s0ᵀφ)_g = Σ_g' Σ_s0(g'→g) φ_g'`.
  - `orpheus.transport.mesh.material_xs_field.MaterialXSField.apply_p0_in_scatter` — `orpheus/transport/mesh/material_xs_field.py:815` — the kernel the operator above delegates to (the single per-material dispatch home).
  - `orpheus.transport.operators.fission.FissionOperator.apply` — `orpheus/transport/operators/fission.py:669` (overloads; class at `:147`) — computes the fission-source term `χ_g·(νΣ_f·φ)`, the group contraction + emission broadcast.
  - `orpheus.derivations.common.eigenvalue._infinite_medium_matrices` — `orpheus/derivations/common/eigenvalue.py:39` — the ORACLE's assembly of the same balance as `(A, F)`; the single site shared by the forward and adjoint references.
  - `orpheus.derivations.common.eigenvalue.kinf_and_spectrum_homogeneous` — `orpheus/derivations/common/eigenvalue.py:98` — solves that balance for `(k, φ)`; this is what every `derive_*` case calls, i.e. the function behind the `equation_labels` claim.
- **implementers — TIER 2** (the same balance in the families that supply 58
  of the 70 claims; each verified to resolve):
  - `orpheus.cp.solver.CPSolver._compute_balance_residual` — `orpheus/cp/solver.py:675` — ⭐ the closest thing in the tree to a literal transcription of this equation: docstring "Compute the neutron balance residual ‖Aφ − B1φ − (1/k)B2φ‖₂", body = collision rate `Σ_t·V·φ` minus (fission + `Σ_sᵀφ` + `2Σ_2ᵀφ`) transported. This is the `balance_residual`-shaped member the brief predicted — it is in CP, not in `Mixture` (see the ⚠ under `Mixture.balance_residual` below).
  - `orpheus.moc.core.MOCSolver.solve_fixed_source` — `orpheus/moc/core.py:121` — assembles `Q[i,:] = (fission + Σ_s0ᵀφ + 2Σ_2ᵀφ)/4π` at `:180-185`, the MoC realisation of the RHS.
  - `orpheus.transport.operators.scattering.LegendreMomentScattering` — `orpheus/transport/operators/scattering.py:115` — the SN in-scatter, ℓ-resolved (the P0 slot is this equation's `Σ_s(g'→g)` sum).
  - `orpheus.transport.operators.scattering.N2NMomentOperator` — `orpheus/transport/operators/scattering.py:300` — the SN `2Σ_2ᵀ` transfer.
  - `orpheus.transport.operators.scattering.ScatteringOperator` — `orpheus/transport/operators/scattering.py:356` — the composed SN scattering source.
- **UNSURE / not enumerated**: the MC family's group-transfer sampling
  (5 claims, `tests/mc/test_monte_carlo.py`). I did not open `orpheus/mc/` —
  a Monte-Carlo realisation of a *balance* is a sampled estimator rather than
  an assembled operator, and deciding whether that counts as "implements" is a
  judgement I am flagging rather than making.
- **⚠ NOT an implementer, despite the name**: `Mixture.balance_residual` /
  `Mixture.assert_balanced` (`orpheus/data/macro_xs/mixture.py:194`, `:223`).
  These compute the **total-cross-section definitional identity**
  `Σ_t = Σ_c + Σ_L + Σ_f + rowsum(Σ_s0) + rowsum(Σ_2n)` — a statement about the
  XS data, with no flux and no `k` in it. The brief's hint ("check for
  `balance_residual`-shaped members") points here by name and is wrong by
  meaning; the right hit is `CPSolver._compute_balance_residual`.
- **confidence**: **high** for Tier 1 (each symbol computes a named term of the
  printed equation, and the derivations tuple corroborates). **medium** for
  Tier 2 — the equation is the same balance but each family adds its own
  transport operator, so a purist could argue those implement the *method*
  page's balance rather than this page's streaming-free one. What would change
  it: a ruling on whether `mg-balance` is (a) the 0-D infinite-medium balance
  of THIS page, in which case Tier 2 is wrong and 58 claims are mis-tagged at
  the marker, or (b) the generic multigroup balance every solver discretises,
  in which case Tier 2 is required. The five `equation_ref` edges — from
  `theory/methods/{collision_probability,method_of_characteristics,monte_carlo,sn/index}`
  — are the tree's own evidence for reading (b).

## `inf-hom-balance`  (12 claims)

- **verdict**: DECLARABLE
- **rationale comment on the page**: none. The authored lead is
  `equation_labels=("one-group-kinf", "inf-hom-balance")` on **`derive_1g`**
  (`orpheus/derivations/continuous/analytical/homogeneous.py:201`) and
  `derive_1g_continuous` (`:550`) — i.e. the page author binds this label to
  the ONE-GROUP reduction, not to the G-group case (which carries
  `mg-balance` instead).
- **what the equation says**: the continuous-energy, angle-integrated,
  streaming-free balance — `Σ_t φ(E) = ∫Σ_s0(E'→E)φ(E')dE' + (χ(E)/k)∫νΣ_f φ dE'`.
  It is `boltzmann` after the three infinite-medium simplifications
  (`∇ψ=0`, position-independent XS, isotropy).
- **implementers** (all verified to resolve):
  - `orpheus.homogeneous.solver._assemble_loss_operator` — `orpheus/homogeneous/solver.py:115` — **the site where the streaming term is dropped**, which is what distinguishes this equation from `boltzmann`: its docstring states "Streaming `L` is identically zero in an infinite medium and dropped". The LHS-minus-in-scatter of exactly this balance.
  - `orpheus.homogeneous.solver.solve_homogeneous_infinite` — `orpheus/homogeneous/solver.py:150` — solves this balance; builds it over a *meshless single cell*, the code-level encoding of "infinite geometry".
  - `orpheus.derivations.common.eigenvalue._infinite_medium_matrices` — `orpheus/derivations/common/eigenvalue.py:39` — the reference assembly of the same streaming-free `(A, F)` pair.
  - `orpheus.derivations.common.eigenvalue.kinf_homogeneous` — `orpheus/derivations/common/eigenvalue.py:67` — the scalar-`k` face of the same solve.
  - `orpheus.derivations.common.eigenvalue.kinf_and_spectrum_homogeneous` — `orpheus/derivations/common/eigenvalue.py:98` — the `(k, φ)` face; the function every `derive_*` case calls.
  - `orpheus.derivations.continuous.analytical.homogeneous.derive_1g` — `orpheus/derivations/continuous/analytical/homogeneous.py:144` — the author-declared realisation: its generated LaTeX writes this balance out as `Σ_a φ = (1/k) νΣ_f φ` (the 1-group reduction, where the two integrals become scalars).
  - `orpheus.derivations.continuous.analytical.homogeneous.derive_1g_continuous` — `orpheus/derivations/continuous/analytical/homogeneous.py:503` — same, in the `ContinuousReferenceSolution` tier.
- **candidate I am NOT confident enough to declare**:
  `orpheus.transport.mesh.material_mesh.MaterialMesh.from_materials`
  (`orpheus/transport/mesh/material_mesh.py`) — the meshless single-cell
  carrier is the object that *encodes* "infinite geometry, no streaming", but
  it is a phase-space constructor with no balance arithmetic in it. Declaring
  it would be declaring the premise, not the equation. Flagging, not claiming.
- **confidence**: **high**. The distinguishing content of this equation
  relative to `mg-balance` is "streaming dropped, angle integrated out", and
  `_assemble_loss_operator` is the one function in the tree that makes that
  decision explicitly.
  What would change it: a ruling that "continuous energy" is the
  distinguishing content instead — ORPHEUS has no continuous-energy solver, so
  under that reading the verdict flips to `NOTHING:canonical-form` (the form
  exhibited to derive the multigroup one). I reject that reading because the
  page author bound the label to `derive_1g`, a *discrete* case.

---

## The `two-group-*` family — one shared preamble

All six carry the same authored lead and the same population, so the shared
facts are stated once here.

- **rationale comment on the page**: none for any of the six. The authored
  lead is `equation_labels` on `derive_2g`
  (`orpheus/derivations/continuous/analytical/homogeneous.py:269-280`) and
  `derive_2g_continuous` (`:594-598`), each of which lists **all six**.
- **claims**: 12 / 12 / 12 / 12 / 8 / 8 — all from `tests/homogeneous/`
  (`test_homogeneous.py` file-level `pytestmark`, `test_continuous_reference.py`
  module-level list). Nothing outside `tests/homogeneous/` claims any of them.
- ⭐ **The find that answers the back half of the family**:
  `orpheus/derivations/continuous/fn_method/origins/k_inf_derivations.py:516`,
  `derive_kinf_mg_matrix_form` — a **fully symbolic G=2** derivation that, in
  one function, builds `A = Σ_t − Σ_s` (2×2), takes `A_inv = A.inv()` (the only
  explicit symbolic 2×2 loss-matrix inverse in the tree), forms
  `M = A_inv·χ·νΣ_f`, takes `M.eigenvals()`, and proves `M` is rank-1 with its
  single non-zero eigenvalue equal to `tr M` — which is *verbatim* the page's
  own paragraph ("the two rows of `M` are proportional — `M` is rank 1 … the
  characteristic polynomial factors as `λ(λ − tr M) = 0` … the second
  eigenvalue is `λ_ = 0`").
  ⚠ Three caveats before declaring it: (i) it uses **Sood's** index convention
  (`Σ_s` already stored TO-row/FROM-column, groups numbered fast = 2), so its
  `A` is the page's `A` under a relabelling, not literally; (ii) it keeps `χ`
  and the upscatter entry **general**, so the page's forms are its
  specialisation at `χ = [1,0]`, `Σ_s(2→1) = 0`; (iii) `[M]` its only consumer
  is `tests/derivations/test_fn_la13511_kinf.py:146`, which carries
  `@pytest.mark.foundation` and **no `verifies` marker** — so none of the 12/8
  claiming tests execute it. Declaring it adds a true implementer that the
  claims cannot reach.
- ⛔ **A doc-drift finding that bears directly on `charpoly` / `roots`**:
  `derive_2g`'s own docstring says *"2-group infinite medium eigenvalue **via
  characteristic polynomial**"* and `derive_2g_continuous`'s `Provenance`
  says *"Solved via characteristic polynomial (closed form for 2x2)"*. `[M]`
  **both are present-tense-false**: the body of each calls
  `kinf_and_spectrum_homogeneous`, which forms `M = np.linalg.solve(A, F)` and
  calls `np.linalg.eig(M)` (`derivations/common/eigenvalue.py:131-133`). No
  characteristic polynomial is ever formed. The same false mechanism is
  restated in `tests/homogeneous/test_homogeneous.py:25-28` as the JUSTIFICATION
  for the `two-group-charpoly` / `two-group-roots` markers ("the analytical
  k_inf is derived symbolically via exactly those equations"). **The markers'
  stated warrant is refuted by the reference's own body** — worth fixing
  independently of any declaration.

## `two-group-A`  (12 claims)

- **verdict**: DECLARABLE
- **what the equation says**: the 2×2 loss matrix for downscatter-only,
  lower-triangular: `diag(Σ_t1 − Σ_s(1→1), Σ_t2 − Σ_s(2→2))` with `−Σ_s(1→2)`
  in the `[1,0]` slot (the transpose is why the downscatter entry sits *below*
  the diagonal).
- **implementers** (all verified to resolve):
  - `orpheus.derivations.continuous.analytical.homogeneous.derive_2g` — `orpheus/derivations/continuous/analytical/homogeneous.py:209` — builds `A_sym` as a literal 2×2 `sp.Matrix` with exactly the printed entries and slot placement (`:228-231`), and renders it into the case's LaTeX. The most specific implementer there is.
  - `orpheus.derivations.common.eigenvalue._infinite_medium_matrices` — `orpheus/derivations/common/eigenvalue.py:39` — `A = diag(σ_t) − (σ_s + 2σ_2)ᵀ`; at G=2 with `σ_2 = 0` this IS the printed matrix, and it is the `A` whose eigenvalue `derive_2g` actually reports.
  - `orpheus.homogeneous.solver._assemble_loss_operator` — `orpheus/homogeneous/solver.py:115` — the production assembly of the same matrix as an operator (`C − K_iso`).
- **secondary** (see preamble caveats): `orpheus.derivations.continuous.fn_method.origins.k_inf_derivations.derive_kinf_mg_matrix_form` — `orpheus/derivations/continuous/fn_method/origins/k_inf_derivations.py:516` — the symbolic 2×2 `A = Σ_t − Σ_s`.
- **confidence**: **high** for the three primary. What would change it: a ruling
  that a general-G assembler must not carry a G=2-specific label (in which case
  only `derive_2g` survives).

## `two-group-F`  (12 claims)

- **verdict**: DECLARABLE
- **what the equation says**: the 2×2 production matrix for `χ = [1,0]` —
  `[[ν₁Σ_f1, ν₂Σ_f2], [0, 0]]`, i.e. the rank-1 dyad `χ ⊗ νΣ_f` with all
  emission into group 1.
- **implementers** (all verified to resolve):
  - `orpheus.derivations.continuous.analytical.homogeneous.derive_2g` — `orpheus/derivations/continuous/analytical/homogeneous.py:209` — builds `F_sym` as the literal 2×2 with the zero second row (`:232-235`).
  - `orpheus.derivations.common.eigenvalue._infinite_medium_matrices` — `orpheus/derivations/common/eigenvalue.py:39` — `F = np.outer(chi, nu_sig_f)`; the G=2, `χ=[1,0]` instance is the printed matrix.
  - `orpheus.transport.operators.fission.FissionOperator` — `orpheus/transport/operators/fission.py:147` — the production dyad `F = χ ⊗ νΣ_f`, materialised densely by `as_matrix` on the homogeneous path. (Already carries an INFERRED edge to the general `fission-matrix`; this is its G=2 face.)
- **confidence**: **high**.

## `two-group-Ainv`  (12 claims)

- **verdict**: DECLARABLE — but weakly, and read the caveat
- **what the equation says**: the closed-form inverse of the lower-triangular
  `A`, entrywise: `[[1/Σ_r1, 0], [Σ_s(1→2)/(Σ_r1 Σ_r2), 1/Σ_r2]]`, with
  `Σ_r,g = Σ_t,g − Σ_s(g→g)` the removal cross section.
- ⚠ **Nothing in ORPHEUS forms `[A⁻¹]` on the production path, by design.**
  The theory page says so in as many words: "the loss matrix is still **solved
  out** of the production, never inverted into an explicit `[A⁻¹]`"
  (`infinite_medium.rst:1078-1080`). `kinf_and_spectrum_homogeneous` uses
  `np.linalg.solve(A, F)`; `MatrixInverseOperator` holds an LU factorisation
  and back-solves. `[M]` there is also **no `Σ_r` symbol anywhere** —
  `grep -rn "removal_xs|Sigma_r|removal cross" orpheus/ --include='*.py'`
  returns only `orpheus/fuel/solver.py` radial *stress* `sig_r`. The nearest
  the tree gets is `Mixture.absorption_xs + Mixture.out_scattering_xs`, which
  is not a named quantity.
- **implementers** (verified to resolve):
  - `orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator` — `orpheus/numerics/matrix_inverse_operator.py:95` — **is** `A⁻¹` for exactly this loss matrix, as an object: eager materialisation of `A` + one LU factorisation, `apply` = the backsolve. It is on the production path and every one of the 12 claiming tests executes it (there is even a dedicated Mode-11 gate asserting its `apply` fires, `test_homogeneous.py:294`).
  - `orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator.as_matrix` — `orpheus/numerics/matrix_inverse_operator.py:251` — ⭐ **the one method in the tree that forms `[A⁻¹]` explicitly**: `lu_solve(lu, I)`, "Materialize `[A⁻¹]` — one batched backsolve on the identity". This is the entrywise-inverse implementer the equation asks for. ⚠ `[M]` it is NOT reached by the homogeneous solve: `solve_homogeneous_infinite` calls `as_matrix` on the *product* `K`, which walks basis columns through `MatrixInverseOperator.apply` (the LU backsolve), never through this override — which is exactly what `test_matrix_inverse_operator_apply_is_on_the_homogeneous_call_path` (`tests/homogeneous/test_homogeneous.py:294`) pins.
  - `orpheus.derivations.continuous.fn_method.origins.k_inf_derivations.derive_kinf_mg_matrix_form` — `orpheus/derivations/continuous/fn_method/origins/k_inf_derivations.py:516` (line `:598`, `A_inv = A.inv()`) — the only place in the tree that forms an explicit symbolic 2×2 loss-matrix inverse. See preamble caveats (i)–(iii).
- **confidence**: **medium**. The declaration is honest for "the inverse of the
  two-group loss matrix"; it is NOT an implementation of the printed *closed
  form's entries*, which no code writes.
  What would change it: a ruling on whether an operator that *represents*
  `A⁻¹` without materialising it implements an equation that *exhibits* `A⁻¹`
  entrywise. If the answer is no, the verdict becomes
  `NOTHING:canonical-form` (a form exhibited to show the triangular structure,
  which the production path deliberately does not take) with
  `derive_kinf_mg_matrix_form` as the sole reference-tier exception.

## `two-group-M`  (12 claims)

- **verdict**: DECLARABLE
- **what the equation says**: the eigenvalue matrix `M = A⁻¹F`, written out
  entrywise for the downscatter-only, `χ=[1,0]` case; its two rows are
  proportional (rank 1) because `F` is a dyad.
- **implementers** (all verified to resolve):
  - `orpheus.derivations.common.eigenvalue.kinf_and_spectrum_homogeneous` — `orpheus/derivations/common/eigenvalue.py:98` — literally `M = np.linalg.solve(A, F)` at `:131`. This is the matrix whose dominant eigenpair every `homo_*` reference reports.
  - `orpheus.homogeneous.solver.solve_homogeneous_infinite` — `orpheus/homogeneous/solver.py:150` — the production spelling: `K = MatrixInverseOperator(loss) @ production`, materialised by `K.as_matrix(basis_shape=(ng,1))` into `[K] = A⁻¹F`.
  - `orpheus.numerics.eigenvalue.direct_eigenvalue` — `orpheus/numerics/eigenvalue.py:594` — the `(A,F)`-posed sibling that forms the same resolvent via `np.linalg.solve`. ⚠ `[M]` it has **zero production consumers** since taxonomy step 5b (the page says so at `:1109-1113`) and is retained as an engine + test oracle — a real implementer that the claiming tests reach only through `test_kinf_matches_direct_eigenvalue_engine_of_the_same_pair`.
- **secondary** (see preamble caveats): `derive_kinf_mg_matrix_form` — `M = A_inv * chi * nuSf` at `k_inf_derivations.py:604`, symbolic.
- **NOT an implementer**: `orpheus.derivations.common.eigenvalue.kinf_and_adjoint_spectrum_homogeneous` (`:157`) — it forms the **daggered** resolvent `(Aᵀ)⁻¹Fᵀ`, explicitly "NOT `(A⁻¹F)ᵀ`" per its own comment. Different object.
- **confidence**: **high**.

## `two-group-charpoly`  (8 claims)

- **verdict**: **UNSURE** — the honest answer is one weak implementer or
  `NOTHING:canonical-form`, and I am not able to settle it without a ruling
- **what the equation says**: the characteristic equation of the 2×2 `M`,
  `λ² − (M₁₁+M₂₂)λ + (M₁₁M₂₂ − M₁₂M₂₁) = 0`.
- **the positive evidence for NOTHING** (I searched; this is not absence of a
  quick hit):
  - `grep -rn "charpoly|char_poly" orpheus/ --include='*.py'` → **2 hits**, both `orpheus/derivations/common/homogenization.py:422-423`, and both a *different* equation: `det((A₀+εδA) − μ(F₀+εδF))` for a synthetic rational **perturbation pencil** in the first-order eigenvalue-shift proof (`derive_first_order_k_shift`). Different matrix, different unknown (`μ = 1/k`), different purpose.
  - `grep -rn "np\.roots|np\.poly[^n]|Poly\(" orpheus/ --include='*.py'` → the only `sp.Poly` is `k_inf_derivations.py:256`, a quadratic in **k** from Sood's `det(M(k)) = 0` (Sood Eq 25) — a different posing from `det(M − λI)`. Every other hit is `np.polynomial.legendre.leggauss`.
  - The production path never forms a polynomial: `dominant_eigenpair` calls `np.linalg.eig` (LAPACK QR), `orpheus/numerics/eigenvalue.py:577`.
- **the one candidate implementer**:
  `orpheus.derivations.continuous.fn_method.origins.k_inf_derivations.derive_kinf_mg_matrix_form`
  — `orpheus/derivations/continuous/fn_method/origins/k_inf_derivations.py:516` — `M.eigenvals()` (`:604`) makes SymPy build `M`'s characteristic polynomial internally, and the function then asserts the rank-1 factorisation `λ(λ − tr M) = 0` explicitly (one zero eigenvalue, the other equal to the trace) — the page's own factored form. But (a) SymPy forms the polynomial, not ORPHEUS code, and (b) none of the 8 claiming tests execute it.
- **if the verdict is NOTHING**, the kind is `canonical-form`: a form exhibited
  to show that `k∞` is the dominant root of a quadratic; the production path
  and the reference path both go straight to an eigensolver.
- **confidence**: **low**, and what would change it is a ruling on whether
  "SymPy computes the char poly inside `eigenvals()`" counts as an
  implementation. Independently of that ruling, the finding that matters is the
  ⛔ doc-drift in the preamble: the marker's stated warrant ("derived
  symbolically via exactly those equations") is false.

## `two-group-roots`  (8 claims)

- **verdict**: DECLARABLE
- **what the equation says**: the two roots of that quadratic,
  `λ± = [(M₁₁+M₂₂) ± √((M₁₁−M₂₂)² + 4M₁₂M₂₁)]/2`; `λ₊` is `k∞` and (for the
  rank-1 `M`) `λ₋ = 0`.
- **why this one IS declarable where `charpoly` is not**: the equation defines
  an OBJECT — the pair of eigenvalues of the 2×2 `M` — and code that computes
  that pair implements it, whatever algorithm it uses. `np.linalg.eig(M)` on a
  2×2 returns exactly `λ±`; the selection of `λ₊` is the next equation
  (`keff-update`). The quadratic-formula *route* is not taken, but the
  quantity is computed, on the production path, by every claiming test.
- **implementers** (all verified to resolve):
  - `orpheus.numerics.eigenvalue.dominant_eigenpair` — `orpheus/numerics/eigenvalue.py:519` — `eigvals, eigvecs = np.linalg.eig(M)` at `:578` computes both roots before selecting the dominant. The shared Perron–Frobenius home; on the production path.
  - `orpheus.derivations.common.eigenvalue.kinf_and_spectrum_homogeneous` — `orpheus/derivations/common/eigenvalue.py:98` — its own `np.linalg.eig(M)` at `:133` (the reference tier, deliberately NOT routed through `dominant_eigenpair` — see `test_kinf_exact`'s independence claim).
- **secondary** (see preamble caveats): `derive_kinf_mg_matrix_form` — `M.eigenvals()` plus the explicit `λ₋ = 0` rank-1 assertion, which is the only place the page's "the second eigenvalue is `λ₋ = 0`" claim is checked.
- **confidence**: **medium-high**. What would change it: if "implements" is
  read as "executes this formula" rather than "computes this quantity", the
  two primary entries fall away and this collapses onto the same
  `UNSURE`/`canonical-form` footing as `two-group-charpoly`. The two labels
  must be ruled on TOGETHER — they are one derivation step split across two
  `.. math::` blocks.

---

## `keff-update`  (8 claims)

- **verdict**: DECLARABLE
- **rationale comment on the page**: none as a `vv-status` comment, but the
  page carries an explicit authored **re-binding note** immediately below the
  equation (`infinite_medium.rst:1140-1155`): *"The labels `fixed-source-solve`
  and `keff-update` historically named the per-iteration fixed-source solve …
  and the production/absorption eigenvalue ratio of the retired power
  iteration. They are retained on the **direct** analogues: … and the
  eigenvalue extraction `k∞ = λ_max(M)`."* That note IS the rationale and it
  names the implementer class unambiguously.
- **what the equation says**: `k∞ = λ_max(M)` and `φ` = the corresponding
  dominant right eigenvector of `M = A⁻¹F`, selected as the eigenpair with the
  largest real eigenvalue and sign-normalised non-negative (Perron–Frobenius).
- **implementers** (all verified to resolve):
  - `orpheus.numerics.eigenvalue.dominant_eigenpair` — `orpheus/numerics/eigenvalue.py:519` — ⭐ **the one home**. Does exactly what the equation says and nothing else: `np.linalg.eig`, `argmax(real)`, complex rejection, sign-normalisation. The page names it as "the shared Perron–Frobenius extraction primitive" and says "This validation has **one home** … every direct spelling delegates to it".
  - `orpheus.derivations.common.eigenvalue.kinf_and_spectrum_homogeneous` — `orpheus/derivations/common/eigenvalue.py:98` — the reference tier's own inline realisation of the same selection (`:134-137` argmax + sign-normalise), deliberately NOT routed through `dominant_eigenpair` so the oracle stays structurally independent of the engine.
  - `orpheus.numerics.eigenvalue.direct_eigenvalue` — `orpheus/numerics/eigenvalue.py:594` — the `(A,F)`-posed engine; forms the resolvent then delegates the extraction. Zero production consumers (see `two-group-M`), retained as engine + oracle.
  - `orpheus.homogeneous.solver.solve_homogeneous_infinite` — `orpheus/homogeneous/solver.py:150` — the call site that *is* the direct analogue the note re-binds the label to (`k_inf, phi = dominant_eigenpair(K.as_matrix(...))`, `:202`). Include it if call sites count; drop it if only the computing symbol does.
- **candidate I deliberately excluded**:
  `orpheus.derivations.common.eigenvalue.kinf_and_adjoint_spectrum_homogeneous`
  (`:157`) — it also extracts `λ_max`, but of the **daggered** resolvent
  `(Aᵀ)⁻¹Fᵀ`. Same `k`, different operator; the equation is stated for `M`.
- **NOT an implementer**: `orpheus.numerics.eigenvalue.power_iteration`
  (`:374`) and the classical `k = (νΣ_f·φ)/(Σ_a·φ)` ratio. The page's own note
  says the iterative ratio is what the label was retired FROM; the
  production/absorption form is carried by `absorption-xs`
  (`Mixture.absorption_xs`, which already has an inferred edge there).
- **confidence**: **high** — this is the cleanest of the twelve. The page names
  the function, the function does only this, and all 8 claiming tests execute it.

## `number-density`  (3 claims)

- **verdict**: DECLARABLE
- **rationale comment on the page**: none, but the CLAIMING TEST names the
  implementer verbatim in its own module docstring:
  `tests/data/test_cross_section_data.py:10` — *"``number-density``  —
  :func:`orpheus.data.macro_xs.recipes._number_density`"*, and again at
  `:364`. That is an authored declaration sitting one file away from the graph.
- **what the equation says**: `N_i = ρ_i /(m_u A_i)`, with the `1e-24`
  cm²/barn factor converting to the library's `1/(barn·cm)` units.
- **implementers** (verified to resolve):
  - `orpheus.data.macro_xs.recipes._number_density` — `orpheus/data/macro_xs/recipes.py:21` — the whole function is this formula: `rho = density_g_cm3 * 1e-24; return rho / _AMU_TO_G / molecular_weight_amu`, with `_AMU_TO_G = 1.660538e-24` at `:18` (the page's `m_u`).
- **complete?** yes — `[M]` `grep -rn "1.660538|number_density" orpheus/` returns
  `recipes.py` only; the eight other hits are all *call sites* of this function
  inside the recipe builders (`aqueous_uranium`, `pwr_like_mix`, and the
  per-material recipes at `:49, :84, :98, :112, :150, :173, :202`). No second
  transcription exists.
- **confidence**: **high**. The only thing that would change it is a ruling
  that the recipe call sites should carry the edge too; they consume the value,
  they do not compute it.

## `normalisation`  (1 claim)

- **verdict**: DECLARABLE
- **rationale comment on the page**: none. The prose above the equation names
  the implementer directly: *"After the eigensolve,
  :func:`~orpheus.homogeneous.solver.solve_homogeneous_infinite` normalises the
  flux so that the **fission** production rate is 100 n/cm³/s"*
  (`infinite_medium.rst:1348-1352`).
- **what the equation says**: rescale the (scale-arbitrary) eigenvector so the
  fission production rate is exactly 100 n/cm³/s — `φ ← φ · 100/(νΣ_f·φ)`,
  with **only** `νΣ_f` in the denominator (no `(n,2n)`, which is a loss-side
  transfer).
- **implementers** (verified to resolve):
  - `orpheus.homogeneous.solver.solve_homogeneous_infinite` — `orpheus/homogeneous/solver.py:150`; the statement is `:214`, `phi = phi * (100.0 / production.evaluate(phi.reshape(ng, 1)))`. `[M]` the **sole** site: `grep -rn "100\.0 /|100 /" orpheus/` returns this line plus `fuel/solver.py` strain-percent conversions and one `matpro.py` correlation — nothing else normalises a flux.
  - `orpheus.transport.reaction_rate_functional.IntegratedReactionRate.evaluate` — `orpheus/transport/reaction_rate_functional.py:169` — computes the denominator `νΣ_f·φ` on the meshless unit-volume cell (the `φ†=1` degenerate of the homogenisation bilinear). Declare it if the denominator counts as part of the equation; it is a general reaction-rate functional, so it is the weaker of the two.
- **confidence**: **high** for `solve_homogeneous_infinite`, **medium** for
  `IntegratedReactionRate.evaluate` (general-purpose functional, not
  normalisation-specific).

## `resolvent-object-gate`  (1 claim)

- **verdict**: DECLARABLE — but it is a **gate**, not a production formula, and
  that changes what the edge should point at
- **rationale comment on the page**: none on this label. But its *sibling*
  `resolvent-similarity`, 36 lines above, carries one, and it is directly about
  this equation (`infinite_medium.rst:1190-1194`): *"Structural linear-algebra
  identity …, NOT a solver claim; explains why the factor-order/transpose
  mutations are spectrally invisible. **Verifiable content is the object-level
  matrix gate ``test_K_operator_as_matrix_is_the_resolvent`` (rtol=1e-12)
  named below.**"* — i.e. the author explicitly routes the verifiable content
  of that pair of equations to THIS one and names the gate.
- **what the equation says**: the materialised operator-algebra resolvent
  `[K]` equals `np.linalg.solve(A, F)` to `rtol = 1e-12` — the object-level
  pin that catches the two spectrally-invisible mutations (factor-order swap,
  transpose), which no `k`-level gate can see.
- **implementers** (all verified to resolve) — the equation is an identity
  between two independently-produced objects, so both sides are implementers:
  - LHS `[K]`: `orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator` — `orpheus/numerics/matrix_inverse_operator.py:95` — the `A⁻¹` factor of the product.
  - LHS `[K]`: `orpheus.numerics.operator.LinearOperator.as_matrix` — `orpheus/numerics/operator.py:952` (class at `:588`; the apply-to-basis materialisation the `OperatorProduct` uses). ⚠ this symbol already carries INFERRED edges to `removal-matrix`, `fission-matrix` and `matrix-eigenvalue` — pure "matrix"-token noise there; here the edge would be genuine, which is worth noting when the noise gets cleaned.
  - LHS `[K]`: `orpheus.numerics.operator.OperatorProduct` — `orpheus/numerics/operator.py:1578` — the `@` composition that IS `K = A⁻¹ @ F`.
  - RHS: `orpheus.numerics.eigenvalue.direct_eigenvalue` — `orpheus/numerics/eigenvalue.py:594` — the tree's production spelling of `np.linalg.solve(A, F)`, i.e. the reference side of the identity.
- **⚠ the alternative reading, which I think is defensible and I am flagging
  rather than choosing**: this equation is the *assertion* of
  `tests/homogeneous/test_homogeneous.py:333`
  `test_K_operator_as_matrix_is_the_resolvent`, which already carries the
  single `tests` edge. Under that reading the equation is
  `NOTHING:canonical-form` (a gate criterion, not a computed quantity), and
  declaring implementers makes the one claim self-corroborating (the test
  executes the code it is being adjudicated against, trivially). That is not
  wrong, but it is weak evidence, and the caller should know the edge buys
  little here.
- **confidence**: **medium**. What would change it: a project ruling on whether
  a labelled equation whose content is a *gate criterion* takes `implements`
  edges at all. There are at least three such on this page
  (`resolvent-object-gate`, and by the same argument the identity halves of
  `resolvent-similarity` and `sigs-in-scatter-transpose`), so the ruling is
  reusable.

---

# Findings beyond the 12 rows (they change how the wider campaign should run)

## F1 ⭐ `equation_labels=` is a large authored declaration surface the graph does not read

`[M]` **58 references across 16 files**, all under `orpheus/derivations/`
(`grep -rn "equation_labels" orpheus/ --include='*.py'`), carried by two types:
`VerificationCase.equation_labels` and
`ContinuousReferenceSolution.equation_labels`. An AST scan of the literal
tuple/list keyword arguments alone (missing the `labels += [...]` idiom used by
`cases/sn.py`, `cases/moc.py`, `flat_source_cp/*`) already yields **44 distinct
equation labels** over 18 call sites.

This is the code-side twin of the `.. (vv-status rationale)` comment the brief
sends you to look for — and on this page it is where the answer actually lives:
it settles 8 of the 12 labels directly and is the ONLY authored statement about
the `two-group-*` family. Any fanout over the remaining ~80 equations should
grep `equation_labels` before anything else.

⚠ It is a CLAIM, not a measurement (see F2), so it is a lead to verify, not an
answer to transcribe.

## F2 ⛔ Present-tense-false mechanism claims, in three places, on one derivation

`derive_2g`'s docstring: *"2-group infinite medium eigenvalue **via
characteristic polynomial**"*. `derive_2g_continuous`'s `Provenance.
derivation_notes`: *"Solved via characteristic polynomial (closed form for
2x2)"*. `tests/homogeneous/test_homogeneous.py:24-28`: *"the analytical k_inf
is derived symbolically via exactly those equations"* — offered as the warrant
for the `two-group-charpoly` / `two-group-roots` markers.

`[M]` all three are false at HEAD: `derive_2g` builds `A_sym`/`F_sym`
symbolically and then calls `kinf_and_spectrum_homogeneous`, which is
`np.linalg.solve(A, F)` + `np.linalg.eig(M)`
(`orpheus/derivations/common/eigenvalue.py:131-133`). No characteristic
polynomial is formed anywhere on that path. Fix the prose independently of any
declaration — the marker's stated warrant is the thing that is wrong, not
necessarily the marker.

## F3 `.. verifies:: <label> :by: <dotted>` mints a `tests` edge, NOT `implements`

`[M]` the corpus has **exactly one** such directive
(`docs/theory/foundations/infinite_medium.rst:404`, `one-group-kinf` `:by:`
`derive_1g`) — `grep -rn "^.. verifies::" docs/theory/` → 1 hit. It produces
`derive_1g --tests--> math:equation:one-group-kinf` with `source=directive`.
So the doc-side declaration channel exists, is used once, and is wired to the
wrong relation for this campaign's purpose: `one-group-kinf` still has **zero**
declared `implements` and two INFERRED ones
(`kinf_and_spectrum_homogeneous`, `kinf_from_cp` — pure "kinf"-token matches).
If a doc-side `implements` declaration is wanted, that directive is the natural
place to extend.

## F4 Inferred noise already sitting on this page's *other* equations

Not in my 12, but adjacent and worth cleaning in the same pass — every one of
these is a name-token guess:

- `fixed-source-solve` ← `orpheus.data.macro_xs.sigma_zeros.solve_sigma_zeros` — a *sigma-zero fixed-point iteration*, nothing to do with the loss-matrix solve. Pure "solve" token.
- `removal-matrix` ← `LinearOperator.as_matrix`, `MatrixInverseOperator` — "matrix" token. The true implementers are `_assemble_loss_operator` and `_infinite_medium_matrices` (both listed above under `two-group-A`).
- `matrix-eigenvalue` ← 9 edges including `power_iteration` and `rayleigh_quotient_iteration`, neither of which is on the homogeneous path.
- `fission-matrix` / `fission-source` ← `FissionOperator` — this one is *correct*, by accident.

## F5 `Mixture.balance_residual` is a decoy for `mg-balance`

`orpheus/data/macro_xs/mixture.py:194` computes the **total-cross-section**
definitional identity `Σ_t = Σ_c + Σ_L + Σ_f + rowsum(Σ_s0) + rowsum(Σ_2n)` —
no flux, no `k`. The brief's "check for `balance_residual`-shaped members"
heuristic lands here by name; the member that actually transcribes the neutron
balance is `orpheus.cp.solver.CPSolver._compute_balance_residual`
(`orpheus/cp/solver.py:675`).

## F6 No `Σ_r` symbol; and a four-way twin path under `mg-balance`

- `[M]` the removal cross section `Σ_r,g = Σ_t,g − Σ_s(g→g)` that
  `two-group-Ainv` / `two-group-M` are written in terms of has **no named
  symbol** in `orpheus/`. The value is spellable as
  `Mixture.absorption_xs + Mixture.out_scattering_xs` and is never spelled.
- `[M]` the in-scatter term `Σ_g' Σ_s(g'→g) φ_g'` has **four independent
  transcriptions**: `IsotropicScattering` (homogeneous + diffusion),
  `LegendreMomentScattering` (SN only — `orpheus/transport/operators/scattering.py`
  is imported by `orpheus/sn/` and nothing else), `cp/solver.py:511`'s
  `_scat_mats[mid].T @ φ`, and `moc/core.py:183`'s `sig_s0[i].T @ phi`.
  A Cardinal-Rule-2 item, and the structural reason `mg-balance`'s 70 claims
  cannot be served by any single implementer.

---

# Resolution check — all 26 distinct symbols named above

`[M]` run against `/Users/rodrigo/git/nuclear/ORPHEUS/.nexus/graph.db`
(2026-08-18), suffix-match on the dotted path, node type restricted to
`{function, method, class, data}`: **26 of 26 resolve, 0 unresolved.**
No module node appears anywhere in this document.

| node id | equations |
|---|---|
| `py:function:orpheus.homogeneous.solver.solve_homogeneous_infinite` | `mg-balance`, `inf-hom-balance`, `two-group-M`, `keff-update`, `normalisation` |
| `py:function:orpheus.homogeneous.solver._assemble_loss_operator` | `mg-balance`, `inf-hom-balance`, `two-group-A` |
| `py:method:orpheus.transport.operators.isotropic_scattering.IsotropicScattering.apply` | `mg-balance` |
| `py:method:orpheus.transport.mesh.material_xs_field.MaterialXSField.apply_p0_in_scatter` | `mg-balance` |
| `py:method:orpheus.transport.operators.fission.FissionOperator.apply` | `mg-balance` |
| `py:class:orpheus.transport.operators.fission.FissionOperator` | `two-group-F` |
| `py:function:orpheus.derivations.common.eigenvalue._infinite_medium_matrices` | `mg-balance`, `inf-hom-balance`, `two-group-A`, `two-group-F` |
| `py:function:orpheus.derivations.common.eigenvalue.kinf_and_spectrum_homogeneous` | `mg-balance`, `inf-hom-balance`, `two-group-M`, `two-group-roots`, `keff-update` |
| `py:function:orpheus.derivations.common.eigenvalue.kinf_homogeneous` | `inf-hom-balance` |
| `py:method:orpheus.cp.solver.CPSolver._compute_balance_residual` | `mg-balance` (tier 2) |
| `py:method:orpheus.moc.core.MOCSolver.solve_fixed_source` | `mg-balance` (tier 2) |
| `py:class:orpheus.transport.operators.scattering.LegendreMomentScattering` | `mg-balance` (tier 2) |
| `py:class:orpheus.transport.operators.scattering.N2NMomentOperator` | `mg-balance` (tier 2) |
| `py:class:orpheus.transport.operators.scattering.ScatteringOperator` | `mg-balance` (tier 2) |
| `py:function:orpheus.derivations.continuous.analytical.homogeneous.derive_1g` | `inf-hom-balance` |
| `py:function:orpheus.derivations.continuous.analytical.homogeneous.derive_1g_continuous` | `inf-hom-balance` |
| `py:function:orpheus.derivations.continuous.analytical.homogeneous.derive_2g` | `two-group-A`, `two-group-F` |
| `py:function:orpheus.derivations.continuous.fn_method.origins.k_inf_derivations.derive_kinf_mg_matrix_form` | `two-group-A/Ainv/M/charpoly/roots` (secondary) |
| `py:class:orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator` | `two-group-Ainv`, `resolvent-object-gate` |
| `py:method:orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator.as_matrix` | `two-group-Ainv` |
| `py:function:orpheus.numerics.eigenvalue.direct_eigenvalue` | `two-group-M`, `keff-update`, `resolvent-object-gate` |
| `py:function:orpheus.numerics.eigenvalue.dominant_eigenpair` | `two-group-roots`, `keff-update` |
| `py:function:orpheus.data.macro_xs.recipes._number_density` | `number-density` |
| `py:method:orpheus.transport.reaction_rate_functional.IntegratedReactionRate.evaluate` | `normalisation` |
| `py:method:orpheus.numerics.operator.LinearOperator.as_matrix` | `resolvent-object-gate` |
| `py:class:orpheus.numerics.operator.OperatorProduct` | `resolvent-object-gate` |

---

# Summary

Ten of the twelve are DECLARABLE with high or medium-high confidence and I
believe the enumerations are complete: `mg-balance`, `inf-hom-balance`,
`two-group-A`, `two-group-F`, `two-group-M`, `two-group-roots`, `keff-update`,
`number-density`, `normalisation`, `resolvent-object-gate`. **`two-group-charpoly`
I could not settle** — the positive evidence is that nothing in ORPHEUS forms
the quadratic (`charpoly` appears twice in the tree, both for a different
perturbation pencil; the only `sp.Poly` is Sood's `det(M(k))=0`, a quadratic in
`k` not in `λ`; the production path goes straight to `np.linalg.eig`), so the
choice is between one weak implementer (`derive_kinf_mg_matrix_form`, where
SymPy — not ORPHEUS — forms the polynomial inside `eigenvals()`, and which none
of the 8 claiming tests executes) and `NOTHING:canonical-form`; it turns on a
project ruling about whether "implements" means *computes the quantity* or
*executes the formula*, and that same ruling also decides whether
`two-group-roots` and `two-group-Ainv` keep their primary implementers, so all
three should be ruled on together. **`resolvent-object-gate` is settled only
conditionally** — it is a gate criterion rather than a computed quantity, so
declaring implementers makes its single claim self-corroborating; the edge is
honest but buys little, and a ruling on gate-shaped equations would apply to at
least two more labels on this same page. Finally, two things about `mg-balance`
should be decided before landing anything: its 70 claims are `[M]` 26 SN / 26
CP / 12 homogeneous / 5 MC / 1 MoC, there is no implementer shared across those
families (four independent in-scatter transcriptions), so a Tier-1-only
declaration will refute 58 of the 70 on execution evidence — and the MC family's
5 claims I did not investigate at all, because whether a sampled estimator
"implements" a balance equation is a judgement I would rather surface than make.
