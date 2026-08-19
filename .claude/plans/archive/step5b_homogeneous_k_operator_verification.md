# Step 5b verification plan — homogeneous solver onto the full K = A⁻¹F operator spelling

**Task #138 · branch `refactor/inverse-as-operator` @ `cb62310` · design-before-implementation.**
Canonical gate: `.venv/bin/python -O -m pytest <paths> -m "not slow" -q -rfE --timeout=300 -p no:xdist -p no:cacheprovider`.

This plan is the verification net for the carve. The design is settled (see the dispatch
brief); the job here is the gate list + cadence + teeth. It is *done* only when every gate
can RED under `-O` (§0.5 standing discipline).

---

## 0. The carve in three sub-changes (SC), and what each activates

| SC | Change | New structural element |
|----|--------|------------------------|
| **SC1** | Extract `dominant_eigenpair(M, *, imag_tol=1e-9) -> (float, ndarray)` in `eigenvalue.py` (eig + argmax-real + complex-reject + real-cast + sign-normalize). `direct_eigenvalue(A, F)` KEEPS its (A,F) shape-validation + `resolvent = np.linalg.solve(A, F)`, then `return dominant_eigenpair(resolvent, imag_tol=imag_tol)`. | A public eig-extraction entry point consumed by the K-path WITHOUT direct_eigenvalue's (A,F) validation. |
| **SC2** | `_assemble_loss_matrix` → `_assemble_loss_operator(mat_xs)` returning the `OperatorSum` `collision − k_iso` UN-densified. | (rename + drop the `.as_matrix()` at the return.) |
| **SC3** | `solve_homogeneous_infinite` spells `K = MatrixInverseOperator(loss, basis_shape=(ng,1)) @ production`; `k_inf, phi = dominant_eigenpair(K.as_matrix(basis_shape=(ng,1)))`. | `OperatorProduct(MatrixInverseOperator, FissionOperator).as_matrix()` — the first production consumer of `MatrixInverseOperator`, and the first `OperatorProduct`-of-(inverse, dyad) materialization. |

**Dependencies.** SC1 ⟂ SC2 (independent, parallelizable). SC3 is the JOIN (needs both
`dominant_eigenpair` and `_assemble_loss_operator`). Land as three sub-commits or one; gates
gate after their SC.

### 0.1 Claim-layer + pillar gate (§1.5) — declared up front

- **`dominant_eigenpair` / `direct_eigenvalue` gates** — pure-linear-algebra machinery.
  Claim layers: *flux-shape* (eigenVECTOR recovery, pinned by the intrinsic law `Mφ = kφ`)
  + *eigenvalue* (dominant-VALUE, closed-form **pillar 1**). References are hand-derived
  closed-form matrices → structural independence **by construction**. The 1-group-degeneracy
  Cardinal Rule does NOT bind (no transport eigenvalue is claimed here — this is
  model-agnostic LA; the same exemption the existing `TestStandardEigenproblem` enjoys).
- **Homogeneous end-to-end gates** — *eigenvalue* claim on the infinite-medium spectrum.
  Pillar: **closed-form** (`test_kinf_exact`, the SymPy `k_inf`), the ONLY structurally-
  independent value anchor. ≥2G discipline satisfied (`homo_2eg`, `homo_4eg`,
  `homo_2eg_n2n`); the asymmetric-n2n case de-vacuums the `2Σ₂ᵀ` loss term.
- **Config-blindness (§0.6) is largely N/A here** and this is worth stating so a future
  reader does not "add a heterogeneous case": the infinite-medium problem has NO
  spatial/angular/redistribution term to null — it reduces to the pure energy-balance
  eigenvalue `Aφ = (1/k)Fφ`. The relevant disciplines are ≥2G (satisfied) and
  asymmetric-n2n (satisfied). Flat-flux / curvilinear / mesh-refinement blindnesses do not
  apply because the terms they guard do not exist in this SUT.

---

## 1. Decisive empirical finding — the equivalence class is PRINCIPLED-EQUIVALENCE, gate at `rtol=1e-12` NOT `==`

Measured on this host (numpy 2.4.4, scipy 1.17.1) for a representative `(A, rank-1 F)`:

- OLD resolvent `np.linalg.solve(A, F)` vs NEW per-column `lu_solve(lu_factor(A), F[:,j])`:
  **bit-identical (max |Δ| = 0.0)**; `k` **bit-identical**.

The numpy/scipy backends share LAPACK here and per-column `getrs (nrhs=1)` is arithmetically
the batched `getrs (nrhs=n)`. **But the gate MUST NOT encode byte-equality.** Reasons
(vv-principles §Bit-identity-vs-principled-equivalence):

1. **Bit-identity is an implementation property, not a math property.** The carve
   deliberately changes the LAPACK call sequence (one batched `gesv` → a held `lu_factor`
   + n single-column `lu_solve` backsolves). A byte contract over-fits *this host's* BLAS.
2. **Portability.** A different build (MKL vs OpenBLAS, a future scipy that blocks `getrs`
   differently, or numpy/scipy linking *different* LAPACKs) can drift the resolvent by
   `κ(A)·ULP`. A `==` gate would RED spuriously on such a host — a false regression.
3. **The three criteria all hold** → the correct contract is a *narrowed* tolerance, not
   byte-equality:
   - **principled** — the resolvent is formed by NAMED operators (`MatrixInverseOperator`
     = A⁻¹, `production` = F, `K = A⁻¹F` = `OperatorProduct`); every intermediate is an
     inspectable operator-algebra object.
   - **structurally-independent reference** — the SymPy anchor `test_kinf_exact` (`k_inf =
     λ_max`, 1e-12), a closed-form pillar, unchanged.
   - **FP-non-associativity, dimensionally bounded** — drift ≤ (reduction-depth ≈ ng) ×
     κ(A) × ULP ≈ 1e-14 for ng ≤ 4, well-conditioned homogeneous A. ≪ 1e-12.

**Contract: keff invariance gates at `rtol=1e-12, atol=0`.** It is satisfied bit-identically
today (Δ=0 < 1e-12) AND survives a κ·ULP-drift host. 1e-12 is ≫ the FP bound (~1e-14) and
≪ any real rewire bug (a wrong composition / dropped factor / transposed resolvent moves
keff O(1e-3)+). Bit-identity is preserved ONLY where the reduction tree is genuinely
unchanged (the rate rerouting, §3.4).

---

## 2. Gate list + per-gate cadence

Legend: **[M]** migrate existing · **[R]** retune existing · **[N]** new · **[=]** stays
green unchanged (regression witness).

### After SC1 — `dominant_eigenpair` extraction · file `tests/numerics/test_eigenvalue.py`

| # | Gate | Kind | Pins | Tol | Mutation that REDs it (teeth) |
|---|------|------|------|-----|-------------------------------|
| **G1** | `TestDominantEigenpair` (new class, binds `_dom(M)`) | [N] | The DIRECT public surface `M → (k,φ)`: hand-value+dominant-selection, complex-dominant→raise (neg), real-from-complex→real dtype (pos), sign-normalize, residual law `Mφ=kφ`. | k `rtol=1e-12`; law `rtol=1e-10` | Reuses the file's `_ref_solve(A=I, …)` teeth apparatus (`TestGatesHaveTeeth` already parametrises argmin / drop_sign / drop_real / invert_swap with A=I — that IS `dominant_eigenpair`'s body). Each documented mutation reddens the matching G1 assertion. |
| **G2** | `test_complex_rejection_is_single_home` (+ `…_sign_normalization_…`) | [N] | The ONE-HOME proof: the complex-reject (and sign-flip) live in `dominant_eigenpair` ALONE; `direct_eigenvalue` inherits by delegation, holds no copy. | — (raises / dtype) | **The mutation IS the test** (§2b). In-process monkeypatch neuters `_eig.dominant_eigenpair`; asserts `direct_eigenvalue`'s rejection VANISHES → no duplicate guard. |
| **G3** | `TestStandard/Generalized/SignNormalisation/ComplexSpectrum/IntrinsicProperties/EdgeContracts/GatesHaveTeeth` | [=] | `direct_eigenvalue` signature+behavior unchanged; TRANSITIVELY exercise `dominant_eigenpair`'s body through the delegation. | unchanged | (regression witnesses — a broken extraction reds them.) |
| **G3b** | `TestRayleighQuotientIteration` | [=] | RQI uses `direct_eigenvalue` as oracle — untouched. | unchanged | — |

### After SC2 — `_assemble_loss_operator` rename · file `tests/homogeneous/test_homogeneous.py`

| # | Gate | Kind | Pins | Tol | Teeth |
|---|------|------|------|-----|-------|
| **G4** | `test_assemble_loss_matrix_matches_fused_oracle` (:119) → spells `_assemble_loss_operator(mat_xs).as_matrix(basis_shape=(ng,1))` | [M] | Operator-composed `A = C − K_iso` (apply-to-basis) == fused `diag(Σt) − (Σs0+2Σ2)ᵀ` on `homo_2eg_n2n`. Same oracle, symbol renamed, `.as_matrix()` added. Keeps `verifies("removal-matrix")`. | `atol=1e-12, rtol=0` | Sign-flip in `k_iso` / dropped `2·` on Σ₂ reds it (asymmetric-n2n → non-vacuous). Shares `mat_xs` data ⟹ NOT structurally independent ⟹ **pairs with `test_kinf_exact`**, never replaces it. |

### After SC3 — the K-spelling rewire · file `tests/homogeneous/test_homogeneous.py`

| # | Gate | Kind | Pins | Tol | Teeth |
|---|------|------|------|-----|-------|
| **G5** | `test_kinf_is_the_direct_eigenvalue_of_the_assembled_pair` (:257) | [R] | CROSS-ENGINE: solver `k_inf` (K-path, `lu_solve` resolvent) vs `direct_eigenvalue(A,F)[0]` (`np.linalg.solve` resolvent). Localizes a REWIRE regression to the resolvent-formation boundary. | `rtol=1e-12, atol=0` (was `==`) | A wrong rewire (factor-swap `F@Minv`, wrong `basis_shape`, transposed resolvent) moves k_inf O(1) → reds. Shares `dominant_eigenpair` ⟹ regression localizer, **pairs with `test_kinf_exact`**. |
| **G6** | `test_direct_eigenvalue_is_on_the_homogeneous_call_path` (:203) → spy `dominant_eigenpair` | [R] | Mode-11 liveness: `solve_homogeneous_infinite` CALLS `dominant_eigenpair` (not a still-live inline eig). **HARD requirement — no silent skip.** | — (counter ≥ 1) | The rewired symbol is on the call graph; a routed-around inline eig leaves the counter 0 → `pytest.fail`. |
| **G7** | `test_matrix_inverse_operator_apply_is_on_the_homogeneous_call_path` | [N] | Mode-11 liveness for the FIRST production consumer: `MatrixInverseOperator.apply` (the LU backsolve) executes during the solve. | — (counter ≥ 1; = ng) | A `K.as_matrix()` routed around the inverse operator leaves the counter 0 → fail. Gold-standard in-process wrap sentinel. |
| **G8** | `test_K_operator_as_matrix_is_the_resolvent` | [N] | The GENUINELY-new structural element: `OperatorProduct(Minv, F).as_matrix()` == `A⁻¹F`. Reference `np.linalg.solve(A_fused, F_fused)` (a DIFFERENT primitive than the operator path). | `rtol=1e-12, atol=0` | Factor-swap (`production @ Minv` → `F·A⁻¹ ≠ A⁻¹F`) reds; a dropped/transposed factor reds. Operator-boundary localizer (faster than the end-to-end eig). |
| **G9** | `test_kinf_exact` (:62), `test_post_solve_production_rate_is_100` (:80), `test_rates_via_integrated_reaction_rate_are_bit_identical` (:303), `test_kinf_gate_executes_the_bare_multiplication_arm` (:144), `test_eg_block_*` (:339,:364) | [=] | End-to-end regression witnesses (see §3.4 for why each stays green under φ's ULP drift). | unchanged | — |
| **G10** | Comment block :179–194 ("k_inf — BIT-IDENTICAL") | [doc] | Amend the equivalence-class note: P4-D BIT-IDENTICAL → step-5b **PRINCIPLED-EQUIVALENCE** (§1 three criteria). | — | — |

**Explicit expectation for G9 / `test_kinf_gate_executes_the_bare_multiplication_arm` (:144).**
This Mode-11 gate perturbs `MultiplicationOperator.apply`'s bare-ndarray arm ×1.5 and expects
k_inf to move O(1). Post-carve the collision arm `C = M[Σt]` executes INSIDE
`MatrixInverseOperator.__init__` (the ctor calls `loss.as_matrix(basis_shape=(ng,1))`, which
drives `collision.apply(e_j)` per basis column) → the perturbation flows into the LU
factorization → wrong resolvent → wrong k_inf → **gate stays green (teeth intact)**. CONFIRM
this in the SC3 run: the perturbation site moved from the old inline `_assemble_loss_matrix`
call to the ctor materialization, but the arm is still on the call graph.

---

## 3. Migration / retune spec (exact spellings)

### 3.1 G4 — `test_assemble_loss_matrix_matches_fused_oracle` (:119) [MIGRATE]

Symbol rename + `.as_matrix()` at the call site. Body otherwise verbatim (same oracle, same
`homo_2eg_n2n`, same `atol=1e-12`). Preserve `@pytest.mark.verifies("removal-matrix")`.

```python
from orpheus.homogeneous.solver import _assemble_loss_operator   # was _assemble_loss_matrix
...
A = _assemble_loss_operator(mat_xs).as_matrix(basis_shape=(mix.ng, 1))   # was _assemble_loss_matrix(mat_xs)
```

Consider renaming the TEST to `test_assemble_loss_operator_matches_fused_oracle` for
greppability with the production symbol (`feedback_naming_consistency_greppable`). Non-load-
bearing; do it in the same commit if cheap.

### 3.2 G6 — `test_direct_eigenvalue_is_on_the_homogeneous_call_path` (:203) [RETUNE — the Mode-11 hazard]

**THE headline hazard.** Left as-is, `target = next((n for n in ("DirectEigenvalue",
"direct_eigenvalue") if hasattr(hsolver, n)), None)` returns `None` post-carve (both symbols
leave the solver namespace) → `pytest.skip(...)` → **silent vacuous green** (Mode 11). Rewire
in the SAME commit as SC3:

- Spy `dominant_eigenpair` in BOTH `hsolver` (the solver's own binding — the rewire imports it
  at module level) AND `_eig` (definition site, belt-and-suspenders for an in-function import).
- Convert the pre-impl skip to a **hard requirement**: the carve has landed, so symbol absence
  is a REGRESSION, not a pre-impl state. Keep `_require` (fires under `-O`).

```python
def test_dominant_eigenpair_is_on_the_homogeneous_call_path(monkeypatch):
    """Mode-11: solve_homogeneous_infinite CALLS dominant_eigenpair (the K-path
    eigenvalue primitive), not a routed-around inline eig. HARD requirement —
    symbol absence is a regression post-carve, never a silent skip."""
    import orpheus.homogeneous.solver as hsolver
    _require(hasattr(hsolver, "dominant_eigenpair"),
             "solver.py does not import dominant_eigenpair — the K-spelling rewire "
             "regressed (or the solver still binds direct_eigenvalue).")
    calls: list[int] = []
    def _wrap(ns, name):
        original = getattr(ns, name)
        def spy(*a, **k):
            calls.append(1); return original(*a, **k)
        monkeypatch.setattr(ns, name, spy)
    _wrap(hsolver, "dominant_eigenpair")
    if hasattr(_eig, "dominant_eigenpair"):
        _wrap(_eig, "dominant_eigenpair")
    case = get("homo_2eg_n2n"); mix = next(iter(case.materials.values()))
    result = solve_homogeneous_infinite(mix)
    _require(len(calls) >= 1,
             "solve_homogeneous_infinite did NOT call dominant_eigenpair — the K-path "
             "eigenvalue primitive is not on the call graph (Mode-11 vacuous green).")
    np.testing.assert_allclose(result.k_inf, case.k_inf, atol=1e-12, rtol=0)
```

Rename the function (`direct_eigenvalue` → `dominant_eigenpair`) so a `grep` for the old
symbol does not resurrect it.

### 3.3 G5 — `test_kinf_is_the_direct_eigenvalue_of_the_assembled_pair` (:257) [RETUNE — byte → rtol]

Post-carve `result.k_inf` (K-path, `lu_solve` resolvent) and `direct_eigenvalue(A,F)[0]`
(`np.linalg.solve` resolvent) are the SAME math by two engines. Drop the byte contract; the
`DE`/`fn` pre-impl branch collapses to `direct_eigenvalue` (which now exists and delegates).

```python
def test_kinf_matches_direct_eigenvalue_engine_of_the_assembled_pair():
    """CROSS-ENGINE equivalence (principled, NOT byte): solver k_inf (K-path —
    MatrixInverseOperator lu_solve resolvent) vs direct_eigenvalue(A,F) (np.linalg.solve
    resolvent). Both then extract via the SAME dominant_eigenpair, so this localizes a
    REWIRE regression to the resolvent-formation boundary; it is NOT structurally
    independent on the eig side and PAIRS with test_kinf_exact (the SymPy anchor)."""
    from orpheus.homogeneous.solver import _assemble_loss_operator
    from orpheus.transport.mesh.material_mesh import MaterialMesh
    from orpheus.transport.operators.fission import FissionOperator
    case = get("homo_2eg_n2n"); mix = next(iter(case.materials.values()))
    mat_xs = MaterialMesh.from_materials({0: mix}).material_xs_field()
    A = _assemble_loss_operator(mat_xs).as_matrix(basis_shape=(mix.ng, 1))
    F = FissionOperator.from_solver_data(mat_xs=mat_xs).as_matrix(basis_shape=(mix.ng, 1))
    k_engine = _eig.direct_eigenvalue(A, F)[0]
    result = solve_homogeneous_infinite(mix)
    np.testing.assert_allclose(result.k_inf, k_engine, rtol=1e-12, atol=0)   # was ==
```

Tolerance justification: §1 (bit-identical here, κ·ULP-portable, ≪ any rewire bug). It keeps
localizing a rewire regression to the engine boundary and pairs with the SymPy anchor.

### 3.4 Why the [=] regression witnesses stay green under φ's ULP drift

Post-carve `phi` is the eigenvector of the (possibly κ·ULP-drifted) resolvent → `phi` is NOT
byte-identical to the pre-carve `phi`. Confirm each witness absorbs that:

- **`test_kinf_exact` (:62)** — `|k_inf − case.k_inf| < 1e-12`. The docstring already
  anticipates FP-reduction-tree changes ("principled-equivalence, ~1 ULP"); 1e-12 absorbs the
  κ·ULP drift. **Stays green unchanged** (do not touch — it is the structural anchor).
- **`test_post_solve_production_rate_is_100` (:80)** — φ is normalized to production=100 by
  construction; the 1e-9 band absorbs drift.
- **`test_rates_via_integrated_reaction_rate_are_bit_identical` (:303)** — compares two rate
  computations on the SAME φ (whatever the solver produced). Byte-identical regardless of φ's
  value. Unaffected — the bit-identity is between `IntegratedReactionRate.evaluate(φ)` and
  `νΣf @ φ`, a reduction tree the carve does not touch.
- **`test_eg_block_*` (:339, :364)** — energy-grid diagnostics; φ-densities gate at `rtol`.

---

## 4. New-gate designs

### 4a. `TestDominantEigenpair` — the direct public surface (deliverable 2a)

Placement: a new class in `tests/numerics/test_eigenvalue.py` (the primitive's home). Add a
binding adapter beside `_solve`:

```python
def _dom(M: np.ndarray) -> tuple[float, np.ndarray]:
    fn = getattr(_eig, "dominant_eigenpair", None)
    if fn is None:
        pytest.skip("step-5b PRE-IMPL: dominant_eigenpair not on eigenvalue.py yet.")
    return fn(M)
```

**Lean on existing, do NOT duplicate.** `dominant_eigenpair(M)` is exactly `_solve(np.eye(n),
M)`'s inner body (the resolvent of `A=I` is `M`). So REUSE the file's hand-derived matrices,
hand values, and `_assert_generalized_eigenpair(np.eye(n), M, k, φ)` law helper verbatim —
only the binding surface (`_dom(M)`) is new. The class is a FOCUSED set proving the *standalone
public contract* (calling `dominant_eigenpair` DIRECTLY applies the guards), NOT a full mirror
of `TestStandard/Generalized/Intrinsic` (those cover the shared body transitively through
`direct_eigenvalue`). Six tests, each a distinct claim:

| Test | Matrix (reused) | Claim | Anti-pattern-#11 role |
|------|-----------------|-------|-----------------------|
| `test_hand_value_and_dominant_selection` | `diag([2,5,1])` | k=5 (dominant in the MIDDLE), φ=e₁ | positive |
| `test_symmetric_hand_pair` | `[[2,1],[1,2]]` | k=3, φ∝[1,1] + residual law | positive |
| `test_complex_dominant_raises` | `[[3,-2,0],[2,3,0],[0,0,2]]` (3±2i, 2) | `ValueError` | **negative** (guard raises) |
| `test_real_dominant_from_complex_spectrum` | `blockdiag(5, rot(1,3))` | k=5 real `float`, φ real dtype | **positive** (guard does NOT over-raise) |
| `test_sign_normalization` | `TestSignNormalisation._M` | φ.sum() ≥ 0 | positive |
| `test_residual_law` | `TestIntrinsicProperties`-style generic M | `Mφ = kφ` | intrinsic |

The complex-reject pair (`test_complex_dominant_raises` + `test_real_dominant_from_complex_
spectrum`) is the vv-principles anti-pattern-#11 positive+negative pair — a `check_X` method
tested with BOTH a raising and a non-raising instance.

**Teeth for G1.** The file's `TestGatesHaveTeeth` already carries `_ref_solve(A, F, mut=…)`
parametrised by argmin / drop_sign / drop_real_k / drop_real_phi / invert_swap and proves each
mutation reddens the matching gate — and every teeth case uses **A = I**, so `_ref_solve(I, M,
mut)` IS `dominant_eigenpair(M)` mutated. G1 inherits its teeth from that class with ZERO new
mutation apparatus; add one line to the teeth docstring noting it now doubles as
`dominant_eigenpair`'s teeth.

### 4b. The ONE-HOME proof — G2 (deliverable 2b)

The discriminator between *relocated* (single home) and *duplicated*: does mutating
`dominant_eigenpair` change `direct_eigenvalue`'s behavior? Single home → YES (delegation
propagates the mutation); duplicate → NO (direct_eigenvalue's own copy is unaffected). The
mutation IS the test. In-process monkeypatch (NEVER `git checkout` — the file carries
uncommitted edits, per process-discipline.md); `-O`-safe (`pytest.raises` / `_require`).

```python
def test_complex_rejection_lives_only_in_dominant_eigenpair(monkeypatch):
    """One-home proof: the complex-dominant rejection is in dominant_eigenpair ALONE;
    direct_eigenvalue inherits it by delegation and holds NO duplicate copy."""
    if getattr(_eig, "dominant_eigenpair", None) is None:
        pytest.skip("step-5b PRE-IMPL.")
    A = np.eye(3)
    F = np.array([[3.0, -2.0, 0.0], [2.0, 3.0, 0.0], [0.0, 0.0, 2.0]])  # spectrum 3±2i, 2
    # BASELINE — both surfaces reject (the one home has teeth; the delegated surface inherits):
    with pytest.raises(ValueError):
        _eig.dominant_eigenpair(F)          # direct surface
    with pytest.raises(ValueError):
        _eig.direct_eigenvalue(A, F)        # delegated surface (transitive contract test)
    # NEUTER the one home (drop complex-reject + real-cast):
    def _neutered(M, *, imag_tol=1e-9):
        vals, vecs = np.linalg.eig(np.asarray(M, dtype=float))
        d = int(np.argmax(np.real(vals)))
        return vals[d], vecs[:, d]          # NO raise, NO real()
    monkeypatch.setattr(_eig, "dominant_eigenpair", _neutered)
    # direct_eigenvalue's rejection now VANISHES → it had no own copy (single home):
    k_direct = _eig.direct_eigenvalue(A, F)[0]
    _require(np.iscomplexobj(np.asarray(k_direct)) or np.imag(k_direct) != 0.0,
             "direct_eigenvalue STILL rejected after neutering dominant_eigenpair — it owns "
             "a DUPLICATE complex-guard, not the single home (validation was copied, not "
             "relocated).")
```

`test_sign_normalization_lives_only_in_dominant_eigenpair` mirrors this: neuter the sign-flip,
assert `direct_eigenvalue(I, TestSignNormalisation._M)` returns a NEGATIVE-sum φ (baseline
positive). The two together prove BOTH relocated validations are single-home.

Why the baseline `pytest.raises` on `direct_eigenvalue` matters: it is exactly the
"transitively REDs the existing direct_eigenvalue contract tests" half the brief asks for —
`TestEdgeContracts::test_complex_dominant_contract` and `TestComplexSpectrum` reach the same
raise through the same delegation; G2 exhibits it locally and adds the neuter discriminator.

### 4c. Mode-11 liveness for the NEW call path (deliverable 2c) — TWO spies, both warranted

**G6 (`dominant_eigenpair` spy)** — the eigenvalue primitive. Warranted: the K-spelling calls
`dominant_eigenpair`, NOT `direct_eigenvalue` (which is bypassed entirely). Without the rewire
the :203 test silently skips. See §3.2.

**G7 (`MatrixInverseOperator.apply` spy) — YES, warranted.** Justification (the brief's exact
question):

- The migrated fused-oracle (G4) pins the LOSS-operator *assembly* (`C − K_iso`), NOT the
  inverse/backsolve. The cross-engine gate (G5) and SymPy anchor pin the *value* of the
  resolvent — but a value gate is structurally BLIND to WHICH code produced it. If some other
  path computed the resolvent, G5/anchor stay green. Only a sentinel proves the LU backsolve is
  the genuine producer.
- `MatrixInverseOperator` is the **first production consumer** (per the step-5 docstring); the
  gold-standard §0.5 / sentinel-harness discipline says: for a NEW production line named as
  evidence, sentinel-instrument that exact line — do NOT trust that a green value gate
  exercised it. `OperatorProduct` does not override `as_matrix`, so `K.as_matrix()`
  *structurally must* route `K.apply(e_j) = Minv.apply(F.apply(e_j))` — but "structurally must"
  is exactly the assumption Mode-11 says to VERIFY, not assume.

```python
def test_matrix_inverse_operator_apply_is_on_the_homogeneous_call_path(monkeypatch):
    """Mode-11: the FIRST production consumer of MatrixInverseOperator — its LU backsolve
    apply — executes during solve_homogeneous_infinite (once per resolvent column)."""
    from orpheus.numerics.matrix_inverse_operator import MatrixInverseOperator
    calls: list[int] = []
    raw = MatrixInverseOperator.apply
    def spy(self, x, /, *, initial_guess=None):
        calls.append(1); return raw(self, x, initial_guess=initial_guess)
    monkeypatch.setattr(MatrixInverseOperator, "apply", spy)
    case = get("homo_2eg_n2n"); mix = next(iter(case.materials.values()))
    result = solve_homogeneous_infinite(mix)
    _require(len(calls) >= 1,
             "MatrixInverseOperator.apply never fired during the solve — K.as_matrix() routed "
             "around the inverse operator (Mode-11 vacuous green for the first-consumer claim).")
    np.testing.assert_allclose(result.k_inf, case.k_inf, atol=1e-12, rtol=0)
```

The counter will read `ng` (one backsolve per resolvent column); `>= 1` is the honest liveness
floor (do not over-pin the exact count — a future `as_matrix` override could batch it, and the
liveness claim is "it executes", not "exactly ng times"). Do NOT sentinel with a bare `assert`
(Mode 8: `-O` strips it) — the `_require`/counter pattern fires under `-O`.

### 4d. First-production-consumer intrinsic gate (deliverable 2d) — the call, with the coverage ledger

**Coverage inventory (verified in the numerics M-gates):**

- `MatrixInverseOperator` on an **OperatorSum inner** is ALREADY covered:
  `test_matrix_inverts_what_green_refuses` (`test_matrix_inverse_operator.py:494`) constructs
  `MatrixInverseOperator(sum_op, basis_shape=(4,))` (an `OperatorSum`) and pins its apply vs a
  hand matrix at 1e-12.
- The EXACT homogeneous loss `C − K_iso` **materialization** is ALREADY covered:
  `test_as_matrix_energy_leaves_vs_storage_oracle` (`…:205`) asserts
  `(C − (IsotropicScattering + IsotropicN2N)).as_matrix(basis_shape=(ng,1))` == the fused
  storage oracle on an asymmetric-2G mixture.

So **`MatrixInverseOperator`-with-an-OperatorSum-inner needs NO new intrinsic gate** — it is a
composition of two already-pinned facts: (a) `OperatorSum.as_matrix` correctness (numerics
`:205` + homogeneous G4), (b) LU-factor-and-backsolve of a materialized matrix (numerics M-gates
`:344`–`:518`). No new behavior emerges.

**The genuinely-new structural element** is the `OperatorProduct(Minv, F).as_matrix()`
composition — the `@` between an inverse operator and a rank-1 dyad, materialized to the
resolvent. NO existing test exercises it. **RECOMMEND ONE new intrinsic gate** (G8):

```python
def test_K_operator_as_matrix_is_the_resolvent():
    """The genuinely-new step-5b structural element: OperatorProduct(A⁻¹, F).as_matrix()
    == A⁻¹F. Reference np.linalg.solve(A_fused, F_fused) — a DIFFERENT primitive than the
    operator-algebra path (procedurally independent), pinning the composition at the
    operator boundary (faster than the end-to-end eig)."""
    from orpheus.homogeneous.solver import _assemble_loss_operator
    from orpheus.numerics.matrix_inverse_operator import MatrixInverseOperator
    from orpheus.transport.mesh.material_mesh import MaterialMesh
    from orpheus.transport.operators.fission import FissionOperator
    case = get("homo_2eg_n2n"); mix = next(iter(case.materials.values())); ng = mix.ng
    mat_xs = MaterialMesh.from_materials({0: mix}).material_xs_field()
    loss = _assemble_loss_operator(mat_xs)
    production = FissionOperator.from_solver_data(mat_xs=mat_xs)
    K = MatrixInverseOperator(loss, basis_shape=(ng, 1)) @ production
    A = loss.as_matrix(basis_shape=(ng, 1))
    F = production.as_matrix(basis_shape=(ng, 1))
    np.testing.assert_allclose(
        K.as_matrix(basis_shape=(ng, 1)), np.linalg.solve(A, F), rtol=1e-12, atol=0)
```

Teeth: swap the factor order (`production @ MatrixInverseOperator(loss, …)` → `F·A⁻¹ ≠ A⁻¹F`)
reds; a transposed/dropped factor reds. It also confirms the `basis_shape=(ng,1)` threading
through the composition (the meshless operators carry no `domain`, so an explicit `basis_shape`
is mandatory — a threading bug reds here at the operator boundary, not obscurely at the eig).

This is the intrinsic-property gate the project standard asks for (`feedback_test_intrinsic_
properties`): the K-operator's DEFINING law is "its matrix is the resolvent". Placement:
`tests/homogeneous/test_homogeneous.py` (it consumes homogeneous XS + the solver's own
`_assemble_loss_operator`), grouped with the other assembly gates.

---

## 5. Coverage call — failure modes + honest scope (deliverable 4)

**In play:**

- **Mode 8 (`-O` strips bare `assert`).** Every gate-critical assertion routes through
  `_require` (a `pytest.fail`) or `np.testing.*` — both fire under `-O`. The file already
  carries `_require`. NEVER a bare `assert` in a gate the canonical `-O` invocation runs.
  (The `assert not S_ao.is_invertible` fixture-verification lines in the numerics M-gates are
  a separate concern — those are non-gate fixture checks; not touched here.)
- **Mode 11 (vacuous-green via silent skip).** THE headline hazard: G6 (:203) silently skips
  post-carve if left targeting `direct_eigenvalue`. Resolved by the §3.2 rewire (spy the NEW
  symbol + hard requirement). G7 closes the first-consumer liveness gap. The `_dom`/`_solve`
  pre-impl skips are transient scaffolding — REMOVE them (or leave as hard `_require`) when the
  carve lands in the same commit; a pre-impl skip that survives the landing is a Mode-11 relapse.
- **Bit-identity vs principled-equivalence (the core re-baseline).** §1: keff moves from
  BIT-IDENTICAL (P4-D) to PRINCIPLED-EQUIVALENCE; gate `rtol=1e-12`, document the three
  criteria in G10 (the :179–194 comment amendment). Bit-identity preserved ONLY at §3.4's rate
  rerouting (unchanged reduction tree).

**Not in play (state it, so a reviewer does not manufacture a gate):**

- **1-group degeneracy** — N/A for the LA-machinery gates (no transport eigenvalue claimed);
  satisfied for the end-to-end gates (≥2G via `homo_2eg/4eg/2eg_n2n`).
- **Flat-flux / curvilinear / redistribution / mesh-refinement (§0.6)** — N/A: the
  infinite-medium SUT has no such term to null (§0.1).
- **Mode 7/10 (MMS bias / activated-unconstrained)** — N/A: this is a rewire, not an MMS; all
  terms are constrained by the SymPy anchor + fused oracle.

**Honest scope — what step 5b does NOT verify:**

1. It does NOT re-verify `MatrixInverseOperator`'s M-materialise / M-direct invariants — those
   are the step-5 M-gates' job (inherited, verified on synthetic + OperatorSum inners).
2. It does NOT verify any spatial/angular transport operator — the infinite-medium problem has
   none (it is the energy-balance eigenvalue only).
3. G5 and G4 are NOT structurally independent (G5 shares `dominant_eigenpair`; G4 shares
   `mat_xs`) — they are REGRESSION LOCALIZERS that PAIR with the SymPy anchor `test_kinf_exact`,
   which is the ONE structurally-independent value reference. The net's structural independence
   rests entirely on `test_kinf_exact` (closed-form pillar) + `TestDominantEigenpair` (hand-
   derived LA); everything else localizes regressions faster but inherits its independence.
4. The `rtol=1e-12` on G5/G8 is a PRINCIPLED tolerance, not a proof of correctness — it bounds
   the FP-drift, and correctness of the value flows from the SymPy anchor.

**Structural-independence ledger (the net at a glance):**

| Gate | Reference | Structurally independent? | Role |
|------|-----------|---------------------------|------|
| `test_kinf_exact` | SymPy `k_inf` (closed form) | **YES** (pillar 1) | THE value anchor |
| `TestDominantEigenpair` | hand-derived LA matrices | **YES** (by construction) | eig-extraction correctness |
| G4 fused-oracle | `diag(Σt)−(Σs0+2Σ2)ᵀ` from same `mat_xs` | NO (shared data) | A-assembly localizer |
| G5 cross-engine | `direct_eigenvalue` (shares `dominant_eigenpair`) | NO (shared eig) | rewire localizer |
| G8 K-resolvent | `np.linalg.solve(A,F)` | procedurally (diff primitive) | composition localizer |

---

## 6. Cadence summary (per-SC gate order)

```
SC1 (dominant_eigenpair) ─▶ G1 TestDominantEigenpair · G2 one-home ·
                            G3/G3b existing eigenvalue suite [=]
SC2 (_assemble_loss_operator) ─▶ G4 migrated fused-oracle
SC3 (K-spelling) ─▶ G5 cross-engine [rtol] · G6 dominant_eigenpair spy [hard] ·
                    G7 MatrixInverseOperator.apply spy · G8 K-resolvent intrinsic ·
                    G9 end-to-end witnesses [=] · G10 comment amendment [doc]
```

Full-suite gate after the join: `.venv/bin/python -O -m pytest tests/numerics/test_eigenvalue.py
tests/numerics/test_matrix_inverse_operator.py tests/homogeneous/test_homogeneous.py -q -rfE
--timeout=300 -p no:xdist -p no:cacheprovider` — expect all green, with G1/G2/G6/G7/G8 as the
new/rewired rows and G3/G9 unchanged.

---

## 7. Design objections / watch-items for the implementer

1. **G6 rename is mandatory, not cosmetic.** If the function keeps the name
   `test_direct_eigenvalue_is_on_the_homogeneous_call_path` while spying `dominant_eigenpair`,
   a future grep for the retired symbol misleads. Rename to `…dominant_eigenpair…`.
2. **Do NOT gate G5/G8 on `==`** even though this host is bit-identical (§1). A byte contract
   is a fragile over-fit; the review will (correctly) read it as an unprincipled bit-identity
   claim on a deliberately-changed LAPACK call sequence.
3. **Confirm the `OperatorProduct` space guard stays None-tolerant** for the meshless path:
   `MatrixInverseOperator(loss).domain == loss.codomain`; the meshless `loss` carries no
   `FunctionSpace` (domain/codomain None), and `production.codomain` is None
   (`from_solver_data` without `full_field_space`) → the guard skips. G8's successful
   construction of `K` IS the assertion of this; if `K = Minv @ production` raised
   `IncompatibleOperatorComposition`, a space leaked in — investigate before proceeding
   (a new space name on either operand is the blocker, per
   `operator_space_guard_only_bites_operatorsum`).
4. **`FissionOperator` (b-factor) is not invertible** (rank-1 dyad → `is_invertible=False`), so
   `K.is_invertible=False`. That is CORRECT and harmless — the K-path only `.as_matrix()`-es and
   `.apply()`-es K; it never inverts it. Do not "fix" it.
```

---

## 8. ✅ AS-BUILT STAMP (2026-07-02 — COMMITTED @ `394d8c0`, 10 files +860/−306)

Every gate landed and is green. Deltas vs this plan:

- **G5 teeth row CORRECTED (empirical, T1–T5 at job-tmp `step5b_teeth.py`):** §2's claim
  "factor-swap … moves k_inf O(1) → reds [G5]" is WRONG for the swap and the transpose —
  `F·A⁻¹` is SIMILAR to `A⁻¹F` (`A·(A⁻¹F)·A⁻¹ = F·A⁻¹`) and `eig(Mᵀ) = eig(M)`, so
  **|Δk| = 0.0 exactly** under both (measured): every k-level gate (G5, the SymPy anchor)
  is **designed-green** on them, and **G8 is the committed catcher** (both move the
  MATRIX O(1): swap 1.456, transpose 1.434 on `homo_2eg_n2n`). G5's real teeth are the
  k-moving mutations: dropped factor (|Δk| = 4.926 measured), wrong `basis_shape`
  (ctor-raises). The φ change under the swap (φ→Aφ) is pinned by G8's column equality,
  not by any spectrum gate. Full derivation: homogeneous.rst `spectral-invisibility` §.
- **G1 tooth count:** shipped 7 tests (the 6 planned + `test_non_square_raises` for the
  boundary guard). `_dom`/`_solve` bind DIRECTLY (no pre-impl skip — SUT+tests land in
  one commit); the P4-D-era `DirectEigenvalue` class arm + skip scaffolding in `_solve`
  RETIRED.
- **AUG1 executed:** `_sign_normalised` folded into SC1 (RQI's inline flip + the
  extraction's copy — byte-identical, both callers pinned).
- **Enforcer post-review folds (SHIP, zero code defects):** test-module docstring swept
  (`DirectEigenvalue`→`direct_eigenvalue` + `dominant_eigenpair` named); the
  severed-consumer parenthetical in `direct_eigenvalue`'s docstring repointed (it no
  longer implies the homogeneous solver routes through it); the "one LU backsolve per
  column" narration trimmed from the solver comment (falsifiable by a future batched
  `as_matrix` override — the liveness gate's own caveat).
- **Numbers:** trio 77/0 · blast radius `tests/numerics/ + tests/homogeneous/` 848/0
  (+1 pre-existing unrelated skip) · ratchet EXACTLY 145 (21/99/25) · Sphinx `-E -W`
  exit 0 (archivist pass: homogeneous.rst K-story + spectral-invisibility §,
  api/numerics.rst, operator_algebra.rst consumer-loop, SN Dev-history 5b row, matrix
  regenerated 30→39 / 12→14).