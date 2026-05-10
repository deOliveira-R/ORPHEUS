# Issue #168 Phase C — Spatial-closure alignment plan (PAUSED)

**Status (2026-05-10)**: **PAUSED.** Phase A and Phase B shipped on
`refactor/sn-operator-algebra` (commits `d73ef68`, `5df90df`,
`ef377cd`). Research dispatches for Phase C completed and surfaced a
**structural blocker**: the user's chosen absolute-value cross-check
ground (F_N method on Sood La13511 cases) **does not cover** the 6
deleted curvilinear regression snapshots. Phase C is deferred until
the F_N infrastructure gap is closed; the architectural map and
verification plan from this session are saved below so a fresh
session can resume without redoing the research.

**User decision (this session)**:

1. **Snapshot strategy** — *"Stop working on this and save the plan.
   If there are no cylindrical cases for us to base on, we shall
   make them first, then we come back to the plan."* The user is
   choosing to **build the F_N infrastructure** (cylinder F_N +
   2G F_N) **before** Phase C resumes, rather than pivoting Phase C
   to a weaker absolute-value ground (k_∞ analytical for the 4
   homogeneous-reflective cases plus regression-contract-only for
   the 2 heterogeneous 3-region cases).

2. **Default-flip sequencing** — *"Implement WDD first, probe
   empirically."* Once Phase C resumes: implement the WDD spatial
   closure rewrite first; run Gate 1.1 (per-ordinate flat-flux
   residual) with `MorelMontryAngularSweep` as the candidate
   default; if Gate 1.1 passes, flip the default and regen
   snapshots; if it fails, keep `LegacyTauSymmetricInterpolation`
   as default and ship Phase C as "spatial-closure aligned;
   angular default flip deferred to Phase D" with a follow-up
   issue carrying the empirical evidence.

---

## Why Phase C is paused — the F_N coverage gap

Phase C was originally scoped to regenerate the 6 deleted
curvilinear regression snapshots and verify each k_eff via the F_N
method on the closest Sood La13511 case
(`orpheus/derivations/continuous/sood_registry/`). The explorer
agent's architectural map (§4 below) found:

1. **F_N cylinder solver does not exist.** The package
   [orpheus/derivations/continuous/fn_method/cylinder/](orpheus/derivations/continuous/fn_method/cylinder/__init__.py)
   is a placeholder with a docstring at
   [cylinder/__init__.py:14](orpheus/derivations/continuous/fn_method/cylinder/__init__.py#L14)
   reading *"TODO: populate after literature-researcher delivers
   Westfall–Metcalf 1973 + Sanchez–Ganapol 1983."* The
   `WIDE_SLICE_STUBS` cylinder cases are pytest-skipped explicitly.

2. **2G F_N for sphere/cylinder is not implemented.** The Siewert–
   Thomas 1986 2G F_N machinery is missing. Existing pytest at
   [test_sood_registry_wide_bare_critical.py:185](tests/derivations/test_sood_registry_wide_bare_critical.py#L185)
   literally skips with reason *"pending 2G FN."*

3. **No direct Sood analogue exists for any of the 6 snapshots.**
   The deleted snapshots use ORPHEUS A/B synthetic mixtures with
   **reflective BC**; F_N solves **vacuum BC bare-critical**. From
   the snapshot factory [tests/sn/regression/_generate_snapshots.py](tests/sn/regression/_generate_snapshots.py):

   | Snapshot | Closest Sood case |
   |---|---|
   | `sphere_2g_homogeneous_dd_n20` | None (2G F_N pending) |
   | `sphere_2g_3reg_dd_n40` | None (multi-region not in Sood corpus) |
   | `sphere_2g_p1_aniso_dd_n20` | UD2Oa/Ob/Oc-1-1-SP cousins are 1G P1, materials differ |
   | `cyl_1g_homogeneous_LS4_dd_n20` | `Ua-1-0-CY_STUB` (cylinder F_N missing) |
   | `cyl_1g_homogeneous_product_dd_n20` | Same as above |
   | `cyl_2g_3reg_LS4_dd_n40` | None (2G + multi-region + cylinder F_N triple-blocker) |

The user's decision: rather than weaken Phase C's verification
chain, build the F_N infrastructure first.

## Prerequisite work tracked as separate issues

Phase C resumption is gated on:

- **GH [#170](https://github.com/deOliveira-R/ORPHEUS/issues/170)**:
  implement F_N cylinder solver (Westfall–Metcalf 1973 +
  Sanchez–Ganapol 1983). Module:
  `orpheus/derivations/continuous/fn_method/cylinder/`. Ground via
  `Ua-1-0-CY_STUB`, `PUb-1-0-CY_STUB`, `UD2O-1-0-CY_STUB` from the
  Sood La13511 wide-slice corpus.

- **GH [#171](https://github.com/deOliveira-R/ORPHEUS/issues/171)**:
  implement 2G F_N for sphere/cylinder (Siewert–Thomas 1986).
  Module: `orpheus/derivations/continuous/fn_method/multi_group/`.
  Ground via the 10 Sood `WIDE_SLICE_STUBS` 2G bare-critical cases.

(Both deliverables are `algebra-of-record` Branch 1 / Branch 2
constructions — SymPy reduction + numpy production — and would
naturally be dispatched to the `method-implementer` agent after a
`literature-researcher` memo on each paper.)

When both ship, the Sood-FN coverage matrix becomes:

| Snapshot | F_N ground (post-prereq) |
|---|---|
| `sphere_2g_homogeneous_dd_n20` | Sood 2G sphere F_N (post-Siewert–Thomas) |
| `sphere_2g_3reg_dd_n40` | Still none (multi-region; needs separate analytical ground) |
| `sphere_2g_p1_aniso_dd_n20` | Sood 2G P1 sphere F_N (post-Siewert–Thomas) |
| `cyl_1g_homogeneous_LS4_dd_n20` | Sood 1G cylinder F_N (post-Westfall–Metcalf) |
| `cyl_1g_homogeneous_product_dd_n20` | Same |
| `cyl_2g_3reg_LS4_dd_n40` | Still none (multi-region; needs separate analytical ground) |

The two heterogeneous 3-region cases (`sphere_2g_3reg`,
`cyl_2g_3reg`) remain without a direct F_N analogue even after the
prereq work; Phase C resumption may decide to: (a) replace them
with single-region Sood-grounded analogues, or (b) ship them as
regression-contract-only with the WDD-aligned sweep as the
structurally-independent reference (transitive verification via
the sweep's own L1 chain). **Decision deferred to Phase C
resumption.**

---

## Phase C resumption sequencing (when prereqs ship)

1. **WDD spatial closure rewrite** (the math change itself)
   - Spherical: rewrite the per-ordinate face-flux closure in
     [transport_operator_matvec_spherical](orpheus/sn/operator.py#L554)
     at lines 694, 698-700, 706, 709, 711.
   - Cylindrical: same change at
     [transport_operator_matvec_cylindrical](orpheus/sn/operator.py#L737)
     lines 837, 841-843, 847, 850, 852.
   - `SNStreamingOperator.apply` is a pure dispatcher — needs no
     parallel change.
   - `SNStreamingOperator.apply_transpose` is constructed by **dense
     probing of `apply`** (no parallel transpose code path —
     reciprocity is automatic for any linear apply). See §4.4
     below.

2. **Gate 1.4 → 1.3 → 1.1 → 1.2 verification cycle**:
   - 1.4 — apply linearity (precondition for 1.3)
   - 1.3 — `apply` ↔ `apply_transpose` reciprocity (rtol=1e-12)
   - 1.1 — per-ordinate flat-flux residual (the canonical
     curvilinear bug class probe)
   - 1.2 — apply ↔ sweep face-flux propagation identity (the
     structural-equivalence test)

3. **Empirical default-flip decision** (per user's sequencing):
   - Run Gate 1.1 with `MorelMontryAngularSweep` as the candidate
     default. If pass → flip
     `SNMesh.pole_angular_closure` default
     `LegacyTauSymmetricInterpolation` → `MorelMontryAngularSweep`;
     flip curvilinear `solve_sn_fixed_source` default
     `"source_iteration"` → `"krylov"`. If fail → keep Legacy as
     default; file follow-up issue with empirical evidence; Phase C
     ships as spatial-closure aligned only.

4. **Snapshot regeneration** (with whichever default ships):
   - 6 curvilinear snapshots regenerated.
   - Cross-check via F_N (post-prereq) per the coverage matrix
     above.
   - Heterogeneous 3-region snapshots: decision per
     "(a)/(b) deferred" above.

5. **Marker removal** (Gate 5.1):
   - 4 `@pytest.mark.xfail(strict=True)` ERR-026 tripwires
     (locations enumerated in §4.4 below).

6. **Documentation**:
   - `error_catalog.md` ERR-026: PARTIAL → CLOSED (or "spatial-
     closure aligned" if default flip deferred).
   - `sn_reshape.md` plan annotation.
   - GitHub Issue #168: closure comment.
   - `tests._harness.audit` clean.

---

## §4 — Architectural map (from explorer dispatch, 2026-05-10)

### §4.1 Apply matvec — face-flux sites

[transport_operator_matvec_spherical](orpheus/sn/operator.py#L554) and
[transport_operator_matvec_cylindrical](orpheus/sn/operator.py#L737) share
the same loop pattern. The face-flux sites Phase C touches:

**Spherical** ([:693-711](orpheus/sn/operator.py#L693)):

| Site | Code | Lines | Closure |
|---|---|---|---|
| Interior right face | `psi_right = 0.5 * (fi[:, n, i, 0] + fi[:, n, i + 1, 0])` | [694](orpheus/sn/operator.py#L694) | symmetric arithmetic average |
| Outer right μ>0 | `psi_right = boundary_face_flux_closure(fi[:, n, :, 0], i)` | [698-700](orpheus/sn/operator.py#L698) | Phase A `BoundaryFaceFlux` (DD extrapolation default) |
| Outer right μ<0 | `psi_right = bf_flux[:, n, 1]` | [706](orpheus/sn/operator.py#L706) | BC-resolved |
| Interior left face | `psi_left = 0.5 * (fi[:, n, i - 1, 0] + fi[:, n, i, 0])` | [709](orpheus/sn/operator.py#L709) | symmetric arithmetic average |
| Pole left face (i=0) | `psi_left = 0.0` | [711](orpheus/sn/operator.py#L711) | hard-coded zero (× A[0]=0) |

**Cylindrical** is byte-parallel at lines [837, 841-843, 847, 850, 852](orpheus/sn/operator.py#L837).

`solution_to_angular_flux_cylindrical = solution_to_angular_flux_spherical`
([:734](orpheus/sn/operator.py#L734)) — full reuse.

The sweep's WDD asymmetric form to mirror:
[orpheus/sn/sweep.py:866-867](orpheus/sn/sweep.py#L866) —
`psi_x[n, ii + ix_out, jj, :] = 2.0 * psi_avg - psi_in_x`. Also at
[orpheus/sn/spatial/diamond.py:520-521](orpheus/sn/spatial/diamond.py#L520) —
the curvilinear sweep's analogous WDD recurrence.

### §4.2 Equation map ordering

[build_equation_map_spherical](orpheus/sn/operator.py#L386):

```python
for ix in range(nx):              # outer = spatial
    for n in range(quad.N):       # inner = ordinate
        if ix == nx - 1 and mu_x[n] < -1e-15:
            continue              # skip incoming at outer
        ords.append(n); ixs.append(ix); iys.append(0)
```

**Ordering is `(spatial outer, ordinate inner)`** — ordinate-
interleaved within spatial blocks. **A per-ordinate WDD sweep
needs cells in sweep direction** (`i: 0 → nx-1` for μ>0,
`i: nx-1 → 0` for μ<0), so Phase C must either:

- (a) iterate cells/ordinates in sweep order outside `eq_map.n_eq`
  and scatter results back through `eq_map` indexing into `lhs`,
  or
- (b) rebuild `eq_map` with `(ordinate outer, spatial inner)`
  ordering (bit-identity-breaking for any consumer of the packed
  vector layout — needs survey).

**Recommendation**: path (a). Less invasive; preserves the
packed-vector contract.

### §4.3 Pole face-flux interaction

The pole face flux at i=0 is **hard-coded zero** at
[:711](orpheus/sn/operator.py#L711) (and [:852](orpheus/sn/operator.py#L852)):

```python
if i > 0:
    psi_left = 0.5 * (fi[:, n, i - 1, 0] + fi[:, n, i, 0])
else:
    psi_left = 0.0   # A[0] = 0 makes the spatial term vanish
```

The streaming term is `mu[n] * (A[i+1]·psi_right - A[i]·psi_left) / V[i]`.
At i=0, A[0]=0, so the value of `psi_left` is multiplied by zero
and never enters the result. **Phase C's WDD rewrite preserves
this**: in the rightward sweep (μ>0) the inflow at i=0 is `0.0`
(consistent); in the leftward sweep (μ<0) the inflow at i=nx is
`bf_flux[:, n, 1]` (BC-resolved), and the WDD propagation walks
the chain inward producing the outflow at i=0 which is then
multiplied by A[0]=0.

`bf_flux` array structure ([:527](orpheus/sn/operator.py#L527)):
- `bf_flux[:, n, 0]` — pole face placeholder (cell-centre at i=0)
- `bf_flux[:, n, 1]` — outer-face BC-resolved value for inward μ<0

### §4.4 `apply_transpose` is built by dense probing — no parallel code path

**Critical finding**: [_ensure_dense_matrix](orpheus/sn/operator.py#L1169)
constructs the operator's dense matrix by probing `apply` on each
canonical basis vector:

```python
def _ensure_dense_matrix(self) -> np.ndarray:
    if self._dense_matrix is None:
        n = self.n_unknowns
        mat = np.empty((n, n), dtype=float)
        basis = np.zeros(n, dtype=float)
        for j in range(n):
            basis[j] = 1.0
            mat[:, j] = self.apply(basis)            # [:1194]
            basis[j] = 0.0
        self._dense_matrix = mat
        self._dense_matrix_T = mat.T.copy()
    return self._dense_matrix

def apply_transpose(self, psi: np.ndarray) -> np.ndarray:
    self._ensure_dense_matrix()
    return self._dense_matrix_T @ psi                # [:1230]
```

**Implication**: there is **no parallel transpose code path to
maintain**. Reciprocity holds **by construction** for any linear
`apply`. Phase C touches only `apply`; `apply_transpose`
automatically tracks via dense probing on next call. Class
docstring at [:1217-1216](orpheus/sn/operator.py#L1217) makes this
explicit.

### §4.5 ERR-026 xfail tripwires (4 tests)

| # | File | Line | Test | Reason summary |
|---|---|---|---|---|
| 1 | [tests/sn/test_mms_curvilinear.py](tests/sn/test_mms_curvilinear.py) | 42-58 | `test_sn_spherical_mms_converges_second_order` | "ERR-026 — same root cause as the anisotropic curvilinear MMS" |
| 2 | [tests/sn/test_mms_curvilinear.py](tests/sn/test_mms_curvilinear.py) | 94-106 | `test_sn_cylindrical_mms_converges_second_order` | "ERR-026 — same root cause as the spherical isotropic MMS" |
| 3 | [tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py](tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py) | 70-96 | `test_sn_spherical_aniso_mms_converges_second_order` | "curvilinear sweep WDD wrong fixed point + curvilinear FD operator boundary face-flux first-order on non-constant" |
| 4 | [tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py](tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py) | 135-149 | `test_sn_cylindrical_aniso_mms_converges_second_order` | "Same root cause as the spherical anisotropic case" |

All 4 are `strict=True` → xpass without removal becomes strict-fail
signalling the bug is fixed. The xfail reasons are **precisely
matched to ERR-026**, not generic tripwires; once Phase C closes
the spatial-closure alignment correctly, strict-mode forces their
removal.

### §4.6 Regression-contract mechanism

[tests/sn/regression/test_dd_regression.py:46-51](tests/sn/regression/test_dd_regression.py#L46):

```python
snapshot_file = SNAPSHOT_DIR / f"{case.name}.npz"
if not snapshot_file.exists():
    pytest.skip(f"snapshot {snapshot_file.name} not yet generated; ...")
```

Tolerance: `rtol=1e-12, atol=1e-13`
([:62-70](tests/sn/regression/test_dd_regression.py#L62)).

Snapshot regeneration: `python -m tests.sn.regression._generate_snapshots`
([_generate_snapshots.py:392-422](tests/sn/regression/_generate_snapshots.py#L392)).

---

## §5 — Verification plan (from test-architect dispatch, 2026-05-10)

13 gates across 5 sets. Levels: L0 (term verification), L1 (equation
verification, MMS, k_∞ closed-form), foundation (software
invariants), regression-contract (bit-identity inheritance).

### Gate Set 1 — WDD spatial closure correctness

**Gate 1.1 (L0) — `test_apply_curvilinear_per_ordinate_flat_flux_residual_under_wdd`**

For ψ constant in space and per-ordinate, on a reflective-BC
homogeneous sphere/cylinder with Σ_t=0 the apply output is
bit-identical to zero (atol=1e-13 per packed-vector entry). With
Σ_t≠0 uniform, output equals `Σ_t · ψ` per ordinate (rtol=1e-12).
**Catches**: #2 (variable swap), #3 (missing factor), #4 (wrong
recursion), Signature 1 (curvilinear sweep divergence under
refinement).

Decorators: `@pytest.mark.l0`, `@pytest.mark.verifies("dd-curvilinear-scalar")`,
`@pytest.mark.catches("ERR-026")`.

Parametrisation: spherical + cylindrical; n_ordinates ∈ {4 GL, 8
GL, LS-4, Product 2x4}; n_groups ∈ {2, 4}; nx ∈ {4, 10, 20}; ψ ≡
const broadcast over the unknown vector.

**Gate 1.2 (foundation) — `test_apply_face_fluxes_match_sweep_recurrence_curvilinear`**

For arbitrary input ψ_cells, the per-ordinate face-flux values
implicitly propagated by `apply` (extracted by inverting the
cell-balance equation per cell) match those produced by the
sweep's WDD recurrence at iteration zero (sweep with same Σ_t,
same BC, source `q = apply(ψ_cells)`). Tolerance: `np.array_equal`
(bit-identical). **NOT** a code-to-code agreement claim — this is
a structural identity claim: the WDD diamond relation
`ψ_face_out = 2·ψ_cell − ψ_face_in` is one closed-form rule
implemented identically by both paths.

Decorators: `@pytest.mark.foundation`,
`@pytest.mark.verifies("phase-c-apply-sweep-equivalence")` (NEW
label), `@pytest.mark.catches("ERR-026")`.

**Gate 1.3 (foundation) — `test_apply_apply_transpose_reciprocity_under_wdd`**

`<L·ψ, φ> = <ψ, Lᵀ·φ>` to round-off (rtol=1e-12, atol=1e-13) on
random ψ, φ. **Free** post-Phase-C if the rewrite preserves
linearity (the dense-matrix probe constructs the transpose).

Decorators: `@pytest.mark.foundation`,
`@pytest.mark.verifies("sn-streaming-reciprocity")`.

**Gate 1.4 (foundation) — `test_apply_linearity_under_wdd_curvilinear`**

`apply(α·ψ + β·φ) = α·apply(ψ) + β·apply(φ)` to rtol=1e-13. Defence
in depth for Gate 1.3.

### Gate Set 2 — Bit-identity preservation

- **2.1**: 5 Cartesian regression snapshots pass bit-identical
  (rtol=1e-12, atol=1e-13). Cartesian path uses upwind FD with
  separate dx/dy denominators — the WDD interior-closure rewrite
  is curvilinear-only.
- **2.2**: Phase A 21 foundation tests
  (`tests/sn/spatial/test_boundary_face_flux.py`) pass.
- **2.3**: Phase B 28 foundation tests
  (`tests/sn/spatial/test_pole_angular_closure.py`) pass.
- **2.4**: Phase B 5 L1 flat-flux-identity tests
  (`tests/sn/l1_analytical/test_pole_closure_flat_flux_identity.py`) pass.

### Gate Set 3 — Convergence (L1, MMS)

**Gate 3.1 — Spatial MMS convergence, spherical, with linear-μ
ansatz**

Trial: `ψ_chosen = (A(r) + B(r)·μ)/W` with the linear-in-μ term
**mandatory** to activate the redistribution path (Mode 7 / ERR-026
hide-the-bug pattern explicitly avoided).

**Note from test-architect**: the originally-proposed
`A(r)=sin(πr/R)`, `B(r)=r·(R−r)/R²` **vanishes at r=0** (both A
and B). Use instead `A(r)=cos(πr/(2R))` (=1 at pole), `B(r)=r/R`
(zero at pole, unit at outer). Confirm with literature (Lewis &
Miller §6.4 or Adams–Larsen §III) before implementation.

Refine nx ∈ {10, 20, 40, 80, 160} at fixed n_ordinates=8 GL. Assert
`min(orders[-3:]) > 1.9` for spatial convergence rate.

Decorators: `@pytest.mark.l1`, `@pytest.mark.slow`,
`@pytest.mark.verifies("sn-mms-spherical-aniso-spatial-convergence")` (NEW),
`@pytest.mark.catches("ERR-026")`.

**Gate 3.2 — Cylindrical MMS convergence, with linear-μ_x ansatz**

Cylindrical analogue of 3.1. ψ_chosen = (A(r) + B(r)·μ_x)/W on
cylindrical mesh r ∈ [0.01, 2.0]. Two quadrature families
parametrised: LS-4 and Product 2x4 (flag any quadrature-dependent
constant — Signature 4 / ERR-004).

Decorators: `@pytest.mark.l1`, `@pytest.mark.slow`,
`@pytest.mark.verifies("sn-mms-cylindrical-aniso-spatial-convergence")` (NEW),
`@pytest.mark.catches("ERR-026")`.

**Gate 3.3 — Angular convergence at fixed mesh**

At nx=80 fixed, vary n_ordinates ∈ {2, 4, 8, 16, 32}. Assert
monotone decrease + saturation to spatial floor.

### Gate Set 4 — Eigenvalue cross-check (BLOCKED ON F_N PREREQ)

**Gate 4.1 — Homogeneous-reflective k_∞ closed-form recovery**

The 2 multi-group homogeneous-reflective curvilinear cases
(`sphere_2g_homogeneous_dd_n20`, `sphere_2g_p1_aniso_dd_n20`)
recover `k_∞ = ν Σ_f / Σ_a` (computed via
`kinf_and_spectrum_homogeneous` on the same Mixture) to ≤ 5e-4
absolute. The 1G cylinder cases are **excluded** from absolute-
value gating per Cardinal Rule 6 (1G eigenvalue is flux-shape
independent — degenerate). Decorators:
`@pytest.mark.l1`, `@pytest.mark.verifies("sn-curvilinear-homogeneous-kinf-recovery")` (NEW).

**Gate 4.2 — Sood/F_N cross-check (BLOCKED)**

**Status**: deferred until F_N cylinder + 2G F_N infrastructure
ships (per user decision this session). When unblocked, the gate
spec is:

For each curvilinear snapshot with a Sood analogue, k_eff matches
the F_N-method semi-analytical computation on the closest Sood
La13511 case to ≤ 5e-4. The 2 heterogeneous 3-region cases
(`sphere_2g_3reg`, `cyl_2g_3reg`) remain without a direct F_N
analogue even after the prereq work; resumption decision deferred
between (a) replace with single-region Sood-grounded analogues, or
(b) regression-contract-only with WDD-aligned sweep as
structurally-independent reference.

### Gate Set 5 — V&V audit + ERR-026 closure

**5.1**: Remove the 4 ERR-026 xfail-strict tripwires (locations in
§4.5). Marker MUST be **removed** entirely (not relaxed to
`strict=False`).

**5.2**: `error_catalog.md` ERR-026 status PARTIAL → CLOSED (or
"spatial-closure aligned; angular default flip deferred to Phase D"
per the empirical-probe outcome).

**5.3**: `tests._harness.audit` clean: new labels
(`phase-c-apply-sweep-equivalence`,
`sn-mms-spherical-aniso-spatial-convergence`,
`sn-mms-cylindrical-aniso-spatial-convergence`,
`sn-curvilinear-homogeneous-kinf-recovery`) gain ≥ 1 incoming
`tests` edge. Phase B's `pole-mm-recurrence` orphan gains a
`tests` edge transitively from Gate 1.2 (since the canonical M-M
form becomes tested under WDD-aligned apply).

### Sphinx labels created by Phase C

1. `phase-c-apply-sweep-equivalence` — WDD diamond identity shared
   between apply's interior face-flux closure and sweep's WDD
   recurrence.
2. `sn-mms-spherical-aniso-spatial-convergence`
3. `sn-mms-cylindrical-aniso-spatial-convergence`
4. `sn-curvilinear-homogeneous-kinf-recovery`
5. `sn-curvilinear-3region-fn-crosscheck` (conditional on Gate 4.2
   resolution)

All 5 to be added to `docs/theory/discrete_ordinates.rst` in the
Phase C narrative subsection by the archivist.

---

## §6 — Open questions (resolved + remaining)

### Resolved by user (this session)

- **Q1: Gate 4.2 fixture path** — RESOLVED: build F_N cylinder +
  2G F_N infrastructure first; Phase C deferred until prereqs ship.
- **Q2: Default-flip sequencing** — RESOLVED: implement WDD spatial
  closure first; probe Gate 1.1 empirically; flip only if Gate 1.1
  passes.

### Resolved by test-architect (architectural review)

- **Q3: 1G permission for MMS Gates 3.1/3.2** — RESOLVED: 1G is
  acceptable for flux-shape claims (Cardinal Rule 6 prohibits 1G
  *eigenvalue* claims, not 1G fixed-source convergence). Citing
  `vv-principles` §1-group degeneracy.
- **Q4: Gate 1.2 face-flux extraction probe stability** — RESOLVED:
  well-defined for the WDD diamond at i=0 since A[0]=0 makes the
  pole face well-conditioned.

### Remaining (for Phase C resumption)

- **Q5: MMS trial function final form**. The originally-proposed
  `A(r)=sin(πr/R)`, `B(r)=r·(R−r)/R²` vanishes at r=0
  (test-architect caught this). Alternative `A(r)=cos(πr/(2R))`,
  `B(r)=r/R` is non-vanishing at the pole. **Recommended**: dispatch
  literature-researcher on Lewis & Miller §6.4 / Adams–Larsen §III
  to pick the canonical curvilinear MMS trial set before
  implementing Gates 3.1-3.2.

- **Q6: Heterogeneous 3-region snapshot strategy** —
  `sphere_2g_3reg` and `cyl_2g_3reg` have NO F_N analogue even
  post-prereq. Decision deferred between (a) replace with
  Sood-grounded analogues (likely needs Phase C resumption to
  rescope), or (b) regression-contract-only.

- **Q7: Empirical Gate 1.1 outcome** — unknowable without
  implementation; gates the default-flip decision.

---

## §7 — Resumption checklist (for the agent / session that picks
this up)

When the F_N prereqs ship:

1. **Re-read this file.**
2. **Re-read** `.claude/agent-memory/method-implementer/issue_168_phase_b_closeout.md`
   for Phase B's empirical context.
3. **Verify the assumed apply matvec call sites are unchanged.**
   Re-grep
   `0.5 \* \(fi\[:, n, i, 0\] \+ fi\[:, n, i [-+] 1, 0\]\)` in
   `orpheus/sn/operator.py`.
4. **Confirm `apply_transpose` is still built by dense probing.**
   Re-read [_ensure_dense_matrix](orpheus/sn/operator.py#L1169).
5. **Confirm the 4 ERR-026 xfail-strict markers are still in place.**
   Re-grep `xfail.*strict=True.*ERR-026` in `tests/sn/`.
6. **Confirm prereq F_N coverage**: at least F_N cylinder solver and
   2G F_N for sphere should exist. If still partial, reassess Phase C
   coverage matrix.
7. **Resolve Q5** (MMS trial function) via literature-researcher
   dispatch before Gate 3 implementation.
8. **Implement Gate 1.4 → 1.3 → 1.1 → 1.2** (foundation/L0 first).
9. **Empirical default-flip decision** based on Gate 1.1 outcome.
10. **Implement Gate 3 (MMS) → Gate 4.1 (k_∞) → Gate 4.2 (F_N/Sood)**.
11. **Resolve Q6** (heterogeneous 3-region strategy).
12. **Snapshot regen + marker removal + ERR-026 closure**.

The rewrite path itself is small (the change is concentrated at 5
sites per geometry × 2 geometries = ~10 sites + the per-ordinate
sweep restructuring of the matvec loop). The verification scaffolding
is the larger deliverable.
