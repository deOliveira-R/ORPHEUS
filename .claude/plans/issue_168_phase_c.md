# Issue #168 Phase C — Sweep-frame apply matvec rewrite (RESUMED 2026-05-12)

**Status (2026-05-12)**: **RESUMED with revised architecture.** Two
post-pause developments dissolve the original blocker and reframe the
scope:

1. The trajectory_resolvent (Peierls Variant α Green's-function)
   campaign (commits `37e3e29` + `cf662a6` + `604f380` + `e10c33c`,
   2026-05-12) shipped cylinder MR Phase 1b, **explicitly built to
   close the cylinder-2G ERR-026 gap**. trajectory_resolvent now
   covers 5 of the 6 deleted curvilinear regression snapshots at
   machine-precision-class precision; the 6th (P1 anisotropic) routes
   to k_∞ closed-form which is shape-independent.

2. The cross-domain-attacker analysis on 2026-05-12 identified a
   **structural inconsistency between the apply matvec's boundary
   closure and the §16A.3 affine BC contract** (`γ_-ψ = R·G·γ_+ψ + q`
   per [eq affine-bc-form](docs/theory/boundary_conditions.rst#L143)).
   The user — flagging this via their "ghost cells for higher-order
   boundary closure" intuition — confirmed the architectural
   reformulation: the apply matvec is rewritten as **one sweep
   iteration**, with the BC trace law owning the boundary edge and the
   `InflowTraceSpace` vector playing the role of the user's "ghost
   cell" (value defined by the realised BC, not extrapolated from
   interior cell-centres).

Phase A and Phase B shipped on `refactor/sn-operator-algebra` (commits
`d73ef68`, `5df90df`, `ef377cd`). Phase A's `BoundaryFaceFlux`
Protocol is now identified as a patch on top of the wrong architecture
and **retires** as part of Phase C.

The F_N method route (issues #170 / #171) is no longer a Phase C
prerequisite; both issues remain open as future work for an
independent F_N-pillar V&V coverage (different mathematical framework
— singular-eigenfunction Fredholm — useful for cross-pillar agreement
but not gating).

---

## §1 — Why Phase C resumes now

### Unblocker 1: trajectory_resolvent coverage of the 6 snapshots

The trajectory_resolvent module ships `solve_greens_function_*` and
`solve_greens_function_*_mr` entry points covering sphere + cylinder
+ slab in homogeneous + multi-region forms. Coverage matrix for the
6 deleted curvilinear regression snapshots:

| # | Snapshot | trajectory_resolvent entry point | Precision (measured) | Coverage |
|---|---|---|---|---|
| 1 | `sphere_2g_homogeneous_dd_n20` | `solve_greens_function_sphere_mg(R=2.0, sigma_t, sigma_s, nu_sigma_f, chi, alpha=1.0)` | k_∞ exact via V_α1 identity (rtol≤1e-10) | ✓ FULL |
| 2 | `sphere_2g_3reg_dd_n40` | `solve_greens_function_sphere_mr(radii=(0.5,1.5,2.0), sigma_t=(3,2), sigma_s=(3,2,2), nu_sigma_f=(3,2), chi, alpha=1.0)` | MR↔MG homogeneous-reduction rtol=1e-9 (2G) | ✓ FULL |
| 3 | `sphere_2g_p1_aniso_dd_n20` | **trajectory_resolvent is P0-only** | — | Routes to Gate 4.1 (k_∞ closed-form, shape-independent of anisotropy) |
| 4 | `cyl_1g_homogeneous_LS4_dd_n20` | `solve_greens_function_cylinder(R=2.0, sigma_t, sigma_s, nu_sigma_f, alpha=1.0)` | k_∞ exact via V_α1_cyl identity; Sood Ua-1-O-CY vacuum at 8.5e-6 | ✓ FULL |
| 5 | `cyl_1g_homogeneous_product_dd_n20` | Same as #4 (different SN quadrature, same continuous reference) | Same | ✓ FULL |
| 6 | `cyl_2g_3reg_LS4_dd_n40` | `solve_greens_function_cylinder_mr(radii=(0.5,1.5,2.0), sigma_t=(3,2), sigma_s=(3,2,2), nu_sigma_f=(3,2), chi, alpha=1.0)` | MR↔MG K=3 2G rtol=1e-9; k_∞ rel_err<1e-6 at α=1 | ✓ FULL |

5/6 covered at machine-precision-class precision via trajectory_resolvent;
snapshot #3 (P1 anisotropic) covered by Gate 4.1 (k_∞ analytical).

**Verification wiring** (per user directive — no shortcuts): the
Gate 4.2 test harness calls the **bare** `solve_greens_function_*_mr`
entry points directly for snapshots 2 and 6. The
`Billiard._infer_geometry_kind` facade does not currently route to MR
variants — that micro-capability gap is tracked at
[#190](https://github.com/deOliveira-R/ORPHEUS/issues/190) and is NOT
a Phase C work item; #190 will be implemented properly in a dedicated
session. Phase C's verification harness consumes the bare functions
explicitly without working around the facade.

### Unblocker 2: SN BC trace-law architecture

Issues #186 / #176 / #188 + Waves 0-12 installed the unified BC
architecture. The relevant pieces for Phase C:

- **`InflowTraceSpace` / `OutflowTraceSpace`** ([trace_space.py](orpheus/numerics/trace_space.py))
  — typed function spaces with per-face directional masks covering
  every Mesh1D coord system (Cartesian / spherical / cylindrical)
  post-Issue #188.
- **`BoundaryTraceLaw`** ([orpheus/geometry/boundary/](orpheus/geometry/boundary/))
  — pure descriptor (no `apply` method post Issue #186/B3). Carries
  `geometry_map`, `response_kernel`, `source` per affine form
  `γ_-ψ = R·G·γ_+ψ + q`.
- **`SNBoundaryRealizer.realize(law, method_space)`** ([orpheus/sn/boundary_realizer.py:122](orpheus/sn/boundary_realizer.py#L122))
  — produces a 1-arg `LinearOperator`. `apply(psi_face_outflow) → psi_face_inflow`.

This is the substrate Phase C consumes. The realised BC operator
applied at the boundary edge gives the proper inflow trace value
without any algebraic extrapolation — that's the user's "ghost cell"
idiom realised as a typed `InflowTraceSpace` vector.

### Unblocker 3: the architectural inconsistency, confirmed

[operator.py:533](orpheus/sn/operator.py#L533):
```python
outgoing = fi[:, :, -1, 0].T   # CELL-CENTRES — violates §16A.3 contract
incoming = bc_outer.apply(outgoing)
```
[boundary_face_flux.py:336](orpheus/sn/spatial/boundary_face_flux.py#L336):
```python
psi_face_out = 1.5 * fi[:, n, N-1, 0] - 0.5 * fi[:, n, N-2, 0]
# Algebraic extrapolation; ignores BC entirely
```

**Magnitude**: silent for specular(α=1) because the reflection
permutation commutes with cell-centre fills; surfaces for vacuum
(α=0), albedo (0<α<1), prescribed inflow, white BC — the regimes
where higher-order behavior should appear.

---

## §2 — Revised architecture: sweep-frame matvec

The apply matvec becomes **one sweep iteration semantically**. The
critical architectural principles:

1. **Vectorise across the ordinate dimension.** Per the user's
   directive (and `[orpheus/sn/angular_operator.py:183-186](orpheus/sn/angular_operator.py#L183)`
   precedent), use `outgoing_mask = quad.mu_x > +eps` and
   `incoming_mask = quad.mu_x < -eps` to operate on whole ordinate
   subsets simultaneously. **No new `for n in range(quad.N)` loops in
   the Phase C matvec** — the systemic cleanup of the 14 existing
   sites is tracked separately at
   [#191](https://github.com/deOliveira-R/ORPHEUS/issues/191).

2. **Use the cell-visit graph by direction, not by ordinate.** The
   cell-visit DAG already encodes direction (the sweep planner's
   per-quadrant `_diag_cache` is keyed by `(sign(mx), sign(my))`, not
   by ordinate). Phase C adds a small API method
   `SNMesh.iter_cells_by_direction(direction_sign[, mu_level_idx])`
   that surfaces what the graph already knows. The existing
   `iter_cell_visits(ordinate_idx=n)` stays for callers that need
   per-ordinate streaming terms; the new method is the direction-only
   variant the vectorised matvec consumes.

3. **BC trace law owns the boundary edge.** The realised BC operator
   `bc_outer.apply(ψ_face_outflow)` is called ONCE per direction at
   the boundary face, consuming the WDD-propagated outflow face flux
   (not cell-centre values) and producing the inflow face flux for
   the partner direction. This is the architectural fix to the §16A.3
   violation.

4. **WDD diamond uniform interior + boundary.** The closure
   `ψ_face_out = 2·ψ_cell − ψ_face_in` applies per cell. The
   boundary face value is determined by either (a) the BC trace law
   for inflow, or (b) the WDD propagation chain for outflow — never
   by algebraic extrapolation through interior cell-centres.

### The matvec body (spherical; cylindrical is structurally identical with per-level redistribution)

```python
def transport_operator_matvec_spherical(
    solution, eq_map, quad, sig_t, nx, ng,
    face_areas, volumes, alpha_half, redist_dAw, tau_mm,
    *,
    bc_outer,
    pole_angular_closure=None,
):
    from .spatial.pole_angular_closure import LegacyTauSymmetricInterpolation
    if pole_angular_closure is None:
        pole_angular_closure = LegacyTauSymmetricInterpolation()

    fi = solution_to_angular_flux_spherical(solution, eq_map, quad, nx, ng)
    A = face_areas
    V = volumes[:, 0]
    eps = 1e-15

    outgoing_mask = quad.mu_x > +eps
    incoming_mask = quad.mu_x < -eps
    mu_out = quad.mu_x[outgoing_mask]
    mu_in  = quad.mu_x[incoming_mask]
    n_out = mu_out.size
    n_in  = mu_in.size

    redist_full = pole_angular_closure(
        fi[..., 0], alpha_half, redist_dAw, tau_mm, V,
        level_indices=None,
    )
    # redist_full shape (ng, N, nx). We slice by mask in each phase.

    # Pre-allocate output, plus the boundary outflow buffer.
    lhs = np.empty((ng, eq_map.n_eq))
    outflow_at_boundary = np.zeros((ng, quad.N))   # full angular vector

    # ── Phase 1: outgoing ordinates (μ > 0), sweep i = 0 → nx-1 ─────
    # Pole face: ψ_face = 0 by symmetry (also multiplied by A[0]=0).
    psi_face_in = np.zeros((ng, n_out))
    for visit in sn_mesh.iter_cells_by_direction(+1):
        i = visit.cell_idx
        psi_cell = fi[:, outgoing_mask, i, 0]                  # (ng, n_out)
        psi_face_out = 2.0 * psi_cell - psi_face_in
        streaming = (
            mu_out[None, :]
            * (A[i + 1] * psi_face_out - A[i] * psi_face_in)
            / V[i]
        )
        redistribution = redist_full[:, outgoing_mask, i]
        collision = sig_t[i, 0, :, None] * psi_cell
        lhs[:, eq_map.unknowns_at_cell_for_mask(i, outgoing_mask)] = (
            streaming + redistribution + collision
        )
        psi_face_in = psi_face_out
    outflow_at_boundary[:, outgoing_mask] = psi_face_out

    # ── BC trace law at the boundary ──────────────────────────────
    # bc_outer.apply consumes the FACE flux (γ_+ψ), per §16A.3.
    # Returns the inflow trace for the same face's incoming ordinates.
    inflow_full = bc_outer.apply(outflow_at_boundary)
    inflow_at_boundary = inflow_full[:, incoming_mask]   # (ng, n_in)

    # ── Phase 2: incoming ordinates (μ < 0), sweep i = nx-1 → 0 ───
    psi_face_in = inflow_at_boundary
    for visit in sn_mesh.iter_cells_by_direction(-1):
        i = visit.cell_idx
        psi_cell = fi[:, incoming_mask, i, 0]
        psi_face_out = 2.0 * psi_cell - psi_face_in
        # For μ<0, ψ_face_in is at A[i+1] (right face = inflow) and
        # ψ_face_out is at A[i] (left face = outflow downstream).
        streaming = (
            mu_in[None, :]
            * (A[i + 1] * psi_face_in - A[i] * psi_face_out)
            / V[i]
        )
        redistribution = redist_full[:, incoming_mask, i]
        collision = sig_t[i, 0, :, None] * psi_cell
        lhs[:, eq_map.unknowns_at_cell_for_mask(i, incoming_mask)] = (
            streaming + redistribution + collision
        )
        psi_face_in = psi_face_out

    return lhs.ravel(order='F')
```

Notes on the pseudocode:

- `eq_map.unknowns_at_cell_for_mask(i, mask)` is a new helper that
  returns the unknown indices `k` for the unknowns at cell `i` whose
  ordinates satisfy `mask`. Built from `eq_map.ordinate` and
  `eq_map.ix` as a precomputed lookup. Replaces the per-equation
  `eq_map.unknown_index_for(n, i)` scatter pattern.
- The `iter_cells_by_direction(±1)` API yields `CellVisit`-like packets
  carrying `cell_idx` only. The packet does NOT need `streaming_terms`
  or `face_area_downstream` because the matvec computes those itself
  from intrinsic cell geometry (`A`, `V`) and ordinate metadata
  (`mu_out`, `mu_in`).
- `bc_outer.apply(outflow_at_boundary)` operates on the full `(ng, N)`
  angular vector at the boundary face. The realised BC operator
  internally consumes the outflow mask (its `OutflowTraceSpace`
  metadata) and writes only the inflow slots — the outflow slots in
  the output are unspecified (likely zero or unchanged). The matvec
  reads back only `inflow_full[:, incoming_mask]`.
- The `redist_full` precomputation handles the pole angular closure
  uniformly. This stays unchanged from Phase B.
- The cylindrical case is structurally identical with the modification
  that `iter_cells_by_direction` takes a `mu_level_idx` second
  argument (per-level cell ordering for cylindrical's azimuthal-DAG
  topology), and `redist_full` is computed per-level via
  `quad.level_indices`.

### What closes by construction

| Defect | Pre-Phase-C status | How Phase C closes it |
|---|---|---|
| 1: outer-face truncation | PARTIAL (Phase A `DDExtrapolation`) | WDD propagation gives outflow face flux at full 2nd order |
| 2: BC-fill cell-centre contamination | PARTIAL (`bf_flux[:, n, 1]` storage split) | BC operator consumes face values directly; no separate face storage needed |
| 3: Bailey ΔA/w pole closure | PARTIAL (Phase B 3-strategy PoleAngularClosure with Legacy as default) | Canonical M-M angular closure becomes consistent (pending Gate 1.1 empirical confirmation) because spatial closure is now WDD throughout |

### What stays

- **`PoleAngularClosure` Protocol** — sphere center is intrinsic
  geometry (coordinate-system singularity), not external BC.
  Three-strategy split (Legacy / BFF / MMS) stays. Phase C's empirical
  Gate 1.1 probe decides the default.
- **`SNStreamingOperator.apply_transpose`** via dense-probe
  construction. Linearity of the rewritten apply guarantees the
  transpose is correctly tracked.
- **`SNMethodSpace` / `SNBoundaryRealizer` / `_BoundBoundaryOperator`**
  — the Wave 8/9/10/11/12 + Issue #186/#176/#188 substrate.

### What retires

- **`BoundaryFaceFlux` Protocol** retires:
  - Delete `orpheus/sn/spatial/boundary_face_flux.py`
  - Delete `tests/sn/spatial/test_boundary_face_flux.py` (21 foundation tests)
  - Remove `SNMesh.boundary_face_flux` constructor field from
    [orpheus/sn/geometry.py](orpheus/sn/geometry.py)
  - Remove `boundary_face_flux_closure` kwarg from
    `transport_operator_matvec_*` and `SNStreamingOperator.apply`
    plumbing
  - Update Sphinx narrative in `docs/theory/discrete_ordinates.rst`
    to document the retirement and the sweep-frame architectural fix

### What gets added

- **`SNMesh.iter_cells_by_direction(direction_sign[, mu_level_idx])`**:
  - For Cartesian / spherical 1D: yields cells in DAG order for the
    given direction sign (`+1` for outward, `-1` for inward). No
    ordinate parameter — the cell-visit graph already encodes
    direction.
  - For cylindrical: requires `mu_level_idx` (cylindrical's
    azimuthal-DAG topology is per-level). Yields cells in DAG order
    for (direction_sign, mu_level_idx).
  - Yields lightweight `CellVisit` packets carrying `cell_idx`. The
    packet does NOT carry per-ordinate `streaming_terms` — the
    vectorised matvec computes streaming terms itself from `(mu_array,
    A, V)`.
  - Foundation test: confirms the yielded `cell_idx` sequence matches
    `iter_cell_visits(ordinate_idx=any_representative_n_for_direction)`
    for every representative ordinate. Pinned bit-identical via
    `np.array_equal`.

- **`EquationMap.unknowns_at_cell_for_mask(i, ordinate_mask) -> np.ndarray`**:
  - Helper that returns the indices `k` for unknowns at cell `i` whose
    ordinates satisfy the mask. Built from `eq_map.ordinate[]` and
    `eq_map.ix[]` as a precomputed lookup at `eq_map` construction
    time. O(1) lookup at matvec call time.

---

## §3 — Implementation sequencing

1. **Add `SNMesh.iter_cells_by_direction` + `EquationMap.unknowns_at_cell_for_mask`**
   with foundation tests pinning equivalence to existing iteration
   patterns.

2. **Implement Gate 1.4 (apply linearity)** FIRST as a pure software
   contract — without this, all downstream gates are invalid.

3. **Implement the sweep-frame matvec rewrite**:
   - `transport_operator_matvec_spherical` per §2 pseudocode
   - `transport_operator_matvec_cylindrical` analogously (per-level
     `redist_full` + `mu_level_idx` to `iter_cells_by_direction`)
   - Simplify `solution_to_angular_flux_spherical`: returns only `fi`
     (no `boundary_face_flux` companion array; BC apply moves into
     matvec body)
   - Update `SNStreamingOperator.apply` dispatcher to drop the
     `boundary_face_flux_closure` plumbing

4. **Retire `BoundaryFaceFlux` Protocol**:
   - Delete `orpheus/sn/spatial/boundary_face_flux.py`
   - Delete `tests/sn/spatial/test_boundary_face_flux.py`
   - Remove `SNMesh.boundary_face_flux` constructor field
   - Update Sphinx narrative

5. **Implement Gate 1.3, 1.1, 1.2, 1.5** (foundation + L0). These are
   tight-loop iteration: Gate 1.1's per-ordinate flat-flux residual
   probe is the canonical curvilinear bug-class diagnostic.

6. **Empirical default-flip decision (per user 2026-05-10 sequencing)**:
   - Run Gate 1.1 with `pole_angular_closure=MorelMontryAngularSweep`.
   - If pass: flip `SNMesh.pole_angular_closure` default
     `LegacyTauSymmetricInterpolation` → `MorelMontryAngularSweep`;
     flip curvilinear `solve_sn_fixed_source` default
     `"source_iteration"` → `"krylov"`.
   - If fail: keep Legacy as default; file follow-up issue with
     empirical evidence; Phase C ships as "sweep-frame matvec aligned;
     angular default flip deferred to Phase D".

7. **Implement Gate 2** preservation tests:
   - 5 Cartesian regression snapshots pass bit-identical (rtol=1e-12).
   - Phase B 28 foundation tests + 5 L1 flat-flux-identity tests pass.

8. **Implement Gate 3 MMS**:
   - Spherical + cylindrical spatial convergence with the **angularly-
     non-trivial ansatz** ψ_chosen = (A(r) + B(r)·μ)/W with A(r) =
     cos(πr/(2R)) (=1 at pole) and B(r) = r/R (zero at pole, unit at
     outer). Linear-μ MANDATORY (Mode 7 / ERR-026 hide-the-bug pattern
     explicitly avoided).
   - Angular convergence at fixed nx=80, varying n_ordinates.

9. **Implement Gate 4 absolute-value cross-check**:
   - Gate 4.1: k_∞ closed-form for snapshots 1 + 3 (multi-group
     homogeneous-reflective, P0 and P1 anisotropic).
   - Gate 4.2: trajectory_resolvent cross-check for snapshots 1, 2,
     4, 5, 6 via bare `solve_greens_function_*` and
     `solve_greens_function_*_mr` entry points per the §1 coverage
     matrix. The verification harness calls bare functions
     explicitly; #190 facade routing is separate work.

10. **Regenerate the 6 deleted curvilinear snapshots** via
    `python -m tests.sn.regression._generate_snapshots`.

11. **Remove the 4 ERR-026 xfail-strict markers** at §4.6 locations.
    Markers must be **removed** (not relaxed to `strict=False`).

12. **Documentation** (handoff to archivist):
    - Sphinx narrative subsection in
      `docs/theory/discrete_ordinates.rst` describing sweep-frame
      matvec architecture with the trace-law BC integration and the
      `iter_cells_by_direction` API
    - Update `docs/theory/boundary_conditions.rst` §16A.3 to note the
      SN apply matvec now honors the trace contract
    - `error_catalog.md` ERR-026 closure narrative (PARTIAL → CLOSED
      or "spatial-closure aligned; angular default flip deferred to
      Phase D" per Gate 1.1 outcome)
    - 6 new `:label:` entries (see §5 below)

13. **Update `sn_reshape.md` Wave H row** + close GH Issue #168.

14. **Submit closeout memo at
    `.claude/agent-memory/method-implementer/issue_168_phase_c_closeout.md`**.

---

## §4 — Architectural map (current branch state, 2026-05-12)

### §4.1 Apply matvec — current face-flux sites (to be replaced)

**Spherical** ([transport_operator_matvec_spherical](orpheus/sn/operator.py#L557)):

| Site | Line | Phase C action |
|---|---|---|
| Interior right face `0.5*(fi[i]+fi[i+1])` | [697](orpheus/sn/operator.py#L697) | REPLACE: WDD `psi_face_out = 2·ψ_cell − ψ_face_in` per cell-visit |
| Outer right μ>0 (BoundaryFaceFlux) | [701-703](orpheus/sn/operator.py#L701) | RETIRE: WDD propagation gives outflow at full order |
| Outer right μ<0 `bf_flux[:, n, 1]` | [709](orpheus/sn/operator.py#L709) | REPLACE: read from BC trace law applied to WDD-propagated outflow at boundary |
| Interior left face `0.5*(fi[i-1]+fi[i])` | [712](orpheus/sn/operator.py#L712) | REPLACE: WDD asymmetric form |
| Pole left face (i=0) `psi_left = 0.0` | [714](orpheus/sn/operator.py#L714) | KEEP (intrinsic; multiplied by A[0]=0) |

**Cylindrical** ([transport_operator_matvec_cylindrical](orpheus/sn/operator.py#L740)):

| Site | Line | Phase C action |
|---|---|---|
| Interior right | [840](orpheus/sn/operator.py#L840) | Same as spherical |
| Outer right μ>0 | [844-846](orpheus/sn/operator.py#L844) | Retire |
| Outer right μ<0 | [850](orpheus/sn/operator.py#L850) | Replace with BC trace law |
| Interior left | [853](orpheus/sn/operator.py#L853) | WDD asymmetric |
| Pole left | [855](orpheus/sn/operator.py#L855) | Keep |

### §4.2 BC apply current site (to be relocated)

[solution_to_angular_flux_spherical:530-552](orpheus/sn/operator.py#L530):
the cell-centre-based BC apply moves out of this function. Function
returns only `fi` (cell-centres). BC apply happens inside the matvec
body, operating on the WDD-propagated outflow face values.

### §4.3 `apply_transpose` mechanism (unchanged)

[_ensure_dense_matrix](orpheus/sn/operator.py#L1172) at lines
1172-1201 builds the transpose by dense-probing `apply`. Phase C
preserves `apply` linearity (Gate 1.4) → transpose tracks
automatically.

### §4.4 Equation map ordering (unchanged; one helper added)

[build_equation_map_spherical](orpheus/sn/operator.py#L389) iterates
`(spatial outer, ordinate inner)`. The sweep-frame matvec scatters
results back via the new `eq_map.unknowns_at_cell_for_mask(i, mask)`
helper, preserving the packed-vector contract.

### §4.5 `SNMesh.iter_cell_visits` (unchanged; sibling method added)

[orpheus/sn/geometry.py:433-540](orpheus/sn/geometry.py#L433) — the
existing per-ordinate iterator stays. Phase C adds the sibling
`iter_cells_by_direction(direction_sign[, mu_level_idx])` that exposes
what the cell-visit graph already encodes via direction sign.

### §4.6 ERR-026 xfail-strict tripwires (4 tests, unchanged)

| # | File | Line | Test |
|---|---|---|---|
| 1 | [tests/sn/test_mms_curvilinear.py](tests/sn/test_mms_curvilinear.py) | 42-53 | `test_sn_spherical_mms_converges_second_order` |
| 2 | [tests/sn/test_mms_curvilinear.py](tests/sn/test_mms_curvilinear.py) | 94-101 | `test_sn_cylindrical_mms_converges_second_order` |
| 3 | [tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py](tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py) | 70-90 | `test_sn_spherical_aniso_mms_converges_second_order` |
| 4 | [tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py](tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py) | 135-143 | `test_sn_cylindrical_aniso_mms_converges_second_order` |

All `@pytest.mark.xfail(strict=True)` → xpass forces marker removal as
Phase C closes the underlying bugs.

### §4.7 Regression-contract mechanism (unchanged)

[tests/sn/regression/test_dd_regression.py:46-70](tests/sn/regression/test_dd_regression.py#L46)
— missing snapshots trigger `pytest.skip`. Tolerance: `rtol=1e-12,
atol=1e-13`. Snapshot regen via
`python -m tests.sn.regression._generate_snapshots`.

---

## §5 — Verification plan (13 gates)

### Gate Set 1 — Sweep-frame matvec correctness

**Gate 1.1 (L0) — `test_apply_curvilinear_per_ordinate_flat_flux_residual`**

For ψ constant in space and per-ordinate, on a reflective-BC
homogeneous sphere/cylinder:
- Σ_t=0 → bit-zero output (atol=1e-13)
- Σ_t≠0 uniform → `Σ_t · ψ` per ordinate (rtol=1e-12)

Parametrisation: spherical + cylindrical; n_ordinates ∈ {4 GL, 8 GL,
LS-4, Product 2x4}; n_groups ∈ {2, 4}; nx ∈ {4, 10, 20}; Σ_t ∈ {0,
0.5}; pole_angular_closure ∈ {Legacy, BFF, MMS}.

**Empirical decision point**: pole_angular_closure=MMS pass/fail
determines default flip.

Decorators: `@pytest.mark.l0`,
`@pytest.mark.verifies("dd-curvilinear-scalar")`,
`@pytest.mark.catches("ERR-026")`.

**Gate 1.2 (foundation) — `test_apply_face_fluxes_match_sweep_recurrence`**

Per-ordinate face-flux propagation identity: implicit face fluxes
in apply (extracted by inverting cell-balance) match sweep's WDD
recurrence on `q = apply(ψ_cells)`. **Structural by construction**
post-Phase-C because both paths consume `iter_cells_by_direction` +
same WDD diamond + same BC trace law.

Tolerance: `np.array_equal` (bit-identical).

Decorators: `@pytest.mark.foundation`,
`@pytest.mark.verifies("phase-c-apply-sweep-equivalence")` (NEW),
`@pytest.mark.catches("ERR-026")`.

**Gate 1.3 (foundation) — `test_apply_apply_transpose_reciprocity_under_sweep_frame`**

`<L·ψ, φ> = <ψ, Lᵀ·φ>` to rtol=1e-12, atol=1e-13 on random ψ, φ.
Free if Gate 1.4 (linearity) passes.

**Gate 1.4 (foundation) — `test_apply_linearity_under_sweep_frame`**

`apply(α·ψ + β·φ) = α·apply(ψ) + β·apply(φ)` to rtol=1e-13.
**Precondition** for Gates 1.2 + 1.3.

**Gate 1.5 (foundation, NEW) — `test_bc_trace_contract_respected_by_matvec`**

For each `BoundaryTraceLaw` concrete (Vacuum, Reflective, White,
Albedo, PrescribedInflow): the apply matvec's BC integration consumes
the WDD-propagated outflow face value (not cell-centres) and produces
the inflow face value consistent with `law.realize().apply(outflow)`.

The architectural-fix anchor. Pins the §16A.3 contract.

Parametrisation: per BC kind × geometry × ordinate count × 5 random
ψ_cells inputs. Compare:
- `inflow_via_matvec = capture_internal_bc_apply_input_output(apply(ψ_cells))`
- `inflow_independent = bc.realize().apply(outflow_at_boundary_for(ψ_cells))`
- Assert bit-identical.

Decorators: `@pytest.mark.foundation`,
`@pytest.mark.verifies("bc-trace-contract-respected-by-matvec")` (NEW).

### Gate Set 2 — Bit-identity preservation

- **2.1**: 5 Cartesian regression snapshots pass bit-identical
  (rtol=1e-12). Cartesian path uses upwind FD; Phase C is curvilinear-
  only.
- **2.2**: Phase B 28 foundation tests pass (PoleAngularClosure suite).
- **2.3**: Phase B 5 L1 flat-flux-identity tests pass.
- **2.4**: 21 Phase A `boundary_face_flux` tests **DELETED** with the
  Protocol retirement. Document the deletion in commit message and
  Sphinx narrative.

### Gate Set 3 — MMS convergence (L1)

**Gate 3.1 — Spatial MMS convergence, spherical, linear-μ ansatz**

ψ_chosen = (A(r) + B(r)·μ)/W with A(r) = cos(πr/(2R)) (=1 at pole),
B(r) = r/R (zero at pole, unit at outer). Linear-μ MANDATORY.

Refine nx ∈ {10, 20, 40, 80, 160} at fixed n_ordinates=8 GL. Assert
`min(orders[-3:]) > 1.9`.

Decorators: `@pytest.mark.l1`, `@pytest.mark.slow`,
`@pytest.mark.verifies("sn-mms-spherical-aniso-spatial-convergence")` (NEW),
`@pytest.mark.catches("ERR-026")`.

**Gate 3.2 — Cylindrical MMS, linear-μ_x ansatz, 2 quadrature families**

Cylindrical analogue. LS-4 and Product 2x4 parametrised (flag
quadrature-dependent constants — Signature 4 / ERR-004).

Decorators: `@pytest.mark.l1`, `@pytest.mark.slow`,
`@pytest.mark.verifies("sn-mms-cylindrical-aniso-spatial-convergence")` (NEW),
`@pytest.mark.catches("ERR-026")`.

**Gate 3.3 — Angular convergence at fixed mesh**

nx=80, n_ordinates ∈ {2, 4, 8, 16, 32}. Assert monotone decrease +
saturation to spatial floor.

### Gate Set 4 — Absolute-value cross-check

**Gate 4.1 — Homogeneous-reflective k_∞ recovery**

Snapshots 1 (`sphere_2g_homogeneous`) and 3 (`sphere_2g_p1_aniso`)
recover `k_∞ = ν Σ_f / Σ_a` to ≤ 5e-4 absolute (target ≤ 1e-6 for P0;
P1 may have larger spatial-discretisation correction).

Decorators: `@pytest.mark.l1`,
`@pytest.mark.verifies("sn-curvilinear-homogeneous-kinf-recovery")` (NEW).

**Gate 4.2 — trajectory_resolvent cross-check (5 P0 snapshots)**

Verification harness calls **bare** `solve_greens_function_*` and
`solve_greens_function_*_mr` entry points (NOT the `Billiard` facade
— #190 tracks facade MR routing as separate work).

| Snapshot | Variant α call | Tolerance |
|---|---|---|
| 1 `sphere_2g_homogeneous` | `solve_greens_function_sphere_mg(R=2.0, sigma_t, sigma_s, nu_sigma_f, chi, alpha=1.0, n_r=24, n_mu=20, n_traj_quad=96)` | rtol ≤ 1e-9 (V_α1 identity floor; relax to ≤ 5e-4 if SN nx-discretisation dominates) |
| 2 `sphere_2g_3reg` | `solve_greens_function_sphere_mr(radii=(R1, R2, R3), sigma_t=(3,2), sigma_s=(3,2,2), nu_sigma_f=(3,2), chi=(3,2), alpha=1.0, ...)` | rtol ≤ 1e-9 (MR↔MG inheritance) |
| 4 `cyl_1g_homogeneous_LS4` | `solve_greens_function_cylinder(R=2.0, sigma_t, sigma_s, nu_sigma_f, alpha=1.0, ...)` | rtol ≤ 1e-9 (V_α1_cyl identity) |
| 5 `cyl_1g_homogeneous_product` | Same as #4 (different SN quadrature, same continuous reference) | Same |
| 6 `cyl_2g_3reg` | `solve_greens_function_cylinder_mr(radii=(R1, R2, R3), sigma_t=(3,2), sigma_s=(3,2,2), nu_sigma_f=(3,2), chi=(3,2), alpha=1.0, ...)` | rtol ≤ 1e-6 (MR k_∞ recovery floor) |

Decorators: `@pytest.mark.l1`, `@pytest.mark.slow`,
`@pytest.mark.verifies("sn-curvilinear-trajectory-resolvent-crosscheck")` (NEW).

### Gate Set 5 — V&V audit + ERR-026 closure

- **5.1**: Remove 4 ERR-026 xfail-strict markers (§4.6 locations).
  Markers REMOVED (not relaxed).
- **5.2**: `error_catalog.md` ERR-026 PARTIAL → CLOSED. Closure
  narrative cites the trajectory_resolvent verification chain.
- **5.3**: `tests._harness.audit` clean. 6 NEW labels gain ≥1 `tests`
  edge. Phase B's `pole-mm-recurrence` orphan gains a `tests` edge
  transitively from Gate 1.5 (canonical M-M tested under sweep-frame
  apply when it becomes default).

### Sphinx labels created (6 new)

1. `phase-c-apply-sweep-equivalence` — structural identity of apply and
   sweep face-flux propagation under WDD + BC trace law
2. `bc-trace-contract-respected-by-matvec` — pins §16A.3 contract
   under sweep-frame matvec
3. `sn-mms-spherical-aniso-spatial-convergence`
4. `sn-mms-cylindrical-aniso-spatial-convergence`
5. `sn-curvilinear-homogeneous-kinf-recovery`
6. `sn-curvilinear-trajectory-resolvent-crosscheck`

All 6 added to `docs/theory/discrete_ordinates.rst` Phase C subsection.

---

## §6 — Open questions

### Resolved this session

- **Q1 — Gate 4.2 fixture path**: trajectory_resolvent unblocker (§1).
  Bare functions, not Billiard facade.
- **Q2 — Architectural inconsistency**: confirmed via cross-domain-
  attacker analysis. Sweep-frame reformulation accepted.
- **Q3 — MMS trial function**: A(r) = cos(πr/(2R)), B(r) = r/R.
- **Q4 — Snapshot strategy for 3-region cases**: trajectory_resolvent
  MR covers; bare-function call in verification harness.
- **Q5 — Ordinate iteration anti-pattern**: vectorise in Phase C; #191
  tracks systemic cleanup.
- **Q6 — `iter_cells_by_direction` vs `iter_cell_visits(ordinate_idx=n)`**:
  add the direction-keyed sibling method as part of Phase C; the graph
  already knows direction.

### Remaining for implementation

- **Q7 — Empirical Gate 1.1 outcome**: gates the default-flip decision.
  Unknowable until implementation. Sweep-frame architecture more
  likely to make MMS angular closure viable (because spatial closure
  is now WDD throughout, matching what MMS expects), but the empirical
  probe at Gate 1.1 is the gating step.

### Tracked as separate work (not blocking)

- **#190**: `Billiard._infer_geometry_kind` MR routing — facade
  cleanup. Phase C uses bare functions.
- **#191**: Systemic `for n in range(quad.N)` cleanup — 14 sites.
  Phase C does not propagate the anti-pattern but does not refactor
  the existing 14.
- **#170**: F_N cylinder solver — future work for independent F_N
  pillar coverage.
- **#171**: 2G F_N for sphere/cylinder — future work.

---

## §7 — Implementation checklist (for method-implementer dispatch)

The dispatch brief points at this file. The implementer:

1. **Reads** this entire file plus:
   - `.claude/agent-memory/method-implementer/issue_168_phase_a_closeout.md`
   - `.claude/agent-memory/method-implementer/issue_168_phase_b_closeout.md`
   - `.claude/agent-memory/cross-domain-attacker/issue_168_phase_c_sweep_frame.md`
     (NEW 2026-05-12)
   - `docs/theory/boundary_conditions.rst` §16A.3 + §bc-realizer-layer
   - `docs/theory/discrete_ordinates.rst` (Phase B sections to keep
     consistent with Phase C additions)

2. **Implements in order**:
   - [ ] `SNMesh.iter_cells_by_direction(direction_sign[, mu_level_idx])`
     + foundation test (equivalence to `iter_cell_visits(ordinate_idx=
     representative_n)`)
   - [ ] `EquationMap.unknowns_at_cell_for_mask(i, mask)` + foundation test
   - [ ] Gate 1.4 (linearity probe) FIRST as a precondition
   - [ ] `transport_operator_matvec_spherical` rewrite per §2 pseudocode
   - [ ] `transport_operator_matvec_cylindrical` rewrite (per-level)
   - [ ] `solution_to_angular_flux_spherical` simplification (BC apply removed)
   - [ ] `SNStreamingOperator.apply` dispatcher cleanup
   - [ ] Retire `BoundaryFaceFlux` Protocol (delete file + tests +
     SNMesh field + matvec kwarg)
   - [ ] Gates 1.3, 1.1, 1.2, 1.5 (foundation + L0)
   - [ ] Empirical default-flip decision (Gate 1.1 with MMS)
   - [ ] Gates 2.x preservation tests
   - [ ] Gates 3.x MMS (with the angularly-non-trivial ansatz)
   - [ ] Gate 4.1 (k_∞) + Gate 4.2 (trajectory_resolvent)
   - [ ] Regenerate 6 curvilinear snapshots
   - [ ] Remove 4 ERR-026 xfail-strict markers
   - [ ] Sphinx narrative + label additions (handoff to archivist)
   - [ ] `error_catalog.md` ERR-026 closure narrative
   - [ ] Closeout memo at
     `.claude/agent-memory/method-implementer/issue_168_phase_c_closeout.md`

3. **Constraints**:
   - **NEVER** introduce a new `for n in range(quad.N)` loop in the
     matvec or any code touched. Vectorise across ordinates.
   - **NEVER** shortcut the verification wiring with workarounds
     where a proper API exists or is being added. Phase C adds
     `iter_cells_by_direction`; use it. #190 will fix `Billiard` MR
     routing separately; do not work around it.
   - **NEVER** approximate the BC trace law's input with cell-centre
     values. The §16A.3 contract requires the boundary face trace.
   - **NEVER** rely on Phase A's `BoundaryFaceFlux` Protocol — it
     retires in this rewrite.
   - **MUST** preserve `apply` linearity (Gate 1.4 first; Gate 1.3
     and `apply_transpose` correctness depend on it).
   - **MUST** preserve the 5 Cartesian regression snapshot bit-
     identity (Gate 2.1) — Phase C is curvilinear-only.
   - **MUST** sequence the default flip per the empirical Gate 1.1
     probe outcome (per user 2026-05-10 decision).

---

## §8 — Cross-references

### Memory artifacts

- `.claude/agent-memory/method-implementer/issue_168_phase_a_closeout.md`
- `.claude/agent-memory/method-implementer/issue_168_phase_b_closeout.md`
- `.claude/agent-memory/literature-researcher/sphere_sn_pole_closure_canonical.md`
- `.claude/agent-memory/cross-domain-attacker/issue_168_phase_c_sweep_frame.md`
  (NEW 2026-05-12 — the structural-frame precedent for this rewrite)

### GitHub issues

- [#168](https://github.com/deOliveira-R/ORPHEUS/issues/168) — main
  issue (Phase C is the closing wave)
- [#170](https://github.com/deOliveira-R/ORPHEUS/issues/170) — F_N
  cylinder (future work; not blocking)
- [#171](https://github.com/deOliveira-R/ORPHEUS/issues/171) — 2G F_N
  (future work; not blocking)
- [#190](https://github.com/deOliveira-R/ORPHEUS/issues/190) —
  Billiard MR routing (dedicated session; Phase C uses bare functions)
- [#191](https://github.com/deOliveira-R/ORPHEUS/issues/191) —
  systemic `for n in range(quad.N)` cleanup (Phase C does not propagate
  but does not refactor existing 14 sites)

### Sphinx pages

- `docs/theory/boundary_conditions.rst` — Layer 2/3 BC architecture
- `docs/theory/discrete_ordinates.rst` — SN spatial + angular theory;
  Phase B pole-closure section + Phase C sweep-frame matvec subsection
  (to be added)
- `docs/theory/operator_algebra.rst` — `LinearOperator` /
  `OperatorSum` / `ScaledOperator` algebra

### Plan files

- `.claude/plans/sn_reshape.md` — campaign plan; Wave H Phase C row to
  be updated to RESUMED status
- `.claude/plans/issue_168_design.md` — original Phase B design memo
