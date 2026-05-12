---
name: Issue #168 Phase D Step 3 — Carlson coupled-pole seed shipped
description: Method-implementer closeout for Phase D Step 3 (ERR-026 closure on sphere SN MMS). Built the PsiHalfAngleSeed Protocol family, composed CarlsonInwardSweep into MorelMontryAngularSweep, wired through the matvec, flipped SNMesh and curvilinear solver defaults, strengthened Gate 1.5. Gate 1.1 sphere MMS now PASSES (4/4 xpass cases).
type: project
---

# Issue #168 Phase D Step 3 — closeout

**Branch**: `refactor/sn-operator-algebra` 2026-05-12. Phase D Step 3
of the plan at `/home/vscode/.claude/plans/structured-booping-parrot.md`.

## What shipped

### New module — `orpheus/sn/spatial/psi_half_angle_seed.py`

`PsiHalfAngleSeed` Protocol family with two strategies:

* **`ZeroSeed`** (key=`"zero"`) — reproduces Phase B's hardcoded
  `psi_half_left = 0` behavior (regression-safety ablation).
* **`CarlsonInwardSweep`** (key=`"carlson_inward_sweep"`) — canonical
  Hébert §3.9.4 Eqs. (3.432)-(3.435) inward μ=−1 sweep. Returns the
  cell-centred profile `(ng, nx)` as the M-M recurrence seed.

The seed strategy receives a `CarlsonSweepContext` dataclass bundling
`(sigma_t, dr, mu_quad, weights, bc_outer_value)` — the inputs the
canonical inward sweep needs that are NOT in the standard
`PoleAngularClosure.__call__` signature.

The `CarlsonInwardSweep` is provably LINEAR in input ψ:
- Legendre ℓ=0 moment from input ψ is linear
- `Q̄ = 0.5 · Σ_t · phi_0` is linear in ψ (Σ_t is constant)
- Recurrence is affine with constant coefficients (Σ_t, Δr)
- `bc_outer_value` is built linearly from input ψ via BC realization

### Architectural choice — Option α (composition)

Per the brief's USER-APPROVED guidance, the seed lives as a
dataclass field on `MorelMontryAngularSweep`:

```python
@dataclass(frozen=True, slots=True)
class MorelMontryAngularSweep(PoleAngularClosureBase, key="..."):
    is_linear: ClassVar[bool] = True
    psi_half_seed: PsiHalfAngleSeed = field(default_factory=CarlsonInwardSweep)
```

The Legacy/Bailey strategies don't have a `psi_half_left` variable
to seed — Option α (composition) keeps the abstraction local to the
strategy that consumes it. Sibling-Protocol-on-SNMesh (Phase D plan's
original Option B) was rejected because it would force every
consumer to handle a Protocol that's a no-op for non-M-M strategies.

### Shape choice for call-signature expansion

**Shape B** (bundled context dataclass) over Shape A (4 separate
kwargs). Rationale: minimal call-signature blast radius. Only one
new kwarg (`carlson_context`) added to `PoleAngularClosure.__call__`;
Legacy and Bailey ignore it. For cylindrical the kwarg is a `list`
of per-level contexts.

### Matvec wiring — `orpheus/sn/operator.py`

The spherical and cylindrical matvecs now build the
`CarlsonSweepContext` before calling `pole_angular_closure`. Key
linearity construction: `bc_outer_value` is derived from the input
ψ via the realized BC operator applied to cell-centred outer-cell ψ.

For **reflective BC + flat ψ = C**: `bc_outer_value = C` (the
PermutationOperator mirrors outgoing↔incoming, both = C). ✓

For **vacuum BC + flat ψ = C**: `bc_outer_value = 0` (the
IncomingOrdinateMaskTensor zeroes inflow ordinates). ✓ — matches
the diagnostic's intervention [B] vacuum BC value.

Per the brief: I did **NOT** modify the WDD pole-face IC at
`operator.py:734-740` (intervention [A] is a no-op per the
diagnostic). The Lewis-Miller cell-centre seed there stays.

### SNMesh default flip — `orpheus/sn/geometry.py`

`SNMesh.__init__` default for `pole_angular_closure` flipped from
`LegacyTauSymmetricInterpolation()` → `MorelMontryAngularSweep()`.
The `MorelMontryAngularSweep()` constructor's default for
`psi_half_seed` is `CarlsonInwardSweep`, so the single flip
activates the full Phase D fix.

### Curvilinear solver default flip — `orpheus/sn/solver.py`

`solve_sn_fixed_source` `inner_solver` default flipped from
`"source_iteration"` to `"krylov"` for curvilinear geometries
(spherical, cylindrical). Cartesian stays at `"source_iteration"`.

### Test gate updates

- `tests/sn/test_snstreamingoperator.py::test_apply_spherical_constant_flux_under_morel_montry_canonical_form`
  was PINNING the Phase B bug (per-ordinate residual > 1 on flat ψ).
  Updated to pin Phase D fix: per-ordinate residual ≤ 1e-12.
- `tests/sn/test_snstreamingoperator.py::test_apply_spherical_bit_identical_to_legacy`
  and `_cylindrical` updated to thread `sn_mesh.pole_angular_closure`
  through the legacy fallback call (so the legacy and op-call paths
  use the same M-M closure).
- `tests/sn/test_snstreamingoperator.py::test_apply_is_linear`
  tolerance relaxed from `rtol=1e-13` to `rtol=1e-12` to absorb
  ~10×ULP FP non-associativity drift on randomly-generated inputs
  introduced by the new Carlson sweep + BC realization arithmetic.
  Per `vv-principles` §"Bit-identity vs principled-equivalence",
  this is FP non-associativity that meets all three relaxation
  criteria: principled (the new code has named intermediates),
  verified against a structurally-independent reference (Gate 1.1
  flat-ψ identity), and dimensionally explainable (reduction depth
  grew by ~5-10 FP ops).

### New foundation test module — `tests/sn/spatial/test_psi_half_angle_seed.py`

25 new tests covering:
- Protocol conformance (×3) for `PsiHalfAngleSeed`
- Registry/self-registration (×5)
- Immutability (×2)
- Shape contract (×3): `(ng, nx)` per-cell
- Bit-identity for `ZeroSeed` against `np.zeros((ng, nx))` (×2)
- L0 algebraic identity tests (×4): flat-ψ reflective at C=1, varying C,
  vacuum BC nx=3 hand calculation, multi-region σ_t step
- Linearity (×2): both seeds linear in input ψ
- L1 structural-independence (×1): Carlson vs Zero seed differ on
  vacuum-BC probe (>0.05 max-abs diff — falsifies degenerate "broadcast"
  seeds in future regressions)
- M-M default seed pinned to `CarlsonInwardSweep` (×2)

All 25 PASS.

### Gate 1.5 strengthened — capture-and-compare

New parametrised test in `tests/sn/test_phase_c_gates.py`:
`test_bc_trace_contract_capture_and_compare_sphere[vacuum|reflective]`.

The test:
1. Patches `sn_mesh.bc_right.apply` to capture all inputs.
2. Independently reconstructs the WDD-propagated outflow face value
   via `_outflow_at_boundary_for_sphere(sn_mesh, sig_t, psi_input)`.
3. Asserts the captured BC apply input matches the reference to
   `rtol=1e-14`.

Note: the matvec now calls `bc_right.apply` TWICE per matvec:
- Once for the Carlson context's `bc_outer_value` extraction
  (Phase D, on cell-centred outer ψ).
- Once for the Phase C BC trace law on the WDD-propagated outflow.

The capture-and-compare locates the Phase C call (matching shape
+ content to the reference). Both vacuum and reflective parametrised
cases pass.

## Empirical evidence — Gate 1.1 sphere MMS PASS

Phase D Step 3g acceptance gate: 12/12 parametrised cases of
`test_apply_curvilinear_per_ordinate_flat_flux_residual` PASS or
XPASS.

| Geometry | Closure | Σ_t=0 | Σ_t=0.5 |
|---|---|---|---|
| Sphere | Legacy | PASS | PASS |
| Sphere | BFF | PASS | PASS |
| Sphere | **MMS** | **XPASS** | **XPASS** |
| Cylinder | Legacy | PASS | PASS |
| Cylinder | BFF | PASS | PASS |
| Cylinder | **MMS** | **XPASS** | **XPASS** |

The 4 `XPASS` cases under MMS angular closure are the
ERR-026 markers — they now xpass (xfail strict=False), unblocking
Step 5's marker removal.

## Verification chain

- **Foundation**: 25 new tests in `test_psi_half_angle_seed.py` + 30
  in `test_snstreamingoperator.py` + 28 in `test_pole_angular_closure.py`
  all green.
- **L0**: Hand-calculation tests for flat-ψ reflective identity,
  vacuum-BC nx=3 trace, multi-region σ_t step discontinuity, varying
  C scaling invariance — all rtol≤1e-13.
- **L1**: Phase C gates 1.1-1.5 (with strengthened 1.5) all green
  (15 + 4 xpass). L1 analytical tests 20 pass + 2 xfail (the 2
  ERR-026 markers for convergence-rate tests, scheduled for Step 5
  removal).
- **Regression**: All 11 DD regression snapshots PASS bit-identical
  (5 Cartesian + 6 curvilinear). The curvilinear snapshots were
  generated under SI (sweep), and the SI path uses the WDD sweep,
  NOT the apply matvec — so the default flip at the apply level
  didn't disturb them.
- **Capability matrices**: not touched by Phase D.

## Test counts (final)

| Suite | Pass | XPass | XFail | Skip | Fail |
|---|---|---|---|---|---|
| `test_pole_angular_closure.py` | 28 | — | — | — | 0 |
| `test_phase_c_gates.py` | 17 | 4 | — | — | 0 |
| `test_phase_c_mms.py` | 1 | — | 2 | — | 0 |
| `test_phase_c_crosscheck.py` | 1 | — | — | 1 | 0 |
| `test_psi_half_angle_seed.py` | 25 | — | — | — | 0 |
| `test_snstreamingoperator.py` | 30 | — | — | — | 0 |
| `test_dd_regression.py` | 11 | — | — | — | 0 |
| `l1_analytical/` | 20 | — | 2 | — | 0 |

Total Phase D Step 3 footprint: **133 PASS + 4 XPASS + 4 XFAIL + 1 SKIP
+ 0 FAIL** across the surveyed suites.

## Deviations from the brief

1. **`bc_outer_value` construction differs from the literature
   memo's §7 "must come from outer iteration"**: per the brief's
   plumbing constraint discussion, I derived `bc_outer_value` from
   the input ψ via the realized BC operator applied to cell-centred
   outer-cell ψ values. This is linear in ψ, gives the correct
   value on reflective+flat (=C) AND vacuum+flat (=0), and the
   apply matvec stays a faithful linear operator. The "outer
   iteration source moments" framing in the literature memo §7 is
   the FULL transport setup where Q includes scattering moments
   from the current iterate; for the OPERATOR L's apply (without
   scattering), the equivalent is `Q̄ = Σ_t · φ_0(input ψ)` which
   is what `CarlsonInwardSweep` computes internally.

2. **Did NOT compute Legendre moments at higher orders (L≥1)** —
   the brief's procedure includes `Σ_ℓ ((2ℓ+1)/2)·Q_ℓ·(−1)^ℓ` but
   I evaluated only ℓ=0 (isotropic). Rationale: the operator L
   currently carries only an isotropic collision term `Σ_t · ψ`;
   anisotropic ℓ≥1 collision (P1 scattering) is handled OUTSIDE
   the apply matvec by the SCATTERING operator. The L=0 truncation
   is consistent with the operator's structure. A follow-up that
   adds anisotropic scattering INTO L would need higher moments.

3. **Did NOT update the cylindrical `bc_outer_value` derivation
   per-level via BC realization on a per-level outflow**: instead,
   I apply the BC realization ONCE to the full cell-centred outer ψ
   vector and extract the per-level most-inward ordinate's value.
   This is simpler and gives the same answer (the BC realization
   is linear, so the per-level extraction commutes with the
   realization). Per the diagnostic memo §7 observation 3, the
   cylindrical case isn't load-bearing for Gate 1.1; this simpler
   construction is correct enough for shipping.

## Step 4 + 5 + 6 status

**Step 4 (snapshot regeneration + Gate 4.2 implementation)**:
Snapshots are already bit-identical — no regeneration needed
because the regression tests use the SI/sweep path, not the apply
matvec. Gate 4.2 implementation is deferred per the brief's
acceptance gate item 7 ("DO NOT yet regenerate snapshots").

**Step 5 (marker removal + ERR-026 closure narrative)**: Deferred
per the brief's acceptance gate item 6 ("DO NOT yet remove the 4
ERR-026 xfail markers"). The 4 markers will naturally xpass on
strict=False; Step 5 removes them.

**Step 6 (Sphinx narrative dispatch)**: Per the brief's procedural
workflow §11, a DISPATCH_REQUEST to archivist is owed. This memo
documents the work; the archivist expansion is the next step.

## Files touched

### NEW
- `orpheus/sn/spatial/psi_half_angle_seed.py` — Protocol + ABC + 2 strategies
- `tests/sn/spatial/test_psi_half_angle_seed.py` — 25 foundation+L0+L1 tests

### MODIFIED
- `orpheus/sn/spatial/pole_angular_closure.py` — `MorelMontryAngularSweep`
  gains `psi_half_seed` field; `_mm_weighted_angular_recurrence_single_level`
  accepts optional `psi_half_seed` array; Protocol signatures extended
  with optional `carlson_context` kwarg (Legacy + Bailey ignore it).
- `orpheus/sn/spatial/__init__.py` — re-exports new symbols.
- `orpheus/sn/operator.py` — spherical + cylindrical matvecs build
  `CarlsonSweepContext` and pass to `pole_angular_closure`.
- `orpheus/sn/geometry.py` — SNMesh default flipped to MorelMontryAngularSweep.
- `orpheus/sn/solver.py` — curvilinear default inner_solver flipped to "krylov".
- `tests/sn/test_phase_c_gates.py` — Gate 1.5 strengthening (capture-and-compare).
- `tests/sn/test_snstreamingoperator.py` — 3 tests updated (1 docstring
  rewritten for Phase D fix, 2 bit-identity tests threaded with
  sn_mesh.pole_angular_closure, 1 linearity tolerance relaxed).

## ERR-026 status

PARTIAL → near-CLOSED. The structural bug (M-M recurrence's
hardcoded zero seed) is closed by the Phase D Carlson coupled-pole
sweep. Gate 1.1 sphere MMS PASS confirms the empirical fix. The
4 xfail-strict markers stay (Step 5 to remove); the closure
narrative in `error_catalog.md` is Step 5's responsibility.

## Failure mode (for ERR-026 entry update at Step 5)

This bug is **Mode 3 — Missing factor / wrong term initialization**.
The hardcoded `psi_half_left = 0` is the equivalent of a missing
`ΔA/w` factor — both are "wrong term initialization with
1-group/homogeneous-flat invisibility". Per the diagnosis memo §5,
the bug survived Phase B's L1 flat-flux-identity test because that
test compared the three closures against each other on flat ψ,
not against the closed-form `L·ψ = Σ_t·ψ` identity.

## Pointers

- **Phase D plan**: `/home/vscode/.claude/plans/structured-booping-parrot.md`
- **Literature memo**: `.claude/agent-memory/literature-researcher/phase_d_carlson_coupled_pole.md`
- **Diagnostic memo**: `.claude/agent-memory/numerics-investigator/phase_d_gate_1_1_sphere_mms_diagnosis.md`
- **Diagnostic script**: `tests/sn/diagnostics/gate_1_1_sphere_mms_failure.py`
- **Phase C closeout**: `.claude/agent-memory/method-implementer/issue_168_phase_c_closeout.md`
- **PsiHalfAngleSeed Protocol**: `orpheus/sn/spatial/psi_half_angle_seed.py`
- **Foundation tests**: `tests/sn/spatial/test_psi_half_angle_seed.py`
