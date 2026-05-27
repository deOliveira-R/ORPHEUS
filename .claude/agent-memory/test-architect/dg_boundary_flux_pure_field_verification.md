---
name: dg-boundary-flux-pure-field-verification
description: D-G verification spec — BoundaryFlux pure-Field carve + SweepScratch split. Pre-D-G regression snapshot, three new test specs (transport.fields, sweep-scratch-split, immutability), failure-mode ranking, principled-equivalence judgment with ULP bound.
metadata:
  type: project
---

# D-G Verification Plan — BoundaryFlux pure Field + SweepScratch split

**Step**: D-G of Depth B (per `.claude/plans/depth_b_field_on_function_space.md`)
**Carve type**: operator-algebra subsystem-crossing — `feedback_no_method_implementer_for_surgical_carves` applies; main agent + turn-by-turn user steering
**Authoring agent**: test-architect (this memo)
**Decision required before implementation**: D-G is PRINCIPLED, not bit-identical; the proof is structured below

This memo is the **verification specification** that MUST be in place BEFORE D-G implementation lands. It is the test-architect's deliverable per the `subagent-handoff-protocol` proactive-trigger row "operator-algebra carve crossing subsystem boundaries".

The carve crosses THREE subsystem boundaries simultaneously:

1. **Mutability boundary** — pre-D-G mutable write-through cache → post-D-G immutable Field with sweep-side scratch.
2. **Storage boundary** — pre-D-G structured-bundle (4 named buffers, two of which are `None` per geometry) → post-D-G flat-buffer `(layout.total_size,)` with `FaceLayout` slice views.
3. **Ownership boundary** — pre-D-G interior wavefront cache conflated INSIDE BoundaryFlux's 2-D `xmin_xmax_buf`/`ymin_ymax_buf` → post-D-G boundary trace lives on BoundaryFlux; interior wavefront cache lives on sweep-private `SweepScratch`.

Per `vv-principles` "Bit-identity vs principled-equivalence", this is the textbook PRINCIPLED case: the reduction tree changes (face slice arithmetic decouples from interior cell arithmetic), the named intermediates become principled (boundary trace is named; wavefront cache is named), and the underlying numerical algorithm (per-cell DD update) is unchanged.

---

## 1. The claim taxonomy — what each test row proves

Per `vv-principles` §"Hierarchical claim taxonomy", each D-G test row is gated on the claim layer it proves AND the pillar it draws evidence from.

| Claim layer | Test row | Pillar | Reference |
|---|---|---|---|
| Convergence-order | `test_sn_p1_aniso_mms_converges_second_order` (existing L1 MMS) | MMS | `build_p1_aniso_mms_case` source-driven |
| Convergence-order | `test_mms_curvilinear_aniso_dd_convergence` (existing L1) | MMS | curvilinear MMS source-driven |
| Flux-shape | `tests/sn/regression/test_dd_regression.py` (existing) | Closed-form via cross-check | DD snapshot was originally pinned against analytical + Lewis-Miller benchmarks |
| Flux-shape | 2-D octant equivalence snapshots (existing 6 cases) | Closed-form fixed-source + DD self-consistency | snapshots pin geometry-resolved expected fluxes |
| Eigenvalue | `tests/sn/l1_analytical/test_kinf_homogeneous.py` (existing) | Closed-form | `k_inf = λ_max(A⁻¹F)`, transfer matrix |
| Bit-identity (algebra) | `tests/transport/fields/test_boundary_flux.py::test_arithmetic_round_trip` (NEW) | Self-consistency | Field ABC dunder math |
| Equivalence (sweep) | `tests/sn/test_sweep_scratch_split.py` (NEW) | Pre-D-G regression baseline | captured-output comparison |
| Type-discipline | `tests/sn/test_boundary_flux_immutability_invariant.py` (NEW) | Software invariant | FrozenInstanceError gate |

**Pillar selection notes:**

- MMS is the correct pillar for convergence-order claims because the source is structurally independent of the production code. The L1 MMS gates ARE structurally independent for the BoundaryFlux carve — the manufactured source machinery in `orpheus.derivations.continuous.mms.sn` does NOT touch boundary buffers.
- Eigenvalue claims use the transfer-matrix closed form (per `vv-principles` "MMS does NOT prove eigenvalues"); these are existing tests and must stay green.
- Bit-identity is the gate for the FIELD ALGEBRA layer of D-G (arithmetic on the flat buffer must agree with arithmetic on the legacy named buffers when reduced to the same face values). It is NOT a bit-identity gate on the SWEEP output — see §5 for that judgment.

---

## 2. Pre-D-G regression snapshot script — the bit-identity gate

The "10 pre-existing DD-regression failures stay AT the same failure set" invariant (plan §11.1) requires a captured baseline BEFORE step D-G's first commit. The baseline IS the regression gate.

### 2.1 Capture-baseline invocation

Run ONCE on the pre-D-G tree (the HEAD that precedes the first D-G commit). The output is committed to the worktree as the gate file.

```bash
# From repo root with .venv active. -O is canonical (per
# feedback_default_test_mode_is_optimize).
mkdir -p tests/sn/regression
python -O -m pytest tests/ \
    --tb=no -q \
    --no-header \
    -o "console_output_style=count" \
    > tests/sn/regression/_pre_dg_baseline.txt 2>&1 \
    || true   # capture failures; do NOT exit on non-zero

# Capture pass/fail set with EXACT test node IDs and result codes.
# This is the comparison artifact.
python -O -m pytest tests/ \
    --co -q                            > tests/sn/regression/_pre_dg_collected.txt
python -O -m pytest tests/ \
    --tb=no -rN -q --no-header        2>&1 \
    | grep -E "^(PASSED|FAILED|ERROR|SKIPPED) " \
    | sort -u                          > tests/sn/regression/_pre_dg_pass_fail.txt \
    || true
```

The file `tests/sn/regression/_pre_dg_pass_fail.txt` is the LOAD-BEARING artifact. It MUST be committed in the prep commit IMMEDIATELY BEFORE the first D-G code change.

### 2.2 The re-runnable gate

At every D-G sub-commit, the gate runs as:

```bash
# Re-runnable: at every D-G sub-commit (and at HEAD before committing).
python -O -m pytest tests/ \
    --tb=no -rN -q --no-header        2>&1 \
    | grep -E "^(PASSED|FAILED|ERROR|SKIPPED) " \
    | sort -u                          > /tmp/_post_dg_pass_fail.txt \
    || true

# The gate: diff must be EMPTY.
diff tests/sn/regression/_pre_dg_baseline_pass_fail.txt /tmp/_post_dg_pass_fail.txt
GATE_STATUS=$?
if [ $GATE_STATUS -ne 0 ]; then
    echo "D-G regression-set DRIFT detected. Pre-existing failures changed."
    exit 1
fi
```

### 2.3 Targeted MMS sub-gate (load-bearing)

The MMS gates are L1 numerical-correctness pins that MUST be green under both pre-D-G and post-D-G. Run as a focused gate at every D-G sub-commit:

```bash
python -O -m pytest \
    tests/sn/test_mms_aniso.py \
    tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py \
    tests/sn/l1_analytical/test_kinf_homogeneous.py \
    tests/sn/l1_analytical/test_kinf_homogeneous_tolerance.py \
    tests/sn/regression/test_dd_regression.py \
    tests/sn/test_2d_octant_sweep_equivalence.py \
    tests/sn/test_boundary_flux_arithmetic.py \
    tests/sn/test_angular_flux_with_boundary.py \
    tests/sn/test_invertible_operator.py \
    -v --tb=short
```

**This sub-gate MUST stay green at every D-G sub-commit.** If any test fails:

1. Confirm the failure was not in the pre-D-G baseline (`_pre_dg_pass_fail.txt`).
2. If it's a NEW failure, the D-G commit MUST roll back; the carve broke a load-bearing test.

### 2.4 Pre-existing 10 DD-regression failures — captured semantics

Per plan §11.1, 10 pre-existing failures exist in the regression suite at the D-G entry point. These are pre-existing technical debt; D-G is NOT required to fix them. The gate is `same failure set`, not `zero failures`.

The pre-D-G baseline captures the exact failure set; the gate diff catches:

- Any NEW failure (a regression D-G introduced) — BLOCKING.
- Any test that flipped from FAILED → PASSED — investigate (could indicate D-G accidentally affected a regression; verify it's a legitimate fix, document, then accept).
- Any test that flipped from PASSED → FAILED — BLOCKING regression.

---

## 3. New test specifications (pytest skeletons, not implementations)

Below are the THREE new test modules D-G requires. Each is specified as a pytest skeleton — the main agent implements the bodies during D-G turn-by-turn.

### 3.1 `tests/transport/fields/test_boundary_flux.py` — pure Field algebra + flat-buffer + FaceLayout

```python
"""L0/L2 — BoundaryFlux pure-Field algebra, FaceLayout slice views, and
flat-buffer round-trips.

The post-D-G BoundaryFlux inherits Field. This module pins:

* Field algebra inherited correctly (no hand-coded dunders survive).
* The flat buffer + FaceLayout descriptor produces correct per-face
  slice views with NO COPIES (memory-shared with the parent buffer).
* Mesh-binding rejection survives (cross-mesh BoundaryFlux arithmetic
  raises).
* Round-trip equivalence: legacy face-slot writes match the new
  flat-buffer layout when slice views are reshaped to the legacy shape.

References
----------

* Depth B plan §3.4 (Option Ω flat-buffer storage).
* `vv-principles` §"Bit-identity vs principled-equivalence" — this
  test is the bit-identity gate at the FIELD-ALGEBRA layer (per
  §1 of this memo).
"""
from __future__ import annotations

import numpy as np
import pytest
from dataclasses import FrozenInstanceError

from orpheus.geometry import BC, Mesh1D, Mesh2D, Region, RegionMesh, StructuredGeometry
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.geometry import SNMesh
from orpheus.transport.fields.boundary_flux import BoundaryFlux
from tests.sn._test_helpers import placeholder_materials


pytestmark = [pytest.mark.foundation]


# ───────────────────────────────────────────────────────────────────────
# Geometry fixtures (reused from legacy test_boundary_flux_arithmetic)
# ───────────────────────────────────────────────────────────────────────

@pytest.fixture
def slab_mesh() -> SNMesh:
    """1-D slab — two faces (xmin, xmax), each (N, ng)."""
    ...


@pytest.fixture
def sphere_mesh() -> SNMesh:
    """1-D sphere — single outer face (xmax)."""
    ...


@pytest.fixture
def cartesian_2d_mesh() -> SNMesh:
    """2-D Cartesian — four faces (xmin, xmax, ymin, ymax). nx=3, ny=2."""
    ...


# ───────────────────────────────────────────────────────────────────────
# A. Field-ABC inheritance and algebra
# ───────────────────────────────────────────────────────────────────────

class TestFieldAlgebraInherited:
    """BoundaryFlux uses Field's __add__/__sub__/__mul__/etc — NO hand-coded
    dunder body survives in the post-D-G implementation."""

    def test_add_returns_new_instance_slab(self, slab_mesh):
        """L0: `bf1 + bf2` returns a fresh BoundaryFlux; originals unchanged.

        Catches: hand-coded dunder accidentally still in module (would
        cause attribute-lookup failures in flat-buffer mode).
        """
        ...

    def test_sub_returns_new_instance_sphere(self, sphere_mesh):
        """L0: same shape for sphere (single face)."""
        ...

    def test_scalar_mul_2d(self, cartesian_2d_mesh):
        """L0: scalar multiplication propagates through flat buffer in ONE numpy call."""
        ...

    def test_scalar_div_propagates(self, slab_mesh):
        ...

    def test_neg_2d(self, cartesian_2d_mesh):
        ...

    def test_distributive_property_2d(self, cartesian_2d_mesh):
        """L1: (bf1 + bf2) * c == bf1*c + bf2*c on the flat buffer.

        Catches: a bug where flat-buffer arithmetic loses associativity
        because the FaceLayout reshape is non-contiguous.
        """
        ...

    def test_linf_via_inherited_property(self, cartesian_2d_mesh):
        """Inherited Field.linf — verifies the property hasn't been
        accidentally shadowed."""
        ...

    def test_l2_via_inherited_property(self, slab_mesh):
        """Inherited Field.l2 — must use FaceLayout-aware metric (if any)
        OR fall back to default Euclidean. Pin the convention."""
        ...

    def test_inner_product_with_self_is_l2_squared(self, sphere_mesh):
        """<bf, bf> == ||bf||^2_2. Algebraic identity."""
        ...


# ───────────────────────────────────────────────────────────────────────
# B. FaceLayout — slice views, no copies, all geometries
# ───────────────────────────────────────────────────────────────────────

class TestFaceLayoutSliceViews:
    """The FaceLayout descriptor maps face names to slices into the flat
    backing buffer. Per-face access via `bf.faces[name]` returns a
    BoundaryFaceFlux that is a VIEW (no copy).
    """

    def test_slab_layout_has_xmin_and_xmax(self, slab_mesh):
        """L0: slab BoundaryFlux exposes exactly two faces, names 'xmin', 'xmax'."""
        ...

    def test_sphere_layout_has_only_xmax(self, sphere_mesh):
        """L0: sphere BoundaryFlux exposes ONE face, name 'xmax'.

        Catches: the geometry-conditional zeros_for_sn_mesh factory
        accidentally allocating an xmin face on a curvilinear mesh.
        """
        ...

    def test_cartesian_2d_layout_has_four_faces(self, cartesian_2d_mesh):
        """L0: 2-D BoundaryFlux exposes four faces: xmin, xmax, ymin, ymax."""
        ...

    def test_face_view_is_memory_shared_slab(self, slab_mesh):
        """L0: bf.faces['xmin'].values shares memory with bf.values.

        Catches: a bug where FaceLayout.slice_view returns a COPY
        instead of a view — sweep writes would be lost.

        Probe: write through the face view; the flat buffer reflects it.
        """
        # bf = BoundaryFlux.zeros_for_sn_mesh(slab_mesh)
        # face = bf.faces['xmin']
        # assert face.values.base is bf.values  # memory-shared
        ...

    def test_face_view_is_memory_shared_2d_xmin(self, cartesian_2d_mesh):
        """L0: 2-D xmin face is a view into the flat buffer at offset 0,
        reshape (N, ng, ny). Memory-shared, not a copy."""
        ...

    def test_face_shapes_match_per_geometry(self, slab_mesh, sphere_mesh, cartesian_2d_mesh):
        """L0: each face's shape from FaceLayout matches the per-geometry
        expected shape per plan §3.4:

        * 1-D slab xmin: (N, ng); same for xmax.
        * 1-D sphere xmax: (N, ng); no xmin.
        * 2-D xmin: (N, ng, ny); xmax: (N, ng, ny); ymin: (N, ng, nx); ymax: (N, ng, nx).

        Catches: a layout bug where 2-D face slices have the wrong
        spatial axis.
        """
        ...

    def test_total_size_consistent_with_face_sizes(self, cartesian_2d_mesh):
        """L0: layout.total_size == sum(slot.flat_size for slot in faces.values()).

        Catches: face overlap (sweep would write the same cell from two
        faces) or face gaps (a slot allocated but never sliced).
        """
        ...

    def test_no_interior_cells_in_layout(self, cartesian_2d_mesh):
        """L0/L2 — CRITICAL: post-D-G 2-D BoundaryFlux contains ONLY face
        slices. The pre-D-G `(N, ng, nx+1, ny)` interior-buffer shape
        is GONE.

        Catches: the most load-bearing carve invariant — interior
        wavefront cache must NOT live in BoundaryFlux. If
        layout.total_size includes interior cells, the split was
        incomplete.

        Concretely (2-D, nx=3, ny=2):

            xmin: N*ng*ny       =  N*ng*2     # ny = 2
            xmax: N*ng*ny       =  N*ng*2
            ymin: N*ng*nx       =  N*ng*3     # nx = 3
            ymax: N*ng*nx       =  N*ng*3
            ───────────────────────────────
            total: N*ng*(2*ny + 2*nx) = N*ng*10

        NOT N*ng*((nx+1)*ny + nx*(ny+1)) = N*ng*(8+9) = N*ng*17
        (the legacy buffer size including interior face columns).
        """
        ...


# ───────────────────────────────────────────────────────────────────────
# C. Mesh-binding rejection (preserved from legacy)
# ───────────────────────────────────────────────────────────────────────

class TestMeshBindingRejection:
    def test_cross_mesh_add_rejected(self, slab_mesh):
        """L0: BoundaryFlux on distinct SNMesh instances cannot be added
        even when structurally identical (preserved from legacy)."""
        ...

    def test_wrong_type_add_rejected(self, slab_mesh):
        """L0: bf + 42 raises (Field._check_partner class identity gate)."""
        ...


# ───────────────────────────────────────────────────────────────────────
# D. Flat-buffer round-trip (the new invariant)
# ───────────────────────────────────────────────────────────────────────

class TestFlatBufferRoundTrip:
    """The flat backing buffer is the source of truth; faces are views.
    Writes through faces propagate to the flat buffer; the converse
    (writes to the flat buffer reflect in the face views) also holds.

    Field is frozen (immutable from the OUTSIDE), but the underlying
    ndarray data IS writable for low-level testing. See immutability
    test module for the public-API frozen-ness verification.
    """

    def test_face_write_reflected_in_flat_buffer(self, slab_mesh):
        """L0: writing to a face view reflects in the flat buffer
        (memory-shared). Probe: bf.faces['xmin'].values[...] = 1.0;
        assert flat region for xmin is all-1.0."""
        ...

    def test_flat_write_reflected_in_face(self, cartesian_2d_mesh):
        """L0: writing the flat buffer slice reflects in the face view
        (the inverse direction)."""
        ...

    def test_round_trip_arithmetic(self, slab_mesh):
        """L1: after `out = bf1 + bf2`, every face of `out` equals the
        sum of the corresponding faces from bf1 and bf2.

        Bit-identity at the FACE-VALUE level. Catches: a FaceLayout
        bug where slice offsets are misaligned."""
        ...


# ───────────────────────────────────────────────────────────────────────
# E. Construction factory and inheritance from Field
# ───────────────────────────────────────────────────────────────────────

class TestConstruction:
    def test_zeros_for_sn_mesh_slab(self, slab_mesh):
        """L0: BoundaryFlux.zeros_for_sn_mesh(slab) returns a Field-typed
        instance with all-zero values."""
        ...

    def test_zeros_for_sn_mesh_sphere(self, sphere_mesh):
        ...

    def test_zeros_for_sn_mesh_2d(self, cartesian_2d_mesh):
        ...

    def test_inherits_field(self):
        """L0: isinstance(BoundaryFlux.zeros_for_sn_mesh(mesh), Field) is True.

        Catches: a refactor that detaches BoundaryFlux from the L1 Field
        hierarchy (would unsubscribe it from operator algebra)."""
        ...

    def test_post_init_validates_shape(self, slab_mesh):
        """L0: constructing with values.shape != (layout.total_size,)
        raises ValueError at __post_init__.

        Catches: shape-validation regression (per Field's contract)."""
        ...
```

### 3.2 `tests/sn/test_sweep_scratch_split.py` — sweep output equivalence after the buffer split

```python
"""L1 — Sweep output bit-identity (or principled-equivalent) AFTER
the BoundaryFlux / SweepScratch buffer split.

The pre-D-G 2-D `xmin_xmax_buf` / `ymin_ymax_buf` buffers carry BOTH
boundary face slices AND the wavefront interior face-flux cache cells
(positions `[:, :, 1, :]` ... `[:, :, nx-1, :]` and `[:, :, :, 1]`
... `[:, :, :, ny-1]`). Post-D-G:

* BoundaryFlux carries ONLY positions 0 and nx/ny (boundary trace).
* SweepScratch carries positions 1..nx-1 / 1..ny-1 (interior cache).

This split MUST preserve the sweep's reflective-BC contract and
the angular-flux output. The 2-D case is the LOAD-BEARING regime
because the conflation lives there.

References
----------

* Depth B plan §3.4 (storage split).
* `vv-principles` §"Bit-identity vs principled-equivalence" — D-G is
  classified PRINCIPLED; ULP-level drift bounded by
  `(reduction depth) × ULP` is acceptable.
"""
from __future__ import annotations

import numpy as np
import pytest

# Tag as foundation — software invariant on the sweep's output
# consistency under the storage refactor.
pytestmark = [pytest.mark.foundation]


# ───────────────────────────────────────────────────────────────────────
# A. Reflective BC reads fresh face state across iterations (CRITICAL)
# ───────────────────────────────────────────────────────────────────────

class TestReflectiveBCStatePersistence:
    """The reflective BC reads outgoing-face flux from the previous
    sweep iteration to compute the partner's incoming flux. Post-D-G,
    BoundaryFlux is IMMUTABLE — the sweep MUST construct a fresh
    BoundaryFlux at each iteration carrying the latest outgoing-face
    state, and the next iteration MUST read from THAT BoundaryFlux.

    THIS IS THE HIGHEST-RISK FAILURE MODE (see §4 of the memo).
    """

    def test_2d_reflective_bc_two_iterations_match_legacy(self):
        """L1: run a 2-D reflective-BC fixed-source problem with manual
        SI driver (two iterations). Capture the angular_flux and
        boundary_flux at the end of each iteration. Compare against
        the pre-D-G regression snapshot.

        Catches: a bug where the new sweep loses the previous-iteration
        outgoing face state (reflective BC reads stale zeros instead
        of the previous iter's outgoing).

        Snapshot: `tests/sn/regression/snapshots/dg_2d_reflective_2iter.npz`
        captured pre-D-G as the bit-identity gate.
        """
        ...

    def test_2d_reflective_bc_converged_solution_match(self):
        """L1: full SI-converged 2-D reflective-BC problem. Compare
        converged scalar flux against pre-D-G snapshot at rtol bounded
        by `(iter_count × ULP)`.

        Reference: `02_reflective_1g_homog_uniformQ_LS4` snapshot from
        the existing regression suite. The flux MUST agree to within
        rtol=1e-13 (principled equivalence; tighter than typical
        regression because the only change is buffer layout, NOT
        algorithm).
        """
        ...

    def test_1d_slab_reflective_match(self):
        """L1: 1-D slab reflective-BC case. NOT a buffer-split case
        (1-D xmin_face/xmax_face are already pure face buffers
        pre-D-G), but tests that the 1-D regime is unaffected by the
        D-G refactor.

        Reference: `slab_2g_3reg_dd_n40.npz` or similar existing snapshot.
        """
        ...


# ───────────────────────────────────────────────────────────────────────
# B. Interior wavefront cache lives on SweepScratch (the split invariant)
# ───────────────────────────────────────────────────────────────────────

class TestInteriorWavefrontCacheLocation:
    """Post-D-G, the interior face-flux cache (cells 1..nx-1 / 1..ny-1
    of the legacy 2-D buffers) lives on a sweep-private SweepScratch
    type, NOT on BoundaryFlux. This contract is checked by inspecting
    BoundaryFlux's flat-buffer total_size.
    """

    def test_2d_boundary_flux_size_is_face_only(self):
        """L0/L2: BoundaryFlux for a 2-D mesh has
        `values.size == N*ng*(2*ny + 2*nx)`, NOT
        `N*ng*((nx+1)*ny + nx*(ny+1))`.

        Catches: the buffer split was botched; interior cells are
        still inside BoundaryFlux.
        """
        ...

    def test_sweep_scratch_carries_interior_cache(self):
        """L0: SweepScratch is initialised by the sweep on first call
        and carries `(N, ng, nx+1, ny)` and `(N, ng, nx, ny+1)`
        arrays — the FULL legacy buffer shape, but now owned by the
        sweep operator.

        Catches: SweepScratch shape regression (sweep would silently
        recompute interior values on every iteration).
        """
        ...

    def test_sweep_scratch_persists_across_iterations(self):
        """L1: SweepScratch is the same Python object across multiple
        sweep calls on the same operator. Pre-D-G the cache was
        carried via BoundaryFlux's write-through; post-D-G it's
        operator-private.

        Catches: a regression where SweepScratch is re-allocated
        every sweep (memory churn the production hot path cannot
        afford — exactly the constraint the pre-D-G mutability
        docstring called out at lines 41-43).
        """
        ...


# ───────────────────────────────────────────────────────────────────────
# C. Output bit-identity / principled-equivalence vs pre-D-G snapshots
# ───────────────────────────────────────────────────────────────────────

class TestOutputEquivalenceVsPreDG:
    """Post-D-G sweep output MUST agree with pre-D-G to within
    `(iter_count × ULP)` — the principled-equivalence bound. The
    underlying algorithm hasn't changed; only the storage location of
    the wavefront cache.
    """

    @pytest.mark.parametrize("snapshot_name", [
        "2d_octant_equivalence_01_smoke_vacuum_1g_homog_uniformQ_LS4",
        "2d_octant_equivalence_02_reflective_1g_homog_uniformQ_LS4",
        "2d_octant_equivalence_03_l7_trap_mixedBC_2g_het_LS4",
        "2d_octant_equivalence_04_vacuum_2g_het_gradientQ_LS6",
        "2d_octant_equivalence_05_qaniso_mixedBC_2g_het_LS4",
        "2d_octant_equivalence_06_purez_vacuum_1g_Lebedev5",
    ])
    def test_2d_snapshot_match(self, snapshot_name):
        """L1: every 2-D regression snapshot reproduces post-D-G.

        Tolerance: assert_array_equal (strict bit-identity) where the
        sweep is a single call; assert_allclose(rtol=1e-14) where SI
        iteration is involved (principled bound).

        Catches: numerical drift introduced by the buffer split —
        e.g. a reshape order error that subtly reorders summation
        and changes ULP-level outputs.
        """
        ...


# ───────────────────────────────────────────────────────────────────────
# D. Per-ordinate flat-flux residual (the curvilinear-bug detector)
# ───────────────────────────────────────────────────────────────────────

class TestPerOrdinateFlatFluxResidual:
    """Inherited from `vv-principles` Signature 1 / ERR-006 / ERR-026:
    per-ordinate flat-flux residual is the canonical curvilinear
    sweep correctness check. D-G is 2-D Cartesian (not curvilinear),
    but the same principle catches angular-redistribution bugs that
    happen to be masked by the storage refactor.
    """

    def test_2d_cartesian_per_ordinate_flat_flux_unchanged(self):
        """L1: set Q uniform, sig_t uniform, vacuum BC, and verify
        ψ_n(r) = Q/Σ_t for every ordinate after one sweep.

        Catches: a buffer-split bug that affects a single ordinate's
        in-plane streaming differently from the others.
        """
        ...


# ───────────────────────────────────────────────────────────────────────
# E. MMS gates (load-bearing) — fresh runs after D-G
# ───────────────────────────────────────────────────────────────────────

class TestMMSGatesPostDG:
    """The L1 MMS tests live in tests/sn/test_mms_aniso.py and
    tests/sn/l1_analytical/. They are not re-implemented here, but
    THIS module records the contract that those gates must stay green
    at every D-G commit (the targeted MMS sub-gate per §2.3 of the
    memo).

    This class contains placeholder tests that ASSERT the L1 gates
    are runnable and unaffected — they import and run a tight
    sub-set as a smoke check.
    """

    def test_p1_aniso_mms_smoke(self):
        """L1 smoke: run the P1 aniso MMS at a single nx and verify
        the L2 error is below the expected threshold for that nx.

        This is not a convergence rate check (that's in the actual
        L1 file); it's a smoke check that the MMS pipeline is
        operational under the D-G refactor.
        """
        ...
```

### 3.3 `tests/sn/test_boundary_flux_immutability_invariant.py` — frozen-instance gate

```python
"""L0 — BoundaryFlux is functionally immutable post-D-G.

Pre-D-G BoundaryFlux was mutable by design (write-through cache for
sweep persistent BC state, per the pre-D-G module docstring at lines
34-43). Post-D-G:

* `@dataclass(frozen=True)` makes attribute assignment raise
  FrozenInstanceError.
* The sweep's "write-through" pattern is replaced by sweep-side
  SweepScratch + fresh-BoundaryFlux construction on each iteration.

This module pins the immutability contract. Per `vv-principles`
anti-pattern #4 / Field's design contract, an immutable Field with
a documented frozen=True is a structural guarantee the sweep
cannot accidentally violate.

The test is NEGATIVE: every mutation attempt MUST raise.

References
----------

* Depth B plan §3.4: "BoundaryFlux is `frozen=True` and `kw_only=True`."
* Plan §3.4: "The pre-D-G mutability (write-through cache for sweep
  persistent BC state) is replaced by sweep-side cache + functional
  reconstruction at iteration boundaries."
"""
from __future__ import annotations

import numpy as np
import pytest
from dataclasses import FrozenInstanceError

pytestmark = [pytest.mark.foundation]


# ───────────────────────────────────────────────────────────────────────
# A. Frozen-dataclass enforcement (Python's structural gate)
# ───────────────────────────────────────────────────────────────────────

class TestFrozenInstance:
    """Direct attribute assignment on a BoundaryFlux instance MUST raise
    FrozenInstanceError per the @dataclass(frozen=True) contract."""

    def test_assign_values_raises(self, slab_mesh):
        """L0: bf.values = np.zeros(...) raises FrozenInstanceError."""
        # bf = BoundaryFlux.zeros_for_sn_mesh(slab_mesh)
        # with pytest.raises(FrozenInstanceError):
        #     bf.values = np.zeros(bf.values.shape)
        ...

    def test_assign_layout_raises(self, slab_mesh):
        """L0: bf.layout = other_layout raises."""
        ...

    def test_assign_space_raises(self, slab_mesh):
        """L0: bf.space = other_space raises (inherited from Field's
        frozen=True)."""
        ...

    def test_assign_mesh_raises(self, slab_mesh):
        """L0: bf.mesh = other_mesh raises."""
        ...


# ───────────────────────────────────────────────────────────────────────
# B. NEGATIVE — the legacy mutability path is GONE
# ───────────────────────────────────────────────────────────────────────

class TestLegacyMutabilityRetired:
    """Pre-D-G, BoundaryFlux.zeros(mesh).xmin_face[...] = 1.0 was the
    canonical write path. Post-D-G:

    * `xmin_face` no longer exists as a top-level dataclass field
      (it's a SLICE-VIEW property derived from the flat buffer).
    * Writing to `bf.faces['xmin'].values[...]` writes through to the
      flat buffer (memory-shared view). This IS still writable at
      the array level — Python ndarray data is not made read-only
      by the dataclass frozen=True flag.

    The contract: the PUBLIC API of BoundaryFlux is frozen
    (attribute assignment fails); the underlying ndarray data is
    writable, which is required for the sweep to construct a fresh
    BoundaryFlux from outgoing-face writes. The sweep does:

        outgoing_face_buf = np.empty(face_shape)
        # ... sweep writes to outgoing_face_buf ...
        new_bf = BoundaryFlux(values=flat_concat(outgoing_face_buf, ...),
                              space=..., layout=..., mesh=...)

    NOT:

        bf.faces['xmax'].values[...] = outgoing  # write-through

    The second form WOULD work (memory-shared view), but is forbidden
    by convention — the sweep MUST construct a fresh BoundaryFlux on
    each iteration. This is verified by the SweepScratch split test
    module (test_sweep_scratch_split.py).
    """

    def test_xmin_face_attribute_no_longer_exists(self, slab_mesh):
        """L0: the legacy bf.xmin_face attribute is GONE post-D-G.

        Accessing bf.xmin_face raises AttributeError. To get the
        xmin slice, use bf.faces['xmin'].

        Catches: a half-completed migration where the legacy attribute
        names linger as @property fallbacks.
        """
        # bf = BoundaryFlux.zeros_for_sn_mesh(slab_mesh)
        # with pytest.raises(AttributeError):
        #     _ = bf.xmin_face
        ...

    def test_xmin_xmax_buf_attribute_no_longer_exists(self, cartesian_2d_mesh):
        """L0: bf.xmin_xmax_buf is GONE. This was the load-bearing
        pre-D-G buffer that conflated face + interior cells.

        Catches: the buffer split was not actually done — the
        legacy 2-D buffer still lives on the dataclass."""
        ...

    def test_ymin_ymax_buf_attribute_no_longer_exists(self, cartesian_2d_mesh):
        """L0: bf.ymin_ymax_buf is GONE."""
        ...


# ───────────────────────────────────────────────────────────────────────
# C. The functional path (construction-only mutation)
# ───────────────────────────────────────────────────────────────────────

class TestFunctionalConstructionPath:
    """The post-D-G mutation pattern is `dataclasses.replace(bf,
    values=new_flat_buf)`. This is the principled path for the sweep
    to produce a new BoundaryFlux from outgoing face writes."""

    def test_replace_creates_new_instance(self, slab_mesh):
        """L0: `dataclasses.replace(bf, values=new_arr)` returns a
        fresh BoundaryFlux with the updated values; the original is
        unmodified."""
        ...

    def test_arithmetic_preserves_frozen_contract(self, slab_mesh):
        """L0: bf1 + bf2 returns a frozen instance (assignment on it
        still raises). The new instance is constructed via the
        Field-inherited dunder, which uses type(self)(...) — same
        frozen path."""
        ...
```

---

## 4. Highest-risk single failure modes — ranked

Per `vv-principles` failure-mode taxonomy and `numerical-bug-signatures` recognition catalog, the D-G carve has the following ranked failure modes. Each lists detection strategy and which of the three new test modules catches it.

### Rank 1 (HIGHEST) — Reflective BC reads stale face state from a non-updated BoundaryFlux

**Failure mode**: post-D-G, the sweep produces a new BoundaryFlux at each iteration carrying the outgoing-face writes. The reflective BC partner reads from this BoundaryFlux. If the sweep does NOT construct the fresh BoundaryFlux before the next iteration reads the partner face, the BC reads stale zeros (or previous-iteration values) instead of the current iteration's outgoing flux.

**Mechanism**: the pre-D-G mutable write-through pattern coupled the partner reads directly to the sweep writes via shared mutable ndarray. Post-D-G immutable BoundaryFlux requires explicit fresh-construction between sweep calls.

**Detection**: `tests/sn/test_sweep_scratch_split.py::TestReflectiveBCStatePersistence::test_2d_reflective_bc_two_iterations_match_legacy` AND `::test_2d_reflective_bc_converged_solution_match`.

**Why it hides**: a single sweep call produces the correct answer; the bug only manifests after the SECOND sweep call (where the reflective partner reads stale state). Single-iteration regression tests would pass.

**Catalog signature**: combines failure mode #2 (variable swap — confusion between "stale BC buffer" and "fresh BC buffer") and failure mode #4 (wrong recursion — the iteration loop omits the fresh-construction step).

### Rank 2 — FaceLayout slice offset / reshape order mismatch in 2-D

**Failure mode**: 2-D `xmin` face values are written to a slot in the flat buffer that the FaceLayout maps to `ymin`. The reflective BC partner sees the WRONG face's values.

**Mechanism**: 2-D has four faces with different per-face shapes — `(N, ng, ny)` for x-faces, `(N, ng, nx)` for y-faces. A FaceLayout slot whose `shape` field is `(N, ng, nx)` but whose `flat_size` is computed as `N*ng*ny` (a typo / convention drift) silently mis-maps.

**Detection**: `tests/transport/fields/test_boundary_flux.py::TestFaceLayoutSliceViews::test_face_shapes_match_per_geometry` AND `::test_total_size_consistent_with_face_sizes`. The 2-D snapshot match in `TestOutputEquivalenceVsPreDG::test_2d_snapshot_match` is the integration-level catch.

**Why it hides**: 1-D tests pass (only two faces, both shape `(N, ng)` — a slot mis-map is structurally impossible); the bug is gated on dimensionality.

**Catalog signature**: failure mode #5 (index error — wrong offset into the flat buffer) AND failure mode #6 (convention drift — slot shape vs flat size).

### Rank 3 — SweepScratch re-allocated on every sweep (perf regression masquerading as correctness)

**Failure mode**: the SweepScratch type is constructed fresh inside the sweep function body (NOT held by the sweep operator across calls). Memory churn becomes the production bottleneck — exactly the constraint the pre-D-G mutability docstring (lines 41-43) called out.

**Mechanism**: the sweep operator's call signature is `sweep(Q, sig_t, sn_mesh, boundary_flux)` pre-D-G. Post-D-G it needs to be `sweep(Q, sig_t, sn_mesh, boundary_flux, sweep_scratch)` where `sweep_scratch` is owned by the SNSweepOperator instance (or by the InvertibleOperator that wraps it). If the carve writes `sweep_scratch = SweepScratch.zeros_for_mesh(sn_mesh)` INSIDE the sweep function body, the cache is rebuilt every call.

**Detection**: `tests/sn/test_sweep_scratch_split.py::TestInteriorWavefrontCacheLocation::test_sweep_scratch_persists_across_iterations`.

**Why it hides**: correctness tests pass (the values are correct; only allocation count is wrong). A correctness-only gate would not catch this — the gate must inspect object identity across iterations.

### Rank 4 — Cross-mesh BoundaryFlux arithmetic accidentally permitted

**Failure mode**: post-D-G `Field._check_partner` enforces same-class + same-space; SAME-SPACE is `(name, shape)` + units dimensionality. Two BoundaryFlux instances built from DIFFERENT SNMesh instances but matching shapes would have the SAME `FunctionSpace` (because the space is "sn_boundary_flat" + shape `(total_size,)`). The Field-level gate would PERMIT cross-mesh arithmetic.

**Mechanism**: the legacy `_validate_partner` had a `self.mesh is not other.mesh` check; this MUST be preserved as an override in the new BoundaryFlux class, OR the space construction must encode the mesh identity (e.g. space name including a mesh hash).

**Detection**: `tests/transport/fields/test_boundary_flux.py::TestMeshBindingRejection::test_cross_mesh_add_rejected`.

**Why it hides**: same-mesh tests pass; cross-mesh is a programmer error that only fires in multi-mesh contexts (rare in unit tests, common in adjoint / sensitivity workflows).

**Catalog signature**: failure mode #6 (convention drift — Field's space-equality semantics differ from BoundaryFlux's mesh-binding semantics; the override must close the gap).

### Rank 5 — Interior cells still in BoundaryFlux (the split was incomplete)

**Failure mode**: the FaceLayout for the 2-D mesh allocates `total_size = N*ng*((nx+1)*ny + nx*(ny+1))` instead of `total_size = N*ng*(2*ny + 2*nx)` — i.e. the carve only renamed the buffers but did NOT split the interior cells out.

**Mechanism**: the `boundary_face_layout` property on SNMesh is implemented as a thin rename of the legacy `xmin_xmax_buf`/`ymin_ymax_buf` shapes. The interior cells live INSIDE BoundaryFlux because the FaceLayout includes them.

**Detection**: `tests/transport/fields/test_boundary_flux.py::TestFaceLayoutSliceViews::test_no_interior_cells_in_layout` AND `tests/sn/test_sweep_scratch_split.py::TestInteriorWavefrontCacheLocation::test_2d_boundary_flux_size_is_face_only`.

**Why it hides**: every sweep test passes (the buffers are just renamed; the math is identical). The carve appears complete but actually delivered NONE of the memory-churn benefit.

### Rank 6 — `dataclasses.replace` round-trip drops critical fields

**Failure mode**: the sweep's fresh-BoundaryFlux construction uses `dataclasses.replace(bf, values=new_flat)`. If `BoundaryFlux` has any non-default fields not explicitly passed to replace (e.g. `layout` or `mesh`), they could revert to defaults. The replace pattern is principled in plan §3.4 but error-prone.

**Mechanism**: dataclass replace COPIES all unmodified fields by default, so this should not happen. But if a subclass overrides `__init__` non-standardly, the replace can silently drop fields.

**Detection**: `tests/sn/test_boundary_flux_immutability_invariant.py::TestFunctionalConstructionPath::test_replace_creates_new_instance` (asserts mesh, layout preserved).

### Rank 7 — Field-ABC l2 norm uses unit-uniform metric, ignores FaceLayout

**Failure mode**: `bf.l2` calls `self.space.norm(self.values)`. The inherited Field uses `FunctionSpace.norm` which uses `inner_product_weights`. For BoundaryFlux on a flat buffer, the metric depends on which face is which — quadrature weights apply on the angular axis, area measure on the spatial axis. If the FunctionSpace's weights are NOT FaceLayout-aware, `bf.l2` returns a Euclidean norm rather than the L²-of-trace norm.

**Detection**: `tests/transport/fields/test_boundary_flux.py::TestFieldAlgebraInherited::test_l2_via_inherited_property`.

**Mitigation**: D-G may ship with a Euclidean fallback (the legacy code did not expose `.l2`); the convention is documented in the test docstring. If the FaceLayout-aware metric is in scope for D-G, the test pins the formula.

---

## 5. Bit-identity vs principled-equivalence judgment for D-G

**Verdict: D-G is PRINCIPLED, NOT bit-identical.** Confirm the plan §7.3 row's classification. The three criteria from `vv-principles` §"Bit-identity vs principled-equivalence":

### Criterion 1 — principled at every step

**SATISFIED.** Every intermediate is named:

- The **boundary trace** (`BoundaryFlux`) is a typed Field with `values + space + layout + mesh`. Each face is a named `BoundaryFaceFlux` accessible via `bf.faces[name]`.
- The **interior wavefront cache** (`SweepScratch`) is a named typed object owned by the sweep operator. Pre-D-G the cache lived in `xmin_xmax_buf[:, :, 1:nx, :]` — a slice of an unnamed-purpose ndarray. Post-D-G it lives as `sweep_scratch.x_cache[:, :, 1:nx, :]` — same data, named-intermediate location.
- The **flat-buffer arithmetic** is principled because `FaceLayout` is a named descriptor; `bf1 + bf2` is one numpy call on a contiguous buffer (one named reduction).

The refactor is a clear application of `coding-elegance` Pattern 3 (named intermediates with domain semantics) — boundary trace and interior cache are separated by ownership AND by name, where pre-D-G they shared a buffer with no naming separation.

### Criterion 2 — verified against structurally-independent reference

**SATISFIED, with caveats.** The reference structure:

- **L1 MMS gates** (`tests/sn/test_mms_aniso.py`, `tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py`) ARE structurally independent of the BoundaryFlux storage — the manufactured-source machinery in `orpheus.derivations.continuous.mms.sn` constructs Q^ext without touching boundary buffers. These are the load-bearing references for D-G.
- **k_inf homogeneous eigenvalue** (`tests/sn/l1_analytical/test_kinf_homogeneous.py`) — closed-form transfer-matrix reference, structurally independent.
- **Pre-D-G regression snapshots** — captured at the pre-D-G HEAD. They are NOT structurally independent of the pre-D-G implementation (they ARE the pre-D-G implementation's output), but their EXISTENCE establishes the "old value" against which the "new value" is compared at the ULP level. The L1 MMS gates establish that the pre-D-G output WAS correct in the first place; the snapshot match establishes that the post-D-G output reproduces it.

The combined evidence chain — L1 MMS proves the pre-D-G output is correct; snapshots prove the post-D-G output matches the pre-D-G output to ULP — closes the verification loop.

### Criterion 3 — drift dimensionally explainable

**SATISFIED.** The drift bound:

- **Within a single sweep**: reduction depth is the per-octant per-anti-diagonal cell update. The flat-buffer + slice-view storage produces the SAME reduction order as the pre-D-G named-buffer storage, because `psi_x[:, :, i, :]` (pre-D-G) and `flat_buf[xmin_offset + i*ny*ng*N : ...].reshape(N, ng, ny)` (post-D-G) point to the SAME memory layout if the FaceLayout offsets are constructed to match the legacy strided layout. **Expected drift per sweep: 0 ULP.** This is the strongest evidence — the cells ARE the same arithmetic on the same memory.
- **Across iterations (SI driver)**: each iteration constructs a fresh BoundaryFlux from outgoing face writes. The `dataclasses.replace` + reshape introduces one extra Python-level construction per iteration, but the underlying ndarray data is unchanged. **Expected drift per iteration: 0 ULP.**
- **WORST CASE bound**: `(reduction depth) × ULP × iter_count` per `vv-principles`. For typical SI runs (50-200 iterations, reduction depth ~10), the bound is `~1e-14 × iter_count = ~1e-12`. The regression-snapshot tolerance is `rtol=1e-13` (assert_array_equal where possible, allclose with strict tolerance where SI is involved).

**Concrete recommended regression tolerance**:

```python
# For 1-D snapshots where the sweep is a single call (no SI):
np.testing.assert_array_equal(out, expected)   # strict bit-identity

# For 2-D snapshots where SI runs to convergence:
np.testing.assert_allclose(out, expected, rtol=1e-13, atol=1e-14)
# Bound: iter_count × ULP ~ 100 × 2.22e-16 ~ 2.2e-14; rtol=1e-13 leaves a 5x safety margin.

# For converged 2-D reflective-BC cases:
np.testing.assert_allclose(out, expected, rtol=1e-12, atol=1e-13)
# Looser by 10x because reflective BC adds an inner partner-face read per iter,
# doubling the reduction depth.
```

### Anti-pattern check

Per `vv-principles` Anti-pattern #6 ("NEVER trust a reference that has not been traced back to a structurally-independent analytical or symbolic ground"):

- The reference chain is: L1 MMS gates (structurally independent, source-driven) → pre-D-G output (verified by MMS) → pre-D-G regression snapshots (frozen pre-D-G output) → post-D-G output (compared bit/ULP against pre-D-G snapshots).
- The chain terminates in MMS, which terminates in SymPy-derived source assembly (structurally independent of the production sweep). NOT reference-contaminated.

**Verdict CONFIRMED**: D-G is PRINCIPLED with `iter_count × ULP` drift bound; pre-D-G regression snapshots are the bit-identity gate at single-sweep granularity; L1 MMS gates are the load-bearing structurally-independent references.

---

## 6. Load-bearing L1 MMS gates for D-G specifically

These tests MUST stay green at every D-G sub-commit. They are the structurally-independent verification chain for the carve.

| Test file | Why load-bearing for D-G |
|---|---|
| `tests/sn/test_mms_aniso.py::test_sn_p1_aniso_mms_converges_second_order` | P1 anisotropic MMS on slab — exercises angular-coupling through reflective BC; load-bearing for 1-D BoundaryFlux refactor (xmin/xmax face slots). |
| `tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py` | Curvilinear MMS — exercises sphere/cylinder BoundaryFlux (xmax-only face). |
| `tests/sn/l1_analytical/test_kinf_homogeneous.py` | k_inf eigenvalue — closed-form transfer matrix; eigenvalue claim layer (not provable by MMS). |
| `tests/sn/l1_analytical/test_kinf_homogeneous_tolerance.py` | Tolerance pin on k_inf — catches regression in iteration-level numerics. |
| `tests/sn/regression/test_dd_regression.py` | Frozen DD snapshots (slab/sphere/cylinder, 1G/2G, homog/3-region) — flux-shape claims. |
| `tests/sn/test_2d_octant_sweep_equivalence.py` | 2-D octant equivalence — THE load-bearing test for the 2-D buffer split (since the conflation lives in 2-D). All 6 snapshots must stay green. |
| `tests/sn/test_boundary_flux_arithmetic.py` | Pre-D-G BoundaryFlux algebra tests — these will migrate to `tests/transport/fields/test_boundary_flux.py` during D-G but the pre-D-G versions must stay green until the migration commit. |
| `tests/sn/test_angular_flux_with_boundary.py` | AngularFlux carrying BoundaryFlux — D-G preserves the legacy `AngularFlux.boundary` field (D-H retires it). Must stay green. |
| `tests/sn/test_invertible_operator.py` | `(L+C).solve` consumes/produces AngularFlux with `.boundary` — exercises the full sweep with new BoundaryFlux. |

---

## 7. Tests that legitimately MIGRATE (not retire) during D-G

Per `feedback_retirement_means_test_migration`:

- `tests/sn/test_boundary_flux_arithmetic.py` — migrates to `tests/transport/fields/test_boundary_flux.py`. The MIGRATION re-writes the legacy mutable-write tests (`bf.xmin_face[...] = 1.0`) to the new flat-buffer functional path (`dataclasses.replace(bf, values=updated_flat)`). The migration commit retires the legacy test file AS the new tests land green.

- The legacy file's tests for `__add__`, `__sub__`, `__mul__`, etc. become foundation tests on the inherited Field algebra (the algebra is correct by inheritance from Field; the migrated tests verify the inheritance is wired).

- The pre-D-G test `test_add_2d_cartesian_propagates_to_both_buffers` (line 189 of the legacy file) migrates to `test_arithmetic_via_flat_buffer_2d` — the legacy "two named buffers" pattern is gone; the new test asserts flat-buffer arithmetic.

---

## 8. Decision matrix — which existing tests pin which legacy behavior

| Legacy behavior | Pinned by | Migration disposition |
|---|---|---|
| Mutable write-through (`bf.xmin_face[...] = ...`) | `tests/sn/test_boundary_flux_arithmetic.py` (writes to `xmin_face`, `xmax_face`) | RETIRE pre-D-G tests; replace with frozen-instance tests in `test_boundary_flux_immutability_invariant.py`. |
| Per-face named-buffer accessor (`bf.xmin_face`, `bf.xmax_face`) | `tests/sn/test_boundary_flux_arithmetic.py` (asserts on `out.xmin_face`) | MIGRATE — accessor becomes `bf.faces['xmin'].values`. New test asserts equivalent shapes/values. |
| 2-D combined buffer (`bf.xmin_xmax_buf`) | `tests/sn/test_boundary_flux_arithmetic.py::test_add_2d_cartesian_propagates_to_both_buffers` AND sweep tests | RETIRE the legacy attribute access; verify the data lives in flat buffer + SweepScratch via the SweepScratch split tests. |
| Geometry-conditional `zeros(mesh)` factory | `tests/sn/test_boundary_flux_arithmetic.py::_slab_mesh / _sphere_mesh / _cylinder_mesh fixtures` | MIGRATE — factory renames to `BoundaryFlux.zeros_for_sn_mesh(mesh)` (or stays as `zeros` — verify the new API in `test_boundary_flux.py::TestConstruction`). |
| Cross-mesh rejection | `tests/sn/test_boundary_flux_arithmetic.py::test_cross_mesh_add_rejected` | PRESERVE — copy verbatim into `tests/transport/fields/test_boundary_flux.py::TestMeshBindingRejection`. |
| Algebra dunder distributivity | `tests/sn/test_boundary_flux_arithmetic.py::test_distributive_property` | PRESERVE — copy into new test file; behavior must be unchanged. |
| SI-driver convergence with reflective BC | `tests/sn/test_2d_octant_sweep_equivalence.py::*_reflective_*` snapshots | PRESERVE — these are the load-bearing snapshots for the reflective BC partner read AT THE ITERATION LEVEL. |
| Reflective BC face read after sweep | `tests/sn/test_2d_octant_sweep_equivalence.py` AND DD-regression snapshots | PRESERVE + ADD targeted test `test_2d_reflective_bc_two_iterations_match_legacy` in `test_sweep_scratch_split.py`. |
| Persistent buffer survives across `(L+C).solve` calls | `tests/sn/test_invertible_operator.py` | PRESERVE — the InvertibleOperator's `_copy_boundary_face_state` helper at `operator.py:2680-2697` is the load-bearing seed-from-rhs path; D-G changes the helper's body but the test contract stays. |

---

## 9. Summary and exit criteria

D-G ships when:

1. ✅ Pre-D-G baseline pass/fail captured to `tests/sn/regression/_pre_dg_baseline_pass_fail.txt` (committed prep step).
2. ✅ `tests/transport/fields/test_boundary_flux.py` green at every D-G sub-commit.
3. ✅ `tests/sn/test_sweep_scratch_split.py` green at every D-G sub-commit.
4. ✅ `tests/sn/test_boundary_flux_immutability_invariant.py` green at every D-G sub-commit.
5. ✅ Re-runnable bit-identity gate (§2.2 diff) produces empty diff at every D-G sub-commit.
6. ✅ Targeted L1 MMS sub-gate (§2.3) green at every D-G sub-commit.
7. ✅ Pre-existing 10 DD-regression failures stay at the same failure set.
8. ✅ `tests/sn/test_boundary_flux_arithmetic.py` retired in the final D-G commit (the new test file supersedes it).

---

## 10. Open questions / decisions deferred to implementation

1. **Should `BoundaryFlux.zeros_for_sn_mesh` be named `zeros` (matching legacy) or `zeros_for_sn_mesh` (per plan §3.4)?** Plan says the latter; legacy says the former. The user has the call.

2. **Does the L²-of-trace metric on `bf.l2` matter for D-G, or defer to D-H/D-I?** The metric ambiguity (Euclidean vs FaceLayout-aware quadrature weighting) is real but does not block D-G. Recommend: D-G ships with Euclidean fallback documented in the test; promote to FaceLayout-aware metric in a follow-up if a consumer needs it.

3. **`InvertibleOperator._copy_boundary_face_state` at `operator.py:2680-2697` — does it survive D-G or get rewritten?** The helper currently does mutable write-through (`dst.xmin_face[...] = src.xmin_face`). Post-D-G immutable BoundaryFlux requires `dst = dataclasses.replace(dst, values=src.values)` or just `dst = src`. The helper becomes a one-line copy or vanishes. Recommend: vanish, replace call site with direct assignment.

4. **Should the migration test file path be `tests/transport/fields/test_boundary_flux.py` (per plan) or `tests/transport/fields/test_boundary_face_flux.py` for the BoundaryFaceFlux per-face wrapper?** Plan §3.4 distinguishes `BoundaryFlux` (over-all-faces) from `BoundaryFaceFlux` (one face). Recommend: single test file `test_boundary_flux.py` covering both; if `BoundaryFaceFlux` gets significant API surface, split later.

---

## References

- Depth B plan: `.claude/plans/depth_b_field_on_function_space.md` §3.4 (pure Field + FaceLayout), §6 step D-G, §7 verification, §11.1 invariants.
- Pre-D-G implementation: `orpheus/sn/boundary_flux.py` (lines 34-43 mutability rationale; lines 96-98 2-D buffer conflation).
- 2-D sweep: `orpheus/sn/sweep.py:697-858` (`_sweep_2d_wavefront`) — the BC apply at lines 819-831 and the persistent buffer scatter at lines 854-855.
- InvertibleOperator boundary copy: `orpheus/sn/operator.py:2680-2697`.
- L1 MMS gates: `tests/sn/test_mms_aniso.py`, `tests/sn/l1_analytical/`.
- Regression snapshots: `tests/sn/regression/snapshots/` (6 2-D + 6 1-D snapshots).
- `vv-principles` skill: claim taxonomy, three pillars, bit-identity vs principled-equivalence, ERR-006/026 curvilinear sweep bugs.
- `numerical-bug-signatures` skill: Signature 1 (curvilinear refinement), recognition catalog for reflective-BC partner reads.

Cross-link: `[[phase1_verification_plan]]`, `[[phase3_verification_plan]]`, `[[feedback_vv_tagging]]`, `[[feedback_diagnostic_promotion]]`.
