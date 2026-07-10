r"""#206 Phase 0 — the A-NEW pre-carve baseline gate (the SHARP Phase-A gate).

This is the one gate the #206 verification spec flagged as "not yet
existing" (test-architect memo
``.claude/agent-memory/test-architect/issue_206_cellupdate_seam_verification.md``,
block "New-gate recommendation A-NEW"). The existing
``test_bc_extraction_matvec.py::TestVacuumMatvecBitIdentity`` pins the
matvec on a random ψ but at **vacuum** BC and through the ``(L+C)`` public
surface; it does NOT pin a single raw **sweep** output, and it runs on the
degenerate vacuum/flat-redistribution regime.  A-NEW is sharper on three
axes that the #206 carve actually touches:

1. **Heterogeneous σ_t** (per-group, per-cell random) — drives the
   curvilinear redistribution denom term ``dA_w · c_out`` out of
   cancellation.  Uniform σ_t + flat ψ NULL it (vv §H2), so a flat-ψ /
   uniform-σ_t gate would pass GREEN while a real curvilinear closure bug
   hides.  Phase A's whole point is routing BOTH the sweep and the matvec
   through the SAME ``scheme`` cell math (#206 Phase A: the closure
   methods; #158 Increment B refactored group-3 to the coefficient model —
   ``affine_scan_coefficients`` → ``(a, inverse_denom, w)`` consumed by the
   generic base reconstruction staticmethods on ``DiscretizationSchemeBase``, the
   per-scheme closure methods retired); the redistribution path is exactly
   where a routing slip would show.

2. **Non-flat random ψ / Q** (≥2 groups) — fixed-seed random bulk +
   boundary trace (matvec) and a fixed-seed random per-ordinate source
   (sweep).  Flat ψ telescopes per-ordinate balance by construction
   (vv §H3) and nulls redistribution (vv §H2); the random non-flat field
   is the discriminating stressor.  1-group is degenerate for eigenvalue
   claims (the 1-group-degeneracy rule) — here the claim is a raw single-application
   bit-identity, but ≥2G is kept regardless because group coupling is one
   of the convention axes the carve crosses (L17 crosswalk).

3. **A bare single-SWEEP output** (``angular_flux`` AND ``scalar_flux``),
   captured directly from ``transport_sweep`` — the production sweep entry
   the carve relocates.  This is the leg the existing matvec-only gate
   does not cover.

For slab + sphere + cylinder, ≥2 groups, the gate captures (under
``--capture-baseline``, run ONCE BEFORE Phase A) and otherwise asserts:

* **Sweep leg:** ONE ``transport_sweep`` driven by a fixed-seed random
  per-ordinate source ``Q`` + heterogeneous σ_t — both the returned
  ``angular_flux`` ``(N, ng, nx)`` and ``scalar_flux`` ``(ng, nx)``.
* **Matvec leg:** ONE ``(L + C).apply(ψ).interior`` on a fixed-seed random
  NON-flat ψ (random bulk + random boundary trace) + heterogeneous σ_t
  — the canonical :func:`tests.sn._test_helpers.het_operands` operands,
  which are tuned precisely so every loss-action term is activated.

Tolerance is the single-source-of-truth
:func:`tests.sn.regression._regression_assert.assert_regression` with
``kind="direct"``: a single sweep / single matvec has no outer iteration,
so the only drift is FP-non-associativity over the reduction tree, bounded
by ``reduction_depth × ULP``.  Phase A/B (pure relocation + closure
single-sourcing) should be **bit-identical** — the
:class:`~tests.sn.regression._regression_assert.DriftWarning` tripwire
fires on ANY ULP move, and ``pytest -W error::DriftWarning`` escalates it
to a hard strict-bit-identity gate.  If Phase A genuinely reorders the
reduction tree (the density-form ``(denom·ψ − numer)/V`` vs scan-form
``b = 2·QV·inv_denom`` regrouping the audit flagged as principled-
equivalent at nULP), ``kind="direct"`` admits the bounded ULP drift while
still hard-failing an algorithmic change masquerading as FP noise.

``-O`` safety
-------------
Every assertion routes through ``assert_regression`` (``np.testing.*`` /
explicit ``raise``) or ``pytest.fail`` — NO bare ``assert`` carries the
load.  The gate fires under the canonical ``python -O -m pytest``
(vv Mode 8).  The two ``assert path.exists(...)`` lines are
collection-time guards (a missing-baseline operator error), NOT the
numerical gate; even were they stripped under ``-O``, ``np.load`` on a
missing file raises ``FileNotFoundError`` — the gate cannot vacuously
pass.

Tagging
-------
``@pytest.mark.foundation`` — software invariant ("the carve did not move
a bit / moved only bounded FP-non-associativity ULP").  No theory-page
``:label:``; the claim is a relocation-fidelity contract, not an L0/L1
equation claim.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.loss_representation import transport_sweep
from orpheus.sn.operators.streaming import StreamingOperator
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.source_sinks import AngularSourceSink
from tests.sn._test_helpers import (
    SN_TESTS_ROOT,
    het_operands,
    placeholder_materials,
)
from tests.sn.regression._regression_assert import assert_regression


pytestmark = [pytest.mark.regression, pytest.mark.foundation]

# ─────────────────────────────────────────────────────────────────────
# Baseline snapshot store + the --capture-baseline flag (ROOT conftest).
# ─────────────────────────────────────────────────────────────────────
#
# ``--capture-baseline`` is declared in the ROOT ``tests/conftest.py``
# (``pytest_addoption`` only fires there). When present the gate WRITES
# the pre-carve sweep/matvec output and skips the assert; absent (normal
# runs, incl. the post-carve gate) it READS the committed snapshots and
# asserts principled agreement via ``assert_regression``.

_BASELINE_DIR = SN_TESTS_ROOT / "_data" / "affine_carve_baseline"

# ≥2 groups (the 1-group-degeneracy rule: 1G is degenerate). The matvec/sweep
# bit-identity is size-independent; small N keeps each run sub-second.
_NG = 2
_N_CELLS = 8
_N_ORD = 4

# slab / sphere / cylinder — the three 1-D geometries the carve relocates.
_GEOMS_1D = ("SLB", "SPH", "CYL")


def _capturing(request) -> bool:
    return bool(request.config.getoption("--capture-baseline", default=False))


def _build_sn_mesh(geometry: str) -> SNMesh:
    """Small ≥2G SNMesh for slab / sphere / cylinder.

    Outer BC vacuum; the curvilinear inner edge (r=0) is the regularity
    pole, declared ``reflective`` per the project convention (NOT a BC —
    a regularity condition).  Cylinder uses ``level_symmetric`` (the
    curvilinear cubature the matvec serves); slab/sphere use
    ``gauss_legendre`` — mirrors ``test_bc_extraction_matvec._build_sn_mesh``
    so the A-NEW operands match the existing gate's geometry exactly.
    """
    mats = placeholder_materials(ng=_NG)
    if geometry == "SLB":
        mesh = Mesh1D(
            edges=np.linspace(0.0, 2.0, _N_CELLS + 1),
            mat_ids=np.zeros(_N_CELLS, dtype=int),
            coord=CoordSystem.CARTESIAN,
            bc_left=BC("vacuum"),
            bc_right=BC("vacuum"),
        )
        quad = Quadrature.gauss_legendre(n_ordinates=_N_ORD)
    elif geometry == "SPH":
        mesh = Mesh1D(
            edges=np.linspace(0.0, 2.0, _N_CELLS + 1),
            mat_ids=np.zeros(_N_CELLS, dtype=int),
            coord=CoordSystem.SPHERICAL,
            bc_left=BC("reflective"),  # r=0 pole — regularity, not a BC
            bc_right=BC("vacuum"),
        )
        quad = Quadrature.gauss_legendre(n_ordinates=_N_ORD)
    elif geometry == "CYL":
        mesh = Mesh1D(
            edges=np.linspace(0.01, 2.0, _N_CELLS + 1),
            mat_ids=np.zeros(_N_CELLS, dtype=int),
            coord=CoordSystem.CYLINDRICAL,
            bc_left=BC("reflective"),
            bc_right=BC("vacuum"),
        )
        quad = Quadrature.level_symmetric(sn_order=_N_ORD)
    else:  # pragma: no cover - guarded by parametrize
        raise ValueError(geometry)
    return SNMesh(mesh, quad, mats)


def _capture_or_assert(
    request,
    key: str,
    actual: np.ndarray,
    *,
    reduction_depth: int,
    quantity: str,
) -> bool:
    """WRITE the baseline under ``--capture-baseline``, else READ + assert.

    Returns ``True`` if it WROTE (capture mode), ``False`` if it READ +
    asserted.  The caller skips ONCE at the end of the test when any leg
    was captured — calling ``pytest.skip`` inside this helper would
    short-circuit the test body and drop the remaining legs (the bug a
    multi-array test must avoid).

    The principled gate is ``assert_regression(kind="direct", ...)`` — a
    single application has no outer iteration, so the only legitimate
    drift is FP-non-associativity over the reduction tree, bounded by
    ``reduction_depth × ULP``.  The :class:`DriftWarning` tripwire fires
    on ANY ULP move; ``-W error::DriftWarning`` makes the gate strict.
    """
    path = _BASELINE_DIR / f"{key}.npy"
    if _capturing(request):
        _BASELINE_DIR.mkdir(parents=True, exist_ok=True)
        np.save(path, actual)
        return True
    if not path.exists():
        pytest.fail(
            f"missing baseline snapshot {path}; run with --capture-baseline "
            f"BEFORE Phase A (at the pre-carve HEAD) to write it."
        )
    expected = np.load(path)
    assert_regression(
        actual,
        expected,
        conv_tol=0.0,  # ignored for kind="direct"
        case_name=key,
        kind="direct",
        reduction_depth=reduction_depth,
        quantity=quantity,
    )
    return False


# ═════════════════════════════════════════════════════════════════════
# SWEEP LEG — one raw transport_sweep on a fixed-seed random Q + het σ_t.
# ═════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
class TestAffineCarveSweepBaseline:
    """One raw ``transport_sweep`` output (angular + scalar) is unmoved.

    Catches: a Phase-A closure-routing slip or a Phase-B relocation that
    perturbs the raw forward sweep on a heterogeneous, non-flat
    configuration — the regime where the curvilinear redistribution term
    is ACTIVE (vv §H2), which the existing vacuum/flat matvec gate cannot
    see.
    """

    @pytest.mark.parametrize("geometry", _GEOMS_1D)
    def test_sweep_angular_and_scalar_unmoved(self, geometry, request):
        """[foundation] One transport_sweep: angular_flux AND scalar_flux.

        Fixed-seed random per-ordinate source ``Q`` (≥2G) + heterogeneous
        σ_t.  Both returned arrays are pinned: ``angular_flux`` reduces
        over ``nx`` cells (the forward scan) → ``reduction_depth = nx``;
        ``scalar_flux`` adds the ``N``-ordinate weight sum on top of the
        ``nx`` scan → ``reduction_depth = nx + N``.
        """
        sn_mesh = _build_sn_mesh(geometry)
        N = sn_mesh.quad.N
        ng = sn_mesh.ng
        nx = sn_mesh.nx

        # Heterogeneous σ_t (per-group, per-cell) and a fixed-seed random
        # per-ordinate source.  Distinct seed from the matvec leg so the
        # two legs are independent stressors.
        rng = np.random.default_rng(20260614)
        sig_t = rng.uniform(0.3, 3.0, size=(ng, *sn_mesh.spatial_shape))
        q_values = rng.standard_normal((N, ng, *sn_mesh.spatial_shape))
        source = AngularSourceSink.from_mesh(q_values, sn_mesh)
        boundary = AngularBoundaryFlux.zeros_on(sn_mesh)

        angular_flux, scalar_flux = transport_sweep(
            source, sig_t, sn_mesh, boundary,
        )

        captured = _capture_or_assert(
            request,
            f"sweep_angular_{geometry}",
            np.asarray(angular_flux, dtype=np.float64),
            reduction_depth=nx,
            quantity="angular_flux",
        )
        captured |= _capture_or_assert(
            request,
            f"sweep_scalar_{geometry}",
            np.asarray(scalar_flux, dtype=np.float64),
            reduction_depth=nx + N,
            quantity="scalar_flux",
        )
        if captured:
            pytest.skip(f"captured sweep baselines for {geometry}")


# ═════════════════════════════════════════════════════════════════════
# MATVEC LEG — one (L+C).apply(ψ).interior on a fixed-seed random non-flat ψ.
# ═════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
class TestAffineCarveMatvecBaseline:
    """One ``(L+C).apply(ψ).interior`` on a het σ_t + non-flat random ψ is
    unmoved.

    Sharper than ``TestVacuumMatvecBitIdentity`` (vacuum BC, flat-
    redistribution regime): the canonical
    :func:`tests.sn._test_helpers.het_operands` builds heterogeneous σ_t
    AND a random bulk + random boundary trace, so every loss-action term
    — including the curvilinear ``dA_w · c_out`` redistribution path the
    carve routes through ``scheme`` — is ACTIVE.  A flat-ψ gate
    would pass while a redistribution-routing bug hides (vv §H2).
    """

    @pytest.mark.parametrize("geometry", _GEOMS_1D)
    def test_matvec_bulk_unmoved(self, geometry, request):
        """[foundation] One (L+C).apply(ψ).interior on het σ_t + non-flat ψ.

        The bulk residual ``(N, ng, nx)`` reduces over ``nx`` cells (the
        per-cell WDD recurrence ``psi_face_in = 2·psi_cell − psi_face_in``)
        → ``reduction_depth = nx``.
        """
        sn_mesh = _build_sn_mesh(geometry)
        nx = sn_mesh.nx

        # het_operands: heterogeneous σ_t (fixed seed) + non-flat random
        # bulk AND boundary trace (every term activated).  A reflective
        # outer BC drives the boundary trace into the matvec (vs vacuum).
        sig_t, psi = het_operands(sn_mesh)
        L = StreamingOperator(sn_mesh)
        C = MultiplicationOperator.from_mesh(sig_t, sn_mesh)
        out = (L + C).apply(psi)

        captured = _capture_or_assert(
            request,
            f"matvec_bulk_{geometry}",
            np.asarray(out.interior.values, dtype=np.float64),
            reduction_depth=nx,
            quantity="matvec_bulk",
        )
        if captured:
            pytest.skip(f"captured matvec baseline for {geometry}")
