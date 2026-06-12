r"""Phase 4 / O.2b 4.6 — prescribed-inflow splitting-invariance (vv Mode 9).

**The claim (a software invariant, NOT an equation/eigenvalue claim).**
A non-zero prescribed inflow injected as the boundary *source* slot
``q.boundary`` is honoured **identically** by SI-Jacobi, SI-Gauss-Seidel,
and Krylov.  These three are different *operator splittings* of the same
affine within-group fixed point

.. math::

   (L + C - S - B)\,\psi = q,
   \qquad q = q_{\rm ext} + (\text{prescribed inflow in } q.\text{boundary}),

so they MUST reach the same :math:`\psi` (``vv-principles`` Mode 9 —
verify operator splittings reach the same fixed point under anisotropic
/ B≠0 stressing).  This is the *floor* the 4.6 non-vacuum MMS rows
(``test_mms_prescribed_inflow.py``) trust before they accept the SI
default path on a non-vacuum inflow.

**Why a foundation test (not L1).**  No theory-page ``:label:`` is being
verified — this pins that three reduction trees of one affine operator
agree on one RHS.  Per the V&V level taxonomy, ``foundation`` NEVER
carries ``verifies()``.

**Promotes the design probe** ``derivations/diagnostics/
diag_p46_prescribed_inflow_source.py`` (the affine-BC-``q``→RHS
hypothesis, verified this session ≤5.6e-13).  Once this test is green
the probe MAY be deleted (PROMOTE then DELETE,
``tests/derivations/_promotion_policy.md``).

**Anti-latent-dud preconditions (the G-4-dud lesson, ``si-gauss-seidel-
rate-recovery-verification-spec``).**  A consistency test where all three
trivially agree on the *zero* solution proves nothing.  The
precondition-asserts make the test self-falsifying if the config
degenerates: the inflow slot must actually be written (>0), the inflow
must non-trivially drive the interior (max|ψ|>1e-3), and the 2-D row
must actually run the B-folding ``_GaussSeidelResolvent`` (not silently
fall back to Jacobi) with an EXPLICIT reflective-y ``Mesh2D`` BC (not a
string kwarg an explicit mesh would override).

See:
- ``.claude/skills/vv-principles/SKILL.md`` (Mode 9 — splitting invariance).
- :func:`orpheus.sn.solver._within_group_triple` / ``_select_si_resolvent``
  / ``_within_group_krylov`` — the operator triple + splittings under test.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.mms.sn import _make_1g_mixture
from orpheus.geometry import BC, Mesh1D, Mesh2D
from orpheus.geometry.coord import CoordSystem
from orpheus.numerics.iteration import SourceIteration
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.geometry import SNMesh
from orpheus.sn.solver import (
    SNSolver,
    _select_si_resolvent,
    _within_group_krylov,
    _within_group_triple,
)
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.boundary_flux import BoundaryFlux
from orpheus.transport.source_sinks import AngularSourceSink, BoundarySourceSink
from orpheus.transport.timed_full_field import TimedFullField


# ── shared helpers (the probe's primitives) ──────────────────────────


def _flux_zero(sn: SNMesh) -> TimedFullField:
    """The iterate template — a FLUX composite (NOT the source space).

    B.5.2 (``iteration.py:685``): the iterate ψ + the returned solution
    live in the SOLUTION (flux) space; ``q_ext`` lives in SOURCE space.
    SI/Krylov need this flux ``initial_guess`` to template the matvec /
    unravel, else ``S.apply`` hits ``AngularSourceSink`` (no
    ``integrate_angular``).
    """
    return TimedFullField.zeros(bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn)


def _prescribed_inflow_source(
    sn: SNMesh, *, psi_in: float, face: str,
) -> TimedFullField:
    """``q_ext`` = composite with the prescribed (isotropic) inflow on ``face``.

    The inflow ordinate slots of ``face`` carry ``psi_in``; every other
    slot is zero (the affine-BC inhomogeneous term ``q`` pushed to the
    RHS).  Bulk source is zero — the inflow alone drives the interior.
    """
    N, ng = sn.quad.N, sn.ng
    bulk_vals = np.zeros((N, ng, *sn.spatial_shape))
    bss = BoundarySourceSink.zeros_on(sn)
    inflow = sn.trace.inflow_indices_for_face(face)
    bss.face_view(face)[inflow] = psi_in
    return TimedFullField(
        bulk=AngularSourceSink.from_mesh(bulk_vals, sn),
        boundary=bss,
        _history=(),
        history_depth=2,
    )


def _reldiff(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a)
    b = np.asarray(b)
    return float(np.abs(a - b).max() / (np.abs(a).max() + np.abs(b).max() + 1e-300))


def _require(condition: bool, message: str) -> None:
    """``-O``-survivable assertion.

    The canonical ORPHEUS test invocation is ``python -O`` (bare
    ``assert`` statements are stripped — a canary that cannot die).  This
    test's load-bearing checks therefore raise explicitly instead of
    via ``assert`` so they FAIL under ``-O`` as intended.
    """
    if not condition:
        raise AssertionError(message)


# ── the test (parametrised over the two discriminating configs) ──────


@pytest.mark.foundation
@pytest.mark.parametrize(
    "config",
    [
        "slab_1d",
        pytest.param("cart2d_reflective_y", marks=pytest.mark.slow),
    ],
)
def test_prescribed_inflow_consistency_si_jacobi_gs_krylov(config: str):
    r"""vv Mode 9 — a non-zero ``q.boundary`` is honoured IDENTICALLY by
    SI-Jacobi, SI-Gauss-Seidel, and Krylov.

    ``slab_1d``: 1-D slab, vacuum both faces + prescribed inflow on xmin.
    SI is always Jacobi in 1-D, so the discriminating pair is **SI ≡
    Krylov** (≤ ~1e-13).  Catches a Krylov matvec that drops
    ``q.boundary`` or an SI that mis-seeds the inflow.

    ``cart2d_reflective_y``: 2-D Cartesian, vacuum-x + REFLECTIVE-y (B≠0)
    + prescribed inflow on xmin.  B≠0 is what makes **SI-Jacobi vs
    SI-Gauss-Seidel** distinct (G-S folds B into ``_GaussSeidelResolvent``;
    Jacobi keeps B as a lagged gain).  Three-way ≡ (J ≡ GS ≡ K, ≤ ~5.6e-13).
    The B≠0 + prescribed-inflow combination is the only config where the
    G-S-folds-B path runs WITH a non-zero boundary source — the ERR-056
    family neighbourhood.

    Pairwise reldiff ceiling 1e-11 leaves headroom for the FP-non-
    associativity of three different reduction trees (bounded by
    ``iter_count × ULP`` per ``vv-principles`` bit-identity §3; the probe
    measured 1.3e-13 … 5.6e-13).
    """
    psi_in = 1.0
    sigma_t, sigma_s = 1.0, 0.5

    if config == "slab_1d":
        quad = Quadrature.gauss_legendre(n_ordinates=8)
        mesh = Mesh1D(
            edges=np.linspace(0.0, 2.0, 21),
            mat_ids=np.zeros(20, dtype=int),
            coord=CoordSystem.CARTESIAN,
            bc_left=BC("vacuum"),
            bc_right=BC("vacuum"),
        )
        run_gs = False
    elif config == "cart2d_reflective_y":
        quad = Quadrature.level_symmetric(4)
        n = 8
        # EXPLICIT reflective-y on the Mesh2D constructor (NOT a
        # boundary_condition='reflective' string kwarg an explicit mesh
        # BC would silently override — the G-4 latent-dud parallel).
        mesh = Mesh2D(
            edges_x=np.linspace(0.0, 2.0, n + 1),
            edges_y=np.linspace(0.0, 2.0, n + 1),
            mat_map=np.zeros((n, n), dtype=int),
            coord=CoordSystem.CARTESIAN,
            bc_xmin=BC("vacuum"),
            bc_xmax=BC("vacuum"),
            bc_ymin=BC("reflective"),  # B≠0 → G-S folds it
            bc_ymax=BC("reflective"),
        )
        run_gs = True
    else:  # pragma: no cover - guarded by parametrize
        raise AssertionError(config)

    sn = SNMesh(mesh, quad, {0: _make_1g_mixture(sigma_t, sigma_s)})
    solver = SNSolver(sn, max_inner=500, inner_tol=1e-13)
    LC, S, B = _within_group_triple(solver)
    n_dof = quad.N * sn.ng * int(np.prod(sn.spatial_shape)) + int(sn.trace.layout.total_size)

    face = "xmin"
    q_ext = _prescribed_inflow_source(sn, psi_in=psi_in, face=face)
    inflow = sn.trace.inflow_indices_for_face(face)

    # PRECONDITION 1 — the inflow slot is actually written (guards the
    # silently-vacuum dud where all three trivially agree on ψ≡0).
    _require(psi_in != 0.0, "psi_in degenerated to zero")
    _require(
        np.abs(q_ext.boundary.face_view(face)[inflow]).max() > 0.0,
        "prescribed inflow slot is empty — q.boundary degenerated to vacuum",
    )

    # SI-Jacobi: resolvent (L+C), gains (S, B).
    rJ, gJ = _select_si_resolvent(LC, S, B, sn, "jacobi")
    psi_j, _ = SourceIteration(rJ, *gJ, max_iter=500, tol=1e-13).solve(
        q_ext, initial_guess=_flux_zero(sn),
    )

    # Krylov: matvec (L+C−S−B), identity precond, full restart.
    kr = _within_group_krylov(LC, S, B, n_dof=n_dof, max_iter=500, tol=1e-13)
    psi_k, _ = kr.solve(q_ext, initial_guess=_flux_zero(sn))

    # PRECONDITION 2 — the inflow non-trivially drives the interior
    # (else J≡GS≡K could hold vacuously at ψ≡0).
    _require(
        np.abs(psi_j.bulk.values).max() > 1e-3,
        "prescribed inflow did not drive the interior — vacuous consistency",
    )

    # SI ≡ Krylov (the 1-D discriminating pair; also a 2-D leg).
    _require(
        _reldiff(psi_j.bulk.values, psi_k.bulk.values) < 1e-11,
        f"SI-Jacobi vs Krylov reldiff = "
        f"{_reldiff(psi_j.bulk.values, psi_k.bulk.values):.2e} >= 1e-11",
    )

    if run_gs:
        rG, gG = _select_si_resolvent(LC, S, B, sn, "gauss_seidel")
        # PRECONDITION 3 — the G-S path is the B-folding resolvent, not a
        # silent Jacobi fall-back (guards the _select_si_resolvent dispatch).
        _require(
            type(rG).__name__ == "_GaussSeidelResolvent",
            f"gauss_seidel did not select the B-folding resolvent: got "
            f"{type(rG).__name__}",
        )
        # ... and reflective-y MUST yield a non-trivial B (else the
        # J-vs-GS distinction collapses to a no-op).
        from orpheus.numerics.operator import CAP_APPLY

        _require(B is not None, "B (SNBoundaryOperator) is None")
        _require(CAP_APPLY in B.capabilities, "B does not advertise CAP_APPLY")

        psi_g, _ = SourceIteration(rG, *gG, max_iter=500, tol=1e-13).solve(
            q_ext, initial_guess=_flux_zero(sn),
        )

        # Three-way ≡ (J ≡ GS ≡ K).
        _require(
            _reldiff(psi_j.bulk.values, psi_g.bulk.values) < 1e-11,
            f"SI-Jacobi vs SI-Gauss-Seidel reldiff = "
            f"{_reldiff(psi_j.bulk.values, psi_g.bulk.values):.2e} >= 1e-11",
        )
        _require(
            _reldiff(psi_g.bulk.values, psi_k.bulk.values) < 1e-11,
            f"SI-Gauss-Seidel vs Krylov reldiff = "
            f"{_reldiff(psi_g.bulk.values, psi_k.bulk.values):.2e} >= 1e-11",
        )

        # The inflow is honoured identically across all three.
        got_j = psi_j.boundary.face_view(face)[inflow]
        got_g = psi_g.boundary.face_view(face)[inflow]
        got_k = psi_k.boundary.face_view(face)[inflow]
        np.testing.assert_allclose(got_j, got_g, rtol=1e-11, atol=1e-13)
        np.testing.assert_allclose(got_j, got_k, rtol=1e-11, atol=1e-13)
