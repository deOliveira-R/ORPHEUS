r"""Piece 2 (#208) — SI convergence diagnostics on the iteration record.

The source-iteration loop measures the iterate increment
:math:`\Delta\psi = \psi^{(i)} - \psi^{(i-1)}` (space norm of the principal
bulk leaf) each pass and records the trajectory as
:attr:`~orpheus.numerics.convergence.IterationRecord.increment_norms`, from
which the record DERIVES the Banach contraction factor :math:`\rho \approx
\lVert\Delta\psi^{(i)}\rVert / \lVert\Delta\psi^{(i-1)}\rVert`
(``record.contraction_ratios``) and the :math:`c\to 1` geometric-tail
estimate (``record.true_error_estimate()``). Relocated off the retired typed
displacement surface at campaign 1 CS3 (2026-08-19). On a homogeneous slab the SI
contraction factor :math:`\rho \to c = \Sigma_s/\Sigma_t` (Adams–Larsen 2002),
so the recorded ratios track :math:`c`.

Marks ``l1`` — a convergence-rate claim, pillar = closed-form (:math:`\rho=c`).
1-group is acceptable HERE: this is a contraction-RATE diagnostic, NOT an
eigenvalue / flux-shape claim; :math:`\rho=c` is flux-shape-independent by
design (the ≥2G mandate guards eigenvalue/flux-shape claims, which this makes
neither). Builds the SI via the production single-source builders
(``build_within_group_system`` / ``_within_group_si``) — no construction duplicated.
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.geometry import BC, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.coupled_system import build_within_group_system
from orpheus.sn.solver import SNSolver, _within_group_si
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.source_sinks import AngularSourceSink, AngularBoundarySourceSink
from orpheus.transport.timed_full_field import TimedFullField


pytestmark = pytest.mark.l1


def _homogeneous_slab_solver(c: float, *, sigma_t: float = 1.0,
                             nx: int = 40, width: float = 10.0, n_ord: int = 8):
    r"""Homogeneous slab SNSolver: Σt=sigma_t, Σs=c·Σt, Σa=(1−c)·Σt, no fission.

    The within-group SI contraction factor is the spectral radius of
    :math:`(L+C)^{-1} S`, which → :math:`c` as the slab thickens.
    """
    sig_t = float(sigma_t)
    sig_s = c * sig_t
    z = np.zeros(1)
    mat = Mixture(
        SigC=np.array([(1.0 - c) * sig_t]),  # absorption: Σt = Σc + Σs (no fission)
        SigL=z.copy(), SigF=z.copy(), SigP=z.copy(),
        SigT=np.array([sig_t]),
        SigS=[csr_matrix(np.array([[sig_s]]))],  # P0 within-group scatter
        Sig2=csr_matrix((1, 1)), chi=z.copy(),
    )
    mesh = Mesh1D(
        edges=np.linspace(0.0, width, nx + 1), mat_ids=np.zeros(nx, dtype=int),
        bc_left=BC("vacuum"), bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=n_ord)
    sn_mesh = SNMesh(mesh, quad, {0: mat})
    return SNSolver(sn_mesh, inner_solver="source_iteration", scattering_order=0)


def _run_si(c: float, **kw):
    """Run the within-group SI on a homogeneous slab; return its
    :class:`IterationRecord` (with ``increment_norms`` populated)."""
    solver = _homogeneous_slab_solver(c, **kw)
    sn_mesh = solver.sn_mesh
    system = build_within_group_system(
        sn_mesh, solver.mat_xs, scattering_op=solver.scattering_op,
    )
    si, _base, _gains, windowed = _within_group_si(
        system, sn_mesh, inner_schedule=solver.inner_schedule,
        max_iter=600, tol=1e-12,
    )
    # 1-D slab never windows (windowing is 2-D Cartesian) → interior=AngularFlux.
    if windowed:
        raise AssertionError("homogeneous slab should not window the SI iterate")
    q_ext = TimedFullField(
        interior=AngularSourceSink.from_isotropic(
            np.full((sn_mesh.ng, *sn_mesh.spatial_shape), 1.0), sn_mesh,
        ),
        boundary=AngularBoundarySourceSink.zeros_on(sn_mesh),
        _history=(), history_depth=2,
    )
    ig = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, space=sn_mesh.full_field_space)
    _, record = si.solve(q_ext, initial_guess=ig)
    return record


def _asymptotic_rho(record) -> float:
    ratios = record.contraction_ratios
    if len(ratios) < 4:
        raise AssertionError(
            f"too few contraction ratios recorded ({len(ratios)}) — the SI loop "
            f"is not populating the record's increment-norm diagnostics."
        )
    return float(np.mean(ratios[-3:]))


@pytest.mark.parametrize("c,rho_lo,rho_hi", [(0.5, 0.40, 0.56), (0.9, 0.80, 0.92)])
def test_contraction_ratio_tracks_scattering_ratio(c, rho_lo, rho_hi):
    r"""``IterationRecord.contraction_ratios`` → ρ ≈ c on a homogeneous slab.

    The record's increment-norm trajectory is the ONLY surface that knows
    "previous"/"step". Proves the derived diagnostic tracks the analytical
    Banach factor ρ ≈ max(Σs/Σt) = c — turning a ρ-blind ‖Δψ‖ reading honest.
    Rate claim (1-group acceptable)."""
    rho = _asymptotic_rho(_run_si(c))
    if not (rho_lo <= rho <= rho_hi):
        raise AssertionError(f"c={c}: measured ρ={rho:.4f} not in [{rho_lo}, {rho_hi}]")


def test_contraction_ratio_increases_with_c():
    r"""ρ(c=0.9) > ρ(c=0.5) by a clear margin — the ratio TRACKS c (the
    discriminating cross-check that ρ is the contraction factor, not a constant)."""
    rho_lo = _asymptotic_rho(_run_si(0.5))
    rho_hi = _asymptotic_rho(_run_si(0.9))
    if not (rho_hi > rho_lo + 0.2):
        raise AssertionError(
            f"ρ did not track c: ρ(0.5)={rho_lo:.4f}, ρ(0.9)={rho_hi:.4f}")


def test_record_feeds_true_error_estimate():
    r"""``record.true_error_estimate()`` == ‖Δψ‖/(1−ρ) of the SAME solve —
    the c→1 false-convergence fix, derived from the recorded trajectory.
    (The per-entry convergence map is now
    :meth:`~orpheus.numerics.field.Field.where_largest`, a property of any
    field — unit-tested in :mod:`tests.numerics.test_field`.)"""
    record = _run_si(0.9)
    if not record.increment_norms:
        raise AssertionError("increment_norms not recorded")
    rho = record.contraction_ratios[-1]
    np.testing.assert_allclose(
        record.true_error_estimate(),
        record.increment_norms[-1] / (1.0 - rho), rtol=1e-12,
    )
