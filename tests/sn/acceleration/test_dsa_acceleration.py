r"""The 3b acceleration gates — FP-invariance (D3/D4), correction→0 (D6), teeth.

The correction→0 partition (verification spec §0): the accelerator's
machinery is value-safe BY CONSTRUCTION (its correction vanishes with
the displacement), so FP-invariance gates prove the SAFETY property
while the object gates + rate tier carry the machinery's correctness
weight. Mode-9 discipline: the invariance configs BREAK the degenerate
coincidences — heterogeneous, ≥2 groups, ANISOTROPIC (ℓ ≥ 1, the
mixtures' real P1 data active via ``scattering_order=1``), vacuum AND
reflective walls — never the isotropic-reflective box.

Teeth (in-process monkeypatch, never a file mutation):

* a SIGN-FLIPPED correction must break the D3 convergence witness —
  the FP gate is constrained, not decorative;
* a trace-arm-ZEROED correction must break the REFLECTIVE case — the
  regression pin of the 3b root-cause: the lagged reflective gain must
  read the corrected outflow (Larsen's reflecting row models a
  current-inflow error equation), and the wall-edge trace arm is the
  load-bearing carrier of that consistency (Mode-10: the arm is
  exercised AND constrained).

The rate tier (ρ vs A&L (3.65), S2 exactness, thickness/c sweeps) is
3c's D11–D13; here the rate appears only as the coarse
``n_inner(dsa) < n_inner(plain)`` acceleration witness.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC
from orpheus.geometry.mesh import Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.acceleration.dsa import DSACorrection
from orpheus.sn.solver import Solution, solve_sn_fixed_source

pytestmark = pytest.mark.l2


_TOL = 1e-11
#: FP-identity band: both runs stop on the same equation-residual tol;
#: the iterates agree to the un-accelerated contraction's error constant
#: (·1/(1−ρ) ≈ 40 at c ≈ 0.975) — 1e-7 is ~3 decades of headroom above
#: that floor and ~4 below any accelerator-off-by-one signal.
_FP_RTOL = 1e-7


def _solve(
    bc: tuple[str, str],
    *,
    acceleration: str | None = None,
    inner_solver: str | None = None,
    max_inner: int = 3000,
) -> Solution:
    """Het 2-zone, 2G, ℓ ≥ 1 slab, S4 — the Mode-9-honest config."""
    mats = {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")}
    mesh = Mesh1D(
        edges=np.linspace(0.0, 20.0, 41),
        mat_ids=np.array([0] * 20 + [1] * 20),
        bc_left=BC(bc[0]),
        bc_right=BC(bc[1]),
    )
    return solve_sn_fixed_source(
        materials=mats,
        mesh=mesh,
        quadrature=Quadrature.gauss_legendre(n_ordinates=4),
        external_source=np.ones((4, 2, 40)),
        scattering_order=1,  # the mixtures' real P1 rows active (ℓ ≥ 1)
        inner_tol=_TOL,
        max_inner=max_inner,
        inner_solver=inner_solver,
        acceleration=acceleration,
    )


def _phi(solution: Solution) -> np.ndarray:
    v = solution.scalar_flux
    return np.asarray(getattr(v, "values", v), dtype=float)


def _inners(solution: Solution) -> int:
    """The iteration count, loudly narrowed (history is Optional on the
    Solution contract; the fixed-source paths always populate it)."""
    history = solution.history
    if history is None or history.n_inner is None:
        pytest.fail("the fixed-source Solution must carry n_inner history")
    return int(history.n_inner)


class TestD3FixedPointInvarianceSI:
    """D3 — SI+DSA converges to the PLAIN-SI fixed point (anisotropic,
    heterogeneous, 2G; vacuum and reflective)."""

    @pytest.mark.parametrize(
        "bc", [("vacuum", "vacuum"), ("reflective", "vacuum")]
    )
    def test_si_dsa_shares_the_unaccelerated_fixed_point(self, bc):
        plain = _solve(bc)
        dsa = _solve(bc, acceleration="dsa")
        np.testing.assert_allclose(
            _phi(dsa), _phi(plain), rtol=_FP_RTOL, atol=0,
            err_msg="SI+DSA must share the plain-SI fixed point",
        )
        if not _inners(dsa) < _inners(plain):
            pytest.fail(
                f"DSA must accelerate ({bc}): "
                f"{_inners(dsa)} !< {_inners(plain)}"
            )


class TestD4FixedPointInvarianceKrylov:
    """D4 — Krylov-preconditioned-by-DSA converges to the
    unpreconditioned GMRES fixed point (safety; effectiveness = D13)."""

    @pytest.mark.parametrize(
        "bc", [("vacuum", "vacuum"), ("reflective", "vacuum")]
    )
    def test_krylov_dsa_shares_the_unpreconditioned_fixed_point(self, bc):
        plain = _solve(bc, inner_solver="krylov")
        dsa = _solve(bc, inner_solver="krylov", acceleration="dsa")
        np.testing.assert_allclose(
            _phi(dsa), _phi(plain), rtol=_FP_RTOL, atol=0,
            err_msg="Krylov+DSA must share the unpreconditioned FP",
        )
        if not _inners(dsa) < _inners(plain):
            pytest.fail(
                f"the DSA preconditioner must reduce GMRES iterations "
                f"({bc}): {_inners(dsa)} !< {_inners(plain)}"
            )


class TestD6CorrectionVanishes:
    """D6 — the correction→0 safety property at the matrix level:
    a zero displacement maps to an exactly-zero correction (G·0 = 0,
    A⁻¹·0 = 0, P·0 = 0 — exact arithmetic, no tolerance)."""

    @pytest.mark.verifies("sn-dsa-correction-vanishes")
    def test_zero_displacement_maps_to_exact_zero(self):
        from orpheus.sn.mesh.augmented_mesh import SNMesh
        from orpheus.transport.fields.angular_flux import AngularFlux
        from orpheus.transport.fields.angular_boundary_flux import (
            AngularBoundaryFlux,
        )
        from orpheus.transport.timed_full_field import TimedFullField

        mesh = Mesh1D(
            edges=np.linspace(0.0, 20.0, 41),
            mat_ids=np.array([0] * 20 + [1] * 20),
            bc_left=BC("reflective"),
            bc_right=BC("vacuum"),
        )
        sn_mesh = SNMesh(
            mesh,
            Quadrature.gauss_legendre(n_ordinates=4),
            {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")},
        )
        corrector = DSACorrection.from_sn_mesh(sn_mesh)
        psi = TimedFullField(
            interior=AngularFlux.zeros(sn_mesh.angular_bulk_space),
            boundary=AngularBoundaryFlux.zeros(sn_mesh.angular_trace),
            _history=(), history_depth=0,
        )
        correction = corrector.apply(psi - psi)
        np.testing.assert_array_equal(
            correction.interior.values,
            np.zeros_like(correction.interior.values),
        )
        np.testing.assert_array_equal(
            correction.boundary.values,
            np.zeros_like(correction.boundary.values),
        )


class TestTeeth:
    """The gates are constrained — seeded corruptions red them."""

    def test_sign_flipped_correction_breaks_convergence(self, monkeypatch):
        """A sign-flipped correction amplifies instead of deflates: the
        D3 convergence witness must break."""
        original = DSACorrection.apply

        def flipped(self, displacement):
            return original(self, displacement) * (-1.0)

        monkeypatch.setattr(DSACorrection, "apply", flipped)
        mutated = _solve(("vacuum", "vacuum"), acceleration="dsa",
                         max_inner=200)
        history = mutated.history
        if history is None or not history.flux_residuals:
            pytest.fail("the mutated run must carry a residual history")
        last_res = history.flux_residuals[-1]
        # The residual VALUE is the witness (an iteration-count bar is
        # off-by-one-prone: a diverging run reports max_inner−1). The
        # healthy accelerator reaches 1e-11 in ~33 inners.
        if last_res < 1e-6:
            pytest.fail(
                "a sign-flipped correction must not converge "
                f"(last residual {last_res:.2e}) — the D3 witness has "
                "no teeth"
            )

    def test_zeroed_trace_arm_breaks_the_reflective_case(self, monkeypatch):
        """The 3b root-cause regression pin: withOUT the wall-edge trace
        correction, the reflective case diverges (measured ρ > 1 — the
        lagged reflective gain reads an uncorrected outflow while the
        low-order's reflecting row assumes the corrected state)."""
        original = DSACorrection.apply

        def bulk_only(self, displacement):
            out = original(self, displacement)
            return out._recombine(
                interior=out.interior, boundary=out.boundary * 0.0,
            )

        monkeypatch.setattr(DSACorrection, "apply", bulk_only)
        mutated = _solve(("reflective", "vacuum"), acceleration="dsa",
                         max_inner=120)
        history = mutated.history
        if history is None or not history.flux_residuals:
            pytest.fail("the mutated run must carry a residual history")
        last_res = history.flux_residuals[-1]
        # Divergence witness on the residual VALUE (measured ~1e+35 at
        # 120 inners across every reflective regime; the healthy trace
        # arm converges in ~33).
        if last_res < 1.0:
            pytest.fail(
                "the bulk-only mutation must DIVERGE on the reflective "
                f"case (last residual {last_res:.2e}) — the trace arm "
                "would be unconstrained"
            )
