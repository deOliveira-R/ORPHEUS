r"""The System-B CYLINDER arm (Q5.6): the seed machinery on folded arcs.

The σ_y-folded cylinder (``Quadrature.folded_product``) is the first
cylindrical configuration whose every level CARRIES an independent ψ½
seed — the machinery built for the sphere activates on a second
geometry, and its landing battery caught it sphere-hardcoded in three
places (each fixed here, each with its own gate):

1. **The q̄½ source fold** applied the sphere's Legendre
   analysis-synthesis AT μ = ±1 — a point outside the arc — to cylinder
   levels.  `[M]` on the L0 flat-flux problem the seed read the two GL
   level families' garbage values (+3.7217 / −0.5801, a NEGATIVE angular
   flux on a positive-source pure absorber) and the solution was 82 %
   off.  The arc's fold is the same analysis-synthesis shape in the
   arc's OWN Gauss family: the staggered arc is Gauss–Chebyshev-1 in
   x = cos ω (T25), so ``T_k(±1) = (±1)^k`` replaces
   ``P_ℓ(±1) = (±1)^ℓ`` and the GC1 discrete orthogonality makes the
   analysis weight-free.
2. **The DD march** baked the sphere's diameter ray (|μ| = 1) into its
   optical step.  The general step is the PATH length ``Δr/|η_start|``
   (level p's ξ = 0 ray traverses Δr of radius over Δr/sinθ_p of path).
   Flat-flux-INVISIBLE (the DD fixed point q̄/σ is cosine-independent),
   so the L0 equilibria cannot see it — the linear-ψ½ gate below is its
   catcher (DD is exact for a linear flux, and the manufactured source
   carries the ``−s_p·B`` streaming term a cosine-less march drops).
3. **The isotropic Emission fold** hard-coded the sphere's ℓ = 0 weight
   ``½ = 1/Σw(GL)``.  On the folded cylinder (Σw = 4π) that
   over-injects the ray scattering source 2π-fold — `[M]` the c = 0.4
   flat-flux equilibrium read 158 % off; the constant function's
   reproducing weight ``1/Σw`` is the geometry-generic identity.

All sphere paths are byte-identical (every generalization divides by
the sphere's 1.0 or multiplies by its ½ through the same spelling).

`[M]` at the fix (folded_product(4,8), 10-cell reflective/reflective
cylinder, serial ``-O``): c = 0 equilibrium 82 % → 2.8e-13; c = 0.4
equilibrium 158 % → 2.2e-13.
"""

from __future__ import annotations

import numpy as np
import pytest

from tests.sn._test_helpers import rc_march

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn import solve_sn_fixed_source
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operators.radial_characteristic import (
    RadialCharacteristicOperator,
)
from orpheus.transport.fields.cross_section_field import CrossSectionField
from orpheus.transport.radial_characteristic_field import (
    RadialCharacteristicField,
)

pytestmark = pytest.mark.l0

_SIG_T = 1.0
_NX = 10
_RADIUS = 2.0


def _mixture(c: float):
    ss = np.array([[c * _SIG_T]])
    return make_mixture(
        sig_t=np.array([_SIG_T]), sig_c=np.array([_SIG_T * (1.0 - c)]),
        sig_f=np.zeros(1), nu=np.zeros(1), chi=np.zeros(1), sig_s=ss,
    )


def _folded_cylinder(c: float = 0.0) -> tuple[SNMesh, Mesh1D]:
    mesh = Mesh1D(
        edges=np.linspace(0.0, _RADIUS, _NX + 1),
        mat_ids=np.zeros(_NX, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_left=BC("reflective"),
        bc_right=BC("reflective"),
    )
    sn = SNMesh(mesh, Quadrature.folded_product(4, 8), {0: _mixture(c)})
    return sn, mesh


# ── the L0 equilibria — per ordinate, never balance-only (vv #8) ────────

@pytest.mark.parametrize("c", [0.0, 0.4], ids=["pure_absorber", "c04"])
def test_flat_flux_equilibrium_per_ordinate(c: float):
    """ψ_n = q/(σ_t(1−c)) per ordinate, and the ψ½ seed carries the same
    constant — the folded cylinder's streaming-equilibrium identity.

    The c = 0 row pins the SOURCE fold (defect 1: it read 82 % off);
    the c = 0.4 row additionally pins the EMISSION weight (defect 3:
    158 % off).  Both are per-ordinate assertions — a scalar-flux
    balance would telescope over exactly the redistribution terms the
    seed feeds (vv anti-pattern #8).
    """
    sn, mesh = _folded_cylinder(c)
    q = np.ones((sn.quad.N, 1, _NX))
    res = solve_sn_fixed_source(
        {0: _mixture(c)}, mesh, sn.quad, q, max_inner=800, inner_tol=1e-13,
    )
    psi_exact = 1.0 / (_SIG_T * (1.0 - c))
    psi = np.asarray(res.angular_flux.interior.values)     # (N, ng, nx)
    np.testing.assert_allclose(
        psi, psi_exact, rtol=5e-12, atol=0.0,
        err_msg=f"per-ordinate equilibrium broken at c={c}",
    )
    seed = np.asarray(res.radial_characteristic.interior.values)
    np.testing.assert_allclose(
        seed, psi_exact, rtol=5e-12, atol=0.0,
        err_msg=f"the ψ½ seed left the flat equilibrium at c={c}",
    )


# ── defect 1's term gate: the fold IS the arc-endpoint synthesis ────────

def test_source_fold_is_the_arc_endpoint_synthesis():
    """A per-level Chebyshev-polynomial source folds to its exact endpoint
    values: q = a + b·T₁(x) + c·T₂(x) ⟹ q̄½(∓) = a ∓ b + c.

    ``x = η/sinθ_p`` are the level's GC1 nodes; the synthesis at the
    start ray ω = π is x = −1 with T_k(−1) = (−1)^k.  The sphere's
    Legendre fold run on these levels produced the +3.72/−0.58 garbage
    this gate's family replaces.
    """
    sn, _ = _folded_cylinder()
    rng = np.random.default_rng(20260807)
    a = rng.uniform(0.5, 2.0, _NX)
    b = rng.uniform(-0.5, 0.5, _NX)
    c2 = rng.uniform(-0.3, 0.3, _NX)

    vals = np.empty((sn.quad.N, 1, _NX))
    for p, ords in enumerate(sn.quad.level_indices):
        ords = np.asarray(ords)
        mu_z0 = float(sn.quad.mu_z[ords[0]])
        x = sn.quad.mu_x[ords] / float(np.sqrt(1.0 - mu_z0**2))
        t2 = 2.0 * x * x - 1.0
        vals[ords, 0, :] = (
            a[None, :] + b[None, :] * x[:, None] + c2[None, :] * t2[:, None]
        )

    seed = RadialCharacteristicField.source_from_angular(vals, sn)
    assert seed is not None
    for p in seed.interior.space.levels:
        np.testing.assert_allclose(
            seed.interior.cells(p, -1)[0], a - b + c2, rtol=0, atol=1e-13,
            err_msg=f"level {p}: the ω = π endpoint synthesis is wrong",
        )
        np.testing.assert_allclose(
            seed.interior.cells(p, +1)[0], a + b + c2, rtol=0, atol=1e-13,
            err_msg=f"level {p}: the ω = 0 endpoint synthesis is wrong",
        )


# ── defect 2's term gate: the march rides the level's own ray ───────────

def test_the_march_rides_the_level_ray_exactly_for_linear_flux():
    """DD is exact for ψ½ linear in r — with the source manufactured ON
    level p's ξ = 0 ray: q = ∓sinθ_p·B + σ(A + B·r).

    The ``∓sinθ_p·B`` streaming term is exactly what a cosine-less
    march (the sphere's |μ| = 1 baked in) mis-weights, so this gate is
    defect 2's designed catcher — the flat equilibria are blind to it
    (the DD fixed point q̄/σ is cosine-independent).
    """
    A, B = 0.7, 0.3
    sn, mesh = _folded_cylinder()
    r = np.asarray(mesh.centers)
    sigma_field = CrossSectionField(values=np.full((1, _NX), _SIG_T), space=sn.bulk_space)
    op = rc_march(sn, sigma_field)

    src = RadialCharacteristicField.source_zeros(sn.radial_characteristic_field_space)
    psi_lin = A + B * r                                     # (nx,)
    for p in src.interior.space.levels:
        ords = np.asarray(sn.quad.level_indices[p])
        mu_z0 = float(sn.quad.mu_z[ords[0]])
        s_p = float(np.sqrt(1.0 - mu_z0**2))
        src.interior.cells(p, -1)[...] = -s_p * B + _SIG_T * psi_lin
        src.interior.cells(p, +1)[...] = +s_p * B + _SIG_T * psi_lin
        src.boundary.corner(p, -1)[...] = A + B * _RADIUS

    flux = op.solve(src)
    for p in flux.interior.space.levels:
        np.testing.assert_allclose(
            flux.interior.cells(p, -1)[0], psi_lin, rtol=1e-12,
            err_msg=(
                f"level {p}: the inward leg does not ride its own ray — "
                f"the march's optical step is not Δr·σ/sinθ_p"
            ),
        )
        np.testing.assert_allclose(
            flux.interior.cells(p, +1)[0], psi_lin, rtol=1e-12,
            err_msg=f"level {p}: the outward leg does not ride its own ray",
        )
        np.testing.assert_allclose(
            flux.boundary.corner(p, +1)[0], A + B * _RADIUS, rtol=1e-12,
            err_msg=f"level {p}: the outflow corner left the linear flux",
        )


# ── the resolvent identity on the new arm ───────────────────────────────

def test_apply_solve_round_trip_closes():
    """``A_BB(A_BB⁻¹ q) = q`` on the folded cylinder — the forward and the
    march must spell the SAME per-level path widths (dr/sinθ_p), or the
    round-trip drifts.

    The contract's corner shape (the sphere's documented property, now on
    the second geometry): the (−1) inflow corner is the identity row; the
    (+1) source corner is an UNUSED slot (the q½ fold never writes it),
    and apply's (+1) row is the streamed−stored self-consistency defect —
    exactly 0.0 for a solved state.
    """
    sn, _ = _folded_cylinder()
    sigma_field = CrossSectionField(values=np.full((1, _NX), _SIG_T), space=sn.bulk_space)
    op = rc_march(sn, sigma_field)

    rng = np.random.default_rng(20260808)
    src = RadialCharacteristicField.source_zeros(sn.radial_characteristic_field_space)
    src.interior.values[...] = rng.uniform(0.5, 2.0, src.interior.values.shape)
    for p in src.interior.space.levels:
        src.boundary.corner(p, -1)[...] = rng.uniform(0.5, 2.0, size=1)

    round_trip = op.apply(op.solve(src))
    np.testing.assert_allclose(
        round_trip.interior.values, src.interior.values, rtol=1e-11,
        err_msg="A_BB ∘ A_BB⁻¹ drifted on the interior cells",
    )
    for p in src.interior.space.levels:
        np.testing.assert_allclose(
            round_trip.boundary.corner(p, -1),
            src.boundary.corner(p, -1),
            rtol=1e-12,
            err_msg=f"level {p}: the inflow-corner identity row drifted",
        )
        np.testing.assert_array_equal(
            round_trip.boundary.corner(p, +1),
            np.zeros_like(round_trip.boundary.corner(p, +1)),
            err_msg=(
                f"level {p}: the outflow self-consistency defect is "
                f"nonzero — the forward does not replay the march's "
                f"path-length arithmetic"
            ),
        )
