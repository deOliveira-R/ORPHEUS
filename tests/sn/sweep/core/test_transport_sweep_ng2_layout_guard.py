r"""ng≥2 layout guard for ``transport_sweep`` — the ERR-055 loud-failure gate.

ERR-055 (2026-06-01): six curvilinear sweep-regression tests fed
``sig_t`` / ``Q`` in the obsolete ``(nx, ng, ny)`` layout after the
production contract flipped to PR-INDEX-5 ``(ng, nx, ny)`` /
``(N, ng, nx, ny)``. The bug only crashed at ng=1, where
``(nx, ng, ny)`` and ``(ng, nx, ny)`` *alias* in element count and
shape-rank — the axis swap declares itself only at a slice/index that
distinguishes the two axes (``CollisionCache.from_geometry``).

**This guard exists so the convention can never drift silently again.**
It runs a 1-D *multigroup* (ng=2) sweep through the SAME production
producer the solver uses (``AngularSourceSink.from_isotropic`` +
``transport_sweep``). With ng=2 the layout aliasing is broken: a
``(nx, ng, ny)`` mismatch makes axis 1 size 2 ≠ nx, so any future
producer/consumer ``(ng, nx, ny)`` drift fails LOUDLY here (axis-size
mismatch) instead of hiding behind the ng=1 degeneracy.

It asserts the canonical PR-INDEX-5 output shapes explicitly:
``ang.shape == (N, ng, nx, ny)`` and ``phi.shape == (ng, nx, ny)`` —
the single most direct pin on the layout contract.

cross-cutting hygiene rule H1 (1-group degeneracy) operating at the
data-layout level: ALWAYS exercise ≥2 groups so axis-swap drift
surfaces. coding-elegance Pattern 7: route through the producer so the
test obeys the SAME convention production does.
"""

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, CoordSystem, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.geometry import SNMesh
from orpheus.sn.sweep import transport_sweep
from orpheus.transport.source_sinks import AngularSourceSink
from orpheus.transport.fields.boundary_flux import BoundaryFlux

pytestmark = pytest.mark.foundation


def _slab_1d(n_cells: int, ng_mix) -> Mesh1D:
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC.reflective, BC.reflective),
    )
    return Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=n_cells),))


def test_transport_sweep_ng2_layout_shapes():
    """1-D ng=2 sweep: output shapes MUST be the PR-INDEX-5 contract.

    Routes a 2-group source + ``sig_t`` through the production producer
    (``AngularSourceSink.from_isotropic`` + ``transport_sweep``) on a
    10-cell slab and asserts ``(N, ng, nx, ny)`` / ``(ng, nx, ny)``. The
    nx=10 ≠ ng=2 asymmetry is the load-bearing choice: it makes a
    ``(nx, ng, ny)`` axis-swap a SIZE mismatch (loud) rather than the
    silent ng=1 alias that hid ERR-055.
    """
    mix = get_mixture("A", "2g")
    ng = 2
    nx = 10
    mesh = _slab_1d(nx, mix)
    quad = Quadrature.gauss_legendre(8)
    sn_mesh = SNMesh(mesh, quad, {0: mix})

    # Principled PR-INDEX-5 inputs: sig_t (ng, nx, ny), isotropic
    # scalar source (ng, nx, ny). A producer/consumer drift to the
    # obsolete (nx, ng, ny) layout would mismatch axis 1 (ng=2 vs nx=10)
    # at CollisionCache.from_geometry and crash loudly here.
    sig_t = np.broadcast_to(mix.SigT[:, None, None], (ng, nx, 1)).copy()
    Q_iso = np.ones((ng, nx, 1))
    source = AngularSourceSink.from_isotropic(Q_iso, sn_mesh)
    boundary_flux = BoundaryFlux.zeros_on(sn_mesh)

    ang, phi = transport_sweep(source, sig_t, sn_mesh, boundary_flux)

    assert ang.shape == (quad.N, ng, nx, 1), (
        f"angular flux layout drift: got {ang.shape}, "
        f"expected (N={quad.N}, ng={ng}, nx={nx}, ny=1)"
    )
    assert phi.shape == (ng, nx, 1), (
        f"scalar flux layout drift: got {phi.shape}, "
        f"expected (ng={ng}, nx={nx}, ny=1)"
    )
    assert np.all(np.isfinite(ang)), "Non-finite angular flux"
    assert np.all(np.isfinite(phi)), "Non-finite scalar flux"


def test_transport_sweep_ng2_per_group_distinct():
    """The two groups must carry DISTINCT flux (no axis collapse/aliasing).

    Group 0 and group 1 of mixture A have different Σ_t, so an
    equilibrium sweep (uniform Q, many sweeps) lands each group at
    φ_g = Q_g/Σ_t,g — distinct per group. If a layout drift collapsed
    or transposed the group axis, the per-group fluxes would coincide
    or scramble; this asserts they are genuinely group-resolved.
    """
    mix = get_mixture("A", "2g")
    ng = 2
    nx = 4
    mesh = _slab_1d(nx, mix)
    quad = Quadrature.gauss_legendre(8)
    sn_mesh = SNMesh(mesh, quad, {0: mix})

    sig_t = np.broadcast_to(mix.SigT[:, None, None], (ng, nx, 1)).copy()
    Q_iso = np.ones((ng, nx, 1))
    source = AngularSourceSink.from_isotropic(Q_iso, sn_mesh)
    boundary_flux = BoundaryFlux.zeros_on(sn_mesh)

    phi = None
    for _ in range(200):
        _, phi = transport_sweep(source, sig_t, sn_mesh, boundary_flux)

    # Equilibrium per group: φ_g = Q_g / Σ_t,g (pure-streaming sweep with
    # no scatter coupling — transport_sweep is the within-group sweep).
    expected = Q_iso / sig_t  # (ng, nx, 1)
    np.testing.assert_allclose(
        phi, expected, rtol=1e-6,
        err_msg="per-group equilibrium φ ≠ Q_g/Σ_t,g — group axis drift?",
    )
    # The two groups MUST differ (mixture A has Σ_t,0 ≠ Σ_t,1).
    assert not np.allclose(phi[0], phi[1]), (
        "group-0 and group-1 flux coincide — the ng axis collapsed "
        "(layout drift would alias the groups)"
    )
