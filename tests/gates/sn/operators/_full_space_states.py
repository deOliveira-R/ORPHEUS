r"""Four small within-group problems and residuals that populate every block of the full typed state.

The fixtures of the composite sweep-inverse identity
(``test_sweep_inverse_identity.py``, ERR-071 and ERR-078), shared with the
gates of the sweep-preconditioned Krylov solve (``tests/gates/sn/solve/
test_krylov_sweep_preconditioner.py``, issue #200): both ask what an operator
does to a residual with EVERY block live (the bulk, the inflow and outflow
trace rows, and on a carrying mesh the ψ½ System-B blocks), which no physical
right-hand side populates.

* :data:`MESHES` / :data:`GEOMS`: a five-region slab with a vacuum or a
  reflective left face, a 2-D Cartesian box (two materials, non-uniform edges,
  reflective and vacuum faces), a folded-product cylinder whose rule carries
  bit-exact pure-azimuthal ordinates, and a Gauss-Legendre sphere. Every trace
  row is live except the cylinder's tangential slots; the two curvilinear
  meshes carry System B.
* :func:`random_state`: a random SOURCE-role residual in an implicit
  operator's domain, every block populated.
* :func:`random_source_composite`: a random SOURCE-role right-hand side on a
  seedless (slab) mesh, every block populated: what source iteration's
  lagged sum adds its gains to.
* :func:`system_a`: the bulk (System-A) member of a possibly-coupled field.
"""
from __future__ import annotations

from typing import Any

import numpy as np

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry import BC, CoordSystem, StructuredGeometry
from orpheus.mesh import CellEdges, CellsByCount, Mesh2D, Mesher
from orpheus.numerics.coupled_system import CoupledField, CoupledOperator
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.problem import SNProblem
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.full_field import FullField
from orpheus.transport.radial_characteristic_field import RadialCharacteristicField


def mixtures():
    mix_a = make_mixture(
        sig_t=np.array([1.0]), sig_c=np.array([0.1]),
        sig_f=np.array([0.0]), nu=np.array([0.0]),
        chi=np.array([0.0]), sig_s=np.array([[0.9]]),
    )
    mix_b = make_mixture(
        sig_t=np.array([2.0]), sig_c=np.array([1.0]),
        sig_f=np.array([0.0]), nu=np.array([0.0]),
        chi=np.array([0.0]), sig_s=np.array([[1.0]]),
    )
    return {0: mix_a, 1: mix_b}


def mesh_slab(left: str) -> SNProblem:
    geom = StructuredGeometry.slab(
        (0.0, 0.5, 3.0, 5.0, 6.0, 8.0), (0, 1, 0, 1, 0),
        left=BC(left), right=BC("vacuum"),
    )
    return SNProblem(
        Mesher(geom).partition((
            CellEdges(np.array([0.0, 0.5])),
            CellEdges(np.array([0.5, 1.5, 3.0])),
            CellEdges(np.array([3.0, 5.0])),
            CellEdges(np.array([5.0, 6.0])),
            CellEdges(np.array([6.0, 8.0])),
        )).mesh,
        Quadrature.gauss_legendre(n_ordinates=4),
        mixtures(),
    )


def mesh_cart2d() -> SNProblem:
    # The 2-D Cartesian row (qa's #200 review, item 4): a 3 x 2 box of two
    # materials on non-uniform edges, reflective on xmin and ymax and vacuum
    # on xmax and ymin (reflection is a coupling GAIN, never inside the bare
    # (L+C), so both kinds of face carry live trace rows). The level-symmetric
    # S4 rule has no ordinate tangent to a face: ``[M]`` 2026-10-05, 12 inflow
    # and 12 outflow slots on each of the four faces, 24 ordinates (a 2-D
    # ``product`` rule leaves half of each face's slots tangential). Seedless,
    # so the implicit operator is the bare (L+C).
    return SNProblem(
        Mesh2D(
            np.array([0.0, 0.4, 1.5, 2.0]), np.array([0.0, 0.7, 2.0]),
            np.array([[0, 1], [1, 0], [0, 0]]),
            face_laws={"xmin": BC("reflective"), "xmax": BC("vacuum"),
                       "ymin": BC("vacuum"), "ymax": BC("reflective")},
            coord=CoordSystem.CARTESIAN,
        ),
        Quadrature.level_symmetric(4),
        mixtures(),
    )


def mesh_cyl() -> SNProblem:
    # The #280 MANDATORY cylinder config, re-posed at the 6.3 flip onto
    # the admitted family: ``folded_product(4, 6)`` — the staggered
    # parent at n_φ ≡ 2 (mod 4) places φ = π/2 exactly, and the
    # roots-of-unity circle (E3) makes those ordinates' μ_r = 0.0
    # BIT-EXACT — so the rule carries degenerate pure-azimuthal
    # ordinates AND (like every admitted cylinder rule) a live ψ½
    # System B.  xmax-only trace layout — exercises the restore's
    # per-face membership loop on the curvilinear face set.
    return SNProblem(
        Mesher(StructuredGeometry.cylinder(
            (0.0, 0.3, 0.8, 1.0), (0, 1, 0), outer=BC("vacuum"),
        )).partition(CellsByCount.uniform_width(1)).mesh,
        Quadrature.folded_product(n_mu=4, n_phi=6),
        mixtures(),
    )


def mesh_sphere() -> SNProblem:
    # The carrying SPHERE row — GL-4, one seed level; the same coupled
    # round-trip as the cylinder row.  The sphere's seed-carrying
    # inverse reciprocity was explicitly deferred pre-6.3 (the "#29
    # domain" note in test_loss_transpose_solve.G3); ERR-078's fix
    # covers both curvilinear arms, so both are gated here.
    return SNProblem(
        Mesher(StructuredGeometry.sphere(
            (0.0, 0.3, 0.8, 1.0), (0, 1, 0), outer=BC("vacuum"),
        )).partition(CellsByCount.uniform_width(1)).mesh,
        Quadrature.gauss_legendre(n_ordinates=4),
        mixtures(),
    )


MESHES = {
    # vacuum walls — every trace row live: inflow identities, outflow defects
    "slab_vacuum": lambda: mesh_slab("vacuum"),
    # the identity is bc-INDEPENDENT: B is a coupling GAIN, never inside
    # the bare (L+C) — a reflective wall must not change the round-trip
    "slab_reflective": lambda: mesh_slab("reflective"),
    "cart2d": mesh_cart2d,
    "cyl_folded": mesh_cyl,
    "sphere_gl": mesh_sphere,
}

GEOMS = list(MESHES)


def zero_source_composite(problem: SNProblem) -> FullField:
    """A zero SOURCE-role System-A carrier for the coupled arm's rhs.

    Role-honest member algebra: a solve's rhs is a SOURCE, and the
    substitution computes ``q_A − Seeding·ψ_B`` — cross-role arithmetic
    is forbidden by the typed fields, so a flux-role rhs raises at the
    block boundary (#289-F2)."""
    from orpheus.transport.source_sinks import (
        AngularBoundarySourceSink,
        AngularSourceSink,
    )

    return FullField(
        interior=AngularSourceSink(values=np.zeros((problem.quad.N, problem.ng, *problem.spatial_shape)), space=problem.angular_bulk_space),
        boundary=AngularBoundarySourceSink.zeros(problem.angular_trace),
    )


def random_state(problem: SNProblem, lc, seed: int):
    """A random rhs in ``lc``'s domain — the bare composite, or the
    coupled (bulk ⊕ ψ½) SOURCE-role state with EVERY member block
    populated (randomized through the coupled ``from_flat``, so the
    ψ½ member's interior ⊕ boundary blocks are live too)."""
    if not isinstance(lc, CoupledOperator):
        return random_composite(problem, seed)
    template = CoupledField(systems=(
        zero_source_composite(problem),
        RadialCharacteristicField.source_zeros(problem.radial_characteristic_field_space),
    ))
    flat = np.asarray(template.to_flat())
    rng = np.random.default_rng(seed + 1)
    return CoupledField.from_flat(rng.normal(size=flat.shape), template)


def system_a(x: Any) -> Any:
    """Project the bulk (System-A) member of a possibly-coupled field."""
    return x.systems[0] if isinstance(x, CoupledField) else x


def random_composite(problem: SNProblem, seed: int) -> FullField:
    """Every block populated — bulk, inflow-trace, AND the outflow-trace
    rows the old sweep dropped — with shapes read off the mesh (so the
    same builder serves slab and the xmax-only curvilinear layout)."""
    rng = np.random.default_rng(seed)
    interior = AngularFlux(values=rng.normal(size=(problem.quad.N, problem.ng, *problem.spatial_shape)), space=problem.angular_bulk_space)
    boundary = AngularBoundaryFlux.zeros(problem.angular_trace)
    for face in boundary.layout.faces:
        view = boundary.face_view(face)
        view[...] = rng.normal(size=view.shape)
    return FullField(interior=interior, boundary=boundary)


def random_source_composite(problem: SNProblem, seed: int) -> FullField:
    """The SOURCE-role sibling of :func:`random_composite`: a random right-hand side with the bulk and every trace
    row populated, in the role a fixed-point step's lagged source ``q + N ψ`` carries (the typed fields refuse a
    flux-role ``q`` there)."""
    from orpheus.transport.source_sinks import AngularBoundarySourceSink, AngularSourceSink

    rng = np.random.default_rng(seed)
    interior = AngularSourceSink(
        values=rng.normal(size=(problem.quad.N, problem.ng, *problem.spatial_shape)), space=problem.angular_bulk_space,
    )
    boundary = AngularBoundarySourceSink.zeros(problem.angular_trace)
    for face in boundary.layout.faces:
        view = boundary.face_view(face)
        view[...] = rng.normal(size=view.shape)
    return FullField(interior=interior, boundary=boundary)
