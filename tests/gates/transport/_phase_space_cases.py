r"""The material meshes the transport gates pose their spaces on (layer DAG Phase 4).

A transport gate holds no method problem: it builds a :class:`MaterialMesh` here, picks a quadrature and a scheme
itself, and builds every space it reads with :mod:`orpheus.transport.phase_spaces` (or the scheme's
``trial_space``, ``FullFieldSpace.from_blocks``, ``HarmonicFrame.moment_space_on``). Nothing here owns a space: a
builder returns a mesh, and the gate spells each space at the call site, so a gate cannot read a space through a
duck-typed carrier (the F1 entries refuse a ``MaterialMesh`` by name).

Two populations:

* the GEOMETRY population (``slab``, ``stretched_slab``, ``sphere``, ``hollow_sphere``, ``cylinder``,
  ``hollow_cylinder``, ``cart2d``): one material, ``ng`` groups of placeholder data, uniform cells; the chart is
  what varies;
* the MONOMORPHIC population (``monomorphic()``): the four two-region, two-group, fissile and anisotropic meshes of
  ``tests/gates/sn/architecture/test_monomorphic_leaves.py`` (non-uniform cells; the sphere's cell volumes span
  three decades), moved here so the transport and numerics rows that read their spaces build the SAME meshes; the
  SN file builds its problems on these.
"""
from __future__ import annotations

from collections.abc import Callable, Mapping
from functools import cache
from typing import NamedTuple

import numpy as np

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.geometry import CoordSystem, StructuredGeometry
from orpheus.geometry.deck import MirrorDeck
from orpheus.mesh import CellEdges, CellsByCount, Mesh1D, Mesh2D, Mesher
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.spaces.full_field_space import FullFieldSpace
from orpheus.transport.mesh.material_mesh import MaterialMesh
from orpheus.transport.phase_spaces import angular_bulk_space, angular_trace_space, scalar_trace_space
from orpheus.transport.spatial import DiamondDifference
from orpheus.transport.spatial.scheme import DiscretizationSchemeBase, StreamingTerms
from tests.gates.transport._materials import anisotropic_mixture, placeholder_materials


class Posed(NamedTuple):
    """The three ingredients a transport test poses its spaces from: a material mesh, a quadrature, a closure.

    It holds ingredients only and exposes no space: every space is spelled at the call site from the raw
    constructors, so a carrier-keyed entry handed a ``Posed`` refuses it (layer DAG Phase 4, spec section 2.0)."""

    mesh: MaterialMesh
    quadrature: Quadrature
    scheme: DiscretizationSchemeBase


# ── The composites (one spelling each of "interior + trace", from the raw constructors) ──


def angular_full_field_space(
    mesh: MaterialMesh, quadrature: Quadrature, scheme: DiscretizationSchemeBase | None = None,
) -> FullFieldSpace:
    """The SN full-field space ``trial (+) angular trace`` of ``mesh`` under ``quadrature`` and ``scheme`` (DD by default)."""
    closure = DiamondDifference() if scheme is None else scheme
    trial = closure.trial_space(angular_bulk_space(mesh, quadrature), mesh.axes)
    return FullFieldSpace.from_blocks(trial, angular_trace_space(mesh, quadrature, closure))


def moment_space(
    mesh: MaterialMesh, quadrature: Quadrature, L: int, scheme: DiscretizationSchemeBase | None = None,
):
    """The harmonic moment space at order ``L`` of ``mesh`` under ``quadrature``: the frame's
    ``HarmonicFrame.for_space(V, L).moment_space_on(V)`` of the angular bulk ``V`` (widened by ``scheme.trial_space``
    when a multi-moment ``scheme`` is given). The second owner of the moment space (ruling O-5); the SN problem's keyed
    ``moment_space`` is the first, and ``tests/gates/sn/mesh/test_hub_and_frame_agree_on_the_moment_space.py`` pins
    the two equal."""
    from orpheus.transport.frames import HarmonicFrame

    space = angular_bulk_space(mesh, quadrature)
    if scheme is not None:
        space = scheme.trial_space(space, mesh.axes)
    return HarmonicFrame.for_space(space, L).moment_space_on(space)


def scalar_full_field_space(mesh: MaterialMesh) -> FullFieldSpace:
    """The diffusion full-field space ``scalar bulk (+) scalar trace`` of ``mesh``."""
    return FullFieldSpace.from_blocks(mesh.bulk_space, scalar_trace_space(mesh))


def slab_streaming_terms(mesh: Mesh1D, quadrature: Quadrature, cell: int, ordinate: int) -> StreamingTerms:
    """One slab cell's streaming terms for one ordinate, by hand: the curvature-neutral faces (area 1 on both sides),
    the cell's width as its volume, ``|mu|`` of the ordinate. The SN slab operator feeds the same terms
    (``tests/gates/sn/mesh/test_reduced_operator.py::test_the_slab_streaming_terms_are_the_neutral_hand_terms``)."""
    edges = np.asarray(mesh.edges, dtype=float)
    return StreamingTerms(
        face_area_inner=1.0, face_area_outer=1.0, volume=float(edges[cell + 1] - edges[cell]),
        abs_mu=abs(float(quadrature.mu_x[ordinate])),
    )


def posing_quadrature(mesh: MaterialMesh, n: int = 4) -> Quadrature:
    """The quadrature a transport gate poses ``mesh`` with, chosen by its CHART: Gauss-Legendre(``n``) on a slab or a
    sphere, the folded product (``n`` polar x ``2n`` azimuthal, the cylinder's azimuthal levels) on a cylinder,
    level-symmetric S(``n``) on a 2-D Cartesian mesh (no ``mu_y == 0`` ordinate). ONE instance per (chart, n): the
    frames intern per rule instance, so a gate asserting interning needs the same rule on every call."""
    if mesh.ndim == 2:
        return _rule("level_symmetric", n)
    if mesh.coord is CoordSystem.CYLINDRICAL:
        return _rule("folded_product", n)
    return _rule("gauss_legendre", n)


@cache
def _rule(family: str, n: int) -> Quadrature:
    match family:
        case "level_symmetric":
            return Quadrature.level_symmetric(sn_order=n)
        case "folded_product":
            return Quadrature.folded_product(n_mu=n, n_phi=2 * n)
    return Quadrature.gauss_legendre(n_ordinates=n)


def posing_scheme(mesh: MaterialMesh) -> DiscretizationSchemeBase:
    """The closure a transport gate poses a mesh with unless it varies it: Diamond Difference (slopeless)."""
    return DiamondDifference()


def slab_mesh(nx: int = 4, ng: int = 2) -> MaterialMesh:
    """A 1-D slab of width 1, ``nx`` equal cells, its two faces tagged ``0``, ``ng`` groups of placeholder data."""
    mesh = Mesher(StructuredGeometry.from_homogeneous(1.0, 0)).partition(CellsByCount.uniform_width(nx)).mesh
    return MaterialMesh(mesh, placeholder_materials(ng=ng))


def _one_d(geometry: StructuredGeometry, nx: int, ng: int) -> MaterialMesh:
    mesh = Mesher(geometry).partition(CellsByCount.uniform_width(nx)).mesh
    return MaterialMesh(mesh, placeholder_materials(ng=ng))


def slab(nx: int = 4, ng: int = 2, width: float = 1.0) -> MaterialMesh:
    """A mirror-bounded slab of ``nx`` equal cells."""
    return _one_d(StructuredGeometry.slab((0.0, width), (0,), left=MirrorDeck(), right=MirrorDeck()), nx, ng)


def stretched_slab(nx: int = 4, ng: int = 2) -> MaterialMesh:
    """``slab`` at twice the width: equal shapes, unequal cell volumes (a content discriminator)."""
    return slab(nx, ng, width=2.0)


def sphere(nx: int = 4, ng: int = 2, inner: float = 0.0) -> MaterialMesh:
    roles = {"outer": 0} if inner == 0.0 else {"outer": 0, "inner": 0}
    return _one_d(StructuredGeometry.sphere((inner, 1.0), (0,), **roles), nx, ng)


def hollow_sphere(nx: int = 4, ng: int = 2) -> MaterialMesh:
    return sphere(nx, ng, inner=0.5)


def cylinder(nx: int = 4, ng: int = 2, inner: float = 0.0) -> MaterialMesh:
    roles = {"outer": 0} if inner == 0.0 else {"outer": 0, "inner": 0}
    return _one_d(StructuredGeometry.cylinder((inner, 1.0), (0,), **roles), nx, ng)


def hollow_cylinder(nx: int = 4, ng: int = 2) -> MaterialMesh:
    return cylinder(nx, ng, inner=0.5)


def cart2d(nx: int = 3, ny: int = 2, ng: int = 2) -> MaterialMesh:
    """A non-square (``nx != ny``) mirror-bounded rectangle: a transposed face shape is visible."""
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 1.0, nx + 1), edges_y=np.linspace(0.0, 1.0, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
        face_roles={face: MirrorDeck() for face in ("xmin", "xmax", "ymin", "ymax")},
    )
    return MaterialMesh(mesh, placeholder_materials(ng=ng))


def legacy_mesh(material_mesh: MaterialMesh) -> "Mesh1D | Mesh2D":
    """The ``Mesh1D`` / ``Mesh2D`` a material mesh of this module was built on (every d <= 2 mesh carries one)."""
    mesh = material_mesh.mesh
    if mesh is None:
        raise ValueError("a d <= 2 material mesh always carries its legacy mesh")
    return mesh


GEOMETRIES: Mapping[str, Callable[[], MaterialMesh]] = {
    "slab": slab, "sphere": sphere, "hollow_sphere": hollow_sphere, "cylinder": cylinder,
    "hollow_cylinder": hollow_cylinder, "cart2d": cart2d,
}

#: The quadrature each geometry is posed with (``posing_quadrature``).
QUADRATURES: Mapping[str, Callable[[], Quadrature]] = {
    name: (lambda name=name: posing_quadrature(GEOMETRIES[name]())) for name in GEOMETRIES
}


# ── The monomorphic population (moved from tests/gates/sn/architecture/test_monomorphic_leaves.py) ──

#: Deliberately non-uniform and geometrically graded: on a sphere ``V_cell`` spans 3.36e3 (MEASURED), so the bulk
#: metric ``V w`` is far from the constant that would make a reciprocity row a Euclidean-transpose check.
NONUNIFORM_EDGES = np.array([0.0, 0.12, 0.35, 0.80, 1.30, 2.00])


def fissile_anisotropic_2g(
    *, sig_t: "list[float]", sig_s0: "list[list[float]]", sig_s1: "list[list[float]]",
    sig_f: "list[float]", chi: "list[float]", nu: float = 2.6,
) -> Mixture:
    r"""The monomorphic population's fissile row of the ONE mixture builder (:func:`~tests.gates.transport._materials.anisotropic_mixture`).

    What this row adds over a within-group fixture is the channels the SN reciprocity gates (G1.4) need ACTIVATED
    and ``make_mixture`` nulls:

    * ``sig_f``/``chi``/``nu`` — producing, so ``F.apply(psi) != 0`` and the ``F`` reciprocity rows are not the
      tautology ``0 == 0``;
    * ``SigL != 0`` — the (n,alpha) channel;
    * an **asymmetric** ``Sig2 != 0`` — ``S``'s (n,2n) term, a separate ell=0 block outside the Legendre fold.

    The non-fissile default of ``anisotropic_mixture`` is correct for the *within-group* gates (fission enters as
    ``q_ext`` at the outer), which is exactly why this row states its extra channels rather than inherit them.
    """
    return anisotropic_mixture(
        sig_t, sig_s0, sig_s1, sig_f=sig_f, chi=chi, nu=nu,
        sig_l=[0.004, 0.011], sig_2=[[0.0, 0.03], [0.01, 0.0]],
    )


def two_region_fissile() -> "dict[int, Mixture]":
    """Two materials: heterogeneous nulls no spatial-distribution bug."""
    return {
        0: fissile_anisotropic_2g(sig_t=[1.1, 2.3], sig_s0=[[0.38, 0.10], [0.05, 0.90]],
                                  sig_s1=[[0.02, 0.01], [0.00, 0.04]], sig_f=[0.02, 0.31], chi=[0.95, 0.05]),
        1: fissile_anisotropic_2g(sig_t=[0.9, 2.6], sig_s0=[[0.22, 0.03], [0.12, 1.10]],
                                  sig_s1=[[0.05, 0.02], [0.01, 0.03]], sig_f=[0.05, 0.12], chi=[0.80, 0.20]),
    }


def two_region_mesh(coord: CoordSystem, edges: np.ndarray, boundaries: tuple) -> Mesh1D:
    """Materials 0 | 1 with the interface at ``edges[2]``: two cells, then three."""
    geometry = StructuredGeometry(
        coord=coord, breakpoints=(edges[0], edges[2], edges[-1]), mat_ids=(0, 1), face_roles=boundaries,
    )
    return Mesher(geometry).partition((CellEdges(edges[:3]), CellEdges(edges[2:]))).mesh


def cart2d_two_region_mesh() -> Mesh2D:
    """Non-square (4 x 3), non-uniform on both axes, materials 0 | 1 across x."""
    mat_map = np.zeros((4, 3), dtype=int)
    mat_map[2:, :] = 1
    return Mesh2D(
        edges_x=np.array([0.0, 0.1, 0.35, 0.9, 1.6]), edges_y=np.array([0.0, 0.2, 0.5, 1.4]),
        mat_map=mat_map, coord=CoordSystem.CARTESIAN,
        face_roles={"xmin": MirrorDeck(), "xmax": 0, "ymin": MirrorDeck(), "ymax": 0},
    )


def monomorphic() -> "dict[str, tuple[MaterialMesh, Quadrature]]":
    """The four monomorphic meshes with their quadratures (GL4 slab and sphere, folded product cylinder, LS4 rectangle)."""
    mats = two_region_fissile()
    return {
        "slab": (MaterialMesh(two_region_mesh(CoordSystem.CARTESIAN, NONUNIFORM_EDGES, (MirrorDeck(), 0)), mats),
                 Quadrature.gauss_legendre(n_ordinates=4)),
        "sphere": (MaterialMesh(two_region_mesh(CoordSystem.SPHERICAL, NONUNIFORM_EDGES, (1,)), mats),
                   Quadrature.gauss_legendre(n_ordinates=4)),
        "cylinder": (MaterialMesh(two_region_mesh(CoordSystem.CYLINDRICAL, NONUNIFORM_EDGES, (1,)), mats),
                     Quadrature.folded_product(n_mu=4, n_phi=8)),
        "cart2d": (MaterialMesh(cart2d_two_region_mesh(), mats), Quadrature.level_symmetric(sn_order=4)),
    }
