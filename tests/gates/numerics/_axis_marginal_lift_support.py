"""Fixtures and independent evaluators for ``test_axis_marginal_lift.py``.

The fixtures are built from production factories only. The evaluators never
call ORPHEUS's ``basis.evaluate`` (X4): real spherical harmonics in ORPHEUS's
no-prefactor convention (``orpheus.numerics.basis.spherical_harmonic_basis``'s
module docstring: polar axis mu_x, cos(phi) = mu_y / sin(theta), sin(phi) =
mu_z / sin(theta); slot [l, l + m]; m > 0 -> cos(m phi), m < 0 -> sin(|m| phi);
Schmidt semi-normalised, no Condon-Shortley phase), built on
``scipy.special.lpmv`` with its (-1)^m phase removed; Legendre from
``numpy.polynomial.legendre``; a 40 x 80 product rule on the sphere built here.
"""
from __future__ import annotations


from functools import cache
from math import factorial

import numpy as np
from numpy.polynomial import legendre as npleg
from scipy.special import lpmv

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, CoordSystem, StructuredGeometry
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics import frame_faces
from orpheus.numerics.axis import Declared, Axis, BasisKind, EnergyAxis
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.space import FunctionSpace


def _mats():
    return {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")}


@cache
def problems():
    from orpheus.sn.solver import _as_problem
    from orpheus.mesh.structured import Mesh2D
    from orpheus.transport.spatial import LinearDiscontinuous

    slab = Mesher(StructuredGeometry.slab(
        (0.0, 0.37, 1.9, 4.13), (0, 1, 0), left=BC.vacuum, right=BC.reflective,
    )).partition((CellsByCount.uniform_width(3), CellsByCount.uniform_width(5), CellsByCount.uniform_width(7))).mesh
    sph = Mesher(StructuredGeometry.sphere((0.0, 1.5, 3.0), (0, 1), outer=BC("vacuum"))).partition(
        CellsByCount.uniform_width(3)).mesh
    cyl = Mesher(StructuredGeometry.cylinder((0.0, 1.0, 2.5), (0, 1), outer=BC.vacuum)).partition(
        (CellsByCount.uniform_width(3), CellsByCount.uniform_width(4))).mesh
    nx, ny = 4, 5
    mat_map = np.zeros((nx, ny), dtype=int)
    mat_map[nx // 2:, :] = 1
    xy = Mesh2D(
        edges_x=np.linspace(0.0, 2.0, nx + 1), edges_y=np.linspace(0.0, 3.0, ny + 1),
        mat_map=mat_map, coord=CoordSystem.CARTESIAN,
        face_laws={"xmin": BC("reflective"), "xmax": BC("vacuum"), "ymin": BC("reflective"), "ymax": BC("vacuum")},
    )
    m = _mats()
    return {
        "slab GL8": _as_problem(slab, Quadrature.gauss_legendre(n_ordinates=8), m),
        "slab GL8 LD": _as_problem(slab, Quadrature.gauss_legendre(n_ordinates=8), m, scheme=LinearDiscontinuous()),
        "sphere GL8": _as_problem(sph, Quadrature.gauss_legendre(n_ordinates=8), m),
        "cylinder folded(4,8)": _as_problem(cyl, Quadrature.folded_product(4, 8), m),
        "xy LS4": _as_problem(xy, Quadrature.level_symmetric(sn_order=4), m),
        "xy product(4,8)": _as_problem(xy, Quadrature.product(4, 8), m),
    }


_W_SPA = np.array([0.2, 0.3, 0.4, 0.7, 1.4])


def _cells():
    return (EnergyAxis.synthetic(2), Axis("spatial", (5,), weights=_W_SPA, origin=Declared(BasisKind.NODAL)))


def head_products():
    """Frame-dressed and continuum heads, each times energy(2) x spatial(5)."""
    from orpheus.numerics.spaces.legendre_space import LegendreSpace
    from orpheus.numerics.spaces.spherical_harmonic_space import SphericalHarmonicSpace

    cells = FunctionSpace.of_axes(*_cells())
    heads = {
        "folded(4,4) frame L2": frame_faces.basis_space(frame_faces.angular_frame(Quadrature.folded_product(4, 4), 2)),
        "folded(4,4) frame L4": frame_faces.basis_space(frame_faces.angular_frame(Quadrature.folded_product(4, 4), 4)),
        "product(4,8) frame L4": frame_faces.basis_space(frame_faces.angular_frame(Quadrature.product(4, 8), 4)),
        "lebedev7 frame L3": frame_faces.basis_space(frame_faces.angular_frame(Quadrature.lebedev(order=7), 3)),
        "GL8 frame L3": frame_faces.basis_space(frame_faces.angular_frame(Quadrature.gauss_legendre(n_ordinates=8), 3)),
        "continuum SH L3": SphericalHarmonicSpace.from_L(3),
        "continuum Legendre L3 z": LegendreSpace.from_L(3, "z"),
    }
    return {k: h * cells for k, h in heads.items()} | {f"{k} [head alone]": h for k, h in heads.items()}


def synthetic():
    return {
        "synthetic nodal": FunctionSpace.of_axes(
            Axis("angular", (4,), weights=np.array([0.3, 0.7, 0.5, 0.5]), origin=Declared(BasisKind.NODAL)),
            Axis("energy", (2,), origin=Declared(BasisKind.NODAL)),
            Axis("spatial", (5,), weights=_W_SPA, origin=Declared(BasisKind.NODAL)),
        ),
        "zero-mass signed": FunctionSpace.of_axes(
            Axis("angular", (2,), weights=np.array([1.0, -1.0]), origin=Declared(BasisKind.NODAL)),
            Axis("spatial", (3,), origin=Declared(BasisKind.NODAL)),
        ),
    }


def zero_weight_nodal():
    """A nodal axis with one zero weight (legal: signed and zero weights are admitted)."""
    return FunctionSpace.of_axes(
        Axis("angular", (3,), weights=np.array([0.5, 0.0, 1.5]), origin=Declared(BasisKind.NODAL)),
        Axis("spatial", (4,), weights=np.array([1.0, 2.0, 0.5, 0.25]), origin=Declared(BasisKind.NODAL)),
    )


def _spd(n, seed):
    g = np.random.default_rng(seed).standard_normal((n, n))
    return g @ g.T + 3.0 * np.eye(n)


def dense_form_spaces():
    """D1, D2, D3 of tests/gates/numerics/test_retraction_adjoint_is_the_pullback.py, verbatim."""
    from orpheus.numerics.metric import DenseMetric, FactoredMetric

    ang = Axis("angular", (4,), weights=np.array([0.3, 0.9, 0.5, 0.3]), origin=Declared(BasisKind.NODAL))
    spa = Axis("spatial", (5,), weights=np.array([1.0, 2.0, 0.5, 0.25, 3.0]), origin=Declared(BasisKind.NODAL))
    mom = Axis("moment", (3,), origin=Declared(BasisKind.NODAL))
    d1 = FunctionSpace("d1", (4, 5, 3), axes=(ang, spa, mom),
                       metric=FactoredMetric((((4,), None), ((5,), None), ((3,), DenseMetric(_spd(3, 1))))))
    d2 = FunctionSpace("d2", (4, 5), axes=(Axis("angular", (4,), origin=Declared(BasisKind.NODAL)), spa),
                       metric=FactoredMetric((((4,), DenseMetric(_spd(4, 2))), ((5,), None))))
    head = FunctionSpace("h", (3,), axes=(mom,), metric=FactoredMetric((((3,), DenseMetric(_spd(3, 3))),)))
    d3 = FunctionSpace.of_axes(ang, spa) * head
    return {"D1": d1, "D2": d2, "D3": d3}


@cache
def registry():
    """name -> space, every member of the S1 population (resolved by NAME at run time, L100)."""
    reg = {}
    for pn, p in problems().items():
        reg[f"{pn} angular_bulk"] = p.angular_bulk_space
        for L in (1, 2):
            reg[f"{pn} moment L{L}"] = p.moment_space(L)
    reg["slab GL8 LD angular_trial"] = problems()["slab GL8 LD"].angular_trial_space
    reg["slab GL8 bulk"] = problems()["slab GL8"].bulk_space
    reg["xy LS4 bulk"] = problems()["xy LS4"].bulk_space
    reg |= {k: v for k, v in head_products().items() if "[head alone]" not in k}
    reg |= synthetic()
    reg["zero-weight nodal"] = zero_weight_nodal()
    reg |= dense_form_spaces()
    return reg


#: (space name, axis label) rows over a NODAL axis.
NODAL_ROWS = (
    [(f"{p} angular_bulk", "angular") for p in ("slab GL8", "slab GL8 LD", "sphere GL8", "cylinder folded(4,8)", "xy LS4", "xy product(4,8)")]
    + [("slab GL8 LD angular_trial", "angular"), ("slab GL8 bulk", "spatial"), ("xy LS4 bulk", "spatial"),
       ("synthetic nodal", "angular"), ("synthetic nodal", "energy"), ("synthetic nodal", "spatial"),
       ("zero-weight nodal", "angular"), ("D1", "angular"), ("D2", "angular"), ("D3", "angular"),
       ("folded(4,4) frame L4", "spatial")]
)
#: Production moment spaces: the head is a frame-dressed MOMENT head (covariant coordinates).
MOMENT_SPACE_ROWS = [(f"{p} moment L{L}", None) for p in ("slab GL8", "slab GL8 LD", "sphere GL8", "cylinder folded(4,8)", "xy LS4", "xy product(4,8)") for L in (1, 2)]
#: Frame-dressed heads x cells (dense positioned forms on three of them).
FRAME_HEAD_ROWS = [(k, None) for k in ("folded(4,4) frame L2", "folded(4,4) frame L4", "product(4,8) frame L4", "lebedev7 frame L3", "GL8 frame L3")]
#: Basis-minted (continuum) heads x cells: COEFFICIENT heads (contravariant coordinates).
COEFFICIENT_HEAD_ROWS = [(k, None) for k in ("continuum SH L3", "continuum Legendre L3 z")]
MODAL_ROWS = MOMENT_SPACE_ROWS + FRAME_HEAD_ROWS + COEFFICIENT_HEAD_ROWS


def resolve(name, label):
    V = registry()[name]
    return V, (V.axes[0].label if label is None else label)




def real_sh_table(L: int, directions: np.ndarray) -> np.ndarray:
    d = np.asarray(directions, float)
    mx, my, mz = d[:, 0], d[:, 1], d[:, 2]
    st = np.sqrt(np.maximum(1.0 - mx * mx, 0.0))
    phi = np.arctan2(mz, my)  # azimuth about mu_x, measured from mu_y toward mu_z
    out = np.zeros((d.shape[0], L + 1, 2 * L + 1))
    for l in range(L + 1):
        for m in range(-l, l + 1):
            am = abs(m)
            norm = 1.0 if m == 0 else np.sqrt(2.0 * factorial(l - am) / factorial(l + am))
            plm = ((-1.0) ** am) * lpmv(am, l, mx)  # strip Condon-Shortley
            ang = 1.0 if m == 0 else (np.cos(am * phi) if m > 0 else np.sin(am * phi))
            out[:, l, l + m] = norm * plm * ang
    return out


def legendre_table(L: int, mu: np.ndarray) -> np.ndarray:
    return npleg.legvander(np.asarray(mu, float), L)


def sphere_rule(n_mu: int = 40, n_phi: int = 80):
    """A high-order product rule on S^2 built here (GL in mu_x x uniform phi)."""
    x, wx = npleg.leggauss(n_mu)
    ph = (np.arange(n_phi) + 0.5) * 2 * np.pi / n_phi
    X, P = np.meshgrid(x, ph, indexing="ij")
    st = np.sqrt(1 - X ** 2)
    dirs = np.stack([X.ravel(), (st * np.cos(P)).ravel(), (st * np.sin(P)).ravel()], axis=1)
    w = (wx[:, None] * np.full(n_phi, 2 * np.pi / n_phi)[None, :]).ravel()
    return dirs, w
