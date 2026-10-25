r"""The LD spatial-moment axis is the cell basis's continuum head (S3, #423).

Ruling (the user, 2026-10-25, option (b)): the scheme's spatial-moment axis is
``Axis.of_basis(spatial_moment_basis(p, d))``, origin ``FromBasis``. Its weights
are the basis's ``coefficient_weights``, the continuum Gram
:math:`G = \bigotimes_a \mathrm{diag}(1, 1/3)`, because the LD moments are
expansion COEFFICIENTS (:math:`\psi \approx \bar\psi + \xi\hat\psi`,
contravariant), and its unit is :math:`e_0`. No frame and no measure enter the
axis (the frame head, ``basis_space``, would install :math:`G^{-1}`: the
covariant-moment metric, wrong on the slope rows).

Kinds
-----
* REFERENCE rows compare the head with tables typed by hand from
  :math:`\theta = 1/3` (``1.0 / 3.0``, the spelling every LD consumer has used
  since #310), so they stay independent after S3 makes the scheme and the
  basis one object (``retirement-audit`` D.14/D.16: the old
  ``weights == moment_mass_diagonal`` rows DEMOTE to a self-comparison once
  the mass reads the Gram; these rows are their surviving pin).
* ROUTE rows (``test-architect`` lessons, the route gate): swap the basis's
  Gram for a decoy and require the scheme's θ and mass to MOVE.
* First reds on today's tree (``[M]`` 2026-10-25, HEAD ``57f2d777`` plus the
  uncommitted S1 code): ``test_the_ld_moment_axis_is_the_cell_basis_head``
  (origin ``Declared``), ``test_theta_and_the_mass_are_read_from_the_basis_gram``
  (θ a class literal), ``test_a_moment_basis_of_the_wrong_dimension_is_refused``
  (the float root accepts it).
"""

from __future__ import annotations

from itertools import product

import numpy as np
import pytest

from orpheus.geometry.deck import MirrorDeck
from orpheus.geometry.manifold import Interval, Product
from orpheus.mesh import AxisMesh
from orpheus.numerics.axis import Axis, Declared, FromBasis, Modal
from orpheus.numerics.basis import IntervalLegendreBasis, TensorProductBasis, spatial_moment_basis
from orpheus.transport.spatial import LinearDiscontinuous
from orpheus.transport.spatial.diamond import DiamondDifference

pytestmark = [pytest.mark.foundation]

_T = 1.0 / 3.0
#: The LD per-volume mass, typed by hand from theta = 1/3 in np.kron's product order (axis 0 outermost).
_HAND_MASS = {
    1: np.array([1.0, _T]),
    2: np.array([1.0, _T, _T, _T * _T]),
    3: np.array([1.0, _T, _T, _T * _T, _T, _T * _T, _T * _T, (_T * _T) * _T]),
}


def _ld_1d():
    # Inlined (not imported from test_spatial_moment_tail_is_the_schemes_axis.py, which imports moment_layout and
    # goes uncollectable mid-S4).
    from orpheus.geometry import StructuredGeometry
    from orpheus.mesh import CellsByCount, Mesher
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.sn.problem import SNProblem
    from orpheus.transport.boundary_conditions import WallCondition
    from tests.gates._exterior import conditions_for
    from tests.gates.sn._test_helpers import placeholder_materials

    mesh = Mesher(StructuredGeometry.slab((0.0, 1.0), (0,), left=0, right=0)).partition(CellsByCount.uniform_width(5)).mesh
    return SNProblem(mesh, Quadrature.gauss_legendre(4), placeholder_materials(ng=2), scheme=LinearDiscontinuous(),
                     boundary_conditions=conditions_for(mesh, {0: WallCondition.vacuum()}))


def _ld_2d():
    from orpheus.geometry import CoordSystem
    from orpheus.mesh import Mesh2D
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.sn.problem import SNProblem
    from orpheus.transport.boundary_conditions import WallCondition
    from tests.gates._exterior import conditions_for
    from tests.gates.sn._test_helpers import placeholder_materials

    mesh = Mesh2D(edges_x=np.linspace(0.0, 1.0, 4), edges_y=np.linspace(0.0, 1.0, 3), mat_map=np.zeros((3, 2), dtype=int),
                  coord=CoordSystem.CARTESIAN, face_roles={f: MirrorDeck() for f in ("xmin", "xmax", "ymin", "ymax")})
    return SNProblem(mesh, Quadrature.level_symmetric(4), placeholder_materials(ng=2), scheme=LinearDiscontinuous(),
                     boundary_conditions=conditions_for(mesh, {0: WallCondition.vacuum()}))


def _cart(d: int):
    return tuple(AxisMesh(edges=np.array([0.0, 1.0]), role_low=MirrorDeck(), role_high=MirrorDeck()) for _ in range(d))


@pytest.mark.parametrize("d", (1, 2, 3))
def test_the_ld_moment_axis_is_the_cell_basis_head(d: int) -> None:
    r"""Origin ``FromBasis(spatial_moment_basis(2, d))``, and the axis IS ``Axis.of_basis`` of it.

    First red: today's mint is ``Declared(Modal(e_0))`` (``scheme.py:1726-1740``)
    with ``support`` and ``layout`` both ``None``. The equality leg is content
    identity (label, shape, weights, kind, support, layout), so a head minted
    with the right origin and a hand-spelled support is still refused.
    """
    ax = LinearDiscontinuous().moment_axis(_cart(d))
    assert not isinstance(ax.origin, Declared), "the LD moment axis is still Declared: S3 has not landed"
    assert isinstance(ax.origin, FromBasis)
    assert ax.origin.basis == spatial_moment_basis(2, d)
    assert ax == Axis.of_basis(spatial_moment_basis(2, d))
    expected_support = Interval(-1.0, 1.0)
    for _ in range(d - 1):
        expected_support = Product(expected_support, Interval(-1.0, 1.0))
    assert ax.support == expected_support
    assert ax.layout is not None
    assert ax.layout == spatial_moment_basis(2, d).layout and ax.layout.label == "spatial_moment"


@pytest.mark.parametrize("d", (1, 2, 3))
def test_the_moment_axis_metric_is_the_gram_not_its_inverse(d: int) -> None:
    r"""``weights`` equal the hand-typed LD mass bit for bit, and are not the Parseval inverse.

    Bit-exact (``[M]`` 2026-10-25, ``scratch/layer_dag/moment_basis/ta/probe_gram_bits.out``):
    ``1.0/(2k+1)``, ``float(Fraction(1, 2k+1))``, the literal ``1.0/3.0`` and
    SymPy's exact Gram rounded once all give the same bits, kron'd in either
    accumulator order, for d = 1, 2, 3 (the per-mode exact product rounded once
    also agrees at p = 2). A quadrature-built Gram does NOT: Gauss-Legendre 4
    and 8 product rules give 5 to 22 ulp at d = 3, which is why the head reads
    the continuum Gram and no frame. Tooth: ``PROTO_ARM=inverse_gram``.
    """
    ax = LinearDiscontinuous().moment_axis(_cart(d))
    assert ax.weights is not None
    np.testing.assert_array_equal(ax.weights, _HAND_MASS[d])
    assert not np.allclose(ax.weights, 1.0 / _HAND_MASS[d]), "the head carries the inverse Gram (the frame's covariant metric)"


@pytest.mark.parametrize("d", (1, 2, 3))
def test_the_moment_axis_unit_is_the_constant_mode(d: int) -> None:
    r"""The unit is :math:`e_0` bitwise and the kind is ``Modal`` (``AVERAGE_MOMENT = 0`` as a basis fact)."""
    ax = LinearDiscontinuous().moment_axis(_cart(d))
    assert isinstance(ax.kind, Modal)
    np.testing.assert_array_equal(ax.unit, np.eye(2 ** d)[0])


def test_the_slopeless_closures_mint_no_moment_axis() -> None:
    r"""DD keeps its refusal: the width-1 head is a state no composed space admits (#240 D5b)."""
    for scheme in (DiamondDifference(),):
        assert scheme.spatial_basis_per_axis == 1
        with pytest.raises(NotImplementedError, match="no moment axis"):
            scheme.moment_axis(_cart(1))


def test_the_per_axis_count_is_the_basis_count() -> None:
    r"""``spatial_basis_per_axis`` is the count of the cell basis's factor (2 for LD)."""
    ax = LinearDiscontinuous().moment_axis(_cart(2))
    assert isinstance(ax.origin, FromBasis) and isinstance(ax.origin.basis, TensorProductBasis)
    assert all(f.n == LinearDiscontinuous.spatial_basis_per_axis for f in ax.origin.basis.factors)


def test_theta_and_the_mass_are_read_from_the_basis_gram(monkeypatch: pytest.MonkeyPatch) -> None:
    r"""ROUTE: a decoy interval Gram ``[1, 1/4]`` must move θ and the LD mass; the honest Gram restores them.

    Ruling 3 of the plan: θ is read from the Gram, one definition. The decoy is
    installed on the class, so every basis instance minted after the patch
    reads it; the activation leg (the honest values before the patch) is what
    makes a green reading mean "routed" rather than "cached". First red today:
    ``LinearDiscontinuous.theta`` is the class literal ``1.0 / 3.0``.
    ⚠ The row presumes θ is read at USE and the basis is minted per call. If S3
    computes θ once at class creation, or caches ``spatial_moment_basis``
    (a module constant, ``lru_cache``, ``instance_cache``), the class-level decoy
    never reaches the cached object and the row reds for the wrong reason
    (lessons L114): then install the decoy before the cache fills, or keep the
    retirement draft's theta-literal row as θ's only witness and say so here.
    """
    honest_mass = LinearDiscontinuous().moment_mass_diagonal(_cart(2))
    np.testing.assert_array_equal(honest_mass, _HAND_MASS[2])
    assert LinearDiscontinuous().theta == _T

    decoy = np.array([1.0, 0.25])
    decoy.setflags(write=False)
    monkeypatch.setattr(IntervalLegendreBasis, "gram_per_degree", property(lambda self: decoy[: self.n]))
    assert IntervalLegendreBasis(2).coefficient_weights[1] == 0.25, "precondition: the decoy did not install"
    assert LinearDiscontinuous().theta == 0.25, "theta is not read from the interval basis's Gram"
    np.testing.assert_array_equal(LinearDiscontinuous().moment_mass_diagonal(_cart(2)), np.kron(decoy, decoy))


@pytest.mark.parametrize("d", (1, 2))
def test_a_moment_basis_of_the_wrong_dimension_is_refused(d: int) -> None:
    r"""A tail axis whose basis has a factor count other than the mesh dimension is refused, not root-inverted.

    ``fields/_bases.py:533`` recovers the per-axis width as
    ``round(n ** (1/ndim))``. A one-factor basis of 4 polynomials on a 2-D mesh
    passes it (``4 ** (1/2) = 2``) and is read as LD's bilinear basis. After S3
    the width is read from the layout (``factor_shape``), so the mismatch is
    visible and must be refused. First red today: the root returns 2.
    """
    from orpheus.numerics.space import FunctionSpace
    from orpheus.transport.fields._bases import BulkField
    sn = {1: _ld_1d, 2: _ld_2d}[d]()
    # d = 2: one factor of 4 polynomials (the root reads 2); d = 1: two factors of 2 (the root reads 4).
    wrong = {2: TensorProductBasis((IntervalLegendreBasis(4),), "spatial_moment"),
             1: TensorProductBasis((IntervalLegendreBasis(2),) * 2, "spatial_moment")}[d]
    assert wrong.layout.shape == (4,) and len(wrong.factors) != d
    assert sn.bulk_space.axes is not None
    space = FunctionSpace.of_axes(*sn.bulk_space.axes, Axis.of_basis(wrong))
    with pytest.raises(ValueError, match="factor"):
        BulkField.spatial_moments_per_axis_of(space)


@pytest.mark.parametrize("ndim", (1, 2))
def test_the_face_transmission_matrix_is_sized_by_the_face_basis(ndim: int) -> None:
    r"""``_face_transmission_matrix`` is square of side ``ndim * dim(face basis)`` (``scheme.py:518``).

    The face basis is ``spatial_moment_basis(2, ndim).face(a)``, dimension
    ``2**(ndim-1)``; at ndim = 1 the face is a point and its count is 1.
    Tooth: ``SITE_ARM=face_count_518`` (the cell count where the face count belongs).
    ndim = 3 is out of reach: ``assemble_inflow_axis`` refuses the interior axis
    (``NotImplementedError``, deferred S2), ``[M]`` 2026-10-25.
    """
    from orpheus.transport.spatial.scheme import _face_transmission_matrix

    face_dim = 1 if ndim == 1 else spatial_moment_basis(2, ndim).face(0).layout.shape[0]
    sigma = _face_transmission_matrix(LinearDiscontinuous(), ndim, (0.7, 1.3, 0.4)[:ndim], 0.9)
    assert sigma.shape == (ndim * face_dim, ndim * face_dim)
