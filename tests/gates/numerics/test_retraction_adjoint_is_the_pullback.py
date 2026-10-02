r"""The retraction's Hilbert adjoint is the pullback (#405 P1 step 6, S6.1-S6.4).

The retraction :math:`R = \pi_*` integrates over one axis,
:math:`(R\psi)(\cdot) = \sum_n w_n \psi(n, \cdot)`. When the full space's
metric is the product (the axis measure) ⊗ (the marginal metric), its
Hilbert adjoint is the pullback :math:`\pi^*`, the plain broadcast, with no
metric in it:
:math:`\langle R\psi, \phi\rangle = \sum_n w_n \langle\psi_n, \phi\rangle
= \langle\psi, \pi^*\phi\rangle`.

Until step 6, ``R.H`` was the generic metric sandwich ♯∘Rᵀ∘♭, the same
operator in exact arithmetic, but rounding 1-3 ULP away from the broadcast
(``[M]`` 2026-10-02: ``array_equal`` failed 1 178 of 1 200 draws over the
six spaces below). The mint now returns the closed form, and it admits only
spaces on which the product holds: it refuses an axis carrying a positioned
form, and the marginal keeps the forms of the axes that remain (before step
6 the marginal dropped them, the spaces D1 and D3 below).

The rows:

- S6.1: ``R.H`` is the pullback type, ``R.H.H is R``, the spaces swap, and
  ``R.H.apply`` is ``array_equal`` to the broadcast.
- S6.2: the product condition and reciprocity, the licence of the closed
  form, with the weights read from the raw axes, never from ``R``.
- S6.3: the mint refuses the one space family on which the product fails
  (a form on the collapsed axis, D2), so ``.H`` needs no guard.
- S6.4: the X4 witness: the closed form agrees with an explicitly built
  sandwich ``AdjointOperator(R)`` to ``nulp=4`` (measured worst 3 ULP).

The specification is ``.claude/plans/reference_p1_spec.md`` §1.6.
"""

from __future__ import annotations

from functools import cache

import numpy as np
import numpy.testing as npt
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, StructuredGeometry
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.axis import Axis, BasisKind
from orpheus.numerics.metric import DenseMetric, FactoredMetric
from orpheus.numerics.operator import (
    AdjointOperator,
    AxisPullbackOperator,
    AxisRetractionOperator,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.space import FunctionSpace
from orpheus.sn.solver import _as_problem
from orpheus.transport.spatial import LinearDiscontinuous

pytestmark = pytest.mark.foundation

_DRAWS = 200


def _spd(n: int, seed: int) -> np.ndarray:
    g = np.random.default_rng(seed).standard_normal((n, n))
    return g @ g.T + 3.0 * np.eye(n)


_W_ANG = np.array([0.3, 0.9, 0.5, 0.3])
_W_SPA = np.array([1.0, 2.0, 0.5, 0.25, 3.0])
_ANG = Axis("angular", (4,), weights=_W_ANG, kind=BasisKind.NODAL)
_SPA = Axis("spatial", (5,), weights=_W_SPA, kind=BasisKind.NODAL)
_MOM = Axis("moment", (3,), kind=BasisKind.NODAL)


def _d1() -> FunctionSpace:
    """A dense form on a NON-collapsed, measure-less axis."""
    return FunctionSpace(
        "d1", (4, 5, 3), axes=(_ANG, _SPA, _MOM),
        metric=FactoredMetric((((4,), None), ((5,), None), ((3,), DenseMetric(_spd(3, 1))))),
    )


def _d2() -> FunctionSpace:
    """A dense form on the COLLAPSED axis itself (a measure-less angular axis)."""
    return FunctionSpace(
        "d2", (4, 5), axes=(Axis("angular", (4,), kind=BasisKind.NODAL), _SPA),
        metric=FactoredMetric((((4,), DenseMetric(_spd(4, 2))), ((5,), None))),
    )


def _d3() -> FunctionSpace:
    """The production composition route: ``(angular ⊗ spatial) * dense head``."""
    head = FunctionSpace("h", (3,), axes=(_MOM,), metric=FactoredMetric((((3,), DenseMetric(_spd(3, 3))),)))
    return FunctionSpace.of_axes(_ANG, _SPA) * head


@cache
def _population() -> dict[str, FunctionSpace]:
    """P_π: the synthetic product, the SN bulk spaces of a 3-region non-uniform
    slab at Gauss-Legendre 5, 8 and 16, that slab's linear-discontinuous trial
    space (it carries the moment axis), a 2-region cylinder's bulk space at
    ``folded_product(4, 8)``, and the dense-form spaces D1 and D3."""
    mats = {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")}
    slab = Mesher(StructuredGeometry.slab(
        (0.0, 0.37, 1.9, 4.13), (0, 1, 0), left=BC.vacuum, right=BC.reflective,
    )).partition((CellsByCount.uniform_width(3), CellsByCount.uniform_width(5), CellsByCount.uniform_width(7))).mesh
    cyl = Mesher(StructuredGeometry.cylinder((0.0, 1.0, 2.5), (0, 1), outer=BC.vacuum)).partition(
        (CellsByCount.uniform_width(3), CellsByCount.uniform_width(4)),
    ).mesh
    spaces = {
        "synthetic": FunctionSpace.of_axes(
            Axis("angular", (4,), weights=np.array([0.2, 0.8, 0.8, 0.2]), kind=BasisKind.NODAL),
            Axis("energy", (2,), kind=BasisKind.NODAL),
            Axis("spatial", (5,), weights=np.array([0.2, 0.3, 0.4, 0.7, 1.4]), kind=BasisKind.NODAL),
        ),
    }
    for n in (5, 8, 16):
        spaces[f"slab GL{n} bulk"] = _as_problem(slab, Quadrature.gauss_legendre(n_ordinates=n), mats).angular_bulk_space
    spaces["slab GL8 LD trial"] = _as_problem(
        slab, Quadrature.gauss_legendre(n_ordinates=8), mats, scheme=LinearDiscontinuous(),
    ).angular_trial_space
    spaces["cylinder folded(4,8) bulk"] = _as_problem(cyl, Quadrature.folded_product(4, 8), mats).angular_bulk_space
    spaces["D1 dense form on a kept axis"] = _d1()
    spaces["D3 (angular x spatial) * dense head"] = _d3()
    return spaces


_NAMES = (
    "synthetic", "slab GL5 bulk", "slab GL8 bulk", "slab GL16 bulk", "slab GL8 LD trial",
    "cylinder folded(4,8) bulk", "D1 dense form on a kept axis", "D3 (angular x spatial) * dense head",
)


def _angular_weights(space: FunctionSpace) -> np.ndarray:
    """The collapsed axis's measure, read from the raw axis (never from R)."""
    assert space.axes is not None  # narrowing only: every member is axis-built
    (axis,) = (ax for ax in space.axes if ax.label == "angular")
    return np.ones(int(np.prod(axis.shape))) if axis.weights is None else np.asarray(axis.weights, float).ravel()


def _angular_dim(space: FunctionSpace) -> int:
    """The ndarray dim the angular axis occupies (every member's is rank 1)."""
    assert space.axes is not None  # narrowing only
    labels = [ax.label for ax in space.axes]
    return sum(len(ax.shape) for ax in space.axes[: labels.index("angular")])


def _phis(R: AxisRetractionOperator, seed: int):
    rng = np.random.default_rng(seed)
    for _ in range(_DRAWS):
        yield rng.standard_normal(R.codomain.shape) * 10.0 ** rng.integers(-6, 7)


@pytest.mark.parametrize("name", _NAMES)
def test_s6_1_the_adjoint_is_the_pullback(name: str) -> None:
    """S6.1. First red (2026-10-02): ``R.H`` was ``AdjointOperator(R)``, and
    ``array_equal`` with the broadcast failed 1 178 of 1 200 draws."""
    V = _population()[name]
    R = V.retraction("angular")
    P = R.H
    if not isinstance(P, AxisPullbackOperator):
        pytest.fail(f"{name}: R.H is a {type(P).__name__}, not the pullback")
    if P.H is not R:
        pytest.fail(f"{name}: R.H.H is not R (the involution is an object identity)")
    if P.domain is not R.codomain or P.codomain is not R.domain:
        pytest.fail(f"{name}: the pullback does not swap the retraction's spaces")
    angular_dim = _angular_dim(V)
    for phi in _phis(R, seed=11):
        broadcast = np.broadcast_to(np.expand_dims(phi, angular_dim), V.shape)
        npt.assert_array_equal(P.apply(phi), broadcast)


@pytest.mark.parametrize("name", _NAMES)
def test_s6_2_the_product_metric_licenses_the_closed_form(name: str) -> None:
    """S6.2. The full metric is (the axis measure) ⊗ (the marginal metric),
    and reciprocity ⟨Rψ, φ⟩ = ⟨ψ, R.Hφ⟩ holds. D1 and D3 failed the product
    by 0.80 and 0.79 (relative) before the mint kept the marginal's forms."""
    V = _population()[name]
    R = V.retraction("angular")
    M = R.codomain
    w = _angular_weights(V)
    rng = np.random.default_rng(29)
    x = rng.standard_normal(V.shape)
    full = np.asarray(V.apply_metric(x))
    angular_dim = _angular_dim(V)
    moved = np.moveaxis(x, angular_dim, 0)
    product = np.moveaxis(
        np.stack([w_n * np.asarray(M.apply_metric(moved[n])) for n, w_n in enumerate(w)]), 0, angular_dim,
    )
    npt.assert_allclose(full, product, rtol=1e-14, atol=1e-14 * np.max(np.abs(full)))
    psi, phi = rng.standard_normal(V.shape), rng.standard_normal(M.shape)
    npt.assert_allclose(
        M.inner_product(np.asarray(R.apply(psi)), phi), V.inner_product(psi, np.asarray(R.H.apply(phi))), rtol=1e-13,
    )


def test_s6_3_the_mint_refuses_a_form_on_the_collapsed_axis() -> None:
    """S6.3. On D2 the retraction integrates with the axis's counting measure
    while the space's metric on that block is a dense form: two measures on
    one block. Before step 6 the mint admitted it and ``R.H`` was the
    sandwich; the pullback would violate reciprocity there by 0.19."""
    with pytest.raises(ValueError, match="two measures on one block"):
        _d2().retraction("angular")


def test_s6_3_positive_leg_the_marginal_keeps_the_kept_forms() -> None:
    """The marginal of D1 carries D1's dense form on the moment block."""
    M = _d1().retraction("angular").codomain
    if not isinstance(M.metric, FactoredMetric):
        pytest.fail("the marginal dropped the kept axes' positioned forms")
    forms = [form for _, form in M.metric.entries]
    if not (forms[0] is None and isinstance(forms[1], DenseMetric)):
        pytest.fail(f"the marginal's forms are {forms}, not (None, DenseMetric)")


def _dense_condition(space: FunctionSpace) -> float:
    """The worst condition number of the space's positioned dense forms, or
    1.0 when it carries none."""
    if not isinstance(space.metric, FactoredMetric):
        return 1.0
    return max(
        (float(np.linalg.cond(form.matrix)) for _, form in space.metric.entries if isinstance(form, DenseMetric)),
        default=1.0,
    )


@pytest.mark.parametrize("name", _NAMES)
def test_s6_4_the_sandwich_witnesses_the_closed_form(name: str) -> None:
    """S6.4, the X4 witness: π* agrees with the generic sandwich ♯∘Rᵀ∘♭,
    built explicitly here (production no longer routes there).

    On a diagonal metric the two differ only by the re-association
    (w·G·φ)/(w·G): ``nulp=4``, measured worst 3 ULP on the LD trial space and
    2 on the others. On a dense form the sandwich's ♯ applies the form's
    pseudo-inverse, whose rounding is κ(G)·ε relative to the result (κ the
    form's condition number); the bound there is 8·κ·ε, derived from κ,
    which the row computes from the form itself. Designed green: π* := Σw·E
    is the same operator and agrees to 0-1 ULP, so this row cannot tell them
    apart (S6.1 can)."""
    V = _population()[name]
    R = V.retraction("angular")
    sandwich = AdjointOperator(R)  # pyright: ignore[reportArgumentType]
    kappa = _dense_condition(V)
    for phi in _phis(R, seed=13):
        closed, generic = np.asarray(R.H.apply(phi)), np.asarray(sandwich.apply(phi))
        if kappa == 1.0:
            npt.assert_array_almost_equal_nulp(closed, generic, nulp=4)
        else:
            npt.assert_allclose(closed, generic, rtol=0, atol=8 * kappa * np.finfo(float).eps * np.max(np.abs(closed)))
