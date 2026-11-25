r"""The retraction's Hilbert adjoint is the pullback (#405 P1 step 6, S6.1-S6.4).

The retraction :math:`R = \pi_*` integrates over one axis against the
covector :math:`m = G_k u` (the axis's metric block applied to its unit; on
the nodal axes here :math:`u = 1`, so :math:`m` is the weights, or the
positioned form's row sums on D2). When the full space's metric is the
product :math:`G_k \otimes G_{\rm rest}`, its Hilbert adjoint is the
pullback :math:`\pi^* = u \otimes \mathrm{id}`, the plain broadcast here,
with no metric in it:
:math:`\langle R\psi, \phi\rangle = \langle\psi, \pi^*\phi\rangle`.

Until step 6, ``R.H`` was the generic metric sandwich ♯∘Rᵀ∘♭, the same
operator in exact arithmetic, but rounding 1-3 ULP away from the broadcast
(``[M]`` 2026-10-02: ``array_equal`` failed 1 178 of 1 200 draws over the
six spaces below). The mint now returns the closed form. Every axis-built
space carries the product metric, and the marginal keeps the forms of the
axes that remain (before step 6 it dropped them, the spaces D1 and D3
below). Until S1 the mint refused a positioned form on the collapsed axis
(D2, "two measures on one block"); since S1 that form IS :math:`G_k`, the
covector is built with it, and D2 is a member of the population.

The rows:

- S6.1: ``R.H`` is the pullback type, ``R.H.H is R``, the spaces swap, and
  ``R.H.apply`` is ``array_equal`` to the broadcast.
- S6.2: the product condition and reciprocity, the licence of the closed
  form, with the axis's metric block read from the raw axes and forms,
  never from ``R``.
- S6.3: the marginal keeps the kept axes' forms.
- S6.4: the X4 witness: the closed form agrees with an explicitly built
  sandwich ``AdjointOperator(R)`` to ``nulp=4`` (measured worst 3 ULP).

The specification is ``.claude/plans/reference_p1_spec.md`` §1.6.
"""

from __future__ import annotations

from functools import cache

import numpy as np
import numpy.testing as npt
import pytest

from orpheus.numerics.axis import Declared, Axis, BasisKind, Modal
from orpheus.numerics.metric import DenseMetric, Face, FactoredMetric, NoForm
from orpheus.numerics.axis_marginal import AxisPullbackOperator, AxisRetractionOperator
from orpheus.numerics.operator import AdjointOperator
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.space import FunctionSpace
from tests.gates.numerics._phase_space_twins import spatial_axis_1d, twin
from orpheus.numerics.factor import Factor, Kind, Subquotient, Variance, VarianceChart
from orpheus.numerics.space import TensorSpace
from orpheus.numerics.metric import DiagonalMetric
from orpheus.numerics.declaration import Coord


pytestmark = pytest.mark.foundation

_DRAWS = 200


def _spd(n: int, seed: int) -> np.ndarray:
    g = np.random.default_rng(seed).standard_normal((n, n))
    return g @ g.T + 3.0 * np.eye(n)


_W_ANG = np.array([0.3, 0.9, 0.5, 0.3])
_W_SPA = np.array([1.0, 2.0, 0.5, 0.25, 3.0])
_ANG = Axis("angular", (4,), origin=Declared(BasisKind.NODAL, DiagonalMetric(_W_ANG)))
_SPA = Axis("spatial", (5,), origin=Declared(BasisKind.NODAL, DiagonalMetric(_W_SPA)))
_MD = Variance(Kind.MEASURE, VarianceChart.DENSITY)


def _dense(n: int, seed: int, label: str) -> Axis:
    """A coefficient (Modal, unit e_0) axis whose form is a dense SPD Gram."""
    return Axis(label, (n,), origin=Declared(Modal(np.eye(n)[0]), DenseMetric(_spd(n, seed), kernel=Coord(n, ()), face=Face.METRIC)))


def _whole(axis: Axis) -> Factor:
    return Factor(axis, _MD, Subquotient.whole(int(np.prod(axis.shape))))


def _d1() -> FunctionSpace:
    """A dense form on a NON-collapsed coefficient (Modal) axis."""
    return TensorSpace.of_axes(_whole(_ANG), _whole(_SPA), _whole(_dense(3, 1, "moment")))


def _d2() -> FunctionSpace:
    """A dense form on the COLLAPSED axis itself (a coefficient axis, Modal with unit e_0: a form lives only there)."""
    return TensorSpace.of_axes(_whole(_dense(4, 2, "angular")), _whole(_SPA))


def _d3() -> FunctionSpace:
    """The production composition route: ``(angular ⊗ spatial) * dense head``."""
    head = TensorSpace.of_axes(_whole(_dense(3, 3, "moment")))
    return TensorSpace.of_axes(Factor.whole_on(_ANG, Kind.MEASURE, VarianceChart.DENSITY), Factor.whole_on(_SPA, Kind.MEASURE, VarianceChart.DENSITY)) * head


@cache
def _population() -> dict[str, FunctionSpace]:
    """P_π: the synthetic product, the angular bulk spaces of a 3-region non-uniform
    slab at Gauss-Legendre 5, 8 and 16, that slab's linear-discontinuous trial
    space (it carries the moment axis), a 2-region cylinder's bulk space at
    ``folded_product(4, 8)``, and the dense-form spaces D1, D2 and D3."""
    # The production-shaped spaces as numerics twins (layer DAG Phase 4; until then read off ``_as_problem`` SN
    # carriers): the same 3-region slab (3 + 5 + 7 cells on 0 | 0.37 | 1.9 | 4.13) and 2-region cylinder (3 + 4 on
    # 0 | 1 | 2.5), two groups.
    slab = spatial_axis_1d("slab", np.concatenate(
        [np.linspace(0.0, 0.37, 4), np.linspace(0.37, 1.9, 6)[1:], np.linspace(1.9, 4.13, 8)[1:]],
    ))
    cyl = spatial_axis_1d("cylinder", np.concatenate([np.linspace(0.0, 1.0, 4), np.linspace(1.0, 2.5, 5)[1:]]))
    spaces = {
        "synthetic": TensorSpace.of_axes(Factor.whole_on(Axis("angular", (4,), origin=Declared(BasisKind.NODAL, DiagonalMetric(np.array([0.2, 0.8, 0.8, 0.2])))), Kind.MEASURE, VarianceChart.DENSITY), Factor.whole_on(Axis("energy", (2,), origin=Declared(BasisKind.NODAL)), Kind.MEASURE, VarianceChart.DENSITY), Factor.whole_on(Axis("spatial", (5,), origin=Declared(BasisKind.NODAL, DiagonalMetric(np.array([0.2, 0.3, 0.4, 0.7, 1.4])))), Kind.MEASURE, VarianceChart.DENSITY)),
    }
    for n in (5, 8, 16):
        spaces[f"slab GL{n} bulk"] = twin(Quadrature.gauss_legendre(n_ordinates=n), slab).angular_bulk
    spaces["slab GL8 LD trial"] = twin(Quadrature.gauss_legendre(n_ordinates=8), slab).ld_trial()
    spaces["cylinder folded(4,8) bulk"] = twin(Quadrature.folded_product(4, 8), cyl).angular_bulk
    spaces["D1 dense form on a kept axis"] = _d1()
    spaces["D2 dense form on the collapsed axis"] = _d2()
    spaces["D3 (angular x spatial) * dense head"] = _d3()
    return spaces


_NAMES = (
    "synthetic", "slab GL5 bulk", "slab GL8 bulk", "slab GL16 bulk", "slab GL8 LD trial",
    "cylinder folded(4,8) bulk", "D1 dense form on a kept axis", "D2 dense form on the collapsed axis",
    "D3 (angular x spatial) * dense head",
)


def _angular_block(space: FunctionSpace) -> np.ndarray:
    """The collapsed axis's metric block G_k as a matrix, read from the raw
    axis and the space's positioned forms (never from R): the form on the
    axis's block if there is one, else diag(weights), else the identity."""
    labels = [ax.label for ax in space.axes]
    k = labels.index("angular")
    axis = space.axes[k]
    form = space.factors[k].form
    if isinstance(form, DenseMetric):
        return np.asarray(form.matrix, float)
    n = int(np.prod(axis.shape))
    return np.eye(n) if isinstance(axis.form, NoForm) else np.diag(np.asarray(axis.form.weights, float).ravel())


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
    ``array_equal`` with the broadcast failed 1 178 of 1 200 draws. The oracle
    is ``u ⊗ φ``, which is the plain broadcast on every nodal row."""
    V = _population()[name]
    R = AxisRetractionOperator.of(V, "angular")
    P = R.H
    if not isinstance(P, AxisPullbackOperator):
        pytest.fail(f"{name}: R.H is a {type(P).__name__}, not the pullback")
    if P.H is not R:
        pytest.fail(f"{name}: R.H.H is not R (the involution is an object identity)")
    if P.domain is not R.codomain or P.codomain is not R.domain:
        pytest.fail(f"{name}: the pullback does not swap the retraction's spaces")
    angular_dim = _angular_dim(V)
    # The pullback is u ⊗ φ, with u the collapsed axis's unit: ones on a nodal
    # axis (the broadcast, bit for bit, since x * 1.0 == x), e_0 on D2's
    # coefficient axis (Modal since U2, plan space_membership_covector_gauge.md).
    unit = np.asarray(V.axis("angular").unit)
    unit_along = unit.reshape([n if d == angular_dim else 1 for d, n in enumerate(V.shape)])
    for phi in _phis(R, seed=11):
        broadcast = np.broadcast_to(np.expand_dims(phi, angular_dim), V.shape) * unit_along
        npt.assert_array_equal(P.apply(phi), broadcast)


@pytest.mark.parametrize("name", _NAMES)
def test_s6_2_the_product_metric_licenses_the_closed_form(name: str) -> None:
    """S6.2. The full metric is G_k ⊗ (the marginal metric), and reciprocity
    ⟨Rψ, φ⟩ = ⟨ψ, R.Hφ⟩ holds. D1 and D3 failed the product by 0.80 and 0.79
    (relative) before the mint kept the marginal's forms."""
    V = _population()[name]
    R = AxisRetractionOperator.of(V, "angular")
    M = R.codomain
    G_k = _angular_block(V)
    rng = np.random.default_rng(29)
    x = rng.standard_normal(V.shape)
    full = np.asarray(V.apply_metric(x))
    angular_dim = _angular_dim(V)
    moved = np.moveaxis(x, angular_dim, 0)
    marginal_metric = np.stack([np.asarray(M.apply_metric(moved[n])) for n in range(G_k.shape[0])])
    product = np.moveaxis(np.einsum("nk,k...->n...", G_k, marginal_metric), 0, angular_dim)
    npt.assert_allclose(full, product, rtol=1e-14, atol=1e-14 * np.max(np.abs(full)))
    psi, phi = rng.standard_normal(V.shape), rng.standard_normal(M.shape)
    npt.assert_allclose(
        M.inner_product(np.asarray(R.apply(psi)), phi), V.inner_product(psi, np.asarray(R.H.apply(phi))), rtol=1e-13,
    )


def test_s6_3_positive_leg_the_marginal_keeps_the_kept_forms() -> None:
    """The marginal of D1 carries D1's dense form on the moment block."""
    M = AxisRetractionOperator.of(_d1(), "angular").codomain
    forms = [f.form for f in M.factors]
    if not (isinstance(forms[0], DiagonalMetric) and isinstance(forms[1], DenseMetric)):
        pytest.fail(f"the marginal's forms are {forms}, not (DiagonalMetric, DenseMetric)")


def _dense_condition(space: FunctionSpace) -> float:
    """The worst condition number of the space's positioned dense forms, or
    1.0 when it carries none."""
    return max(
        (float(np.linalg.cond(f.form.matrix)) for f in space.factors if isinstance(f.form, DenseMetric)),
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
    form's condition number, computed by the row from the form itself); the
    bound there is 8·κ·ε, where the factor 8 is chosen headroom over the
    measured worst, 2.35·κ·ε on D1 (2026-10-02). Designed green: π* := Σw·E
    is the same operator and agrees to 0-1 ULP, so this row cannot tell them
    apart (S6.1 can)."""
    V = _population()[name]
    R = AxisRetractionOperator.of(V, "angular")
    sandwich = AdjointOperator(R)  # pyright: ignore[reportArgumentType]
    kappa = _dense_condition(V)
    for phi in _phis(R, seed=13):
        closed, generic = np.asarray(R.H.apply(phi)), np.asarray(sandwich.apply(phi))
        if kappa == 1.0:
            npt.assert_array_almost_equal_nulp(closed, generic, nulp=4)
        else:
            npt.assert_allclose(closed, generic, rtol=0, atol=8 * kappa * np.finfo(float).eps * np.max(np.abs(closed)))


def test_s6_1_a_second_pullback_cannot_be_minted() -> None:
    """The pullback is minted by its retraction, once, so ``R.H.H is R``
    holds for every pullback that exists."""
    R = AxisRetractionOperator.of(_population()["synthetic"], "angular")
    with pytest.raises(TypeError, match="minted by its retraction, once"):
        AxisPullbackOperator(R)


def test_the_axis_verbs_keep_the_positioned_forms() -> None:
    """``without_axis`` and ``with_axis`` are the one route by which a
    product loses or gains an axis; both keep every other axis's positioned
    form (the mint and the moment-tail composer use them). ``of_axes`` over
    the same axes drops the form, which is the defect they replace."""
    D1 = _d1()
    without = D1.without_axis("angular")
    if not isinstance(without.factor("moment").form, DenseMetric):
        pytest.fail("without_axis dropped the kept axis's dense form")
    tail = Factor.whole_on(Axis("tail", (2,), origin=Declared(BasisKind.NODAL)), Kind.MEASURE, VarianceChart.DENSITY)
    widened = without.with_axis(tail)
    forms = [f.form for f in widened.factors]
    if not (len(forms) == 3 and isinstance(forms[1], DenseMetric) and isinstance(forms[2], DiagonalMetric)):
        pytest.fail(f"with_axis lost a form: {forms}")
