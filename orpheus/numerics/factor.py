r"""The factor: one slot of a tensor product, typed by its axis, its variance and the subquotient it spans.

A :class:`Factor` is an :class:`~orpheus.numerics.axis.Axis` (the index set, its form and its origin) together with a
**variance** (how the slot transforms: a function or a measure, a density or an integrated chart) and a
**subquotient** :math:`S/K` (the subspace of the axis's flattened coordinates :math:`\mathbb{R}^n` the slot ranges over, modulo the
kernel :math:`K` it is quotiented by). Every construction states its variance and its subquotient; there is no default.

The two variance types
======================

A factor on a **formed** axis (an axis whose form is a :class:`~orpheus.numerics.metric.DiagonalMetric` or a
:class:`~orpheus.numerics.metric.DenseMetric`) carries a :class:`Variance`, the Klein four-group torsor of ``(kind, chart)``
with its two involutions: :meth:`Variance.dual` flips both coordinates and :meth:`Variance.riesz` flips the chart, the change
of chart a form induces. A factor on a **formless** axis (the form is :class:`~orpheus.numerics.metric.NoForm`) carries a
:class:`FormlessVariance`, two values, ``FUNCTION`` (a density) and ``MEASURE`` (integrated), with ``dual`` only: a chart
change needs a form to pass through, so a formless axis has no ``riesz``, and a function that is integrated on a formless axis is
not a value of either type. :attr:`Factor.kind` and :attr:`Factor.chart` read the slot's :class:`Kind` and
:class:`VarianceChart` off the variance of either type.

The form a subquotient induces
==============================

On the density end the form of :math:`S/K` is the coset minimum of the axis form: the quadratic form minimised over each coset
of :math:`K`, the Schur complement of the Gram on the quotient; on a diagonal form it is the weights on :math:`S \setminus K` and
zero elsewhere.
The form is computed once, at construction, so a constructed factor always has its form. On the integrated end it is the inverse face of its density partner: the density form with its ``face`` flipped
on a dense form (:meth:`~orpheus.numerics.metric.DenseMetric.flipped`, one stated matrix and no second inversion), the dual
metric of the density form on a diagonal form.

Constructibility
================

The factor refuses, at construction, a slot on a singular form that is not quotiented by the form's kernel: a density factor
needs :math:`S \cap \ker G \subseteq K`, an integrated factor needs its density partner (its ``riesz`` image) to satisfy it. It
refuses an indefinite diagonal form when a weight anywhere in :math:`S` of the density partner is negative, on :math:`S \setminus K`
and on :math:`K` alike: the coset minimum :math:`\min_{k \in K} (x + k)^{\mathsf T} G (x + k)` of a diagonal form with a negative
weight at an index of :math:`K` is minus infinity, so no induced form exists. A negative weight outside :math:`S` does not enter the
restriction and does not refuse. Positivity is required here, where the Hilbert structure begins, and the Krein arm of the form is
not built (issue #706).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum, unique

import numpy as np

from orpheus.values.content import ContentIdentity

from .axis import Axis, FromFrame
from .declaration import Ann, Coord, Declaration, Spanned, contains, intersect
from .metric import DenseMetric, DiagonalMetric, Face, HilbertMetric, NoForm, _masked_reciprocal

__all__ = ["Factor", "FormlessVariance", "Kind", "Subquotient", "Variance", "VarianceChart"]


@unique
class Kind(Enum):  # concept: factor | role: component | kind: set
    r"""Whether a slot holds the values of a function or the weights of a measure."""

    FUNCTION = "function"
    MEASURE = "measure"

    def flip(self) -> "Kind":
        r"""The other kind."""
        return Kind.MEASURE if self is Kind.FUNCTION else Kind.FUNCTION


@unique
class VarianceChart(Enum):  # concept: factor | role: component | kind: set
    r"""Whether a slot is read in the density chart or the integrated chart."""

    DENSITY = "density"
    INTEGRATED = "integrated"

    def flip(self) -> "VarianceChart":
        r"""The other chart."""
        return VarianceChart.INTEGRATED if self is VarianceChart.DENSITY else VarianceChart.DENSITY


@unique
class FormlessVariance(Enum):  # concept: factor | role: component | kind: set
    r"""The variance of a factor on a formless axis: one bit.

    A function is a density and a measure is integrated; the chart is implied by the kind, so there is no chart change and no
    ``riesz``.
    """

    FUNCTION = "function"
    MEASURE = "measure"

    def dual(self) -> "FormlessVariance":
        r"""The other value."""
        return FormlessVariance.MEASURE if self is FormlessVariance.FUNCTION else FormlessVariance.FUNCTION

    @property
    def kind(self) -> Kind:
        r"""The slot's :class:`Kind`."""
        return Kind.FUNCTION if self is FormlessVariance.FUNCTION else Kind.MEASURE

    @property
    def chart(self) -> VarianceChart:
        r"""The slot's :class:`VarianceChart`, implied by the kind."""
        return VarianceChart.DENSITY if self is FormlessVariance.FUNCTION else VarianceChart.INTEGRATED


@dataclass(frozen=True)
class Variance:  # concept: factor | role: component | kind: set
    r"""The variance of a factor on a formed axis: the pair ``(kind, chart)``, all four values legal (the Klein four-group torsor)."""

    kind: Kind
    chart: VarianceChart

    def __post_init__(self) -> None:
        if not isinstance(self.kind, Kind):
            raise TypeError(f"Variance kind must be a Kind, got {type(self.kind).__name__}")
        if not isinstance(self.chart, VarianceChart):
            raise TypeError(f"Variance chart must be a VarianceChart, got {type(self.chart).__name__}")

    def dual(self) -> "Variance":
        r"""Both coordinates flipped: the variance of the dual slot."""
        return Variance(self.kind.flip(), self.chart.flip())

    def riesz(self) -> "Variance":
        r"""The chart flipped, the kind kept: the change of chart a form induces."""
        return Variance(self.kind, self.chart.flip())


@dataclass(frozen=True, eq=False)
class Subquotient(ContentIdentity):  # concept: factor | role: component | kind: space
    r"""The subspace ``S`` of :math:`\mathbb{R}^n` a slot ranges over, modulo the kernel ``K`` it is quotiented by, with
    :math:`K \subseteq S` decided.

    ``(in ∪ tang)/tang`` and ``Sub(in)`` are isomorphic and distinct (their duals differ) and compare unequal.
    """

    S: Declaration
    K: Declaration

    def __post_init__(self) -> None:
        for name, d in (("S", self.S), ("K", self.K)):
            if not isinstance(d, (Coord, Spanned, Ann)):
                raise TypeError(f"Subquotient {name} must be a Declaration (Coord, Spanned or Ann), got {type(d).__name__}")
        if self.S.n != self.K.n:
            raise ValueError(f"subquotient: S is declared on R^{self.S.n} and K on R^{self.K.n}")
        if not contains(self.S, self.K):  # an undecided inclusion raises "subspace inclusion is undecided"
            raise ValueError(f"subquotient kernel is not declared inside its subspace: K = {self.K!r}, S = {self.S!r}")

    @classmethod
    def whole(cls, n: int) -> "Subquotient":
        r"""The whole space modulo nothing."""
        return cls(Coord(n, tuple(range(n))), Coord(n, ()))

    @property
    def n(self) -> int:
        r"""The dimension of the ambient space."""
        return self.S.n

    @property
    def live(self) -> tuple[int, ...] | None:
        r"""The live set ``S`` minus ``K`` as sorted coordinate indices when both are ``Coord``, else ``None`` (a declared
        subspace has no enumerable coordinates).

        The meet ``S`` with the annihilator of ``K``, spelled in coordinates: the one place the set difference is written, until a
        non-``Coord`` arm of the declaration algebra enumerates its own.
        """
        if isinstance(self.S, Coord) and isinstance(self.K, Coord):
            return tuple(sorted(set(self.S.indices) - set(self.K.indices)))
        return None

    @property
    def is_whole(self) -> bool:
        r"""True for the whole space modulo nothing."""
        return self == Subquotient.whole(self.n)

    @property
    def dim(self) -> int:
        r"""``dim S - dim K``."""
        return self.S.dim - self.K.dim

    def dual(self) -> "Subquotient":
        r"""The subquotient of the dual: ``(Ann K) / (Ann S)``."""
        return Subquotient(self.K.ann(), self.S.ann())


def _flat_weights(form: DiagonalMetric) -> np.ndarray:
    r"""The weights of a diagonal form on the flattened index set (row-major)."""
    return np.asarray(form.weights).ravel()


@dataclass(frozen=True, eq=False)
class Factor(ContentIdentity):  # concept: factor | role: head | kind: space
    r"""One slot of a tensor product: an axis, the variance the slot is read in and the subquotient it spans.

    Refused at construction: a subquotient on another number of coordinates than the axis (``ValueError``); a variance of the
    wrong type for the axis's form (``TypeError``); a slot on a singular form that is not quotiented by its kernel; an indefinite
    diagonal form when a weight anywhere in ``S`` of the density partner is negative (the Krein arm, issue #706); the
    zero-dimensional state on a dense form; a state whose induced form the table of :meth:`_density_form` does not build (a ``Spanned``
    kernel not the form's own, an annihilator that is not the form's). The induced form is built once, here, so no read of a constructed
    factor refuses.
    """

    axis: Axis
    variance: FormlessVariance | Variance
    sk: Subquotient
    #: The density form of the density partner, computed once at construction (``None`` on a formless axis): a constructed factor has
    #: its form, so construction refuses exactly the states the build cannot make.
    _induced: DiagonalMetric | DenseMetric | None = field(init=False, compare=False, repr=False)

    def __post_init__(self) -> None:
        if not isinstance(self.axis, Axis):
            raise TypeError(f"Factor axis must be an Axis, got {type(self.axis).__name__}")
        if not isinstance(self.variance, (FormlessVariance, Variance)):
            raise TypeError(f"Factor variance must be a FormlessVariance or a Variance, got {type(self.variance).__name__}")
        if not isinstance(self.sk, Subquotient):
            raise TypeError(f"Factor sk must be a Subquotient, got {type(self.sk).__name__}")
        n = int(np.prod(self.axis.shape)) if self.axis.shape else 1
        if self.sk.n != n:
            raise ValueError(f"factor: the subquotient lives in R^{self.sk.n}, the axis {self.axis.label!r} has {n} coordinates")
        formless = isinstance(self.axis.form, NoForm)
        # ELEGANCE-DEBT[guard] #757 the variance type must agree with the axis's form at run time; retires when the factor splits
        # by the axis's form into a formless and a formed class, so that the disagreement is unspellable.
        if isinstance(self.variance, FormlessVariance) != formless:
            raise TypeError(
                f"factor: axis {self.axis.label!r} is {'formless' if formless else 'formed'}, so its factor takes a "
                f"{'FormlessVariance' if formless else 'Variance'}, got {type(self.variance).__name__}"
            )
        if formless:
            object.__setattr__(self, "_induced", None)
            return
        self._require_constructible()
        if isinstance(self.axis.form, DenseMetric) and self.sk.dim == 0:
            # SCOPE-BOUNDARY[guard]
            # machinery: a zero-form arm of Form, which a DenseMetric of dimension zero would need
            # ruling: ON03-2 (scratch/ontology/w1/open/N03/rulings.md), the orchestrator, 2026-11-24
            # revisit: W8a
            raise ValueError("the induced dense form is the zero form, which DenseMetric does not spell")
        axis_form = self.axis.form
        assert not isinstance(axis_form, NoForm)  # narrowing only: the formless case returned above
        object.__setattr__(self, "_induced", self._density_form(self._density_partner(), axis_form))

    def _density_partner(self) -> Subquotient:
        r"""The subquotient of the density chart this slot's form is induced from: this one on the density end, its dual on the
        integrated end (the ``riesz`` image)."""
        return self.sk if self.variance.chart is VarianceChart.DENSITY else self.sk.dual()

    def _require_constructible(self) -> None:
        sk, form = self._density_partner(), self.axis.form
        assert not isinstance(form, NoForm)  # narrowing only: the formless case returned before
        if isinstance(form, DiagonalMetric):
            weights = _flat_weights(form)
            # the weights the restriction to S reads: S itself, K included (a coset minimum over K of a negative weight is minus infinity)
            in_s = sk.S.indices if isinstance(sk.S, Coord) else range(weights.size)
            if bool((weights[list(in_s)] < 0.0).any()):
                raise ValueError(
                    "factor: the form has a negative weight on the subspace the slot ranges over, so it has no Hilbert chart; the Krein "
                    "arm of the form is not built (issue #706)"
                )
        kernel = form.kernel
        meet = intersect(sk.S, kernel)  # an undecided meet raises "subspace intersection is undecided"
        if not contains(sk.K, meet):
            raise ValueError(f"factor: a slot on a singular form needs the quotient by its kernel (quotient by {kernel!r})")

    @property
    def form(self) -> HilbertMetric | NoForm:
        r"""The form of the slot: :class:`NoForm` on a formless axis; else the density form of the density partner, and on the
        integrated end its inverse face (the dense form with its ``face`` flipped, the dual metric of a diagonal form)."""
        density = self._induced
        if density is None:
            return self.axis.form  # a NoForm
        if self.variance.chart is VarianceChart.DENSITY:
            return density
        if isinstance(density, DenseMetric):
            return density.flipped()
        if isinstance(self.axis.origin, FromFrame):
            # ELEGANCE-DEBT[guard] #756 the diagonal frame head's integrated form is the masked reciprocal of its Gram's diagonal,
            # the arithmetic the head has always used; retires when N10 gives the diagonal integrated end one spelling.
            return DiagonalMetric(_masked_reciprocal(_flat_weights(density).reshape(self.axis.shape)))
        return density.dual()

    def _density_form(self, sk: Subquotient, axis_form: DiagonalMetric | DenseMetric) -> DiagonalMetric | DenseMetric:
        r"""The coset minimum of ``axis_form`` on ``sk``, stored on the axis's full index set."""
        if sk.is_whole:
            return axis_form
        if isinstance(axis_form, DenseMetric) and isinstance(sk.S, Coord) and sk.S.dim == sk.S.n and isinstance(sk.K, Spanned):
            if sk.K != axis_form.kernel:  # the whole space modulo a declared subspace is the form's own only modulo ITS kernel
                raise ValueError("subquotient kernel is not the form's declared kernel")
            return axis_form
        if isinstance(axis_form, DiagonalMetric) and isinstance(sk.S, Coord) and isinstance(sk.K, Coord):
            keep = np.zeros(sk.n, dtype=bool)
            keep[list(sk.live or ())] = True
            return DiagonalMetric(np.where(keep, _flat_weights(axis_form), 0.0).reshape(self.axis.shape))
        if isinstance(axis_form, DenseMetric) and isinstance(sk.S, Coord) and isinstance(sk.K, Coord):
            live = list(sk.live or ())
            # the form is exactly zero on its own Coord kernel and couples to nothing there: eliminate only the rest of K
            own_kernel = set(axis_form.kernel.indices) if isinstance(axis_form.kernel, Coord) else set()
            quotiented = [i for i in sk.K.indices if i not in own_kernel]
            gram = np.asarray(axis_form.matrix)
            out = np.zeros_like(gram)
            if live:
                block = gram[np.ix_(live, live)]
                if quotiented:
                    coupling = gram[np.ix_(live, quotiented)]
                    block = block - coupling @ np.linalg.solve(gram[np.ix_(quotiented, quotiented)], coupling.T)
                out[np.ix_(live, live)] = block
            dead = tuple(sorted(set(sk.K.indices) | (set(range(sk.n)) - set(sk.S.indices))))
            return DenseMetric(out, kernel=Coord(sk.n, dead), face=Face.METRIC)
        if isinstance(axis_form, DenseMetric) and isinstance(sk.S, Ann):
            kernel = axis_form.kernel
            if not (isinstance(kernel, Spanned) and sk.S == kernel.ann()):
                raise ValueError("subquotient subspace is not the annihilator of the form's declared kernel")
            if sk.K != Coord(sk.n, ()):
                raise ValueError("subquotient kernel is not the form's declared kernel")
            return axis_form
        raise ValueError("a declared subspace has no enumerable coordinates: the induced form of this (form, S, K) is not built")

    @property
    def kind(self) -> Kind:
        r"""The slot's :class:`Kind`, read off the variance of either type."""
        return self.variance.kind

    @property
    def chart(self) -> VarianceChart:
        r"""The slot's :class:`VarianceChart`, read off the variance of either type."""
        return self.variance.chart

    @property
    def stored_shape(self) -> tuple[int, ...]:
        r"""The axis's index shape: what the slot stores."""
        return self.axis.shape

    @property
    def dim(self) -> int:
        r"""``dim S - dim K``."""
        return self.sk.dim

    def dual(self) -> "Factor":
        r"""The dual slot: both variance coordinates flipped, the subquotient dualised."""
        return Factor(self.axis, self.variance.dual(), self.sk.dual())

    def riesz(self) -> "Factor":
        r"""The slot's image under the form's change of chart.

        A formless factor has no ``riesz``: its :class:`FormlessVariance` has none, so the call raises ``AttributeError``.
        """
        # a FormlessVariance has no riesz: the AttributeError of the formless case is raised by the variance (#757)
        return Factor(self.axis, self.variance.riesz(), self.sk.dual())  # pyright: ignore[reportAttributeAccessIssue]
