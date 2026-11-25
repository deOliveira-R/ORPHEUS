r"""Function spaces for the operator-algebra framework.

A *function space* is the domain or codomain of a linear operator. In
the matrix-free transport algebra, every operator :math:`A : V \to W`
acts on discrete flux distributions that live in a definite
:class:`FunctionSpace` — angular-flux space, scalar-flux space, a
boundary-trace space, and so on. Tagging operators with their
:attr:`domain` and :attr:`range` lets the composition machinery in
:mod:`orpheus.numerics.operator` reject ill-formed compositions at
construction time (raising :class:`IncompatibleOperatorComposition`),
preventing the harmful-stub anti-pattern where a downstream Krylov
consumer hits a shape mismatch mid-iteration.

The *Hilbert-adjoint* refinement (:meth:`LinearOperator.adjoint`,
``A.H``) further requires the domain and range to carry an inner
product. :class:`FunctionSpace` stores the inner-product weights as
metadata: the L² inner product
:math:`\langle x, y \rangle_w = \sum_i w_i \, x_i \, y_i` reduces to
the Euclidean :math:`\sum_i x_i \, y_i` when the weights are absent,
which is the natural default for spaces that have no canonical
quadrature (cell-flux, region-flux). For angular-flux and
boundary-trace spaces the quadrature weights ARE the canonical
inner-product weights, so the adjoint identity
:math:`\langle A x, y \rangle_W = \langle x, A^* y \rangle_V`
becomes a non-trivial consistency check (see test
``tests/gates/numerics/test_operator.py::test_hilbert_adjoint_weighted_identity``).

Where the anticipated specialisations landed
============================================

The Grand Report v3 (§5.3, §6.1) anticipated a richer Space ontology as
subclasses layered on :class:`FunctionSpace`. It landed as AXES and as
whole-boundary spaces instead, under the ruling that a space is the
ordered product of its axes (:class:`~orpheus.numerics.axis.Axis`):

* ``MeshFunctionSpace`` — the spatial axis, minted by the hub from the
  cell volumes (``MaterialMesh.bulk_space``).
* **AngularTraceSpace** — shipped as ONE whole-boundary space
  (:mod:`orpheus.numerics.spaces.angular_trace_space`); inflow / outflow
  are selectors over its signed :math:`\Omega\cdot\hat n`, not
  directional tags.
* ``RegionSpace`` — a coarser stage of the same filtration, reached by
  the axis marginal (:mod:`orpheus.numerics.axis_marginal`: retraction
  and section), not by a class.
* ``EnergyGroupSpace`` — :class:`~orpheus.numerics.axis.EnergyAxis`, a
  one-dimensional mesh in energy.
* ``DiscreteAngularSpace`` —
  :meth:`~orpheus.numerics.axis.Axis.of_measure` on the
  quadrature's measure, which carries the invariance group.

Of the compositional dunders anticipated with them, ``S * T`` (the tensor
product, :class:`TensorSpace`, the space of factors) and ``S.dual()`` (the dual space)
ship on :class:`TensorSpace`; the direct sum is not
a dunder but the typed composites
:class:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace` and
:class:`~orpheus.numerics.coupled_system.CoupledSpace`.

References
----------

* Trefethen, L.N. & Bau, D. (1997). *Numerical Linear Algebra*. SIAM.
  §1 — vector spaces, inner products, the Hilbert adjoint vs.
  representation transpose distinction.
* Reed, M. & Simon, B. (1980). *Methods of Modern Mathematical
  Physics I: Functional Analysis*, §III.6 (Hilbert spaces, Riesz
  representation, the adjoint operator).
"""

from __future__ import annotations

from collections.abc import Mapping
from functools import cached_property
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import TYPE_CHECKING, Any, Callable, Generic, Optional, TypeVar, cast

import numpy as np
from numpy.typing import NDArray

from .axis import Axis, BasisKind, Declared
from .coordinate_space import coordinates, leading_coordinates
from .factor import Factor, FormlessVariance, Kind, Subquotient, Variance, VarianceChart
from orpheus.values.content import content_digest, name_digest
from .metric import (
    DenseMetric,
    DiagonalMetric,
    DualMetric,
    FactoredMetric,
    HilbertMetric,
    NoForm,
    _require_positive,
)

if TYPE_CHECKING:
    from .basis.base import Basis, TruncatedBasis
    from .measure import DiscreteMeasure


__all__ = [
    "CarrierSpec",
    "Covector",
    "DualSpace",
    "FunctionSpace",
    "Member",
    "TensorSpace",
    "act_on_direct_sum_coordinates",
    "act_on_flat_coordinates",
    "admit_matrix",
    "coordinate_space_of",
    "factor_tail",
]

#: The space's element type (what a member of :math:`V` IS): a bare
#: ``NDArray`` for the leaf spaces, a bulk ⊕ boundary composite field for
#: :class:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace`. The
#: PEP-696 default is ``Any`` because a bare ``FunctionSpace`` annotation
#: (the operator layer's ``domain`` / ``codomain`` slots) genuinely holds
#: EITHER realization today — specializing those slots to the operator's
#: own carrier (``FunctionSpace[Domain]``, the #65 two-param discipline
#: extended to the space layer) is the follow-on step, not this one.
Carrier = TypeVar("Carrier", default=Any)


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class Member:  # concept: member | role: head | kind: element.member
    r"""A typed element of a space: its coordinates and the space they belong to.

    The nominal type the space's membership check reads. A member is
    compared with the space by content identity on every run
    (:meth:`FunctionSpace.coordinates_of`), which is the only check that
    separates two spaces of equal shape and different metric; a bare array
    is the raw path. Its subclasses are
    :class:`~orpheus.numerics.field.Field`, every production quantity on
    :math:`V` with its physical role and units, and :class:`Covector`, an
    element of :math:`V^{*}`. A bare ``Member`` is an element of :math:`V`
    with no physical role: the constant function :math:`1` that
    :attr:`FunctionSpace.integral` lowers is one. This base
    lives here, with the space that asks the question, because ``field.py``
    imports this module and the concept order runs space before field.

    Parameters
    ----------
    values : NDArray
        The coordinates; ``values.shape`` must equal ``space.shape``.
    space : FunctionSpace
        The space the coordinates belong to.

    Notes
    -----
    ``eq=False``: an ndarray field makes value equality element-wise, so a
    member is identified by object identity.
    """

    values: NDArray
    space: "FunctionSpace"

    def __post_init__(self) -> None:
        coordinates(self.values, self.space.shape, owner=type(self).__name__)


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class Covector(Member):  # concept: covector | role: head | kind: element.covector
    r"""A co-vector: a member of the dual :math:`V^{*}`, paired canonically with :math:`V`.

    The functional :math:`\langle w, \cdot\rangle_V` of a density :math:`w`
    on :math:`V` is the element :math:`\flat_V w = G\,w` of :math:`V^{*}`,
    and its value on a member :math:`x` of :math:`V` is the canonical,
    metric-free pairing

    .. math::

       c(x) = \sum_i c_i\, x_i , \qquad c = G\,w ,
       \qquad\text{so}\qquad c(x) = \langle w, x\rangle_V .

    The metric is read once, at the mint (:meth:`FunctionSpace.flat`), and
    never again: :meth:`evaluate` contracts coordinates. Its ``space`` is
    :math:`V^{*}` (``V.dual()``), so ``evaluate`` asks the dual's primal
    whether ``x`` is a member, and a member of an equal-shape space with
    another metric is refused.

    The specialised constructor is :meth:`FunctionSpace.flat`; the raw one,
    ``Covector(values=…, space=V.dual())``, is for a dependency-free test
    and must be handed the coordinates the mint would compute.

    Conforms to :class:`~orpheus.numerics.functional.Functional`
    structurally (it exposes ``evaluate``); it is declared here, beside
    :class:`DualSpace`, because :meth:`FunctionSpace.flat` mints it and the
    concept order runs space before functional.
    """

    def evaluate(self, x: Any, /) -> float:
        r"""Return :math:`c(x) = \sum_i c_i x_i` for a member ``x`` of :math:`V`."""
        return float(np.sum(self.values * self.space.dual().coordinates_of(x)))

    def __repr__(self) -> str:
        return f"Covector(shape={self.values.shape}, space={self.space.name!r})"


#: The subclasses of :class:`FunctionSpace`: module-qualified class name to the unit that removes the arm. :class:`TensorSpace` is the arm that
#: stays; the nine legacy arms leave with the unit named (N05 the direct sums, N06 the trace spaces, N07 the dual). The class-creation hook
#: of :class:`FunctionSpace` refuses every other subclass at any depth. A module-level constant, so that the names it holds are data and
#: not reference edges of the concept table.
FUNCTION_SPACE_ARMS: Mapping[str, str | None] = MappingProxyType({
    "orpheus.numerics.space.TensorSpace": None,
    "orpheus.numerics.coupled_system.CoupledSpace": "N05",
    "orpheus.numerics.spaces.full_field_space.FullFieldSpace": "N05",
    "orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace": "N06",
    "orpheus.numerics.spaces.angular_trace_space.AngularFaceTraceSpace": "N06",
    "orpheus.numerics.spaces.scalar_trace_space.ScalarTraceSpace": "N06",
    "orpheus.numerics.spaces.radial_characteristic_space._RadialCharacteristicSubSpace": "N06",
    "orpheus.numerics.spaces.radial_characteristic_space.RadialCharacteristicInteriorSpace": "N06",
    "orpheus.numerics.spaces.radial_characteristic_space.RadialCharacteristicBoundarySpace": "N06",
    "orpheus.numerics.space.DualSpace": "N07",
})


@dataclass(frozen=True)
class CarrierSpec:  # concept: tensor | role: component | kind: set
    r"""The element classes of a space and of its dual: the pair ``(member, dual_member)``, :class:`Member` and :class:`Covector` by default.

    A space mints its members as ``carrier.member``; the dual of a space carries the swapped pair (:meth:`dual`), so the members of
    ``V.dual()`` are covectors, the data side, and those of ``V.dual().dual()`` are members again.

    SCOPE-BOUNDARY[guard]
    machinery: members that follow their space by type, so that the member class is a function of the space's content (W3)
    ruling: issue #707 (second comment, 2026-11-11, confirmed 2026-11-12) and section 6.5 of .claude/plans/ontology_w1_spec.md
    revisit: W3, the known defect: two equal spaces carrying different carriers are equal and mint different member classes; the space's own verbs reach it, since partial duals composed to the all-slot dual carry the default pair (witness: the strict xfail in tests/gates/numerics/w1/test_n04_review_repairs.py)
    """

    member: type[Member] = Member
    dual_member: type[Member] = Covector

    def __post_init__(self) -> None:
        for name, cls in (("member", self.member), ("dual_member", self.dual_member)):
            if not (isinstance(cls, type) and issubclass(cls, Member)):
                raise TypeError(f"CarrierSpec {name} must be a Member class, got {cls!r}")

    def dual(self) -> "CarrierSpec":
        r"""The pair of the dual space: the two classes swapped."""
        return CarrierSpec(self.dual_member, self.member)


@dataclass(frozen=True, eq=False)
class FunctionSpace(Generic[Carrier]):  # concept: tensor | role: weld | kind: space | welded: tensor, fmetric, factor
    r"""The common type of the spaces of discrete fields: a :class:`TensorSpace` over factors, and the legacy leaves and sums that leave with N05 to N07.

    A leaf is stored by its ``name``, its ``shape`` and its metric, which is either diagonal weights (``inner_product_weights``, broadcastable
    to the shape, leading-aligned) or a metric object (``metric``), or neither, which is the Euclidean inner product. The axis-built arm is
    :class:`TensorSpace`: a space is the ordered product of the factors of its slots, each with its form, its variance and its subquotient, and
    the metric is the Kronecker product of the forms in slot order. The subclasses are an allow-list
    (:data:`FUNCTION_SPACE_ARMS`): a class that is not on it is refused when it is defined.

    A space encodes pure *geometry*: the degrees of freedom (``shape``), the inner product and the composition algebra (``*``, :meth:`dual`).
    Units do NOT live on the space; they are a property of the quantity, carried by the :class:`~orpheus.numerics.field.Field` role leaf
    (issues #205, #207, #208).

    **Identity** is :meth:`identity_key`, the one definition from which ``==`` and ``hash`` derive: a leaf is its class, name, shape and the
    content digests of its metric; a tensor space is its factors; a dual, a direct sum and a coupled system are their members' keys, tagged,
    so covariance travels into a composite. Two spaces with the same ``name`` and ``shape`` need not be equal.

    Parameters
    ----------
    name : str
        Human-readable identifier, used by :meth:`__repr__` and by the composition guard's messages.
    shape : tuple[int, ...]
        Tensor shape of the elements.
    inner_product_weights : NDArray | None, default None
        Diagonal weights of the L² inner product :math:`\langle x, y \rangle = \sum_i w_i x_i y_i`, broadcast-compatible with ``shape``.
    metric : HilbertMetric | None, keyword-only, default None
        A metric object for a form no diagonal array can spell; validated by :meth:`HilbertMetric.validate_for`. One metric source only.
    """

    name: str
    shape: tuple[int, ...]
    # Metadata, NOT identity: ``compare=False`` keeps an array out of the content of any dataclass that would read the fields (a compared
    # array makes the generated ``__eq__`` of a subclass return an ndarray and ``hash`` raise).
    inner_product_weights: Optional[NDArray] = field(
        default=None, repr=False, compare=False,
    )
    metric: Optional[HilbertMetric] = field(
        default=None, repr=False, compare=False, kw_only=True,
    )

    def __init_subclass__(cls, **kwargs: Any) -> None:
        super().__init_subclass__(**kwargs)
        key = f"{cls.__module__}.{cls.__qualname__}"
        if key not in FUNCTION_SPACE_ARMS:
            raise TypeError(
                f"{key}: arm is not on the FunctionSpace allow-list; a space is a TensorSpace over factors "
                f"(the arms are {sorted(FUNCTION_SPACE_ARMS)})"
            )

    def __post_init__(self) -> None:
        # ELEGANCE-DEBT[guard] #707 legacy-leaf: the leaf arm validates its stored metric; retires at N07, when the leaf fields leave the base.
        if self.metric is not None:
            # The object knows its own admission (a DenseMetric must span the flattened leading block).
            self.metric.validate_for(self.shape)
        if self.inner_product_weights is not None and self.metric is not None:
            raise ValueError(
                f"space {self.name!r} carries BOTH dense "
                f"inner_product_weights and a metric object — one metric "
                f"source only (install the object OR the array, never both)"
            )
        if self.inner_product_weights is not None:
            DiagonalMetric(self.inner_product_weights).validate_for(self.shape)

    def identity_key(self) -> tuple[Any, ...]:
        r"""The content that IS this space: the one definition of its identity.

        ``==`` and ``hash`` derive from it here, in the base class only (#631, finding F7 of the S6a review). A tagged tuple per kind of space:

        * a leaf: ``("named", class, name, shape, weights, metric)``, the class keeping a hand-built space wearing a factory leaf's name
          from passing for one, and the content digests of the weights and of the metric keeping two spaces of one name and shape with
          different metrics apart;
        * a tensor space, a dual, a direct sum and a coupled system override this with their members' keys, tagged.

        ELEGANCE-DEBT[guard] #707 legacy-leaf: the leaf key is the identity of the five legacy leaves; retires at N07, when they are tensor spaces.
        """
        return (
            "named", type(self), self.name, self.shape,
            None if self.inner_product_weights is None else content_digest(np.asarray(self.inner_product_weights)),
            None if self.metric is None else content_digest(self.metric),
        )

    @cached_property
    def _identity(self) -> tuple[Any, ...]:
        return self.identity_key()

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, FunctionSpace):
            return NotImplemented
        return other is self or self._identity == other._identity

    def __hash__(self) -> int:
        return hash(self._identity)

    def __repr__(self) -> str:
        return f"FunctionSpace({self.name!r}, shape={self.shape})"

    # ------------------------------------------------------------------
    # The leaf arm of the axis verbs: a leaf has no named factors
    # ------------------------------------------------------------------

    @property
    def has_coordinate_cone(self) -> bool | None:
        r"""Whether per-component positivity is meaningful on this space.

        Three-valued: ``True`` on a tensor space whose factors are all ``NODAL`` (components are point or cell values, the coordinate cone
        :math:`K = \{x \ge 0\}` is the physical positive cone); ``False`` where any factor is ``MODAL`` (components are expansion coefficients,
        so a per-component sign test is meaningless and is refused, not answered); ``None`` on a leaf, whose structure cannot answer.

        ELEGANCE-DEBT[guard] #707 legacy-leaf: a leaf answers ``None``; retires at N07, when no space is a leaf.
        """
        return None

    def _axis_index(self, label: str) -> int:
        r"""Position of the unique factor labeled ``label``: the one by-label resolver, which a leaf refuses (it has no named factors).

        ELEGANCE-DEBT[guard] #707 legacy-leaf: the refusal of a space that is not a tensor space; retires at N07.
        """
        raise TypeError(
            f"axis lookup {label!r}: {self!r} is a leaf space with no named factors. Compose the space with TensorSpace.of_axes."
        )

    def axis(self, label: str) -> Axis:
        r"""The axis of the factor labeled ``label``: :class:`TensorSpace` answers; a leaf has no named factors and refuses it."""
        raise TypeError(
            f"axis lookup {label!r}: {self!r} is a leaf space with no named factors. Compose the space with TensorSpace.of_axes."
        )

    def _axis_block(self, k: int) -> DiagonalMetric | DenseMetric | DualMetric | None:
        r"""The metric block of factor ``k``, which a leaf refuses (it has no named factors).

        ELEGANCE-DEBT[guard] #707 legacy-leaf: the refusal of a space that is not a tensor space; retires at N07.
        """
        raise TypeError(f"{self!r} is a leaf space with no named factors, so it has no metric block {k}.")

    def without_axis(self, label: str) -> "FunctionSpace":
        r"""The product of the OTHER factors; a leaf refuses it (it has no named factors).

        ELEGANCE-DEBT[guard] #707 legacy-leaf: the refusal of a space that is not a tensor space; retires at N07.
        """
        raise TypeError(f"without_axis {label!r}: {self!r} is a leaf space with no named factors.")

    def with_axis(self, factor: Factor) -> "FunctionSpace":
        r"""This product with ``factor`` appended last; a leaf refuses it (it has no named factors).

        ELEGANCE-DEBT[guard] #707 legacy-leaf: the refusal of a space that is not a tensor space; retires at N07.
        """
        raise TypeError(f"with_axis: {self!r} is a leaf space with no named factors.")

    def positions_spanned_by(self, sub: "FunctionSpace") -> tuple[int, ...]:
        r"""The array positions of this space that the factors of ``sub`` occupy, read by factor label.

        A leaf cannot say which of its positions a factor occupies, and is refused rather than guessed (the user's ruling, 2026-10-20).

        ELEGANCE-DEBT[guard] #707 legacy-leaf: the refusal of a space that is not a tensor space; retires at N07.
        """
        raise TypeError(
            f"positions_spanned_by: the space {self!r} is not built from "
            f"labelled factors, so the positions {sub!r} spans cannot be read off it. Compose the space with TensorSpace.of_axes."
        )

    # ------------------------------------------------------------------
    # Inner product / norm
    # ------------------------------------------------------------------

    def coordinates_of(self, x: Any) -> Any:
        r"""Return the coordinates of ``x`` as a member of this space, or refuse it.

        Membership is the space's verb, asked by every verb that takes an
        element, and it is a parse with two legs, both checked on every run
        (``-O`` included):

        * **the identity leg.** A typed member, a :class:`Member` (a
          :class:`~orpheus.numerics.field.Field`, a :class:`Covector` on a
          dual space, or the bare constant function), is compared with this
          space by content identity, and its coordinates are returned. This
          is the only check that separates two spaces of equal shape and
          different metric.
        * **the shape leg.** A bare array must be an element of the
          coordinate space :math:`\mathbb R^{\text{shape}}`: its shape equals
          :attr:`shape` exactly
          (:func:`~orpheus.numerics.coordinate_space.coordinates`). A scalar
          verb's elements are members, and a trailing batch has no inner
          product (a pairing over a payload is posed on the explicit tensor
          product ``V * P``).

        The shape leg is what keeps numpy's broadcast closure out of every
        verb: handed a ``(2,)`` array where the space is ``(2, 1)``, the
        product would silently be an outer product (ERR-122). A bare array
        passes only the shape leg; an equal-shape array of another space is
        indistinguishable from a member, which is why production hands typed
        members (#664 tracks the operator algebra's bare iterate).
        """
        if isinstance(x, Member):
            if x.space != self:
                raise ValueError(
                    f"{type(x).__name__} on {x.space!r} is not a member of {self!r}."
                )
            return x.values
        return coordinates(x, self.shape, owner=repr(self))

    def leading_coordinates_of(self, x: Any) -> Any:
        r"""Return the coordinates of ``x`` for an ACTION on this space, or refuse it.

        An action (a metric leg, an operator's matvec) acts on the space's
        leading axes, so a trailing payload rides through unchanged: the
        action is :math:`A \otimes I_P` on :math:`V \otimes P`. A typed member
        is compared by content identity (see :meth:`coordinates_of`); a bare
        array's leading shape must equal :attr:`shape`
        (:func:`~orpheus.numerics.coordinate_space.leading_coordinates`). Both
        legs run on every run.
        """
        if isinstance(x, Member):
            return self.coordinates_of(x)
        return leading_coordinates(x, self.shape, owner=repr(self))

    def inner_product(self, x: Carrier, y: Carrier) -> float:
        r"""Return :math:`\langle x, y \rangle`.

        The **surface** is carrier-generic: a space defines the inner
        product on its own element type. The base **realization**
        (:meth:`_diagonal_inner_product`) is the diagonal-weight form for
        the bare-array carrier; a space with a structured carrier
        overrides this surface with its own realization
        (:class:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace`
        dispatches per direct-sum block).

        Parameters
        ----------
        x, y : Carrier
            Members of this space: typed members (compared by content
            identity), or bare arrays of exactly :attr:`shape` on the raw
            path (see :meth:`coordinates_of`).

        Returns
        -------
        float
            Scalar inner product.
        """
        return self._diagonal_inner_product(self.coordinates_of(x), self.coordinates_of(y))

    def _diagonal_inner_product(self, x: Any, y: Any) -> float:
        r"""The bare-array realization of :meth:`inner_product`.

        With diagonal weights ``w`` the inner product is the weighted
        sum :math:`\sum_i w_i \, x_i \, y_i`. Without weights it
        reduces to the Euclidean :math:`\sum_i x_i \, y_i`. The
        weights array is broadcast against ``x * y`` through the resolved
        :class:`~orpheus.numerics.metric.DiagonalMetric`
        (:func:`orpheus.numerics.metric._broadcast_leading` — the SAME
        leading-axis convention :meth:`_diagonal_apply_metric` uses, one
        home) — so a 1-D weight vector
        along (say) the ordinate axis acts on the full
        ``(n_ordinates, n_groups, *spatial)`` tensor without manual
        reshaping. Valid only for an ``NDArray``-carried space — the
        ``Any`` parameters are the realization/surface seam, not an
        open contract.

        Notes
        -----
        ⛔ **This routed through numpy's default (TRAILING) broadcast
        until 2026-08-04, while :meth:`_diagonal_apply_metric` used the
        LEADING convention — the same metric applied along different
        axes by two methods of one space.** The divergence is invisible
        whenever ``w.ndim >= x.ndim`` (the leading broadcast is then a
        no-op), which is every case the tree exercised, so it shipped.
        It bites the moment a space carries a leading-axis metric over
        an element with trailing axes:

        * non-square shapes — :meth:`inner_product` raised
          ``ValueError`` while :meth:`apply_metric` worked. `[M]`
          the degree-3 spherical-harmonic head on its PRODUCTION layout
          ``(L+1, 2L+1, ng, *spatial)`` did exactly this.
        * square shapes — both succeeded and **silently disagreed**
          (`[M]` 456 vs 552 on a ``(3, 3)`` probe), which is the
          dangerous half.

        Either way :math:`\langle Ax, y\rangle = \langle x, A^\dagger
        y\rangle` is false by construction, since ``AdjointOperator``
        builds :math:`A^\dagger = G_V^{-1}A^{\mathsf T}G_W` from
        :meth:`apply_metric` while the pairing that judges it came from
        here — the ERR-067 family, one layer down. The fix is
        bit-identical wherever the old path did not raise: it only ever
        pads ``w`` with trailing 1s, which is a no-op when ``w`` already
        spans every axis.
        """
        m = self._resolved_metric
        if m is None:
            return float(np.sum(x * y))
        # Σ (G⊙x)·y through the realization — for the resolved diagonal
        # case this is np.sum((w_b*x)*y), the SAME reduction tree as the
        # legacy np.sum(w*x*y) (left-to-right), so the reroute is
        # bit-identical; the matmul spelling y@(diag(w)@x) is NOT ([M] up
        # to 1792 ULP at n=15) and is deliberately unspellable here.
        return m.pairing(x, y)

    def norm(self, x: Carrier) -> float:
        r"""Return the induced :math:`L^2` norm
        :math:`\sqrt{\langle x, x \rangle}`.

        Carrier-generic through :meth:`inner_product` — valid unchanged
        for a structured-carrier space that overrides the inner product.

        Raises ``ValueError`` on a space that reads a negative direction (:attr:`negative`, summed over the members of a composite), in every
        build: the square root of a signed form is not a norm, and the Krein arm of the form is not built (issue #706). The reading is exact on a
        diagonal form (and on its dual and its tensor products) and declared as 0 on a dense one, whose positivity rests on the formal path; a
        raw indefinite dense form under ``python -O`` is caught by the one exact check left, a negative :math:`\langle x, x\rangle`, with the
        same message.
        """
        _require_positive(self, "FunctionSpace.norm")
        square = self.inner_product(x, x)
        if square < 0.0:
            raise ValueError(
                f"FunctionSpace.norm needs a positive form: <x, x> = {square} is negative, and the Krein arm of the form is not built (issue #706)"
            )
        return float(np.sqrt(square))

    @property
    def negative(self) -> int:
        r"""The number of negative directions of the space's form: its resolved metric's, 0 on the Euclidean space; a composite sums its members'."""
        m = self._resolved_metric
        return 0 if m is None else m.negative

    @property
    def positive(self) -> bool:
        r"""Whether the space's form has no negative direction: the predicate :func:`~orpheus.numerics.metric._require_positive` asks of a space
        and of a form alike."""
        return self.negative == 0

    def flat(self, density: Any) -> Covector:
        r"""Mint the co-vector :math:`\flat_V w = G\,w`, the functional :math:`\langle w, \cdot\rangle_V`.

        The one specialised constructor of :class:`Covector`: ``density``
        is a member of this space (a typed member compared by content, or a
        bare array of the exact shape), and the result
        is a member of :math:`V^{*}` whose value on a member :math:`x` is
        :math:`\langle w, x\rangle_V`. The metric is applied here, once,
        through :meth:`apply_metric`, the same arithmetic as the lowering
        leg :class:`~orpheus.numerics.operator.RieszLowerOperator`; the
        number one is :attr:`unit`, never ``np.ones``.
        """
        return Covector(
            values=np.asarray(self.apply_metric(self.coordinates_of(density))),
            space=self.dual(),
        )

    @property
    def unit(self) -> NDArray:
        r"""The coordinates of the constant function :math:`1` on this space: ones on a leaf, whose coordinates are nodal by construction.

        :class:`TensorSpace` overrides it with the tensor product of its axes' units. The integration co-vector
        :math:`\langle 1, \cdot\rangle_V` is :attr:`integral`.

        ELEGANCE-DEBT[guard] #707 legacy-leaf: the leaf unit is ones; retires at N07, when no space is a leaf.
        """
        return np.ones(self.shape)

    @property
    def integral(self) -> Covector:
        r"""The integration co-vector :math:`\int \cdot\,dV = \langle 1, \cdot\rangle_V = \flat_V 1`.

        Minted through :meth:`flat` from the constant function as a typed
        member of this space (its coordinates are :attr:`unit`), so the
        metric is read once, at the mint, and the number one is never
        spelled ``np.ones``: on a frame-dressed head the constant function's
        coordinates are its moments, and ``flat(ones)`` is not the integral.
        """
        return self.flat(Member(values=self.unit, space=self))

    # ------------------------------------------------------------------
    # Metric application (Hilbert adjoint building blocks, Wave O / O.2b)
    # ------------------------------------------------------------------

    def metric_diagonal(self) -> NDArray | None:
        r"""The diagonal of the metric over this space's coordinates, or ``None`` if it is not diagonal.

        :math:`G\mathbf 1`, of :attr:`shape`: all ones on a Euclidean space,
        the resolved metric's :meth:`~orpheus.numerics.metric.HilbertMetric.diagonal`
        otherwise. A consumer that needs node weights (a frame's measure on a
        subset of the coordinates, the kernel gauge's summands) reads them
        here. Composite spaces override it with the concatenation of their
        blocks' diagonals in their flat layout.
        """
        metric = self._resolved_metric
        return np.ones(self.shape) if metric is None else metric.diagonal(self.shape)

    def apply_metric(self, x: Carrier) -> Carrier:
        r"""Apply the Hilbert metric :math:`G\odot x` (identity if Euclidean).

        Carrier-generic surface; the base realization
        (:meth:`_diagonal_apply_metric`) delegates to the space's resolved
        :class:`~orpheus.numerics.metric.HilbertMetric` — a diagonal
        weight broadcast against the leading axes of a bare-array ``x``,
        or a dense matrix on its flattened leading block. This is the
        building block :class:`~orpheus.numerics.operator.AdjointOperator`
        applies to the codomain before the transpose. Composite spaces
        (bulk :math:`\oplus` trace) OVERRIDE this to apply a per-block metric
        to a structured field (the Wave-O direct-sum adjoint).
        """
        return self._diagonal_apply_metric(self.leading_coordinates_of(x))

    def _diagonal_apply_metric(self, x: Any) -> Any:
        r"""The bare-array realization of :meth:`apply_metric`."""
        m = self._resolved_metric
        if m is None:
            return x
        return m.apply(x)

    def apply_inverse_metric(self, x: Carrier) -> Carrier:
        r"""Apply the Moore–Penrose pseudo-inverse metric :math:`G^{+}\odot x`.

        Carrier-generic surface (see :meth:`apply_metric`); the base
        realization is :meth:`_diagonal_apply_inverse_metric`, delegating
        to the resolved :class:`~orpheus.numerics.metric.HilbertMetric`.
        For the diagonal realization this is
        ``(1/G)⊙x`` where ``G ≠ 0``, and ``0`` on the metric's null space
        (``G = 0`` — e.g. the tangential partial-current trace slots where
        ``|Ω·n| = 0``). The pseudo-inverse is exact for the Hilbert adjoint:
        the null-space components carry zero ``⟨·,·⟩_G`` weight and are zero
        in any matvec output by construction. Identity if Euclidean. Applied
        to the domain after the transpose. For a strictly-positive metric
        (e.g. the angular quadrature weights ``w_n``) this is plain ``x/G``,
        so the spherical-harmonic adjoint path stays bit-identical.
        """
        return self._diagonal_apply_inverse_metric(self.leading_coordinates_of(x))

    def _diagonal_apply_inverse_metric(self, x: Any) -> Any:
        r"""The bare-array realization of :meth:`apply_inverse_metric`."""
        m = self._resolved_metric
        if m is None:
            return x
        return m.apply_inverse(x)

    @property
    def n_coordinates(self) -> int:
        r"""The length of the space's flat chart, :math:`\dim V = \prod` :attr:`shape`.

        A member's coordinates are an array of :attr:`shape`; their C-order
        ravel is a vector of :math:`\mathbb R^{n}` with this :math:`n`. A
        matrix realising a map :math:`V \to W` is
        ``(W.n_coordinates, V.n_coordinates)`` (:func:`admit_matrix`).
        """
        return int(np.prod(self.shape))

    @property
    def rank(self) -> int:
        r"""The dimension of :math:`V / \ker G`, on which this space's inner
        product is non-degenerate (numerics_untangling.md ruling 74): the
        element size when the metric is non-degenerate (or Euclidean), less
        when it has a kernel (the tangential trace slots, dead moment-head
        slots, a rank-deficient Parseval head). The adjoint is posed on the
        quotient, so its laws hold modulo the kernel."""
        metric = self._resolved_metric
        if metric is None:
            return self.n_coordinates
        return metric.rank(self.shape)

    def kernel_basis(self) -> NDArray:
        r"""An orthonormal basis of :math:`\ker G`, the ``(n, n − rank)`` columns
        in the flat C-order chart (:meth:`HilbertMetric.kernel_basis`; none on a
        Euclidean space). A composite space overrides it block-wise, in its flat
        chart's block order, as it overrides :attr:`rank`."""
        metric = self._resolved_metric
        if metric is None:
            return np.zeros((self.n_coordinates, 0))
        return metric.kernel_basis(self.shape)

    def member(self, coordinates: NDArray) -> Carrier:
        r"""The member whose coordinates are ``coordinates``: the inverse of the
        coordinate chart, so ``coordinates_of(member(c)) == c`` and
        ``member(c).space == self`` (the same object on a leaf and its dual; a
        composite field rebuilds an equal space per instance).

        The class follows the space (#707 ruling R5): on a primal leaf it is a
        bare :class:`Member`, an element of :math:`V` with no physical role;
        :meth:`DualSpace.member` mints a :class:`Covector`, the integrated data
        side. A space whose members are typed above this layer (a composite)
        overrides this with its builder-wired member class.
        """
        return cast(Carrier, Member(values=self._chart_values(coordinates), space=self))

    def _chart_values(self, coordinates: NDArray) -> NDArray:
        """The coordinates as a float array of :attr:`shape` (flat or shaped)."""
        return np.asarray(coordinates, dtype=float).reshape(self.shape)

    @property
    def _resolved_metric(self) -> Optional[HilbertMetric]:
        r"""The metric SOURCE of a leaf resolved to its realization, PER CALL: an explicit ``metric`` object, else the ``inner_product_weights``
        array as a :class:`~orpheus.numerics.metric.DiagonalMetric`, else ``None`` (Euclidean: the verbs short-circuit without an object).

        ELEGANCE-DEBT[guard] #707 legacy-leaf: the leaf resolution; retires at N07, when the leaf fields leave the base.

        ⛔ Deliberately NOT a ``cached_property``: the mutation-battery idiom mutates a frozen space's weight FIELD in place
        (``object.__setattr__``), and the adjoint-certification propagation probes assert the mutation reaches the metric surface at the
        next read (`[M]` 2026-08-30, P7 exit gate).
        """
        if self.metric is not None:
            return self.metric
        if self.inner_product_weights is not None:
            return DiagonalMetric(self.inner_product_weights)
        return None

    def dual(self) -> "FunctionSpace":
        r"""Return the dual space :math:`V^*`, the space of covectors on this one.

        :math:`V^*` has this space's coordinates and the inverse Riesz form
        :math:`G^{+}` as its metric (numerics_untangling.md ruling 26, issue
        #617), and it is a distinct space: ``V.dual() != V`` (ruling 75), so
        the composition guard refuses a covector where a vector is expected.
        A leaf's dual is a :class:`DualSpace`; a composite's is the composite
        of its members' duals (biproduct duality, :meth:`_make_dual`).

        Memoised per space, since the dual is a function of the space alone:
        ``V.dual() is V.dual()``, and ``V.dual().dual() is V`` (the
        involution :math:`V^{**} = V`, by object identity). `[M]` 2026-10-21,
        the S6a spec: one slab adjoint solve minted 5 645 DualSpaces (10 965
        at tight tolerances) from the dual arrows each matvec rebuilds; the
        memo bounds the count by the distinct spaces asked for.

        The return type is :class:`FunctionSpace`, not :class:`DualSpace`,
        because ``dual`` is reflexive: applied to a :class:`DualSpace` it
        returns the primal, which is any space, and a composite's dual is a
        composite.
        """
        return self._dual

    @cached_property
    def _dual(self) -> "FunctionSpace":
        dual = self._make_dual()
        # The involution as object identity: the dual's own dual is this
        # space, seeded on the dual unless it already knows its dual (a
        # DualSpace's dual is its primal, never overwritten).
        dual.__dict__.setdefault("_dual", self)
        return dual

    def _make_dual(self) -> "FunctionSpace":
        r"""Build :math:`V^*` (called once per space, by :meth:`dual`)."""
        return DualSpace(name=f"{self.name}*", shape=self.shape, primal=self)



@dataclass(frozen=True, eq=False)
class DualSpace(FunctionSpace):  # concept: dual | role: head | kind: space
    r"""The dual :math:`V^*` of a leaf :class:`FunctionSpace`: its covectors.

    :math:`V^*` has the primal's coordinates and carries the inverse Riesz
    form :math:`G^{+}` as its metric (numerics_untangling.md ruling 26,
    issue #617): its metric is the primal's metric read the other way round
    (:meth:`~orpheus.numerics.metric.HilbertMetric.dual`), never a copy of
    :math:`G`, so the Riesz legs are posed on it and
    :math:`\langle \phi, \chi\rangle_{V^*} = \langle \sharp\phi, \sharp\chi\rangle_V`.
    Its content IS its primal: it copies no weights, axes or metric object,
    and it equals another dual iff their primals are equal, and never equals
    a primal space (ruling 75).

    Built only through :meth:`FunctionSpace.dual` (memoised per space), whose
    :meth:`of` constructor it keeps; ``V.dual().dual() is V``.

    Parameters
    ----------
    primal : FunctionSpace
        The primal space :math:`V` of which this is the dual.
    """

    # Required (a DualSpace WITHOUT its primal is an illegal state);
    # kw_only sidesteps the inherited-defaults field-ordering rule.
    primal: "FunctionSpace" = field(kw_only=True, compare=False, repr=False)

    def identity_key(self) -> tuple[Any, ...]:
        """Covariance is identity (ruling 75): a dual is tagged and keyed by
        its primal, so it equals another dual iff their primals are equal and
        never equals a primal space."""
        return ("dual", self.primal.identity_key())

    def __post_init__(self) -> None:
        # Its content IS its primal: a dual carries the primal's shape and no
        # metric source of its own (the inverse form is read through the
        # primal), so a dual that disagrees with its primal is an illegal
        # state.
        if self.shape != self.primal.shape:
            raise ValueError(
                f"DualSpace {self.name!r}: shape {self.shape} is not its "
                f"primal's {self.primal.shape}."
            )
        if self.metric is not None or self.inner_product_weights is not None:
            raise ValueError(
                f"DualSpace {self.name!r} carries a metric source of its own; "
                f"a dual's metric is its primal's inverse form, read through "
                f"the primal."
            )
        super().__post_init__()

    @classmethod
    def of(cls, primal: "FunctionSpace") -> "FunctionSpace":
        r"""The dual of ``primal``: :meth:`FunctionSpace.dual`, memoised.

        Idempotent: ``of(of(V)) is V``. A composite's dual is the composite
        of its members' duals, never a :class:`DualSpace` over it.
        """
        return primal.dual()

    def _make_dual(self) -> "FunctionSpace":
        return self.primal  # V** = V

    def member(self, coordinates: NDArray) -> Covector:
        r"""The co-vector with these integrated coordinates: the data side's
        member class (#707 ruling R5). The inverse of :meth:`coordinates_of`,
        the raw mint of :meth:`FunctionSpace.flat`'s output."""
        return Covector(values=self._chart_values(coordinates), space=self)

    @property
    def _resolved_metric(self) -> Optional[HilbertMetric]:
        r"""The inverse Riesz form: the primal's resolved metric, dualised.

        Read per call through the primal, so a battery that mutates the
        primal's weight field in place reaches the dual too (the reason
        :meth:`FunctionSpace._resolved_metric` is not cached).
        """
        metric = self.primal._resolved_metric
        return None if metric is None else metric.dual()

    def __repr__(self) -> str:
        return f"DualSpace({self.name!r}, shape={self.shape})"



@dataclass(frozen=True, eq=False)
class TensorSpace(FunctionSpace):  # concept: tensor | role: head | kind: space
    r"""The tensor space of factors :math:`V = F_1 \otimes \dots \otimes F_k`: one value of one type, for every product of slots.

    A slot is a :class:`~orpheus.numerics.factor.Factor`: an axis (the index set, its form and its origin), a variance (a function or a
    measure, a density or an integrated chart) and a subquotient :math:`S/K`. The space is its ordered tuple of factors, so ``of_axes(a, b)``
    puts ``a`` first, the shape is the concatenation of the axes' shapes, and the metric is the Kronecker product of the factors' forms in
    slot order (:meth:`flat_on`, :meth:`sharp_on`): the one place the order of the application is spelled. A slot whose form is
    :class:`~orpheus.numerics.metric.NoForm` is skipped by the walk and counts in :meth:`inner_product` and :meth:`norm` (the energy slot
    of a flux until N11); a counting form is applied, as a product with one.

    **The Klein action per slot.** :meth:`dual` flips the variance of every slot (kind and chart, and dualises the subquotient), :meth:`riesz`
    flips the chart only, and each takes ``on=`` a set of labels, so a space may be dual on one slot and primal on another: the two
    involutions commute and ``dual(on=L).dual(on=M) == dual(on=L ^ M)``. ``riesz`` is refused on a formless slot: a change of chart needs a
    form to pass through.

    **The member class** follows the space (issue #707): a member of ``V`` is a :class:`Member`, a member of ``V.dual()`` a :class:`Covector`
    (the data side). The pair is the space's :class:`CarrierSpec`, not part of its identity.

    **Constructors.** ``of_axes(*factors)`` takes factors, so every call states its variances; ``of_measure``, ``of_basis`` and ``declared``
    are one-factor constructors whose body is ``of_axes(Factor(axis, variance, Subquotient.whole(n)))``; ``a * b`` concatenates the slot
    lists; ``sub(*labels)`` is the sub-product on those labels in slot order, of which ``tail``, ``without_axis`` and the head read are calls.
    The labels are unique (``ValueError``), and the empty product is refused.

    Parameters
    ----------
    factors : tuple[Factor, ...]
        The slots, in order.
    carrier : CarrierSpec | None, keyword-only, default None
        The member classes; ``None`` is the default pair ``(Member, Covector)``.
    """

    # The leaf fields of the base are fixed here, so that a stored shape or a second metric cannot be passed: ``name`` and ``shape`` are
    # derived from the factors and the other two stay ``None``.
    # ELEGANCE-DEBT[guard] #707 legacy-leaf: the four leaf fields are re-declared closed until N07 removes them from the base.
    name: str = field(init=False)
    shape: tuple[int, ...] = field(init=False)
    inner_product_weights: Optional[NDArray] = field(init=False, default=None, repr=False, compare=False)
    metric: Optional[HilbertMetric] = field(init=False, default=None, repr=False, compare=False, kw_only=True)
    factors: tuple[Factor, ...]
    carrier: Optional[CarrierSpec] = field(default=None, repr=False, compare=False, kw_only=True)

    def __post_init__(self) -> None:
        factors = tuple(self.factors)
        if not factors:
            raise ValueError("a tensor space needs at least one factor: the empty product has no index set")
        for f in factors:
            if not isinstance(f, Factor):
                raise TypeError(
                    f"a tensor space is built from Factors (state the variance of each slot), got {type(f).__name__}"
                )
        if self.carrier is not None and not isinstance(self.carrier, CarrierSpec):
            raise TypeError(f"carrier must be a CarrierSpec or None, got {type(self.carrier).__name__}")
        labels = [f.axis.label for f in factors]
        repeated = sorted({l for l in labels if labels.count(l) > 1})
        if repeated:
            raise ValueError(f"a tensor space has unique labels; {repeated} repeat in {labels}")
        shape: tuple[int, ...] = ()
        for f in factors:
            shape = shape + f.axis.shape
        readable = "*".join(f"{f.axis.label}{f.axis.shape}" for f in factors)
        object.__setattr__(self, "factors", factors)
        object.__setattr__(self, "shape", shape)
        object.__setattr__(self, "name", f"{readable}#{name_digest(factors)}")

    def identity_key(self) -> tuple[Any, ...]:
        r"""The space is its ordered factors: ``("tensor", factors)``, so the variance, the subquotient and the slot order are identity."""
        return ("tensor", self.factors)

    def __repr__(self) -> str:
        return f"TensorSpace({self.name!r}, shape={self.shape})"

    # ------------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------------

    @classmethod
    def of_axes(cls, *factors: Factor) -> "TensorSpace":
        r"""The ordered product of ``factors``: the one composition. It takes :class:`Factor` s (an :class:`Axis` is a ``TypeError``)."""
        return cls(factors=tuple(factors))

    @classmethod
    def of_measure(cls, measure: "DiscreteMeasure", label: str, variance: FormlessVariance | Variance) -> "TensorSpace":
        r"""The discrete :math:`L^2(\mu)` space a measure samples: one factor over :meth:`Axis.of_measure
        <orpheus.numerics.axis.Axis.of_measure>`, whole, in the stated variance."""
        axis = Axis.of_measure(measure, label)
        return cls.of_axes(Factor(axis, variance, Subquotient.whole(int(np.prod(axis.shape)))))

    @classmethod
    def of_basis(cls, basis: "Basis | TruncatedBasis", variance: FormlessVariance | Variance) -> "TensorSpace":
        r"""The coefficient space a basis spans: one factor over :meth:`Axis.of_basis <orpheus.numerics.axis.Axis.of_basis>`, whole, in the
        stated variance."""
        axis = Axis.of_basis(basis)
        return cls.of_axes(Factor(axis, variance, Subquotient.whole(int(np.prod(axis.shape)))))

    @classmethod
    def declared(
        cls, name: str, shape: tuple[int, ...], variance: FormlessVariance | Variance, *, form: DiagonalMetric | DenseMetric | NoForm | None = None,
    ) -> "TensorSpace":
        r"""The space of a declared index set no generating object minted: one nodal factor, whole, in the stated variance.

        ``form`` is the factor's form, carried by the axis's ``Declared`` origin; ``None`` is the counting form over ``shape``.
        """
        axis = Axis(name, tuple(int(d) for d in shape), origin=Declared(BasisKind.NODAL, form))
        return cls.of_axes(Factor(axis, variance, Subquotient.whole(int(np.prod(axis.shape)))))

    # ------------------------------------------------------------------
    # The slots
    # ------------------------------------------------------------------

    @property
    def labels(self) -> tuple[str, ...]:
        r"""The labels of the slots, in slot order."""
        return tuple(f.axis.label for f in self.factors)

    @property
    def axes(self) -> tuple[Axis, ...]:
        r"""The factors' axes, in slot order (a restricted axis is returned only where the space was built with one)."""
        return tuple(f.axis for f in self.factors)

    @property
    def formed_labels(self) -> tuple[str, ...]:
        r"""The labels of the slots that carry a form (every slot but the formless ones), in slot order."""
        return tuple(f.axis.label for f in self.factors if not isinstance(f.variance, FormlessVariance))

    def _axis_index(self, label: str) -> int:
        labels = self.labels
        hits = [i for i, l in enumerate(labels) if l == label]
        if len(hits) != 1:
            raise ValueError(f"axis lookup: label {label!r} names {len(hits)} axes of {self!r} (have {list(labels)}) — need exactly one.")
        return hits[0]

    def position(self, label: str) -> int:
        r"""The slot index of the factor labeled ``label``."""
        return self._axis_index(label)

    def factor(self, label: str) -> Factor:
        r"""The factor labeled ``label``."""
        return self.factors[self._axis_index(label)]

    def axis(self, label: str) -> Axis:
        r"""The axis of the factor labeled ``label``: ``space.axis("energy").shape[0]`` is the group count."""
        return self.factor(label).axis

    def _axis_block(self, k: int) -> DiagonalMetric | DenseMetric | DualMetric | None:
        form = self.factors[k].form
        return None if isinstance(form, NoForm) else cast("DiagonalMetric | DenseMetric | DualMetric", form)

    def _labels(self, on: Any) -> frozenset[str]:
        r"""``on`` parsed as a collection of labels of this space (a ``str`` is refused: it would read as its characters)."""
        if isinstance(on, str):
            raise TypeError(f"on= takes a collection of labels, not the string {on!r}")
        labels = frozenset(on)
        unknown = sorted(labels - set(self.labels))
        if unknown:
            raise KeyError(f"unknown labels {unknown}; the space has {list(self.labels)}")
        return labels

    def _again(self, factors: tuple[Factor, ...], carrier: Optional[CarrierSpec]) -> "TensorSpace":
        return TensorSpace(factors=factors, carrier=carrier)

    def with_axis(self, factor: Factor) -> "TensorSpace":
        r"""This product with ``factor`` appended last, every other factor unchanged. It takes a :class:`Factor`, so a variance is stated."""
        return self._again(self.factors + (factor,), self.carrier)

    def sub(self, *labels: str) -> "TensorSpace":
        r"""The sub-product on ``labels``, in slot order: the one verb of which ``tail``, ``without_axis`` and the head read are calls.

        An unknown label is a ``KeyError``, a label list not in slot order (or repeating a label, or empty) a ``ValueError``.
        """
        known = self._labels(labels)
        positions = [self._axis_index(l) for l in labels]
        if not labels or positions != sorted(set(positions)) or len(known) != len(labels):
            raise ValueError(f"sub: the labels {list(labels)} are not a non-empty sub-product, in slot order, of {list(self.labels)}")
        return self._again(tuple(self.factors[i] for i in positions), self.carrier)

    def tail(self, k: int) -> "TensorSpace":
        r"""The tensor space of the last ``k`` factors, with their variances."""
        if isinstance(k, bool) or not isinstance(k, int) or not 1 <= k <= len(self.factors):
            raise ValueError(f"tail takes 1 <= k <= {len(self.factors)} (an int), got {k!r}")
        return self.sub(*self.labels[len(self.factors) - k:])

    def without_axis(self, label: str) -> "TensorSpace":
        r"""The product of the OTHER factors, each keeping its variance, form and subquotient."""
        self._axis_index(label)
        return self.sub(*(l for l in self.labels if l != label))

    def __mul__(self, other: "FunctionSpace") -> "TensorSpace":
        r"""The product :math:`V \otimes W`: the concatenation of the slot lists, associative; ``NotImplemented`` against a leaf.

        The product's carrier is the factors' when they agree and the default pair when they differ.
        """
        if not isinstance(other, TensorSpace):
            return NotImplemented
        return self._again(self.factors + other.factors, self.carrier if self.carrier == other.carrier else None)

    def positions_spanned_by(self, sub: "FunctionSpace") -> tuple[int, ...]:
        r"""The array positions of this space that the factors of ``sub`` occupy, read by label.

        ``sub`` is a sub-product of this space's factors, in this space's order; an operator acting on it acts on this space as
        :math:`A \otimes I` on the other positions, and these are the positions :math:`A` acts on. A factor may span several array
        positions (a 2-D spatial mesh), so the answer is in array positions, not factor indices.
        """
        if not isinstance(sub, TensorSpace):
            raise TypeError(
                f"positions_spanned_by: the space {sub!r} is not built from labelled factors, so the positions it spans cannot be read off {self!r}."
            )
        labels = sub.labels
        ordered = tuple(l for l in self.labels if l in labels)
        if ordered != labels:
            raise ValueError(
                f"positions_spanned_by: {sub!r} (axes {labels}) is not a sub-product, in order, of {self!r} (axes {self.labels})."
            )
        positions: list[int] = []
        dimension = 0
        for axis in self.axes:
            rank = len(axis.shape)
            if axis.label in labels:
                positions.extend(range(dimension, dimension + rank))
            dimension += rank
        return tuple(positions)

    @property
    def has_coordinate_cone(self) -> bool:
        return all(ax.kind is BasisKind.NODAL for ax in self.axes)

    @property
    def unit(self) -> NDArray:
        r"""The coordinates of the constant function :math:`1`: the tensor product of the axes' units, :math:`u = u_1 \otimes \dots \otimes u_k`
        (ones on a nodal axis, the carried unit on a modal one, which on a frame head is the moments :math:`\Phi^{\mathsf T} w`)."""
        unit = np.ones(())
        for ax in self.axes:
            unit = np.multiply.outer(unit, np.asarray(ax.unit))
        return unit.reshape(self.shape)

    # ------------------------------------------------------------------
    # The metric: the Kronecker product of the forms, in slot order
    # ------------------------------------------------------------------

    @cached_property
    def _metrics(self) -> dict[frozenset[str], FactoredMetric]:
        return {}

    def _metric_on(self, labels: frozenset[str]) -> FactoredMetric:
        r"""The Kronecker product of the forms of the slots in ``labels``, in slot order; the other slots carry no form (``NoForm``)."""
        metric = self._metrics.get(labels)
        if metric is None:
            metric = FactoredMetric(tuple(
                (f.axis.shape, f.form if f.axis.label in labels else NoForm()) for f in self.factors
            ))
            self._metrics[labels] = metric
        return metric

    def flat_on(self, x: Any, labels: Any) -> Any:
        r"""The product of the forms of the slots in ``labels`` applied to ``x``, one slot at a time in slot order (2.4); a ``NoForm`` slot is skipped.

        The one body of the forward walk: ``flat``, ``apply_metric``, ``integral`` and the pairing are calls of it.
        """
        return self._metric_on(self._labels(labels)).apply(x)

    def sharp_on(self, x: Any, labels: Any) -> Any:
        r"""The inverse of :meth:`flat_on`: the inverse faces of the forms of the slots in ``labels``, in slot order (the masked reciprocal on a
        diagonal form)."""
        return self._metric_on(self._labels(labels)).apply_inverse(x)

    @property
    def _resolved_metric(self) -> HilbertMetric:
        return self._metric_on(frozenset(self.labels))

    def apply_metric(self, x: Any) -> Any:
        r"""The metric :math:`G x`: :meth:`flat_on` over every slot (a ``NoForm`` slot is skipped by the walk)."""
        return self.flat_on(self.leading_coordinates_of(x), self.labels)

    def apply_inverse_metric(self, x: Any) -> Any:
        r"""The pseudo-inverse metric :math:`G^{+} x`: :meth:`sharp_on` over every slot (a ``NoForm`` slot is skipped by the walk)."""
        return self.sharp_on(self.leading_coordinates_of(x), self.labels)

    def flat(self, density: Any) -> Covector:
        r"""The covector :math:`\flat_V w = G w`: :meth:`flat_on` over every slot, a member of the dual space."""
        return Covector(values=np.asarray(self.flat_on(self.coordinates_of(density), self.labels)), space=self.dual())

    @property
    def integral(self) -> Covector:
        r"""The integration covector :math:`\langle 1, \cdot\rangle_V = \flat_V 1`, minted through :meth:`flat` from the constant function."""
        return self.flat(Member(values=self.unit, space=self))

    # ------------------------------------------------------------------
    # The two involutions, per slot
    # ------------------------------------------------------------------

    def dual(self, on: Any = None) -> "TensorSpace":
        r"""The dual on the slots ``on`` (every slot by default): each of those factors dualised, the others unchanged.

        ``dual(on=L).dual(on=M) == dual(on=L ^ M)`` (law L4) and ``dual(on=all) is dual()``, memoised (so ``V.dual().dual() is V``);
        ``on=()`` is the space itself; an unknown label is a ``KeyError``. The member class is a function of the factors: a flip that lands on the
        space itself or on its all-slot dual returns that memoised object (so composed partial duals reach ``V.dual()`` with its carrier), and a
        genuinely partial dual carries the default pair, whichever route reached it.
        """
        return self._flip(on, Factor.dual)

    def riesz(self, on: Any = None) -> "TensorSpace":
        r"""The Riesz image on the slots ``on`` (every slot by default): the chart of each flipped through its form (law L5).

        Refused, "riesz is not defined: ", where a slot of ``on`` is formless: a change of chart needs a form to pass through. The member class
        follows the factors, as for :meth:`dual`.
        """
        labels = frozenset(self.labels) if on is None else self._labels(on)
        formless = [f.axis.label for f in self.factors if f.axis.label in labels and isinstance(f.variance, FormlessVariance)]
        if formless:
            raise ValueError(f"riesz is not defined: the slots {formless} have no form to change the chart through")
        return self._flip(on, Factor.riesz)

    def _flip(self, on: Any, per_factor: Callable[[Factor], Factor]) -> "TensorSpace":
        r"""The one core of :meth:`dual` and :meth:`riesz`: ``per_factor`` on each slot of ``on``, the others unchanged, in the carrier that
        the resulting factors name (the space itself, its memoised all-slot dual, or the default pair)."""
        labels = frozenset(self.labels) if on is None else self._labels(on)
        factors = tuple(per_factor(f) if f.axis.label in labels else f for f in self.factors)
        if factors == self.factors:
            return self
        flipped = cast("TensorSpace", self._dual)
        return flipped if factors == flipped.factors else self._again(factors, None)

    def _make_dual(self) -> "TensorSpace":
        carrier = (self.carrier if self.carrier is not None else CarrierSpec()).dual()
        return self._again(tuple(f.dual() for f in self.factors), carrier)

    def member(self, coordinates: NDArray) -> Any:
        r"""The member whose coordinates are ``coordinates``, of the class its carrier names (a :class:`Member` on a primal, a
        :class:`Covector` on the all-slot dual)."""
        carrier = self.carrier if self.carrier is not None else CarrierSpec()
        return carrier.member(values=self._chart_values(coordinates), space=self)

    # ------------------------------------------------------------------
    # The completions
    # ------------------------------------------------------------------

    def completed(self, on: Any) -> "TensorSpace":
        r"""This space with the forms of the slots ``on`` replaced by the counting form: the one core of :meth:`euclidean` and :meth:`counting`.

        A replaced slot is the counting axis over the same shape, whole, and keeps its kind and chart (a formless slot takes the
        :class:`Variance` its :class:`FormlessVariance` implies); the counting form is applied by the walk, as a product with one.
        """
        labels = self._labels(on)
        if not labels:
            return self
        return self._again(tuple(_counting_factor(f) if f.axis.label in labels else f for f in self.factors), self.carrier)

    def euclidean(self) -> "TensorSpace":
        r"""The completion of the formless slots by the counting form (the energy slot); the angular and spatial forms are kept."""
        return self.completed(on=frozenset(l for l in self.labels if l not in self.formed_labels))

    def counting(self) -> "TensorSpace":
        r"""The completion of every slot by the counting form."""
        return self.completed(on=self.labels)


def _counting_factor(f: Factor) -> Factor:
    r"""The slot ``f`` with its axis replaced by the counting axis over the same shape, whole, in the variance it had.

    A formless variance is re-stated as the literal :class:`Variance` it implies (the production-mint gate requires every mint to state its variance
    literally, so the map :attr:`FormlessVariance.kind` and :attr:`FormlessVariance.chart` define is written out here, not read).
    """
    shape = f.axis.shape
    axis = Axis(f.axis.label, shape, origin=Declared(BasisKind.NODAL, DiagonalMetric.counting(shape)))
    whole = Subquotient.whole(int(np.prod(shape)))
    if isinstance(f.variance, FormlessVariance):
        if f.variance is FormlessVariance.FUNCTION:
            return Factor(axis, Variance(Kind.FUNCTION, VarianceChart.DENSITY), whole)
        return Factor(axis, Variance(Kind.MEASURE, VarianceChart.INTEGRATED), whole)
    return Factor(axis, f.variance, whole)

# ---------------------------------------------------------------------------
# Pre-populated common-space factories
# ---------------------------------------------------------------------------
#
# These are *factory functions*, not module-level singletons, because
# the shape of an angular-flux space depends on the mesh dimensions —
# every solver instance carries its own (n_cells, n_ordinates, n_groups)
# triple. Factories give the caller named, well-typed construction
# sites without forcing premature commitment to a global instance.


def act_on_direct_sum_coordinates(
    space: FunctionSpace,
    members: "tuple[FunctionSpace, ...]",
    x: Any,
    leg: "Callable[[FunctionSpace, NDArray], NDArray]",
) -> NDArray:
    r"""A direct sum's action on its FLAT coordinates: each member's ``leg`` on its own slice.

    The coordinates of :math:`V_1 \oplus \dots \oplus V_N` are the members'
    coordinates concatenated in member order (the layout the composite carriers'
    ``to_flat`` packs). An action that is block-diagonal on the sum (a metric
    leg, :math:`G = \operatorname{diag}(G_1, \dots, G_N)`) acts on that chart by
    cutting it at the member sizes, handing each member its slice in the
    member's own shape, and concatenating the results.

    The operand is parsed at ``space``'s ACTION verb
    (:meth:`FunctionSpace.leading_coordinates_of`), so a trailing payload rides
    through as :math:`G \otimes I_P`. This is the chart the operator algebra
    hands a composite's metric legs (an operator's ``.H`` is
    :math:`\sharp \circ \mathrm{dual} \circ \flat` on flat coordinates); the
    composite's typed carrier is the other chart of the same element.
    """
    coordinates = np.asarray(space.leading_coordinates_of(x))
    payload = coordinates.shape[1:]
    pieces = []
    start = 0
    for member in members:
        size = member.n_coordinates
        block = coordinates[start:start + size].reshape(tuple(member.shape) + payload)
        pieces.append(np.asarray(leg(member, block)).reshape((size,) + payload))
        start += size
    return np.concatenate(pieces, axis=0)


def coordinate_space_of(shape: tuple[int, ...]) -> FunctionSpace:
    r"""The coordinate space :math:`\mathbb R^{\text{shape}}` (Euclidean): what an UNPOSED matrix carrier acts on.

    A carrier posed on no space still acts between coordinate spaces; this is
    that space, so a posed and an unposed carrier go through one
    :func:`act_on_flat_coordinates`. Minted only by :func:`admit_matrix`.
    """
    return TensorSpace.declared(
        "R^" + "x".join(str(d) for d in shape), tuple(int(d) for d in shape), Variance(Kind.FUNCTION, VarianceChart.DENSITY),
    )


def act_on_flat_coordinates(
    source: FunctionSpace,
    target: FunctionSpace,
    x: Any,
    action: "Callable[[NDArray], NDArray]",
) -> NDArray:
    r"""A map :math:`V \to W` held as a matrix, acting on a member: ``x`` in, an element of ``target`` out.

    A matrix :math:`M` realises a linear map between two spaces in their flat
    charts (the C-order ravel of the coordinates, :attr:`FunctionSpace.n_coordinates`
    long). Acting on a member is three steps, and this is their one
    definition: parse ``x`` at ``source``'s ACTION verb
    (:meth:`FunctionSpace.leading_coordinates_of`, so a trailing payload
    :math:`P` rides through as :math:`M \otimes I_P`), ravel the leading axes
    into the chart, act, and return the result in ``target``'s shape. The
    result is an element of ``target`` because it takes ``target``'s exact
    shape (#683).

    A matrix carrier's transpose is the SAME call with the ends swapped,
    :math:`M^{\mathsf T}` from ``target`` to ``source``, so the forward and the
    transpose arm cannot disagree on which space each end is.
    """
    coordinates = np.asarray(source.leading_coordinates_of(x))
    payload = coordinates.shape[len(source.shape):]
    if not payload:
        return np.asarray(action(coordinates.reshape(source.n_coordinates))).reshape(target.shape)
    # The payload folds into ONE column axis: a matrix action is M @ X on an
    # (n, p) block, and a higher-rank block is not that (LAPACK's lu_solve
    # batches its leading axes, so an (n, p, q) block is solved wrongly).
    columns = coordinates.reshape(source.n_coordinates, int(np.prod(payload)))
    return np.asarray(action(columns)).reshape(tuple(target.shape) + payload)


def admit_matrix(
    shape: tuple[int, int],
    domain: "FunctionSpace | None",
    codomain: "FunctionSpace | None",
    *,
    owner: str,
    unposed: "tuple[tuple[int, ...], tuple[int, ...]] | None" = None,
) -> "tuple[FunctionSpace, FunctionSpace]":
    r"""Parse a matrix and its ends into the two charts it acts between, ``(source, target)``.

    The one admission of every carrier that holds a matrix. A matrix carrier
    is posed on both ends or on neither, as the adjoint is
    (:func:`~orpheus.numerics.operator.posed`):

    * **posed**: the charts are ``domain`` and ``codomain``, and the matrix
      must be ``(codomain.n_coordinates, domain.n_coordinates)``;
    * **unposed**: the charts are the coordinate spaces
      (:func:`coordinate_space_of`) of the ``unposed`` shapes the carrier
      declares, ``(source shape, target shape)``, by default the matrix's own
      :math:`\mathbb R^{\text{cols}}` and :math:`\mathbb R^{\text{rows}}`;
      the matrix must map one into the other.

    The charts are what :func:`act_on_flat_coordinates` acts between, and they
    are read off the same :attr:`FunctionSpace.n_coordinates`, so the
    admission and the action cannot disagree.

    Raises
    ------
    ValueError
        One end posed and the other not, or a matrix that does not map the
        source chart into the target chart.
    """
    if (domain is None) != (codomain is None):
        raise ValueError(
            f"{owner}: a matrix carrier is posed on both ends or on neither; got domain={domain!r}, "
            f"codomain={codomain!r}."
        )
    if domain is not None and codomain is not None:
        source, target = domain, codomain
    else:
        source_shape, target_shape = unposed if unposed is not None else ((shape[1],), (shape[0],))
        source, target = coordinate_space_of(source_shape), coordinate_space_of(target_shape)
    expected = (target.n_coordinates, source.n_coordinates)
    if tuple(shape) != expected:
        raise ValueError(
            f"{owner}: a {shape[0]}×{shape[1]} matrix does not map {source!r} into {target!r}, "
            f"which needs {expected[0]}×{expected[1]} (codomain.n_coordinates × domain.n_coordinates)."
        )
    return source, target


def factor_tail(width: int) -> tuple[int, ...]:
    r"""The shape suffix a factor of ``width`` coefficients adds to a composed space.

    ``()`` for a one-dimensional factor, ``(width,)`` otherwise. A
    one-dimensional factor is :math:`\mathbb R`, and
    :math:`V \otimes \mathbb R \cong V`, so composing with it leaves the space
    (and every buffer on it) unchanged: a cell-average closure's bulk and face
    buffers stay byte-identical to the moment-free ones (#240 D5b). The one rule
    behind every trailing spatial-moment axis, on cells and on faces alike.
    """
    return () if width == 1 else (width,)

