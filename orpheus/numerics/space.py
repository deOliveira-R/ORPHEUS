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
product, :class:`TensorProductSpace`) and ``S.dual()`` (the dual space,
:class:`DualSpace`) ship on :class:`FunctionSpace`; the direct sum is not
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

from functools import cached_property
from dataclasses import dataclass, field, replace
from typing import TYPE_CHECKING, Any, Generic, Optional, TypeVar

import numpy as np
from numpy.typing import NDArray

from .axis import Axis, BasisKind
from .content import content_digest, name_digest
from .metric import (
    DenseMetric,
    DiagonalMetric,
    DualMetric,
    FactoredMetric,
    HilbertMetric,
    PositionedMetric,
)

if TYPE_CHECKING:
    from .basis.base import Basis, TruncatedBasis
    from .measure import DiscreteMeasure


__all__ = [
    "DualSpace",
    "FunctionSpace",
    "TensorProductSpace",
    "angular_flux_space",
    "scalar_flux_space",
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


@dataclass(frozen=True, eq=False)
class FunctionSpace(Generic[Carrier]):
    r"""A finite-dimensional vector space of discrete fields.

    Parameters
    ----------
    name : str
        Human-readable identifier. Used by :meth:`__repr__` and by
        :class:`IncompatibleOperatorComposition` error messages. Two
        spaces with the same ``name`` and ``shape`` compare equal even
        if they were constructed via different factory functions —
        ``name`` is the **identity** of the space, not a description
        of its contents.
    shape : tuple[int, ...]
        Tensor shape of the elements. A ``shape=(n_cells,
        n_ordinates, n_groups)`` space contains 3-D arrays; a
        ``shape=(n_cells,)`` space contains 1-D arrays.
    inner_product_weights : NDArray | None, default None
        Diagonal weights for the L² inner product
        :math:`\langle x, y \rangle = \sum_i w_i \, x_i \, y_i`. The
        weights array MUST be broadcast-compatible with ``shape``;
        most commonly a 1-D array along one axis (e.g. quadrature
        weights along the ordinate axis of an angular-flux space).
        ``None`` selects the Euclidean inner product
        :math:`\sum_i x_i \, y_i` — unless a ``metric`` object is
        installed (below), which this slot then does NOT describe: a
        dense-metric space reads ``inner_product_weights is None``
        while carrying a real, non-Euclidean metric.
    metric : HilbertMetric | None, keyword-only, default None
        The metric **object** (campaign 1, P7) — a
        :class:`~orpheus.numerics.metric.HilbertMetric` realization for
        forms no diagonal array can spell (a
        :class:`~orpheus.numerics.metric.DenseMetric` Gram inverse being
        the founding occupant). On an AXIS-BUILT space it is the
        positioned OVERLAY of forms (CS4c step 6 item 6.2c-i, ruling
        R-6.2c-2): a :class:`~orpheus.numerics.metric.FactoredMetric` with
        one entry per axis, a form only on a block whose axis carries NO
        measure (a dense-Gram head), ``None`` elsewhere — the space's
        metric is DERIVED from its axes with that overlay merged in, so
        every axis keeps supplying its own diagonal block. Resolution
        order: the axes-derived object (overlay merged) > an explicit
        object on an axes-less space > ``inner_product_weights`` >
        Euclidean. Validated at construction via
        :meth:`HilbertMetric.validate_for` and, beside axes, by the
        positioning guard.
    axes : tuple[Axis, ...] | None, default None
        The generator record of an axis-composed space (campaign 1,
        CS1): the ordered tensor factors this space is the product of.
        ``None`` means *legacy / not axis-built* — every pre-CS1
        construction path. Populated by :meth:`of_axes` (and threaded
        through ``*``); when present, the factor measures live PER AXIS
        (``inner_product_weights`` stays ``None`` — never both, enforced
        at construction) and the metric is DERIVED from them as one
        :class:`~orpheus.numerics.metric.FactoredMetric` (item 6.2c-i).
        Since the identity flip (CS4c step 6,
        2026-09-07) it IS the identity of an axis-built space —
        ``__eq__``/``__hash__`` read it directly; the dataclass field
        stays ``compare=False`` only so a subclass's generated ``__eq__``
        never reaches an ndarray (see the identity paragraph below).

    Notes
    -----
    The class is **frozen**. A :class:`FunctionSpace` encodes pure
    *geometry* — the discrete degrees of freedom (``shape``), the
    inner-product metric, and the composition algebra (``*`` /
    :meth:`dual`). It is **role- and dimension-agnostic** (the
    "View-G" decision, issues #205 / #207): a flux ``ψ`` and a
    reaction-rate density ``Lψ`` live on the *same* geometric space
    even though they carry different units. **Units do NOT live on the
    space** — they are a property of the *quantity*, carried by the
    :class:`~orpheus.numerics.field.Field` role-leaf (as a class
    constant ``UNITS``) and, for maps, by the operator's unit-gain
    (issue #208). This keeps one space per grid (no ``flux_space`` vs
    ``ratedensity_space`` duplication) and lets ``L`` / ``L⁻¹`` type as
    geometric endomorphisms on the bulk grid with a dimensional gain.

    **Identity — the chartered doctrine, realized (campaign 1, item 0.8;
    the identity flip landed at CS4c step 6, 2026-09-07).** The doctrine
    (ruled 2026-08-19/20): a space's identity is the **structural content
    of its axes** plus its tags — and since the metric lives on the axes,
    *metric differences imply space differences* (a quotient point with
    unit weight and a genuine one-cell mesh with :math:`V \neq 1` are
    DIFFERENT spaces; the old reading "two copies of :math:`\mathbb{R}^n`
    are 'the same' space regardless of which inner product is installed"
    is OVERTURNED). The realization is ONE key, :meth:`identity_key`, from
    which ``==`` and ``hash`` derive in this class only (S7 of the numerics
    untangling, #631; until then twelve subclasses carried their own pair):

    * an **axis-built** space is its ``axes`` tuple (``Axis.__eq__`` reads
      the axis content, support and atoms included, and never the origin;
      ruling 71) together with the dense forms POSITIONED over it, so an
      all-axes ``*`` product and the :meth:`of_axes` mint of the same axes
      are one space, a dense-dressed head and the same axes under another
      Gram are two, and an axis-built space is never equal to a hand-named
      one wearing its label;
    * an **axes-less leaf** is its class, name and shape: the trace spaces
      and the radial-characteristic sub-spaces fold their content into a
      digest name (CS4b S3), so that triple IS content identity there, and
      the class keeps a hand-built space wearing such a name from passing
      for one; a hand-built legacy space carries whatever its author wrote;
    * a **dual**, a non-axis **tensor product**, a **direct sum**
      (:class:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace`) and
      a **coupled system** are their members' keys, tagged, so covariance
      travels into a composite (F7 of the S6a review: until S7 a composite
      folded its members' NAMES, and ``from_blocks(a, b).dual()`` equalled
      ``from_blocks(a*, b*)``).

    Until the identity flip (CS4c step 6) the doctrine flowed through a
    BRIDGE: identity was ``(name, shape)`` for every space and
    :meth:`of_axes` derived the name injectively from the axes' content. The
    derived name survives as the readable label and as what keeps ``repr``
    and guard messages content-distinguishing.
    """

    name: str
    shape: tuple[int, ...]
    # Metadata, NOT identity (see class docstring): on an axes-less leaf the
    # key is (class, name, shape), so two such spaces are equal regardless of
    # the installed dense weights (the leaf classes derive their names from
    # content); on an axis-built space the measures live ON the axes, so they
    # enter the identity through ``axes`` and this slot is ``None``.
    # ``compare=False`` keeps an array out of the content of any dataclass
    # that would read the fields.
    inner_product_weights: Optional[NDArray] = field(
        default=None, repr=False, compare=False,
    )
    # The axis-composition record (campaign 1, CS1) — the identity of an
    # axis-built space (with the positioned overlay; :meth:`identity_key`)
    # since the identity flip (CS4c step 6) — and
    # LOAD-BEARING for the metric: when present, the factor measures live
    # per axis and the metric is DERIVED from them as one FactoredMetric
    # (never densified; item 6.2c-i retired the inline per-axis loop).
    axes: Optional[tuple[Axis, ...]] = field(
        default=None, repr=False, compare=False,
    )
    # The metric OBJECT (campaign 1, P7): the third metric source — a
    # first-class HilbertMetric realization (today: DenseMetric, the
    # non-Hadamard case no array-or-axes source can spell). Resolution
    # order in _resolved_metric: axes-derived (with this object merged in
    # as the positioned OVERLAY of forms) > an object on an axes-less
    # space > dense weights > None; beside axes the object must be
    # POSITIONED over them (one entry per axis, a form only where the
    # axis carries no measure) — the guard checks it.
    # ``compare=False`` is structurally MANDATORY, not taste: [M] a
    # compared metric field makes the dataclass-generated ``__eq__`` of
    # subclasses return an ndarray and ``hash()`` raise — the same
    # mechanism recorded for the weights above, re-measured for P7.
    # ``kw_only`` so subclass positional fields keep their order.
    metric: Optional[HilbertMetric] = field(
        default=None, repr=False, compare=False, kw_only=True,
    )

    def __post_init__(self) -> None:
        if self.metric is not None:
            # The object knows its own admission (a DenseMetric must span
            # the flattened leading block; a DiagonalMetric is as
            # permissive as the legacy array slot).
            self.metric.validate_for(self.shape)
        # ONE metric source, enforced pairwise over all three (P7 S2 —
        # illegal states unrepresentable; these are construction bugs,
        # not user input). ⚠ Until P7 the check lived under the
        # axes-early-return below, so the (dense, metric) arm was
        # structurally unreachable; each arm now has its own witness
        # (battery arms M10a/b/c).
        if self.axes is not None and self.inner_product_weights is not None:
            raise ValueError(
                f"space {self.name!r} carries BOTH per-axis measures and "
                f"dense inner_product_weights — one metric source only "
                f"(the axes own the measure on an axis-built space)"
            )
        if self.axes is not None and self.metric is not None:
            # CS4c step 6 item 6.2c-i: an axis-built space's metric is a
            # DERIVED object over its axes; an explicit object is admitted
            # only as the positioned OVERLAY of forms (one entry per axis,
            # in order, a form only on a measure-less axis's block) — the
            # home of a dense Gram on an axis-built head (ruling
            # R-6.2c-2). Anything else is two metric sources.
            _require_positioned_over(self.metric, self.axes, name=self.name)
        if self.inner_product_weights is not None and self.metric is not None:
            raise ValueError(
                f"space {self.name!r} carries BOTH dense "
                f"inner_product_weights and a metric object — one metric "
                f"source only (install the object OR the array, never both)"
            )
        if self.axes is None:
            return
        # The remaining axis-built guard: a shape that IS the axes'
        # concatenation.
        concat: tuple[int, ...] = ()
        for ax in self.axes:
            concat = concat + ax.shape
        if concat != self.shape:
            raise ValueError(
                f"space {self.name!r} shape {self.shape} does not equal "
                f"the concatenation {concat} of its axes' shapes"
            )

    def identity_key(self) -> tuple[Any, ...]:
        r"""The content that IS this space: the one definition of its identity.

        ``==`` and ``hash`` derive from it here, in the base class only, so no
        member of the family re-states equality (#631, finding F7 of the S6a
        review). A tagged tuple per kind of space:

        * axis-built — ``("axes", axes, overlay)``: the axis tuple (its
          content, ruling 71) and the content digests of the dense forms
          POSITIONED over it (``None`` where an axis carries none), so a
          dense-dressed head and the same axes with another Gram are two
          spaces;
        * axes-less leaf — ``("named", class, name, shape)``: the
          factory-built leaves (the two trace spaces, the two
          radial-characteristic spaces) fold their content into a digest
          name, and the class keeps a hand-built space wearing such a name
          from passing for one;
        * a dual, a non-axis tensor product, a direct sum and a coupled
          system override this with their members' keys, tagged.
        """
        if self.axes is not None:
            return ("axes", self.axes, tuple(
                None if form is None else content_digest(form) for form in self._overlay_forms()
            ))
        return ("named", type(self), self.name, self.shape)

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
    # Axis composition (campaign 1, CS1)
    # ------------------------------------------------------------------

    @classmethod
    def of_measure(cls, measure: DiscreteMeasure, label: str) -> "FunctionSpace":
        r"""The discrete :math:`L^2(\mu)` space a measure samples: one axis.

        A discrete measure :math:`\mu = \sum_n w_n \delta_{x_n}` induces the
        space of its sampled values, per-node arrays under
        :math:`\langle f, g\rangle_\mu = \sum_n w_n f_n g_n` (the weights ARE
        the metric). It is the space over :meth:`Axis.of_measure
        <orpheus.numerics.axis.Axis.of_measure>`, so its identity is the
        axis's content (support, atoms, weights, label; ruling 71) and the
        measure space of a frame over a quadrature IS the SN angular space
        over the same rule (ruling 63). The label is the caller's: the role
        the factor plays, which a caller holding only the measure reads off
        :attr:`~orpheus.numerics.measure.DiscreteMeasure.axis_label`.
        """
        return cls.of_axes(Axis.of_measure(measure, label))

    @classmethod
    def of_basis(cls, basis: Basis | TruncatedBasis) -> "FunctionSpace":
        r"""The coefficient space a basis spans: one axis, indexed by its layout.

        The space over :meth:`Axis.of_basis <orpheus.numerics.axis.Axis.of_basis>`
        (ruling 81): the layout, the continuum Gram and the domain are the
        basis's, read with no space import on the basis side. Its identity is
        that axis's content, so two families on two domains (the harmonics on
        :math:`S^2` and their σ-even restriction on :math:`S^2/\sigma_y`) and
        two partitions with different edges span two spaces (#631 pairs 3
        and 4). Replaced each basis's own ``space`` property at S7.
        """
        return cls.of_axes(Axis.of_basis(basis))

    @classmethod
    def of_axes(cls, *axes: Axis) -> "FunctionSpace":
        r"""Compose a space as the ordered product of its axes.

        The ONE composition mechanism of the axis doctrine: a space IS
        its axis tuple. Shape is the concatenation of the axes' shapes;
        the factor measures stay PER AXIS (``inner_product_weights`` is
        never populated — **no densification**, structurally: composing
        two 2000-point weighted axes stores 2 × 2000 weights, never the
        4 000 000-entry outer product); the metric machinery routes
        through the per-axis path.

        **The derived name — the readable label, and the carrier the
        axes-less composites fold.** The name is derived
        DETERMINISTICALLY and INJECTIVELY from the axes' structural
        content (label, shape, kind, measure bytes, subclass identity —
        via a content digest, never ``hash()``, so it is stable across
        processes). Until the identity flip (CS4c step 6, 2026-09-07)
        that injectivity was the identity BRIDGE: space identity was
        ``(name, shape)`` for every space, so different axis tuples had
        to mint different names to be different spaces. Since the flip
        an axis-built space compares by its ``axes`` tuple directly
        (:meth:`__eq__`), and the digest's injectivity is load-bearing
        only where a derived name is folded into an axes-less composite's
        own digest
        (:meth:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace.from_blocks`)
        and for what ``repr`` and the guard messages can tell apart. Two
        same-``ng`` energy axes with different partitions, or two
        same-shape spatial axes with different measures, compose into
        UNEQUAL spaces either way — the chartered "metric differences
        imply space differences" doctrine, now directly.

        Always returns a plain :class:`FunctionSpace` — an axis product
        is not a different *kind* of space (ruled Q-T4); invoking
        through a subclass does not change the return type.

        **The second composition, ``V * W``** (the ``*`` dunder →
        :class:`TensorProductSpace`), threads ``axes`` when every factor
        is axis-built — then it IS this composition, with ``factors``
        metadata on top — and otherwise positions the factor measures in
        a lazy :class:`~orpheus.numerics.metric.FactoredMetric`, per
        axis for an axis-built factor. ⚠ Until CS4c step 6 item 6.2a
        (2026-09-07) that second arm DENSIFIED the metric into an
        outer-product ``inner_product_weights`` tensor, bridging
        axis-borne measures through a densifier so no value was lost —
        `[M]` 18 state-sized tensors per 2-D windowed solve. The
        densifier retired with the dense arm; what keeps ``*`` alive is
        the axes-less angular head of the moment products, which item
        6.2c makes axis-built. New axis-aware code composes with
        ``of_axes``.
        """
        if not axes:
            raise ValueError(
                "of_axes needs at least one axis — a space with no "
                "factors has no index set"
            )
        shape: tuple[int, ...] = ()
        for ax in axes:
            shape = shape + ax.shape
        digest = name_digest(tuple(axes))
        readable = "*".join(f"{ax.label}{ax.shape}" for ax in axes)
        return FunctionSpace(
            name=f"{readable}#{digest}",
            shape=shape,
            inner_product_weights=None,
            axes=tuple(axes),
        )

    @property
    def has_coordinate_cone(self) -> bool | None:
        r"""Whether per-component positivity is meaningful on this space.

        Three-valued, deliberately:

        * ``True`` — axis-built, ALL factors ``NODAL``: components are
          point/cell values, the coordinate cone :math:`K = \{x \ge 0\}`
          is the physical positive cone, and per-component sign tests
          (``Field.cone_violations``) are meaningful.
        * ``False`` — axis-built, ANY factor ``MODAL``: components are
          expansion coefficients; a positive function may have negative
          coefficients, so a per-component sign test is MEANINGLESS and
          must be refused, not answered.
        * ``None`` — ``axes is None`` (legacy / not migrated): the
          question cannot be answered structurally; consumers keep their
          pre-CS1 behavior. Collapsing ``None`` into ``False`` would
          make the cone refusal fire on every legacy space in the tree.
        """
        if self.axes is None:
            return None
        return all(ax.kind is BasisKind.NODAL for ax in self.axes)

    # ------------------------------------------------------------------
    # Axis lookup — the public by-label accessor (un-weld arc S-1)
    # ------------------------------------------------------------------

    def _axis_index(self, label: str) -> int:
        """Position of the unique axis labeled ``label`` in :attr:`axes`.

        The one home of by-label axis resolution: the public
        :meth:`axis` accessor and the axis-marginal mint
        (:mod:`orpheus.numerics.axis_marginal`) both route through here,
        so the refusal vocabulary cannot drift between them.
        """
        if self.axes is None:
            raise TypeError(
                f"axis lookup {label!r}: {self!r} is not axis-built "
                f"(axes is None) — no named factors to look up. Compose "
                f"the space with FunctionSpace.of_axes."
            )
        hits = [i for i, ax in enumerate(self.axes) if ax.label == label]
        if len(hits) != 1:
            raise ValueError(
                f"axis lookup: label {label!r} names {len(hits)} axes "
                f"of {self!r} (have {[ax.label for ax in self.axes]}) — "
                f"need exactly one."
            )
        return hits[0]

    def axis(self, label: str) -> Axis:
        r"""The unique axis labeled ``label`` — by-label factor access.

        The public spelling of reads that used to route through carrier
        shape metadata: the space carries its axes, so the axis IS the
        metadata's home — ``space.axis("energy").shape[0]`` is the group
        count, ``space.axis("spatial").shape`` the spatial shape.
        Refuses a legacy name-built space (``axes is None``, TypeError)
        and a label naming zero or several axes (ValueError, naming the
        inventory) — the same vocabulary the axis marginal refuses in
        (shared resolver).
        """
        k = self._axis_index(label)
        axes = self.axes
        assert axes is not None  # narrowing only: _axis_index refused the None case
        return axes[k]

    def _overlay_forms(self) -> tuple[DiagonalMetric | DenseMetric | None, ...]:
        r"""The positioned overlay of forms, one entry per axis (``None`` where
        the axis carries no form, and on every axis of a space with no overlay).

        The one reader of the overlay for the product composition
        (:func:`_positioned_forms`) and the axis verbs :meth:`without_axis` and
        :meth:`with_axis`; :func:`_axis_entries` reads the factored metric's
        ``entries`` directly instead.
        """
        assert self.axes is not None  # narrowing only: callers are axis-built
        if isinstance(self.metric, FactoredMetric):
            return _overlay_leaf_forms(self.metric)
        return (None,) * len(self.axes)

    def _axis_block(self, k: int) -> DiagonalMetric | DenseMetric | None:
        r"""The metric block :math:`G_k` of axis ``k``: the positioned form on
        its block if there is one, else its measure, else ``None`` (the
        counting measure). The same entry the space's derived metric
        composes (:func:`_axis_entries`), so the axis marginal's covector
        :math:`m = G_k u` and the space's metric cannot disagree.
        """
        assert self.axes is not None  # narrowing only: callers are axis-built
        return _axis_entries(self.axes, self.metric)[k][1]

    def _with_overlay(
        self, forms: tuple[DiagonalMetric | DenseMetric | None, ...],
    ) -> "FunctionSpace":
        """This axis-built space with ``forms`` positioned over its axes (none: no object)."""
        assert self.axes is not None  # narrowing only: callers are axis-built
        if all(form is None for form in forms):
            return self
        return replace(self, metric=FactoredMetric(tuple(
            (ax.shape, form) for ax, form in zip(self.axes, forms)
        )))

    def without_axis(self, label: str) -> "FunctionSpace":
        r"""The product of the OTHER axes, each keeping its positioned form.

        The marginal space of a fibre integration over ``label``: its metric
        is this space's metric with the axis's block removed, so the full
        metric is exactly (the axis measure) ⊗ (this marginal's metric)
        whenever the dropped axis carries no form itself.
        """
        k = self._axis_index(label)
        assert self.axes is not None  # narrowing only: _axis_index refused the None case
        forms = self._overlay_forms()
        others = tuple(ax for i, ax in enumerate(self.axes) if i != k)
        return FunctionSpace.of_axes(*others)._with_overlay(
            tuple(form for i, form in enumerate(forms) if i != k)
        )

    def with_axis(self, axis: Axis) -> "FunctionSpace":
        r"""This product with ``axis`` appended last, every positioned form kept."""
        if self.axes is None:
            raise TypeError(f"with_axis: {self!r} is not axis-built (axes is None)")
        return FunctionSpace.of_axes(*self.axes, axis)._with_overlay(self._overlay_forms() + (None,))

    # ------------------------------------------------------------------
    # Inner product / norm
    # ------------------------------------------------------------------

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
            Elements of this space (arrays of shape :attr:`shape` for
            the default realization).

        Returns
        -------
        float
            Scalar inner product.
        """
        return self._diagonal_inner_product(x, y)

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
        """
        return float(np.sqrt(self.inner_product(x, x)))

    # ------------------------------------------------------------------
    # Metric application (Hilbert adjoint building blocks, Wave O / O.2b)
    # ------------------------------------------------------------------

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
        return self._diagonal_apply_metric(x)

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
        return self._diagonal_apply_inverse_metric(x)

    def _diagonal_apply_inverse_metric(self, x: Any) -> Any:
        r"""The bare-array realization of :meth:`apply_inverse_metric`."""
        m = self._resolved_metric
        if m is None:
            return x
        return m.apply_inverse(x)

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
            return int(np.prod(self.shape)) if self.shape else 1
        return metric.rank(self.shape)

    def project_live(self, x: Carrier) -> Carrier:
        r"""The G-orthogonal projector onto the live part, :math:`P = G^{+}G`.

        :math:`P` is the orthogonal projector onto :math:`\operatorname{range} G
        = (\ker G)^{\perp}`; it is the identity on a non-degenerate space (to
        rounding: `[M]` at most 2.7e-16 relative, the two metric faces applied
        in turn) and zeroes the kernel coordinates of a degenerate one. Two elements are
        the same element of the quotient :math:`V / \ker G` iff their live
        projections agree, which is how an adjoint law that holds modulo the
        kernel is stated and gated (ruling 74; e.g. ``A.H.H`` and ``A`` agree
        after the codomain's projection). One definition for every space: the
        composites' block-wise metric verbs make it block-wise.
        """
        return self.apply_inverse_metric(self.apply_metric(x))

    @property
    def _resolved_metric(self) -> Optional[HilbertMetric]:
        r"""The metric SOURCE resolved to its realization, PER CALL.

        Resolution order (CS4c step 6 item 6.2c-i, ruling R-6.2c-2): an
        axis-built space's metric is DERIVED from its axes as a
        :class:`~orpheus.numerics.metric.FactoredMetric` — one
        :class:`~orpheus.numerics.metric.DiagonalMetric` per weighted axis
        on that axis's block, ``None`` on a counting-measure axis — with
        the space's positioned OVERLAY of forms merged in (a dense Gram's
        :class:`~orpheus.numerics.metric.DenseMetric` on the block of an
        axis that carries no measure; the guard admitted it), no object at
        all when every block is Euclidean (:func:`_axes_metric`); else an
        explicit ``metric`` object on an axes-less space; else a legacy
        ``inner_product_weights`` array resolves to a ``DiagonalMetric``;
        ``None`` means Euclidean and the verbs short-circuit without an
        object. ONE realization for every source:
        :meth:`DiagonalMetric.apply_block` is operation-for-operation the
        per-axis reshape-and-multiply the retired ``_apply_axes_weights``
        performed inline (`[M]` bit-identical on every axis-built space —
        gated in ``tests/gates/numerics/test_axis_metric_is_a_derived_object.py``),
        and the legacy diagonal arm is operation-for-operation the arms
        that used to live inline in the three ``_diagonal_*`` realizations.

        ⛔ Deliberately NOT a ``cached_property``, and the reason is a
        measured red: the mutation-battery idiom mutates a frozen
        space's weight FIELD in place (``object.__setattr__``), and the
        adjoint-certification propagation probes assert the mutation
        reaches the metric surface at the next read. `[M]` 2026-08-30
        (P7 exit gate): a cached resolution served the stale
        ``DiagonalMetric`` and
        ``test_gsd_metric_drop_is_k_blind_but_vector_red``'s own
        ``after == before`` probe reddened — the pre-P7 semantics read
        the field per call, and every such battery relies on it. The
        wrapper is a tiny frozen object; the numpy work dominates.
        """
        if self.axes is not None:
            return _axes_metric(self.axes, self.metric)
        if self.metric is not None:
            return self.metric
        if self.inner_product_weights is not None:
            return DiagonalMetric(self.inner_product_weights)
        return None

    # ------------------------------------------------------------------
    # Space algebra (Depth B step D-B)
    # ------------------------------------------------------------------

    def __mul__(self, other: "FunctionSpace") -> "TensorProductSpace":
        r"""Return the tensor product :math:`V \otimes W` of this space
        with ``other``.

        Implements ``V * W`` per grand-report §6.1: the resulting
        :class:`TensorProductSpace` carries the concatenated shape
        ``self.shape + other.shape``, the factor measures positioned per
        block (concatenated axes when every factor is axis-built, a lazy
        :class:`~orpheus.numerics.metric.FactoredMetric` otherwise —
        never a stored outer product, since CS4c step 6 item 6.2a), and
        the multiplied units. Associative on its inputs:
        ``(A * B) * C`` and ``A * (B * C)`` both produce a flat
        3-factor :class:`TensorProductSpace`.

        Loadbearing for the Wave T tensor-network rewires per the
        grand report §15.1 (streaming as
        :math:`L = \sum_{\text{axis}} D_{\text{axis}} \otimes \Omega_{\text{axis}} \otimes I_g`)
        and §16A.10 (boundary as
        :math:`B = G_{\text{patch}} \otimes K_\omega \otimes K_g`).
        See ``.claude/plans/depth_b_field_on_function_space.md`` §6
        step D-B for the design.
        """
        if not isinstance(other, FunctionSpace):
            return NotImplemented
        self_factors = (
            self.factors if isinstance(self, TensorProductSpace) else (self,)
        )
        other_factors = (
            other.factors if isinstance(other, TensorProductSpace) else (other,)
        )
        return TensorProductSpace.from_factors(self_factors + other_factors)

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


# ---------------------------------------------------------------------------
# Tensor-product and dual space (Depth B step D-B)
# ---------------------------------------------------------------------------


def _overlay_leaf_forms(
    overlay: FactoredMetric,
) -> tuple[DiagonalMetric | DenseMetric | None, ...]:
    r"""The forms of an axis-built space's positioned OVERLAY: leaf metrics
    only (the overlay guard admits a diagonal or dense form per axis, never a
    dual's inverse form, which lives on a :class:`DualSpace` that has no
    axes)."""
    forms: list[DiagonalMetric | DenseMetric | None] = []
    for _, form in overlay.entries:
        if form is None or isinstance(form, (DiagonalMetric, DenseMetric)):
            forms.append(form)
        else:
            raise TypeError(
                f"an axis-built space's overlay holds diagonal or dense forms "
                f"only, got {type(form).__name__}"
            )
    return tuple(forms)


def _axis_entries(
    axes: tuple[Axis, ...],
    positioned: Optional[HilbertMetric] = None,
) -> list[tuple[tuple[int, ...], DiagonalMetric | DenseMetric | None]]:
    r"""One positioned metric entry per AXIS — the block form of the axis
    doctrine's "the measure lives per axis" — with the space's positioned
    OVERLAY of forms merged in.

    A weighted axis contributes a
    :class:`~orpheus.numerics.metric.DiagonalMetric` carrying its measure on
    its own index block; a counting-measure axis (``weights is None``, the
    canonical spelling) contributes ``None`` — no pass over that block —
    unless ``positioned`` (the space's admitted overlay: one entry per
    axis, a form only where the axis carries no measure) places a FORM on
    that block, which then IS the entry: the dense Gram no diagonal
    measure can spell (ruling R-6.2c-2, 2026-09-08 — *each axis supplies
    its diagonal block from its measure; a dense-Gram head supplies a
    positioned DenseMetric on its block*). The guard makes shadowing
    unspellable, so the merge never has to choose.
    :meth:`DiagonalMetric.apply_block` performs exactly the explicit
    reshape-and-multiply (leading 1s for the preceding axes' ranks, the
    axis shape, trailing 1s for the rest — exact for interior axes, never
    materializing the outer product; the ERR-067-family divergence is
    unspellable here) that ``FunctionSpace._apply_axes_weights`` performed
    inline until CS4c step 6 item 6.2c-i retired it into this one
    realization.
    """
    forms: tuple[DiagonalMetric | DenseMetric | None, ...] = (
        (None,) * len(axes)
        if not isinstance(positioned, FactoredMetric)
        else _overlay_leaf_forms(positioned)
    )
    return [
        (
            ax.shape,
            form if form is not None
            else (None if ax.weights is None else DiagonalMetric(ax.weights)),
        )
        for ax, form in zip(axes, forms)
    ]


def _axes_metric(
    axes: tuple[Axis, ...], positioned: Optional[HilbertMetric] = None,
) -> Optional[FactoredMetric]:
    r"""The metric an axis-built space DERIVES from its axes (CS4c step 6
    item 6.2c-i): a lazy :class:`~orpheus.numerics.metric.FactoredMetric`
    over :func:`_axis_entries` (the positioned overlay merged in), or
    ``None`` when every block is Euclidean (the default — no object, no
    allocation).
    Rebuilt per resolution (a handful of frozen wrappers around the axes'
    own arrays; nothing is copied), so an in-place weight mutation reaches
    the metric surface at the next read exactly as the dense-slot arm's
    does.
    """
    entries = _axis_entries(axes, positioned)
    if all(entry is None for _, entry in entries):
        return None
    return FactoredMetric(tuple(entries))


def _require_positioned_over(
    metric: HilbertMetric, axes: tuple[Axis, ...], *, name: str,
) -> None:
    r"""Admit an explicit metric object on an AXIS-BUILT space only as the
    positioned OVERLAY of forms — a
    :class:`~orpheus.numerics.metric.FactoredMetric` with exactly one entry
    per axis, in order, each on that axis's block, a form ONLY on a block
    whose axis carries no measure, and never a diagonal one.

    That is the one shape in which a form no diagonal measure can spell —
    a dense Gram's :class:`~orpheus.numerics.metric.DenseMetric` on a
    moment head's block — lives beside per-axis measures without becoming
    a second metric source (ruling R-6.2c-2, 2026-09-08: *a Gram is a FORM,
    never on* ``Axis.weights``; the space's metric is the derived object
    and a dense block is positioned ON it). Everything else is refused by
    name, each a second spelling of one measure: a diagonal array spanning
    the whole shape or a factored object whose blocks do not follow the
    axes (the derived object would silently disagree with it — `[M]` a
    relaxed guard alone left the per-axis short-circuit reading the axes
    and ignoring the object by ``max|Δ| = 14.6`` on a dense-Gram head); a
    positioned :class:`~orpheus.numerics.metric.DiagonalMetric` (a diagonal
    measure is the axis's own to carry); a form on an axis that ALSO
    carries a measure (two sources on one block — the overlay merge would
    have to pick one and the loser would be dead data).
    """
    if not isinstance(metric, FactoredMetric):
        raise ValueError(
            f"space {name!r} carries per-axis measures and a "
            f"{type(metric).__name__} — an axis-built space's metric is "
            f"DERIVED from its axes; an explicit object must be a "
            f"FactoredMetric positioned over them (one entry per axis)."
        )
    blocks = tuple(block for block, _ in metric.entries)
    axis_blocks = tuple(ax.shape for ax in axes)
    if blocks != axis_blocks:
        raise ValueError(
            f"space {name!r}: the metric object's blocks {blocks} do not "
            f"follow the axes' blocks {axis_blocks} — a positioned object "
            f"carries exactly one entry per axis, in order."
        )
    for k, ((_, form), ax) in enumerate(zip(metric.entries, axes)):
        if form is None:
            continue
        if isinstance(form, DiagonalMetric):
            raise ValueError(
                f"space {name!r}: entry {k} (axis {ax.label!r}) positions a "
                f"DiagonalMetric — a diagonal measure is the axis's own to "
                f"carry (Axis.weights); only a FORM no measure can spell (a "
                f"dense Gram) is positioned on the derived object."
            )
        if ax.weights is not None:
            raise ValueError(
                f"space {name!r}: axis {ax.label!r} carries a measure AND "
                f"entry {k} positions a {type(form).__name__} on its block — "
                f"two metric sources on one block; a form REPLACES the "
                f"measure (mint that axis with weights=None)."
            )


def _positioned_forms(
    factors: tuple["FunctionSpace", ...],
) -> Optional[FactoredMetric]:
    r"""The positioned OVERLAY an all-axis-built product carries beside its
    concatenated axes: each factor's admitted overlay (one entry per axis)
    in order, ``None`` on every axis of a factor that carries none — or no
    object at all when no factor carries a form. The product then derives
    exactly the entries its factors derive (item 6.2c-i).
    """
    entries: list[
        tuple[tuple[int, ...], PositionedMetric | None]
    ] = []
    for f in factors:
        assert f.axes is not None  # narrowing only: the caller's arm
        entries.extend((ax.shape, form) for ax, form in zip(f.axes, f._overlay_forms()))
    if all(form is None for _, form in entries):
        return None
    return FactoredMetric(tuple(entries))


def _tensor_product_factored_metric(
    factors: tuple["FunctionSpace", ...],
) -> Optional[FactoredMetric]:
    r"""Assemble a tensor product's metric as a lazy per-BLOCK product —
    the one non-axis arm of ``*`` (CS4c step 6 item 6.2a, 2026-09-07).

    One positioned entry per block, and the Kronecker product is never
    materialized:

    * an **axis-built** factor contributes one entry per AXIS — a
      :class:`~orpheus.numerics.metric.DiagonalMetric` carrying the
      axis measure on that axis's block, ``None`` for a counting-measure
      axis — so the product applies exactly the entries the factor
      DERIVES for itself (:func:`_axis_entries`, CS4c step 6 item 6.2c-i;
      until then the factor's inline ``_apply_axes_weights`` loop spelled
      the same arithmetic — :meth:`~orpheus.numerics.metric.DiagonalMetric.apply_block`
      is operation-for-operation that reshape-and-multiply: bit-identical
      on an axis block by construction);
    * a factor carrying a metric **object** rides verbatim — a nested
      :class:`~orpheus.numerics.metric.FactoredMetric` flattens (a
      product of products is one product), a
      :class:`~orpheus.numerics.metric.DenseMetric` positions on the
      factor's block (the Parseval-dressed head of a dense-Gram frame);
    * a **dense-slot leaf** (``inner_product_weights`` set, no axes — the
      spherical-harmonic / Legendre heads, the trace spaces, hand-built
      legacy spaces) becomes a ``DiagonalMetric`` on its block;
    * a **Euclidean** factor contributes ``None`` (no pass over its block).

    Returns ``None`` when every block is Euclidean, preserving the
    Euclidean default (no object, no allocation).

    **History.** Until 6.2a ``*`` had a DENSE arm that formed
    :math:`w_1 \otimes w_2 \otimes \cdots` as one weights tensor on every
    mint with an axes-less factor — `[M]` 18 ``(L+1, 2L+1, ng, nx, ny)``
    tensors per 2-D windowed solve at ``max_inner = 6``,
    ``2·max_inner + 6`` in general — and this builder (the P7
    dense-factor arm) was unreachable on every SN path (no shipped factor
    carried a metric object) *and* itself densified an axis-built factor's
    measure. The factored application agrees with the retired outer
    product to `[M]` 2 ULP worst over 1600 draws (200 seeds × 8
    (geometry × L) rows — one extra rounding, reduction depth +1; the
    ``vv-principles`` bit-identity criterion 3), gated at ``nulp = 4`` in
    ``tests/gates/numerics/test_tensor_product_metric_is_factored.py``.
    """
    entries: list[
        tuple[tuple[int, ...], "PositionedMetric | None"]
    ] = []
    for f in factors:
        # ONE metric resolution: the factor's own (an axis-built factor's
        # per-axis entries with its overlay; a metric object; a dense-slot
        # leaf's DiagonalMetric; a dual space's inverse form), never its raw
        # fields read a second time here (Pattern 7: until 2026-10-21 this
        # loop re-read ``axes``/``metric``/``inner_product_weights`` and so
        # read a dual factor, which carries none of them, as Euclidean).
        m = f._resolved_metric
        if m is None:
            entries.append((f.shape, None))
        elif isinstance(m, FactoredMetric):
            # A product of products is one product.
            entries.extend(m.entries)
        elif isinstance(m, (DiagonalMetric, DenseMetric, DualMetric)):
            entries.append((f.shape, m))
        else:
            # Type-narrowing refusal: the positioned application is
            # defined on the leaf realizations; anything else would
            # die one call later with a worse message.
            raise TypeError(
                f"a tensor product can position only diagonal/dense "
                f"factor metrics or their duals, got {type(m).__name__}"
            )
    if all(entry is None for _, entry in entries):
        return None
    return FactoredMetric(tuple(entries))


@dataclass(frozen=True, eq=False)
class TensorProductSpace(FunctionSpace):
    r"""A function space that decomposes as
    :math:`V = V_1 \otimes V_2 \otimes \cdots \otimes V_k`.

    The tensor-product structure makes algebraic identities of operators
    on this space (adjoint distributivity, composition distributivity,
    representation polymorphism) checkable at the type level. See grand-
    report §5.3, §15, §32.4 for the L1 motivation; see
    ``.claude/plans/wave_t_tensor_network.md`` for the production
    consumers being wired in Wave T.

    Construction
    ------------
    Two equivalent paths:

    * **Operator-algebra dispatch** — ``A * B`` where ``A`` and ``B``
      are :class:`FunctionSpace` instances returns a
      :class:`TensorProductSpace`. The dunder is associative on its
      inputs (``(A * B) * C`` flattens to a 3-factor product, never
      nests).
    * **Explicit factory** — :meth:`from_factors` with a tuple of
      factor spaces.

    Notes
    -----
    The class is **frozen**. Identity is :meth:`identity_key`: an
    all-axes product is its axes (one space with the :meth:`of_axes` mint of
    them); otherwise it is the ordered keys of its ``factors``. Two
    :class:`TensorProductSpace` instances with the same factor sequence
    compare equal even if reached via different composition paths.

    The ``factors`` field also supports introspection (e.g.,
    operator-algebra factor matching for ``(A & B) ∘ (C & D)`` rewriting).

    Parameters
    ----------
    factors : tuple[FunctionSpace, ...]
        The factor spaces, in order. Should have ``len >= 2`` for a
        meaningful tensor product; trivial 1-factor TensorProductSpaces
        are permitted by the dataclass but produced only via the
        :meth:`__mul__` flattening edge cases.
    """

    # Required (a non-axis product IS its factors' keys; S7): kw_only for the field order.
    factors: tuple["FunctionSpace", ...] = field(kw_only=True, compare=False, repr=False)

    def identity_key(self) -> tuple[Any, ...]:
        """An axis-built product is its axes (the ``of_axes`` mint of them is
        one space); otherwise the ordered keys of its factors."""
        if self.axes is not None:
            return super().identity_key()
        return ("tensor", tuple(f.identity_key() for f in self.factors))

    @classmethod
    def from_factors(
        cls, factors: tuple["FunctionSpace", ...],
    ) -> "TensorProductSpace":
        r"""Construct a :class:`TensorProductSpace` from a tuple of factor
        spaces.

        Derives:
        * ``name`` from ``" ⊗ ".join(f.name for f in factors)``
        * ``shape`` from concatenated factor shapes
        * the metric: concatenated ``axes`` when every factor is
          axis-built, else a lazy
          :class:`~orpheus.numerics.metric.FactoredMetric` positioned
          per block (``None`` if every block is Euclidean) —
          ``inner_product_weights`` is never populated by a product
        """
        if len(factors) < 2:
            raise ValueError(
                f"TensorProductSpace.from_factors requires at least 2 "
                f"factors; got {len(factors)}"
            )
        name = " ⊗ ".join(f.name for f in factors)
        shape: tuple[int, ...] = ()
        for f in factors:
            shape = shape + f.shape
        # Two arms, and neither materializes a weights tensor (CS4c step
        # 6 item 6.2a retired the dense arm and its densifier bridge):
        # * axis threading (CS1, gate B7) — when EVERY factor carries an
        #   axes record, the product's record is the concatenation and
        #   the measure rides the per-axis path; a factor's positioned
        #   OVERLAY of forms (a dense-Gram head, 6.2c-i) rides beside it,
        #   concatenated per axis, so the product derives the same
        #   entries its factors do;
        # * the factored arm — an axes-less factor on either side leaves
        #   ``axes=None`` (never fabricate an axis for a space that did
        #   not declare one) and the product carries a lazy
        #   FactoredMetric positioned per block: per AXIS for an
        #   axis-built factor, per factor for a dense-slot leaf, verbatim
        #   for a metric object (P7 S2: dropping an object would silently
        #   revert the product to Euclidean on that factor, a VALUE bug —
        #   [M] 33.0 where G ⊗ w gives 109.0, on the harmonic-frame mint
        #   path).
        factor_axes = [f.axes for f in factors]
        metric: Optional[HilbertMetric] = None
        if all(fa is not None for fa in factor_axes):
            axes: Optional[tuple[Axis, ...]] = tuple(
                ax for fa in factor_axes if fa is not None for ax in fa
            )
            metric = _positioned_forms(factors)
        else:
            axes = None
            metric = _tensor_product_factored_metric(factors)
        return cls(
            name=name,
            shape=shape,
            inner_product_weights=None,
            axes=axes,
            metric=metric,
            factors=factors,
        )

    def _make_dual(self) -> "FunctionSpace":
        r"""The dual of a tensor product is the product of the duals,
        :math:`(V \otimes W)^* = V^* \otimes W^*`, each factor carrying its
        own inverse form (numerics_untangling.md rulings 26 and 75)."""
        return TensorProductSpace.from_factors(tuple(f.dual() for f in self.factors))

    def __repr__(self) -> str:
        return f"TensorProductSpace({self.name!r}, shape={self.shape})"


@dataclass(frozen=True, eq=False)
class DualSpace(FunctionSpace):
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
        if (
            self.axes is not None
            or self.metric is not None
            or self.inner_product_weights is not None
        ):
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


# ---------------------------------------------------------------------------
# Pre-populated common-space factories
# ---------------------------------------------------------------------------
#
# These are *factory functions*, not module-level singletons, because
# the shape of an angular-flux space depends on the mesh dimensions —
# every solver instance carries its own (n_cells, n_ordinates, n_groups)
# triple. Factories give the caller named, well-typed construction
# sites without forcing premature commitment to a global instance.


def angular_flux_space(
    n_cells: int,
    n_ordinates: int,
    n_groups: int,
    *,
    quadrature_weights: NDArray | None = None,
) -> FunctionSpace:
    r"""Construct the angular-flux space for an SN solve.

    Shape is ``(n_cells, n_ordinates, n_groups)``. When
    ``quadrature_weights`` is provided, it is broadcast along the
    ordinate axis as the inner-product metadata so the canonical
    angular inner product
    :math:`\langle \psi, \varphi \rangle_\Omega = \sum_n w_n \,
    \psi_n \, \varphi_n` (summed over cells / groups too) becomes
    :meth:`FunctionSpace.inner_product`.

    Parameters
    ----------
    n_cells, n_ordinates, n_groups : int
        Tensor dimensions.
    quadrature_weights : NDArray, optional
        Shape ``(n_ordinates,)`` quadrature weights along the
        ordinate axis. Reshaped to ``(1, n_ordinates, 1)`` so it
        broadcasts against the full angular-flux tensor.
    """
    weights: NDArray | None = None
    if quadrature_weights is not None:
        w = np.asarray(quadrature_weights)
        if w.shape != (n_ordinates,):
            raise ValueError(
                f"quadrature_weights must have shape ({n_ordinates},), "
                f"got {w.shape}"
            )
        weights = w.reshape(1, n_ordinates, 1)
    return FunctionSpace(
        name="angular_flux",
        shape=(n_cells, n_ordinates, n_groups),
        inner_product_weights=weights,
    )


def scalar_flux_space(n_cells: int, n_groups: int) -> FunctionSpace:
    r"""Construct the scalar-flux space.

    Shape is ``(n_cells, n_groups)``. No inner-product weights — the
    canonical inner product on scalar-flux space is the Euclidean
    sum (or, equivalently, the volume-weighted L² inner product when
    the cell volumes are absorbed into the operator that produced
    :math:`\phi`; see :meth:`Mesh1D.volume_measure` for the
    volume-weighted variant).
    """
    return FunctionSpace(
        name="scalar_flux",
        shape=(n_cells, n_groups),
    )


