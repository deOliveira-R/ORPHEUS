r"""The direct sum of function spaces: one type for the full field and the coupled system (concept ``dsum``).

A :class:`DirectSumSpace` is the biproduct :math:`V = \bigoplus_l V_l` of labelled summand spaces. Its coordinates are the summands'
coordinates concatenated in label order, and the layout that cuts the concatenation (:class:`~orpheus.numerics.coproduct_layout.CoproductLayout`)
is the one place the offsets are computed (law L13). The injection :math:`\iota_l` and the projection :math:`\pi_l` are :meth:`pack` and
:meth:`unpack` by label, and they satisfy the biproduct laws: :math:`\pi_l \iota_m = \delta_{lm}`, :math:`\sum_l \iota_l \pi_l = \mathrm{id}`, and
:math:`\langle \pi_l x, y\rangle_{V_l} = \langle x, \iota_l y\rangle_V` under the block-diagonal metric.

**Identity** is the ordered ``(label, summand key)`` pairs: two sums over isomorphic summands with different labels are unequal, and so are two
sums over the same labelled summands in two orders (the offsets differ). A nested sum is not equal to its flattening.

**The typed element.** A sum reads and rebuilds its typed element through three accesses of the carrier class and nothing else: ``blocks``, the tuple
of the summand elements in label order, ``labels``, the summand labels the element is admitted by (an object without them is not a carrier of a sum), and
``with_blocks(blocks)``, the same concrete class rebuilt with those blocks. A block whose summand is a
sum is handed to that summand's verb, which returns the block's own class; a block whose summand is not a sum is a
:class:`~orpheus.numerics.space.Member` whose values the summand's verb acts on and which :func:`with_values` rebuilds, in one place.

**The carrier.** A sum cannot mint its own typed element (the element class above ``numerics`` is not a function of the summand spaces), so the
layer that builds the sum installs the zero factories in its :class:`SumCarrier`; the dual of a sum carries the swapped pair.

**The metric** is block diagonal and every metric verb acts summand by summand through the summand's own verb: ``apply_metric`` and
``apply_inverse_metric`` take the typed carrier (returned as the same class) and the flat coordinates with a trailing payload kept
(:math:`G \otimes I_P`); ``inner_product`` and ``norm`` take the typed carrier only. A sum has no ``axes``.

What was tried and refused (the amendment record of N05): a second partition type beside the layout (the layout is the partition: a second type
would hold the same index sets twice); the sum stored in its carriers (a twin datum: the carrier derives its space from its blocks and its
labels); the label stored on the carrier with a class default (a default label is a positional guess); a flat-only sum (it deletes the typed chart
production uses); the typed element read by the space (it welds the space to member classes, which are W3's); a derived zero (the element class of
a full field is not a function of its summand spaces); one four-field carrier pair shared with the tensor space (a tensor space cannot hold a
factory and a sum cannot hold a member class).
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from functools import cached_property
from typing import Any, Callable, Optional, TypeIs

import numpy as np
from numpy.typing import NDArray
from scipy import linalg as scipy_linalg

from orpheus.numerics.coproduct_layout import CoproductLayout
from orpheus.numerics.metric import HilbertMetric
from orpheus.numerics.space import (
    FunctionSpace,
    Member,
    TensorSpace,
    _refused_legacy_metric_field,
)
from orpheus.values.content import name_digest

__all__ = [
    "BOUNDARY",
    "BULK",
    "INTERIOR",
    "TRACE",
    "DirectSumSpace",
    "SumCarrier",
    "bulk_and_trace",
    "with_values",
]

#: The summand labels of the full field: the bulk (volumetric) summand and the trace summand.
BULK = "bulk"
TRACE = "trace"
#: The summand labels of System B's ray instance: the interior summand and the boundary corner.
INTERIOR = "interior"
BOUNDARY = "boundary"


@dataclass(frozen=True)
class SumCarrier:  # concept: dsum | role: component | kind: set
    r"""The pair of zero factories of a sum: ``zeros`` mints the zero element of the sum and ``dual_zeros`` the zero element of its dual.

    A pair that was never installed stays ``None`` and its use raises. :meth:`dual` swaps the two, so the dual of a sum carries the swapped
    pair. The carrier is no part of a sum's identity (``compare=False`` on the sum); it is the member seam of the layer that builds the sum.

    SCOPE-BOUNDARY[guard]
    machinery: member classes that follow their space by type, so that the element class of a sum is a function of its space's content (W3)
    ruling: issue #707 (second comment, 2026-11-11), [ON04-2] and [ON05-4] of the amendment record in .claude/plans/ontology_w1_spec.md
    revisit: W3, the known defect: two equal sums carrying different pairs are equal and mint different element classes (the dual of a full field has no zero factory)
    """

    zeros: Optional[Callable[[], Any]] = None
    dual_zeros: Optional[Callable[[], Any]] = None

    def dual(self) -> "SumCarrier":
        r"""The pair of the dual space: the two factories swapped."""
        return SumCarrier(zeros=self.dual_zeros, dual_zeros=self.zeros)


def with_values(block: Any, values: NDArray) -> Any:
    r"""The leaf block ``block`` rebuilt with new ``values``: the one place a block that is not itself a sum is rebuilt, whether by a sum's verb or by a
    composite's ``from_flat`` (law L13: the leaf rebuild is spelled once).

    SCOPE-BOUNDARY[guard]
    machinery: members typed by their space, so that a leaf summand mints its own typed member (W3, #707)
    ruling: [ON05-16] of the amendment record in .claude/plans/ontology_w1_spec.md
    revisit: W3, when a leaf summand's ``member`` returns the typed member and this helper retires
    """
    return replace(block, values=values)


def _is_sum(space: Any) -> TypeIs["DirectSumSpace"]:
    return isinstance(space, DirectSumSpace)


@dataclass(frozen=True, eq=False)
class DirectSumSpace(FunctionSpace):  # concept: dsum | role: head | kind: space
    r"""The direct sum :math:`V = \bigoplus_l V_l` of labelled summand spaces, with the block-diagonal metric.

    Parameters
    ----------
    summands : tuple[tuple[str, FunctionSpace], ...]
        The ordered ``(label, space)`` pairs. At least one; the labels are non-empty strings, pairwise distinct. A summand that is itself a sum
        nests. The summands are typed :class:`~orpheus.numerics.space.FunctionSpace`; no runtime check refuses another object (a misuse is the
        caller's, and the formal path is :meth:`of_pairs`).
    carrier : SumCarrier, keyword-only, default ``SumCarrier()``
        The zero factories (``compare=False``, no part of the identity). The default installs none, and the doors :meth:`zero` and :meth:`dual_zero`
        are the one refusal of a sum that carries none; a consumer asks :attr:`has_zero` and calls the door, and never reads the carrier's fields.

    Attributes
    ----------
    layout : CoproductLayout[str]
        The labelled slices of the flat chart, derived from the summands' sizes (``init=False``).
    """

    # The leaf fields of the base are closed (``name`` and ``shape`` are derived, the two metric fields are refused on read), as on a tensor space.
    # ELEGANCE-DEBT[guard] #759 legacy-leaf: the four leaf fields are re-declared closed until N07 removes them from the base.
    name: str = field(init=False)
    shape: tuple[int, ...] = field(init=False)
    inner_product_weights: Optional[NDArray] = field(init=False, default=None, repr=False, compare=False)
    metric: Optional[HilbertMetric] = field(init=False, default=None, repr=False, compare=False, kw_only=True)
    summands: tuple[tuple[str, FunctionSpace], ...]
    layout: CoproductLayout[str] = field(init=False, repr=False, compare=False)
    carrier: SumCarrier = field(default_factory=SumCarrier, repr=False, compare=False, kw_only=True)

    def __post_init__(self) -> None:
        pairs = tuple(self.summands)
        if not pairs:
            raise ValueError("a direct sum needs at least one member space (summand): the empty sum has no coordinates")
        for pair in pairs:
            if not (isinstance(pair, tuple) and len(pair) == 2):
                raise TypeError(f"a direct sum is built from (label, space) pairs; got {pair!r}")
            if not (isinstance(pair[0], str) and pair[0]):
                raise TypeError(f"a summand label is a non-empty string; got {pair[0]!r}")
        labels = [label for label, _ in pairs]
        repeated = sorted({l for l in labels if labels.count(l) > 1})
        if repeated:
            raise ValueError(f"the labels of a direct sum are pairwise distinct; {repeated} repeat in {labels}")
        if not isinstance(self.carrier, SumCarrier):
            raise TypeError(f"carrier must be a SumCarrier, got {type(self.carrier).__name__}")
        layout = CoproductLayout(tuple((label, tuple(s.shape)) for label, s in pairs))
        if layout.sizes != tuple(int(s.n_coordinates) for _, s in pairs):
            raise ValueError(
                f"the layout sizes {layout.sizes} are not the summands' coordinate counts {tuple(s.n_coordinates for _, s in pairs)}"
            )
        digest = name_digest(tuple((label, s.name, s.shape) for label, s in pairs))
        object.__setattr__(self, "summands", pairs)
        object.__setattr__(self, "layout", layout)
        object.__setattr__(self, "shape", (layout.total_size,))
        object.__setattr__(self, "name", f"dsum({' ⊕ '.join(labels)})#{digest}")

    # ------------------------------------------------------------------
    # Construction and identity
    # ------------------------------------------------------------------

    @classmethod
    def of(cls, **summands: FunctionSpace) -> "DirectSumSpace":
        r"""The sum of the keyword summands, in keyword order: ``of(x=A, y=B)``."""
        return cls.of_pairs(tuple(summands.items()))

    @classmethod
    def of_pairs(
        cls, pairs: "tuple[tuple[str, FunctionSpace], ...]", *, carrier: Optional[SumCarrier] = None,
    ) -> "DirectSumSpace":
        r"""The sum of ``(label, space)`` pairs, in order. The carrier is keyword-only on this door alone, so no summand label can collide with it."""
        return cls(summands=tuple(pairs), carrier=SumCarrier() if carrier is None else carrier)

    def identity_key(self) -> tuple[Any, ...]:
        r"""The ordered ``(label, summand key)`` pairs, tagged: the labels and the order are identity."""
        return ("dsum", tuple((label, s.identity_key()) for label, s in self.summands))

    def __repr__(self) -> str:
        return f"DirectSumSpace({self.name!r}, shape={self.shape})"

    @property
    def labels(self) -> tuple[str, ...]:
        r"""The summand labels, in order."""
        return tuple(label for label, _ in self.summands)

    @cached_property
    def _by_label(self) -> dict[str, FunctionSpace]:
        return dict(self.summands)

    def summand(self, label: str) -> FunctionSpace:
        r"""The summand labelled ``label`` (a constant-time lookup); a ``KeyError`` listing the labels otherwise."""
        try:
            return self._by_label[label]
        except KeyError:
            raise KeyError(f"{label!r} is not a summand label; the labels are {list(self.labels)}") from None

    def leaf(self, *path: str) -> FunctionSpace:
        r"""The summand at a nested label path: ``leaf("a", "b")`` is ``summand("a").summand("b")``."""
        found: FunctionSpace = self
        for label in path:
            if not _is_sum(found):
                raise TypeError(f"leaf{path}: {found!r} is not a direct sum, so it has no summand {label!r}")
            found = found.summand(label)
        return found

    def with_summand(self, label: str, space: FunctionSpace) -> "DirectSumSpace":
        r"""This sum with the summand ``label`` replaced by ``space`` (the order and the labels kept; the carrier is not carried over, since its
        factories mint elements of the old summand)."""
        self.summand(label)
        return DirectSumSpace.of_pairs(tuple((l, space if l == label else s) for l, s in self.summands))

    def slice_of(self, label: str) -> slice:
        r"""The slice of the flat chart that holds summand ``label``."""
        return self.layout.slice_of(label)

    def axis_index_of_coordinates(self, label: str) -> NDArray:
        r"""The index along the axis ``label`` of every flat coordinate: each summand's answer, placed through the layout (the sum has no axes of its own)."""
        out = np.zeros(self.n_coordinates, dtype=np.intp)
        for name, s in self.summands:
            out[self.layout.slice_of(name)] = s.axis_index_of_coordinates(label)
        return out

    def pack(self, parts: Any) -> NDArray:
        r"""The injections summed: the flat chart of the parts ``{label: array of the summand's shape (and the common payload)}``."""
        return self.layout.pack(parts)

    def unpack(self, flat: Any) -> dict[str, NDArray]:
        r"""The projections: ``{label: the summand's part of the flat chart}``, as views."""
        return self.layout.unpack(flat)

    # ------------------------------------------------------------------
    # The duals: summand-wise, and by label for the factors
    # ------------------------------------------------------------------

    def _make_dual(self) -> "DirectSumSpace":
        return DirectSumSpace.of_pairs(tuple((label, s.dual()) for label, s in self.summands), carrier=self.carrier.dual())

    def _held(self) -> frozenset[str]:
        r"""The factor labels held by the summands of this sum, or a refusal where a summand is a legacy leaf.

        ELEGANCE-DEBT[guard] #759 legacy-leaf: a flip addressed by label needs factors, which a legacy leaf does not hold; retires at N07, when no space is a leaf.
        """
        held: set[str] = set()
        for label, s in self.summands:
            # ELEGANCE-DEBT[guard] #759 legacy-leaf: the door tells a tensor space from a legacy leaf and nothing else;
            # it retires at N07, when no space is a legacy leaf (G04o pins it).
            if isinstance(s, TensorSpace):
                held |= set(s.labels)
            elif _is_sum(s):
                held |= s._held()
            else:
                raise TypeError(
                    f"a legacy leaf holds no factor labels: the summand {label!r} ({s!r}) is not built from factors, so a flip addressed to factor labels "
                    f"cannot be applied to this sum; flip it whole with dual()"
                )
        return frozenset(held)

    def _on(self, on: Any) -> frozenset[str]:
        if isinstance(on, str):
            raise TypeError(f"on= takes a collection of factor labels, not the string {on!r}")
        labels = frozenset(on)
        if not labels:
            return labels
        held = self._held()
        unknown = sorted(labels - held)
        if unknown:
            raise KeyError(f"unknown factor labels {unknown}; the summands hold {sorted(held)} (summand labels are not factor labels)")
        return labels

    def _flip(self, on: Any, verb: Callable[[FunctionSpace, Any], FunctionSpace]) -> "DirectSumSpace":
        labels = self._on(on)
        if not labels:
            return self
        pairs = []
        for label, s in self.summands:
            held = s._held() if _is_sum(s) else frozenset(s.labels)  # type: ignore[attr-defined]
            mine = labels & held
            pairs.append((label, verb(s, mine) if mine else s))
        return DirectSumSpace.of_pairs(tuple(pairs))

    def dual(self, on: Any = None) -> "DirectSumSpace":
        r"""The dual: summand-wise, memoised (``V.dual().dual() is V``); with ``on=`` a set of factor labels, the flip of the factors with those labels
        in the summands that hold them (law L6). An unknown label is a ``KeyError``; a sum with a legacy-leaf summand refuses a non-empty ``on=``.
        The partial flip (``on=``), like :meth:`riesz`, returns a sum with no carrier: the swapped carrier pair of the whole dual does not mint the
        elements of a half-flipped space, as :meth:`with_summand` also drops it."""
        if on is None:
            return self._dual  # type: ignore[return-value]
        return self._flip(on, lambda s, mine: s.dual(on=mine))  # type: ignore[call-arg]

    def riesz(self, on: Any = None) -> "DirectSumSpace":
        r"""The Riesz image: summand-wise (the chart flipped through each form), or with ``on=`` the factors with those labels (law L6).

        Refused where a flipped slot is formless (the summand's own refusal) and where a summand is a legacy leaf, whole: nothing half applied is returned.
        The result carries no carrier (a Riesz image is neither the dual nor the primal space, so no installed pair mints its elements).
        """
        if on is None:
            self._held()
            return DirectSumSpace.of_pairs(tuple((label, s.riesz()) for label, s in self.summands))  # type: ignore[attr-defined]
        return self._flip(on, lambda s, mine: s.riesz(on=mine))  # type: ignore[attr-defined]

    # ------------------------------------------------------------------
    # The typed element and the two charts
    # ------------------------------------------------------------------

    def _blocks_of(self, x: Any, verb: str) -> tuple[Any, ...]:
        r"""The blocks of a typed element: the one door through which every typed verb reaches them. A refusal names what ``verb`` takes; the blocks are
        admitted by arity and by label, and an object without labels is not a carrier of a sum.

        ELEGANCE-DEBT[guard] #666: the summand of each block is not compared with the block's own space, so the verbs hand the summands' verbs bare
        coordinates. A composite on a DUAL space still labels its blocks with the primal spaces (the role of a carrier under a dual arrow is #633), so a typed
        block would be refused there. Retires when a composite member is compared with this space by content and the summand verbs take typed blocks.
        """
        if not (hasattr(x, "blocks") and hasattr(x, "with_blocks")):
            raise AttributeError(
                f"DirectSumSpace.{verb} takes the typed element of the sum (a carrier with blocks, with_blocks and labels); got {type(x).__name__}"
            )
        blocks = tuple(x.blocks)
        if len(blocks) != len(self.summands):
            raise ValueError(
                f"DirectSumSpace.{verb}: the element carries {len(blocks)} blocks but the sum has {len(self.summands)} summands {list(self.labels)} — "
                f"the field/space pairing is inconsistent."
            )
        labels = x.labels
        if tuple(labels) != self.labels:
            raise ValueError(
                f"DirectSumSpace.{verb}: the element is labelled {list(labels)} but the sum's summands are {list(self.labels)}; "
                f"the blocks are admitted by label, never by position."
            )
        return blocks

    def _on_blocks(self, x: Any, verb: str, leg: Callable[[FunctionSpace, Any], Any]) -> Any:
        r"""A block-diagonal action on the typed element, returned as the element's own class: each block goes to its summand's ``leg``."""
        blocks = self._blocks_of(x, verb)
        out = []
        for (_, s), block in zip(self.summands, blocks):
            out.append(leg(s, block) if _is_sum(s) else with_values(block, np.asarray(leg(s, block.values))))
        return x.with_blocks(tuple(out))

    def _on_flat(self, x: Any, leg: Callable[[FunctionSpace, Any], Any]) -> NDArray:
        r"""A block-diagonal action on the flat coordinates, with a trailing payload kept (:math:`G \otimes I_P`)."""
        coordinates = np.asarray(self.leading_coordinates_of(x))
        payload = tuple(coordinates.shape[1:])
        parts = self.layout.unpack(coordinates)
        return self.layout.pack({
            label: np.asarray(leg(s, parts[label])).reshape(tuple(s.shape) + payload) for label, s in self.summands
        })

    def apply_metric(self, x: Any) -> Any:
        r"""Apply the block-diagonal Hilbert metric :math:`G\odot x`, summand by summand, to the typed element (returned as its own class) or to the flat coordinates."""
        leg = lambda s, v: s.apply_metric(v)  # noqa: E731
        return self._on_flat(x, leg) if _is_flat(x) else self._on_blocks(x, "apply_metric", leg)

    def apply_inverse_metric(self, x: Any) -> Any:
        r"""Apply the block-diagonal pseudo-inverse metric :math:`G^{+}\odot x`, summand by summand, in either chart."""
        leg = lambda s, v: s.apply_inverse_metric(v)  # noqa: E731
        return self._on_flat(x, leg) if _is_flat(x) else self._on_blocks(x, "apply_inverse_metric", leg)

    def inner_product(self, x: Any, y: Any) -> float:
        r"""The direct-sum inner product :math:`\sum_l \langle x_l, y_l\rangle_{G_l}` of two typed elements (the blocks through :meth:`_blocks_of`)."""
        total = 0.0
        for (_, s), bx, by in zip(self.summands, self._blocks_of(x, "inner_product"), self._blocks_of(y, "inner_product")):
            total += s.inner_product(bx, by) if _is_sum(s) else s.inner_product(bx.values, by.values)
        return float(total)

    def metric_diagonal(self) -> Optional[NDArray]:
        r"""The summands' diagonals concatenated in label order, or ``None`` if any summand is not diagonal."""
        diagonals = tuple(s.metric_diagonal() for _, s in self.summands)
        if any(d is None for d in diagonals):
            return None
        return np.concatenate([np.asarray(d).reshape(-1) for d in diagonals])

    @property
    def n_coordinates(self) -> int:
        return int(self.layout.total_size)

    @property
    def rank(self) -> int:
        r"""The sum of the summands' ranks."""
        return int(sum(s.rank for _, s in self.summands))

    def kernel_basis(self) -> NDArray:
        r"""The summands' kernel bases, block diagonally, in the flat chart's label order."""
        return np.asarray(scipy_linalg.block_diag(*(s.kernel_basis() for _, s in self.summands)))

    @property
    def negative(self) -> int:
        r"""The negative directions of the sum: the sum over the summands', so that a signed summand refuses :meth:`norm` (issue #706)."""
        return int(sum(s.negative for _, s in self.summands))

    @property
    def positive(self) -> bool:
        r"""Whether every summand's form has no negative direction (what ``_require_positive`` asks of a space)."""
        return all(s.positive for _, s in self.summands)

    @property
    def unit(self) -> NDArray:
        r"""The coordinates of the constant function 1: the summands' units concatenated in label order."""
        return np.concatenate([np.asarray(s.unit).reshape(-1) for _, s in self.summands])

    @property
    def has_coordinate_cone(self) -> bool | None:
        r"""The three-valued conjunction of the summands' answers: False if any summand answers False, otherwise None if any answers None, otherwise True."""
        answers = [s.has_coordinate_cone for _, s in self.summands]
        if any(a is False for a in answers):
            return False
        if any(a is None for a in answers):
            return None
        return True

    # ------------------------------------------------------------------
    # The element: the inverse of the coordinate chart
    # ------------------------------------------------------------------

    @property
    def has_zero(self) -> bool:
        r"""Whether the layer that built this sum installed its zero-element factory."""
        return self.carrier.zeros is not None

    @property
    def has_dual_zero(self) -> bool:
        r"""Whether the layer that built this sum installed the zero-element factory of its dual."""
        return self.carrier.dual_zeros is not None

    def zero(self) -> Any:
        r"""The zero typed element of the sum, minted by the factory its builder installed; the one refusal of a sum that carries none.

        SCOPE-BOUNDARY[guard]
        machinery: member classes that follow their space by type, so that the element class of a sum is a function of its space's content (W3)
        ruling: issue #707 (second comment, 2026-11-11) and [ON05-4] of the amendment record in .claude/plans/ontology_w1_spec.md
        revisit: W3, when the zero element is derived from the summands and no factory is installed
        """
        if self.carrier.zeros is None:
            raise RuntimeError(
                f"{self!r} carries no zero-element factory: the layer that builds the sum installs its SumCarrier "
                f"(the member seam, #705; the typed elements live above the numerics layer)."
            )
        return self.carrier.zeros()

    def dual_zero(self) -> Any:
        r"""The zero typed element of the dual of the sum (the cotangent-side zero), minted by the factory its builder installed; refused where none is.

        SCOPE-BOUNDARY[guard]
        machinery: member classes that follow their space by type, so that the element class of a sum is a function of its space's content (W3)
        ruling: issue #707 (second comment, 2026-11-11) and [ON05-4] of the amendment record in .claude/plans/ontology_w1_spec.md
        revisit: W3, when the zero element of the dual is derived from the summands and no factory is installed
        """
        if self.carrier.dual_zeros is None:
            raise RuntimeError(
                f"{self!r} carries no DUAL zero-element factory: the layer that builds the sum installs its SumCarrier "
                f"(the member seam, #705; the typed elements live above the numerics layer)."
            )
        return self.carrier.dual_zeros()

    def member(self, coordinates: NDArray) -> Any:
        r"""The typed element with these flat coordinates, on the zero element's class: the inverse of :meth:`coordinates_of`."""
        return self.with_coordinates(self.zero(), coordinates)

    def with_coordinates(self, x: Any, coordinates: NDArray) -> Any:
        r"""The typed element ``x``, its class and everything it carries kept, rebuilt with these flat coordinates.

        Each block is rebuilt from its part of the coordinates: a block whose summand is a sum by that summand's own :meth:`with_coordinates`, any other
        block through :func:`with_values`. This is the one place a flat chart is cut into blocks (law L13): the cut is the layout's.
        """
        flat = np.asarray(coordinates, dtype=float).ravel()
        if flat.size != self.n_coordinates:
            raise ValueError(
                f"DirectSumSpace.with_coordinates: {flat.size} coordinates given, the sum has {self.n_coordinates} ({list(self.labels)} of sizes {self.layout.sizes})"
            )
        parts = self.layout.unpack(flat)
        blocks = []
        for (label, s), block in zip(self.summands, self._blocks_of(x, "with_coordinates")):
            blocks.append(s.with_coordinates(block, parts[label]) if _is_sum(s) else with_values(block, parts[label]))
        return x.with_blocks(tuple(blocks))

    def coordinates_of(self, x: Any) -> Any:
        r"""The flat coordinates of ``x``: a typed element (its ``to_flat``, summands in label order) after the identity leg (its space equals this one),
        a plain :class:`~orpheus.numerics.space.Member` of the sum, or a bare array (the base's shape leg)."""
        if isinstance(x, (np.ndarray, Member)) or not hasattr(x, "blocks"):
            return super().coordinates_of(x)
        self._blocks_of(x, "coordinates_of")
        space = x.space
        if space != self:
            raise ValueError(f"{type(x).__name__} on {space!r} is not a member of {self!r}.")
        return np.asarray(x.to_flat(), dtype=float)

    # ------------------------------------------------------------------
    # The leaf-arm axis verbs: a sum has no named factors of its own
    # ------------------------------------------------------------------

    def _no_factors(self, verb: str) -> TypeError:
        r"""The refusal of an axis verb: the factors are the summands', addressed through ``summand(label)``.

        ELEGANCE-DEBT[guard] #759 legacy-leaf: retires at N07, when the axis verbs leave the base for the tensor space and this refusal is the absence of the verb.
        """
        return TypeError(f"{verb}: {self!r} is a direct sum with no named factors of its own; ask a summand, summand(label), one of {list(self.labels)}")

    def _axis_index(self, label: str) -> int:
        raise self._no_factors(f"axis lookup {label!r}")

    def axis(self, label: str) -> Any:
        raise self._no_factors(f"axis {label!r}")

    def _axis_block(self, k: int) -> Any:
        raise self._no_factors(f"axis block {k}")

    def without_axis(self, label: str) -> Any:
        raise self._no_factors(f"without_axis {label!r}")

    def with_axis(self, factor: Any) -> Any:
        raise self._no_factors("with_axis")

    def positions_spanned_by(self, sub: Any) -> Any:
        raise self._no_factors("positions_spanned_by")


def _is_flat(x: Any) -> bool:
    r"""Whether ``x`` is read in the flat chart: anything that is not a typed element (a carrier with ``blocks`` and ``with_blocks``)."""
    return not (hasattr(x, "blocks") and hasattr(x, "with_blocks"))


# ELEGANCE-DEBT[guard] #759 legacy-leaf: the two metric fields of the leaf arm are refused on a sum, as on a tensor space; retires at N07, when the fields leave the base.
for _field_name in ("inner_product_weights", "metric"):
    setattr(DirectSumSpace, _field_name, _refused_legacy_metric_field(_field_name, "DirectSumSpace"))
del _field_name


def bulk_and_trace(space: Any, *, owner: str) -> tuple[FunctionSpace, FunctionSpace]:
    r"""The bulk and the trace summands of a full-field sum, or a refusal naming ``owner``.

    Accepts by the labels (``bulk``, ``trace``) and never by arity: a ray instance (``interior``, ``boundary``) and a coupled sum are refused.

    SCOPE-BOUNDARY[guard]
    machinery: operator ends typed by their space, so that an operator's domain says which sum it takes (W2 of the chart, operators on true ends)
    ruling: [ON05-3] and [ON05-19] of the amendment record in .claude/plans/ontology_w1_spec.md
    revisit: W2
    """
    if isinstance(space, DirectSumSpace) and space.labels == (BULK, TRACE):
        return space.summand(BULK), space.summand(TRACE)
    labels = list(space.labels) if isinstance(space, DirectSumSpace) else None
    raise TypeError(
        f"{owner} needs a composite DirectSumSpace labelled ({BULK!r}, {TRACE!r}) (the full field); got {space!r}"
        + ("" if labels is None else f" with labels {labels}")
    )

