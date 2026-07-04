r"""The composite full-field space :math:`V = V_{\rm bulk} \oplus V_{\rm trace}`.

The carrier space of the FULL SN operator :math:`L` and of every composite
that couples bulk and boundary — :math:`(L + C - S - F - B)` and its
within-group sub-sums. Its elements are bulk :math:`\oplus` boundary
*composite fields* (an
:class:`~orpheus.transport.timed_full_field.TimedFullField`: a bulk
:class:`~orpheus.transport.fields._bases.BulkField` paired with a boundary
:class:`~orpheus.transport.fields._bases.AngularBoundaryField`).

Why a direct-sum space (Wave O / O.2b R5)
=========================================

The Hilbert adjoint ``A.H`` is the **G-adjoint** :math:`A^\dagger =
G^{-1} A^{\mathsf T} G`, NOT the plain Euclidean transpose. For an SN
loss operator the metric :math:`G` is **block-diagonal** on the direct
sum:

* **bulk block** :math:`G_{\rm bulk} = V_{\rm cell}\,w_n` — the full
  phase-space measure :math:`\mathrm{d}V\,\mathrm{d}\Omega` (cell volume
  times angular quadrature weight);
* **trace block** :math:`G_{\rm trace} = |\Omega\cdot\hat n_f|\,w_n` —
  the partial-current surface measure (the cosine-weighted angular
  quadrature under which the boundary operators are self-adjoint; see
  :mod:`orpheus.numerics.spaces.angular_trace_space`).

Both carry :math:`w_n`; they differ in the spatial factor (volume vs.
oriented surface). :class:`~orpheus.numerics.operator._AdjointOperator`
realises ``A.H`` as
:math:`G^{-1}\!\odot\,A^{\mathsf T}\!\big(G\odot y\big)` by calling
:meth:`apply_metric` on the codomain (here) before the transpose and
:meth:`apply_inverse_metric` on the domain (here) after — so this space
is exactly the object that teaches the otherwise metric-agnostic adjoint
wrapper *which* per-block metric to apply.

Composition over a new metric kernel (Pattern 2 / Pattern 5)
============================================================

:class:`FullFieldSpace` owns ONLY the direct-sum *structure*: it splits a
composite field into its ``bulk`` / ``boundary`` blocks and routes each to
the matching leaf :class:`~orpheus.numerics.space.FunctionSpace`. Each leaf
space already implements the metric primitives
(:meth:`~orpheus.numerics.space.FunctionSpace.apply_metric` /
:meth:`~orpheus.numerics.space.FunctionSpace.apply_inverse_metric` /
:meth:`~orpheus.numerics.space.FunctionSpace.inner_product`) — including
the Moore–Penrose pseudo-inverse the singular partial-current trace metric
needs (``G_trace = 0`` on the tangential ordinates where
:math:`|\Omega\cdot\hat n| = 0`). The composite therefore introduces **no
new metric arithmetic**; it is a pure direct-sum dispatcher. The
pseudo-inverse is exact for the adjoint: the tangential trace slots carry
zero :math:`\langle\cdot,\cdot\rangle_G` weight and are identically zero in
every matvec output (the inflow / outflow selectors exclude tangential
ordinates), so they never appear on the adjoint's null space in practice.

Composite-field contract (duck-typed — no transport import)
===========================================================

This space lives in the ``numerics`` layer and must not import the
``transport`` layer. It operates on the composite field structurally: the
field is a frozen dataclass exposing ``.bulk`` and ``.boundary`` leaf
fields, each itself a frozen dataclass with a ``.values`` ndarray. The
metric methods rebuild the composite with :func:`dataclasses.replace` —
no concrete type import is needed. This mirrors the duck-typed contract
:class:`~orpheus.numerics.operator._AdjointOperator` already relies on
when it calls ``space.apply_metric(y)`` on a composite field.

Identity
========

Identity is the inherited ``(name, shape)`` tuple, with ``name =
"full_field"`` and ``shape = (n_bulk + n_trace,)`` (the flat direct-sum
dimension). The name is method-agnostic (P4.5 W-D de-SN-ified it from the
former ``"sn_full_field"``: the operators that advertise it — :math:`C`,
:math:`S`, :math:`F` — are cross-method, not SN-specific). The block spaces
are ``compare=False`` leaf metadata — two composites over meshes of the same
total dimension compare equal, so the
:class:`~orpheus.numerics.operator.OperatorSum` composition guard accepts the
full within-group loss ``L + C - S - B``: every bulk operand (:math:`L`,
:math:`C`, :math:`S`) and the boundary :math:`B` reports the SAME composite
domain. (P4.5 W-D gave :math:`C`/:math:`S`/:math:`F` real spaces; before that
only :math:`L`/:math:`B` did, and the guard silently skipped the
``None``-spaced summands.)

References
----------

* Lewis, E.E. & Miller, W.F. (1993). *Computational Methods of Neutron
  Transport*. ANS. §3.7 (boundary partial-current inner product), §6
  (curvilinear angular redistribution).
* Trefethen, L.N. & Bau, D. (1997). *Numerical Linear Algebra*. SIAM. §1
  (the Hilbert adjoint :math:`G^{-1}A^{\mathsf T}G` vs. the representation
  transpose).
* Issue #208 / ``.claude/plans/glimmering-launching-lantern.md`` — the
  Wave-O operator-algebra G-adjoint design (R5).
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import TYPE_CHECKING, Any, ClassVar, Optional, Protocol

import numpy as np

from orpheus.numerics.space import FunctionSpace

if TYPE_CHECKING:
    from orpheus.numerics.spaces.angular_trace_space import AngularTraceSpace


__all__ = ["CompositeField", "FullFieldSpace"]


class _CompositeLeaf(Protocol):
    """One block of a composite field: a frozen dataclass carrying ``values``.

    The ``__dataclass_fields__`` ClassVar makes the protocol satisfy
    typeshed's ``DataclassInstance``, so :func:`dataclasses.replace`
    accepts a leaf without a concrete transport import.
    """

    __dataclass_fields__: ClassVar[dict[str, Any]]

    @property
    def values(self) -> np.ndarray: ...


class CompositeField(Protocol):
    """The composite-field contract (see module docstring) made static.

    The carrier of :class:`FullFieldSpace` — a bulk ⊕ boundary
    (⊕ optional starting-direction, #282 route (a)) block tuple with
    the polymorphic ``_recombine`` rebuild hook. Structural, so the
    transport carriers (``FullField`` / ``TimedFullField``) satisfy it
    without an import edge out of the numerics layer.
    """

    @property
    def bulk(self) -> _CompositeLeaf: ...

    @property
    def boundary(self) -> _CompositeLeaf: ...

    @property
    def starting_direction(self) -> Optional[_CompositeLeaf]: ...

    def _recombine(
        self,
        *,
        bulk: _CompositeLeaf,
        boundary: _CompositeLeaf,
        starting_direction: Optional[_CompositeLeaf],
    ) -> "CompositeField": ...


@dataclass(frozen=True)
class FullFieldSpace(FunctionSpace[CompositeField]):
    r"""Direct sum :math:`V_{\rm bulk} \oplus V_{\rm trace}` with a per-block metric.

    Parameters
    ----------
    name, shape, inner_product_weights
        Inherited from :class:`~orpheus.numerics.space.FunctionSpace`.
        ``name`` is ``"full_field"`` and ``shape`` is the flat
        direct-sum dimension ``(n_bulk + n_trace,)``. The composite's own
        ``inner_product_weights`` stays ``None`` — the metric is carried
        per block by :attr:`bulk_space` / :attr:`trace_space`, and the
        block-aware :meth:`apply_metric` / :meth:`apply_inverse_metric` /
        :meth:`inner_product` overrides never read the base slot.
    bulk_space : FunctionSpace
        The bulk leaf space, carrying the phase-space metric
        ``G_bulk = V_cell · w_n`` as its ``inner_product_weights`` (shape
        broadcast against the ``(N, ng, nx, ny)`` bulk tensor). ``compare
        =False`` leaf metadata (not part of the ``(name, shape)`` identity).
    trace_space : FunctionSpace
        The boundary leaf space — the angular
        :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`
        (partial-current metric ``G_trace = |Ω·n|·w_n``) for SN, the
        :class:`~orpheus.numerics.spaces.scalar_trace_space.ScalarTraceSpace`
        (face-area metric over ``(J⁺, J⁻)`` pairs) for diffusion / CP
        (#290 P2). The composite reads only the carrier-generic
        FunctionSpace metric surface, so the block dispatch is
        family-blind. ``compare=False`` leaf metadata.
    starting_direction_space : FunctionSpace, optional
        The OPTIONAL third block (#282 route (a), 2.5d): the
        :class:`~orpheus.numerics.spaces.starting_direction_space.StartingDirectionSpace`
        carrying the per-level ψ½ layout and the ALL-ZERO ghost metric
        (``G_sd = 0`` — the starting-direction rays carry no quadrature
        measure; the masked pseudo-inverse and the zero inner-product
        contribution follow from the leaf space's zero weights, no
        composite-level arithmetic). ``None`` on meshes with no
        seed-carrying level (R12a) and for non-SN methods.
        ``compare=False`` leaf metadata.
    """

    bulk_space: Optional[FunctionSpace] = field(
        default=None, repr=False, compare=False,
    )
    trace_space: Optional[FunctionSpace] = field(
        default=None, repr=False, compare=False,
    )
    starting_direction_space: Optional[FunctionSpace] = field(
        default=None, repr=False, compare=False,
    )

    # ── Equality / hashing inherited from FunctionSpace ───────────────
    #
    # Same rationale as :class:`AngularTraceSpace` / :class:`TensorProductSpace`:
    # the @dataclass(frozen=True) decorator would otherwise auto-generate
    # an __eq__ / __hash__ over the fields; explicit delegation restores
    # the (name, shape) identity convention. The block fields are already
    # excluded via ``compare=False``.

    def __eq__(self, other: object) -> bool:
        return FunctionSpace.__eq__(self, other)

    def __hash__(self) -> int:
        return FunctionSpace.__hash__(self)

    def __repr__(self) -> str:
        return f"FullFieldSpace({self.name!r}, shape={self.shape})"

    @classmethod
    def from_blocks(
        cls,
        bulk_space: FunctionSpace,
        trace_space: FunctionSpace,
        starting_direction_space: Optional[FunctionSpace] = None,
    ) -> "FullFieldSpace":
        r"""Build the composite from its bulk and trace (and optional seed) leaf spaces.

        Derives the flat direct-sum ``shape`` from the leaf shapes:
        ``(prod(bulk) + prod(trace) [+ prod(starting_direction)],)`` —
        ``prod`` on every block so the identity dimension is robust to a
        future multi-axis trace (today ``trace_space.shape ==
        (total_size,)``). The starting-direction block is present
        exactly when the mesh predicate supplies a space (R12a — see
        :attr:`starting_direction_space`).
        """
        n_bulk = int(np.prod(bulk_space.shape))
        n_trace = int(np.prod(trace_space.shape))
        n_seed = (
            0
            if starting_direction_space is None
            else int(np.prod(starting_direction_space.shape))
        )
        return cls(
            name="full_field",
            shape=(n_bulk + n_trace + n_seed,),
            bulk_space=bulk_space,
            trace_space=trace_space,
            starting_direction_space=starting_direction_space,
        )

    # ------------------------------------------------------------------
    # Direct-sum metric dispatch (per-block, on a composite field)
    # ------------------------------------------------------------------

    def _require_blocks(self) -> tuple[FunctionSpace, FunctionSpace]:
        r"""Return the ``(bulk_space, trace_space)`` pair, guarding the bare-constructor footgun.

        The ``bulk_space`` / ``trace_space`` fields default to ``None`` (the
        ``compare=False`` dataclass-field convention), so a composite built via
        the bare constructor instead of :meth:`from_blocks` would
        ``AttributeError`` deep inside the adjoint path. Fail at the boundary
        with intent instead (parse-don't-validate), and return the narrowed
        non-``None`` pair so callers bind locals that type cleanly.
        """
        if self.bulk_space is None or self.trace_space is None:
            raise RuntimeError(
                "FullFieldSpace has no block spaces; build it via "
                "FullFieldSpace.from_blocks(bulk_space, trace_space) (or "
                "SNMesh.full_field_space), not the bare dataclass constructor."
            )
        return self.bulk_space, self.trace_space

    def _seed_space_for(self, x: CompositeField) -> Optional[FunctionSpace]:
        r"""The starting-direction leaf space for ``x``'s seed block, or ``None``.

        Field-driven dispatch (2.5d transitional discipline): the seed
        arithmetic runs exactly when the FIELD carries the block. A
        field WITHOUT the block on a space WITH one is legal (the space
        grew at d1; production fields adopt the block at d3) — the seed
        paths are simply not exercised. The REVERSE — a field carrying a
        ψ½ block through a composite space that has no seed leaf space —
        is an illegal pairing and raises.
        """
        if getattr(x, "starting_direction", None) is None:
            return None
        if self.starting_direction_space is None:
            raise RuntimeError(
                "FullFieldSpace: the composite field carries a "
                "starting_direction block but this space has no "
                "starting_direction_space — the field/space pairing is "
                "inconsistent (build the space via "
                "FullFieldSpace.from_blocks(..., starting_direction_space=...) "
                "/ SNMesh.full_field_space)."
            )
        return self.starting_direction_space

    @staticmethod
    def _rebuild(
        x: CompositeField,
        bulk_values: np.ndarray,
        boundary_values: np.ndarray,
        starting_direction_values: Optional[np.ndarray] = None,
    ) -> CompositeField:
        r"""Return a copy of composite field ``x`` with new block ``values``.

        Rebuilds the frozen leaves (``x.bulk`` / ``x.boundary`` / the
        optional ``x.starting_direction``) via
        :func:`dataclasses.replace` — preserving each leaf's ``space`` /
        ``mesh`` — then routes the recombined blocks through the composite's
        polymorphic :meth:`_recombine` hook.  The hook gives the correct
        concrete type for either carrier: a timeless
        :class:`~orpheus.transport.full_field.FullField` rebuilds a
        ``FullField`` (no history concept), a timed
        :class:`~orpheus.transport.timed_full_field.TimedFullField` rebuilds a
        ``TimedFullField`` with its ``history_depth`` and an EMPTY history
        (a metric application is algebra, not a state step).  #257 S8a — since
        operator matvec outputs are now the TIMELESS ``FullField`` (base
        arrows), the G-adjoint metric path receives either carrier; routing
        through ``_recombine`` (not a hardcoded ``_history=()``) handles both.
        No concrete type import (duck-typed on the ``_recombine`` hook).
        ``starting_direction_values`` must be supplied exactly when ``x``
        carries the seed block (the callers thread it from
        :meth:`_seed_space_for`'s field-driven dispatch).
        """
        seed = x.starting_direction
        if seed is not None and starting_direction_values is None:
            raise RuntimeError(
                "FullFieldSpace._rebuild: the composite carries a "
                "starting_direction block but no rebuilt seed values were "
                "supplied — a metric path silently dropping the ψ½ block "
                "is the silent-drop bug class (§16.A A3)."
            )
        return x._recombine(
            bulk=replace(x.bulk, values=bulk_values),
            boundary=replace(x.boundary, values=boundary_values),
            starting_direction=(
                None
                if seed is None
                else replace(seed, values=starting_direction_values)
            ),
        )

    def apply_metric(self, x: CompositeField) -> CompositeField:
        r"""Apply the block-diagonal Hilbert metric :math:`G\odot x`.

        ``x`` is a composite field; the bulk block is weighted by
        ``G_bulk = V·w_n`` and the trace block by ``G_trace = |Ω·n|·w_n``,
        each delegated to the matching leaf space. A carried
        starting-direction block is weighted by its ALL-ZERO ghost
        metric — the output seed block is identically zero (the
        Mode-12 metric-blindness mechanism, documented on
        :mod:`orpheus.numerics.spaces.starting_direction_space`).
        """
        bulk_space, trace_space = self._require_blocks()
        seed_space = self._seed_space_for(x)
        return self._rebuild(
            x,
            bulk_space.apply_metric(x.bulk.values),
            trace_space.apply_metric(x.boundary.values),
            None
            if seed_space is None
            else seed_space.apply_metric(
                x.starting_direction.values  # type: ignore[union-attr]
            ),
        )

    def apply_inverse_metric(self, x: CompositeField) -> CompositeField:
        r"""Apply the block-diagonal pseudo-inverse metric :math:`G^{+}\odot x`.

        Plain ``1/G`` on the strictly-positive bulk block; the
        Moore–Penrose masked inverse on the trace block (zero on the
        tangential null space ``|Ω·n| = 0``) and on the ENTIRE
        starting-direction block (its ghost metric is all-zero, so the
        whole block is metric null space). Delegated per block.
        """
        bulk_space, trace_space = self._require_blocks()
        seed_space = self._seed_space_for(x)
        return self._rebuild(
            x,
            bulk_space.apply_inverse_metric(x.bulk.values),
            trace_space.apply_inverse_metric(x.boundary.values),
            None
            if seed_space is None
            else seed_space.apply_inverse_metric(
                x.starting_direction.values  # type: ignore[union-attr]
            ),
        )

    def inner_product(self, x: CompositeField, y: CompositeField) -> float:
        r"""Return the direct-sum inner product
        :math:`\langle x, y\rangle_G = \langle x_{\rm bulk}, y_{\rm bulk}\rangle_{G_{\rm bulk}}
        + \langle x_{\rm trace}, y_{\rm trace}\rangle_{G_{\rm trace}}
        \;(+\;0)`.

        A carried starting-direction block contributes EXACTLY zero
        (all-zero ghost metric) — the term is still evaluated through
        the leaf space so the accounting is visible, and presence must
        MATCH between ``x`` and ``y`` (mixed presence is the silent-drop
        bug class).
        """
        bulk_space, trace_space = self._require_blocks()
        x_seed = getattr(x, "starting_direction", None)
        y_seed = getattr(y, "starting_direction", None)
        if (x_seed is None) != (y_seed is None):
            raise ValueError(
                "FullFieldSpace.inner_product: mixed starting-direction "
                "presence between the two composite operands — on a "
                "seed-carrying mesh every composite must carry the block."
            )
        seed_term = 0.0
        if x_seed is not None:
            seed_space = self._seed_space_for(x)
            assert seed_space is not None  # narrowed by _seed_space_for
            seed_term = seed_space.inner_product(x_seed.values, y_seed.values)  # type: ignore[union-attr]
        return (
            bulk_space.inner_product(x.bulk.values, y.bulk.values)
            + trace_space.inner_product(
                x.boundary.values, y.boundary.values
            )
            + seed_term
        )
