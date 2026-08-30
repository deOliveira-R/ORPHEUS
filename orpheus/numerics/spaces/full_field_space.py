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
oriented surface). :class:`~orpheus.numerics.operator.AdjointOperator`
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
field is a frozen dataclass exposing ``.interior`` and ``.boundary`` leaf
fields, each itself a frozen dataclass with a ``.values`` ndarray. The
metric methods rebuild the composite with :func:`dataclasses.replace` —
no concrete type import is needed. This mirrors the duck-typed contract
:class:`~orpheus.numerics.operator.AdjointOperator` already relies on
when it calls ``space.apply_metric(y)`` on a composite field.

Identity
========

Identity is the inherited ``(name, shape)`` tuple, with ``name =
"full_field"`` and ``shape = (n_interior + n_trace,)`` (the flat direct-sum
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

import hashlib
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

    The carrier of :class:`FullFieldSpace` — a bulk ⊕ boundary block pair
    with the polymorphic ``_recombine`` rebuild hook. Structural, so the
    transport carriers (``FullField`` / ``TimedFullField`` / the System-B
    ``RadialCharacteristicField``) satisfy it without an import edge
    out of the numerics layer.
    """

    @property
    def interior(self) -> _CompositeLeaf: ...

    @property
    def boundary(self) -> _CompositeLeaf: ...

    def _recombine(
        self,
        *,
        interior: _CompositeLeaf,
        boundary: _CompositeLeaf,
    ) -> "CompositeField": ...


@dataclass(frozen=True)
class FullFieldSpace(FunctionSpace[CompositeField]):
    r"""Direct sum :math:`V_{\rm bulk} \oplus V_{\rm trace}` with a per-block metric.

    Parameters
    ----------
    name, shape, inner_product_weights
        Inherited from :class:`~orpheus.numerics.space.FunctionSpace`.
        ``name`` is ``"full_field"`` and ``shape`` is the flat
        direct-sum dimension ``(n_interior + n_trace,)``. The composite's own
        ``inner_product_weights`` stays ``None`` — the metric is carried
        per block by :attr:`interior_space` / :attr:`trace_space`, and the
        block-aware :meth:`apply_metric` / :meth:`apply_inverse_metric` /
        :meth:`inner_product` overrides never read the base slot.
    interior_space : FunctionSpace
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

    Notes
    -----
    A curvilinear carrying mesh's ψ½ ray is **System B** — its own
    2-block instance of THIS class (``name="radial_characteristic"``,
    members the split interior / boundary ray spaces; B.2b DP1: one
    family-blind composite-space class, instances differ in members).
    The transient third ``radial_characteristic_space`` slot the 2.5d
    interim carried here was evicted at Phase B.2d — the coupled DOF
    count is the honest two-system sum, never a padded third block.
    """

    interior_space: Optional[FunctionSpace] = field(
        default=None, repr=False, compare=False,
    )
    trace_space: Optional[FunctionSpace] = field(
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
        interior_space: FunctionSpace,
        trace_space: FunctionSpace,
    ) -> "FullFieldSpace":
        r"""Build the composite from its bulk and trace leaf spaces.

        Derives the flat direct-sum ``shape`` from the leaf shapes:
        ``(prod(bulk) + prod(trace),)`` — ``prod`` on every block so the
        identity dimension is robust to a future multi-axis trace (today
        ``trace_space.shape == (total_size,)``).

        **The name is derived from the members (CS4b S4 — the F4 ruling),
        exactly as ``of_axes`` derives a product's name from its factors:**
        ``full_field#<digest>`` folds each member's ``(name, shape)`` pair,
        so the composite's inherited ``(name, shape)`` identity is
        precisely as content-keyed as its members' — a direct sum of
        content-digest-named mints compares by content (twin carriers
        EQUAL, volumes/quadrature/family differences UNEQUAL), and a
        direct sum of hand-built raw spaces compares by whatever their
        names carry, the family convention everywhere else. This is what
        retired the R2 block-blindness (a bare ``"full_field"`` name made
        any two composites of equal flat dimension identical).

        The retired ``name=`` role parameter is deliberate, not lost:
        role is CLASS identity (G2.3) — System B's ψ½ composite is
        ``RadialCharacteristicField`` the *field class*, instantiating
        this SAME family-blind space class (B.2b DP1); tagging the space
        ``"radial_characteristic"`` was the role-flavoured naming the S2a
        ``_SPACE_NAME`` retirement removed from the leaves. Every
        instance is somebody's "full field" (interior ⊕ boundary of ONE
        system).
        """
        n_interior = int(np.prod(interior_space.shape))
        n_trace = int(np.prod(trace_space.shape))
        payload = "|".join(
            f"{s.name}:{s.shape}" for s in (interior_space, trace_space)
        ).encode()
        digest = hashlib.blake2b(payload, digest_size=8).hexdigest()
        return cls(
            name=f"full_field#{digest}",
            shape=(n_interior + n_trace,),
            interior_space=interior_space,
            trace_space=trace_space,
        )

    # ------------------------------------------------------------------
    # Direct-sum metric dispatch (per-block, on a composite field)
    # ------------------------------------------------------------------

    def _require_blocks(self) -> tuple[FunctionSpace, FunctionSpace]:
        r"""Return the ``(interior_space, trace_space)`` pair, guarding the bare-constructor footgun.

        The ``interior_space`` / ``trace_space`` fields default to ``None`` (the
        ``compare=False`` dataclass-field convention), so a composite built via
        the bare constructor instead of :meth:`from_blocks` would
        ``AttributeError`` deep inside the adjoint path. Fail at the boundary
        with intent instead (parse-don't-validate), and return the narrowed
        non-``None`` pair so callers bind locals that type cleanly.
        """
        if self.interior_space is None or self.trace_space is None:
            raise RuntimeError(
                "FullFieldSpace has no block spaces; build it via "
                "FullFieldSpace.from_blocks(interior_space, trace_space) (or "
                "SNMesh.full_field_space), not the bare dataclass constructor."
            )
        return self.interior_space, self.trace_space

    @staticmethod
    def _rebuild(
        x: CompositeField,
        bulk_values: np.ndarray,
        boundary_values: np.ndarray,
    ) -> CompositeField:
        r"""Return a copy of composite field ``x`` with new block ``values``.

        Rebuilds the frozen leaves (``x.interior`` / ``x.boundary``) via
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
        """
        return x._recombine(
            interior=replace(x.interior, values=bulk_values),
            boundary=replace(x.boundary, values=boundary_values),
        )

    def apply_metric(self, x: CompositeField) -> CompositeField:
        r"""Apply the block-diagonal Hilbert metric :math:`G\odot x`.

        ``x`` is a composite field; the bulk block is weighted by
        ``G_bulk = V·w_n`` and the trace block by ``G_trace = |Ω·n|·w_n``,
        each delegated to the matching leaf space. (For System B's ψ½
        instance the same two-block dispatch scales the interior cells by
        the SPD state metric ``G_sd = V_cell`` and the corner block by its
        trace metric — the family-blind dispatch needs no third arm.)
        """
        interior_space, trace_space = self._require_blocks()
        return self._rebuild(
            x,
            interior_space.apply_metric(x.interior.values),
            trace_space.apply_metric(x.boundary.values),
        )

    def apply_inverse_metric(self, x: CompositeField) -> CompositeField:
        r"""Apply the block-diagonal pseudo-inverse metric :math:`G^{+}\odot x`.

        Plain ``1/G`` on the strictly-positive bulk block; the
        Moore–Penrose masked inverse is needed only on the trace block
        (zero on the tangential null space ``|Ω·n| = 0``). Delegated per
        block.
        """
        interior_space, trace_space = self._require_blocks()
        return self._rebuild(
            x,
            interior_space.apply_inverse_metric(x.interior.values),
            trace_space.apply_inverse_metric(x.boundary.values),
        )

    def inner_product(self, x: CompositeField, y: CompositeField) -> float:
        r"""Return the direct-sum inner product
        :math:`\langle x, y\rangle_G = \langle x_{\rm bulk}, y_{\rm bulk}\rangle_{G_{\rm bulk}}
        + \langle x_{\rm trace}, y_{\rm trace}\rangle_{G_{\rm trace}}`.
        """
        interior_space, trace_space = self._require_blocks()
        return (
            interior_space.inner_product(x.interior.values, y.interior.values)
            + trace_space.inner_product(
                x.boundary.values, y.boundary.values
            )
        )
