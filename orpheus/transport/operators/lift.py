r"""The lift — how a BULK action enters the composite ``bulk ⊕ trace`` carrier.

Every volumetric operator of the transport algebra (the gains :math:`S`,
:math:`N_{2n}`, :math:`F`; the multiplier :math:`M[f]`; a diffusion energy
binding lifted onto its scalar composite) acts on the bulk block alone and
emits **nothing** on the trace: *extension by zero*. Until CS4c step 5 that
one fact was spelled **nine times** across four modules (`[M]` 2026-09-04,
the step-5 census: ``transfer.py`` 2, ``fission.py`` 2,
``isotropic_transfer.py`` 1, ``multiplication_operator.py`` 4), once per
operator per carrier family, each spelling naming the trace's zero class by
hand. This module is the ONE home (``coding-elegance`` Pattern 2), and the
ruling that made it possible is R-2 of the step-5 design round
(``.claude/plans/cs4c_binding_design.md`` §19.5): the zero-trace emission
is **blessed** as the bulk operators' composite action — not a transitional
shim — and the carriers declare their role partners
(:class:`~orpheus.transport.fields._bases.RolePair`), so the verb names the
output ROLE and the operand's own class supplies the leaf.

Three verbs and one operator:

* :func:`admit_composite` / :func:`admit_array` — the family's ONLY two
  carrier parses. *Each binding acts through the body its ends select*
  (the step-5 outcome): a composite-bound operator admits exactly the
  :class:`~orpheus.transport.full_field.FullField` whose interior rides the
  bound end's interior space; a plain-bound operator admits exactly the
  bare array of the bound end's shape. Everything else is a typed refusal
  naming the operator — the ``singledispatchmethod`` tables and their
  ``isinstance`` carrier arms (`[M]` 3 tables / 13 arms / 12 carrier
  parses at ``f90f7914``) are retired by these two.
* :func:`lift_bulk_action` — *a bulk action enters the composite by
  extension-by-zero on the trace*: run the interior body, emit the zero
  trace of the operand's boundary class in the output role.
* :func:`embed_bulk_assembly` — the assembly-mode twin: a bulk operator's
  sparse emission embedded in the composite flat layout ``[bulk C-ravel |
  trace]`` (index-identity on the bulk block, zero trace rows/columns).
* :class:`BulkLift` — a bulk-born operator (bound on the interior spaces)
  as an operator on the composite: ``apply``/``apply_transpose`` through
  the verb, ``assemble`` through the embed. This is what a scalar-composite
  consumer (diffusion) builds over its plain-bound energy bindings, so the
  energy bindings themselves stay plain-bound (R-4: the array carriers are
  theirs) and the composite is the lift's — one composite action for the
  whole family instead of a ``FullField`` arm on each energy operator.

The angular gains' lift (the frame's :math:`\ell = 0` conjugation plus the
:math:`\ell \ge 1` redistribution, selected at construction from which end
of the retained analysis face the domain's interior is) is
:class:`~orpheus.transport.operators.angular_lift.AngularLift`, built on
these verbs.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal, cast

import numpy as np

from orpheus.numerics.field import Field
from orpheus.numerics.operator import (
    BlockRole,
    LinearOperator,
    MissingAdjoint,
    MissingAssembly,
    adjointable,
    assemblable,
)
from orpheus.numerics.spaces.full_field_space import FullFieldSpace
from orpheus.transport.fields._bases import BoundaryField, BulkField, FieldRole
from orpheus.transport.full_field import FullField
from orpheus.transport.operators.bound_operator import BoundOperator

if TYPE_CHECKING:
    from orpheus.numerics.assembled_operator import SparseAssembledOperator
    from orpheus.numerics.space import FunctionSpace

__all__ = [
    "BulkLift",
    "admit_array",
    "admit_composite",
    "embed_bulk_assembly",
    "interior_space_of",
    "lift_bulk_action",
]

#: Which of a binding's two ends an admission reads.
End = Literal["domain", "codomain"]


def interior_space_of(space: "FunctionSpace", *, owner: str) -> "FunctionSpace":
    r"""The bulk (interior) block of a composite end — the ONE composite-space
    parse of the operator tier.

    A composite binding's end is a block-bearing
    :class:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace`; its
    interior sizes every bulk action. A plain space carries no interior
    and refuses loudly, naming the owner.
    """
    if not isinstance(space, FullFieldSpace) or space.interior_space is None:
        raise TypeError(
            f"{owner} binds the composite FullFieldSpace (bulk ⊕ trace) — "
            f"this end carries no interior to size the bulk action; got "
            f"{type(space).__name__}."
        )
    return space.interior_space


def admit_composite(
    op: "LinearOperator", x: object, *, end: End = "domain",
) -> FullField:
    r"""Admit ``x`` as the composite element of ``op``'s bound ``end``.

    The carrier is the :class:`~orpheus.transport.full_field.FullField`
    whose interior rides the end's interior space — content equality
    (``==``, the ``(name, shape)`` identity every carrier-minted space
    carries), never object identity. A bare array, a typed bulk field
    outside its composite, or a composite on ANOTHER interior (the moment
    iterate handed to an angular-bound gain — the D16 non-endomorphism the
    step-5 census measured at 143 feeds per windowed solve) is refused with
    the operator's name and both spaces, so the caller can bind the
    operator on the end it meant (``on_moment_domain()`` for the moment
    iterate).
    """
    owner = type(op).__name__
    space = cast("FunctionSpace", getattr(op, end))
    interior = interior_space_of(space, owner=owner)
    if not isinstance(x, FullField):
        raise TypeError(
            f"{owner}: this binding acts on the composite FullField "
            f"(bulk ⊕ trace) of its bound {end}; got {type(x).__name__}. A "
            f"typed bulk field rides inside a FullField; a bare array is "
            f"the PLAIN binding's carrier (the ends select the carrier)."
        )
    if x.interior.space != interior:
        raise TypeError(
            f"{owner}: the composite's interior rides {x.interior.space!r} "
            f"(a {type(x.interior).__name__}) but this binding's {end} "
            f"interior is {interior!r} — every binding acts through the "
            f"body its ends select; bind the operator on the space the "
            f"iterate lives on (on_moment_domain() for a moment iterate)."
        )
    return x


def admit_array(
    op: "LinearOperator", x: object, *, end: End = "domain",
) -> np.ndarray:
    r"""Admit ``x`` as the bare array element of ``op``'s plain bound ``end``.

    A plain binding (an energy binding on the scalar space, a multiplier
    on a bulk space) acts on bare arrays of its end's shape — the
    model-portable contract every solver family feeds. A composite is
    refused naming :class:`BulkLift` (the composite action is the lift's);
    a typed field is refused naming ``.values`` (the ends select the
    carrier: a typed field is the composite binding's operand, not the
    plain binding's); a wrong shape is refused with both shapes.
    """
    owner = type(op).__name__
    space = cast("FunctionSpace", getattr(op, end))
    if isinstance(x, FullField):
        raise TypeError(
            f"{owner}: bound on the plain space {space!r}, fed the composite "
            f"FullField — lift the binding onto the composite with "
            f"BulkLift({owner}(...), domain=composite, codomain=composite)."
        )
    if isinstance(x, Field):
        raise TypeError(
            f"{owner}: the plain binding acts on bare arrays — pass "
            f"{type(x).__name__}.values (the ends select the carrier: a "
            f"typed field is the COMPOSITE binding's operand)."
        )
    if not isinstance(x, np.ndarray):
        raise TypeError(
            f"{owner}: the plain binding acts on bare numpy arrays of shape "
            f"{tuple(space.shape)}; got {type(x).__name__}."
        )
    if tuple(x.shape) != tuple(space.shape):
        raise TypeError(
            f"{owner}: the plain binding's {end} is {space!r} (shape "
            f"{tuple(space.shape)}); got an array of shape {tuple(x.shape)}."
        )
    return x


def lift_bulk_action(
    x: FullField,
    interior_action: Callable[[BulkField], Field],
    *,
    trace_role: FieldRole,
) -> FullField:
    r"""A bulk action enters the composite by extension-by-zero on the trace.

    ``interior_action`` maps the operand's bulk block to the output bulk
    block (the operator's own body, typed on its own terms); the output
    trace is the ZERO field of the operand's boundary class in
    ``trace_role`` — the same class for an emission (``SOURCE_SINK``) as
    the nine retired hand spellings named, and the flux class for an
    inverse (``FLUX``). The zero rides the operand's trace space (CS4b
    S4: output blocks ride the operand's blocks), allocated through the
    ONE :meth:`~orpheus.numerics.field.Field.zeros` primitive.
    """
    trace = x.boundary
    partner = cast("type[BoundaryField]", type(trace).role_partner(trace_role))
    return FullField(
        interior=cast("BulkField", interior_action(x.interior)),
        boundary=partner.zeros(trace.space),
    )


def embed_bulk_assembly(
    bulk: "SparseAssembledOperator",
    *,
    domain: "FunctionSpace",
    codomain: "FunctionSpace",
) -> "SparseAssembledOperator":
    r"""Embed a bulk operator's sparse emission in the composite flat layout.

    The composite flat vector is ``[bulk C-ravel | trace]``
    (:meth:`~orpheus.transport.full_field.Composite.to_flat`), so a bulk
    operator's ``n_bulk × n_bulk`` emission is index-identity on the
    leading block and zero on every trace row and column — the assembly
    twin of :func:`lift_bulk_action`. The entries are carried verbatim
    (no arithmetic), so the embedded action reproduces the bulk action
    bit-for-bit.
    """
    from scipy import sparse

    from orpheus.numerics.assembled_operator import SparseAssembledOperator

    n_rows_bulk = int(np.prod(interior_space_of(codomain, owner="embed_bulk_assembly").shape))
    n_cols_bulk = int(np.prod(interior_space_of(domain, owner="embed_bulk_assembly").shape))
    rows, cols = bulk.shape
    if (rows, cols) != (n_rows_bulk, n_cols_bulk):
        raise ValueError(
            f"embed_bulk_assembly: the bulk emission is {rows}×{cols} but "
            f"the composite interiors are {n_rows_bulk} (codomain) × "
            f"{n_cols_bulk} (domain) — the lifted operator is not bound on "
            f"the composite's interiors."
        )
    coo = bulk.matrix.tocoo()
    matrix = sparse.coo_array(
        (coo.data, (coo.row, coo.col)),
        shape=(int(codomain.shape[0]), int(domain.shape[0])),
    )
    return SparseAssembledOperator(matrix, domain=domain, codomain=codomain)


@dataclass(eq=False)
class BulkLift(BoundOperator["FullField"]):
    r"""A bulk-born operator lifted onto the composite ``bulk ⊕ trace``.

    ``inner`` is bound on the composite's interior spaces (its ``domain``
    / ``codomain`` content-equal the composite ends' interiors — admitted
    at construction, so a lift of an operator bound elsewhere is
    unspellable). The lift acts by :func:`lift_bulk_action` in both
    directions and assembles by :func:`embed_bulk_assembly`; its
    capability axes are the inner's (adjointable iff the inner is,
    assemblable iff the inner is). ``block_role`` is BULK by definition —
    this class IS the bulk-role composite action.

    The diffusion solver's binding (CS4c step 5, R-4): the energy
    bindings :class:`~orpheus.transport.operators.isotropic_transfer.IsotropicScattering`
    + :class:`~orpheus.transport.operators.isotropic_transfer.IsotropicN2N`
    on the mesh's scalar ``bulk_space``, lifted here onto the scalar
    composite ``full_field_space`` so the loss ``L + C − lift(K_iso) − B``
    composes under the ``OperatorSum`` guard and assembles for the exact
    resolvent.

    Parameters
    ----------
    inner : LinearOperator
        The bulk operator — bare-array in, bare-array out on the interior
        spaces.
    domain, codomain : FullFieldSpace
        The composite ends (kw-only, write-once — the base).
    """

    inner: "LinearOperator"

    block_role = BlockRole.BULK

    def __post_init__(self) -> None:
        owner = type(self).__name__
        dom = interior_space_of(self.domain, owner=owner)
        cod = interior_space_of(self.codomain, owner=owner)
        if self.inner.domain != dom or self.inner.codomain != cod:
            raise TypeError(
                f"{owner}: the lifted {type(self.inner).__name__} is bound "
                f"on {self.inner.domain!r} → {self.inner.codomain!r}, but "
                f"the composite's interiors are {dom!r} → {cod!r} — bind "
                f"the inner operator on the composite's interior spaces."
            )

    @property
    def _domain_interior(self) -> "FunctionSpace":
        return interior_space_of(self.domain, owner=type(self).__name__)

    @property
    def _codomain_interior(self) -> "FunctionSpace":
        return interior_space_of(self.codomain, owner=type(self).__name__)

    def apply(self, x: FullField, /) -> FullField:
        r"""``inner`` on the bulk block, the zero source/sink on the trace."""
        psi = admit_composite(self, x, end="domain")
        codomain_interior = self._codomain_interior

        def emit(bulk: BulkField) -> Field:
            return bulk.into_role(
                FieldRole.SOURCE_SINK, self.inner.apply(bulk.values),
                space=codomain_interior,
            )

        return lift_bulk_action(psi, emit, trace_role=FieldRole.SOURCE_SINK)

    def apply_transpose(self, x: FullField, /) -> FullField:
        r"""``innerᵀ`` on the bulk cotangent, the zero source/sink on the trace."""
        chi = admit_composite(self, x, end="codomain")
        domain_interior = self._domain_interior
        inner = self.inner
        if not adjointable(inner):
            raise MissingAdjoint(
                f"BulkLift.apply_transpose requires an adjointable inner "
                f"operator; {type(inner).__name__}.is_adjointable is False."
            )

        def emit(bulk: BulkField) -> Field:
            return bulk.into_role(
                FieldRole.SOURCE_SINK, inner.apply_transpose(bulk.values),
                space=domain_interior,
            )

        return lift_bulk_action(chi, emit, trace_role=FieldRole.SOURCE_SINK)

    @property
    def is_adjointable(self) -> bool:
        return self.inner.is_adjointable

    @property
    def is_assemblable(self) -> bool:
        return self.inner.is_assemblable

    def assemble(self) -> "SparseAssembledOperator":
        r"""The inner's emission in the composite flat layout."""
        inner = self.inner
        if not assemblable(inner):
            raise MissingAssembly(
                f"BulkLift.assemble requires an assemblable inner operator; "
                f"{type(inner).__name__}.is_assemblable is False."
            )
        return embed_bulk_assembly(
            inner.assemble(), domain=self.domain, codomain=self.codomain,
        )

    def __repr__(self) -> str:
        return f"BulkLift({self.inner!r})"
