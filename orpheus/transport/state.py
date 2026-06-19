r"""The composite carrier of the transport operator algebra.

The operator algebra in :mod:`orpheus.numerics.operator` reduces every
transport solve to compositions of a small set of linear operators
acting on a discrete state :math:`\psi`:

.. math::

    (L + C - S - F - B)\,\psi \;=\; q .

:mod:`orpheus.numerics.vector` names the *bare* contract that solve
requires — :class:`~orpheus.numerics.vector.Vector`, the four
vector-space operations (``+``, ``-``, ``scalar *``, ``/ scalar``).
That contract is deliberately blind to *what kind* of carrier flows
through: a flat ``np.ndarray`` at the scipy boundary satisfies it, and
so does a single bulk :class:`~orpheus.numerics.field.Field` leaf. Both
are honest :class:`~orpheus.numerics.vector.Vector`\ s, neither is a
transport *state*.

This module names the stronger contract the production transport
operators actually read off the carrier. A :class:`TransportState` is a
:class:`~orpheus.numerics.vector.Vector` (so it carries the four
vector-space ops the iteration drivers compose) **plus** the composite
accessors the #208 block-role dispatch needs — the ``bulk`` field, its
``boundary`` partner, and the ``history_depth`` the algebra threads
through every result:

.. math::

    \psi \;=\; \underbrace{\psi|_{\rm bulk}}_{\text{volumetric}}
        \;\oplus\;
        \underbrace{\psi|_\Gamma}_{\text{boundary trace}}
        \;\oplus\;
        \underbrace{(\psi^{n-1}, \dots)}_{\text{history}} .

The SN operator leaves read ``psi.bulk.values`` (the per-ordinate grid)
and ``psi.bulk.mesh`` (the geometry the codomain field is rebuilt on),
emit a boundary partner, and propagate ``history_depth`` into the
result — so ``bulk`` / ``boundary`` / ``history_depth`` are precisely
the composite structure the algebra speaks, and nothing more.

Why a *structural* (duck-typed) Protocol, not an ABC
----------------------------------------------------

:class:`~orpheus.transport.timed_full_field.TimedFullField` is the one
concrete carrier today, and it satisfies :class:`TransportState`
**without inheritance** — it already exposes ``bulk`` / ``boundary`` /
``history_depth`` and the four vector-space dunders. A nominal ABC
would force ``TimedFullField`` (and every future CP / MoC / diffusion
carrier) to inherit from a transport base purely to claim the contract,
buying nothing the structural Protocol does not already give. Matching
the :class:`~orpheus.numerics.operator.LinearOperator` /
:class:`~orpheus.numerics.vector.Vector` style ("any object providing
these members") keeps the algebra duck-typed end to end.

The discriminating check is the payoff. ``np.ndarray`` IS a
:class:`~orpheus.numerics.vector.Vector` (it has the four dunders) but
is NOT a :class:`TransportState` — it has no ``bulk`` / ``boundary`` /
``history_depth``. A bare bulk leaf such as
:class:`~orpheus.transport.fields.angular_flux.AngularFlux` is likewise
a :class:`~orpheus.numerics.vector.Vector` (and a
:class:`~orpheus.numerics.field.Field`) but NOT a
:class:`TransportState` — it is the ``bulk`` member, not the composite
that carries one. So :class:`TransportState` is *strictly stronger*
than :class:`~orpheus.numerics.vector.Vector`: every transport state is
a vector, but the flat serialization wire format and the individual
field leaves are vectors that are not transport states. Typing an SN
operator leaf ``apply(x: TransportState) -> TransportState`` makes
"pass a bare ``np.ndarray`` as a transport state" a *type error* rather
than an ``AttributeError`` discovered three frames deep at
``x.bulk.values`` (``coding-elegance`` Pattern 4 — illegal states
unrepresentable; Pattern 1 — the signature reads like the domain).

Layering note
-------------

:mod:`orpheus.numerics` sits *below* :mod:`orpheus.transport`, so the
numerics :class:`~orpheus.numerics.vector.Vector` cannot name
``BulkField`` / ``BoundaryField`` (those are transport types). This
contract lives in :mod:`orpheus.transport` precisely because it names
them: :class:`TransportState` is the transport-level *refinement* of
the numerics-level :class:`~orpheus.numerics.vector.Vector`, adding the
composite-accessor structure that only transport can see. The numerics
algebra stays generic over ``Generic[V: Vector]`` (the flat
``np.ndarray`` path and the individual leaves still satisfy it); the
transport-level SN operators narrow ``V`` to :class:`TransportState`,
reading like the domain.
"""

from __future__ import annotations

from typing import Protocol, runtime_checkable

from orpheus.numerics.vector import Vector
from orpheus.transport.fields._bases import BoundaryField, BulkField

__all__ = ["TransportState"]


@runtime_checkable
class TransportState(Vector, Protocol):
    r"""Structural contract for the composite transport carrier.

    Refines :class:`~orpheus.numerics.vector.Vector`: a
    :class:`TransportState` carries the four vector-space operations
    (``+``, ``-``, ``scalar *``, ``/ scalar`` — inherited from
    :class:`~orpheus.numerics.vector.Vector`) **and** the composite
    accessors the production transport operators read — the volumetric
    ``bulk`` field, its ``boundary`` partner, and the ``history_depth``
    that the algebra propagates into every result.

    The contract is intentionally minimal — exactly the members the
    operator algebra reads off the carrier, mirroring
    :class:`~orpheus.numerics.vector.Vector`'s minimalism:

    * ``bulk`` — the volumetric / bulk field
      (:class:`~orpheus.transport.fields._bases.BulkField`). The SN
      leaves read ``psi.bulk.values`` (the per-ordinate grid) and
      ``psi.bulk.mesh`` (the geometry the codomain is rebuilt on).
    * ``boundary`` — the boundary partner field
      (:class:`~orpheus.transport.fields._bases.BoundaryField`), the
      trace the boundary-action operators emit.
    * ``history_depth`` — the iteration / time-derivative history depth
      the algebra threads through composed results.

    Carrier utilities — ``to_flat`` / ``from_flat`` (the scipy
    serialization boundary), ``copy``, ``advance`` / ``at_lag`` (the
    history shift-register) — are deliberately NOT part of the contract;
    they belong to the concrete leaf
    (:class:`~orpheus.transport.timed_full_field.TimedFullField`), not
    to this structural type, exactly as
    :class:`~orpheus.numerics.vector.Vector` excludes ``.copy()`` /
    ``.to_flat()``.

    Satisfied — without inheritance — by
    :class:`~orpheus.transport.timed_full_field.TimedFullField`.
    Deliberately NOT satisfied by ``np.ndarray`` (a
    :class:`~orpheus.numerics.vector.Vector` with no ``bulk`` /
    ``boundary`` / ``history_depth``) nor by a bare bulk leaf such as
    :class:`~orpheus.transport.fields.angular_flux.AngularFlux` (the
    ``bulk`` member, not the composite that carries one) — the
    discriminating check that makes :class:`TransportState` strictly
    stronger than :class:`~orpheus.numerics.vector.Vector`.

    The three members are contracted as **read-only** properties, not
    writable attributes. A writable Protocol attribute is *invariant*
    (the implementing type's attribute must be writable too), which
    would reject the frozen
    :class:`~orpheus.transport.timed_full_field.TimedFullField` whose
    fields are immutable. Read-only is also the *honest* contract: the
    operator algebra only ever READS ``bulk`` / ``boundary`` /
    ``history_depth`` off the carrier (``psi.bulk.values``,
    ``psi.history_depth``) and rebuilds a fresh composite — it never
    assigns them. A read-only property getter is satisfied by both a
    frozen dataclass field and a plain attribute.
    """

    @property
    def bulk(self) -> BulkField:
        """The volumetric / bulk field (the SN leaves read ``.values`` / ``.mesh``)."""
        ...

    @property
    def boundary(self) -> BoundaryField:
        """The boundary partner field (the trace the boundary action emits)."""
        ...

    @property
    def history_depth(self) -> int:
        """The iteration / time-derivative history depth threaded through results."""
        ...
