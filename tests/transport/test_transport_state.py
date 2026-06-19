r"""Foundation / L0 tests for the :class:`TransportState` Protocol (#257 S4).

The §5.5 "honest generic" — the named bulk :math:`\oplus` boundary
:math:`\oplus` history carrier that the #208 block-role dispatch reads,
the transport-level *refinement* of the numerics
:class:`~orpheus.numerics.vector.Vector` (it carries the four
vector-space ops PLUS the ``bulk`` / ``boundary`` / ``history_depth``
composite accessors).

This file pins the DISCRIMINATING type-check the issue calls for — the
intrinsic property that makes :class:`TransportState` *strictly
stronger* than :class:`~orpheus.numerics.vector.Vector`:

1. :class:`~orpheus.transport.timed_full_field.TimedFullField` IS a
   :class:`TransportState` (and a :class:`~orpheus.numerics.vector.Vector`)
   — satisfied structurally, without inheritance.
2. ``np.ndarray`` IS a :class:`~orpheus.numerics.vector.Vector` (the four
   dunders) but is NOT a :class:`TransportState` (no ``bulk`` /
   ``boundary`` / ``history_depth``) — the discriminating check.
3. A bare bulk leaf
   (:class:`~orpheus.transport.fields.angular_flux.AngularFlux`) is a
   :class:`~orpheus.numerics.vector.Vector` / :class:`Field` but is NOT
   a :class:`TransportState` — it is the ``bulk`` member, not the
   composite that carries one.

``foundation`` — a software-invariant type-check (the Protocol
contract), not a theory-page equation claim.

References
----------

* ``.claude/plans/issue_257_coefficient_field_promotion.md`` — S4.
* :mod:`orpheus.transport.state` — the Protocol (module docstring
  carries the WHY).
* Grand Report v3 §5.5 (the composite carrier / honest generic).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.vector import Vector
from orpheus.sn.geometry import SNMesh
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.boundary_flux import BoundaryFlux
from orpheus.transport.state import TransportState
from orpheus.transport.timed_full_field import TimedFullField
from tests.sn._test_helpers import placeholder_materials

pytestmark = pytest.mark.foundation


# vv Mode-8: the canonical ORPHEUS invocation is ``python -O``, which
# strips bare ``assert`` to a NO-OP. These structural gates (isinstance /
# Protocol membership) MUST fire under -O, so they route through
# ``pytest.fail`` (a function call) rather than a bare ``assert`` (a
# false-green tripwire under -O).
def _require(condition: bool, message: str) -> None:
    """A -O-firing assertion (NOT a bare ``assert``)."""
    if not condition:
        pytest.fail(message)


def _is_a(candidate: object, protocol: type) -> bool:
    """Runtime ``isinstance`` against a runtime-checkable Protocol.

    The ``candidate: object`` parameter is load-bearing: it is the
    DELIBERATE point of this whole file — a structural check of a
    concrete runtime object against a Protocol. Passing ``candidate``
    through an ``object``-typed boundary keeps the runtime semantics
    IDENTICAL (the ``isinstance`` still fires) while telling the type
    checker "this is an intentional runtime structural probe, not a
    static narrowing" — so pyright does not run its
    structural-overlap analysis on a concrete literal (which would warn
    that ``np.ndarray`` overlaps ``Vector`` "unsafely", which is
    EXACTLY the discriminating fact the test asserts, not a hazard).
    """
    return isinstance(candidate, protocol)


def _slab_mesh(nx: int = 4, ng: int = 2) -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


# ═══════════════════════════════════════════════════════════════════════
# The discriminating type-check: TransportState ⊋ Vector.
# ═══════════════════════════════════════════════════════════════════════


class TestTransportStateProtocol:
    """``TransportState`` is strictly stronger than ``Vector``."""

    def test_timed_full_field_is_a_transport_state(self):
        """The concrete carrier satisfies the Protocol — without inheritance."""
        sn = _slab_mesh()
        tff = TimedFullField.zeros(
            bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn,
        )
        _require(
            _is_a(tff, TransportState),
            "TimedFullField must satisfy the TransportState Protocol "
            "(it exposes bulk / boundary / history_depth + the four "
            "vector-space dunders), structurally and without inheritance.",
        )

    def test_timed_full_field_is_also_a_vector(self):
        """Every transport state is a vector (TransportState refines Vector)."""
        sn = _slab_mesh()
        tff = TimedFullField.zeros(
            bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn,
        )
        _require(
            _is_a(tff, Vector),
            "TimedFullField must also satisfy Vector — TransportState "
            "refines Vector, so a transport state IS a vector.",
        )

    def test_ndarray_is_a_vector_but_not_a_transport_state(self):
        """``np.ndarray`` is a Vector (4 dunders) but NOT a TransportState.

        The discriminating check: the flat scipy-serialization wire
        format has ``+`` / ``-`` / ``scalar *`` / ``/ scalar`` (so it is
        an honest :class:`Vector`) but carries no ``bulk`` / ``boundary``
        / ``history_depth`` — it is not a composite transport state.
        This pins that ``TransportState`` is *strictly* stronger than
        ``Vector``.
        """
        arr = np.zeros((4, 2, 3))
        _require(
            _is_a(arr, Vector),
            "np.ndarray must satisfy Vector (it has +/-/scalar */ /scalar).",
        )
        _require(
            not _is_a(arr, TransportState),
            "np.ndarray must NOT satisfy TransportState — it has no "
            "bulk / boundary / history_depth. This is the discriminating "
            "check (TransportState ⊋ Vector). Passing a bare ndarray as a "
            "transport state is now a type error, not a deep AttributeError.",
        )

    def test_bare_bulk_leaf_is_a_vector_but_not_a_transport_state(self):
        """A bare ``AngularFlux`` is a Vector / Field but NOT a TransportState.

        The bulk leaf is the ``bulk`` MEMBER of a transport state, not
        the composite that carries one — it has ``.values`` and the
        field dunders (so it is a :class:`Vector`) but no ``.bulk`` /
        ``.boundary`` / ``.history_depth``.
        """
        sn = _slab_mesh()
        flux = AngularFlux.zeros_on(sn)
        _require(
            _is_a(flux, Vector),
            "AngularFlux must satisfy Vector (Field dunders).",
        )
        _require(
            not _is_a(flux, TransportState),
            "A bare AngularFlux must NOT satisfy TransportState — it is "
            "the bulk member, not the composite carrier (no bulk / "
            "boundary / history_depth).",
        )
