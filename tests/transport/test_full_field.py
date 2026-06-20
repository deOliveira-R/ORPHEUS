r"""Foundation / L0 tests for the :class:`FullField` base (#257 S4.5 / #217).

The #217 timeless-``FullField`` extraction: the named ``bulk ⊕
boundary`` carrier the operator algebra acts on, with
:class:`~orpheus.transport.timed_full_field.TimedFullField` as the
history-bearing subclass (the cofree comonad
``Cofree(FullField, depth=d)``).

This file pins the DISCRIMINATING checks the extraction calls for.
Where the retired S4 ``TransportState`` Protocol made the check
*structural* (membership of a runtime-checkable Protocol), the
concrete-base split makes it *nominal* (``isinstance`` against the
:class:`FullField` base) plus a timeless-vs-timed distinction:

1. :class:`~orpheus.transport.timed_full_field.TimedFullField` IS a
   :class:`FullField` (by inheritance) AND a
   :class:`~orpheus.numerics.vector.Vector`.
2. ``np.ndarray`` is a :class:`~orpheus.numerics.vector.Vector` (the
   four dunders) but is NOT a :class:`FullField` — the discriminating
   check.
3. A bare bulk leaf
   (:class:`~orpheus.transport.fields.angular_flux.AngularFlux`) is a
   :class:`~orpheus.numerics.vector.Vector` / :class:`Field` but is NOT
   a :class:`FullField` — it is the ``bulk`` member, not the composite.
4. A plain timeless :class:`FullField` IS a :class:`FullField` but is
   NOT a :class:`TimedFullField` (no history) — the timeless/timed
   distinction.
5. The polymorphic-recombine teeth: ``FullField + FullField →
   FullField``; ``TimedFullField + TimedFullField → TimedFullField``
   (type preserved, empty history, ``history_depth`` preserved). This
   pins the :meth:`_recombine` correctness — the load-bearing new
   guarantee of the extraction (the algebra is defined ONCE on the base
   yet returns the right concrete type per subclass).

``foundation`` — a software-invariant type-check (the class hierarchy +
the recombine contract), not a theory-page equation claim.

References
----------

* GH **issue #217** — the timeless-``FullField`` extraction.
* ``.claude/plans/issue_257_coefficient_field_promotion.md`` — S4.5.
* :mod:`orpheus.transport.full_field` — the base (module docstring
  carries the cofree-comonad WHY).
* Grand Report v3 §5.5 (the composite carrier).
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
from orpheus.transport.full_field import FullField
from orpheus.transport.timed_full_field import TimedFullField
from tests.sn._test_helpers import placeholder_materials

pytestmark = pytest.mark.foundation


# vv Mode-8: the canonical ORPHEUS invocation is ``python -O``, which
# strips bare ``assert`` to a NO-OP. These structural gates (isinstance /
# class-hierarchy membership) MUST fire under -O, so they route through
# ``pytest.fail`` (a function call) rather than a bare ``assert`` (a
# false-green tripwire under -O).
def _require(condition: bool, message: str) -> None:
    """A -O-firing assertion (NOT a bare ``assert``)."""
    if not condition:
        pytest.fail(message)


def _is_a(candidate: object, kind: type) -> bool:
    """Runtime ``isinstance`` against a concrete class OR runtime-checkable Protocol.

    The ``candidate: object`` parameter is load-bearing: it is the
    DELIBERATE point of this whole file — a structural check of a
    concrete runtime object against a class / Protocol. Passing
    ``candidate`` through an ``object``-typed boundary keeps the runtime
    semantics IDENTICAL (the ``isinstance`` still fires) while telling
    the type checker "this is an intentional runtime structural probe,
    not a static narrowing" — so pyright does not run its
    structural-overlap analysis on a concrete literal (which would warn
    that ``np.ndarray`` overlaps ``Vector`` "unsafely", which is EXACTLY
    the discriminating fact the test asserts, not a hazard).
    """
    return isinstance(candidate, kind)


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


def _timeless_full_field(m: SNMesh) -> FullField:
    """A plain timeless ``FullField`` (NOT a ``TimedFullField``)."""
    return FullField(
        bulk=AngularFlux.zeros_on(m),
        boundary=BoundaryFlux.zeros_on(m),
    )


# ═══════════════════════════════════════════════════════════════════════
# The discriminating type-check: a FullField is bulk ⊕ boundary.
# ═══════════════════════════════════════════════════════════════════════


class TestFullFieldMembership:
    """``FullField`` is the composite; ndarray / bare leaves are NOT."""

    def test_timed_full_field_is_a_full_field(self):
        """The history-bearing carrier IS a FullField — by inheritance."""
        sn = _slab_mesh()
        tff = TimedFullField.zeros(
            bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn,
        )
        _require(
            _is_a(tff, FullField),
            "TimedFullField must be a FullField (it inherits the timeless "
            "bulk ⊕ boundary composite base).",
        )

    def test_timed_full_field_is_also_a_vector(self):
        """Every full field is a Vector (FullField carries the four dunders)."""
        sn = _slab_mesh()
        tff = TimedFullField.zeros(
            bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn,
        )
        _require(
            _is_a(tff, Vector),
            "TimedFullField must also satisfy Vector — FullField carries "
            "the four vector-space dunders, so a full field IS a vector.",
        )

    def test_timeless_full_field_is_a_full_field_and_a_vector(self):
        """A plain timeless FullField is a FullField and a Vector."""
        sn = _slab_mesh()
        ff = _timeless_full_field(sn)
        _require(_is_a(ff, FullField), "A plain FullField must be a FullField.")
        _require(
            _is_a(ff, Vector),
            "A plain FullField must satisfy Vector (the four dunders).",
        )

    def test_ndarray_is_a_vector_but_not_a_full_field(self):
        """``np.ndarray`` is a Vector (4 dunders) but NOT a FullField.

        The discriminating check: the flat scipy-serialization wire
        format has ``+`` / ``-`` / ``scalar *`` / ``/ scalar`` (so it is
        an honest :class:`Vector`) but it is not the composite ``bulk ⊕
        boundary`` carrier. Passing a bare ndarray as a full field is a
        type error, not a deep AttributeError.
        """
        arr = np.zeros((4, 2, 3))
        _require(
            _is_a(arr, Vector),
            "np.ndarray must satisfy Vector (it has +/-/scalar */ /scalar).",
        )
        _require(
            not _is_a(arr, FullField),
            "np.ndarray must NOT be a FullField — it is the flat wire "
            "format, not the composite bulk ⊕ boundary carrier.",
        )

    def test_bare_bulk_leaf_is_a_vector_but_not_a_full_field(self):
        """A bare ``AngularFlux`` is a Vector / Field but NOT a FullField.

        The bulk leaf is the ``bulk`` MEMBER of a full field, not the
        composite that carries one — it has ``.values`` and the field
        dunders (so it is a :class:`Vector`) but it is not the ``bulk ⊕
        boundary`` composite.
        """
        sn = _slab_mesh()
        flux = AngularFlux.zeros_on(sn)
        _require(
            _is_a(flux, Vector),
            "AngularFlux must satisfy Vector (Field dunders).",
        )
        _require(
            not _is_a(flux, FullField),
            "A bare AngularFlux must NOT be a FullField — it is the bulk "
            "member, not the composite carrier.",
        )


# ═══════════════════════════════════════════════════════════════════════
# The timeless / timed distinction.
# ═══════════════════════════════════════════════════════════════════════


class TestTimelessVsTimed:
    """A timeless FullField has no history; a TimedFullField adds it."""

    def test_timeless_full_field_is_not_a_timed_full_field(self):
        """A plain FullField is NOT a TimedFullField (no history)."""
        sn = _slab_mesh()
        ff = _timeless_full_field(sn)
        _require(
            not _is_a(ff, TimedFullField),
            "A plain timeless FullField must NOT be a TimedFullField — it "
            "carries no history (no advance / at_lag / history_depth).",
        )

    def test_timeless_full_field_has_no_history_machinery(self):
        """The timeless base does not expose the history shift-register verbs."""
        sn = _slab_mesh()
        ff = _timeless_full_field(sn)
        _require(
            not hasattr(ff, "advance"),
            "A timeless FullField must NOT expose advance (a source has no "
            "time).",
        )
        _require(
            not hasattr(ff, "history_depth"),
            "A timeless FullField must NOT carry history_depth.",
        )


# ═══════════════════════════════════════════════════════════════════════
# The polymorphic-recombine teeth (the load-bearing new guarantee).
# ═══════════════════════════════════════════════════════════════════════


class TestPolymorphicRecombine:
    """The inherited dunders return the RIGHT concrete type (the _recombine hook).

    The algebra is defined ONCE on :class:`FullField`, routed through the
    polymorphic :meth:`FullField._recombine` hook. The teeth: a
    :class:`TimedFullField` operand must yield a :class:`TimedFullField`
    (type preserved, empty history, ``history_depth`` preserved) — NOT a
    bare :class:`FullField`. A plain :class:`FullField` operand must
    yield a plain :class:`FullField`.
    """

    def test_full_field_plus_full_field_is_a_full_field(self):
        """FullField + FullField → FullField (the base recombine hook)."""
        sn = _slab_mesh()
        a = _timeless_full_field(sn)
        b = _timeless_full_field(sn)
        # ψ + ψ is forbidden by the #208 affine gate at the leaf; the
        # torsor ψ + (ψ' ⊖ ψ) recombines. Use sub then add-displacement.
        d = b - a  # composite displacement
        out = a + d
        _require(
            type(out) is FullField,
            f"FullField + displacement must return a bare FullField; got "
            f"{type(out).__name__}.",
        )

    def test_full_field_sub_is_a_full_field(self):
        """FullField - FullField → FullField (a composite displacement)."""
        sn = _slab_mesh()
        a = _timeless_full_field(sn)
        b = _timeless_full_field(sn)
        out = a - b
        _require(
            type(out) is FullField,
            f"FullField - FullField must return a bare FullField; got "
            f"{type(out).__name__}.",
        )

    def test_full_field_scalar_mul_is_a_full_field(self):
        """scalar * FullField → FullField."""
        sn = _slab_mesh()
        a = _timeless_full_field(sn)
        out = 2.0 * a
        _require(
            type(out) is FullField,
            f"scalar * FullField must return a bare FullField; got "
            f"{type(out).__name__}.",
        )

    def test_timed_plus_timed_is_a_timed_full_field_empty_history(self):
        """TimedFullField + TimedFullField → TimedFullField (type, empty history).

        The load-bearing recombine tooth: the algebra lives on the base
        but a ``TimedFullField`` operand must round-trip to a
        ``TimedFullField`` (NOT a bare FullField), with empty history and
        preserved ``history_depth`` (#217: algebra results carry empty
        history).
        """
        sn = _slab_mesh()
        a = TimedFullField.zeros(
            bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn, history_depth=5,
        )
        b = TimedFullField.zeros(
            bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn, history_depth=5,
        )
        # Build a non-trivial history on ``a`` so the empty-history claim
        # is a real assertion (not a tautology of zeros()).
        nb = AngularFlux.zeros_on(sn)
        nf = BoundaryFlux.zeros_on(sn)
        a = a.advance(nb, nf)
        _require(a.history_length == 1, "precondition: a must carry history.")

        d = b - a  # composite displacement
        out = a + d
        _require(
            type(out) is TimedFullField,
            f"TimedFullField + displacement must return a TimedFullField, "
            f"NOT a bare FullField; got {type(out).__name__}. This is the "
            f"_recombine correctness tooth.",
        )
        _require(
            out._history == (),
            "Algebra results must carry EMPTY history (#217: history is "
            "iteration metadata, not algebraic state).",
        )
        _require(
            out.history_depth == 5,
            f"Algebra results must preserve history_depth; got "
            f"{out.history_depth} (expected 5).",
        )

    def test_timed_scalar_mul_preserves_type_and_depth(self):
        """scalar * TimedFullField → TimedFullField (type, depth preserved)."""
        sn = _slab_mesh()
        a = TimedFullField.zeros(
            bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn, history_depth=3,
        )
        out = 2.0 * a
        _require(
            type(out) is TimedFullField,
            f"scalar * TimedFullField must return a TimedFullField; got "
            f"{type(out).__name__}.",
        )
        _require(
            out.history_depth == 3 and out._history == (),
            "scalar * TimedFullField must preserve history_depth and carry "
            "empty history.",
        )

    def test_timed_neg_preserves_type(self):
        """-TimedFullField → TimedFullField."""
        sn = _slab_mesh()
        a = TimedFullField.zeros(
            bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn,
        )
        out = -a
        _require(
            type(out) is TimedFullField,
            f"-TimedFullField must return a TimedFullField; got "
            f"{type(out).__name__}.",
        )
